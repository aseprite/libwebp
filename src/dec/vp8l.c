// Copyright 2012 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// main entry for the decoder
//
// Author: Vikas Arora (vikaas.arora@gmail.com)
//         jyrki@google.com (Jyrki Alakuijala)

#include <stdlib.h>

#include "./vp8li.h"
#include "../utils/huffman.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#define LOSSLESS_MAGIC_BYTE 0x64

static const size_t kHeaderBytes = 5;
static const uint32_t kImageSizeBits = 14;

static const int kCodeLengthLiterals = 16;
static const int kCodeLengthRepeatCode = 16;
static const int kCodeLengthExtraBits[3] = { 2, 3, 7 };
static const int kCodeLengthRepeatOffsets[3] = { 3, 3, 11 };

#define ARGB_BLACK 0xff000000
#define NUM_LENGTH_CODES    24
#define NUM_DISTANCE_CODES  40
#define NUM_CODES_PER_BYTE 256
// -----------------------------------------------------------------------------
//  Five Huffman codes are used at each meta code:
//  1. green + length prefix codes + color cache codes,
//  2. alpha,
//  3. red,
//  4. blue, and,
//  5. distance prefix codes.
typedef enum {
  GREEN = 0,
  RED   = 1,
  BLUE  = 2,
  ALPHA = 3,
  DIST  = 4,
} HuffIndex;

static const uint16_t kAlphabetSize[HUFFMAN_CODES_PER_META_CODE] = {
  NUM_CODES_PER_BYTE + NUM_LENGTH_CODES,
  NUM_CODES_PER_BYTE, NUM_CODES_PER_BYTE, NUM_CODES_PER_BYTE,
  NUM_DISTANCE_CODES};


#define NUM_CODE_LENGTH_CODES       19
static const uint8_t kCodeLengthCodeOrder[NUM_CODE_LENGTH_CODES] = {
  17, 18, 0, 1, 2, 3, 4, 5, 16, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
};

#define CODE_TO_PLANE_CODES        120
static const uint8_t code_to_plane_lut[CODE_TO_PLANE_CODES] = {
   0x18, 0x07, 0x17, 0x19, 0x28, 0x06, 0x27, 0x29, 0x16, 0x1a,
   0x26, 0x2a, 0x38, 0x05, 0x37, 0x39, 0x15, 0x1b, 0x36, 0x3a,
   0x25, 0x2b, 0x48, 0x04, 0x47, 0x49, 0x14, 0x1c, 0x35, 0x3b,
   0x46, 0x4a, 0x24, 0x2c, 0x58, 0x45, 0x4b, 0x34, 0x3c, 0x03,
   0x57, 0x59, 0x13, 0x1d, 0x56, 0x5a, 0x23, 0x2d, 0x44, 0x4c,
   0x55, 0x5b, 0x33, 0x3d, 0x68, 0x02, 0x67, 0x69, 0x12, 0x1e,
   0x66, 0x6a, 0x22, 0x2e, 0x54, 0x5c, 0x43, 0x4d, 0x65, 0x6b,
   0x32, 0x3e, 0x78, 0x01, 0x77, 0x79, 0x53, 0x5d, 0x11, 0x1f,
   0x64, 0x6c, 0x42, 0x4e, 0x76, 0x7a, 0x21, 0x2f, 0x75, 0x7b,
   0x31, 0x3f, 0x63, 0x6d, 0x52, 0x5e, 0x00, 0x74, 0x7c, 0x41,
   0x4f, 0x10, 0x20, 0x62, 0x6e, 0x30, 0x73, 0x7d, 0x51, 0x5f,
   0x40, 0x72, 0x7e, 0x61, 0x6f, 0x50, 0x71, 0x7f, 0x60, 0x70,
};


//------------------------------------------------------------------------------

static int DecodeImageStream(uint32_t xsize, uint32_t ysize,
                             VP8LDecoder* const dec,
                             argb_t** const decoded_data);

static int ReadImageSize(BitReader* const br,
                         uint32_t* const width, uint32_t* const height) {
  const int signature = VP8LReadBits(br, 8);
  if (signature != LOSSLESS_MAGIC_BYTE) return 0;
  *width = VP8LReadBits(br, kImageSizeBits) + 1;
  *height = VP8LReadBits(br, kImageSizeBits) + 1;
  return 1;
}

int VP8LGetInfo(const uint8_t* data, uint32_t data_size,
                int* width, int* height) {
  if (data_size >= kHeaderBytes) {
    uint32_t w, h;
    BitReader br;
    VP8LInitBitReader(&br, data, kHeaderBytes);
    if (ReadImageSize(&br, &w, &h)) {
      *width = w;
      *height = h;
      return 1;
    } else {
      return 0;  // Could not locate signature.
    }
  } else {
    return 0;         // not enough data
  }
}

static WEBP_INLINE uint32_t SubSampleSize(uint32_t size,
                                          uint32_t sampling_bits) {
  return (size + (1 << sampling_bits) - 1) >> sampling_bits;
}

static WEBP_INLINE uint32_t GetCopyDistance(uint32_t distance_symbol,
                                            BitReader* const br) {
  uint32_t extra_bits, offset;
  if (distance_symbol < 4) {
    return distance_symbol + 1;
  }
  extra_bits = (distance_symbol - 2) >> 1;
  offset = (2 + (distance_symbol & 1)) << extra_bits;
  return offset + VP8LReadBits(br, extra_bits) + 1;
}

static WEBP_INLINE uint32_t GetCopyLength(uint32_t length_symbol,
                                          BitReader* const br) {
  // Length and distance prefixes are encoded the same way.
  return GetCopyDistance(length_symbol, br);
}

static WEBP_INLINE int PlaneCodeToDistance(int xsize, uint32_t plane_code) {
  int dist_code, yoffset, xoffset;
  if (plane_code > CODE_TO_PLANE_CODES) {
    return plane_code - CODE_TO_PLANE_CODES;
  }
  dist_code = code_to_plane_lut[plane_code - 1];
  yoffset = dist_code >> 4;
  xoffset = 8 - (dist_code & 0xf);
  return yoffset * xsize + xoffset;
}

// Decodes the next Huffman code from bit-stream.
// FillBitWindow(br) needs to be called at minimum every second call
// to ReadSymbol.
static WEBP_INLINE int ReadSymbol(const HuffmanTree* tree,
                                  BitReader* const br) {
  const HuffmanTreeNode* node = &tree->nodes_[0];
  while (!HuffmanTreeNodeIsLeaf(node)) {
    node = node->child_[VP8LReadOneBitUnsafe(br)];
  }
  assert(node);
  return node->symbol_;
}

static WEBP_INLINE int ReadSymbolSafe(const HuffmanTree* tree,
                                      BitReader* const br) {
  const HuffmanTreeNode* node = &tree->nodes_[0];
  const int read_safe = (br->pos_ > br->len_ - 8);
  if (read_safe) {
    while (!HuffmanTreeNodeIsLeaf(node)) {
      node = node->child_[VP8LReadOneBit(br)];
    }
  } else {
    while (!HuffmanTreeNodeIsLeaf(node)) {
      node = node->child_[VP8LReadOneBitUnsafe(br)];
    }
  }
  assert(node);

  return node->symbol_;
}

static int ReadHuffmanCodeLengths(
    VP8LDecoder* const dec, const uint32_t* const code_length_code_lengths,
    uint32_t num_codes, uint32_t num_symbols, uint32_t** const code_lengths) {
  int ok = 0;
  BitReader* const br = &dec->br_;
  uint32_t max_length = 0;
  uint32_t code_idx, sym_cnt;
  uint32_t* code_lengths_lcl = *code_lengths;
  int code_len, prev_code_len;
  const int use_length = VP8LReadBits(br, 1);
  HuffmanTree tree;

  if (!HuffmanTreeBuild(code_length_code_lengths, NULL, NULL, num_codes,
                        &tree)) {
    dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
    return 0;
  }

  if (use_length) {
    const int length_nbits = (VP8LReadBits(br, 3) + 1) * 2;
    max_length = VP8LReadBits(br, length_nbits) + 2;
  }
  sym_cnt = 0;
  prev_code_len = 8;
  for (code_idx = 0; code_idx < num_symbols; ++code_idx) {
    if (use_length && ++sym_cnt > max_length) break;
    VP8LFillBitWindow(br);
    code_len = ReadSymbol(&tree, br);
    if (code_len < kCodeLengthLiterals) {
      code_lengths_lcl[code_idx] = code_len;
      if (code_len != 0) prev_code_len = code_len;
    } else {
      const int repeat_ix = code_len - kCodeLengthLiterals;
      const int extra_bits = kCodeLengthExtraBits[repeat_ix];
      const int repeat_offset = kCodeLengthRepeatOffsets[repeat_ix];
      const int rep_cnt = VP8LReadBits(br, extra_bits) + repeat_offset;
      const int use_prev = (code_len == kCodeLengthRepeatCode);
      int rep_iter;
      if (code_idx + rep_cnt > num_symbols) {
        dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
        goto End;
      }
      for (rep_iter = 0; rep_iter < rep_cnt; ++rep_iter) {
        code_lengths_lcl[code_idx + rep_iter] = use_prev ? prev_code_len : 0;
      }
      code_idx += (rep_cnt - 1);
    }
  }
  ok = 1;

 End:
  HuffmanTreeRelease(&tree);
  return ok;
}

static int ReadHuffmanCode(int alphabet_size, VP8LDecoder* const dec,
                           HuffmanTree* const tree) {
  int* symbols = NULL;
  int* codes = NULL;
  uint32_t* code_lengths = NULL;
  int num_symbols;
  int ok = 1;
  BitReader* const br = &dec->br_;
  const int simple_code = VP8LReadBits(br, 1);

  if (simple_code) {
    // Read symbols, codes & code lengths directly.
    int nbits;
    nbits = VP8LReadBits(br, 3);
    if (nbits == 0) {
      num_symbols = 1;
    } else {
      num_symbols = VP8LReadBits(br, 1) + 1;
    }

    symbols = (int*)malloc(num_symbols * sizeof(*symbols));
    codes = (int*)malloc(num_symbols * sizeof(*codes));
    code_lengths = (uint32_t*)malloc(num_symbols * sizeof(*code_lengths));
    if (symbols == NULL || codes == NULL || code_lengths == NULL) {
      free(symbols); free(codes); free(code_lengths);
      dec->status_ = VP8_STATUS_OUT_OF_MEMORY;
      return 0;
    }

    if (nbits == 0) {
      symbols[0] = 0; codes[0] = 0; code_lengths[0] = 0;
    } else {
      const int num_bits = (nbits - 1) * 2 + 4;
      int i;
      for (i = 0; i < num_symbols; ++i) {
        symbols[i] = VP8LReadBits(br, num_bits);
        codes[i] = i;
        code_lengths[i] = num_symbols - 1;
      }
    }
  } else {
    // Decode Huffman-coded code lengths.
    uint32_t code_idx;
    uint32_t code_length_code_lengths[NUM_CODE_LENGTH_CODES] = { 0 };
    const uint32_t num_codes = VP8LReadBits(br, 4) + 4;

    if (num_codes > NUM_CODE_LENGTH_CODES) {
      dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
      return 0;
    }

    code_lengths = (uint32_t*)calloc(alphabet_size, sizeof(code_lengths[0]));
    if (code_lengths == NULL) {
      dec->status_ = VP8_STATUS_OUT_OF_MEMORY;
      return 0;
    }

    for (code_idx = 0; code_idx < num_codes; ++code_idx) {
      code_length_code_lengths[kCodeLengthCodeOrder[code_idx]] =
          VP8LReadBits(br, 3);
    }
    ok = ReadHuffmanCodeLengths(dec, code_length_code_lengths,
                                NUM_CODE_LENGTH_CODES,
                                alphabet_size, &code_lengths);
    num_symbols = alphabet_size;
  }

  // Build Huffman tree.
  ok = ok && HuffmanTreeBuild(code_lengths, codes, symbols, num_symbols, tree);
  free(symbols); free(codes); free(code_lengths);

  if (br->error_) {
    ok = 0;
    dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
  }

  return ok;
}

static int ReadHuffmanCodes(
    int xsize, int ysize, VP8LDecoder* const dec,
    int* const color_cache_bits, int* const color_cache_x_subsample_bits,
    uint32_t** const huffman_image, uint32_t* const huffman_subsample_bits,
    uint32_t** const meta_codes, uint32_t* const meta_code_size,
    HuffmanTree** htrees, uint32_t* const num_huffman_trees) {
  int ok = 0;
  int use_color_cache, color_cache_size;
  uint32_t tree_idx;

  uint32_t* huffman_image_lcl = NULL;
  HuffmanTree* htrees_lcl = NULL;
  uint32_t* meta_codes_lcl = NULL;

  BitReader* const br = &dec->br_;
  const int use_meta_huff_codes = VP8LReadBits(br, 1);

  uint32_t num_huffman_trees_lcl = HUFFMAN_CODES_PER_META_CODE;
  if (use_meta_huff_codes) {
    uint32_t hpc, mc, hc, mci;
    uint32_t meta_codes_nbits, num_meta_codes, nbits;
    uint32_t huffman_xsize, huffman_ysize, huffman_pixs;

    *huffman_subsample_bits = VP8LReadBits(br, 4);
    huffman_xsize = SubSampleSize(xsize, *huffman_subsample_bits);
    huffman_ysize = SubSampleSize(ysize, *huffman_subsample_bits);
    huffman_pixs = huffman_xsize * huffman_ysize;
    if (!DecodeImageStream(huffman_xsize, huffman_ysize, dec,
                           &huffman_image_lcl)) goto Error;
    for (hpc = 0; hpc < huffman_pixs; ++hpc) {
      // The huffman data is stored in red and green bytes.
      huffman_image_lcl[hpc] >>= 8;
      huffman_image_lcl[hpc] &= 0xffff;
    }
    meta_codes_nbits = VP8LReadBits(br, 4);
    num_meta_codes = VP8LReadBits(br, meta_codes_nbits) + 2;
    nbits = VP8LReadBits(br, 4);
    *meta_code_size = num_meta_codes * HUFFMAN_CODES_PER_META_CODE;
    meta_codes_lcl = (uint32_t*)calloc(
        *meta_code_size, sizeof(meta_codes_lcl[0]));
    if (meta_codes_lcl == NULL) {
      dec->status_ = VP8_STATUS_OUT_OF_MEMORY;
      goto Error;
    }

    mci = 0;
    for (mc = 0; mc < num_meta_codes; ++mc) {
      for (hc = 0; hc < HUFFMAN_CODES_PER_META_CODE; ++hc) {
        const uint32_t tree_index = VP8LReadBits(br, nbits);
        meta_codes_lcl[mci] = tree_index;
        if (num_huffman_trees_lcl < tree_index + 1) {
          num_huffman_trees_lcl = tree_index + 1;
        }
        ++mci;
      }
    }
  } else {
    int hc;
    *meta_code_size = HUFFMAN_CODES_PER_META_CODE;
    meta_codes_lcl = (uint32_t*)calloc(
        *meta_code_size, sizeof(meta_codes_lcl[0]));
    if (meta_codes_lcl == NULL) {
      dec->status_ = VP8_STATUS_OUT_OF_MEMORY;
      goto Error;
    }

    for (hc = 0; hc < HUFFMAN_CODES_PER_META_CODE; ++hc) {
      meta_codes_lcl[hc] = hc;
    }
  }

  use_color_cache = VP8LReadBits(br, 1);
  if (use_color_cache) {
    *color_cache_x_subsample_bits = VP8LReadBits(br, 4);
    *color_cache_bits = VP8LReadBits(br, 4);
    color_cache_size = 1 << *color_cache_bits;
  } else {
    *color_cache_x_subsample_bits = 0;
    *color_cache_bits = 0;
    color_cache_size = 0;
  }

  htrees_lcl =
      (HuffmanTree*)calloc(num_huffman_trees_lcl, sizeof(htrees_lcl[0]));
  if (htrees_lcl == NULL) {
    free(meta_codes_lcl);
    dec->status_ = VP8_STATUS_OUT_OF_MEMORY;
    goto Error;
  }
  if (br->error_) goto Error;

  ok = 1;
  for (tree_idx = 0; ok && tree_idx < num_huffman_trees_lcl; ++tree_idx) {
    const int tree_type = tree_idx % HUFFMAN_CODES_PER_META_CODE;
    int alphabet_size = kAlphabetSize[tree_type];
    if (tree_type == 0) {
      alphabet_size += color_cache_size;
    }
    ok = ReadHuffmanCode(alphabet_size, dec, &htrees_lcl[tree_idx]);
    ok = ok & !br->error_;
  }
  if (ok) {
    *huffman_image = huffman_image_lcl;
    *num_huffman_trees = num_huffman_trees_lcl;
    *meta_codes = meta_codes_lcl;
    *htrees = htrees_lcl;
    return ok;
  }

 Error:
  {
    uint32_t i;
    free(huffman_image_lcl);
    free(meta_codes_lcl);
    if (htrees_lcl) {
      for (i = 0; i < num_huffman_trees_lcl; ++i) {
        HuffmanTreeRelease(&htrees_lcl[i]);
      }
      free(htrees_lcl);
    }
    dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
    return 0;
  }
}

static WEBP_INLINE int GetMetaIndex(
    const uint32_t* const image, uint32_t xsize, uint32_t bits, int x, int y) {
  if (bits == 0) return 0;
  return image[xsize * (y >> bits) + (x >> bits)];
}

static WEBP_INLINE void UpdateHtreeForPos(VP8LDecoder* const dec,
                                          int32_t x, uint32_t y) {
  VP8LMetadata* const hdr = &dec->hdr_;
  HuffmanTree* const htrees = hdr->htrees_;
  const uint32_t* const meta_codes = hdr->meta_codes_;
  const int meta_index = HUFFMAN_CODES_PER_META_CODE *
      GetMetaIndex(hdr->huffman_image_, hdr->huffman_xsize_,
                   hdr->huffman_subsample_bits_, x, y);

  if (hdr->meta_index_ != meta_index) {
    hdr->meta_htrees_[GREEN] = &htrees[meta_codes[meta_index + GREEN]];
    hdr->meta_htrees_[RED] = &htrees[meta_codes[meta_index + RED]];
    hdr->meta_htrees_[BLUE] = &htrees[meta_codes[meta_index + BLUE]];
    hdr->meta_htrees_[ALPHA] = &htrees[meta_codes[meta_index + ALPHA]];
    hdr->meta_htrees_[DIST] = &htrees[meta_codes[meta_index + DIST]];
    hdr->meta_index_ = meta_index;
  }
}

int VP8LDecodePixels(VP8LDecoder* const dec, argb_t** const decoded_data) {
  int ok = 1;
  int red, green, blue;
  int alpha = ARGB_BLACK;
  argb_t* data = *decoded_data;
  uint32_t col = 0;
  uint32_t row = dec->row_;
  uint32_t xsize = dec->xsize_;
  uint32_t pix_ix = dec->pix_;
  size_t num_pixs = xsize * dec->ysize_;
  BitReader* const br = &dec->br_;
  VP8LMetadata* const hdr = &dec->hdr_;
  VP8LColorCache* const color_cache = hdr->color_cache_;
  // Values in range [NUM_CODES_PER_BYTE .. color_cache_limit[ are
  // color cache codes.
  const int color_cache_limit =
      NUM_CODES_PER_BYTE + hdr->color_cache_size_;

  while (!br->eos_ && pix_ix < num_pixs) {
    VP8LFillBitWindow(br);

    // Only update the huffman code when moving from one block to the next.
    if ((col & hdr->huffman_mask_) == 0) {
      UpdateHtreeForPos(dec, col, row);
    }

    green = ReadSymbolSafe(hdr->meta_htrees_[GREEN], br);
    if (green < NUM_CODES_PER_BYTE) {
      // Literal.
      // Decode and save this pixel.
      red = ReadSymbolSafe(hdr->meta_htrees_[RED], br);
      VP8LFillBitWindow(br);
      blue = ReadSymbolSafe(hdr->meta_htrees_[BLUE], br);
      alpha = ReadSymbolSafe(hdr->meta_htrees_[ALPHA], br);

      data[pix_ix] = (alpha << 24) + (red << 16) + (green << 8) + blue;
      if (color_cache) VP8LColorCacheInsert(color_cache, col, data[pix_ix]);

      // Update location pointers.
      ++pix_ix; ++col;
      if (col == xsize) {
        ++row; col = 0;
        VP8LCloneBitReader(&dec->br_check_point_, &dec->br_);
        dec->pix_ = pix_ix;
      }
    } else if (green < color_cache_limit) {
      // Color cache.
      // Decode and save this pixel.
      const int color_cache_key = green - NUM_CODES_PER_BYTE;
      argb_t argb;
      ok = VP8LColorCacheLookup(color_cache, col, color_cache_key, &argb);
      if (!ok) goto Error;
      data[pix_ix] = argb;
      VP8LColorCacheInsert(color_cache, col, argb);

      // Update location pointers.
      ++pix_ix; ++col;
      if (col == xsize) {
        ++row; col = 0;
        VP8LCloneBitReader(&dec->br_check_point_, &dec->br_);
        dec->pix_ = pix_ix;
      }
    } else if (green - color_cache_limit < NUM_LENGTH_CODES) {
      // Backward reference
      int dist_symbol;
      uint32_t i, dist_code, dist;
      int clone_br = 0;
      int length_sym = green - color_cache_limit;
      const uint32_t length = GetCopyLength(length_sym, br);
      dist_symbol = ReadSymbolSafe(hdr->meta_htrees_[DIST], br);
      VP8LFillBitWindow(br);
      dist_code = GetCopyDistance(dist_symbol, br);
      dist = PlaneCodeToDistance(xsize, dist_code);
      if ((dist > pix_ix) || (pix_ix + length > num_pixs)) {
        ok = 0;
        goto Error;
      }

      // Fill data for specified (backward-ref) length and update location.
      if (color_cache) {
        for (i = 0; i < length; ++i) {
          data[pix_ix] = data[pix_ix - dist];
          VP8LColorCacheInsert(color_cache, col, data[pix_ix]);
          ++pix_ix; ++col;
          if (col == xsize) {
            ++row; col = 0;
            clone_br = 1;
          }
        }
      } else {
        for (i = 0; i < length; ++i) {
          data[pix_ix] = data[pix_ix - dist];
          ++pix_ix;
        }
        col += length;
        while(col >= xsize) {
          col -= xsize;
          ++row;
          clone_br = 1;
        }
      }
      if (clone_br) {
        VP8LCloneBitReader(&dec->br_check_point_, &dec->br_);
        dec->pix_ = pix_ix;
      }

      if (br->error_) goto Error;

      if (!br->eos_ && pix_ix < num_pixs) {
        UpdateHtreeForPos(dec, col, row);
      }
    } else {
      // Code flow should not come here.
      ok = 0;
      goto Error;
    }
  }
  dec->row_ = row;
  dec->pix_ = pix_ix;

 Error:
  if (br->error_ || !ok) {
    ok = 0;
    dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
  } else if (pix_ix == num_pixs && dec->level_ == 1) {
    dec->state_ = READ_DATA;
  }

  return ok;
}

static WEBP_INLINE uint32_t Average2(uint32_t a0, uint32_t a1) {
  return (((a0 ^ a1) & 0xfefefefeL) >> 1) + (a0 & a1);
}

static WEBP_INLINE uint32_t Average3(uint32_t a0, uint32_t a1, uint32_t a2) {
  return Average2(Average2(a0, a2), a1);
}

static WEBP_INLINE uint32_t Average4(uint32_t a0, uint32_t a1,
                                     uint32_t a2, uint32_t a3) {
  return Average2(Average2(a0, a1), Average2(a2, a3));
}

static WEBP_INLINE uint32_t Add(uint32_t a, uint32_t b) {
  // This computes the sum of each component with mod 256.
  const uint32_t alpha_and_green = (a & 0xff00ff00) + (b & 0xff00ff00);
  const uint32_t red_and_blue = (a & 0x00ff00ff) + (b & 0x00ff00ff);
  return (alpha_and_green & 0xff00ff00) | (red_and_blue & 0x00ff00ff);
}

static WEBP_INLINE uint32_t Clip255(uint32_t a) {
  if (a < NUM_CODES_PER_BYTE) {
    return a;
  }
  // return 0, when a is a negative integer.
  // return 255, when a is positive.
  return ~a >> 24;
}

static WEBP_INLINE int AddSubtractComponentFull(int a, int b, int c) {
  return Clip255(a + b - c);
}

static WEBP_INLINE uint32_t ClampedAddSubtractFull(uint32_t c0, uint32_t c1,
                                                   uint32_t c2) {
  const int a = AddSubtractComponentFull(c0 >> 24, c1 >> 24, c2 >> 24);
  const int r = AddSubtractComponentFull((c0 >> 16) & 0xff,
                                         (c1 >> 16) & 0xff,
                                         (c2 >> 16) & 0xff);
  const int g = AddSubtractComponentFull((c0 >> 8) & 0xff,
                                         (c1 >> 8) & 0xff,
                                         (c2 >> 8) & 0xff);
  const int b = AddSubtractComponentFull(c0 & 0xff, c1 & 0xff, c2 & 0xff);
  return (a << 24) | (r << 16) | (g << 8) | b;
}

static WEBP_INLINE int AddSubtractComponentHalf(int a, int b) {
  return Clip255(a + (a - b) / 2);
}

static WEBP_INLINE uint32_t ClampedAddSubtractHalf(uint32_t c0, uint32_t c1,
                                                   uint32_t c2) {
  const uint32_t ave = Average2(c0, c1);
  const int a = AddSubtractComponentHalf(ave >> 24, c2 >> 24);
  const int r = AddSubtractComponentHalf((ave >> 16) & 0xff, (c2 >> 16) & 0xff);
  const int g = AddSubtractComponentHalf((ave >> 8) & 0xff, (c2 >> 8) & 0xff);
  const int b = AddSubtractComponentHalf((ave >> 0) & 0xff, (c2 >> 0) & 0xff);
  return (a << 24) | (r << 16) | (g << 8) | b;
}

static WEBP_INLINE uint32_t Select(uint32_t a, uint32_t b, uint32_t c) {
  const int p0 = (int)(a >> 24) + (int)(b >> 24) - (int)(c >> 24);
  const int p1 = (int)((a >> 16) & 0xff) + (int)((b >> 16) & 0xff) -
      (int)((c >> 16) & 0xff);
  const int p2 = (int)((a >> 8) & 0xff) + (int)((b >> 8) & 0xff) -
      (int)((c >> 8) & 0xff);
  const int p3 = (int)(a & 0xff) + (int)(b & 0xff) - (int)(c & 0xff);
  const int pa = abs(p0 - (a >> 24)) +
      abs(p1 - ((a >> 16) & 0xff)) +
      abs(p2 - ((a >> 8) & 0xff)) +
      abs(p3 - (a & 0xff));
  const int pb = abs(p0 - (b >> 24)) +
      abs(p1 - ((b >> 16) & 0xff)) +
      abs(p2 - ((b >> 8) & 0xff)) +
      abs(p3 - (b & 0xff));
  if (pa <= pb) {
    return a;
  } else {
    return b;
  }
}

static WEBP_INLINE argb_t PredictValue(uint32_t pred_mode, int xsize,
                                       const argb_t* const argb) {
  switch(pred_mode) {
    case 0: return ARGB_BLACK;
    case 1: return argb[-1];
    case 2: return argb[-xsize];
    case 3: return argb[-xsize + 1];
    case 4: return argb[-xsize - 1];
    case 5: return Average3(argb[-1], argb[-xsize], argb[-xsize + 1]);
    case 6: return Average2(argb[-1], argb[-xsize - 1]);
    case 7: return Average2(argb[-1], argb[-xsize]);
    case 8: return Average2(argb[-xsize - 1], argb[-xsize]);
    case 9: return Average2(argb[-xsize], argb[-xsize + 1]);
    case 10: return Average4(argb[-1], argb[-xsize - 1],
                             argb[-xsize], argb[-xsize + 1]);
    case 11: return Select(argb[-xsize], argb[-1], argb[-xsize - 1]);
    case 12:
      return ClampedAddSubtractFull(argb[-1], argb[-xsize], argb[-xsize - 1]);
    case 13:
      return ClampedAddSubtractHalf(argb[-1], argb[-xsize], argb[-xsize - 1]);
  }
  return 0;
}

static void PredictorInverseTransform(const VP8LTransform* const transform,
                                      argb_t* const decoded_data) {
  size_t row, col, col_start, tile_offset;
  uint32_t pix_ix = 0;
  const uint32_t tile_size = 1 << transform->bits_;
  const uint32_t tiles_per_row = SubSampleSize(transform->xsize_,
                                               transform->bits_);
  // First Row follows the L (mode=1) mode.
  decoded_data[0] = Add(decoded_data[0], ARGB_BLACK);
  for (col = 1; col < transform->xsize_; ++col) {
    decoded_data[col] = Add(decoded_data[col], decoded_data[col - 1]);
  }
  pix_ix += transform->xsize_;

  for (row = 1; row < transform->ysize_; ++row) {
    const uint32_t tile_base_ix = tiles_per_row * (row >> transform->bits_);
    for (tile_offset = 0, col_start = 0; tile_offset < tiles_per_row;
         ++tile_offset, col_start += tile_size) {
      argb_t pred;
      // Pick the appropriate predictor mode (at start of every tile).
      const uint32_t pred_mode =
          (transform->data_[tile_base_ix + tile_offset] >> 8) & 0xff;
      uint32_t col_end = col_start + tile_size;
      if (col_end > transform->xsize_) col_end = transform->xsize_;

      // First col follows the T (mode=2) mode.
      pred = (col_start == 0) ? decoded_data[pix_ix - transform->xsize_] :
          PredictValue(pred_mode, transform->xsize_, decoded_data + pix_ix);
      decoded_data[pix_ix] = Add(decoded_data[pix_ix], pred);
      ++pix_ix;

      // Subsequent columns.
      for (col = col_start + 1; col < col_end; ++col, ++pix_ix) {
        pred = PredictValue(pred_mode, transform->xsize_,
                            decoded_data + pix_ix);
        decoded_data[pix_ix] = Add(decoded_data[pix_ix], pred);
      }
    }
  }
}

static void AddGreenToBlueAndRed(const VP8LTransform* const transform,
                                 argb_t* const decoded_data) {
  size_t i;
  const size_t num_pixs = transform->xsize_ * transform->ysize_;
  for (i = 0; i < num_pixs; ++i) {
    const argb_t argb = decoded_data[i];
    argb_t green = (argb >> 8) & 0xff;
    argb_t red_blue = argb & 0x00ff00ff;
    red_blue += ((green << 16) + green);
    red_blue &= 0x00ff00ff;
    decoded_data[i] = (argb & 0xff00ff00) + red_blue;
  }
}

typedef struct {
  int green_to_red_;
  int green_to_blue_;
  int red_to_blue_;
} ColorTransformElem;

static WEBP_INLINE argb_t ColorTransformDelta(signed char color_pred,
                                              signed char color) {
  return (argb_t)((int)(color_pred) * color) >> 5;
}

static WEBP_INLINE void ColorTransformElemInitFromCode(ColorTransformElem* elem,
                                                       argb_t color_pred) {
  elem->green_to_red_ = color_pred & 0xff;
  elem->green_to_blue_ = (color_pred >> 8) & 0xff;
  elem->red_to_blue_ = (color_pred >> 16) & 0xff;
}

static WEBP_INLINE argb_t TransformColor(const ColorTransformElem* const elem,
                                         argb_t argb) {
  const argb_t green = argb >> 8;
  const argb_t red = argb >> 16;
  argb_t new_red = red;
  argb_t new_blue = argb;

  new_red += ColorTransformDelta(elem->green_to_red_, green);
  new_red &= 0xff;
  new_blue += ColorTransformDelta(elem->green_to_blue_, green);
  new_blue += ColorTransformDelta(elem->red_to_blue_, new_red);
  new_blue &= 0xff;
  return (argb & 0xff00ff00) | (new_red << 16) | (new_blue);
}

static void ColorSpaceInverseTransform(const VP8LTransform* const transform,
                                       argb_t* const decoded_data) {
  size_t row, col, col_start, tile_offset;
  uint32_t pix_ix = 0;
  const uint32_t tile_size = 1 << transform->bits_;
  const uint32_t tiles_per_row = SubSampleSize(transform->xsize_,
                                               transform->bits_);
  for (row = 0; row < transform->ysize_; ++row) {
    const uint32_t tile_base_ix = tiles_per_row * (row >> transform->bits_);
    for (tile_offset = 0, col_start = 0; tile_offset < tiles_per_row;
         ++tile_offset, col_start += tile_size) {
      ColorTransformElem color_pred;
      uint32_t col_end = col_start + tile_size;
      if (col_end > transform->xsize_) col_end = transform->xsize_;
      // Pick the appropriate color predictor mode (at start of every tile).
      ColorTransformElemInitFromCode(
          &color_pred, transform->data_[tile_base_ix + tile_offset]);
      for (col = col_start; col < col_end; ++col, ++pix_ix) {
        decoded_data[pix_ix] = TransformColor(&color_pred,
                                              decoded_data[pix_ix]);
      }
    }
  }
}

static void ColorIndexingInverseTransform(const VP8LTransform* const transform,
                                          argb_t* const decoded_data) {
  size_t i;
  const size_t num_pixs = transform->xsize_ * transform->ysize_;
  for (i = 0; i < num_pixs; ++i) {
    decoded_data[i] = transform->data_[(decoded_data[i] >> 8) & 0xff];
  }
}

static int PixelBundleInverseTransform(const VP8LTransform* const transform,
                                       argb_t** const decoded_data) {
  uint32_t row, col, tile_x;
  const uint32_t bit_depth = 8 >> transform->bits_;
  const uint32_t num_cols = 1 << transform->bits_;
  const uint32_t bit_mask = (1 << bit_depth) - 1;
  const uint32_t xs = SubSampleSize(transform->xsize_, transform->bits_);

  uint32_t* tmp = (uint32_t*)calloc(
      transform->xsize_ * transform->ysize_, sizeof(*tmp));
  if (tmp == NULL) return 0;

  for (row = 0; row < transform->ysize_; ++row) {
    for (tile_x = 0; tile_x < xs; ++tile_x) {
      const argb_t* const rowp = (*decoded_data) + row * xs;
      uint32_t tile_code = (rowp[tile_x] >> 8) & 0xff;
      for (col = 0; col < num_cols; ++col) {
        const uint32_t x_all = tile_x * num_cols + col;
        if (x_all < transform->xsize_) {
          const uint32_t ix = row * transform->xsize_ + x_all;
          const uint32_t green = tile_code & bit_mask;
          tmp[ix] = ARGB_BLACK | (green << 8);
          tile_code >>= bit_depth;
        }
      }
    }
  }
  free(*decoded_data);
  *decoded_data = tmp;

  return 1;
}

static int ApplyInverseTransforms(VP8LDecoder* const dec, int start_idx,
                                  argb_t** const decoded_data) {
  int ok = 1;
  int n = dec->next_transform_;
  assert(start_idx >= 0);
  while(ok && n > start_idx) {
    VP8LTransform* const transform = &(dec->transforms_[--n]);
    switch (transform->type_) {
      case SUBTRACT_GREEN:
        AddGreenToBlueAndRed(transform, *decoded_data);
        break;
      case PREDICTOR_TRANSFORM:
        PredictorInverseTransform(transform, *decoded_data);
        break;
      case CROSS_COLOR_TRANSFORM:
        ColorSpaceInverseTransform(transform, *decoded_data);
        break;
      case COLOR_INDEXING_TRANSFORM:
        ColorIndexingInverseTransform(transform, *decoded_data);
        break;
      case PIXEL_BUNDLE_TRANSFORM:
        ok = PixelBundleInverseTransform(transform, decoded_data);
        if (!ok) {
          dec->status_ = VP8_STATUS_OUT_OF_MEMORY;
        }
        break;
      default:
        dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
        ok = 0;
        break;
    }
    free(transform->data_);
    transform->data_ = NULL;
  }
  dec->next_transform_ = n;

  return ok;
}

static int ReadTransform(int* const xsize, int* const ysize,
                         VP8LDecoder* const dec) {
  int ok = 1;
  BitReader* const br = &dec->br_;
  VP8LTransform* transform = &(dec->transforms_[dec->next_transform_]);
  const VP8LImageTransformType type =
      (VP8LImageTransformType)VP8LReadBits(br, 3);

  transform->type_ = type;
  transform->xsize_ = *xsize;
  transform->ysize_ = *ysize;
  transform->data_ = NULL;

  switch (type) {
    case PREDICTOR_TRANSFORM:
    case CROSS_COLOR_TRANSFORM:
      transform->bits_ = VP8LReadBits(br, 4);
      ok = DecodeImageStream(SubSampleSize(transform->xsize_, transform->bits_),
                             SubSampleSize(transform->ysize_, transform->bits_),
                             dec, &transform->data_);
      break;
    case COLOR_INDEXING_TRANSFORM:
      {
        const int num_colors = VP8LReadBits(br, 8) + 1;
        ok = DecodeImageStream(num_colors, 1, dec, &transform->data_);
        if (ok) {
          int i;
          for (i = 1; i < num_colors; ++i) {
            transform->data_[i] = Add(transform->data_[i] ,
                                      transform->data_[i - 1]);
          }
        }
      }
      break;
    case PIXEL_BUNDLE_TRANSFORM:
      transform->bits_ = VP8LReadBits(br, 2);
      transform->data_ = NULL;
      *xsize = SubSampleSize(transform->xsize_, transform->bits_);
      break;
    case SUBTRACT_GREEN:
      break;
    default:
      dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
      ok = 0;
  }

  if (ok) ++dec->next_transform_;

  return ok;
}

void VP8LInitDecoder(VP8LDecoder* const dec) {
  VP8LMetadata* const hdr = &dec->hdr_;
  dec->xsize_ = 0;
  dec->ysize_ = 0;
  dec->row_ = 0;
  dec->pix_ = 0;
  dec->next_transform_ = 0;
  dec->argb_ = NULL;
  dec->level_ = 0;

  hdr->meta_index_ = -1;
  hdr->color_cache_ = NULL;
  hdr->color_cache_size_ = 0;
  hdr->huffman_image_ = NULL;
  hdr->huffman_subsample_bits_ = 0;
  hdr->meta_codes_ = NULL;
  hdr->htrees_ = NULL;
  hdr->num_huffman_trees_ = 0;
  hdr->huffman_xsize_ = 0;
  hdr->huffman_mask_ = 0;
}

static void VP8LClearMetadata(VP8LDecoder* const dec) {
  VP8LMetadata* hdr;
  if (dec == NULL) return;
  hdr = &dec->hdr_;

  free(hdr->huffman_image_);
  hdr->huffman_image_ = NULL;
  free(hdr->meta_codes_);
  hdr->meta_codes_ = NULL;

  if (hdr->htrees_ != NULL) {
    int i;
    for (i = 0; i < hdr->num_huffman_trees_; ++i) {
      HuffmanTreeRelease(&(hdr->htrees_[i]));
    }
    free(hdr->htrees_);
    hdr->htrees_ = NULL;
    hdr->num_huffman_trees_ = 0;
    hdr->huffman_subsample_bits_ = 0;
    hdr->huffman_xsize_ = 0;
    hdr->huffman_mask_ = 0;
  }

  if (hdr->color_cache_) {
    VP8LColorCacheRelease(hdr->color_cache_);
    free(hdr->color_cache_);
    hdr->color_cache_ = NULL;
    hdr->color_cache_size_ = 0;
  }
  hdr->meta_index_ = -1;
}

static void UpdateDecoder(
    uint32_t xsize, uint32_t ysize,
    VP8LColorCache* const color_cache, int color_cache_size,
    uint32_t* const huffman_image, uint32_t huffman_subsample_bits,
    uint32_t* const meta_codes, HuffmanTree* htrees, int num_huffman_trees_,
    VP8LDecoder* const dec) {
  VP8LMetadata* const hdr = &dec->hdr_;
  dec->xsize_ = xsize;
  dec->ysize_ = ysize;
  dec->row_ = 0;
  dec->pix_ = 0;

  hdr->meta_index_ = -1;
  hdr->color_cache_ = color_cache;
  hdr->color_cache_size_ = color_cache_size;
  hdr->huffman_image_ = huffman_image;
  hdr->huffman_subsample_bits_ = huffman_subsample_bits;
  hdr->meta_codes_ = meta_codes;
  hdr->htrees_ = htrees;
  hdr->num_huffman_trees_ = num_huffman_trees_;
  hdr->huffman_xsize_ = SubSampleSize(xsize, huffman_subsample_bits);
  hdr->huffman_mask_ = (huffman_subsample_bits == 0) ?
      ~0 : (1 << huffman_subsample_bits) - 1;
}

static int DecodeImageStream(uint32_t xsize, uint32_t ysize,
                             VP8LDecoder* const dec,
                             argb_t** const decoded_data) {
  int ok = 1;
  int transform_xsize = xsize;
  int transform_ysize = ysize;

  argb_t* data = NULL;
  uint32_t* huffman_image = NULL;
  HuffmanTree* htrees = NULL;
  uint32_t num_huffman_trees = 0;
  // TODO: htrees & num_huffman_trees should be part of a struct.
  uint32_t huffman_subsample_bits = 0;
  uint32_t* meta_codes = NULL;
  uint32_t num_meta_codes = 0;
  VP8LColorCache* color_cache = NULL;
  int color_cache_bits = 0;
  int color_cache_size = 0;
  int color_cache_x_subsample_bits = 0;

  BitReader* const br = &dec->br_;
  int transform_start_idx = dec->next_transform_;
  ++dec->level_;

  // Step#1: Read the transforms.
  while(ok && VP8LReadBits(br, 1)) {
    ok = ReadTransform(&transform_xsize, &transform_ysize, dec);
  }

  // Step#2: Read the Huffman codes.
  ok = ok && ReadHuffmanCodes(transform_xsize, transform_ysize, dec,
                              &color_cache_bits, &color_cache_x_subsample_bits,
                              &huffman_image, &huffman_subsample_bits,
                              &meta_codes, &num_meta_codes,
                              &htrees, &num_huffman_trees);

  if (!ok) goto End;

  if (dec->level_ == 1) {
    dec->state_ = READ_HDR;
    if (dec->action_ == READ_HDR) goto End;
  }

  if (color_cache_bits > 0) {
    color_cache = (VP8LColorCache*)malloc(sizeof(*color_cache));
    if (color_cache == NULL) {
      dec->status_ = VP8_STATUS_OUT_OF_MEMORY;
      ok = 0;
      goto End;
    }
    color_cache_size = 1 << color_cache_bits;
    VP8LColorCacheInit(color_cache, transform_xsize,
                       color_cache_x_subsample_bits, color_cache_bits);
  }

  data = (argb_t*)calloc(transform_xsize * transform_ysize, sizeof(argb_t));
  if (data == NULL) {
    dec->status_ = VP8_STATUS_OUT_OF_MEMORY;
    ok = 0;
    goto End;
  }

  UpdateDecoder(transform_xsize, transform_ysize,
                color_cache, color_cache_size, huffman_image,
                huffman_subsample_bits, meta_codes, htrees, num_huffman_trees,
                dec);
  // Step#3: Use the Huffman trees to decode the LZ77 encoded data.
  ok = VP8LDecodePixels(dec, &data);

  // Step#4: Apply transforms on the decoded data.
  ok = ok && ApplyInverseTransforms(dec, transform_start_idx, &data);
  ok = ok & !br->error_;

 End:
  UpdateDecoder(transform_xsize, transform_ysize,
                color_cache, color_cache_size, huffman_image,
                huffman_subsample_bits, meta_codes, htrees, num_huffman_trees,
                dec);
  if (dec->level_ > 1) VP8LClearMetadata(dec);

  if (!ok) {
    free(data);
    dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
  } else {
    *decoded_data = data;
  }
  --dec->level_;
  return ok;
}

int VP8LDecodeImage(VP8LDecoder* const dec, VP8Io* const io, uint32_t offset) {
  uint32_t width, height;
  argb_t* decoded_data = NULL;

  if (dec == NULL) return 0;
  if (io == NULL || offset > io->data_size) {
    dec->status_ = VP8_STATUS_INVALID_PARAM;
    return 0;
  }

  dec->status_ = VP8_STATUS_OK;
  VP8LInitBitReader(&dec->br_, io->data + offset, io->data_size - offset);
  if (!ReadImageSize(&dec->br_, &width, &height)) {
    dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
    return 0;
  }
  VP8LInitDecoder(dec);
  dec->state_ = READ_DIM;
  dec->width_ = width;
  dec->height_ = height;
  VP8LCloneBitReader(&dec->br_check_point_, &dec->br_);

  if (DecodeImageStream(width, height, dec, &decoded_data)) {
    dec->argb_ = decoded_data;
    if (dec->state_ == READ_DATA) VP8LClearMetadata(dec);
  } else {
    VP8LClear(dec);
    if (dec->status_ == VP8_STATUS_OK) {
      dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
    }
    return 0;
  }

  return 1;
}

void VP8LClear(VP8LDecoder* const dec) {
  if (dec == NULL) return;
  VP8LClearMetadata(dec);

  free(dec->argb_);
  dec->argb_ = NULL;
  while (dec->next_transform_ > 0) {
    free(dec->transforms_[--dec->next_transform_].data_);
    dec->transforms_[dec->next_transform_].data_ = NULL;
  }
}

//------------------------------------------------------------------------------

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif
