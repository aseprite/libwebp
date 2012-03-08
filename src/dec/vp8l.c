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

#include <stdio.h>
#include <stdlib.h>
#include "./vp8li.h"
#include "../dsp/lossless.h"
#include "../utils/huffman.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#define LOSSLESS_MAGIC_BYTE 0x64

static const int kHeaderBytes = 5;
static const uint32_t kImageSizeBits = 14;

static const int kCodeLengthLiterals = 16;
static const int kCodeLengthRepeatCode = 16;
static const int kCodeLengthExtraBits[3] = { 2, 3, 7 };
static const int kCodeLengthRepeatOffsets[3] = { 3, 3, 11 };

#define NUM_LENGTH_CODES    24
#define NUM_DISTANCE_CODES  40
#define DEFAULT_CODE_LENGTH 8

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
  NUM_LITERAL_CODES + NUM_LENGTH_CODES,
  NUM_LITERAL_CODES, NUM_LITERAL_CODES, NUM_LITERAL_CODES,
  NUM_DISTANCE_CODES
};


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
   0x40, 0x72, 0x7e, 0x61, 0x6f, 0x50, 0x71, 0x7f, 0x60, 0x70
};

static int DecodeImageStream(int xsize, int ysize,
                             VP8LDecoder* const dec,
                             argb_t** const decoded_data);

//------------------------------------------------------------------------------

static int ReadImageSize(BitReader* const br,
                         int* const width, int* const height) {
  const int signature = VP8LReadBits(br, 8);
  if (signature != LOSSLESS_MAGIC_BYTE) return 0;
  *width = VP8LReadBits(br, kImageSizeBits) + 1;
  *height = VP8LReadBits(br, kImageSizeBits) + 1;
  return 1;
}

int VP8LGetInfo(const uint8_t* data, int data_size,
                int* width, int* height) {
  if (data_size < kHeaderBytes) {
    return 0;         // not enough data
  } else {
    int w, h;
    BitReader br;
    VP8LInitBitReader(&br, data, data_size);
    if (!ReadImageSize(&br, &w, &h)) {
      return 0;
    }
    *width = w;
    *height = h;
    return 1;
  }
}

static WEBP_INLINE int GetCopyDistance(int distance_symbol,
                                       BitReader* const br) {
  int extra_bits, offset;
  if (distance_symbol < 4) {
    return distance_symbol + 1;
  }
  extra_bits = (distance_symbol - 2) >> 1;
  offset = (2 + (distance_symbol & 1)) << extra_bits;
  return offset + VP8LReadBits(br, extra_bits) + 1;
}

static WEBP_INLINE int GetCopyLength(int length_symbol,
                                     BitReader* const br) {
  // Length and distance prefixes are encoded the same way.
  return GetCopyDistance(length_symbol, br);
}

static WEBP_INLINE int PlaneCodeToDistance(int xsize, int plane_code) {
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
static int ReadSymbol(const HuffmanTree* tree, BitReader* const br) {
  const HuffmanTreeNode* node = &tree->nodes_[0];
  while (!HuffmanTreeNodeIsLeaf(node)) {
    node = node->child_[VP8LReadOneBitUnsafe(br)];
  }
  assert(node);
  return node->symbol_;
}

static WEBP_INLINE int ReadSymbolSafe(const HuffmanTree* tree,
                                      BitReader* const br) {
  const int read_safe = (br->pos_ > br->len_ - 8);
  if (!read_safe) {
    return ReadSymbol(tree, br);
  } else {
    const HuffmanTreeNode* node = &tree->nodes_[0];
    while (!HuffmanTreeNodeIsLeaf(node)) {
      node = node->child_[VP8LReadOneBit(br)];
    }
    assert(node);
    return node->symbol_;
  }
}

static int ReadHuffmanCodeLengths(
    VP8LDecoder* const dec, const int* const code_length_code_lengths,
    int num_codes, int num_symbols, int* const code_lengths) {
  int ok = 0;
  BitReader* const br = &dec->br_;
  int symbol;
  int max_symbol;
  int prev_code_len = DEFAULT_CODE_LENGTH;
  HuffmanTree tree;

  if (!HuffmanTreeBuild(code_length_code_lengths, NULL, NULL, num_codes,
                        &tree)) {
    dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
    return 0;
  }

  if (VP8LReadBits(br, 1)) {    // use length
    const int length_nbits = 2 + 2 * VP8LReadBits(br, 3);
    max_symbol = 2 + VP8LReadBits(br, length_nbits);
    if (max_symbol > num_symbols) {
      dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
      goto End;
    }
  } else {
    max_symbol = num_symbols;
  }

  symbol = 0;
  while (symbol < num_symbols) {
    int code_len;
    if (max_symbol-- == 0) break;
    VP8LFillBitWindow(br);
    code_len = ReadSymbol(&tree, br);
    if (code_len < kCodeLengthLiterals) {
      code_lengths[symbol++] = code_len;
      if (code_len != 0) prev_code_len = code_len;
    } else {
      const int use_prev = (code_len == kCodeLengthRepeatCode);
      const int slot = code_len - kCodeLengthLiterals;
      const int extra_bits = kCodeLengthExtraBits[slot];
      const int repeat_offset = kCodeLengthRepeatOffsets[slot];
      int repeat = VP8LReadBits(br, extra_bits) + repeat_offset;
      if (symbol + repeat > num_symbols) {
        dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
        goto End;
      } else {
        const int length = use_prev ? prev_code_len : 0;
        while (repeat-- > 0) code_lengths[symbol++] = length;
      }
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
  int* code_lengths = NULL;
  int num_symbols;
  int ok = 0;
  BitReader* const br = &dec->br_;
  const int simple_code = VP8LReadBits(br, 1);

  if (simple_code) {  // Read symbols, codes & code lengths directly.
    const int nbits = VP8LReadBits(br, 3);
    num_symbols = 1 + ((nbits == 0) ? 0 : VP8LReadBits(br, 1));

    symbols = (int*)malloc(num_symbols * sizeof(*symbols));
    codes = (int*)malloc(num_symbols * sizeof(*codes));
    code_lengths = (int*)malloc(num_symbols * sizeof(*code_lengths));
    ok = (symbols != NULL && codes != NULL && code_lengths != NULL);
    if (!ok) {
      dec->status_ = VP8_STATUS_OUT_OF_MEMORY;
      goto End;
    }

    if (nbits == 0) {
      symbols[0] = 0;
      codes[0] = 0;
      code_lengths[0] = 0;
    } else {
      const int num_bits = (nbits - 1) * 2 + 4;
      int i;
      for (i = 0; i < num_symbols; ++i) {
        symbols[i] = VP8LReadBits(br, num_bits);
        codes[i] = i;
        code_lengths[i] = num_symbols - 1;
      }
    }
  } else {  // Decode Huffman-coded code lengths.
    int i;
    int code_length_code_lengths[NUM_CODE_LENGTH_CODES] = { 0 };
    const int num_codes = VP8LReadBits(br, 4) + 4;
    if (num_codes > NUM_CODE_LENGTH_CODES) {
      dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
      goto End;
    }

    code_lengths = (int*)calloc(alphabet_size, sizeof(*code_lengths));
    if (code_lengths == NULL) {
      dec->status_ = VP8_STATUS_OUT_OF_MEMORY;
      goto End;
    }

    for (i = 0; i < num_codes; ++i) {
      code_length_code_lengths[kCodeLengthCodeOrder[i]] =
          VP8LReadBits(br, 3);
    }
    ok = ReadHuffmanCodeLengths(dec, code_length_code_lengths,
                                NUM_CODE_LENGTH_CODES,
                                alphabet_size, code_lengths);
    num_symbols = alphabet_size;
  }

  // Build Huffman tree.
  ok = ok &&
       HuffmanTreeBuild(code_lengths, codes, symbols, num_symbols, tree) &&
       !br->error_;
  if (!ok) {
    dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
    goto End;
  }

 End:
  free(symbols);
  free(codes);
  free(code_lengths);
  return ok;
}

static int ReadHuffmanCodes(
    int xsize, int ysize, VP8LDecoder* const dec,
    int* const color_cache_bits_ptr,
    int* const color_cache_x_subsample_bits_ptr,
    uint32_t** const huffman_image_ptr, int* const huffman_precision_ptr,
    uint32_t** const meta_codes_ptr, int* const meta_codes_size_ptr,
    HuffmanTree** htrees_ptr, int* const num_huffman_trees_ptr) {
  int ok = 0;
  int i;
  BitReader* const br = &dec->br_;

  uint32_t* huffman_image = NULL;
  HuffmanTree* htrees = NULL;
  int num_htrees = HUFFMAN_CODES_PER_META_CODE;
  uint32_t* meta_codes = NULL;
  int meta_codes_size;
  int color_cache_size;

  if (VP8LReadBits(br, 1)) {      // use meta Huffman codes
    int meta_codes_nbits, num_meta_codes, nbits;
    const int huffman_precision = VP8LReadBits(br, 4);
    const int huffman_xsize = VP8LSubSampleSize(xsize, huffman_precision);
    const int huffman_ysize = VP8LSubSampleSize(ysize, huffman_precision);
    const int huffman_pixs = huffman_xsize * huffman_ysize;
    if (!DecodeImageStream(huffman_xsize, huffman_ysize, dec,
                           &huffman_image)) {
      dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
      goto Error;
    }
    *huffman_precision_ptr = huffman_precision;
    for (i = 0; i < huffman_pixs; ++i) {
      // The huffman data is stored in red and green bytes.
      huffman_image[i] = (huffman_image[i] >> 8) & 0xffff;
    }

    meta_codes_nbits = VP8LReadBits(br, 4);
    num_meta_codes = 2 + VP8LReadBits(br, meta_codes_nbits);
    nbits = VP8LReadBits(br, 4);
    meta_codes_size = num_meta_codes * HUFFMAN_CODES_PER_META_CODE;
    meta_codes = (uint32_t*)calloc(meta_codes_size, sizeof(*meta_codes));
    if (meta_codes == NULL) {
      dec->status_ = VP8_STATUS_OUT_OF_MEMORY;
      goto Error;
    }

    for (i = 0; i < meta_codes_size; ++i) {
      const int tree_index = VP8LReadBits(br, nbits);
      meta_codes[i] = tree_index;
      if (num_htrees <= tree_index) {
        num_htrees = tree_index + 1;
      }
    }
  } else {
    meta_codes_size = HUFFMAN_CODES_PER_META_CODE;
    meta_codes = (uint32_t*)malloc(meta_codes_size * sizeof(*meta_codes));
    if (meta_codes == NULL) {
      dec->status_ = VP8_STATUS_OUT_OF_MEMORY;
      goto Error;
    }
    for (i = 0; i < meta_codes_size; ++i) meta_codes[i] = i;
  }

  if (VP8LReadBits(br, 1)) {    // use color cache
    *color_cache_x_subsample_bits_ptr = VP8LReadBits(br, 4);
    *color_cache_bits_ptr = VP8LReadBits(br, 4);
    color_cache_size = 1 << *color_cache_bits_ptr;
  } else {
    *color_cache_x_subsample_bits_ptr = 0;
    *color_cache_bits_ptr = 0;
    color_cache_size = 0;
  }

  htrees = (HuffmanTree*)calloc(num_htrees, sizeof(*htrees));
  if (htrees == NULL) {
    dec->status_ = VP8_STATUS_OUT_OF_MEMORY;
    goto Error;
  }

  ok = !br->error_;
  for (i = 0; ok && i < num_htrees; ++i) {
    const int tree_type = i % HUFFMAN_CODES_PER_META_CODE;
    int alphabet_size = kAlphabetSize[tree_type];
    if (tree_type == 0) {
      alphabet_size += color_cache_size;
    }
    ok = ReadHuffmanCode(alphabet_size, dec, &htrees[i]);
    ok = ok & !br->error_;
  }

  if (ok) {   // finalize pointers and return
    *huffman_image_ptr     = huffman_image;
    *meta_codes_size_ptr   = meta_codes_size;
    *meta_codes_ptr        = meta_codes;
    *num_huffman_trees_ptr = num_htrees;
    *htrees_ptr            = htrees;
    return ok;
  }

 Error:
  free(huffman_image);
  free(meta_codes);
  if (htrees) {
    for (i = 0; i < num_htrees; ++i) {
      HuffmanTreeRelease(&htrees[i]);
    }
    free(htrees);
  }
  return 0;
}

static WEBP_INLINE int GetMetaIndex(
    const uint32_t* const image, int xsize, int bits, int x, int y) {
  if (bits == 0) return 0;
  return image[xsize * (y >> bits) + (x >> bits)];
}

static WEBP_INLINE void UpdateHtreeForPos(VP8LDecoder* const dec,
                                          int x, int y) {
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

static int DecodePixels(VP8LDecoder* const dec, argb_t* const data) {
  int ok = 1;
  int col = 0;
  int pix_ix = 0;
  int row = dec->row_;
  int xsize = dec->xsize_;
  const int num_pixs = xsize * dec->ysize_;
  BitReader* const br = &dec->br_;
  VP8LMetadata* const hdr = &dec->hdr_;
  VP8LColorCache* const color_cache = hdr->color_cache_;
  // Values in range [NUM_CODES_PER_BYTE .. color_cache_limit[ are
  // color cache codes.
  const int color_cache_limit =
      NUM_LITERAL_CODES + hdr->color_cache_size_;

  assert(hdr->htrees_ != NULL);
  assert(hdr->meta_codes_ != NULL);

  while (!br->eos_ && pix_ix < num_pixs) {
    int code;
    VP8LFillBitWindow(br);

    // Only update the huffman code when moving from one block to the next.
    if ((col & hdr->huffman_mask_) == 0) {
      UpdateHtreeForPos(dec, col, row);
    }

    code = ReadSymbolSafe(hdr->meta_htrees_[GREEN], br);
    if (code < NUM_LITERAL_CODES) {
      // Literal.
      // Decode and save this pixel.
      int red, green, blue, alpha;
      red = ReadSymbolSafe(hdr->meta_htrees_[RED], br);
      green = code;
      VP8LFillBitWindow(br);
      blue = ReadSymbolSafe(hdr->meta_htrees_[BLUE], br);
      alpha = ReadSymbolSafe(hdr->meta_htrees_[ALPHA], br);

      data[pix_ix] = (alpha << 24) + (red << 16) + (green << 8) + blue;
      if (color_cache) VP8LColorCacheInsert(color_cache, col, data[pix_ix]);

      // Update location pointers.
      ++pix_ix;
      ++col;
      if (col == xsize) {
        ++row;
        col = 0;
      }
    } else if (code < color_cache_limit) {
      // Color cache.
      // Decode and save this pixel.
      const int color_cache_key = code - NUM_LITERAL_CODES;
      argb_t argb;
      ok = VP8LColorCacheLookup(color_cache, col, color_cache_key, &argb);
      if (!ok) goto Error;
      data[pix_ix] = argb;
      VP8LColorCacheInsert(color_cache, col, argb);

      // Update location pointers.
      ++pix_ix;
      ++col;
      if (col == xsize) {
        ++row;
        col = 0;
      }
    } else if (code - color_cache_limit < NUM_LENGTH_CODES) {
      // Backward reference
      int dist_code, dist;
      int length_sym = code - color_cache_limit;
      const int length = GetCopyLength(length_sym, br);
      const int dist_symbol = ReadSymbolSafe(hdr->meta_htrees_[DIST], br);
      VP8LFillBitWindow(br);
      dist_code = GetCopyDistance(dist_symbol, br);
      dist = PlaneCodeToDistance(xsize, dist_code);
      if ((dist > pix_ix) || (pix_ix + length > num_pixs)) {
        ok = 0;
        goto Error;
      }

      // Fill data for specified (backward-ref) length and update location.
      if (color_cache) {
        int i;
        for (i = 0; i < length; ++i) {
          data[pix_ix] = data[pix_ix - dist];
          VP8LColorCacheInsert(color_cache, col, data[pix_ix]);
          ++pix_ix;
          ++col;
          if (col == xsize) {
            ++row;
            col = 0;
          }
        }
      } else {
        int i;
        for (i = 0; i < length; ++i) {
          data[pix_ix] = data[pix_ix - dist];
          ++pix_ix;
        }
        col += length;
        while (col >= xsize) {
          col -= xsize;
          ++row;
        }
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

 Error:
  if (br->error_ || !ok) {
    ok = 0;
    dec->status_ = (!br->eos_) ?
        VP8_STATUS_BITSTREAM_ERROR : VP8_STATUS_SUSPENDED;
  } else if (pix_ix == num_pixs && dec->level_ == 1) {
    dec->state_ = READ_DATA;
  }

  return ok;
}

static void ClearTransform(VP8LTransform* const transform) {
  free(transform->data_);
  transform->data_ = NULL;
}

static int ApplyInverseTransforms(VP8LDecoder* const dec, int start_idx,
                                  argb_t** const decoded_data) {
  int ok = 1;
  int n = dec->next_transform_;
  assert(start_idx >= 0);
  while (ok && n-- > start_idx) {
    VP8LTransform* const transform = &dec->transforms_[n];
    dec->status_ = VP8LInverseTransform(transform, 0, transform->ysize_,
                                        decoded_data);
    ok = (dec->status_ == VP8_STATUS_OK);
    ClearTransform(transform);
  }
  dec->next_transform_ = start_idx;
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
      ok = DecodeImageStream(VP8LSubSampleSize(transform->xsize_,
                                               transform->bits_),
                             VP8LSubSampleSize(transform->ysize_,
                                               transform->bits_),
                             dec, &transform->data_);
      break;
    case COLOR_INDEXING_TRANSFORM:
      {
        const int num_colors = VP8LReadBits(br, 8) + 1;
        ok = DecodeImageStream(num_colors, 1, dec, &transform->data_);
        if (ok) {
          int i;
          for (i = 1; i < num_colors; ++i) {
            transform->data_[i] = VP8LAddPixels(transform->data_[i] ,
                                                transform->data_[i - 1]);
          }
        }
      }
      break;
    case PIXEL_BUNDLE_TRANSFORM:
      transform->bits_ = VP8LReadBits(br, 2);
      *xsize = VP8LSubSampleSize(transform->xsize_, transform->bits_);
      break;
    case SUBTRACT_GREEN:
      break;
    default:
      dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
      ok = 0;
      break;
  }

  if (ok) ++dec->next_transform_;

  return ok;
}

void VP8LInitDecoder(VP8LDecoder* const dec) {
  VP8LMetadata* const hdr = &dec->hdr_;
  dec->xsize_ = 0;
  dec->ysize_ = 0;
  dec->row_ = 0;
  dec->next_transform_ = 0;
  dec->argb_ = NULL;
  dec->level_ = 0;
  dec->decoded_data_ = NULL;
  dec->output_colorspace_ = MODE_BGRA;

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

  hdr->meta_index_ = -1;
  hdr->color_cache_ = color_cache;
  hdr->color_cache_size_ = color_cache_size;
  hdr->huffman_image_ = huffman_image;
  hdr->huffman_subsample_bits_ = huffman_subsample_bits;
  hdr->meta_codes_ = meta_codes;
  hdr->htrees_ = htrees;
  hdr->num_huffman_trees_ = num_huffman_trees_;
  hdr->huffman_xsize_ = VP8LSubSampleSize(xsize, huffman_subsample_bits);
  hdr->huffman_mask_ = (huffman_subsample_bits == 0) ?
      ~0 : (1 << huffman_subsample_bits) - 1;
}

static int DecodeImageStream(int xsize, int ysize,
                             VP8LDecoder* const dec,
                             argb_t** const decoded_data) {
  int ok = 1;
  int transform_xsize = xsize;
  int transform_ysize = ysize;

  argb_t* data = NULL;
  uint32_t* huffman_image = NULL;
  HuffmanTree* htrees = NULL;
  int num_huffman_trees = 0;
  // TODO: htrees & num_huffman_trees should be part of a struct.
  int huffman_subsample_bits = 0;
  uint32_t* meta_codes = NULL;
  int num_meta_codes = 0;
  VP8LColorCache* color_cache = NULL;
  int color_cache_bits = 0;
  int color_cache_size = 0;
  int color_cache_x_subsample_bits = 0;

  BitReader* const br = &dec->br_;
  int transform_start_idx = dec->next_transform_;
  ++dec->level_;

  // Step#1: Read the transforms.
  while (ok && VP8LReadBits(br, 1)) {
    ok = ReadTransform(&transform_xsize, &transform_ysize, dec);
  }

  // Step#2: Read the Huffman codes.
  ok = ok && ReadHuffmanCodes(transform_xsize, transform_ysize, dec,
                              &color_cache_bits, &color_cache_x_subsample_bits,
                              &huffman_image, &huffman_subsample_bits,
                              &meta_codes, &num_meta_codes,
                              &htrees, &num_huffman_trees);

  if (!ok) {
    dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
    goto End;
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

  UpdateDecoder(transform_xsize, transform_ysize,
                color_cache, color_cache_size, huffman_image,
                huffman_subsample_bits, meta_codes, htrees, num_huffman_trees,
                dec);

  if (dec->level_ == 1) {
    dec->state_ = READ_HDR;
    if (dec->action_ == READ_HDR) goto End;
  }

  data = (argb_t*)calloc(transform_xsize * transform_ysize, sizeof(argb_t));
  if (data == NULL) {
    dec->status_ = VP8_STATUS_OUT_OF_MEMORY;
    ok = 0;
    goto End;
  }

  // Step#3: Use the Huffman trees to decode the LZ77 encoded data.
  ok = DecodePixels(dec, data);

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
    // If not enough data (br.eos_) resulted in BIT_STREAM_ERROR, update the
    // status appropriately.
    if (dec->status_ == VP8_STATUS_BITSTREAM_ERROR && dec->br_.eos_) {
      dec->status_ = VP8_STATUS_SUSPENDED;
    }
  } else {
    *decoded_data = data;
  }
  --dec->level_;
  return ok;
}

int VP8LDecodeHeader(VP8LDecoder* const dec, VP8Io* const io) {
  int width, height;
  argb_t* decoded_data = NULL;

  if (dec == NULL) return 0;
  if (io == NULL || dec->br_offset_ > io->data_size) {
    dec->status_ = VP8_STATUS_INVALID_PARAM;
    return 0;
  }

  dec->status_ = VP8_STATUS_OK;
  VP8LInitBitReader(&dec->br_, io->data + dec->br_offset_,
                    io->data_size - dec->br_offset_);
  if (!ReadImageSize(&dec->br_, &width, &height)) {
    dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
    return 0;
  }
  dec->state_ = READ_DIM;
  dec->width_ = width;
  dec->height_ = height;

  dec->decoded_data_ = NULL;
  dec->action_ = READ_HDR;
  if (!DecodeImageStream(width, height, dec, &decoded_data)) {
    free(decoded_data);
    VP8LClear(dec);
    assert(dec->status_ != VP8_STATUS_OK);
    return 0;
  }
  return 1;
}

int VP8LDecodeImage(VP8LDecoder* const dec) {
  argb_t* data = NULL;
  int xsize, ysize, num_pixels;
  if (dec == NULL) return 0;

  if (dec->next_transform_ > 0) {
    xsize = dec->transforms_[dec->next_transform_ - 1].xsize_;
    ysize = dec->transforms_[dec->next_transform_ - 1].ysize_;
  } else {
    xsize = dec->width_;
    ysize = dec->height_;
  }
  num_pixels = xsize * ysize;

  data = (argb_t*)calloc(num_pixels, sizeof(argb_t));
  if (data == NULL) {
    dec->status_ = VP8_STATUS_OUT_OF_MEMORY;
    goto Err;
  }

  dec->action_ = READ_DATA;
  if (!DecodePixels(dec, data)) {
    goto Err;
  }

  if (!ApplyInverseTransforms(dec, 0, &data)) {
    dec->status_ = VP8_STATUS_BITSTREAM_ERROR;
    goto Err;
  }

  if(!VP8LConvertFromBGRA(data, num_pixels, dec->output_colorspace_,
                          &dec->decoded_data_)) {
    dec->status_ = VP8_STATUS_INVALID_PARAM;
    goto Err;
  }

  dec->argb_ = data;
  VP8LClearMetadata(dec);
  return 1;

 Err:
  free(data);
  VP8LClear(dec);
  assert(dec->status_ != VP8_STATUS_OK);
  return 0;
}

void VP8LClear(VP8LDecoder* const dec) {
  int i;
  if (dec == NULL) return;
  VP8LClearMetadata(dec);

  free(dec->argb_);
  dec->argb_ = NULL;
  free(dec->decoded_data_);
  dec->decoded_data_ = NULL;
  for (i = 0; i < dec->next_transform_; ++i) {
    ClearTransform(&dec->transforms_[i]);
  }
  dec->next_transform_ = 0;
}

//------------------------------------------------------------------------------

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif
