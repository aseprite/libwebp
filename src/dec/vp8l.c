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

#include <stdlib.h>

#include "./vp8li.h"
#include "../utils/bit_reader.h"
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

#define NUM_LENGTH_CODES    24
#define NUM_DISTANCE_CODES  40
// -----------------------------------------------------------------------------
//  Five Huffman codes are used at each meta code:
//  1. green + length prefix codes + palette codes,
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
  HUFFMAN_CODES_PER_META_CODE = 5
} HuffIndex;

static const uint16_t kAlphabetSize[HUFFMAN_CODES_PER_META_CODE] = {
  256 + NUM_LENGTH_CODES, 256, 256, 256, NUM_DISTANCE_CODES};


#define NUM_CODE_LENGTH_CODES       19
static const uint8_t kCodeLengthCodeOrder[NUM_CODE_LENGTH_CODES] = {
  17, 18, 0, 1, 2, 3, 4, 5, 16, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
};

#define CODE_TO_PLANE_CODES        120
static const unsigned char code_to_plane_lut[CODE_TO_PLANE_CODES] = {
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
                             VP8LDecoder* const dec, uint32_t** decoded_data);

static int ReadImageSize(BitReader* const br,
                         uint32_t* width, uint32_t* height) {
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

static uint32_t SubSampleSize(uint32_t size, uint32_t sampling_bits) {
  return (size + (1 << sampling_bits) - 1) >> sampling_bits;
}

static inline uint32_t GetCopyDistance(uint32_t distance_symbol,
                                       BitReader* const br) {
  uint32_t extra_bits, offset;
  if (distance_symbol < 4) {
    return distance_symbol + 1;
  }
  extra_bits = (distance_symbol - 2) >> 1;
  offset = (2 + (distance_symbol & 1)) << extra_bits;
  return offset + VP8LReadBits(br, extra_bits) + 1;
}

static inline uint32_t GetCopyLength(uint32_t length_symbol,
                                     BitReader* const br) {
  // Length and distance prefixes are encoded the same way.
  return GetCopyDistance(length_symbol, br);
}

static int PlaneCodeToDistance(int xsize, uint32_t plane_code) {
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
static inline int ReadSymbol(const HuffmanTreeNode* root, BitReader* const br) {
  const HuffmanTreeNode* node = root;
  while (!HuffmanTreeNodeIsLeaf(node)) {
    node = node->child_[VP8LReadOneBitUnsafe(br)];
  }
  assert(node);
  assert(node->symbol_ >= 0);

  return node->symbol_;
}

static int ReadHuffmanCodeLengths(
    VP8LDecoder* const dec, const uint32_t* const code_length_code_lengths,
    uint32_t num_codes, uint32_t num_symbols, uint32_t** const code_lengths) {
  int ok = 1;
  BitReader* const br = dec->br_;
  uint32_t max_length = 0;
  uint32_t code_idx, sym_cnt;
  uint32_t* code_lengths_lcl = *code_lengths;
  int code_len, prev_code_len;
  const int use_length = VP8LReadBits(br, 1);
  HuffmanTreeNode* root = NULL;

  root = HuffmanTreeNodeNew();
  if (root == NULL) return 0;
  if (!HuffmanTreeBuild(root, code_length_code_lengths, num_codes)) return 0;
  if (!HuffmanTreeIsFull(root)) return 0;

  if (use_length) {
    const int length_nbits = (VP8LReadBits(br, 3) + 1) * 2;
    max_length = VP8LReadBits(br, length_nbits) + 2;
  }
  sym_cnt = 0;
  prev_code_len = 8;
  for (code_idx = 0; code_idx < num_symbols; ++code_idx) {
    if (use_length && ++sym_cnt > max_length) break;
    VP8LFillBitWindow(br);
    code_len = ReadSymbol(root, br);
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
      assert(code_idx + rep_cnt <= num_symbols);
      for (rep_iter = 0; rep_iter < rep_cnt; ++rep_iter) {
        code_lengths_lcl[code_idx + rep_iter] = use_prev ? prev_code_len : 0;
      }
      code_idx += (rep_cnt - 1);
    }
  }
  HuffmanTreeRelease(root);

  return ok;
}

static int ReadHuffmanCode(int num_symbols, VP8LDecoder* const dec,
                           HuffmanTreeNode* const root) {
  int ok = 1;
  BitReader* const br = dec->br_;
  uint32_t* code_lengths = NULL;
  const int simple_code = VP8LReadBits(br, 1);
  if (simple_code) {
    int sym_idx;
    int nbits, num_bits;
    num_symbols = VP8LReadBits(br, 1) + 1;
    nbits = VP8LReadBits(br, 3);
    num_bits = (nbits > 0) ? ((nbits - 1) * 2 + 4) : 0;
    code_lengths = (uint32_t*)calloc(num_symbols, sizeof(code_lengths[0]));
    if (code_lengths == NULL) return 0;
    for (sym_idx = 0; ok && sym_idx < num_symbols; ++sym_idx) {
      code_lengths[sym_idx] = VP8LReadBits(br, num_bits);
    }
  } else {
    uint32_t code_idx;
    uint32_t code_length_code_lengths[NUM_CODE_LENGTH_CODES] = { 0 };
    const uint32_t num_codes = VP8LReadBits(br, 4) + 4;

    if (num_codes > NUM_CODE_LENGTH_CODES) return 0;

    code_lengths = (uint32_t*)calloc(num_symbols, sizeof(code_lengths[0]));
    if (code_lengths == NULL) {
      return 0;
    }

    for (code_idx = 0; code_idx < num_codes; ++code_idx) {
      code_length_code_lengths[kCodeLengthCodeOrder[code_idx]] =
          VP8LReadBits(br, 3);
    }
    ok = ReadHuffmanCodeLengths(dec, code_length_code_lengths,
                                NUM_CODE_LENGTH_CODES,
                                num_symbols, &code_lengths);
  }

  HuffmanTreeNodeInit(root);
  ok = HuffmanTreeBuild(root, code_lengths, num_symbols);

  free(code_lengths);

  return ok;
}

static int ReadHuffmanCodes(
    int xsize, int ysize, VP8LDecoder* const dec,
    uint32_t* palette_size,
    uint32_t** huffman_image, uint32_t* huffman_subsample_bits,
    uint32_t** meta_codes, uint32_t* meta_code_size,
    HuffmanTreeNode** htrees, uint32_t* num_huffman_trees) {
  int ok = 1;
  int use_palette, palette_x_subsample_bits, palette_code_bits;
  uint32_t tree_idx;

  uint32_t* huffman_image_lcl = NULL;
  HuffmanTreeNode* htrees_lcl = NULL;
  uint32_t* meta_codes_lcl = NULL;

  BitReader* const br = dec->br_;
  const int use_meta_huff_codes = VP8LReadBits(br, 1);

  *num_huffman_trees = HUFFMAN_CODES_PER_META_CODE;
  if (use_meta_huff_codes) {
    uint32_t hpc, mc, hc, mci;
    uint32_t meta_codes_nbits, num_meta_codes, nbits;
    uint32_t huffman_xsize, huffman_ysize, huffman_pixs;

    *huffman_subsample_bits = VP8LReadBits(br, 4);
    huffman_xsize = SubSampleSize(xsize, *huffman_subsample_bits);
    huffman_ysize = SubSampleSize(ysize, *huffman_subsample_bits);
    huffman_pixs = huffman_xsize * huffman_ysize;
    DecodeImageStream(huffman_xsize, huffman_ysize, dec, &huffman_image_lcl);
    for (hpc = 0; hpc < huffman_pixs; ++hpc) {
      // The huffman data is stored in R & G bytes.
      huffman_image_lcl[hpc] >>= 8;
      huffman_image_lcl[hpc] &= 0xffff;
    }
    meta_codes_nbits = VP8LReadBits(br, 4);
    num_meta_codes = VP8LReadBits(br, meta_codes_nbits) + 2;
    nbits = VP8LReadBits(br, 4);
    *meta_code_size = num_meta_codes * HUFFMAN_CODES_PER_META_CODE;
    meta_codes_lcl = (uint32_t*)calloc(
        *meta_code_size, sizeof(meta_codes_lcl[0]));
    if (meta_codes_lcl == NULL) return 0;

    mci = 0;
    for (mc = 0; mc < num_meta_codes; ++mc) {
      for (hc = 0; hc < HUFFMAN_CODES_PER_META_CODE; ++hc) {
        const uint32_t tree_index = VP8LReadBits(br, nbits);
        meta_codes_lcl[mci] = tree_index;
        if (*num_huffman_trees < tree_index + 1) {
          *num_huffman_trees = tree_index + 1;
        }
        ++mci;
      }
    }
  } else {
    int hc;
    *meta_code_size = HUFFMAN_CODES_PER_META_CODE;
    meta_codes_lcl = (uint32_t*)calloc(
        *meta_code_size, sizeof(meta_codes_lcl[0]));
    if (meta_codes_lcl == NULL) return 0;

    for (hc = 0; hc < HUFFMAN_CODES_PER_META_CODE; ++hc) {
      meta_codes_lcl[hc] = hc;
    }
  }

  use_palette = VP8LReadBits(br, 1);
  if (use_palette) {
    palette_x_subsample_bits = VP8LReadBits(br, 4);
    palette_code_bits = VP8LReadBits(br, 4);
    *palette_size = 1 << palette_code_bits;
  } else {
    *palette_size = 0;
  }

  htrees_lcl = (HuffmanTreeNode *)calloc(
      *num_huffman_trees, sizeof(htrees_lcl[0]));
  if (htrees_lcl == NULL) {
    free(meta_codes_lcl);
    return 0;
  }

  for (tree_idx = 0; ok && tree_idx < *num_huffman_trees; ++tree_idx) {
    const int tree_type = tree_idx % HUFFMAN_CODES_PER_META_CODE;
    int alphabet_size = kAlphabetSize[tree_type];
    if (tree_type == 0) {
      alphabet_size += *palette_size;
    }
    ok = ReadHuffmanCode(alphabet_size, dec, &htrees_lcl[tree_idx]);
  }

  *huffman_image = huffman_image_lcl;
  *meta_codes = meta_codes_lcl;
  *htrees = htrees_lcl;

  return ok;
}

static inline int GetMetaIndex(
    const uint32_t* image, uint32_t xsize, uint32_t bits, int x, int y) {
  if (bits == 0) return 0;
  return image[xsize * (y >> bits) + (x >> bits)];
}

typedef HuffmanTreeNode* HuffmanTreeNodeArray[HUFFMAN_CODES_PER_META_CODE];

static void UpdateHuffmanSet(
    const uint32_t* const huffman_image, const uint32_t* const meta_codes,
    HuffmanTreeNode* const htrees, uint32_t huffman_xsize,
    uint32_t huffman_subsample_bits, uint32_t x, uint32_t y,
    int* orig_meta_ix, HuffmanTreeNodeArray* huffs) {
  const int meta_index = HUFFMAN_CODES_PER_META_CODE *
      GetMetaIndex(huffman_image, huffman_xsize, huffman_subsample_bits, x, y);

  if (*orig_meta_ix != meta_index) {
    HuffmanTreeNode** huffs_lcl = *huffs;
    huffs_lcl[GREEN] = &htrees[meta_codes[meta_index + GREEN]];
    huffs_lcl[RED] = &htrees[meta_codes[meta_index + RED]];
    huffs_lcl[BLUE] = &htrees[meta_codes[meta_index + BLUE]];
    huffs_lcl[ALPHA] = &htrees[meta_codes[meta_index + ALPHA]];
    huffs_lcl[DIST] = &htrees[meta_codes[meta_index + DIST]];
    *orig_meta_ix = meta_index;
  }
}

static int DecodeBackwardRefs(
    VP8LDecoder* const dec,
    uint32_t xsize, uint32_t ysize, uint32_t palette_size,
    const uint32_t* const huffman_image, uint32_t huffman_subsample_bits,
    const uint32_t* const meta_codes, HuffmanTreeNode* htrees,
    uint32_t** decoded_data) {
  int ok = 1;
  int red, green, blue;
  int alpha = 0xff000000;
  int meta_ix = -1;
  uint32_t pos;
  uint32_t x = 0;
  uint32_t y = 0;
  uint32_t* data = *decoded_data;
  BitReader* const br = dec->br_;
  HuffmanTreeNode* huffs[HUFFMAN_CODES_PER_META_CODE] = { 0 };

  // Green values >= 256 but < palette_limit are from the palette.
  const int palette_limit = 256 + palette_size;
  const int huffman_mask = (huffman_subsample_bits == 0) ?
      ~0 : (1 << huffman_subsample_bits) - 1;
  const uint32_t huffman_xsize = SubSampleSize(xsize, huffman_subsample_bits);

  data = (uint32_t*)calloc(xsize * ysize,sizeof(uint32_t));
  if (data == NULL) return 0;

  for (pos = 0; pos < xsize * ysize; ) {
    int length_sym;
    VP8LFillBitWindow(br);

    // Only update the huffman code when moving from one block to the next.
    if ((x & huffman_mask) == 0) {
      UpdateHuffmanSet(huffman_image, meta_codes, htrees, huffman_xsize,
                       huffman_subsample_bits, x, y, &meta_ix, &huffs);
    }

    green = ReadSymbol(huffs[GREEN], br);
    // Literal
    if (green < 256) {
      red = ReadSymbol(huffs[RED], br);
      VP8LFillBitWindow(br);
      blue = ReadSymbol(huffs[BLUE], br);
      alpha = ReadSymbol(huffs[ALPHA], br);

      data[pos] = (alpha << 24) + (red << 16) + (green << 8) + blue;
      ++x;
      if (x >= xsize) {
        x = 0;
        ++y;
      }
      ++pos;
      continue;
    }
    // Backward reference
    length_sym = green - palette_limit;
    if (length_sym < NUM_LENGTH_CODES) {
      int dist_symbol;
      uint32_t i, dist;
      const uint32_t length = GetCopyLength(length_sym, br);
      // Here, we have read the length code prefix + extra bits for the length,
      // so reading the next 15 bits can exhaust the bit window.
      // We must fill the window before the next read.
      VP8LFillBitWindow(br);
      dist_symbol = ReadSymbol(huffs[DIST], br);
      dist = GetCopyDistance(dist_symbol, br);
      dist = PlaneCodeToDistance(xsize, dist);
      assert(dist <= pos);
      assert(pos + length <= xsize * ysize);

      for (i = 0; i < length; ++i) {
        data[pos] = data[pos - dist];
        ++pos;
      }
      x += length;
      while(x >= xsize) {
        x -= xsize;
        ++y;
      }

      if (pos == xsize * ysize) {
        break;
      }

      UpdateHuffmanSet(huffman_image, meta_codes, htrees, huffman_xsize,
                       huffman_subsample_bits, x, y, &meta_ix, &huffs);

      continue;
    }
    assert(0);
  }
  return ok;
}

static int ApplyInverseImageTransform(VP8LTransform* transform) {
  int ok = 1;
  (void)transform;

  return ok;
}

static int ApplyInverseTransforms(VP8LDecoder* const dec, int start_idx) {
  int ok = 1;

  if (dec->next_transform_ == 0) return ok;
  while(ok && dec->next_transform_ > start_idx) {
    VP8LTransform* transform = &(dec->transforms_[dec->next_transform_ - 1]);
    ok = ApplyInverseImageTransform(transform);
    --dec->next_transform_;
  }

  return ok;
}

static int ReadTransform(int* xsize, int* ysize, VP8LDecoder* const dec) {
  int ok = 1;
  BitReader* const br = dec->br_;
  VP8LTransform* transform = &(dec->transforms_[dec->next_transform_]);
  const ImageTransformType type = (ImageTransformType)VP8LReadBits(br, 3);

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
      ok = 0;
  }

  if (ok) ++dec->next_transform_;

  return ok;
}

static int DecodeImageStream(uint32_t xsize, uint32_t ysize,
                             VP8LDecoder* const dec, uint32_t** decoded_data) {
  int transform_xsize = xsize;
  int transform_ysize = ysize;
  int ok = 1;
  uint32_t tree_idx;
  uint32_t* huffman_image = NULL;
  uint32_t* meta_codes = NULL;
  uint32_t huffman_subsample_bits = 0;
  uint32_t num_meta_codes = 0;
  uint32_t num_huffman_trees = 0;
  uint32_t palette_size = 0;

  HuffmanTreeNode* htrees = NULL;
  BitReader* const br = dec->br_;
  int transform_start_idx = dec->next_transform_;

  // Step#1: Read the transforms.
  while(ok && VP8LReadBits(br, 1)) {
    ok = ReadTransform(&transform_xsize, &transform_ysize, dec);
  }

  // Step#2: Read the Huffman codes.
  // TODO: Return Huffman codes (trees) from ReadHuffmanCodes.
  ok = ReadHuffmanCodes(transform_xsize, transform_ysize, dec,
                        &palette_size, &huffman_image, &huffman_subsample_bits,
                        &meta_codes, &num_meta_codes,
                        &htrees, &num_huffman_trees);

  // Step#3: Use the Huffman trees to decode the LZ77 encoded data.
  // TODO: Implementation of this method.
  ok = DecodeBackwardRefs(dec, xsize, ysize, palette_size,
                          huffman_image, huffman_subsample_bits,
                          meta_codes, htrees, decoded_data);

  free(huffman_image);
  free(meta_codes);
  for (tree_idx = 0; tree_idx < num_huffman_trees; ++tree_idx) {
    HuffmanTreeRelease(&htrees[tree_idx]);
  }
  free(htrees);

  // Step#4: Appply transforms on the decoded data.
  ok = ApplyInverseTransforms(dec, transform_start_idx);

  return ok;
}

int VP8LDecodeImage(VP8LDecoder* const dec, VP8Io* const io, uint32_t offset) {
  uint32_t width, height;
  argb_t** decoded_data = NULL;
  BitReader br;
  (void)io;
  assert(dec);

  if (offset > io->data_size) return 0;

  VP8LInitBitReader(&br, io->data + offset, io->data_size - offset);
  if (!ReadImageSize(&br, &width, &height)) return 0;
  dec->br_ = &br;
  dec->next_transform_ = 0;
  if (!DecodeImageStream(width, height, dec, decoded_data)) return 0;
  dec->argb_ = *decoded_data;

  return 1;
}
//------------------------------------------------------------------------------

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif
