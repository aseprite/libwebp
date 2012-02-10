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
#define HUFFMAN_CODES_PER_META_CODE  5
static const uint16_t kAlphabetSize[HUFFMAN_CODES_PER_META_CODE] = {
  256 + NUM_LENGTH_CODES, 256, 256, 256, NUM_DISTANCE_CODES};


#define NUM_CODE_LENGTH_CODES       19
static const uint8_t kCodeLengthCodeOrder[NUM_CODE_LENGTH_CODES] = {
  17, 18, 0, 1, 2, 3, 4, 5, 16, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
};

//------------------------------------------------------------------------------

static int DecodeImageStream(int xsize, int ysize, VP8LDecoder* const dec,
                             uint32_t** data);

static int ReadImageSize(BitReader* const br,
                         int* width, int* height) {
  const int signature = VP8LReadBits(br, 8);
  if (signature != LOSSLESS_MAGIC_BYTE) return 0;
  *width = VP8LReadBits(br, kImageSizeBits) + 1;
  *height = VP8LReadBits(br, kImageSizeBits) + 1;
  return 1;
}

int VP8LGetInfo(const uint8_t* data, uint32_t data_size,
                int* width, int* height) {
  if (data_size >= kHeaderBytes) {
    int w, h;
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

static int SubSampleSize(int size, int sampling_bits) {
  return (size + (1 << sampling_bits) - 1) >> sampling_bits;
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

static int ReadHuffmanCodes(int xsize, int ysize, VP8LDecoder* const dec,
                            uint32_t** meta_codes, uint32_t* meta_code_size,
                            HuffmanTreeNode** htrees,
                            uint32_t* num_huffman_trees) {
  int ok = 1;
  BitReader* const br = dec->br_;
  int use_palette, palette_x_subsample_bits, palette_code_bits, palette_size;
  uint32_t tree_idx;
  uint32_t* huffman_data;
  HuffmanTreeNode* htrees_lcl = *htrees;
  uint32_t* meta_codes_lcl = *meta_codes;
  const int use_meta_huff_codes = VP8LReadBits(br, 1);

  *num_huffman_trees = HUFFMAN_CODES_PER_META_CODE;
  if (use_meta_huff_codes) {
    int meta_codes_nbits, num_meta_codes, nbits, mc, hc, mci;
    const int huffman_subsample_bits = VP8LReadBits(br, 4);
    DecodeImageStream(SubSampleSize(xsize, huffman_subsample_bits),
                      SubSampleSize(ysize, huffman_subsample_bits),
                      dec, &huffman_data);
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
    palette_size = 1 << palette_code_bits;
  } else {
    palette_size = 0;
  }

  htrees_lcl = (HuffmanTreeNode *)calloc(
      *num_huffman_trees, sizeof(htrees_lcl[0]));
  if (htrees_lcl == NULL) {
    free(meta_codes_lcl);
    meta_codes_lcl = NULL;
    return 0;
  }

  for (tree_idx = 0; ok && tree_idx < *num_huffman_trees; ++tree_idx) {
    const int tree_type = tree_idx % HUFFMAN_CODES_PER_META_CODE;
    int alphabet_size = kAlphabetSize[tree_type];
    if (tree_type == 0) {
      alphabet_size += palette_size;
    }
    ok = ReadHuffmanCode(alphabet_size, dec, &htrees_lcl[tree_idx]);
  }

  return ok;
}

static int DecodeBackwardRefs(VP8LDecoder* const dec,
                              uint32_t* meta_codes, uint32_t num_meta_codes,
                              HuffmanTreeNode* htrees,
                              uint32_t num_huffman_trees) {
  int ok = 1;
  (void)dec;
  (void)meta_codes;
  (void)num_meta_codes;
  (void)htrees;
  (void)num_huffman_trees;

  return ok;
}

static int ApplyTransforms(VP8LDecoder* const dec) {
  int ok = 1;
  (void)dec;

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

static int DecodeImageStream(int xsize, int ysize, VP8LDecoder* const dec,
                             uint32_t** data) {
  int transform_xsize = xsize;
  int transform_ysize = ysize;
  int ok = 1;
  uint32_t tree_idx;
  uint32_t* meta_codes = NULL;
  uint32_t num_meta_codes = 0;
  uint32_t num_huffman_trees = 0;
  HuffmanTreeNode* htrees = NULL;
  BitReader* const br = dec->br_;
  (void)data;

  // Step#1: Read the transforms.
  while(ok && VP8LReadBits(br, 1)) {
    ok = ReadTransform(&transform_xsize, &transform_ysize, dec);
  }

  // Step#2: Read the Huffman codes.
  // TODO: Return Huffman codes (trees) from ReadHuffmanCodes.
  ok = ReadHuffmanCodes(transform_xsize, transform_ysize, dec,
                        &meta_codes, &num_meta_codes,
                        &htrees, &num_huffman_trees);

  // Step#3: Use the Huffman trees to decode the LZ77 encoded data.
  // TODO: Implementation of this method.
  ok = DecodeBackwardRefs(dec, meta_codes, num_meta_codes,
                          htrees, num_huffman_trees);

  // Step#4: Appply transforms on the decoded data.
  // TODO: Implementation of this method.
  ok = ApplyTransforms(dec);

  free(meta_codes);
  for (tree_idx = 0; tree_idx < num_huffman_trees; ++tree_idx) {
    HuffmanTreeRelease(&htrees[tree_idx]);
  }
  free(htrees);

  return ok;
}

int VP8LDecodeImage(VP8LDecoder* const dec, VP8Io* const io, uint32_t offset) {
  int width, height;
  argb_t** data = NULL;
  BitReader br;
  (void)io;
  assert(dec);

  if (offset > io->data_size) return 0;

  VP8LInitBitReader(&br, io->data + offset, io->data_size - offset);
  if (!ReadImageSize(&br, &width, &height)) return 0;
  dec->br_ = &br;
  dec->next_transform_ = 0;
  if (!DecodeImageStream(width, height, dec, data)) return 0;
  dec->argb_ = *data;

  return 1;
}
//------------------------------------------------------------------------------

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif
