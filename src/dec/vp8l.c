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

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#define LOSSLESS_MAGIC_BYTE 0x64

static const size_t kHeaderBytes = 5;
static const uint32_t kImageSizeBits = 14;

// -----------------------------------------------------------------------------
//  Five Huffman codes are used at each meta code:
//  1. green + length prefix codes + palette codes,
//  2. alpha,
//  3. red,
//  4. blue, and,
//  5. distance prefix codes.
//  For these 5 Huffman indices, there are three tree-types.
//  0 is for green + length prefix codes + palette codes,
//  1 is for alpha, red and blue,
//  2 is for distance prefix codes.
#define HUFFMAN_CODES_PER_META_CODE  5
static const uint8_t kTreeType[HUFFMAN_CODES_PER_META_CODE] = { 0, 1, 1, 1, 2 };
static const uint16_t kAlphabetSize[HUFFMAN_CODES_PER_META_CODE] = {
  280, 256, 256, 256, 24 };

static const int kNumLengthSymbols = 24;
static const int kCodeLengthLiterals = 16;
static const int kCodeLengthRepeatCode = 16;
static const int kCodeLengthExtraBits[3] = { 2, 3, 7 };
static const int kCodeLengthRepeatOffsets[3] = { 3, 3, 11 };

#define CODE_LENGTH_CODES           19
static const uint8_t kCodeLengthCodeOrder[CODE_LENGTH_CODES] = {
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

static int BuildHuffmanTree(int* code_lengths, int num_symbols) {
  int ok = 1;
  (void)code_lengths;
  (void)num_symbols;
  // TODO: Integrate Huffman BuildTree from array of code_lengths.

  return ok;
}

static int ReadHuffmanCodeLengths(VP8LDecoder* const dec,
                                  int* code_length_code_lengths, int num_codes,
                                  int num_symbols, int** code_lengths) {
  int ok = 1;
  BitReader* const br = dec->br_;
  int max_length = 0;
  int code_idx, sym_cnt;
  int* code_lengths_lcl = *code_lengths;
  int code_len, prev_code_len;
  const int use_length = VP8LReadBits(br, 1);
  (void)code_length_code_lengths;
  (void)num_codes;

  if (use_length) {
    const int length_nbits = (VP8LReadBits(br, 3) + 1) * 2;
    max_length = VP8LReadBits(br, length_nbits) + 2;
  }
  sym_cnt = 0;
  prev_code_len = 8;
  for (code_idx = 0; code_idx < num_symbols; ++code_idx) {
    if (use_length && ++sym_cnt > max_length) break;
    VP8LFillBitWindow(br);
    // TODO: Construct huffman tree from code_length_code_lengths and read
    // symbol from this tree instead of reading 4 bits.
    code_len = VP8LReadBits(br, 4);
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
      for (rep_iter = 0; rep_iter < rep_cnt; ++rep_iter) {
        code_lengths_lcl[code_idx + rep_iter] = use_prev ? prev_code_len : 0;
      }
      code_idx += (rep_cnt - 1);
    }
  }

  return ok;
}

static int ReadHuffmanCode(int num_symbols, VP8LDecoder* const dec) {
  int ok = 1;
  BitReader* const br = dec->br_;
  int* code_lengths = NULL;
  const int simple_code = VP8LReadBits(br, 1);
  if (simple_code) {
    int sym_idx;
    int nbits, num_bits;
    num_symbols = VP8LReadBits(br, 1) + 1;
    nbits = VP8LReadBits(br, 3);
    num_bits = (nbits > 0) ? ((nbits - 1) * 2 + 4) : 0;
    code_lengths = (int*)malloc(num_symbols * sizeof(code_lengths[0]));
    memset(code_lengths, 0, num_symbols * sizeof(code_lengths[0]));
    if (code_lengths == NULL) return 0;
    for (sym_idx = 0; ok && sym_idx < num_symbols; ++sym_idx) {
      code_lengths[sym_idx] = VP8LReadBits(br, num_bits);
    }
  } else {
    int code_idx;
    const int num_codes = VP8LReadBits(br, 4) + 4;
    int* code_length_code_lengths = (int*)malloc(
        num_codes * sizeof(code_length_code_lengths[0]));
    if (code_length_code_lengths == NULL) return 0;

    code_lengths = (int*)malloc(num_symbols * sizeof(code_lengths[0]));
    memset(code_lengths, 0, num_symbols * sizeof(code_lengths[0]));
    if (code_lengths == NULL) {
      free(code_length_code_lengths);
      return 0;
    }

    for (code_idx = 0; code_idx < num_codes; ++code_idx) {
      code_length_code_lengths[kCodeLengthCodeOrder[code_idx]] =
          VP8LReadBits(br, 3);
    }
    ok = ReadHuffmanCodeLengths(dec, code_length_code_lengths, num_codes,
                                num_symbols, &code_lengths);
    free(code_length_code_lengths);
  }

  // TODO: This function will generate thhe root of HuffmanTree.
  ok = BuildHuffmanTree(code_lengths, num_symbols);

  free(code_lengths);

  return ok;
}

static int ReadHuffmanCodes(int xsize, int ysize, VP8LDecoder* const dec,
                            int** meta_codes, int* meta_code_size) {
  int ok = 1;
  BitReader* const br = dec->br_;
  int use_palette, palette_x_subsample_bits, palette_code_bits, palette_size;
  uint32_t* huffman_data;
  int tree_idx;
  int num_huffman_trees = HUFFMAN_CODES_PER_META_CODE;
  const int use_meta_huff_codes = VP8LReadBits(br, 1);

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
    *meta_codes = (int*)malloc(*meta_code_size * sizeof((*meta_codes)[0]));
    if (*meta_codes == NULL) return 0;

    mci = 0;
    for (mc = 0; mc < num_meta_codes; ++mc) {
      for (hc = 0; hc < HUFFMAN_CODES_PER_META_CODE; ++hc) {
        const int tree_index = VP8LReadBits(br, nbits);
        (*meta_codes)[mci] = tree_index;
        if (num_huffman_trees < tree_index + 1) {
          num_huffman_trees = tree_index + 1;
        }
        ++mci;
      }
    }
  } else {
    int hc;
    *meta_code_size = HUFFMAN_CODES_PER_META_CODE;
    *meta_codes = (int*)malloc(*meta_code_size * sizeof((*meta_codes)[0]));
    if (*meta_codes == NULL) return 0;

    for (hc = 0; hc < HUFFMAN_CODES_PER_META_CODE; ++hc) {
      (*meta_codes)[hc] = hc;
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

  for (tree_idx = 0; ok && tree_idx < num_huffman_trees; ++tree_idx) {
    const int tree_type = tree_idx % HUFFMAN_CODES_PER_META_CODE;
    int alphabet_size = kAlphabetSize[tree_type];
    if (tree_type == 0) {
      alphabet_size += palette_size;
    }
    ok = ReadHuffmanCode(alphabet_size, dec);
  }

  return ok;
}

static int DecodeBackwardRefs(int* meta_codes, int num_meta_codes,
                              VP8LDecoder* const dec) {
  int ok = 1;
  (void)meta_codes;
  (void)num_meta_codes;
  (void)dec;

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
  int* meta_codes = NULL;
  int num_meta_codes = 0;
  BitReader* const br = dec->br_;
  (void)data;

  // Step#1: Read the transforms.
  while(ok && VP8LReadBits(br, 1)) {
    ok = ReadTransform(&transform_xsize, &transform_ysize, dec);
  }

  // Step#2: Read the Huffman codes.
  // TODO: Return Huffman codes (trees) from ReadHuffmanCodes.
  ok = ReadHuffmanCodes(transform_xsize, transform_ysize, dec,
                        &meta_codes, &num_meta_codes);

  // Step#3: Use the Huffman trees to decode the LZ77 encoded data.
  // TODO: Implementation of this method.
  ok = DecodeBackwardRefs(meta_codes, num_meta_codes, dec);

  // Step#4: Appply transforms on the decoded data.
  // TODO: Implementation of this method.
  ok = ApplyTransforms(dec);

  free(meta_codes);
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
  if (!DecodeImageStream(width, height, dec, data)) return 0;
  dec->argb_ = *data;

  return 1;
}
//------------------------------------------------------------------------------

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif
