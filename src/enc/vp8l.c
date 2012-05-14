// Copyright 2012 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// main entry for the lossless encoder.
//
// Author: Vikas Arora (vikaas.arora@gmail.com)
//

#ifdef USE_LOSSLESS_ENCODER

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "./backward_references.h"
#include "./vp8enci.h"
#include "./vp8li.h"
#include "../dsp/lossless.h"
#include "../utils/bit_writer.h"
#include "../utils/huffman_encode.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#define MAX_HUFF_IMAGE_SIZE       (32 * 1024 * 1024)

// TODO(vikas): find a common place between enc and dec for these:
#define PREDICTOR_TRANSFORM      0
#define CROSS_COLOR_TRANSFORM    1
#define SUBTRACT_GREEN           2
#define COLOR_INDEXING_TRANSFORM 3
#define TRANSFORM_PRESENT 1

#define IMAGE_SIZE_BITS 14

// -----------------------------------------------------------------------------
// Palette

static int CompareColors(const void* p1, const void* p2) {
  const uint32_t a = *(const uint32_t*)p1;
  const uint32_t b = *(const uint32_t*)p2;
  return (a < b) ? -1 : (a > b) ? 1 : 0;
}

// If number of colors in the image is less than or equal to MAX_PALETTE_SIZE,
// creates a palette and returns true, else returns false.
static int AnalyzeAndCreatePalette(const uint32_t* const argb, int num_pix,
                                   uint32_t palette[MAX_PALETTE_SIZE],
                                   int* const palette_size) {
  int i, key;
  int num_colors = 0;
  uint8_t in_use[MAX_PALETTE_SIZE * 4] = { 0 };
  uint32_t colors[MAX_PALETTE_SIZE * 4];
  static const uint32_t kHashMul = 0x1e35a7bd;

  key = (kHashMul * argb[0]) >> PALETTE_KEY_RIGHT_SHIFT;
  colors[key] = argb[0];
  in_use[key] = 1;
  ++num_colors;

  for (i = 1; i < num_pix; ++i) {
    if (argb[i] == argb[i - 1]) {
      continue;
    }
    key = (kHashMul * argb[i]) >> PALETTE_KEY_RIGHT_SHIFT;
    while (1) {
      if (!in_use[key]) {
        colors[key] = argb[i];
        in_use[key] = 1;
        ++num_colors;
        if (num_colors > MAX_PALETTE_SIZE) {
          return 0;
        }
        break;
      } else if (colors[key] == argb[i]) {
        // The color is already there.
        break;
      } else {
        // Some other color sits there.
        // Do linear conflict resolution.
        ++key;
        key &= (MAX_PALETTE_SIZE * 4 - 1);  // key mask for 1K buffer.
      }
    }
  }

  num_colors = 0;
  for (i = 0; i < (int)(sizeof(in_use) / sizeof(in_use[0])); ++i) {
    if (in_use[i]) {
      palette[num_colors] = colors[i];
      ++num_colors;
    }
  }

  qsort(palette, num_colors, sizeof(*palette), CompareColors);
  *palette_size = num_colors;
  return 1;
}

static int AnalyzeEntropy(const uint32_t const *argb, int xsize, int ysize,
                          double* const nonpredicted_bits,
                          double* const predicted_bits) {
  int i;
  VP8LHistogram* nonpredicted = NULL;
  VP8LHistogram* predicted = (VP8LHistogram*)malloc(2 * sizeof(*predicted));
  if (predicted == NULL) return 0;
  nonpredicted = predicted + 1;

  VP8LHistogramInit(predicted, 0);
  VP8LHistogramInit(nonpredicted, 0);
  for (i = 1; i < xsize * ysize; ++i) {
    const uint32_t pix = argb[i];
    const uint32_t pix_diff = VP8LSubPixels(pix, argb[i - 1]);
    if (pix_diff == 0) continue;
    if (i >= xsize && pix == argb[i - xsize]) {
      continue;
    }
    {
      const PixOrCopy pix_token = PixOrCopyCreateLiteral(pix);
      const PixOrCopy pix_diff_token = PixOrCopyCreateLiteral(pix_diff);
      VP8LHistogramAddSinglePixOrCopy(nonpredicted, &pix_token);
      VP8LHistogramAddSinglePixOrCopy(predicted, &pix_diff_token);
    }
  }
  *nonpredicted_bits = VP8LHistogramEstimateBitsBulk(nonpredicted);
  *predicted_bits = VP8LHistogramEstimateBitsBulk(predicted);
  free(predicted);
  return 1;
}

static int VP8LEncAnalyze(VP8LEncoder* const enc) {
  const WebPPicture* const pic = enc->pic_;
  assert(pic != NULL && pic->argb != NULL);

  enc->use_palette_ =
        AnalyzeAndCreatePalette(pic->argb, pic->width * pic->height,
                                enc->palette_, &enc->palette_size_);
  if (!enc->use_palette_) {
    double non_pred_entropy, pred_entropy;
    if (!AnalyzeEntropy(pic->argb, pic->width, pic->height,
                        &non_pred_entropy, &pred_entropy)) {
      return 0;
    }

    if (pred_entropy < 0.95 * non_pred_entropy) {
      enc->use_predict_ = 1;
      enc->use_cross_color_ = 1;
    }
  }
  return 1;
}

// -----------------------------------------------------------------------------

// Heuristics for selecting the stride ranges to collapse.
static int ValuesShouldBeCollapsedToStrideAverage(int a, int b) {
  return abs(a - b) < 4;
}

// Change the population counts in a way that the consequent
// Hufmann tree compression, especially its rle-part will be more
// likely to compress this data more efficiently.
//
// length contains the size of the histogram.
// data contains the population counts.
static int OptimizeHuffmanForRle(int length, int* counts) {
  int stride;
  int limit;
  int sum;
  uint8_t* good_for_rle;
  // 1) Let's make the Huffman code more compatible with rle encoding.
  int i;
  for (; length >= 0; --length) {
    if (length == 0) {
      return 1;  // All zeros.
    }
    if (counts[length - 1] != 0) {
      // Now counts[0..length - 1] does not have trailing zeros.
      break;
    }
  }
  // 2) Let's mark all population counts that already can be encoded
  // with an rle code.
  good_for_rle = (uint8_t*)calloc(length, 1);
  if (good_for_rle == NULL) {
    return 0;
  }
  {
    // Let's not spoil any of the existing good rle codes.
    // Mark any seq of 0's that is longer as 5 as a good_for_rle.
    // Mark any seq of non-0's that is longer as 7 as a good_for_rle.
    int symbol = counts[0];
    int stride = 0;
    for (i = 0; i < length + 1; ++i) {
      if (i == length || counts[i] != symbol) {
        if ((symbol == 0 && stride >= 5) ||
            (symbol != 0 && stride >= 7)) {
          int k;
          for (k = 0; k < stride; ++k) {
            good_for_rle[i - k - 1] = 1;
          }
        }
        stride = 1;
        if (i != length) {
          symbol = counts[i];
        }
      } else {
        ++stride;
      }
    }
  }
  // 3) Let's replace those population counts that lead to more rle codes.
  stride = 0;
  limit = counts[0];
  sum = 0;
  for (i = 0; i < length + 1; ++i) {
    if (i == length || good_for_rle[i] ||
        (i != 0 && good_for_rle[i - 1]) ||
        !ValuesShouldBeCollapsedToStrideAverage(counts[i], limit)) {
      if (stride >= 4 || (stride >= 3 && sum == 0)) {
        int k;
        // The stride must end, collapse what we have, if we have enough (4).
        int count = (sum + stride / 2) / stride;
        if (count < 1) {
          count = 1;
        }
        if (sum == 0) {
          // Don't make an all zeros stride to be upgraded to ones.
          count = 0;
        }
        for (k = 0; k < stride; ++k) {
          // We don't want to change value at counts[i],
          // that is already belonging to the next stride. Thus - 1.
          counts[i - k - 1] = count;
        }
      }
      stride = 0;
      sum = 0;
      if (i < length - 3) {
        // All interesting strides have a count of at least 4,
        // at least when non-zeros.
        limit = (counts[i] + counts[i + 1] +
                 counts[i + 2] + counts[i + 3] + 2) / 4;
      } else if (i < length) {
        limit = counts[i];
      } else {
        limit = 0;
      }
    }
    ++stride;
    if (i != length) {
      sum += counts[i];
      if (stride >= 4) {
        limit = (sum + stride / 2) / stride;
      }
    }
  }
  free(good_for_rle);
  return 1;
}

// TODO(vikasa): Wrap bit_codes and bit_lengths in a Struct.
static int GetHuffBitLengthsAndCodes(
    const VP8LHistogramSet* const histogram_image,
    int use_color_cache, int* const bit_length_sizes,
    uint16_t** const bit_codes, uint8_t** const bit_lengths) {
  int i, k;
  int ok = 1;
  int total_length_size = 0;
  uint8_t* mem_buf = NULL;
  const int histogram_image_size = histogram_image->size;

  // Iterate over all histograms and get the aggregate number of codes used.
  for (i = 0; i < histogram_image_size; ++i) {
    VP8LHistogram* const histo = histogram_image->histograms[i];
    const int num_literals = VP8LHistogramNumCodes(histo);
    k = 0;
    bit_length_sizes[5 * i] = num_literals;
    total_length_size += num_literals;
    for (k = 1; k < 5; ++k) {
      const int val = (k == 4) ? DISTANCE_CODES_MAX : 256;
      bit_length_sizes[5 * i + k] = val;
      total_length_size += val;
    }
  }

  // Allocate and Set bit_lengths and bit_codes.
  {
    uint16_t* codes;
    uint8_t* lengths;
    const size_t total_buf_size = total_length_size * sizeof(*lengths)
                                + total_length_size * sizeof(*codes);
    mem_buf = (uint8_t*)calloc(total_buf_size, 1);
    if (mem_buf == NULL) {
      ok = 0;
      goto End;
    }
    codes = (uint16_t*)mem_buf;
    lengths = (uint8_t*)&codes[total_length_size];
    for (i = 0; i < 5 * histogram_image_size; ++i) {
      const int bit_length = bit_length_sizes[i];
      bit_codes[i] = codes;
      bit_lengths[i] = lengths;
      codes += bit_length;
      lengths += bit_length;
    }
  }

  // Create Huffman trees.
  for (i = 0; i < histogram_image_size; ++i) {
    const int ix = 5 * i;
    VP8LHistogram* const histo = histogram_image->histograms[i];
    const int num_literals = bit_length_sizes[ix];
    // For each component, optimize histogram for Huffman with RLE compression.
    ok = ok && OptimizeHuffmanForRle(num_literals, histo->literal_);
    if (!use_color_cache) {
      // Implies that palette_bits == 0,
      // and so number of palette entries = (1 << 0) = 1.
      // Optimization might have smeared population count in this single
      // palette entry, so zero it out.
      histo->literal_[256 + kLengthCodes] = 0;
    }
    // Create a Huffman tree (in the form of bit lengths) for each component.
    ok = ok && VP8LCreateHuffmanTree(histo->literal_, num_literals, 15,
                                     bit_lengths[ix]);

    ok = ok && OptimizeHuffmanForRle(256, histo->red_);
    ok = ok && VP8LCreateHuffmanTree(histo->red_, 256, 15,
                                     bit_lengths[ix + 1]);

    ok = ok && OptimizeHuffmanForRle(256, histo->blue_);
    ok = ok && VP8LCreateHuffmanTree(histo->blue_, 256, 15,
                                     bit_lengths[ix + 2]);

    ok = ok && OptimizeHuffmanForRle(256, histo->alpha_);
    ok = ok && VP8LCreateHuffmanTree(histo->alpha_, 256, 15,
                                     bit_lengths[ix + 3]);

    ok = ok && OptimizeHuffmanForRle(DISTANCE_CODES_MAX, histo->distance_);
    ok = ok && VP8LCreateHuffmanTree(histo->distance_, DISTANCE_CODES_MAX, 15,
                                     bit_lengths[ix + 4]);

    // Create the actual bit codes for the bit lengths.
    for (k = 0; k < 5; ++k) {
      VP8LConvertBitDepthsToSymbols(bit_lengths[ix + k],
                                    bit_length_sizes[ix + k],
                                    bit_codes[ix + k]);
    }
  }

 End:
  if (!ok) free(mem_buf);
  return ok;
}

static void ClearHuffmanTreeIfOnlyOneSymbol(const int num_symbols,
                                            uint8_t* lengths,
                                            uint16_t* symbols) {
  int k;
  int count = 0;
  for (k = 0; k < num_symbols; ++k) {
    if (lengths[k] != 0) {
      ++count;
      if (count > 1) return;
    }
  }
  for (k = 0; k < num_symbols; ++k) {
    lengths[k] = 0;
    symbols[k] = 0;
  }
}

static void StoreHuffmanTreeOfHuffmanTreeToBitMask(
    VP8LBitWriter* const bw, const uint8_t* code_length_bitdepth) {
  // RFC 1951 will calm you down if you are worried about this funny sequence.
  // This sequence is tuned from that, but more weighted for lower symbol count,
  // and more spiking histograms.
  static const uint8_t kStorageOrder[CODE_LENGTH_CODES] = {
    17, 18, 0, 1, 2, 3, 4, 5, 16, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
  };
  int i;
  // Throw away trailing zeros:
  int codes_to_store = CODE_LENGTH_CODES;
  for (; codes_to_store > 4; --codes_to_store) {
    if (code_length_bitdepth[kStorageOrder[codes_to_store - 1]] != 0) {
      break;
    }
  }
  VP8LWriteBits(bw, 4, codes_to_store - 4);
  for (i = 0; i < codes_to_store; ++i) {
    VP8LWriteBits(bw, 3, code_length_bitdepth[kStorageOrder[i]]);
  }
}

static void StoreHuffmanTreeToBitMask(
    VP8LBitWriter* const bw,
    const HuffmanTreeToken* const tokens,
    const int num_tokens,
    const uint8_t* code_length_bitdepth,
    const uint16_t* code_length_bitdepth_symbols) {
  int i;
  for (i = 0; i < num_tokens; ++i) {
    const int ix = tokens[i].code;
    const int extra_bits = tokens[i].extra_bits;
    VP8LWriteBits(bw, code_length_bitdepth[ix],
                  code_length_bitdepth_symbols[ix]);
    switch (ix) {
      case 16:
        VP8LWriteBits(bw, 2, extra_bits);
        break;
      case 17:
        VP8LWriteBits(bw, 3, extra_bits);
        break;
      case 18:
        VP8LWriteBits(bw, 7, extra_bits);
        break;
    }
  }
}

static int StoreFullHuffmanCode(VP8LBitWriter* const bw,
                                const uint8_t* const bit_lengths,
                                int bit_lengths_size) {
  int ok = 0;
  uint8_t code_length_bitdepth[CODE_LENGTH_CODES] = { 0 };
  uint16_t code_length_bitdepth_symbols[CODE_LENGTH_CODES] = { 0 };
  int num_tokens;
  HuffmanTreeToken* const tokens =
      (HuffmanTreeToken*)malloc(bit_lengths_size * sizeof(*tokens));
  if (tokens == NULL) return 0;

  VP8LWriteBits(bw, 1, 0);
  num_tokens = VP8LCreateCompressedHuffmanTree(bit_lengths, bit_lengths_size,
                                               tokens, bit_lengths_size);
  {
    int histogram[CODE_LENGTH_CODES] = { 0 };
    int i;
    for (i = 0; i < num_tokens; ++i) {
      ++histogram[tokens[i].code];
    }

    if (!VP8LCreateHuffmanTree(histogram, CODE_LENGTH_CODES,
                               7, code_length_bitdepth)) {
      goto End;
    }
  }
  VP8LConvertBitDepthsToSymbols(code_length_bitdepth, CODE_LENGTH_CODES,
                                code_length_bitdepth_symbols);
  StoreHuffmanTreeOfHuffmanTreeToBitMask(bw, code_length_bitdepth);
  ClearHuffmanTreeIfOnlyOneSymbol(CODE_LENGTH_CODES,
                                  code_length_bitdepth,
                                  code_length_bitdepth_symbols);
  {
    int trailing_zero_bits = 0;
    int trimmed_length = num_tokens;
    int write_trimmed_length;
    int length;
    int i = num_tokens;
    while (i-- > 0) {
      const int ix = tokens[i].code;
      if (ix == 0 || ix == 17 || ix == 18) {
        --trimmed_length;   // discount trailing zeros
        trailing_zero_bits += code_length_bitdepth[ix];
        if (ix == 17) {
          trailing_zero_bits += 3;
        } else if (ix == 18) {
          trailing_zero_bits += 7;
        }
      } else {
        break;
      }
    }
    write_trimmed_length = (trimmed_length > 1 && trailing_zero_bits > 12);
    length = write_trimmed_length ? trimmed_length : num_tokens;
    VP8LWriteBits(bw, 1, write_trimmed_length);
    if (write_trimmed_length) {
      const int nbits = VP8LBitsLog2Ceiling(trimmed_length - 1);
      const int nbitpairs = (nbits == 0) ? 1 : (nbits + 1) / 2;
      VP8LWriteBits(bw, 3, nbitpairs - 1);
      VP8LWriteBits(bw, nbitpairs * 2, trimmed_length - 2);
    }
    StoreHuffmanTreeToBitMask(bw, tokens,
                              length, code_length_bitdepth,
                              code_length_bitdepth_symbols);
  }
  ok = 1;
 End:
  free(tokens);
  return ok;
}

static int StoreHuffmanCode(VP8LBitWriter* const bw,
                            const uint8_t* const bit_lengths,
                            int bit_lengths_size) {
  int i;
  int count = 0;
  int symbols[2] = { 0, 0 };
  const int kMaxBits = 8;
  const int kMaxSymbol = 1 << kMaxBits;

  // Check whether it's a small tree.
  for (i = 0; i < bit_lengths_size && count < 3; ++i) {
    if (bit_lengths[i] != 0) {
      if (count < 2) symbols[count] = i;
      ++count;
    }
  }

  if (count == 0) {   // emit minimal tree for empty cases
    // bits: small tree marker: 1, count-1: 0, large 8-bit code: 0, code: 0
    VP8LWriteBits(bw, 4, 0x01);
    return 1;
  } else if (count <= 2 && symbols[0] < kMaxSymbol && symbols[1] < kMaxSymbol) {
    VP8LWriteBits(bw, 1, 1);  // Small tree marker to encode 1 or 2 symbols.
    VP8LWriteBits(bw, 1, count - 1);
    if (symbols[0] <= 1) {
      VP8LWriteBits(bw, 1, 0);  // Code bit for small (1 bit) symbol value.
      VP8LWriteBits(bw, 1, symbols[0]);
    } else {
      VP8LWriteBits(bw, 1, 1);
      VP8LWriteBits(bw, 8, symbols[0]);
    }
    if (count == 2) {
      VP8LWriteBits(bw, 8, symbols[1]);
    }
    return 1;
  } else {
    return StoreFullHuffmanCode(bw, bit_lengths, bit_lengths_size);
  }
}

static void StoreImageToBitMask(
    VP8LBitWriter* const bw, int width, int histo_bits,
    const VP8LBackwardRefs* const refs,
    const uint16_t* histogram_symbols,
    uint8_t** const bitdepths, uint16_t** const bit_symbols) {
  // x and y trace the position in the image.
  int x = 0;
  int y = 0;
  const int histo_xsize = histo_bits ? VP8LSubSampleSize(width, histo_bits) : 1;
  int i;
  for (i = 0; i < refs->size; ++i) {
    const PixOrCopy* const v = &refs->refs[i];
    const int histogram_ix = histogram_symbols[histo_bits ?
                                               (y >> histo_bits) * histo_xsize +
                                               (x >> histo_bits) : 0];
    if (PixOrCopyIsCacheIdx(v)) {
      const int code = PixOrCopyCacheIdx(v);
      const int literal_ix = 256 + kLengthCodes + code;
      VP8LWriteBits(bw, bitdepths[5 * histogram_ix][literal_ix],
                    bit_symbols[5 * histogram_ix][literal_ix]);
    } else if (PixOrCopyIsLiteral(v)) {
      static const int order[] = { 1, 2, 0, 3 };
      int k;
      for (k = 0; k < 4; ++k) {
        const int code = PixOrCopyLiteral(v, order[k]);
        VP8LWriteBits(bw, bitdepths[5 * histogram_ix + k][code],
                      bit_symbols[5 * histogram_ix + k][code]);
      }
    } else {
      int bits, n_bits;
      int code, distance;
      int len_ix;
      PrefixEncode(v->len, &code, &n_bits, &bits);
      len_ix = 256 + code;
      VP8LWriteBits(bw, bitdepths[5 * histogram_ix][len_ix],
                    bit_symbols[5 * histogram_ix][len_ix]);
      VP8LWriteBits(bw, n_bits, bits);

      distance = PixOrCopyDistance(v);
      PrefixEncode(distance, &code, &n_bits, &bits);
      VP8LWriteBits(bw, bitdepths[5 * histogram_ix + 4][code],
                    bit_symbols[5 * histogram_ix + 4][code]);
      VP8LWriteBits(bw, n_bits, bits);
    }
    x += PixOrCopyLength(v);
    while (x >= width) {
      x -= width;
      ++y;
    }
  }
}

static int EncodeImageInternal(VP8LBitWriter* const bw,
                               const uint32_t* const argb,
                               int width, int height, int quality,
                               int cache_bits, int histogram_bits) {
  int i;
  int ok = 0;
  int write_histogram_image;
  int* bit_lengths_sizes = NULL;
  uint8_t** bit_lengths = NULL;
  uint16_t** bit_codes = NULL;
  const int use_2d_locality = 1;
  const int use_color_cache = (cache_bits > 0);
  const int color_cache_size = use_color_cache ? (1 << cache_bits) : 0;
  const int histogram_image_xysize =
      VP8LSubSampleSize(width, histogram_bits) *
      VP8LSubSampleSize(height, histogram_bits);
  VP8LHistogramSet* histogram_image =
      VP8LAllocateHistogramSet(histogram_image_xysize, 0);
  int histogram_image_size = 0;
  int bit_array_size = 0;
  VP8LBackwardRefs refs;
  uint16_t* const histogram_symbols =
      (uint16_t*)malloc(histogram_image_xysize * sizeof(*histogram_symbols));

  if (histogram_image == NULL || histogram_symbols == NULL) goto Error;

  // Calculate backward references from ARGB image.
  if (!VP8LGetBackwardReferences(width, height, argb, quality, cache_bits,
                                 use_2d_locality, &refs)) {
    goto Error;
  }
  // Build histogram image & symbols from backward references.
  if (!VP8LGetHistoImageSymbols(width, height, &refs,
                                quality, histogram_bits, cache_bits,
                                histogram_image,
                                histogram_symbols)) {
    goto Error;
  }
  // Create Huffman bit lengths & codes for each histogram image.
  histogram_image_size = histogram_image->size;
  bit_array_size = 5 * histogram_image_size;
  bit_lengths_sizes = (int*)calloc(bit_array_size, sizeof(*bit_lengths_sizes));
  bit_lengths = (uint8_t**)calloc(bit_array_size, sizeof(*bit_lengths));
  bit_codes = (uint16_t**)calloc(bit_array_size, sizeof(*bit_codes));
  if (bit_lengths_sizes == NULL || bit_lengths == NULL || bit_codes == NULL ||
      !GetHuffBitLengthsAndCodes(histogram_image, use_color_cache,
                                 bit_lengths_sizes, bit_codes, bit_lengths)) {
    goto Error;
  }

  // Color Cache parameters.
  VP8LWriteBits(bw, 1, use_color_cache);
  if (use_color_cache) {
    VP8LWriteBits(bw, 4, cache_bits);
  }

  // Huffman image + meta huffman.
  write_histogram_image = (histogram_image_size > 1);
  VP8LWriteBits(bw, 1, write_histogram_image);
  if (write_histogram_image) {
    uint32_t* const histogram_argb =
        (uint32_t*)malloc(histogram_image_xysize * sizeof(*histogram_argb));
    int max_index = 0;
    if (histogram_argb == NULL) goto Error;
    for (i = 0; i < histogram_image_xysize; ++i) {
      const int index = histogram_symbols[i] & 0xffff;
      histogram_argb[i] = 0xff000000 | (index << 8);
      if (index >= max_index) {
        max_index = index + 1;
      }
    }
    histogram_image_size = max_index;

    VP8LWriteBits(bw, 4, histogram_bits);
    ok = EncodeImageInternal(bw, histogram_argb,
                             VP8LSubSampleSize(width, histogram_bits),
                             VP8LSubSampleSize(height, histogram_bits),
                             quality, 0, 0);
    free(histogram_argb);
    if (!ok) goto Error;
  }

  // Store Huffman codes.
  for (i = 0; i < histogram_image_size; ++i) {
    int k;
    for (k = 0; k < 5; ++k) {
      const uint8_t* const cur_bit_lengths =  bit_lengths[5 * i + k];
      const int cur_bit_lengths_size = (k == 0) ?
                   256 + kLengthCodes + color_cache_size :
                   bit_lengths_sizes[5 * i + k];
      if (!StoreHuffmanCode(bw, cur_bit_lengths, cur_bit_lengths_size)) {
        goto Error;
      }
    }
  }

  // Free combined histograms.
  free(histogram_image);
  histogram_image = NULL;

  // Emit no bits if there is only one symbol in the histogram.
  // This gives better compression for some images.
  for (i = 0; i < 5 * histogram_image_size; ++i) {
    ClearHuffmanTreeIfOnlyOneSymbol(bit_lengths_sizes[i], bit_lengths[i],
                                    bit_codes[i]);
  }
  // Store actual literals.
  StoreImageToBitMask(bw, width, histogram_bits, &refs,
                      histogram_symbols, bit_lengths, bit_codes);
  ok = 1;

 Error:
  if (!ok) free(histogram_image);

  VP8LClearBackwardRefs(&refs);
  free(bit_lengths_sizes);
  free(bit_lengths);
  free(*bit_codes);
  free(bit_codes);
  free(histogram_symbols);
  return ok;
}

// -----------------------------------------------------------------------------
// Transforms

// Check if it would be a good idea to subtract green from red and blue. We
// only impact entropy in red/blue components, don't bother to look at others.
static int EvalAndApplySubtractGreen(const VP8LEncoder* const enc,
                                     int width, int height,
                                     VP8LBitWriter* const bw) {
  if (!enc->use_palette_) {
    int i;
    const uint32_t* const argb = enc->argb_;
    double bit_cost_before, bit_cost_after;
    VP8LHistogram* const histo = (VP8LHistogram*)malloc(sizeof(*histo));
    if (histo == NULL) return 0;

    VP8LHistogramInit(histo, 1);
    for (i = 0; i < width * height; ++i) {
      const uint32_t c = argb[i];
      ++histo->red_[(c >> 16) & 0xff];
      ++histo->blue_[(c >> 0) & 0xff];
    }
    bit_cost_before = VP8LHistogramEstimateBits(histo);

    VP8LHistogramInit(histo, 1);
    for (i = 0; i < width * height; ++i) {
      const uint32_t c = argb[i];
      const int green = (c >> 8) & 0xff;
      ++histo->red_[((c >> 16) - green) & 0xff];
      ++histo->blue_[((c >> 0) - green) & 0xff];
    }
    bit_cost_after = VP8LHistogramEstimateBits(histo);
    free(histo);

    // Check if subtracting green yields low entropy.
    if (bit_cost_after < bit_cost_before) {
      VP8LWriteBits(bw, 1, TRANSFORM_PRESENT);
      VP8LWriteBits(bw, 2, SUBTRACT_GREEN);
      VP8LSubtractGreenFromBlueAndRed(enc->argb_, width * height);
    }
  }
  return 1;
}

static int ApplyPredictFilter(const VP8LEncoder* const enc,
                              int width, int height, int quality,
                              VP8LBitWriter* const bw) {
  const int pred_bits = enc->transform_bits_;
  const int transform_width = VP8LSubSampleSize(width, pred_bits);
  const int transform_height = VP8LSubSampleSize(height, pred_bits);

  VP8LResidualImage(width, height, pred_bits, enc->argb_, enc->argb_scratch_,
                    enc->transform_data_);
  VP8LWriteBits(bw, 1, TRANSFORM_PRESENT);
  VP8LWriteBits(bw, 2, PREDICTOR_TRANSFORM);
  VP8LWriteBits(bw, 4, pred_bits);
  if (!EncodeImageInternal(bw, enc->transform_data_,
                           transform_width, transform_height, quality, 0, 0)) {
    return 0;
  }
  return 1;
}

static int ApplyCrossColorFilter(const VP8LEncoder* const enc,
                                 int width, int height, int quality,
                                 VP8LBitWriter* const bw) {
  const int ccolor_transform_bits = enc->transform_bits_;
  const int transform_width = VP8LSubSampleSize(width, ccolor_transform_bits);
  const int transform_height = VP8LSubSampleSize(height, ccolor_transform_bits);
  const int step = (quality == 0) ? 32 : 8;

  VP8LColorSpaceTransform(width, height, ccolor_transform_bits, step,
                          enc->argb_, enc->transform_data_);
  VP8LWriteBits(bw, 1, TRANSFORM_PRESENT);
  VP8LWriteBits(bw, 2, CROSS_COLOR_TRANSFORM);
  VP8LWriteBits(bw, 4, ccolor_transform_bits);
  if (!EncodeImageInternal(bw, enc->transform_data_,
                           transform_width, transform_height, quality, 0, 0)) {
    return 0;
  }
  return 1;
}

// -----------------------------------------------------------------------------

static void PutLE32(uint8_t* const data, uint32_t val) {
  data[0] = (val >>  0) & 0xff;
  data[1] = (val >>  8) & 0xff;
  data[2] = (val >> 16) & 0xff;
  data[3] = (val >> 24) & 0xff;
}

static WebPEncodingError WriteRiffHeader(const VP8LEncoder* const enc,
                                         size_t riff_size, size_t vp8l_size) {
  const WebPPicture* const pic = enc->pic_;
  uint8_t riff[HEADER_SIZE + SIGNATURE_SIZE] = {
    'R', 'I', 'F', 'F', 0, 0, 0, 0, 'W', 'E', 'B', 'P',
    'V', 'P', '8', 'L', 0, 0, 0, 0, LOSSLESS_MAGIC_BYTE,
  };
  if (riff_size < (vp8l_size + TAG_SIZE + CHUNK_HEADER_SIZE)) {
    return VP8_ENC_ERROR_INVALID_CONFIGURATION;
  }
  PutLE32(riff + TAG_SIZE, (uint32_t)riff_size);
  PutLE32(riff + RIFF_HEADER_SIZE + TAG_SIZE, (uint32_t)vp8l_size);
  if (!pic->writer(riff, sizeof(riff), pic)) {
    return VP8_ENC_ERROR_BAD_WRITE;
  }
  return VP8_ENC_OK;
}

static void WriteImageSize(VP8LEncoder* const enc, VP8LBitWriter* const bw) {
  WebPPicture* const pic = enc->pic_;
  const int width = pic->width - 1;
  const int height = pic->height -1;
  assert(width < WEBP_MAX_DIMENSION && height < WEBP_MAX_DIMENSION);

  VP8LWriteBits(bw, IMAGE_SIZE_BITS, width);
  VP8LWriteBits(bw, IMAGE_SIZE_BITS, height);
}

static WebPEncodingError WriteImage(const VP8LEncoder* const enc,
                                    VP8LBitWriter* const bw,
                                    size_t* const coded_size) {
  size_t riff_size, vp8l_size, webpll_size, pad;
  const WebPPicture* const pic = enc->pic_;
  WebPEncodingError err = VP8_ENC_OK;
  const uint8_t* const webpll_data = VP8LBitWriterFinish(bw);

  webpll_size = VP8LBitWriterNumBytes(bw);
  vp8l_size = SIGNATURE_SIZE + webpll_size;
  pad = vp8l_size & 1;
  vp8l_size += pad;

  riff_size = TAG_SIZE + CHUNK_HEADER_SIZE + vp8l_size;
  err = WriteRiffHeader(enc, riff_size, vp8l_size);
  if (err != VP8_ENC_OK) goto Error;

  if (!pic->writer(webpll_data, webpll_size, pic)) {
    err = VP8_ENC_ERROR_BAD_WRITE;
    goto Error;
  }

  if (pad) {
    const uint8_t pad_byte[1] = { 0 };
    if (!pic->writer(pad_byte, 1, pic)) {
      err = VP8_ENC_ERROR_BAD_WRITE;
      goto Error;
    }
  }
  *coded_size = vp8l_size;
  return VP8_ENC_OK;

 Error:
  return err;
}

// -----------------------------------------------------------------------------

// Allocates the memory for argb (W x H) buffer, 2 rows of context for
// prediction and transform data.
static WebPEncodingError AllocateTransformBuffer(VP8LEncoder* const enc,
                                                 int width, int height) {
  WebPEncodingError err = VP8_ENC_OK;
  const size_t tile_size = 1 << enc->transform_bits_;
  const size_t image_size = height * width;
  const size_t argb_scratch_size = (tile_size + 1) * width;
  const size_t transform_data_size =
      VP8LSubSampleSize(height, enc->transform_bits_) *
      VP8LSubSampleSize(width, enc->transform_bits_);
  const size_t total_size =
      image_size + argb_scratch_size + transform_data_size;
  uint32_t* mem = (uint32_t*)malloc(total_size * sizeof(*mem));
  if (mem == NULL) {
    err = VP8_ENC_ERROR_OUT_OF_MEMORY;
    goto Error;
  }
  enc->argb_ = mem;
  mem += image_size;
  enc->argb_scratch_ = mem;
  mem += argb_scratch_size;
  enc->transform_data_ = mem;
  enc->current_width_ = width;

 Error:
  return err;
}

// Bundles multiple (2, 4 or 8) pixels into a single pixel.
// Returns the new xsize.
static void BundleColorMap(const uint32_t* const argb,
                           int width, int height, int xbits,
                           uint32_t* bundled_argb, int xs) {
  int x, y;
  const int bit_depth = 1 << (3 - xbits);
  uint32_t code = 0;

  for (y = 0; y < height; ++y) {
    for (x = 0; x < width; ++x) {
      const int mask = (1 << xbits) - 1;
      const int xsub = x & mask;
      if (xsub == 0) {
        code = 0;
      }
      // TODO(vikasa): simplify the bundling logic.
      code |= (argb[y * width + x] & 0xff00) << (bit_depth * xsub);
      bundled_argb[y * xs + (x >> xbits)] = 0xff000000 | code;
    }
  }
}

// Note: Expects "enc->palette_" to be set properly.
// Also, "enc->palette_" will be modified after this call and should not be used
// later.
static WebPEncodingError ApplyPalette(VP8LBitWriter* const bw,
                                      VP8LEncoder* const enc,
                                      int width, int height, int quality) {
  WebPEncodingError err = VP8_ENC_OK;
  int i;
  uint32_t* const argb = enc->pic_->argb;
  uint32_t* const palette = enc->palette_;
  const int palette_size = enc->palette_size_;

  // Replace each input pixel by corresponding palette index.
  for (i = 0; i < width * height; ++i) {
    int k;
    for (k = 0; k < palette_size; ++k) {
      const uint32_t pix = argb[i];
      if (pix == palette[k]) {
        argb[i] = 0xff000000u | (k << 8);
        break;
      }
    }
  }

  // Save palette to bitstream.
  VP8LWriteBits(bw, 1, TRANSFORM_PRESENT);
  VP8LWriteBits(bw, 2, COLOR_INDEXING_TRANSFORM);
  VP8LWriteBits(bw, 8, palette_size - 1);
  for (i = palette_size - 1; i >= 1; --i) {
    palette[i] = VP8LSubPixels(palette[i], palette[i - 1]);
  }
  if (!EncodeImageInternal(bw, palette, palette_size, 1, quality, 0, 0)) {
    err = VP8_ENC_ERROR_INVALID_CONFIGURATION;
    goto Error;
  }

  if (palette_size <= 16) {
    // Image can be packed (multiple pixels per uint32_t).
    int xbits = 1;
    if (palette_size <= 2) {
      xbits = 3;
    } else if (palette_size <= 4) {
      xbits = 2;
    }
    err = AllocateTransformBuffer(enc, VP8LSubSampleSize(width, xbits), height);
    if (err != VP8_ENC_OK) goto Error;
    BundleColorMap(argb, width, height, xbits, enc->argb_, enc->current_width_);
  }

 Error:
  return err;
}

// -----------------------------------------------------------------------------

static int GetHistoBits(const WebPConfig* const config,
                        const WebPPicture* const pic) {
  const int width = pic->width;
  const int height = pic->height;
  const size_t hist_size = sizeof(VP8LHistogram);
  int histo_bits = 9 - (int)(config->quality / 16.f + .5f);
  while (1) {
    const size_t huff_image_size = VP8LSubSampleSize(width, histo_bits) *
                                   VP8LSubSampleSize(height, histo_bits) *
                                   hist_size;
    if (huff_image_size <= MAX_HUFF_IMAGE_SIZE) break;
    ++histo_bits;
  }
  return (histo_bits < 3) ? 3 : (histo_bits > 10) ? 10 : histo_bits;
}

static void InitEncParams(VP8LEncoder* const enc) {
  const WebPConfig* const config = enc->config_;
  const WebPPicture* const picture = enc->pic_;
  const int method = config->method;
  const float quality = config->quality;
  enc->transform_bits_ = (method < 4) ? 5 : (method > 4) ? 3 : 4;
  enc->histo_bits_ = GetHistoBits(config, picture);
  enc->cache_bits_ = (quality <= 25.f) ? 0 : 7;
}

// -----------------------------------------------------------------------------
// VP8LEncoder

static VP8LEncoder* NewVP8LEncoder(const WebPConfig* const config,
                                   WebPPicture* const picture) {
  VP8LEncoder* const enc = (VP8LEncoder*)calloc(1, sizeof(*enc));
  if (enc == NULL) {
    WebPEncodingSetError(picture, VP8_ENC_ERROR_OUT_OF_MEMORY);
    return NULL;
  }
  enc->config_ = config;
  enc->pic_ = picture;
  return enc;
}

static void DeleteVP8LEncoder(VP8LEncoder* enc) {
  free(enc->argb_);
  free(enc);
}

// -----------------------------------------------------------------------------
// Main call

int VP8LEncodeImage(const WebPConfig* const config,
                    WebPPicture* const picture) {
  int ok = 0;
  int width, height, quality;
  size_t coded_size;
  VP8LEncoder* enc = NULL;
  WebPEncodingError err = VP8_ENC_OK;
  VP8LBitWriter bw;

  if (config == NULL || picture == NULL) return 0;

  if (picture->argb == NULL) {
    err = VP8_ENC_ERROR_NULL_PARAMETER;
    goto Error;
  }

  enc = NewVP8LEncoder(config, picture);
  if (enc == NULL) {
    err = VP8_ENC_ERROR_OUT_OF_MEMORY;
    goto Error;
  }
  width = picture->width;
  height = picture->height;
  quality = config->quality;

  InitEncParams(enc);

  // ---------------------------------------------------------------------------
  // Analyze image (entropy, num_palettes etc)

  if (!VP8LEncAnalyze(enc)) {
    err = VP8_ENC_ERROR_OUT_OF_MEMORY;
    goto Error;
  }

  // Write image size.
  VP8LBitWriterInit(&bw, (width * height) >> 1);
  WriteImageSize(enc, &bw);

  if (enc->use_palette_) {
    err = ApplyPalette(&bw, enc, width, height, quality);
    if (err != VP8_ENC_OK) goto Error;
    enc->cache_bits_ = 0;
  }

  // In case image is not packed.
  if (enc->argb_ == NULL) {
    const size_t image_size = height * width;
    err = AllocateTransformBuffer(enc, width, height);
    if (err != VP8_ENC_OK) goto Error;
    memcpy(enc->argb_, picture->argb, image_size * sizeof(*enc->argb_));
    enc->current_width_ = width;
  }

  // ---------------------------------------------------------------------------
  // Apply transforms and write transform data.

  if (!EvalAndApplySubtractGreen(enc, enc->current_width_, height, &bw)) {
    err = VP8_ENC_ERROR_OUT_OF_MEMORY;
    goto Error;
  }

  if (enc->use_predict_) {
    if (!ApplyPredictFilter(enc, enc->current_width_, height, quality, &bw)) {
      err = VP8_ENC_ERROR_INVALID_CONFIGURATION;
      goto Error;
    }
  }

  if (enc->use_cross_color_) {
    if (!ApplyCrossColorFilter(enc, enc->current_width_, height, quality,
                               &bw)) {
      err = VP8_ENC_ERROR_INVALID_CONFIGURATION;
      goto Error;
    }
  }

  VP8LWriteBits(&bw, 1, !TRANSFORM_PRESENT);  // No more transforms.

  // ---------------------------------------------------------------------------
  // Estimate the color cache size.

  if (enc->cache_bits_ > 0) {
    if (!VP8LCalculateEstimateForCacheSize(enc->argb_, enc->current_width_,
                                           height, &enc->cache_bits_)) {
      err = VP8_ENC_ERROR_INVALID_CONFIGURATION;
      goto Error;
    }
  }

  // ---------------------------------------------------------------------------
  // Encode and write the transformed image.

  ok = EncodeImageInternal(&bw, enc->argb_, enc->current_width_, height,
                           quality, enc->cache_bits_, enc->histo_bits_);
  if (!ok) goto Error;

  err = WriteImage(enc, &bw, &coded_size);
  if (err != VP8_ENC_OK) {
    ok = 0;
    goto Error;
  }

  if (picture->stats != NULL) {
    WebPAuxStats* const stats = picture->stats;
    memset(stats, 0, sizeof(*stats));
    stats->PSNR[0] = 99.;
    stats->PSNR[1] = 99.;
    stats->PSNR[2] = 99.;
    stats->PSNR[3] = 99.;
    stats->coded_size = coded_size;
  }

  if (picture->extra_info != NULL) {
    const int mb_w = (width + 15) >> 4;
    const int mb_h = (height + 15) >> 4;
    memset(picture->extra_info, 0, mb_w * mb_h * sizeof(*picture->extra_info));
  }

 Error:
  VP8LBitWriterDestroy(&bw);
  DeleteVP8LEncoder(enc);
  if (!ok) {
    assert(err != VP8_ENC_OK);
    WebPEncodingSetError(picture, err);
  }
  return ok;
}

//------------------------------------------------------------------------------

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif
