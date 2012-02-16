// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)
//
// This module is mostly for actually writing the bits, but has also
// some algorithmic code and some higher level code.

#include "./encode.h"

#include <stdlib.h>
#include <string.h>

#include "../common/integral_types.h"
#include "../common/predictor.h"
#include "./backward_references.h"
#include "./bit_writer.h"
#include "./entropy_encode.h"
#include "./histogram_image.h"
#include "./predictor.h"

static void PutLE32(uint8_t* const data, uint32_t val) {
  data[0] = (val >>  0) & 0xff;
  data[1] = (val >>  8) & 0xff;
  data[2] = (val >> 16) & 0xff;
  data[3] = (val >> 24) & 0xff;
}

static void PutRiffHeader(uint8_t *bytes, size_t webpll_size) {
  uint8_t riff[HEADER_SIZE + SIGNATURE_SIZE] = {
    'R', 'I', 'F', 'F', 0, 0, 0, 0, 'W', 'E', 'B', 'P',
    'V', 'P', '8', 'L', 0, 0, 0, 0, 0x64
  };
  const size_t riff_size = TAG_SIZE + CHUNK_HEADER_SIZE + webpll_size;
  PutLE32(riff + TAG_SIZE, riff_size);
  PutLE32(riff + RIFF_HEADER_SIZE + TAG_SIZE, webpll_size);
  memcpy(bytes, riff, sizeof(riff));
}

// Heuristics for selecting the stride ranges to collapse.
static int ValuesShouldBeCollapsedToStrideAverage(int a, int b) {
  return abs(a - b) < 4;
}

// Change the population counts in a way that the consequent
// Hufmann tree compression, especially its rle-part will be more
// likely to compress this data more efficiently.
//
// length containts the size of the histogram.
// data contains the population counts.
static int OptimizeHuffmanForRle(int length, int *counts) {
  int stride;
  int limit;
  int sum;
  char *good_for_rle;
  // Let's make the Huffman code more compatible with rle encoding.
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
  good_for_rle = (char *)malloc(length);
  if (!good_for_rle) {
    return 0;
  }
  memset(good_for_rle, 0, length);
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

static void ClearHuffmanTreeIfOnlyOneSymbol(const int num_symbols,
                                            uint8_t* lengths,
                                            uint16_t* symbols) {
  int count = 0;
  int k;
  for (k = 0; k < num_symbols; ++k) {
    if (lengths[k] != 0) ++count;
  }
  if (count == 1) {
    for (k = 0; k < num_symbols; ++k) {
      lengths[k] = 0;
      symbols[k] = 0;
    }
  }
}

static int PackLiteralBitLengths(const uint8_t* bit_lengths,
                                 int palette_bits, int use_palette,
                                 int *new_length_size,
                                 uint8_t** new_lengths) {
  int i;
  int length_code_offset = 256;
  *new_length_size = 256 + kLengthCodes;
  if (use_palette) {
    *new_length_size += 1 << palette_bits;
  }
  *new_lengths = (uint8_t*)malloc(*new_length_size);
  if (*new_lengths == NULL) {
    return 0;
  }
  for (i = 0; i < 256; ++i) {
    (*new_lengths)[i] = bit_lengths[i];
  }
  if (use_palette) {
    for (i = 0; i < (1 << palette_bits); ++i) {
      (*new_lengths)[256 + i] = bit_lengths[256 + i];
    }
    length_code_offset += 1 << palette_bits;
  }
  for (i = 0; i < kLengthCodes; ++i) {
    (*new_lengths)[length_code_offset + i] =
        bit_lengths[256 + (1 << palette_bits) + i];
  }
  return 1;
}

static void StoreHuffmanTreeOfHuffmanTreeToBitMask(
    const uint8_t *code_length_bitdepth, BitWriter* const bw) {
  // RFC 1951 will calm you down if you are worried about this funny sequence.
  // This sequence is tuned from that, but more weighted for lower symbol count,
  // and more spiking histograms.
  int i;
  static const uint8_t kStorageOrder[CODE_LENGTH_CODES] = {
    17, 18, 0, 1, 2, 3, 4, 5, 16, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
  };
  // Throw away trailing zeros:
  int codes_to_store = sizeof(kStorageOrder);
  for (; codes_to_store > 4; --codes_to_store) {
    if (code_length_bitdepth[kStorageOrder[codes_to_store - 1]] != 0) {
      break;
    }
  }
  // How many code length codes we write above the first four (see RFC 1951).
  WriteBits(4, codes_to_store - 4, bw);
  for (i = 0; i < codes_to_store; ++i) {
    WriteBits(3, code_length_bitdepth[kStorageOrder[i]], bw);
  }
}

static void StoreHuffmanTreeToBitMask(
    const uint8_t *huffman_tree,
    const uint8_t *huffman_tree_extra_bits,
    const int num_symbols,
    const uint8_t *code_length_bitdepth,
    const uint16_t *code_length_bitdepth_symbols,
    BitWriter* const bw) {
  int i;
  for (i = 0; i < num_symbols; ++i) {
    int ix = huffman_tree[i];
    WriteBits(code_length_bitdepth[ix], code_length_bitdepth_symbols[ix], bw);
    switch (ix) {
      case 16:
        WriteBits(2, huffman_tree_extra_bits[i], bw);
        break;
      case 17:
        WriteBits(3, huffman_tree_extra_bits[i], bw);
        break;
      case 18:
        WriteBits(7, huffman_tree_extra_bits[i], bw);
        break;
    }
  }
}

static int StoreHuffmanCode(uint8_t *bit_lengths, int bit_lengths_size,
                            BitWriter* const bw) {
  int ok = 1;
  int count = 0;
  int symbols[2] = { 0, 0 };
  int i;
  uint8_t* huffman_tree = (uint8_t*)malloc(bit_lengths_size);
  uint8_t* huffman_tree_extra_bits = (uint8_t*)malloc(bit_lengths_size);
  int huffman_tree_size = 0;
  uint8_t code_length_bitdepth[CODE_LENGTH_CODES];
  uint16_t code_length_bitdepth_symbols[CODE_LENGTH_CODES];
  int huffman_tree_histogram[CODE_LENGTH_CODES];
  if (huffman_tree == NULL ||
      huffman_tree_extra_bits == NULL) {
    ok = 0;
    goto exit_label;
  }
  for (i = 0; i < bit_lengths_size; ++i) {
    if (bit_lengths[i] != 0) {
      if (count < 2) symbols[count] = i;
      ++count;
    }
  }
  if (count <= 2) {
    int num_bits = 4;
    // 0, 1 or 2 symbols to encode.
    WriteBits(1, 1, bw);
    if (count == 0) {
      WriteBits(4, 0, bw);
      return 1;
    }
    WriteBits(1, count - 1, bw);
    while (symbols[count - 1] >= (1 << num_bits)) num_bits += 2;
    WriteBits(3, (num_bits - 4) / 2 + 1, bw);
    for (i = 0; i < count; ++i) {
      WriteBits(num_bits, symbols[i], bw);
    }
    return 1;
  }
  WriteBits(1, 0, bw);
  CreateCompressedHuffmanTree(bit_lengths,
                              bit_lengths_size,
                              &huffman_tree_size,
                              &huffman_tree[0],
                              &huffman_tree_extra_bits[0]);
  memset(&huffman_tree_histogram[0], 0, sizeof(huffman_tree_histogram));
  for (i = 0; i < huffman_tree_size; ++i) {
    ++huffman_tree_histogram[huffman_tree[i]];
  }
  memset(&code_length_bitdepth[0], 0, sizeof(code_length_bitdepth));
  memset(&code_length_bitdepth_symbols[0], 0,
         sizeof(code_length_bitdepth_symbols));
  ok = ok &&
      CreateHuffmanTree(&huffman_tree_histogram[0], CODE_LENGTH_CODES,
                        7, &code_length_bitdepth[0]);
  if (!ok) {
    goto exit_label;
  }
  ConvertBitDepthsToSymbols(&code_length_bitdepth[0], CODE_LENGTH_CODES,
                            &code_length_bitdepth_symbols[0]);
  StoreHuffmanTreeOfHuffmanTreeToBitMask(code_length_bitdepth,
                                         bw);
  ClearHuffmanTreeIfOnlyOneSymbol(CODE_LENGTH_CODES,
                                  &code_length_bitdepth[0],
                                  &code_length_bitdepth_symbols[0]);
  {
    int num_trailing_zeros = 0;
    int trailing_zero_bits = 0;
    int trimmed_length;
    int write_length;
    int length;
    for (i = huffman_tree_size; i > 0; --i) {
      int ix = huffman_tree[i - 1];
      if (ix == 0 || ix == 17 || ix == 18) {
        ++num_trailing_zeros;
        trailing_zero_bits += code_length_bitdepth[ix];
        if (ix == 17) trailing_zero_bits += 3;
        if (ix == 18) trailing_zero_bits += 7;
      } else {
        break;
      }
    }
    trimmed_length = huffman_tree_size - num_trailing_zeros;
    write_length = (trimmed_length > 1 && trailing_zero_bits > 12);
    length = write_length ? trimmed_length : huffman_tree_size;
    WriteBits(1, write_length, bw);
    if (write_length) {
      const int nbits = BitsLog2Ceiling(trimmed_length - 1);
      const int nbitpairs = nbits == 0 ? 1 : (nbits + 1) / 2;
      WriteBits(3, nbitpairs - 1, bw);
      WriteBits(nbitpairs * 2, trimmed_length - 2, bw);
    }
    StoreHuffmanTreeToBitMask(&huffman_tree[0],
                              &huffman_tree_extra_bits[0],
                              length,
                              &code_length_bitdepth[0],
                              &code_length_bitdepth_symbols[0],
                              bw);
  }
exit_label:
  if (huffman_tree) free(huffman_tree);
  if (huffman_tree_extra_bits) free(huffman_tree_extra_bits);
  return ok;
}

static inline int MetaSize(int size, int bits) {
  return (size + (1 << bits) - 1) >> bits;
}

static void StoreImageToBitMask(
    int xsize, int histo_bits, int palette_bits,
    const PixOrCopy *literals, int literals_size,
    const uint32_t *histogram_symbols,
    uint8_t** const bitdepths, uint16_t** const bit_symbols,
    BitWriter* const bw) {
  const int histo_xsize = histo_bits ? MetaSize(xsize, histo_bits) : 1;
  // x and y trace the position in the image.
  int x = 0;
  int y = 0;
  int i;
  for (i = 0; i < literals_size; ++i) {
    const PixOrCopy v = literals[i];
    int histogram_ix = histogram_symbols[histo_bits ?
                                        (y >> histo_bits) * histo_xsize +
                                        (x >> histo_bits) : 0];
    if (PixOrCopy_IsPaletteIx(&v)) {
      const int code = PixOrCopy_PaletteIx(&v);
      int literal_ix = 256 + code;
      WriteBits(bitdepths[5 * histogram_ix][literal_ix],
                bit_symbols[5 * histogram_ix][literal_ix], bw);
    } else if (PixOrCopy_IsLiteral(&v)) {
      static const int order[] = {1, 2, 0, 3};
      int k;
      for (k = 0; k < 4; ++k) {
        const int code = PixOrCopy_Literal(&v, order[k]);
        WriteBits(bitdepths[5 * histogram_ix + k][code],
                  bit_symbols[5 * histogram_ix + k][code], bw);
      }
    } else {
      int code;
      int n_bits;
      int bits;
      int distance;
      int len_ix;
      PixOrCopy_LengthCodeAndBits(&v, &code, &n_bits, &bits);
      len_ix = 256 + (1 << palette_bits) + code;
      WriteBits(bitdepths[5 * histogram_ix][len_ix],
                bit_symbols[5 * histogram_ix][len_ix], bw);
      WriteBits(n_bits, bits, bw);

      distance = PixOrCopy_Distance(&v);
      PrefixEncode(distance, &code, &n_bits, &bits);
      WriteBits(bitdepths[5 * histogram_ix + 4][code],
                bit_symbols[5 * histogram_ix + 4][code], bw);
      WriteBits(n_bits, bits, bw);
    }
    x += PixOrCopy_Length(&v);
    while (x >= xsize) {
      x -= xsize;
      ++y;
    }
  }
}

static void ShiftHistogramImage(int image_size, uint32_t* image) {
  int i;
  for (i = 0; i < image_size; ++i) {
    image[i] <<= 8;
    image[i] |= 0xff000000;
  }
}

static int Uint32Order(const void *p1, const void *p2) {
  uint32_t a = *(const uint32_t *)p1;
  uint32_t b = *(const uint32_t *)p2;
  if (a < b) {
    return -1;
  }
  if (a == b) {
    return 0;
  }
  return 1;
}

static int CreatePalette256(int n, const uint32_t *argb,
                            int* palette_size,
                            uint32_t* palette) {
  int i;
  const int max_palette_size = 256;
  int current_size = 0;
  uint8_t in_use[1024];
  uint32_t colors[1024];
  static const uint32_t kHashMul = 0x1e35a7bd;
  memset(&in_use[0], 0, sizeof(in_use));
  for (i = 0; i < n; ++i) {
    int addr;
    if (i != 0 && argb[i - 1] == argb[i]) {
      continue;
    }
    addr = (kHashMul * argb[i]) >> 22;
    for(;;) {
      if (!in_use[addr]) {
        colors[addr] = argb[i];
        in_use[addr] = 1;
        ++current_size;
        if (current_size > max_palette_size) {
          return 0;
        }
        break;
      } else if (colors[addr] == argb[i]) {
        // The color is already there.
        break;
      } else {
        // Some other color sits there.
        // Do linear conflict resolution.
        ++addr;
        addr &= 0x3ff;
      }
    }
  }
  *palette_size = 0;
  for (i = 0; i < (int)sizeof(in_use); ++i) {
    if (in_use[i]) {
      palette[*palette_size] = colors[i];
      ++(*palette_size);
    }
  }
  qsort(&palette[0], *palette_size, sizeof(palette[0]), Uint32Order);
  return 1;
}

// Bundles multiple (2, 4 or 8) pixels into a single pixel.
// Returns the new xsize.
static void BundlePixels(int xsize, int ysize, int xbits,
                         uint32_t *argb, int *new_xsize) {
  const int bit_depth = 1 << (3 - xbits);
  int xs = MetaSize(xsize, xbits);
  int y;
  int x;
  for (y = 0; y < ysize; ++y) {
    uint32_t code;
    for (x = 0; x < xsize; ++x) {
      const int xsub = x & ((1 << xbits) - 1);
      if (xsub == 0) {
        code = 0;
      }
      code |= (argb[y * xsize + x] & 0xff00) << (bit_depth * xsub);
      argb[y * xs + (x >> xbits)] = 0xff000000 | code;
    }
  }
  *new_xsize = xs;
}

static void DeleteHistograms(int size, Histogram** histograms) {
  int i;
  for (i = 0; i < size; ++i) {
    free(histograms[i]);
  }
}

static int GetBackwardReferences(int xsize, int ysize,
                                 const uint32_t* argb,
                                 int quality, int use_palette,
                                 int palette_bits, int use_2d_locality,
                                 int *backward_refs_size,
                                 PixOrCopy** backward_refs) {
  // Backward Reference using LZ77.
  int lz77_is_useful;
  int backward_refs_rle_only_size;
  Histogram *histo_rle = (Histogram *)malloc(sizeof(*histo_rle));
  Histogram *histo_lz77 = (Histogram *)malloc(sizeof(*histo_lz77));
  PixOrCopy *backward_refs_lz77 = (PixOrCopy *)
      malloc(xsize * ysize * sizeof(*backward_refs_lz77));
  PixOrCopy *backward_refs_rle_only = (PixOrCopy *)
      malloc(xsize * ysize * sizeof(*backward_refs_rle_only));
  int backward_refs_lz77_size;
  int ok = 1;
  *backward_refs = NULL;
  if (histo_rle == NULL || histo_lz77 == NULL ||
      backward_refs_lz77 == NULL || backward_refs_rle_only == NULL) {
    ok = 0;
    goto exit_label;
  }
  BackwardReferencesHashChain(xsize, ysize, use_palette, &argb[0], palette_bits,
                              &backward_refs_lz77[0], &backward_refs_lz77_size);
  Histogram_Init(histo_lz77, palette_bits);
  Histogram_Build(histo_lz77, &backward_refs_lz77[0], backward_refs_lz77_size);

  // Backward Reference using RLE only.
  BackwardReferencesRle(xsize, ysize, &argb[0], &backward_refs_rle_only[0],
                        &backward_refs_rle_only_size);

  Histogram_Init(histo_rle, palette_bits);
  Histogram_Build(histo_rle,
                  &backward_refs_rle_only[0], backward_refs_rle_only_size);

  // Check if LZ77 is useful.
  lz77_is_useful =
      (Histogram_EstimateBits(histo_rle) > Histogram_EstimateBits(histo_lz77));

  // Choose appropriate backward reference.
  if (quality >= 50 && lz77_is_useful) {
    const int recursion_level = (xsize * ysize < 320 * 200) ? 1 : 0;
    free(backward_refs_rle_only);
    free(backward_refs_lz77);
    *backward_refs = (PixOrCopy *)
        malloc(xsize * ysize * sizeof(*backward_refs));
    BackwardReferencesTraceBackwards(xsize, ysize, recursion_level, use_palette,
                                     &argb[0], palette_bits, *backward_refs,
                                     backward_refs_size);
  } else {
    if (lz77_is_useful) {
      *backward_refs = backward_refs_lz77;
      *backward_refs_size = backward_refs_lz77_size;
      free(backward_refs_rle_only);
    } else {
      *backward_refs = backward_refs_rle_only;
      *backward_refs_size = backward_refs_rle_only_size;
      free(backward_refs_lz77);
    }
  }

  // Verify.
  VERIFY(VerifyBackwardReferences(&argb[0], xsize, ysize, palette_bits,
                                  *backward_refs, *backward_refs_size));

  if (use_2d_locality) {
    // Use backward reference with 2D locality.
    BackwardReferences2DLocality(xsize, *backward_refs_size, *backward_refs);
  }
exit_label:
  if (histo_rle) free(histo_rle);
  if (histo_lz77) free(histo_lz77);
  if (!ok) {
    if (*backward_refs) free(*backward_refs);
    *backward_refs = NULL;
  }
  return ok;
}

static void GetHistImageSymbols(int xsize, int ysize,
                                PixOrCopy* backward_refs,
                                int backward_refs_size,
                                int quality, int histogram_bits,
                                int palette_bits,
                                int *histogram_image_size,
                                Histogram*** histogram_image,
                                uint32_t* histogram_symbols) {
  // Build histogram image.
  const int max_refinement_iters = 1;
  Histogram** histogram_image_raw;
  int histogram_image_raw_size;
  int i;
  BuildHistogramImage(xsize, ysize, histogram_bits, palette_bits,
                      backward_refs, backward_refs_size,
                      &histogram_image_raw,
                      &histogram_image_raw_size);
  // Collapse similar histograms.
  CombineHistogramImage(histogram_image_raw, histogram_image_raw_size,
                        quality, palette_bits,
                        histogram_image,
                        histogram_image_size);
  // Refine histogram image.
  for (i = 0; i < histogram_image_raw_size; ++i) {
    histogram_symbols[i] = -1;
  }
  for (i = 0; i < max_refinement_iters; ++i) {
    RefineHistogramImage(histogram_image_raw, histogram_image_raw_size,
                         histogram_symbols,
                         *histogram_image_size,
                         *histogram_image);
    if (quality < 30) {
      break;
    }
  }
  DeleteHistograms(histogram_image_raw_size, histogram_image_raw);
  free(histogram_image_raw);
}

// Returns 0 when no error has occured.
static int GetHuffBitLengthsAndCodes(
    int histogram_image_size,
    Histogram** histogram_image,
    int use_palette,
    int **bit_length_sizes,
    uint8_t*** bit_lengths,
    uint16_t*** bit_codes) {
  int i;
  int k;
  int ok = 1;
  for (i = 0; i < histogram_image_size; ++i) {
    const int lit_len = Histogram_NumPixOrCopyCodes(histogram_image[i]);
    (*bit_length_sizes)[5 * i] = lit_len;
    (*bit_lengths)[5 * i] = (uint8_t *)calloc(lit_len, 1);
    (*bit_codes)[5 * i] = (uint16_t *)
        malloc(lit_len * sizeof(*(*bit_codes)[5 * i]));

    // For each component, optimize histogram for Huffman with RLE compression.
    ok = ok && OptimizeHuffmanForRle(lit_len, &histogram_image[i]->literal_[0]);
    if (!use_palette) {
      // Implies that palette_bits == 0,
      // and so number of palette entries = (1 << 0) = 1.
      // Optimization might have smeared population count in this single
      // palette entry, so zero it out.
      histogram_image[i]->literal_[256] = 0;
    }
    ok = ok && OptimizeHuffmanForRle(256, &histogram_image[i]->red_[0]);
    ok = ok && OptimizeHuffmanForRle(256, &histogram_image[i]->blue_[0]);
    ok = ok && OptimizeHuffmanForRle(256, &histogram_image[i]->alpha_[0]);
    ok = ok && OptimizeHuffmanForRle(DISTANCE_CODES_MAX,
                                     &histogram_image[i]->distance_[0]);

    // Create a Huffman tree (in the form of bit lengths) for each component.
    ok = ok && CreateHuffmanTree(histogram_image[i]->literal_, lit_len, 15,
                                 (*bit_lengths)[5 * i]);
    for (k = 1; k < 5; ++k) {
      int val = 256;
      if (k == 4) {
        val = DISTANCE_CODES_MAX;
      }
      (*bit_length_sizes)[5 * i + k] = val;
      (*bit_lengths)[5 * i + k] = (uint8_t *)calloc(val, 1);
      (*bit_codes)[5 * i + k] = (uint16_t *)calloc(val, sizeof(bit_codes[0]));
    }
    ok = ok && CreateHuffmanTree(histogram_image[i]->red_, 256, 15,
                                 (*bit_lengths)[5 * i + 1]) &&
        CreateHuffmanTree(histogram_image[i]->blue_, 256, 15,
                          (*bit_lengths)[5 * i + 2]) &&
        CreateHuffmanTree(histogram_image[i]->alpha_, 256, 15,
                          (*bit_lengths)[5 * i + 3]) &&
        CreateHuffmanTree(histogram_image[i]->distance_,
                          DISTANCE_CODES_MAX, 15,
                          (*bit_lengths)[5 * i + 4]);
    // Create the actual bit codes for the bit lengths.
    for (k = 0; k < 5; ++k) {
      int ix = 5 * i + k;
      ConvertBitDepthsToSymbols((*bit_lengths)[ix], (*bit_length_sizes)[ix],
                                (*bit_codes)[ix]);
    }
  }
  return ok;
}

static void EncodeImageInternal(int xsize, int ysize,
                                const uint32_t *argb, int quality,
                                int palette_bits, int histogram_bits,
                                int use_2d_locality, BitWriter *bw) {
  int histogram_image_size;
  Histogram **histogram_image;

  int* bit_lengths_sizes;
  uint8_t** bit_lengths;
  uint16_t** bit_codes;

  int write_histogram_image;
  int i;
  const int use_palette = palette_bits ? 1 : 0;

  int backward_refs_size;
  PixOrCopy* backward_refs;
  const int histogram_image_xysize = MetaSize(xsize, histogram_bits) *
      MetaSize(ysize, histogram_bits);
  uint32_t *histogram_symbols = (uint32_t *)
      calloc(histogram_image_xysize, sizeof(*histogram_symbols));
  // Calculate backward references from ARGB image.
  GetBackwardReferences(xsize, ysize, argb, quality, use_palette, palette_bits,
                        use_2d_locality, &backward_refs_size, &backward_refs);

  // Build histogram image & symbols from backward references.
  GetHistImageSymbols(xsize, ysize, backward_refs, backward_refs_size,
                      quality, histogram_bits,
                      palette_bits,
                      &histogram_image_size,
                      &histogram_image,
                      histogram_symbols);

  // Create Huffman bit lengths & codes for each histogram image.
  bit_lengths_sizes = (int *)calloc(5 * histogram_image_size,
                                    sizeof(*bit_lengths_sizes));
  bit_lengths = (uint8_t**)calloc(5 * histogram_image_size,
                                  sizeof(*bit_lengths));
  bit_codes = (uint16_t**)calloc(5 * histogram_image_size,
                                 sizeof(*bit_codes));
  GetHuffBitLengthsAndCodes(histogram_image_size, histogram_image,
                            use_palette, &bit_lengths_sizes,
                            &bit_lengths, &bit_codes);

  // No transforms.
  WriteBits(1, 0, bw);

  // Huffman image + meta huffman.
  write_histogram_image = (histogram_image_size > 1);
  WriteBits(1, write_histogram_image, bw);
  if (write_histogram_image) {
    int nbits;
    int image_size_bits;
    int num_histograms;
    uint32_t* histogram_argb = (uint32_t*)
        malloc(histogram_image_xysize * sizeof(*histogram_argb));
    memcpy(histogram_argb, histogram_symbols,
           histogram_image_xysize * sizeof(*histogram_argb));

    ShiftHistogramImage(histogram_image_xysize, &histogram_argb[0]);
    WriteBits(4, histogram_bits, bw);
    EncodeImageInternal(MetaSize(xsize, histogram_bits),
                        MetaSize(ysize, histogram_bits),
                        &histogram_argb[0],
                        quality,
                        0,
                        0,      // no histogram bits
                        1,   // use_2d_locality
                        bw);
    image_size_bits = BitsLog2Ceiling(histogram_image_size - 1);
    WriteBits(4, image_size_bits, bw);
    WriteBits(image_size_bits, histogram_image_size - 2, bw);
    num_histograms = 5 * histogram_image_size;
    nbits = BitsLog2Ceiling(num_histograms);
    WriteBits(4, nbits, bw);
    for (i = 0; i < num_histograms; ++i) {
      WriteBits(nbits, i, bw);
    }
  }

  // Palette parameters.
  WriteBits(1, use_palette, bw);
  if (use_palette) {
    WriteBits(4, kRowHasherXSubsampling, bw);
    WriteBits(4, palette_bits, bw);
  }

  // Store Huffman codes.
  for (i = 0; i < histogram_image_size; ++i) {
    int k;
    int literal_lengths_size;
    uint8_t* literal_lengths;
    PackLiteralBitLengths(&bit_lengths[5 * i][0], palette_bits, use_palette,
                          &literal_lengths_size, &literal_lengths);
    StoreHuffmanCode(literal_lengths, literal_lengths_size, bw);
    free(literal_lengths);
    for (k = 1; k < 5; ++k) {
      StoreHuffmanCode(&bit_lengths[5 * i + k][0],
                       bit_lengths_sizes[5 * i + k], bw);
    }
  }
  // free combined histograms
  DeleteHistograms(histogram_image_size, &histogram_image[0]);
  free(histogram_image);

  // Emit no bits if there is only one symbol in the histogram.
  // This gives better compression for some images.
  for (i = 0; i < 5 * histogram_image_size; ++i) {
    ClearHuffmanTreeIfOnlyOneSymbol(bit_lengths_sizes[i], &bit_lengths[i][0],
                                    &bit_codes[i][0]);
  }
  // Store actual literals
  StoreImageToBitMask(xsize, histogram_bits, palette_bits,
                      backward_refs, backward_refs_size,
                      histogram_symbols, bit_lengths, bit_codes, bw);
  for (i = 0; i < 5 * histogram_image_size; ++i) {
    free(bit_lengths[i]);
    free(bit_codes[i]);
  }
  free(bit_lengths_sizes);
  free(bit_lengths);
  free(bit_codes);
}

static inline void WriteImageSize(uint32_t xsize, uint32_t ysize,
                                  BitWriter *bw) {
  --xsize;
  --ysize;
  VERIFY(xsize < 0x4000 && ysize < 0x4000);
  WriteBits(14, xsize, bw);
  WriteBits(14, ysize, bw);
}

// Returns 1 on success.
int EncodeWebpLLImage(int xsize, int ysize, const uint32_t *argb_orig,
                      EncodingStrategy *strategy,
                      size_t *num_bytes, uint8_t **bytes) {
  int i;
  const int quality = strategy->quality;
  //  const int use_lz77 = strategy->use_lz77;
  int palette_bits = strategy->palette_bits;
  const int use_small_palette = strategy->use_small_palette;
  const int predict = strategy->predict;
  const int predict_bits = strategy->predict_bits;
  const int histogram_bits = strategy->histogram_bits;
  const int cross_color_transform = strategy->cross_color_transform;
  const int cross_color_transform_bits = strategy->cross_color_transform_bits;
  int use_emerging_palette = 1;
  uint32_t* argb = (uint32_t *)malloc(xsize * ysize * sizeof(*argb));
  const int kEstmatedEncodeSize = 0.5 * xsize * ysize;  // 4 bpp.
  int palette_size;
  uint32_t palette[256];

  BitWriter bw;
  if (!BitWriterInit(&bw, kEstmatedEncodeSize)) return 0;

  memcpy(argb, argb_orig, xsize * ysize * sizeof(*argb));

  // Write image size.
  WriteImageSize(xsize, ysize, &bw);

  if (!use_small_palette) {
    Histogram *after;
    Histogram *before = (Histogram *)malloc(sizeof(*before));
    Histogram_Init(before, 1);
    for (i = 0; i < xsize * ysize; ++i) {
      Histogram_AddSinglePixOrCopy(before, PixOrCopy_CreateLiteral(argb[i]));
    }
    SubtractGreenFromBlueAndRed(xsize * ysize, &argb[0]);
    after = (Histogram *)malloc(sizeof(*after));
    Histogram_Init(after, 1);
    for (i = 0; i < xsize * ysize; ++i) {
      Histogram_AddSinglePixOrCopy(after, PixOrCopy_CreateLiteral(argb[i]));
    }
    if (Histogram_EstimateBits(after) < Histogram_EstimateBits(before)) {
      WriteBits(1, 1, &bw);
      WriteBits(3, 2, &bw);
    } else {
      // Undo subtract green from blue and red -- rewrite with original data.
      memcpy(argb, argb_orig, xsize * ysize * sizeof(argb[0]));
    }
    free(before);
    free(after);
  }

  if (use_small_palette &&
      CreatePalette256(xsize * ysize, &argb[0], &palette_size, &palette[0])) {
    for (i = 0; i < xsize * ysize; ++i) {
      int k;
      for (k = 0; k < palette_size; ++k) {
        if (argb[i] == palette[k]) {
          argb[i] = 0xff000000 | (k << 8);
          break;
        }
      }
    }
    WriteBits(1, 1, &bw);
    WriteBits(3, 3, &bw);
    WriteBits(8, palette_size - 1, &bw);
    for (i = palette_size - 1; i >= 1; --i) {
      palette[i] = Subtract(palette[i], palette[i - 1]);
    }
    EncodeImageInternal(palette_size, 1, &palette[0], quality, 0, 0, 1, &bw);
    use_emerging_palette = 0;
    if (palette_size <= 16) {
      int xbits = 1;
      if (palette_size <= 2) {
        xbits = 3;
      } else if (palette_size <= 4) {
        xbits = 2;
      }
      BundlePixels(xsize, ysize, xbits, argb, &xsize);
      WriteBits(1, 1, &bw);
      WriteBits(3, 4, &bw);
      WriteBits(2, xbits, &bw);
    }
  }

  if (predict) {
    const int predictor_image_size =
        MetaSize(xsize, predict_bits) *
        MetaSize(ysize, predict_bits);
    uint32_t* predictor_image = (uint32_t*)
        malloc(predictor_image_size * sizeof(*predictor_image));
    uint32_t* from_argb = (uint32_t*)malloc(xsize * ysize * sizeof(*from_argb));
    memcpy(from_argb, &argb[0], xsize * ysize * sizeof(*argb));
    PredictorImage(xsize, ysize, predict_bits,
                   &from_argb[0],
                   &argb[0],
                   &predictor_image[0]);
    WriteBits(1, 1, &bw);
    WriteBits(3, 0, &bw);
    WriteBits(4, predict_bits, &bw);
    EncodeImageInternal(MetaSize(xsize, predict_bits),
                        MetaSize(ysize, predict_bits),
                        &predictor_image[0],
                        quality,
                        0,
                        0,      // no histogram bits
                        1,   // use_2d_locality
                        &bw);
    free(from_argb);
    free(predictor_image);
  }

  if (cross_color_transform) {
    const int color_space_image_size =
        MetaSize(xsize, cross_color_transform_bits) *
        MetaSize(ysize, cross_color_transform_bits);
    uint32_t* color_space_image = (uint32_t *)
        malloc(color_space_image_size * sizeof(*color_space_image));
    ColorSpaceTransform(xsize, ysize, cross_color_transform_bits,
                        &argb[0],
                        quality,
                        &argb[0],
                        color_space_image);
    WriteBits(1, 1, &bw);
    WriteBits(3, 1, &bw);
    WriteBits(4, cross_color_transform_bits, &bw);
    EncodeImageInternal(MetaSize(xsize, cross_color_transform_bits),
                        MetaSize(ysize, cross_color_transform_bits),
                        color_space_image,
                        quality,
                        0,
                        0,      // no histogram bits
                        1,   // use_2d_locality
                        &bw);
    free(color_space_image);
  }

  palette_bits = 0;
  if (use_emerging_palette) {
    palette_bits = (quality == 0) ?
        7 :
        CalculateEstimateForPaletteSize(&argb[0], xsize, ysize);
  }

  EncodeImageInternal(xsize,
                      ysize,
                      &argb[0],
                      quality,
                      palette_bits,
                      histogram_bits,
                      1,   // use_2d_locality
                      &bw);

  free(argb);

  {
    const size_t webpll_size = BitWriterNumBytes(&bw);
    uint8_t *webpll_data = BitWriterFinish(&bw);
    *num_bytes = HEADER_SIZE + SIGNATURE_SIZE + webpll_size;
    // TODO(vikasa): This memory allocation can be avoided, if RIFF header is
    // writen to BitWriter in the begining and BitWriter's buffer can be owned
    // and passed back instead.
    *bytes = (uint8_t *)malloc(*num_bytes);
    if (*bytes == NULL) return 0;
    PutRiffHeader(*bytes, webpll_size);
    memcpy(*bytes + HEADER_SIZE + SIGNATURE_SIZE, webpll_data, webpll_size);
    BitWriterDestroy(&bw);
  }
  return 1;
}
