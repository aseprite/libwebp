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

#include "webp_bit_stream.h"

#include <string.h>

#include <vector>

#include "../common/integral_types.h"
#include "../common/predictor.h"
#include "./backward_distance.h"
#include "./backward_references.h"
#include "./bit_writer.h"
#include "./encode.h"
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
bool ValuesShouldBeCollapsedToStrideAverage(int stride, int a, int b) {
  return abs(a - b) < 4;
}

// Change the population counts in a way that the consequent
// Hufmann tree compression, especially its rle-part will be more
// likely to compress this data more efficiently.
//
// length containts the size of the histogram.
// data contains the population counts.
void OptimizeHuffmanForRle(int length, int *counts) {
  // Let's make the Huffman code more compatible with rle encoding.
  for (; length >= 0; --length) {
    if (length == 0) {
      return;
    }
    if (counts[length - 1] != 0) {
      // Now counts[0..length - 1] does not have trailing zeros.
      break;
    }
  }
  // 2) Let's mark all population counts that already can be encoded
  // with an rle code.
  char* good_for_rle = (char *)malloc(length);
  memset(good_for_rle, 0, length);
  {
    // Let's not spoil any of the existing good rle codes.
    // Mark any seq of 0's that is longer as 5 as a good_for_rle.
    // Mark any seq of non-0's that is longer as 7 as a good_for_rle.
    int symbol = counts[0];
    int stride = 0;
    for (int i = 0; i < length + 1; ++i) {
      if (i == length || counts[i] != symbol) {
        if ((symbol == 0 && stride >= 5) ||
            (symbol != 0 && stride >= 7)) {
          for (int k = 0; k < stride; ++k) {
            good_for_rle[i - k - 1] = true;
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
  int stride = 0;
  int limit = counts[0];
  int sum = 0;
  for (int i = 0; i < length + 1; ++i) {
    if (i == length || good_for_rle[i] ||
        (i != 0 && good_for_rle[i - 1]) ||
        !ValuesShouldBeCollapsedToStrideAverage(stride, counts[i], limit)) {
      if (stride >= 4 || (stride >= 3 && sum == 0)) {
        // The stride must end, collapse what we have, if we have enough (4).
        int count = (sum + stride / 2) / stride;
        if (count < 1) {
          count = 1;
        }
        if (sum == 0) {
          // Don't make an all zeros stride to be upgraded to ones.
          count = 0;
        }
        for (int k = 0; k < stride; ++k) {
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
}

void ClearHuffmanTreeIfOnlyOneSymbol(const int num_symbols,
                                     uint8_t* lengths,
                                     std::vector<uint16_t>* symbols) {
  int count = 0;
  for (int k = 0; k < num_symbols; ++k) {
    if (lengths[k] != 0) ++count;
  }
  if (count == 1) {
    for (int k = 0; k < num_symbols; ++k) {
      lengths[k] = 0;
      (*symbols)[k] = 0;
    }
  }
}

void PackLiteralBitLengths(const std::vector<uint8_t> bit_lengths,
                           int palette_bits, bool use_palette,
                           std::vector<uint8_t>* new_lengths) {
  for (int i = 0; i < 256; ++i) {
    new_lengths->push_back(bit_lengths[i]);
  }
  if (use_palette) {
    for (int i = 0; i < (1 << palette_bits); ++i) {
      new_lengths->push_back(bit_lengths[256 + i]);
    }
  }
  for (int i = 0; i < kLengthCodes; ++i) {
    new_lengths->push_back(bit_lengths[256 + (1 << palette_bits) + i]);
  }
}

void StoreHuffmanTreeOfHuffmanTreeToBitMask(
    const uint8_t *code_length_bitdepth, BitWriter* const bw) {
  // RFC 1951 will calm you down if you are worried about this funny sequence.
  // This sequence is tuned from that, but more weighted for lower symbol count,
  // and more spiking histograms.
  static const uint8_t kStorageOrder[kCodeLengthCodes] = {
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
  for (uint32_t i = 0; i < codes_to_store; ++i) {
    WriteBits(3, code_length_bitdepth[kStorageOrder[i]], bw);
  }
}

void StoreHuffmanTreeToBitMask(
    const std::vector<uint8_t> &huffman_tree,
    const std::vector<uint8_t> &huffman_tree_extra_bits,
    const int num_symbols,
    const uint8_t *code_length_bitdepth,
    const std::vector<uint16_t> &code_length_bitdepth_symbols,
    BitWriter* const bw) {
  for (uint32_t i = 0; i < num_symbols; ++i) {
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

void StoreHuffmanCode(const std::vector<uint8_t> &bit_lengths,
                      BitWriter* const bw) {
  int count = 0;
  int symbols[2] = { 0 };
  for (int i = 0; i < bit_lengths.size(); ++i) {
    if (bit_lengths[i] != 0) {
      if (count < 2) symbols[count] = i;
      ++count;
    }
  }
  if (count <= 2) {
    // 0, 1 or 2 symbols to encode.
    WriteBits(1, 1, bw);
    if (count == 0) {
      WriteBits(4, 0, bw);
      return;
    }
    WriteBits(1, count - 1, bw);
    int num_bits = 4;
    while (symbols[count - 1] >= (1 << num_bits)) num_bits += 2;
    WriteBits(3, (num_bits - 4) / 2 + 1, bw);
    for (int i = 0; i < count; ++i) {
      WriteBits(num_bits, symbols[i], bw);
    }
    return;
  }
  WriteBits(1, 0, bw);
  std::vector<uint8_t> huffman_tree(bit_lengths.size());
  std::vector<uint8_t> huffman_tree_extra_bits(bit_lengths.size());
  int huffman_tree_size = 0;
  CreateCompressedHuffmanTree(&bit_lengths[0],
                              bit_lengths.size(),
                              &huffman_tree_size,
                              &huffman_tree[0],
                              &huffman_tree_extra_bits[0]);
  huffman_tree.resize(huffman_tree_size);
  huffman_tree_extra_bits.resize(huffman_tree_size);

  int huffman_tree_histogram[kCodeLengthCodes] = { 0 };
  for (int i = 0; i < huffman_tree.size(); ++i) {
    ++huffman_tree_histogram[huffman_tree[i]];
  }
  uint8_t code_length_bitdepth[kCodeLengthCodes] = { 0 };
  std::vector<uint16_t> code_length_bitdepth_symbols(kCodeLengthCodes);
  CreateHuffmanTree(&huffman_tree_histogram[0], kCodeLengthCodes,
                    7, &code_length_bitdepth[0]);
  ConvertBitDepthsToSymbols(&code_length_bitdepth[0], kCodeLengthCodes,
                            &code_length_bitdepth_symbols[0]);
  StoreHuffmanTreeOfHuffmanTreeToBitMask(code_length_bitdepth,
                                         bw);
  ClearHuffmanTreeIfOnlyOneSymbol(kCodeLengthCodes,
                                  code_length_bitdepth,
                                  &code_length_bitdepth_symbols);
  int num_trailing_zeros = 0;
  int trailing_zero_bits = 0;
  for (int i = huffman_tree.size(); i > 0; --i) {
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
  const int trimmed_length = huffman_tree.size() - num_trailing_zeros;
  const bool write_length = (trimmed_length > 1 && trailing_zero_bits > 12);
  const int length = write_length ? trimmed_length : huffman_tree.size();
  WriteBits(1, write_length, bw);
  if (write_length) {
    const int nbits = BitsLog2Ceiling(trimmed_length - 1);
    const int nbitpairs = nbits == 0 ? 1 : (nbits + 1) / 2;
    WriteBits(3, nbitpairs - 1, bw);
    WriteBits(nbitpairs * 2, trimmed_length - 2, bw);
  }
  StoreHuffmanTreeToBitMask(huffman_tree,
                            huffman_tree_extra_bits,
                            length,
                            &code_length_bitdepth[0],
                            code_length_bitdepth_symbols,
                            bw);
}

void StoreImageToBitMask(
    const int xsize,
    const int ysize,
    const int histobits,
    const int palette_bits,
    const LiteralOrCopy *literal,
    const int n_literal_and_length,
    const std::vector<uint32_t> &histogram_symbol,
    const std::vector< std::vector<uint8_t> > &bitdepth,
    const std::vector< std::vector<uint16_t> > &bit_symbols,
    BitWriter *bw) {
  int histo_xsize = histobits ? (xsize + (1 << histobits) - 1) >> histobits : 1;
  // x and y trace the position in the image.
  int x = 0;
  int y = 0;
  for (int i = 0; i < n_literal_and_length; ++i) {
    const LiteralOrCopy &v = literal[i];
    int histogram_ix = histogram_symbol[histobits ?
                                        (y >> histobits) * histo_xsize +
                                        (x >> histobits) : 0];
    if (v.IsPaletteIx()) {
      const int code = v.PaletteIx();
      int literal_ix = 256 + code;
      WriteBits(bitdepth[5 * histogram_ix][literal_ix],
                bit_symbols[5 * histogram_ix][literal_ix], bw);
    } else if (v.IsLiteral()) {
      int order[] = {1, 2, 0, 3};
      for (int i = 0; i < 4; ++i) {
        const int code = v.Literal(order[i]);
        WriteBits(bitdepth[5 * histogram_ix + i][code],
                  bit_symbols[5 * histogram_ix + i][code], bw);
      }
    } else {
      int code;
      int n_bits;
      int bits;
      v.LengthCodeAndBits(&code, &n_bits, &bits);
      int len_ix = 256 + (1 << palette_bits) + code;
      WriteBits(bitdepth[5 * histogram_ix][len_ix],
                bit_symbols[5 * histogram_ix][len_ix], bw);
      WriteBits(n_bits, bits, bw);

      const int distance = v.Distance();
      BackwardDistance::Encode(distance, &code, &n_bits, &bits);
      WriteBits(bitdepth[5 * histogram_ix + 4][code],
                bit_symbols[5 * histogram_ix + 4][code], bw);
      WriteBits(n_bits, bits, bw);
    }
    x += v.Length();
    while (x >= xsize) {
      x -= xsize;
      ++y;
    }
  }
}

void ShiftHistogramImage(std::vector<uint32_t>* image) {
  for (int i = 0; i < image->size(); ++i) {
    (*image)[i] <<= 8;
    (*image)[i] |= 0xff000000;
  }
}

int NumberOfUsedRBAComponents(const std::vector<uint32_t>& argb) {
  int num = 0;
  for (int i = 0; i < argb.size() && num < 3; ++i) {
    if (((argb[i] >> 24) & 0xff) != 0xff) return 3;
    if ((argb[i] & 0xff) != 0 && num < 2) num = 2;
    if (((argb[i] >> 16) & 0xff) != 0 && num < 1) num = 1;
  }
  return num;
}

int MetaSize(int size, int bits) {
  return (size + (1 << bits) - 1) >> bits;
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

bool CreatePalette256(int n, const uint32_t *argb,
                      std::vector<uint32_t>* palette) {
  int i;
  const int max_palette_size = 256;
  int current_size = 0;
  uint8_t in_use[1024];
  uint32_t colors[1024];
  static const uint32_t kHashMul = 0x1e35a7bd;
  memset(&in_use[0], 0, sizeof(in_use));
  for (i = 0; i < n; ++i) {
    if (i != 0 && argb[i - 1] == argb[i]) {
      continue;
    }
    int addr = (kHashMul * argb[i]) >> 22;
    for(;;) {
      if (!in_use[addr]) {
        colors[addr] = argb[i];
        in_use[addr] = 1;
        ++current_size;
        if (current_size > max_palette_size) {
          return false;
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
  for (i = 0; i < sizeof(in_use); ++i) {
    if (in_use[i]) {
      palette->push_back(colors[i]);
    }
  }
  qsort(&(*palette)[0], palette->size(), sizeof((*palette)[0]), Uint32Order);
  return true;
}

// Bundles multiple (2, 4 or 8) pixels into a single pixel.
// Returns the new xsize.
int BundlePixels(int xsize, int ysize, int xbits,
                 const std::vector<uint32_t>& from_argb,
                 std::vector<uint32_t>* to_argb) {
  const int bit_depth = 1 << (3 - xbits);
  const int xs = MetaSize(xsize, xbits);
  VERIFY(xsize * ysize == from_argb.size());
  to_argb->resize(xs * ysize);
  for (int y = 0; y < ysize; ++y) {
    uint32_t code;
    for (int x = 0; x < xsize; ++x) {
      const int xsub = x & ((1 << xbits) - 1);
      if (xsub == 0) {
        code = 0;
      }
      const uint32_t green = from_argb[y * xsize + x] & 0xff00;
      code |= green << (bit_depth * xsub);
      (*to_argb)[y * xs + (x >> xbits)] = 0xff000000 | code;
    }
  }
  return xs;
}

static void DeleteHistograms(std::vector<Histogram *> &histograms) {
  int i;
  for (i = 0; i < histograms.size(); ++i) {
    delete histograms[i];
  }
}

static void GetBackwardReferences(int xsize, int ysize,
                                  const std::vector<uint32_t>& argb,
                                  int quality, int use_palette,
                                  int palette_bits, bool use_2d_locality,
                                  std::vector<LiteralOrCopy>& backward_refs) {
  // Backward Reference using LZ77.
  std::vector<LiteralOrCopy> backward_refs_lz77;
  BackwardReferencesHashChain(xsize, ysize, use_palette, &argb[0], palette_bits,
                              &backward_refs_lz77);
  Histogram *histo_lz77 = new Histogram(palette_bits);
  histo_lz77->Build(&backward_refs_lz77[0], backward_refs_lz77.size());

  // Backward Reference using RLE only.
  std::vector<LiteralOrCopy> backward_refs_rle_only;
  BackwardReferencesRle(xsize, ysize, &argb[0], &backward_refs_rle_only);
  Histogram *histo_rle = new Histogram(palette_bits);
  histo_rle->Build(&backward_refs_rle_only[0], backward_refs_rle_only.size());

  // Check if LZ77 is useful.
  const bool lz77_is_useful =
      (histo_rle->EstimateBits() > histo_lz77->EstimateBits());
  if (lz77_is_useful) {
    backward_refs_rle_only.clear();
  } else {
    backward_refs_lz77.clear();
  }
  delete histo_rle;
  delete histo_lz77;

  // Choose appropriate backward reference.
  if (quality >= 50 && lz77_is_useful) {
    const int recursion_level = (xsize * ysize < 320 * 200) ? 1 : 0;
    BackwardReferencesTraceBackwards(xsize, ysize, recursion_level, use_palette,
                                     &argb[0], palette_bits, &backward_refs);
  } else {
    if (lz77_is_useful) {
      backward_refs.swap(backward_refs_lz77);
    } else {
      backward_refs.swap(backward_refs_rle_only);
    }
  }

  // Verify.
  VERIFY(VerifyBackwardReferences(&argb[0], xsize, ysize, palette_bits,
                                  &backward_refs[0], backward_refs.size()));

  if (use_2d_locality) {
    // Use backward reference with 2D locality.
    BackwardReferences2DLocality(xsize, ysize, backward_refs.size(),
                                 &backward_refs[0]);
  }
}

static void GetHistImageSymbols(int xsize, int ysize,
                                const std::vector<LiteralOrCopy>& backward_refs,
                                int quality, int histogram_bits,
                                int palette_bits, bool use_2d_locality,
                                std::vector<Histogram*>& histogram_image,
                                std::vector<uint32_t>& histogram_symbols) {
  // Build histogram image.
  Histogram** histogram_image_raw;
  int histogram_image_raw_size;
  BuildHistogramImage(xsize, ysize, histogram_bits, palette_bits,
                      &backward_refs[0], backward_refs.size(),
                      &histogram_image_raw,
                      &histogram_image_raw_size);
  // Collapse similar histograms.
  histogram_symbols.clear();
  histogram_symbols.resize(histogram_image_raw_size, -1);
  {
    // TODO(jyrki): remove these once this function does not use a vector<>.
    Histogram **no_vec;
    int count;
    CombineHistogramImage(histogram_image_raw, histogram_image_raw_size,
                          quality, palette_bits,
                          &no_vec,
                          &count);
    histogram_image.resize(count);
    for (int i = 0; i < count; ++i) {
      histogram_image[i] = no_vec[i];
    }
    free(no_vec);
  }
  // Refine histogram image.
  const int max_refinement_iters = 1;
  histogram_symbols.resize(histogram_image_raw_size);
  for (int iter = 0; iter < max_refinement_iters; ++iter) {
    RefineHistogramImage(histogram_image_raw, histogram_image_raw_size,
                         &histogram_symbols[0],
                         histogram_image.size(),
                         &histogram_image[0]);

    if (quality < 30) {
      break;
    }
  }
  // TODO(jyrki): Restore the use of DeleteHistograms(histogram_image_raw);
  {
    int i;
    for (i = 0; i < histogram_image_raw_size; ++i) {
      delete histogram_image_raw[i];
    }
    free(histogram_image_raw);
  }
}

static void GetHuffBitLengthsAndCodes(
    const std::vector<Histogram*>& histogram_image,
    int use_palette,
    std::vector< std::vector<uint8_t> >& bit_lengths,
    std::vector< std::vector<uint16_t> >& bit_codes) {
  for (int i = 0; i < histogram_image.size(); ++i) {
    bit_lengths[5 * i].resize(histogram_image[i]->NumLiteralOrCopyCodes());

    // For each component, optimize histogram for Huffman with RLE compression.
    OptimizeHuffmanForRle(histogram_image[i]->NumLiteralOrCopyCodes(),
                          &histogram_image[i]->literal_[0]);
    if (!use_palette) {
      // Implies that palette_bits == 0,
      // and so number of palette entries = (1 << 0) = 1.
      // Optimization might have smeared population count in this single
      // palette entry, so zero it out.
      histogram_image[i]->literal_[256] = 0;
    }
    OptimizeHuffmanForRle(256, &histogram_image[i]->red_[0]);
    OptimizeHuffmanForRle(256, &histogram_image[i]->blue_[0]);
    OptimizeHuffmanForRle(256, &histogram_image[i]->alpha_[0]);
    OptimizeHuffmanForRle(kDistanceCodes, &histogram_image[i]->distance_[0]);

    // Create a Huffman tree (in the form of bit lengths) for each component.
    CreateHuffmanTree(histogram_image[i]->literal_,
                      histogram_image[i]->NumLiteralOrCopyCodes(), 15,
                      &bit_lengths[5 * i][0]);
    bit_lengths[5 * i + 1].resize(256);
    CreateHuffmanTree(histogram_image[i]->red_, 256, 15,
                      &bit_lengths[5 * i + 1][0]);

    bit_lengths[5 * i + 2].resize(256);
    CreateHuffmanTree(histogram_image[i]->blue_, 256, 15,
                      &bit_lengths[5 * i + 2][0]);
    bit_lengths[5 * i + 3].resize(256);
    CreateHuffmanTree(histogram_image[i]->alpha_, 256, 15,
                      &bit_lengths[5 * i + 3][0]);
    bit_lengths[5 * i + 4].resize(kDistanceCodes);
    CreateHuffmanTree(histogram_image[i]->distance_, kDistanceCodes, 15,
                      &bit_lengths[5 * i + 4][0]);
  }

  // Create the actual bit codes for the bit lengths.
  bit_codes.resize(bit_lengths.size());
  for (int i = 0; i < bit_codes.size(); ++i) {
    bit_codes[i].resize(bit_lengths[i].size());
    ConvertBitDepthsToSymbols(&bit_lengths[i][0], bit_lengths[i].size(),
                              &bit_codes[i][0]);
  }
}

static void EncodeImageInternal(int xsize, int ysize,
                                const std::vector<uint32_t>& argb, int quality,
                                int palette_bits, int histogram_bits,
                                bool use_2d_locality, BitWriter *bw) {
  const int use_palette = palette_bits ? 1 : 0;

  // Calculate backward references from ARGB image.
  std::vector<LiteralOrCopy> backward_refs;
  GetBackwardReferences(xsize, ysize, argb, quality, use_palette, palette_bits,
                        use_2d_locality, backward_refs);

  // Build histogram image & symbols from backward references.
  std::vector<Histogram*> histogram_image;
  std::vector<uint32_t> histogram_symbols;
  GetHistImageSymbols(xsize, ysize, backward_refs, quality, histogram_bits,
                      palette_bits, use_2d_locality, histogram_image,
                      histogram_symbols);

  // Create Huffman bit lengths & codes for each histogram image.
  std::vector< std::vector<uint8_t> > bit_lengths(5 * histogram_image.size());
  std::vector< std::vector<uint16_t> > bit_codes;
  GetHuffBitLengthsAndCodes(histogram_image, use_palette,
                            bit_lengths, bit_codes);

  // No transforms.
  WriteBits(1, 0, bw);

  // Huffman image + meta huffman.
  const bool write_histogram_image = (histogram_image.size() > 1);
  WriteBits(1, write_histogram_image, bw);
  if (write_histogram_image) {
    std::vector<uint32_t> histogram_argb(histogram_symbols.begin(),
                                       histogram_symbols.end());
    ShiftHistogramImage(&histogram_argb);
    WriteBits(4, histogram_bits, bw);
    EncodeImageInternal(MetaSize(xsize, histogram_bits),
                        MetaSize(ysize, histogram_bits),
                        histogram_argb,
                        quality,
                        0,
                        0,      // no histogram bits
                        true,   // use_2d_locality
                        bw);
    const int image_size_bits = BitsLog2Ceiling(histogram_image.size() - 1);
    WriteBits(4, image_size_bits, bw);
    WriteBits(image_size_bits, histogram_image.size() - 2, bw);
    const int num_histograms = 5 * histogram_image.size();
    int nbits = BitsLog2Ceiling(num_histograms);
    WriteBits(4, nbits, bw);
    for (int i = 0; i < num_histograms; ++i) {
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
  for (int i = 0; i < histogram_image.size(); ++i) {
    std::vector<uint8_t> literal_lengths;
    PackLiteralBitLengths(bit_lengths[5 * i], palette_bits, use_palette,
                          &literal_lengths);
    StoreHuffmanCode(literal_lengths, bw);
    for (int k = 1; k < 5; ++k) {
      StoreHuffmanCode(bit_lengths[5 * i + k], bw);
    }
  }
  DeleteHistograms(histogram_image);  // free combined histograms

  // Emit no bits if there is only one symbol in the histogram.
  // This gives better compression for some images.
  for (int i = 0; i < bit_lengths.size(); ++i) {
    ClearHuffmanTreeIfOnlyOneSymbol(bit_lengths[i].size(), &(bit_lengths[i])[0],
                                    &bit_codes[i]);
  }

  // Store actual literals
  StoreImageToBitMask(xsize, ysize, histogram_bits, palette_bits,
                      &backward_refs[0], backward_refs.size(),
                      histogram_symbols, bit_lengths, bit_codes, bw);
}

inline void WriteImageSize(uint32_t xsize, uint32_t ysize, BitWriter *bw) {
  --xsize;
  --ysize;
  VERIFY(xsize < 0x4000 && ysize < 0x4000);
  WriteBits(14, xsize, bw);
  WriteBits(14, ysize, bw);
}

int EncodeWebpLLImage(int xsize, int ysize, const uint32_t *argb_orig,
                      EncodingStrategy *strategy,
                      size_t *num_bytes, uint8_t **bytes) {
  const int quality = strategy->quality;
  //  const int use_lz77 = strategy->use_lz77;
  int palette_bits = strategy->palette_bits;
  const int use_small_palette = strategy->use_small_palette;
  const int predict = strategy->predict;
  const int predict_bits = strategy->predict_bits;
  const int histogram_bits = strategy->histogram_bits;
  const int cross_color_transform = strategy->cross_color_transform;
  const int cross_color_transform_bits = strategy->cross_color_transform_bits;
  bool use_emerging_palette = true;
  std::vector<uint32_t> argb(argb_orig, argb_orig + xsize * ysize);
  const int kEstmatedEncodeSize = 0.5 * xsize * ysize;  // 4 bpp.

  BitWriter bw;
  if (!BitWriterInit(&bw, kEstmatedEncodeSize)) return false;

  // Write image size.
  WriteImageSize(xsize, ysize, &bw);

  if (!use_small_palette) {
    Histogram *before = new Histogram(1);
    for (int i = 0; i < xsize * ysize; ++i) {
      before->AddSingleLiteralOrCopy(LiteralOrCopy::CreateLiteral(argb[i]));
    }
    SubtractGreenFromBlueAndRed(xsize * ysize, &argb[0]);
    Histogram *after = new Histogram(1);
    for (int i = 0; i < xsize * ysize; ++i) {
      after->AddSingleLiteralOrCopy(LiteralOrCopy::CreateLiteral(argb[i]));
    }
    if (after->EstimateBits() < before->EstimateBits()) {
      WriteBits(1, 1, &bw);
      WriteBits(3, 2, &bw);
    } else {
      // Undo subtract green from blue and red -- rewrite with original data.
      argb = std::vector<uint32_t>(argb_orig, argb_orig + xsize * ysize);
    }
    delete before;
    delete after;
  }

  std::vector<uint32_t> palette;
  if (use_small_palette &&
      CreatePalette256(xsize * ysize, &argb[0], &palette)) {
    for (int i = 0; i < xsize * ysize; ++i) {
      for (int k = 0; k < palette.size(); ++k) {
        if (argb[i] == palette[k]) {
          argb[i] = 0xff000000 | (k << 8);
          break;
        }
      }
    }
    WriteBits(1, 1, &bw);
    WriteBits(3, 3, &bw);
    WriteBits(8, palette.size() - 1, &bw);
    for (int i = int(palette.size()) - 1; i >= 1; --i) {
      palette[i] = Subtract(palette[i], palette[i - 1]);
    }
    EncodeImageInternal(palette.size(), 1, palette, quality, 0, 0, true, &bw);
    use_emerging_palette = false;
    if (palette.size() <= 16) {
      int xbits = 1;
      if (palette.size() <= 2) {
        xbits = 3;
      } else if (palette.size() <= 4) {
        xbits = 2;
      }
      std::vector<uint32_t> from_argb(argb.begin(), argb.end());
      xsize = BundlePixels(xsize, ysize, xbits, from_argb, &argb);
      WriteBits(1, 1, &bw);
      WriteBits(3, 4, &bw);
      WriteBits(2, xbits, &bw);
    }
  }

  if (predict) {
    const int predictor_image_size =
        MetaSize(xsize, predict_bits) *
        MetaSize(ysize, predict_bits);
    std::vector<uint32_t> predictor_image(predictor_image_size);
    std::vector<uint32_t> from_argb(argb.begin(), argb.end());
    PredictorImage(xsize, ysize, predict_bits,
                   &from_argb[0],
                   &argb[0],
                   &predictor_image[0]);
    WriteBits(1, 1, &bw);
    WriteBits(3, 0, &bw);
    WriteBits(4, predict_bits, &bw);
    EncodeImageInternal(MetaSize(xsize, predict_bits),
                        MetaSize(ysize, predict_bits),
                        predictor_image,
                        quality,
                        0,
                        0,      // no histogram bits
                        true,   // use_2d_locality
                        &bw);
  }

  if (cross_color_transform) {
    const int color_space_image_size =
        MetaSize(xsize, cross_color_transform_bits) *
        MetaSize(ysize, cross_color_transform_bits);
    std::vector<uint32_t> color_space_image(color_space_image_size);
    ColorSpaceTransform(xsize, ysize, cross_color_transform_bits,
                        &argb[0],
                        quality,
                        &argb[0],
                        &color_space_image[0]);
    WriteBits(1, 1, &bw);
    WriteBits(3, 1, &bw);
    WriteBits(4, cross_color_transform_bits, &bw);
    EncodeImageInternal(MetaSize(xsize, cross_color_transform_bits),
                        MetaSize(ysize, cross_color_transform_bits),
                        color_space_image,
                        quality,
                        0,
                        0,      // no histogram bits
                        true,   // use_2d_locality
                        &bw);
  }

  palette_bits = 0;
  if (use_emerging_palette) {
    palette_bits = (quality == 0) ?
        7 :
        CalculateEstimateForPaletteSize(&argb[0], xsize, ysize);
  }

  EncodeImageInternal(xsize,
                      ysize,
                      argb,
                      quality,
                      palette_bits,
                      histogram_bits,
                      true,   // use_2d_locality
                      &bw);

  const size_t webpll_size = BitWriterNumBytes(&bw);
  uint8_t *webpll_data = BitWriterFinish(&bw);

  *num_bytes = HEADER_SIZE + SIGNATURE_SIZE + webpll_size;
  // TODO(vikasa): This memory allocation can be avoided, if RIFF header is
  // writen to BitWriter in the begining and BitWriter's buffer can be owned
  // and passed back instead.
  *bytes = (uint8_t *)malloc(*num_bytes);
  if (*bytes == NULL) return false;
  PutRiffHeader(*bytes, webpll_size);
  memcpy(*bytes + HEADER_SIZE + SIGNATURE_SIZE, webpll_data, webpll_size);
  BitWriterDestroy(&bw);

  return true;
}
