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

#include <algorithm>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "../common/integral_types.h"
#include "../common/predictor.h"
#include "backward_distance.h"
#include "backward_references.h"
#include "bit_stream.h"
#include "encode.h"
#include "entropy_encode.h"
#include "histogram_image.h"
#include "predictor.h"

static const uint8 kMagicByte = 0xa3;

static void PutLE32(uint8* const data, uint32 val) {
  data[0] = (val >>  0) & 0xff;
  data[1] = (val >>  8) & 0xff;
  data[2] = (val >> 16) & 0xff;
  data[3] = (val >> 24) & 0xff;
}

static void PutRiffHeader(uint8 *bytes, size_t webpll_size) {
  uint8 riff[HEADER_SIZE + SIGNATURE_SIZE] = {
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
        int count = std::max(1, (sum + stride / 2) / stride);
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
                                     uint8* lengths,
                                     std::vector<uint16>* symbols) {
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

void PackGreenBitLengths(const std::vector<uint8> bit_lengths,
                         const int green_bits,
                         const int palette_bits,
                         const bool use_palette,
                         std::vector<uint8>* new_lengths) {
  for (int i = 0; i < (1 << green_bits); ++i) {
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
    const uint8 *code_length_bitdepth,
    int *storage_ix,
    uint8 *storage) {
  // RFC 1951 will calm you down if you are worried about this funny sequence.
  // This sequence is tuned from that, but more weighted for lower symbol count,
  // and more spiking histograms.
  static const uint8 kStorageOrder[kCodeLengthCodes] = {
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
  WriteBits(4, codes_to_store - 4, storage_ix, storage);
  for (uint32 i = 0; i < codes_to_store; ++i) {
    WriteBits(3, code_length_bitdepth[kStorageOrder[i]], storage_ix, storage);
  }
}

void StoreHuffmanTreeToBitMask(
    const std::vector<uint8> &huffman_tree,
    const std::vector<uint8> &huffman_tree_extra_bits,
    const int num_symbols,
    const uint8 *code_length_bitdepth,
    const std::vector<uint16> &code_length_bitdepth_symbols,
    int * __restrict storage_ix,
    uint8 * __restrict storage) {
  for (uint32 i = 0; i < num_symbols; ++i) {
    int ix = huffman_tree[i];
    WriteBits(code_length_bitdepth[ix],
              code_length_bitdepth_symbols[ix],
              storage_ix,
              storage);
    switch (ix) {
      case 16:
        WriteBits(2, huffman_tree_extra_bits[i], storage_ix, storage);
        break;
      case 17:
        WriteBits(3, huffman_tree_extra_bits[i], storage_ix, storage);
        break;
      case 18:
        WriteBits(7, huffman_tree_extra_bits[i], storage_ix, storage);
        break;
    }
  }
}

void StoreHuffmanCode(const std::vector<uint8> &bit_lengths,
                      int* storage_ix,
                      uint8* storage) {
  int count = 0;
  int symbols[2] = { 0 };
  for (int i = 0; i < bit_lengths.size(); ++i) {
    if (bit_lengths[i] != 0) {
      if (count < 2) symbols[count] = i;
      ++count;
    }
  }
  if (count <= 2 && (count == 0 || symbols[count - 1] < (1 << 12))) {
    WriteBits(1, 1, storage_ix, storage);
    if (count == 0) {
      WriteBits(4, 0, storage_ix, storage);
      return;
    }
    WriteBits(1, count - 1, storage_ix, storage);
    int num_bits = 4;
    while (symbols[count - 1] >= (1 << num_bits)) num_bits += 2;
    WriteBits(3, (num_bits - 4) / 2 + 1, storage_ix, storage);
    for (int i = 0; i < count; ++i) {
      WriteBits(num_bits, symbols[i], storage_ix, storage);
    }
    return;
  }
  WriteBits(1, 0, storage_ix, storage);
  std::vector<uint8> huffman_tree(bit_lengths.size());
  std::vector<uint8> huffman_tree_extra_bits(bit_lengths.size());
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
  uint8 code_length_bitdepth[kCodeLengthCodes] = { 0 };
  std::vector<uint16> code_length_bitdepth_symbols(kCodeLengthCodes);
  CreateHuffmanTree(&huffman_tree_histogram[0], kCodeLengthCodes,
                    7, &code_length_bitdepth[0]);
  ConvertBitDepthsToSymbols(&code_length_bitdepth[0], kCodeLengthCodes,
                            &code_length_bitdepth_symbols[0]);
  StoreHuffmanTreeOfHuffmanTreeToBitMask(code_length_bitdepth,
                                         storage_ix, storage);
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
  WriteBits(1, write_length, storage_ix, storage);
  if (write_length) {
    const int nbits = BitsLog2Ceiling(trimmed_length - 1);
    const int nbitpairs = nbits == 0 ? 1 : (nbits + 1) / 2;
    WriteBits(3, nbitpairs - 1, storage_ix, storage);
    WriteBits(nbitpairs * 2, trimmed_length - 2, storage_ix, storage);
  }
  StoreHuffmanTreeToBitMask(huffman_tree,
                            huffman_tree_extra_bits,
                            length,
                            &code_length_bitdepth[0],
                            code_length_bitdepth_symbols,
                            storage_ix, storage);
}

void StoreImageToBitMask(
    const int xsize,
    const int ysize,
    const int histobits,
    const int palette_bits,
    const LiteralOrCopy *literal,
    const int n_literal_and_length,
    const std::vector<uint32> &histogram_symbol,
    const std::vector< std::vector<uint8> > &bitdepth,
    const std::vector< std::vector<uint16> > &bit_symbols,
    int *storage_ix,
    uint8 *storage) {
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
                bit_symbols[5 * histogram_ix][literal_ix],
                storage_ix, storage);
    } else if (v.IsLiteral()) {
      int order[] = {1, 2, 0, 3};
      for (int i = 0; i < 4; ++i) {
        const int code = v.Literal(order[i]);
        WriteBits(bitdepth[5 * histogram_ix + i][code],
                  bit_symbols[5 * histogram_ix + i][code],
                  storage_ix, storage);
      }
    } else {
      int code;
      int n_bits;
      int bits;
      v.LengthCodeAndBits(&code, &n_bits, &bits);
      int len_ix = 256 + (1 << palette_bits) + code;
      WriteBits(bitdepth[5 * histogram_ix][len_ix],
                bit_symbols[5 * histogram_ix][len_ix],
                storage_ix, storage);
      WriteBits(n_bits, bits, storage_ix, storage);

      const int distance = v.Distance();
      BackwardDistance::Encode(distance, &code, &n_bits, &bits);
      WriteBits(bitdepth[5 * histogram_ix + 4][code],
                bit_symbols[5 * histogram_ix + 4][code],
                storage_ix, storage);
      WriteBits(n_bits, bits, storage_ix, storage);
    }
    x += v.Length();
    while (x >= xsize) {
      x -= xsize;
      ++y;
    }
  }
}

void ShiftHistogramImage(std::vector<uint32>* image) {
  for (int i = 0; i < image->size(); ++i) {
    (*image)[i] <<= 8;
    (*image)[i] |= 0xff000000;
  }
}

int NumberOfUsedRBAComponents(const std::vector<uint32>& argb) {
  int num = 0;
  for (int i = 0; i < argb.size() && num < 3; ++i) {
    if (((argb[i] >> 24) & 0xff) != 0xff) return 3;
    if ((argb[i] & 0xff) != 0 && num < 2) num = 2;
    if (((argb[i] >> 16) & 0xff) != 0 && num < 1) num = 1;
  }
  return num;
}

int GreenBitDepth(const std::vector<uint32>& argb) {
  int bits = 1;
  for (int i = 0; i < argb.size(); ++i) {
    while (((argb[i] >> 8) & 0xff) >= (1 << bits)) ++bits;
  }
  return bits;
}

int MetaSize(int size, int bits) {
  return (size + (1 << bits) - 1) >> bits;
}

bool CreatePalette(int n, const uint32 *argb,
                   int max_palette_size,
                   std::vector<uint32>* palette) {
  if (max_palette_size == 0 || n == 0) return false;
  std::set<uint32> k;
  k.insert(argb[0]);
  for (int i = 1; i < n; ++i) {
    if (argb[i - 1] == argb[i]) {
      continue;
    }
    if (k.find(argb[i]) == k.end()) {
      k.insert(argb[i]);
      if (k.size() > max_palette_size) {
        return false;
      }
    }
  }
  for (std::set<uint32>::const_iterator it = k.begin(); it != k.end(); ++it) {
    palette->push_back(*it);
  }
  return true;
}

void BundlePixels(int* xsize, int* ysize, int xbits, int ybits, int nbits,
                  const std::vector<uint32>& from_argb,
                  std::vector<uint32>* to_argb) {
  VERIFY((*xsize) * (*ysize) == from_argb.size());
  int xs = MetaSize(*xsize, xbits);
  int ys = MetaSize(*ysize, ybits);
  to_argb->resize(xs * ys);
  for (int tile_y = 0; tile_y < ys; ++tile_y) {
    for (int tile_x = 0; tile_x < xs; ++tile_x) {
      int tile_ix = tile_y * xs + tile_x;
      uint32 tile_code = 0;
      for (int y = 0; y < (1 << ybits); ++y) {
        int all_y = tile_y * (1 << ybits) + y;
        if (all_y >= (*ysize)) continue;
        for (int x = 0; x < (1 << xbits); ++x) {
          int all_x = tile_x * (1 << xbits) + x;
          if (all_x >= (*xsize)) continue;
          int all_ix = all_y * (*xsize) + all_x;
          int bit_position = (y * (1 << xbits) + x) * nbits;
          uint32 mask = (1 << nbits) - 1;
          if (bit_position < 16) {
            bit_position += 8;
          } else if (bit_position < 24) {
            bit_position -= 16;
          }
          tile_code |= ((from_argb[all_ix] >> 8) & mask) << bit_position;
        }
      }
      (*to_argb)[tile_ix] = tile_code;
      if ((1 << xbits) * (1 << ybits) * nbits <= 24) {
        (*to_argb)[tile_ix] = (*to_argb)[tile_ix] | 0xff000000;
      }
    }
  }
  *xsize = xs;
  *ysize = ys;
}

static void DeleteHistograms(std::vector<Histogram *> &histograms) {
  int i;
  for (i = 0; i < histograms.size(); ++i) {
    delete histograms[i];
  }
}

void EncodeImageInternal(const int xsize,
                         const int ysize,
                         const std::vector<uint32>& argb,
                         const int quality,
                         const int palette_bits,
                         const int histogram_bits,
                         const bool use_2d_locality,
                         const bool write_error_detection_bits,
                         int *storage_ix,
                         uint8 *storage) {
  const int use_palette = palette_bits ? 1 : 0;
  // First, check if it is at all a good idea to use LZ77
  bool lz77_is_useful = false;
  std::vector<LiteralOrCopy> backward_refs_lz77;
  std::vector<LiteralOrCopy> backward_refs_rle_only;
  {
    BackwardReferencesHashChain(
        xsize,
        ysize,
        use_palette,
        &argb[0],
        palette_bits,
        &backward_refs_lz77);
    Histogram *histo_lz77 = new Histogram(palette_bits);
    histo_lz77->Build(&backward_refs_lz77[0], backward_refs_lz77.size());

    BackwardReferencesRle(
        xsize,
        ysize,
        &argb[0],
        &backward_refs_rle_only);
    Histogram *histo_rle = new Histogram(palette_bits);
    histo_rle->Build(&backward_refs_rle_only[0], backward_refs_rle_only.size());

    lz77_is_useful = histo_rle->EstimateBits() > histo_lz77->EstimateBits();
    if (lz77_is_useful) {
      backward_refs_rle_only.clear();
    } else {
      backward_refs_lz77.clear();
    }
    delete histo_rle;
    delete histo_lz77;
  }

  std::vector<LiteralOrCopy> backward_refs;
  if (quality >= 50 && lz77_is_useful) {
    int recursion_level = 0;
    if (xsize * ysize < 320 * 200) {
      recursion_level = 1;
    }
    BackwardReferencesTraceBackwards(
        xsize,
        ysize,
        recursion_level,
        use_palette,
        &argb[0],
        palette_bits,
        &backward_refs);
  } else {
    if (lz77_is_useful) {
      backward_refs.swap(backward_refs_lz77);
    } else {
      backward_refs.swap(backward_refs_rle_only);
    }
  }
  VERIFY(VerifyBackwardReferences(&argb[0], xsize, ysize,
                                 palette_bits,
                                 backward_refs));

  if (use_2d_locality) {
    BackwardReferences2DLocality(xsize,
                                 ysize,
                                 backward_refs.size(),
                                 &backward_refs[0]);
  }

  std::vector<Histogram *> histogram_image;
  std::vector<Histogram *> histogram_image_raw;
  BuildHistogramImage(xsize,
                      ysize,
                      histogram_bits,
                      palette_bits,
                      backward_refs,
                      &histogram_image_raw);
  std::vector<uint32> histogram_symbols(histogram_image_raw.size(), -1);
  CombineHistogramImage(histogram_image_raw,
                        quality,
                        palette_bits,
                        &histogram_image);
  const int max_refinement_iters = 1;
  for (int iter = 0; iter < max_refinement_iters; ++iter) {
    RefineHistogramImage(histogram_image_raw, &histogram_symbols,
                         &histogram_image);
    if (quality < 30) {
      break;
    }
  }

  DeleteHistograms(histogram_image_raw);  // free raw histograms

  // Create bit lengths for each histogram code.
  std::vector< std::vector<uint8> > bit_length(5 * histogram_image.size());
  for (int i = 0; i < histogram_image.size(); ++i) {
    bit_length[5 * i].resize(histogram_image[i]->NumLiteralOrCopyCodes());
    OptimizeHuffmanForRle(histogram_image[i]->NumLiteralOrCopyCodes(),
                          &histogram_image[i]->literal_[0]);
    if (!use_palette) {
      // Optimization might have smeared population counts to palette entries,
      // so zero them out afterwards.
      for (int k = 0; k < (1 << palette_bits); ++k) {
        histogram_image[i]->literal_[256 + k] = 0;
      }
    }
    OptimizeHuffmanForRle(256, &histogram_image[i]->red_[0]);
    OptimizeHuffmanForRle(256, &histogram_image[i]->blue_[0]);
    OptimizeHuffmanForRle(256, &histogram_image[i]->alpha_[0]);
    OptimizeHuffmanForRle(kDistanceCodes, &histogram_image[i]->distance_[0]);

    CreateHuffmanTree(histogram_image[i]->literal_,
                      histogram_image[i]->NumLiteralOrCopyCodes(), 15,
                      &bit_length[5 * i][0]);
    bit_length[5 * i + 1].resize(256);
    CreateHuffmanTree(histogram_image[i]->red_, 256, 15,
                      &bit_length[5 * i + 1][0]);

    bit_length[5 * i + 2].resize(256);
    CreateHuffmanTree(histogram_image[i]->blue_, 256, 15,
                      &bit_length[5 * i + 2][0]);
    bit_length[5 * i + 3].resize(256);
    CreateHuffmanTree(histogram_image[i]->alpha_, 256, 15,
                      &bit_length[5 * i + 3][0]);
    bit_length[5 * i + 4].resize(kDistanceCodes);
    CreateHuffmanTree(histogram_image[i]->distance_, kDistanceCodes, 15,
                      &bit_length[5 * i + 4][0]);
  }
  // We have all Huffman trees (modeled by their bit lengths).

  // Now create the actual bit codes for the bit lengths.
  std::vector< std::vector<uint16> > bit_codes(bit_length.size());
  for (int i = 0; i < bit_codes.size(); ++i) {
    // Create actual bit codes
    bit_codes[i].resize(bit_length[i].size());
    ConvertBitDepthsToSymbols(&bit_length[i][0], bit_length[i].size(),
                              &bit_codes[i][0]);
  }

  // No transforms.
  WriteBits(1, 0, storage_ix, storage);

  WriteBits(1, write_error_detection_bits, storage_ix, storage);
  if (write_error_detection_bits) {
    WriteBits(8, kMagicByte, storage_ix, storage);
  }

  // This many of the red, blue, alpha components we have in addition to green
  // (in this order).
  const int num_rba = NumberOfUsedRBAComponents(argb);
  WriteBits(2, num_rba, storage_ix, storage);

  //
  // Huffman image + meta huffman
  //
  const bool write_histogram_image = (histogram_image.size() > 1);
  WriteBits(1, write_histogram_image, storage_ix, storage);
  if (write_histogram_image) {
    std::vector<uint32> histogram_argb(histogram_symbols.begin(),
                                  histogram_symbols.end());
    ShiftHistogramImage(&histogram_argb);
    WriteBits(4, histogram_bits, storage_ix, storage);
    EncodeImageInternal(MetaSize(xsize, histogram_bits),
                        MetaSize(ysize, histogram_bits),
                        histogram_argb,
                        quality,
                        0,
                        0,      // no histogram bits
                        true,   // use_2d_locality
                        write_error_detection_bits,
                        storage_ix, storage);
    const int image_size_bits = BitsLog2Ceiling(histogram_image.size() - 1);
    WriteBits(4, image_size_bits, storage_ix, storage);
    WriteBits(image_size_bits, histogram_image.size() - 2, storage_ix, storage);
    const int num_histograms = (num_rba + 2) * histogram_image.size();
    int nbits = BitsLog2Ceiling(num_histograms);
    WriteBits(4, nbits, storage_ix, storage);
    for (int i = 0; i < num_histograms; ++i) {
      WriteBits(nbits, i, storage_ix, storage);
    }
    if (write_error_detection_bits) {
      WriteBits(8, kMagicByte, storage_ix, storage);
    }
  }

  // Palette parameters
  WriteBits(1, use_palette, storage_ix, storage);
  if (use_palette) {
    WriteBits(4, kRowHasherXSubsampling, storage_ix, storage);
    WriteBits(4, palette_bits, storage_ix, storage);
  }

  // Green component bit depth, with the property that the green component
  // of all pixels is less than 1 << green_bit_depth.
  const int green_bit_depth = GreenBitDepth(argb);
  VERIFY(green_bit_depth >= 1);
  VERIFY(green_bit_depth <= 8);
  WriteBits(3, green_bit_depth - 1, storage_ix, storage);

  //
  // Huffman codes
  //
  for (int i = 0; i < histogram_image.size(); ++i) {
    std::vector<uint8> green_lengths;
    PackGreenBitLengths(bit_length[5 * i], green_bit_depth,
                        palette_bits, use_palette, &green_lengths);
    StoreHuffmanCode(green_lengths, storage_ix, storage);
    for (int k = 1; k <= num_rba; ++k) {
      StoreHuffmanCode(bit_length[5 * i + k], storage_ix, storage);
    }
    StoreHuffmanCode(bit_length[5 * i + 4], storage_ix, storage);
  }
  DeleteHistograms(histogram_image);  // free combined histograms

  if (write_error_detection_bits) {
    WriteBits(8, kMagicByte, storage_ix, storage);
  }

  // Emit no bits if there is only one symbol in the histogram.
  // This gives ~5 % for lena.png.
  for (int i = 0; i < bit_length.size(); ++i) {
    ClearHuffmanTreeIfOnlyOneSymbol(bit_length[i].size(),
                                    &(bit_length[i])[0],
                                    &bit_codes[i]);
  }

  // Store actual literals
  StoreImageToBitMask(
      xsize, ysize,
      histogram_bits,
      palette_bits,
      &backward_refs[0],
      backward_refs.size(),
      histogram_symbols,
      bit_length,
      bit_codes,
      storage_ix,
      storage);

  if (write_error_detection_bits) {
    WriteBits(8, kMagicByte, storage_ix, storage);
  }
}

int EncodeWebpLLImage(int xsize, int ysize, const uint32 *argb_orig,
                      EncodingStrategy *strategy,
                      size_t *num_bytes, uint8 **bytes) {
  const int quality = strategy->quality;
  //  const int use_lz77 = strategy->use_lz77;
  int palette_bits = strategy->palette_bits;
  const int use_small_palette = strategy->use_small_palette;
  const int predict = strategy->predict;
  const int predict_bits = strategy->predict_bits;
  const int histogram_bits = strategy->histogram_bits;
  const int cross_color_transform = strategy->cross_color_transform;
  const int cross_color_transform_bits = strategy->cross_color_transform_bits;
  const int write_error_detection_bits = false;
  bool use_emerging_palette = true;
  const int kHeaderSize = 2048;
  std::vector<uint32> argb(argb_orig, argb_orig + xsize * ysize);
  // TODO: Come up with a good estimate w.r.t 'storage' size.
  std::vector<uint8> storage(xsize * ysize * 8 + kHeaderSize);
  int storage_ix = 0;

  WriteBitsPrepareStorage(storage_ix, &storage[0]);

  // Write image size.
  WriteImageSize(xsize, ysize, &storage_ix, &storage[0]);

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
      WriteBits(1, 1, &storage_ix, &storage[0]);
      WriteBits(3, 2, &storage_ix, &storage[0]);
    } else {
      // Undo subtract green from blue and red -- rewrite with original data.
      argb = std::vector<uint32>(argb_orig, argb_orig + xsize * ysize);
    }
    delete before;
    delete after;
  }

  std::vector<uint32> palette;
  const int max_palette_size = 256;
  if (use_small_palette &&
      CreatePalette(xsize * ysize, &argb[0], max_palette_size, &palette)) {
    for (int i = 0; i < xsize * ysize; ++i) {
      for (int k = 0; k < palette.size(); ++k) {
        if (argb[i] == palette[k]) {
          argb[i] = 0xff000000 | (k << 8);
          break;
        }
      }
    }
    WriteBits(1, 1, &storage_ix, &storage[0]);
    WriteBits(3, 3, &storage_ix, &storage[0]);
    WriteBits(8, palette.size() - 1, &storage_ix, &storage[0]);
    for (int i = int(palette.size()) - 1; i >= 1; --i) {
      palette[i] = Subtract(palette[i], palette[i - 1]);
    }
    EncodeImageInternal(palette.size(), 1, palette, quality, 0, 0, true,
                        write_error_detection_bits, &storage_ix, &storage[0]);
    use_emerging_palette = false;
    int ybits = 0;
    int bit_depth = 8;
    int xbits = 0;
    if (palette.size() <= 2) {
      bit_depth = 1;
      xbits = 3;
    } else if (palette.size() <= 4) {
      bit_depth = 2;
      xbits = 2;
    } else if (palette.size() <= 16) {
      bit_depth = 4;
      xbits = 1;
    }
    std::vector<uint32> from_argb(argb.begin(), argb.end());
    BundlePixels(&xsize, &ysize, xbits, ybits, bit_depth,
                 from_argb, &argb);
    WriteBits(1, 1, &storage_ix, &storage[0]);
    WriteBits(3, 4, &storage_ix, &storage[0]);
    WriteBits(3, xbits, &storage_ix, &storage[0]);
    WriteBits(3, ybits, &storage_ix, &storage[0]);
    WriteBits(3, bit_depth - 1, &storage_ix, &storage[0]);
  }

  if (predict) {
    const int predictor_image_size =
        MetaSize(xsize, predict_bits) *
        MetaSize(ysize, predict_bits);
    std::vector<uint32> predictor_image(predictor_image_size);
    std::vector<uint32> from_argb(argb.begin(), argb.end());
    PredictorImage(xsize, ysize, predict_bits,
                   &from_argb[0],
                   &argb[0],
                   &predictor_image[0]);
    WriteBits(1, 1, &storage_ix, &storage[0]);
    WriteBits(3, 0, &storage_ix, &storage[0]);
    WriteBits(4, predict_bits, &storage_ix, &storage[0]);
    EncodeImageInternal(MetaSize(xsize, predict_bits),
                        MetaSize(ysize, predict_bits),
                        predictor_image,
                        quality,
                        0,
                        0,      // no histogram bits
                        true,   // use_2d_locality
                        write_error_detection_bits,
                        &storage_ix, &storage[0]);
  }

  if (cross_color_transform) {
    const int color_space_image_size =
        MetaSize(xsize, cross_color_transform_bits) *
        MetaSize(ysize, cross_color_transform_bits);
    std::vector<uint32> color_space_image(color_space_image_size);
    ColorSpaceTransform(xsize, ysize, cross_color_transform_bits,
                        &argb[0],
                        quality,
                        &argb[0],
                        &color_space_image[0]);
    WriteBits(1, 1, &storage_ix, &storage[0]);
    WriteBits(3, 1, &storage_ix, &storage[0]);
    WriteBits(4, cross_color_transform_bits, &storage_ix, &storage[0]);
    EncodeImageInternal(MetaSize(xsize, cross_color_transform_bits),
                        MetaSize(ysize, cross_color_transform_bits),
                        color_space_image,
                        quality,
                        0,
                        0,      // no histogram bits
                        true,   // use_2d_locality
                        write_error_detection_bits,
                        &storage_ix, &storage[0]);
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
                      write_error_detection_bits,
                      &storage_ix, &storage[0]);

  const size_t webpll_size = (storage_ix + 7) >> 3;
  uint8 *webpll_data = &storage[0];

  *num_bytes = HEADER_SIZE + SIGNATURE_SIZE + webpll_size;
  *bytes = (uint8 *)malloc(*num_bytes);
  if (*bytes == NULL) return false;
  PutRiffHeader(*bytes, webpll_size);
  memcpy(*bytes + HEADER_SIZE + SIGNATURE_SIZE, webpll_data, webpll_size);

  return true;
}
