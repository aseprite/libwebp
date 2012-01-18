// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)
//
// Models the histograms of literal and distance codes.

#ifndef WEBP_HISTOGRAM_H_
#define WEBP_HISTOGRAM_H_

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../common/integral_types.h"
#include "backward_references.h"

class LiteralOrCopy;

// A simple container for histograms of data.
// This data structure models five histograms:
// 1) green literal, palette-code and copy-length-prefix histogram
// 2) red histogram
// 3) blue histogram
// 4) alpha histogram
// 5) distance-prefix histogram.
struct Histogram {
  explicit Histogram(int palette_code_bits)
      : palette_code_bits_(palette_code_bits) {
    Clear();
  }

  // Create the histogram.
  //
  // The input data is the LiteralOrCopy data, which models the
  // literals, stop codes and backward references (both distances and lengths)
  void Build(const LiteralOrCopy *literal_and_length,
             int n_literal_and_length);
  void Clear() {
    memset(&literal_[0], 0, sizeof(literal_));
    memset(&red_[0], 0, sizeof(red_));
    memset(&blue_[0], 0, sizeof(blue_));
    memset(&alpha_[0], 0, sizeof(alpha_));
    memset(&distance_[0], 0, sizeof(distance_));
  }
  void AddSingleLiteralOrCopy(const LiteralOrCopy &v);

  // Estimate how many bits the combined entropy of literals and distance
  // approximately maps to.
  double EstimateBits() const;

  // This function estimates the Huffman dictionary + other block overhead
  // size for creating a new deflate block.
  double EstimateBitsHeader() const;

  // This function estimates the cost in bits excluding the bits needed to
  // represent the entropy code itself.
  double EstimateBitsBulk() const;

  void Add(const Histogram &a) {
    for (int i = 0; i < kLiteralOrCopyCodesMax; ++i) {
      literal_[i] += a.literal_[i];
    }
    for (int i = 0; i < kDistanceCodes; ++i) {
      distance_[i] += a.distance_[i];
    }
    for (int i = 0; i < 256; ++i) {
      red_[i] += a.red_[i];
      blue_[i] += a.blue_[i];
      alpha_[i] += a.alpha_[i];
    }
  }

  void Remove(const Histogram &a) {
    for (int i = 0; i < kLiteralOrCopyCodesMax; ++i) {
      literal_[i] -= a.literal_[i];
      assert(literal_[i] >= 0);
    }
    for (int i = 0; i < kDistanceCodes; ++i) {
      distance_[i] -= a.distance_[i];
      assert(distance_[i] >= 0);
    }
    for (int i = 0; i < 256; ++i) {
      red_[i] -= a.red_[i];
      blue_[i] -= a.blue_[i];
      alpha_[i] -= a.alpha_[i];
      assert(red_[i] >= 0);
      assert(blue_[i] >= 0);
      assert(alpha_[i] >= 0);
    }
  }

  int NumLiteralOrCopyCodes() const {
    return 256 + kLengthCodes + (1 << palette_code_bits_);
  }

  // green or length-prefix code histogram.
  int literal_[kLiteralOrCopyCodesMax];
  int red_[256];
  int blue_[256];
  int alpha_[256];
  // Backward reference prefix-code histogram.
  int distance_[kDistanceCodes];
  int palette_code_bits_;
 private:
  Histogram();
};

void ConvertPopulationCountTableToBitEstimates(
    int n, const int *population_counts,
    double *output);

double BitsEntropy(const int *array, int n);

#endif  // WEBP_HISTOGRAM_H_
