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
typedef struct {
  // literal_ contains green literal, palette-code and
  // copy-length-prefix histogram
  int literal_[kLiteralOrCopyCodesMax];
  int red_[256];
  int blue_[256];
  int alpha_[256];
  // Backward reference prefix-code histogram.
  int distance_[kDistanceCodes];
  int palette_code_bits_;
} Histogram;

static inline void Histogram_Clear(Histogram *p) {
  memset(&p->literal_[0], 0, sizeof(p->literal_));
  memset(&p->red_[0], 0, sizeof(p->red_));
  memset(&p->blue_[0], 0, sizeof(p->blue_));
  memset(&p->alpha_[0], 0, sizeof(p->alpha_));
  memset(&p->distance_[0], 0, sizeof(p->distance_));
}

static inline void Histogram_Init(Histogram *p, int palette_code_bits) {
  p->palette_code_bits_ = palette_code_bits;
  Histogram_Clear(p);
}

// Create the histogram.
//
// The input data is the LiteralOrCopy data, which models the
// literals, stop codes and backward references (both distances and lengths)
void Histogram_Build(Histogram *p,
                     const LiteralOrCopy *literal_and_length,
                     int n_literal_and_length);
void Histogram_AddSingleLiteralOrCopy(Histogram *p,
                                      const LiteralOrCopy &v);

// Estimate how many bits the combined entropy of literals and distance
// approximately maps to.
double Histogram_EstimateBits(const Histogram * const p);

// This function estimates the Huffman dictionary + other block overhead
// size for creating a new deflate block.
double Histogram_EstimateBitsHeader(const Histogram * const p);

// This function estimates the cost in bits excluding the bits needed to
// represent the entropy code itself.
double Histogram_EstimateBitsBulk(const Histogram * const p);

static inline void Histogram_Add(Histogram *p, const Histogram *a) {
  int i;
  for (i = 0; i < kLiteralOrCopyCodesMax; ++i) {
    p->literal_[i] += a->literal_[i];
  }
  for (i = 0; i < kDistanceCodes; ++i) {
    p->distance_[i] += a->distance_[i];
  }
  for (i = 0; i < 256; ++i) {
    p->red_[i] += a->red_[i];
    p->blue_[i] += a->blue_[i];
    p->alpha_[i] += a->alpha_[i];
  }
}

static inline void Histogram_Remove(Histogram *p, const Histogram *a) {
  int i;
  for (i = 0; i < kLiteralOrCopyCodesMax; ++i) {
    p->literal_[i] -= a->literal_[i];
    assert(p->literal_[i] >= 0);
  }
  for (i = 0; i < kDistanceCodes; ++i) {
    p->distance_[i] -= a->distance_[i];
    assert(p->distance_[i] >= 0);
  }
  for (i = 0; i < 256; ++i) {
    p->red_[i] -= a->red_[i];
    p->blue_[i] -= a->blue_[i];
    p->alpha_[i] -= a->alpha_[i];
    assert(p->red_[i] >= 0);
    assert(p->blue_[i] >= 0);
    assert(p->alpha_[i] >= 0);
  }
}

static inline int Histogram_NumLiteralOrCopyCodes(const Histogram *p) {
  return 256 + kLengthCodes + (1 << p->palette_code_bits_);
}

void ConvertPopulationCountTableToBitEstimates(
    int n, const int *population_counts,
    double *output);

double BitsEntropy(const int *array, int n);

#endif  // WEBP_HISTOGRAM_H_
