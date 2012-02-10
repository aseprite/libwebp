// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)
//
// Models a 2d image of histograms.

#ifndef WEBP_HISTOGRAM_IMAGE_H_
#define WEBP_HISTOGRAM_IMAGE_H_

#include "backward_references.h"
#include "histogram.h"
#include "../common/integral_types.h"

void BuildHistogramImage(int xsize, int ysize,
                         int histobits,
                         int palette_bits,
                         const PixOrCopy *backward_refs,
                         int backward_refs_size,
                         Histogram ***image,
                         int *histogram_size);

// Combines several histograms into fewer histograms.
void CombineHistogramImage(Histogram **in,
                           int in_size,
                           int quality,
                           int palette_bits,
                           Histogram ***out,
                           int *out_size);

void RefineHistogramImage(Histogram **raw,
                          int raw_size,
                          uint32_t *symbols,
                          int out_size,
                          Histogram **out);


#endif  // WEBP_HISTOGRAM_IMAGE_H_
