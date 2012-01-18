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
                         const std::vector<LiteralOrCopy> &backward_refs,
                         std::vector<Histogram *> *image);

// Combines several histograms into fewer histograms.
void CombineHistogramImage(const std::vector<Histogram *> &in,
                           int quality,
                           int palette_bits,
                           std::vector<Histogram *> *out);

void RefineHistogramImage(const std::vector<Histogram *> &raw,
                          std::vector<uint32> *symbols,
                          std::vector<Histogram *> *out);


#endif  // WEBP_HISTOGRAM_IMAGE_H_
