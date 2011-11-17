// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)
//
// Simple predictor for webp lossless.

#ifndef WEBP_PREDICTOR_H_
#define WEBP_PREDICTOR_H_

#include "../common/integral_types.h"

void PredictorImage(int xsize, int ysize, int bits,
                    const uint32 *original_argb,
                    uint32 *to_argb,
                    uint32 *image);

// Get in the image, produces a subresolution transform image.
// The transform minimize locally the entropy of red and blue by
// finding cross-component correlation to green.
void ColorSpaceTransform(int xsize, int ysize, int bits,
                         const uint32 *original_argb,
                         int quality,
                         uint32 *to_argb,
                         uint32 *image);

void SubtractGreenFromBlueAndRed(int n, uint32 *argb_array);

#endif  // WEBP_PREDICTOR_H_
