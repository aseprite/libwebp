// Copyright 2011 Google Inc. All Rights Reserved.
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

#include <stdint.h>

void PredictorImage(int xsize, int ysize, int bits,
                    const uint32_t *original_argb,
                    uint32_t *to_argb,
                    uint32_t *image);

// Get in the image, produces a subresolution transform image.
// The transform minimize locally the entropy of red and blue by
// finding cross-component correlation to green.
void ColorSpaceTransform(int xsize, int ysize, int bits,
                         const uint32_t *original_argb,
                         int quality,
                         uint32_t *to_argb,
                         uint32_t *image);

void SubtractGreenFromBlueAndRed(int n, uint32_t *argb_array);

#endif  // WEBP_PREDICTOR_H_
