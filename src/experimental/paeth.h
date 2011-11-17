// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Spatial prediction using Paeth filter
//
// Author: Urvang (urvang@google.com)

#ifndef WEBP_UTILS_PAETH_H_
#define WEBP_UTILS_PAETH_H_

#include "../webp/types.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

// Filter the given data using Paeth predictor.
// 'data' corresponds to a 2-dimensional pixel array of size (stride * height)
// in raster order.
// 'bpp' is number of bytes per pixel, and
// 'stride' is number of bytes per scan line (with possible padding).
// 'filtered_data' should be pre-allocated.
void PaethFilter(const uint8_t* data, int width, int height, int bpp,
                 int stride, uint8_t* filtered_data);

// Reconstruct the original data from the given filtered data.
// 'recon_data' should be pre-allocated.
void PaethReconstruct(const uint8_t* data, int width, int height, int bpp,
                      int stride, uint8_t* recon_data);

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif  /* WEBP_UTILS_PAETH_H_ */
