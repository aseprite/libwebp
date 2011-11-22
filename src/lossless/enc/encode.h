// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: somnath@google.com (Somnath Banerjee)

#ifndef WEBP_ENC_ENCODE_H_
#define WEBP_ENC_ENCODE_H_

#include <stdlib.h>
#include "../common/integral_types.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

// Main encoding call for webp lossless compression.
// Input is an array of argb pixel values. Output "stream" is the compressed
// webp lossless image.
// Returns false in case of error, otherwise true.
int EncodeWebpLLImage(const int xsize,
                      const int ysize,
                      const uint32 *argb,
                      const int quality,
                      const int use_small_palette,
                      const int predict,
                      const int predict_bits,
                      const int histogram_bits,
                      const int cross_color_transform,
                      const int color_transform_bits,
                      const int write_error_detection_bits,
                      size_t *num_bytes,
                      uint8 **bytes);

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif  // WEBP_ENC_ENCODE_H_
