// Copyright 2011 Google Inc. All Rights Reserved.
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

// A struct to control the encoding strategy.
typedef struct {
  int quality;
  int use_lz77;
  int use_small_palette;
  int palette_bits;
  int predict;
  int predict_bits;
  int histogram_bits;
  int cross_color_transform;
  int cross_color_transform_bits;
} EncodingStrategy;

// Main encoding call for webp lossless compression.
// Input is an array of argb pixel values. Output "stream" is the compressed
// webp lossless image.
// Returns false in case of error, otherwise true.
int EncodeWebpLLImage(const int xsize,
                      const int ysize,
                      const uint32 *argb,
                      EncodingStrategy *encoding_strategy,
                      size_t *num_bytes,
                      uint8 **bytes);

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif  // WEBP_ENC_ENCODE_H_
