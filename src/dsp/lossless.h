// Copyright 2012 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Image transforms and color space conversion methods for lossless decoder.
//
// Author: Vikas Arora (vikaas.arora@gmail.com)
//         jyrki@google.com (Jyrki Alakuijala)

#ifndef WEBP_DSP_LOSSLESS_H_
#define WEBP_DSP_LOSSLESS_H_

#include "../webp/types.h"
#include "../webp/decode.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

//------------------------------------------------------------------------------
// Inverse image transforms.

struct VP8LTransform;  // Defined in dec/vp8li.h.

// Performs inverse transform of data given transform information.
// Returns appropriate status code.
VP8StatusCode VP8LInverseTransform(const struct VP8LTransform* const transform,
                                   argb_t** const data);

//------------------------------------------------------------------------------
// Color space conversion.

// Converts from BGRA to other color spaces.
// Returns true on success.
int VP8LConvertColorSpaceFromBGRA(const uint8_t* const in_data,
                                  size_t num_pixels,
                                  WEBP_CSP_MODE out_colorspace,
                                  uint8_t** const output_data);

//------------------------------------------------------------------------------
// Misc methods.

// Computes the sum of each component with mod 256.
static WEBP_INLINE uint32_t VP8LAddPixels(uint32_t a, uint32_t b) {
  const uint32_t alpha_and_green = (a & 0xff00ff00) + (b & 0xff00ff00);
  const uint32_t red_and_blue = (a & 0x00ff00ff) + (b & 0x00ff00ff);
  return (alpha_and_green & 0xff00ff00) | (red_and_blue & 0x00ff00ff);
}

// Computes sampled size of 'size' when sampling using 'sampling bits'.
static WEBP_INLINE uint32_t VP8LSubSampleSize(uint32_t size,
                                              uint32_t sampling_bits) {
  return (size + (1 << sampling_bits) - 1) >> sampling_bits;
}

//------------------------------------------------------------------------------

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif  // WEBP_DSP_LOSSLESS_H_
