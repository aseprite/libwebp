// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Alpha plane encoding and decoding library.
//
// Author: vikasa@google.com (Vikas Arora)

#ifndef WEBP_EXPERIMENTAL_ALPHA_H_
#define WEBP_EXPERIMENTAL_ALPHA_H_

#include <stdlib.h>

#include "../webp/types.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

// Experimental method for encoding/decoding Alpha plane.
// The pre-processing (Quantization) is performed if 'quality' is less than 100.
// For such cases, the encoding is lossy.
// 'method' can be:
//   0: - No compression;
//   1: - zlib;
//   2: - Use WebPLL

int EncodeAlphaExperimental(const uint8_t* data, int width, int height,
                            int stride, int quality, int method,
                            uint8_t** output, size_t* output_size);
// Decode the alpha data with the ad-hoc method, and fills the 'output' plane,
// which must be allocated and of size 'height x stride'.
int DecodeAlphaExperimental(const uint8_t* data, size_t data_size,
                            int width, int height, int stride, uint8_t* output);

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif  /* WEBP_EXPERIMENTAL_ALPHA_H_ */
