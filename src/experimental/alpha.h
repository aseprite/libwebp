// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Alpha plane encoding and decoding library.
//
// Author: vikasa@google.com (Vikas Arora)

#ifndef WEBP_UTILS_ALPHA_H_
#define WEBP_UTILS_ALPHA_H_

#include <stdlib.h>

#include "../webp/types.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

// Encodes the given Alpha data 'data' of size 'stride'x'height' via specified
// compression method 'method'. The pre-processing (Quantization) is
// performed if 'quality' is less than 100. For such cases, the encoding is
// lossy. Valid ranges for 'quality' is [0, 100] and 'method' is [0, 2]:
//   'method = 0' - No compression;
//   'method = 1' - zlib;
//   'method = 2' - Run Length Encoding of Symbols and Arithmetic Encoding of
//                  <symbol, lenght> pair.
// 'output' corresponds to the buffer containing compressed Alpha data.
//          This buffer is allocated by this method and caller should call
//          free(*output) when done.
// 'output_size' corresponds to size of this compressed Alpha buffer.
//
// Returns 1 on successfully encoding the Alpha and
//         0 if either:
//           data, output or output_size is NULL, or
//           inappropriate width, height or stride, or
//           invalid quality or method, or
//           Memory allocation for the compressed data fails.

int EncodeAlpha(const uint8_t* data, int width, int height, int stride,
                int quality, int method,
                uint8_t** output, size_t* output_size);

// Decodes the compressed data 'data' of size 'data_size' into the 'output'.
// The 'output' buffer should be pre-alloacated and must be of the same
// dimension 'height'x'stride', as that of the image.
//
// Returns 1 on successfully decoding the compressed Alpha and
//         0 if either:
//           data or output is NULL, or
//           Error in bit-stream header (invalid compression mode or qbits), or
//           Error returned by approppriate compression method.
int DecodeAlpha(const uint8_t* data, size_t data_size,
                int width, int height, int stride, uint8_t* output);

// Replace the input 'data' of size 'width'x'height' with 'num-levels'
// quantized values. If not NULL, 'mse' will contain the mean-squared error.
// Valid range for 'num_levels' is [2, 256].
// Returns false in case of error (data is NULL, or parameters are invalid).
int QuantizeLevels(uint8_t* data, int width, int height, int num_levels,
                   float* mse);

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif  /* WEBP_UTILS_ALPHA_H_ */
