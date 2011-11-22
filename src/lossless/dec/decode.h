// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: szabadka@google.com (Zoltan Szabadka)

#ifndef WEBP_DEC_DECODE_H_
#define WEBP_DEC_DECODE_H_

#include <stdlib.h>
#include "../common/integral_types.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

// Main decoding function for webp lossless image.
// Input is the webp lossless image and output is the width ("xsize"),
// height ("ysize") of the image and the array of argb pixel values.
//
// Return false in case error, otherwise true.
int DecodeWebpLLImage(size_t encoded_image_size,
                      const uint8* const encoded_image,
                      int* xsize,
                      int* ysize,
                      uint32** argb_image);

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif  // WEBP_DEC_DECODE_H_
