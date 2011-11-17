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

#include "../common/integral_types.h"

bool DecodeWebpLLImage(int encoded_image_size,
                       const uint8* const encoded_image,
                       int* xsize,
                       int* ysize,
                       uint32** argb_image);

#endif  // WEBP_DEC_DECODE_H_
