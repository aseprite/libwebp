// Copyright 2012 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// main entry for the decoder
//
// Author: Vikas Arora (vikaas.arora@gmail.com)

#include <stdlib.h>

#include "./vp8li.h"
#include "../utils/bit_reader.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#define LOSSLESS_MAGIC_BYTE 0x64

static const size_t kHeaderBytes = 5;
static const uint32_t kImageSizeBits = 14;

//------------------------------------------------------------------------------

int VP8LGetInfo(const uint8_t* data, uint32_t data_size,
                int* width, int* height) {
  if (data_size >= kHeaderBytes) {
    BitReader br;
    int signature;

    InitSimpleBitReader(&br, data, kHeaderBytes);
    signature = ReadBits(&br, 8);
    if (signature != LOSSLESS_MAGIC_BYTE) return 0;
    *width = ReadBits(&br, kImageSizeBits) + 1;
    *height = ReadBits(&br, kImageSizeBits) + 1;
    return 1;
  } else {
    return 0;         // not enough data
  }
}

//------------------------------------------------------------------------------

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif
