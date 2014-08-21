// Copyright 2014 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------

#ifndef WEBP_DSP_IO_H_
#define WEBP_DSP_IO_H_

#include "../dec/decode_vp8.h"
#include "../dec/webpi.h"

static int GetAlphaSourceRow(const VP8Io* const io,
                             const uint8_t** alpha, int* const num_rows) {
  int start_y = io->mb_y;
  *num_rows = io->mb_h;

  // Compensate for the 1-line delay of the fancy upscaler.
  // This is similar to EmitFancyRGB().
  if (io->fancy_upsampling) {
    if (start_y == 0) {
      // We don't process the last row yet. It'll be done during the next call.
      --*num_rows;
    } else {
      --start_y;
      // Fortunately, *alpha data is persistent, so we can go back
      // one row and finish alpha blending, now that the fancy upscaler
      // completed the YUV->RGB interpolation.
      *alpha -= io->width;
    }
    if (io->crop_top + io->mb_y + io->mb_h == io->crop_bottom) {
      // If it's the very last call, we process all the remaining rows!
      *num_rows = io->crop_bottom - io->crop_top - start_y;
    }
  }
  return start_y;
}

typedef int (*VP8EmitAlphaFunc)(const VP8Io* const io, WebPDecParams* const p);

extern VP8EmitAlphaFunc VP8EmitAlphaRGB;

void VP8IOInit(void);

#endif  // WEBP_DSP_IO_H_
