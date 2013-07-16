// Copyright 2013 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
// Utilities for premultiplying or remultiplying ARGB values
//
// Author: Pascal Massimino (skal@google.com)

#ifndef WEBP_UTILS_ALPHA_MULTIPLY_H_
#define WEBP_UTILS_ALPHA_MULTIPLY_H_

#include "../webp/types.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

// Pre-Multiply or Remultiply argb values in a row.
void WebPMultARGBRow(uint32_t* const ptr, int width, int do_remultiply);

// Same a WebPMultARGBRow(), but for several rows.
void WebPMultARGBRows(uint8_t* ptr, int stride, int width, int height,
                      int do_remultiply);

// same for a row of single values, with side alpha values.
void WebPMultRow(uint8_t* const ptr, const uint8_t* const alpha,
                 int width, int do_remultiply);

// Same a WebPMultRow(), but for several 'num_rows' rows.
void WebPMultRows(uint8_t* ptr, int stride,
                  const uint8_t* alpha, int alpha_stride,
                  int width, int num_rows, int do_remultiply);

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif  // WEBP_UTILS_ALPHA_MULTIPLY_H_
