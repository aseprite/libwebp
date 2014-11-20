// Copyright 2014 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
//   ARGB making functions.
//
// Author: Djordje Pesut (djordje.pesut@imgtec.com)

#ifndef WEBP_DSP_ARGB_H_
#define WEBP_DSP_ARGB_H_

#include "../webp/types.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void (*VP8PackARGB)(const uint8_t* a, const uint8_t* r,
                           const uint8_t* g, const uint8_t* b, int len,
                           int step, int offset, uint32_t* out);
extern void (*VP8PackRGB)(const uint8_t* r, const uint8_t* g, const uint8_t* b,
                          int len, int step, int offset, uint32_t* out);

void VP8EncDspARGBInit(void);

#ifdef __cplusplus
}    // extern "C"
#endif

#endif  /* WEBP_DSP_ARGB_H_ */
