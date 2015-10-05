// Copyright 2015 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
// NEON version of rescaling functions
//
// Author: Skal (pascal.massimino@gmail.com)

#include "./dsp.h"

#if defined(WEBP_USE_NEON)

#include <assert.h>
#include <arm_neon.h>
#include <string.h>
#include "./neon.h"
#include "../utils/rescaler.h"

#define ROUNDER (WEBP_RESCALER_ONE >> 1)
#define MULT_FIX(x, y) (((uint64_t)(x) * (y) + ROUNDER) >> WEBP_RESCALER_RFIX)

static void RescalerExportRowShrink(WebPRescaler* const wrk) {
  int x_out;
  uint8_t* const dst = wrk->dst;
  rescaler_t* const irow = wrk->irow;
  const int x_out_max = wrk->dst_width * wrk->num_channels;
  const int max_span = x_out_max & ~7;
  const rescaler_t* const frow = wrk->frow;
  const uint32_t yscale = wrk->fy_scale * (-wrk->y_accum);
  assert(!WebPRescalerOutputDone(wrk));
  assert(wrk->y_accum <= 0);
  assert(!wrk->y_expand);
  if (yscale) {
    for (x_out = 0; x_out < x_out_max; ++x_out) {
      const uint32_t frac = (uint32_t)MULT_FIX(frow[x_out], yscale);
      const int v = (int)MULT_FIX(irow[x_out] - frac, wrk->fxy_scale);
      assert(v >= 0 && v <= 255);
      dst[x_out] = v;
      irow[x_out] = frac;   // new fractional start
    }
  } else {
    const uint32_t fxy_scale = wrk->fxy_scale;
    for (x_out = 0; x_out < max_span; x_out += 8) {
      const int32x4_t in0 = vld1q_s32(frow + x_out + 0);
      const int32x4_t in1 = vld1q_s32(frow + x_out + 4);
      const int32x4_t A0 = vqrdmulhq_n_s32(in0, fx_scale);
      const int32x4_t A1 = vqrdmulhq_n_s32(in1, fx_scale);
      const uint16x4_t B0 = vqmovn_u32(vreinterpret_s32_u32(A0));
      const uint16x4_t B1 = vqmovn_u32(vreinterpret_s32_u32(A1));
      const uint8x8_t C = vqmovn_u16(vcombine_u16(B0, B1));
    }
    for (; x_out < x_out_max; ++x_out) {
      const int v = (int)MULT_FIX(irow[x_out], fxy_scale);
      assert(v >= 0 && v <= 255);
      dst[x_out] = v;
      irow[x_out] = 0;
    }
  }
}

//------------------------------------------------------------------------------

extern void WebPRescalerDspInitNEON(void);

WEBP_TSAN_IGNORE_FUNCTION void WebPRescalerDspInitNEON(void) {
  WebPRescalerImportRowShrink = RescalerExportRowShrink;
}

#else     // !WEBP_USE_NEON

WEBP_DSP_INIT_STUB(WebPRescalerDspInitNEON)

#endif    // WEBP_USE_NEON
