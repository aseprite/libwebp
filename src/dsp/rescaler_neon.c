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
  const int32_t yscale = wrk->fy_scale * (-wrk->y_accum);
  const int32_t fxy_scale = wrk->fxy_scale;
  const int32x4_t zero = vdupq_n_s32(0);
  assert(!WebPRescalerOutputDone(wrk));
  assert(wrk->y_accum <= 0);
  assert(!wrk->y_expand);
  if (yscale) {
    for (x_out = 0; x_out < max_span; x_out += 8) {
      const int32x4_t in0 = vld1q_s32(frow + x_out + 0);
      const int32x4_t in1 = vld1q_s32(frow + x_out + 4);
      const int32x4_t A0 = vqrdmulhq_n_s32(in0, yscale);
      const int32x4_t A1 = vqrdmulhq_n_s32(in1, yscale);
      const int32x4_t in2 = vld1q_s32(irow + x_out + 0);
      const int32x4_t in3 = vld1q_s32(irow + x_out + 4);
      const int32x4_t B0 = vqsubq_s32(in2, A0);
      const int32x4_t B1 = vqsubq_s32(in3, A1);
      const int32x4_t C0 = vqrdmulhq_n_s32(B0, fxy_scale);
      const int32x4_t C1 = vqrdmulhq_n_s32(B1, fxy_scale);
      const int16x4_t D0 = vqmovn_s32(C0);
      const int16x4_t D1 = vqmovn_s32(C1);
      const uint8x8_t E = vqmovun_s16(vcombine_s16(D0, D1));
      vst1_u8(dst + x_out, E);
      vst1q_s32(irow + x_out + 0, A0);
      vst1q_s32(irow + x_out + 4, A1);
    }
    for (; x_out < x_out_max; ++x_out) {
      const uint32_t frac = (uint32_t)MULT_FIX(frow[x_out], yscale);
      const int v = (int)MULT_FIX(irow[x_out] - frac, wrk->fxy_scale);
      assert(v >= 0 && v <= 255);
      dst[x_out] = v;
      irow[x_out] = frac;   // new fractional start
    }
  } else {
    for (x_out = 0; x_out < max_span; x_out += 8) {
      const int32x4_t in0 = vld1q_s32(irow + x_out + 0);
      const int32x4_t in1 = vld1q_s32(irow + x_out + 4);
      const int32x4_t A0 = vqrdmulhq_n_s32(in0, fxy_scale);
      const int32x4_t A1 = vqrdmulhq_n_s32(in1, fxy_scale);
      const int16x4_t B0 = vqmovn_s32(A0);
      const int16x4_t B1 = vqmovn_s32(A1);
      const uint8x8_t C = vqmovun_s16(vcombine_s16(B0, B1));
      vst1_u8(dst + x_out, C);
      vst1q_s32(irow + x_out + 0, zero);
      vst1q_s32(irow + x_out + 4, zero);
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
  WebPRescalerExportRowShrink = RescalerExportRowShrink;
  (void)RescalerExportRowShrink;
}

#else     // !WEBP_USE_NEON

WEBP_DSP_INIT_STUB(WebPRescalerDspInitNEON)

#endif    // WEBP_USE_NEON
