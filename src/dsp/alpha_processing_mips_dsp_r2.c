// Copyright 2014 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
// Utilities for processing transparent channel.
//
// Author(s): Branimir Vasic (branimir.vasic@imgtec.com)
//            Djordje Pesut  (djordje.pesut@imgtec.com)

#include "./dsp.h"

#if defined(WEBP_USE_MIPS_DSP_R2)

static uint32_t FindAlphaMask(uint8_t* const base_rgba, const int alpha_first,
                              const int num_rows, const int mb_w,
                              const uint8_t* alpha, const int io_width,
                              const int buf_stride) {
  uint8_t* dst = base_rgba + (alpha_first ? 0 : 3);
  uint32_t alpha_mask = 0xffffffff;
  int i, j, temp0;

  for (j = 0; j < num_rows; ++j) {
    uint8_t* pdst = dst;
    uint8_t* palpha = (uint8_t*)alpha;
    for (i = 0; i < (mb_w >> 2); ++i) {
      int temp1, temp2, temp3;

      __asm__ volatile (
        "ulw    %[temp0],      0(%[palpha])                \n\t"
        "addiu  %[palpha],     %[palpha],     4            \n\t"
        "addiu  %[pdst],       %[pdst],       16           \n\t"
        "srl    %[temp1],      %[temp0],      8            \n\t"
        "srl    %[temp2],      %[temp0],      16           \n\t"
        "srl    %[temp3],      %[temp0],      24           \n\t"
        "and    %[alpha_mask], %[alpha_mask], %[temp0]     \n\t"
        "sb     %[temp0],      -16(%[pdst])                \n\t"
        "sb     %[temp1],      -12(%[pdst])                \n\t"
        "sb     %[temp2],      -8(%[pdst])                 \n\t"
        "sb     %[temp3],      -4(%[pdst])                 \n\t"
        : [temp0]"=&r"(temp0), [temp1]"=&r"(temp1), [temp2]"=&r"(temp2),
          [temp3]"=&r"(temp3), [palpha]"+r"(palpha), [pdst]"+r"(pdst),
          [alpha_mask]"+r"(alpha_mask)
        :
        : "memory"
      );
    }

    for (i = 0; i < (mb_w & 3); ++i) {
      __asm__ volatile (
        "lbu    %[temp0],      0(%[palpha])                \n\t"
        "addiu  %[palpha],     %[palpha],     1            \n\t"
        "sb     %[temp0],      0(%[pdst])                  \n\t"
        "and    %[alpha_mask], %[alpha_mask], %[temp0]     \n\t"
        "addiu  %[pdst],       %[pdst],       4            \n\t"
        : [temp0]"=&r"(temp0), [palpha]"+r"(palpha), [pdst]"+r"(pdst),
          [alpha_mask]"+r"(alpha_mask)
        :
        : "memory"
      );
    }
    alpha += io_width;
    dst += buf_stride;
  }

  __asm__ volatile (
    "ext    %[temp0],      %[alpha_mask], 0, 16            \n\t"
    "srl    %[alpha_mask], %[alpha_mask], 16               \n\t"
    "and    %[alpha_mask], %[alpha_mask], %[temp0]         \n\t"
    "ext    %[temp0],      %[alpha_mask], 0, 8             \n\t"
    "srl    %[alpha_mask], %[alpha_mask], 8                \n\t"
    "and    %[alpha_mask], %[alpha_mask], %[temp0]         \n\t"
    : [temp0]"=&r"(temp0), [alpha_mask]"+r"(alpha_mask)
    :
  );

  return alpha_mask;
}

#endif  // WEBP_USE_MIPS_DSP_R2

//------------------------------------------------------------------------------
// Init function

extern void WebPInitAlphaProcessingMIPSdspR2(void);

void WebPInitAlphaProcessingMIPSdspR2(void) {
#if defined(WEBP_USE_MIPS_DSP_R2)
  WebPFindAlphaMask = FindAlphaMask;
#endif
}
