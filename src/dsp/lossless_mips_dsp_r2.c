// Copyright 2014 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
// Image transforms and color space conversion methods for lossless decoder.
//
// Author(s):  Djordje Pesut    (djordje.pesut@imgtec.com)
//             Jovan Zelincevic (jovan.zelincevic@imgtec.com)

#include "./dsp.h"

#if defined(WEBP_USE_MIPS_DSP_R2)

#include "./lossless.h"

#define MAP_COLOR_FUNCS(FUNC_NAME, TYPE, GET_INDEX, GET_VALUE)                 \
static void FUNC_NAME(const TYPE* src,                                         \
                      const uint32_t* const color_map,                         \
                      TYPE* dst, int y_start, int y_end,                       \
                      int width) {                                             \
  int y;                                                                       \
  for (y = y_start; y < y_end; ++y) {                                          \
    int x;                                                                     \
    for (x = 0; x < (width >> 2); ++x) {                                       \
      int tmp1, tmp2, tmp3, tmp4;                                              \
      __asm__ volatile (                                                       \
      ".ifc        "#TYPE",  uint8_t                    \n\t"                  \
        "lbu       %[tmp1],  0(%[src])                  \n\t"                  \
        "lbu       %[tmp2],  1(%[src])                  \n\t"                  \
        "lbu       %[tmp3],  2(%[src])                  \n\t"                  \
        "lbu       %[tmp4],  3(%[src])                  \n\t"                  \
        "addiu     %[src],   %[src],      4             \n\t"                  \
      ".endif                                           \n\t"                  \
      ".ifc        "#TYPE",  uint32_t                   \n\t"                  \
        "lw        %[tmp1],  0(%[src])                  \n\t"                  \
        "lw        %[tmp2],  4(%[src])                  \n\t"                  \
        "lw        %[tmp3],  8(%[src])                  \n\t"                  \
        "lw        %[tmp4],  12(%[src])                 \n\t"                  \
        "ext       %[tmp1],  %[tmp1],     8,        8   \n\t"                  \
        "ext       %[tmp2],  %[tmp2],     8,        8   \n\t"                  \
        "ext       %[tmp3],  %[tmp3],     8,        8   \n\t"                  \
        "ext       %[tmp4],  %[tmp4],     8,        8   \n\t"                  \
        "addiu     %[src],   %[src],      16            \n\t"                  \
      ".endif                                           \n\t"                  \
        "sll       %[tmp1],  %[tmp1],     2             \n\t"                  \
        "sll       %[tmp2],  %[tmp2],     2             \n\t"                  \
        "sll       %[tmp3],  %[tmp3],     2             \n\t"                  \
        "sll       %[tmp4],  %[tmp4],     2             \n\t"                  \
        "lwx       %[tmp1],  %[tmp1](%[color_map])      \n\t"                  \
        "lwx       %[tmp2],  %[tmp2](%[color_map])      \n\t"                  \
        "lwx       %[tmp3],  %[tmp3](%[color_map])      \n\t"                  \
        "lwx       %[tmp4],  %[tmp4](%[color_map])      \n\t"                  \
      ".ifc        "#TYPE",  uint8_t                    \n\t"                  \
        "ext       %[tmp1],  %[tmp1],     8,        8   \n\t"                  \
        "ext       %[tmp2],  %[tmp2],     8,        8   \n\t"                  \
        "ext       %[tmp3],  %[tmp3],     8,        8   \n\t"                  \
        "ext       %[tmp4],  %[tmp4],     8,        8   \n\t"                  \
        "sb        %[tmp1],  0(%[dst])                  \n\t"                  \
        "sb        %[tmp2],  1(%[dst])                  \n\t"                  \
        "sb        %[tmp3],  2(%[dst])                  \n\t"                  \
        "sb        %[tmp4],  3(%[dst])                  \n\t"                  \
        "addiu     %[dst],   %[dst],      4             \n\t"                  \
      ".endif                                           \n\t"                  \
      ".ifc        "#TYPE",  uint32_t                   \n\t"                  \
        "sw        %[tmp1],  0(%[dst])                  \n\t"                  \
        "sw        %[tmp2],  4(%[dst])                  \n\t"                  \
        "sw        %[tmp3],  8(%[dst])                  \n\t"                  \
        "sw        %[tmp4],  12(%[dst])                 \n\t"                  \
        "addiu     %[dst],   %[dst],      16            \n\t"                  \
      ".endif                                           \n\t"                  \
        : [tmp1]"=&r"(tmp1), [tmp2]"=&r"(tmp2), [tmp3]"=&r"(tmp3),             \
          [tmp4]"=&r"(tmp4), [src]"+&r"(src), [dst]"+r"(dst)                   \
        : [color_map]"r"(color_map)                                            \
        : "memory"                                                             \
      );                                                                       \
    }                                                                          \
    for (x = 0; x < (width & 3); ++x) {                                        \
      *dst++ = GET_VALUE(color_map[GET_INDEX(*src++)]);                        \
    }                                                                          \
  }                                                                            \
}

MAP_COLOR_FUNCS(MapARGB, uint32_t, VP8GetARGBIndex, VP8GetARGBValue)
MAP_COLOR_FUNCS(MapAlpha, uint8_t, VP8GetAlphaIndex, VP8GetAlphaValue)

#undef MAP_COLOR_FUNCS

static WEBP_INLINE uint32_t ClampedAddSubtractFull(uint32_t c0, uint32_t c1,
                                                   uint32_t c2) {
  int temp0, temp1, temp2, temp3, temp4, temp5;
  __asm__ volatile (
    "preceu.ph.qbr   %[temp1],   %[c0]                 \n\t"
    "preceu.ph.qbl   %[temp2],   %[c0]                 \n\t"
    "preceu.ph.qbr   %[temp3],   %[c1]                 \n\t"
    "preceu.ph.qbl   %[temp4],   %[c1]                 \n\t"
    "preceu.ph.qbr   %[temp5],   %[c2]                 \n\t"
    "preceu.ph.qbl   %[temp0],   %[c2]                 \n\t"
    "subq.ph         %[temp3],   %[temp3],   %[temp5]  \n\t"
    "subq.ph         %[temp4],   %[temp4],   %[temp0]  \n\t"
    "addq.ph         %[temp1],   %[temp1],   %[temp3]  \n\t"
    "addq.ph         %[temp2],   %[temp2],   %[temp4]  \n\t"
    "shll_s.ph       %[temp1],   %[temp1],   7         \n\t"
    "shll_s.ph       %[temp2],   %[temp2],   7         \n\t"
    "precrqu_s.qb.ph %[temp2],   %[temp2],   %[temp1]  \n\t"
    : [temp0]"=r"(temp0), [temp1]"=&r"(temp1), [temp2]"=&r"(temp2),
      [temp3]"=&r"(temp3), [temp4]"=&r"(temp4), [temp5]"=&r"(temp5)
    : [c0]"r"(c0), [c1]"r"(c1), [c2]"r"(c2)
    : "memory"
  );
  return temp2;
}

static WEBP_INLINE uint32_t ClampedAddSubtractHalf(uint32_t c0, uint32_t c1,
                                                   uint32_t c2) {
  int temp0, temp1, temp2, temp3, temp4, temp5;
  __asm__ volatile (
    "adduh.qb         %[temp5],   %[c0],      %[c1]       \n\t"
    "preceu.ph.qbr    %[temp3],   %[c2]                   \n\t"
    "preceu.ph.qbr    %[temp1],   %[temp5]                \n\t"
    "preceu.ph.qbl    %[temp2],   %[temp5]                \n\t"
    "preceu.ph.qbl    %[temp4],   %[c2]                   \n\t"
    "subq.ph          %[temp3],   %[temp1],   %[temp3]    \n\t"
    "subq.ph          %[temp4],   %[temp2],   %[temp4]    \n\t"
    "shrl.ph          %[temp5],   %[temp3],   15          \n\t"
    "shrl.ph          %[temp0],   %[temp4],   15          \n\t"
    "addq.ph          %[temp3],   %[temp3],   %[temp5]    \n\t"
    "addq.ph          %[temp4],   %[temp0],   %[temp4]    \n\t"
    "shra.ph          %[temp3],   %[temp3],   1           \n\t"
    "shra.ph          %[temp4],   %[temp4],   1           \n\t"
    "addq.ph          %[temp1],   %[temp1],   %[temp3]    \n\t"
    "addq.ph          %[temp2],   %[temp2],   %[temp4]    \n\t"
    "shll_s.ph        %[temp1],   %[temp1],   7           \n\t"
    "shll_s.ph        %[temp2],   %[temp2],   7           \n\t"
    "precrqu_s.qb.ph  %[temp1],   %[temp2],   %[temp1]    \n\t"
    : [temp0]"=r"(temp0), [temp1]"=&r"(temp1), [temp2]"=&r"(temp2),
      [temp3]"=&r"(temp3), [temp4]"=r"(temp4), [temp5]"=&r"(temp5)
    : [c0]"r"(c0), [c1]"r"(c1), [c2]"r"(c2)
    : "memory"
  );
  return temp1;
}

static void SubtractGreenFromBlueAndRed(uint32_t* argb_data,
                                        int num_pixels) {
  uint32_t temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7;
  uint32_t* const p_loop1_end = argb_data + (num_pixels & ~3);
  uint32_t* const p_loop2_end = p_loop1_end + (num_pixels & 3);
  __asm__ volatile (
    ".set       push                                          \n\t"
    ".set       noreorder                                     \n\t"
    "beq        %[argb_data],    %[p_loop1_end],     3f       \n\t"
    " nop                                                     \n\t"
  "0:                                                         \n\t"
    "lw         %[temp0],        0(%[argb_data])              \n\t"
    "lw         %[temp1],        4(%[argb_data])              \n\t"
    "lw         %[temp2],        8(%[argb_data])              \n\t"
    "lw         %[temp3],        12(%[argb_data])             \n\t"
    "ext        %[temp4],        %[temp0],           8,    8  \n\t"
    "ext        %[temp5],        %[temp1],           8,    8  \n\t"
    "ext        %[temp6],        %[temp2],           8,    8  \n\t"
    "ext        %[temp7],        %[temp3],           8,    8  \n\t"
    "addiu      %[argb_data],    %[argb_data],       16       \n\t"
    "replv.ph   %[temp4],        %[temp4]                     \n\t"
    "replv.ph   %[temp5],        %[temp5]                     \n\t"
    "replv.ph   %[temp6],        %[temp6]                     \n\t"
    "replv.ph   %[temp7],        %[temp7]                     \n\t"
    "subu.qb    %[temp0],        %[temp0],           %[temp4] \n\t"
    "subu.qb    %[temp1],        %[temp1],           %[temp5] \n\t"
    "subu.qb    %[temp2],        %[temp2],           %[temp6] \n\t"
    "subu.qb    %[temp3],        %[temp3],           %[temp7] \n\t"
    "sw         %[temp0],        -16(%[argb_data])            \n\t"
    "sw         %[temp1],        -12(%[argb_data])            \n\t"
    "sw         %[temp2],        -8(%[argb_data])             \n\t"
    "bne        %[argb_data],    %[p_loop1_end],     0b       \n\t"
    " sw        %[temp3],        -4(%[argb_data])             \n\t"
  "3:                                                         \n\t"
    "beq        %[argb_data],    %[p_loop2_end],     2f       \n\t"
    " nop                                                     \n\t"
  "1:                                                         \n\t"
    "lw         %[temp0],        0(%[argb_data])              \n\t"
    "addiu      %[argb_data],    %[argb_data],       4        \n\t"
    "ext        %[temp4],        %[temp0],           8,    8  \n\t"
    "replv.ph   %[temp4],        %[temp4]                     \n\t"
    "subu.qb    %[temp0],        %[temp0],           %[temp4] \n\t"
    "bne        %[argb_data],    %[p_loop2_end],     1b       \n\t"
    " sw        %[temp0],        -4(%[argb_data])             \n\t"
  "2:                                                         \n\t"
    ".set       pop                                           \n\t"
    : [argb_data]"+&r"(argb_data), [temp0]"=&r"(temp0),
      [temp1]"=&r"(temp1), [temp2]"=&r"(temp2), [temp3]"=&r"(temp3),
      [temp4]"=&r"(temp4), [temp5]"=&r"(temp5), [temp6]"=&r"(temp6),
      [temp7]"=&r"(temp7)
    : [p_loop1_end]"r"(p_loop1_end), [p_loop2_end]"r"(p_loop2_end)
    : "memory"
  );
}

static WEBP_INLINE uint32_t Select(uint32_t a, uint32_t b, uint32_t c) {
  int temp0, temp1, temp2, temp3, temp4, temp5;
  __asm__ volatile (
    "cmpgdu.lt.qb %[temp1], %[c],     %[b]             \n\t"
    "pick.qb      %[temp1], %[b],     %[c]             \n\t"
    "pick.qb      %[temp2], %[c],     %[b]             \n\t"
    "cmpgdu.lt.qb %[temp4], %[c],     %[a]             \n\t"
    "pick.qb      %[temp4], %[a],     %[c]             \n\t"
    "pick.qb      %[temp5], %[c],     %[a]             \n\t"
    "subu.qb      %[temp3], %[temp1], %[temp2]         \n\t"
    "subu.qb      %[temp0], %[temp4], %[temp5]         \n\t"
    "raddu.w.qb   %[temp3], %[temp3]                   \n\t"
    "raddu.w.qb   %[temp0], %[temp0]                   \n\t"
    "subu         %[temp3], %[temp3], %[temp0]         \n\t"
    "slti         %[temp0], %[temp3], 0x1              \n\t"
    "movz         %[a],     %[b],     %[temp0]         \n\t"
    : [temp1]"=&r"(temp1), [temp2]"=&r"(temp2), [temp3]"=&r"(temp3),
      [temp4]"=&r"(temp4), [temp5]"=&r"(temp5), [temp0]"=&r"(temp0),
      [a]"+&r"(a)
    : [b]"r"(b), [c]"r"(c)
  );
  return a;
}

static uint32_t Predictor11(uint32_t left, const uint32_t* const top) {
  return Select(top[0], left, top[-1]);
}

static uint32_t Predictor12(uint32_t left, const uint32_t* const top) {
  return ClampedAddSubtractFull(left, top[0], top[-1]);
}

static uint32_t Predictor13(uint32_t left, const uint32_t* const top) {
  return ClampedAddSubtractHalf(left, top[0], top[-1]);
}

static void TransformColor(const VP8LMultipliers* const m, uint32_t* data,
                           int num_pixels) {
  int temp0, temp1, temp2, temp3, temp4, temp5;
  uint32_t argb, argb1, new_red, new_red1;
  uint32_t G_to_R = m->green_to_red_;
  const uint32_t G_to_B = m->green_to_blue_;
  uint32_t R_to_B = m->red_to_blue_;
  uint32_t* const p_loop_end = data + (num_pixels & ~1);
  const uint32_t rest = num_pixels & 1;
  __asm__ volatile (
    ".set            push                                    \n\t"
    ".set            noreorder                               \n\t"
    "beq             %[data],      %[p_loop_end],  1f        \n\t"
    " nop                                                    \n\t"
    "replv.ph        %[temp0],     %[G_to_R]                 \n\t"
    "replv.ph        %[temp1],     %[G_to_B]                 \n\t"
    "replv.ph        %[temp2],     %[R_to_B]                 \n\t"
    "shll.ph         %[temp0],     %[temp0],       8         \n\t"
    "shll.ph         %[temp1],     %[temp1],       8         \n\t"
    "shll.ph         %[temp2],     %[temp2],       8         \n\t"
    "shra.ph         %[temp0],     %[temp0],       8         \n\t"
    "shra.ph         %[temp1],     %[temp1],       8         \n\t"
    "shra.ph         %[temp2],     %[temp2],       8         \n\t"
  "0:                                                        \n\t"
    "lw              %[argb],      0(%[data])                \n\t"
    "lw              %[argb1],     4(%[data])                \n\t"
    "lhu             %[new_red],   2(%[data])                \n\t"
    "lhu             %[new_red1],  6(%[data])                \n\t"
    "precrq.qb.ph    %[temp3],     %[argb],        %[argb1]  \n\t"
    "precr.qb.ph     %[temp4],     %[argb],        %[argb1]  \n\t"
    "preceu.ph.qbra  %[temp3],     %[temp3]                  \n\t"
    "preceu.ph.qbla  %[temp4],     %[temp4]                  \n\t"
    "shll.ph         %[temp3],     %[temp3],       8         \n\t"
    "shll.ph         %[temp4],     %[temp4],       8         \n\t"
    "shra.ph         %[temp3],     %[temp3],       8         \n\t"
    "shra.ph         %[temp4],     %[temp4],       8         \n\t"
    "mul.ph          %[temp5],     %[temp3],       %[temp0]  \n\t"
    "mul.ph          %[temp3],     %[temp3],       %[temp1]  \n\t"
    "mul.ph          %[temp4],     %[temp4],       %[temp2]  \n\t"
    "addiu           %[data],      %[data],        8         \n\t"
    "ins             %[new_red1],  %[new_red],     16,   16  \n\t"
    "ins             %[argb1],     %[argb],        16,   16  \n\t"
    "shra.ph         %[temp5],     %[temp5],       5         \n\t"
    "shra.ph         %[temp3],     %[temp3],       5         \n\t"
    "shra.ph         %[temp4],     %[temp4],       5         \n\t"
    "subu.ph         %[new_red1],  %[new_red1],    %[temp5]  \n\t"
    "subu.ph         %[argb1],     %[argb1],       %[temp3]  \n\t"
    "preceu.ph.qbra  %[temp5],     %[new_red1]               \n\t"
    "subu.ph         %[argb1],     %[argb1],       %[temp4]  \n\t"
    "preceu.ph.qbra  %[temp3],     %[argb1]                  \n\t"
    "sb              %[temp5],     -2(%[data])               \n\t"
    "sb              %[temp3],     -4(%[data])               \n\t"
    "sra             %[temp5],     %[temp5],       16        \n\t"
    "sra             %[temp3],     %[temp3],       16        \n\t"
    "sb              %[temp5],     -6(%[data])               \n\t"
    "bne             %[data],      %[p_loop_end],  0b        \n\t"
    " sb             %[temp3],     -8(%[data])               \n\t"
  "1:                                                        \n\t"
    "beqz            %[rest],      2f                        \n\t"
    " nop                                                    \n\t"
    "lw              %[argb],      0(%[data])                \n\t"
    "lhu             %[new_red],   2(%[data])                \n\t"
    "lb              %[temp3],     1(%[data])                \n\t"
    "ins             %[G_to_R],    %[G_to_B],      16,   16  \n\t"
    "replv.ph        %[temp3],     %[temp3]                  \n\t"
    "lb              %[temp0],     2(%[data])                \n\t"
    "shll.ph         %[G_to_R],    %[G_to_R],      8         \n\t"
    "shra.ph         %[G_to_R],    %[G_to_R],      8         \n\t"
    "shll.ph         %[R_to_B],    %[R_to_B],      8         \n\t"
    "shra.ph         %[R_to_B],    %[R_to_B],      8         \n\t"
    "mul.ph          %[temp5],     %[temp3],       %[G_to_R] \n\t"
    "mul.ph          %[temp0],     %[temp0],       %[R_to_B] \n\t"
    "ins             %[new_red],   %[argb],        16,   16  \n\t"
    "shra.ph         %[temp5],     %[temp5],       5         \n\t"
    "shra.ph         %[temp0],     %[temp0],       5         \n\t"
    "subu.ph         %[temp4],     %[new_red],     %[temp5]  \n\t"
    "sra             %[temp5],     %[temp4],       16        \n\t"
    "subu.ph         %[temp5],     %[temp5],       %[temp0]  \n\t"
    "sb              %[temp4],     2(%[data])                \n\t"
    "sb              %[temp5],     0(%[data])                \n\t"
  "2:                                                        \n\t"
    ".set            pop                                     \n\t"
    : [temp0]"=&r"(temp0), [temp1]"=&r"(temp1), [temp2]"=&r"(temp2),
      [temp3]"=&r"(temp3), [temp4]"=&r"(temp4), [temp5]"=&r"(temp5),
      [new_red1]"=&r"(new_red1), [new_red]"=&r"(new_red),
      [argb]"=&r"(argb), [argb1]"=&r"(argb1), [data]"+&r"(data),
      [G_to_R]"+r"(G_to_R), [R_to_B]"+r"(R_to_B)
    : [G_to_B]"r"(G_to_B), [p_loop_end]"r"(p_loop_end),
      [rest]"r"(rest)
    : "memory", "hi", "lo"
  );
}

#endif  // WEBP_USE_MIPS_DSP_R2

//------------------------------------------------------------------------------

extern void VP8LDspInitMIPSdspR2(void);

void VP8LDspInitMIPSdspR2(void) {
#if defined(WEBP_USE_MIPS_DSP_R2)
  VP8LMapColor32b = MapARGB;
  VP8LMapColor8b = MapAlpha;
  VP8LPredictors[11] = Predictor11;
  VP8LPredictors[12] = Predictor12;
  VP8LPredictors[13] = Predictor13;
  VP8LSubtractGreenFromBlueAndRed = SubtractGreenFromBlueAndRed;
  VP8LTransformColor = TransformColor;
#endif  // WEBP_USE_MIPS_DSP_R2
}

//------------------------------------------------------------------------------
