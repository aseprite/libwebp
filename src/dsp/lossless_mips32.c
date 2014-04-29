// Copyright 2014 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
// MIPS version of lossless functions
//
// Author(s):  Djordje Pesut    (djordje.pesut@imgtec.com)
//             Jovan Zelincevic (jovan.zelincevic@imgtec.com)

#include "./dsp.h"
#include "./lossless.h"

#if defined(WEBP_USE_MIPS32)

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define APPROX_LOG_WITH_CORRECTION_MAX  65536
#define APPROX_LOG_MAX                   4096
#define LOG_2_RECIPROCAL 1.44269504088896338700465094007086

static float FastSLog2Slow(uint32_t v) {
  assert(v >= LOG_LOOKUP_IDX_MAX);
  if (v < APPROX_LOG_WITH_CORRECTION_MAX) {
    uint32_t log_cnt, y, correction;
    const int c24 = 24;
    const float v_f = (float)v;
    uint32_t temp;

    // Xf = 256 = 2^8
    // log_cnt is index of leading one in upper 24 bits
    __asm__ volatile(
      "clz      %[log_cnt], %[v]                      \n\t"
      "addiu    %[y],       $zero,        1           \n\t"
      "subu     %[log_cnt], %[c24],       %[log_cnt]  \n\t"
      "sllv     %[y],       %[y],         %[log_cnt]  \n\t"
      "srlv     %[temp],    %[v],         %[log_cnt]  \n\t"
      : [log_cnt]"=&r"(log_cnt), [y]"=&r"(y),
        [temp]"=r"(temp)
      : [c24]"r"(c24), [v]"r"(v)
    );

    // vf = (2^log_cnt) * Xf; where y = 2^log_cnt and Xf < 256
    // Xf = floor(Xf) * (1 + (v % y) / v)
    // log2(Xf) = log2(floor(Xf)) + log2(1 + (v % y) / v)
    // The correction factor: log(1 + d) ~ d; for very small d values, so
    // log2(1 + (v % y) / v) ~ LOG_2_RECIPROCAL * (v % y)/v
    // LOG_2_RECIPROCAL ~ 23/16

    // (v % y) = (v % 2^log_cnt) = v & (2^log_cnt - 1)
    correction = (23 * (v & (y - 1))) >> 4;
    return v_f * (kLog2Table[temp] + log_cnt) + correction;
  } else {
    return (float)(LOG_2_RECIPROCAL * v * log((double)v));
  }
}

static float FastLog2Slow(uint32_t v) {
  assert(v >= LOG_LOOKUP_IDX_MAX);
  if (v < APPROX_LOG_WITH_CORRECTION_MAX) {
    uint32_t log_cnt, y;
    const int c24 = 24;
    double log_2;
    uint32_t temp;

    __asm__ volatile(
      "clz      %[log_cnt], %[v]                      \n\t"
      "addiu    %[y],       $zero,        1           \n\t"
      "subu     %[log_cnt], %[c24],       %[log_cnt]  \n\t"
      "sllv     %[y],       %[y],         %[log_cnt]  \n\t"
      "srlv     %[temp],    %[v],         %[log_cnt]  \n\t"
      : [log_cnt]"=&r"(log_cnt), [y]"=&r"(y),
        [temp]"=r"(temp)
      : [c24]"r"(c24), [v]"r"(v)
    );

    log_2 = kLog2Table[temp] + log_cnt;
    if (v >= APPROX_LOG_MAX) {
      // Since the division is still expensive, add this correction factor only
      // for large values of 'v'.

      const uint32_t correction = (23 * (v & (y - 1))) >> 4;
      log_2 += (double)correction / v;
    }
    return (float)log_2;
  } else {
    return (float)(LOG_2_RECIPROCAL * log((double)v));
  }
}

// C version of this function:
//   int i = 0;
//   int64_t cost = 0;
//   const uint32_t* pop = &population[4];
//   const uint32_t* LoopEnd = &population[length];
//   while (pop != LoopEnd) {
//     ++i;
//     cost += i * *pop;
//     cost += i * *(pop + 1);
//     pop += 2;
//   }
//   return (double)cost;
static double ExtraCost(const uint32_t* const population, int length) {
  int i, temp0, temp1;
  const uint32_t* pop = &population[4];
  const uint32_t* const LoopEnd = &population[length];

  __asm__ volatile(
    "mult   $zero,    $zero                  \n\t"
    "xor    %[i],     %[i],       %[i]       \n\t"
    "beq    %[pop],   %[LoopEnd], 2f         \n\t"
  "1:                                        \n\t"
    "lw     %[temp0], 0(%[pop])              \n\t"
    "lw     %[temp1], 4(%[pop])              \n\t"
    "addiu  %[i],     %[i],       1          \n\t"
    "addiu  %[pop],   %[pop],     8          \n\t"
    "madd   %[i],     %[temp0]               \n\t"
    "madd   %[i],     %[temp1]               \n\t"
    "bne    %[pop],   %[LoopEnd], 1b         \n\t"
  "2:                                        \n\t"
    "mfhi   %[temp0]                         \n\t"
    "mflo   %[temp1]                         \n\t"
    : [temp0]"=&r"(temp0), [temp1]"=&r"(temp1),
      [i]"=&r"(i), [pop]"+r"(pop)
    : [LoopEnd]"r"(LoopEnd)
    : "memory", "hi", "lo"
  );

  return (double)((int64_t)temp0 << 32 | temp1);
}

// C version of this function:
//   int i = 0;
//   int64_t cost = 0;
//   const uint32_t* pX = &X[4];
//   const uint32_t* pY = &Y[4];
//   const uint32_t* LoopEnd = &X[length];
//   while (pX != LoopEnd) {
//     const uint32_t xy0 = *pX + *pY;
//     const uint32_t xy1 = *(pX + 1) + *(pY + 1);
//     ++i;
//     cost += i * xy0;
//     cost += i * xy1;
//     pX += 2;
//     pY += 2;
//   }
//   return (double)cost;
static double ExtraCostCombined(const uint32_t* const X,
                                const uint32_t* const Y, int length) {
  int i, temp0, temp1, temp2, temp3;
  const uint32_t* pX = &X[4];
  const uint32_t* pY = &Y[4];
  const uint32_t* const LoopEnd = &X[length];

  __asm__ volatile(
    "mult   $zero,    $zero                  \n\t"
    "xor    %[i],     %[i],       %[i]       \n\t"
    "beq    %[pX],    %[LoopEnd], 2f         \n\t"
  "1:                                        \n\t"
    "lw     %[temp0], 0(%[pX])               \n\t"
    "lw     %[temp1], 0(%[pY])               \n\t"
    "lw     %[temp2], 4(%[pX])               \n\t"
    "lw     %[temp3], 4(%[pY])               \n\t"
    "addiu  %[i],     %[i],       1          \n\t"
    "addu   %[temp0], %[temp0],   %[temp1]   \n\t"
    "addu   %[temp2], %[temp2],   %[temp3]   \n\t"
    "addiu  %[pX],    %[pX],      8          \n\t"
    "addiu  %[pY],    %[pY],      8          \n\t"
    "madd   %[i],     %[temp0]               \n\t"
    "madd   %[i],     %[temp2]               \n\t"
    "bne    %[pX],    %[LoopEnd], 1b         \n\t"
  "2:                                        \n\t"
    "mfhi   %[temp0]                         \n\t"
    "mflo   %[temp1]                         \n\t"
    : [temp0]"=&r"(temp0), [temp1]"=&r"(temp1),
      [temp2]"=&r"(temp2), [temp3]"=&r"(temp3),
      [i]"=&r"(i), [pX]"+r"(pX), [pY]"+r"(pY)
    : [LoopEnd]"r"(LoopEnd)
    : "memory", "hi", "lo"
  );

  return (double)((int64_t)temp0 << 32 | temp1);
}

#define HUFFMAN_COST_PASS                                 \
  __asm__ volatile(                                       \
    "sll   %[temp1],  %[temp0],    3           \n\t"      \
    "addiu %[temp3],  %[streak],   -3          \n\t"      \
    "addu  %[temp2],  %[pstreaks], %[temp1]    \n\t"      \
    "blez  %[temp3],  1f                       \n\t"      \
    "srl   %[temp1],  %[temp1],    1           \n\t"      \
    "addu  %[temp3],  %[pcnts],    %[temp1]    \n\t"      \
    "lw    %[temp0],  4(%[temp2])              \n\t"      \
    "lw    %[temp1],  0(%[temp3])              \n\t"      \
    "addu  %[temp0],  %[temp0],    %[streak]   \n\t"      \
    "addiu %[temp1],  %[temp1],    1           \n\t"      \
    "sw    %[temp0],  4(%[temp2])              \n\t"      \
    "sw    %[temp1],  0(%[temp3])              \n\t"      \
    "b     2f                                  \n\t"      \
  "1:                                          \n\t"      \
    "lw    %[temp0],  0(%[temp2])              \n\t"      \
    "addu  %[temp0],  %[temp0],    %[streak]   \n\t"      \
    "sw    %[temp0],  0(%[temp2])              \n\t"      \
  "2:                                          \n\t"      \
    : [temp1]"=&r"(temp1), [temp2]"=&r"(temp2),           \
      [temp3]"=&r"(temp3), [temp0]"+r"(temp0)             \
    : [pstreaks]"r"(pstreaks), [pcnts]"r"(pcnts),         \
      [streak]"r"(streak)                                 \
    : "memory"                                            \
  );

// Returns the various RLE counts
static VP8LStreaks HuffmanCostCount(const uint32_t* population, int length) {
  int i;
  int streak = 0;
  VP8LStreaks stats;
  int* const pstreaks = &stats.streaks[0][0];
  int* const pcnts = &stats.counts[0];
  int temp0, temp1, temp2, temp3;
  memset(&stats, 0, sizeof(stats));
  for (i = 0; i < length - 1; ++i) {
    ++streak;
    if (population[i] == population[i + 1]) {
      continue;
    }
    temp0 = (population[i] != 0);
    HUFFMAN_COST_PASS
    streak = 0;
  }
  ++streak;
  temp0 = (population[i] != 0);
  HUFFMAN_COST_PASS

  return stats;
}

static VP8LStreaks HuffmanCostCombinedCount(const uint32_t* X,
                                            const uint32_t* Y, int length) {
  int i;
  int streak = 0;
  VP8LStreaks stats;
  int* const pstreaks = &stats.streaks[0][0];
  int* const pcnts = &stats.counts[0];
  int temp0, temp1, temp2, temp3;
  memset(&stats, 0, sizeof(stats));
  for (i = 0; i < length - 1; ++i) {
    const uint32_t xy = X[i] + Y[i];
    const uint32_t xy_next = X[i + 1] + Y[i + 1];
    ++streak;
    if (xy == xy_next) {
      continue;
    }
    temp0 = (xy != 0);
    HUFFMAN_COST_PASS
    streak = 0;
  }
  {
    const uint32_t xy = X[i] + Y[i];
    ++streak;
    temp0 = (xy != 0);
    HUFFMAN_COST_PASS
  }

  return stats;
}

#define ASM_START                                       \
  __asm__ volatile(                                     \
    ".set   push                            \n\t"       \
    ".set   at                              \n\t"       \
    ".set   macro                           \n\t"       \
  "1:                                       \n\t"

// P2 = P0 + P1
// A..D - offsets
// E and F - temp variables to tell macro
//           if pointers should be incremented
// literal_ and successive histograms could be unaligned
// so we must use ulw and usw
#define ADD_TO_OUT(A, B, C, D, E, F, P0, P1, P2)        \
    "ulw    %[temp0], "#A"(%["#P0"])        \n\t"       \
    "ulw    %[temp1], "#B"(%["#P0"])        \n\t"       \
    "ulw    %[temp2], "#C"(%["#P0"])        \n\t"       \
    "ulw    %[temp3], "#D"(%["#P0"])        \n\t"       \
    "ulw    %[temp4], "#A"(%["#P1"])        \n\t"       \
    "ulw    %[temp5], "#B"(%["#P1"])        \n\t"       \
    "ulw    %[temp6], "#C"(%["#P1"])        \n\t"       \
    "ulw    %[temp7], "#D"(%["#P1"])        \n\t"       \
    "addu   %[temp4], %[temp4],   %[temp0]  \n\t"       \
    "addu   %[temp5], %[temp5],   %[temp1]  \n\t"       \
    "addu   %[temp6], %[temp6],   %[temp2]  \n\t"       \
    "addu   %[temp7], %[temp7],   %[temp3]  \n\t"       \
  ".if "#E" == 1                            \n\t"       \
    "addiu  %["#P0"],  %["#P0"],  16        \n\t"       \
  ".endif                                   \n\t"       \
  ".if "#F" == 1                            \n\t"       \
    "addiu  %["#P1"],  %["#P1"],  16        \n\t"       \
  ".endif                                   \n\t"       \
    "usw    %[temp4], "#A"(%["#P2"])        \n\t"       \
    "usw    %[temp5], "#B"(%["#P2"])        \n\t"       \
    "usw    %[temp6], "#C"(%["#P2"])        \n\t"       \
    "usw    %[temp7], "#D"(%["#P2"])        \n\t"       \
  ".if "#E" == 1                            \n\t"       \
    "addiu  %["#P2"], %["#P2"],   16        \n\t"       \
    "bne    %["#P0"], %[LoopEnd], 1b        \n\t"       \
    ".set   pop                             \n\t"       \
  ".endif                                   \n\t"

#define ASM_END_0                                       \
    : [temp0]"=&r"(temp0), [temp1]"=&r"(temp1),         \
      [temp2]"=&r"(temp2), [temp3]"=&r"(temp3),         \
      [temp4]"=&r"(temp4), [temp5]"=&r"(temp5),         \
      [temp6]"=&r"(temp6), [temp7]"=&r"(temp7),         \
      [pa]"+r"(pa), [pb]"+r"(pb), [pout]"+r"(pout)      \
    : [LoopEnd]"r"(LoopEnd)                             \
    : "memory", "at"                                    \
  );

#define ASM_END_1                                       \
    : [temp0]"=&r"(temp0), [temp1]"=&r"(temp1),         \
      [temp2]"=&r"(temp2), [temp3]"=&r"(temp3),         \
      [temp4]"=&r"(temp4), [temp5]"=&r"(temp5),         \
      [temp6]"=&r"(temp6), [temp7]"=&r"(temp7),         \
      [pa]"+r"(pa), [pout]"+r"(pout)                    \
    : [LoopEnd]"r"(LoopEnd)                             \
    : "memory", "at"                                    \
  );

static void HistogramAdd(const VP8LHistogram* const a,
                         const VP8LHistogram* const b,
                         VP8LHistogram* const out) {
  uint32_t* pout;
  uint32_t* pa;
  uint32_t* pb;
  uint32_t* LoopEnd;
  uint32_t temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, i;

  assert(a->palette_code_bits_ == b->palette_code_bits_);

  if (b != out) {
    pa = (uint32_t*)a->literal_;
    pb = (uint32_t*)b->literal_;
    pout = out->literal_;
    LoopEnd = pa + NUM_LITERAL_CODES + NUM_LENGTH_CODES;
    // (NUM_LITERAL_CODES + NUM_LENGTH_CODES) % 4 = 0
    ASM_START
    ADD_TO_OUT(0, 4, 8, 12, 1, 1, pa, pb, pout)
    ASM_END_0
    // remaining literals, if there is any
    temp0 = (a->palette_code_bits_ > 0) ? (1 << a->palette_code_bits_) : 0;
    for (i = 0; i < temp0; ++i) {
      *pout++ = *pa++ + *pb++;
    }

    pa = (uint32_t*)a->distance_;
    pb = (uint32_t*)b->distance_;
    pout = out->distance_;
    LoopEnd = pa + NUM_DISTANCE_CODES;
    // NUM_DISTANCE_CODES % 4 = 0
    ASM_START
    ADD_TO_OUT(0, 4, 8, 12, 1, 1, pa, pb, pout)
    ASM_END_0

    pa = (uint32_t*)a->red_;
    pb = (uint32_t*)b->red_;
    pout = out->red_;
    LoopEnd = pa + NUM_LITERAL_CODES;
    // works only if 'uint32_t red_[256]', 'uint32_t blue_[256]'
    // and 'uint32_t alpha_[256]' are successive in memory
    // &blue_[0] == &red_[256]
    // &alpha_[0] == &red_[512]
    // NUM_LITERAL_CODES % 4 = 0
    ASM_START
    ADD_TO_OUT(   0,    4,    8,   12, 0, 0, pa, pb, pout)
    ADD_TO_OUT(1024, 1028, 1032, 1036, 0, 0, pa, pb, pout)
    ADD_TO_OUT(2048, 2052, 2056, 2060, 1, 1, pa, pb, pout)
    ASM_END_0
  } else {
    pa = (uint32_t*)a->literal_;
    pout = out->literal_;
    LoopEnd = pa + NUM_LITERAL_CODES + NUM_LENGTH_CODES;
    // (NUM_LITERAL_CODES + NUM_LENGTH_CODES) % 4 = 0
    ASM_START
    ADD_TO_OUT(0, 4, 8, 12, 1, 0, pa, pout, pout)
    ASM_END_1
    // remaining literals, if there is any
    temp0 = (a->palette_code_bits_ > 0) ? (1 << a->palette_code_bits_) : 0;
    for (i = 0; i < temp0; ++i) {
      *pout++ += *pa++;
    }

    pa = (uint32_t*)a->distance_;
    pout = out->distance_;
    LoopEnd = pa + NUM_DISTANCE_CODES;
    // NUM_DISTANCE_CODES % 4 = 0
    ASM_START
    ADD_TO_OUT(0, 4, 8, 12, 1, 0, pa, pout, pout)
    ASM_END_1

    pa = (uint32_t*)a->red_;
    pout = out->red_;
    LoopEnd = pa + NUM_LITERAL_CODES;
    // works only if 'int red_[256]', 'int blue_[256]' and 'int alpha_[256]'
    // are successive in memory
    // &blue_[0] == &red_[256]
    // &alpha_[0] == &red_[512]
    // NUM_LITERAL_CODES % 4 = 0
    ASM_START
    ADD_TO_OUT(   0,    4,    8,   12, 0, 0, pa, pout, pout)
    ADD_TO_OUT(1024, 1028, 1032, 1036, 0, 0, pa, pout, pout)
    ADD_TO_OUT(2048, 2052, 2056, 2060, 1, 0, pa, pout, pout)
    ASM_END_1
  }
}

#undef ASM_END_1
#undef ASM_END_0
#undef ADD_TO_OUT
#undef ASM_START

#endif  // WEBP_USE_MIPS32

//------------------------------------------------------------------------------
// Entry point

extern void VP8LDspInitMIPS32(void);

void VP8LDspInitMIPS32(void) {
#if defined(WEBP_USE_MIPS32)
  VP8LFastSLog2Slow = FastSLog2Slow;
  VP8LFastLog2Slow = FastLog2Slow;
  VP8LExtraCost = ExtraCost;
  VP8LExtraCostCombined = ExtraCostCombined;
  VP8LHuffmanCostCount = HuffmanCostCount;
  VP8LHuffmanCostCombinedCount = HuffmanCostCombinedCount;
  VP8LHistogramAdd = HistogramAdd;
#endif  // WEBP_USE_MIPS32
}
