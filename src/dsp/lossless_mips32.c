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
// Author(s):  Jovan Zelincevic (jovan.zelincevic@imgtec.com)
//             Djordje Pesut    (djordje.pesut@imgtec.com)

#include "./dsp.h"
#include "./lossless.h"
#include "../enc/histogram.h"

#if defined(WEBP_USE_MIPS32)

#include <math.h>
#include <stdlib.h>
#include <assert.h>

#define APPROX_LOG_WITH_CORRECTION_MAX  65536
#define APPROX_LOG_MAX                   4096
#define LOG_2_RECIPROCAL 1.44269504088896338700465094007086

static float FastSLog2SlowMIPS32(int v) {
  assert(v >= LOG_LOOKUP_IDX_MAX);
  if (v < APPROX_LOG_WITH_CORRECTION_MAX) {
    int log_cnt, y, correction;
    const int c24 = 24;
    const float v_f = (float)v;
    int temp;

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

static float FastLog2SlowMIPS32(int v) {
  assert(v >= LOG_LOOKUP_IDX_MAX);
  if (v < APPROX_LOG_WITH_CORRECTION_MAX) {
    int log_cnt, y;
    const int c24 = 24;
    double log_2;
    int temp;

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

      const int correction = (23 * (v & (y - 1))) >> 4;
      log_2 += (double)correction / v;
    }
    return (float)log_2;
  } else {
    return (float)(LOG_2_RECIPROCAL * log((double)v));
  }
}

static WEBP_INLINE double InitialHuffmanCost(void) {
  // Small bias because Huffman code length is typically not stored in
  // full length.
  static const int kHuffmanCodeOfHuffmanCodeSize = CODE_LENGTH_CODES * 3;
  static const double kSmallBias = 9.1;
  return kHuffmanCodeOfHuffmanCodeSize - kSmallBias;
}

// Returns the cost encode the rle-encoded entropy code.
// The constants in this function are experimental.
static double HuffmanCostMIPS32(const int* const population, int length) {
  double retval = InitialHuffmanCost();
  int streak = 0;
  int i = 0;
  int streak1 = 0, streak2 = 0, streak3 = 0, streak4 = 0;
  int cnt1 = 0, cnt2 = 0;
  int temp0;
  for (; i < length - 1; ++i) {
    ++streak;
    if (population[i] == population[i + 1]) {
      continue;
    }
    // population[i] points now to the symbol in the streak of same values.
    temp0 = population[i];
    if (streak > 3) {
      int temp1, temp2, temp3, temp4;
      __asm__ volatile(
        "addu  %[temp1],   %[streak1], %[streak]    \n\t"
        "addu  %[temp2],   %[streak2], %[streak]    \n\t"
        "addiu %[temp3],   %[cnt1],    1            \n\t"
        "addiu %[temp4],   %[cnt2],    1            \n\t"
        "movz  %[streak1], %[temp1],   %[temp0]     \n\t"
        "movn  %[streak2], %[temp2],   %[temp0]     \n\t"
        "movz  %[cnt1],    %[temp3],   %[temp0]     \n\t"
        "movn  %[cnt2],    %[temp4],   %[temp0]     \n\t"
        : [streak1]"+r"(streak1), [streak2]"+r"(streak2),
          [temp1]"=&r"(temp1), [temp2]"=&r"(temp2),
          [temp3]"=&r"(temp3), [temp4]"=&r"(temp4),
          [cnt1]"+r"(cnt1), [cnt2]"+r"(cnt2)
        : [streak]"r"(streak), [temp0]"r"(temp0)
      );
    } else {
      int temp1, temp2;
      __asm__ volatile(
        "addu  %[temp1],   %[streak3], %[streak]    \n\t"
        "addu  %[temp2],   %[streak4], %[streak]    \n\t"
        "movz  %[streak3], %[temp1],   %[temp0]     \n\t"
        "movn  %[streak4], %[temp2],   %[temp0]     \n\t"
        : [streak3]"+r"(streak3), [streak4]"+r"(streak4),
          [temp1]"=&r"(temp1), [temp2]"=&r"(temp2)
        : [streak]"r"(streak), [temp0]"r"(temp0)
      );
    }
    streak = 0;
  }

  ++streak;
  temp0 = population[i];
  if (streak > 3) {
    int temp1, temp2, temp3, temp4;
    __asm__ volatile(
      "addu  %[temp1],   %[streak1], %[streak]    \n\t"
      "addu  %[temp2],   %[streak2], %[streak]    \n\t"
      "addiu %[temp3],   %[cnt1],    1            \n\t"
      "addiu %[temp4],   %[cnt2],    1            \n\t"
      "movz  %[streak1], %[temp1],   %[temp0]     \n\t"
      "movn  %[streak2], %[temp2],   %[temp0]     \n\t"
      "movz  %[cnt1],    %[temp3],   %[temp0]     \n\t"
      "movn  %[cnt2],    %[temp4],   %[temp0]     \n\t"
      : [streak1]"+r"(streak1), [streak2]"+r"(streak2),
        [temp1]"=&r"(temp1), [temp2]"=&r"(temp2),
        [temp3]"=&r"(temp3), [temp4]"=&r"(temp4),
        [cnt1]"+r"(cnt1), [cnt2]"+r"(cnt2)
      : [streak]"r"(streak), [temp0]"r"(temp0)
    );
  } else {
    int temp1, temp2;
    __asm__ volatile(
      "addu  %[temp1],   %[streak3], %[streak]    \n\t"
      "addu  %[temp2],   %[streak4], %[streak]    \n\t"
      "movz  %[streak3], %[temp1],   %[temp0]     \n\t"
      "movn  %[streak4], %[temp2],   %[temp0]     \n\t"
      : [streak3]"+r"(streak3), [streak4]"+r"(streak4),
        [temp1]"=&r"(temp1), [temp2]"=&r"(temp2)
      : [streak]"r"(streak), [temp0]"r"(temp0)
    );
  }

  // float/double operations moved out of loop
  retval += cnt1 * 1.5625 + 0.234375 * streak1;
  retval += cnt2 * 2.578125 + 0.703125 * streak2;
  retval += 1.796875 * streak3;
  retval += 3.28125 * streak4;

  return retval;
}

static double HuffmanCostCombinedMIPS32(const int* const X, const int* const Y,
                                        int length) {
  double retval = InitialHuffmanCost();
  int streak = 0;
  int i = 0;
  int streak1 = 0, streak2 = 0, streak3 = 0, streak4 = 0;
  int cnt1 = 0, cnt2 = 0;
  int temp0;

  for (; i < length - 1; ++i) {
    const int xy = X[i] + Y[i];
    const int xy_next = X[i + 1] + Y[i + 1];
    ++streak;
    if (xy == xy_next) {
      continue;
    }
    temp0 = xy;
    if (streak > 3) {
      int temp1, temp2, temp3, temp4;
      __asm__ volatile(
        "addu  %[temp1],   %[streak1], %[streak]    \n\t"
        "addu  %[temp2],   %[streak2], %[streak]    \n\t"
        "addiu %[temp3],   %[cnt1],    1            \n\t"
        "addiu %[temp4],   %[cnt2],    1            \n\t"
        "movz  %[streak1], %[temp1],   %[temp0]     \n\t"
        "movn  %[streak2], %[temp2],   %[temp0]     \n\t"
        "movz  %[cnt1],    %[temp3],   %[temp0]     \n\t"
        "movn  %[cnt2],    %[temp4],   %[temp0]     \n\t"
        : [streak1]"+r"(streak1), [streak2]"+r"(streak2),
          [temp1]"=&r"(temp1), [temp2]"=&r"(temp2),
          [temp3]"=&r"(temp3), [temp4]"=&r"(temp4),
          [cnt1]"+r"(cnt1), [cnt2]"+r"(cnt2)
        : [streak]"r"(streak), [temp0]"r"(temp0)
      );
    } else {
      int temp1, temp2;
      __asm__ volatile(
        "addu  %[temp1],   %[streak3], %[streak]    \n\t"
        "addu  %[temp2],   %[streak4], %[streak]    \n\t"
        "movz  %[streak3], %[temp1],   %[temp0]     \n\t"
        "movn  %[streak4], %[temp2],   %[temp0]     \n\t"
        : [streak3]"+r"(streak3), [streak4]"+r"(streak4),
          [temp1]"=&r"(temp1), [temp2]"=&r"(temp2)
        : [streak]"r"(streak), [temp0]"r"(temp0)
      );
    }
    streak = 0;
  }

  ++streak;
  temp0 = X[i] + Y[i];
  if (streak > 3) {
    int temp1, temp2, temp3, temp4;
    __asm__ volatile(
      "addu  %[temp1],   %[streak1], %[streak]    \n\t"
      "addu  %[temp2],   %[streak2], %[streak]    \n\t"
      "addiu %[temp3],   %[cnt1],    1            \n\t"
      "addiu %[temp4],   %[cnt2],    1            \n\t"
      "movz  %[streak1], %[temp1],   %[temp0]     \n\t"
      "movn  %[streak2], %[temp2],   %[temp0]     \n\t"
      "movz  %[cnt1],    %[temp3],   %[temp0]     \n\t"
      "movn  %[cnt2],    %[temp4],   %[temp0]     \n\t"
      : [streak1]"+r"(streak1), [streak2]"+r"(streak2),
        [temp1]"=&r"(temp1), [temp2]"=&r"(temp2),
        [temp3]"=&r"(temp3), [temp4]"=&r"(temp4),
        [cnt1]"+r"(cnt1), [cnt2]"+r"(cnt2)
      : [streak]"r"(streak), [temp0]"r"(temp0)
    );
  } else {
    int temp1, temp2;
    __asm__ volatile(
      "addu  %[temp1],   %[streak3], %[streak]    \n\t"
      "addu  %[temp2],   %[streak4], %[streak]    \n\t"
      "movz  %[streak3], %[temp1],   %[temp0]     \n\t"
      "movn  %[streak4], %[temp2],   %[temp0]     \n\t"
      : [streak3]"+r"(streak3), [streak4]"+r"(streak4),
        [temp1]"=&r"(temp1), [temp2]"=&r"(temp2)
      : [streak]"r"(streak), [temp0]"r"(temp0)
    );
  }

  // float/double operations moved out of loop
  retval += cnt1 * 1.5625 + 0.234375 * streak1;
  retval += cnt2 * 2.578125 + 0.703125 * streak2;
  retval += 1.796875 * streak3;
  retval += 3.28125 * streak4;

  return retval;
}

#endif  // WEBP_USE_MIPS32

//------------------------------------------------------------------------------
// Entry point

extern void VP8LDspInitMIPS32(void);
extern void HistogramInitMIPS32(void);

void VP8LDspInitMIPS32(void) {
#if defined(WEBP_USE_MIPS32)
  VP8LFastSLog2Slow = FastSLog2SlowMIPS32;
  VP8LFastLog2Slow = FastLog2SlowMIPS32;
#endif  // WEBP_USE_MIPS32
}

void HistogramInitMIPS32(void) {
#if defined(WEBP_USE_MIPS32)
  HuffmanCost = HuffmanCostMIPS32;
  HuffmanCostCombined = HuffmanCostCombinedMIPS32;
#endif  // WEBP_USE_MIPS32
}
