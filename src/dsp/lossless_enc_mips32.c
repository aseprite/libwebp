// Copyright 2015 Google Inc. All Rights Reserved.
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

static WEBP_INLINE double BitsEntropyRefine(int nonzeros, int sum, int max_val,
                                            double retval) {
  double mix;
  if (nonzeros < 5) {
    if (nonzeros <= 1) {
      return 0;
    }
    // Two symbols, they will be 0 and 1 in a Huffman code.
    // Let's mix in a bit of entropy to favor good clustering when
    // distributions of these are combined.
    if (nonzeros == 2) {
      return 0.99 * sum + 0.01 * retval;
    }
    // No matter what the entropy says, we cannot be better than min_limit
    // with Huffman coding. I am mixing a bit of entropy into the
    // min_limit since it produces much better (~0.5 %) compression results
    // perhaps because of better entropy clustering.
    if (nonzeros == 3) {
      mix = 0.95;
    } else {
      mix = 0.7;  // nonzeros == 4.
    }
  } else {
    mix = 0.627;
  }

  {
    double min_limit = 2 * sum - max_val;
    min_limit = mix * min_limit + (1.0 - mix) * retval;
    return (retval < min_limit) ? min_limit : retval;
  }
}

static double InitialHuffmanCost(void) {
  // Small bias because Huffman code length is typically not stored in
  // full length.
  static const int kHuffmanCodeOfHuffmanCodeSize = CODE_LENGTH_CODES * 3;
  static const double kSmallBias = 9.1;
  return kHuffmanCodeOfHuffmanCodeSize - kSmallBias;
}

// Finalize the Huffman cost based on streak numbers and length type (<3 or >=3)
static double FinalHuffmanCost(const VP8LStreaks* const stats) {
  double retval = InitialHuffmanCost();
  retval += stats->counts[0] * 1.5625 + 0.234375 * stats->streaks[0][1];
  retval += stats->counts[1] * 2.578125 + 0.703125 * stats->streaks[1][1];
  retval += 1.796875 * stats->streaks[0][0];
  retval += 3.28125 * stats->streaks[1][0];
  return retval;
}

static double GetCombinedEntropy(const uint32_t* const X,
                                 const uint32_t* const Y, int length) {
  double bits_entropy_combined;
  double huffman_cost_combined;
  int i;

  // Bit entropy variables.
  double retval = 0.;
  int sum = 0;
  int nonzeros = 0;
  uint32_t max_val = 0;
  int i_prev;
  uint32_t xy;

  // Huffman cost variables.
  int streak = 0;
  uint32_t xy_prev;
  VP8LStreaks stats;

  // Macro variables.
  int* const pstreaks = &stats.streaks[0][0];
  int* const pcnts = &stats.counts[0];
  int temp0, temp1, temp2, temp3;

  memset(&stats, 0, sizeof(stats));

  // Treat the first value for the huffman cost: this is keeping the original
  // behavior, even though there is no first streak.
  // TODO(vrabaud): study proper behavior
  xy = X[0] + Y[0];
  ++stats.streaks[xy != 0][0];
  xy_prev = xy;
  i_prev = 0;

  for (i = 1; i < length; ++i) {
    xy = X[i] + Y[i];

    // Process data by streaks for both bit entropy and huffman cost.
    if (xy != xy_prev) {
      streak = i - i_prev;

      // Gather info for the bit entropy.
      if (xy_prev != 0) {
        sum += xy_prev * streak;
        nonzeros += streak;
        retval -= VP8LFastSLog2(xy_prev) * streak;
        if (max_val < xy_prev) {
          max_val = xy_prev;
        }
      }

      // Gather info for the huffman cost.
      temp0 = (xy != 0);
      HUFFMAN_COST_PASS

      xy_prev = xy;
      i_prev = i;
    }
  }

  // Finish off the last streak for bit entropy.
  if (xy != 0) {
    streak = i - i_prev;
    sum += xy * streak;
    nonzeros += streak;
    retval -= VP8LFastSLog2(xy) * streak;
    if (max_val < xy) {
      max_val = xy;
    }
  }
  // Huffman cost is not updated with the last streak to keep original behavior.
  // TODO(vrabaud): study proper behavior

  retval += VP8LFastSLog2(sum);
  bits_entropy_combined = BitsEntropyRefine(nonzeros, sum, max_val, retval);

  huffman_cost_combined = FinalHuffmanCost(&stats);

  return bits_entropy_combined + huffman_cost_combined;
}

#define ASM_START                                       \
  __asm__ volatile(                                     \
    ".set   push                            \n\t"       \
    ".set   at                              \n\t"       \
    ".set   macro                           \n\t"       \
  "1:                                       \n\t"

// P2 = P0 + P1
// A..D - offsets
// E - temp variable to tell macro
//     if pointer should be incremented
// literal_ and successive histograms could be unaligned
// so we must use ulw and usw
#define ADD_TO_OUT(A, B, C, D, E, P0, P1, P2)           \
    "ulw    %[temp0], " #A "(%[" #P0 "])    \n\t"       \
    "ulw    %[temp1], " #B "(%[" #P0 "])    \n\t"       \
    "ulw    %[temp2], " #C "(%[" #P0 "])    \n\t"       \
    "ulw    %[temp3], " #D "(%[" #P0 "])    \n\t"       \
    "ulw    %[temp4], " #A "(%[" #P1 "])    \n\t"       \
    "ulw    %[temp5], " #B "(%[" #P1 "])    \n\t"       \
    "ulw    %[temp6], " #C "(%[" #P1 "])    \n\t"       \
    "ulw    %[temp7], " #D "(%[" #P1 "])    \n\t"       \
    "addu   %[temp4], %[temp4],   %[temp0]  \n\t"       \
    "addu   %[temp5], %[temp5],   %[temp1]  \n\t"       \
    "addu   %[temp6], %[temp6],   %[temp2]  \n\t"       \
    "addu   %[temp7], %[temp7],   %[temp3]  \n\t"       \
    "addiu  %[" #P0 "],  %[" #P0 "],  16    \n\t"       \
  ".if " #E " == 1                          \n\t"       \
    "addiu  %[" #P1 "],  %[" #P1 "],  16    \n\t"       \
  ".endif                                   \n\t"       \
    "usw    %[temp4], " #A "(%[" #P2 "])    \n\t"       \
    "usw    %[temp5], " #B "(%[" #P2 "])    \n\t"       \
    "usw    %[temp6], " #C "(%[" #P2 "])    \n\t"       \
    "usw    %[temp7], " #D "(%[" #P2 "])    \n\t"       \
    "addiu  %[" #P2 "], %[" #P2 "],   16    \n\t"       \
    "bne    %[" #P0 "], %[LoopEnd], 1b      \n\t"       \
    ".set   pop                             \n\t"       \

#define ASM_END_COMMON_0                                \
    : [temp0]"=&r"(temp0), [temp1]"=&r"(temp1),         \
      [temp2]"=&r"(temp2), [temp3]"=&r"(temp3),         \
      [temp4]"=&r"(temp4), [temp5]"=&r"(temp5),         \
      [temp6]"=&r"(temp6), [temp7]"=&r"(temp7),         \
      [pa]"+r"(pa), [pout]"+r"(pout)

#define ASM_END_COMMON_1                                \
    : [LoopEnd]"r"(LoopEnd)                             \
    : "memory", "at"                                    \
  );

#define ASM_END_0                                       \
    ASM_END_COMMON_0                                    \
      , [pb]"+r"(pb)                                    \
    ASM_END_COMMON_1

#define ASM_END_1                                       \
    ASM_END_COMMON_0                                    \
    ASM_END_COMMON_1

#define ADD_VECTOR(A, B, OUT, SIZE, EXTRA_SIZE)  do {   \
  const uint32_t* pa = (const uint32_t*)(A);            \
  const uint32_t* pb = (const uint32_t*)(B);            \
  uint32_t* pout = (uint32_t*)(OUT);                    \
  const uint32_t* const LoopEnd = pa + (SIZE);          \
  assert((SIZE) % 4 == 0);                              \
  ASM_START                                             \
  ADD_TO_OUT(0, 4, 8, 12, 1, pa, pb, pout)              \
  ASM_END_0                                             \
  if ((EXTRA_SIZE) > 0) {                               \
    const int last = (EXTRA_SIZE);                      \
    int i;                                              \
    for (i = 0; i < last; ++i) pout[i] = pa[i] + pb[i]; \
  }                                                     \
} while (0)

#define ADD_VECTOR_EQ(A, OUT, SIZE, EXTRA_SIZE)  do {   \
  const uint32_t* pa = (const uint32_t*)(A);            \
  uint32_t* pout = (uint32_t*)(OUT);                    \
  const uint32_t* const LoopEnd = pa + (SIZE);          \
  assert((SIZE) % 4 == 0);                              \
  ASM_START                                             \
  ADD_TO_OUT(0, 4, 8, 12, 0, pa, pout, pout)            \
  ASM_END_1                                             \
  if ((EXTRA_SIZE) > 0) {                               \
    const int last = (EXTRA_SIZE);                      \
    int i;                                              \
    for (i = 0; i < last; ++i) pout[i] += pa[i];        \
  }                                                     \
} while (0)

static void HistogramAdd(const VP8LHistogram* const a,
                         const VP8LHistogram* const b,
                         VP8LHistogram* const out) {
  uint32_t temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7;
  const int extra_cache_size = VP8LHistogramNumCodes(a->palette_code_bits_)
                             - (NUM_LITERAL_CODES + NUM_LENGTH_CODES);
  assert(a->palette_code_bits_ == b->palette_code_bits_);

  if (b != out) {
    ADD_VECTOR(a->literal_, b->literal_, out->literal_,
               NUM_LITERAL_CODES + NUM_LENGTH_CODES, extra_cache_size);
    ADD_VECTOR(a->distance_, b->distance_, out->distance_,
               NUM_DISTANCE_CODES, 0);
    ADD_VECTOR(a->red_, b->red_, out->red_, NUM_LITERAL_CODES, 0);
    ADD_VECTOR(a->blue_, b->blue_, out->blue_, NUM_LITERAL_CODES, 0);
    ADD_VECTOR(a->alpha_, b->alpha_, out->alpha_, NUM_LITERAL_CODES, 0);
  } else {
    ADD_VECTOR_EQ(a->literal_, out->literal_,
                  NUM_LITERAL_CODES + NUM_LENGTH_CODES, extra_cache_size);
    ADD_VECTOR_EQ(a->distance_, out->distance_, NUM_DISTANCE_CODES, 0);
    ADD_VECTOR_EQ(a->red_, out->red_, NUM_LITERAL_CODES, 0);
    ADD_VECTOR_EQ(a->blue_, out->blue_, NUM_LITERAL_CODES, 0);
    ADD_VECTOR_EQ(a->alpha_, out->alpha_, NUM_LITERAL_CODES, 0);
  }
}

#undef ADD_VECTOR_EQ
#undef ADD_VECTOR
#undef ASM_END_1
#undef ASM_END_0
#undef ASM_END_COMMON_1
#undef ASM_END_COMMON_0
#undef ADD_TO_OUT
#undef ASM_START

//------------------------------------------------------------------------------
// Entry point

extern void VP8LEncDspInitMIPS32(void);

WEBP_TSAN_IGNORE_FUNCTION void VP8LEncDspInitMIPS32(void) {
  VP8LFastSLog2Slow = FastSLog2Slow;
  VP8LFastLog2Slow = FastLog2Slow;
  VP8LExtraCost = ExtraCost;
  VP8LExtraCostCombined = ExtraCostCombined;
  VP8LHuffmanCostCount = HuffmanCostCount;
  VP8LGetCombinedEntropy = GetCombinedEntropy;
  VP8LHistogramAdd = HistogramAdd;
}

#else  // !WEBP_USE_MIPS32

WEBP_DSP_INIT_STUB(VP8LEncDspInitMIPS32)

#endif  // WEBP_USE_MIPS32
