// Copyright 2014 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
// MIPS version of speed-critical encoding functions.
//
// Author(s): Djordje Pesut    (djordje.pesut@imgtec.com)
//            Jovan Zelincevic (jovan.zelincevic@imgtec.com)

#include "./dsp.h"

#if defined(WEBP_USE_MIPS32)

#include "../enc/vp8enci.h"

static const int kC1 = 20091 + (1 << 16);
static const int kC2 = 35468;

// macro for one vertical pass in ITransformOne
// MUL macro inlined
// temp0..temp15 holds tmp[0]..tmp[15]
// A..D - offsets in bytes to load from in buffer
// TEMP0..TEMP3 - registers for corresponding tmp elements
// TEMP4..TEMP5 - temporary registers
#define VERTICAL_PASS(A, B, C, D, TEMP4, TEMP5, TEMP0, TEMP1, TEMP2, TEMP3) \
  "lh      %[temp16],      "#A"(%[temp20])                 \n\t"            \
  "lh      %[temp18],      "#B"(%[temp20])                 \n\t"            \
  "lh      %[temp17],      "#C"(%[temp20])                 \n\t"            \
  "lh      %[temp19],      "#D"(%[temp20])                 \n\t"            \
  "addu    %["#TEMP4"],    %[temp16],      %[temp18]       \n\t"            \
  "subu    %["#TEMP5"],    %[temp16],      %[temp18]       \n\t"            \
  "mul     %[temp16],      %[temp17],      %[kC2]          \n\t"            \
  "mul     %[temp18],      %[temp19],      %[kC1]          \n\t"            \
  "mul     %[temp17],      %[temp17],      %[kC1]          \n\t"            \
  "mul     %[temp19],      %[temp19],      %[kC2]          \n\t"            \
  "sra     %[temp16],      %[temp16],      16              \n\n"            \
  "sra     %[temp18],      %[temp18],      16              \n\n"            \
  "sra     %[temp17],      %[temp17],      16              \n\n"            \
  "sra     %[temp19],      %[temp19],      16              \n\n"            \
  "subu    %["#TEMP2"],    %[temp16],      %[temp18]       \n\t"            \
  "addu    %["#TEMP3"],    %[temp17],      %[temp19]       \n\t"            \
  "addu    %["#TEMP0"],    %["#TEMP4"],    %["#TEMP3"]     \n\t"            \
  "addu    %["#TEMP1"],    %["#TEMP5"],    %["#TEMP2"]     \n\t"            \
  "subu    %["#TEMP2"],    %["#TEMP5"],    %["#TEMP2"]     \n\t"            \
  "subu    %["#TEMP3"],    %["#TEMP4"],    %["#TEMP3"]     \n\t"

// macro for one horizontal pass in ITransformOne
// MUL and STORE macros inlined
// a = clip_8b(a) is replaced with: a = max(a, 0); a = min(a, 255)
// temp0..temp15 holds tmp[0]..tmp[15]
// A..D - offsets in bytes to load from ref and store to dst buffer
// TEMP0, TEMP4, TEMP8 and TEMP12 - registers for corresponding tmp elements
#define HORIZONTAL_PASS(A, B, C, D, TEMP0, TEMP4, TEMP8, TEMP12)            \
  "addiu   %["#TEMP0"],    %["#TEMP0"],    4               \n\t"            \
  "addu    %[temp16],      %["#TEMP0"],    %["#TEMP8"]     \n\t"            \
  "subu    %[temp17],      %["#TEMP0"],    %["#TEMP8"]     \n\t"            \
  "mul     %["#TEMP0"],    %["#TEMP4"],    %[kC2]          \n\t"            \
  "mul     %["#TEMP8"],    %["#TEMP12"],   %[kC1]          \n\t"            \
  "mul     %["#TEMP4"],    %["#TEMP4"],    %[kC1]          \n\t"            \
  "mul     %["#TEMP12"],   %["#TEMP12"],   %[kC2]          \n\t"            \
  "sra     %["#TEMP0"],    %["#TEMP0"],    16              \n\t"            \
  "sra     %["#TEMP8"],    %["#TEMP8"],    16              \n\t"            \
  "sra     %["#TEMP4"],    %["#TEMP4"],    16              \n\t"            \
  "sra     %["#TEMP12"],   %["#TEMP12"],   16              \n\t"            \
  "subu    %[temp18],      %["#TEMP0"],    %["#TEMP8"]     \n\t"            \
  "addu    %[temp19],      %["#TEMP4"],    %["#TEMP12"]    \n\t"            \
  "addu    %["#TEMP0"],    %[temp16],      %[temp19]       \n\t"            \
  "addu    %["#TEMP4"],    %[temp17],      %[temp18]       \n\t"            \
  "subu    %["#TEMP8"],    %[temp17],      %[temp18]       \n\t"            \
  "subu    %["#TEMP12"],   %[temp16],      %[temp19]       \n\t"            \
  "lw      %[temp20],      0(%[args])                      \n\t"            \
  "sra     %["#TEMP0"],    %["#TEMP0"],    3               \n\t"            \
  "sra     %["#TEMP4"],    %["#TEMP4"],    3               \n\t"            \
  "sra     %["#TEMP8"],    %["#TEMP8"],    3               \n\t"            \
  "sra     %["#TEMP12"],   %["#TEMP12"],   3               \n\t"            \
  "lbu     %[temp16],      "#A"(%[temp20])                 \n\t"            \
  "lbu     %[temp17],      "#B"(%[temp20])                 \n\t"            \
  "lbu     %[temp18],      "#C"(%[temp20])                 \n\t"            \
  "lbu     %[temp19],      "#D"(%[temp20])                 \n\t"            \
  "addu    %["#TEMP0"],    %[temp16],      %["#TEMP0"]     \n\t"            \
  "addu    %["#TEMP4"],    %[temp17],      %["#TEMP4"]     \n\t"            \
  "addu    %["#TEMP8"],    %[temp18],      %["#TEMP8"]     \n\t"            \
  "addu    %["#TEMP12"],   %[temp19],      %["#TEMP12"]    \n\t"            \
  "slt     %[temp16],      %["#TEMP0"],    $zero           \n\t"            \
  "slt     %[temp17],      %["#TEMP4"],    $zero           \n\t"            \
  "slt     %[temp18],      %["#TEMP8"],    $zero           \n\t"            \
  "slt     %[temp19],      %["#TEMP12"],   $zero           \n\t"            \
  "movn    %["#TEMP0"],    $zero,          %[temp16]       \n\t"            \
  "movn    %["#TEMP4"],    $zero,          %[temp17]       \n\t"            \
  "movn    %["#TEMP8"],    $zero,          %[temp18]       \n\t"            \
  "movn    %["#TEMP12"],   $zero,          %[temp19]       \n\t"            \
  "slt     %[temp16],      %["#TEMP0"],    %[temp21]       \n\t"            \
  "slt     %[temp17],      %["#TEMP4"],    %[temp21]       \n\t"            \
  "slt     %[temp18],      %["#TEMP8"],    %[temp21]       \n\t"            \
  "slt     %[temp19],      %["#TEMP12"],   %[temp21]       \n\t"            \
  "lw      %[temp20],      8(%[args])                      \n\t"            \
  "movz    %["#TEMP0"],    %[temp21],      %[temp16]       \n\t"            \
  "movz    %["#TEMP4"],    %[temp21],      %[temp17]       \n\t"            \
  "movz    %["#TEMP8"],    %[temp21],      %[temp18]       \n\t"            \
  "movz    %["#TEMP12"],   %[temp21],      %[temp19]       \n\t"            \
  "sb      %["#TEMP0"],    "#A"(%[temp20])                 \n\t"            \
  "sb      %["#TEMP4"],    "#B"(%[temp20])                 \n\t"            \
  "sb      %["#TEMP8"],    "#C"(%[temp20])                 \n\t"            \
  "sb      %["#TEMP12"],   "#D"(%[temp20])                 \n\t"

// Does one or two inverse transforms.
static WEBP_INLINE void ITransformOneMIPS32(const uint8_t* ref,
                                            const int16_t* in,
                                            uint8_t* dst) {
  int temp0, temp1, temp2, temp3, temp4, temp5, temp6;
  int temp7, temp8, temp9, temp10, temp11, temp12, temp13;
  int temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21;
  const int* args[3] = {(const int*)ref, (const int*)in, (const int*)dst};

  __asm__ volatile(
    "lw      %[temp20],      4(%[args])                      \n\t"
    VERTICAL_PASS(0, 16,  8, 24, temp4,  temp5,  temp0,  temp1,  temp2,  temp3)
    VERTICAL_PASS(2, 18, 10, 26, temp8,  temp9,  temp4,  temp5,  temp6,  temp7)
    VERTICAL_PASS(4, 20, 12, 28, temp12, temp13, temp8,  temp9,  temp10, temp11)
    VERTICAL_PASS(6, 22, 14, 30, temp20, temp21, temp12, temp13, temp14, temp15)
    "addiu   %[temp21],      $zero,          255             \n\t"
    HORIZONTAL_PASS( 0,  1,  2,  3, temp0, temp4, temp8,  temp12)
    HORIZONTAL_PASS(16, 17, 18, 19, temp1, temp5, temp9,  temp13)
    HORIZONTAL_PASS(32, 33, 34, 35, temp2, temp6, temp10, temp14)
    HORIZONTAL_PASS(48, 49, 50, 51, temp3, temp7, temp11, temp15)

    : [temp0]"=&r"(temp0), [temp1]"=&r"(temp1), [temp2]"=&r"(temp2),
      [temp3]"=&r"(temp3), [temp4]"=&r"(temp4), [temp5]"=&r"(temp5),
      [temp6]"=&r"(temp6), [temp7]"=&r"(temp7), [temp8]"=&r"(temp8),
      [temp9]"=&r"(temp9), [temp10]"=&r"(temp10), [temp11]"=&r"(temp11),
      [temp12]"=&r"(temp12), [temp13]"=&r"(temp13), [temp14]"=&r"(temp14),
      [temp15]"=&r"(temp15), [temp16]"=&r"(temp16), [temp17]"=&r"(temp17),
      [temp18]"=&r"(temp18), [temp19]"=&r"(temp19), [temp20]"=&r"(temp20),
      [temp21]"=&r"(temp21)
    : [args]"r"(args), [kC1]"r"(kC1), [kC2]"r"(kC2)
    : "memory", "hi", "lo"
  );
}

static void ITransformMIPS32(const uint8_t* ref, const int16_t* in,
                             uint8_t* dst, int do_two) {
  ITransformOneMIPS32(ref, in, dst);
  if (do_two) {
    ITransformOneMIPS32(ref + 4, in + 16, dst + 4);
  }
}

#undef VERTICAL_PASS
#undef HORIZONTAL_PASS

// macro for one pass through for loop in QuantizeBlock
// QUANTDIV macro inlined
// J - offset in bytes (kZigzag[n] * 2)
// K - offset in bytes (kZigzag[n] * 4)
// N - offset in bytes (n * 2)
#define QUANTIZE_ONE(J, K, N)                                               \
  "lh           %[temp0],       "#J"(%[ppin])                       \n\t"   \
  "lhu          %[temp1],       "#J"(%[ppsharpen])                  \n\t"   \
  "lw           %[temp2],       "#K"(%[ppzthresh])                  \n\t"   \
  "sra          %[sign],        %[temp0],           15              \n\t"   \
  "xor          %[coeff],       %[temp0],           %[sign]         \n\t"   \
  "subu         %[coeff],       %[coeff],           %[sign]         \n\t"   \
  "addu         %[coeff],       %[coeff],           %[temp1]        \n\t"   \
  "slt          %[temp4],       %[temp2],           %[coeff]        \n\t"   \
  "addiu        %[temp5],       $zero,              0               \n\t"   \
  "addiu        %[level],       $zero,              0               \n\t"   \
  "beqz         %[temp4],       2f                                  \n\t"   \
  "lhu          %[temp1],       "#J"(%[ppiq])                       \n\t"   \
  "lw           %[temp2],       "#K"(%[ppbias])                     \n\t"   \
  "lhu          %[temp3],       "#J"(%[ppq])                        \n\t"   \
  "mul          %[level],       %[coeff],           %[temp1]        \n\t"   \
  "addu         %[level],       %[level],           %[temp2]        \n\t"   \
  "sra          %[level],       %[level],           17              \n\t"   \
  "slt          %[temp4],       %[max_level],       %[level]        \n\t"   \
  "movn         %[level],       %[max_level],       %[temp4]        \n\t"   \
  "xor          %[level],       %[level],           %[sign]         \n\t"   \
  "subu         %[level],       %[level],           %[sign]         \n\t"   \
  "mul          %[temp5],       %[level],           %[temp3]        \n\t"   \
"2:                                                                 \n\t"   \
  "sh           %[temp5],       "#J"(%[ppin])                       \n\t"   \
  "sh           %[level],       "#N"(%[pout])                       \n\t"

static int QuantizeBlockMIPS32(int16_t in[16], int16_t out[16],
                               const VP8Matrix* const mtx) {
  int temp0, temp1, temp2, temp3, temp4, temp5;
  int sign, coeff, level, i;
  int max_level = MAX_LEVEL;

  int16_t* ppin             = &in[0];
  int16_t* pout             = &out[0];
  const uint16_t* ppsharpen = &mtx->sharpen_[0];
  const uint32_t* ppzthresh = &mtx->zthresh_[0];
  const uint16_t* ppq       = &mtx->q_[0];
  const uint16_t* ppiq      = &mtx->iq_[0];
  const uint32_t* ppbias    = &mtx->bias_[0];

  __asm__ volatile(
    QUANTIZE_ONE( 0,  0,  0)
    QUANTIZE_ONE( 2,  4,  2)
    QUANTIZE_ONE( 8, 16,  4)
    QUANTIZE_ONE(16, 32,  6)
    QUANTIZE_ONE(10, 20,  8)
    QUANTIZE_ONE( 4,  8, 10)
    QUANTIZE_ONE( 6, 12, 12)
    QUANTIZE_ONE(12, 24, 14)
    QUANTIZE_ONE(18, 36, 16)
    QUANTIZE_ONE(24, 48, 18)
    QUANTIZE_ONE(26, 52, 20)
    QUANTIZE_ONE(20, 40, 22)
    QUANTIZE_ONE(14, 28, 24)
    QUANTIZE_ONE(22, 44, 26)
    QUANTIZE_ONE(28, 56, 28)
    QUANTIZE_ONE(30, 60, 30)

    : [temp0]"=&r"(temp0), [temp1]"=&r"(temp1),
      [temp2]"=&r"(temp2), [temp3]"=&r"(temp3),
      [temp4]"=&r"(temp4), [temp5]"=&r"(temp5),
      [sign]"=&r"(sign), [coeff]"=&r"(coeff),
      [level]"=&r"(level)
    : [pout]"r"(pout), [ppin]"r"(ppin),
      [ppiq]"r"(ppiq), [max_level]"r"(max_level),
      [ppbias]"r"(ppbias), [ppzthresh]"r"(ppzthresh),
      [ppsharpen]"r"(ppsharpen), [ppq]"r"(ppq)
    : "memory", "hi", "lo"
  );

  // moved out from macro to increase possibility for earlier breaking
  for (i = 15; i >= 0; i--) {
    if (out[i]) return 1;
  }
  return 0;
}

#undef QUANTIZE_ONE

#endif  // WEBP_USE_MIPS32

//------------------------------------------------------------------------------
// Entry point

extern void VP8EncDspInitMIPS32(void);

void VP8EncDspInitMIPS32(void) {
#if defined(WEBP_USE_MIPS32)
  VP8ITransform = ITransformMIPS32;
  VP8EncQuantizeBlock = QuantizeBlockMIPS32;
#endif  // WEBP_USE_MIPS32
}
