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
// Author(s): Darko Laus (darko.laus@imgtec.com)
//            Mirko Raus (mirko.raus@imgtec.com)

#include "./dsp.h"

#if defined(WEBP_USE_MIPS_DSP_R2)

#include "./mips_macro.h"
#include "../enc/cost.h"
#include "../enc/vp8enci.h"

static const int kC1 = 20091 + (1 << 16);
static const int kC2 = 35468;

// O - output
// I - input (macro doesn't change it)
#define ADD_SUB_HALVES_X4(O0, O1, O2, O3, O4, O5, O6, O7,                      \
                          I0, I1, I2, I3, I4, I5, I6, I7)                      \
  "addq.ph          %["#O0"],   %["#I0"],  %["#I1"]           \n\t"            \
  "subq.ph          %["#O1"],   %["#I0"],  %["#I1"]           \n\t"            \
  "addq.ph          %["#O2"],   %["#I2"],  %["#I3"]           \n\t"            \
  "subq.ph          %["#O3"],   %["#I2"],  %["#I3"]           \n\t"            \
  "addq.ph          %["#O4"],   %["#I4"],  %["#I5"]           \n\t"            \
  "subq.ph          %["#O5"],   %["#I4"],  %["#I5"]           \n\t"            \
  "addq.ph          %["#O6"],   %["#I6"],  %["#I7"]           \n\t"            \
  "subq.ph          %["#O7"],   %["#I6"],  %["#I7"]           \n\t"

// IO - input/output
#define ABS_X8(IO0, IO1, IO2, IO3, IO4, IO5, IO6, IO7)                         \
  "absq_s.ph        %["#IO0"],   %["#IO0"]                    \n\t"            \
  "absq_s.ph        %["#IO1"],   %["#IO1"]                    \n\t"            \
  "absq_s.ph        %["#IO2"],   %["#IO2"]                    \n\t"            \
  "absq_s.ph        %["#IO3"],   %["#IO3"]                    \n\t"            \
  "absq_s.ph        %["#IO4"],   %["#IO4"]                    \n\t"            \
  "absq_s.ph        %["#IO5"],   %["#IO5"]                    \n\t"            \
  "absq_s.ph        %["#IO6"],   %["#IO6"]                    \n\t"            \
  "absq_s.ph        %["#IO7"],   %["#IO7"]                    \n\t"

// dpa.w.ph $ac0 temp0 ,temp1
//  $ac += temp0[31..16] * temp1[31..16] + temp0[15..0] * temp1[15..0]
// dpax.w.ph $ac0 temp0 ,temp1
//  $ac += temp0[31..16] * temp1[15..0] + temp0[15..0] * temp1[31..16]
// O - output
// I - input (macro doesn't change it)
#define MUL_HALF(O0, I0, I1, I2, I3, I4, I5, I6, I7,                           \
                 I8, I9, I10, I11, I12, I13, I14, I15)                         \
    "mult            $ac0,      $zero,     $zero              \n\t"            \
    "dpa.w.ph        $ac0,      %["#I2"],  %["#I0"]           \n\t"            \
    "dpax.w.ph       $ac0,      %["#I5"],  %["#I6"]           \n\t"            \
    "dpa.w.ph        $ac0,      %["#I8"],  %["#I9"]           \n\t"            \
    "dpax.w.ph       $ac0,      %["#I11"], %["#I4"]           \n\t"            \
    "dpa.w.ph        $ac0,      %["#I12"], %["#I7"]           \n\t"            \
    "dpax.w.ph       $ac0,      %["#I13"], %["#I1"]           \n\t"            \
    "dpa.w.ph        $ac0,      %["#I14"], %["#I3"]           \n\t"            \
    "dpax.w.ph       $ac0,      %["#I15"], %["#I10"]          \n\t"            \
    "mflo            %["#O0"],  $ac0                          \n\t"

#define OUTPUT_EARLY_CLOBBER_REGS_17()                                         \
  OUTPUT_EARLY_CLOBBER_REGS_10(),                                              \
  [temp11]"=&r"(temp11), [temp12]"=&r"(temp12), [temp13]"=&r"(temp13),         \
  [temp14]"=&r"(temp14), [temp15]"=&r"(temp15), [temp16]"=&r"(temp16),         \
  [temp17]"=&r"(temp17)

// macro for one horizontal pass in FTransform
// temp0..temp15 holds tmp[0]..tmp[15]
// A - offset in bytes to load from src and ref buffers
// TEMP0..TEMP3 - registers for corresponding tmp elements
#define HORIZONTAL_PASS(A, TEMP0, TEMP1, TEMP2, TEMP3)                         \
  "lw              %["#TEMP0"],   0(%[args])                        \n\t"      \
  "lw              %["#TEMP1"],   4(%[args])                        \n\t"      \
  "lw              %["#TEMP2"],   "XSTR(BPS)"*"#A"(%["#TEMP0"])     \n\t"      \
  "lw              %["#TEMP3"],   "XSTR(BPS)"*"#A"(%["#TEMP1"])     \n\t"      \
  "preceu.ph.qbl   %["#TEMP0"],   %["#TEMP2"]                       \n\t"      \
  "preceu.ph.qbl   %["#TEMP1"],   %["#TEMP3"]                       \n\t"      \
  "preceu.ph.qbr   %["#TEMP2"],   %["#TEMP2"]                       \n\t"      \
  "preceu.ph.qbr   %["#TEMP3"],   %["#TEMP3"]                       \n\t"      \
  "subq.ph         %["#TEMP0"],   %["#TEMP0"],   %["#TEMP1"]        \n\t"      \
  "subq.ph         %["#TEMP2"],   %["#TEMP2"],   %["#TEMP3"]        \n\t"      \
  "rotr            %["#TEMP0"],   %["#TEMP0"],   16                 \n\t"      \
  "addq.ph         %["#TEMP1"],   %["#TEMP2"],   %["#TEMP0"]        \n\t"      \
  "subq.ph         %["#TEMP3"],   %["#TEMP2"],   %["#TEMP0"]        \n\t"      \
  "seh             %["#TEMP0"],   %["#TEMP1"]                       \n\t"      \
  "sra             %[temp16],     %["#TEMP1"],   16                 \n\t"      \
  "seh             %[temp19],     %["#TEMP3"]                       \n\t"      \
  "sra             %["#TEMP3"],   %["#TEMP3"],   16                 \n\t"      \
  "subu            %["#TEMP2"],   %["#TEMP0"],   %[temp16]          \n\t"      \
  "addu            %["#TEMP0"],   %["#TEMP0"],   %[temp16]          \n\t"      \
  "mul             %[temp17],     %[temp19],     %[c2217]           \n\t"      \
  "mul             %[temp18],     %["#TEMP3"],   %[c5352]           \n\t"      \
  "mul             %["#TEMP1"],   %[temp19],     %[c5352]           \n\t"      \
  "mul             %[temp16],     %["#TEMP3"],   %[c2217]           \n\t"      \
  "sll             %["#TEMP2"],   %["#TEMP2"],   3                  \n\t"      \
  "sll             %["#TEMP0"],   %["#TEMP0"],   3                  \n\t"      \
  "subu            %["#TEMP3"],   %[temp17],     %[temp18]          \n\t"      \
  "addu            %["#TEMP1"],   %[temp16],     %["#TEMP1"]        \n\t"      \
  "addiu           %["#TEMP3"],   %["#TEMP3"],   937                \n\t"      \
  "addiu           %["#TEMP1"],   %["#TEMP1"],   1812               \n\t"      \
  "sra             %["#TEMP3"],   %["#TEMP3"],   9                  \n\t"      \
  "sra             %["#TEMP1"],   %["#TEMP1"],   9                  \n\t"

// macro for one vertical pass in FTransform
// temp0..temp15 holds tmp[0]..tmp[15]
// A..D - offsets in bytes to store to out buffer
// TEMP0, TEMP4, TEMP8 and TEMP12 - registers for corresponding tmp elements
#define VERTICAL_PASS(A, B, C, D, TEMP0, TEMP4, TEMP8, TEMP12)                 \
  "addu            %[temp16],     %["#TEMP0"],   %["#TEMP12"] \n\t"            \
  "subu            %[temp19],     %["#TEMP0"],   %["#TEMP12"] \n\t"            \
  "addu            %[temp17],     %["#TEMP4"],   %["#TEMP8"]  \n\t"            \
  "subu            %[temp18],     %["#TEMP4"],   %["#TEMP8"]  \n\t"            \
  "mul             %["#TEMP8"],   %[temp19],     %[c2217]     \n\t"            \
  "mul             %["#TEMP12"],  %[temp18],     %[c2217]     \n\t"            \
  "mul             %["#TEMP4"],   %[temp19],     %[c5352]     \n\t"            \
  "mul             %[temp18],     %[temp18],     %[c5352]     \n\t"            \
  "addiu           %[temp16],     %[temp16],     7            \n\t"            \
  "addu            %["#TEMP0"],   %[temp16],     %[temp17]    \n\t"            \
  "sra             %["#TEMP0"],   %["#TEMP0"],   4            \n\t"            \
  "addu            %["#TEMP12"],  %["#TEMP12"],  %["#TEMP4"]  \n\t"            \
  "subu            %["#TEMP4"],   %[temp16],     %[temp17]    \n\t"            \
  "sra             %["#TEMP4"],   %["#TEMP4"],   4            \n\t"            \
  "addiu           %["#TEMP8"],   %["#TEMP8"],   30000        \n\t"            \
  "addiu           %["#TEMP12"],  %["#TEMP12"],  12000        \n\t"            \
  "addiu           %["#TEMP8"],   %["#TEMP8"],   21000        \n\t"            \
  "subu            %["#TEMP8"],   %["#TEMP8"],   %[temp18]    \n\t"            \
  "sra             %["#TEMP12"],  %["#TEMP12"],  16           \n\t"            \
  "sra             %["#TEMP8"],   %["#TEMP8"],   16           \n\t"            \
  "addiu           %[temp16],     %["#TEMP12"],  1            \n\t"            \
  "movn            %["#TEMP12"],  %[temp16],     %[temp19]    \n\t"            \
  "sh              %["#TEMP0"],   "#A"(%[temp20])             \n\t"            \
  "sh              %["#TEMP4"],   "#C"(%[temp20])             \n\t"            \
  "sh              %["#TEMP8"],   "#D"(%[temp20])             \n\t"            \
  "sh              %["#TEMP12"],  "#B"(%[temp20])             \n\t"

static void FTransform(const uint8_t* src, const uint8_t* ref, int16_t* out) {
  const int c2217 = 2217;
  const int c5352 = 5352;
  int temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8;
  int temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16;
  int temp17, temp18, temp19, temp20;
  const int* const args[3] =
      { (const int*)src, (const int*)ref, (const int*)out };

  __asm__ volatile (
    HORIZONTAL_PASS(0, temp0,  temp1,  temp2,  temp3)
    HORIZONTAL_PASS(1, temp4,  temp5,  temp6,  temp7)
    HORIZONTAL_PASS(2, temp8,  temp9,  temp10, temp11)
    HORIZONTAL_PASS(3, temp12, temp13, temp14, temp15)
    "lw            %[temp20],     8(%[args])                  \n\t"
    VERTICAL_PASS(0,  8, 16, 24, temp0, temp4, temp8,  temp12)
    VERTICAL_PASS(2, 10, 18, 26, temp1, temp5, temp9,  temp13)
    VERTICAL_PASS(4, 12, 20, 28, temp2, temp6, temp10, temp14)
    VERTICAL_PASS(6, 14, 22, 30, temp3, temp7, temp11, temp15)
    OUTPUT_EARLY_CLOBBER_REGS_18(),
      [temp0]"=&r"(temp0), [temp19]"=&r"(temp19), [temp20]"=&r"(temp20)
    : [args]"r"(args), [c2217]"r"(c2217), [c5352]"r"(c5352)
    : "memory", "hi", "lo"
  );
}

#undef VERTICAL_PASS
#undef HORIZONTAL_PASS

static WEBP_INLINE void ITransformOne(const uint8_t* ref, const int16_t* in,
                                      uint8_t* dst) {
  int temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9;
  int temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18;

  __asm__ volatile (
    "ulw              %[temp1],   0(%[in])                 \n\t"
    "ulw              %[temp2],   16(%[in])                \n\t"
    LOAD_IN_X2(temp5, temp6, 24, 26)
    ADD_SUB_HALVES(temp3, temp4, temp1, temp2)
    LOAD_IN_X2(temp1, temp2, 8, 10)
    MUL_SHIFT_SUM(temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14,
                  temp10, temp8, temp9, temp7, temp1, temp2, temp5, temp6,
                  temp13, temp11, temp14, temp12)
    INSERT_HALF_X2(temp8, temp7, temp10, temp9)
    "ulw              %[temp17],  4(%[in])                 \n\t"
    "ulw              %[temp18],  20(%[in])                \n\t"
    ADD_SUB_HALVES(temp1, temp2, temp3, temp8)
    ADD_SUB_HALVES(temp5, temp6, temp4, temp7)
    ADD_SUB_HALVES(temp7, temp8, temp17, temp18)
    LOAD_IN_X2(temp17, temp18, 12, 14)
    LOAD_IN_X2(temp9, temp10, 28, 30)
    MUL_SHIFT_SUM(temp11, temp12, temp13, temp14, temp15, temp16, temp4, temp17,
                  temp12, temp14, temp11, temp13, temp17, temp18, temp9, temp10,
                  temp15, temp4, temp16, temp17)
    INSERT_HALF_X2(temp11, temp12, temp13, temp14)
    ADD_SUB_HALVES(temp17, temp8, temp8, temp11)
    ADD_SUB_HALVES(temp3, temp4, temp7, temp12)

    // horizontal
    SRA_16(temp9, temp10, temp11, temp12, temp1, temp2, temp5, temp6)
    INSERT_HALF_X2(temp1, temp6, temp5, temp2)
    SRA_16(temp13, temp14, temp15, temp16, temp3, temp4, temp17, temp8)
    "repl.ph          %[temp2],   0x4                      \n\t"
    INSERT_HALF_X2(temp3, temp8, temp17, temp4)
    "addq.ph          %[temp1],   %[temp1],  %[temp2]      \n\t"
    "addq.ph          %[temp6],   %[temp6],  %[temp2]      \n\t"
    ADD_SUB_HALVES(temp2, temp4, temp1, temp3)
    ADD_SUB_HALVES(temp5, temp7, temp6, temp8)
    MUL_SHIFT_SUM(temp1, temp3, temp6, temp8, temp9, temp13, temp17, temp18,
                  temp3, temp13, temp1, temp9, temp9, temp13, temp11, temp15,
                  temp6, temp17, temp8, temp18)
    MUL_SHIFT_SUM(temp6, temp8, temp18, temp17, temp11, temp15, temp12, temp16,
                  temp8, temp15, temp6, temp11, temp12, temp16, temp10, temp14,
                  temp18, temp12, temp17, temp16)
    INSERT_HALF_X2(temp1, temp3, temp9, temp13)
    INSERT_HALF_X2(temp6, temp8, temp11, temp15)
    SHIFT_R_SUM_X2(temp9, temp10, temp11, temp12, temp13, temp14, temp15,
                   temp16, temp2, temp4, temp5, temp7, temp3, temp1, temp8,
                   temp6)
    PACK_2_HALVES_TO_WORD(temp1, temp2, temp3, temp4, temp9, temp12, temp13,
                          temp16, temp11, temp10, temp15, temp14)
    LOAD_WITH_OFFSET_X4(temp10, temp11, temp14, temp15, ref,
                        0, 0, 0, 0,
                        0, 1, 2, 3,
                        BPS)
    CONVERT_2_BYTES_TO_HALF(temp5, temp6, temp7, temp8, temp17, temp18, temp10,
                            temp11, temp10, temp11, temp14, temp15)
    STORE_SAT_SUM_X2(temp5, temp6, temp7, temp8, temp17, temp18, temp10, temp11,
                     temp9, temp12, temp1, temp2, temp13, temp16, temp3, temp4,
                     dst, 0, 1, 2, 3, BPS)

    OUTPUT_EARLY_CLOBBER_REGS_18()
    : [dst]"r"(dst), [in]"r"(in), [kC1]"r"(kC1), [kC2]"r"(kC2), [ref]"r"(ref)
    : "memory", "hi", "lo"
  );
}

static void ITransform(const uint8_t* ref, const int16_t* in, uint8_t* dst,
                       int do_two) {
  ITransformOne(ref, in, dst);
  if (do_two) {
    ITransformOne(ref + 4, in + 16, dst + 4);
  }
}

static int Disto4x4(const uint8_t* const a, const uint8_t* const b,
                    const uint16_t* const w) {
  int temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9;
  int temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17;

  __asm__ volatile (
    LOAD_WITH_OFFSET_X4(temp1, temp2, temp3, temp4, a,
                        0, 0, 0, 0,
                        0, 1, 2, 3,
                        BPS)
    CONVERT_2_BYTES_TO_HALF(temp5, temp6, temp7, temp8, temp9,temp10, temp11,
                            temp12, temp1, temp2, temp3, temp4)
    ADD_SUB_HALVES_X4(temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8,
                      temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12)
    PACK_2_HALVES_TO_WORD(temp9, temp10, temp11, temp12, temp1, temp3, temp5,
                          temp7, temp2, temp4, temp6, temp8)
    ADD_SUB_HALVES_X4(temp2, temp4, temp6, temp8, temp9, temp1, temp3, temp10,
                      temp1, temp9, temp3, temp10, temp5, temp11, temp7, temp12)
    ADD_SUB_HALVES_X4(temp5, temp11, temp7, temp2, temp9, temp3, temp6, temp12,
                      temp2, temp9, temp6, temp3, temp4, temp1, temp8, temp10)
    ADD_SUB_HALVES_X4(temp1, temp4, temp10, temp8, temp7, temp11, temp5, temp2,
                      temp5, temp7, temp11, temp2, temp9, temp6, temp3, temp12)
    ABS_X8(temp1, temp4, temp10, temp8, temp7, temp11, temp5, temp2)
    LOAD_WITH_OFFSET_X4(temp3, temp6, temp9, temp12, w,
                        0, 4, 8, 12,
                        0, 0, 0, 0,
                        0)
    LOAD_WITH_OFFSET_X4(temp13, temp14, temp15, temp16, w,
                        0, 4, 8, 12,
                        1, 1, 1, 1,
                        16)
    MUL_HALF(temp17, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8,
             temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16)
    LOAD_WITH_OFFSET_X4(temp1, temp2, temp3, temp4, b,
                        0, 0, 0, 0,
                        0, 1, 2, 3,
                        BPS)
    CONVERT_2_BYTES_TO_HALF(temp5,temp6, temp7, temp8, temp9,temp10, temp11,
                            temp12, temp1, temp2, temp3, temp4)
    ADD_SUB_HALVES_X4(temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8,
                      temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12)
    PACK_2_HALVES_TO_WORD(temp9, temp10, temp11, temp12, temp1, temp3, temp5,
                          temp7, temp2, temp4, temp6, temp8)
    ADD_SUB_HALVES_X4(temp2, temp4, temp6, temp8, temp9, temp1, temp3, temp10,
                      temp1, temp9, temp3, temp10, temp5, temp11, temp7, temp12)
    ADD_SUB_HALVES_X4(temp5, temp11, temp7, temp2, temp9, temp3, temp6, temp12,
                      temp2, temp9, temp6, temp3, temp4, temp1, temp8, temp10)
    ADD_SUB_HALVES_X4(temp1, temp4, temp10, temp8, temp7, temp11, temp5, temp2,
                      temp5, temp7, temp11, temp2, temp9, temp6, temp3, temp12)
    ABS_X8(temp1, temp4, temp10, temp8, temp7, temp11, temp5, temp2)
    LOAD_WITH_OFFSET_X4(temp3, temp6, temp9, temp12, w,
                        0, 4, 8, 12,
                        0, 0, 0, 0,
                        0)
    LOAD_WITH_OFFSET_X4(temp13, temp14, temp15, temp16, w,
                        0, 4, 8, 12,
                        1, 1, 1, 1,
                        16)
    MUL_HALF(temp3, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8,
             temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16)
    OUTPUT_EARLY_CLOBBER_REGS_17()
    : [a]"r"(a), [b]"r"(b), [w]"r"(w)
    : "memory", "hi", "lo"
  );
  return abs(temp3 - temp17) >> 5;
}

static int Disto16x16(const uint8_t* const a, const uint8_t* const b,
                      const uint16_t* const w) {
  int D = 0;
  int x, y;
  for (y = 0; y < 16 * BPS; y += 4 * BPS) {
    for (x = 0; x < 16; x += 4) {
      D += Disto4x4(a + x + y, b + x + y, w);
    }
  }
  return D;
}

//------------------------------------------------------------------------------
// Intra predictions

static WEBP_INLINE uint8_t clip_8b(int v) {
  return (!(v & ~0xff)) ? v : (v < 0) ? 0 : 255;
}

static uint8_t clip1[255 + 510 + 1];    // clips [-255,510] to [0,255]

static volatile int tables_ok = 0;

static void InitTables(void) {
  if (!tables_ok) {
    int i;
    for (i = -255; i <= 255 + 255; ++i) {
      clip1[255 + i] = clip_8b(i);
    }
    tables_ok = 1;
  }
}

#define FILL_PART(J, SIZE)                                          \
    "usw        %[value],  0+"#J"*"XSTR(BPS)"(%[dst])    \n\t"      \
    "usw        %[value],  4+"#J"*"XSTR(BPS)"(%[dst])    \n\t"      \
  ".if "#SIZE" == 16                                     \n\t"      \
    "usw        %[value],  8+"#J"*"XSTR(BPS)"(%[dst])    \n\t"      \
    "usw        %[value], 12+"#J"*"XSTR(BPS)"(%[dst])    \n\t"      \
  ".endif                                                \n\t"

#define FILL_8_OR_16(DST, VALUE, SIZE) do {                         \
  int value = (VALUE);                                              \
  __asm__ volatile (                                                \
    "replv.qb   %[value],  %[value]                      \n\t"      \
    FILL_PART( 0, SIZE)                                             \
    FILL_PART( 1, SIZE)                                             \
    FILL_PART( 2, SIZE)                                             \
    FILL_PART( 3, SIZE)                                             \
    FILL_PART( 4, SIZE)                                             \
    FILL_PART( 5, SIZE)                                             \
    FILL_PART( 6, SIZE)                                             \
    FILL_PART( 7, SIZE)                                             \
  ".if "#SIZE" == 16                                     \n\t"      \
    FILL_PART( 8, 16)                                               \
    FILL_PART( 9, 16)                                               \
    FILL_PART(10, 16)                                               \
    FILL_PART(11, 16)                                               \
    FILL_PART(12, 16)                                               \
    FILL_PART(13, 16)                                               \
    FILL_PART(14, 16)                                               \
    FILL_PART(15, 16)                                               \
  ".endif                                                \n\t"      \
    : [value]"+&r"(value)                                           \
    : [dst]"r"((DST))                                               \
    : "memory"                                                      \
  );                                                                \
} while (0)

#define VERTICAL_PRED(DST, TOP, SIZE)                                          \
static WEBP_INLINE void VerticalPred##SIZE(uint8_t* (DST),                     \
                                           const uint8_t* (TOP)) {             \
  int j;                                                                       \
  if ((TOP)) {                                                                 \
    for (j = 0; j < (SIZE); ++j) memcpy((DST) + j * BPS, (TOP), (SIZE));       \
  } else {                                                                     \
    FILL_8_OR_16((DST), 127, (SIZE));                                          \
  }                                                                            \
}

VERTICAL_PRED(dst, top, 8)
VERTICAL_PRED(dst, top, 16)

#undef VERTICAL_PRED

#define HORIZONTAL_PRED(DST, LEFT, SIZE)                                       \
static WEBP_INLINE void HorizontalPred##SIZE(uint8_t* (DST),                   \
                                             const uint8_t* (LEFT)) {          \
  if (LEFT) {                                                                  \
    int j;                                                                     \
    for (j = 0; j < (SIZE); ++j) {                                             \
      memset((DST) + j * BPS, (LEFT)[j], (SIZE));                              \
    }                                                                          \
  } else {                                                                     \
    FILL_8_OR_16((DST), 129, (SIZE));                                          \
  }                                                                            \
}

HORIZONTAL_PRED(dst, left, 8)
HORIZONTAL_PRED(dst, left, 16)

#undef HORIZONTAL_PRED

#define CLIPPING()                                                             \
  "preceu.ph.qbl   %[temp2],   %[temp0]                  \n\t"                 \
  "preceu.ph.qbr   %[temp0],   %[temp0]                  \n\t"                 \
  "preceu.ph.qbl   %[temp3],   %[temp1]                  \n\t"                 \
  "preceu.ph.qbr   %[temp1],   %[temp1]                  \n\t"                 \
  "addu.ph         %[temp2],   %[temp2],   %[leftY_1]    \n\t"                 \
  "addu.ph         %[temp0],   %[temp0],   %[leftY_1]    \n\t"                 \
  "addu.ph         %[temp3],   %[temp3],   %[leftY_1]    \n\t"                 \
  "addu.ph         %[temp1],   %[temp1],   %[leftY_1]    \n\t"                 \
  "shll_s.ph       %[temp2],   %[temp2],   7             \n\t"                 \
  "shll_s.ph       %[temp0],   %[temp0],   7             \n\t"                 \
  "shll_s.ph       %[temp3],   %[temp3],   7             \n\t"                 \
  "shll_s.ph       %[temp1],   %[temp1],   7             \n\t"                 \
  "precrqu_s.qb.ph %[temp0],   %[temp2],   %[temp0]      \n\t"                 \
  "precrqu_s.qb.ph %[temp1],   %[temp3],   %[temp1]      \n\t"

#define CLIP_8B_TO_DST(DST, LEFT, TOP, SIZE) do {                              \
  int leftY_1 = ((int)(LEFT)[y] << 16) + (LEFT)[y];                            \
  int temp0, temp1, temp2, temp3;                                              \
  __asm__ volatile (                                                           \
    "replv.ph        %[leftY_1], %[leftY_1]              \n\t"                 \
    "ulw             %[temp0],   0(%[top])               \n\t"                 \
    "ulw             %[temp1],   4(%[top])               \n\t"                 \
    "subu.ph         %[leftY_1], %[leftY_1], %[left_1]   \n\t"                 \
    CLIPPING()                                                                 \
    "usw             %[temp0],   0(%[dst])               \n\t"                 \
    "usw             %[temp1],   4(%[dst])               \n\t"                 \
  ".if "#SIZE" == 16                                     \n\t"                 \
    "ulw             %[temp0],   8(%[top])               \n\t"                 \
    "ulw             %[temp1],   12(%[top])              \n\t"                 \
    CLIPPING()                                                                 \
    "usw             %[temp0],   8(%[dst])               \n\t"                 \
    "usw             %[temp1],   12(%[dst])              \n\t"                 \
  ".endif                                                \n\t"                 \
    : [leftY_1]"+&r"(leftY_1), [temp0]"=&r"(temp0), [temp1]"=&r"(temp1),       \
      [temp2]"=&r"(temp2), [temp3]"=&r"(temp3)                                 \
    : [left_1]"r"(left_1), [top]"r"((TOP)), [dst]"r"((DST))                    \
    : "memory"                                                                 \
  );                                                                           \
} while (0)

#define CLIP_TO_DST(DST, LEFT, TOP, SIZE) do {                                 \
  int y;                                                                       \
  const int left_1 = ((int)(LEFT)[-1] << 16) + (LEFT)[-1];                     \
  for (y = 0; y < (SIZE); ++y) {                                               \
    CLIP_8B_TO_DST((DST), (LEFT), (TOP), (SIZE));                              \
    (DST) += BPS;                                                              \
  }                                                                            \
} while (0)

#define TRUE_MOTION(DST, LEFT, TOP, SIZE)                                      \
static WEBP_INLINE void TrueMotion##SIZE(uint8_t* (DST), const uint8_t* (LEFT),\
                                         const uint8_t* (TOP)) {               \
  if (LEFT) {                                                                  \
    if (TOP) {                                                                 \
      CLIP_TO_DST((DST), (LEFT), (TOP), (SIZE));                               \
    } else {                                                                   \
      HorizontalPred##SIZE((DST), (LEFT));                                     \
    }                                                                          \
  } else {                                                                     \
    /* true motion without left samples (hence: with default 129 value)    */  \
    /* is equivalent to VE prediction where you just copy the top samples. */  \
    /* Note that if top samples are not available, the default value is    */  \
    /* then 129, and not 127 as in the VerticalPred case.                  */  \
    if (TOP) {                                                                 \
      VerticalPred##SIZE((DST), (TOP));                                        \
    } else {                                                                   \
      FILL_8_OR_16((DST), 129, (SIZE));                                        \
    }                                                                          \
  }                                                                            \
}

TRUE_MOTION(dst, left, top, 8)
TRUE_MOTION(dst, left, top, 16)

#undef TRUE_MOTION
#undef CLIP_TO_DST
#undef CLIP_8B_TO_DST
#undef CLIPPING

static WEBP_INLINE void DCMode16(uint8_t* dst, const uint8_t* left,
                                 const uint8_t* top) {
  int DC, DC1;
  int temp0, temp1, temp2, temp3;

  __asm__ volatile(
    "beqz        %[top],   2f                  \n\t"
    LOAD_WITH_OFFSET_X4(temp0, temp1, temp2, temp3, top,
                        0, 4, 8, 12,
                        0, 0, 0, 0,
                        0)
    "raddu.w.qb  %[temp0], %[temp0]            \n\t"
    "raddu.w.qb  %[temp1], %[temp1]            \n\t"
    "raddu.w.qb  %[temp2], %[temp2]            \n\t"
    "raddu.w.qb  %[temp3], %[temp3]            \n\t"
    "addu        %[temp0], %[temp0], %[temp1]  \n\t"
    "addu        %[temp2], %[temp2], %[temp3]  \n\t"
    "addu        %[DC],    %[temp0], %[temp2]  \n\t"
    "move        %[DC1],   %[DC]               \n\t"
    "beqz        %[left],  1f                  \n\t"
    LOAD_WITH_OFFSET_X4(temp0, temp1, temp2, temp3, left,
                        0, 4, 8, 12,
                        0, 0, 0, 0,
                        0)
    "raddu.w.qb  %[temp0], %[temp0]            \n\t"
    "raddu.w.qb  %[temp1], %[temp1]            \n\t"
    "raddu.w.qb  %[temp2], %[temp2]            \n\t"
    "raddu.w.qb  %[temp3], %[temp3]            \n\t"
    "addu        %[temp0], %[temp0], %[temp1]  \n\t"
    "addu        %[temp2], %[temp2], %[temp3]  \n\t"
    "addu        %[DC1],   %[temp0], %[temp2]  \n\t"
  "1:                                          \n\t"
    "addu        %[DC],   %[DC],     %[DC1]    \n\t"
    "j           3f                            \n\t"
  "2:                                          \n\t"
    "beqz        %[left],  4f                  \n\t"
    LOAD_WITH_OFFSET_X4(temp0, temp1, temp2, temp3, left,
                        0, 4, 8, 12,
                        0, 0, 0, 0,
                        0)
    "raddu.w.qb  %[temp0], %[temp0]            \n\t"
    "raddu.w.qb  %[temp1], %[temp1]            \n\t"
    "raddu.w.qb  %[temp2], %[temp2]            \n\t"
    "raddu.w.qb  %[temp3], %[temp3]            \n\t"
    "addu        %[temp0], %[temp0], %[temp1]  \n\t"
    "addu        %[temp2], %[temp2], %[temp3]  \n\t"
    "addu        %[DC],    %[temp0], %[temp2]  \n\t"
    "addu        %[DC],    %[DC],    %[DC]     \n\t"
  "3:                                          \n\t"
    "shra_r.w    %[DC],    %[DC],    5         \n\t"
    "j           5f                            \n\t"
  "4:                                          \n\t"
    "li          %[DC],    0x80                \n\t"
  "5:                                          \n\t"
    : [temp0]"=&r"(temp0), [temp1]"=&r"(temp1), [DC]"=&r"(DC),
      [temp2]"=&r"(temp2), [temp3]"=&r"(temp3), [DC1]"=&r"(DC1)
    : [left]"r"(left), [top]"r"(top)
    : "memory"
  );

  FILL_8_OR_16(dst, DC, 16);
}

static WEBP_INLINE void DCMode8(uint8_t* dst, const uint8_t* left,
                                const uint8_t* top) {
  int DC, DC1;
  int temp0, temp1, temp2, temp3;

  __asm__ volatile(
    "beqz        %[top],   2f                  \n\t"
    "ulw         %[temp0], 0(%[top])           \n\t"
    "ulw         %[temp1], 4(%[top])           \n\t"
    "raddu.w.qb  %[temp0], %[temp0]            \n\t"
    "raddu.w.qb  %[temp1], %[temp1]            \n\t"
    "addu        %[DC],    %[temp0], %[temp1]  \n\t"
    "move        %[DC1],   %[DC]               \n\t"
    "beqz        %[left],  1f                  \n\t"
    "ulw         %[temp2], 0(%[left])          \n\t"
    "ulw         %[temp3], 4(%[left])          \n\t"
    "raddu.w.qb  %[temp2], %[temp2]            \n\t"
    "raddu.w.qb  %[temp3], %[temp3]            \n\t"
    "addu        %[DC1],   %[temp2], %[temp3]  \n\t"
  "1:                                          \n\t"
    "addu        %[DC],    %[DC],    %[DC1]    \n\t"
    "j           3f                            \n\t"
  "2:                                          \n\t"
    "beqz        %[left],  4f                  \n\t"
    "ulw         %[temp2], 0(%[left])          \n\t"
    "ulw         %[temp3], 4(%[left])          \n\t"
    "raddu.w.qb  %[temp2], %[temp2]            \n\t"
    "raddu.w.qb  %[temp3], %[temp3]            \n\t"
    "addu        %[DC],    %[temp2], %[temp3]  \n\t"
    "addu        %[DC],    %[DC],    %[DC]     \n\t"
  "3:                                          \n\t"
    "shra_r.w    %[DC], %[DC], 4               \n\t"
    "j           5f                            \n\t"
  "4:                                          \n\t"
    "li          %[DC], 0x80                   \n\t"
  "5:                                          \n\t"
    : [temp0]"=&r"(temp0), [temp1]"=&r"(temp1), [DC]"=&r"(DC),
      [temp2]"=&r"(temp2), [temp3]"=&r"(temp3), [DC1]"=&r"(DC1)
    : [left]"r"(left), [top]"r"(top)
    : "memory"
  );

  FILL_8_OR_16(dst, DC, 8);
}

//------------------------------------------------------------------------------
// Chroma 8x8 prediction (paragraph 12.2)

static void IntraChromaPreds(uint8_t* dst, const uint8_t* left,
                             const uint8_t* top) {
  // U block
  DCMode8(C8DC8 + dst, left, top);
  VerticalPred8(C8VE8 + dst, top);
  HorizontalPred8(C8HE8 + dst, left);
  TrueMotion8(C8TM8 + dst, left, top);
  // V block
  dst += 8;
  if (top) top += 8;
  if (left) left += 16;
  DCMode8(C8DC8 + dst, left, top);
  VerticalPred8(C8VE8 + dst, top);
  HorizontalPred8(C8HE8 + dst, left);
  TrueMotion8(C8TM8 + dst, left, top);
}

//------------------------------------------------------------------------------
// luma 16x16 prediction (paragraph 12.3)

static void Intra16Preds(uint8_t* dst,
                         const uint8_t* left, const uint8_t* top) {
  DCMode16(I16DC16 + dst, left, top);
  VerticalPred16(I16VE16 + dst, top);
  HorizontalPred16(I16HE16 + dst, left);
  TrueMotion16(I16TM16 + dst, left, top);
}

#undef FILL_8_OR_16
#undef FILL_PART
#undef OUTPUT_EARLY_CLOBBER_REGS_17
#undef MUL_HALF
#undef ABS_X8
#undef ADD_SUB_HALVES_X4

#endif  // WEBP_USE_MIPS_DSP_R2

//------------------------------------------------------------------------------
// Entry point

extern WEBP_TSAN_IGNORE_FUNCTION void VP8EncDspInitMIPSdspR2(void);

WEBP_TSAN_IGNORE_FUNCTION void VP8EncDspInitMIPSdspR2(void) {
#if defined(WEBP_USE_MIPS_DSP_R2)
  InitTables();
  VP8FTransform = FTransform;
  VP8ITransform = ITransform;
  VP8TDisto4x4 = Disto4x4;
  VP8TDisto16x16 = Disto16x16;
  VP8EncPredLuma16 = Intra16Preds;
  VP8EncPredChroma8 = IntraChromaPreds;
#endif  // WEBP_USE_MIPS_DSP_R2
}
