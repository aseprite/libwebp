// Copyright 2014 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
// MIPS version of dsp functions
//
// Author(s):  Djordje Pesut    (djordje.pesut@imgtec.com)
//             Jovan Zelincevic (jovan.zelincevic@imgtec.com)

#include "./dsp.h"

#if defined(WEBP_USE_MIPS_DSP_R2)

static const int kC1 = 20091 + (1 << 16);
static const int kC2 = 35468;

#define MUL(a, b) (((a) * (b)) >> 16)

// temp0[31..16 | 15..0] = temp0[31..16 | 15..0] + temp8[31..16 | 15..0]
// temp0[31..16 | 15..0] = temp0[31..16 <<(s) 7 | 15..0 <<(s) 7]
// temp1..temp7 same as temp0
// precrqu_s.qb.ph temp0, temp1, temp0:
//   temp0 = temp1[31..24] | temp1[15..8] | temp0[31..24] | temp0[15..8]
// store temp0 to dst
// IO - input/output
// I - input (macro doesn't change it)
#define STORE_SAT_SUM_X2(IO0, IO1, IO2, IO3, IO4, IO5, IO6, IO7,               \
                         I0, I1, I2, I3, I4, I5, I6, I7)                       \
  "addq.ph          %["#IO0"],  %["#IO0"],  %["#I0"]          \n\t"            \
  "addq.ph          %["#IO1"],  %["#IO1"],  %["#I1"]          \n\t"            \
  "addq.ph          %["#IO2"],  %["#IO2"],  %["#I2"]          \n\t"            \
  "addq.ph          %["#IO3"],  %["#IO3"],  %["#I3"]          \n\t"            \
  "addq.ph          %["#IO4"],  %["#IO4"],  %["#I4"]          \n\t"            \
  "addq.ph          %["#IO5"],  %["#IO5"],  %["#I5"]          \n\t"            \
  "addq.ph          %["#IO6"],  %["#IO6"],  %["#I6"]          \n\t"            \
  "addq.ph          %["#IO7"],  %["#IO7"],  %["#I7"]          \n\t"            \
  "shll_s.ph        %["#IO0"],  %["#IO0"],  7                 \n\t"            \
  "shll_s.ph        %["#IO1"],  %["#IO1"],  7                 \n\t"            \
  "shll_s.ph        %["#IO2"],  %["#IO2"],  7                 \n\t"            \
  "shll_s.ph        %["#IO3"],  %["#IO3"],  7                 \n\t"            \
  "shll_s.ph        %["#IO4"],  %["#IO4"],  7                 \n\t"            \
  "shll_s.ph        %["#IO5"],  %["#IO5"],  7                 \n\t"            \
  "shll_s.ph        %["#IO6"],  %["#IO6"],  7                 \n\t"            \
  "shll_s.ph        %["#IO7"],  %["#IO7"],  7                 \n\t"            \
  "precrqu_s.qb.ph  %["#IO0"],  %["#IO1"],  %["#IO0"]         \n\t"            \
  "precrqu_s.qb.ph  %["#IO2"],  %["#IO3"],  %["#IO2"]         \n\t"            \
  "precrqu_s.qb.ph  %["#IO4"],  %["#IO5"],  %["#IO4"]         \n\t"            \
  "precrqu_s.qb.ph  %["#IO6"],  %["#IO7"],  %["#IO6"]         \n\t"            \
  "usw              %["#IO0"],  0(%[dst])                     \n\t"            \
  "usw              %["#IO2"],  32(%[dst])                    \n\t"            \
  "usw              %["#IO4"],  64(%[dst])                    \n\t"            \
  "usw              %["#IO6"],  96(%[dst])                    \n\t"

// temp0[31..16 | 15..0] = temp8[31..16 | 15..0] + temp12[31..16 | 15..0]
// temp1[31..16 | 15..0] = temp8[31..16 | 15..0] - temp12[31..16 | 15..0]
// temp0[31..16 | 15..0] = temp0[31..16 >> 3 | 15..0 >> 3]
// temp1[31..16 | 15..0] = temp1[31..16 >> 3 | 15..0 >> 3]
// O - output
// I - input (macro doesn't change it)
#define SHIFT_R_SUM_X2(O0, O1, O2, O3, O4, O5, O6, O7,                         \
                       I0, I1, I2, I3, I4, I5, I6, I7)                         \
  "addq.ph          %["#O0"],   %["#I0"],   %["#I4"]          \n\t"            \
  "subq.ph          %["#O1"],   %["#I0"],   %["#I4"]          \n\t"            \
  "addq.ph          %["#O2"],   %["#I1"],   %["#I5"]          \n\t"            \
  "subq.ph          %["#O3"],   %["#I1"],   %["#I5"]          \n\t"            \
  "addq.ph          %["#O4"],   %["#I2"],   %["#I6"]          \n\t"            \
  "subq.ph          %["#O5"],   %["#I2"],   %["#I6"]          \n\t"            \
  "addq.ph          %["#O6"],   %["#I3"],   %["#I7"]          \n\t"            \
  "subq.ph          %["#O7"],   %["#I3"],   %["#I7"]          \n\t"            \
  "shra.ph          %["#O0"],   %["#O0"],   3                 \n\t"            \
  "shra.ph          %["#O1"],   %["#O1"],   3                 \n\t"            \
  "shra.ph          %["#O2"],   %["#O2"],   3                 \n\t"            \
  "shra.ph          %["#O3"],   %["#O3"],   3                 \n\t"            \
  "shra.ph          %["#O4"],   %["#O4"],   3                 \n\t"            \
  "shra.ph          %["#O5"],   %["#O5"],   3                 \n\t"            \
  "shra.ph          %["#O6"],   %["#O6"],   3                 \n\t"            \
  "shra.ph          %["#O7"],   %["#O7"],   3                 \n\t"

// preceu.ph.qbr temp0, temp8
//   temp0 = 0 | 0 | temp8[23..16] | temp8[7..0]
// preceu.ph.qbl temp1, temp8
//   temp1 = temp8[23..16] | temp8[7..0] | 0 | 0
// O - output
// I - input (macro doesn't change it)
#define CONVERT_2_BYTES_TO_HALF(O0, O1, O2, O3, O4, O5, O6, O7,                \
                                I0, I1, I2, I3)                                \
  "preceu.ph.qbr    %["#O0"],   %["#I0"]                      \n\t"            \
  "preceu.ph.qbl    %["#O1"],   %["#I0"]                      \n\t"            \
  "preceu.ph.qbr    %["#O2"],   %["#I1"]                      \n\t"            \
  "preceu.ph.qbl    %["#O3"],   %["#I1"]                      \n\t"            \
  "preceu.ph.qbr    %["#O4"],   %["#I2"]                      \n\t"            \
  "preceu.ph.qbl    %["#O5"],   %["#I2"]                      \n\t"            \
  "preceu.ph.qbr    %["#O6"],   %["#I3"]                      \n\t"            \
  "preceu.ph.qbl    %["#O7"],   %["#I3"]                      \n\t"

// O - output
#define LOAD_DST(O0, O1, O2, O3)                                               \
  "ulw              %["#O0"],  0(%[dst])                      \n\t"            \
  "ulw              %["#O1"],  32(%[dst])                     \n\t"            \
  "ulw              %["#O2"],  64(%[dst])                     \n\t"            \
  "ulw              %["#O3"],  96(%[dst])                     \n\t"

// precrq.ph.w temp0, temp8, temp2
//   temp0 = temp8[31..16] | temp2[31..16]
// ins temp2, temp8, 16, 16
//   temp2 = temp8[31..16] | temp2[15..0]
// O - output
// IO - input/output
// I - input (macro doesn't change it)
#define PACK_2_HALVES_TO_WORD(O0, O1, O2, O3,                                  \
                              IO0, IO1, IO2, IO3,                              \
                              I0, I1, I2, I3)                                  \
  "precrq.ph.w      %["#O0"],    %["#I0"],  %["#IO0"]         \n\t"            \
  "precrq.ph.w      %["#O1"],    %["#I1"],  %["#IO1"]         \n\t"            \
  "ins              %["#IO0"],   %["#I0"],  16,    16         \n\t"            \
  "ins              %["#IO1"],   %["#I1"],  16,    16         \n\t"            \
  "precrq.ph.w      %["#O2"],    %["#I2"],  %["#IO2"]         \n\t"            \
  "precrq.ph.w      %["#O3"],    %["#I3"],  %["#IO3"]         \n\t"            \
  "ins              %["#IO2"],   %["#I2"],  16,    16         \n\t"            \
  "ins              %["#IO3"],   %["#I3"],  16,    16         \n\t"

// O - output
// IO - input/output
// I - input (macro doesn't change it)
#define MUL_SHIFT_SUM(O0, O1, O2, O3, O4, O5, O6, O7,                          \
                      IO0, IO1, IO2, IO3,                                      \
                      I0, I1, I2, I3, I4, I5, I6, I7)                          \
  "mul              %["#O0"],   %["#I0"],   %[kC2]            \n\t"            \
  "mul              %["#O1"],   %["#I0"],   %[kC1]            \n\t"            \
  "mul              %["#O2"],   %["#I1"],   %[kC2]            \n\t"            \
  "mul              %["#O3"],   %["#I1"],   %[kC1]            \n\t"            \
  "mul              %["#O4"],   %["#I2"],   %[kC2]            \n\t"            \
  "mul              %["#O5"],   %["#I2"],   %[kC1]            \n\t"            \
  "mul              %["#O6"],   %["#I3"],   %[kC2]            \n\t"            \
  "mul              %["#O7"],   %["#I3"],   %[kC1]            \n\t"            \
  "sra              %["#O0"],   %["#O0"],   16                \n\t"            \
  "sra              %["#O1"],   %["#O1"],   16                \n\t"            \
  "sra              %["#O2"],   %["#O2"],   16                \n\t"            \
  "sra              %["#O3"],   %["#O3"],   16                \n\t"            \
  "sra              %["#O4"],   %["#O4"],   16                \n\t"            \
  "sra              %["#O5"],   %["#O5"],   16                \n\t"            \
  "sra              %["#O6"],   %["#O6"],   16                \n\t"            \
  "sra              %["#O7"],   %["#O7"],   16                \n\t"            \
  "addu             %["#IO0"],  %["#IO0"],  %["#I4"]          \n\t"            \
  "addu             %["#IO1"],  %["#IO1"],  %["#I5"]          \n\t"            \
  "subu             %["#IO2"],  %["#IO2"],  %["#I6"]          \n\t"            \
  "subu             %["#IO3"],  %["#IO3"],  %["#I7"]          \n\t"

// O - output
// I - input (macro doesn't change it)
#define ADD_SUB_HALVES(O0, O1,                                                 \
                       I0, I1)                                                 \
  "addq.ph          %["#O0"],   %["#I0"],  %["#I1"]           \n\t"            \
  "subq.ph          %["#O1"],   %["#I0"],  %["#I1"]           \n\t"

// O - output
// I - input (macro doesn't change it)
// I[0/1] - offset in bytes
#define LOAD_IN_X2(O0, O1,                                                     \
                   I0, I1)                                                     \
  "lh               %["#O0"],   "#I0"(%[in])                  \n\t"            \
  "lh               %["#O1"],   "#I1"(%[in])                  \n\t"

// O - output
// I - input (macro doesn't change it)
#define SRA_16(O0, O1, O2, O3,                                                 \
               I0, I1, I2, I3)                                                 \
  "sra              %["#O0"],  %["#I0"],  16                  \n\t"            \
  "sra              %["#O1"],  %["#I1"],  16                  \n\t"            \
  "sra              %["#O2"],  %["#I2"],  16                  \n\t"            \
  "sra              %["#O3"],  %["#I3"],  16                  \n\t"

// O - output
// I - input (macro doesn't change it)
#define INSERT_HALF_X2(O0, O1,                                                 \
                       I0, I1)                                                 \
  "ins              %["#O0"],   %["#I0"], 16,    16           \n\t"            \
  "ins              %["#O1"],   %["#I1"], 16,    16           \n\t"

#define OUTPUT_EARLY_CLOBBER_REGS_10()                                         \
  : [temp1]"=&r"(temp1), [temp2]"=&r"(temp2), [temp3]"=&r"(temp3),             \
    [temp4]"=&r"(temp4), [temp5]"=&r"(temp5), [temp6]"=&r"(temp6),             \
    [temp7]"=&r"(temp7), [temp8]"=&r"(temp8), [temp9]"=&r"(temp9),             \
    [temp10]"=&r"(temp10)

#define OUTPUT_EARLY_CLOBBER_REGS_18()                                         \
  OUTPUT_EARLY_CLOBBER_REGS_10(),                                              \
  [temp11]"=&r"(temp11), [temp12]"=&r"(temp12), [temp13]"=&r"(temp13),         \
  [temp14]"=&r"(temp14), [temp15]"=&r"(temp15), [temp16]"=&r"(temp16),         \
  [temp17]"=&r"(temp17), [temp18]"=&r"(temp18)

static void TransformDC(const int16_t* in, uint8_t* dst) {
  int temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10;

  __asm__ volatile (
    LOAD_DST(temp1, temp2, temp3, temp4)
    "lh               %[temp5],  0(%[in])               \n\t"
    "addiu            %[temp5],  %[temp5],  4           \n\t"
    "ins              %[temp5],  %[temp5],  16, 16      \n\t"
    "shra.ph          %[temp5],  %[temp5],  3           \n\t"
    CONVERT_2_BYTES_TO_HALF(temp6, temp7, temp8, temp9, temp10, temp1, temp2,
                            temp3, temp1, temp2, temp3, temp4)
    STORE_SAT_SUM_X2(temp6, temp7, temp8, temp9, temp10, temp1, temp2, temp3,
                     temp5, temp5, temp5, temp5, temp5, temp5, temp5, temp5)

    OUTPUT_EARLY_CLOBBER_REGS_10()
    : [in]"r"(in), [dst]"r"(dst)
    : "memory"
  );
}

static void TransformAC3(const int16_t* in, uint8_t* dst) {
  const int a = in[0] + 4;
  int c4 = MUL(in[4], kC2);
  const int d4 = MUL(in[4], kC1);
  const int c1 = MUL(in[1], kC2);
  const int d1 = MUL(in[1], kC1);
  int temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9;
  int temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18;

  __asm__ volatile (
    "ins              %[c4],      %[d4],     16,       16    \n\t"
    "replv.ph         %[temp1],   %[a]                       \n\t"
    "replv.ph         %[temp4],   %[d1]                      \n\t"
    ADD_SUB_HALVES(temp2, temp3, temp1, c4)
    "replv.ph         %[temp5],   %[c1]                      \n\t"
    SHIFT_R_SUM_X2(temp1, temp6, temp7, temp8, temp2, temp9, temp10, temp4,
                   temp2, temp2, temp3, temp3, temp4, temp5, temp4, temp5)
    LOAD_DST(temp3, temp5, temp11, temp12)
    CONVERT_2_BYTES_TO_HALF(temp13, temp14, temp3, temp15, temp5, temp16,
                            temp11, temp17, temp3, temp5, temp11, temp12)
    PACK_2_HALVES_TO_WORD(temp12, temp18, temp7, temp6, temp1, temp8, temp2,
                          temp4, temp7, temp6, temp10, temp9)
    STORE_SAT_SUM_X2(temp13, temp14, temp3, temp15, temp5, temp16, temp11,
                     temp17, temp12, temp18, temp1, temp8, temp2, temp4,
                     temp7, temp6)

    OUTPUT_EARLY_CLOBBER_REGS_18(),
      [c4]"+&r"(c4)
    : [dst]"r"(dst), [a]"r"(a), [d1]"r"(d1), [d4]"r"(d4), [c1]"r"(c1)
    : "memory"
  );
}

static void TransformOne(const int16_t* in, uint8_t* dst) {
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
    LOAD_DST(temp10, temp11, temp14, temp15)
    CONVERT_2_BYTES_TO_HALF(temp5, temp6, temp7, temp8, temp17, temp18, temp10,
                            temp11, temp10, temp11, temp14, temp15)
    STORE_SAT_SUM_X2(temp5, temp6, temp7, temp8, temp17, temp18, temp10, temp11,
                     temp9, temp12, temp1, temp2, temp13, temp16, temp3, temp4)

    OUTPUT_EARLY_CLOBBER_REGS_18()
    : [dst]"r"(dst), [in]"r"(in), [kC1]"r"(kC1), [kC2]"r"(kC2)
    : "memory", "hi", "lo"
  );
}

static void TransformTwo(const int16_t* in, uint8_t* dst, int do_two) {
  TransformOne(in, dst);
  if (do_two) {
    TransformOne(in + 16, dst + 4);
  }
}

static WEBP_INLINE void FilterLoop26(uint8_t* p,
                                     int hstride, int vstride, int size,
                                     int thresh, int ithresh, int hev_thresh) {
  const int thresh2 = 2 * thresh + 1;
  int temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9;
  int temp10, temp11, temp12, temp13, temp14, temp15;

  __asm__ volatile (
    ".set      push                                      \n\t"
    ".set      noreorder                                 \n\t"
  "1:                                                    \n\t"
    "negu      %[temp1],  %[hstride]                     \n\t"
    "addiu     %[size],   %[size],        -1             \n\t"
    "sll       %[temp2],  %[hstride],     1              \n\t"
    "sll       %[temp3],  %[temp1],       1              \n\t"
    "addu      %[temp4],  %[temp2],       %[hstride]     \n\t"
    "addu      %[temp5],  %[temp3],       %[temp1]       \n\t"
    "lbu       %[temp7],  0(%[p])                        \n\t"
    "sll       %[temp6],  %[temp3],       1              \n\t"
    "lbux      %[temp8],  %[temp5](%[p])                 \n\t"
    "lbux      %[temp9],  %[temp3](%[p])                 \n\t"
    "lbux      %[temp10], %[temp1](%[p])                 \n\t"
    "lbux      %[temp11], %[temp6](%[p])                 \n\t"
    "lbux      %[temp12], %[hstride](%[p])               \n\t"
    "lbux      %[temp13], %[temp2](%[p])                 \n\t"
    "lbux      %[temp14], %[temp4](%[p])                 \n\t"
    "subu      %[temp1],  %[temp10],      %[temp7]       \n\t"
    "subu      %[temp2],  %[temp9],       %[temp12]      \n\t"
    "absq_s.w  %[temp3],  %[temp1]                       \n\t"
    "absq_s.w  %[temp4],  %[temp2]                       \n\t"
    "negu      %[temp1],  %[temp1]                       \n\t"
    "sll       %[temp3],  %[temp3],       2              \n\t"
    "addu      %[temp15], %[temp3],       %[temp4]       \n\t"
    "subu      %[temp3],  %[temp15],      %[thresh2]     \n\t"
    "sll       %[temp6],  %[temp1],       1              \n\t"
    "bgtz      %[temp3],  3f                             \n\t"
    " subu     %[temp4],  %[temp11],      %[temp8]       \n\t"
    "absq_s.w  %[temp4],  %[temp4]                       \n\t"
    "shll_s.w  %[temp2],  %[temp2],       23             \n\t"
    "subu      %[temp4],  %[temp4],       %[ithresh]     \n\t"
    "bgtz      %[temp4],  3f                             \n\t"
    " subu     %[temp3],  %[temp8],       %[temp9]       \n\t"
    "absq_s.w  %[temp3],  %[temp3]                       \n\t"
    "subu      %[temp3],  %[temp3],       %[ithresh]     \n\t"
    "bgtz      %[temp3],  3f                             \n\t"
    " subu     %[temp5],  %[temp9],       %[temp10]      \n\t"
    "absq_s.w  %[temp3],  %[temp5]                       \n\t"
    "absq_s.w  %[temp5],  %[temp5]                       \n\t"
    "subu      %[temp3],  %[temp3],       %[ithresh]     \n\t"
    "bgtz      %[temp3],  3f                             \n\t"
    " subu     %[temp3],  %[temp14],      %[temp13]      \n\t"
    "absq_s.w  %[temp3],  %[temp3]                       \n\t"
    "slt       %[temp5],  %[hev_thresh],  %[temp5]       \n\t"
    "subu      %[temp3],  %[temp3],       %[ithresh]     \n\t"
    "bgtz      %[temp3],  3f                             \n\t"
    " subu     %[temp3],  %[temp13],      %[temp12]      \n\t"
    "absq_s.w  %[temp3],  %[temp3]                       \n\t"
    "sra       %[temp4],  %[temp2],       23             \n\t"
    "subu      %[temp3],  %[temp3],       %[ithresh]     \n\t"
    "bgtz      %[temp3],  3f                             \n\t"
    " subu     %[temp15], %[temp12],      %[temp7]       \n\t"
    "absq_s.w  %[temp3],  %[temp15]                      \n\t"
    "absq_s.w  %[temp15], %[temp15]                      \n\t"
    "subu      %[temp3],  %[temp3],       %[ithresh]     \n\t"
    "bgtz      %[temp3],  3f                             \n\t"
    " slt      %[temp15], %[hev_thresh],  %[temp15]      \n\t"
    "addu      %[temp3],  %[temp6],       %[temp1]       \n\t"
    "or        %[temp2],  %[temp5],       %[temp15]      \n\t"
    "addu      %[temp5],  %[temp4],       %[temp3]       \n\t"
    "beqz      %[temp2],  4f                             \n\t"
    " shra_r.w %[temp1],  %[temp5],       3              \n\t"
    "addiu     %[temp2],  %[temp5],       3              \n\t"
    "sra       %[temp2],  %[temp2],       3              \n\t"
    "shll_s.w  %[temp1],  %[temp1],       27             \n\t"
    "shll_s.w  %[temp2],  %[temp2],       27             \n\t"
    "subu      %[temp3],  %[p],           %[hstride]     \n\t"
    "sra       %[temp1],  %[temp1],       27             \n\t"
    "sra       %[temp2],  %[temp2],       27             \n\t"
    "subu      %[temp1],  %[temp7],       %[temp1]       \n\t"
    "addu      %[temp2],  %[temp10],      %[temp2]       \n\t"
    "lbux      %[temp2],  %[temp2](%[VP8kclip1])         \n\t"
    "lbux      %[temp1],  %[temp1](%[VP8kclip1])         \n\t"
    "sb        %[temp2],  0(%[temp3])                    \n\t"
    "j         3f                                        \n\t"
    " sb       %[temp1],  0(%[p])                        \n\t"
  "4:                                                    \n\t"
    "shll_s.w  %[temp5],  %[temp5],       23             \n\t"
    "subu      %[temp14], %[p],           %[hstride]     \n\t"
    "subu      %[temp11], %[temp14],      %[hstride]     \n\t"
    "sra       %[temp6],  %[temp5],       23             \n\t"
    "sll       %[temp1],  %[temp6],       3              \n\t"
    "subu      %[temp15], %[temp11],      %[hstride]     \n\t"
    "addu      %[temp2],  %[temp6],       %[temp1]       \n\t"
    "sll       %[temp3],  %[temp2],       1              \n\t"
    "addu      %[temp4],  %[temp3],       %[temp2]       \n\t"
    "addiu     %[temp2],  %[temp2],       63             \n\t"
    "addiu     %[temp3],  %[temp3],       63             \n\t"
    "addiu     %[temp4],  %[temp4],       63             \n\t"
    "sra       %[temp2],  %[temp2],       7              \n\t"
    "sra       %[temp3],  %[temp3],       7              \n\t"
    "sra       %[temp4],  %[temp4],       7              \n\t"
    "addu      %[temp1],  %[temp8],       %[temp2]       \n\t"
    "addu      %[temp5],  %[temp9],       %[temp3]       \n\t"
    "addu      %[temp6],  %[temp10],      %[temp4]       \n\t"
    "subu      %[temp8],  %[temp7],       %[temp4]       \n\t"
    "subu      %[temp7],  %[temp12],      %[temp3]       \n\t"
    "addu      %[temp10], %[p],           %[hstride]     \n\t"
    "subu      %[temp9],  %[temp13],      %[temp2]       \n\t"
    "addu      %[temp12], %[temp10],      %[hstride]     \n\t"
    "lbux      %[temp2],  %[temp1](%[VP8kclip1])         \n\t"
    "lbux      %[temp3],  %[temp5](%[VP8kclip1])         \n\t"
    "lbux      %[temp4],  %[temp6](%[VP8kclip1])         \n\t"
    "lbux      %[temp5],  %[temp8](%[VP8kclip1])         \n\t"
    "lbux      %[temp6],  %[temp7](%[VP8kclip1])         \n\t"
    "lbux      %[temp8],  %[temp9](%[VP8kclip1])         \n\t"
    "sb        %[temp2],  0(%[temp15])                   \n\t"
    "sb        %[temp3],  0(%[temp11])                   \n\t"
    "sb        %[temp4],  0(%[temp14])                   \n\t"
    "sb        %[temp5],  0(%[p])                        \n\t"
    "sb        %[temp6],  0(%[temp10])                   \n\t"
    "sb        %[temp8],  0(%[temp12])                   \n\t"
  "3:                                                    \n\t"
    "bgtz      %[size],   1b                             \n\t"
    " addu     %[p],      %[p],           %[vstride]     \n\t"
    ".set      pop                                       \n\t"
    : [temp1]"=&r"(temp1), [temp2]"=&r"(temp2),[temp3]"=&r"(temp3),
      [temp4]"=&r"(temp4), [temp5]"=&r"(temp5), [temp6]"=&r"(temp6),
      [temp7]"=&r"(temp7),[temp8]"=&r"(temp8),[temp9]"=&r"(temp9),
      [temp10]"=&r"(temp10),[temp11]"=&r"(temp11),[temp12]"=&r"(temp12),
      [temp13]"=&r"(temp13),[temp14]"=&r"(temp14),[temp15]"=&r"(temp15),
      [size]"+&r"(size), [p]"+&r"(p)
    : [hstride]"r"(hstride), [thresh2]"r"(thresh2),
      [ithresh]"r"(ithresh),[vstride]"r"(vstride), [hev_thresh]"r"(hev_thresh),
      [VP8kclip1]"r"(VP8kclip1)
    : "memory"
  );
}

static WEBP_INLINE void FilterLoop24(uint8_t* p,
                                     int hstride, int vstride, int size,
                                     int thresh, int ithresh, int hev_thresh) {
  int p0, q0, p1, q1, p2, q2, p3, q3;
  int step1, step2, temp1, temp2, temp3, temp4;
  uint8_t* p_pom;
  uint8_t* p_pom2;
  const int thresh2 = 2 * thresh + 1;

  __asm__ volatile (
    ".set      push                                   \n\t"
    ".set      noreorder                              \n\t"
    "bltz      %[size],    3f                         \n\t"
    " nop                                             \n\t"
  "2:                                                 \n\t"
    "negu      %[step1],   %[hstride]                 \n\t"
    "lbu       %[q0],      0(%[p])                    \n\t"
    "lbux      %[p0],      %[step1](%[p])             \n\t"
    "subu      %[step1],   %[step1],      %[hstride]  \n\t"
    "lbux      %[q1],      %[hstride](%[p])           \n\t"
    "subu      %[temp1],   %[p0],         %[q0]       \n\t"
    "lbux      %[p1],      %[step1](%[p])             \n\t"
    "addu      %[step2],   %[hstride],    %[hstride]  \n\t"
    "absq_s.w  %[temp2],   %[temp1]                   \n\t"
    "subu      %[temp3],   %[p1],         %[q1]       \n\t"
    "absq_s.w  %[temp4],   %[temp3]                   \n\t"
    "sll       %[temp2],   %[temp2],      2           \n\t"
    "addu      %[temp2],   %[temp2],      %[temp4]    \n\t"
    "subu      %[temp4],   %[temp2],      %[thresh2]  \n\t"
    "subu      %[step1],   %[step1],      %[hstride]  \n\t"
    "bgtz      %[temp4],   0f                         \n\t"
    " lbux     %[p2],      %[step1](%[p])             \n\t"
    "subu      %[step1],   %[step1],      %[hstride]  \n\t"
    "lbux      %[q2],      %[step2](%[p])             \n\t"
    "lbux      %[p3],      %[step1](%[p])             \n\t"
    "subu      %[temp4],   %[p2],         %[p1]       \n\t"
    "addu      %[step2],   %[step2],      %[hstride]  \n\t"
    "subu      %[temp2],   %[p3],         %[p2]       \n\t"
    "absq_s.w  %[temp4],   %[temp4]                   \n\t"
    "absq_s.w  %[temp2],   %[temp2]                   \n\t"
    "lbux      %[q3],      %[step2](%[p])             \n\t"
    "subu      %[temp4],   %[temp4],      %[ithresh]  \n\t"
    "negu      %[temp1],   %[temp1]                   \n\t"
    "bgtz      %[temp4],   0f                         \n\t"
    " subu     %[temp2],   %[temp2],      %[ithresh]  \n\t"
    "subu      %[p3],      %[p1],         %[p0]       \n\t"
    "bgtz      %[temp2],   0f                         \n\t"
    " absq_s.w %[p3],      %[p3]                      \n\t"
    "subu      %[temp4],   %[q3],         %[q2]       \n\t"
    "subu      %[p_pom],   %[p],          %[hstride]  \n\t"
    "absq_s.w  %[temp4],   %[temp4]                   \n\t"
    "subu      %[temp2],   %[p3],         %[ithresh]  \n\t"
    "sll       %[step1],   %[temp1],      1           \n\t"
    "bgtz      %[temp2],   0f                         \n\t"
    " subu     %[temp4],   %[temp4],      %[ithresh]  \n\t"
    "subu      %[temp2],   %[q2],         %[q1]       \n\t"
    "bgtz      %[temp4],   0f                         \n\t"
    " absq_s.w %[temp2],   %[temp2]                   \n\t"
    "subu      %[q3],      %[q1],         %[q0]       \n\t"
    "absq_s.w  %[q3],      %[q3]                      \n\t"
    "subu      %[temp2],   %[temp2],      %[ithresh]  \n\t"
    "addu      %[temp1],   %[temp1],      %[step1]    \n\t"
    "bgtz      %[temp2],   0f                         \n\t"
    " subu     %[temp4],   %[q3],         %[ithresh]  \n\t"
    "slt       %[p3],      %[hev_thresh], %[p3]       \n\t"
    "bgtz      %[temp4],   0f                         \n\t"
    " slt      %[q3],      %[hev_thresh], %[q3]       \n\t"
    "or        %[q3],      %[q3],         %[p3]       \n\t"
    "bgtz      %[q3],      1f                         \n\t"
    " shra_r.w %[temp2],   %[temp1],      3           \n\t"
    "addiu     %[temp1],   %[temp1],      3           \n\t"
    "sra       %[temp1],   %[temp1],      3           \n\t"
    "shll_s.w  %[temp2],   %[temp2],      27          \n\t"
    "shll_s.w  %[temp1],   %[temp1],      27          \n\t"
    "addu      %[p_pom2],  %[p],          %[hstride]  \n\t"
    "sra       %[temp2],   %[temp2],      27          \n\t"
    "sra       %[temp1],   %[temp1],      27          \n\t"
    "addiu     %[step1],   %[temp2],      1           \n\t"
    "sra       %[step1],   %[step1],      1           \n\t"
    "addu      %[p0],      %[p0],         %[temp1]    \n\t"
    "addu      %[p1],      %[p1],         %[step1]    \n\t"
    "subu      %[q0],      %[q0],         %[temp2]    \n\t"
    "subu      %[q1],      %[q1],         %[step1]    \n\t"
    "lbux      %[temp2],   %[p0](%[VP8kclip1])        \n\t"
    "lbux      %[temp3],   %[q0](%[VP8kclip1])        \n\t"
    "lbux      %[temp4],   %[q1](%[VP8kclip1])        \n\t"
    "sb        %[temp2],   0(%[p_pom])                \n\t"
    "lbux      %[temp1],   %[p1](%[VP8kclip1])        \n\t"
    "subu      %[p_pom],   %[p_pom], %[hstride]       \n\t"
    "sb        %[temp3],   0(%[p])                    \n\t"
    "sb        %[temp4],   0(%[p_pom2])               \n\t"
    "j         0f                                     \n\t"
    " sb       %[temp1],   0(%[p_pom])                \n\t"
  "1:                                                 \n\t"
    "shll_s.w  %[temp3],   %[temp3],      24          \n\t"
    "sra       %[temp3],   %[temp3],      24          \n\t"
    "addu      %[temp1],   %[temp1],      %[temp3]    \n\t"
    "shra_r.w  %[temp2],   %[temp1],      3           \n\t"
    "addiu     %[temp1],   %[temp1],      3           \n\t"
    "shll_s.w  %[temp2],   %[temp2],      27          \n\t"
    "sra       %[temp1],   %[temp1],      3           \n\t"
    "shll_s.w  %[temp1],   %[temp1],      27          \n\t"
    "sra       %[temp2],   %[temp2],      27          \n\t"
    "sra       %[temp1],   %[temp1],      27          \n\t"
    "addu      %[p0],      %[p0],         %[temp1]    \n\t"
    "subu      %[q0],      %[q0],         %[temp2]    \n\t"
    "lbux      %[temp1],   %[p0](%[VP8kclip1])        \n\t"
    "lbux      %[temp2],   %[q0](%[VP8kclip1])        \n\t"
    "sb        %[temp2],   0(%[p])                    \n\t"
    "sb        %[temp1],   0(%[p_pom])                \n\t"
  "0:                                                 \n\t"
    "subu      %[size],    %[size],       1           \n\t"
    "bgtz      %[size],    2b                         \n\t"
    " addu     %[p],       %[p],          %[vstride]  \n\t"
  "3:                                                 \n\t"
    ".set      pop                                    \n\t"
    : [p0]"=&r"(p0), [q0]"=&r"(q0), [p1]"=&r"(p1), [q1]"=&r"(q1),
      [p2]"=&r"(p2), [q2]"=&r"(q2), [p3]"=&r"(p3), [q3]"=&r"(q3),
      [step2]"=&r"(step2), [step1]"=&r"(step1), [temp1]"=&r"(temp1),
      [temp2]"=&r"(temp2), [temp3]"=&r"(temp3), [temp4]"=&r"(temp4),
      [p_pom]"=&r"(p_pom), [p_pom2]"=&r"(p_pom2), [p]"+&r"(p),
      [size]"+&r"(size)
    : [vstride]"r"(vstride), [ithresh]"r"(ithresh),
      [hev_thresh]"r"(hev_thresh), [hstride]"r"(hstride),
      [VP8kclip1]"r"(VP8kclip1), [thresh2]"r"(thresh2)
    : "memory"
  );
}

// on macroblock edges
static void VFilter16(uint8_t* p, int stride,
                      int thresh, int ithresh, int hev_thresh) {
  FilterLoop26(p, stride, 1, 16, thresh, ithresh, hev_thresh);
}

static void HFilter16(uint8_t* p, int stride,
                      int thresh, int ithresh, int hev_thresh) {
  FilterLoop26(p, 1, stride, 16, thresh, ithresh, hev_thresh);
}

// 8-pixels wide variant, for chroma filtering
static void VFilter8(uint8_t* u, uint8_t* v, int stride,
                     int thresh, int ithresh, int hev_thresh) {
  FilterLoop26(u, stride, 1, 8, thresh, ithresh, hev_thresh);
  FilterLoop26(v, stride, 1, 8, thresh, ithresh, hev_thresh);
}

static void HFilter8(uint8_t* u, uint8_t* v, int stride,
                     int thresh, int ithresh, int hev_thresh) {
  FilterLoop26(u, 1, stride, 8, thresh, ithresh, hev_thresh);
  FilterLoop26(v, 1, stride, 8, thresh, ithresh, hev_thresh);
}

// on three inner edges
static void VFilter16i(uint8_t* p, int stride,
                       int thresh, int ithresh, int hev_thresh) {
  int k;
  for (k = 3; k > 0; --k) {
    p += 4 * stride;
    FilterLoop24(p, stride, 1, 16, thresh, ithresh, hev_thresh);
  }
}

static void HFilter16i(uint8_t* p, int stride,
                       int thresh, int ithresh, int hev_thresh) {
  int k;
  for (k = 3; k > 0; --k) {
    p += 4;
    FilterLoop24(p, 1, stride, 16, thresh, ithresh, hev_thresh);
  }
}

static void VFilter8i(uint8_t* u, uint8_t* v, int stride,
                      int thresh, int ithresh, int hev_thresh) {
  FilterLoop24(u + 4 * stride, stride, 1, 8, thresh, ithresh, hev_thresh);
  FilterLoop24(v + 4 * stride, stride, 1, 8, thresh, ithresh, hev_thresh);
}

static void HFilter8i(uint8_t* u, uint8_t* v, int stride,
                      int thresh, int ithresh, int hev_thresh) {
  FilterLoop24(u + 4, 1, stride, 8, thresh, ithresh, hev_thresh);
  FilterLoop24(v + 4, 1, stride, 8, thresh, ithresh, hev_thresh);
}

#undef OUTPUT_EARLY_CLOBBER_REGS_18
#undef OUTPUT_EARLY_CLOBBER_REGS_10
#undef INSERT_HALF_X2
#undef SRA_16
#undef LOAD_IN_X2
#undef ADD_SUB_HALVES
#undef MUL_SHIFT_SUM
#undef PACK_2_HALVES_TO_WORD
#undef LOAD_DST
#undef CONVERT_2_BYTES_TO_HALF
#undef SHIFT_R_SUM_X2
#undef STORE_SAT_SUM_X2
#undef MUL

#endif  // WEBP_USE_MIPS_DSP_R2

//------------------------------------------------------------------------------
// Entry point

extern void VP8DspInitMIPSdspR2(void);

void VP8DspInitMIPSdspR2(void) {
#if defined(WEBP_USE_MIPS_DSP_R2)
  VP8TransformDC = TransformDC;
  VP8TransformAC3 = TransformAC3;
  VP8Transform = TransformTwo;
  VP8VFilter16 = VFilter16;
  VP8HFilter16 = HFilter16;
  VP8VFilter8 = VFilter8;
  VP8HFilter8 = HFilter8;
  VP8VFilter16i = VFilter16i;
  VP8HFilter16i = HFilter16i;
  VP8VFilter8i = VFilter8i;
  VP8HFilter8i = HFilter8i;
#endif
}
