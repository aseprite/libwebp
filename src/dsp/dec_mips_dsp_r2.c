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


#define STORE_SAT_SUM_X2(TEMP0, TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6,      \
                         TEMP7, TEMP8, TEMP9, TEMP10, TEMP11, TEMP12, TEMP13,  \
                         TEMP14, TEMP15)                                       \
  "addq.ph          %["#TEMP0"],  %["#TEMP0"],  %["#TEMP8"]     \n\t"          \
  "addq.ph          %["#TEMP1"],  %["#TEMP1"],  %["#TEMP9"]     \n\t"          \
  "addq.ph          %["#TEMP2"],  %["#TEMP2"],  %["#TEMP10"]    \n\t"          \
  "addq.ph          %["#TEMP3"],  %["#TEMP3"],  %["#TEMP11"]    \n\t"          \
  "addq.ph          %["#TEMP4"],  %["#TEMP4"],  %["#TEMP12"]    \n\t"          \
  "addq.ph          %["#TEMP5"],  %["#TEMP5"],  %["#TEMP13"]    \n\t"          \
  "addq.ph          %["#TEMP6"],  %["#TEMP6"],  %["#TEMP14"]    \n\t"          \
  "addq.ph          %["#TEMP7"],  %["#TEMP7"],  %["#TEMP15"]    \n\t"          \
  "shll_s.ph        %["#TEMP0"],  %["#TEMP0"],  7               \n\t"          \
  "shll_s.ph        %["#TEMP1"],  %["#TEMP1"],  7               \n\t"          \
  "shll_s.ph        %["#TEMP2"],  %["#TEMP2"],  7               \n\t"          \
  "shll_s.ph        %["#TEMP3"],  %["#TEMP3"],  7               \n\t"          \
  "shll_s.ph        %["#TEMP4"],  %["#TEMP4"],  7               \n\t"          \
  "shll_s.ph        %["#TEMP5"],  %["#TEMP5"],  7               \n\t"          \
  "shll_s.ph        %["#TEMP6"],  %["#TEMP6"],  7               \n\t"          \
  "shll_s.ph        %["#TEMP7"],  %["#TEMP7"],  7               \n\t"          \
  "precrqu_s.qb.ph  %["#TEMP0"],  %["#TEMP1"],  %["#TEMP0"]     \n\t"          \
  "precrqu_s.qb.ph  %["#TEMP2"],  %["#TEMP3"],  %["#TEMP2"]     \n\t"          \
  "precrqu_s.qb.ph  %["#TEMP4"],  %["#TEMP5"],  %["#TEMP4"]     \n\t"          \
  "precrqu_s.qb.ph  %["#TEMP6"],  %["#TEMP7"],  %["#TEMP6"]     \n\t"          \
  "usw              %["#TEMP0"],  0(%[dst])                     \n\t"          \
  "usw              %["#TEMP2"],  32(%[dst])                    \n\t"          \
  "usw              %["#TEMP4"],  64(%[dst])                    \n\t"          \
  "usw              %["#TEMP6"],  96(%[dst])                    \n\t"

#define SHIFT_R_SUM_X2(TEMP0, TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6, TEMP7, \
                       TEMP8, TEMP9, TEMP10, TEMP11, TEMP12, TEMP13, TEMP14,   \
                       TEMP15)                                                 \
  "addq.ph          %["#TEMP0"],   %["#TEMP8"],   %["#TEMP12"]      \n\t"      \
  "subq.ph          %["#TEMP1"],   %["#TEMP8"],   %["#TEMP12"]      \n\t"      \
  "addq.ph          %["#TEMP2"],   %["#TEMP9"],   %["#TEMP13"]      \n\t"      \
  "subq.ph          %["#TEMP3"],   %["#TEMP9"],   %["#TEMP13"]      \n\t"      \
  "addq.ph          %["#TEMP4"],   %["#TEMP10"],  %["#TEMP14"]      \n\t"      \
  "subq.ph          %["#TEMP5"],   %["#TEMP10"],  %["#TEMP14"]      \n\t"      \
  "addq.ph          %["#TEMP6"],   %["#TEMP11"],  %["#TEMP15"]      \n\t"      \
  "subq.ph          %["#TEMP7"],   %["#TEMP11"],  %["#TEMP15"]      \n\t"      \
  "shra.ph          %["#TEMP0"],   %["#TEMP0"],   3                 \n\t"      \
  "shra.ph          %["#TEMP1"],   %["#TEMP1"],   3                 \n\t"      \
  "shra.ph          %["#TEMP2"],   %["#TEMP2"],   3                 \n\t"      \
  "shra.ph          %["#TEMP3"],   %["#TEMP3"],   3                 \n\t"      \
  "shra.ph          %["#TEMP4"],   %["#TEMP4"],   3                 \n\t"      \
  "shra.ph          %["#TEMP5"],   %["#TEMP5"],   3                 \n\t"      \
  "shra.ph          %["#TEMP6"],   %["#TEMP6"],   3                 \n\t"      \
  "shra.ph          %["#TEMP7"],   %["#TEMP7"],   3                 \n\t"

#define PRECEU_RL(TEMP0, TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6,             \
                  TEMP7, TEMP8, TEMP9, TEMP10, TEMP11)                         \
  "preceu.ph.qbr    %["#TEMP0"],   %["#TEMP8"]                 \n\t"           \
  "preceu.ph.qbl    %["#TEMP1"],   %["#TEMP8"]                 \n\t"           \
  "preceu.ph.qbr    %["#TEMP2"],   %["#TEMP9"]                 \n\t"           \
  "preceu.ph.qbl    %["#TEMP3"],   %["#TEMP9"]                 \n\t"           \
  "preceu.ph.qbr    %["#TEMP4"],   %["#TEMP10"]                \n\t"           \
  "preceu.ph.qbl    %["#TEMP5"],   %["#TEMP10"]                \n\t"           \
  "preceu.ph.qbr    %["#TEMP6"],   %["#TEMP11"]                \n\t"           \
  "preceu.ph.qbl    %["#TEMP7"],   %["#TEMP11"]                \n\t"

#define LOAD_DST(TEMP0, TEMP1, TEMP2, TEMP3)                                   \
  "ulw              %["#TEMP0"],  0(%[dst])                \n\t"               \
  "ulw              %["#TEMP1"],  32(%[dst])               \n\t"               \
  "ulw              %["#TEMP2"],  64(%[dst])               \n\t"               \
  "ulw              %["#TEMP3"],  96(%[dst])               \n\t"

#define PRECRQ_X2(TEMP0, TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6,             \
                  TEMP7, TEMP8, TEMP9, TEMP10, TEMP11)                         \
  "precrq.ph.w      %["#TEMP0"],   %["#TEMP8"],  %["#TEMP2"]      \n\t"        \
  "precrq.ph.w      %["#TEMP1"],   %["#TEMP9"],  %["#TEMP3"]      \n\t"        \
  "ins              %["#TEMP2"],   %["#TEMP8"],  16,    16        \n\t"        \
  "ins              %["#TEMP3"],   %["#TEMP9"],  16,    16        \n\t"        \
  "precrq.ph.w      %["#TEMP4"],   %["#TEMP10"], %["#TEMP6"]      \n\t"        \
  "precrq.ph.w      %["#TEMP5"],   %["#TEMP11"], %["#TEMP7"]      \n\t"        \
  "ins              %["#TEMP6"],   %["#TEMP10"], 16,    16        \n\t"        \
  "ins              %["#TEMP7"],   %["#TEMP11"], 16,    16        \n\t"

#define MUL_SHIFT_SUM(TEMP0, TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6,         \
                      TEMP7, TEMP8, TEMP9, TEMP10, TEMP11, TEMP12, TEMP13,     \
                      TEMP14, TEMP15, TEMP16, TEMP17, TEMP18, TEMP19)          \
  "mul              %["#TEMP0"],   %["#TEMP8"],   %[kC2]          \n\t"        \
  "mul              %["#TEMP1"],   %["#TEMP8"],   %[kC1]          \n\t"        \
  "mul              %["#TEMP2"],   %["#TEMP9"],   %[kC2]          \n\t"        \
  "mul              %["#TEMP3"],   %["#TEMP9"],   %[kC1]          \n\t"        \
  "mul              %["#TEMP4"],   %["#TEMP10"],  %[kC2]          \n\t"        \
  "mul              %["#TEMP5"],   %["#TEMP10"],  %[kC1]          \n\t"        \
  "mul              %["#TEMP6"],   %["#TEMP11"],  %[kC2]          \n\t"        \
  "mul              %["#TEMP7"],   %["#TEMP11"],  %[kC1]          \n\t"        \
  "sra              %["#TEMP0"],   %["#TEMP0"],   16              \n\t"        \
  "sra              %["#TEMP1"],   %["#TEMP1"],   16              \n\t"        \
  "sra              %["#TEMP2"],   %["#TEMP2"],   16              \n\t"        \
  "sra              %["#TEMP3"],   %["#TEMP3"],   16              \n\t"        \
  "sra              %["#TEMP4"],   %["#TEMP4"],   16              \n\t"        \
  "sra              %["#TEMP5"],   %["#TEMP5"],   16              \n\t"        \
  "sra              %["#TEMP6"],   %["#TEMP6"],   16              \n\t"        \
  "sra              %["#TEMP7"],   %["#TEMP7"],   16              \n\t"        \
  "addu             %["#TEMP12"],  %["#TEMP12"],  %["#TEMP16"]    \n\t"        \
  "addu             %["#TEMP13"],  %["#TEMP13"],  %["#TEMP17"]    \n\t"        \
  "subu             %["#TEMP14"],  %["#TEMP14"],  %["#TEMP18"]    \n\t"        \
  "subu             %["#TEMP15"],  %["#TEMP15"],  %["#TEMP19"]    \n\t"

#define ADD_SUB_HALVES(TEMP0, TEMP1, TEMP2, TEMP3)                             \
  "addq.ph          %["#TEMP0"],   %["#TEMP2"],  %["#TEMP3"]      \n\t"        \
  "subq.ph          %["#TEMP1"],   %["#TEMP2"],  %["#TEMP3"]      \n\t"

#define LOAD_IN_X2(TEMP0, TEMP1, A, B)                                         \
  "lh               %["#TEMP0"],   "#A"(%[in])                \n\t"            \
  "lh               %["#TEMP1"],   "#B"(%[in])                \n\t"

#define SRA_16(TEMP0, TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6, TEMP7)         \
  "sra              %["#TEMP0"],  %["#TEMP4"],  16            \n\t"            \
  "sra              %["#TEMP1"],  %["#TEMP5"],  16            \n\t"            \
  "sra              %["#TEMP2"],  %["#TEMP6"],  16            \n\t"            \
  "sra              %["#TEMP3"],  %["#TEMP7"],  16            \n\t"

#define INSERT_HALF_X2(TEMP0, TEMP1, TEMP2, TEMP3)                             \
  "ins              %["#TEMP0"],   %["#TEMP2"], 16,    16     \n\t"            \
  "ins              %["#TEMP1"],   %["#TEMP3"], 16,    16     \n\t"

static void TransformDC(const int16_t* in, uint8_t* dst) {
  int temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10;

  __asm__ volatile (
    LOAD_DST(temp1, temp2, temp3, temp4)
    "lh               %[temp5],  0(%[in])               \n\t"
    "addiu            %[temp5],  %[temp5],  4           \n\t"
    "ins              %[temp5],  %[temp5],  16, 16      \n\t"
    "shra.ph          %[temp5],  %[temp5],  3           \n\t"
    PRECEU_RL(temp6, temp7, temp8, temp9, temp10, temp1, temp2, temp3,
              temp1, temp2, temp3, temp4)
    STORE_SAT_SUM_X2(temp6, temp7, temp8, temp9, temp10, temp1, temp2, temp3,
                     temp5, temp5, temp5, temp5, temp5, temp5, temp5, temp5)

    : [temp1]"=&r"(temp1), [temp2]"=&r"(temp2), [temp3]"=&r"(temp3),
      [temp4]"=&r"(temp4), [temp5]"=&r"(temp5), [temp6]"=&r"(temp6),
      [temp7]"=&r"(temp7), [temp8]"=&r"(temp8), [temp9]"=&r"(temp9),
      [temp10]"=&r"(temp10)
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
    PRECEU_RL(temp13, temp14, temp3, temp15, temp5, temp16, temp11, temp17,
              temp3, temp5, temp11, temp12)
    PRECRQ_X2(temp12, temp18, temp1, temp8, temp7, temp6, temp2, temp4,
              temp7, temp6, temp10, temp9)
    STORE_SAT_SUM_X2(temp13, temp14, temp3, temp15, temp5, temp16, temp11,
                     temp17, temp12, temp18, temp1, temp8, temp2, temp4,
                     temp7, temp6)

    : [temp1]"=&r"(temp1), [temp2]"=&r"(temp2), [temp3]"=&r"(temp3),
      [temp4]"=&r"(temp4), [temp5]"=&r"(temp5), [temp6]"=&r"(temp6),
      [temp7]"=&r"(temp7), [temp8]"=&r"(temp8), [temp9]"=&r"(temp9),
      [temp10]"=&r"(temp10), [temp11]"=&r"(temp11), [temp12]"=&r"(temp12),
      [temp13]"=&r"(temp13), [temp14]"=&r"(temp14), [temp15]"=&r"(temp15),
      [temp16]"=&r"(temp16), [temp17]"=&r"(temp17), [temp18]"=&r"(temp18),
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
                  temp1, temp2, temp5, temp6, temp10, temp8, temp9, temp7,
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
                  temp17, temp18, temp9, temp10, temp12, temp14, temp11, temp13,
                  temp15, temp4, temp16, temp17)
    INSERT_HALF_X2(temp11, temp12, temp13, temp14)
    ADD_SUB_HALVES(temp17, temp8, temp8, temp11)
    ADD_SUB_HALVES(temp3, temp4, temp7, temp12)

    // horizintal
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
                  temp9, temp13, temp11, temp15, temp3, temp13, temp1, temp9,
                  temp6, temp17, temp8, temp18)
    MUL_SHIFT_SUM(temp6, temp8, temp18, temp17, temp11, temp15, temp12, temp16,
                  temp12, temp16, temp10, temp14, temp8, temp15, temp6, temp11,
                  temp18, temp12, temp17, temp16)
    INSERT_HALF_X2(temp1, temp3, temp9, temp13)
    INSERT_HALF_X2(temp6, temp8, temp11, temp15)
    SHIFT_R_SUM_X2(temp9, temp10, temp11, temp12, temp13, temp14, temp15,
                   temp16, temp2, temp4, temp5, temp7, temp3, temp1, temp8,
                   temp6)
    PRECRQ_X2(temp1, temp2, temp9, temp12, temp3, temp4, temp13, temp16,
              temp11, temp10, temp15, temp14)
    LOAD_DST(temp10, temp11, temp14, temp15)
    PRECEU_RL(temp5, temp6, temp7, temp8, temp17, temp18, temp10, temp11,
              temp10, temp11, temp14, temp15)
    STORE_SAT_SUM_X2(temp5, temp6, temp7, temp8, temp17, temp18, temp10, temp11,
                     temp9, temp12, temp1, temp2, temp13, temp16, temp3, temp4)

    : [temp1]"=&r"(temp1), [temp2]"=&r"(temp2), [temp3]"=&r"(temp3),
      [temp4]"=&r"(temp4), [temp5]"=&r"(temp5), [temp6]"=&r"(temp6),
      [temp7]"=&r"(temp7), [temp8]"=&r"(temp8), [temp9]"=&r"(temp9),
      [temp10]"=&r"(temp10), [temp11]"=&r"(temp11), [temp12]"=&r"(temp12),
      [temp13]"=&r"(temp13), [temp14]"=&r"(temp14), [temp15]"=&r"(temp15),
      [temp16]"=&r"(temp16), [temp17]"=&r"(temp17), [temp18]"=&r"(temp18)
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

#undef INSERT_HALF_X2
#undef SRA_16
#undef LOAD_IN_X2
#undef ADD_SUB_HALVES
#undef MUL_SHIFT_SUM
#undef PRECRQ_X2
#undef LOAD_DST
#undef PRECEU_RL
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
#endif
}
