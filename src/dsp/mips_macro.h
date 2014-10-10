#ifndef WEBP_DSP_MIPS_MACRO_H_
#define WEBP_DSP_MIPS_MACRO_H_

#define ADD_SUB_HALVES(O0, O1,                                                 \
                       I0, I1)                                                 \
  "addq.ph          %["#O0"],   %["#I0"],  %["#I1"]           \n\t"            \
  "subq.ph          %["#O1"],   %["#I0"],  %["#I1"]           \n\t"

#define LOAD_IN_X2(O0, O1,                                                     \
                   I0, I1)                                                     \
  "lh               %["#O0"],   "#I0"(%[in])                  \n\t"            \
  "lh               %["#O1"],   "#I1"(%[in])                  \n\t"

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

 #define INSERT_HALF_X2(O0, O1,                                                \
                        I0, I1)                                                \
  "ins              %["#O0"],   %["#I0"], 16,    16           \n\t"            \
  "ins              %["#O1"],   %["#I1"], 16,    16           \n\t"

#define SRA_16(O0, O1, O2, O3,                                                 \
               I0, I1, I2, I3)                                                 \
  "sra              %["#O0"],  %["#I0"],  16                  \n\t"            \
  "sra              %["#O1"],  %["#I1"],  16                  \n\t"            \
  "sra              %["#O2"],  %["#I2"],  16                  \n\t"            \
  "sra              %["#O3"],  %["#I3"],  16                  \n\t"

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

#endif  // WEBP_DSP_MIPS_MACRO_H_