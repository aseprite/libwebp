/*
 * Copyright (c) 2013
 *      Imagination Technologies Limited and/or its affiliated group
 * companies ("Imagination").
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of Imagination, nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY IMAGINATION ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL IMAGINATION BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 * Author(s):  Djordje Pesut (djordje.pesut@imgtec.com)
 */

#include "./dsp.h"

#if defined(WEBP_USE_MIPS32)

static const int kC1 = 20091 + (1 << 16);
static const int kC2 = 35468;

static void TransformOneMIPS32(const int16_t* in, uint8_t* dst) {
  int temp0, temp1, temp2, temp3, temp4;
  int temp5, temp6, temp7, temp8, temp9;
  int temp10, temp11, temp12, temp13, temp14;
  int temp15, temp16, temp17, temp18;

  int16_t *p_in = (int16_t *)in;

  __asm__ volatile(
    "lh       %[temp0],  0(%[in])                      \n\t"
    "lh       %[temp8],  16(%[in])                     \n\t"
    "lh       %[temp4],  8(%[in])                      \n\t"
    "lh       %[temp12], 24(%[in])                     \n\t"
    "addu     %[temp16], %[temp0],  %[temp8]           \n\t"
    "subu     %[temp0],  %[temp0],  %[temp8]           \n\t"
    "mul      %[temp8],  %[temp4],  %[kC2]             \n\t"
    "mul      %[temp17], %[temp12], %[kC1]             \n\t"
    "mul      %[temp4],  %[temp4],  %[kC1]             \n\t"
    "mul      %[temp12], %[temp12], %[kC2]             \n\t"
    "lh       %[temp1],  2(%[in])                      \n\t"
    "lh       %[temp5],  10(%[in])                     \n\t"
    "lh       %[temp9],  18(%[in])                     \n\t"
    "lh       %[temp13], 26(%[in])                     \n\t"
    "sra      %[temp8],  %[temp8],  16                 \n\t"
    "sra      %[temp17], %[temp17], 16                 \n\t"
    "sra      %[temp4],  %[temp4],  16                 \n\t"
    "sra      %[temp12], %[temp12], 16                 \n\t"
    "lh       %[temp2],  4(%[in])                      \n\t"
    "lh       %[temp6],  12(%[in])                     \n\t"
    "lh       %[temp10], 20(%[in])                     \n\t"
    "lh       %[temp14], 28(%[in])                     \n\t"
    "subu     %[temp17], %[temp8],  %[temp17]          \n\t"
    "addu     %[temp4],  %[temp4],  %[temp12]          \n\t"
    "addu     %[temp8],  %[temp16], %[temp4]           \n\t"
    "subu     %[temp4],  %[temp16], %[temp4]           \n\t"
    "addu     %[temp16], %[temp1],  %[temp9]           \n\t"
    "subu     %[temp1],  %[temp1],  %[temp9]           \n\t"
    "lh       %[temp3],  6(%[in])                      \n\t"
    "lh       %[temp7],  14(%[in])                     \n\t"
    "lh       %[temp11], 22(%[in])                     \n\t"
    "lh       %[temp15], 30(%[in])                     \n\t"
    "addu     %[temp12], %[temp0],  %[temp17]          \n\t"
    "subu     %[temp0],  %[temp0],  %[temp17]          \n\t"
    "mul      %[temp9],  %[temp5],  %[kC2]             \n\t"
    "mul      %[temp17], %[temp13], %[kC1]             \n\t"
    "mul      %[temp5],  %[temp5],  %[kC1]             \n\t"
    "mul      %[temp13], %[temp13], %[kC2]             \n\t"
    "sra      %[temp9],  %[temp9],  16                 \n\t"
    "sra      %[temp17], %[temp17], 16                 \n\t"
    "subu     %[temp17], %[temp9],  %[temp17]          \n\t"
    "sra      %[temp5],  %[temp5],  16                 \n\t"
    "sra      %[temp13], %[temp13], 16                 \n\t"
    "addu     %[temp5],  %[temp5],  %[temp13]          \n\t"
    "addu     %[temp13], %[temp1],  %[temp17]          \n\t"
    "subu     %[temp1],  %[temp1],  %[temp17]          \n\t"
    "mul      %[temp17], %[temp14], %[kC1]             \n\t"
    "mul      %[temp14], %[temp14], %[kC2]             \n\t"
    "addu     %[temp9],  %[temp16], %[temp5]           \n\t"
    "subu     %[temp5],  %[temp16], %[temp5]           \n\t"
    "addu     %[temp16], %[temp2],  %[temp10]          \n\t"
    "subu     %[temp2],  %[temp2],  %[temp10]          \n\t"
    "mul      %[temp10], %[temp6],  %[kC2]             \n\t"
    "mul      %[temp6],  %[temp6],  %[kC1]             \n\t"
    "sra      %[temp17], %[temp17], 16                 \n\t"
    "sra      %[temp14], %[temp14], 16                 \n\t"
    "sra      %[temp10], %[temp10], 16                 \n\t"
    "sra      %[temp6],  %[temp6],  16                 \n\t"
    "subu     %[temp17], %[temp10], %[temp17]          \n\t"
    "addu     %[temp6],  %[temp6],  %[temp14]          \n\t"
    "addu     %[temp10], %[temp16], %[temp6]           \n\t"
    "subu     %[temp6],  %[temp16], %[temp6]           \n\t"
    "addu     %[temp14], %[temp2],  %[temp17]          \n\t"
    "subu     %[temp2],  %[temp2],  %[temp17]          \n\t"
    "mul      %[temp17], %[temp15], %[kC1]             \n\t"
    "mul      %[temp15], %[temp15], %[kC2]             \n\t"
    "addu     %[temp16], %[temp3],  %[temp11]          \n\t"
    "subu     %[temp3],  %[temp3],  %[temp11]          \n\t"
    "mul      %[temp11], %[temp7],  %[kC2]             \n\t"
    "mul      %[temp7],  %[temp7],  %[kC1]             \n\t"
    "addiu    %[temp8],  %[temp8],  4                  \n\t"
    "addiu    %[temp12], %[temp12], 4                  \n\t"
    "addiu    %[temp0],  %[temp0],  4                  \n\t"
    "addiu    %[temp4],  %[temp4],  4                  \n\t"
    "sra      %[temp17], %[temp17], 16                 \n\t"
    "sra      %[temp15], %[temp15], 16                 \n\t"
    "sra      %[temp11], %[temp11], 16                 \n\t"
    "sra      %[temp7],  %[temp7],  16                 \n\t"
    "subu     %[temp17], %[temp11], %[temp17]          \n\t"
    "addu     %[temp7],  %[temp7],  %[temp15]          \n\t"
    "addu     %[temp15], %[temp3],  %[temp17]          \n\t"
    "subu     %[temp3],  %[temp3],  %[temp17]          \n\t"
    "addu     %[temp11], %[temp16], %[temp7]           \n\t"
    "subu     %[temp7],  %[temp16], %[temp7]           \n\t"
    "addu     %[temp16], %[temp8],  %[temp10]          \n\t"
    "subu     %[temp8],  %[temp8],  %[temp10]          \n\t"
    "mul      %[temp10], %[temp9],  %[kC2]             \n\t"
    "mul      %[temp17], %[temp11], %[kC1]             \n\t"
    "mul      %[temp9],  %[temp9],  %[kC1]             \n\t"
    "mul      %[temp11], %[temp11], %[kC2]             \n\t"
    "sra      %[temp10], %[temp10], 16                 \n\t"
    "sra      %[temp17], %[temp17], 16                 \n\t"
    "sra      %[temp9],  %[temp9],  16                 \n\t"
    "sra      %[temp11], %[temp11], 16                 \n\t"
    "subu     %[temp17], %[temp10], %[temp17]          \n\t"
    "addu     %[temp11], %[temp9],  %[temp11]          \n\t"
    "addu     %[temp10], %[temp12], %[temp14]          \n\t"
    "subu     %[temp12], %[temp12], %[temp14]          \n\t"
    "mul      %[temp14], %[temp13], %[kC2]             \n\t"
    "mul      %[temp9],  %[temp15], %[kC1]             \n\t"
    "mul      %[temp13], %[temp13], %[kC1]             \n\t"
    "mul      %[temp15], %[temp15], %[kC2]             \n\t"
    "sra      %[temp14], %[temp14], 16                 \n\t"
    "sra      %[temp9],  %[temp9],  16                 \n\t"
    "sra      %[temp13], %[temp13], 16                 \n\t"
    "sra      %[temp15], %[temp15], 16                 \n\t"
    "subu     %[temp9],  %[temp14], %[temp9]           \n\t"
    "addu     %[temp15], %[temp13], %[temp15]          \n\t"
    "addu     %[temp14], %[temp0],  %[temp2]           \n\t"
    "subu     %[temp0],  %[temp0],  %[temp2]           \n\t"
    "mul      %[temp2],  %[temp1],  %[kC2]             \n\t"
    "mul      %[temp13], %[temp3],  %[kC1]             \n\t"
    "mul      %[temp1],  %[temp1],  %[kC1]             \n\t"
    "mul      %[temp3],  %[temp3],  %[kC2]             \n\t"
    "sra      %[temp2],  %[temp2],  16                 \n\t"
    "sra      %[temp13], %[temp13], 16                 \n\t"
    "sra      %[temp1],  %[temp1],  16                 \n\t"
    "sra      %[temp3],  %[temp3],  16                 \n\t"
    "subu     %[temp13], %[temp2],  %[temp13]          \n\t"
    "addu     %[temp3],  %[temp1],  %[temp3]           \n\t"
    "addu     %[temp2],  %[temp4],  %[temp6]           \n\t"
    "subu     %[temp4],  %[temp4],  %[temp6]           \n\t"
    "mul      %[temp6],  %[temp5],  %[kC2]             \n\t"
    "mul      %[temp1],  %[temp7],  %[kC1]             \n\t"
    "mul      %[temp5],  %[temp5],  %[kC1]             \n\t"
    "mul      %[temp7],  %[temp7],  %[kC2]             \n\t"
    "sra      %[temp6],  %[temp6],  16                 \n\t"
    "sra      %[temp1],  %[temp1],  16                 \n\t"
    "sra      %[temp5],  %[temp5],  16                 \n\t"
    "sra      %[temp7],  %[temp7],  16                 \n\t"
    "subu     %[temp1],  %[temp6],  %[temp1]           \n\t"
    "addu     %[temp7],  %[temp5],  %[temp7]           \n\t"
    "addu     %[temp5],  %[temp16], %[temp11]          \n\t"
    "subu     %[temp16], %[temp16], %[temp11]          \n\t"
    "addu     %[temp11], %[temp8],  %[temp17]          \n\t"
    "subu     %[temp8],  %[temp8],  %[temp17]          \n\t"
    "sra      %[temp5],  %[temp5],  3                  \n\t"
    "sra      %[temp16], %[temp16], 3                  \n\t"
    "sra      %[temp11], %[temp11], 3                  \n\t"
    "sra      %[temp8],  %[temp8],  3                  \n\t"
    "addu     %[temp17], %[temp10], %[temp15]          \n\t"
    "subu     %[temp10], %[temp10], %[temp15]          \n\t"
    "addu     %[temp15], %[temp12], %[temp9]           \n\t"
    "subu     %[temp12], %[temp12], %[temp9]           \n\t"
    "sra      %[temp17], %[temp17], 3                  \n\t"
    "sra      %[temp10], %[temp10], 3                  \n\t"
    "sra      %[temp15], %[temp15], 3                  \n\t"
    "sra      %[temp12], %[temp12], 3                  \n\t"
    "addu     %[temp9],  %[temp14], %[temp3]           \n\t"
    "subu     %[temp14], %[temp14], %[temp3]           \n\t"
    "addu     %[temp3],  %[temp0],  %[temp13]          \n\t"
    "subu     %[temp0],  %[temp0],  %[temp13]          \n\t"
    "sra      %[temp9],  %[temp9],  3                  \n\t"
    "sra      %[temp14], %[temp14], 3                  \n\t"
    "sra      %[temp3],  %[temp3],  3                  \n\t"
    "sra      %[temp0],  %[temp0],  3                  \n\t"
    "addu     %[temp13], %[temp2],  %[temp7]           \n\t"
    "subu     %[temp2],  %[temp2],  %[temp7]           \n\t"
    "addu     %[temp7],  %[temp4],  %[temp1]           \n\t"
    "subu     %[temp4],  %[temp4],  %[temp1]           \n\t"
    "sra      %[temp13], %[temp13], 3                  \n\t"
    "sra      %[temp2],  %[temp2],  3                  \n\t"
    "sra      %[temp7],  %[temp7],  3                  \n\t"
    "sra      %[temp4],  %[temp4],  3                  \n\t"
    "addiu    %[temp6],  $zero,     255                \n\t"
    "lbu      %[temp1],  0(%[dst])                     \n\t"
    "addu     %[temp1],  %[temp1],  %[temp5]           \n\t"
    "sra      %[temp5],  %[temp1],  8                  \n\t"
    "sra      %[temp18], %[temp1],  31                 \n\t"
    "beqz     %[temp5],  1f                            \n\t"
    "xor      %[temp1],  %[temp1],  %[temp1]           \n\t"
    "movz     %[temp1],  %[temp6],  %[temp18]          \n\t"
  "1:                                                  \n\t"
    "lbu      %[temp18], 1(%[dst])                     \n\t"
    "sb       %[temp1],  0(%[dst])                     \n\t"
    "addu     %[temp18], %[temp18], %[temp11]          \n\t"
    "sra      %[temp11], %[temp18], 8                  \n\t"
    "sra      %[temp1],  %[temp18], 31                 \n\t"
    "beqz     %[temp11], 2f                            \n\t"
    "xor      %[temp18], %[temp18], %[temp18]          \n\t"
    "movz     %[temp18], %[temp6],  %[temp1]           \n\t"
  "2:                                                  \n\t"
    "lbu      %[temp1],  2(%[dst])                     \n\t"
    "sb       %[temp18], 1(%[dst])                     \n\t"
    "addu     %[temp1],  %[temp1],  %[temp8]           \n\t"
    "sra      %[temp8],  %[temp1],  8                  \n\t"
    "sra      %[temp18], %[temp1],  31                 \n\t"
    "beqz     %[temp8],  3f                            \n\t"
    "xor      %[temp1],  %[temp1],  %[temp1]           \n\t"
    "movz     %[temp1],  %[temp6],  %[temp18]          \n\t"
  "3:                                                  \n\t"
    "lbu      %[temp18], 3(%[dst])                     \n\t"
    "sb       %[temp1],  2(%[dst])                     \n\t"
    "addu     %[temp18], %[temp18], %[temp16]          \n\t"
    "sra      %[temp16], %[temp18], 8                  \n\t"
    "sra      %[temp1],  %[temp18], 31                 \n\t"
    "beqz     %[temp16], 4f                            \n\t"
    "xor      %[temp18], %[temp18], %[temp18]          \n\t"
    "movz     %[temp18], %[temp6],  %[temp1]           \n\t"
  "4:                                                  \n\t"
    "sb       %[temp18], 3(%[dst])                     \n\t"
    "lbu      %[temp5],  32(%[dst])                    \n\t"
    "lbu      %[temp8],  33(%[dst])                    \n\t"
    "lbu      %[temp11], 34(%[dst])                    \n\t"
    "lbu      %[temp16], 35(%[dst])                    \n\t"
    "addu     %[temp5],  %[temp5],  %[temp17]          \n\t"
    "addu     %[temp8],  %[temp8],  %[temp15]          \n\t"
    "addu     %[temp11], %[temp11], %[temp12]          \n\t"
    "addu     %[temp16], %[temp16], %[temp10]          \n\t"
    "sra      %[temp18], %[temp5],  8                  \n\t"
    "sra      %[temp1],  %[temp5],  31                 \n\t"
    "beqz     %[temp18], 5f                            \n\t"
    "xor      %[temp5],  %[temp5],  %[temp5]           \n\t"
    "movz     %[temp5],  %[temp6],  %[temp1]           \n\t"
  "5:                                                  \n\t"
    "sra      %[temp18], %[temp8],  8                  \n\t"
    "sra      %[temp1],  %[temp8],  31                 \n\t"
    "beqz     %[temp18], 6f                            \n\t"
    "xor      %[temp8],  %[temp8],  %[temp8]           \n\t"
    "movz     %[temp8],  %[temp6],  %[temp1]           \n\t"
  "6:                                                  \n\t"
    "sra      %[temp18], %[temp11], 8                  \n\t"
    "sra      %[temp1],  %[temp11], 31                 \n\t"
    "sra      %[temp17], %[temp16], 8                  \n\t"
    "sra      %[temp15], %[temp16], 31                 \n\t"
    "beqz     %[temp18], 7f                            \n\t"
    "xor      %[temp11], %[temp11], %[temp11]          \n\t"
    "movz     %[temp11], %[temp6],  %[temp1]           \n\t"
  "7:                                                  \n\t"
    "beqz     %[temp17], 8f                            \n\t"
    "xor      %[temp16], %[temp16], %[temp16]          \n\t"
    "movz     %[temp16], %[temp6],  %[temp15]          \n\t"
  "8:                                                  \n\t"
    "sb       %[temp5],  32(%[dst])                    \n\t"
    "sb       %[temp8],  33(%[dst])                    \n\t"
    "sb       %[temp11], 34(%[dst])                    \n\t"
    "sb       %[temp16], 35(%[dst])                    \n\t"
    "lbu      %[temp5],  64(%[dst])                    \n\t"
    "lbu      %[temp8],  65(%[dst])                    \n\t"
    "lbu      %[temp11], 66(%[dst])                    \n\t"
    "lbu      %[temp16], 67(%[dst])                    \n\t"
    "addu     %[temp5],  %[temp5],  %[temp9]           \n\t"
    "addu     %[temp8],  %[temp8],  %[temp3]           \n\t"
    "addu     %[temp11], %[temp11], %[temp0]           \n\t"
    "addu     %[temp16], %[temp16], %[temp14]          \n\t"
    "sra      %[temp18], %[temp5],  8                  \n\t"
    "sra      %[temp1],  %[temp5],  31                 \n\t"
    "sra      %[temp17], %[temp8],  8                  \n\t"
    "sra      %[temp15], %[temp8],  31                 \n\t"
    "sra      %[temp12], %[temp11], 8                  \n\t"
    "sra      %[temp10], %[temp11], 31                 \n\t"
    "sra      %[temp9],  %[temp16], 8                  \n\t"
    "sra      %[temp3],  %[temp16], 31                 \n\t"
    "beqz     %[temp18], 9f                            \n\t"
    "xor      %[temp5],  %[temp5],  %[temp5]           \n\t"
    "movz     %[temp5],  %[temp6],  %[temp1]           \n\t"
  "9:                                                  \n\t"
    "beqz     %[temp17], 10f                           \n\t"
    "xor      %[temp8],  %[temp8],  %[temp8]           \n\t"
    "movz     %[temp8],  %[temp6],  %[temp15]          \n\t"
  "10:                                                 \n\t"
    "beqz     %[temp12], 11f                           \n\t"
    "xor      %[temp11], %[temp11], %[temp11]          \n\t"
    "movz     %[temp11], %[temp6],  %[temp10]          \n\t"
  "11:                                                 \n\t"
    "beqz     %[temp9],  12f                           \n\t"
    "xor      %[temp16], %[temp16], %[temp16]          \n\t"
    "movz     %[temp16], %[temp6],  %[temp3]           \n\t"
  "12:                                                 \n\t"
    "sb       %[temp5],  64(%[dst])                    \n\t"
    "sb       %[temp8],  65(%[dst])                    \n\t"
    "sb       %[temp11], 66(%[dst])                    \n\t"
    "sb       %[temp16], 67(%[dst])                    \n\t"
    "lbu      %[temp5],  96(%[dst])                    \n\t"
    "lbu      %[temp8],  97(%[dst])                    \n\t"
    "lbu      %[temp11], 98(%[dst])                    \n\t"
    "lbu      %[temp16], 99(%[dst])                    \n\t"
    "addu     %[temp5],  %[temp5],  %[temp13]          \n\t"
    "addu     %[temp8],  %[temp8],  %[temp7]           \n\t"
    "addu     %[temp11], %[temp11], %[temp4]           \n\t"
    "addu     %[temp16], %[temp16], %[temp2]           \n\t"
    "sra      %[temp18], %[temp5],  8                  \n\t"
    "sra      %[temp1],  %[temp5],  31                 \n\t"
    "sra      %[temp17], %[temp8],  8                  \n\t"
    "sra      %[temp15], %[temp8],  31                 \n\t"
    "sra      %[temp12], %[temp11], 8                  \n\t"
    "sra      %[temp10], %[temp11], 31                 \n\t"
    "sra      %[temp9],  %[temp16], 8                  \n\t"
    "sra      %[temp3],  %[temp16], 31                 \n\t"
    "beqz     %[temp18], 13f                           \n\t"
    "xor      %[temp5],  %[temp5],  %[temp5]           \n\t"
    "movz     %[temp5],  %[temp6],  %[temp1]           \n\t"
  "13:                                                 \n\t"
    "beqz     %[temp17], 14f                           \n\t"
    "xor      %[temp8],  %[temp8],  %[temp8]           \n\t"
    "movz     %[temp8],  %[temp6],  %[temp15]          \n\t"
  "14:                                                 \n\t"
    "beqz     %[temp12], 15f                           \n\t"
    "xor      %[temp11], %[temp11], %[temp11]          \n\t"
    "movz     %[temp11], %[temp6],  %[temp10]          \n\t"
  "15:                                                 \n\t"
    "beqz     %[temp9],  16f                           \n\t"
    "xor      %[temp16], %[temp16], %[temp16]          \n\t"
    "movz     %[temp16], %[temp6],  %[temp3]           \n\t"
  "16:                                                 \n\t"
    "sb       %[temp5],  96(%[dst])                    \n\t"
    "sb       %[temp8],  97(%[dst])                    \n\t"
    "sb       %[temp11], 98(%[dst])                    \n\t"
    "sb       %[temp16], 99(%[dst])                    \n\t"

    : [temp0]"=&r"(temp0), [temp1]"=&r"(temp1), [temp2]"=&r"(temp2),
      [temp3]"=&r"(temp3), [temp4]"=&r"(temp4), [temp5]"=&r"(temp5),
      [temp6]"=&r"(temp6), [temp7]"=&r"(temp7), [temp8]"=&r"(temp8),
      [temp9]"=&r"(temp9), [temp10]"=&r"(temp10), [temp11]"=&r"(temp11),
      [temp12]"=&r"(temp12), [temp13]"=&r"(temp13), [temp14]"=&r"(temp14),
      [temp15]"=&r"(temp15), [temp16]"=&r"(temp16), [temp17]"=&r"(temp17),
      [temp18]"=&r"(temp18)
    : [in]"r"(p_in), [kC1]"r"(kC1), [kC2]"r"(kC2), [dst]"r"(dst)
    : "memory", "hi", "lo"
  );
}

static void TransformTwoMIPS32(const int16_t* in, uint8_t* dst, int do_two) {
  TransformOneMIPS32(in, dst);
  if (do_two) {
    TransformOneMIPS32(in + 16, dst + 4);
  }
}

#endif /* WEBP_USE_MIPS32 */

/*------------------------------------------------------------------------------
 Entry point */

extern void VP8DspInitMIPS32(void);

void VP8DspInitMIPS32(void) {
#if defined(WEBP_USE_MIPS32)
  VP8Transform = TransformTwoMIPS32;
#endif /* WEBP_USE_MIPS32 */
}
