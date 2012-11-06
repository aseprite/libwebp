// Copyright 2012 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// ARM NEON version of speed-critical encoding functions.
//
// adapted from libvpx (http://www.webmproject.org/code/)

#include "./dsp.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#if defined(WEBP_USE_NEON)

#include "../enc/vp8enci.h"

//------------------------------------------------------------------------------
// Transforms (Paragraph 14.4)

// Inverse transform.
// This code is pretty much the same as TransformOneNEON in the decoder, except
// for subtraction to *ref. See the comments there for algorithmic explanations.
static void ITransformOne(const uint8_t* ref,
                          const int16_t* in, uint8_t* dst) {
  const int kBPS = BPS;
  const int16_t kC1C2[] = { 20091, 17734, 0, 0 };  // kC1 / (kC2 >> 1) / 0 / 0

  __asm__ volatile (
    "vld1.16         {q1, q2}, [%[in]]           \n"
    "vld1.16         {d0}, [%[kC1C2]]            \n"

    // d2: in[0]
    // d3: in[8]
    // d4: in[4]
    // d5: in[12]
    "vswp            d3, d4                      \n"

    // q8 = {in[4], in[12]} * kC1 * 2 >> 16
    // q9 = {in[4], in[12]} * kC2 >> 16
    "vqdmulh.s16     q8, q2, d0[0]               \n"
    "vqdmulh.s16     q9, q2, d0[1]               \n"

    // d22 = a = in[0] + in[8]
    // d23 = b = in[0] - in[8]
    "vqadd.s16       d22, d2, d3                 \n"
    "vqsub.s16       d23, d2, d3                 \n"

    //  q8 = in[4]/[12] * kC1 >> 16
    "vshr.s16        q8, q8, #1                  \n"

    // Add {in[4], in[12]} back after the multiplication.
    "vqadd.s16       q8, q2, q8                  \n"

    // d20 = c = in[4]*kC2 - in[12]*kC1
    // d21 = d = in[4]*kC1 + in[12]*kC2
    "vqsub.s16       d20, d18, d17               \n"
    "vqadd.s16       d21, d19, d16               \n"

    // d2 = tmp[0] = a + d
    // d3 = tmp[1] = b + c
    // d4 = tmp[2] = b - c
    // d5 = tmp[3] = a - d
    "vqadd.s16       d2, d22, d21                \n"
    "vqadd.s16       d3, d23, d20                \n"
    "vqsub.s16       d4, d23, d20                \n"
    "vqsub.s16       d5, d22, d21                \n"

    "vzip.16         q1, q2                      \n"
    "vzip.16         q1, q2                      \n"

    "vswp            d3, d4                      \n"

    // q8 = {tmp[4], tmp[12]} * kC1 * 2 >> 16
    // q9 = {tmp[4], tmp[12]} * kC2 >> 16
    "vqdmulh.s16     q8, q2, d0[0]               \n"
    "vqdmulh.s16     q9, q2, d0[1]               \n"

    // d22 = a = tmp[0] + tmp[8]
    // d23 = b = tmp[0] - tmp[8]
    "vqadd.s16       d22, d2, d3                 \n"
    "vqsub.s16       d23, d2, d3                 \n"

    "vshr.s16        q8, q8, #1                  \n"
    "vqadd.s16       q8, q2, q8                  \n"

    // d20 = c = in[4]*kC2 - in[12]*kC1
    // d21 = d = in[4]*kC1 + in[12]*kC2
    "vqsub.s16       d20, d18, d17               \n"
    "vqadd.s16       d21, d19, d16               \n"

    // d2 = tmp[0] = a + d
    // d3 = tmp[1] = b + c
    // d4 = tmp[2] = b - c
    // d5 = tmp[3] = a - d
    "vqadd.s16       d2, d22, d21                \n"
    "vqadd.s16       d3, d23, d20                \n"
    "vqsub.s16       d4, d23, d20                \n"
    "vqsub.s16       d5, d22, d21                \n"

    "vld1.32         d6[0], [%[ref]], %[kBPS]    \n"
    "vld1.32         d6[1], [%[ref]], %[kBPS]    \n"
    "vld1.32         d7[0], [%[ref]], %[kBPS]    \n"
    "vld1.32         d7[1], [%[ref]], %[kBPS]    \n"

    "sub         %[ref], %[ref], %[kBPS], lsl #2 \n"

    // (val) + 4 >> 3
    "vrshr.s16       d2, d2, #3                  \n"
    "vrshr.s16       d3, d3, #3                  \n"
    "vrshr.s16       d4, d4, #3                  \n"
    "vrshr.s16       d5, d5, #3                  \n"

    "vzip.16         q1, q2                      \n"
    "vzip.16         q1, q2                      \n"

    // Must accumulate before saturating
    "vmovl.u8        q8, d6                      \n"
    "vmovl.u8        q9, d7                      \n"

    "vqadd.s16       q1, q1, q8                  \n"
    "vqadd.s16       q2, q2, q9                  \n"

    "vqmovun.s16     d0, q1                      \n"
    "vqmovun.s16     d1, q2                      \n"

    "vst1.32         d0[0], [%[dst]], %[kBPS]    \n"
    "vst1.32         d0[1], [%[dst]], %[kBPS]    \n"
    "vst1.32         d1[0], [%[dst]], %[kBPS]    \n"
    "vst1.32         d1[1], [%[dst]]             \n"

    : [in] "+r"(in), [dst] "+r"(dst)               // modified registers
    : [kBPS] "r"(kBPS), [kC1C2] "r"(kC1C2), [ref] "r"(ref)  // constants
    : "memory", "q0", "q1", "q2", "q8", "q9", "q10", "q11"  // clobbered
  );
}

static void ITransform(const uint8_t* ref,
                       const int16_t* in, uint8_t* dst, int do_two) {
  ITransformOne(ref, in, dst);
  if (do_two) {
    ITransformOne(ref + 4, in + 16, dst + 4);
  }
}

// Same code as dec_neon.c
static void ITransformWHT(const int16_t* in, int16_t* out) {
  const int kStep = 32; // The store is only incrementing the pointer as if we
                        // had stored a single byte.
  __asm__ volatile (
    // part 1
    // load data into q0, q1
    "vld1.16         {q0, q1}, [%[in]]           \n"

    "vaddl.s16       q2, d0, d3                  \n" // a0 = in[0] + in[12]
    "vaddl.s16       q3, d1, d2                  \n" // a1 = in[4] + in[8]
    "vsubl.s16       q4, d1, d2                  \n" // a2 = in[4] - in[8]
    "vsubl.s16       q5, d0, d3                  \n" // a3 = in[0] - in[12]

    "vadd.s32        q0, q2, q3                  \n" // tmp[0] = a0 + a1
    "vsub.s32        q2, q2, q3                  \n" // tmp[8] = a0 - a1
    "vadd.s32        q1, q5, q4                  \n" // tmp[4] = a3 + a2
    "vsub.s32        q3, q5, q4                  \n" // tmp[12] = a3 - a2

    // Transpose
    // q0 = tmp[0, 4, 8, 12], q1 = tmp[2, 6, 10, 14]
    // q2 = tmp[1, 5, 9, 13], q3 = tmp[3, 7, 11, 15]
    "vswp            d1, d4                      \n" // vtrn.64 q0, q2
    "vswp            d3, d6                      \n" // vtrn.64 q1, q3
    "vtrn.32         q0, q1                      \n"
    "vtrn.32         q2, q3                      \n"

    "vmov.s32        q4, #3                      \n" // dc = 3
    "vadd.s32        q0, q0, q4                  \n" // dc = tmp[0] + 3
    "vadd.s32        q6, q0, q3                  \n" // a0 = dc + tmp[3]
    "vadd.s32        q7, q1, q2                  \n" // a1 = tmp[1] + tmp[2]
    "vsub.s32        q8, q1, q2                  \n" // a2 = tmp[1] - tmp[2]
    "vsub.s32        q9, q0, q3                  \n" // a3 = dc - tmp[3]

    "vadd.s32        q0, q6, q7                  \n"
    "vshrn.s32       d0, q0, #3                  \n" // (a0 + a1) >> 3
    "vadd.s32        q1, q9, q8                  \n"
    "vshrn.s32       d1, q1, #3                  \n" // (a3 + a2) >> 3
    "vsub.s32        q2, q6, q7                  \n"
    "vshrn.s32       d2, q2, #3                  \n" // (a0 - a1) >> 3
    "vsub.s32        q3, q9, q8                  \n"
    "vshrn.s32       d3, q3, #3                  \n" // (a3 - a2) >> 3

    // set the results to output
    "vst1.16         d0[0], [%[out]], %[kStep]      \n"
    "vst1.16         d1[0], [%[out]], %[kStep]      \n"
    "vst1.16         d2[0], [%[out]], %[kStep]      \n"
    "vst1.16         d3[0], [%[out]], %[kStep]      \n"
    "vst1.16         d0[1], [%[out]], %[kStep]      \n"
    "vst1.16         d1[1], [%[out]], %[kStep]      \n"
    "vst1.16         d2[1], [%[out]], %[kStep]      \n"
    "vst1.16         d3[1], [%[out]], %[kStep]      \n"
    "vst1.16         d0[2], [%[out]], %[kStep]      \n"
    "vst1.16         d1[2], [%[out]], %[kStep]      \n"
    "vst1.16         d2[2], [%[out]], %[kStep]      \n"
    "vst1.16         d3[2], [%[out]], %[kStep]      \n"
    "vst1.16         d0[3], [%[out]], %[kStep]      \n"
    "vst1.16         d1[3], [%[out]], %[kStep]      \n"
    "vst1.16         d2[3], [%[out]], %[kStep]      \n"
    "vst1.16         d3[3], [%[out]], %[kStep]      \n"

    : [out] "+r"(out)  // modified registers
    : [in] "r"(in), [kStep] "r"(kStep)  // constants
    : "memory", "q0", "q1", "q2", "q3", "q4",
      "q5", "q6", "q7", "q8", "q9" // clobbered
  );
}

// Forward transform.

// adapted from vp8/encoder/arm/neon/shortfdct_neon.asm
static const int16_t kCoeff16[] = {
  5352,  5352,  5352, 5352, 2217,  2217,  2217, 2217
};
static const int32_t kCoeff32[] = {
  14500, 14500, 14500, 14500,
  7500,   7500,  7500,  7500,
  12000, 12000, 12000, 12000,
  51000, 51000, 51000, 51000
};

static void FTransform(const uint8_t* src, const uint8_t* ref,
                       int16_t* out) {
 const int kBPS = BPS;
  const uint8_t* src_ptr = src;
  const uint8_t* ref_ptr = ref;
  int16_t* coeff16 = kCoeff16;
  int32_t* coeff32 = kCoeff32;

  __asm__ volatile (
    // load src into q4, q5 in high half
    "vld1.8 {d8},  [%[src_ptr]], %[kBPS]      \n"
    "vld1.8 {d10}, [%[src_ptr]], %[kBPS]      \n"
    "vld1.8 {d9},  [%[src_ptr]], %[kBPS]      \n"
    "vld1.8 {d11}, [%[src_ptr]]               \n"

    // load ref into q6, q7 in high half
    "vld1.8 {d12}, [%[ref_ptr]], %[kBPS]      \n"
    "vld1.8 {d14}, [%[ref_ptr]], %[kBPS]      \n"
    "vld1.8 {d13}, [%[ref_ptr]], %[kBPS]      \n"
    "vld1.8 {d15}, [%[ref_ptr]]               \n"

    // Pack the high values in to q4 and q6
    "vtrn.32     q4, q5                       \n"
    "vtrn.32     q6, q7                       \n"

    // {q0, q1} = q4 - q6
    "vsubl.u8    q0, d8, d12                  \n"
    "vsubl.u8    q1, d9, d13                  \n"

    // load coeff16 into q8(d16=5352, d17=2217)
    "vld1.16     {q8}, [%[coeff16]]           \n"

    // load coeff32 high half into q9 = 14500, q10 = 7500
    "vld1.32     {q9, q10}, [%[coeff32]]!     \n"

    // load coeff32 low half into q11=12000, q12=51000
    "vld1.32     {q11,q12}, [%[coeff32]]      \n"

    // part 1
    // transpose d0=ip[0], d1=ip[1], d2=ip[2], d3=ip[3]
    "vtrn.32         d0, d2                   \n"
    "vtrn.32         d1, d3                   \n"
    "vtrn.16         d0, d1                   \n"
    "vtrn.16         d2, d3                   \n"

    "vadd.s16        d4, d0, d3               \n" // a1 = ip[0] + ip[3]
    "vadd.s16        d5, d1, d2               \n" // b1 = ip[1] + ip[2]
    "vsub.s16        d6, d1, d2               \n" // c1 = ip[1] - ip[2]
    "vsub.s16        d7, d0, d3               \n" // d1 = ip[0] - ip[3]

    "vshl.s16        q2, q2, #3               \n" // (a1, b1) << 3
    "vshl.s16        q3, q3, #3               \n" // (c1, d1) << 3

    "vadd.s16        d0, d4, d5               \n" // op[0] = a1 + b1
    "vsub.s16        d2, d4, d5               \n" // op[2] = a1 - b1
    "vmlal.s16       q9, d7, d16              \n" // d1*5352 + 14500
    "vmlal.s16       q10, d7, d17             \n" // d1*2217 + 7500
    "vmlal.s16       q9, d6, d17              \n" // c1*2217 + d1*5352 + 14500
    "vmlsl.s16       q10, d6, d16             \n" // d1*2217 - c1*5352 + 7500

    // op[1] = (c1*2217 + d1*5352 + 14500) >> 12
    // op[3] = (d1*2217 - c1*5352 +  7500) >> 12
    "vshrn.s32       d1, q9, #12              \n"
    "vshrn.s32       d3, q10, #12             \n"

    // part 2
    // transpose d0=ip[0], d1=ip[4], d2=ip[8], d3=ip[12]
    "vtrn.32         d0, d2                   \n"
    "vtrn.32         d1, d3                   \n"
    "vtrn.16         d0, d1                   \n"
    "vtrn.16         d2, d3                   \n"

    "vmov.s16        d26, #7                  \n"

    "vadd.s16        d4, d0, d3               \n" // a1 = ip[0] + ip[12]
    "vadd.s16        d5, d1, d2               \n" // b1 = ip[4] + ip[8]
    "vsub.s16        d6, d1, d2               \n" // c1 = ip[4] - ip[8]
    "vadd.s16        d4, d4, d26              \n" // a1 + 7
    "vsub.s16        d7, d0, d3               \n" // d1 = ip[0] - ip[12]

    "vadd.s16        d0, d4, d5               \n" // op[0] = a1 + b1 + 7
    "vsub.s16        d2, d4, d5               \n" // op[8] = a1 - b1 + 7

    "vmlal.s16       q11, d7, d16             \n" // d1*5352 + 12000
    "vmlal.s16       q12, d7, d17             \n" // d1*2217 + 51000

    "vceq.s16        d4, d7, #0               \n"

    "vshr.s16        d0, d0, #4               \n"
    "vshr.s16        d2, d2, #4               \n"

    "vmlal.s16       q11, d6, d17             \n" // c1*2217 + d1*5352 + 12000
    "vmlsl.s16       q12, d6, d16             \n" // d1*2217 - c1*5352 + 51000

    "vmvn.s16        d4, d4                   \n"
    // op[4] = (c1*2217 + d1*5352 + 12000)>>16
    "vshrn.s32       d1, q11, #16             \n"
    // op[4] += (d1!=0)
    "vsub.s16        d1, d1, d4               \n"
    // op[12]= (d1*2217 - c1*5352 + 51000)>>16
    "vshrn.s32       d3, q12, #16             \n"

    // set result to out array
    "vst1.16         {q0, q1}, [%[out]]   \n"
    : [src_ptr] "+r"(src_ptr), [ref_ptr] "+r"(ref_ptr),
      [coeff32] "+r"(coeff32)          // modified registers
    : [kBPS] "r"(kBPS), [coeff16] "r"(coeff16),
      [out] "r"(out)                   // constants
    : "memory", "q0", "q1", "q2", "q3", "q4", "q5", "q6", "q7", "q8", "q9",
      "q10", "q11", "q12", "q13"       // clobbered
  );
}

static void FTransformWHT(const int16_t* in, int16_t* out) {
  const int kStep = 32;
  __asm__ volatile (
    // d0 = in[0 * 16] , d1 = in[1 * 16]
    // d2 = in[2 * 16] , d3 = in[3 * 16]
    "vld1.16         d0[0], [%[in]], %[kStep]   \n"
    "vld1.16         d1[0], [%[in]], %[kStep]   \n"
    "vld1.16         d2[0], [%[in]], %[kStep]   \n"
    "vld1.16         d3[0], [%[in]], %[kStep]   \n"
    "vld1.16         d0[1], [%[in]], %[kStep]   \n"
    "vld1.16         d1[1], [%[in]], %[kStep]   \n"
    "vld1.16         d2[1], [%[in]], %[kStep]   \n"
    "vld1.16         d3[1], [%[in]], %[kStep]   \n"
    "vld1.16         d0[2], [%[in]], %[kStep]   \n"
    "vld1.16         d1[2], [%[in]], %[kStep]   \n"
    "vld1.16         d2[2], [%[in]], %[kStep]   \n"
    "vld1.16         d3[2], [%[in]], %[kStep]   \n"
    "vld1.16         d0[3], [%[in]], %[kStep]   \n"
    "vld1.16         d1[3], [%[in]], %[kStep]   \n"
    "vld1.16         d2[3], [%[in]], %[kStep]   \n"
    "vld1.16         d3[3], [%[in]], %[kStep]   \n"

    "vaddl.s16       q2, d0, d2                 \n"
    "vshl.s32        q2, q2, #2                 \n" // a0=(in[0*16]+in[2*16])<<2
    "vaddl.s16       q3, d1, d3                 \n"
    "vshl.s32        q3, q3, #2                 \n" // a1=(in[1*16]+in[3*16])<<2
    "vsubl.s16       q4, d1, d3                 \n"
    "vshl.s32        q4, q4, #2                 \n" // a2=(in[1*16]-in[3*16])<<2
    "vsubl.s16       q5, d0, d2                 \n"
    "vshl.s32        q5, q5, #2                 \n" // a3=(in[0*16]-in[2*16])<<2

    "vceq.s32        q10, q2, #0               \n" // (a0 != 0)
    "vadd.s32        q6, q2, q3                \n" // (a0 + a1)
    "vsub.s32        q6, q6, q10               \n" // (a0 + a1) + (a0 != 0)
    "vadd.s32        q7, q5, q4                \n" // a3 + a2
    "vsub.s32        q8, q5, q4                \n" // a3 - a2
    "vsub.s32        q9, q2, q3                \n" // a0 - a1

    // Transponse
    // q6 = tmp[0, 1,  2,  3] ; q7 = tmp[ 4,  5,  6,  7]
    // q8 = tmp[8, 9, 10, 11] ; q9 = tmp[12, 13, 14, 15]
    "vtrn.32         q6, q7                     \n"
    "vtrn.32         q8, q9                     \n"
    "vswp            d16, d13                   \n"
    "vswp            d18, d15                   \n"

    "vqadd.s32        q0, q6, q8                \n" // a0 = tmp[0] + tmp[8]
    "vqadd.s32        q1, q7, q9                \n" // a1 = tmp[4] + tmp[12]
    "vqsub.s32        q2, q7, q9                \n" // a2 = tmp[4] - tmp[12]
    "vqsub.s32        q3, q6, q8                \n" // a3 = tmp[0] - tmp[8]

    "vqadd.s32        q4, q0, q1                \n" // b0 = a0 + a1
    "vqadd.s32        q5, q3, q2                \n" // b1 = a3 + a2
    "vqsub.s32        q6, q3, q2                \n" // b2 = a3 - a2
    "vqsub.s32        q7, q0, q1                \n" // b3 = a0 - a1

    "mov             r10, #3                    \n"
    "vdup.32         q0, r10                    \n" // q0 = 3

    "vceq.s32        q1, q4, #0                 \n" // (b0 != 0)
    "vqsub.s32        q2, q4, q1                \n" // (b0 + (b0 != 0))
    "vqadd.s32        q3, q2, q0                \n" // (b0 + (b0 != 0) + 3)
    "vrshrn.s32       d18, q3, #3               \n" // (b0 + (b0 != 0) + 3) >> 3

    "vceq.s32        q1, q5, #0                 \n" // (b1 != 0)
    "vqsub.s32        q2, q5, q1                \n" // (b1 + (b1 != 0))
    "vqadd.s32        q3, q2, q0                \n" // (b1 + (b1 != 0) + 3)
    "vrshrn.s32       d19, q3, #3               \n" // (b1 + (b1 != 0) + 3) >> 3

    "vceq.s32        q1, q6, #0                 \n" // (b2 != 0)
    "vqsub.s32        q2, q6, q1                \n" // (b2 + (b2 != 0))
    "vqadd.s32        q3, q2, q0                \n" // (b2 + (b2 != 0) + 3)
    "vrshrn.s32       d20, q3, #3               \n" // (b2 + (b2 != 0) + 3) >> 3

    "vceq.s32        q1, q7, #0                 \n" // (b3 != 0)
    "vqsub.s32        q2, q7, q1                \n" // (b3 + (b3 != 0))
    "vqadd.s32        q3, q2, q0                \n" // (b3 + (b3 != 0) + 3)
    "vrshrn.s32       d21, q3, #3               \n" // (b3 + (b3 != 0) + 3) >> 3

    "vst1.16         {q9, q10}, [%[out]]        \n"

    : [in] "+r"(in), [out] "+r"(out)
    : [kStep] "r"(kStep)
    : "memory", "r10", "q0", "q1", "q2", "q3", "q4", "q5",
      "q6", "q7", "q8", "q9", "q10"       // clobbered
  ) ;
}

#endif   // WEBP_USE_NEON

//------------------------------------------------------------------------------
// Entry point

extern void VP8EncDspInitNEON(void);

void VP8EncDspInitNEON(void) {
#if defined(WEBP_USE_NEON)
  VP8ITransform = ITransform;
  VP8FTransform = FTransform;

  VP8ITransformWHT = ITransformWHT;
  VP8FTransformWHT = FTransformWHT;
#endif   // WEBP_USE_NEON
}

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif
