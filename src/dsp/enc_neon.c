// Copyright 2012 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
// ARM NEON version of speed-critical encoding functions.
//
// adapted from libvpx (http://www.webmproject.org/code/)

#include "./dsp.h"

#if defined(WEBP_USE_NEON)

#define USE_INTRINSICS   // use intrinsics when possible

#include <arm_neon.h>

#include "../enc/vp8enci.h"

//------------------------------------------------------------------------------
// Transforms (Paragraph 14.4)

// Inverse transform.
// This code is pretty much the same as TransformOne in the dec_neon.c, except
// for subtraction to *ref. See the comments there for algorithmic explanations.

static const int16_t kC1 = 20091;
static const int16_t kC2 = 17734;  // half of kC2, actually. See comment above.

// This code works but is *slower* than the inlined-asm version below
// (with gcc-4.6). So we disable it for now. Later, it'll be conditional to
// USE_INTRINSICS define.
// With gcc-4.8, it's a little faster speed than inlined-assembly.
#if 0  // defined(USE_INTRINSICS)

// Treats 'v' as an uint8x8_t and zero extends to an int16x8_t.
static WEBP_INLINE int16x8_t ConvertU8ToS16(uint32x2_t v) {
  return vreinterpretq_s16_u16(vmovl_u8(vreinterpret_u8_u32(v)));
}

// Performs unsigned 8b saturation on 'dst01' and 'dst23' storing the result
// to the corresponding rows of 'dst'.
static WEBP_INLINE void SaturateAndStore4x4(uint8_t* const dst,
                                            const int16x8_t dst01,
                                            const int16x8_t dst23) {
  // Unsigned saturate to 8b.
  const uint8x8_t dst01_u8 = vqmovun_s16(dst01);
  const uint8x8_t dst23_u8 = vqmovun_s16(dst23);

  // Store the results.
  vst1_lane_u32((uint32_t*)(dst + 0 * BPS), vreinterpret_u32_u8(dst01_u8), 0);
  vst1_lane_u32((uint32_t*)(dst + 1 * BPS), vreinterpret_u32_u8(dst01_u8), 1);
  vst1_lane_u32((uint32_t*)(dst + 2 * BPS), vreinterpret_u32_u8(dst23_u8), 0);
  vst1_lane_u32((uint32_t*)(dst + 3 * BPS), vreinterpret_u32_u8(dst23_u8), 1);
}

static WEBP_INLINE void Add4x4(const int16x8_t row01, const int16x8_t row23,
                               const uint8_t* const ref, uint8_t* const dst) {
  uint32x2_t dst01 = {0, 0};
  uint32x2_t dst23 = {0, 0};

  // Load the source pixels.
  dst01 = vld1_lane_u32((uint32_t*)(ref + 0 * BPS), dst01, 0);
  dst23 = vld1_lane_u32((uint32_t*)(ref + 2 * BPS), dst23, 0);
  dst01 = vld1_lane_u32((uint32_t*)(ref + 1 * BPS), dst01, 1);
  dst23 = vld1_lane_u32((uint32_t*)(ref + 3 * BPS), dst23, 1);

  {
    // Convert to 16b.
    const int16x8_t dst01_s16 = ConvertU8ToS16(dst01);
    const int16x8_t dst23_s16 = ConvertU8ToS16(dst23);

    // Descale with rounding.
    const int16x8_t out01 = vrsraq_n_s16(dst01_s16, row01, 3);
    const int16x8_t out23 = vrsraq_n_s16(dst23_s16, row23, 3);
    // Add the inverse transform.
    SaturateAndStore4x4(dst, out01, out23);
  }
}

static WEBP_INLINE void Transpose8x2(const int16x8_t in0, const int16x8_t in1,
                                     int16x8x2_t* const out) {
  // a0 a1 a2 a3 | b0 b1 b2 b3   => a0 b0 c0 d0 | a1 b1 c1 d1
  // c0 c1 c2 c3 | d0 d1 d2 d3      a2 b2 c2 d2 | a3 b3 c3 d3
  const int16x8x2_t tmp0 = vzipq_s16(in0, in1);   // a0 c0 a1 c1 a2 c2 ...
                                                  // b0 d0 b1 d1 b2 d2 ...
  *out = vzipq_s16(tmp0.val[0], tmp0.val[1]);
}

static WEBP_INLINE void TransformPass(int16x8x2_t* const rows) {
  // {rows} = in0 | in4
  //          in8 | in12
  // B1 = in4 | in12
  const int16x8_t B1 =
      vcombine_s16(vget_high_s16(rows->val[0]), vget_high_s16(rows->val[1]));
  // C0 = kC1 * in4 | kC1 * in12
  // C1 = kC2 * in4 | kC2 * in12
  const int16x8_t C0 = vsraq_n_s16(B1, vqdmulhq_n_s16(B1, kC1), 1);
  const int16x8_t C1 = vqdmulhq_n_s16(B1, kC2);
  const int16x4_t a = vqadd_s16(vget_low_s16(rows->val[0]),
                                vget_low_s16(rows->val[1]));   // in0 + in8
  const int16x4_t b = vqsub_s16(vget_low_s16(rows->val[0]),
                                vget_low_s16(rows->val[1]));   // in0 - in8
  // c = kC2 * in4 - kC1 * in12
  // d = kC1 * in4 + kC2 * in12
  const int16x4_t c = vqsub_s16(vget_low_s16(C1), vget_high_s16(C0));
  const int16x4_t d = vqadd_s16(vget_low_s16(C0), vget_high_s16(C1));
  const int16x8_t D0 = vcombine_s16(a, b);      // D0 = a | b
  const int16x8_t D1 = vcombine_s16(d, c);      // D1 = d | c
  const int16x8_t E0 = vqaddq_s16(D0, D1);      // a+d | b+c
  const int16x8_t E_tmp = vqsubq_s16(D0, D1);   // a-d | b-c
  const int16x8_t E1 = vcombine_s16(vget_high_s16(E_tmp), vget_low_s16(E_tmp));
  Transpose8x2(E0, E1, rows);
}

static void ITransformOne(const uint8_t* ref,
                          const int16_t* in, uint8_t* dst) {
  int16x8x2_t rows = {{ vld1q_s16(in + 0), vld1q_s16(in + 8) }};
  TransformPass(&rows);
  TransformPass(&rows);
  Add4x4(rows.val[0], rows.val[1], ref, dst);
}

#else

static void ITransformOne(const uint8_t* ref,
                          const int16_t* in, uint8_t* dst) {
  const int kBPS = BPS;
  const int16_t kC1C2[] = { kC1, kC2, 0, 0 };

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

#endif    // USE_INTRINSICS

static void ITransform(const uint8_t* ref,
                       const int16_t* in, uint8_t* dst, int do_two) {
  ITransformOne(ref, in, dst);
  if (do_two) {
    ITransformOne(ref + 4, in + 16, dst + 4);
  }
}

// Forward transform.

// adapted from vp8/encoder/arm/neon/shortfdct_neon.asm
static const int16_t kCoeff16[] = {
  5352,  5352,  5352, 5352, 2217,  2217,  2217, 2217
};
static const int32_t kCoeff32[] = {
   1812,  1812,  1812,  1812,
    937,   937,   937,   937,
  12000, 12000, 12000, 12000,
  51000, 51000, 51000, 51000
};

static void FTransform(const uint8_t* src, const uint8_t* ref,
                       int16_t* out) {
  const int kBPS = BPS;
  const uint8_t* src_ptr = src;
  const uint8_t* ref_ptr = ref;
  const int16_t* coeff16 = kCoeff16;
  const int32_t* coeff32 = kCoeff32;

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

    // d[0-3] = src - ref
    "vsubl.u8    q0, d8, d12                  \n"
    "vsubl.u8    q1, d9, d13                  \n"

    // load coeff16 into q8(d16=5352, d17=2217)
    "vld1.16     {q8}, [%[coeff16]]           \n"

    // load coeff32 high half into q9 = 1812, q10 = 937
    "vld1.32     {q9, q10}, [%[coeff32]]!     \n"

    // load coeff32 low half into q11=12000, q12=51000
    "vld1.32     {q11,q12}, [%[coeff32]]      \n"

    // part 1
    // Transpose. Register dN is the same as dN in C
    "vtrn.32         d0, d2                   \n"
    "vtrn.32         d1, d3                   \n"
    "vtrn.16         d0, d1                   \n"
    "vtrn.16         d2, d3                   \n"

    "vadd.s16        d4, d0, d3               \n" // a0 = d0 + d3
    "vadd.s16        d5, d1, d2               \n" // a1 = d1 + d2
    "vsub.s16        d6, d1, d2               \n" // a2 = d1 - d2
    "vsub.s16        d7, d0, d3               \n" // a3 = d0 - d3

    "vadd.s16        d0, d4, d5               \n" // a0 + a1
    "vshl.s16        d0, d0, #3               \n" // temp[0+i*4] = (a0+a1) << 3
    "vsub.s16        d2, d4, d5               \n" // a0 - a1
    "vshl.s16        d2, d2, #3               \n" // (temp[2+i*4] = (a0-a1) << 3

    "vmlal.s16       q9, d7, d16              \n" // a3*5352 + 1812
    "vmlal.s16       q10, d7, d17             \n" // a3*2217 + 937
    "vmlal.s16       q9, d6, d17              \n" // a2*2217 + a3*5352 + 1812
    "vmlsl.s16       q10, d6, d16             \n" // a3*2217 + 937 - a2*5352

    // temp[1+i*4] = (d2*2217 + d3*5352 + 1812) >> 9
    // temp[3+i*4] = (d3*2217 + 937 - d2*5352) >> 9
    "vshrn.s32       d1, q9, #9               \n"
    "vshrn.s32       d3, q10, #9              \n"

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

    "vmvn            d4, d4                   \n" // !(d1 == 0)
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

    "vaddl.s16       q2, d0, d2                 \n" // a0=(in[0*16]+in[2*16])
    "vaddl.s16       q3, d1, d3                 \n" // a1=(in[1*16]+in[3*16])
    "vsubl.s16       q4, d1, d3                 \n" // a2=(in[1*16]-in[3*16])
    "vsubl.s16       q5, d0, d2                 \n" // a3=(in[0*16]-in[2*16])

    "vqadd.s32       q6, q2, q3                 \n" // a0 + a1
    "vqadd.s32       q7, q5, q4                 \n" // a3 + a2
    "vqsub.s32       q8, q5, q4                 \n" // a3 - a2
    "vqsub.s32       q9, q2, q3                 \n" // a0 - a1

    // Transpose
    // q6 = tmp[0, 1,  2,  3] ; q7 = tmp[ 4,  5,  6,  7]
    // q8 = tmp[8, 9, 10, 11] ; q9 = tmp[12, 13, 14, 15]
    "vswp            d13, d16                   \n" // vtrn.64 q0, q2
    "vswp            d15, d18                   \n" // vtrn.64 q1, q3
    "vtrn.32         q6, q7                     \n"
    "vtrn.32         q8, q9                     \n"

    "vqadd.s32       q0, q6, q8                 \n" // a0 = tmp[0] + tmp[8]
    "vqadd.s32       q1, q7, q9                 \n" // a1 = tmp[4] + tmp[12]
    "vqsub.s32       q2, q7, q9                 \n" // a2 = tmp[4] - tmp[12]
    "vqsub.s32       q3, q6, q8                 \n" // a3 = tmp[0] - tmp[8]

    "vqadd.s32       q4, q0, q1                 \n" // b0 = a0 + a1
    "vqadd.s32       q5, q3, q2                 \n" // b1 = a3 + a2
    "vqsub.s32       q6, q3, q2                 \n" // b2 = a3 - a2
    "vqsub.s32       q7, q0, q1                 \n" // b3 = a0 - a1

    "vshrn.s32       d18, q4, #1                \n" // b0 >> 1
    "vshrn.s32       d19, q5, #1                \n" // b1 >> 1
    "vshrn.s32       d20, q6, #1                \n" // b2 >> 1
    "vshrn.s32       d21, q7, #1                \n" // b3 >> 1

    "vst1.16         {q9, q10}, [%[out]]        \n"

    : [in] "+r"(in)
    : [kStep] "r"(kStep), [out] "r"(out)
    : "memory", "q0", "q1", "q2", "q3", "q4", "q5",
      "q6", "q7", "q8", "q9", "q10"       // clobbered
  ) ;
}

//------------------------------------------------------------------------------
// Texture distortion
//
// We try to match the spectral content (weighted) between source and
// reconstructed samples.

// Hadamard transform
// Returns the weighted sum of the absolute value of transformed coefficients.
// This uses a TTransform helper function in C
static int Disto4x4(const uint8_t* const a, const uint8_t* const b,
                    const uint16_t* const w) {
  const int kBPS = BPS;
  const uint8_t* A = a;
  const uint8_t* B = b;
  const uint16_t* W = w;
  int sum;
  __asm__ volatile (
    "vld1.32         d0[0], [%[a]], %[kBPS]   \n"
    "vld1.32         d0[1], [%[a]], %[kBPS]   \n"
    "vld1.32         d2[0], [%[a]], %[kBPS]   \n"
    "vld1.32         d2[1], [%[a]]            \n"

    "vld1.32         d1[0], [%[b]], %[kBPS]   \n"
    "vld1.32         d1[1], [%[b]], %[kBPS]   \n"
    "vld1.32         d3[0], [%[b]], %[kBPS]   \n"
    "vld1.32         d3[1], [%[b]]            \n"

    // a d0/d2, b d1/d3
    // d0/d1: 01 01 01 01
    // d2/d3: 23 23 23 23
    // But: it goes 01 45 23 67
    // Notice the middle values are transposed
    "vtrn.16         q0, q1                   \n"

    // {a0, a1} = {in[0] + in[2], in[1] + in[3]}
    "vaddl.u8        q2, d0, d2               \n"
    "vaddl.u8        q10, d1, d3              \n"
    // {a3, a2} = {in[0] - in[2], in[1] - in[3]}
    "vsubl.u8        q3, d0, d2               \n"
    "vsubl.u8        q11, d1, d3              \n"

    // tmp[0] = a0 + a1
    "vpaddl.s16      q0, q2                   \n"
    "vpaddl.s16      q8, q10                  \n"

    // tmp[1] = a3 + a2
    "vpaddl.s16      q1, q3                   \n"
    "vpaddl.s16      q9, q11                  \n"

    // No pair subtract
    // q2 = {a0, a3}
    // q3 = {a1, a2}
    "vtrn.16         q2, q3                   \n"
    "vtrn.16         q10, q11                 \n"

    // {tmp[3], tmp[2]} = {a0 - a1, a3 - a2}
    "vsubl.s16       q12, d4, d6              \n"
    "vsubl.s16       q13, d5, d7              \n"
    "vsubl.s16       q14, d20, d22            \n"
    "vsubl.s16       q15, d21, d23            \n"

    // separate tmp[3] and tmp[2]
    // q12 = tmp[3]
    // q13 = tmp[2]
    "vtrn.32         q12, q13                 \n"
    "vtrn.32         q14, q15                 \n"

    // Transpose tmp for a
    "vswp            d1, d26                  \n" // vtrn.64
    "vswp            d3, d24                  \n" // vtrn.64
    "vtrn.32         q0, q1                   \n"
    "vtrn.32         q13, q12                 \n"

    // Transpose tmp for b
    "vswp            d17, d30                 \n" // vtrn.64
    "vswp            d19, d28                 \n" // vtrn.64
    "vtrn.32         q8, q9                   \n"
    "vtrn.32         q15, q14                 \n"

    // The first Q register is a, the second b.
    // q0/8 tmp[0-3]
    // q13/15 tmp[4-7]
    // q1/9 tmp[8-11]
    // q12/14 tmp[12-15]

    // These are still in 01 45 23 67 order. We fix it easily in the addition
    // case but the subtraction propagates them.
    "vswp            d3, d27                  \n"
    "vswp            d19, d31                 \n"

    // a0 = tmp[0] + tmp[8]
    "vadd.s32        q2, q0, q1               \n"
    "vadd.s32        q3, q8, q9               \n"

    // a1 = tmp[4] + tmp[12]
    "vadd.s32        q10, q13, q12            \n"
    "vadd.s32        q11, q15, q14            \n"

    // a2 = tmp[4] - tmp[12]
    "vsub.s32        q13, q13, q12            \n"
    "vsub.s32        q15, q15, q14            \n"

    // a3 = tmp[0] - tmp[8]
    "vsub.s32        q0, q0, q1               \n"
    "vsub.s32        q8, q8, q9               \n"

    // b0 = a0 + a1
    "vadd.s32        q1, q2, q10              \n"
    "vadd.s32        q9, q3, q11              \n"

    // b1 = a3 + a2
    "vadd.s32        q12, q0, q13             \n"
    "vadd.s32        q14, q8, q15             \n"

    // b2 = a3 - a2
    "vsub.s32        q0, q0, q13              \n"
    "vsub.s32        q8, q8, q15              \n"

    // b3 = a0 - a1
    "vsub.s32        q2, q2, q10              \n"
    "vsub.s32        q3, q3, q11              \n"

    "vld1.64         {q10, q11}, [%[w]]       \n"

    // abs(b0)
    "vabs.s32        q1, q1                   \n"
    "vabs.s32        q9, q9                   \n"
    // abs(b1)
    "vabs.s32        q12, q12                 \n"
    "vabs.s32        q14, q14                 \n"
    // abs(b2)
    "vabs.s32        q0, q0                   \n"
    "vabs.s32        q8, q8                   \n"
    // abs(b3)
    "vabs.s32        q2, q2                   \n"
    "vabs.s32        q3, q3                   \n"

    // expand w before using.
    "vmovl.u16       q13, d20                 \n"
    "vmovl.u16       q15, d21                 \n"

    // w[0] * abs(b0)
    "vmul.u32        q1, q1, q13              \n"
    "vmul.u32        q9, q9, q13              \n"

    // w[4] * abs(b1)
    "vmla.u32        q1, q12, q15             \n"
    "vmla.u32        q9, q14, q15             \n"

    // expand w before using.
    "vmovl.u16       q13, d22                 \n"
    "vmovl.u16       q15, d23                 \n"

    // w[8] * abs(b1)
    "vmla.u32        q1, q0, q13              \n"
    "vmla.u32        q9, q8, q13              \n"

    // w[12] * abs(b1)
    "vmla.u32        q1, q2, q15              \n"
    "vmla.u32        q9, q3, q15              \n"

    // Sum the arrays
    "vpaddl.u32      q1, q1                   \n"
    "vpaddl.u32      q9, q9                   \n"
    "vadd.u64        d2, d3                   \n"
    "vadd.u64        d18, d19                 \n"

    // Hadamard transform needs 4 bits of extra precision (2 bits in each
    // direction) for dynamic raw. Weights w[] are 16bits at max, so the maximum
    // precision for coeff is 8bit of input + 4bits of Hadamard transform +
    // 16bits for w[] + 2 bits of abs() summation.
    //
    // This uses a maximum of 31 bits (signed). Discarding the top 32 bits is
    // A-OK.

    // sum2 - sum1
    "vsub.u32        d0, d2, d18              \n"
    // abs(sum2 - sum1)
    "vabs.s32        d0, d0                   \n"
    // abs(sum2 - sum1) >> 5
    "vshr.u32        d0, #5                   \n"

    // It would be better to move the value straight into r0 but I'm not
    // entirely sure how this works with inline assembly.
    "vmov.32         %[sum], d0[0]            \n"

    : [sum] "=r"(sum), [a] "+r"(A), [b] "+r"(B), [w] "+r"(W)
    : [kBPS] "r"(kBPS)
    : "memory", "q0", "q1", "q2", "q3", "q4", "q5", "q6", "q7", "q8", "q9",
      "q10", "q11", "q12", "q13", "q14", "q15"  // clobbered
  ) ;

  return sum;
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

static WEBP_INLINE void AccumulateSSE16(const uint8_t* const a,
                                        const uint8_t* const b,
                                        uint32x4_t* const sum ) {
  const uint8x16_t a0 = vld1q_u8(a);
  const uint8x16_t b0 = vld1q_u8(b);
  const uint8x16_t abs_diff = vabdq_u8(a0, b0);
  const uint16x8_t prod_l = vmull_u8(vget_low_u8(abs_diff),
                                     vget_low_u8(abs_diff));
  const uint16x8_t prod_h = vmull_u8(vget_high_u8(abs_diff),
                                     vget_high_u8(abs_diff));
  const uint16x8_t prod = vaddq_u16(prod_h, prod_l);
  *sum = vpadalq_u16(*sum, prod);      // pair-wise multiply and accumulate
}

// Horizontal sum of all four uint32_t values in 'sum'.
static int SumToInt(uint32x4_t sum) {
  const uint64x2_t sum2 = vpaddlq_u32(sum);
  const uint64_t sum3 = vgetq_lane_u64(sum2, 0) + vgetq_lane_u64(sum2, 1);
  return (int)sum3;
}

static int SSE16x16(const uint8_t* a, const uint8_t* b) {
  uint32x4_t sum = { 0, 0, 0, 0 };
  int y;
  for (y = 0; y < 16; ++y) {
    AccumulateSSE16(a + y * BPS, b + y * BPS, &sum);
  }
  return SumToInt(sum);
}

static int SSE16x8(const uint8_t* a, const uint8_t* b) {
  uint32x4_t sum = { 0, 0, 0, 0 };
  int y;
  for (y = 0; y < 8; ++y) {
    AccumulateSSE16(a + y * BPS, b + y * BPS, &sum);
  }
  return SumToInt(sum);
}

static int SSE8x8(const uint8_t* a, const uint8_t* b) {
  uint32x4_t sum = { 0, 0, 0, 0 };
  int y;
  for (y = 0; y < 8; ++y) {
    const uint8x8_t a0 = vld1_u8(a + y * BPS);
    const uint8x8_t b0 = vld1_u8(b + y * BPS);
    const uint8x8_t abs_diff = vabd_u8(a0, b0);
    const uint16x8_t prod = vmull_u8(abs_diff, abs_diff);
    sum = vpadalq_u16(sum, prod);
  }
  return SumToInt(sum);
}

#define LOAD_LANE_32b(src, VALUE, LANE) \
    (VALUE) = vld1q_lane_u32((const uint32_t*)(src), (VALUE), (LANE))

static int SSE4x4(const uint8_t* a, const uint8_t* b) {
  uint32x4_t a0 = { 0, 0, 0, 0 };
  uint32x4_t b0 = { 0, 0, 0, 0 };
  // Load all 4x4 pixels into a single uin32x4_t variable.
  LOAD_LANE_32b(a + 0 * BPS, a0, 0);
  LOAD_LANE_32b(a + 1 * BPS, a0, 1);
  LOAD_LANE_32b(a + 2 * BPS, a0, 2);
  LOAD_LANE_32b(a + 3 * BPS, a0, 3);
  LOAD_LANE_32b(b + 0 * BPS, b0, 0);
  LOAD_LANE_32b(b + 1 * BPS, b0, 1);
  LOAD_LANE_32b(b + 2 * BPS, b0, 2);
  LOAD_LANE_32b(b + 3 * BPS, b0, 3);
  {
    const uint8x16_t abs_diff = vabdq_u8(vreinterpretq_u8_u32(a0),
                                         vreinterpretq_u8_u32(b0));
    const uint16x8_t prod_l = vmull_u8(vget_low_u8(abs_diff),
                                       vget_low_u8(abs_diff));
    const uint16x8_t prod_h = vmull_u8(vget_high_u8(abs_diff),
                                       vget_high_u8(abs_diff));
    const uint32x4_t sum = vpaddlq_u16(vaddq_u16(prod_h, prod_l));
    return SumToInt(sum);
  }
}

#undef LOAD_LANE_32b

#endif   // WEBP_USE_NEON

//------------------------------------------------------------------------------
// Entry point

extern void VP8EncDspInitNEON(void);

void VP8EncDspInitNEON(void) {
#if defined(WEBP_USE_NEON)
  VP8ITransform = ITransform;
  VP8FTransform = FTransform;

  VP8FTransformWHT = FTransformWHT;

  VP8TDisto4x4 = Disto4x4;
  VP8TDisto16x16 = Disto16x16;
#if defined(USE_INTRINSICS)
  VP8SSE16x16 = SSE16x16;
  VP8SSE16x8 = SSE16x8;
  VP8SSE8x8 = SSE8x8;
  VP8SSE4x4 = SSE4x4;
#endif
#endif   // WEBP_USE_NEON
}
