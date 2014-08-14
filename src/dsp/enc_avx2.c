// Copyright 2014 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
// AVX2 version of speed-critical encoding functions.

#include "./dsp.h"

#if defined(WEBP_USE_AVX2)
#include "../enc/vp8enci.h"
#include <immintrin.h>

#define BROADCAST64(mem) \
  _mm256_castpd_si256(_mm256_broadcast_sd((double const *)mem));

#define BROADCAST128(mem) \
  _mm256_castpd_si256(_mm256_broadcast_pd((__m128d const *)mem));


WEBP_ALIGNED_DECL(32, static const uint8_t, kFilter4[32]) = {
0, 128, 1, 128, 2, 128, 3, 128, 128, 128, 128, 128, 128, 128, 128, 128,
4, 128, 5, 128, 6, 128, 7, 128, 128, 128, 128, 128, 128, 128, 128, 128};

WEBP_ALIGNED_DECL(32, static const uint8_t, kFilter8[32]) = {
0, 128, 1, 128, 2, 128, 3, 128, 4, 128, 5, 128, 6, 128, 7, 128,
8, 128, 9, 128, 10, 128, 11, 128, 12, 128, 13, 128, 14, 128, 15, 128};

static void SSE4x4x2(const uint8_t* a, const uint8_t* b, int64_t *D1, int64_t *D2) {
  const __m256i zero = _mm256_setzero_si256();
  __m256i a0, a1, a2, a3;
  __m256i b0, b1, b2, b3;
  __m256i sum2;
  int32_t tmp[8];

  // Load values. Note that we read 8 pixels instead of 4,
  // but the a/b buffers are over-allocated to that effect.
  a0 = BROADCAST64(&a[BPS * 0]);
  a1 = BROADCAST64(&a[BPS * 1]);
  a2 = BROADCAST64(&a[BPS * 2]);
  a3 = BROADCAST64(&a[BPS * 3]);
  b0 = _mm256_load_si256((__m256i const *)(b + 0));
  b1 = _mm256_load_si256((__m256i const *)(b + 32));
  b2 = _mm256_load_si256((__m256i const *)(b + 32*2));
  b3 = _mm256_load_si256((__m256i const *)(b + 32*3));

  {
    // Combine pair of lines and convert to 16b.
    const __m256i a01 = _mm256_unpacklo_epi32(a0, a1);
    const __m256i a23 = _mm256_unpacklo_epi32(a2, a3);
    const __m256i b01 = _mm256_unpacklo_epi32(b0, b1);
    const __m256i b23 = _mm256_unpacklo_epi32(b2, b3);
    const __m256i a01s = _mm256_unpacklo_epi8(a01, zero);
    const __m256i a23s = _mm256_unpacklo_epi8(a23, zero);
    const __m256i b01s = _mm256_unpacklo_epi8(b01, zero);
    const __m256i b23s = _mm256_unpacklo_epi8(b23, zero);

    // Compute differences; (a-b)^2 = (abs(a-b))^2 = (sat8(a-b) + sat8(b-a))^2
    const __m256i d0 = _mm256_subs_epu8(a01s, b01s);
    const __m256i d1 = _mm256_subs_epu8(b01s, a01s);
    const __m256i d2 = _mm256_subs_epu8(a23s, b23s);
    const __m256i d3 = _mm256_subs_epu8(b23s, a23s);

    // Square and add them all together.
    const __m256i madd0 = _mm256_madd_epi16(d0, d0);
    const __m256i madd1 = _mm256_madd_epi16(d1, d1);
    const __m256i madd2 = _mm256_madd_epi16(d2, d2);
    const __m256i madd3 = _mm256_madd_epi16(d3, d3);
    const __m256i sum0 = _mm256_add_epi32(madd0, madd1);
    const __m256i sum1 = _mm256_add_epi32(madd2, madd3);
    sum2 = _mm256_add_epi32(sum0, sum1);
  }
  _mm256_storeu_si256((__m256i*)tmp, sum2);
  *D1 = tmp[3] + tmp[2] + tmp[1] + tmp[0];
  *D2 = tmp[7] + tmp[6] + tmp[5] + tmp[4];
}

//------------------------------------------------------------------------------
// Transforms (Paragraph 14.4)
// Does two inverse transforms.
static void ITransform4x2(const uint8_t* ref, const int16_t* in, uint8_t* dst,

                           int do_two) {
  // This implementation makes use of 16-bit fixed point versions of two
  // multiply constants:
  //    K1 = sqrt(2) * cos (pi/8) ~= 85627 / 2^16
  //    K2 = sqrt(2) * sin (pi/8) ~= 35468 / 2^16
  //
  // To be able to use signed 16-bit integers, we use the following trick to
  // have constants within range:
  // - Associated constants are obtained by subtracting the 16-bit fixed point
  //   version of one:
  //      k = K - (1 << 16)  =>  K = k + (1 << 16)
  //      K1 = 85267  =>  k1 =  20091
  //      K2 = 35468  =>  k2 = -30068
  // - The multiplication of a variable by a constant become the sum of the
  //   variable and the multiplication of that variable by the associated
  //   constant:
  //      (x * K) >> 16 = (x * (k + (1 << 16))) >> 16 = ((x * k ) >> 16) + x
  const __m256i k1 = _mm256_set1_epi16(20091);
  const __m256i k2 = _mm256_set1_epi16(-30068);
  __m256i T0, T1, T2, T3;

  // Load and concatenate the transform coefficients (we'll do two inverse
  // transforms in parallel). In the case of only one inverse transform, the
  // second half of the vectors will just contain random value we'll never
  // use nor store.
  __m256i in0, in1, in2, in3;
  {
    in0 = _mm256_loadu_si256((__m256i*)&in[0]);
    in1 = _mm256_loadu_si256((__m256i*)&in[4]);
    in2 = _mm256_loadu_si256((__m256i*)&in[16]);
    in3 = _mm256_loadu_si256((__m256i*)&in[20]);
    // a00 a10 a20 a30   x x x x
    // a01 a11 a21 a31   x x x x
    // a02 a12 a22 a32   x x x x
    // a03 a13 a23 a33   x x x x
    if (do_two) {
      const __m256i inB0 = _mm256_loadu_si256((__m256i*)&in[24]);
      const __m256i inB1 = _mm256_loadu_si256((__m256i*)&in[28]);
      const __m256i inB2 = _mm256_loadu_si256((__m256i*)&in[32]);
      const __m256i inB3 = _mm256_loadu_si256((__m256i*)&in[36]);
      in0 = _mm256_unpacklo_epi64(in0, inB0);
      in1 = _mm256_unpacklo_epi64(in1, inB1);
      in2 = _mm256_unpacklo_epi64(in2, inB2);
      in3 = _mm256_unpacklo_epi64(in3, inB3);
      // a00 a10 a20 a30   b00 b10 b20 b30
      // a01 a11 a21 a31   b01 b11 b21 b31
      // a02 a12 a22 a32   b02 b12 b22 b32
      // a03 a13 a23 a33   b03 b13 b23 b33
    }
  }

  // Vertical pass and subsequent transpose.
  {
    // First pass, c and d calculations are longer because of the "trick"
    // multiplications.
    const __m256i a = _mm256_add_epi16(in0, in2);
    const __m256i b = _mm256_sub_epi16(in0, in2);
    // c = MUL(in1, K2) - MUL(in3, K1) = MUL(in1, k2) - MUL(in3, k1) + in1 - in3
    const __m256i c1 = _mm256_mulhi_epi16(in1, k2);
    const __m256i c2 = _mm256_mulhi_epi16(in3, k1);
    const __m256i c3 = _mm256_sub_epi16(in1, in3);
    const __m256i c4 = _mm256_sub_epi16(c1, c2);
    const __m256i c = _mm256_add_epi16(c3, c4);
    // d = MUL(in1, K1) + MUL(in3, K2) = MUL(in1, k1) + MUL(in3, k2) + in1 + in3
    const __m256i d1 = _mm256_mulhi_epi16(in1, k1);
    const __m256i d2 = _mm256_mulhi_epi16(in3, k2);
    const __m256i d3 = _mm256_add_epi16(in1, in3);
    const __m256i d4 = _mm256_add_epi16(d1, d2);
    const __m256i d = _mm256_add_epi16(d3, d4);

    // Second pass.
    const __m256i tmp0 = _mm256_add_epi16(a, d);
    const __m256i tmp1 = _mm256_add_epi16(b, c);
    const __m256i tmp2 = _mm256_sub_epi16(b, c);
    const __m256i tmp3 = _mm256_sub_epi16(a, d);

    // Transpose the two 4x4.
    // a00 a01 a02 a03   b00 b01 b02 b03
    // a10 a11 a12 a13   b10 b11 b12 b13
    // a20 a21 a22 a23   b20 b21 b22 b23
    // a30 a31 a32 a33   b30 b31 b32 b33
    const __m256i transpose0_0 = _mm256_unpacklo_epi16(tmp0, tmp1);
    const __m256i transpose0_1 = _mm256_unpacklo_epi16(tmp2, tmp3);
    const __m256i transpose0_2 = _mm256_unpackhi_epi16(tmp0, tmp1);
    const __m256i transpose0_3 = _mm256_unpackhi_epi16(tmp2, tmp3);
    // a00 a10 a01 a11   a02 a12 a03 a13
    // a20 a30 a21 a31   a22 a32 a23 a33
    // b00 b10 b01 b11   b02 b12 b03 b13
    // b20 b30 b21 b31   b22 b32 b23 b33
    const __m256i transpose1_0 = _mm256_unpacklo_epi32(transpose0_0, transpose0_1);
    const __m256i transpose1_1 = _mm256_unpacklo_epi32(transpose0_2, transpose0_3);
    const __m256i transpose1_2 = _mm256_unpackhi_epi32(transpose0_0, transpose0_1);
    const __m256i transpose1_3 = _mm256_unpackhi_epi32(transpose0_2, transpose0_3);
    // a00 a10 a20 a30 a01 a11 a21 a31
    // b00 b10 b20 b30 b01 b11 b21 b31
    // a02 a12 a22 a32 a03 a13 a23 a33
    // b02 b12 a22 b32 b03 b13 b23 b33
    T0 = _mm256_unpacklo_epi64(transpose1_0, transpose1_1);
    T1 = _mm256_unpackhi_epi64(transpose1_0, transpose1_1);
    T2 = _mm256_unpacklo_epi64(transpose1_2, transpose1_3);
    T3 = _mm256_unpackhi_epi64(transpose1_2, transpose1_3);
    // a00 a10 a20 a30   b00 b10 b20 b30
    // a01 a11 a21 a31   b01 b11 b21 b31
    // a02 a12 a22 a32   b02 b12 b22 b32
    // a03 a13 a23 a33   b03 b13 b23 b33
  }

  // Horizontal pass and subsequent transpose.
  {
    // First pass, c and d calculations are longer because of the "trick"
    // multiplications.
    const __m256i four = _mm256_set1_epi16(4);
    const __m256i dc = _mm256_add_epi16(T0, four);
    const __m256i a =  _mm256_add_epi16(dc, T2);
    const __m256i b =  _mm256_sub_epi16(dc, T2);
    // c = MUL(T1, K2) - MUL(T3, K1) = MUL(T1, k2) - MUL(T3, k1) + T1 - T3
    const __m256i c1 = _mm256_mulhi_epi16(T1, k2);
    const __m256i c2 = _mm256_mulhi_epi16(T3, k1);
    const __m256i c3 = _mm256_sub_epi16(T1, T3);
    const __m256i c4 = _mm256_sub_epi16(c1, c2);
    const __m256i c = _mm256_add_epi16(c3, c4);
    // d = MUL(T1, K1) + MUL(T3, K2) = MUL(T1, k1) + MUL(T3, k2) + T1 + T3
    const __m256i d1 = _mm256_mulhi_epi16(T1, k1);
    const __m256i d2 = _mm256_mulhi_epi16(T3, k2);
    const __m256i d3 = _mm256_add_epi16(T1, T3);
    const __m256i d4 = _mm256_add_epi16(d1, d2);
    const __m256i d = _mm256_add_epi16(d3, d4);

    // Second pass.
    const __m256i tmp0 = _mm256_add_epi16(a, d);
    const __m256i tmp1 = _mm256_add_epi16(b, c);
    const __m256i tmp2 = _mm256_sub_epi16(b, c);
    const __m256i tmp3 = _mm256_sub_epi16(a, d);
    const __m256i shifted0 = _mm256_srai_epi16(tmp0, 3);
    const __m256i shifted1 = _mm256_srai_epi16(tmp1, 3);
    const __m256i shifted2 = _mm256_srai_epi16(tmp2, 3);
    const __m256i shifted3 = _mm256_srai_epi16(tmp3, 3);

    // Transpose the two 4x4.
    // a00 a01 a02 a03   b00 b01 b02 b03
    // a10 a11 a12 a13   b10 b11 b12 b13
    // a20 a21 a22 a23   b20 b21 b22 b23
    // a30 a31 a32 a33   b30 b31 b32 b33
    const __m256i transpose0_0 = _mm256_unpacklo_epi16(shifted0, shifted1);
    const __m256i transpose0_1 = _mm256_unpacklo_epi16(shifted2, shifted3);
    const __m256i transpose0_2 = _mm256_unpackhi_epi16(shifted0, shifted1);
    const __m256i transpose0_3 = _mm256_unpackhi_epi16(shifted2, shifted3);
    // a00 a10 a01 a11   a02 a12 a03 a13
    // a20 a30 a21 a31   a22 a32 a23 a33
    // b00 b10 b01 b11   b02 b12 b03 b13
    // b20 b30 b21 b31   b22 b32 b23 b33
    const __m256i transpose1_0 = _mm256_unpacklo_epi32(transpose0_0, transpose0_1);
    const __m256i transpose1_1 = _mm256_unpacklo_epi32(transpose0_2, transpose0_3);
    const __m256i transpose1_2 = _mm256_unpackhi_epi32(transpose0_0, transpose0_1);
    const __m256i transpose1_3 = _mm256_unpackhi_epi32(transpose0_2, transpose0_3);
    // a00 a10 a20 a30 a01 a11 a21 a31
    // b00 b10 b20 b30 b01 b11 b21 b31
    // a02 a12 a22 a32 a03 a13 a23 a33
    // b02 b12 a22 b32 b03 b13 b23 b33
    T0 = _mm256_unpacklo_epi64(transpose1_0, transpose1_1);
    T1 = _mm256_unpackhi_epi64(transpose1_0, transpose1_1);
    T2 = _mm256_unpacklo_epi64(transpose1_2, transpose1_3);
    T3 = _mm256_unpackhi_epi64(transpose1_2, transpose1_3);
    // a00 a10 a20 a30   b00 b10 b20 b30
    // a01 a11 a21 a31   b01 b11 b21 b31
    // a02 a12 a22 a32   b02 b12 b22 b32
    // a03 a13 a23 a33   b03 b13 b23 b33
  }

  // Add inverse transform to 'ref' and store.
  {
    // Load the reference(s).
    __m256i ref0, ref1, ref2, ref3, flt;

    if (do_two) {
      // Load 16 bytes/pixels per line.
      ref0 = BROADCAST128(&ref[0 * BPS]);
      ref1 = BROADCAST128(&ref[1 * BPS]);
      ref2 = BROADCAST128(&ref[2 * BPS]);
      ref3 = BROADCAST128(&ref[3 * BPS]);
      flt = _mm256_load_si256((__m256i const *)kFilter8);
    } else {
      // Load 8 bytes/pixels per line.
      ref0 = BROADCAST64(&ref[0 * BPS]);
      ref1 = BROADCAST64(&ref[1 * BPS]);
      ref2 = BROADCAST64(&ref[2 * BPS]);
      ref3 = BROADCAST64(&ref[3 * BPS]);
      flt = _mm256_load_si256((__m256i const *)kFilter4);
    }
    // Convert to 16b.
    ref0 = _mm256_shuffle_epi8(ref0, flt);
    ref1 = _mm256_shuffle_epi8(ref1, flt);
    ref2 = _mm256_shuffle_epi8(ref2, flt);
    ref3 = _mm256_shuffle_epi8(ref3, flt);
    // Add the inverse transform(s).
    ref0 = _mm256_add_epi16(ref0, T0);
    ref1 = _mm256_add_epi16(ref1, T1);
    ref2 = _mm256_add_epi16(ref2, T2);
    ref3 = _mm256_add_epi16(ref3, T3);
    // Unsigned saturate to 8b.
    // Store the results.
    _mm256_store_si256((__m256i *)(dst + 32*0), _mm256_packus_epi16(ref0, ref0));
    _mm256_store_si256((__m256i *)(dst + 32*1), _mm256_packus_epi16(ref1, ref1));
    _mm256_store_si256((__m256i *)(dst + 32*2), _mm256_packus_epi16(ref2, ref2));
    _mm256_store_si256((__m256i *)(dst + 32*3), _mm256_packus_epi16(ref3, ref3));
  }
}

static void FTransform4x2(const uint8_t* src, const uint8_t* ref, int16_t* out) {
  const __m256i zero = _mm256_setzero_si256();
  const __m256i seven = _mm256_set1_epi16(7);
  const __m256i k937 = _mm256_set1_epi32(937);
  const __m256i k1812 = _mm256_set1_epi32(1812);
  const __m256i k51000 = _mm256_set1_epi32(51000);
  const __m256i k12000_plus_one = _mm256_set1_epi32(12000 + (1 << 16));
  const __m256i k5352_2217 = _mm256_set_epi16(5352,  2217, 5352,  2217,
                                              5352,  2217, 5352,  2217,
                                              5352,  2217, 5352,  2217,
                                              5352,  2217, 5352,  2217);
  const __m256i k2217_5352 = _mm256_set_epi16(2217, -5352, 2217, -5352,
                                              2217, -5352, 2217, -5352,
                                              2217, -5352, 2217, -5352,
                                              2217, -5352, 2217, -5352);
  const __m256i k88p = _mm256_set_epi16(8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8);
  const __m256i k88m = _mm256_set_epi16(-8, 8, -8, 8, -8, 8, -8, 8, -8, 8, -8, 8, -8, 8, -8, 8);
  const __m256i k5352_2217p = _mm256_set_epi16(2217, 5352, 2217, 5352,
                                               2217, 5352, 2217, 5352,
                                               2217, 5352, 2217, 5352,
                                               2217, 5352, 2217, 5352);
  const __m256i k5352_2217m = _mm256_set_epi16(-5352, 2217, -5352, 2217,
                                               -5352, 2217, -5352, 2217,
                                               -5352, 2217, -5352, 2217,
                                               -5352, 2217, -5352, 2217);
  __m256i v01, v32;

  // Difference between src and ref and initial transpose.
  {
    __m256i src0, src1, src2, src3;
    __m256i ref0, ref1, ref2, ref3;

    // Load src
    // for processing two iterations at a time a 256 bit register is used.
    // the data of the src is being broadcast to the first lane and the second lane
    // of the 256 bit register. A floating point broadcast is being used here because
    // versus to "integer broadcast" this broadcast is done from the memeory
    // and not from a register (1 cycle vs 3 cycles)
    src0 = BROADCAST64(&src[0 * BPS]);
    src1 = BROADCAST64(&src[1 * BPS]);
    src2 = BROADCAST64(&src[2 * BPS]);
    src3 = BROADCAST64(&src[3 * BPS]);

    // Load ref
    // The first 8 element are being broadcast to the first lane and second lane of
    // the 256 bit register.
    ref0 = BROADCAST64(&ref[0 * BPS]);
    ref1 = BROADCAST64(&ref[1 * BPS]);
    ref2 = BROADCAST64(&ref[2 * BPS]);
    ref3 = BROADCAST64(&ref[3 * BPS]);
    {
      // convert src to 16b.
      const __m256i filter4 = _mm256_load_si256((__m256i const *)kFilter4);
      const __m256i src_0 = _mm256_unpacklo_epi8(src0, zero);
      const __m256i src_1 = _mm256_unpacklo_epi8(src1, zero);
      const __m256i src_2 = _mm256_unpacklo_epi8(src2, zero);
      const __m256i src_3 = _mm256_unpacklo_epi8(src3, zero);

      // The data in the register is now converted to 16b
      // The first 4 byte reside in the first lane and the next 4 bytes resides in the second lane
      const __m256i ref_0 = _mm256_shuffle_epi8(ref0, filter4);
      const __m256i ref_1 = _mm256_shuffle_epi8(ref1, filter4);
      const __m256i ref_2 = _mm256_shuffle_epi8(ref2, filter4);
      const __m256i ref_3 = _mm256_shuffle_epi8(ref3, filter4);

      // Compute difference. -> 00 01 02 03 00 00 00 00
      const __m256i diff0 = _mm256_sub_epi16(src_0, ref_0);
      const __m256i diff1 = _mm256_sub_epi16(src_1, ref_1);
      const __m256i diff2 = _mm256_sub_epi16(src_2, ref_2);
      const __m256i diff3 = _mm256_sub_epi16(src_3, ref_3);

      // Unpack and shuffle
      // 00 01 02 03   0 0 0 0
      // 10 11 12 13   0 0 0 0
      // 20 21 22 23   0 0 0 0
      // 30 31 32 33   0 0 0 0
      const __m256i shuf01 = _mm256_unpacklo_epi32(diff0, diff1);
      const __m256i shuf23 = _mm256_unpacklo_epi32(diff2, diff3);

      // 00 01 10 11 02 03 12 13
      // 20 21 30 31 22 23 32 33
      const __m256i shuf01_p =
          _mm256_shufflehi_epi16(shuf01, _MM_SHUFFLE(2, 3, 0, 1));
      const __m256i shuf23_p =
          _mm256_shufflehi_epi16(shuf23, _MM_SHUFFLE(2, 3, 0, 1));
      // 00 01 10 11 03 02 13 12
      // 20 21 30 31 23 22 33 32
      const __m256i s01 = _mm256_unpacklo_epi64(shuf01_p, shuf23_p);
      const __m256i s32 = _mm256_unpackhi_epi64(shuf01_p, shuf23_p);
      // 00 01 10 11 20 21 30 31
      // 03 02 13 12 23 22 33 32
      const __m256i a01 = _mm256_add_epi16(s01, s32);
      const __m256i a32 = _mm256_sub_epi16(s01, s32);
      // [d0 + d3 | d1 + d2 | ...] = [a0 a1 | a0' a1' | ... ]
      // [d0 - d3 | d1 - d2 | ...] = [a3 a2 | a3' a2' | ... ]

      const __m256i tmp0 = _mm256_madd_epi16(a01, k88p);  // [ (a0 + a1) << 3, ... ]
      const __m256i tmp2 = _mm256_madd_epi16(a01, k88m);  // [ (a0 - a1) << 3, ... ]
      const __m256i tmp1_1 = _mm256_madd_epi16(a32, k5352_2217p);
      const __m256i tmp3_1 = _mm256_madd_epi16(a32, k5352_2217m);
      const __m256i tmp1_2 = _mm256_add_epi32(tmp1_1, k1812);
      const __m256i tmp3_2 = _mm256_add_epi32(tmp3_1, k937);
      const __m256i tmp1   = _mm256_srai_epi32(tmp1_2, 9);
      const __m256i tmp3   = _mm256_srai_epi32(tmp3_2, 9);
      const __m256i s03 = _mm256_packs_epi32(tmp0, tmp2);
      const __m256i s12 = _mm256_packs_epi32(tmp1, tmp3);
      const __m256i  s_lo = _mm256_unpacklo_epi16(s03, s12);   // 0 1 0 1 0 1...
      const __m256i  s_hi = _mm256_unpackhi_epi16(s03, s12);   // 2 3 2 3 2 3
      const __m256i  v23 = _mm256_unpackhi_epi32(s_lo, s_hi);
      v01 = _mm256_unpacklo_epi32(s_lo, s_hi);
      v32 = _mm256_shuffle_epi32(v23, _MM_SHUFFLE(1, 0, 3, 2));  // 3 2 3 2 3 2..
    }
  }

  // Second pass
  {
    // Same operations are done on the (0,3) and (1,2) pairs.
    // a0 = v0 + v3
    // a1 = v1 + v2
    // a3 = v0 - v3
    // a2 = v1 - v2
    const __m256i a01 = _mm256_add_epi16(v01, v32);
    const __m256i a32 = _mm256_sub_epi16(v01, v32);
    const __m256i a11 = _mm256_unpackhi_epi64(a01, a01);
    const __m256i a22 = _mm256_unpackhi_epi64(a32, a32);
    const __m256i a01_plus_7 = _mm256_add_epi16(a01, seven);

    // d0 = (a0 + a1 + 7) >> 4;
    // d2 = (a0 - a1 + 7) >> 4;
    const __m256i c0 = _mm256_add_epi16(a01_plus_7, a11);
    const __m256i c2 = _mm256_sub_epi16(a01_plus_7, a11);
    const __m256i d0 = _mm256_srai_epi16(c0, 4);
    const __m256i d2 = _mm256_srai_epi16(c2, 4);

    // f1 = ((b3 * 5352 + b2 * 2217 + 12000) >> 16)
    // f3 = ((b3 * 2217 - b2 * 5352 + 51000) >> 16)
    const __m256i b23 = _mm256_unpacklo_epi16(a22, a32);
    const __m256i c1 = _mm256_madd_epi16(b23, k5352_2217);
    const __m256i c3 = _mm256_madd_epi16(b23, k2217_5352);
    const __m256i d1 = _mm256_add_epi32(c1, k12000_plus_one);
    const __m256i d3 = _mm256_add_epi32(c3, k51000);
    const __m256i e1 = _mm256_srai_epi32(d1, 16);
    const __m256i e3 = _mm256_srai_epi32(d3, 16);
    const __m256i f1 = _mm256_packs_epi32(e1, e1);
    const __m256i f3 = _mm256_packs_epi32(e3, e3);
    // f1 = f1 + (a3 != 0);
    // The compare will return (0xffff, 0) for (==0, !=0). To turn that into the
    // desired (0, 1), we add one earlier through k12000_plus_one.
    // -> f1 = f1 + 1 - (a3 == 0)
    const __m256i g1 = _mm256_add_epi16(f1, _mm256_cmpeq_epi16(a32, zero));

    const __m256i d0_g1 = _mm256_unpacklo_epi64(d0, g1);
    const __m256i d2_f3 = _mm256_unpacklo_epi64(d2, f3);
    _mm256_store_si256((__m256i*)&out[0], d0_g1);
    _mm256_store_si256((__m256i*)&out[16], d2_f3);
  }
}

//------------------------------------------------------------------------------
// Texture distortion
//
// We try to match the spectral content (weighted) between source and
// reconstructed samples.

// Hadamard transform
// Returns the difference between the weighted sum of the absolute value of
// transformed coefficients.
static void TTransform4x2(const uint8_t* inA,const uint8_t* inB,
                      const uint16_t* const w,
                      int64_t* SD1,
                      int64_t* SD2) {
  int32_t sum[8];
  __m256i tmp_0, tmp_1, tmp_2, tmp_3;
  const __m256i zero = _mm256_setzero_si256();

  // Load, combine and transpose inputs.
  {
    __m256i inA_0, inA_1, inA_2, inA_3; 
    // A floating point broadcast is being used here because versus to
    // "integer broadcast" this broadcast is done from the memeory
    // and not from a register (1 cycle vs 3 cycles)
    inA_0 = BROADCAST64(&inA[BPS * 0]);
    inA_1 = BROADCAST64(&inA[BPS * 1]);
    inA_2 = BROADCAST64(&inA[BPS * 2]);
    inA_3 = BROADCAST64(&inA[BPS * 3]);
    {
      const __m256i inB_0 = _mm256_load_si256((__m256i const *)(inB + 32*0));
      const __m256i inB_1 = _mm256_load_si256((__m256i const *)(inB + 32*1));
      const __m256i inB_2 = _mm256_load_si256((__m256i const *)(inB + 32*2));
      const __m256i inB_3 = _mm256_load_si256((__m256i const *)(inB + 32*3));

      // Combine inA and inB (we'll do two transforms in parallel).
      const __m256i inAB_0 = _mm256_unpacklo_epi8(inA_0, inB_0);
      const __m256i inAB_1 = _mm256_unpacklo_epi8(inA_1, inB_1);
      const __m256i inAB_2 = _mm256_unpacklo_epi8(inA_2, inB_2);
      const __m256i inAB_3 = _mm256_unpacklo_epi8(inA_3, inB_3);
      // a00 b00 a01 b01 a02 b03 a03 b03   0 0 0 0 0 0 0 0
      // a10 b10 a11 b11 a12 b12 a13 b13   0 0 0 0 0 0 0 0
      // a20 b20 a21 b21 a22 b22 a23 b23   0 0 0 0 0 0 0 0
      // a30 b30 a31 b31 a32 b32 a33 b33   0 0 0 0 0 0 0 0

      // Transpose the two 4x4, discarding the filling zeroes.
      const __m256i transpose0_0 = _mm256_unpacklo_epi8(inAB_0, inAB_2);
      const __m256i transpose0_1 = _mm256_unpacklo_epi8(inAB_1, inAB_3);
      // a00 a20  b00 b20  a01 a21  b01 b21  a02 a22  b02 b22  a03 a23  b03 b23
      // a10 a30  b10 b30  a11 a31  b11 b31  a12 a32  b12 b32  a13 a33  b13 b33
      const __m256i transpose1_0 = _mm256_unpacklo_epi8(transpose0_0, transpose0_1);
      const __m256i transpose1_1 = _mm256_unpackhi_epi8(transpose0_0, transpose0_1);
      // a00 a10 a20 a30  b00 b10 b20 b30  a01 a11 a21 a31  b01 b11 b21 b31
      // a02 a12 a22 a32  b02 b12 b22 b32  a03 a13 a23 a33  b03 b13 b23 b33

      // Convert to 16b.
      tmp_0 = _mm256_unpacklo_epi8(transpose1_0, zero);
      tmp_1 = _mm256_unpackhi_epi8(transpose1_0, zero);
      tmp_2 = _mm256_unpacklo_epi8(transpose1_1, zero);
      tmp_3 = _mm256_unpackhi_epi8(transpose1_1, zero);
      // a00 a10 a20 a30   b00 b10 b20 b30
      // a01 a11 a21 a31   b01 b11 b21 b31
      // a02 a12 a22 a32   b02 b12 b22 b32
      // a03 a13 a23 a33   b03 b13 b23 b33
    }
  }
  // Horizontal pass and subsequent transpose.
  {
    // Calculate a and b (two 4x4 at once).
    const __m256i a0 = _mm256_add_epi16(tmp_0, tmp_2);
    const __m256i a1 = _mm256_add_epi16(tmp_1, tmp_3);
    const __m256i a2 = _mm256_sub_epi16(tmp_1, tmp_3);
    const __m256i a3 = _mm256_sub_epi16(tmp_0, tmp_2);
    const __m256i b0 = _mm256_add_epi16(a0, a1);
    const __m256i b1 = _mm256_add_epi16(a3, a2);
    const __m256i b2 = _mm256_sub_epi16(a3, a2);
    const __m256i b3 = _mm256_sub_epi16(a0, a1);
    // a00 a01 a02 a03   b00 b01 b02 b03
    // a10 a11 a12 a13   b10 b11 b12 b13
    // a20 a21 a22 a23   b20 b21 b22 b23
    // a30 a31 a32 a33   b30 b31 b32 b33

    // Transpose the two 4x4.
    const __m256i transpose0_0 = _mm256_unpacklo_epi16(b0, b1);
    const __m256i transpose0_1 = _mm256_unpacklo_epi16(b2, b3);
    const __m256i transpose0_2 = _mm256_unpackhi_epi16(b0, b1);
    const __m256i transpose0_3 = _mm256_unpackhi_epi16(b2, b3);
    // a00 a10 a01 a11   a02 a12 a03 a13
    // a20 a30 a21 a31   a22 a32 a23 a33
    // b00 b10 b01 b11   b02 b12 b03 b13
    // b20 b30 b21 b31   b22 b32 b23 b33
    const __m256i transpose1_0 = _mm256_unpacklo_epi32(transpose0_0, transpose0_1);
    const __m256i transpose1_1 = _mm256_unpacklo_epi32(transpose0_2, transpose0_3);
    const __m256i transpose1_2 = _mm256_unpackhi_epi32(transpose0_0, transpose0_1);
    const __m256i transpose1_3 = _mm256_unpackhi_epi32(transpose0_2, transpose0_3);
    // a00 a10 a20 a30 a01 a11 a21 a31
    // b00 b10 b20 b30 b01 b11 b21 b31
    // a02 a12 a22 a32 a03 a13 a23 a33
    // b02 b12 a22 b32 b03 b13 b23 b33
    tmp_0 = _mm256_unpacklo_epi64(transpose1_0, transpose1_1);
    tmp_1 = _mm256_unpackhi_epi64(transpose1_0, transpose1_1);
    tmp_2 = _mm256_unpacklo_epi64(transpose1_2, transpose1_3);
    tmp_3 = _mm256_unpackhi_epi64(transpose1_2, transpose1_3);
    // a00 a10 a20 a30   b00 b10 b20 b30
    // a01 a11 a21 a31   b01 b11 b21 b31
    // a02 a12 a22 a32   b02 b12 b22 b32
    // a03 a13 a23 a33   b03 b13 b23 b33
  }

  // Vertical pass and difference of weighted sums.
  {
    __m256i w_0, w_8;
    // Load all inputs.
    w_0 = BROADCAST128(&w[0]);
    w_8 = BROADCAST128(&w[8]);
    {
      // Calculate a and b (two 4x4 at once).
      const __m256i a0 = _mm256_add_epi16(tmp_0, tmp_2);
      const __m256i a1 = _mm256_add_epi16(tmp_1, tmp_3);
      const __m256i a2 = _mm256_sub_epi16(tmp_1, tmp_3);
      const __m256i a3 = _mm256_sub_epi16(tmp_0, tmp_2);
      const __m256i b0 = _mm256_add_epi16(a0, a1);
      const __m256i b1 = _mm256_add_epi16(a3, a2);
      const __m256i b2 = _mm256_sub_epi16(a3, a2);
      const __m256i b3 = _mm256_sub_epi16(a0, a1);

      // Separate the transforms of inA and inB.
      __m256i A_b0 = _mm256_unpacklo_epi64(b0, b1);
      __m256i A_b2 = _mm256_unpacklo_epi64(b2, b3);
      __m256i B_b0 = _mm256_unpackhi_epi64(b0, b1);
      __m256i B_b2 = _mm256_unpackhi_epi64(b2, b3);

      {
        // sign(b) = b >> 15  (0x0000 if positive, 0xffff if negative)
        const __m256i sign_A_b0 = _mm256_srai_epi16(A_b0, 15);
        const __m256i sign_A_b2 = _mm256_srai_epi16(A_b2, 15);
        const __m256i sign_B_b0 = _mm256_srai_epi16(B_b0, 15);
        const __m256i sign_B_b2 = _mm256_srai_epi16(B_b2, 15);

        // b = abs(b) = (b ^ sign) - sign
        A_b0 = _mm256_xor_si256(A_b0, sign_A_b0);
        A_b2 = _mm256_xor_si256(A_b2, sign_A_b2);
        B_b0 = _mm256_xor_si256(B_b0, sign_B_b0);
        B_b2 = _mm256_xor_si256(B_b2, sign_B_b2);
        A_b0 = _mm256_sub_epi16(A_b0, sign_A_b0);
        A_b2 = _mm256_sub_epi16(A_b2, sign_A_b2);
        B_b0 = _mm256_sub_epi16(B_b0, sign_B_b0);
        B_b2 = _mm256_sub_epi16(B_b2, sign_B_b2);
      }

      // weighted sums
      A_b0 = _mm256_madd_epi16(A_b0, w_0);
      A_b2 = _mm256_madd_epi16(A_b2, w_8);
      B_b0 = _mm256_madd_epi16(B_b0, w_0);
      B_b2 = _mm256_madd_epi16(B_b2, w_8);
      A_b0 = _mm256_add_epi32(A_b0, A_b2);
      B_b0 = _mm256_add_epi32(B_b0, B_b2);

      // difference of weighted sums
      A_b0 = _mm256_sub_epi32(A_b0, B_b0);
      _mm256_storeu_si256((__m256i*)&sum[0], A_b0);
    }
  }
  *SD1 =  sum[0] + sum[1] + sum[2] + sum[3];
  *SD2 =  sum[4] + sum[5] + sum[6] + sum[7];
}

static void Disto4x4x2(const uint8_t* const a, const uint8_t* b,
                          const uint16_t* const w,
                          int64_t* SD1,
                          int64_t* SD2) {
  TTransform4x2(a, b, w, SD1, SD2);
  *SD1 = abs(*SD1) >> 5;
  *SD2 = abs(*SD2) >> 5;
}

#define QFIX2 0
static WEBP_INLINE void DoQuantizeBlock4x2(int16_t in[32], int16_t out[32],
                                           int shift,
                                           const uint16_t* const sharpen,
                                           const VP8Matrix* const mtx,
                                           uint32_t* nz1,
                                           uint32_t* nz2) {
  const __m256i max_coeff_2047 = _mm256_set1_epi16(MAX_LEVEL);
  const __m256i zero = _mm256_setzero_si256();
  __m256i coeff0, coeff8;
  __m256i out0, out8;
  __m256i packed_out;
  __m256i iq0, iq8, q0, q8;
  uint32_t res;

  // Load all inputs.
  __m256i in0 = _mm256_load_si256((__m256i*)&in[0]);
  __m256i in8 = _mm256_load_si256((__m256i*)&in[16]);

  // extract sign(in)  (0x0000 if positive, 0xffff if negative)
  const __m256i sign0 = _mm256_cmpgt_epi16(zero, in0);
  const __m256i sign8 = _mm256_cmpgt_epi16(zero, in8);

  iq0 = BROADCAST128(&mtx->iq_[0]);
  iq8 = BROADCAST128(&mtx->iq_[8]);
  q0 = BROADCAST128(&mtx->q_[0]);
  q8 = BROADCAST128(&mtx->q_[8]);

  // coeff = abs(in) = (in ^ sign) - sign
  coeff0 = _mm256_xor_si256(in0, sign0);
  coeff8 = _mm256_xor_si256(in8, sign8);
  coeff0 = _mm256_sub_epi16(coeff0, sign0);
  coeff8 = _mm256_sub_epi16(coeff8, sign8);

  // coeff = abs(in) + sharpen
  if (sharpen != NULL) {
    __m256i sharpen0, sharpen8; 
    sharpen0 = BROADCAST128(&sharpen[0]);
    sharpen8 = BROADCAST128(&sharpen[8]);
    coeff0 = _mm256_add_epi16(coeff0, sharpen0);
    coeff8 = _mm256_add_epi16(coeff8, sharpen8);
  }

  // out = (coeff * iQ + B) >> (QFIX + QFIX2 - shift)
  {
    __m256i bias_00, bias_04, bias_08, bias_12;
    // doing calculations with 32b precision (QFIX=17)
    // out = (coeff * iQ)
    const __m256i coeff_iQ0H = _mm256_mulhi_epu16(coeff0, iq0);
    const __m256i coeff_iQ0L = _mm256_mullo_epi16(coeff0, iq0);
    const __m256i coeff_iQ8H = _mm256_mulhi_epu16(coeff8, iq8);
    const __m256i coeff_iQ8L = _mm256_mullo_epi16(coeff8, iq8);
    __m256i out_00 = _mm256_unpacklo_epi16(coeff_iQ0L, coeff_iQ0H);
    __m256i out_04 = _mm256_unpackhi_epi16(coeff_iQ0L, coeff_iQ0H);
    __m256i out_08 = _mm256_unpacklo_epi16(coeff_iQ8L, coeff_iQ8H);
    __m256i out_12 = _mm256_unpackhi_epi16(coeff_iQ8L, coeff_iQ8H);
    // out = (coeff * iQ + B)
    bias_00 = BROADCAST128(&mtx->bias_[0]);
    bias_04 = BROADCAST128(&mtx->bias_[4]);
    bias_08 = BROADCAST128(&mtx->bias_[8]);
    bias_12 = BROADCAST128(&mtx->bias_[12]);
    out_00 = _mm256_add_epi32(out_00, bias_00);
    out_04 = _mm256_add_epi32(out_04, bias_04);
    out_08 = _mm256_add_epi32(out_08, bias_08);
    out_12 = _mm256_add_epi32(out_12, bias_12);
    // out = QUANTDIV(coeff, iQ, B, QFIX + QFIX2 - shift)
    out_00 = _mm256_srai_epi32(out_00, QFIX + QFIX2 - shift);
    out_04 = _mm256_srai_epi32(out_04, QFIX + QFIX2 - shift);
    out_08 = _mm256_srai_epi32(out_08, QFIX + QFIX2 - shift);
    out_12 = _mm256_srai_epi32(out_12, QFIX + QFIX2 - shift);

    // pack result as 16b
    out0 = _mm256_packs_epi32(out_00, out_04);
    out8 = _mm256_packs_epi32(out_08, out_12);

    // if (coeff > 2047) coeff = 2047
    out0 = _mm256_min_epi16(out0, max_coeff_2047);
    out8 = _mm256_min_epi16(out8, max_coeff_2047);
  }

  // get sign back (if (sign[j]) out_n = -out_n)
  out0 = _mm256_xor_si256(out0, sign0);
  out8 = _mm256_xor_si256(out8, sign8);
  out0 = _mm256_sub_epi16(out0, sign0);
  out8 = _mm256_sub_epi16(out8, sign8);

  // in = out * Q
  in0 = _mm256_mullo_epi16(out0, q0);
  in8 = _mm256_mullo_epi16(out8, q8);

  _mm256_storeu_si256((__m256i*)&in[0], in0);
  _mm256_storeu_si256((__m256i*)&in[16], in8);

  // zigzag the output before storing it.
  //
  // The zigzag pattern can almost be reproduced with a small sequence of
  // shuffles. After it, we only need to swap the 7th (ending up in third
  // position instead of twelfth) and 8th values.
  {
    __m256i outZ0, outZ8 , replace12, replace3;
    outZ0 = _mm256_shufflehi_epi16(out0,  _MM_SHUFFLE(2, 1, 3, 0));
    outZ0 = _mm256_shuffle_epi32  (outZ0, _MM_SHUFFLE(3, 1, 2, 0));
    outZ0 = _mm256_shufflehi_epi16(outZ0, _MM_SHUFFLE(3, 1, 0, 2));
    outZ8 = _mm256_shufflelo_epi16(out8,  _MM_SHUFFLE(3, 0, 2, 1));
    outZ8 = _mm256_shuffle_epi32  (outZ8, _MM_SHUFFLE(3, 1, 2, 0));
    outZ8 = _mm256_shufflelo_epi16(outZ8, _MM_SHUFFLE(1, 3, 2, 0));
    //switch between the third element to the 12th element
    replace12 = _mm256_srli_si256(outZ8, 2);
    replace3 = _mm256_slli_si256(outZ0 , 2);
    outZ0 = _mm256_blend_epi16(outZ0, replace12, 8);
    outZ8 = _mm256_blend_epi16(outZ8, replace3, 16);
    _mm_storeu_si128((__m128i*)&out[0], _mm256_castsi256_si128(outZ0));
    _mm_storeu_si128((__m128i*)&out[8], _mm256_castsi256_si128(outZ8));
    _mm_storeu_si128((__m128i*)&out[16], _mm256_extracti128_si256(outZ0, 1));
    _mm_storeu_si128((__m128i*)&out[24], _mm256_extracti128_si256(outZ8, 1));
    packed_out = _mm256_packs_epi16(outZ0, outZ8);
  }

  // detect if all 'out' values are zeroes or not
  res = _mm256_movemask_epi8(_mm256_cmpeq_epi8(packed_out, zero));
  *nz1 = ((res & 0xffff) != 0xffff);
  *nz2 = ((res >> 16) != 0xffff);
}


static void QuantizeBlock4x2(int16_t in[32], int16_t out[32],
                             const VP8Matrix* const mtx,
                             uint32_t* nz1, uint32_t* nz2) {
  DoQuantizeBlock4x2(in, out, 0, &mtx->sharpen_[0], mtx, nz1, nz2);
}

extern int PickBestIntra4x2(struct VP8EncIterator* const it, struct VP8ModeScore* const rd);
#endif  // WEBP_USE_AVX2

//------------------------------------------------------------------------------
// Entry point

extern void VP8EncDspInitAVX2(void);

void VP8EncDspInitAVX2(void) {
#if defined(WEBP_USE_AVX2)
  VP8EncQuantizeBlock4x2 = QuantizeBlock4x2;
  VP8ITransform4x2 = ITransform4x2;
  VP8FTransform4x2 = FTransform4x2;
  VP8SSE4x4x2 = SSE4x4x2;
  VP8TDisto4x4x2 = Disto4x4x2;
  VP8EncPickBestIntra4 = PickBestIntra4x2;
#endif
}
