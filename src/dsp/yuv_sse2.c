// Copyright 2014 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
// YUV->RGB conversion functions
//
// Author: Skal (pascal.massimino@gmail.com)

#include "./yuv.h"

#if defined(WEBP_USE_SSE2)

#include <emmintrin.h>
#include <string.h>   // for memcpy

//-----------------------------------------------------------------------------
// Convert spans of 32 pixels to various RGB formats for the fancy upsampler.

// These constants are 14b fixed-point version of ITU-R BT.601 constants.
// R = (19077 * y             + 26149 * v - 14234) >> 6
// G = (19077 * y -  6419 * u - 13320 * v +  8708) >> 6
// B = (19077 * y + 33050 * u             - 17685) >> 6
static void ConvertYUV444ToRGB(const __m128i* const Y0,
                               const __m128i* const U0,
                               const __m128i* const V0,
                               __m128i* const R,
                               __m128i* const G,
                               __m128i* const B) {
  const __m128i k19077 = _mm_set1_epi16(19077);
  const __m128i k26149 = _mm_set1_epi16(26149);
  const __m128i k14234 = _mm_set1_epi16(14234);
  const __m128i k33050 = _mm_set1_epi16(33050);
  const __m128i k17685 = _mm_set1_epi16(17685);
  const __m128i k6419  = _mm_set1_epi16(6419);
  const __m128i k13320 = _mm_set1_epi16(13320);
  const __m128i k8708  = _mm_set1_epi16(8708);

  const __m128i Y1 = _mm_mulhi_epu16(*Y0, k19077);

  const __m128i R0 = _mm_mulhi_epu16(*V0, k26149);
  const __m128i R1 = _mm_sub_epi16(Y1, k14234);
  const __m128i R2 = _mm_add_epi16(R1, R0);

  const __m128i G0 = _mm_mulhi_epu16(*U0, k6419);
  const __m128i G1 = _mm_mulhi_epu16(*V0, k13320);
  const __m128i G2 = _mm_add_epi16(Y1, k8708);
  const __m128i G3 = _mm_add_epi16(G0, G1);
  const __m128i G4 = _mm_sub_epi16(G2, G3);

  // be careful with the saturated *unsigned* arithmetic here!
  const __m128i B0 = _mm_mulhi_epu16(*U0, k33050);
  const __m128i B1 = _mm_adds_epu16(B0, Y1);
  const __m128i B2 = _mm_subs_epu16(B1, k17685);

  // use logical shift for B2, which can be larger than 32767
  *R = _mm_srai_epi16(R2, 6);   // range: [-14234, 30815]
  *G = _mm_srai_epi16(G4, 6);   // range: [-10953, 27710]
  *B = _mm_srli_epi16(B2, 6);   // range: [0, 34238]
}

// Load the bytes into the *upper* part of 16b words. That's "<< 8", basically.
static WEBP_INLINE __m128i Load_HI_16(const uint8_t* src) {
  const __m128i zero = _mm_setzero_si128();
  return _mm_unpacklo_epi8(zero, _mm_loadl_epi64((const __m128i*)src));
}

// Load and replicate the U/V samples
static WEBP_INLINE __m128i Load_UV_HI_8(const uint8_t* src) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i tmp0 = _mm_cvtsi32_si128(*(const uint32_t*)src);
  const __m128i tmp1 = _mm_unpacklo_epi8(zero, tmp0);
  return _mm_unpacklo_epi16(tmp1, tmp1);   // replicate samples
}

// Convert 32 samples of YUV444 to R/G/B
static void YUV444ToRGB(const uint8_t* const y,
                        const uint8_t* const u,
                        const uint8_t* const v,
                        __m128i* const R, __m128i* const G, __m128i* const B) {
  const __m128i Y0 = Load_HI_16(y), U0 = Load_HI_16(u), V0 = Load_HI_16(v);
  ConvertYUV444ToRGB(&Y0, &U0, &V0, R, G, B);
}

// Convert 32 samples of YUV420 to R/G/B
static void YUV420ToRGB(const uint8_t* const y,
                        const uint8_t* const u,
                        const uint8_t* const v,
                        __m128i* const R, __m128i* const G, __m128i* const B) {
  const __m128i Y0 = Load_HI_16(y), U0 = Load_UV_HI_8(u), V0 = Load_UV_HI_8(v);
  ConvertYUV444ToRGB(&Y0, &U0, &V0, R, G, B);
}

// Pack R/G/B/A results into 32b output.
static WEBP_INLINE void PackAndStore4(const __m128i* const R,
                                      const __m128i* const G,
                                      const __m128i* const B,
                                      const __m128i* const A,
                                      uint8_t* const dst) {
  const __m128i rb = _mm_packus_epi16(*R, *B);
  const __m128i ga = _mm_packus_epi16(*G, *A);
  const __m128i rg = _mm_unpacklo_epi8(rb, ga);
  const __m128i ba = _mm_unpackhi_epi8(rb, ga);
  const __m128i RGBA_lo = _mm_unpacklo_epi16(rg, ba);
  const __m128i RGBA_hi = _mm_unpackhi_epi16(rg, ba);
  _mm_storeu_si128((__m128i*)(dst +  0), RGBA_lo);
  _mm_storeu_si128((__m128i*)(dst + 16), RGBA_hi);
}

// Pack R/G/B results into 24b output.
static WEBP_INLINE void PackAndStore3(const __m128i* const R,
                                      const __m128i* const G,
                                      const __m128i* const B,
                                      uint8_t* const dst) {
  const __m128i tmp0 = _mm_packus_epi16(*R, *R);
  const __m128i tmp1 = _mm_packus_epi16(*G, *G);
  const __m128i tmp2 = _mm_packus_epi16(*B, *B);
  _mm_storel_epi64((__m128i*)(dst + 0 * 32), tmp0);
  _mm_storel_epi64((__m128i*)(dst + 1 * 32), tmp1);
  _mm_storel_epi64((__m128i*)(dst + 2 * 32), tmp2);
}

// Converts 32 samples in src[3][32] to interleaved RGB24 in dst[]
static WEBP_INLINE void PlanarTo24b(const uint8_t* src, uint8_t* dst) {
#if 1
  // This code is faster than the version below (left there for reference).
  // It's also endian-dependent but we're only targeting x86.
  const uint8_t* const end = src + 32;
  for (; src < end; src += 4, dst += 12) {
    const uint32_t A = (src[0 + 0 * 32] <<  0) | (src[0 + 1 * 32] <<  8)
                     | (src[0 + 2 * 32] << 16) | (src[1 + 0 * 32] << 24);
    const uint32_t B = (src[1 + 1 * 32] <<  0) | (src[1 + 2 * 32] <<  8)
                     | (src[2 + 0 * 32] << 16) | (src[2 + 1 * 32] << 24);
    const uint32_t C = (src[2 + 2 * 32] <<  0) | (src[3 + 0 * 32] <<  8)
                     | (src[3 + 1 * 32] << 16) | (src[3 + 2 * 32] << 24);
    *(uint32_t*)(dst + 0) = A;
    *(uint32_t*)(dst + 4) = B;
    *(uint32_t*)(dst + 8) = C;
  }
#else
  int n;
  for (n = 0; n < 32; ++n) {
    *dst++ = src[0 * 32 + n];
    *dst++ = src[1 * 32 + n];
    *dst++ = src[2 * 32 + n];
    }
#endif
}

void VP8YuvToRgba32(const uint8_t* y, const uint8_t* u, const uint8_t* v,
                    uint8_t* dst) {
  const __m128i kAlpha = _mm_set1_epi16(255);
  int n;
  for (n = 0; n < 32; n += 8, dst += 32) {
    __m128i R, G, B;
    YUV444ToRGB(y + n, u + n, v + n, &R, &G, &B);
    PackAndStore4(&R, &G, &B, &kAlpha, dst);
  }
}

void VP8YuvToBgra32(const uint8_t* y, const uint8_t* u, const uint8_t* v,
                    uint8_t* dst) {
  const __m128i kAlpha = _mm_set1_epi16(255);
  int n;
  for (n = 0; n < 32; n += 8, dst += 32) {
    __m128i R, G, B;
    YUV444ToRGB(y + n, u + n, v + n, &R, &G, &B);
    PackAndStore4(&B, &G, &R, &kAlpha, dst);
  }
}

void VP8YuvToRgb32(const uint8_t* y, const uint8_t* u, const uint8_t* v,
                   uint8_t* dst) {
  int n;
  uint8_t tmp[32 * 3];
  for (n = 0; n < 32; n += 8) {
    __m128i R, G, B;
    YUV444ToRGB(y + n, u + n, v + n, &R, &G, &B);
    PackAndStore3(&R, &G, &B, tmp + n);
  }
  PlanarTo24b(tmp, dst);
}

void VP8YuvToBgr32(const uint8_t* y, const uint8_t* u, const uint8_t* v,
                   uint8_t* dst) {
  int n;
  uint8_t tmp[32 * 3];
  for (n = 0; n < 32; n += 8) {
    __m128i R, G, B;
    YUV444ToRGB(y + n, u + n, v + n, &R, &G, &B);
    PackAndStore3(&B, &G, &R, tmp + n);
  }
  PlanarTo24b(tmp, dst);
}

//-----------------------------------------------------------------------------
// Arbitrary-length row conversion functions

static void YuvToRgbaRow(const uint8_t* y, const uint8_t* u, const uint8_t* v,
                         uint8_t* dst, int len) {
  const __m128i kAlpha = _mm_set1_epi16(255);
  int n;
  for (n = 0; n + 8 <= len; n += 8, dst += 32) {
    __m128i R, G, B;
    YUV420ToRGB(y, u, v, &R, &G, &B);
    PackAndStore4(&R, &G, &B, &kAlpha, dst);
    y += 8;
    u += 4;
    v += 4;
  }
  for (; n < len; ++n) {   // Finish off
    VP8YuvToRgba(y[0], u[0], v[0], dst);
    dst += 4;
    y += 1;
    u += (n & 1);
    v += (n & 1);
  }
}

static void YuvToBgraRow(const uint8_t* y, const uint8_t* u, const uint8_t* v,
                         uint8_t* dst, int len) {
  const __m128i kAlpha = _mm_set1_epi16(255);
  int n;
  for (n = 0; n + 8 <= len; n += 8, dst += 32) {
    __m128i R, G, B;
    YUV420ToRGB(y, u, v, &R, &G, &B);
    PackAndStore4(&B, &G, &R, &kAlpha, dst);
    y += 8;
    u += 4;
    v += 4;
  }
  for (; n < len; ++n) {   // Finish off
    VP8YuvToBgra(y[0], u[0], v[0], dst);
    dst += 4;
    y += 1;
    u += (n & 1);
    v += (n & 1);
  }
}

static void YuvToArgbRow(const uint8_t* y, const uint8_t* u, const uint8_t* v,
                         uint8_t* dst, int len) {
  const __m128i kAlpha = _mm_set1_epi16(255);
  int n;
  for (n = 0; n + 8 <= len; n += 8, dst += 32) {
    __m128i R, G, B;
    YUV420ToRGB(y, u, v, &R, &G, &B);
    PackAndStore4(&kAlpha, &R, &G, &B, dst);
    y += 8;
    u += 4;
    v += 4;
  }
  for (; n < len; ++n) {   // Finish off
    VP8YuvToArgb(y[0], u[0], v[0], dst);
    dst += 4;
    y += 1;
    u += (n & 1);
    v += (n & 1);
  }
}

static void YuvToRgbRow(const uint8_t* y, const uint8_t* u, const uint8_t* v,
                        uint8_t* dst, int len) {
  int n;
  for (n = 0; n + 32 <= len; n += 32, dst += 32 * 3) {
    uint8_t tmp[32 * 3];
    __m128i R, G, B;
    YUV420ToRGB(y +  0, u +  0, v +  0, &R, &G, &B);
    PackAndStore3(&R, &G, &B, tmp +  0);
    YUV420ToRGB(y +  8, u +  4, v +  4, &R, &G, &B);
    PackAndStore3(&R, &G, &B, tmp +  8);
    YUV420ToRGB(y + 16, u +  8, v +  8, &R, &G, &B);
    PackAndStore3(&R, &G, &B, tmp + 16);
    YUV420ToRGB(y + 24, u + 12, v + 12, &R, &G, &B);
    PackAndStore3(&R, &G, &B, tmp + 24);
    PlanarTo24b(tmp, dst);
    y += 32;
    u += 16;
    v += 16;
  }
  for (; n < len; ++n) {   // Finish off
    VP8YuvToRgb(y[0], u[0], v[0], dst);
    dst += 3;
    y += 1;
    u += (n & 1);
    v += (n & 1);
  }
}

static void YuvToBgrRow(const uint8_t* y, const uint8_t* u, const uint8_t* v,
                        uint8_t* dst, int len) {
  int n;
  for (n = 0; n + 32 <= len; n += 32, dst += 32 * 3) {
    uint8_t tmp[32 * 3];
    __m128i R, G, B;
    YUV420ToRGB(y +  0, u +  0, v +  0, &R, &G, &B);
    PackAndStore3(&B, &G, &R, tmp +  0);
    YUV420ToRGB(y +  8, u +  4, v +  4, &R, &G, &B);
    PackAndStore3(&B, &G, &R, tmp +  8);
    YUV420ToRGB(y + 16, u +  8, v +  8, &R, &G, &B);
    PackAndStore3(&B, &G, &R, tmp + 16);
    YUV420ToRGB(y + 24, u + 12, v + 12, &R, &G, &B);
    PackAndStore3(&B, &G, &R, tmp + 24);
    PlanarTo24b(tmp, dst);
    y += 32;
    u += 16;
    v += 16;
  }
  for (; n < len; ++n) {   // Finish off
    VP8YuvToBgr(y[0], u[0], v[0], dst);
    dst += 3;
    y += 1;
    u += (n & 1);
    v += (n & 1);
  }
}

//------------------------------------------------------------------------------
// Entry point

extern void WebPInitSamplersSSE2(void);

WEBP_TSAN_IGNORE_FUNCTION void WebPInitSamplersSSE2(void) {
  WebPSamplers[MODE_RGB]  = YuvToRgbRow;
  WebPSamplers[MODE_RGBA] = YuvToRgbaRow;
  WebPSamplers[MODE_BGR]  = YuvToBgrRow;
  WebPSamplers[MODE_BGRA] = YuvToBgraRow;
  WebPSamplers[MODE_ARGB] = YuvToArgbRow;
}

//------------------------------------------------------------------------------
// RGB24/32 -> YUV converters

// Load eight 16b-words from *src.
#define LOAD_16(src) _mm_loadu_si128((const __m128i*)(src))
// Store either 16b-words into *dst
#define STORE_16(V, dst) _mm_storeu_si128((__m128i*)(dst), (V))

// Convert 8 packed RGB or BGR samples to r[], g[], b[]
static WEBP_INLINE void RGB24PackedToPlanar(const uint8_t* const rgb,
                                            __m128i* const r,
                                            __m128i* const g,
                                            __m128i* const b,
                                            int input_is_bgr) {
  const __m128i zero = _mm_setzero_si128();
  // in0: r0 g0 b0 r1 | g1 b1 r2 g2 | b2 r3 g3 b3 | r4 g4 b4 r5
  // in1: b2 r3 g3 b3 | r4 g4 b4 r5 | g5 b5 r6 g6 | b6 r7 g7 b7
  const __m128i in0 = LOAD_16(rgb + 0);
  const __m128i in1 = LOAD_16(rgb + 8);
  // A0: | r2 g2 b2 r3 | g3 b3 r4 g4 | b4 r5 ...
  // A1:                   ... b2 r3 | g3 b3 r4 g4 | b4 r5 g5 b5 |
  const __m128i A0 = _mm_srli_si128(in0, 6);
  const __m128i A1 = _mm_slli_si128(in1, 6);
  // B0: r0 r2 g0 g2 | b0 b2 r1 r3 | g1 g3 b1 b3 | r2 r4 b2 b4
  // B1: g3 g5 b3 b5 | r4 r6 g4 g6 | b4 b6 r5 r7 | g5 g7 b5 b7
  const __m128i B0 = _mm_unpacklo_epi8(in0, A0);
  const __m128i B1 = _mm_unpackhi_epi8(A1, in1);
  // C0: r1 r3 g1 g3 | b1 b3 r2 r4 | b2 b4 ...
  // C1:                 ... g3 g5 | b3 b5 r4 r6 | g4 g6 b4 b6
  const __m128i C0 = _mm_srli_si128(B0, 6);
  const __m128i C1 = _mm_slli_si128(B1, 6);
  // D0: r0 r1 r2 r3 | g0 g1 g2 g3 | b0 b1 b2 b3 | r1 r2 r3 r4
  // D1: b3 b4 b5 b6 | r4 r5 r6 r7 | g4 g5 g6 g7 | b4 b5 b6 b7 |
  const __m128i D0 = _mm_unpacklo_epi8(B0, C0);
  const __m128i D1 = _mm_unpackhi_epi8(C1, B1);
  // r4 r5 r6 r7 | g4 g5 g6 g7 | b4 b5 b6 b7 | 0
  const __m128i D2 = _mm_srli_si128(D1, 4);
  // r0 r1 r2 r3 | r4 r5 r6 r7 | g0 g1 g2 g3 | g4 g5 g6 g7
  const __m128i E0 = _mm_unpacklo_epi32(D0, D2);
  // b0 b1 b2 b3 | b4 b5 b6 b7 | r1 r2 r3 r4 | 0
  const __m128i E1 = _mm_unpackhi_epi32(D0, D2);
  // g0 g1 g2 g3 | g4 g5 g6 g7 | 0
  const __m128i E2 = _mm_srli_si128(E0, 8);
  const __m128i F0 = _mm_unpacklo_epi8(E0, zero);  // -> R
  const __m128i F1 = _mm_unpacklo_epi8(E1, zero);  // -> B
  const __m128i F2 = _mm_unpacklo_epi8(E2, zero);  // -> G
  *g = F2;
  if (input_is_bgr) {
    *r = F1;
    *b = F0;
  } else {
    *r = F0;
    *b = F1;
  }
}

// Convert 8 packed ARGB to r[], g[], b[]
static WEBP_INLINE void RGB32PackedToPlanar(const uint32_t* const argb,
                                            __m128i* const r,
                                            __m128i* const g,
                                            __m128i* const b) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i in0 = LOAD_16(argb + 0);    // argb3 | argb2 | argb1 | argb0
  const __m128i in1 = LOAD_16(argb + 4);    // argb7 | argb6 | argb5 | argb4
  // column-wise transpose
  const __m128i A0 = _mm_unpacklo_epi8(in0, in1);
  const __m128i A1 = _mm_unpackhi_epi8(in0, in1);
  const __m128i B0 = _mm_unpacklo_epi8(A0, A1);
  const __m128i B1 = _mm_unpackhi_epi8(A0, A1);
  // C0 = g7 g6 ... g1 g0 | b7 b6 ... b1 b0
  // C1 = a7 a6 ... a1 a0 | r7 r6 ... r1 r0
  const __m128i C0 = _mm_unpacklo_epi8(B0, B1);
  const __m128i C1 = _mm_unpackhi_epi8(B0, B1);
  // store 16b
  *r = _mm_unpacklo_epi8(C1, zero);
  *g = _mm_unpackhi_epi8(C0, zero);
  *b = _mm_unpacklo_epi8(C0, zero);
}

// This macro computes (RG * MULT_RG + GB * MULT_GB + ROUNDER) >> DESCALE_FIX
// It's a macro and not a function because we need to use immediate values with
// srai_epi32, e.g.
#define TRANSFORM(RG_LO, RG_HI, GB_LO, GB_HI, MULT_RG, MULT_GB, \
                  ROUNDER, DESCALE_FIX, OUT) do {               \
  const __m128i V0_lo = _mm_madd_epi16(RG_LO, MULT_RG);         \
  const __m128i V0_hi = _mm_madd_epi16(RG_HI, MULT_RG);         \
  const __m128i V1_lo = _mm_madd_epi16(GB_LO, MULT_GB);         \
  const __m128i V1_hi = _mm_madd_epi16(GB_HI, MULT_GB);         \
  const __m128i V2_lo = _mm_add_epi32(V0_lo, V1_lo);            \
  const __m128i V2_hi = _mm_add_epi32(V0_hi, V1_hi);            \
  const __m128i V3_lo = _mm_add_epi32(V2_lo, ROUNDER);          \
  const __m128i V3_hi = _mm_add_epi32(V2_hi, ROUNDER);          \
  const __m128i V5_lo = _mm_srai_epi32(V3_lo, DESCALE_FIX);     \
  const __m128i V5_hi = _mm_srai_epi32(V3_hi, DESCALE_FIX);     \
  (OUT) = _mm_packs_epi32(V5_lo, V5_hi);                        \
} while (0)

#define MK_CST_16(A, B) _mm_set_epi16((B), (A), (B), (A), (B), (A), (B), (A))
static WEBP_INLINE void ConvertRGBToY(const __m128i* const R,
                                      const __m128i* const G,
                                      const __m128i* const B,
                                      __m128i* const Y) {
  const __m128i kRG_y = MK_CST_16(16839, 33059 - 16384);
  const __m128i kGB_y = MK_CST_16(16384, 6420);
  const __m128i kHALF_Y = _mm_set1_epi32((16 << YUV_FIX) + YUV_HALF);

  const __m128i RG_lo = _mm_unpacklo_epi16(*R, *G);
  const __m128i RG_hi = _mm_unpackhi_epi16(*R, *G);
  const __m128i GB_lo = _mm_unpacklo_epi16(*G, *B);
  const __m128i GB_hi = _mm_unpackhi_epi16(*G, *B);
  TRANSFORM(RG_lo, RG_hi, GB_lo, GB_hi, kRG_y, kGB_y, kHALF_Y, YUV_FIX, *Y);
}

static WEBP_INLINE void ConvertRGBToUV(const __m128i* const R,
                                       const __m128i* const G,
                                       const __m128i* const B,
                                       __m128i* const U, __m128i* const V) {
  const __m128i kRG_u = MK_CST_16(-9719, -19081);
  const __m128i kGB_u = MK_CST_16(0, 28800);
  const __m128i kRG_v = MK_CST_16(28800, 0);
  const __m128i kGB_v = MK_CST_16(-24116, -4684);
  const __m128i kHALF_UV = _mm_set1_epi32(((128 << YUV_FIX) + YUV_HALF) << 2);

  const __m128i RG_lo = _mm_unpacklo_epi16(*R, *G);
  const __m128i RG_hi = _mm_unpackhi_epi16(*R, *G);
  const __m128i GB_lo = _mm_unpacklo_epi16(*G, *B);
  const __m128i GB_hi = _mm_unpackhi_epi16(*G, *B);
  TRANSFORM(RG_lo, RG_hi, GB_lo, GB_hi, kRG_u, kGB_u,
            kHALF_UV, YUV_FIX + 2, *U);
  TRANSFORM(RG_lo, RG_hi, GB_lo, GB_hi, kRG_v, kGB_v,
            kHALF_UV, YUV_FIX + 2, *V);
}

#undef MK_CST_16
#undef TRANSFORM

static void ConvertRGB24ToY(const uint8_t* rgb, uint8_t* y, int width) {
  const int max_width = width & ~15;
  int i;
  for (i = 0; i < max_width; i += 16, rgb += 3 * 16) {
    __m128i r, g, b, Y0, Y1;
    RGB24PackedToPlanar(rgb + 0 * 8, &r, &g, &b, 0);
    ConvertRGBToY(&r, &g, &b, &Y0);
    RGB24PackedToPlanar(rgb + 3 * 8, &r, &g, &b, 0);
    ConvertRGBToY(&r, &g, &b, &Y1);
    STORE_16(_mm_packus_epi16(Y0, Y1), y + i);
  }
  for (; i < width; ++i, rgb += 3) {   // left-over
    y[i] = VP8RGBToY(rgb[0], rgb[1], rgb[2], YUV_HALF);
  }
}

static void ConvertBGR24ToY(const uint8_t* bgr, uint8_t* y, int width) {
  int i;
  const int max_width = width & ~15;
  for (i = 0; i < max_width; i += 16, bgr += 3 * 16) {
    __m128i r, g, b, Y0, Y1;
    RGB24PackedToPlanar(bgr + 0 * 8, &r, &g, &b, 1);
    ConvertRGBToY(&r, &g, &b, &Y0);
    RGB24PackedToPlanar(bgr + 3 * 8, &r, &g, &b, 1);
    ConvertRGBToY(&r, &g, &b, &Y1);
    STORE_16(_mm_packus_epi16(Y0, Y1), y + i);
  }
  for (; i < width; ++i, bgr += 3) {  // left-over
    y[i] = VP8RGBToY(bgr[2], bgr[1], bgr[0], YUV_HALF);
  }
}

static void ConvertARGBToY(const uint32_t* argb, uint8_t* y, int width) {
  const int max_width = width & ~15;
  int i;
  for (i = 0; i < max_width; i += 16) {
    __m128i r, g, b, Y0, Y1;
    RGB32PackedToPlanar(&argb[i + 0], &r, &g, &b);
    ConvertRGBToY(&r, &g, &b, &Y0);
    RGB32PackedToPlanar(&argb[i + 8], &r, &g, &b);
    ConvertRGBToY(&r, &g, &b, &Y1);
    STORE_16(_mm_packus_epi16(Y0, Y1), y + i);
  }
  for (; i < width; ++i) {   // left-over
    const uint32_t p = argb[i];
    y[i] = VP8RGBToY((p >> 16) & 0xff, (p >> 8) & 0xff, (p >>  0) & 0xff,
                     YUV_HALF);
  }
}

// Horizontal add (doubled) of two 16b values, result is 16b.
// in: A | B | C | D | ... -> out: 2*(A+B) | 2*(C+D) | ...
static void HorizontalAddPack(const __m128i* const A, const __m128i* const B,
                              __m128i* const out) {
  const __m128i k2 = _mm_set1_epi16(2);
  const __m128i C = _mm_madd_epi16(*A, k2);
  const __m128i D = _mm_madd_epi16(*B, k2);
  *out = _mm_packs_epi32(C, D);
}

static void ConvertARGBToUV(const uint32_t* argb, uint8_t* u, uint8_t* v,
                            int src_width, int do_store) {
  const int max_width = src_width & ~31;
  int i;
  for (i = 0; i < max_width; i += 32, u += 16, v += 16) {
    __m128i r0, g0, b0, r1, g1, b1, U0, V0, U1, V1;
    RGB32PackedToPlanar(&argb[i +  0], &r0, &g0, &b0);
    RGB32PackedToPlanar(&argb[i +  8], &r1, &g1, &b1);
    HorizontalAddPack(&r0, &r1, &r0);
    HorizontalAddPack(&g0, &g1, &g0);
    HorizontalAddPack(&b0, &b1, &b0);
    ConvertRGBToUV(&r0, &g0, &b0, &U0, &V0);

    RGB32PackedToPlanar(&argb[i + 16], &r0, &g0, &b0);
    RGB32PackedToPlanar(&argb[i + 24], &r1, &g1, &b1);
    HorizontalAddPack(&r0, &r1, &r0);
    HorizontalAddPack(&g0, &g1, &g0);
    HorizontalAddPack(&b0, &b1, &b0);
    ConvertRGBToUV(&r0, &g0, &b0, &U1, &V1);

    U0 = _mm_packus_epi16(U0, U1);
    V0 = _mm_packus_epi16(V0, V1);
    if (!do_store) {
      const __m128i prev_u = LOAD_16(u);
      const __m128i prev_v = LOAD_16(v);
      U0 = _mm_avg_epu8(U0, prev_u);
      V0 = _mm_avg_epu8(V0, prev_v);
    }
    STORE_16(U0, u);
    STORE_16(V0, v);
  }
  if (i < src_width) {  // left-over
    WebPConvertARGBToUV_C(argb + i, u, v, src_width - i, do_store);
  }
}

// Convert 16 packed ARGB 16b-values to r[], g[], b[]
static WEBP_INLINE void RGBA32PackedToPlanar_16b(const uint16_t* const rgbx,
                                                 __m128i* const r,
                                                 __m128i* const g,
                                                 __m128i* const b) {
  const __m128i in0 = LOAD_16(rgbx +  0);  // r0 | g0 | b0 |x| r1 | g1 | b1 |x
  const __m128i in1 = LOAD_16(rgbx +  8);  // r2 | g2 | b2 |x| r3 | g3 | b3 |x
  const __m128i in2 = LOAD_16(rgbx + 16);  // r4 | ...
  const __m128i in3 = LOAD_16(rgbx + 24);  // r6 | ...
  // column-wise transpose
  const __m128i A0 = _mm_unpacklo_epi16(in0, in1);
  const __m128i A1 = _mm_unpackhi_epi16(in0, in1);
  const __m128i A2 = _mm_unpacklo_epi16(in2, in3);
  const __m128i A3 = _mm_unpackhi_epi16(in2, in3);
  const __m128i B0 = _mm_unpacklo_epi16(A0, A1);  // r0 r1 r2 r3 | g0 g1 ..
  const __m128i B1 = _mm_unpackhi_epi16(A0, A1);  // b0 b1 b2 b3 | x x x x
  const __m128i B2 = _mm_unpacklo_epi16(A2, A3);  // r4 r5 r6 r7 | g4 g5 ..
  const __m128i B3 = _mm_unpackhi_epi16(A2, A3);  // b4 b5 b6 b7 | x x x x
  *r = _mm_unpacklo_epi64(B0, B2);
  *g = _mm_unpackhi_epi64(B0, B2);
  *b = _mm_unpacklo_epi64(B1, B3);
}

static void ConvertRGBA32ToUV(const uint16_t* rgb,
                              uint8_t* u, uint8_t* v, int width) {
  const int max_width = width & ~15;
  const uint16_t* const last_rgb = rgb + 4 * max_width;
  while (rgb < last_rgb) {
    __m128i r, g, b, U0, V0, U1, V1;
    RGBA32PackedToPlanar_16b(rgb +  0, &r, &g, &b);
    ConvertRGBToUV(&r, &g, &b, &U0, &V0);
    RGBA32PackedToPlanar_16b(rgb + 32, &r, &g, &b);
    ConvertRGBToUV(&r, &g, &b, &U1, &V1);
    STORE_16(_mm_packus_epi16(U0, U1), u);
    STORE_16(_mm_packus_epi16(V0, V1), v);
    u += 16;
    v += 16;
    rgb += 2 * 32;
  }
  if (max_width < width) {  // left-over
    WebPConvertRGBA32ToUV_C(rgb, u, v, width - max_width);
  }
}

//------------------------------------------------------------------------------

extern void WebPInitConvertARGBToYUVSSE2(void);

WEBP_TSAN_IGNORE_FUNCTION void WebPInitConvertARGBToYUVSSE2(void) {
  WebPConvertARGBToY = ConvertARGBToY;
  WebPConvertARGBToUV = ConvertARGBToUV;

  WebPConvertRGB24ToY = ConvertRGB24ToY;
  WebPConvertBGR24ToY = ConvertBGR24ToY;

  WebPConvertRGBA32ToUV = ConvertRGBA32ToUV;
}

#else  // !WEBP_USE_SSE2

WEBP_DSP_INIT_STUB(WebPInitSamplersSSE2)
WEBP_DSP_INIT_STUB(WebPInitConvertARGBToYUVSSE2)

#endif  // WEBP_USE_SSE2
