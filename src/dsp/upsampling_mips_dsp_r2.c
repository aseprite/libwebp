// Copyright 2014 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
// YUV to RGB upsampling functions.
//
// Author(s): Branimir Vasic (branimir.vasic@imgtec.com)
//            Djordje Pesut  (djordje.pesut@imgtec.com)

#include "./dsp.h"

#if defined(WEBP_USE_MIPS_DSP_R2)

#include "./yuv.h"
#include <assert.h>

#if !defined(WEBP_YUV_USE_TABLE)

#define VP8YUV_FUNCTIONS(Y, U, V) do {                                         \
    t1 = kYScale * Y;                                                          \
    t2 = kVToR * V;                                                            \
    t3 = kUToG * U;                                                            \
    t4 = kUToB * U;                                                            \
    t5 = kVToG * V;                                                            \
    t2 = t1 + t2;                                                              \
    t3 = t1 - t3;                                                              \
    t4 = t1 + t4;                                                              \
    t2 = t2 + kRCst;                                                           \
    t3 = t3 - t5 + kGCst;                                                      \
    t4 = t4 + kBCst;                                                           \
    __asm__ volatile (                                                         \
      "shll_s.w         %[t2],      %[t2],        9                 \n\t"      \
      "shll_s.w         %[t3],      %[t3],        9                 \n\t"      \
      "shll_s.w         %[t4],      %[t4],        9                 \n\t"      \
      "precrqu_s.qb.ph  %[t2],      %[t2],        $zero             \n\t"      \
      "precrqu_s.qb.ph  %[t3],      %[t3],        $zero             \n\t"      \
      "precrqu_s.qb.ph  %[t4],      %[t4],        $zero             \n\t"      \
      "srl              %[t2],      %[t2],        24                \n\t"      \
      "srl              %[t3],      %[t3],        24                \n\t"      \
      "srl              %[t4],      %[t4],        24                \n\t"      \
      : [t2]"+r"(t2), [t3]"+r"(t3), [t4]"+r"(t4)                               \
      :                                                                        \
    );                                                                         \
  } while (0)

static WEBP_INLINE void VP8YuvToRgbMips(int y, int u, int v,
                                        uint8_t* const rgb) {
  int t1, t2, t3, t4, t5;
  VP8YUV_FUNCTIONS(y,u,v);
  rgb[0] = t2;
  rgb[1] = t3;
  rgb[2] = t4;
}
static WEBP_INLINE void VP8YuvToBgrMips(int y, int u, int v,
                                        uint8_t* const bgr) {
  int t1, t2, t3, t4, t5;
  VP8YUV_FUNCTIONS(y,u,v);
  bgr[0] = t4;
  bgr[1] = t3;
  bgr[2] = t2;
}
static WEBP_INLINE void VP8YuvToRgb565Mips(int y, int u, int v,
                                           uint8_t* const rgb) {
  int t1, t2, t3, t4, t5;
  VP8YUV_FUNCTIONS(y,u,v);
  {
    const int rg = (t2 & 0xf8) | (t3 >> 5);
    const int gb = ((t3 << 3) & 0xe0) | (t4 >> 3);
#ifdef WEBP_SWAP_16BIT_CSP
    rgb[0] = gb;
    rgb[1] = rg;
#else
    rgb[0] = rg;
    rgb[1] = gb;
#endif
  }
}
static WEBP_INLINE void VP8YuvToRgba4444Mips(int y, int u, int v,
                                             uint8_t* const argb) {
  int t1, t2, t3, t4, t5;
  VP8YUV_FUNCTIONS(y,u,v);
  {
    const int rg = (t2 & 0xf0) | (t3 >> 4);
    const int ba = (t4 & 0xf0) | 0x0f;     // overwrite the lower 4 bits
#ifdef WEBP_SWAP_16BIT_CSP
    argb[0] = ba;
    argb[1] = rg;
#else
    argb[0] = rg;
    argb[1] = ba;
#endif
   }
}
#endif  // WEBP_YUV_USE_TABLE

//-----------------------------------------------------------------------------
// Alpha handling variants

static WEBP_INLINE void VP8YuvToArgbMips(uint8_t y, uint8_t u, uint8_t v,
                                         uint8_t* const argb) {
  int t1, t2, t3, t4, t5;
  VP8YUV_FUNCTIONS(y,u,v);
  argb[0] = 0xff;
  argb[1] =(t2);
  argb[2] =(t3);
  argb[3] =(t4);
}
static WEBP_INLINE void VP8YuvToBgraMips(uint8_t y, uint8_t u, uint8_t v,
                                         uint8_t* const bgra) {
  int t1, t2, t3, t4, t5;
  VP8YUV_FUNCTIONS(y,u,v);
  bgra[0] =(t4);
  bgra[1] =(t3);
  bgra[2] =(t2);
  bgra[3] = 0xff;
}
static WEBP_INLINE void VP8YuvToRgbaMips(uint8_t y, uint8_t u, uint8_t v,
                                         uint8_t* const rgba) {
  int t1, t2, t3, t4, t5;
  VP8YUV_FUNCTIONS(y,u,v);
  rgba[0] =(t2);
  rgba[1] =(t3);
  rgba[2] =(t4);
  rgba[3] = 0xff;
}

//------------------------------------------------------------------------------
// Fancy upsampler

#ifdef FANCY_UPSAMPLING

// Given samples laid out in a square as:
//  [a b]
//  [c d]
// we interpolate u/v as:
//  ([9*a + 3*b + 3*c +   d    3*a + 9*b + 3*c +   d] + [8 8]) / 16
//  ([3*a +   b + 9*c + 3*d      a + 3*b + 3*c + 9*d]   [8 8]) / 16

// We process u and v together stashed into 32bit (16bit each).
#define LOAD_UV(u, v) ((u) | ((v) << 16))

#define UPSAMPLE_FUNC(FUNC_NAME, FUNC, XSTEP)                                  \
static void FUNC_NAME(const uint8_t* top_y, const uint8_t* bottom_y,           \
                      const uint8_t* top_u, const uint8_t* top_v,              \
                      const uint8_t* cur_u, const uint8_t* cur_v,              \
                      uint8_t* top_dst, uint8_t* bottom_dst, int len) {        \
  int x;                                                                       \
  const int last_pixel_pair = (len - 1) >> 1;                                  \
  uint32_t tl_uv = LOAD_UV(top_u[0], top_v[0]);   /* top-left sample */        \
  uint32_t l_uv  = LOAD_UV(cur_u[0], cur_v[0]);   /* left-sample */            \
  assert(top_y != NULL);                                                       \
  {                                                                            \
    const uint32_t uv0 = (3 * tl_uv + l_uv + 0x00020002u) >> 2;                \
    FUNC(top_y[0], uv0 & 0xff, (uv0 >> 16), top_dst);                          \
  }                                                                            \
  if (bottom_y != NULL) {                                                      \
    const uint32_t uv0 = (3 * l_uv + tl_uv + 0x00020002u) >> 2;                \
    FUNC(bottom_y[0], uv0 & 0xff, (uv0 >> 16), bottom_dst);                    \
  }                                                                            \
  for (x = 1; x <= last_pixel_pair; ++x) {                                     \
    const uint32_t t_uv = LOAD_UV(top_u[x], top_v[x]);  /* top sample */       \
    const uint32_t uv   = LOAD_UV(cur_u[x], cur_v[x]);  /* sample */           \
    /* precompute invariant values associated with first and second diagonals*/\
    const uint32_t avg = tl_uv + t_uv + l_uv + uv + 0x00080008u;               \
    const uint32_t diag_12 = (avg + 2 * (t_uv + l_uv)) >> 3;                   \
    const uint32_t diag_03 = (avg + 2 * (tl_uv + uv)) >> 3;                    \
    {                                                                          \
      const uint32_t uv0 = (diag_12 + tl_uv) >> 1;                             \
      const uint32_t uv1 = (diag_03 + t_uv) >> 1;                              \
      FUNC(top_y[2 * x - 1], uv0 & 0xff, (uv0 >> 16),                          \
           top_dst + (2 * x - 1) * XSTEP);                                     \
      FUNC(top_y[2 * x - 0], uv1 & 0xff, (uv1 >> 16),                          \
           top_dst + (2 * x - 0) * XSTEP);                                     \
    }                                                                          \
    if (bottom_y != NULL) {                                                    \
      const uint32_t uv0 = (diag_03 + l_uv) >> 1;                              \
      const uint32_t uv1 = (diag_12 + uv) >> 1;                                \
      FUNC(bottom_y[2 * x - 1], uv0 & 0xff, (uv0 >> 16),                       \
           bottom_dst + (2 * x - 1) * XSTEP);                                  \
      FUNC(bottom_y[2 * x + 0], uv1 & 0xff, (uv1 >> 16),                       \
           bottom_dst + (2 * x + 0) * XSTEP);                                  \
    }                                                                          \
    tl_uv = t_uv;                                                              \
    l_uv = uv;                                                                 \
  }                                                                            \
  if (!(len & 1)) {                                                            \
    {                                                                          \
      const uint32_t uv0 = (3 * tl_uv + l_uv + 0x00020002u) >> 2;              \
      FUNC(top_y[len - 1], uv0 & 0xff, (uv0 >> 16),                            \
           top_dst + (len - 1) * XSTEP);                                       \
    }                                                                          \
    if (bottom_y != NULL) {                                                    \
      const uint32_t uv0 = (3 * l_uv + tl_uv + 0x00020002u) >> 2;              \
      FUNC(bottom_y[len - 1], uv0 & 0xff, (uv0 >> 16),                         \
           bottom_dst + (len - 1) * XSTEP);                                    \
    }                                                                          \
  }                                                                            \
}

// All variants implemented.
UPSAMPLE_FUNC(UpsampleRgbLinePair,  VP8YuvToRgbMips,  3)
UPSAMPLE_FUNC(UpsampleBgrLinePair,  VP8YuvToBgrMips,  3)
UPSAMPLE_FUNC(UpsampleRgbaLinePair, VP8YuvToRgbaMips, 4)
UPSAMPLE_FUNC(UpsampleBgraLinePair, VP8YuvToBgraMips, 4)
UPSAMPLE_FUNC(UpsampleArgbLinePair, VP8YuvToArgbMips, 4)
UPSAMPLE_FUNC(UpsampleRgba4444LinePair, VP8YuvToRgba4444Mips, 2)
UPSAMPLE_FUNC(UpsampleRgb565LinePair,  VP8YuvToRgb565Mips,  2)

#undef LOAD_UV
#undef UPSAMPLE_FUNC

#endif  // FANCY_UPSAMPLING

//------------------------------------------------------------------------------
// YUV444 converter

#define YUV444_FUNC(FUNC_NAME, FUNC, XSTEP)                                    \
static void FUNC_NAME(const uint8_t* y, const uint8_t* u, const uint8_t* v,    \
                      uint8_t* dst, int len) {                                 \
  int i;                                                                       \
  for (i = 0; i < len; ++i) FUNC(y[i], u[i], v[i], &dst[i * XSTEP]);           \
}

YUV444_FUNC(Yuv444ToRgb,      VP8YuvToRgbMips,  3)
YUV444_FUNC(Yuv444ToBgr,      VP8YuvToBgrMips,  3)
YUV444_FUNC(Yuv444ToRgba,     VP8YuvToRgbaMips, 4)
YUV444_FUNC(Yuv444ToBgra,     VP8YuvToBgraMips, 4)
YUV444_FUNC(Yuv444ToArgb,     VP8YuvToArgbMips, 4)
YUV444_FUNC(Yuv444ToRgba4444, VP8YuvToRgba4444Mips, 2)
YUV444_FUNC(Yuv444ToRgb565,   VP8YuvToRgb565Mips, 2)

#undef YUV444_FUNC

#endif  // WEBP_USE_MIPS_DSP_R2

//------------------------------------------------------------------------------

extern void WebPInitUpsamplersMIPSdspR2(void);

void WebPInitUpsamplersMIPSdspR2(void) {
#if defined(WEBP_USE_MIPS_DSP_R2)
#ifdef FANCY_UPSAMPLING
  WebPUpsamplers[MODE_RGB]       = UpsampleRgbLinePair;
  WebPUpsamplers[MODE_RGBA]      = UpsampleRgbaLinePair;
  WebPUpsamplers[MODE_BGR]       = UpsampleBgrLinePair;
  WebPUpsamplers[MODE_BGRA]      = UpsampleBgraLinePair;
  WebPUpsamplers[MODE_ARGB]      = UpsampleArgbLinePair;
  WebPUpsamplers[MODE_RGBA_4444] = UpsampleRgba4444LinePair;
  WebPUpsamplers[MODE_RGB_565]   = UpsampleRgb565LinePair;
  WebPUpsamplers[MODE_rgbA]      = UpsampleRgbaLinePair;
  WebPUpsamplers[MODE_bgrA]      = UpsampleBgraLinePair;
  WebPUpsamplers[MODE_Argb]      = UpsampleArgbLinePair;
  WebPUpsamplers[MODE_rgbA_4444] = UpsampleRgba4444LinePair;
#endif  // FANCY_UPSAMPLING
  WebPYUV444Converters[MODE_RGB]       = Yuv444ToRgb;
  WebPYUV444Converters[MODE_RGBA]      = Yuv444ToRgba;
  WebPYUV444Converters[MODE_BGR]       = Yuv444ToBgr;
  WebPYUV444Converters[MODE_BGRA]      = Yuv444ToBgra;
  WebPYUV444Converters[MODE_ARGB]      = Yuv444ToArgb;
  WebPYUV444Converters[MODE_RGBA_4444] = Yuv444ToRgba4444;
  WebPYUV444Converters[MODE_RGB_565]   = Yuv444ToRgb565;
  WebPYUV444Converters[MODE_rgbA]      = Yuv444ToRgba;
  WebPYUV444Converters[MODE_bgrA]      = Yuv444ToBgra;
  WebPYUV444Converters[MODE_Argb]      = Yuv444ToArgb;
  WebPYUV444Converters[MODE_rgbA_4444] = Yuv444ToRgba4444;
#endif  // WEBP_USE_MIPS_DSP_R2
}

//------------------------------------------------------------------------------
