// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// NEON version of YUV to RGB upsampling functions.
//
// Author: mans@mansr.com (Mans Rullgard)
// Based on SSE code by: somnath@google.com (Somnath Banerjee)

#include "./dsp.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#if defined(WEBP_USE_NEON)

#include <assert.h>
#include <arm_neon.h>
#include <string.h>
#include "./yuv.h"

#ifdef FANCY_UPSAMPLING

// Loads 9 pixels each from rows r1 and r2 and generates 16 pixels.
#define UPSAMPLE_16PIXELS(r1, r2, out) {                                \
  uint8x8_t a = vld1_u8(r1);                                            \
  uint8x8_t b = vld1_u8(r1 + 1);                                        \
  uint8x8_t c = vld1_u8(r2);                                            \
  uint8x8_t d = vld1_u8(r2 + 1);                                        \
                                                                        \
  uint16x8_t al = vshll_n_u8(a, 1);                                     \
  uint16x8_t bl = vshll_n_u8(b, 1);                                     \
  uint16x8_t cl = vshll_n_u8(c, 1);                                     \
  uint16x8_t dl = vshll_n_u8(d, 1);                                     \
                                                                        \
  uint8x8_t diag1, diag2;                                               \
  uint16x8_t sl;                                                        \
                                                                        \
  /* a + b + c + d */                                                   \
  sl = vaddl_u8(a,  b);                                                 \
  sl = vaddw_u8(sl, c);                                                 \
  sl = vaddw_u8(sl, d);                                                 \
                                                                        \
  al = vaddq_u16(sl, al); /* 3a +  b +  c +  d */                       \
  bl = vaddq_u16(sl, bl); /*  a + 3b +  c +  d */                       \
                                                                        \
  al = vaddq_u16(al, dl); /* 3a +  b +  c + 3d */                       \
  bl = vaddq_u16(bl, cl); /*  a + 3b + 3c +  d */                       \
                                                                        \
  diag2 = vshrn_n_u16(al, 3);                                           \
  diag1 = vshrn_n_u16(bl, 3);                                           \
                                                                        \
  a = vrhadd_u8(a, diag1);                                              \
  b = vrhadd_u8(b, diag2);                                              \
  c = vrhadd_u8(c, diag2);                                              \
  d = vrhadd_u8(d, diag1);                                              \
                                                                        \
  vst2_u8(out,      (uint8x8x2_t){{ a, b }});                           \
  vst2_u8(out + 32, (uint8x8x2_t){{ c, d }});                           \
}

// Turn the macro into a function for reducing code-size when non-critical
static void Upsample16Pixels(const uint8_t *r1, const uint8_t *r2,
                             uint8_t *out) {
  UPSAMPLE_16PIXELS(r1, r2, out);
}

#define UPSAMPLE_LAST_BLOCK(tb, bb, num_pixels, out) {                  \
  uint8_t r1[9], r2[9];                                                 \
  memcpy(r1, (tb), (num_pixels));                                       \
  memcpy(r2, (bb), (num_pixels));                                       \
  /* replicate last byte */                                             \
  memset(r1 + (num_pixels), r1[(num_pixels) - 1], 9 - (num_pixels));    \
  memset(r2 + (num_pixels), r2[(num_pixels) - 1], 9 - (num_pixels));    \
  Upsample16Pixels(r1, r2, out);                                        \
}

#define CY  (76276  / 4)
#define CVR (104560 / 4)
#define CUG (25660  / 4)
#define CVG (53256  / 4)
#define CUB (132156 / 4)

static inline int clip(int x) {
  return x < 0 ? 0 : x > 255 ? 255 : x;
}

static inline int yuv2r(int y, int u, int v) {
  return clip((CY * y + CVR * v + 8192) >> 14);
}

static inline int yuv2g(int y, int u, int v) {
  return clip((CY * y - CUG * u - CVG * v + 8192) >> 14);
}

static inline int yuv2b(int y, int u, int v) {
  return clip((CY * y + CUB * u + 8192) >> 14);
}

static inline void yuv2rgb(uint8_t *out, int y, int u, int v) {
  out[0] = yuv2r(y, u, v);
  out[1] = yuv2g(y, u, v);
  out[2] = yuv2b(y, u, v);
}

static inline void yuv2bgr(uint8_t *out, int y, int u, int v) {
  out[0] = yuv2b(y, u, v);
  out[1] = yuv2g(y, u, v);
  out[2] = yuv2r(y, u, v);
}

static inline void yuv2rgba(uint8_t *out, int y, int u, int v) {
  yuv2rgb(out, y, u, v);
  out[3] = 255;
}

static inline void yuv2bgra(uint8_t *out, int y, int u, int v) {
  yuv2bgr(out, y, u, v);
  out[3] = 255;
}

static const int16_t coef[] = { CY, CVR, CUG, CVG };

#define v255 vmov_n_u8(255)

#define STR_rgb(out, r, g, b)  vst3_u8(out, (uint8x8x3_t){{ r, g, b }})
#define STR_bgr(out, r, g, b)  vst3_u8(out, (uint8x8x3_t){{ b, g, r }})
#define STR_rgba(out, r, g, b) vst4_u8(out, (uint8x8x4_t){{ r, g, b, v255 }})
#define STR_bgra(out, r, g, b) vst4_u8(out, (uint8x8x4_t){{ b, g, r, v255 }})

#define CONVERT8(FMT, XSTEP, N, src_y, src_uv, out, cur_x) {            \
  int i;                                                                \
  for (i = 0; i < N; i += 8) {                                          \
    int off = ((cur_x) + i) * XSTEP;                                    \
    uint8x8_t y  = vld1_u8(src_y + (cur_x)  + i);                       \
    uint8x8_t u  = vld1_u8((src_uv) + i);                               \
    uint8x8_t v  = vld1_u8((src_uv) + i + 16);                          \
    int16x8_t yy = vreinterpretq_s16_u16(vsubl_u8(y, u16));             \
    int16x8_t uu = vreinterpretq_s16_u16(vsubl_u8(u, u128));            \
    int16x8_t vv = vreinterpretq_s16_u16(vsubl_u8(v, u128));            \
    int32x4_t yl = vmull_lane_s16(vget_low_s16(yy),  cf16, 0);          \
    int32x4_t yh = vmull_lane_s16(vget_high_s16(yy), cf16, 0);          \
    int32x4_t rl = vmlal_lane_s16(yl, vget_low_s16(vv),  cf16, 1);      \
    int32x4_t rh = vmlal_lane_s16(yh, vget_high_s16(vv), cf16, 1);      \
    int32x4_t gl = vmlsl_lane_s16(yl, vget_low_s16(uu),  cf16, 2);      \
    int32x4_t gh = vmlsl_lane_s16(yh, vget_high_s16(uu), cf16, 2);      \
    int32x4_t bl = vmovl_s16(vget_low_s16(uu));                         \
    int32x4_t bh = vmovl_s16(vget_high_s16(uu));                        \
    gl = vmlsl_lane_s16(gl, vget_low_s16(vv),  cf16, 3);                \
    gh = vmlsl_lane_s16(gh, vget_high_s16(vv), cf16, 3);                \
    yl = vmlaq_lane_s32(yl, bl, cf32, 0);                               \
    yh = vmlaq_lane_s32(yh, bh, cf32, 0);                               \
    y = vqmovun_s16(vcombine_s16(vrshrn_n_s32(rl, 14),                  \
                                 vrshrn_n_s32(rh, 14)));                \
    u = vqmovun_s16(vcombine_s16(vrshrn_n_s32(gl, 14),                  \
                                 vrshrn_n_s32(gh, 14)));                \
    v = vqmovun_s16(vcombine_s16(vrshrn_n_s32(yl, 14),                  \
                                 vrshrn_n_s32(yh, 14)));                \
    STR_ ## FMT(out + off, y, u, v);                                    \
  }                                                                     \
}

#define CONVERT1(FUNC, XSTEP, N, src_y, src_uv, rgb, cur_x) {           \
  int i;                                                                \
  for (i = 0; i < N; i++) {                                             \
    int off = ((cur_x) + i) * XSTEP;                                    \
    int y = src_y[(cur_x) + i] - 16;                                    \
    int u = (src_uv)[i]        - 128;                                   \
    int v = (src_uv)[i + 16]   - 128;                                   \
    FUNC(rgb + off, y, u, v);                                           \
  }                                                                     \
}

#define CONVERT2RGB_8(FMT, XSTEP, top_y, bottom_y, uv,                  \
                      top_dst, bottom_dst, cur_x, len) {                \
  if (top_y) {                                                          \
    CONVERT8(FMT, XSTEP, len, top_y, uv, top_dst, cur_x)                \
  }                                                                     \
  if (bottom_y) {                                                       \
    CONVERT8(FMT, XSTEP, len, bottom_y, (uv) + 32, bottom_dst, cur_x)   \
  }                                                                     \
}

#define CONVERT2RGB_1(FUNC, XSTEP, top_y, bottom_y, uv,                 \
                      top_dst, bottom_dst, cur_x, len) {                \
  if (top_y) {                                                          \
    CONVERT1(FUNC, XSTEP, len, top_y, uv, top_dst, cur_x);              \
  }                                                                     \
  if (bottom_y) {                                                       \
    CONVERT1(FUNC, XSTEP, len, bottom_y, (uv) + 32, bottom_dst, cur_x); \
  }                                                                     \
}

#define NEON_UPSAMPLE_FUNC(FUNC_NAME, FMT, XSTEP)                       \
static void FUNC_NAME(const uint8_t *top_y, const uint8_t *bottom_y,    \
                      const uint8_t *top_u, const uint8_t *top_v,       \
                      const uint8_t *cur_u, const uint8_t *cur_v,       \
                      uint8_t *top_dst, uint8_t *bottom_dst, int len) { \
  int b;                                                                \
  /* 16 byte aligned array to cache reconstructed u and v */            \
  uint8_t uv_buf[2 * 32 + 15];                                          \
  uint8_t *const r_uv = (uint8_t*)((uintptr_t)(uv_buf + 15) & ~15);     \
  const int uv_len = (len + 1) >> 1;                                    \
  /* 9 pixels must be read-able for each block */                       \
  const int num_blocks = (uv_len - 1) >> 3;                             \
  const int leftover = uv_len - num_blocks * 8;                         \
  const int last_pos = 1 + 16 * num_blocks;                             \
                                                                        \
  const int u_diag = ((top_u[0] + cur_u[0]) >> 1) + 1;                  \
  const int v_diag = ((top_v[0] + cur_v[0]) >> 1) + 1;                  \
                                                                        \
  const int16x4_t cf16 = vld1_s16(coef);                                \
  const int32x2_t cf32 = vmov_n_s32(CUB);                               \
  const uint8x8_t u16  = vmov_n_u8(16);                                 \
  const uint8x8_t u128 = vmov_n_u8(128);                                \
                                                                        \
  /* Treat the first pixel in regular way */                            \
  if (top_y) {                                                          \
    const int u0 = (top_u[0] + u_diag) >> 1;                            \
    const int v0 = (top_v[0] + v_diag) >> 1;                            \
    yuv2 ## FMT(top_dst, top_y[0] - 16, u0 - 128, v0 - 128);            \
  }                                                                     \
  if (bottom_y) {                                                       \
    const int u0 = (cur_u[0] + u_diag) >> 1;                            \
    const int v0 = (cur_v[0] + v_diag) >> 1;                            \
    yuv2 ## FMT(bottom_dst, bottom_y[0] - 16, u0 - 128, v0 - 128);      \
  }                                                                     \
                                                                        \
  for (b = 0; b < num_blocks; ++b) {                                    \
    UPSAMPLE_16PIXELS(top_u, cur_u, r_uv);                              \
    UPSAMPLE_16PIXELS(top_v, cur_v, r_uv + 16);                         \
    CONVERT2RGB_8(FMT, XSTEP, top_y, bottom_y, r_uv,                    \
                  top_dst, bottom_dst, 16 * b + 1, 16);                 \
    top_u += 8;                                                         \
    cur_u += 8;                                                         \
    top_v += 8;                                                         \
    cur_v += 8;                                                         \
  }                                                                     \
                                                                        \
  UPSAMPLE_LAST_BLOCK(top_u, cur_u, leftover, r_uv);                    \
  UPSAMPLE_LAST_BLOCK(top_v, cur_v, leftover, r_uv + 16);               \
  CONVERT2RGB_1(yuv2 ## FMT, XSTEP, top_y, bottom_y, r_uv,              \
                top_dst, bottom_dst, last_pos, len - last_pos);         \
}

// NEON variants of the fancy upsampler.
NEON_UPSAMPLE_FUNC(UpsampleRgbLinePairNEON,  rgb,  3)
NEON_UPSAMPLE_FUNC(UpsampleBgrLinePairNEON,  bgr,  3)
NEON_UPSAMPLE_FUNC(UpsampleRgbaLinePairNEON, rgba, 4)
NEON_UPSAMPLE_FUNC(UpsampleBgraLinePairNEON, bgra, 4)

#endif  // FANCY_UPSAMPLING

#endif   // WEBP_USE_NEON

//------------------------------------------------------------------------------

extern WebPUpsampleLinePairFunc WebPUpsamplers[/* MODE_LAST */];

void WebPInitUpsamplersNEON(void) {
#if defined(WEBP_USE_NEON)
  WebPUpsamplers[MODE_RGB]  = UpsampleRgbLinePairNEON;
  WebPUpsamplers[MODE_RGBA] = UpsampleRgbaLinePairNEON;
  WebPUpsamplers[MODE_BGR]  = UpsampleBgrLinePairNEON;
  WebPUpsamplers[MODE_BGRA] = UpsampleBgraLinePairNEON;
#endif   // WEBP_USE_NEON
}

void WebPInitPremultiplyNEON(void) {
#if defined(WEBP_USE_NEON)
  WebPUpsamplers[MODE_rgbA] = UpsampleRgbaLinePairNEON;
  WebPUpsamplers[MODE_bgrA] = UpsampleBgraLinePairNEON;
#endif   // WEBP_USE_NEON
}

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif
