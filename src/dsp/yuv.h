// Copyright 2010 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// inline YUV<->RGB conversion function
//
// The exact naming is Y'CbCr, following the ITU-R BT.601 standard.
// More information at: http://en.wikipedia.org/wiki/YCbCr
// Y = 0.2569 * R + 0.5044 * G + 0.0979 * B + 16
// U = -0.1483 * R - 0.2911 * G + 0.4394 * B + 128
// V = 0.4394 * R - 0.3679 * G - 0.0715 * B + 128
// We use 16bit fixed point operations for RGB->YUV conversion.
//
// Author: Skal (pascal.massimino@gmail.com)

#ifndef WEBP_DSP_YUV_H_
#define WEBP_DSP_YUV_H_

#include "../dec/decode_vp8.h"

#if defined(WEBP_EXPERIMENTAL_FEATURES)
// Do NOT activate this feature for real compression. This is only experimental!
// This flag is for comparison purpose against JPEG's "YUVj" natural colorspace.
// This colorspace is close to Rec.601's Y'CbCr model with the notable
// difference of allowing larger range for luma/chroma.
// See http://en.wikipedia.org/wiki/YCbCr#JPEG_conversion paragraph, and its
// difference with http://en.wikipedia.org/wiki/YCbCr#ITU-R_BT.601_conversion
// #define USE_YUVj
#endif

//------------------------------------------------------------------------------
// YUV -> RGB conversion

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

enum { YUV_FIX = 16,                // fixed-point precision
       YUV_RANGE_MIN = -227,        // min value of r/g/b output
       YUV_RANGE_MAX = 256 + 226    // max value of r/g/b output
};
extern int16_t VP8kVToR[256], VP8kUToB[256];
extern int32_t VP8kVToG[256], VP8kUToG[256];
extern uint8_t VP8kClip[YUV_RANGE_MAX - YUV_RANGE_MIN];
extern uint8_t VP8kClip4Bits[YUV_RANGE_MAX - YUV_RANGE_MIN];

static WEBP_INLINE void VP8YuvToRgb(uint8_t y, uint8_t u, uint8_t v,
                                    uint8_t* const rgb) {
  const int r_off = VP8kVToR[v];
  const int g_off = (VP8kVToG[v] + VP8kUToG[u]) >> YUV_FIX;
  const int b_off = VP8kUToB[u];
  rgb[0] = VP8kClip[y + r_off - YUV_RANGE_MIN];
  rgb[1] = VP8kClip[y + g_off - YUV_RANGE_MIN];
  rgb[2] = VP8kClip[y + b_off - YUV_RANGE_MIN];
}

static WEBP_INLINE void VP8YuvToRgb565(uint8_t y, uint8_t u, uint8_t v,
                                       uint8_t* const rgb) {
  const int r_off = VP8kVToR[v];
  const int g_off = (VP8kVToG[v] + VP8kUToG[u]) >> YUV_FIX;
  const int b_off = VP8kUToB[u];
  const uint8_t rg = ((VP8kClip[y + r_off - YUV_RANGE_MIN] & 0xf8) |
                      (VP8kClip[y + g_off - YUV_RANGE_MIN] >> 5));
  const uint8_t gb = (((VP8kClip[y + g_off - YUV_RANGE_MIN] << 3) & 0xe0) |
                      (VP8kClip[y + b_off - YUV_RANGE_MIN] >> 3));
#ifdef WEBP_SWAP_16BIT_CSP
  rgb[0] = gb;
  rgb[1] = rg;
#else
  rgb[0] = rg;
  rgb[1] = gb;
#endif
}

static WEBP_INLINE void VP8YuvToArgb(uint8_t y, uint8_t u, uint8_t v,
                                     uint8_t* const argb) {
  argb[0] = 0xff;
  VP8YuvToRgb(y, u, v, argb + 1);
}

static WEBP_INLINE void VP8YuvToRgba4444(uint8_t y, uint8_t u, uint8_t v,
                                         uint8_t* const argb) {
  const int r_off = VP8kVToR[v];
  const int g_off = (VP8kVToG[v] + VP8kUToG[u]) >> YUV_FIX;
  const int b_off = VP8kUToB[u];
  const uint8_t rg = ((VP8kClip4Bits[y + r_off - YUV_RANGE_MIN] << 4) |
                      VP8kClip4Bits[y + g_off - YUV_RANGE_MIN]);
  const uint8_t ba = (VP8kClip4Bits[y + b_off - YUV_RANGE_MIN] << 4) | 0x0f;
#ifdef WEBP_SWAP_16BIT_CSP
  argb[0] = ba;
  argb[1] = rg;
#else
  argb[0] = rg;
  argb[1] = ba;
#endif
}

static WEBP_INLINE void VP8YuvToBgr(uint8_t y, uint8_t u, uint8_t v,
                                    uint8_t* const bgr) {
  const int r_off = VP8kVToR[v];
  const int g_off = (VP8kVToG[v] + VP8kUToG[u]) >> YUV_FIX;
  const int b_off = VP8kUToB[u];
  bgr[0] = VP8kClip[y + b_off - YUV_RANGE_MIN];
  bgr[1] = VP8kClip[y + g_off - YUV_RANGE_MIN];
  bgr[2] = VP8kClip[y + r_off - YUV_RANGE_MIN];
}

static WEBP_INLINE void VP8YuvToBgra(uint8_t y, uint8_t u, uint8_t v,
                                     uint8_t* const bgra) {
  VP8YuvToBgr(y, u, v, bgra);
  bgra[3] = 0xff;
}

static WEBP_INLINE void VP8YuvToRgba(uint8_t y, uint8_t u, uint8_t v,
                                     uint8_t* const rgba) {
  VP8YuvToRgb(y, u, v, rgba);
  rgba[3] = 0xff;
}

// Must be called before everything, to initialize the tables.
void VP8YUVInit(void);

//------------------------------------------------------------------------------
// RGB -> YUV conversion

static WEBP_INLINE int VP8ClipUV(int v) {
  v = (v + (257 << (YUV_FIX + 2 - 1))) >> (YUV_FIX + 2);
  return ((v & ~0xff) == 0) ? v : (v < 0) ? 0 : 255;
}

#ifndef USE_YUVj

static WEBP_INLINE int VP8RGBToY(int r, int g, int b) {
  const int kRound = (1 << (YUV_FIX - 1)) + (16 << YUV_FIX);
  const int luma = 16839 * r + 33059 * g + 6420 * b;
  return (luma + kRound) >> YUV_FIX;  // no need to clip
}

static WEBP_INLINE int VP8RGBToU(int r, int g, int b) {
  const int u = -9719 * r - 19081 * g + 28800 * b;
  return VP8ClipUV(u);
}

static WEBP_INLINE int VP8RGBToV(int r, int g, int b) {
  const int v = +28800 * r - 24116 * g - 4684 * b;
  return VP8ClipUV(v);
}

#else

// This JPEG-YUV colorspace, only for comparison!
// These are also 16-bit precision coefficients from Rec.601, but with full
// [0..255] output range.
static WEBP_INLINE int VP8RGBToY(int r, int g, int b) {
  const int kRound = (1 << (YUV_FIX - 1));
  const int luma = 19595 * r + 38470 * g + 7471 * b;
  return (luma + kRound) >> YUV_FIX;  // no need to clip
}

static WEBP_INLINE int VP8RGBToU(int r, int g, int b) {
  const int u = -11058 * r - 21710 * g + 32768 * b;
  return VP8ClipUV(u);
}

static WEBP_INLINE int VP8RGBToV(int r, int g, int b) {
  const int v = 32768 * r - 27439 * g - 5329 * b;
  return VP8ClipUV(v);
}

#endif    // USE_YUVj

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif  /* WEBP_DSP_YUV_H_ */
