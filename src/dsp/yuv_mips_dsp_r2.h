// Copyright 2014 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
// MIPS dspr2 version of some functions
//
// Author(s): Branimir Vasic (branimir.vasic@imgtec.com)
//            Djordje Pesut  (djordje.pesut@imgtec.com)

#ifndef WEBP_DSP_YUV_MIPS_DSP_R2_H_
#define WEBP_DSP_YUV_MIPS_DSP_R2_H_

#if !defined(WEBP_YUV_USE_TABLE)

#if defined(WEBP_USE_MIPS32)
static WEBP_INLINE int VP8Clip8(int v) {
  return ((v & ~YUV_MASK2) == 0) ? (v >> YUV_FIX2) : (v < 0) ? 0 : 255;
}
#endif  // WEBP_USE_MIPS32

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

static WEBP_INLINE void VP8YuvToRgb(int y, int u, int v,
                                    uint8_t* const rgb) {
  int t1, t2, t3, t4, t5;
  VP8YUV_FUNCTIONS(y,u,v);
  rgb[0] = t2;
  rgb[1] = t3;
  rgb[2] = t4;
}
static WEBP_INLINE void VP8YuvToBgr(int y, int u, int v,
                                    uint8_t* const bgr) {
  int t1, t2, t3, t4, t5;
  VP8YUV_FUNCTIONS(y,u,v);
  bgr[0] = t4;
  bgr[1] = t3;
  bgr[2] = t2;
}
static WEBP_INLINE void VP8YuvToRgb565(int y, int u, int v,
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
static WEBP_INLINE void VP8YuvToRgba4444(int y, int u, int v,
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

static WEBP_INLINE void VP8YuvToArgb(uint8_t y, uint8_t u, uint8_t v,
                                     uint8_t* const argb) {
  int t1, t2, t3, t4, t5;
  VP8YUV_FUNCTIONS(y,u,v);
  argb[0] = 0xff;
  argb[1] =(t2);
  argb[2] =(t3);
  argb[3] =(t4);
}
static WEBP_INLINE void VP8YuvToBgra(uint8_t y, uint8_t u, uint8_t v,
                                     uint8_t* const bgra) {
  int t1, t2, t3, t4, t5;
  VP8YUV_FUNCTIONS(y,u,v);
  bgra[0] =(t4);
  bgra[1] =(t3);
  bgra[2] =(t2);
  bgra[3] = 0xff;
}
static WEBP_INLINE void VP8YuvToRgba(uint8_t y, uint8_t u, uint8_t v,
                                     uint8_t* const rgba) {
  int t1, t2, t3, t4, t5;
  VP8YUV_FUNCTIONS(y,u,v);
  rgba[0] =(t2);
  rgba[1] =(t3);
  rgba[2] =(t4);
  rgba[3] = 0xff;
}

#endif  // WEBP_DSP_YUV_MIPS_DSP_R2_H_
