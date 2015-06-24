// Copyright 2015 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
// SSE2 variant of methods for lossless encoder
//
// Author: Skal (pascal.massimino@gmail.com)

#include "./dsp.h"

#if defined(WEBP_USE_SSE2)
#include <assert.h>
#include <emmintrin.h>
#include "./lossless.h"

// For sign-extended multiplying constants, pre-shifted by 5:
#define CST_5b(X)  (((int16_t)(X << 8)) >> 5)

//------------------------------------------------------------------------------
// Subtract-Green Transform

static void SubtractGreenFromBlueAndRed(uint32_t* argb_data, int num_pixels) {
  int i;
  for (i = 0; i + 4 <= num_pixels; i += 4) {
    const __m128i in = _mm_loadu_si128((__m128i*)&argb_data[i]); // argb
    const __m128i A = _mm_srli_epi16(in, 8);     // 0 a 0 g
    const __m128i B = _mm_shufflelo_epi16(A, _MM_SHUFFLE(2, 2, 0, 0));
    const __m128i C = _mm_shufflehi_epi16(B, _MM_SHUFFLE(2, 2, 0, 0));  // 0g0g
    const __m128i out = _mm_sub_epi8(in, C);
    _mm_storeu_si128((__m128i*)&argb_data[i], out);
  }
  // fallthrough and finish off with plain-C
  VP8LSubtractGreenFromBlueAndRed_C(argb_data + i, num_pixels - i);
}

//------------------------------------------------------------------------------
// Color Transform

static WEBP_INLINE void TransformColor(const VP8LMultipliers* const m,
                                       uint32_t* argb_data, int num_pixels) {
  const __m128i mults_rb = _mm_set_epi16(
      CST_5b(m->green_to_red_), CST_5b(m->green_to_blue_),
      CST_5b(m->green_to_red_), CST_5b(m->green_to_blue_),
      CST_5b(m->green_to_red_), CST_5b(m->green_to_blue_),
      CST_5b(m->green_to_red_), CST_5b(m->green_to_blue_));
  const __m128i mults_b2 = _mm_set_epi16(
      CST_5b(m->red_to_blue_), 0, CST_5b(m->red_to_blue_), 0,
      CST_5b(m->red_to_blue_), 0, CST_5b(m->red_to_blue_), 0);
  const __m128i mask_ag = _mm_set1_epi32(0xff00ff00);  // alpha-green masks
  const __m128i mask_rb = _mm_set1_epi32(0x00ff00ff);  // red-blue masks
  int i;
  for (i = 0; i + 4 <= num_pixels; i += 4) {
    const __m128i in = _mm_loadu_si128((__m128i*)&argb_data[i]); // argb
    const __m128i A = _mm_and_si128(in, mask_ag);     // a   0   g   0
    const __m128i B = _mm_shufflelo_epi16(A, _MM_SHUFFLE(2, 2, 0, 0));
    const __m128i C = _mm_shufflehi_epi16(B, _MM_SHUFFLE(2, 2, 0, 0));  // g0g0
    const __m128i D = _mm_mulhi_epi16(C, mults_rb);    // x dr  x db1
    const __m128i E = _mm_slli_epi16(in, 8);           // r 0   b   0
    const __m128i F = _mm_mulhi_epi16(E, mults_b2);    // x db2 0   0
    const __m128i G = _mm_srli_epi32(F, 16);           // 0 0   x db2
    const __m128i H = _mm_add_epi8(G, D);              // x dr  x  db
    const __m128i I = _mm_and_si128(H, mask_rb);       // 0 dr  0  db
    const __m128i out = _mm_sub_epi8(in, I);
    _mm_storeu_si128((__m128i*)&argb_data[i], out);
  }
  // fallthrough and finish off with plain-C
  VP8LTransformColor_C(m, argb_data + i, num_pixels - i);
}

static void CollectColorBlueTransforms(const uint32_t* argb, int stride,
                                       int tile_width, int tile_height,
                                       int green_to_blue, int red_to_blue,
                                       int histo[]) {
  const __m128i mults_r = _mm_set_epi16(
      CST_5b(red_to_blue), 0, CST_5b(red_to_blue), 0,
      CST_5b(red_to_blue), 0, CST_5b(red_to_blue), 0);
  const __m128i mults_g = _mm_set_epi16(
      0, CST_5b(green_to_blue), 0, CST_5b(green_to_blue),
      0, CST_5b(green_to_blue), 0, CST_5b(green_to_blue));
  const __m128i mask_g = _mm_set1_epi32(0x00ff00);  // green mask
  const __m128i mask_b = _mm_set1_epi32(0x0000ff);  // blue mask
  const int span = 4;
  int y;
  for (y = 0; y < tile_height; ++y) {
    const uint32_t* const src = argb + y * stride;
    int i, x;
    for (x = 0; x + span <= tile_width; x += span) {
      uint32_t values[span];
      const __m128i in = _mm_loadu_si128((__m128i*)&src[x]);  // argb
      const __m128i A = _mm_slli_epi16(in, 8);        // r 0  | b 0
      const __m128i B = _mm_and_si128(in, mask_g);    // 0 0  | g 0
      const __m128i C = _mm_mulhi_epi16(A, mults_r);  // x db | 0 0
      const __m128i D = _mm_mulhi_epi16(B, mults_g);  // 0 0  | x db
      const __m128i E = _mm_sub_epi8(in, D);          // x x  | x b'
      const __m128i F = _mm_srli_epi32(C, 16);        // 0 0  | x db
      const __m128i G = _mm_sub_epi8(E, F);           // 0 0  | x b'
      const __m128i H = _mm_and_si128(G, mask_b);     // 0 0  | 0 b
      _mm_storeu_si128((__m128i*)values, H);
      for (i = 0; i < span; ++i) ++histo[values[i]];
    }
  }
  {
    const int left_over = tile_width & (span - 1);
    if (left_over > 0) {
      VP8LCollectColorBlueTransforms_C(argb + tile_width - left_over, stride,
                                       left_over, tile_height,
                                       green_to_blue, red_to_blue, histo);
    }
  }
}

static void CollectColorRedTransforms(const uint32_t* argb, int stride,
                                      int tile_width, int tile_height,
                                      int green_to_red, int histo[]) {
  const __m128i mults_g = _mm_set_epi16(
      0, CST_5b(green_to_red), 0, CST_5b(green_to_red),
      0, CST_5b(green_to_red), 0, CST_5b(green_to_red));
  const __m128i mask_g = _mm_set1_epi32(0x00ff00);  // green mask
  const __m128i mask = _mm_set1_epi32(0xff);
  const int span = 4;
  int y;
  for (y = 0; y < tile_height; ++y) {
    const uint32_t* const src = argb + y * stride;
    int i, x;
    for (x = 0; x + span <= tile_width; x += span) {
      uint32_t values[span];
      const __m128i in = _mm_loadu_si128((__m128i*)&src[x]);  // argb
      const __m128i A = _mm_and_si128(in, mask_g);    // 0 0  | g 0
      const __m128i B = _mm_srli_epi32(in, 16);       // 0 0  | x r
      const __m128i C = _mm_mulhi_epi16(A, mults_g);  // 0 0  | x dr
      const __m128i E = _mm_sub_epi8(B, C);           // x x  | x r'
      const __m128i F = _mm_and_si128(E, mask);       // 0 0  | 0 r'
      _mm_storeu_si128((__m128i*)values, F);
      for (i = 0; i < span; ++i) ++histo[values[i]];
    }
  }
  {
    const int left_over = tile_width & (span - 1);
    if (left_over > 0) {
      VP8LCollectColorRedTransforms_C(argb + tile_width - left_over, stride,
                                       left_over, tile_height,
                                       green_to_red, histo);
    }
  }
}

//------------------------------------------------------------------------------

#define LINE_SIZE 16    // 8 or 16
static void AddVector(const uint32_t* a, const uint32_t* b, uint32_t* out,
                      int size) {
  int i;
  assert(size % LINE_SIZE == 0);
  for (i = 0; i < size; i += LINE_SIZE) {
    const __m128i a0 = _mm_loadu_si128((const __m128i*)&a[i +  0]);
    const __m128i a1 = _mm_loadu_si128((const __m128i*)&a[i +  4]);
#if (LINE_SIZE == 16)
    const __m128i a2 = _mm_loadu_si128((const __m128i*)&a[i +  8]);
    const __m128i a3 = _mm_loadu_si128((const __m128i*)&a[i + 12]);
#endif
    const __m128i b0 = _mm_loadu_si128((const __m128i*)&b[i +  0]);
    const __m128i b1 = _mm_loadu_si128((const __m128i*)&b[i +  4]);
#if (LINE_SIZE == 16)
    const __m128i b2 = _mm_loadu_si128((const __m128i*)&b[i +  8]);
    const __m128i b3 = _mm_loadu_si128((const __m128i*)&b[i + 12]);
#endif
    _mm_storeu_si128((__m128i*)&out[i +  0], _mm_add_epi32(a0, b0));
    _mm_storeu_si128((__m128i*)&out[i +  4], _mm_add_epi32(a1, b1));
#if (LINE_SIZE == 16)
    _mm_storeu_si128((__m128i*)&out[i +  8], _mm_add_epi32(a2, b2));
    _mm_storeu_si128((__m128i*)&out[i + 12], _mm_add_epi32(a3, b3));
#endif
  }
}

static void AddVectorEq(const uint32_t* a, uint32_t* out, int size) {
  int i;
  assert(size % LINE_SIZE == 0);
  for (i = 0; i < size; i += LINE_SIZE) {
    const __m128i a0 = _mm_loadu_si128((const __m128i*)&a[i +  0]);
    const __m128i a1 = _mm_loadu_si128((const __m128i*)&a[i +  4]);
#if (LINE_SIZE == 16)
    const __m128i a2 = _mm_loadu_si128((const __m128i*)&a[i +  8]);
    const __m128i a3 = _mm_loadu_si128((const __m128i*)&a[i + 12]);
#endif
    const __m128i b0 = _mm_loadu_si128((const __m128i*)&out[i +  0]);
    const __m128i b1 = _mm_loadu_si128((const __m128i*)&out[i +  4]);
#if (LINE_SIZE == 16)
    const __m128i b2 = _mm_loadu_si128((const __m128i*)&out[i +  8]);
    const __m128i b3 = _mm_loadu_si128((const __m128i*)&out[i + 12]);
#endif
    _mm_storeu_si128((__m128i*)&out[i +  0], _mm_add_epi32(a0, b0));
    _mm_storeu_si128((__m128i*)&out[i +  4], _mm_add_epi32(a1, b1));
#if (LINE_SIZE == 16)
    _mm_storeu_si128((__m128i*)&out[i +  8], _mm_add_epi32(a2, b2));
    _mm_storeu_si128((__m128i*)&out[i + 12], _mm_add_epi32(a3, b3));
#endif
  }
}
#undef LINE_SIZE

// Note we are adding uint32_t's as *signed* int32's (using _mm_add_epi32). But
// that's ok since the histogram values are less than 1<<28 (max picture size).
static void HistogramAdd(const VP8LHistogram* const a,
                         const VP8LHistogram* const b,
                         VP8LHistogram* const out) {
  int i;
  const int literal_size = VP8LHistogramNumCodes(a->palette_code_bits_);
  assert(a->palette_code_bits_ == b->palette_code_bits_);
  if (b != out) {
    AddVector(a->literal_, b->literal_, out->literal_, NUM_LITERAL_CODES);
    AddVector(a->red_, b->red_, out->red_, NUM_LITERAL_CODES);
    AddVector(a->blue_, b->blue_, out->blue_, NUM_LITERAL_CODES);
    AddVector(a->alpha_, b->alpha_, out->alpha_, NUM_LITERAL_CODES);
  } else {
    AddVectorEq(a->literal_, out->literal_, NUM_LITERAL_CODES);
    AddVectorEq(a->red_, out->red_, NUM_LITERAL_CODES);
    AddVectorEq(a->blue_, out->blue_, NUM_LITERAL_CODES);
    AddVectorEq(a->alpha_, out->alpha_, NUM_LITERAL_CODES);
  }
  for (i = NUM_LITERAL_CODES; i < literal_size; ++i) {
    out->literal_[i] = a->literal_[i] + b->literal_[i];
  }
  for (i = 0; i < NUM_DISTANCE_CODES; ++i) {
    out->distance_[i] = a->distance_[i] + b->distance_[i];
  }
}

//------------------------------------------------------------------------------
// Entry point

extern void VP8LEncDspInitSSE2(void);

WEBP_TSAN_IGNORE_FUNCTION void VP8LEncDspInitSSE2(void) {
  VP8LSubtractGreenFromBlueAndRed = SubtractGreenFromBlueAndRed;
  VP8LTransformColor = TransformColor;
  VP8LCollectColorBlueTransforms = CollectColorBlueTransforms;
  VP8LCollectColorRedTransforms = CollectColorRedTransforms;
  VP8LHistogramAdd = HistogramAdd;
}

#else  // !WEBP_USE_SSE2

WEBP_DSP_INIT_STUB(VP8LEncDspInitSSE2)

#endif  // WEBP_USE_SSE2
