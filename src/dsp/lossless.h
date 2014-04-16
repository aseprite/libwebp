// Copyright 2012 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
// Image transforms and color space conversion methods for lossless decoder.
//
// Authors: Vikas Arora (vikaas.arora@gmail.com)
//          Jyrki Alakuijala (jyrki@google.com)

#ifndef WEBP_DSP_LOSSLESS_H_
#define WEBP_DSP_LOSSLESS_H_

#include "../webp/types.h"
#include "../webp/decode.h"
#include "../webp/format_constants.h"

#ifdef __cplusplus
extern "C" {
#endif

//------------------------------------------------------------------------------
// Signatures and generic function-pointers

typedef uint32_t (*VP8LPredictorFunc)(uint32_t left, const uint32_t* const top);
extern VP8LPredictorFunc VP8LPredictors[16];

typedef void (*VP8LProcessBlueAndRedFunc)(uint32_t* argb_data, int num_pixels);
extern VP8LProcessBlueAndRedFunc VP8LSubtractGreenFromBlueAndRed;
extern VP8LProcessBlueAndRedFunc VP8LAddGreenToBlueAndRed;

typedef struct {
  // Note: the members are uint8_t, so that any negative values are
  // automatically converted to "mod 256" values.
  uint8_t green_to_red_;
  uint8_t green_to_blue_;
  uint8_t red_to_blue_;
} VP8LMultipliers;
typedef void (*VP8LTransformColorFunc)(const VP8LMultipliers* const m,
                                       uint32_t* argb_data, int num_pixels);
extern VP8LTransformColorFunc VP8LTransformColor;
extern VP8LTransformColorFunc VP8LTransformColorInverse;

typedef void (*VP8LConvertFunc)(const uint32_t* src, int num_pixels,
                                uint8_t* dst);
extern VP8LConvertFunc VP8LConvertBGRAToRGB;
extern VP8LConvertFunc VP8LConvertBGRAToRGBA;
extern VP8LConvertFunc VP8LConvertBGRAToRGBA4444;
extern VP8LConvertFunc VP8LConvertBGRAToRGB565;
extern VP8LConvertFunc VP8LConvertBGRAToBGR;

// Expose some C-only fallback functions
extern void VP8LTransformColor_C(const VP8LMultipliers* const m,
                                 uint32_t* data, int num_pixels);
extern void VP8LTransformColorInverse_C(const VP8LMultipliers* const m,
                                        uint32_t* data, int num_pixels);

extern void VP8LConvertBGRAToRGB_C(const uint32_t* src,
                                   int num_pixels, uint8_t* dst);
extern void VP8LConvertBGRAToRGBA_C(const uint32_t* src,
                                    int num_pixels, uint8_t* dst);
extern void VP8LConvertBGRAToRGBA4444_C(const uint32_t* src,
                                        int num_pixels, uint8_t* dst);
extern void VP8LConvertBGRAToRGB565_C(const uint32_t* src,
                                      int num_pixels, uint8_t* dst);
extern void VP8LConvertBGRAToBGR_C(const uint32_t* src,
                                   int num_pixels, uint8_t* dst);
extern void VP8LSubtractGreenFromBlueAndRed_C(uint32_t* argb_data,
                                              int num_pixels);
extern void VP8LAddGreenToBlueAndRed_C(uint32_t* data, int num_pixels);

// Must be called before calling any of the above methods.
void VP8LDspInit(void);

//------------------------------------------------------------------------------
// Image transforms.

struct VP8LTransform;  // Defined in dec/vp8li.h.

// Performs inverse transform of data given transform information, start and end
// rows. Transform will be applied to rows [row_start, row_end[.
// The *in and *out pointers refer to source and destination data respectively
// corresponding to the intermediate row (row_start).
void VP8LInverseTransform(const struct VP8LTransform* const transform,
                          int row_start, int row_end,
                          const uint32_t* const in, uint32_t* const out);

// Similar to the static method ColorIndexInverseTransform() that is part of
// lossless.c, but used only for alpha decoding. It takes uint8_t (rather than
// uint32_t) arguments for 'src' and 'dst'.
void VP8LColorIndexInverseTransformAlpha(
    const struct VP8LTransform* const transform, int y_start, int y_end,
    const uint8_t* src, uint8_t* dst);

void VP8LResidualImage(int width, int height, int bits,
                       uint32_t* const argb, uint32_t* const argb_scratch,
                       uint32_t* const image);

void VP8LColorSpaceTransform(int width, int height, int bits, int quality,
                             uint32_t* const argb, uint32_t* image);

//------------------------------------------------------------------------------
// Color space conversion.

// Converts from BGRA to other color spaces.
void VP8LConvertFromBGRA(const uint32_t* const in_data, int num_pixels,
                         WEBP_CSP_MODE out_colorspace, uint8_t* const rgba);

//------------------------------------------------------------------------------
// Misc methods.

// Computes sampled size of 'size' when sampling using 'sampling bits'.
static WEBP_INLINE uint32_t VP8LSubSampleSize(uint32_t size,
                                              uint32_t sampling_bits) {
  return (size + (1 << sampling_bits) - 1) >> sampling_bits;
}

// -----------------------------------------------------------------------------
// Faster logarithm for integers. Small values use a look-up table.
#define LOG_LOOKUP_IDX_MAX 256
extern const float kLog2Table[LOG_LOOKUP_IDX_MAX];
extern const float kSLog2Table[LOG_LOOKUP_IDX_MAX];
typedef float (*VP8LFastLog2SlowFunc)(int v);

extern VP8LFastLog2SlowFunc VP8LFastLog2Slow;
extern VP8LFastLog2SlowFunc VP8LFastSLog2Slow;

static WEBP_INLINE float VP8LFastLog2(int v) {
  return (v < LOG_LOOKUP_IDX_MAX) ? kLog2Table[v] : VP8LFastLog2Slow(v);
}
// Fast calculation of v * log2(v) for integer input.
static WEBP_INLINE float VP8LFastSLog2(int v) {
  return (v < LOG_LOOKUP_IDX_MAX) ? kSLog2Table[v] : VP8LFastSLog2Slow(v);
}

// -----------------------------------------------------------------------------
// cost and entropy related functions.

static WEBP_INLINE double VP8LBitsEntropyRefine(int nonzeros, int sum,
                                                int max_val, double retval) {
  double mix;
  if (nonzeros < 5) {
    if (nonzeros <= 1) {
      return 0;
    }
    // Two symbols, they will be 0 and 1 in a Huffman code.
    // Let's mix in a bit of entropy to favor good clustering when
    // distributions of these are combined.
    if (nonzeros == 2) {
      return 0.99 * sum + 0.01 * retval;
    }
    // No matter what the entropy says, we cannot be better than min_limit
    // with Huffman coding. I am mixing a bit of entropy into the
    // min_limit since it produces much better (~0.5 %) compression results
    // perhaps because of better entropy clustering.
    if (nonzeros == 3) {
      mix = 0.95;
    } else {
      mix = 0.7;  // nonzeros == 4.
    }
  } else {
    mix = 0.627;
  }

  {
    double min_limit = 2 * sum - max_val;
    min_limit = mix * min_limit + (1.0 - mix) * retval;
    return (retval < min_limit) ? min_limit : retval;
  }
}

static WEBP_INLINE double VP8LBitsEntropy(const int* const array, int n) {
  double retval = 0.;
  int sum = 0;
  int nonzeros = 0;
  int max_val = 0;
  int i;
  for (i = 0; i < n; ++i) {
    if (array[i] != 0) {
      sum += array[i];
      ++nonzeros;
      retval -= VP8LFastSLog2(array[i]);
      if (max_val < array[i]) {
        max_val = array[i];
      }
    }
  }
  retval += VP8LFastSLog2(sum);
  return VP8LBitsEntropyRefine(nonzeros, sum, max_val, retval);
}

static WEBP_INLINE double VP8LBitsEntropyCombined(const int* const X,
                                                  const int* const Y,
                                                  int n) {
  double retval = 0.;
  int sum = 0;
  int nonzeros = 0;
  int max_val = 0;
  int i;
  for (i = 0; i < n; ++i) {
    const int xy = X[i] + Y[i];
    if (xy != 0) {
      sum += xy;
      ++nonzeros;
      retval -= VP8LFastSLog2(xy);
      if (max_val < xy) {
        max_val = xy;
      }
    }
  }
  retval += VP8LFastSLog2(sum);
  return VP8LBitsEntropyRefine(nonzeros, sum, max_val, retval);
}

static WEBP_INLINE double VP8LInitialHuffmanCost(void) {
  // Small bias because Huffman code length is typically not stored in
  // full length.
  static const int kHuffmanCodeOfHuffmanCodeSize = CODE_LENGTH_CODES * 3;
  static const double kSmallBias = 9.1;
  return kHuffmanCodeOfHuffmanCodeSize - kSmallBias;
}

// Used to finalized the Huffman cost:
// cnt_z / cnt_nz: counts the number of 0's and non-0's
// streak_{z,nz}_le3 / streak_{z,nz}_gt3: number of streaks larger than 3
// or less-or-equal than 3.
static WEBP_INLINE double VP8LFinalHuffmanCost(int cnt_z, int streak_z_le3,
                                               int streak_z_gt3, int cnt_nz,
                                               int streak_nz_le3,
                                               int streak_nz_gt3) {
  double retval = VP8LInitialHuffmanCost();
  retval += cnt_z * 1.5625 + 0.234375 * streak_z_gt3;
  retval += cnt_nz * 2.578125 + 0.703125 * streak_nz_gt3;
  retval += 1.796875 * streak_z_le3;
  retval += 3.28125 * streak_nz_le3;
  return retval;
}

typedef double (*VP8LCostFunc)(const int* const population, int length);
typedef double (*VP8LCostCombinedFunc)(const int* const X,
                                       const int* const Y,
                                       int length);

extern VP8LCostFunc VP8LExtraCost;
extern VP8LCostCombinedFunc VP8LExtraCostCombined;
extern VP8LCostFunc VP8LPopulationCost;
extern VP8LCostCombinedFunc VP8LGetCombinedEntropy;

// -----------------------------------------------------------------------------
// PrefixEncode()

// use GNU builtins where available.
#if defined(__GNUC__) && \
    ((__GNUC__ == 3 && __GNUC_MINOR__ >= 4) || __GNUC__ >= 4)
static WEBP_INLINE int BitsLog2Floor(uint32_t n) {
  return 31 ^ __builtin_clz(n);
}
#elif defined(_MSC_VER) && _MSC_VER > 1310 && \
      (defined(_M_X64) || defined(_M_IX86))
#include <intrin.h>
#pragma intrinsic(_BitScanReverse)

static WEBP_INLINE int BitsLog2Floor(uint32_t n) {
  unsigned long first_set_bit;
  _BitScanReverse(&first_set_bit, n);
  return first_set_bit;
}
#else
// Returns (int)floor(log2(n)). n must be > 0.
static WEBP_INLINE int BitsLog2Floor(uint32_t n) {
  int log = 0;
  uint32_t value = n;
  int i;

  for (i = 4; i >= 0; --i) {
    const int shift = (1 << i);
    const uint32_t x = value >> shift;
    if (x != 0) {
      value = x;
      log += shift;
    }
  }
  return log;
}
#endif

static WEBP_INLINE int VP8LBitsLog2Ceiling(uint32_t n) {
  const int log_floor = BitsLog2Floor(n);
  if (n == (n & ~(n - 1)))  // zero or a power of two.
    return log_floor;
  else
    return log_floor + 1;
}

// Splitting of distance and length codes into prefixes and
// extra bits. The prefixes are encoded with an entropy code
// while the extra bits are stored just as normal bits.
static WEBP_INLINE void VP8LPrefixEncodeBitsNoLUT(int distance, int* const code,
                                                  int* const extra_bits) {
  const int highest_bit = BitsLog2Floor(--distance);
  const int second_highest_bit = (distance >> (highest_bit - 1)) & 1;
  *extra_bits = highest_bit - 1;
  *code = 2 * highest_bit + second_highest_bit;
}

static WEBP_INLINE void VP8LPrefixEncodeNoLUT(int distance, int* const code,
                                              int* const extra_bits,
                                              int* const extra_bits_value) {
  const int highest_bit = BitsLog2Floor(--distance);
  const int second_highest_bit = (distance >> (highest_bit - 1)) & 1;
  *extra_bits = highest_bit - 1;
  *extra_bits_value = distance & ((1 << *extra_bits) - 1);
  *code = 2 * highest_bit + second_highest_bit;
}

#define PREFIX_LOOKUP_IDX_MAX   512
typedef struct {
  int8_t code_;
  int8_t extra_bits_;
} VP8LPrefixCode;

// These tables are derived using VP8LPrefixEncodeNoLUT.
extern const VP8LPrefixCode kPrefixEncodeCode[PREFIX_LOOKUP_IDX_MAX];
extern const uint8_t kPrefixEncodeExtraBitsValue[PREFIX_LOOKUP_IDX_MAX];
static WEBP_INLINE void VP8LPrefixEncodeBits(int distance, int* const code,
                                             int* const extra_bits) {
  if (distance < PREFIX_LOOKUP_IDX_MAX) {
    const VP8LPrefixCode prefix_code = kPrefixEncodeCode[distance];
    *code = prefix_code.code_;
    *extra_bits = prefix_code.extra_bits_;
  } else {
    VP8LPrefixEncodeBitsNoLUT(distance, code, extra_bits);
  }
}

static WEBP_INLINE void VP8LPrefixEncode(int distance, int* const code,
                                         int* const extra_bits,
                                         int* const extra_bits_value) {
  if (distance < PREFIX_LOOKUP_IDX_MAX) {
    const VP8LPrefixCode prefix_code = kPrefixEncodeCode[distance];
    *code = prefix_code.code_;
    *extra_bits = prefix_code.extra_bits_;
    *extra_bits_value = kPrefixEncodeExtraBitsValue[distance];
  } else {
    VP8LPrefixEncodeNoLUT(distance, code, extra_bits, extra_bits_value);
  }
}

// In-place difference of each component with mod 256.
static WEBP_INLINE uint32_t VP8LSubPixels(uint32_t a, uint32_t b) {
  const uint32_t alpha_and_green =
      0x00ff00ffu + (a & 0xff00ff00u) - (b & 0xff00ff00u);
  const uint32_t red_and_blue =
      0xff00ff00u + (a & 0x00ff00ffu) - (b & 0x00ff00ffu);
  return (alpha_and_green & 0xff00ff00u) | (red_and_blue & 0x00ff00ffu);
}

void VP8LBundleColorMap(const uint8_t* const row, int width,
                        int xbits, uint32_t* const dst);

//------------------------------------------------------------------------------

#ifdef __cplusplus
}    // extern "C"
#endif

#endif  // WEBP_DSP_LOSSLESS_H_
