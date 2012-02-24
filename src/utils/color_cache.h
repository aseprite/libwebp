// Copyright 2012 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Color Cache for WebP Lossless
//
// Authors: jyrki@google.com (Jyrki Alakuijala)
//          urvang@google.com (Urvang Joshi)

#ifndef WEBP_UTILS_COLOR_CACHE_H_
#define WEBP_UTILS_COLOR_CACHE_H_

#include "../webp/types.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

// Represents cache for colors in a vertical stripe of an image.
typedef struct {
  uint32_t *data_;
} VP8LColorCacheColumn;

// Main color cache struct.
typedef struct {
  VP8LColorCacheColumn* cache_columns_;  // A collection of cache columns.
  int num_cache_columns_;      // Number of cache columns
  int hash_size_;              // Size of a cache column: 1 << hash_bits.
                               // (where hash_bits = number of bits used for
                               // hash addressing in a cache column).
  int hash_shift_;             // Hash shift: 32 - hash_bits.
  int x_downsample_bits_;      // Number of columns:
                               // ceil(image_width/x_downsample_bits_).
} VP8LColorCache;

//------------------------------------------------------------------------------
// Helper methods.

static const uint32_t kHashMul = 0x1e35a7bd;
// kNotInitialized is a special value which can be inserted into the
// ColorCacheColumn, but is never recalled as a value.
static const uint32_t kNotInitialized = 0x1e35a7bd;

static WEBP_INLINE uint32_t ColorCacheColumnGetKey(uint32_t argb,
                                                   int hash_shift) {
  return (kHashMul * argb) >> hash_shift;
}

static WEBP_INLINE int ColorCacheColumnLookup(
    const VP8LColorCacheColumn* const cc, uint32_t key, uint32_t* const argb) {
  if (cc->data_[key] != kNotInitialized) {
    *argb = cc->data_[key];
    return 1;
  }
  return 0;
}

static WEBP_INLINE void ColorCacheColumnInsert(VP8LColorCacheColumn* const cc,
                                               uint32_t argb,
                                               int hash_shift) {
  const uint32_t key = ColorCacheColumnGetKey(argb, hash_shift);
  cc->data_[key] = argb;
}

//------------------------------------------------------------------------------
// Main APIs.

static WEBP_INLINE int ColorCacheGetColumn(
    const VP8LColorCache* const color_cache, uint32_t x_pos) {
  assert(color_cache != NULL);
  return x_pos >> color_cache->x_downsample_bits_;
}

// Initializes the color cache based on the given parameters.
// Note: 'color_cache' must already be allocated.
void VP8LColorCacheInit(VP8LColorCache* const color_cache, uint32_t xsize,
                        int x_downsample_bits, int hash_bits);

// Releases the data owned by color cache.
// Note: It does NOT free 'color_cache' itself.
void VP8LColorCacheRelease(VP8LColorCache* const color_cache);

// TODO(urvang): Use 'argb_t' type for 'argb' once 'argb_t' is part of
// webp/types.h
// Given the ARGB value and x position of a pixel, insert it into the cache.
static WEBP_INLINE void VP8LColorCacheInsert(VP8LColorCache* const color_cache,
                                             uint32_t x_pos, uint32_t argb) {
  const int col = ColorCacheGetColumn(color_cache, x_pos);
  ColorCacheColumnInsert(&color_cache->cache_columns_[col], argb,
                         color_cache->hash_shift_);
}

// Given the key and x position of a pixel, find out its ARGB value.
// Returns 1 if the ARGB value was found.
static WEBP_INLINE int VP8LColorCacheLookup(
    const VP8LColorCache* const color_cache, uint32_t x_pos, uint32_t key,
    uint32_t* const argb) {
  int i, col;
  assert(color_cache != NULL);
  assert(argb != NULL);
  col = ColorCacheGetColumn(color_cache, x_pos);
  if (ColorCacheColumnLookup(&color_cache->cache_columns_[col], key, argb)) {
    return 1;
  }
  for (i = 1; i < color_cache->num_cache_columns_; ++i) {
    if (col - i >= 0 &&
        ColorCacheColumnLookup(&color_cache->cache_columns_[col - i], key,
                               argb)) {
      return 1;
    }
    if (col + i < color_cache->num_cache_columns_ &&
        ColorCacheColumnLookup(&color_cache->cache_columns_[col + i], key,
                               argb)) {
      return 1;
    }
  }
  return 0;
}

// Given the ARGB value and x position of a pixel, find out its key.
// Returns 1 if the key was found.
int VL8LColorCacheInverseLookup(const VP8LColorCache* const color_cache,
                                uint32_t x_pos, uint32_t argb,
                                uint32_t* const key);

//------------------------------------------------------------------------------

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif  // WEBP_UTILS_COLOR_CACHE_H_
