// Copyright 2012 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Color Cache for WebP Lossless
//
// Author: jyrki@google.com (Jyrki Alakuijala)

#include <assert.h>
#include <stdlib.h>
#include "./color_cache.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

//------------------------------------------------------------------------------
// Helper functions.

// TODO: This function is static in vp8l.c as well. Move to a common header.
static WEBP_INLINE uint32_t SubSampleSize(uint32_t size,
                                          uint32_t sampling_bits) {
  return (size + (1 << sampling_bits) - 1) >> sampling_bits;
}

//------------------------------------------------------------------------------
// ColorCacheColumn.

static WEBP_INLINE void ColumnInit(ColorCacheColumn* const cc, int hash_bits) {
  uint32_t i;
  assert(cc != NULL);
  cc->hash_shift_ = 32 - hash_bits;
  cc->hash_size_ = 1 << hash_bits;
  cc->data_ = (uint32_t*)malloc(cc->hash_size_ * sizeof(cc->data_[0]));
  for (i = 0; i < cc->hash_size_; ++i) {
    cc->data_[i] = kNotInitialized;
  }
}

static WEBP_INLINE void ColumnRelease(ColorCacheColumn* const cc) {
  free(cc->data_);
  cc->data_ = NULL;
}

// Returns 1 if the key location in the column contains this ARGB value,
//         0 if the key contains a different ARGB value and
//        -1 if the key is uninitialized.
static WEBP_INLINE int ColumnContains(const ColorCacheColumn* const cc,
                                      uint32_t argb) {
  assert(cc != NULL);
  if (argb == kNotInitialized) {
    return 0;
  } else {
    const uint32_t key = ColorCacheColumnGetKey(argb, cc->hash_shift_);
    if (cc->data_[key] == kNotInitialized) {
      return -1;
    } else {
      return (cc->data_[key] == argb);
    }
  }
}

//------------------------------------------------------------------------------
// Helper functions for VP8LColorCache.

static WEBP_INLINE uint32_t ColorCacheGetIndex(
    const VP8LColorCache* const color_cache, uint32_t argb) {
  assert(color_cache != NULL);
  return ColorCacheColumnGetKey(argb, 32 - color_cache->hash_bits_);
}

static WEBP_INLINE int ColorCacheContains(
    const VP8LColorCache* const color_cache, uint32_t x_pos, uint32_t argb) {
  int i;
  int ret;
  const int col = ColorCacheGetColumn(color_cache, x_pos);
  const ColorCacheColumn* const columns = &color_cache->cache_columns_[col];

  // Search in current column.
  ret = ColumnContains(&columns[0], argb);
  if (ret != -1) return ret;

  // Search in nearby columns one-by-one.
  for (i = 1; i < color_cache->num_cache_columns_; ++i) {
    if (col - i >= 0) {
      ret = ColumnContains(&columns[-i], argb);
      if (ret != -1) return ret;
    }
    if (col + i < color_cache->num_cache_columns_) {
      ret = ColumnContains(&columns[i], argb);
      if (ret != -1) return ret;
    }
  }

  // No match.
  return 0;
}

//------------------------------------------------------------------------------
// VP8LColorCache.

void VP8LColorCacheInit(VP8LColorCache* const color_cache, uint32_t xsize,
                        int x_downsample_bits, int hash_bits) {
  int i;
  assert(color_cache != NULL);
  if (hash_bits == 0) {
    hash_bits = 1;
  }
  color_cache->x_downsample_bits_ = x_downsample_bits;
  color_cache->hash_bits_ = hash_bits;
  color_cache->num_cache_columns_ = SubSampleSize(xsize, x_downsample_bits);
  color_cache->cache_columns_ =
      (ColorCacheColumn*)malloc(color_cache->num_cache_columns_ *
                                sizeof(color_cache->cache_columns_[0]));
  for (i = 0; i < color_cache->num_cache_columns_; ++i) {
    ColumnInit(&color_cache->cache_columns_[i], hash_bits);
  }
}

void VP8LColorCacheRelease(VP8LColorCache* const color_cache) {
  int i;
  assert(color_cache != NULL);
  for (i = 0; i < color_cache->num_cache_columns_; ++i) {
    ColumnRelease(&color_cache->cache_columns_[i]);
  }
  free(color_cache->cache_columns_);
}

int VL8LColorCacheInverseLookup(const VP8LColorCache* const color_cache,
                                uint32_t x_pos, uint32_t argb,
                                uint32_t* const key) {
  assert(color_cache != NULL);
  assert(key != NULL);
  if(ColorCacheContains(color_cache, x_pos, argb)) {
    *key = ColorCacheGetIndex(color_cache, argb);
    return 1;
  }
  return 0;
}

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

