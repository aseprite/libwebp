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

// TODO(urvang): Evaluate if inlining VP8LColorCacheInsert() and/or
// VP8LColorCacheLookup() gives a better performance.

// Represents cache for colors in a vertical stripe of an image.
typedef struct ColorCacheColumn ColorCacheColumn;

// Main color cache struct. It contains multiple cache columns.
typedef struct {
  ColorCacheColumn* cache_columns_;
  int num_cache_columns_;
  int hash_bits_;
  int x_downsample_bits_;
} VP8LColorCache;

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
void VP8LColorCacheInsert(VP8LColorCache* const color_cache, uint32_t x_pos,
                          uint32_t argb);

// Given the key and x position of a pixel, find out its ARGB value.
// Returns 1 if the ARGB value was found.
int VP8LColorCacheLookup(const VP8LColorCache* const color_cache,
                         uint32_t x_pos, uint32_t key, uint32_t* const argb);

// Given the ARGB value and x position of a pixel, find out its key.
// Returns 1 if the key was found.
int VL8LColorCacheInverseLookup(const VP8LColorCache* const color_cache,
                                uint32_t x_pos, uint32_t argb,
                                uint32_t* const key);

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif  // WEBP_UTILS_COLOR_CACHE_H_
