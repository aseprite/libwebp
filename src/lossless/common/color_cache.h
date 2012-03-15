// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)

#ifndef WEBP_COLOR_CACHE_H_
#define WEBP_COLOR_CACHE_H_

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdint.h>

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

static const uint32_t kHashMul = 0x1e35a7bd;

typedef struct {
  uint32_t* data_;
  uint32_t hash_shift_;
} VP8LColorCache;

static inline int VP8LColorCacheInit(VP8LColorCache* p, int hash_bits) {
  if (hash_bits == 0) {
    hash_bits = 1;
  }
  p->hash_shift_ = 32 - hash_bits;
  p->data_ = (uint32_t*)calloc(1 << hash_bits, sizeof(p->data_[0]));
  if (p->data_ == NULL) {
    return 0;
  }
  return 1;
}

static inline void VP8LColorCacheInsert(VP8LColorCache* p, uint32_t argb) {
  const uint32_t key = (kHashMul * argb) >> p->hash_shift_;
  p->data_[key] = argb;
}

static inline void VP8LColorCacheDelete(VP8LColorCache* p) {
  if (p->data_) {
    free(p->data_);
  }
}

static inline int VP8LColorCacheGetIndex(const VP8LColorCache* p,
                                         uint32_t argb) {
  return (kHashMul * argb) >> p->hash_shift_;
}

static inline int VP8LColorCacheContains(const VP8LColorCache* p,
                                         uint32_t argb) {
  uint32_t key = (kHashMul * argb) >> p->hash_shift_;
  return p->data_[key] == argb;
}

static inline uint32_t VP8LColorCacheLookup(VP8LColorCache* p, uint32_t key) {
  return p->data_[key];
}

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif  // WEBP_COLOR_CACHE_H_
