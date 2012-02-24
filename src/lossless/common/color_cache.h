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

// kNotInitialized is a special value which can be inserted into the
// VP8LColorCacheColumn, but is never recalled as a value.
static const uint32_t kNotInitialized = 0x1e35a7bd;

typedef struct {
  uint32_t *data_;
  uint32_t hash_shift_;
  uint32_t hash_size_;
} VP8LColorCacheColumn;

static int VP8LColorCacheColumnInit(VP8LColorCacheColumn *p, int hash_bits) {
  uint32_t i;
  p->hash_shift_ = 32 - hash_bits;
  p->hash_size_ = 1 << hash_bits;
  p->data_ = (uint32_t *)malloc(p->hash_size_ * sizeof(p->data_[0]));
  if (p->data_ == NULL) {
    return 0;
  }
  for (i = 0; i < p->hash_size_; ++i) {
    p->data_[i] = kNotInitialized;
  }
  return 1;
}

static void VP8LColorCacheColumnDelete(VP8LColorCacheColumn *p) {
  if (p->data_) {
    free(p->data_);
  }
}

static void VP8LColorCacheColumnInsert(VP8LColorCacheColumn *p, uint32_t argb) {
  const uint32_t key = (kHashMul * argb) >> p->hash_shift_;
  p->data_[key] = argb;
}

static int VP8LColorCacheColumnIsInitialized(const VP8LColorCacheColumn *p,
                                             uint32_t argb) {
  const uint32_t key = (kHashMul * argb) >> p->hash_shift_;
  return p->data_[key] != kNotInitialized;
}

static int VP8LColorCacheColumnContains(const VP8LColorCacheColumn *p,
                                        uint32_t argb) {
  uint32_t key;
  if (argb == kNotInitialized) {
    return 0;
  }
  key = (kHashMul * argb) >> p->hash_shift_;
  return p->data_[key] == argb;
}

static int VP8LColorCacheColumnLookup(VP8LColorCacheColumn *p,
                                      uint32_t key, uint32_t* argb) {
  assert(key < p->hash_size_);
  if (p->data_[key] != kNotInitialized) {
    *argb = p->data_[key];
    return 1;
  }
  return 0;
}

typedef struct {
  VP8LColorCacheColumn* hashers_;
  int hash_bits_;
  int x_downsample_bits_;
  int hashers_size_;
} VP8LColorCache;

static inline int VP8LColorCacheInit(VP8LColorCache *p,
                                      int xsize,
                                      int x_downsample_bits,
                                      int hash_bits) {
  int i;
  if (hash_bits == 0) {
    hash_bits = 1;
  }
  p->x_downsample_bits_ = x_downsample_bits;
  p->hash_bits_ = hash_bits;
  p->hashers_size_ =
      (xsize + (1 << x_downsample_bits) - 1) >> x_downsample_bits;

  p->hashers_ = (VP8LColorCacheColumn *)
      malloc(p->hashers_size_ * sizeof(p->hashers_[0]));
  if (p->hashers_ == NULL) {
    return 0;
  }
  for (i = 0; i < p->hashers_size_; ++i) {
    VP8LColorCacheColumnInit(&p->hashers_[i], hash_bits);
  }
  return 1;
}

static inline void VP8LColorCacheDelete(VP8LColorCache *p) {
  int i;
  if (p->hashers_) {
    for (i = 0; i < p->hashers_size_; ++i) {
      VP8LColorCacheColumnDelete(&p->hashers_[i]);
    }
    free(p->hashers_);
  }
}

static inline void VP8LColorCacheInsert(VP8LColorCache *p,
                                        int x, uint32_t argb) {
  VP8LColorCacheColumnInsert(&p->hashers_[x >> p->x_downsample_bits_], argb);
}

static inline uint32_t VP8LColorCacheGetIndex(const VP8LColorCache *p,
                                              uint32_t argb) {
  uint32_t val = kHashMul * argb;
  val >>= 32 - p->hash_bits_;
  return val;
}

static inline int VP8LColorCacheContains(VP8LColorCache *p,
                                         int x, uint32_t argb) {
  int i;
  int ix = x >> p->x_downsample_bits_;
  if (VP8LColorCacheColumnContains(&p->hashers_[ix], argb)) {
    return 1;
  }
  if (VP8LColorCacheColumnIsInitialized(&p->hashers_[ix], argb)) {
    return 0;
  }
  for (i = 1; i < p->hashers_size_; ++i) {
    if (ix - i >= 0) {
      if (VP8LColorCacheColumnContains(&p->hashers_[ix - i], argb)) {
        return 1;
      }
      if (VP8LColorCacheColumnIsInitialized(&p->hashers_[ix - i], argb)) {
        return 0;
      }
    }
    if (ix + i < p->hashers_size_) {
      if (VP8LColorCacheColumnContains(&p->hashers_[ix + i], argb)) {
        return 1;
      }
      if (VP8LColorCacheColumnIsInitialized(&p->hashers_[ix + i], argb)) {
        return 0;
      }
    }
  }
  return 0;
}

static inline uint32_t VP8LColorCacheLookup(VP8LColorCache *p,
                                            int x, uint32_t hash) {
  int i;
  int ix = x >> p->x_downsample_bits_;
  uint32_t argb;
  if (VP8LColorCacheColumnLookup(&p->hashers_[ix], hash, &argb)) {
    return argb;
  }
  for (i = 1; i < p->hashers_size_; ++i) {
    if (ix - i >= 0 &&
        VP8LColorCacheColumnLookup(&p->hashers_[ix - i], hash, &argb)) {
      return argb;
    }
    if (ix + i < p->hashers_size_ &&
        VP8LColorCacheColumnLookup(&p->hashers_[ix + i], hash, &argb)) {
      return argb;
    }
  }
  return 0;
}

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif  // WEBP_COLOR_CACHE_H_
