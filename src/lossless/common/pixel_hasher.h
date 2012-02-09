// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)

#ifndef WEBP_PIXEL_HASHER_H_
#define WEBP_PIXEL_HASHER_H_

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdint.h>

static const uint32_t kHashMul = 0x1e35a7bd;

// kNotInitialized is a special value which can be inserted into the
// PixelHasher, but is never recalled as a value.
static const uint32_t kNotInitialized = 0x1e35a7bd;

struct PixelHasher {
  uint32_t *data_;
  uint32_t hash_shift_;
  uint32_t hash_size_;
};

static void VP8LPixelHasherInit(PixelHasher *p, int hash_bits) {
  int i;
  p->hash_shift_ = 32 - hash_bits;
  p->hash_size_ = 1 << hash_bits;
  p->data_ = (uint32_t *)malloc(p->hash_size_ * sizeof(p->data_[0]));
  for (i = 0; i < p->hash_size_; ++i) {
    p->data_[i] = kNotInitialized;
  }
}

static void VP8LPixelHasherDelete(PixelHasher *p) {
  free(p->data_);
}

static void VP8LPixelHasherInsert(PixelHasher *p, uint32_t argb) {
  const uint32_t key = (kHashMul * argb) >> p->hash_shift_;
  p->data_[key] = argb;
}

static bool VP8LPixelHasherIsInitialized(const PixelHasher *p, uint32_t argb) {
  const uint32_t key = (kHashMul * argb) >> p->hash_shift_;
  return p->data_[key] != kNotInitialized;
}

static bool VP8LPixelHasherContains(const PixelHasher *p, uint32_t argb) {
  uint32_t key;
  if (argb == kNotInitialized) {
    return false;
  }
  key = (kHashMul * argb) >> p->hash_shift_;
  return p->data_[key] == argb;
}

static bool VP8LPixelHasherLookup(PixelHasher *p,
                                  uint32_t key, uint32_t* argb) {
  assert(key < p->hash_size_);
  if (p->data_[key] != kNotInitialized) {
    *argb = p->data_[key];
    return true;
  }
  return false;
}

struct PixelHasherLine {
  void Init(int xsize, int x_downsample_bits, int hash_bits) {
    int i;
    if (hash_bits == 0) {
      hash_bits = 1;
    }
    x_downsample_bits_ = x_downsample_bits;
    hash_bits_ = hash_bits;
    hashers_size_ =
        (xsize + (1 << x_downsample_bits) - 1) >> x_downsample_bits;

    hashers_ = (PixelHasher *)malloc(hashers_size_ * sizeof(hashers_[0]));
    for (i = 0; i < hashers_size_; ++i) {
      VP8LPixelHasherInit(&hashers_[i], hash_bits);
    }
  }
  void Delete() {
    int i;
    for (i = 0; i < hashers_size_; ++i) {
      VP8LPixelHasherDelete(&hashers_[i]);
    }
    free(hashers_);
  }
  void Insert(int x, uint32_t argb) {
    VP8LPixelHasherInsert(&hashers_[x >> x_downsample_bits_], argb);
  }
  uint32_t GetIndex(uint32_t argb) const {
    uint32_t val = kHashMul * argb;
    val >>= 32 - hash_bits_;
    return val;
  }
  bool Contains(int x, uint32_t argb) {
    int i;
    int ix = x >> x_downsample_bits_;
    if (VP8LPixelHasherContains(&hashers_[ix], argb)) {
      return true;
    }
    if (VP8LPixelHasherIsInitialized(&hashers_[ix], argb)) {
      return false;
    }
    for (i = 1; i < hashers_size_; ++i) {
      if (ix - i >= 0) {
        if (VP8LPixelHasherContains(&hashers_[ix - i], argb)) {
          return true;
        }
        if (VP8LPixelHasherIsInitialized(&hashers_[ix - i], argb)) {
          return false;
        }
      }
      if (ix + i < hashers_size_) {
        if (VP8LPixelHasherContains(&hashers_[ix + i], argb)) {
          return true;
        }
        if (VP8LPixelHasherIsInitialized(&hashers_[ix + i], argb)) {
          return false;
        }
      }
    }
    return false;
  }
  uint32_t Lookup(int x, uint32_t hash) {
    int i;
    int ix = x >> x_downsample_bits_;
    uint32_t argb;
    if (VP8LPixelHasherLookup(&hashers_[ix], hash, &argb)) {
      return argb;
    }
    for (i = 1; i < hashers_size_; ++i) {
      if (ix - i >= 0 &&
          VP8LPixelHasherLookup(&hashers_[ix - i], hash, &argb)) {
        return argb;
      }
      if (ix + i < hashers_size_ &&
          VP8LPixelHasherLookup(&hashers_[ix + i], hash, &argb)) {
        return argb;
      }
    }
    return 0;
  }
  PixelHasher* hashers_;
  int hash_bits_;
  int x_downsample_bits_;
  int hashers_size_;
};

#endif  // WEBP_PIXEL_HASHER_H_
