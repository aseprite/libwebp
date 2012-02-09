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
static const uint32_t kNotInitialized = kHashMul;

struct PixelHasher {
  void Init(int hash_bits) {
    int i;
    hash_shift_ = 32 - hash_bits;
    hash_size_ = 1 << hash_bits;
    data_ = (uint32_t *)malloc(hash_size_ * sizeof(data_[0]));
    for (i = 0; i < hash_size_; ++i) {
      data_[i] = kNotInitialized;
    }
  }
  void Delete() {
    free(data_);
  }
  void insert(uint32_t argb) {
    const uint32_t key = (kHashMul * argb) >> hash_shift_;
    data_[key] = argb;
  }
  bool is_initialized(uint32_t argb) const {
    const uint32_t key = (kHashMul * argb) >> hash_shift_;
    return data_[key] != kNotInitialized;
  }
  bool contains(uint32_t argb) const {
    uint32_t key;
    if (argb == kNotInitialized) {
      return false;
    }
    key = (kHashMul * argb) >> hash_shift_;
    return data_[key] == argb;
  }
  bool lookup(uint32_t key, uint32_t* argb) {
    assert(key < hash_size_);
    if (data_[key] != kNotInitialized) {
      *argb = data_[key];
      return true;
    }
    return false;
  }

  uint32_t *data_;
  uint32_t hash_shift_;
  uint32_t hash_size_;
};

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
      hashers_[i].Init(hash_bits);
    }
  }
  void Delete() {
    int i;
    for (i = 0; i < hashers_size_; ++i) {
      hashers_[i].Delete();
    }
    free(hashers_);
  }
  void Insert(int x, uint32_t argb) {
    hashers_[x >> x_downsample_bits_].insert(argb);
  }
  uint32_t GetIndex(uint32_t argb) const {
    uint32_t val = kHashMul * argb;
    val >>= 32 - hash_bits_;
    return val;
  }
  bool Contains(int x, uint32_t argb) {
    int i;
    int ix = x >> x_downsample_bits_;
    if (hashers_[ix].contains(argb)) {
      return true;
    }
    if (hashers_[ix].is_initialized(argb)) {
      return false;
    }
    for (i = 1; i < hashers_size_; ++i) {
      if (ix - i >= 0) {
        if (hashers_[ix - i].contains(argb)) {
          return true;
        }
        if (hashers_[ix - i].is_initialized(argb)) {
          return false;
        }
      }
      if (ix + i < hashers_size_) {
        if (hashers_[ix + i].contains(argb)) {
          return true;
        }
        if (hashers_[ix + i].is_initialized(argb)) {
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
    if (hashers_[ix].lookup(hash, &argb)) {
      return argb;
    }
    for (i = 1; i < hashers_size_; ++i) {
      if (ix - i >= 0 &&
          hashers_[ix - i].lookup(hash, &argb)) {
        return argb;
      }
      if (ix + i < hashers_size_ &&
          hashers_[ix + i].lookup(hash, &argb)) {
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
