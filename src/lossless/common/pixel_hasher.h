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

#include <stddef.h>
#include <stdlib.h>

#include "integral_types.h"

static const int kNotInitialized = 0x00e1c300;
static const uint32 kHashMul = 0x1e35a7bd;

class PixelHasher {
 public:
  void Init(int hash_bits) {
    data_ = new uint32[1 << hash_bits];
    hash_bits_ = hash_bits;
    for (int i = 0; i < (1 << hash_bits_); ++i) {
      data_[i] = kNotInitialized;
    }
  }
  void Delete() {
    delete[] data_;
  }
  void insert(uint32 argb) {
    uint32 val = kHashMul * argb;
    val >>= 32 - hash_bits_;
    if (hash_bits_ == 0) {
      val = 0;
    }
    data_[val] = argb;
  }
  bool is_initialized(uint32 argb) const {
    uint32 val = kHashMul * argb;
    val >>= 32 - hash_bits_;
    if (hash_bits_ == 0) {
      val = 0;
    }
    return data_[val] != kNotInitialized;
  }
  bool contains(uint32 argb) const {
    uint32 val = kHashMul * argb;
    val >>= 32 - hash_bits_;
    if (hash_bits_ == 0) {
      val = 0;
    }
    if (argb == kNotInitialized) {
      return false;
    }
    return data_[val] == argb;
  }
  bool lookup(uint32 val, uint32* argb) {
    if (val < (1 << hash_bits_) && data_[val] != kNotInitialized) {
      *argb = data_[val];
      return true;
    }
    return false;
  }

 private:
  uint32 *data_;
  uint32 hash_bits_;
};

class PixelHasherLine {
 public:
  PixelHasherLine(int xsize, int x_downsample_bits, int hash_bits) {
    x_downsample_bits_ = x_downsample_bits;
    hash_bits_ = hash_bits;
    hashers_size_ =
        (xsize + (1 << x_downsample_bits) - 1) >> x_downsample_bits;

    hashers_ = (PixelHasher *)malloc(hashers_size_ * sizeof(hashers_[0]));
    for (int i = 0; i < hashers_size_; ++i) {
      hashers_[i].Init(hash_bits);
    }
    InitDefault();
  }
  ~PixelHasherLine() {
    default_.Delete();
    for (int i = 0; i < hashers_size_; ++i) {
      hashers_[i].Delete();
    }
    free(hashers_);
  }
  void InitDefault() {
    // TODO(jyrki): This is ugly, make this initialization more beautiful.

    default_.Init(hash_bits_);
    // 12 bit palette
    for (int r = 0; r < 256; r += 17) {
      for (int g = 0; g < 256; g += 17) {
        for (int b = 0; b < 256; b += 17) {
          default_.insert(0xff000000 + (r << 16) + (g << 8) + b);
        }
      }
    }
    // 6 bit palette
    for (int r = 0; r < 256; r += 85) {
      for (int g = 0; g < 256; g += 85) {
        for (int b = 0; b < 256; b += 85) {
          default_.insert(0xff000000 + (r << 16) + (g << 8) + b);
        }
      }
    }
    // 3 bit palette
    for (int r = 0; r < 256; r += 255) {
      for (int g = 0; g < 256; g += 255) {
        for (int b = 0; b < 256; b += 255) {
          default_.insert(0xff000000 + (r << 16) + (g << 8) + b);
        }
      }
    }
    // Around zero for predictor mode
    // TODO(jyrki): this should only be added if there is a spatial prediction
    // image.
    for (int i0 = -3; i0 < 4; ++i0) {
      for (int i1 = -3; i1 < 4; ++i1) {
        for (int i2 = -3; i2 < 4; ++i2) {
          default_.insert(((i0 & 0xff) << 24) +
                          ((i1 & 0xff) << 16) +
                          (i2 & 0xff));
        }
      }
    }
    // Green is special.
    default_.insert(0xff001f00);
    default_.insert(0xff003f00);
    default_.insert(0xff005f00);
    default_.insert(0xff007f00);
    default_.insert(0xff009f00);
    default_.insert(0xff00bf00);
    default_.insert(0xff00df00);
    default_.insert(0xff00ff00);
    // Fully transparent, black, and, white.
    default_.insert(0x00000000);
    default_.insert(0xff000000);
    default_.insert(0xffffffff);
  }
  void Insert(int x, uint32 argb) {
    hashers_[x >> x_downsample_bits_].insert(argb);
  }
  uint32 GetIndex(uint32 argb) const {
    uint32 val = kHashMul * argb;
    val >>= 32 - hash_bits_;
    if (hash_bits_ == 0) {
      return 0;
    }
    return val;
  }
  bool Contains(int x, uint32 argb) {
    int ix = x >> x_downsample_bits_;
    if (hashers_[ix].contains(argb)) {
      return true;
    }
    if (hashers_[ix].is_initialized(argb)) {
      return false;
    }
    for (int i = 1; i < hashers_size_; ++i) {
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
    return default_.contains(argb);
  }
  uint32 Lookup(int x, uint32 hash) {
    int ix = x >> x_downsample_bits_;
    uint32 argb;
    if (hashers_[ix].lookup(hash, &argb)) {
      return argb;
    }
    for (int i = 1; i < hashers_size_; ++i) {
      if (ix - i >= 0 &&
          hashers_[ix - i].lookup(hash, &argb)) {
        return argb;
      }
      if (ix + i < hashers_size_ &&
          hashers_[ix + i].lookup(hash, &argb)) {
        return argb;
      }
    }
    if (default_.lookup(hash, &argb)) {
      return argb;
    }
    return 0;
  }
  PixelHasher* hashers_;
  PixelHasher default_;
  int hash_bits_;
  int x_downsample_bits_;
  int hashers_size_;
};

#endif  // WEBP_PIXEL_HASHER_H_
