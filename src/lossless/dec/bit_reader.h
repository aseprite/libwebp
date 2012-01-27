// Copyright 2012 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: vikasa@google.com (Vikas Arora)
//
// Bit Reader for WebP Lossless

#ifndef WEBP_LOSSLESS_DEC_BIT_READER_H_
#define WEBP_LOSSLESS_DEC_BIT_READER_H_

#include <assert.h>

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

// -----------------------------------------------------------------------------
// Bitreader

#define MAX_NUM_BIT_READ 25
const uint32 kBitMask[MAX_NUM_BIT_READ] = {
  0, 1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047, 4095, 8191, 16383, 32767,
  65535, 131071, 262143, 524287, 1048575, 2097151, 4194303, 8388607, 16777215
};

typedef struct {
  uint64       val_;
  const uint8* buf_;
  size_t       len_;
  size_t       pos_;
  int          bit_pos_;
} BitReader;

void InitBitReader(BitReader* const br,
                   const uint8* const start,
                   size_t length) {
  int i;
  assert(br);
  assert(start);

  br->buf_ = start;
  br->len_ = length;
  br->val_ = 0;
  br->pos_ = 0;
  br->bit_pos_ = 0;
  for (i = 0; i < sizeof(br->val_) && i < br->len_; ++i) {
    br->val_ |= ((uint64)br->buf_[br->pos_]) << (8 * i);
    ++br->pos_;
  }
}

static void ShiftBytes(BitReader* const br) {
  while (br->bit_pos_ >= 8 && br->pos_ < br->len_) {
    br->val_ >>= 8;
    br->val_ |= ((uint64)br->buf_[br->pos_]) << 56;
    ++br->pos_;
    br->bit_pos_ -= 8;
  }
}

inline uint32 ReadBits(BitReader* const br, int n_bits) {
  const uint32 val = (br->val_ >> br->bit_pos_) & kBitMask[n_bits];
  assert(n_bits < MAX_NUM_BIT_READ);
  assert(n_bits >= 0);
  br->bit_pos_ += n_bits;
  if (br->bit_pos_ >= 40) {
    if (br->pos_ < br->len_ - 5) {
      br->val_ >>= 40;
      br->val_ |=
          (((uint64)br->buf_[br->pos_ + 0]) << 24) |
          (((uint64)br->buf_[br->pos_ + 1]) << 32) |
          (((uint64)br->buf_[br->pos_ + 2]) << 40) |
          (((uint64)br->buf_[br->pos_ + 3]) << 48) |
          (((uint64)br->buf_[br->pos_ + 4]) << 56);
      br->pos_ += 5;
      br->bit_pos_ -= 40;
    }
    if (br->bit_pos_ >= 8) {
      ShiftBytes(br);
    }
  }
  return val;
}

inline void FillBitWindow(BitReader* const br) {
  if (br->bit_pos_ >= 32) {
#if defined(__x86_64__)
    if (br->pos_ < br->len_ - 8) {
      br->val_ >>= 32;
      // The expression below needs a little-endian arch to work correctly.
      // This gives a large speedup for decoding speed.
      br->val_ |= *(const uint64 *)(br->buf_ + br->pos_) << 32;
      br->pos_ += 4;
      br->bit_pos_ -= 32;
    } else {
      // Slow path.
      ShiftBytes(br);
    }
#else
    // Always the slow path.
    ShiftBytes(br);
#endif
  }
}

// If ReadOneBitUnsafe is much faster than ReadBits(..., 1), but it can be
// called only 32 times after the last FillBitWindow, and consequent calls
// may return invalid data.
inline uint32 ReadOneBitUnsafe(BitReader* const br) {
  const uint32 val = (br->val_ >> br->bit_pos_) & 1;
  assert(br->bit_pos_ < 64);
  ++br->bit_pos_;
  return val;
}

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif  // WEBP_LOSSLESS_DEC_BIT_READER_H_
