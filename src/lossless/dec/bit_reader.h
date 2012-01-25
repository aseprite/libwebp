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
  const uint8* buf_;
  size_t       len_;
  size_t       pos_;
  uint64       val_;
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
    ShiftBytes(br);
  }
  return val;
}

inline uint32 ReadOneBit(BitReader* const br) {
  const uint32 val = (br->val_ >> br->bit_pos_) & 1;
  ++br->bit_pos_;
  if (br->bit_pos_ == 40) {
    ShiftBytes(br);
  }
  return val;
}

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif  // WEBP_LOSSLESS_DEC_BIT_READER_H_
