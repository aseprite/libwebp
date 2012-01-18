// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: szabadka@google.com (Zoltan Szabadka)

#include <assert.h>
#include "bit_stream.h"
#include "../common/integral_types.h"

BitStream::BitStream(int length, const uint8* stream)
    : length_(length), stream_(stream),
      byte_position_(8), bit_position_(0), window_(0) {

  // Read first four bytes to window_
  for (int i = 0; i < 8 && i < length_; ++i) {
    window_ += (uint64)(stream_[i]) << (8 * i);
  }
}

void BitStream::Shift(int num_bytes) {
  for (int i = 0; i < num_bytes; ++i) {
    ShiftOneByte();
  }
}

uint32 BitStream::Read(int num_bits) {
  assert(num_bits <= 25);
  assert(num_bits >= 0);

  uint32 result = (window_ >> bit_position_) & ((1 << num_bits) - 1);

  bit_position_ += num_bits;
  Shift(bit_position_ >> 3);
  bit_position_ &= 7;

  return result;
}
