// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: szabadka@google.com (Zoltan Szabadka)

#include "bit_stream.h"
#include "../common/integral_types.h"

BitStream::BitStream(int length, const uint8* stream)
    : length_(length), stream_(stream),
      byte_position_(0), bit_position_(0), window_(0) {
  // Read first four bytes to window_
  for (int i = 0; i < 4 && i < length_; ++i) {
    window_ += stream_[i] << (8 * i);
  }
}

void BitStream::Shift(int num_bytes) {
  for (int i = 0; i< num_bytes; ++i) {
    ShiftOneByte();
  }
}

void BitStream::ShiftOneByte() {
  window_ = window_ >> 8;
  if (byte_position_ + 4 < length_) {
    window_ += stream_[byte_position_ + 4] << 24;
  }
  ++byte_position_;
}

uint32 BitStream::Read(int num_bits) {
  if (num_bits == 0) { return 0; }

  VERIFY(!end());
  VERIFY(num_bits <= 25);
  VERIFY(num_bits > 0);

  uint32 result = (window_ >> bit_position_) & ((1 << num_bits) - 1);

  bit_position_ += num_bits;
  Shift(bit_position_ >> 3);
  bit_position_ &= 7;

  return result;
}
