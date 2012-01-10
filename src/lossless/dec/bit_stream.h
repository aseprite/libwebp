// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: szabadka@google.com (Zoltan Szabadka)
//
// Helper class to read bits.

#ifndef WEBP_DEC_BIT_STREAM_H_
#define WEBP_DEC_BIT_STREAM_H_

#include "../common/integral_types.h"

class BitStream {
 public:
  BitStream(int length, const uint8* stream);

  bool end() { return byte_position_ >= length_; }

  // Can read at most 25 bits at a time.
  uint32 Read(int num_bits);
  int ReadOneBit();

 private:
  inline void Shift(int num_bytes);
  inline void ShiftOneByte();

  const int length_;
  const uint8* stream_;

  // Current position in the stream.
  int byte_position_;
  int bit_position_;

  // The next four bytes stored.
  uint32 window_;
};

inline void BitStream::ShiftOneByte() {
  window_ = window_ >> 8;
  if (byte_position_ < length_) {
    window_ += stream_[byte_position_] << 24;
  } else {
    VERIFY(byte_position_ < length_ + 4);
  }
  ++byte_position_;
}

inline int BitStream::ReadOneBit() {
  uint32 result = (window_ >> bit_position_) & 1;
  ++bit_position_;
  if (bit_position_ >= 8) {
    ShiftOneByte();
    bit_position_ = 0;
  }
  return result;
}

#endif  // WEBP_DEC_BIT_STREAM_H_
