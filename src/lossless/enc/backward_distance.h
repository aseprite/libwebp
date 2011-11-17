// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)
//
// This class defines the encoding of backward distances into
// prefix codes, the amount of extra bits, and the actual
// values of the extra bits.

#ifndef WEBP_BACKWARD_DISTANCE_H_
#define WEBP_BACKWARD_DISTANCE_H_

#include "../common/integral_types.h"

// use GNU builtins where available
#if defined(__GNUC__) && \
    ((__GNUC__ == 3 && __GNUC_MINOR__ >= 4) || __GNUC__ >= 4)
inline int BitsLog2Floor(uint32 n) {
  return n == 0 ? -1 : 31 ^ __builtin_clz(n);
}
#else
inline int BitsLog2Floor(uint32 n) {
  if (n == 0)
    return -1;
  int log = 0;
  uint32 value = n;
  for (int i = 4; i >= 0; --i) {
    int shift = (1 << i);
    uint32 x = value >> shift;
    if (x != 0) {
      value = x;
      log += shift;
    }
  }
  return log;
}
#endif

inline int BitsLog2Ceiling(uint32 n) {
  int floor = BitsLog2Floor(n);
  if (n == (n &~ (n - 1)))  // zero or a power of two.
    return floor;
  else
    return floor + 1;
}


// Splitting of distance and length codes into prefixes and
// extra bits. The prefixes are encoded with an entropy code
// while the extra bits are stored just as normal bits.
//
// This code allows for large window sizes, up to 32 bits
// (16 GB for 4 byte pixels).
class BackwardDistance {
 public:
  struct DistanceCode {
    uint8 extra_bits;
    uint8 distance_code;
  };
  static inline void Encode(
      int distance,
      int * __restrict code,
      int * __restrict extra_bits_count,
      int * __restrict extra_bits_value) {
    --distance;
    const int highest_bit = BitsLog2Floor(distance);
    // & 0x3f is to make behavior well defined when highest_bit is -1 or 0.
    const int second_highest_bit_value =
        (distance >> ((highest_bit - 1) & 0x3f)) & 1;
    *extra_bits_count =
        distance_code_lut_offset_[2 * highest_bit].extra_bits;
    *extra_bits_value = distance & ((1 << *extra_bits_count) - 1);
    *code =
        distance_code_lut_offset_[2 * highest_bit +
                                  second_highest_bit_value].distance_code;
  }
 private:
  static const DistanceCode distance_code_lut_[66];
  static const DistanceCode *distance_code_lut_offset_;
};

// Same entropy code for length and distance.
typedef BackwardDistance BackwardLength;

#endif  // WEBP_BACKWARD_DISTANCE_H_
