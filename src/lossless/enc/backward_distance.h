// Copyright 2011 Google Inc. All Rights Reserved.
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

#include <stdint.h>

// use GNU builtins where available
#if defined(__GNUC__) && \
    ((__GNUC__ == 3 && __GNUC_MINOR__ >= 4) || __GNUC__ >= 4)
inline int BitsLog2Floor(uint32_t n) {
  return n == 0 ? -1 : 31 ^ __builtin_clz(n);
}
#else
inline int BitsLog2Floor(uint32_t n) {
  int log;
  uint32_t value;
  int i;
  if (n == 0)
    return -1;
  log = 0;
  value = n;
  for (i = 4; i >= 0; --i) {
    int shift = (1 << i);
    uint32_t x = value >> shift;
    if (x != 0) {
      value = x;
      log += shift;
    }
  }
  return log;
}
#endif

inline int BitsLog2Ceiling(uint32_t n) {
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
static inline void PrefixEncode(
    int distance,
    int * __restrict code,
    int * __restrict extra_bits_count,
    int * __restrict extra_bits_value) {
  const int highest_bit = BitsLog2Floor(--distance);
  // & 0x3f is to make behavior well defined when highest_bit is -1 or 0.
  const int second_highest_bit =
      (distance >> ((highest_bit - 1) & 0x3f)) & 1;
  *extra_bits_count = (highest_bit > 0) ? highest_bit - 1 : 0;
  *extra_bits_value = distance & ((1 << *extra_bits_count) - 1);
  *code = (highest_bit > 0) ? 2 * highest_bit + second_highest_bit :
      (highest_bit == 0) ? 1 : 0;
}

#endif  // WEBP_BACKWARD_DISTANCE_H_
