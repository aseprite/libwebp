// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)
//
// Bit writing helpers

#ifndef WEBP_LOSSLESS_ENC_BIT_WRITER_H_
#define WEBP_LOSSLESS_ENC_BIT_WRITER_H_

#ifdef __linux__
#include <endian.h>
#else
#include <machine/endian.h>  // Gross overassumption.
#endif

#include "./backward_distance.h"
#include "../common/integral_types.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#define UNALIGNED_LOAD32(_p) (*reinterpret_cast<const uint32 *>(_p))
#define UNALIGNED_STORE32(_p, _val) (*reinterpret_cast<uint32 *>(_p) = (_val))

#define TAG_SIZE 4
#define CHUNK_HEADER_SIZE 8
#define RIFF_HEADER_SIZE 12
#define HEADER_SIZE      (RIFF_HEADER_SIZE + CHUNK_HEADER_SIZE)
#define SIGNATURE_SIZE   1

typedef struct {
  uint8* buf_;
  size_t bit_pos_;
} BitWriter;

void BitWriterInit(BitWriter* bw, uint8* buf, size_t size) {
  bw->buf_ = buf;
  bw->bit_pos_ = 0;
  UNALIGNED_STORE32(bw->buf_, 0);
}

inline size_t BitWriterNumBytes(BitWriter* bw) {
  return (bw->bit_pos_ + 7) >> 3;
}

// This function writes bits into bytes in increasing addresses, and within
// a byte least-significant-bit first.
//
// The function can write up to 16 bits in one go with WriteBits
// Example: let's assume that 3 bits (Rs below) have been written already:
//
// BYTE-0     BYTE+1       BYTE+2
//
// 0000 0RRR    0000 0000    0000 0000
//
// Now, we could write 5 or less bits in MSB by just sifting by 3
// and OR'ing to BYTE-0.
//
// For n bits, we take the last 5 bytes, OR that with high bits in BYTE-0,
// and locate the rest in BYTE+1 and BYTE+2.
inline void WriteBits(const int n_bits, const uint32 bits, BitWriter* bw) {
#ifdef LITTLE_ENDIAN
  // Technically, this branch of the code can write up to 25 bits at a time,
  // but in deflate, the maximum number of bits written is 16 at a time.
  uint8 *p = &bw->buf_[bw->bit_pos_ >> 3];
  uint32 v = UNALIGNED_LOAD32(p);
  v |= bits << (bw->bit_pos_ & 7);
  UNALIGNED_STORE32(p, v);  // Set some bits.

  // This buffer could be memset before, but it would be expensive to know
  // how much to memset, and would require another passes through the memory.
  //
  // So, it is probably better to be cleared here.
  UNALIGNED_STORE32(p + 4, 0);  // Clear bits that will be or'ed later.
  // The first 4 bytes need to be cleared by the user of WriteBits.

  bw->bit_pos_ += n_bits;
#else
  // implicit & 0xff is assumed for uint8 arithmetics
  uint8 *p = &bw->buf_[bw->bit_pos_ >> 3];
  const int bits_reserved_in_first_byte = (bw->bit_pos_ & 7);
  *p++ |= (bits << bits_reserved_in_first_byte);
  const int bits_left_to_write = n_bits - 8 + bits_reserved_in_first_byte;
  if (bits_left_to_write >= 1) {
    *p++ = bits >> (8 - bits_reserved_in_first_byte);
    if (bits_left_to_write >= 9) {
      *p++ = bits >> (16 - bits_reserved_in_first_byte);
    }
  }
  *p = 0;
  bw->bit_pos_ += n_bits;
#endif
}

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif  // WEBP_LOSSLESS_ENC_BIT_WRITER_H_
