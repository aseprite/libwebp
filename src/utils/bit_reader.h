// Copyright 2010 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
// Boolean decoder
//
// Author: Skal (pascal.massimino@gmail.com)
//         Vikas Arora (vikaas.arora@gmail.com)

#ifndef WEBP_UTILS_BIT_READER_H_
#define WEBP_UTILS_BIT_READER_H_

#include <assert.h>
#ifdef _MSC_VER
#include <stdlib.h>  // _byteswap_ulong
#endif
#include "../webp/types.h"

#ifdef __cplusplus
extern "C" {
#endif

// The Boolean decoder needs to maintain infinite precision on the value_ field.
// But since range_ is only 8bit, we only need an active window of 8 bits
// for value_. Left bits (MSB) gets zeroed and shifted away when value_ falls
// below 128, range_ is updated, and fresh bits read from the bitstream are
// brought in as LSB. To avoid reading the fresh bits one by one (slow), we
// cache BITS of them ahead. The total of (BITS + 8) bits must fit into a
// natural register (with type bit_t). To fetch BITS bits from bitstream we
// use a type lbit_t.
//
// value_ contains 8 active bits left-justified, followed by at most BITS
// fresh unused bits
//
// BITS can be any multiple of 8 from 8 to 56 (inclusive).
// Pick values that fit natural register size.

#if !defined(WEBP_REFERENCE_IMPLEMENTATION)

#if defined(__i386__) || defined(_M_IX86)      // x86 32bit
#define BITS 16
#elif defined(__x86_64__) || defined(_M_X64)   // x86 64bit
#define BITS 56
#elif defined(__arm__) || defined(_M_ARM)      // ARM
#define BITS 24
#elif defined(__mips__)                        // MIPS
#define BITS 24
#else                      // reasonable default
#define BITS 24
#endif

#else     // reference choices

#define BITS 8

#endif

// some endian fix (e.g.: mips-gcc doesn't define __BIG_ENDIAN__)
#if !defined(__BIG_ENDIAN__) && defined(__BYTE_ORDER__) && \
    (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
#define __BIG_ENDIAN__
#endif

//------------------------------------------------------------------------------
// Derived types and constants

// bit_t = natural register type
// lbit_t = natural type for memory I/O

#if (BITS > 32)
typedef uint64_t bit_t;
typedef uint64_t lbit_t;
#elif (BITS == 32)
typedef uint64_t bit_t;
typedef uint32_t lbit_t;
#elif (BITS == 24)
typedef uint32_t bit_t;
typedef uint32_t lbit_t;
#elif (BITS == 16)
typedef uint32_t bit_t;
typedef uint16_t lbit_t;
#else
typedef uint32_t bit_t;
typedef uint8_t lbit_t;
#endif

typedef uint32_t range_t;  // range_ only uses 8bits here. No need for bit_t.

//------------------------------------------------------------------------------
// Bitreader

typedef struct VP8BitReader VP8BitReader;
struct VP8BitReader {
  const uint8_t* buf_;        // next byte to be read
  const uint8_t* buf_end_;    // end of read buffer
  int eof_;                   // true if input is exhausted

  // boolean decoder
  range_t range_;            // current range minus 1. In [127, 254] interval.
  bit_t value_;              // current value
  int bits_;                 // number of valid bits left
};

// Initialize the bit reader and the boolean decoder.
void VP8InitBitReader(VP8BitReader* const br,
                      const uint8_t* const start, const uint8_t* const end);

// Update internal pointers to displace the byte buffer by the
// relative offset 'offset'.
void VP8RemapBitReader(VP8BitReader* const br, ptrdiff_t offset);

// return the next value made of 'num_bits' bits
uint32_t VP8GetValue(VP8BitReader* const br, int num_bits);
static WEBP_INLINE uint32_t VP8Get(VP8BitReader* const br) {
  return VP8GetValue(br, 1);
}

// return the next value with sign-extension.
int32_t VP8GetSignedValue(VP8BitReader* const br, int num_bits);

// Read a bit with proba 'prob'. Speed-critical function!
extern const uint8_t kVP8Log2Range[128];
extern const range_t kVP8NewRange[128];

void VP8LoadFreshBytes(VP8BitReader* const br);  // replenish br->bits_
void VP8LoadFinalBytes(VP8BitReader* const br);  // special case for the tail

static WEBP_INLINE void VP8LoadNewBytes(VP8BitReader* const br) {
  assert(br != NULL && br->buf_ != NULL);
  if (br->bits_ < 0) {
    if (br->buf_ + sizeof(lbit_t) <= br->buf_end_) {
      VP8LoadFreshBytes(br);
    } else {
      VP8LoadFinalBytes(br);
    }
  }
}

static WEBP_INLINE int VP8BitUpdate(VP8BitReader* const br, range_t split) {
  // Make sure we have a least BITS bits in 'value_'
  VP8LoadNewBytes(br);
  {
    const range_t value = (range_t)(br->value_ >> BITS);
    if (value > split) {
      br->range_ -= split + 1;
      br->value_ -= (bit_t)(split + 1) << BITS;
      return 1;
    } else {
      br->range_ = split;
      return 0;
    }
  }
}

static WEBP_INLINE int VP8GetBit(VP8BitReader* const br, int prob) {
  const range_t split = (br->range_ * prob) >> 8;
  const int bit = VP8BitUpdate(br, split);
  // range_ is in [0..127] interval here.
  const range_t idx = br->range_;
  if (idx <= 126) {
    const int shift = kVP8Log2Range[idx];
    br->range_ = kVP8NewRange[idx];
    br->value_ <<= shift;
    br->bits_ -= shift;
  }
  return bit;
}

static WEBP_INLINE int VP8GetSigned(VP8BitReader* const br, int v) {
  VP8LoadNewBytes(br);
  // simplified version of GetBit() for prob=0x80
  {
    const range_t split = br->range_ >> 1;
    const range_t value = (range_t)(br->value_ >> BITS);
    if (value > split) {
      br->range_ -= 1;     // r - (r>>1) - 1 = (r-1) >> 1
      br->value_ -= (bit_t)(split + 1) << BITS;
      v = -v;
    }
    // shift is always 1 here
    br->value_ <<= 1;
    br->bits_ -= 1;
    br->range_ |= 1;
    return v;
  }
}

// -----------------------------------------------------------------------------
// Bitreader for lossless format

// maximum number of bits (inclusive) the bit-reader can handle:
#define VP8L_MAX_NUM_BIT_READ 24

typedef uint64_t vp8l_val_t;  // right now, this bit-reader can only use 64bit.

typedef struct {
  vp8l_val_t     val_;        // pre-fetched bits
  const uint8_t* buf_;        // input byte buffer
  size_t         len_;        // buffer length
  size_t         pos_;        // byte position in buf_
  int            bit_pos_;    // current bit-reading position in val_
  int            eos_;        // bitstream is finished
  int            error_;      // an error occurred (buffer overflow attempt...)
} VP8LBitReader;

void VP8LInitBitReader(VP8LBitReader* const br,
                       const uint8_t* const start,
                       size_t length);

//  Sets a new data buffer.
void VP8LBitReaderSetBuffer(VP8LBitReader* const br,
                            const uint8_t* const buffer, size_t length);

// Reads the specified number of bits from read buffer.
// Flags an error in case end_of_stream or n_bits is more than the allowed limit
// of VP8L_MAX_NUM_BIT_READ (inclusive).
// Flags eos_ if this read attempt is going to cross the read buffer.
uint32_t VP8LReadBits(VP8LBitReader* const br, int n_bits);

// Return the prefetched bits, so they can be looked up.
static WEBP_INLINE uint32_t VP8LPrefetchBits(VP8LBitReader* const br) {
  return (uint32_t)(br->val_ >> br->bit_pos_);
}

// For jumping over a number of bits in the bit stream when accessed with
// VP8LPrefetchBits and VP8LFillBitWindow.
static WEBP_INLINE void VP8LSetBitPos(VP8LBitReader* const br, int val) {
  br->bit_pos_ = val;
}

// Advances the read buffer by 4 bytes to make room for reading next 32 bits.
void VP8LFillBitWindow(VP8LBitReader* const br);

#ifdef __cplusplus
}    // extern "C"
#endif

#endif  /* WEBP_UTILS_BIT_READER_H_ */
