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

#include "./bit_reader.h"

//------------------------------------------------------------------------------
// VP8BitReader

void VP8InitBitReader(VP8BitReader* const br,
                      const uint8_t* const start, const uint8_t* const end) {
  assert(br != NULL);
  assert(start != NULL);
  assert(start <= end);
  br->range_   = 255 - 1;
  br->buf_     = start;
  br->buf_end_ = end;
  br->value_   = 0;
  br->bits_    = -8;   // to load the very first 8bits
  br->eof_     = 0;
  VP8LoadNewBytes(br);
}

void VP8RemapBitReader(VP8BitReader* const br, ptrdiff_t offset) {
  if (br->buf_ != NULL) {
    br->buf_ += offset;
    br->buf_end_ += offset;
  }
}

const uint8_t kVP8Log2Range[128] = {
     7, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  0
};

// range = ((range + 1) << kVP8Log2Range[range]) - 1
const range_t kVP8NewRange[128] = {
  127, 127, 191, 127, 159, 191, 223, 127,
  143, 159, 175, 191, 207, 223, 239, 127,
  135, 143, 151, 159, 167, 175, 183, 191,
  199, 207, 215, 223, 231, 239, 247, 127,
  131, 135, 139, 143, 147, 151, 155, 159,
  163, 167, 171, 175, 179, 183, 187, 191,
  195, 199, 203, 207, 211, 215, 219, 223,
  227, 231, 235, 239, 243, 247, 251, 127,
  129, 131, 133, 135, 137, 139, 141, 143,
  145, 147, 149, 151, 153, 155, 157, 159,
  161, 163, 165, 167, 169, 171, 173, 175,
  177, 179, 181, 183, 185, 187, 189, 191,
  193, 195, 197, 199, 201, 203, 205, 207,
  209, 211, 213, 215, 217, 219, 221, 223,
  225, 227, 229, 231, 233, 235, 237, 239,
  241, 243, 245, 247, 249, 251, 253, 127
};

void VP8LoadFinalBytes(VP8BitReader* const br) {
  assert(br != NULL && br->buf_ != NULL && br->bits_ < 0);
  // Only read 8bits at a time
  if (br->buf_ < br->buf_end_) {
    br->bits_ += 8;
    br->value_ |= (bit_t)(*br->buf_++) << ((BITS) - br->bits_);
  } else if (!br->eof_) {
    br->eof_ = 1;
  }
}

void VP8LoadFreshBytes(VP8BitReader* const br) {
  assert(br != NULL && br->buf_ != NULL && br->bits_ < 0);
  // Read 'BITS' bits at a time if possible.
  assert(br->buf_ + sizeof(lbit_t) <= br->buf_end_);
  {
    // convert memory type to register type (with some zero'ing!)
    bit_t bits;
#if defined(__mips__)                          // MIPS
    // This is needed because of un-aligned read.
    lbit_t in_bits;
    lbit_t* p_buf_ = (lbit_t*)br->buf_;
    __asm__ volatile(
      ".set   push                             \n\t"
      ".set   at                               \n\t"
      ".set   macro                            \n\t"
      "ulw    %[in_bits], 0(%[p_buf_])         \n\t"
      ".set   pop                              \n\t"
      : [in_bits]"=r"(in_bits)
      : [p_buf_]"r"(p_buf_)
      : "memory", "at"
    );
#else
    const lbit_t in_bits = *(const lbit_t*)br->buf_;
#endif
    br->buf_ += (BITS) >> 3;
#if !defined(__BIG_ENDIAN__)
#if (BITS > 32)
// gcc 4.3 has builtin functions for swap32/swap64
#if defined(__GNUC__) && \
           (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 3))
    bits = (bit_t)__builtin_bswap64(in_bits);
#elif defined(_MSC_VER)
    bits = (bit_t)_byteswap_uint64(in_bits);
#elif defined(__x86_64__)
    __asm__ volatile("bswapq %0" : "=r"(bits) : "0"(in_bits));
#else  // generic code for swapping 64-bit values (suggested by bdb@)
    bits = (bit_t)in_bits;
    bits = ((bits & 0xffffffff00000000ull) >> 32) |
           ((bits & 0x00000000ffffffffull) << 32);
    bits = ((bits & 0xffff0000ffff0000ull) >> 16) |
           ((bits & 0x0000ffff0000ffffull) << 16);
    bits = ((bits & 0xff00ff00ff00ff00ull) >> 8) |
           ((bits & 0x00ff00ff00ff00ffull) << 8);
#endif
    bits >>= 64 - BITS;
#elif (BITS >= 24)
#if defined(__i386__) || defined(__x86_64__)
    {
      lbit_t swapped_in_bits;
      __asm__ volatile("bswap %k0" : "=r"(swapped_in_bits) : "0"(in_bits));
      bits = (bit_t)swapped_in_bits;   // 24b/32b -> 32b/64b zero-extension
    }
#elif defined(_MSC_VER)
    bits = (bit_t)_byteswap_ulong(in_bits);
#else
    bits = (bit_t)(in_bits >> 24) | ((in_bits >> 8) & 0xff00)
         | ((in_bits << 8) & 0xff0000)  | (in_bits << 24);
#endif  // x86
    bits >>= 32 - BITS;
#elif (BITS == 16)
    // gcc will recognize a 'rorw $8, ...' here:
    bits = (bit_t)(in_bits >> 8) | ((in_bits & 0xff) << 8);
#else   // BITS == 8
    bits = (bit_t)in_bits;
#endif
#else    // BIG_ENDIAN
    bits = (bit_t)in_bits;
    bits >>= (8 * sizeof(bit_t) - BITS);
#endif
    br->value_ |= bits << (-br->bits_);
    br->bits_ += (BITS);
  }
}

//------------------------------------------------------------------------------
// Higher-level calls

uint32_t VP8GetValue(VP8BitReader* const br, int bits) {
  uint32_t v = 0;
  while (bits-- > 0) {
    v |= VP8GetBit(br, 0x80) << bits;
  }
  return v;
}

int32_t VP8GetSignedValue(VP8BitReader* const br, int bits) {
  const int value = VP8GetValue(br, bits);
  return VP8Get(br) ? -value : value;
}

//------------------------------------------------------------------------------
// VP8LBitReader

#define LBITS 64      // Number of bits prefetched.
#define WBITS 32      // Minimum number of bytes needed after VP8LFillBitWindow.
#define LOG8_WBITS 4  // Number of bytes needed to store WBITS bits.

static const uint32_t kBitMask[VP8L_MAX_NUM_BIT_READ + 1] = {
  0, 1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047, 4095, 8191, 16383, 32767,
  65535, 131071, 262143, 524287, 1048575, 2097151, 4194303, 8388607, 16777215
};

void VP8LInitBitReader(VP8LBitReader* const br,
                       const uint8_t* const start,
                       size_t length) {
  size_t i;
  assert(br != NULL);
  assert(start != NULL);
  assert(length < 0xfffffff8u);   // can't happen with a RIFF chunk.

  br->buf_ = start;
  br->len_ = length;
  br->val_ = 0;
  br->pos_ = 0;
  br->bit_pos_ = 0;
  br->eos_ = 0;
  br->error_ = 0;
  for (i = 0; i < sizeof(br->val_) && i < br->len_; ++i) {
    br->val_ |= ((vp8l_val_t)br->buf_[br->pos_]) << (8 * i);
    ++br->pos_;
  }
}

// Special version that assumes br->pos_ <= br_len_.
static int IsEndOfStreamSpecial(const VP8LBitReader* const br) {
  assert(br->pos_ <= br->len_);
  return br->pos_ == br->len_ && br->bit_pos_ >= LBITS;
}

static int IsEndOfStream(const VP8LBitReader* const br) {
  return (br->pos_ > br->len_) || IsEndOfStreamSpecial(br);
}

void VP8LBitReaderSetBuffer(VP8LBitReader* const br,
                            const uint8_t* const buf, size_t len) {
  assert(br != NULL);
  assert(buf != NULL);
  assert(len < 0xfffffff8u);   // can't happen with a RIFF chunk.
  br->buf_ = buf;
  br->len_ = len;
  br->eos_ = IsEndOfStream(br);
}

// If not at EOS, reload up to LBITS byte-by-byte
static void ShiftBytes(VP8LBitReader* const br) {
  while (br->bit_pos_ >= 8 && br->pos_ < br->len_) {
    br->val_ >>= 8;
    br->val_ |= ((vp8l_val_t)br->buf_[br->pos_]) << (LBITS - 8);
    ++br->pos_;
    br->bit_pos_ -= 8;
  }
}

void VP8LFillBitWindow(VP8LBitReader* const br) {
  if (br->bit_pos_ >= WBITS) {
#if (defined(__x86_64__) || defined(_M_X64))
    if (br->pos_ + sizeof(br->val_) < br->len_) {
      br->val_ >>= WBITS;
      br->bit_pos_ -= WBITS;
      // The expression below needs a little-endian arch to work correctly.
      // This gives a large speedup for decoding speed.
      br->val_ |= *(const vp8l_val_t*)(br->buf_ + br->pos_) << (LBITS - WBITS);
      br->pos_ += LOG8_WBITS;
      return;
    }
#endif
    ShiftBytes(br);       // Slow path.
    br->eos_ = IsEndOfStreamSpecial(br);
  }
}

uint32_t VP8LReadBits(VP8LBitReader* const br, int n_bits) {
  assert(n_bits >= 0);
  // Flag an error if end_of_stream or n_bits is more than allowed limit.
  if (!br->eos_ && n_bits <= VP8L_MAX_NUM_BIT_READ) {
    const uint32_t val =
        (uint32_t)(br->val_ >> br->bit_pos_) & kBitMask[n_bits];
    const int new_bits = br->bit_pos_ + n_bits;
    br->bit_pos_ = new_bits;
    // If this read is going to cross the read buffer, set the eos flag.
    br->eos_ = IsEndOfStreamSpecial(br);
    ShiftBytes(br);
    return val;
  } else {
    br->error_ = 1;
    return 0;
  }
}

//------------------------------------------------------------------------------

