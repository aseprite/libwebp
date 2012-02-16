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

#include <string.h>   // for memcpy()
#include "../common/integral_types.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#define UNALIGNED_LOAD32(_p) (*(const uint32_t *)(_p))
#define UNALIGNED_STORE32(_p, _val) (*(uint32_t *)(_p) = (_val))

#define TAG_SIZE 4
#define CHUNK_HEADER_SIZE 8
#define RIFF_HEADER_SIZE 12
#define HEADER_SIZE      (RIFF_HEADER_SIZE + CHUNK_HEADER_SIZE)
#define SIGNATURE_SIZE   1

typedef struct {
  uint8_t* buf_;
  size_t bit_pos_;
  size_t max_bytes_;

  // After all bits are written, the caller must observe the state of
  // error_. A value of 1 indicates that a memory allocation failure
  // has happened during bit writing.
  char error_;
} BitWriter;

static inline size_t BitWriterNumBytes(BitWriter* const bw) {
  return (bw->bit_pos_ + 7) >> 3;
}

// Returns 1 on success.
static int BitWriterResize(BitWriter* const bw, size_t extra_size) {
  uint8_t* allocated_buf;
  size_t allocated_size;
  const size_t size_required = BitWriterNumBytes(bw) + extra_size;
  if ((bw->max_bytes_ > 0) && (size_required <= bw->max_bytes_)) return 1;
  allocated_size = (3 * bw->max_bytes_) >> 1;
  if (allocated_size < size_required) {
    allocated_size = size_required;
  }
  // Make Allocated size multiple of KBs
  allocated_size = (((allocated_size >> 10) + 1) << 10);
  allocated_buf = (uint8_t*)malloc(allocated_size);
  if (allocated_buf == NULL) return 0;
  memset(allocated_buf, 0, allocated_size);
  if (bw->bit_pos_ > 0) {
    memcpy(allocated_buf, bw->buf_, BitWriterNumBytes(bw));
  }
  free(bw->buf_);
  bw->buf_ = allocated_buf;
  bw->max_bytes_ = allocated_size;
  return 1;
}

// Returns 1 on success.
static int BitWriterInit(BitWriter* const bw, size_t expected_size) {
  memset(bw, 0, sizeof(*bw));
  return BitWriterResize(bw, expected_size);
}

static inline uint8_t* BitWriterFinish(BitWriter* const bw) {
  return bw->buf_;
}

static void BitWriterDestroy(BitWriter* const bw) {
  assert(bw != NULL);
  free(bw->buf_);
  memset(bw, 0, sizeof(*bw));
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
//
// returns 1 on success.
static inline void WriteBits(int n_bits, uint32_t bits, BitWriter* const bw) {
  if (n_bits < 1) return;
#ifdef LITTLE_ENDIAN
  // Technically, this branch of the code can write up to 25 bits at a time,
  // but in deflate, the maximum number of bits written is 16 at a time.
  {
    uint8_t *p = &bw->buf_[bw->bit_pos_ >> 3];
    uint32_t v = UNALIGNED_LOAD32(p);
    v |= bits << (bw->bit_pos_ & 7);
    UNALIGNED_STORE32(p, v);  // Set some bits.
    bw->bit_pos_ += n_bits;
  }
#else
  // implicit & 0xff is assumed for uint8_t arithmetics
  {
    uint8_t *p = &bw->buf_[bw->bit_pos_ >> 3];
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
  }
#endif
  if ((bw->bit_pos_ >> 3) > (bw->max_bytes_ - 8)) {
    const size_t kAdditionalBuffer = 32768 + bw->max_bytes_;
    if (!BitWriterResize(bw, kAdditionalBuffer)) {
      bw->bit_pos_ = 0;
      bw->error_ = 1;
    }
  }
}

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif  // WEBP_LOSSLESS_ENC_BIT_WRITER_H_
