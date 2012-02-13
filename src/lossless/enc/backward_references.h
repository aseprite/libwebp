// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)

#ifndef WEBP_BACKWARD_REFERENCES_H_
#define WEBP_BACKWARD_REFERENCES_H_

#include <assert.h>
#include <stdint.h>

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

// Backward reference distance prefix codes
#define DISTANCE_CODES_MAX 40

// Compression constants
#define CODE_LENGTH_CODES 19
static const int kRowHasherXSubsampling = 7;
static const int kLengthCodes = 24;
static const int kPaletteCodeBitsMax = 11;
#define PIX_OR_COPY_CODES_MAX (256 + 24 + (1 << 11))
static const int kMaxLength = 4096;

// use GNU builtins where available.
#if defined(__GNUC__) && \
    ((__GNUC__ == 3 && __GNUC_MINOR__ >= 4) || __GNUC__ >= 4)
static inline int BitsLog2Floor(uint32_t n) {
  return n == 0 ? -1 : 31 ^ __builtin_clz(n);
}
#else
static inline int BitsLog2Floor(uint32_t n) {
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

static inline int BitsLog2Ceiling(uint32_t n) {
  int floor = BitsLog2Floor(n);
  if (n == (n & ~(n - 1)))  // zero or a power of two.
    return floor;
  else
    return floor + 1;
}

// Splitting of distance and length codes into prefixes and
// extra bits. The prefixes are encoded with an entropy code
// while the extra bits are stored just as normal bits.
static inline void PrefixEncode(
    int distance,
    int *code,
    int *extra_bits_count,
    int *extra_bits_value) {
  // Collect the two most significant bits where the highest bit is 1.
  const int highest_bit = BitsLog2Floor(--distance);
  // & 0x3f is to make behavior well defined when highest_bit
  // does not exist or is the least significant bit.
  const int second_highest_bit =
      (distance >> ((highest_bit - 1) & 0x3f)) & 1;
  *extra_bits_count = (highest_bit > 0) ? highest_bit - 1 : 0;
  *extra_bits_value = distance & ((1 << *extra_bits_count) - 1);
  *code = (highest_bit > 0) ? 2 * highest_bit + second_highest_bit :
      (highest_bit == 0) ? 1 : 0;
}

enum Mode {
  kLiteral,
  kPaletteIx,
  kCopy,
  kNone,
};

typedef struct {
  // mode as uint8_t to make the memory layout to be exactly 8 bytes.
  uint8_t mode;
  uint16_t len;
  uint32_t argb_or_offset;
} PixOrCopy;

static inline PixOrCopy PixOrCopy_CreateCopy(uint32_t offset_arg,
                                             uint16_t len_arg) {
  PixOrCopy retval;
  retval.mode = kCopy;
  retval.argb_or_offset = offset_arg;
  retval.len = len_arg;
  return retval;
}

static inline PixOrCopy PixOrCopy_CreatePaletteIx(int ix) {
  PixOrCopy retval;
  assert(ix >= 0);
  assert(ix < (1 << kPaletteCodeBitsMax));
  retval.mode = kPaletteIx;
  retval.argb_or_offset = ix;
  retval.len = 1;
  return retval;
}

static inline PixOrCopy PixOrCopy_CreateLiteral(uint32_t argb_arg) {
  PixOrCopy retval;
  retval.mode = kLiteral;
  retval.argb_or_offset = argb_arg;
  retval.len = 1;
  return retval;
}

static inline int PixOrCopy_IsLiteral(const PixOrCopy *p) {
  return p->mode == kLiteral;
}

static inline int PixOrCopy_IsPaletteIx(const PixOrCopy *p) {
  return p->mode == kPaletteIx;
}

static inline int PixOrCopy_IsCopy(const PixOrCopy *p) {
  return p->mode == kCopy;
}

static inline uint32_t PixOrCopy_Literal(const PixOrCopy *p, int component) {
  assert(p->mode == kLiteral);
  return (p->argb_or_offset >> (component * 8)) & 0xff;
}

static inline uint32_t PixOrCopy_Length(const PixOrCopy *p) {
  return p->len;
}

static inline uint32_t PixOrCopy_Argb(const PixOrCopy *p) {
  assert(p->mode == kLiteral);
  return p->argb_or_offset;
}

static inline uint32_t PixOrCopy_PaletteIx(const PixOrCopy *p) {
  assert(p->mode == kPaletteIx);
  assert(p->argb_or_offset < (1 << kPaletteCodeBitsMax));
  return p->argb_or_offset;
}

static inline uint32_t PixOrCopy_Distance(const PixOrCopy *p) {
  assert(p->mode == kCopy);
  return p->argb_or_offset;
}

static inline void PixOrCopy_LengthCodeAndBits(
    const PixOrCopy *p, int *code, int *n_bits, int *bits) {
  assert(p->len >= 1 && p->len <= kMaxLength);
  PrefixEncode(p->len, code, n_bits, bits);
}


// Ridiculously simple backward references for images where it is unlikely
// that there are large backward references (photos).
void BackwardReferencesRle(
    int xsize,
    int ysize,
    const uint32_t *argb,
    PixOrCopy *stream,
    int *stream_size);

// This is a simple fast function for obtaining backward references
// based on simple heuristics.
void BackwardReferencesHashChain(
    int xsize,
    int ysize,
    int use_palette,
    const uint32_t *argb,
    int palette_bits,
    PixOrCopy *stream,
    int *stream_size);

// This method looks for a shortest path through the backward reference
// network based on a cost model generated by a first round of compression.
void BackwardReferencesTraceBackwards(
    int xsize,
    int ysize,
    int recursive_cost_model,
    int use_palette,
    const uint32_t *argb,
    int palette_bits,
    PixOrCopy *stream,
    int *stream_size);

// Convert backward references that are of linear distance along
// the image scan lines to have a 2d locality indexing where
// smaller values are used for backward references that are close by.
void BackwardReferences2DLocality(int xsize, int data_size,
                                  PixOrCopy *data);

// Internals of locality transform exposed for testing use.
int DistanceToPlaneCode(int xsize, int distance);

// Returns true if the given backward references actually produce
// the image given in tuple (argb, xsize, ysize).
int VerifyBackwardReferences(const uint32_t* argb,
                             int xsize, int ysize,
                             int palette_bits,
                             const PixOrCopy *lit,
                             int lit_size);

// Produce an estimate for a good emerging palette size for the image.
int CalculateEstimateForPaletteSize(const uint32_t *argb, int xsize, int ysize);

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif  // WEBP_BACKWARD_REFERENCES_H_
