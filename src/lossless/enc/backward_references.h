// Copyright 2011 Google Inc.
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
#include <vector>

#include "backward_distance.h"
#include "../common/integral_types.h"

// Backward reference distance codes, for all 32-bit values.
static const int kDistanceCodes = 40;

// Compression constants
static const int kLengthCodes = 24;
static const int kCodeLengthCodes = 19;
static const int kRowHasherXSubsampling = 7;
static const int kPaletteCodeBitsMax = 11;
static const int kLiteralOrCopyCodesMax =
    256 + kLengthCodes + (1 << kPaletteCodeBitsMax);
static const int kMaxLength = 4096;

struct LiteralOrCopy {
  LiteralOrCopy() {
    mode = kNone;
  }
  LiteralOrCopy(const LiteralOrCopy &r) {
    mode = r.mode;
    argb_or_offset = r.argb_or_offset;
    len = r.len;
  }

  const LiteralOrCopy &operator =(const LiteralOrCopy &r) {
    mode = r.mode;
    argb_or_offset = r.argb_or_offset;
    len = r.len;
    return *this;
  }

  static LiteralOrCopy CreateCopy(uint32 offset_arg, uint16 len_arg) {
    LiteralOrCopy retval;
    retval.mode = kCopy;
    retval.argb_or_offset = offset_arg;
    retval.len = len_arg;
    return retval;
  }
  static LiteralOrCopy CreatePaletteIx(int ix) {
    assert(ix >= 0);
    assert(ix < (1 << kPaletteCodeBitsMax));
    LiteralOrCopy retval;
    retval.mode = kPaletteIx;
    retval.argb_or_offset = ix;
    retval.len = 1;
    return retval;
  }
  static LiteralOrCopy CreateLiteral(uint32 argb_arg) {
    LiteralOrCopy retval;
    retval.mode = kLiteral;
    retval.argb_or_offset = argb_arg;
    retval.len = 1;
    return retval;
  }
  enum Mode {
    kLiteral,
    kPaletteIx,
    kCopy,
    kNone,
  };
  bool IsLiteral() const {
    return mode == kLiteral;
  }
  bool IsPaletteIx() const {
    return mode == kPaletteIx;
  }
  bool IsCopy() const {
    return mode == kCopy;
  }
  uint32 Literal(int component) const {
    assert(mode == kLiteral);
    return (argb_or_offset >> (component * 8)) & 0xff;
  }
  uint32 Length() const {
    return len;
  }
  uint32 Argb() const {
    assert(mode == kLiteral);
    return argb_or_offset;
  }
  uint32 PaletteIx() const {
    assert(mode == kPaletteIx);
    assert(argb_or_offset >= 0);
    assert(argb_or_offset < (1 << kPaletteCodeBitsMax));
    return argb_or_offset;
  }
  uint32 Distance() const {
    assert(mode == kCopy);
    return argb_or_offset;
  }
  inline void LengthCodeAndBits(int *code, int *n_bits, int *bits) const {
    assert(len >= 1 && len <= kMaxLength);
    // Unlike flate, distance and length are encoded the same way.
    BackwardLength::Encode(len, code, n_bits, bits);
  }

  // mode as uint 8, and not as type Mode, to make the memory layout
  // of this class as 8 bytes.
  uint8 mode;
  uint16 len;
  uint32 argb_or_offset;
};

// Ridiculously simple backward references for images where it is unlikely
// that there are large backward references (photos).
void BackwardReferencesRle(
    int xsize,
    int ysize,
    const uint32 *argb,
    std::vector<LiteralOrCopy> *stream);

// This is a simple fast function for obtaining backward references
// based on simple heuristics.
void BackwardReferencesHashChain(
    int xsize,
    int ysize,
    bool use_palette,
    const uint32 *argb,
    int palette_bits,
    std::vector<LiteralOrCopy> *stream);

// This method looks for a shortest path through the backward reference
// network based on a cost model generated by a first round of compression.
void BackwardReferencesTraceBackwards(
    int xsize,
    int ysize,
    int recursive_cost_model,
    bool use_palette,
    const uint32 *argb,
    int palette_bits,
    std::vector<LiteralOrCopy> *stream);


// Convert backward references that are of linear distance along
// the image scan lines to have a 2d locality indexing where
// smaller values are used for backward references that are close by.
void BackwardReferences2DLocality(int xsize,
                                  int ysize,
                                  int data_size,
                                  LiteralOrCopy *data);

// Internals of locality transform exposed for testing use.
int DistanceToPlaneCode(int xsize, int ysize, int distance);

// Returns true if the given backward references actually produce
// the image given in tuple (argb, xsize, ysize).
bool VerifyBackwardReferences(const uint32* argb,
                              int xsize, int ysize,
                              int palette_bits,
                              const std::vector<LiteralOrCopy>& v);

// Produce an estimate for a good emerging palette size for the image.
int CalculateEstimateForPaletteSize(const uint32 *argb, int xsize, int ysize);

#endif  // WEBP_BACKWARD_REFERENCES_H_
