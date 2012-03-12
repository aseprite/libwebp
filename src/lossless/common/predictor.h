// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)
//
// Simple predictor for webp lossless.

#ifndef WEBP_LOSSLESS_COMMON_PREDICTOR_H_
#define WEBP_LOSSLESS_COMMON_PREDICTOR_H_

#include <stdint.h>

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

static inline uint32_t ColorTransformDelta(int8_t transform, int8_t color) {
  return (uint32_t)(transform * color) >> 5;
}

uint32_t PredictValue(int mode, int xsize, const uint32_t* argb);

uint32_t Add(uint32_t a, uint32_t b);
uint32_t Subtract(uint32_t a, uint32_t b);

// from_argb and to_argb represent images of size (xsize, ysize).
// This function copies only the tile at tile_x, tile_y where the tile size is
// (1 << bits, 1 << bits).
void CopyTileWithPrediction(int xsize, int ysize,
                            const uint32_t* from_argb,
                            int tile_x, int tile_y, int bits, int mode,
                            uint32_t* to_argb);

void PredictorInverseTransform(int xsize, int ysize, int bits,
                               const uint32_t* image,
                               const uint32_t* original_argb,
                               uint32_t* to_argb);

// The transform consists of element contributing the green value into
// red and blue, and the red value into blue, using a linear weighting and
// summing.
typedef struct {
  uint8_t green_to_red_;
  uint8_t green_to_blue_;
  uint8_t red_to_blue_;
} ColorSpaceTransformElement;

static inline void ColorSpaceTransformElementClear(
    ColorSpaceTransformElement* p) {
  p->green_to_red_ = 0;
  p->green_to_blue_ = 0;
  p->red_to_blue_ = 0;
}

static inline uint32_t ColorSpaceTransformElementTransform(
    const ColorSpaceTransformElement* p, uint32_t argb) {
  // green, red, new_red and new_blue temporarily contain dirty upper bits,
  // which is not a problem since they are masked off in the end.
  const uint32_t green = argb >> 8;
  const uint32_t red = argb >> 16;
  uint32_t new_red = red;
  uint32_t new_blue = argb;
  new_red -= ColorTransformDelta(p->green_to_red_, green);
  new_blue -= ColorTransformDelta(p->green_to_blue_, green);
  new_blue -= ColorTransformDelta(p->red_to_blue_, red);
  new_red &= 0xff;
  new_blue &= 0xff;
  return (argb & 0xff00ff00u) | (new_red << 16) | (new_blue);
}

static inline uint32_t ColorSpaceTransformElementInverseTransform(
    const ColorSpaceTransformElement* p, uint32_t argb) {
  // green, red, new_red and new_blue temporarily contain dirty upper bits,
  // which is not a problem since they are masked off in the end.
  const uint32_t green = argb >> 8;
  uint32_t red = argb >> 16;
  uint32_t blue = argb;
  red += ColorTransformDelta(p->green_to_red_, green);
  blue += ColorTransformDelta(p->green_to_blue_, green);
  red &= 0xff;
  blue += ColorTransformDelta(p->red_to_blue_, red);
  blue &= 0xff;
  return (argb & 0xff00ff00u) | (red << 16) | (blue);
}

static inline uint32_t ColorSpaceTransformElementCode(
    const ColorSpaceTransformElement* p) {
  return
      0xff000000u |
      ((uint32_t)(p->red_to_blue_) << 16) |
      ((uint32_t)(p->green_to_blue_) << 8) |
      p->green_to_red_;
}

static inline void ColorSpaceTransformElementInitFromCode(
    ColorSpaceTransformElement* p,
    uint32_t code) {
  p->green_to_red_ = code & 0xff;
  p->green_to_blue_ = (code >> 8) & 0xff;
  p->red_to_blue_ = (code >> 16) & 0xff;
}

void CopyTileWithColorTransform(int xsize, int ysize,
                                const uint32_t* from_argb,
                                int tile_x, int tile_y, int bits,
                                ColorSpaceTransformElement color_transform,
                                int inverse,
                                uint32_t* to_argb);

void ColorSpaceInverseTransform(int xsize, int ysize, int bits,
                                const uint32_t* image,
                                const uint32_t* original_argb,
                                uint32_t* to_argb);

void AddGreenToBlueAndRed(int n, uint32_t* argb_array);

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif  // WEBP_LOSSLESS_COMMON_PREDICTOR_H_
