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

#include "integral_types.h"

static inline uint32 ColorTransformDelta(signed char transform,
                                         signed char color) {
  return uint32(int(transform) * color) >> 5;
}

uint32 PredictValue(int mode, int xsize, const uint32 *argb);

uint32 Add(uint32 a, uint32 b);
uint32 Subtract(uint32 a, uint32 b);

// from_argb and to_argb represent images of size (xsize, ysize).
// This function copies only the tile at tile_x, tile_y where the tile size is
// (1 << bits, 1 << bits).
void CopyTileWithPrediction(int xsize, int ysize,
                            const uint32 *from_argb,
                            int tile_x, int tile_y, int bits, int mode,
                            uint32 *to_argb);

void PredictorInverseTransform(int xsize, int ysize, int bits,
                               const uint32* image,
                               const uint32 *original_argb,
                               uint32* to_argb);

// The transform consists of element contributing the green value into
// red and blue, and the red value into blue, using a linear weighting and
// summing.
struct ColorSpaceTransformElement {
  uint8 green_to_red_;
  uint8 green_to_blue_;
  uint8 red_to_blue_;

  void Clear() {
    green_to_red_ = 0;
    green_to_blue_ = 0;
    red_to_blue_ = 0;
  }

  uint32 Transform(uint32 argb) const {
    // green, red, new_red and new_blue temporarily contain dirty upper bits,
    // which is not a problem since they are masked off in the end.
    const uint32 green = argb >> 8;
    const uint32 red = argb >> 16;
    uint32 new_red = red;
    uint32 new_blue = argb;

    new_red -= ColorTransformDelta(green_to_red_, green);
    new_blue -= ColorTransformDelta(green_to_blue_, green);
    new_blue -= ColorTransformDelta(red_to_blue_, red);

    new_red &= 0xff;
    new_blue &= 0xff;
    return (argb & 0xff00ff00) | (new_red << 16) | (new_blue);
  }

  uint32 InverseTransform(uint32 argb) const {
    // green, red, new_red and new_blue temporarily contain dirty upper bits,
    // which is not a problem since they are masked off in the end.
    const uint32 green = argb >> 8;
    const uint32 red = argb >> 16;
    uint32 new_red = red;
    uint32 new_blue = argb;

    new_red += ColorTransformDelta(green_to_red_, green);
    new_blue += ColorTransformDelta(green_to_blue_, green);
    new_red &= 0xff;
    new_blue += ColorTransformDelta(red_to_blue_, new_red);
    new_blue &= 0xff;
    return (argb & 0xff00ff00) | (new_red << 16) | (new_blue);
  }

  uint32 Code() const {
    return
        0xff000000 |
        (static_cast<uint32>(red_to_blue_) << 16) |
        (static_cast<uint32>(green_to_blue_) << 8) |
        green_to_red_;
  }

  void InitFromCode(uint32 code) {
    green_to_red_ = code & 0xff;
    green_to_blue_ = (code >> 8) & 0xff;
    red_to_blue_ = (code >> 16) & 0xff;
  }
};

void CopyTileWithColorTransform(int xsize, int ysize,
                                const uint32 *from_argb,
                                int tile_x, int tile_y, int bits,
                                ColorSpaceTransformElement color_transform,
                                bool inverse,
                                uint32 *to_argb);

void ColorSpaceInverseTransform(int xsize, int ysize, int bits,
                                const uint32* image,
                                const uint32 *original_argb,
                                uint32* to_argb);

void AddGreenToBlueAndRed(int n, uint32 *argb_array);

#endif  // WEBP_LOSSLESS_COMMON_PREDICTOR_H_
