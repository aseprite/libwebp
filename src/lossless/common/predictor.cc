// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)

#include "predictor.h"

#include <assert.h>
#include <stdlib.h>

#include "integral_types.h"

#define ARGB_BLACK 0xff000000

static inline uint32_t Average2(uint32_t a0, uint32_t a1) {
  return (((a0 ^ a1) & 0xfefefefeL) >> 1) + (a0 & a1);
}

static inline uint32_t Average3(uint32_t a0, uint32_t a1, uint32_t a2) {
  return Average2(Average2(a0, a2), a1);
}

static inline uint32_t Average4(uint32_t a0, uint32_t a1,
                                uint32_t a2, uint32_t a3) {
  return Average2(Average2(a0, a1), Average2(a2, a3));
}

uint32_t Add(uint32_t a, uint32_t b) {
  // This computes the sum of each component with mod 256.
  const uint32_t alpha_and_green = (a & 0xff00ff00) + (b & 0xff00ff00);
  const uint32_t red_and_blue = (a & 0x00ff00ff) + (b & 0x00ff00ff);
  return (alpha_and_green & 0xff00ff00) | (red_and_blue & 0x00ff00ff);
}

uint32_t Subtract(uint32_t a, uint32_t b) {
  // This subtracts each component with mod 256.
  const uint32_t alpha_and_green =
      0x00ff00ff + (a & 0xff00ff00) - (b & 0xff00ff00);
  const uint32_t red_and_blue =
      0xff00ff00 + (a & 0x00ff00ff) - (b & 0x00ff00ff);
  return (alpha_and_green & 0xff00ff00) | (red_and_blue & 0x00ff00ff);
}

static uint32_t Select(uint32_t a, uint32_t b, uint32_t c) {
  const int p0 = int(a >> 24) + int(b >> 24) - int(c >> 24);
  const int p1 = int((a >> 16) & 0xff) + int((b >> 16) & 0xff) -
      int((c >> 16) & 0xff);
  const int p2 = int((a >> 8) & 0xff) + int((b >> 8) & 0xff) -
      int((c >> 8) & 0xff);
  const int p3 = int(a & 0xff) + int(b & 0xff) - int(c & 0xff);
  const int pa = abs(p0 - (a >> 24)) +
      abs(p1 - ((a >> 16) & 0xff)) +
      abs(p2 - ((a >> 8) & 0xff)) +
      abs(p3 - (a & 0xff));
  const int pb = abs(p0 - (b >> 24)) +
      abs(p1 - ((b >> 16) & 0xff)) +
      abs(p2 - ((b >> 8) & 0xff)) +
      abs(p3 - (b & 0xff));
  if (pa <= pb) {
    return a;
  } else {
    return b;
  }
}

inline uint32_t Clamp(uint32_t a) {
  if (a <= 255) {
    return a;
  }
  // return 0, when a is a negative integer.
  // return 255, when a is positive.
  return ~a >> 24;
}

inline int AddSubtractComponentFull(int a, int b, int c) {
  return Clamp(a + b - c);
}

static uint32_t ClampedAddSubtractFull(uint32_t c0, uint32_t c1, uint32_t c2) {
  const int a = AddSubtractComponentFull(c0 >> 24, c1 >> 24, c2 >> 24);
  const int r = AddSubtractComponentFull((c0 >> 16) & 0xff,
                                         (c1 >> 16) & 0xff,
                                         (c2 >> 16) & 0xff);
  const int g = AddSubtractComponentFull((c0 >> 8) & 0xff,
                                         (c1 >> 8) & 0xff,
                                         (c2 >> 8) & 0xff);
  const int b = AddSubtractComponentFull(c0 & 0xff,
                                         c1 & 0xff,
                                         c2 & 0xff);
  return (a << 24) | (r << 16) | (g << 8) | b;
}

inline int AddSubtractComponentHalf(int a, int b) {
  return Clamp(a + (a - b) / 2);
}

static uint32_t ClampedAddSubtractHalf(uint32_t c0, uint32_t c1, uint32_t c2) {
  const uint32_t ave = Average2(c0, c1);
  const int a = AddSubtractComponentHalf(ave >> 24, c2 >> 24);
  const int r = AddSubtractComponentHalf((ave >> 16) & 0xff, (c2 >> 16) & 0xff);
  const int g = AddSubtractComponentHalf((ave >> 8) & 0xff, (c2 >> 8) & 0xff);
  const int b = AddSubtractComponentHalf((ave >> 0) & 0xff, (c2 >> 0) & 0xff);
  return (a << 24) | (r << 16) | (g << 8) | b;
}

uint32_t PredictValue(int mode, int xsize, const uint32_t *argb) {
  switch (mode) {
    case 0: return ARGB_BLACK;
    case 1: return argb[-1];
    case 2: return argb[-xsize];
    case 3: return argb[-xsize + 1];
    case 4: return argb[-xsize - 1];
    case 5: return Average3(argb[-1], argb[-xsize], argb[-xsize + 1]);
    case 6: return Average2(argb[-1], argb[-xsize - 1]);
    case 7: return Average2(argb[-1], argb[-xsize]);
    case 8: return Average2(argb[-xsize - 1], argb[-xsize]);
    case 9: return Average2(argb[-xsize], argb[-xsize + 1]);
    case 10: return Average4(argb[-1], argb[-xsize - 1],
                             argb[-xsize], argb[-xsize + 1]);
    case 11: return Select(argb[-xsize], argb[-1], argb[-xsize - 1]);
    case 12:
      return ClampedAddSubtractFull(argb[-1], argb[-xsize], argb[-xsize - 1]);
    case 13:
      return ClampedAddSubtractHalf(argb[-1], argb[-xsize], argb[-xsize - 1]);
  }
  return 0;
}

void CopyTileWithPrediction(int xsize, int ysize,
                            const uint32_t *from_argb,
                            int tile_x, int tile_y, int bits, int mode,
                            uint32_t *to_argb) {
  int ymax = 1 << bits;
  int xmax = 1 << bits;
  int y;
  if (ymax > ysize - (tile_y << bits)) {
    ymax = ysize - (tile_y << bits);
  }
  if (xmax > xsize - (tile_x << bits)) {
    xmax = xsize - (tile_x << bits);
  }
  for (y = 0; y < ymax; ++y) {
    const int all_y = (tile_y << bits) + y;
    int x;
    for (x = 0; x < xmax; ++x) {
      const int all_x = (tile_x << bits) + x;
      const int ix = all_y * xsize + all_x;
      uint32_t predict;
      if (all_y == 0) {
        if (all_x == 0) {
          predict = ARGB_BLACK;
        } else {
          predict = from_argb[ix - 1];
        }
      } else if (all_x == 0) {
        predict = from_argb[ix - xsize];
      } else {
        predict = PredictValue(mode, xsize, from_argb + ix);
      }
      to_argb[ix] = Subtract(from_argb[ix], predict);
    }
  }
}

void PredictorInverseTransform(int xsize, int ysize, int bits,
                               const uint32_t* image,
                               const uint32_t *from_argb,
                               uint32_t* to_argb) {
  int image_y;
  const int tile_xsize = (xsize + (1 << bits) - 1) >> bits;
  {
    // The first row is special.
    int image_x;
    to_argb[0] = from_argb[0] + ARGB_BLACK;
    for (image_x = 1; image_x < xsize; ++image_x) {
      const uint32_t predict = to_argb[image_x - 1];
      to_argb[image_x] = Add(from_argb[image_x], predict);
    }
  }
  for (image_y = 1; image_y < ysize; ++image_y) {
    const int tile_y = image_y >> bits;
    int ix = image_y * xsize;
    int tile_x;
    for (tile_x = 0; tile_x < tile_xsize; ++tile_x) {
      const int tile_ix = tile_y * tile_xsize + tile_x;
      const int mode = ((image[tile_ix] >> 8) & 0xff);
      int image_x = (tile_x << bits);
      int xend = image_x + (1 << bits);
      if (xend > xsize) {
        xend = xsize;
      }
      if (image_x == 0) {
        const uint32_t predict = to_argb[ix + image_x - xsize];
        to_argb[ix + image_x] = Add(from_argb[ix + image_x], predict);
        ++image_x;
      }
      if (mode == 0) {
        for (; image_x < xend; ++image_x) {
          to_argb[ix + image_x] = from_argb[ix + image_x] + ARGB_BLACK;
        }
      } else if (mode == 1) {
        for (; image_x < xend; ++image_x) {
          to_argb[ix + image_x] =
              Add(from_argb[ix + image_x], to_argb[ix + image_x - 1]);
        }
      } else if (mode == 2) {
        for (; image_x < xend; ++image_x) {
          to_argb[ix + image_x] =
              Add(from_argb[ix + image_x], to_argb[ix + image_x - xsize]);
        }
      } else {
        for (; image_x < xend; ++image_x) {
          const uint32_t predict =
              PredictValue(mode, xsize, to_argb + ix + image_x);
          to_argb[ix + image_x] = Add(from_argb[ix + image_x], predict);
        }
      }
    }
  }
}

void CopyTileWithColorTransform(int xsize, int ysize,
                                const uint32_t *from_argb,
                                int tile_x, int tile_y, int bits,
                                ColorSpaceTransformElement color_transform,
                                int inverse,
                                uint32_t *to_argb) {
  int xscan = 1 << bits;
  int yscan = 1 << bits;
  tile_x <<= bits;
  tile_y <<= bits;
  if (xscan > xsize - tile_x) {
    xscan = xsize - tile_x;
  }
  if (yscan > ysize - tile_y) {
    yscan = ysize - tile_y;
  }
  yscan += tile_y;
  if (inverse) {
    int y;
    for (y = tile_y; y < yscan; ++y) {
      int ix = y * xsize + tile_x;
      const int end_ix = ix + xscan;
      for (;ix < end_ix; ++ix) {
        to_argb[ix] = ColorSpaceTransformElementInverseTransform(
            &color_transform, from_argb[ix]);
      }
    }
  } else {
    int y;
    for (y = tile_y; y < yscan; ++y) {
      int ix = y * xsize + tile_x;
      const int end_ix = ix + xscan;
      for (;ix < end_ix; ++ix) {
        to_argb[ix] = ColorSpaceTransformElementTransform(
            &color_transform, from_argb[ix]);
      }
    }
  }
}

void ColorSpaceInverseTransform(int xsize, int ysize, int bits,
                                const uint32_t* image,
                                const uint32_t *from_argb,
                                uint32_t* to_argb) {
  const int tile_xsize = (xsize + (1 << bits) - 1) >> bits;
  const int tile_ysize = (ysize + (1 << bits) - 1) >> bits;
  int tile_ix = 0;
  int tile_y;
  int tile_x;
  for (tile_y = 0; tile_y < tile_ysize; ++tile_y) {
    for (tile_x = 0; tile_x < tile_xsize; ++tile_x, ++tile_ix) {
      ColorSpaceTransformElement color_transform;
      ColorSpaceTransformElementInitFromCode(&color_transform, image[tile_ix]);
      CopyTileWithColorTransform(xsize, ysize, from_argb,
                                 tile_x, tile_y, bits, color_transform,
                                 true,  // inverse transform
                                 to_argb);
    }
  }
}

void AddGreenToBlueAndRed(int n, uint32_t *argb_array) {
  int i;
  for (i = 0; i < n; ++i) {
    const uint32_t argb = argb_array[i];
    uint32_t green = ((argb >> 8) & 0xff);
    green += green << 16;
    argb_array[i] = (argb & 0xff00ff00) +
        (((argb & 0x00ff00ff) + green) & 0x00ff00ff);
  }
}
