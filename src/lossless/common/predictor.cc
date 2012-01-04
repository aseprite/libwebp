// Copyright 2011 Google Inc.
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

static inline uint32 Average2(uint32 a0, uint32 a1) {
  return ( ((((a0) ^ (a1)) & 0xfefefefeL) >> 1) + ((a0) & (a1)) );
}

static inline uint32 Average3(uint32 a0, uint32 a1, uint32 a2) {
  return Average2(Average2(a0, a2), a1);
}

static inline uint32 Average4(uint32 a0, uint32 a1, uint32 a2, uint32 a3) {
  return Average2(Average2(a0, a1), Average2(a2, a3));
}

uint32 Add(uint32 a, uint32 b) {
  // This computes the sum of each component with mod 256.
  uint32 alpha_and_green = (a & 0xff00ff00) + (b & 0xff00ff00);
  uint32 red_and_blue = (a & 0x00ff00ff) + (b & 0x00ff00ff);
  return (alpha_and_green & 0xff00ff00) | (red_and_blue & 0x00ff00ff);
}

uint32 Subtract(uint32 a, uint32 b) {
  // This subtracts each component with mod 256.
  uint32 alpha_and_green = 0x00ff00ff + (a & 0xff00ff00) - (b & 0xff00ff00);
  uint32 red_and_blue = 0xff00ff00 + (a & 0x00ff00ff) - (b & 0x00ff00ff);
  return (alpha_and_green & 0xff00ff00) | (red_and_blue & 0x00ff00ff);
}

uint32 Paeth32(uint32 a, uint32 b, uint32 c) {
  const int p0 = int(a >> 24) + int(b >> 24) - int(c >> 24);
  const int p1 = int((a >> 16) & 0xff) + int((b >> 16) & 0xff) -
      int((c >> 16) & 0xff);
  const int p2 = int((a >> 8) & 0xff) + int((b >> 8) & 0xff) -
      int((c >> 8) & 0xff);
  const int p3 = int(a & 0xff) + int(b & 0xff) - int(c & 0xff);

  int pa = abs(p0 - (a >> 24)) +
      abs(p1 - ((a >> 16) & 0xff)) +
      abs(p2 - ((a >> 8) & 0xff)) +
      abs(p3 - (a & 0xff));
  int pb = abs(p0 - (b >> 24)) +
      abs(p1 - ((b >> 16) & 0xff)) +
      abs(p2 - ((b >> 8) & 0xff)) +
      abs(p3 - (b & 0xff));
  int pc = abs(p0 - (c >> 24)) +
      abs(p1 - ((c >> 16) & 0xff)) +
      abs(p2 - ((c >> 8) & 0xff)) +
      abs(p3 - (c & 0xff));

  if (pa <= pb && pa <= pc) {
    return a;
  } else if (pb <= pc) {
    return b;
  }
  return c;
}

inline int Clamp(int a) {
  if (a < 0) return 0;
  if (a > 255) return 255;
  return a;
}

inline int AddSubtractComponentFull(int a, int b, int c) {
  return Clamp(a + b - c);
}

uint32 ClampedAddSubtractFull(uint32 c0, uint32 c1, uint32 c2) {
  int a = AddSubtractComponentFull(c0 >> 24, c1 >> 24, c2 >> 24);
  int r = AddSubtractComponentFull((c0 >> 16) & 0xff,
                                   (c1 >> 16) & 0xff,
                                   (c2 >> 16) & 0xff);
  int g = AddSubtractComponentFull((c0 >> 8) & 0xff,
                                   (c1 >> 8) & 0xff,
                                   (c2 >> 8) & 0xff);
  int b = AddSubtractComponentFull((c0 >> 0) & 0xff,
                                   (c1 >> 0) & 0xff,
                                   (c2 >> 0) & 0xff);
  return (a << 24) | (r << 16) | (g << 8) | b;
}

inline int AddSubtractComponentHalf(int a, int b) {
  return Clamp(a + (a - b) / 2);
}

uint32 ClampedAddSubtractHalf(uint32 c0, uint32 c1) {
  int a = AddSubtractComponentHalf(c0 >> 24, c1 >> 24);
  int r = AddSubtractComponentHalf((c0 >> 16) & 0xff, (c1 >> 16) & 0xff);
  int g = AddSubtractComponentHalf((c0 >> 8) & 0xff, (c1 >> 8) & 0xff);
  int b = AddSubtractComponentHalf((c0 >> 0) & 0xff, (c1 >> 0) & 0xff);
  return (a << 24) | (r << 16) | (g << 8) | b;
}

uint32 PredictValue(int mode, int xsize, const uint32 *argb) {
  switch (mode) {
    case 0: return 0xff000000;
    case 1: return argb[-1];
    case 2: return argb[-xsize];
    case 3: return argb[-xsize + 1];
    case 4: return argb[-xsize - 1];

    case 5: return Average3(argb[-1],
                            argb[-xsize],
                            argb[-xsize + 1]);
    case 6: return Average2(argb[-1],
                            argb[-xsize - 1]);
    case 7: return Average2(argb[-1],
                            argb[-xsize]);
    case 8: return Average2(argb[-xsize - 1],
                            argb[-xsize]);
    case 9: return Average2(argb[-xsize],
                            argb[-xsize + 1]);
    case 10: return Average4(argb[-1],
                             argb[-xsize - 1],
                             argb[-xsize],
                             argb[-xsize + 1]);
    case 11: return Paeth32(argb[-xsize],
                            argb[-1],
                            argb[-xsize - 1]);
    case 12:
      return ClampedAddSubtractFull(argb[-1],
                                    argb[-xsize],
                                    argb[-xsize - 1]);
    case 13:
      {
        uint32 ave = Average2(argb[-1], argb[-xsize]);
        return ClampedAddSubtractHalf(ave, argb[-xsize - 1]);
      }
  }
  return 0;
}

void CopyTileWithPrediction(int xsize, int ysize,
                            const uint32 *from_argb,
                            int tile_x, int tile_y, int bits, int mode,
                            uint32 *to_argb) {
  int ymax = 1 << bits;
  if (ymax > ysize - (tile_y << bits)) {
    ymax = ysize - (tile_y << bits);
  }
  int xmax = 1 << bits;
  if (xmax > xsize - (tile_x << bits)) {
    xmax = xsize - (tile_x << bits);
  }
  for (int y = 0; y < ymax; ++y) {
    int all_y = (tile_y << bits) + y;
    for (int x = 0; x < xmax; ++x) {
      int all_x = (tile_x << bits) + x;
      const int ix = all_y * xsize + all_x;
      uint32 predict;
      if (all_y == 0) {
        if (all_x == 0) {
          predict = 0xff000000;
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
                               const uint32* image,
                               const uint32 *from_argb,
                               uint32* to_argb) {
  const int tile_xsize = (xsize + (1 << bits) - 1) >> bits;
  {
    // The first row is special.
    to_argb[0] = from_argb[0] + 0xff000000;
    for (int image_x = 1; image_x < xsize; ++image_x) {
      uint32 predict = to_argb[image_x - 1];
      to_argb[image_x] = Add(from_argb[image_x], predict);
    }
  }
  for (int image_y = 1; image_y < ysize; ++image_y) {
    int tile_y = image_y >> bits;
    int ix = image_y * xsize;
    for (int tile_x = 0; tile_x < tile_xsize; ++tile_x) {
      int tile_ix = tile_y * tile_xsize + tile_x;
      int mode = ((image[tile_ix] >> 8) & 0xff);
      int image_x = (tile_x << bits);
      int xend = image_x + (1 << bits);
      if (xend > xsize) {
        xend = xsize;
      }
      if (image_x == 0) {
        uint32 predict = to_argb[ix + image_x - xsize];
        to_argb[ix + image_x] = Add(from_argb[ix + image_x], predict);
        ++image_x;
      }
      if (mode == 0) {
        for (; image_x < xend; ++image_x) {
          to_argb[ix + image_x] = from_argb[ix + image_x] + 0xff000000;
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
          uint32 predict = PredictValue(mode, xsize, to_argb + ix + image_x);
          to_argb[ix + image_x] = Add(from_argb[ix + image_x], predict);
        }
      }
    }
  }
}

void CopyTileWithColorTransform(int xsize, int ysize,
                                const uint32 *from_argb,
                                int tile_x, int tile_y, int bits,
                                ColorSpaceTransformElement color_transform,
                                bool inverse,
                                uint32 *to_argb) {
  tile_x <<= bits;
  tile_y <<= bits;
  int xscan = 1 << bits;
  if (xscan > xsize - tile_x) {
    xscan = xsize - tile_x;
  }
  int yscan = 1 << bits;
  if (yscan > ysize - tile_y) {
    yscan = ysize - tile_y;
  }
  yscan += tile_y;
  if (inverse) {
    for (int y = tile_y; y < yscan; ++y) {
      int ix = y * xsize + tile_x;
      int end_ix = ix + xscan;
      for (;ix < end_ix; ++ix) {
        to_argb[ix] = color_transform.InverseTransform(from_argb[ix]);
      }
    }
  } else {
    for (int y = tile_y; y < yscan; ++y) {
      int ix = y * xsize + tile_x;
      int end_ix = ix + xscan;
      for (;ix < end_ix; ++ix) {
        to_argb[ix] = color_transform.Transform(from_argb[ix]);
      }
    }
  }
}

void ColorSpaceInverseTransform(int xsize, int ysize, int bits,
                                const uint32* image,
                                const uint32 *from_argb,
                                uint32* to_argb) {
  int tile_xsize = (xsize + (1 << bits) - 1) >> bits;
  int tile_ysize = (ysize + (1 << bits) - 1) >> bits;
  int tile_ix = 0;
  for (int tile_y = 0; tile_y < tile_ysize; ++tile_y) {
    for (int tile_x = 0; tile_x < tile_xsize; ++tile_x, ++tile_ix) {
      ColorSpaceTransformElement color_transform;
      color_transform.InitFromCode(image[tile_ix]);
      CopyTileWithColorTransform(xsize, ysize, from_argb,
                                 tile_x, tile_y, bits, color_transform,
                                 true,  // inverse transform
                                 to_argb);
    }
  }
}

void AddGreenToBlueAndRed(int n, uint32 *argb_array) {
  for (int i = 0; i < n; ++i) {
    uint32 argb = argb_array[i];
    uint32 green = ((argb >> 8) & 0xff);
    green += green << 16;
    argb_array[i] = (argb & 0xff00ff00) +
        (((argb & 0x00ff00ff) + green) & 0x00ff00ff);
  }
}
