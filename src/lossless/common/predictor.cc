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

inline uint32 Average2(uint32 a0, uint32 a1) {
  uint32 alpha_and_green =
      ((((a0 >> 8) & 0x00ff00ff) +
        ((a1 >> 8) & 0x00ff00ff) + 0x00010001) >> 1) & 0x00ff00ff;
  uint32 red_and_blue =
      (((a0 & 0x00ff00ff) +
        (a1 & 0x00ff00ff) + 0x00010001) >> 1) & 0x00ff00ff;
  return (alpha_and_green << 8) | red_and_blue;
}

uint32 Average(uint32 a0, uint32 a1, uint32 a2) {
  int alpha = ((a0 >> 24) + (a1 >> 24) + (a2 >> 24) + 2) / 3;
  int red = (((a0 >> 16) & 0xff) +
             ((a1 >> 16) & 0xff) +
             ((a2 >> 16) & 0xff) + 2) / 3;
  int green = (((a0 >> 8) & 0xff) +
               ((a1 >> 8) & 0xff) +
               ((a2 >> 8) & 0xff) + 2) / 3;
  int blue = (((a0 >> 0) & 0xff) +
              ((a1 >> 0) & 0xff) +
              ((a2 >> 0) & 0xff) + 2) / 3;
  return (alpha << 24) | (red << 16) | (green << 8) | blue;
}

uint32 Average(uint32 a0, uint32 a1, uint32 a2, uint32 a3) {
  uint32 alpha_and_green =
      ((((a0 >> 8) & 0x00ff00ff) +
        ((a1 >> 8) & 0x00ff00ff) +
        ((a2 >> 8) & 0x00ff00ff) +
        ((a3 >> 8) & 0x00ff00ff) + 0x00020002) >> 2) & 0x00ff00ff;
  uint32 red_and_blue =
      (((a0 & 0x00ff00ff) +
        (a1 & 0x00ff00ff) +
        (a2 & 0x00ff00ff) +
        (a3 & 0x00ff00ff) + 0x00020002) >> 2) & 0x00ff00ff;
  return (alpha_and_green << 8) | red_and_blue;
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

inline uint8 Paeth8(uint8 a, uint8 b, uint8 c) {
  int p = int(a) + int(b) - int(c);
  int pa = abs(p - a);
  int pb = abs(p - b);
  int pc = abs(p - c);
  if (pa <= pb && pa <= pc) {
    return a;
  } else if (pb <= pc) {
    return b;
  }
  return c;
}

uint32 Paeth32(uint32 a, uint32 b, uint32 c) {
  return
      (Paeth8(a >> 24, b >> 24, c >> 24) << 24) |
      (Paeth8((a >> 16) & 0xff,
              (b >> 16) & 0xff,
              (c >> 16) & 0xff) << 16) |
      (Paeth8((a >> 8) & 0xff,
              (b >> 8) & 0xff,
              (c >> 8) & 0xff) << 8) |
      Paeth8(a & 0xff, b & 0xff, c & 0xff);
}

inline int Clamp(int a) {
  if (a < 0) return 0;
  if (a > 255) return 255;
  return a;
}

inline int AddSubtractComponent(int mul, int a, int b, int c) {
  return Clamp(a + (mul * (b - c)) / 256);
}

uint32 ClampedAddSubtract(int mul, uint32 c0, uint32 c1, uint32 c2) {
  int a = AddSubtractComponent(mul,
                               c0 >> 24, c1 >> 24, c2 >> 24);
  int r = AddSubtractComponent(mul,
                               (c0 >> 16) & 0xff,
                               (c1 >> 16) & 0xff,
                               (c2 >> 16) & 0xff);
  int g = AddSubtractComponent(mul,
                               (c0 >> 8) & 0xff,
                               (c1 >> 8) & 0xff,
                               (c2 >> 8) & 0xff);
  int b = AddSubtractComponent(mul,
                               (c0 >> 0) & 0xff,
                               (c1 >> 0) & 0xff,
                               (c2 >> 0) & 0xff);
  return (a << 24) | (r << 16) | (g << 8) | b;
}


uint32 PredictValue(int mode, int x, int y, int xsize, const uint32 *argb) {
  if (x == 0) {
    if (y == 0) {
      return 0xff000000;
    }
    return argb[(y - 1) * xsize + x];
  }
  if (y == 0) {
    return argb[y * xsize + x - 1];
  }
  if (x == 1 && mode == 11) {
    mode = 1;
  }
  if (x >= xsize - 2 && mode == 12) {
    mode = 3;
  }
  if (x == xsize - 1 && (mode == 3 || mode == 5 || mode == 9 || mode == 10)) {
    mode = 2;
  }
  switch (mode) {
    case 1: return argb[y * xsize + x - 1];
    case 2: return argb[(y - 1) * xsize + x];
    case 3: return argb[(y - 1) * xsize + x + 1];
    case 4: return argb[(y - 1) * xsize + x - 1];

    case 5: return Average(argb[y * xsize + x - 1],
                           argb[(y - 1) * xsize + x],
                           argb[(y - 1) * xsize + x + 1]);
    case 6: return Average2(argb[y * xsize + x - 1],
                            argb[(y - 1) * xsize + x - 1]);
    case 7: return Average2(argb[y * xsize + x - 1],
                            argb[(y - 1) * xsize + x]);
    case 8: return Average2(argb[(y - 1) * xsize + x - 1],
                            argb[(y - 1) * xsize + x]);
    case 9: return Average2(argb[(y - 1) * xsize + x],
                            argb[(y - 1) * xsize + x + 1]);
    case 10: return Average(argb[y * xsize + x - 1],
                            argb[(y - 1) * xsize + x - 1],
                            argb[(y - 1) * xsize + x],
                            argb[(y - 1) * xsize + x + 1]);
    case 11: return argb[(y - 1) * xsize + x - 2];
    case 12: return argb[(y - 1) * xsize + x + 2];
    case 13: return Paeth32(argb[y * xsize + x - 1],
                            argb[(y - 1) * xsize + x],
                            argb[(y - 1) * xsize + x - 1]);
    case 14:
      return ClampedAddSubtract(256,
                                argb[y * xsize + x - 1],
                                argb[(y - 1) * xsize + x],
                                argb[(y - 1) * xsize + x - 1]);
    case 15:
      {
        uint32 ave = Average2(argb[y * xsize + x - 1],
                              argb[(y - 1) * xsize + x]);
        return ClampedAddSubtract(128, ave, ave, argb[(y - 1) * xsize + x - 1]);
      }
  }
  printf("Impossible %d\n", mode);
  abort();
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
      const uint32 predict =
          PredictValue(mode, all_x, all_y, xsize, from_argb);
      to_argb[ix] = Subtract(from_argb[ix], predict);
    }
  }
}

void PredictorInverseTransform(int xsize, int ysize, int bits,
                               const uint32* image,
                               const uint32 *from_argb,
                               uint32* to_argb) {
  const int tile_xsize = (xsize + (1 << bits) - 1) >> bits;
  for (int image_y = 0; image_y < ysize; ++image_y) {
    int tile_y = image_y >> bits;
    int ix = image_y * xsize;
    for (int tile_x = 0; tile_x < tile_xsize; ++tile_x) {
      int tile_ix = tile_y * tile_xsize + tile_x;
      int mode = ((image[tile_ix] >> 8) & 0xff) + 1;
      int image_x = (tile_x << bits);
      int xend = image_x + (1 << bits);
      if (xend > xsize) {
        xend = xsize;
      }
      for (; image_x < xend; ++image_x) {
        const uint32 predict =
            PredictValue(mode, image_x, image_y, xsize, to_argb);
        to_argb[ix + image_x] = Add(from_argb[ix + image_x], predict);
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
  int xscan = 1 << bits;
  if (xscan > xsize - (tile_x << bits)) {
    xscan = xsize - (tile_x << bits);
  }
  if (inverse) {
    for (int y = 0; y < (1 << bits); ++y) {
      int all_y = (tile_y << bits) + y;
      if (all_y >= ysize) {
        break;
      }
      int ix = all_y * xsize + (tile_x << bits);
      for (int x = 0; x < xscan; ++x, ++ix) {
        to_argb[ix] = color_transform.InverseTransform(from_argb[ix]);
      }
    }
  } else {
    for (int y = 0; y < (1 << bits); ++y) {
      int all_y = (tile_y << bits) + y;
      if (all_y >= ysize) {
        break;
      }
      int ix = all_y * xsize + (tile_x << bits);
      for (int x = 0; x < xscan; ++x, ++ix) {
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
  for (int tile_y = 0; tile_y < tile_ysize; ++tile_y) {
    for (int tile_x = 0; tile_x < tile_xsize; ++tile_x) {
      ColorSpaceTransformElement color_transform;
      int tile_ix = tile_y * tile_xsize + tile_x;
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
    argb_array[i] = (argb & 0xff00ff00) |
        (((argb & 0x00ff00ff) + green) & 0x00ff00ff);
  }
}
