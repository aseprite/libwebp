// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)

#include "predictor.h"

#include <stdlib.h>

#include "integral_types.h"

uint32 Average2(uint32 a0, uint32 a1) {
  int alpha = ((a0 >> 24) + (a1 >> 24) + 1) / 2;
  int red = (((a0 >> 16) & 0xff) + ((a1 >> 16) & 0xff) + 1) / 2;
  int green = (((a0 >> 8) & 0xff) + ((a1 >> 8) & 0xff) + 1) / 2;
  int blue = (((a0 >> 0) & 0xff) + ((a1 >> 0) & 0xff) + 1) / 2;
  return (alpha << 24) | (red << 16) | (green << 8) | blue;
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
  int alpha = ((a0 >> 24) + (a1 >> 24) + (a2 >> 24) + (a3 >> 24) + 2) / 4;
  int red = (((a0 >> 16) & 0xff) +
             ((a1 >> 16) & 0xff) +
             ((a2 >> 16) & 0xff) +
             ((a3 >> 16) & 0xff) + 2) / 4;
  int green = (((a0 >> 8) & 0xff) +
               ((a1 >> 8) & 0xff) +
               ((a2 >> 8) & 0xff) +
               ((a3 >> 8) & 0xff) + 2) / 4;
  int blue = (((a0 >> 0) & 0xff) +
              ((a1 >> 0) & 0xff) +
              ((a2 >> 0) & 0xff) +
              ((a3 >> 0) & 0xff) + 2) / 4;
  return (alpha << 24) | (red << 16) | (green << 8) | blue;
}

uint32 Add(uint32 a, uint32 b) {
  int alpha = ((a >> 24) & 0xff) + ((b >> 24) & 0xff);
  int red = ((a >> 16) & 0xff) + ((b >> 16) & 0xff);
  int green = ((a >> 8) & 0xff) + ((b >> 8) & 0xff);
  int blue = ((a >> 0) & 0xff) + ((b >> 0) & 0xff);
  alpha &= 0xff;
  red &= 0xff;
  green &= 0xff;
  blue &= 0xff;
  return (alpha << 24) | (red << 16) | (green << 8) | blue;
}

uint32 Subtract(uint32 a, uint32 b) {
  int alpha = ((a >> 24) & 0xff) - ((b >> 24) & 0xff);
  int red = ((a >> 16) & 0xff) - ((b >> 16) & 0xff);
  int green = ((a >> 8) & 0xff) - ((b >> 8) & 0xff);
  int blue = ((a >> 0) & 0xff) - ((b >> 0) & 0xff);
  alpha &= 0xff;
  red &= 0xff;
  green &= 0xff;
  blue &= 0xff;
  return (alpha << 24) | (red << 16) | (green << 8) | blue;
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
  uint32 retval = 0;
  for (int i = 0; i < 4; ++i) {
    uint32 v = Paeth8((a >> (i * 8)) & 0xff,
                      (b >> (i * 8)) & 0xff,
                      (c >> (i * 8)) & 0xff);
    retval |= v << (i * 8);
  }
  return retval;
}

int Clamp(int a) {
  if (a < 0) return 0;
  if (a > 255) return 255;
  return a;
}

int AddSubtractComponent(int mul, int a, int b, int c) {
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
  if (mode == 0) {
    return 0xff000000;
  }
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
  for (int y = 0; y < (1 << bits); ++y) {
    int all_y = (tile_y << bits) + y;
    if (all_y >= ysize) {
      break;
    }
    for (int x = 0; x < (1 << bits); ++x) {
      int all_x = (tile_x << bits) + x;
      if (all_x >= xsize) {
        break;
      }
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
  int tile_xsize = (xsize + (1 << bits) - 1) >> bits;
  for (int y = 0; y < ysize; ++y) {
    int tile_y = y >> bits;
    for (int x = 0; x < xsize; ++x) {
      int tile_x = x >> bits;
      int tile_ix = tile_y * tile_xsize + tile_x;
      int mode = ((image[tile_ix] >> 8) & 0xff) + 1;
      int ix = y * xsize + x;
      const uint32 predict = PredictValue(mode, x, y, xsize, to_argb);
      to_argb[ix] = Add(from_argb[ix], predict);
    }
  }
}

void CopyTileWithColorTransform(int xsize, int ysize,
                                const uint32 *from_argb,
                                int tile_x, int tile_y, int bits,
                                ColorSpaceTransformElement color_transform,
                                bool inverse,
                                uint32 *to_argb) {
  for (int y = 0; y < (1 << bits); ++y) {
    int all_y = (tile_y << bits) + y;
    if (all_y >= ysize) {
      break;
    }
    for (int x = 0; x < (1 << bits); ++x) {
      int all_x = (tile_x << bits) + x;
      if (all_x >= xsize) {
        break;
      }
      const int ix = all_y * xsize + all_x;
      if (inverse) {
        to_argb[ix] = color_transform.InverseTransform(from_argb[ix]);
      } else {
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
    uint32 green = (argb >> 8) & 0xff;
    uint32 new_r = (((argb >> 16) & 0xff) + green) & 0xff;
    uint32 new_b = ((argb & 0xff) + green) & 0xff;
    argb_array[i] = (argb & 0xff00ff00) | (new_r << 16) | new_b;
  }
}
