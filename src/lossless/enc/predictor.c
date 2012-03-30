// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)

#include "./predictor.h"

#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <math.h>

#include "./backward_references.h"
#include "./histogram.h"
#include "../common/predictor.h"

static double PredictionCostSpatial(const int* counts,
                                    int weight_0, double exp_val) {
  const int significant_symbols = 16;
  const double exp_decay_factor = 0.6;
  double bits = weight_0 * counts[0];
  int i;
  for (i = 1; i < significant_symbols; ++i) {
    bits += exp_val * (counts[i] + counts[256 - i]);
    exp_val *= exp_decay_factor;
  }
  return -0.1 * bits;
}

static double PredictionCostSpatialHistogram(int accumulated[4][256],
                                             int tile[4][256]) {
  int i;
  int k;
  int combo[256];
  double retval = 0;
  for (i = 0; i < 4; ++i) {
    const double exp_val = 0.94;
    retval += PredictionCostSpatial(&tile[i][0], 1, exp_val);
    retval += ShannonEntropy(&tile[i][0], 256);
    for (k = 0; k < 256; ++k) {
      combo[k] = accumulated[i][k] + tile[i][k];
    }
    retval += ShannonEntropy(&combo[0], 256);
  }
  return retval;
}

static int GetBestPredictorForTile(int tile_x, int tile_y, int max_tile_size,
                                   int xsize, int ysize,
                                   int accumulated[4][256],
                                   const uint32_t* argb) {
  const int num_pred_modes = 14;
  const int tile_y_offset = tile_y * max_tile_size;
  const int tile_x_offset = tile_x * max_tile_size;
  double cur_diff;
  double best_diff = 1e99;
  int best_mode = 0;
  int mode;
  int all_x_max = tile_x_offset + max_tile_size;
  int all_y_max = tile_y_offset + max_tile_size;
  //  int mask = 0x827;
  //  int mask = 0x3fff;
  int histo[4][256];
  if (all_x_max > xsize) {
    all_x_max = xsize;
  }
  if (all_y_max > ysize) {
    all_y_max = ysize;
  }
  for (mode = 0; mode < num_pred_modes; ++mode) {
    int all_y;
    // TODO(jyrki): add a mask control here to have more control for fastest
    // compression (with mask >>= 1 in the loop).
    //    if ((mask & 1) == 0) {
    //      continue;
    //    }
    memset(&histo[0][0], 0, sizeof(histo));
    for (all_y = tile_y_offset; all_y < all_y_max; ++all_y) {
      int all_x;
      for (all_x = tile_x_offset; all_x < all_x_max; ++all_x) {
        uint32_t predict;
        uint32_t predict_diff;
        if (all_y == 0) {
          if (all_x == 0) {
            predict = 0xff000000;
          } else {
            predict = argb[all_x - 1];  // Top Row: Pick Left Element.
          }
        } else if (all_x == 0) {
          predict = argb[(all_y - 1) * xsize];  // First Col: Pick Top Element.
        } else {
          predict = PredictValue(mode, xsize, argb + all_y * xsize + all_x);
        }
        predict_diff = Subtract(argb[all_y * xsize + all_x], predict);
        ++histo[0][predict_diff >> 24];
        ++histo[1][((predict_diff >> 16) & 0xff)];
        ++histo[2][((predict_diff >> 8) & 0xff)];
        ++histo[3][(predict_diff & 0xff)];
      }
    }
    cur_diff = PredictionCostSpatialHistogram(accumulated, histo);
    if (cur_diff < best_diff) {
      best_diff = cur_diff;
      best_mode = mode;
    }
  }
  return best_mode;
}

void PredictorImage(int xsize, int ysize, int bits, const uint32_t* from_argb,
                    uint32_t* to_argb, uint32_t* image) {
  const int max_tile_size = 1 << bits;
  const int tile_xsize = (xsize + max_tile_size - 1) >> bits;
  const int tile_ysize = (ysize + max_tile_size - 1) >> bits;
  int tile_y;
  int histo[4][256];
#ifndef NDEBUG
  uint32_t* argb_orig = (uint32_t*)
      malloc(sizeof(from_argb[0]) * xsize * ysize);
  memcpy(argb_orig, from_argb, sizeof(from_argb[0]) * xsize * ysize);
#endif
  memset(histo, 0, sizeof(histo));
  for (tile_y = 0; tile_y < tile_ysize; ++tile_y) {
    const int tile_y_offset = tile_y * max_tile_size;
    int tile_x;
    for (tile_x = 0; tile_x < tile_xsize; ++tile_x) {
      int pred;
      int y;
      const int tile_x_offset = tile_x * max_tile_size;
      int all_x_max = tile_x_offset + max_tile_size;
      if (all_x_max > xsize) {
        all_x_max = xsize;
      }
      pred = GetBestPredictorForTile(tile_x, tile_y, max_tile_size,
                                     xsize, ysize, histo, from_argb);
      image[tile_y * tile_xsize + tile_x] = 0xff000000 | (pred << 8);
      CopyTileWithPrediction(xsize, ysize, from_argb,
                             tile_x, tile_y, bits, pred,
                             to_argb);
      for (y = 0; y < max_tile_size; ++y) {
        int ix;
        int all_x;
        int all_y = tile_y_offset + y;
        if (all_y >= ysize) {
          break;
        }
        ix = all_y * xsize + tile_x_offset;
        for (all_x = tile_x_offset; all_x < all_x_max; ++all_x, ++ix) {
          const uint32_t a = to_argb[ix];
          ++histo[0][a >> 24];
          ++histo[1][((a >> 16) & 0xff)];
          ++histo[2][((a >> 8) & 0xff)];
          ++histo[3][(a & 0xff)];
        }
      }
    }
  }
#ifndef NDEBUG
  {
    int i;
    uint32_t* argb = (uint32_t*)malloc(sizeof(to_argb[0]) * xsize * ysize);
    memcpy(argb, to_argb, sizeof(to_argb[0]) * xsize * ysize);
    PredictorInverseTransform(xsize, ysize, bits, image,
                              &argb[0], &argb[0]);
    for (i = 0; i < xsize * ysize; ++i) {
      assert(argb[i] == argb_orig[i]);
    }
    free(argb);
  }
  free(argb_orig);
#endif
}

static double PredictionCostCrossColor(int* accumulated, int* counts) {
  // Favor low entropy, locally and globally.
  const int length = 256;
  int combo[256];
  int i;
  for (i = 0; i < length; ++i) {
    combo[i] = accumulated[i] + counts[i];
  }
  return ShannonEntropy(&combo[0], length) +
      ShannonEntropy(counts, length) +
      PredictionCostSpatial(counts, 3, 2.4); // Favor small absolute values.
}

static inline int SkipRepeatedPixels(const uint32_t* argb, int ix, int xsize) {
  const uint32_t v = argb[ix];
  if (ix >= xsize + 3) {
    if (v == argb[ix - xsize] &&
        argb[ix - 1] == argb[ix - xsize - 1] &&
        argb[ix - 2] == argb[ix - xsize - 2] &&
        argb[ix - 3] == argb[ix - xsize - 3]) {
      return 1;
    }
    return v == argb[ix - 3] &&
        v == argb[ix - 2] &&
        v == argb[ix - 1];
  } else if (ix >= 3) {
    return v == argb[ix - 3] &&
        v == argb[ix - 2] &&
        v == argb[ix - 1];
  }
  return 0;
}

static ColorSpaceTransformElement GetBestColorTransformForTile(
    int tile_x, int tile_y, int bits,
    ColorSpaceTransformElement prevX,
    ColorSpaceTransformElement prevY,
    int quality, int xsize, int ysize,
    int* accumulated_red_histo,
    int* accumulated_blue_histo,
    const uint32_t* argb) {
  double best_diff = 1e99;
  double cur_diff;
  const int step = (quality == 0) ? 32 : 8;
  const int halfstep = step / 2;
  const int max_tile_size = 1 << bits;
  const int tile_y_offset = tile_y * max_tile_size;
  const int tile_x_offset = tile_x * max_tile_size;
  int green_to_red;
  int green_to_blue;
  int red_to_blue;
  int all_x_max = tile_x_offset + max_tile_size;
  int all_y_max = tile_y_offset + max_tile_size;
  ColorSpaceTransformElement best_tx;
  ColorSpaceTransformElementClear(&best_tx);
  if (all_x_max > xsize) {
    all_x_max = xsize;
  }
  if (all_y_max > ysize) {
    all_y_max = ysize;
  }
  for (green_to_red = -64; green_to_red <= 64; green_to_red += halfstep) {
    int histo[256] = { 0 };
    int all_y;
    ColorSpaceTransformElement tx;
    ColorSpaceTransformElementClear(&tx);
    tx.green_to_red_ = green_to_red & 0xff;

    for (all_y = tile_y_offset; all_y < all_y_max; ++all_y) {
      uint32_t predict;
      int ix = all_y * xsize + tile_x_offset;
      int all_x;
      for (all_x = tile_x_offset; all_x < all_x_max; ++all_x, ++ix) {
        if (SkipRepeatedPixels(argb, ix, xsize)) {
          continue;
        }
        predict = ColorSpaceTransformElementTransform(&tx, argb[ix]);
        ++histo[(predict >> 16) & 0xff];  // red.
      }
    }
    cur_diff = PredictionCostCrossColor(&accumulated_red_histo[0], &histo[0]);
    if (tx.green_to_red_ == prevX.green_to_red_) {
      cur_diff -= 3;  // favor keeping the areas locally similar
    }
    if (tx.green_to_red_ == prevY.green_to_red_) {
      cur_diff -= 3;  // favor keeping the areas locally similar
    }
    if (tx.green_to_red_ == 0) {
      cur_diff -= 3;
    }
    if (cur_diff < best_diff) {
      best_diff = cur_diff;
      best_tx = tx;
    }
  }
  best_diff = 1e99;
  green_to_red = best_tx.green_to_red_;
  for (green_to_blue = -32; green_to_blue <= 32; green_to_blue += step) {
    for (red_to_blue = -32; red_to_blue <= 32; red_to_blue += step) {
      int all_y;
      int histo[256] = { 0 };
      ColorSpaceTransformElement tx;
      tx.green_to_red_ = green_to_red;
      tx.green_to_blue_ = green_to_blue;
      tx.red_to_blue_ = red_to_blue;
      for (all_y = tile_y_offset; all_y < all_y_max; ++all_y) {
        uint32_t predict;
        int all_x;
        int ix = all_y * xsize + tile_x_offset;
        for (all_x = tile_x_offset; all_x < all_x_max; ++all_x, ++ix) {
          if (SkipRepeatedPixels(argb, ix, xsize)) {
            continue;
          }
          predict = ColorSpaceTransformElementTransform(&tx, argb[ix]);
          ++histo[predict & 0xff];  // blue.
        }
      }
      cur_diff =
        PredictionCostCrossColor(&accumulated_blue_histo[0], &histo[0]);
      if (tx.green_to_blue_ == prevX.green_to_blue_) {
        cur_diff -= 3;  // favor keeping the areas locally similar
      }
      if (tx.green_to_blue_ == prevY.green_to_blue_) {
        cur_diff -= 3;  // favor keeping the areas locally similar
      }
      if (tx.red_to_blue_ == prevX.red_to_blue_) {
        cur_diff -= 3;  // favor keeping the areas locally similar
      }
      if (tx.red_to_blue_ == prevY.red_to_blue_) {
        cur_diff -= 3;  // favor keeping the areas locally similar
      }
      if (tx.green_to_blue_ == 0) {
        cur_diff -= 3;
      }
      if (tx.red_to_blue_ == 0) {
        cur_diff -= 3;
      }
      if (cur_diff < best_diff) {
        best_diff = cur_diff;
        best_tx = tx;
      }
    }
  }
  return best_tx;
}

void ColorSpaceTransform(int xsize, int ysize, int bits,
                         const uint32_t* from_argb, int quality,
                         uint32_t* to_argb, uint32_t* image) {
  const int max_tile_size = 1 << bits;
  int tile_xsize = (xsize + max_tile_size - 1) >> bits;
  int tile_ysize = (ysize + max_tile_size - 1) >> bits;
  int accumulated_red_histo[256] = { 0 };
  int accumulated_blue_histo[256] = { 0 };
  int tile_y;
  int tile_x;
  ColorSpaceTransformElement prevX;
  ColorSpaceTransformElement prevY;
#ifndef NDEBUG
  uint32_t* argb_orig = (uint32_t*)malloc(sizeof(from_argb[0]) * xsize * ysize);
  memcpy(argb_orig, from_argb, sizeof(from_argb[0]) * xsize * ysize);
#endif
  ColorSpaceTransformElementClear(&prevY);
  ColorSpaceTransformElementClear(&prevX);
  for (tile_y = 0; tile_y < tile_ysize; ++tile_y) {
    for (tile_x = 0; tile_x < tile_xsize; ++tile_x) {
      ColorSpaceTransformElement color_transform;
      int all_x_max;
      int y;
      const int tile_y_offset = tile_y * max_tile_size;
      const int tile_x_offset = tile_x * max_tile_size;
      if (tile_y != 0) {
        ColorSpaceTransformElementInitFromCode(
            &prevY, image[(tile_y - 1) * tile_xsize + tile_x]);
        ColorSpaceTransformElementInitFromCode(
            &prevX, image[tile_y * tile_xsize + tile_x - 1]);
      } else if (tile_x != 0) {
        ColorSpaceTransformElementInitFromCode(
            &prevX, image[tile_y * tile_xsize + tile_x - 1]);
      }
      color_transform =
          GetBestColorTransformForTile(tile_x, tile_y, bits,
                                       prevX, prevY,
                                       quality, xsize, ysize,
                                       &accumulated_red_histo[0],
                                       &accumulated_blue_histo[0],
                                       from_argb);
      image[tile_y * tile_xsize + tile_x] =
          ColorSpaceTransformElementCode(&color_transform);
      CopyTileWithColorTransform(xsize, ysize, from_argb,
                                 tile_x, tile_y, bits, color_transform,
                                 0,  // forward transform
                                 to_argb);

      // Gather accumulated histogram data.
      all_x_max = tile_x_offset + max_tile_size;
      if (all_x_max > xsize) {
        all_x_max = xsize;
      }
      for (y = 0; y < max_tile_size; ++y) {
        int ix;
        int all_x;
        int all_y = tile_y_offset + y;
        if (all_y >= ysize) {
          break;
        }
        ix = all_y * xsize + tile_x_offset;
        for (all_x = tile_x_offset; all_x < all_x_max; ++all_x, ++ix) {
          if (ix >= 2 &&
              to_argb[ix] == to_argb[ix - 2] &&
              to_argb[ix] == to_argb[ix - 1]) {
            continue;  // repeated pixels are handled by backward references
          }
          if (ix >= xsize + 2 &&
              to_argb[ix - 2] == to_argb[ix - xsize - 2] &&
              to_argb[ix - 1] == to_argb[ix - xsize - 1] &&
              to_argb[ix] == to_argb[ix - xsize]) {
            continue;  // repeated pixels are handled by backward references
          }
          ++accumulated_red_histo[(to_argb[ix] >> 16) & 0xff];
          ++accumulated_blue_histo[to_argb[ix] & 0xff];
        }
      }
    }
  }
#ifndef NDEBUG
  {
    int i;
    uint32_t* argb = (uint32_t*)malloc(sizeof(to_argb[0]) * xsize * ysize);
    memcpy(argb, to_argb, sizeof(to_argb[0]) * xsize * ysize);
    ColorSpaceInverseTransform(xsize, ysize, bits, image, &argb[0], &argb[0]);
    for (i = 0; i < xsize * ysize; ++i) {
      assert(argb[i] == argb_orig[i]);
    }
    free(argb);
  }
  free(argb_orig);
#endif
}

void SubtractGreenFromBlueAndRed(int n, uint32_t* argb_array) {
  int i;
  for (i = 0; i < n; ++i) {
    uint32_t argb = argb_array[i];
    uint32_t green = (argb >> 8) & 0xff;
    uint32_t new_r = (((argb >> 16) & 0xff) - green) & 0xff;
    uint32_t new_b = ((argb & 0xff) - green) & 0xff;
    argb_array[i] = (argb & 0xff00ff00) | (new_r << 16) | new_b;
  }
}
