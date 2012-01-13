// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)

#include "./predictor.h"

#include <limits.h>
#include <stdlib.h>
#include <math.h>

#include <algorithm>

#include "./backward_references.h"
#include "./histogram.h"
#include "../common/integral_types.h"
#include "../common/predictor.h"

static double PredictionCostSpatial(int *counts, int weight_0,
                                    double exp_val) {
  const int significant_symbols = 16;
  const double exp_decay_factor = 0.6;
  double bits = weight_0 * counts[0];
  for (int i = 1; i < significant_symbols; ++i) {
    bits += exp_val * (counts[i] + counts[256 - i]);
    exp_val *= exp_decay_factor;
  }
  return -0.1 * bits;
}

static double PredictionCostSpatialHistogram(Histogram *accumulated,
                                             Histogram *tile) {
  const double exp_val = 0.94;
  Histogram combo(0);
  combo.Add(*accumulated);
  combo.Add(*tile);
  return
      tile->EstimateBitsBulk() +
      combo.EstimateBitsBulk() +
      PredictionCostSpatial(tile->alpha_, 1, exp_val) +
      PredictionCostSpatial(tile->literal_, 1, exp_val) +
      PredictionCostSpatial(tile->red_, 1, exp_val) +
      PredictionCostSpatial(tile->blue_, 1, exp_val);
}

static int GetBestPredictorForTile(int tile_x, int tile_y, int max_tile_size,
                                   int xsize, int ysize,
                                   Histogram *accumulated,
                                   const uint32 *argb) {
  const int num_pred_modes = 14;
  const int tile_y_offset = tile_y * max_tile_size;
  const int tile_x_offset = tile_x * max_tile_size;
  int all_x_max = tile_x_offset + max_tile_size;
  if (all_x_max > xsize) {
    all_x_max = xsize;
  }
  int all_y_max = tile_y_offset + max_tile_size;
  if (all_y_max > ysize) {
    all_y_max = ysize;
  }
  double best_diff = 1e99;
  int best_mode = 0;
  for (int mode = 0; mode < num_pred_modes; ++mode) {
    Histogram histo(0);  // 0 is for only 1 (unused) palette value.
    for (int all_y = tile_y_offset; all_y < all_y_max; ++all_y) {
      for (int all_x = tile_x_offset; all_x < all_x_max; ++all_x) {
        uint32 predict;
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
        const uint32 predict_diff =
            Subtract(argb[all_y * xsize + all_x], predict);
        histo.AddSingleLiteralOrCopy(
            LiteralOrCopy::CreateLiteral(predict_diff));
      }
    }
    const double cur_diff = PredictionCostSpatialHistogram(accumulated, &histo);
    if (cur_diff < best_diff) {
      best_diff = cur_diff;
      best_mode = mode;
    }
  }
  return best_mode;
}

void PredictorImage(int xsize, int ysize, int bits,
                    const uint32 *from_argb, uint32 *to_argb, uint32 *image) {
  uint32 *argb_orig = reinterpret_cast<uint32 *>(
      malloc(sizeof(from_argb[0]) * xsize * ysize));
  memcpy(argb_orig, from_argb, sizeof(from_argb[0]) * xsize * ysize);
  const int max_tile_size = 1 << bits;
  const int tile_xsize = (xsize + max_tile_size - 1) >> bits;
  const int tile_ysize = (ysize + max_tile_size - 1) >> bits;
  Histogram histo(0);
  for (int tile_y = 0; tile_y < tile_ysize; ++tile_y) {
    const int tile_y_offset = tile_y * max_tile_size;
    for (int tile_x = 0; tile_x < tile_xsize; ++tile_x) {
      const int tile_x_offset = tile_x * max_tile_size;
      int all_x_max = tile_x_offset + max_tile_size;
      if (all_x_max > xsize) {
        all_x_max = xsize;
      }
      int pred = GetBestPredictorForTile(tile_x, tile_y, max_tile_size,
                                         xsize, ysize, &histo, from_argb);
      image[tile_y * tile_xsize + tile_x] = 0xff000000 | (pred << 8);
      CopyTileWithPrediction(xsize, ysize, from_argb,
                             tile_x, tile_y, bits, pred,
                             to_argb);
      for (int y = 0; y < max_tile_size; ++y) {
        int all_y = tile_y_offset + y;
        if (all_y >= ysize) {
          break;
        }
        int ix = all_y * xsize + tile_x_offset;
        for (int all_x = tile_x_offset; all_x < all_x_max; ++all_x, ++ix) {
          histo.AddSingleLiteralOrCopy(
              LiteralOrCopy::CreateLiteral(to_argb[ix]));
        }
      }
    }
  }
#ifndef NDEBUG
  uint32 *argb = reinterpret_cast<uint32 *>(
      malloc(sizeof(to_argb[0]) * xsize * ysize));
  memcpy(argb, to_argb, sizeof(to_argb[0]) * xsize * ysize);
  PredictorInverseTransform(xsize, ysize, bits, image,
                            &argb[0], &argb[0]);
  for (int i = 0; i < xsize * ysize; ++i) {
    VERIFY(argb[i] == argb_orig[i]);
  }
  free(argb);
#endif
  free(argb_orig);
}

static double PredictionCostCrossColor(int *accumulated, int *counts) {
  // Favor low entropy, locally and globally.
  const int length = 256;
  int combo[256];
  for (int i = 0; i < length; ++i) {
    combo[i] = accumulated[i] + counts[i];
  }
  double bits = BitsEntropy(&combo[0], length) + BitsEntropy(counts, length);
  // Favor small absolute values.
  bits += PredictionCostSpatial(counts, 3, 2.4);
  return bits;
}

inline bool SkipRepeatedPixels(const uint32 *argb, int ix, int xsize) {
  const uint32 v = argb[ix];
  if (ix >= xsize + 3) {
    if (v == argb[ix - xsize] &&
        argb[ix - 1] == argb[ix - xsize - 1] &&
        argb[ix - 2] == argb[ix - xsize - 2] &&
        argb[ix - 3] == argb[ix - xsize - 3]) {
      return true;
    }
    return v == argb[ix - 3] &&
        v == argb[ix - 2] &&
        v == argb[ix - 1];
  } else if (ix >= 3) {
    return v == argb[ix - 3] &&
        v == argb[ix - 2] &&
        v == argb[ix - 1];
  }
  return false;
}

static ColorSpaceTransformElement GetBestColorTransformForTile(
    int tile_x, int tile_y, int bits,
    ColorSpaceTransformElement prevX,
    ColorSpaceTransformElement prevY,
    int quality, int xsize, int ysize,
    int *accumulated_red_histo,
    int *accumulated_blue_histo,
    const uint32 *argb) {
  double best_diff = 1e99;
  ColorSpaceTransformElement best_tx;
  best_tx.Clear();
  const int step = (quality == 0) ? 16 : 8;
  const int halfstep = step / 2;
  const int max_tile_size = 1 << bits;
  const int tile_y_offset = tile_y * max_tile_size;
  const int tile_x_offset = tile_x * max_tile_size;
  int all_x_max = tile_x_offset + max_tile_size;
  if (all_x_max > xsize) {
    all_x_max = xsize;
  }
  int all_y_max = tile_y_offset + max_tile_size;
  if (all_y_max > xsize) {
    all_y_max = xsize;
  }
  for (int green_to_red = -64; green_to_red <= 64; green_to_red += halfstep) {
    ColorSpaceTransformElement tx;
    tx.Clear();
    tx.green_to_red_ = green_to_red & 0xff;

    int histo[256] = { 0 };
    for (int all_y = tile_y_offset; all_y < all_y_max; ++all_y) {
      int ix = all_y * xsize + tile_x_offset;
      for (int all_x = tile_x_offset; all_x < all_x_max; ++all_x, ++ix) {
        if (SkipRepeatedPixels(argb, ix, xsize)) {
          continue;
        }
        const uint32 predict = tx.Transform(argb[ix]);
        ++histo[(predict >> 16) & 0xff];  // red.
      }
    }
    double cur_diff =
        PredictionCostCrossColor(&accumulated_red_histo[0], &histo[0]);
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
  const int green_to_red = best_tx.green_to_red_;
  for (int green_to_blue = -32; green_to_blue <= 32; green_to_blue += step) {
    for (int red_to_blue = -32; red_to_blue <= 32; red_to_blue += step) {
      ColorSpaceTransformElement tx;
      tx.Clear();
      tx.green_to_red_ = green_to_red;
      tx.green_to_blue_ = green_to_blue;
      tx.red_to_blue_ = red_to_blue;
      int histo[256] = { 0 };
      for (int all_y = tile_y_offset; all_y < all_y_max; ++all_y) {
        int ix = all_y * xsize + tile_x_offset;
        for (int all_x = tile_x_offset; all_x < all_x_max; ++all_x, ++ix) {
          if (SkipRepeatedPixels(argb, ix, xsize)) {
            continue;
          }
          const uint32 predict = tx.Transform(argb[ix]);
          ++histo[predict & 0xff];  // blue.
        }
      }
      double cur_diff =
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
                         const uint32 *from_argb, int quality,
                         uint32 *to_argb, uint32 *image) {
  uint32 *argb_orig = reinterpret_cast<uint32 *>(
      malloc(sizeof(from_argb[0]) * xsize * ysize));
  memcpy(argb_orig, from_argb, sizeof(from_argb[0]) * xsize * ysize);
  const int max_tile_size = 1 << bits;
  int tile_xsize = (xsize + max_tile_size - 1) >> bits;
  int tile_ysize = (ysize + max_tile_size - 1) >> bits;
  int accumulated_red_histo[256] = { 0 };
  int accumulated_blue_histo[256] = { 0 };
  ColorSpaceTransformElement prevX;
  ColorSpaceTransformElement prevY;
  prevY.Clear();
  prevX.Clear();
  for (int tile_y = 0; tile_y < tile_ysize; ++tile_y) {
    for (int tile_x = 0; tile_x < tile_xsize; ++tile_x) {
      const int tile_y_offset = tile_y * max_tile_size;
      const int tile_x_offset = tile_x * max_tile_size;
      if (tile_y != 0) {
        prevY.InitFromCode(image[(tile_y - 1) * tile_xsize + tile_x]);
        prevX.InitFromCode(image[tile_y * tile_xsize + tile_x - 1]);
      } else if (tile_x != 0) {
        prevX.InitFromCode(image[tile_y * tile_xsize + tile_x - 1]);
      }
      const ColorSpaceTransformElement color_transform =
          GetBestColorTransformForTile(tile_x, tile_y, bits,
                                       prevX, prevY,
                                       quality, xsize, ysize,
                                       &accumulated_red_histo[0],
                                       &accumulated_blue_histo[0],
                                       from_argb);
      image[tile_y * tile_xsize + tile_x] = color_transform.Code();
      CopyTileWithColorTransform(xsize, ysize, from_argb,
                                 tile_x, tile_y, bits, color_transform,
                                 false,  // forward transform
                                 to_argb);
      // Gather accumulated histogram data.

      int all_x_max = tile_x_offset + max_tile_size;
      if (all_x_max > xsize) {
        all_x_max = xsize;
      }
      for (int y = 0; y < max_tile_size; ++y) {
        int all_y = tile_y_offset + y;
        if (all_y >= ysize) {
          break;
        }
        int ix = all_y * xsize + tile_x_offset;
        for (int all_x = tile_x_offset; all_x < all_x_max; ++all_x, ++ix) {
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
  uint32 *argb = reinterpret_cast<uint32 *>(
      malloc(sizeof(to_argb[0]) * xsize * ysize));
  memcpy(argb, to_argb, sizeof(to_argb[0]) * xsize * ysize);
  ColorSpaceInverseTransform(xsize, ysize, bits, image, &argb[0], &argb[0]);
  for (int i = 0; i < xsize * ysize; ++i) {
    VERIFY(argb[i] == argb_orig[i]);
  }
  free(argb);
#endif
  free(argb_orig);
}

void SubtractGreenFromBlueAndRed(int n, uint32 *argb_array) {
  for (int i = 0; i < n; ++i) {
    uint32 argb = argb_array[i];
    uint32 green = (argb >> 8) & 0xff;
    uint32 new_r = (((argb >> 16) & 0xff) - green) & 0xff;
    uint32 new_b = ((argb & 0xff) - green) & 0xff;
    argb_array[i] = (argb & 0xff00ff00) | (new_r << 16) | new_b;
  }
}
