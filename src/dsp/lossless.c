// Copyright 2012 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Image transforms and color space conversion methods for lossless decoder.
//
// Author: Vikas Arora (vikaas.arora@gmail.com)
//         jyrki@google.com (Jyrki Alakuijala)

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#include <stdlib.h>
#include "./lossless.h"
#include "../dec/vp8li.h"

static WEBP_INLINE uint32_t Average2(uint32_t a0, uint32_t a1) {
  return (((a0 ^ a1) & 0xfefefefeL) >> 1) + (a0 & a1);
}

static WEBP_INLINE uint32_t Average3(uint32_t a0, uint32_t a1, uint32_t a2) {
  return Average2(Average2(a0, a2), a1);
}

static WEBP_INLINE uint32_t Average4(uint32_t a0, uint32_t a1,
                                     uint32_t a2, uint32_t a3) {
  return Average2(Average2(a0, a1), Average2(a2, a3));
}

static WEBP_INLINE uint32_t Clip255(uint32_t a) {
  if (a < NUM_CODES_PER_BYTE) {
    return a;
  }
  // return 0, when a is a negative integer.
  // return 255, when a is positive.
  return ~a >> 24;
}

static WEBP_INLINE int AddSubtractComponentFull(int a, int b, int c) {
  return Clip255(a + b - c);
}

static WEBP_INLINE uint32_t ClampedAddSubtractFull(uint32_t c0, uint32_t c1,
                                                   uint32_t c2) {
  const int a = AddSubtractComponentFull(c0 >> 24, c1 >> 24, c2 >> 24);
  const int r = AddSubtractComponentFull((c0 >> 16) & 0xff,
                                         (c1 >> 16) & 0xff,
                                         (c2 >> 16) & 0xff);
  const int g = AddSubtractComponentFull((c0 >> 8) & 0xff,
                                         (c1 >> 8) & 0xff,
                                         (c2 >> 8) & 0xff);
  const int b = AddSubtractComponentFull(c0 & 0xff, c1 & 0xff, c2 & 0xff);
  return (a << 24) | (r << 16) | (g << 8) | b;
}

static WEBP_INLINE int AddSubtractComponentHalf(int a, int b) {
  return Clip255(a + (a - b) / 2);
}

static WEBP_INLINE uint32_t ClampedAddSubtractHalf(uint32_t c0, uint32_t c1,
                                                   uint32_t c2) {
  const uint32_t ave = Average2(c0, c1);
  const int a = AddSubtractComponentHalf(ave >> 24, c2 >> 24);
  const int r = AddSubtractComponentHalf((ave >> 16) & 0xff, (c2 >> 16) & 0xff);
  const int g = AddSubtractComponentHalf((ave >> 8) & 0xff, (c2 >> 8) & 0xff);
  const int b = AddSubtractComponentHalf((ave >> 0) & 0xff, (c2 >> 0) & 0xff);
  return (a << 24) | (r << 16) | (g << 8) | b;
}

static WEBP_INLINE uint32_t Select(uint32_t a, uint32_t b, uint32_t c) {
  const int p0 = (int)(a >> 24) + (int)(b >> 24) - (int)(c >> 24);
  const int p1 = (int)((a >> 16) & 0xff) + (int)((b >> 16) & 0xff) -
      (int)((c >> 16) & 0xff);
  const int p2 = (int)((a >> 8) & 0xff) + (int)((b >> 8) & 0xff) -
      (int)((c >> 8) & 0xff);
  const int p3 = (int)(a & 0xff) + (int)(b & 0xff) - (int)(c & 0xff);
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

static WEBP_INLINE argb_t PredictValue(uint32_t pred_mode, int xsize,
                                       const argb_t* const argb) {
  switch(pred_mode) {
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

void VP8LPredictorInverseTransform(const VP8LTransform* const transform,
                                   argb_t* const data) {
  size_t row, col, col_start, tile_offset;
  uint32_t pix_ix = 0;
  const uint32_t tile_size = 1 << transform->bits_;
  const uint32_t tiles_per_row = VP8LSubSampleSize(transform->xsize_,
                                               transform->bits_);
  // First Row follows the L (mode=1) mode.
  data[0] = VP8LAddPixels(data[0], ARGB_BLACK);
  for (col = 1; col < transform->xsize_; ++col) {
    data[col] = VP8LAddPixels(data[col], data[col - 1]);
  }
  pix_ix += transform->xsize_;

  for (row = 1; row < transform->ysize_; ++row) {
    const uint32_t tile_base_ix = tiles_per_row * (row >> transform->bits_);
    for (tile_offset = 0, col_start = 0; tile_offset < tiles_per_row;
         ++tile_offset, col_start += tile_size) {
      argb_t pred;
      // Pick the appropriate predictor mode (at start of every tile).
      const uint32_t pred_mode =
          (transform->data_[tile_base_ix + tile_offset] >> 8) & 0xff;
      uint32_t col_end = col_start + tile_size;
      if (col_end > transform->xsize_) col_end = transform->xsize_;

      // First col follows the T (mode=2) mode.
      pred = (col_start == 0) ? data[pix_ix - transform->xsize_] :
          PredictValue(pred_mode, transform->xsize_, data + pix_ix);
      data[pix_ix] = VP8LAddPixels(data[pix_ix], pred);
      ++pix_ix;

      // Subsequent columns.
      for (col = col_start + 1; col < col_end; ++col, ++pix_ix) {
        pred = PredictValue(pred_mode, transform->xsize_,
                            data + pix_ix);
        data[pix_ix] = VP8LAddPixels(data[pix_ix], pred);
      }
    }
  }
}

void VP8LAddGreenToBlueAndRed(const VP8LTransform* const transform,
                              argb_t* const data) {
  size_t i;
  const size_t num_pixs = transform->xsize_ * transform->ysize_;
  for (i = 0; i < num_pixs; ++i) {
    const argb_t argb = data[i];
    argb_t green = (argb >> 8) & 0xff;
    argb_t red_blue = argb & 0x00ff00ff;
    red_blue += ((green << 16) + green);
    red_blue &= 0x00ff00ff;
    data[i] = (argb & 0xff00ff00) + red_blue;
  }
}

typedef struct {
  int green_to_red_;
  int green_to_blue_;
  int red_to_blue_;
} ColorTransformElem;

static WEBP_INLINE argb_t ColorTransformDelta(signed char color_pred,
                                              signed char color) {
  return (argb_t)((int)(color_pred) * color) >> 5;
}

static WEBP_INLINE void ColorTransformElemInitFromCode(ColorTransformElem* elem,
                                                       argb_t color_pred) {
  elem->green_to_red_ = color_pred & 0xff;
  elem->green_to_blue_ = (color_pred >> 8) & 0xff;
  elem->red_to_blue_ = (color_pred >> 16) & 0xff;
}

static WEBP_INLINE argb_t TransformColor(const ColorTransformElem* const elem,
                                         argb_t argb) {
  const argb_t green = argb >> 8;
  const argb_t red = argb >> 16;
  argb_t new_red = red;
  argb_t new_blue = argb;

  new_red += ColorTransformDelta(elem->green_to_red_, green);
  new_red &= 0xff;
  new_blue += ColorTransformDelta(elem->green_to_blue_, green);
  new_blue += ColorTransformDelta(elem->red_to_blue_, new_red);
  new_blue &= 0xff;
  return (argb & 0xff00ff00) | (new_red << 16) | (new_blue);
}

void VP8LColorSpaceInverseTransform(const VP8LTransform* const transform,
                                    argb_t* const data) {
  size_t row, col, col_start, tile_offset;
  uint32_t pix_ix = 0;
  const uint32_t tile_size = 1 << transform->bits_;
  const uint32_t tiles_per_row = VP8LSubSampleSize(transform->xsize_,
                                                   transform->bits_);
  for (row = 0; row < transform->ysize_; ++row) {
    const uint32_t tile_base_ix = tiles_per_row * (row >> transform->bits_);
    for (tile_offset = 0, col_start = 0; tile_offset < tiles_per_row;
         ++tile_offset, col_start += tile_size) {
      ColorTransformElem color_pred;
      uint32_t col_end = col_start + tile_size;
      if (col_end > transform->xsize_) col_end = transform->xsize_;
      // Pick the appropriate color predictor mode (at start of every tile).
      ColorTransformElemInitFromCode(
          &color_pred, transform->data_[tile_base_ix + tile_offset]);
      for (col = col_start; col < col_end; ++col, ++pix_ix) {
        data[pix_ix] = TransformColor(&color_pred,
                                              data[pix_ix]);
      }
    }
  }
}

void VP8LColorIndexingInverseTransform(const VP8LTransform* const transform,
                                       argb_t* const data) {
  size_t i;
  const size_t num_pixs = transform->xsize_ * transform->ysize_;
  for (i = 0; i < num_pixs; ++i) {
    data[i] = transform->data_[(data[i] >> 8) & 0xff];
  }
}

int VP8LPixelBundleInverseTransform(const VP8LTransform* const transform,
                                    argb_t** const data) {
  uint32_t row, col, tile_x;
  const uint32_t bit_depth = 8 >> transform->bits_;
  const uint32_t num_cols = 1 << transform->bits_;
  const uint32_t bit_mask = (1 << bit_depth) - 1;
  const uint32_t xs = VP8LSubSampleSize(transform->xsize_, transform->bits_);

  uint32_t* tmp = (uint32_t*)calloc(
      transform->xsize_ * transform->ysize_, sizeof(*tmp));
  if (tmp == NULL) return 0;

  for (row = 0; row < transform->ysize_; ++row) {
    for (tile_x = 0; tile_x < xs; ++tile_x) {
      const argb_t* const rowp = (*data) + row * xs;
      uint32_t tile_code = (rowp[tile_x] >> 8) & 0xff;
      for (col = 0; col < num_cols; ++col) {
        const uint32_t x_all = tile_x * num_cols + col;
        if (x_all < transform->xsize_) {
          const uint32_t ix = row * transform->xsize_ + x_all;
          const uint32_t green = tile_code & bit_mask;
          tmp[ix] = ARGB_BLACK | (green << 8);
          tile_code >>= bit_depth;
        }
      }
    }
  }
  free(*data);
  *data = tmp;

  return 1;
}

static int IsAlphaMode(WEBP_CSP_MODE mode) {
  return (mode == MODE_RGBA || mode == MODE_BGRA || mode == MODE_ARGB);
}

// TODO: This function assumes that little-ending byte order is used.
// Need to add logic for big-endian.
int VP8LConvertColorSpaceFromBGRA(const uint8_t* const in_data,
                                  size_t num_pixels,
                                  WEBP_CSP_MODE out_colorspace,
                                  uint8_t** const output_data) {
  // RGBA_4444 & RGB_565 are unsupported for now & YUV modes are invalid.
  if (out_colorspace >= MODE_RGBA_4444) {
    return 0;
  } else {
    // Allocate as per mode.
    const int need_alpha = IsAlphaMode(out_colorspace);
    const int in_pixel_size = 4;
    size_t IN_BLUE  = 0;
    size_t IN_GREEN = 1;
    size_t IN_RED   = 2;
    size_t IN_ALPHA = 3;
    const size_t out_pixel_size = need_alpha ? 4 : 3;
    size_t OUT_BLUE;
    size_t OUT_GREEN;
    size_t OUT_RED;
    size_t OUT_ALPHA;
    size_t i;
    uint8_t* const output_data_lcl =
        (uint8_t*)malloc(num_pixels * out_pixel_size);
    if (output_data_lcl == NULL) return 0;
    switch (out_colorspace) {
      case MODE_RGB:
        OUT_RED   = 0;
        OUT_GREEN = 1;
        OUT_BLUE  = 2;
        break;
      case MODE_RGBA:
        OUT_RED   = 0;
        OUT_GREEN = 1;
        OUT_BLUE  = 2;
        OUT_ALPHA = 3;
        break;
      case MODE_BGR:
        OUT_BLUE  = 0;
        OUT_GREEN = 1;
        OUT_RED   = 2;
        break;
      case MODE_BGRA:
        OUT_BLUE  = 0;
        OUT_GREEN = 1;
        OUT_RED   = 2;
        OUT_ALPHA = 3;
        break;
      case MODE_ARGB:
        OUT_ALPHA = 0;
        OUT_RED   = 1;
        OUT_GREEN = 2;
        OUT_BLUE  = 3;
        break;
      default:
        // Code flow should not reach here.
        assert(0);
    }
    for (i = 0; i < num_pixels; ++i) {
      output_data_lcl[OUT_RED]   = in_data[IN_RED];
      output_data_lcl[OUT_GREEN] = in_data[IN_GREEN];
      output_data_lcl[OUT_BLUE]  = in_data[IN_BLUE];
      IN_RED    += in_pixel_size;
      IN_GREEN  += in_pixel_size;
      IN_BLUE   += in_pixel_size;
      OUT_RED   += out_pixel_size;
      OUT_GREEN += out_pixel_size;
      OUT_BLUE  += out_pixel_size;
      if(need_alpha) {
        output_data_lcl[OUT_ALPHA] = in_data[IN_ALPHA];
        IN_ALPHA  += in_pixel_size;
        OUT_ALPHA += out_pixel_size;
      }
    }
    *output_data = output_data_lcl;
    return 1;
  }
}

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif
