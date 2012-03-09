// Copyright 2012 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Image transforms and color space conversion methods for lossless decoder.
//
// Authors: Vikas Arora (vikaas.arora@gmail.com)
//          jyrki@google.com (Jyrki Alakuijala)
//          Urvang Joshi (urvang@google.com)

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#include <stdlib.h>
#include "./lossless.h"
#include "../dec/vp8li.h"

//------------------------------------------------------------------------------
// Inverse image transforms.

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
  if (a < NUM_LITERAL_CODES) {
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

// Inverse prediction.
static void PredictorInverseTransform(const VP8LTransform* const transform,
                                      size_t row_start, size_t row_end,
                                      argb_t* const data) {
  size_t row, col, col_start, tile_offset;
  uint32_t pix_ix = row_start * transform->xsize_;
  const uint32_t tile_size = 1 << transform->bits_;
  const uint32_t tiles_per_row = VP8LSubSampleSize(transform->xsize_,
                                                   transform->bits_);
  if (row_start == 0) {
    // First Row follows the L (mode=1) mode.
    data[0] = VP8LAddPixels(data[0], ARGB_BLACK);
    for (col = 1; col < transform->xsize_; ++col) {
      data[col] = VP8LAddPixels(data[col], data[col - 1]);
    }
    pix_ix += transform->xsize_;
    ++row_start;
  }

  for (row = row_start; row < row_end; ++row) {
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

// Add Green to Blue and Red channels (i.e. perform the inverse transform of
// 'Subtract Green').
static void AddGreenToBlueAndRed(const VP8LTransform* const transform,
                                 size_t row_start, size_t row_end,
                                 argb_t* const data) {
  size_t i;
  const size_t pix_start = row_start * transform->xsize_;
  const size_t pix_end = row_end * transform->xsize_;
  for (i = pix_start; i < pix_end; ++i) {
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

// Color space inverse transform.
static void ColorSpaceInverseTransform(const VP8LTransform* const transform,
                                       size_t row_start, size_t row_end,
                                       argb_t* const data) {
  size_t row, col, col_start, tile_offset;
  const uint32_t tile_size = 1 << transform->bits_;
  const uint32_t tiles_per_row = VP8LSubSampleSize(transform->xsize_,
                                                   transform->bits_);
  uint32_t pix_ix = row_start * transform->xsize_;

  for (row = row_start; row < row_end; ++row) {
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
        data[pix_ix] = TransformColor(&color_pred, data[pix_ix]);
      }
    }
  }
}

// Recover actual color values of pixels from their color indices.
static void ColorIndexingInverseTransform(const VP8LTransform* const transform,
                                          size_t row_start, size_t row_end,
                                          argb_t* const data) {
  size_t i;
  const size_t pix_start = row_start * transform->xsize_;
  const size_t pix_end = row_end * transform->xsize_;
  for (i = pix_start; i < pix_end; ++i) {
    data[i] = transform->data_[(data[i] >> 8) & 0xff];
  }
}

// Separate out pixels packed together using Pixel bundling.
// TODO: Move the allocation out of this function, and then make it generic,
// so that it can transform only the given rows.
static int PixelBundleInverseTransform(const VP8LTransform* const transform,
                                       argb_t** const data) {
  uint32_t row;
  const uint32_t bits_per_pixel = 8 >> transform->bits_;
  const uint32_t pixels_per_byte = 1 << transform->bits_;
  const uint32_t bit_mask = (1 << bits_per_pixel) - 1;
  const uint32_t num_packed_columns =
      VP8LSubSampleSize(transform->xsize_, transform->bits_);

  uint32_t* tmp = (uint32_t*)calloc(
      transform->xsize_ * transform->ysize_, sizeof(*tmp));
  if (tmp == NULL) return 0;

  for (row = 0; row < transform->ysize_; ++row) {
    const argb_t* const rowp = (*data) + row * num_packed_columns;
    const uint32_t pixels_done_before = row * transform->xsize_;

    // Unpack all packed pixels of this row except the last one.
    uint32_t col;
    for (col = 0; col < num_packed_columns - 1; ++col) {
      const uint32_t pixels_done_in_row = col * pixels_per_byte;
      const uint32_t pixels_done = pixels_done_before + pixels_done_in_row;
      uint32_t packed_pixel = (rowp[col] >> 8) & 0xff;
      // Unpack this 'packed_pixel'.
      uint32_t p;
      for (p = 0; p < pixels_per_byte; ++p) {
        const uint32_t green = packed_pixel & bit_mask;
        tmp[pixels_done + p] = ARGB_BLACK | (green << 8);
        packed_pixel >>= bits_per_pixel;
      }
    }

    // Last packed pixel is a special case, as it can be unpacked to anywhere
    // from 1 to (pixels_per_byte - 1) pixels.
    {
      const uint32_t last_col = num_packed_columns - 1;
      const uint32_t pixels_done_in_row = last_col * pixels_per_byte;
      const uint32_t pixels_done = pixels_done_before + pixels_done_in_row;
      uint32_t packed_pixel = (rowp[last_col] >> 8) & 0xff;
      // Unpack this 'packed_pixel'.
      uint32_t p;
      for (p = 0; p < pixels_per_byte; ++p) {
        if (pixels_done_in_row + p < transform->xsize_) {
          const uint32_t green = packed_pixel & bit_mask;
          tmp[pixels_done + p] = ARGB_BLACK | (green << 8);
          packed_pixel >>= bits_per_pixel;
        } else {
          break;
        }
      }
    }
  }
  free(*data);
  *data = tmp;
  return 1;
}

VP8StatusCode VP8LInverseTransform(const VP8LTransform* const transform,
                                   size_t row_start, size_t row_end,
                                   argb_t** const data) {
  VP8StatusCode status = VP8_STATUS_OK;
  assert(row_start < row_end);
  assert(row_end <= transform->ysize_);
  switch (transform->type_) {
    case SUBTRACT_GREEN:
      AddGreenToBlueAndRed(transform, row_start, row_end, *data);
      break;
    case PREDICTOR_TRANSFORM:
      PredictorInverseTransform(transform, row_start, row_end, *data);
      break;
    case CROSS_COLOR_TRANSFORM:
      ColorSpaceInverseTransform(transform, row_start, row_end, *data);
      break;
    case COLOR_INDEXING_TRANSFORM:
      ColorIndexingInverseTransform(transform, row_start, row_end, *data);
      break;
    case PIXEL_BUNDLE_TRANSFORM:
      if (!PixelBundleInverseTransform(transform, data)) {
        status = VP8_STATUS_OUT_OF_MEMORY;
      }
      break;
    default:
      status = VP8_STATUS_BITSTREAM_ERROR;
      break;
  }
  return status;
}

//------------------------------------------------------------------------------
// Color space conversion.

static void ConvertBGRAToRGB(const argb_t* const in_data, size_t num_pixels,
                             uint8_t* const out_data) {
  const argb_t* in_wordp = in_data;
  uint8_t* out_bytep = out_data;
  size_t i;
  for (i = 0; i < num_pixels; ++i) {
    const argb_t argb = in_wordp[i];
    *out_bytep++ = (argb >> 16) & 0xff;
    *out_bytep++ = (argb >> 8) & 0xff;
    *out_bytep++ = (argb >> 0) & 0xff;
  }
}

static void ConvertBGRAToRGBA(const argb_t* const in_data, size_t num_pixels,
                              uint8_t* const out_data) {
  const argb_t* in_wordp = in_data;
  uint8_t* out_bytep = out_data;
  size_t i;
  for (i = 0; i < num_pixels; ++i) {
    const argb_t argb = in_wordp[i];
    *out_bytep++ = (argb >> 16) & 0xff;
    *out_bytep++ = (argb >> 8) & 0xff;
    *out_bytep++ = (argb >> 0) & 0xff;
    *out_bytep++ = (argb >> 24) & 0xff;
  }
}

static void ConvertBGRAToBGR(const argb_t* const in_data, size_t num_pixels,
                             uint8_t* const out_data) {
  const argb_t* in_wordp = in_data;
  uint8_t* out_bytep = out_data;
  size_t i;
  for (i = 0; i < num_pixels; ++i) {
    const argb_t argb = in_wordp[i];
    *out_bytep++ = (argb >> 0) & 0xff;
    *out_bytep++ = (argb >> 8) & 0xff;
    *out_bytep++ = (argb >> 16) & 0xff;
  }
}

static void ConvertBGRAToARGB(const argb_t* const in_data, size_t num_pixels,
                              uint8_t* const out_data) {
  const argb_t* in_wordp = in_data;
  uint8_t* out_bytep = out_data;
  size_t i;
  for (i = 0; i < num_pixels; ++i) {
    const argb_t argb = in_wordp[i];
    *out_bytep++ = (argb >> 24) & 0xff;
    *out_bytep++ = (argb >> 16) & 0xff;
    *out_bytep++ = (argb >> 8) & 0xff;
    *out_bytep++ = (argb >> 0) & 0xff;
  }
}

void VP8LConvertFromBGRA(const argb_t* const in_data, size_t num_pixels,
                        WEBP_CSP_MODE out_colorspace,
                        uint8_t* const rgba) {
  switch (out_colorspace) {
    case MODE_RGB:
      ConvertBGRAToRGB(in_data, num_pixels, rgba);
      break;
    case MODE_RGBA:
      ConvertBGRAToRGBA(in_data, num_pixels, rgba);
      break;
    case MODE_BGR:
      ConvertBGRAToBGR(in_data, num_pixels, rgba);
      break;
    case MODE_BGRA:
      memcpy(rgba, in_data, num_pixels * 4);
      break;
    case MODE_ARGB:
      ConvertBGRAToARGB(in_data, num_pixels, rgba);
      break;
    default:
      // Code flow should not reach here.
      assert(0);
  }
}

//------------------------------------------------------------------------------

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif
