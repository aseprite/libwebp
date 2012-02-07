// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: szabadka@google.com (Zoltan Szabadka)

#ifndef WEBP_DEC_IMAGE_TRANSFORM_H_
#define WEBP_DEC_IMAGE_TRANSFORM_H_

#include <stdio.h>
#include <stdint.h>

#include "../common/predictor.h"

typedef enum ImageTransformType {
  PREDICTOR_TRANSFORM = 0,
  CROSS_COLOR_TRANSFORM = 1,
  SUBTRACT_GREEN = 2,
  COLOR_INDEXING_TRANSFORM = 3,
  PIXEL_BUNDLE_TRANSFORM = 4,
} ImageTransformType;

struct ImageTransform {
  ImageTransformType type;
  int xsize;
  int ysize;
  void* data;
};

struct PerTileTransformData {
  int bits;
  uint32_t* image;
};

struct PixelBundleTransformData {
  int xbits;
};

void FreeImageTransformData(ImageTransform transform) {
  if (transform.type == PREDICTOR_TRANSFORM ||
      transform.type == CROSS_COLOR_TRANSFORM ||
      transform.type == COLOR_INDEXING_TRANSFORM) {
    PerTileTransformData* data = (PerTileTransformData*)transform.data;
    free(data->image);
  }
  free(transform.data);
}

void PixelBundleInverseTransform(int xsize, int ysize,
                                 PixelBundleTransformData data,
                                 uint32_t** image) {
  const int bit_depth = 1 << (3 - data.xbits);
  const int xs = (xsize + (1 << data.xbits) - 1) >> data.xbits;
  uint32_t* tmp_image = (uint32_t*)malloc(xs * ysize * sizeof(uint32_t));
  memcpy(tmp_image, *image, xs * ysize * sizeof(uint32_t));
  *image = (uint32_t*)realloc(*image, xsize * ysize * sizeof(uint32_t));
  for (int y = 0; y < ysize; ++y) {
    for (int tile_x = 0; tile_x < xs; ++tile_x) {
      uint32_t tile_code = (tmp_image[y * xs + tile_x] >> 8) & 0xff;
      for (int x = 0; x < 1 << data.xbits; ++x) {
        const int x_all = tile_x * (1 << data.xbits) + x;
        if (x_all >= xsize) continue;
        const int ix = y * xsize + x_all;
        const uint32_t g = tile_code & ((1 << bit_depth) - 1);
        (*image)[ix] = 0xff000000 | (g << 8);
        tile_code >>= bit_depth;
      }
    }
  }
  free(tmp_image);
}

void ColorIndexingInverseTransform(int xsize, int ysize,
                                   const uint32_t* palette, uint32_t* image) {
  for (int i = 0; i < xsize * ysize; ++i) {
    image[i] = palette[(image[i] >> 8) & 0xff];
  }
}

void ApplyInverseImageTransform(ImageTransform transform,
                                uint32_t** argb_image) {
  switch (transform.type) {
    case PREDICTOR_TRANSFORM:
      {
        PerTileTransformData* data = (PerTileTransformData*)transform.data;
        int img_size = transform.xsize * transform.ysize;
        uint32_t* tmp_image = (uint32_t*)malloc(img_size * sizeof(uint32_t));
        memcpy(tmp_image, *argb_image, img_size * sizeof(uint32_t));
        PredictorInverseTransform(transform.xsize, transform.ysize,
                                  data->bits, data->image,
                                  tmp_image, *argb_image);
        free(tmp_image);
      }
      break;
    case CROSS_COLOR_TRANSFORM:
      {
        PerTileTransformData* data = (PerTileTransformData*)transform.data;
        ColorSpaceInverseTransform(transform.xsize, transform.ysize,
                                   data->bits, data->image,
                                   *argb_image, *argb_image);
      }
      break;
    case SUBTRACT_GREEN:
      AddGreenToBlueAndRed(transform.xsize * transform.ysize, *argb_image);
      break;
    case COLOR_INDEXING_TRANSFORM:
      {
        PerTileTransformData* data = (PerTileTransformData*)transform.data;
        ColorIndexingInverseTransform(transform.xsize, transform.ysize,
                                      data->image, *argb_image);
      }
      break;
    case PIXEL_BUNDLE_TRANSFORM:
      {
        PixelBundleTransformData* data =
            (PixelBundleTransformData*)transform.data;
        PixelBundleInverseTransform(transform.xsize, transform.ysize,
                                    *data, argb_image);
      }
      break;
    default:
      // Nothing to do if we don't recognize the transform
      break;
  }
}

#endif  // WEBP_DEC_IMAGE_TRANSFORM_H_
