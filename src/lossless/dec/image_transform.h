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

#include "../common/integral_types.h"
#include "../common/predictor.h"

typedef enum ImageTransformType {
  PREDICTOR_TRANSFORM = 0,
  CROSS_COLOR_TRANSFORM = 1,
  SUBTRACT_GREEN = 2,
  COMPONENT_SUBSAMPLING_TRANSFORM = 3,
  COLOR_INDEXING_TRANSFORM = 5,
  PIXEL_BUNDLE_TRANSFORM = 6,
} ImageTransformType;

struct ImageTransform {
  ImageTransformType type;
  int xsize;
  int ysize;
  void* data;
};

struct PerTileTransformData {
  int bits;
  uint32* image;
};

struct PixelBundleTransformData {
  int xbits;
  int ybits;
  int bit_depth;
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
                                 uint32** image) {
  int xs = (xsize + (1 << data.xbits) - 1) >> data.xbits;
  int ys = (ysize + (1 << data.ybits) - 1) >> data.ybits;
  uint32* tmp_image = (uint32*)malloc(xs * ys * sizeof(uint32));
  memcpy(tmp_image, *image, xs * ys * sizeof(uint32));
  *image = (uint32*)realloc(*image, xsize * ysize * sizeof(uint32));
  for (int tile_y = 0; tile_y < ys; ++tile_y) {
    for (int tile_x = 0; tile_x < xs; ++tile_x) {
      uint32 tile_code = tmp_image[tile_y * xs + tile_x];
      for (int y = 0; y < 1 << data.ybits; ++y) {
        int y_all = tile_y * (1 << data.ybits) + y;
        if (y_all >= ysize) continue;
        for (int x = 0; x < 1 << data.xbits; ++x) {
          int x_all = tile_x * (1 << data.xbits) + x;
          if (x_all >= xsize) continue;
          int ix = y_all * xsize + x_all;
          int bit_position = (y * (1 << data.xbits) + x) * data.bit_depth;
          if (bit_position < 16) {
            bit_position += 8;
          } else if (bit_position < 24) {
            bit_position -= 16;
          }
          uint32 g = (tile_code >> bit_position) & ((1 << data.bit_depth) - 1);
          (*image)[ix] = 0xff000000 | (g << 8);
        }
      }
    }
  }
  free(tmp_image);
}

// Compute a lookup table for unpacking packed component dynamics.
//
// bits indicates how many bits are used to store the full range [0..255].
// for the case of bits == 4, the values stored are in the range of [0..15],
// and 15 is expanded to 255, i.e., 14 to 238, i.e., the highest bits are
// used to fill a range of lower bits.
static void ConstructSubsamplingLut(int bits, uint8 *lut) {
  int sbits = 8 - bits;
  for (uint32 i = 0; i < 256; ++i) {
    uint32 c = i << bits;
    for (int nbits = sbits; nbits < 8; nbits += sbits) {
      int shift = 8 - nbits - sbits;
      if (shift >= 0) {
        c |= i << shift;
      } else {
        c |= i >> -shift;
      }
    }
    lut[i] = c;
  }
}

void ComponentSubsamplingInverseTransform(int xsize, int ysize,
                                          const int* bits, uint32* image) {
  uint8 lut_a[256];
  uint8 lut_r[256];
  uint8 lut_g[256];
  uint8 lut_b[256];
  ConstructSubsamplingLut(bits[3], lut_a);
  ConstructSubsamplingLut(bits[2], lut_r);
  ConstructSubsamplingLut(bits[1], lut_g);
  ConstructSubsamplingLut(bits[0], lut_b);
  for (int i = 0; i < xsize * ysize; ++i) {
    uint32 a = lut_a[(image[i] >> 24)] << 24;
    uint32 r = lut_r[(image[i] >> 16) & 0xff] << 16;
    uint32 g = lut_g[(image[i] >> 8) & 0xff] << 8;
    uint32 b = lut_b[image[i] & 0xff];
    image[i] = a | r | g | b;
  }
}

void ColorIndexingInverseTransform(int xsize, int ysize,
                                   const uint32* palette, uint32* image) {
  for (int i = 0; i < xsize * ysize; ++i) {
    image[i] = palette[(image[i] >> 8) & 0xff];
  }
}

void ApplyInverseImageTransform(ImageTransform transform, uint32** argb_image) {
  switch (transform.type) {
    case PREDICTOR_TRANSFORM:
      {
        PerTileTransformData* data = (PerTileTransformData*)transform.data;
        int img_size = transform.xsize * transform.ysize;
        uint32* tmp_image = (uint32*)malloc(img_size * sizeof(uint32));
        memcpy(tmp_image, *argb_image, img_size * sizeof(uint32));
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
    case COMPONENT_SUBSAMPLING_TRANSFORM:
      {
        int* bits = (int*)transform.data;
        ComponentSubsamplingInverseTransform(transform.xsize, transform.ysize,
                                             bits, *argb_image);
      }
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
