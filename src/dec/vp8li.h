// Copyright 2012 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// VP8 decoder: internal header.
//
// Author: Skal (pascal.massimino@gmail.com)
//         Vikas Arora(vikaas.arora@gmail.com)

#ifndef WEBP_DEC_VP8LI_H_
#define WEBP_DEC_VP8LI_H_

#include <string.h>     // for memcpy()
#include "../utils/bit_reader.h"
#include "../utils/color_cache.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

typedef unsigned int argb_t;

#define NUM_TRANSFORMS              16
#define HUFFMAN_CODES_PER_META_CODE  5

struct HuffmanTreeNode;

typedef enum ImageTransformType {
  PREDICTOR_TRANSFORM = 0,
  CROSS_COLOR_TRANSFORM = 1,
  SUBTRACT_GREEN = 2,
  COLOR_INDEXING_TRANSFORM = 3,
  PIXEL_BUNDLE_TRANSFORM = 4
} ImageTransformType;

typedef struct VP8LTransform VP8LTransform;
struct VP8LTransform {
  ImageTransformType    type_;    // transform type.
  int                   bits_;    // subsampling bits defining transform window.
  size_t                xsize_;   // transform window X index.
  size_t                ysize_;   // transform window Y index.
  uint32_t              *data_;   // transform data.
};

typedef struct VP8LDecoder VP8LDecoder;
struct VP8LDecoder {
  VP8StatusCode   status_;
  argb_t          *argb_;
  BitReader       *br_;

  uint32_t        xsize_;
  uint32_t        ysize_;
  uint32_t        x_ix_;
  uint32_t        y_ix_;
  uint32_t        pix_ix_;

  int             color_cache_size_;
  VP8LColorCache  *color_cache_;

  int             meta_index_;
  int             huffman_mask_;
  uint32_t        huffman_subsample_bits_;
  uint32_t        huffman_xsize_;
  uint32_t        *meta_codes_;
  uint32_t        *huffman_image_;
  struct HuffmanTreeNode *htrees_;
  struct HuffmanTreeNode *meta_htrees_[HUFFMAN_CODES_PER_META_CODE];

  int             next_transform_;
  VP8LTransform transforms_[NUM_TRANSFORMS];
};

//------------------------------------------------------------------------------
// internal functions. Not public.

// in vp8l.c

// Validates the VP8L data-header and retrieves basic header information viz
// width and height. Returns 0 in case of formatting error. *width/*height
// can be passed NULL.
int VP8LGetInfo(const uint8_t* data,
                uint32_t data_size,    // data available so far
                int *width, int *height);

// Decode a picture. Returns false in case of error.
int VP8LDecodeImage(VP8LDecoder* const dec, VP8Io* const io, uint32_t offset);

// Resets the decoder in its initial state, reclaiming memory.
void VP8LClear(VP8LDecoder* const dec);
//------------------------------------------------------------------------------

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif  /* WEBP_DEC_VP8LI_H_ */
