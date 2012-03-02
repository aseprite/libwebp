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

#define NUM_TRANSFORMS               8
#define HUFFMAN_CODES_PER_META_CODE  5

struct HuffmanTree;

typedef enum VP8LDecodeState {
  READ_DATA = 0,
  READ_HDR = 1,
  READ_DIM = 2
} VP8LDecodeState;

typedef enum VP8LImageTransformType {
  PREDICTOR_TRANSFORM = 0,
  CROSS_COLOR_TRANSFORM = 1,
  SUBTRACT_GREEN = 2,
  COLOR_INDEXING_TRANSFORM = 3,
  PIXEL_BUNDLE_TRANSFORM = 4
} VP8LImageTransformType;

typedef struct VP8LTransform VP8LTransform;
struct VP8LTransform {
  VP8LImageTransformType type_;   // transform type.
  int                    bits_;   // subsampling bits defining transform window.
  size_t                 xsize_;  // transform window X index.
  size_t                 ysize_;  // transform window Y index.
  uint32_t               *data_;  // transform data.
};

typedef struct VP8LMetadata VP8LMetadata;
struct VP8LMetadata {
  int             color_cache_size_;
  VP8LColorCache  *color_cache_;

  int             meta_index_;
  int             num_huffman_trees_;
  int             huffman_mask_;
  int             huffman_subsample_bits_;
  int             huffman_xsize_;
  uint32_t        *meta_codes_;
  uint32_t        *huffman_image_;
  struct HuffmanTree *htrees_;
  struct HuffmanTree *meta_htrees_[HUFFMAN_CODES_PER_META_CODE];
};

typedef struct VP8LDecoder VP8LDecoder;
struct VP8LDecoder {
  VP8StatusCode    status_;
  VP8LDecodeState  action_;
  VP8LDecodeState  state_;

  argb_t           *argb_;  // Internal data: always in BGRA color mode.
  BitReader        br_;
  BitReader        br_check_point_;

  WEBP_CSP_MODE    output_colorspace_;
  uint8_t*         decoded_data_;  // Output in 'output_colorspace_' color mode.

  int              width_;
  int              height_;

  int              xsize_;
  int              ysize_;
  int              row_;
  int              pix_;
  int              level_;

  VP8LMetadata     hdr_;

  int              next_transform_;
  VP8LTransform    transforms_[NUM_TRANSFORMS];
};

//------------------------------------------------------------------------------
// internal functions. Not public.

// in vp8l.c

// Validates the VP8L data-header and retrieves basic header information viz
// width and height. Returns 0 in case of formatting error. width/height
// can be passed NULL.
int VP8LGetInfo(const uint8_t* data,
                uint32_t data_size,    // data available so far
                int *width, int *height);

// Initializes the decoder object with the given output colorspace.
void VP8LInitDecoder(VP8LDecoder* const dec);

// Decodes the image header. Returns false in case of error.
int VP8LDecodeHeader(VP8LDecoder* const dec, VP8Io* const io, uint32_t offset);

// Decodes a image pixels. It's required to decode the lossless header before
// calling this funcion. Returns false in case of error.
int VP8LDecodeImage(VP8LDecoder* const dec);

// Decode image pixels.
int VP8LDecodePixels(VP8LDecoder* const dec, argb_t* const decoded_data);

// Convert decoded lossless stream to the specified ColorSpace.
int VP8LConvertColorSpaceFromBGRA(uint8_t* const in_data, size_t num_pixels,
                                  WEBP_CSP_MODE out_colorspace,
                                  uint8_t** const output_data);

// Resets the decoder in its initial state, reclaiming memory.
void VP8LClear(VP8LDecoder* const dec);
//------------------------------------------------------------------------------

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif  /* WEBP_DEC_VP8LI_H_ */
