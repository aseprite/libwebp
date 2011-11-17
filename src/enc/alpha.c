// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Alpha-plane compression.
//
// Author: Skal (pascal.massimino@gmail.com)

#include <assert.h>
#include <stdlib.h>
#include "vp8enci.h"

#ifdef WEBP_EXPERIMENTAL_FEATURES
#include "zlib.h"
#endif

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#ifdef WEBP_EXPERIMENTAL_FEATURES

//------------------------------------------------------------------------------

static int CompressAlpha(const uint8_t* data, size_t data_size,
                         uint8_t** output, size_t* output_size,
                         int algo) {
  (void)data;
  (void)data_size;
  (void)output;
  (void)output_size;
  (void)algo;
  // Don't compress Alpha-plane, as it's compressed outside WebPEncode.
  // TODO: Link to Alpha compression lib to compress the alpha data.
  *output = NULL;
  *output_size = 0;
  return 1;
}

#endif    /* WEBP_EXPERIMENTAL_FEATURES */

void VP8EncInitAlpha(VP8Encoder* enc) {
  enc->has_alpha_ = (enc->pic_->a != NULL);
  enc->alpha_data_ = NULL;
  enc->alpha_data_size_ = 0;
}

void VP8EncCodeAlphaBlock(VP8EncIterator* it) {
  (void)it;
  // Nothing for now. We just ZLIB-compress in the end.
}

int VP8EncFinishAlpha(VP8Encoder* enc) {
  if (enc->has_alpha_) {
#ifdef WEBP_EXPERIMENTAL_FEATURES
    const WebPPicture* pic = enc->pic_;
    assert(pic->a);
    if (!CompressAlpha(pic->a, pic->width * pic->height,
                       &enc->alpha_data_, &enc->alpha_data_size_,
                       enc->config_->alpha_compression)) {
      return 0;
    }
#endif
  }
  return 1;
}

void VP8EncDeleteAlpha(VP8Encoder* enc) {
  free(enc->alpha_data_);
  enc->alpha_data_ = NULL;
  enc->alpha_data_size_ = 0;
  enc->has_alpha_ = 0;
}

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif
