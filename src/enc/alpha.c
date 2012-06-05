// Copyright 2011 Google Inc. All Rights Reserved.
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

#include "./vp8enci.h"
#include "../utils/filters.h"
#include "../utils/quant_levels.h"
#include "../webp/format_constants.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

// -----------------------------------------------------------------------------
// int EncodeAlpha(const uint8_t* data, int width, int height, int stride,
//                 int quality, int method, int filter, int effort_level,
//                 uint8_t** output, size_t* output_size)
//
// Encodes the given alpha data 'data' of size 'stride'x'height' via specified
// compression method 'method'. The pre-processing (Quantization) is
// performed if 'quality' is less than 100. For such cases, the encoding is
// lossy. Valid ranges for 'quality' is [0, 100] and 'method' is [0, 1]:
//   'method = 0' - No compression;
//   'method = 1' - Use lossless coder on the alpha plane only
// 'filter' values [0, 4] correspond to prediction modes none, horizontal,
// vertical & gradient filters. The prediction mode 4 will try all the
// prediction modes (0 to 3) and pick the best prediction mode.
// 'effort_level': specifies how much effort must be spent to try and reduce
//  the compressed output size. In range 0 (quick) to 6 (slow).
//
// 'output' corresponds to the buffer containing compressed alpha data.
//          This buffer is allocated by this method and caller should call
//          free(*output) when done.
// 'output_size' corresponds to size of this compressed alpha buffer.
//
// Returns 1 on successfully encoding the alpha and
//         0 if either:
//           invalid quality or method, or
//           memory allocation for the compressed data fails.

#ifdef USE_LOSSLESS_ENCODER

#include "../enc/vp8li.h"

static int EncodeLossless(const uint8_t* data, int width, int height,
                          int effort_level,  // in [0..6] range
                          VP8BitWriter* const bw) {

  int ok = 0;
  WebPConfig config;
  WebPPicture picture;
  VP8LBitWriter tmp_bw;

  WebPPictureInit(&picture);
  picture.width = width;
  picture.height = height;
  picture.use_argb_input = 1;
  if (!WebPPictureAlloc(&picture)) return 0;

  // Transfer the alpha values to the green channel.
  {
    int i, j;
    uint32_t* dst = picture.argb;
    const uint8_t* src = data;
    for (j = 0; j < picture.height; ++j) {
      for (i = 0; i < picture.width; ++i) {
        dst[i] = (src[i] << 8) | 0xff000000u;
      }
      src += width;
      dst += picture.argb_stride;
    }
  }

  WebPConfigInit(&config);
  config.lossless = 1;
  config.method = effort_level;  // impact is very small
  // quality below 50 doesn't change things much (in speed and size).
  // quality above 80 can be very very slow.
  config.quality = 40 + 10.f * effort_level;

  VP8LBitWriterInit(&tmp_bw, (width * height) >> 3);
  ok = (VP8LEncodeStream(&config, &picture, &tmp_bw) == VP8_ENC_OK);
  WebPPictureFree(&picture);
  if (ok) {
    const uint8_t* const data = VP8LBitWriterFinish(&tmp_bw);
    const size_t data_size = VP8LBitWriterNumBytes(&tmp_bw);
    VP8BitWriterAppend(bw, data, data_size);
  }
  VP8LBitWriterDestroy(&tmp_bw);
  return ok && !bw->error_;
}

#endif

// -----------------------------------------------------------------------------

static int EncodeAlphaInternal(const uint8_t* data, int width, int height,
                               int method, int filter, int reduce_levels,
                               int effort_level,  // in [0..6] range
                               uint8_t* tmp_alpha, VP8BitWriter* const bw) {
  int ok = 0;
  const uint8_t* alpha_src;
  WebPFilterFunc filter_func;
  uint8_t header;
  size_t expected_size;
  const size_t data_size = width * height;

#ifndef USE_LOSSLESS_ENCODER
  method = ALPHA_NO_COMPRESSION;
#endif
  assert(filter >= 0 && filter < WEBP_FILTER_LAST);
  assert(method >= ALPHA_NO_COMPRESSION);
  assert(method <= ALPHA_LOSSLESS_COMPRESSION);
  assert(sizeof(header) == ALPHA_HEADER_LEN);
  // TODO(skal): have a common function and #define's to validate alpha params.

  expected_size =
      (method == ALPHA_NO_COMPRESSION) ? (ALPHA_HEADER_LEN + data_size)
                                       : (data_size >> 5);
  header = method | (filter << 2);
  if (reduce_levels) header |= ALPHA_PREPROCESSED_LEVELS << 4;

  VP8BitWriterInit(bw, expected_size);
  VP8BitWriterAppend(bw, &header, ALPHA_HEADER_LEN);

  filter_func = WebPFilters[filter];
  if (filter_func) {
    filter_func(data, width, height, 1, width, tmp_alpha);
    alpha_src = tmp_alpha;
  }  else {
    alpha_src = data;
  }

  if (method == ALPHA_NO_COMPRESSION) {
    ok = VP8BitWriterAppend(bw, alpha_src, width * height);
    ok = ok && !bw->error_;
  } else {
#ifdef USE_LOSSLESS_ENCODER
    ok = EncodeLossless(alpha_src, width, height, effort_level, bw);
    VP8BitWriterFinish(bw);
#else
    (void)effort_level;
    assert(0);  // not reached.
#endif
  }
  return ok;
}

// -----------------------------------------------------------------------------

// TODO(skal): move to dsp/ ?
static void CopyPlane(const uint8_t* src, int src_stride,
                      uint8_t* dst, int dst_stride, int width, int height) {
  while (height-- > 0) {
    memcpy(dst, src, width);
    src += src_stride;
    dst += dst_stride;
  }
}

static int EncodeAlpha(const uint8_t* data, int width, int height, int stride,
                       int quality, int method, int filter,
                       int effort_level,
                       uint8_t** output, size_t* output_size) {
  uint8_t* quant_alpha = NULL;
  const size_t data_size = width * height;
  int ok = 1;
  const int reduce_levels = (quality < 100);

  // quick sanity checks
  assert(data != NULL && output != NULL && output_size != NULL);
  assert(width > 0 && height > 0);
  assert(stride >= width);
  assert(filter >= WEBP_FILTER_NONE && filter <= WEBP_FILTER_FAST);

  if (quality < 0 || quality > 100) {
    return 0;
  }

  if (method < ALPHA_NO_COMPRESSION || method > ALPHA_LOSSLESS_COMPRESSION) {
    return 0;
  }

  quant_alpha = (uint8_t*)malloc(data_size);
  if (quant_alpha == NULL) {
    return 0;
  }

  // Extract alpha data (width x height) from raw_data (stride x height).
  CopyPlane(data, stride, quant_alpha, width, width, height);

  if (reduce_levels) {  // No Quantization required for 'quality = 100'.
    // 16 alpha levels gives quite a low MSE w.r.t original alpha plane hence
    // mapped to moderate quality 70. Hence Quality:[0, 70] -> Levels:[2, 16]
    // and Quality:]70, 100] -> Levels:]16, 256].
    const int alpha_levels = (quality <= 70) ? (2 + quality / 5)
                                             : (16 + (quality - 70) * 8);
    ok = QuantizeLevels(quant_alpha, width, height, alpha_levels, NULL);
  }

  if (ok) {
    VP8BitWriter bw;
    size_t best_score;
    int test_filter;
    uint8_t* filtered_alpha = NULL;

    // We always test WEBP_FILTER_NONE first.
    ok = EncodeAlphaInternal(quant_alpha, width, height,
                             method, WEBP_FILTER_NONE, reduce_levels,
                             effort_level, NULL, &bw);
    if (!ok) {
      VP8BitWriterWipeOut(&bw);
      goto End;
    }
    best_score = VP8BitWriterSize(&bw);

    if (filter == WEBP_FILTER_FAST) {  // Quick estimate of a second candidate?
      filter = EstimateBestFilter(quant_alpha, width, height, width);
    }
    // Stop?
    if (filter == WEBP_FILTER_NONE) {
      goto Ok;
    }

    filtered_alpha = (uint8_t*)malloc(data_size);
    ok = (filtered_alpha != NULL);
    if (!ok) {
      goto End;
    }

    // Try the other mode(s).
    for (test_filter = WEBP_FILTER_HORIZONTAL;
         ok && (test_filter <= WEBP_FILTER_GRADIENT);
         ++test_filter) {
      VP8BitWriter tmp_bw;
      if (filter != WEBP_FILTER_BEST && test_filter != filter) {
        continue;
      }

      ok = EncodeAlphaInternal(quant_alpha, width, height,
                               method, test_filter, reduce_levels,
                               effort_level, filtered_alpha, &tmp_bw);
      if (ok) {
        const size_t score = VP8BitWriterSize(&tmp_bw);
        if (score < best_score) {
          // swap bitwriter objects.
          VP8BitWriter tmp = tmp_bw;
          tmp_bw = bw;
          bw = tmp;
          best_score = score;
        }
      } else {
        VP8BitWriterWipeOut(&bw);
      }
      VP8BitWriterWipeOut(&tmp_bw);
    }
 Ok:
    if (ok) {
      *output_size = VP8BitWriterSize(&bw);
      *output = VP8BitWriterBuf(&bw);
    }
    free(filtered_alpha);
  }
 End:
  free(quant_alpha);
  return ok;
}


//------------------------------------------------------------------------------
// Main calls

void VP8EncInitAlpha(VP8Encoder* enc) {
  enc->has_alpha_ = WebPPictureHasTransparency(enc->pic_);
  enc->alpha_data_ = NULL;
  enc->alpha_data_size_ = 0;
}

int VP8EncFinishAlpha(VP8Encoder* enc) {
  if (enc->has_alpha_) {
    const WebPConfig* config = enc->config_;
    const WebPPicture* pic = enc->pic_;
    uint8_t* tmp_data = NULL;
    size_t tmp_size = 0;
    const int effort_level = config->method;  // maps to [0..6]
    const WEBP_FILTER_TYPE filter =
        (config->alpha_filtering == 0) ? WEBP_FILTER_NONE :
        (config->alpha_filtering == 1) ? WEBP_FILTER_FAST :
                                         WEBP_FILTER_BEST;

    assert(pic->a);
    if (!EncodeAlpha(pic->a, pic->width, pic->height, pic->a_stride,
                     config->alpha_quality, config->alpha_compression,
                     filter, effort_level, &tmp_data, &tmp_size)) {
      return 0;
    }
    if (tmp_size != (uint32_t)tmp_size) {  // Sanity check.
      free(tmp_data);
      return 0;
    }
    enc->alpha_data_size_ = (uint32_t)tmp_size;
    enc->alpha_data_ = tmp_data;
  }
  return WebPReportProgress(enc->pic_, enc->percent_ + 20, &enc->percent_);
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
