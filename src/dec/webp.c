// Copyright 2010 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Main decoding functions for WEBP images.
//
// Author: Skal (pascal.massimino@gmail.com)

#include <stdlib.h>

#include "./vp8i.h"
#include "./vp8li.h"
#include "./webpi.h"
#include "../mux/muxi.h"  // For MAX_CHUNK_PAYLOAD.
#include "../webp/mux.h"  // For 'ALPHA_FLAG'.

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

//------------------------------------------------------------------------------
// RIFF layout is:
//   Offset  tag
//   0...3   "RIFF" 4-byte tag
//   4...7   size of image data (including metadata) starting at offset 8
//   8...11  "WEBP"   our form-type signature
// The RIFF container (12 bytes) is followed by appropriate chunks:
//   12..15  "VP8 ": 4-bytes tags, signaling the use of VP8 video format
//   16..19  size of the raw VP8 image data, starting at offset 20
//   20....  the VP8 bytes
// Or,
//   12..15  "VP8L": 4-bytes tags, signaling the use of VP8L lossless format
//   16..19  size of the raw VP8L image data, starting at offset 20
//   20....  the VP8L bytes
// Or,
//   12..15  "VP8X": 4-bytes tags, describing the extended-VP8 chunk.
//   16..19  size of the VP8X chunk starting at offset 20.
//   20..23  VP8X flags bit-map corresponding to the chunk-types present.
//   24..27  Width of the Canvas Image.
//   28..31  Height of the Canvas Image.
// There can be extra chunks after the "VP8X" chunk (ICCP, TILE, FRM, VP8,
// META  ...)
// All 32-bits sizes are in little-endian order.
// Note: chunk data must be padded to multiple of 2 in size

static WEBP_INLINE uint32_t get_le32(const uint8_t* const data) {
  return data[0] | (data[1] << 8) | (data[2] << 16) | (data[3] << 24);
}

VP8StatusCode WebPParseRIFF(const uint8_t** data, uint32_t* data_size,
                            uint32_t* riff_size) {
  assert(data);
  assert(data_size);
  assert(riff_size);

  if (*data_size >= RIFF_HEADER_SIZE && !memcmp(*data, "RIFF", TAG_SIZE)) {
    if (memcmp(*data + 8, "WEBP", TAG_SIZE)) {
      return VP8_STATUS_BITSTREAM_ERROR;  // Wrong image file signature.
    } else {
      *riff_size = get_le32(*data + TAG_SIZE);
      // Check that we have at least one chunk (i.e "WEBP" + "VP8?nnnn").
      if (*riff_size < TAG_SIZE + CHUNK_HEADER_SIZE) {
        return VP8_STATUS_BITSTREAM_ERROR;
      }
      // We have a RIFF container. Skip it.
      *data += RIFF_HEADER_SIZE;
      *data_size -= RIFF_HEADER_SIZE;
    }
  } else {
    *riff_size = 0;  // Did not get full RIFF header.
  }
  return VP8_STATUS_OK;
}

VP8StatusCode WebPParseVP8X(const uint8_t** data, uint32_t* data_size,
                            uint32_t* bytes_skipped,
                            int* width, int* height, uint32_t* flags) {
  assert(data);
  assert(data_size);
  assert(bytes_skipped);

  *bytes_skipped = 0;

  if (*data_size < CHUNK_HEADER_SIZE + VP8X_CHUNK_SIZE) {
    return VP8_STATUS_NOT_ENOUGH_DATA;  // Insufficient data.
  }

  if (!memcmp(*data, "VP8X", TAG_SIZE)) {
    const uint32_t chunk_size = get_le32(*data + TAG_SIZE);
    if (chunk_size != VP8X_CHUNK_SIZE) {
      return VP8_STATUS_BITSTREAM_ERROR;  // Wrong chunk size.
    }
    if (flags != NULL) {
      *flags = get_le32(*data + 8);
    }
    if (width != NULL) {
      *width = get_le32(*data + 12);
    }
    if (height != NULL) {
      *height = get_le32(*data + 16);
    }
    // Skip over VP8X header bytes.
    *bytes_skipped = CHUNK_HEADER_SIZE + VP8X_CHUNK_SIZE;
    *data += *bytes_skipped;
    *data_size -= *bytes_skipped;
  }
  return VP8_STATUS_OK;
}

VP8StatusCode WebPParseOptionalChunks(const uint8_t** data, uint32_t* data_size,
                                      uint32_t riff_size,
                                      uint32_t* bytes_skipped,
                                      const uint8_t** alpha_data,
                                      uint32_t* alpha_size) {
  const uint8_t* buf;
  uint32_t buf_size;

  assert(data);
  assert(data_size);
  assert(bytes_skipped);
  assert(alpha_data);
  assert(alpha_size);

  buf = *data;
  buf_size = *data_size;
  *bytes_skipped = 0;
  *alpha_data = NULL;
  *alpha_size = 0;

  while (1) {
    uint32_t chunk_size;
    uint32_t cur_skip_size;
    const uint32_t bytes_skipped_header = TAG_SIZE +           // "WEBP".
                                          CHUNK_HEADER_SIZE +  // "VP8Xnnnn".
                                          VP8X_CHUNK_SIZE;     // data.
    *data = buf;
    *data_size = buf_size;

    if (buf_size < CHUNK_HEADER_SIZE) {  // Insufficient data.
      return VP8_STATUS_NOT_ENOUGH_DATA;
    }

    chunk_size = get_le32(buf + TAG_SIZE);
    // For odd-sized chunk-payload, there's one byte padding at the end.
    cur_skip_size = (CHUNK_HEADER_SIZE + chunk_size + 1) & ~1;

    // Check that total bytes skipped along with current chunk size
    // does not exceed riff_size.
    if (riff_size > 0 &&
        (bytes_skipped_header + *bytes_skipped + cur_skip_size > riff_size)) {
      return VP8_STATUS_BITSTREAM_ERROR;  // Not a valid chunk size.
    }

    if (buf_size < cur_skip_size) {  // Insufficient data.
      return VP8_STATUS_NOT_ENOUGH_DATA;
    }

    if (!memcmp(buf, "ALPH", TAG_SIZE)) {  // A valid ALPH header.
      *alpha_data = buf + CHUNK_HEADER_SIZE;
      *alpha_size = chunk_size;
    } else if (!memcmp(buf, "VP8 ", TAG_SIZE)) {  // A valid VP8 header.
      return VP8_STATUS_OK;  // Found.
    }

    // We have a full and valid chunk; skip it.
    buf += cur_skip_size;
    buf_size -= cur_skip_size;
    *bytes_skipped += cur_skip_size;
  }
}

VP8StatusCode WebPParseVP8Header(const uint8_t** data, uint32_t* data_size,
                                 uint32_t riff_size, uint32_t* bytes_skipped,
                                 uint32_t* chunk_size, int* is_lossless) {
  const int is_vp8 = !memcmp(*data, "VP8 ", TAG_SIZE);
  const int is_vp8l = !memcmp(*data, "VP8L", TAG_SIZE);
  assert(data);
  assert(data_size);
  assert(bytes_skipped);
  assert(chunk_size);

  *bytes_skipped = 0;
  *chunk_size = 0;

  if (*data_size < CHUNK_HEADER_SIZE) {
    return VP8_STATUS_NOT_ENOUGH_DATA;  // Insufficient data.
  }

  if (is_vp8 || is_vp8l) {
    if (is_lossless) *is_lossless = is_vp8l;
    *chunk_size = get_le32(*data + TAG_SIZE);
    if ((riff_size >= TAG_SIZE + CHUNK_HEADER_SIZE) &&  // "WEBP" + "VP8 nnnn".
        (*chunk_size > riff_size - (TAG_SIZE + CHUNK_HEADER_SIZE))) {
      return VP8_STATUS_BITSTREAM_ERROR;  // Inconsistent size information.
    }
    // Skip over CHUNK_HEADER_SIZE bytes from VP8 Header.
    *bytes_skipped = CHUNK_HEADER_SIZE;
    *data += *bytes_skipped;
    *data_size -= *bytes_skipped;
  }
  return VP8_STATUS_OK;
}

VP8StatusCode WebPParseHeaders(const uint8_t** data, uint32_t* data_size,
                               uint32_t* vp8_size, uint32_t* bytes_skipped,
                               const uint8_t** alpha_data,
                               uint32_t* alpha_size, int* is_lossless) {
  const uint8_t* buf;
  uint32_t buf_size;
  uint32_t riff_size;
  uint32_t vp8_size_tmp;
  uint32_t optional_data_size;
  uint32_t vp8x_skip_size;
  uint32_t vp8_skip_size;
  VP8StatusCode status;

  assert(data);
  assert(data_size);
  assert(vp8_size);
  assert(bytes_skipped);
  assert(alpha_data);
  assert(alpha_size);

  buf = *data;
  buf_size = *data_size;

  *vp8_size = 0;
  *bytes_skipped = 0;
  *alpha_data = NULL;
  *alpha_size = 0;

  if (buf == NULL || buf_size < RIFF_HEADER_SIZE) {
    return VP8_STATUS_NOT_ENOUGH_DATA;
  }

  // Skip over RIFF header.
  if (WebPParseRIFF(&buf, &buf_size, &riff_size) != VP8_STATUS_OK) {
    return VP8_STATUS_BITSTREAM_ERROR;  // Wrong RIFF header.
  }

  // Skip over VP8X header.
  status = WebPParseVP8X(&buf, &buf_size, &vp8x_skip_size, NULL, NULL, NULL);
  if (status != VP8_STATUS_OK) {
    return status;  // Wrong VP8X chunk / insufficient data.
  }
  if (vp8x_skip_size > 0) {
    // Skip over optional chunks.
    status = WebPParseOptionalChunks(&buf, &buf_size, riff_size,
                                     &optional_data_size,
                                     alpha_data, alpha_size);
    if (status != VP8_STATUS_OK) {
      return status;  // Found an invalid chunk size / insufficient data.
    }
  }

  // Skip over VP8 chunk header.
  status = WebPParseVP8Header(&buf, &buf_size, riff_size, &vp8_skip_size,
                              &vp8_size_tmp, is_lossless);
  if (status != VP8_STATUS_OK) {
    return status;  // Invalid VP8 header / insufficient data.
  }
  if (vp8_skip_size > 0) {
    *vp8_size = vp8_size_tmp;
  }

  *bytes_skipped = (uint32_t)(buf - *data);
  assert(buf - *data < MAX_CHUNK_PAYLOAD);
  assert(*bytes_skipped == *data_size - buf_size);
  *data = buf;
  *data_size = buf_size;
  return VP8_STATUS_OK;
}

//------------------------------------------------------------------------------
// WebPDecParams

void WebPResetDecParams(WebPDecParams* const params) {
  if (params) {
    memset(params, 0, sizeof(*params));
  }
}

//------------------------------------------------------------------------------
// "Into" decoding variants

// Main flow
static VP8StatusCode DecodeInto(const uint8_t* data, uint32_t data_size,
                                WebPDecParams* const params) {
  VP8Decoder* dec = VP8New();
  VP8StatusCode status = VP8_STATUS_OK;
  VP8Io io;

  assert(params);
  if (dec == NULL) {
    return VP8_STATUS_INVALID_PARAM;
  }

  VP8InitIo(&io);
  io.data = data;
  io.data_size = data_size;
  WebPInitCustomIo(params, &io);  // Plug the I/O functions.

#ifdef WEBP_USE_THREAD
  dec->use_threads_ = params->options && (params->options->use_threads > 0);
#else
  dec->use_threads_ = 0;
#endif

  // Decode bitstream header, update io->width/io->height.
  if (!VP8GetHeaders(dec, &io)) {
    status = VP8_STATUS_BITSTREAM_ERROR;
  } else {
    // Allocate/check output buffers.
    status = WebPAllocateDecBuffer(io.width, io.height, params->options,
                                   params->output);
    if (status == VP8_STATUS_OK) {
      // Decode
      if (!dec->is_lossless_) {
        if (!VP8Decode(dec, &io)) {
          status = dec->status_;
        }
      } else {
        VP8LDecoder* const vp8l_decoder = &dec->vp8l_decoder_;
        if (!VP8LDecodeImage(vp8l_decoder)) {
          status = VP8_STATUS_BITSTREAM_ERROR;
        }
      }
    }
  }
  VP8Delete(dec);
  if (status != VP8_STATUS_OK) {
    WebPFreeDecBuffer(params->output);
  }
  return status;
}

// Helpers
static uint8_t* DecodeIntoRGBABuffer(WEBP_CSP_MODE colorspace,
                                     const uint8_t* data, uint32_t data_size,
                                     uint8_t* rgba, int stride, int size) {
  WebPDecParams params;
  WebPDecBuffer buf;
  if (rgba == NULL) {
    return NULL;
  }
  WebPInitDecBuffer(&buf);
  WebPResetDecParams(&params);
  params.output = &buf;
  buf.colorspace    = colorspace;
  buf.u.RGBA.rgba   = rgba;
  buf.u.RGBA.stride = stride;
  buf.u.RGBA.size   = size;
  buf.is_external_memory = 1;
  if (DecodeInto(data, data_size, &params) != VP8_STATUS_OK) {
    return NULL;
  }
  return rgba;
}

uint8_t* WebPDecodeRGBInto(const uint8_t* data, uint32_t data_size,
                           uint8_t* output, int size, int stride) {
  return DecodeIntoRGBABuffer(MODE_RGB, data, data_size, output, stride, size);
}

uint8_t* WebPDecodeRGBAInto(const uint8_t* data, uint32_t data_size,
                            uint8_t* output, int size, int stride) {
  return DecodeIntoRGBABuffer(MODE_RGBA, data, data_size, output, stride, size);
}

uint8_t* WebPDecodeARGBInto(const uint8_t* data, uint32_t data_size,
                            uint8_t* output, int size, int stride) {
  return DecodeIntoRGBABuffer(MODE_ARGB, data, data_size, output, stride, size);
}

uint8_t* WebPDecodeBGRInto(const uint8_t* data, uint32_t data_size,
                           uint8_t* output, int size, int stride) {
  return DecodeIntoRGBABuffer(MODE_BGR, data, data_size, output, stride, size);
}

uint8_t* WebPDecodeBGRAInto(const uint8_t* data, uint32_t data_size,
                            uint8_t* output, int size, int stride) {
  return DecodeIntoRGBABuffer(MODE_BGRA, data, data_size, output, stride, size);
}

uint8_t* WebPDecodeYUVInto(const uint8_t* data, uint32_t data_size,
                           uint8_t* luma, int luma_size, int luma_stride,
                           uint8_t* u, int u_size, int u_stride,
                           uint8_t* v, int v_size, int v_stride) {
  WebPDecParams params;
  WebPDecBuffer output;
  if (luma == NULL) return NULL;
  WebPInitDecBuffer(&output);
  WebPResetDecParams(&params);
  params.output = &output;
  output.colorspace      = MODE_YUV;
  output.u.YUVA.y        = luma;
  output.u.YUVA.y_stride = luma_stride;
  output.u.YUVA.y_size   = luma_size;
  output.u.YUVA.u        = u;
  output.u.YUVA.u_stride = u_stride;
  output.u.YUVA.u_size   = u_size;
  output.u.YUVA.v        = v;
  output.u.YUVA.v_stride = v_stride;
  output.u.YUVA.v_size   = v_size;
  output.is_external_memory = 1;
  if (DecodeInto(data, data_size, &params) != VP8_STATUS_OK) {
    return NULL;
  }
  return luma;
}

//------------------------------------------------------------------------------

static uint8_t* Decode(WEBP_CSP_MODE mode, const uint8_t* data,
                       uint32_t data_size, int* width, int* height,
                       WebPDecBuffer* keep_info) {
  WebPDecParams params;
  WebPDecBuffer output;

  WebPInitDecBuffer(&output);
  WebPResetDecParams(&params);
  params.output = &output;
  output.colorspace = mode;

  // Retrieve (and report back) the required dimensions from bitstream.
  if (!WebPGetInfo(data, data_size, &output.width, &output.height)) {
    return NULL;
  }
  if (width != NULL) *width = output.width;
  if (height != NULL) *height = output.height;

  // Decode
  if (DecodeInto(data, data_size, &params) != VP8_STATUS_OK) {
    return NULL;
  }
  if (keep_info != NULL) {    // keep track of the side-info
    WebPCopyDecBuffer(&output, keep_info);
  }
  // return decoded samples (don't clear 'output'!)
  return (mode >= MODE_YUV) ? output.u.YUVA.y : output.u.RGBA.rgba;
}

uint8_t* WebPDecodeRGB(const uint8_t* data, uint32_t data_size,
                       int* width, int* height) {
  return Decode(MODE_RGB, data, data_size, width, height, NULL);
}

uint8_t* WebPDecodeRGBA(const uint8_t* data, uint32_t data_size,
                        int* width, int* height) {
  return Decode(MODE_RGBA, data, data_size, width, height, NULL);
}

uint8_t* WebPDecodeARGB(const uint8_t* data, uint32_t data_size,
                        int* width, int* height) {
  return Decode(MODE_ARGB, data, data_size, width, height, NULL);
}

uint8_t* WebPDecodeBGR(const uint8_t* data, uint32_t data_size,
                       int* width, int* height) {
  return Decode(MODE_BGR, data, data_size, width, height, NULL);
}

uint8_t* WebPDecodeBGRA(const uint8_t* data, uint32_t data_size,
                        int* width, int* height) {
  return Decode(MODE_BGRA, data, data_size, width, height, NULL);
}

uint8_t* WebPDecodeYUV(const uint8_t* data, uint32_t data_size,
                       int* width, int* height, uint8_t** u, uint8_t** v,
                       int* stride, int* uv_stride) {
  WebPDecBuffer output;   // only to preserve the side-infos
  uint8_t* const out = Decode(MODE_YUV, data, data_size,
                              width, height, &output);

  if (out != NULL) {
    const WebPYUVABuffer* const buf = &output.u.YUVA;
    *u = buf->u;
    *v = buf->v;
    *stride = buf->y_stride;
    *uv_stride = buf->u_stride;
    assert(buf->u_stride == buf->v_stride);
  }
  return out;
}

static void DefaultFeatures(WebPBitstreamFeatures* const features) {
  assert(features);
  memset(features, 0, sizeof(*features));
  features->bitstream_version = 0;
}

static VP8StatusCode GetFeatures(const uint8_t* data, uint32_t data_size,
                                 WebPBitstreamFeatures* const features) {
  uint32_t chunk_size = 0;
  uint32_t riff_size = 0;
  uint32_t flags = 0;
  uint32_t vp8x_skip_size = 0;
  uint32_t vp8_skip_size = 0;
  int is_lossless = 0;
  int* const width = &features->width;
  int* const height = &features->height;
  VP8StatusCode status;

  if (features == NULL || data == NULL) {
    return VP8_STATUS_INVALID_PARAM;
  }
  DefaultFeatures(features);

  // Skip over RIFF header.
  status = WebPParseRIFF(&data, &data_size, &riff_size);
  if (status != VP8_STATUS_OK) {
    return status;   // Wrong RIFF header / insufficient data.
  }

  // Skip over VP8X.
  status = WebPParseVP8X(&data, &data_size, &vp8x_skip_size,
                         width, height, &flags);
  if (status != VP8_STATUS_OK) {
    return status;  // Wrong VP8X / insufficient data.
  }
  features->has_alpha = !!(flags & ALPHA_FLAG);
  if (vp8x_skip_size > 0) {
    return VP8_STATUS_OK;  // Return features from VP8X header.
  }

  // Skip over VP8 header.
  status = WebPParseVP8Header(&data, &data_size, riff_size, &vp8_skip_size,
                              &chunk_size, &is_lossless);
  if (status != VP8_STATUS_OK) {
    return status;  // Wrong VP8 chunk-header / insufficient data.
  }
  if (vp8_skip_size == 0) {
    chunk_size = data_size;  // No VP8 chunk wrapper over raw VP8 data.
  }

  if (!is_lossless) {
    // Validates raw VP8 data.
    if (!VP8GetInfo(data, data_size, chunk_size, width, height)) {
      return VP8_STATUS_BITSTREAM_ERROR;
    }
  } else {
    // Validates raw VP8L data.
    if (!VP8LGetInfo(data, data_size, width, height)) {
      return VP8_STATUS_BITSTREAM_ERROR;
    }
    features->has_alpha = 1;
  }

  return VP8_STATUS_OK;  // Return features from VP8 header.
}

//------------------------------------------------------------------------------
// WebPGetInfo()

int WebPGetInfo(const uint8_t* data, uint32_t data_size,
                int* width, int* height) {
  WebPBitstreamFeatures features;

  if (GetFeatures(data, data_size, &features) != VP8_STATUS_OK) {
    return 0;
  }

  if (width != NULL) {
    *width  = features.width;
  }
  if (height != NULL) {
    *height = features.height;
  }

  return 1;
}

//------------------------------------------------------------------------------
// Advance decoding API

int WebPInitDecoderConfigInternal(WebPDecoderConfig* const config,
                                  int version) {
  if (version != WEBP_DECODER_ABI_VERSION) {
    return 0;   // version mismatch
  }
  if (config == NULL) {
    return 0;
  }
  memset(config, 0, sizeof(*config));
  DefaultFeatures(&config->input);
  WebPInitDecBuffer(&config->output);
  return 1;
}

VP8StatusCode WebPGetFeaturesInternal(const uint8_t* data, uint32_t data_size,
                                      WebPBitstreamFeatures* const features,
                                      int version) {
  VP8StatusCode status;
  if (version != WEBP_DECODER_ABI_VERSION) {
    return VP8_STATUS_INVALID_PARAM;   // version mismatch
  }
  if (features == NULL) {
    return VP8_STATUS_INVALID_PARAM;
  }

  status = GetFeatures(data, data_size, features);
  if (status == VP8_STATUS_NOT_ENOUGH_DATA) {
    return VP8_STATUS_BITSTREAM_ERROR;  // Not-enough-data treated as error.
  }
  return status;
}

VP8StatusCode WebPDecode(const uint8_t* data, uint32_t data_size,
                         WebPDecoderConfig* const config) {
  WebPDecParams params;
  VP8StatusCode status;

  if (config == NULL) {
    return VP8_STATUS_INVALID_PARAM;
  }

  status = GetFeatures(data, data_size, &config->input);
  if (status != VP8_STATUS_OK) {
    if (status == VP8_STATUS_NOT_ENOUGH_DATA) {
      return VP8_STATUS_BITSTREAM_ERROR;  // Not-enough-data treated as error.
    }
    return status;
  }

  WebPResetDecParams(&params);
  params.output = &config->output;
  params.options = &config->options;
  status = DecodeInto(data, data_size, &params);

  return status;
}

//------------------------------------------------------------------------------
// Cropping & rescaling.

int WebPIoInitFromOptions(const WebPDecoderOptions* const options,
                          VP8Io* const io) {
  const int W = io->width;
  const int H = io->height;
  int x = 0, y = 0, w = W, h = H;

  // Cropping
  io->use_cropping = (options != NULL) && (options->use_cropping > 0);
  if (io->use_cropping) {
    w = options->crop_width;
    h = options->crop_height;
    // TODO(skal): take colorspace into account. Don't assume YUV420.
    x = options->crop_left & ~1;
    y = options->crop_top & ~1;
    if (x < 0 || y < 0 || w <= 0 || h <= 0 || x + w > W || y + h > H) {
      return 0;  // out of frame boundary error
    }
  }
  io->crop_left   = x;
  io->crop_top    = y;
  io->crop_right  = x + w;
  io->crop_bottom = y + h;
  io->mb_w = w;
  io->mb_h = h;

  // Scaling
  io->use_scaling = (options != NULL) && (options->use_scaling > 0);
  if (io->use_scaling) {
    if (options->scaled_width <= 0 || options->scaled_height <= 0) {
      return 0;
    }
    io->scaled_width = options->scaled_width;
    io->scaled_height = options->scaled_height;
  }

  // Filter
  io->bypass_filtering = options && options->bypass_filtering;

  // Fancy upsampler
#ifdef FANCY_UPSAMPLING
  io->fancy_upsampling = (options == NULL) || (!options->no_fancy_upsampling);
#endif

  if (io->use_scaling) {
    // disable filter (only for large downscaling ratio).
    io->bypass_filtering = (io->scaled_width < W * 3 / 4) &&
                           (io->scaled_height < H * 3 / 4);
    io->fancy_upsampling = 0;
  }
  return 1;
}

#define RFIX 30
#define MULT_FIX(x,y) (((int64_t)(x) * (y) + (1 << (RFIX - 1))) >> RFIX)

void WebPRescalerInit(WebPRescaler* const wrk, int src_width, int src_height,
                      uint8_t* const dst, int dst_width, int dst_height,
                      int dst_stride, int num_channels, int x_add, int x_sub,
                      int y_add, int y_sub, int32_t* const work) {
  wrk->x_expand = (src_width < dst_width);
  wrk->src_width = src_width;
  wrk->src_height = src_height;
  wrk->dst_width = dst_width;
  wrk->dst_height = dst_height;
  wrk->dst = dst;
  wrk->dst_stride = dst_stride;
  wrk->num_channels = num_channels;
  // for 'x_expand', we use bilinear interpolation
  wrk->x_add = wrk->x_expand ? (x_sub - 1) : x_add - x_sub;
  wrk->x_sub = wrk->x_expand ? (x_add - 1) : x_sub;
  wrk->y_accum = y_add;
  wrk->y_add = y_add;
  wrk->y_sub = y_sub;
  wrk->fx_scale = (1 << RFIX) / x_sub;
  wrk->fy_scale = (1 << RFIX) / y_sub;
  wrk->fxy_scale = wrk->x_expand ?
      ((int64_t)dst_height << RFIX) / (x_sub * src_height) :
      ((int64_t)dst_height << RFIX) / (x_add * src_height);
  wrk->irow = work;
  wrk->frow = work + num_channels * dst_width;
}

void WebPRescalerImportRow(const uint8_t* const src, int channel,
                           WebPRescaler* const wrk) {
  const int x_stride = wrk->num_channels;
  const int x_out_max = wrk->dst_width * wrk->num_channels;
  int x_in = channel;
  int x_out;
  int accum = 0;
  if (!wrk->x_expand) {
    int sum = 0;
    for (x_out = channel; x_out < x_out_max; x_out += x_stride) {
      accum += wrk->x_add;
      for (; accum > 0; accum -= wrk->x_sub) {
        sum += src[x_in];
        x_in += x_stride;
      }
      {        // Emit next horizontal pixel.
        const int32_t base = src[x_in];
        const int32_t frac = base * (-accum);
        x_in += x_stride;
        wrk->frow[x_out] = (sum + base) * wrk->x_sub - frac;
        // fresh fractional start for next pixel
        sum = (int)MULT_FIX(frac, wrk->fx_scale);
      }
    }
  } else {        // simple bilinear interpolation
    int left = src[channel], right = src[channel];
    for (x_out = channel; x_out < x_out_max; x_out += x_stride) {
      if (accum < 0) {
        left = right;
        x_in += x_stride;
        right = src[x_in];
        accum += wrk->x_add;
      }
      wrk->frow[x_out] = right * wrk->x_add + (left - right) * accum;
      accum -= wrk->x_sub;
    }
  }
  // Accumulate the new row's contribution
  for (x_out = channel; x_out < x_out_max; x_out += x_stride) {
    wrk->irow[x_out] += wrk->frow[x_out];
  }
}

void WebPRescalerExportRow(WebPRescaler* const wrk) {
  int x_out;
  const int yscale = wrk->fy_scale * (-wrk->y_accum);
  const int x_out_max = wrk->dst_width * wrk->num_channels;
  assert(wrk->y_accum <= 0);
  for (x_out = 0; x_out < x_out_max; ++x_out) {
    const int frac = (int)MULT_FIX(wrk->frow[x_out], yscale);
    const int v = (int)MULT_FIX(wrk->irow[x_out] - frac, wrk->fxy_scale);
    wrk->dst[x_out] = (!(v & ~0xff)) ? v : (v < 0) ? 0 : 255;
    wrk->irow[x_out] = frac;   // new fractional start
  }
  wrk->y_accum += wrk->y_add;
  wrk->dst += wrk->dst_stride;
}

#undef MULT_FIX
#undef RFIX

//------------------------------------------------------------------------------

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif
