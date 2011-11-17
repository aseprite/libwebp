// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Alpha plane encoding and decoding library.
//
// Author: vikasa@google.com (Vikas Arora)

#include <string.h>  // for memcpy()
#include "./alpha.h"

#include "../utils/bit_reader.h"
#include "../utils/bit_writer.h"
#include "zlib.h"


#define MAX_SYMBOLS      255
#define ZLIB_CHUNK_SIZE  8192
#define ALPHA_HEADER_LEN 2

// -----------------------------------------------------------------------------
// Alpha Encode.

static int EncodeIdent(const uint8_t* data, int width, int height,
                       uint8_t** output, size_t* output_size) {
  const size_t data_size = height * width;
  uint8_t* alpha = NULL;
  assert((output != NULL) && (output_size != NULL));

  if (data == NULL) {
    return 0;
  }

  alpha = (uint8_t*)malloc(data_size);
  if (alpha == NULL) {
    return 0;
  }
  memcpy(alpha, data, data_size);
  *output_size = data_size;
  *output = alpha;
  return 1;
}

// -----------------------------------------------------------------------------

static int EncodeZlib(const uint8_t* data, int width, int height,
                      uint8_t** output, size_t* output_size, int algo) {
  unsigned char chunk[ZLIB_CHUNK_SIZE];
  const size_t data_size = height * width;
  int ret = Z_OK;
  z_stream strm;

  assert((data != NULL) && (output != NULL) && (output_size != NULL));

  *output = NULL;
  *output_size = 0;
  memset(&strm, 0, sizeof(strm));
  if (deflateInit(&strm, algo ? Z_BEST_SPEED : Z_BEST_COMPRESSION) != Z_OK) {
    return 0;
  }

  strm.next_in = (unsigned char*)data;
  strm.avail_in = data_size;
  do {
    size_t size_out;

    strm.next_out = chunk;
    strm.avail_out = ZLIB_CHUNK_SIZE;
    ret = deflate(&strm, Z_FINISH);
    if (ret == Z_STREAM_ERROR) {
      break;
    }
    size_out = ZLIB_CHUNK_SIZE - strm.avail_out;
    if (size_out) {
      size_t new_size = *output_size + size_out;
      uint8_t* new_output = (uint8_t*)realloc(*output, new_size);
      if (new_output == NULL) {
        ret = Z_MEM_ERROR;
        break;
      }
      memcpy(new_output + *output_size, chunk, size_out);
      *output_size = new_size;
      *output = new_output;
    }
  } while (ret != Z_STREAM_END || strm.avail_out == 0);

  deflateEnd(&strm);
  if (ret != Z_STREAM_END) {
    free(*output);
    output_size = 0;
    return 0;
  }
  return 1;
}

// -----------------------------------------------------------------------------

int EncodeAlpha(const uint8_t* data, int width, int height, int stride,
                int quality, int method,
                uint8_t** output, size_t* output_size) {
  const int kMaxImageDim = (1 << 14) - 1;
  uint8_t* compressed_alpha = NULL;
  uint8_t* quant_alpha = NULL;
  uint8_t* out = NULL;
  size_t compressed_size = 0;
  size_t data_size = height * width;
  float mse = 0.0;
  int ok = 0;
  int h;

  if ((data == NULL) || (output == NULL) || (output_size == NULL)) {
    return 0;
  }

  if (width <= 0 || width > kMaxImageDim ||
      height <= 0 || height > kMaxImageDim || stride < width) {
    return 0;
  }

  if (quality < 0 || quality > 100) {
    return 0;
  }

  if (method < 0 || method > 2) {
    return 0;
  }

  quant_alpha = (uint8_t*)malloc(data_size);
  if (quant_alpha == NULL) {
    return 0;
  }

  // Extract the alpha data (WidthXHeight) from raw_data (StrideXHeight).
  for (h = 0; h < height; ++h) {
    memcpy(quant_alpha + h * width, data + h * stride, width * sizeof(*data));
  }

  if (quality < 100) {  // No Quantization required for 'quality = 100'.
    // 16 Alpha levels gives quite a low MSE w.r.t Original Alpha plane hence
    // mapped to moderate quality 70. Hence Quality:[0, 70] -> Levels:[2, 16]
    // and Quality:]70, 100] -> Levels:]16, 256].
    const int alpha_levels = (quality <= 70) ?
                             2 + quality / 5 :
                             16 + (quality - 70) * 8;

    ok = QuantizeLevels(quant_alpha, width, height, alpha_levels, &mse);
    if (!ok) {
      free(quant_alpha);
      return 0;
    }
  }

  if (method == 0) {
    ok = EncodeIdent(quant_alpha, width, height,
                     &compressed_alpha, &compressed_size);
  } else if (method == 1) {
    ok = EncodeZlib(quant_alpha, width, height,
                    &compressed_alpha, &compressed_size, 0);
  }

  free(quant_alpha);
  if (!ok) {
    return 0;
  }

  out = (uint8_t*)malloc(compressed_size + ALPHA_HEADER_LEN);
  if (out == NULL) {
    free(compressed_alpha);
    return 0;
  } else {
    *output = out;
  }

  // Alpha bit-stream Header:
  // Byte0: Compression Method.
  // Byte1: Reserved for later extension.
  out[0] = method & 0xff;
  out[1] = 0;  // Reserved Byte.
  out += ALPHA_HEADER_LEN;
  memcpy(out, compressed_alpha, compressed_size);
  free(compressed_alpha);
  out += compressed_size;

  *output_size = out - *output;

  return 1;
}

// -----------------------------------------------------------------------------
// Alpha Decode.

static int DecodeIdent(const uint8_t* data, size_t data_size,
                       uint8_t* output) {
  assert((data != NULL) && (output != NULL));
  memcpy(output, data, data_size);
  return 1;
}

static int DecodeZlib(const uint8_t* data, size_t data_size,
                      uint8_t* output, size_t output_size) {
  z_stream strm;
  int ret = Z_OK;

  assert((data != NULL) && (output != NULL));

  memset(&strm, 0, sizeof(strm));
  if (inflateInit(&strm) != Z_OK) {
    return 0;
  }

  strm.avail_in = data_size;
  strm.next_in = (unsigned char*)data;
  do {
    strm.avail_out = output_size;
    strm.next_out = output;
    ret = inflate(&strm, Z_NO_FLUSH);
    if (ret == Z_NEED_DICT || ret == Z_DATA_ERROR || ret == Z_MEM_ERROR) {
      break;
    }
  } while (strm.avail_out == 0);

  inflateEnd(&strm);
  if (ret != Z_STREAM_END) {
    return 0;    // error
  }

  return 1;
}

// -----------------------------------------------------------------------------

int DecodeAlpha(const uint8_t* data, size_t data_size,
                int width, int height, int stride,
                uint8_t* output) {
  uint8_t* decoded_data = NULL;
  int ok = 0;
  int method;
  size_t decoded_size = height * width;

  if (data == NULL || output == NULL) {
    return 0;
  }

  if (data_size <= ALPHA_HEADER_LEN) {
    return 0;
  }

  if (width <= 0 || height <= 0 || stride < width) {
    return 0;
  }

  method = data[0];
  if (method > 2) {
    return 0;
  }

  decoded_data = (uint8_t*)malloc(decoded_size);
  if (decoded_data == NULL) {
    return 0;
  }

  data_size -= ALPHA_HEADER_LEN;
  data += ALPHA_HEADER_LEN;

  if (method == 0) {
    ok = DecodeIdent(data, data_size, decoded_data);
  } else if (method == 1) {
    ok = DecodeZlib(data, data_size, decoded_data, decoded_size);
  }
  if (ok) {
    // Construct raw_data (HeightXStride) from the alpha data (HeightXWidth).
    int h;
    for (h = 0; h < height; ++h) {
      memcpy(output + h * stride, decoded_data + h * width,
             width * sizeof(*data));
    }
  }
  free(decoded_data);

  return ok;
}
