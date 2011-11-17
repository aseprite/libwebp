// Copyright 2010 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
//  Command-line tool for converting a WebP image to PNG image.
//  This binary is WebP-Mux bit-stream aware and extracts Alpha and converts it
//  to PNG appropriately.
//
//  Compile with: gcc -o webp2png webp2png.c -lexperimental -lwebpmux -lwebp
//
// Author: Vikas Arora (vikasa@google.com)

#ifdef WEBP_EXPERIMENTAL_FEATURES

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef WEBP_HAVE_PNG
#include <png.h>
#endif

#ifdef HAVE_WINCODEC_H
#ifdef __MINGW32__
#define INITGUID  // Without this GUIDs are declared extern and fail to link
#endif
#define CINTERFACE
#define COBJMACROS
#define _WIN32_IE 0x500  // Workaround bug in shlwapi.h when compiling C++
                         // code with COBJMACROS.
#include <shlwapi.h>
#include <windows.h>
#include <wincodec.h>
#endif

#include "experimental/alpha.h"
#include "webp/decode.h"
#include "webp/mux.h"
#include "stopwatch.h"

static int verbose = 0;

//------------------------------------------------------------------------------

#ifdef HAVE_WINCODEC_H

#define IFS(fn)                \
  do {                         \
     if (SUCCEEDED(hr))        \
     {                         \
        hr = (fn);             \
        if (FAILED(hr) && verbose)           \
          printf(#fn " failed %08x\n", hr);  \
     }                         \
  } while (0)

#ifdef __cplusplus
#define MAKE_REFGUID(x) (x)
#else
#define MAKE_REFGUID(x) &(x)
#endif

static HRESULT CreateOutputStream(const char* out_file_name,
                                  IStream** ppStream) {
  HRESULT hr = S_OK;
  IFS(SHCreateStreamOnFileA(out_file_name, STGM_WRITE | STGM_CREATE, ppStream));
  if (FAILED(hr))
    printf("Error opening output file %s (%08x)\n", out_file_name, hr);
  return hr;
}

static HRESULT WriteUsingWIC(const char* out_file_name, REFGUID container_guid,
                             unsigned char* rgb, int stride,
                             uint32_t width, uint32_t height, int has_alpha) {
  HRESULT hr = S_OK;
  IWICImagingFactory* pFactory = NULL;
  IWICBitmapFrameEncode* pFrame = NULL;
  IWICBitmapEncoder* pEncoder = NULL;
  IStream* pStream = NULL;
  WICPixelFormatGUID pixel_format = has_alpha ? GUID_WICPixelFormat32bppBGRA
                                              : GUID_WICPixelFormat24bppBGR;

  IFS(CoInitialize(NULL));
  IFS(CoCreateInstance(MAKE_REFGUID(CLSID_WICImagingFactory), NULL,
          CLSCTX_INPROC_SERVER, MAKE_REFGUID(IID_IWICImagingFactory),
          (LPVOID*)&pFactory));
  if (hr == REGDB_E_CLASSNOTREG) {
    printf("Couldn't access Windows Imaging Component (are you running \n");
    printf("Windows XP SP3 or newer?). PNG support not available.\n");
    printf("Use -ppm or -pgm for available PPM and PGM formats.\n");
  }
  IFS(CreateOutputStream(out_file_name, &pStream));
  IFS(IWICImagingFactory_CreateEncoder(pFactory, container_guid, NULL,
          &pEncoder));
  IFS(IWICBitmapEncoder_Initialize(pEncoder, pStream,
                                   WICBitmapEncoderNoCache));
  IFS(IWICBitmapEncoder_CreateNewFrame(pEncoder, &pFrame, NULL));
  IFS(IWICBitmapFrameEncode_Initialize(pFrame, NULL));
  IFS(IWICBitmapFrameEncode_SetSize(pFrame, width, height));
  IFS(IWICBitmapFrameEncode_SetPixelFormat(pFrame, &pixel_format));
  IFS(IWICBitmapFrameEncode_WritePixels(pFrame, height, stride,
          height * stride, rgb));
  IFS(IWICBitmapFrameEncode_Commit(pFrame));
  IFS(IWICBitmapEncoder_Commit(pEncoder));

  if (pFrame != NULL) IUnknown_Release(pFrame);
  if (pEncoder != NULL) IUnknown_Release(pEncoder);
  if (pFactory != NULL) IUnknown_Release(pFactory);
  if (pStream != NULL) IUnknown_Release(pStream);
  return hr;
}

static int WritePNG(const char* out_file_name,
                    const WebPDecBuffer* const buffer) {
  const uint32_t width = buffer->width;
  const uint32_t height = buffer->height;
  unsigned char* const rgb = buffer->u.RGBA.rgba;
  const int stride = buffer->u.RGBA.stride;
  const int has_alpha = (buffer->colorspace == MODE_BGRA);

  return SUCCEEDED(WriteUsingWIC(out_file_name,
             MAKE_REFGUID(GUID_ContainerFormatPng), rgb, stride, width,
             height, has_alpha));
}

#elif defined(WEBP_HAVE_PNG)    // !HAVE_WINCODEC_H
static void PNGAPI error_function(png_structp png, png_const_charp dummy) {
  (void)dummy;  // remove variable-unused warning
  longjmp(png_jmpbuf(png), 1);
}

static int WritePNG(FILE* out_file, const WebPDecBuffer* const buffer) {
  const uint32_t width = buffer->width;
  const uint32_t height = buffer->height;
  unsigned char* const rgba = buffer->u.RGBA.rgba;
  const int stride = buffer->u.RGBA.stride;
  const int has_alpha = (buffer->colorspace == MODE_RGBA);
  png_structp png;
  png_infop info;
  png_uint_32 y;

  png = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                NULL, error_function, NULL);
  if (png == NULL) {
    return 0;
  }
  info = png_create_info_struct(png);
  if (info == NULL) {
    png_destroy_write_struct(&png, NULL);
    return 0;
  }
  if (setjmp(png_jmpbuf(png))) {
    png_destroy_write_struct(&png, &info);
    return 0;
  }
  png_init_io(png, out_file);
  png_set_IHDR(png, info, width, height, 8,
               has_alpha ? PNG_COLOR_TYPE_RGBA : PNG_COLOR_TYPE_RGB,
               PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
               PNG_FILTER_TYPE_DEFAULT);
  png_write_info(png, info);
  for (y = 0; y < height; ++y) {
    png_bytep row = rgba + y * stride;
    png_write_rows(png, &row, 1);
  }
  png_write_end(png, info);
  png_destroy_write_struct(&png, &info);
  return 1;
}
#else    // !HAVE_WINCODEC_H && !WEBP_HAVE_PNG

typedef uint32_t png_uint_32;

static int WritePNG(FILE* out_file, const WebPDecBuffer* const buffer) {
  (void)out_file;
  (void)buffer;
  printf("PNG support not compiled. Please install the libpng development "
         "package before building.\n");
  return 0;
}
#endif

static void SaveOutput(const WebPDecBuffer* const buffer,
                       const uint8_t* alpha_data, const char* const out_file) {
  FILE* fout = NULL;
  Stopwatch stop_watch;
  const uint32_t width = buffer->width;
  const uint32_t height = buffer->height;
  const int stride = buffer->u.RGBA.stride;
  unsigned char* const rgb = buffer->u.RGBA.rgba;
  int needs_open_file = 1;
  int ok = 1;
#ifdef HAVE_WINCODEC_H
  const int has_alpha =
      (buffer->colorspace == MODE_BGRA) && (alpha_data != NULL);
#else
  const int has_alpha =
      (buffer->colorspace == MODE_RGBA) && (alpha_data != NULL);
#endif

  if (has_alpha) {
    const uint32_t alpha_stride = width * sizeof(*rgb);
    uint32_t y, x;

    // Add Alpha to decoded RGB data.
    for (y = 0; y < height; ++y) {
      uint8_t* rgb_row = rgb + y * stride;
      const uint8_t* alpha_row = alpha_data + y * alpha_stride;
      for (x = 0; x < width; ++x) {
        rgb_row[4 * x + 3] = alpha_row[x];
      }
    }
  }

  if (verbose)
    StopwatchReadAndReset(&stop_watch);

#ifdef HAVE_WINCODEC_H
  needs_open_file = 0;
#endif
  if (needs_open_file) {
    fout = fopen(out_file, "wb");
    if (!fout) {
      fprintf(stderr, "Error opening output file %s\n", out_file);
      return;
    }
  }

#ifdef HAVE_WINCODEC_H
  ok &= WritePNG(out_file, buffer);
#else
  ok &= WritePNG(fout, buffer);
#endif

  if (fout) {
    fclose(fout);
  }
  if (ok) {
    printf("Saved file %s\n", out_file);
    if (verbose) {
      const double time = StopwatchReadAndReset(&stop_watch);
      printf("Time to write output: %.3fs\n", time);
    }
  } else {
    fprintf(stderr, "Error writing file %s !!\n", out_file);
  }
}

static void Help(void) {
  printf("Usage: webp2png in_file [options] [-o out_file]\n\n"
         "Decodes the WebP image file to PNG format\n"
         " Other options are:\n"
         "  -version  .... print version number and exit.\n"
         "  -h     ....... this help message.\n"
         "  -v     ....... verbose (e.g. print encoding/decoding times)\n"
        );
}

static const char* const kStatusMessages[] = {
  "OK", "OUT_OF_MEMORY", "INVALID_PARAM", "BITSTREAM_ERROR",
  "UNSUPPORTED_FEATURE", "SUSPENDED", "USER_ABORT", "NOT_ENOUGH_DATA"
};

int main(int argc, const char *argv[]) {
  const char *in_file = NULL;
  const char *out_file = NULL;
  uint8_t* alpha_data = NULL;
  int c;
  WebPDecoderConfig config;
  WebPDecBuffer* const output_buffer = &config.output;

  if (!WebPInitDecoderConfig(&config)) {
    fprintf(stderr, "Library version mismatch!\n");
    return -1;
  }

  // TODO: Split this function into smaller ones like ParseCommandLine,
  // ReadInput, Decode and SaveOutput.
  for (c = 1; c < argc; ++c) {
    if (!strcmp(argv[c], "-h") || !strcmp(argv[c], "-help")) {
      Help();
      return 0;
    } else if (!strcmp(argv[c], "-o") && c < argc - 1) {
      out_file = argv[++c];
    } else if (!strcmp(argv[c], "-version")) {
      const int version = WebPGetDecoderVersion();
      printf("%d.%d.%d\n",
        (version >> 16) & 0xff, (version >> 8) & 0xff, version & 0xff);
      return 0;
    } else if (!strcmp(argv[c], "-v")) {
      verbose = 1;
    } else if (argv[c][0] == '-') {
      printf("Unknown option '%s'\n", argv[c]);
      Help();
      return -1;
    } else {
      in_file = argv[c];
    }
  }

  if (in_file == NULL) {
    printf("missing input file!!\n");
    Help();
    return -1;
  }

  {
    Stopwatch stop_watch;
    WebPMux* mux = NULL;
    WebPMuxError err = WEBP_MUX_OK;
    VP8StatusCode status = VP8_STATUS_OK;
    int ok;
    uint32_t data_size = 0;
    const uint8_t* rgb_data = NULL;
    const uint8_t* compressed_alpha = NULL;
    uint32_t rgb_size = 0;
    uint32_t compressed_alpha_size = 0;

    void* data = NULL;
    FILE* const in = fopen(in_file, "rb");

    if (!in) {
      fprintf(stderr, "cannot open input file '%s'\n", in_file);
      return 1;
    }
    fseek(in, 0, SEEK_END);
    data_size = ftell(in);
    fseek(in, 0, SEEK_SET);
    data = malloc(data_size);
    ok = (fread(data, data_size, 1, in) == 1);
    fclose(in);
    if (!ok) {
      fprintf(stderr, "Could not read %d bytes of data from file %s\n",
              data_size, in_file);
      free(data);
      return -1;
    }

    mux = WebPMuxCreate((const uint8_t*)data, data_size, 1);
    if (mux == NULL) goto End;

    err = WebPMuxGetImage(mux, &rgb_data, &rgb_size,
                          &compressed_alpha, &compressed_alpha_size);
    if (err != WEBP_MUX_OK) goto End;

    if (verbose) {
      StopwatchReadAndReset(&stop_watch);
    }

#ifdef HAVE_WINCODEC_H
    output_buffer->colorspace = compressed_alpha ? MODE_BGRA : MODE_BGR;
#else
    output_buffer->colorspace = compressed_alpha ? MODE_RGBA : MODE_RGB;
#endif

    // TODO: Handle Alpha Decoding inside WebPDecode.
    status = WebPDecode(rgb_data, rgb_size, &config);

    if (status == VP8_STATUS_OK) {
      if (compressed_alpha != NULL) {
        const int width = output_buffer->width;
        const int height = output_buffer->height;

        alpha_data = (uint8_t*)malloc(height * width);
        if (alpha_data == NULL) {
          printf("Error: Could not allocate memory\n");
          goto End;
        }

        if (!DecodeAlpha(compressed_alpha, compressed_alpha_size,
                         width, height, width, alpha_data)) {
          goto End;
        }
      }
    }

    if (verbose) {
      const double time = StopwatchReadAndReset(&stop_watch);
      printf("Time to decode picture: %.3fs\n", time);
    }

 End:
    free(data);
    WebPMuxDelete(mux);

    ok = (status == VP8_STATUS_OK);
    if (!ok) {
      fprintf(stderr, "Decoding of %s failed.\n", in_file);
      fprintf(stderr, "Status: %d (%s)\n", status, kStatusMessages[status]);
      return -1;
    }
  }

  if (out_file) {
    printf("Decoded %s. Dimensions: %d x %d%s. Now saving...\n", in_file,
           output_buffer->width, output_buffer->height,
           alpha_data ? " (with alpha)" : "");
    SaveOutput(output_buffer, alpha_data, out_file);
  } else {
    printf("File %s can be decoded (dimensions: %d x %d)%s.\n",
           in_file, output_buffer->width, output_buffer->height,
           alpha_data ? " (with alpha)" : "");
    printf("Nothing written; use -o flag to save the result as e.g. PNG.\n");
  }
  WebPFreeDecBuffer(output_buffer);
  free(alpha_data);

  return 0;
}

//------------------------------------------------------------------------------

#endif  // WEBP_EXPERIMENTAL_FEATURES
