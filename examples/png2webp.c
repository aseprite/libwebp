// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
//  Command line to convert PNG image to WebP, converting the Alpha plane to
//  WebP-container.
//
// Author: Vikas Arora (vikasa@google.com)

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

#ifndef GUID_WICPixelFormat24bppRGB
// From Microsoft SDK 7.0a
DEFINE_GUID(GUID_WICPixelFormat24bppRGB,
    0x6fddc324, 0x4e03, 0x4bfe, 0xb1, 0x85, 0x3d, 0x77, 0x76, 0x8d, 0xc9, 0x0d);
#endif
#ifndef GUID_WICPixelFormat32bppRGBA
DEFINE_GUID(GUID_WICPixelFormat32bppRGBA,
    0xf5c7ad2d, 0x6a8d, 0x43dd, 0xa7, 0xa8, 0xa2, 0x99, 0x35, 0x26, 0x1a, 0xe9);
#endif
#endif  /* HAVE_WINCODEC_H */

#include "utils/alpha.h"
#include "webp/encode.h"
#include "stopwatch.h"

//------------------------------------------------------------------------------

static int verbose = 0;

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

static HRESULT OpenInputStream(const char* filename, IStream** ppStream) {
  HRESULT hr = S_OK;
  IFS(SHCreateStreamOnFileA(filename, STGM_READ, ppStream));
  if (FAILED(hr))
    printf("Error opening input file %s (%08x)\n", filename, hr);
  return hr;
}

static HRESULT ReadPictureWithWIC(const char* filename,
                                  WebPPicture* const pic, int keep_alpha) {
  HRESULT hr = S_OK;
  IWICBitmapFrameDecode* pFrame = NULL;
  IWICFormatConverter* pConverter = NULL;
  IWICImagingFactory* pFactory = NULL;
  IWICBitmapDecoder* pDecoder = NULL;
  IStream* pStream = NULL;
  UINT frameCount = 0;
  UINT width, height = 0;
  BYTE* rgb = NULL;
  WICPixelFormatGUID srcPixelFormat = { 0 };
  GUID srcContainerFormat = { 0 };
  const GUID* alphaContainers[] = {
    &GUID_ContainerFormatBmp,
    &GUID_ContainerFormatPng,
    &GUID_ContainerFormatTiff
  };
  int has_alpha = 0;
  int i, stride;

  IFS(CoInitialize(NULL));
  IFS(CoCreateInstance(MAKE_REFGUID(CLSID_WICImagingFactory), NULL,
          CLSCTX_INPROC_SERVER, MAKE_REFGUID(IID_IWICImagingFactory),
          (LPVOID*)&pFactory));
  if (hr == REGDB_E_CLASSNOTREG) {
    printf("Couldn't access Windows Imaging Component (are you running \n");
    printf("Windows XP SP3 or newer?). Most formats not available.\n");
    printf("Use -s for the available YUV input.\n");
  }
  // Prepare for image decoding.
  IFS(OpenInputStream(filename, &pStream));
  IFS(IWICImagingFactory_CreateDecoderFromStream(pFactory, pStream, NULL,
          WICDecodeMetadataCacheOnDemand, &pDecoder));
  IFS(IWICBitmapDecoder_GetFrameCount(pDecoder, &frameCount));
  if (SUCCEEDED(hr) && frameCount == 0) {
    printf("No frame found in input file.\n");
    hr = E_FAIL;
  }
  IFS(IWICBitmapDecoder_GetFrame(pDecoder, 0, &pFrame));
  IFS(IWICBitmapFrameDecode_GetPixelFormat(pFrame, &srcPixelFormat));
  IFS(IWICBitmapDecoder_GetContainerFormat(pDecoder, &srcContainerFormat));

  has_alpha = keep_alpha;
  for (i = 0;
       has_alpha && i < sizeof(alphaContainers)/sizeof(alphaContainers[0]);
       ++i) {
    if (IsEqualGUID(&srcContainerFormat, alphaContainers[i])) {
      has_alpha =
          IsEqualGUID(&srcPixelFormat, &GUID_WICPixelFormat32bppRGBA) ||
          IsEqualGUID(&srcPixelFormat, &GUID_WICPixelFormat32bppBGRA);
      break;
    }
  }

  // Prepare for pixel format conversion (if necessary).
  IFS(IWICImagingFactory_CreateFormatConverter(pFactory, &pConverter));
  IFS(IWICFormatConverter_Initialize(pConverter, (IWICBitmapSource*)pFrame,
          has_alpha ? MAKE_REFGUID(GUID_WICPixelFormat32bppRGBA)
                    : MAKE_REFGUID(GUID_WICPixelFormat24bppRGB),
          WICBitmapDitherTypeNone,
          NULL, 0.0, WICBitmapPaletteTypeCustom));

  // Decode.
  IFS(IWICFormatConverter_GetSize(pConverter, &width, &height));
  stride = (has_alpha ? 4 : 3) * width * sizeof(*rgb);
  if (SUCCEEDED(hr)) {
    rgb = (BYTE*)malloc(stride * height);
    if (rgb == NULL)
      hr = E_OUTOFMEMORY;
  }
  IFS(IWICFormatConverter_CopyPixels(pConverter, NULL, stride,
          stride * height, rgb));

  // WebP conversion.
  if (SUCCEEDED(hr)) {
    int ok;
    pic->width = width;
    pic->height = height;
    ok = has_alpha ? WebPPictureImportRGBA(pic, rgb, stride)
                   : WebPPictureImportRGB(pic, rgb, stride);
    if (!ok)
      hr = E_FAIL;
  }

  // Cleanup.
  if (pConverter != NULL) IUnknown_Release(pConverter);
  if (pFrame != NULL) IUnknown_Release(pFrame);
  if (pDecoder != NULL) IUnknown_Release(pDecoder);
  if (pFactory != NULL) IUnknown_Release(pFactory);
  if (pStream != NULL) IUnknown_Release(pStream);
  free(rgb);
  return hr;
}

static int ReadPicture(const char* const filename, WebPPicture* const pic,
                       int keep_alpha) {
  int ok = SUCCEEDED(ReadPictureWithWIC(filename, pic, keep_alpha));
  if (!ok) {
    fprintf(stderr, "Error! Could not process file %s\n", filename);
  }
  return ok;
}

#else  // !HAVE_WINCODEC_H

#ifdef WEBP_HAVE_PNG
static void PNGAPI error_function(png_structp png, png_const_charp dummy) {
  (void)dummy;  // remove variable-unused warning
  longjmp(png_jmpbuf(png), 1);
}

static int ReadPNG(FILE* in_file, WebPPicture* const pic, int keep_alpha) {
  png_structp png;
  png_infop info;
  int color_type, bit_depth, interlaced;
  int has_alpha;
  int num_passes;
  int p;
  int ok = 0;
  png_uint_32 width, height, y;
  int stride;
  uint8_t* rgb = NULL;

  png = png_create_read_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
  if (png == NULL) {
    goto End;
  }

  png_set_error_fn(png, 0, error_function, NULL);
  if (setjmp(png_jmpbuf(png))) {
 Error:
    png_destroy_read_struct(&png, NULL, NULL);
    free(rgb);
    goto End;
  }

  info = png_create_info_struct(png);
  if (info == NULL) goto Error;

  png_init_io(png, in_file);
  png_read_info(png, info);
  if (!png_get_IHDR(png, info,
                    &width, &height, &bit_depth, &color_type, &interlaced,
                    NULL, NULL)) goto Error;

  png_set_strip_16(png);
  png_set_packing(png);
  if (color_type == PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png);
  if (color_type == PNG_COLOR_TYPE_GRAY ||
      color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
    if (bit_depth < 8) {
      png_set_expand_gray_1_2_4_to_8(png);
    }
    png_set_gray_to_rgb(png);
  }
  if (png_get_valid(png, info, PNG_INFO_tRNS)) {
    png_set_tRNS_to_alpha(png);
    has_alpha = 1;
  } else {
    has_alpha = !!(color_type & PNG_COLOR_MASK_ALPHA);
  }

  if (!keep_alpha) {
    png_set_strip_alpha(png);
    has_alpha = 0;
  }

  num_passes = png_set_interlace_handling(png);
  png_read_update_info(png, info);

  stride = ((keep_alpha && has_alpha) ? 4 : 3) * width * sizeof(*rgb);
  rgb = (uint8_t*)malloc(stride * height);
  if (rgb == NULL) goto Error;
  for (p = 0; p < num_passes; ++p) {
    for (y = 0; y < height; ++y) {
      png_bytep row = rgb + y * stride;
      png_read_rows(png, &row, NULL, 1);
    }
  }
  png_read_end(png, info);
  png_destroy_read_struct(&png, &info, NULL);

  pic->width = width;
  pic->height = height;
  ok = (keep_alpha && has_alpha) ?
      WebPPictureImportRGBA(pic, rgb, stride) :
      WebPPictureImportRGB(pic, rgb, stride);
  free(rgb);

 End:
  return ok;
}
#else
static int ReadPNG(FILE* in_file, WebPPicture* const pic, int keep_alpha) {
  (void)in_file;
  (void)pic;
  (void)keep_alpha;
  printf("PNG support not compiled. Please install the libpng development "
         "package before building.\n");
  return 0;
}
#endif

typedef enum {
  PNG = 0,
  UNSUPPORTED,
} InputFileFormat;

static InputFileFormat GetImageType(FILE* in_file) {
  InputFileFormat format = UNSUPPORTED;
  unsigned int magic;
  unsigned char buf[4];

  if ((fread(&buf[0], 4, 1, in_file) != 1) ||
      (fseek(in_file, 0, SEEK_SET) != 0)) {
    return format;
  }

  magic = (buf[0] << 24) | (buf[1] << 16) | (buf[2] << 8) | buf[3];
  if (magic == 0x89504E47U) {
    format = PNG;
  }
  return format;
}

static int ReadPicture(const char* const filename, WebPPicture* const pic,
                       int keep_alpha) {
  int ok = 0;
  InputFileFormat format = UNSUPPORTED;
  FILE* in_file = fopen(filename, "rb");
  if (in_file == NULL) {
    fprintf(stderr, "Error! Cannot open input file '%s'\n", filename);
    return ok;
  }

  // Try to decode the specified PNG file.
  format = GetImageType(in_file);
  if (format == PNG) {
    ok = ReadPNG(in_file, pic, keep_alpha);
    if (!ok) {
      fprintf(stderr, "Error! Could not read PNG file %s\n", filename);
    }
  } else {
    fprintf(stderr, "Error! not a PNG file %s\n", filename);
  }

  fclose(in_file);
  return ok;
}

#endif  // !HAVE_WINCODEC_H

//------------------------------------------------------------------------------

typedef struct {
  uint8_t** mem;
  size_t    max_size;
  size_t*   size;
} WebPMemoryWriter;

static void InitMemoryWriter(WebPMemoryWriter* const writer) {
  *writer->mem = NULL;
  *writer->size = 0;
  writer->max_size = 0;
}

static int WebPMemoryWrite(const uint8_t* data, size_t data_size,
                           const WebPPicture* const picture) {
  WebPMemoryWriter* const w = (WebPMemoryWriter*)picture->custom_ptr;
  size_t next_size;
  if (w == NULL) {
    return 1;
  }
  next_size = (*w->size) + data_size;
  if (next_size > w->max_size) {
    uint8_t* new_mem;
    size_t next_max_size = w->max_size * 2;
    if (next_max_size < next_size) next_max_size = next_size;
    if (next_max_size < 8192) next_max_size = 8192;
    new_mem = (uint8_t*)malloc(next_max_size);
    if (new_mem == NULL) {
      return 0;
    }
    if ((*w->size) > 0) {
      memcpy(new_mem, *w->mem, *w->size);
    }
    free(*w->mem);
    *w->mem = new_mem;
    w->max_size = next_max_size;
  }
  if (data_size) {
    memcpy((*w->mem) + (*w->size), data, data_size);
    *w->size += data_size;
  }
  return 1;
}

//------------------------------------------------------------------------------

static void Help(void) {
  printf("Usage:\n");
  printf(" png2webp [options] in_file [-o out_file]\n\n");
#ifdef HAVE_WINCODEC_H
  printf("Windows builds can take as input any of the files handled by WIC\n");
#endif
  printf("options:\n");
  printf("  -h / -help  ............ help\n");
  printf("  -q <float> ............. quality factor (0:small..100:big).\n");
  printf("  -alpha_q <int> ......... Transparency-compression quality "
         "(0..100).\n");
  printf("  -alpha_method <int> .... Transparency-compression method.\n");
  printf("  -noalpha ............... discard any transparency information.\n");
  printf("\n");
  printf("  -quiet ................. don't print anything.\n");
  printf("  -version ............... print version number and exit.\n");
  printf("  -v ..................... verbose, e.g. print encoding/decoding "
         "times.\n");
  printf("\n");
}

//------------------------------------------------------------------------------
// Error messages

static const char* const kErrorMessages[] = {
  "OK",
  "OUT_OF_MEMORY: Out of memory allocating objects",
  "BITSTREAM_OUT_OF_MEMORY: Out of memory re-allocating byte buffer",
  "NULL_PARAMETER: NULL parameter passed to function",
  "INVALID_CONFIGURATION: configuration is invalid",
  "BAD_DIMENSION: Bad picture dimension. Maximum width and height "
  "allowed is 16383 pixels.",
  "PARTITION0_OVERFLOW: Partition #0 is too big to fit 512k.\n"
  "To reduce the size of this partition, try using less segments "
  "with the -segments option, and eventually reduce the number of "
  "header bits using -partition_limit. More details are available "
  "in the manual (`man cwebp`)",
  "PARTITION_OVERFLOW: Partition is too big to fit 16M",
  "BAD_WRITE: Picture writer returned an I/O error"
  "FILE_TOO_BIG: File would be too big to fit in 4G"
};

//------------------------------------------------------------------------------

int main(int argc, const char *argv[]) {
  const char *in_file = NULL, *out_file = NULL;
  FILE *out = NULL;
  uint8_t* compressed_rgb = NULL;
  size_t compressed_rgb_size = 0;
  int keep_alpha = 1;
  int c;
  int quiet = 0;

  WebPPicture picture;
  WebPConfig config;
  WebPAuxStats stats;
  WebPMemoryWriter wrt;
  Stopwatch stop_watch;

  if (!WebPPictureInit(&picture) || !WebPConfigInit(&config)) {
    fprintf(stderr, "Error! Version mismatch!\n");
    goto Error;
  }

  if (argc == 1) {
    Help();
    return 0;
  }

  // Set some high-quality for RGB WEBP compression as default.
  config.quality = 90;

  for (c = 1; c < argc; ++c) {
    if (!strcmp(argv[c], "-h") || !strcmp(argv[c], "-help")) {
      Help();
      return 0;
    } else if (!strcmp(argv[c], "-o") && c < argc - 1) {
      out_file = argv[++c];
    } else if (!strcmp(argv[c], "-q") && c < argc - 1) {
      config.quality = (float)strtod(argv[++c], NULL);
    } else if (!strcmp(argv[c], "-alpha_q") && c < argc - 1) {
      config.alpha_quality = strtol(argv[++c], NULL, 0);
    } else if (!strcmp(argv[c], "-alpha_method") && c < argc - 1) {
      config.alpha_compression = strtol(argv[++c], NULL, 0);
    } else if (!strcmp(argv[c], "-noalpha")) {
      keep_alpha = 0;
    } else if (!strcmp(argv[c], "-version")) {
      const int version = WebPGetEncoderVersion();
      printf("%d.%d.%d\n",
        (version >> 16) & 0xff, (version >> 8) & 0xff, version & 0xff);
      return 0;
    } else if (!strcmp(argv[c], "-quiet")) {
      quiet = 1;
    } else if (!strcmp(argv[c], "-v")) {
      verbose = 1;
    } else if (argv[c][0] == '-') {
      fprintf(stderr, "Error! Unknown option '%s'\n", argv[c]);
      Help();
      return -1;
    } else {
      in_file = argv[c];
    }
  }

  if (!WebPValidateConfig(&config)) {
    fprintf(stderr, "Error! Invalid configuration.\n");
    goto Error;
  }

  // Read the input
  if (verbose) {
    StopwatchReadAndReset(&stop_watch);
  }

  if (!ReadPicture(in_file, &picture, keep_alpha)) {
    fprintf(stderr, "Error! Cannot read input picture\n");
    goto Error;
  }

  picture.writer = WebPMemoryWrite;
  picture.custom_ptr = &wrt;
  wrt.mem = &compressed_rgb;
  wrt.size = &compressed_rgb_size;
  InitMemoryWriter(&wrt);
  picture.stats = &stats;

  if (verbose) {
    const double time = StopwatchReadAndReset(&stop_watch);
    fprintf(stderr, "Time to read input: %.3fs\n", time);
  }

  // Compress
  if (verbose) {
    StopwatchReadAndReset(&stop_watch);
  }

  if (!WebPEncode(&config, &picture)) {
    fprintf(stderr, "Error! Cannot encode picture as WebP\n");
    fprintf(stderr, "Error code: %d (%s)\n",
            picture.error_code, kErrorMessages[picture.error_code]);
    goto Error;
  }

  if (verbose) {
    const double time = StopwatchReadAndReset(&stop_watch);
    fprintf(stderr, "Time to encode picture: %.3fs\n", time);
  }

  if (out_file) {
    out = fopen(out_file, "wb");
    if (out == NULL) {
      fprintf(stderr, "Error! Cannot open output file '%s'\n", out_file);
      goto Error;
    } else {
      if (!quiet) {
        fprintf(stderr, "Saving file '%s'\n", out_file);
      }
    }

    if (fwrite(compressed_rgb, compressed_rgb_size, 1, out) != 1) {
      fprintf(stderr, "Error! Cannot save image '%s'\n", out_file);
      goto Error;
    }
  } else {
    out = NULL;
    if (!quiet) {
      fprintf(stderr, "No output file specified (no -o flag). Encoding will\n");
      fprintf(stderr, "be performed, but its results discarded.\n\n");
    }
  }

 Error:
  free(compressed_rgb);
  WebPPictureFree(&picture);
  if (out != NULL) {
    fclose(out);
  }

  return 0;
}

//------------------------------------------------------------------------------
