// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "../common/integral_types.h"
#include "../dec/decode.h"
#include "png.h"

static std::string ReadFileToString(const char *path) {
  FILE *fp = fopen(path, "rb");
  VERIFY(fp != 0);
  VERIFY(fseek(fp, 0, SEEK_END) == 0);
  size_t len = ftell(fp);
  VERIFY(fseek(fp, 0, SEEK_SET) == 0);
  std::string retval(len, '\0');
  VERIFY(fread(&retval[0], len, 1, fp) == 1);
  VERIFY(fclose(fp) == 0);
  return retval;
}

static void PNGAPI error_function(png_structp png, png_const_charp dummy) {
  (void)dummy;  // remove variable-unused warning
  longjmp(png_jmpbuf(png), 1);
}

static int WritePng(const char *path, int xsize, int ysize, uint32_t *argb) {
  // Convert 32-bit ARGB to 8-bit RGBA for png encoding.
  uint8_t* const rgba_image = (uint8_t *)malloc(xsize * ysize * 4);
  uint8_t* rgba = rgba_image;
  for (int i = 0; i < xsize * ysize; ++i) {
    *rgba++ = (argb[i] >> 16) & 0xff;
    *rgba++ = (argb[i] >> 8) & 0xff;
    *rgba++ = (argb[i] >> 0) & 0xff;
    *rgba++ = (argb[i] >> 24) & 0xff;
  }
  // Write the png.
  png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                            NULL, error_function, NULL);
  if (png == NULL) {
    return 1;
  }
  png_infop info = png_create_info_struct(png);
  if (info == NULL) {
    png_destroy_write_struct(&png, NULL);
    return 1;
  }
  if (setjmp(png_jmpbuf(png))) {
    png_destroy_write_struct(&png, &info);
    return 1;
  }
  FILE *out_file = fopen(path, "wb");
  VERIFY(out_file != NULL);
  png_init_io(png, out_file);
  png_set_IHDR(png, info, xsize, ysize, 8,
               PNG_COLOR_TYPE_RGBA,
               PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
               PNG_FILTER_TYPE_DEFAULT);
  png_write_info(png, info);
  for (int y = 0; y < ysize; ++y) {
    png_bytep row = (png_bytep)(&rgba_image[y * xsize * 4]);
    png_write_rows(png, &row, 1);
  }
  png_write_end(png, info);
  png_destroy_write_struct(&png, &info);
  VERIFY(fclose(out_file) == 0);
  free(rgba_image);
  return 0;
}

int main(int argc, char **argv) {
  char *in_path = NULL;
  char *out_path = NULL;
  int timing_iters = 1;
  int skip_png = 0;
  bool show_usage = argc < 2;
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i], "-o")) {
      if (i < argc - 1) {
        out_path = argv[i + 1];
        ++i;
      } else {
        show_usage = true;
      }
    } else if (!strcmp(argv[i], "-t")) {
      if (i < argc - 1) {
        timing_iters = atoi(argv[i + 1]);
        ++i;
      } else {
        show_usage = true;
      }
    } else if (!strcmp(argv[i], "-s")) {
      skip_png = 1;
    } else if (!strcmp(argv[i], "-h") ||
               !strcmp(argv[i], "--help")) {
      show_usage = true;
      break;
    } else {
      if (argv[i][0] == '-') {
        fprintf(stderr, "%s: unrecognized option -- '%s'\n", argv[0], argv[i]);
        show_usage = true;
        break;
      }
      if (in_path != NULL) {
        fprintf(stderr,
                "only one input file name allowed, two given:"
                "'%s' and '%s'", in_path, argv[i]);
        exit(1);
      }
      in_path = argv[i];
    }
  }
  if (show_usage || in_path == NULL || strlen(in_path) == 0) {
    fprintf(stderr,
            "Usage: %s xyzzy.webpll -o xyzzy.png [-s] [-t 1..n] [-h]\n"
            "  -s for skipping the output phase\n"
            "  -t X for running decoding X times "
            "(useful for timing tests)\n"
            "  -h shows the usage information\n",
            argv[0]);
    exit(1);
  }
  std::string alternate_out_path;
  if (!out_path) {
    alternate_out_path.assign(in_path);
    alternate_out_path.append(".png");
    out_path = &alternate_out_path[0];
  }
  std::string input_data = ReadFileToString(in_path);
  int xsize = 0;
  int ysize = 0;
  uint32_t* argb_image;
  for (int i = 0; i < timing_iters; ++i) {
    if (i != 0) {
      free(argb_image);
    }
    if (!DecodeWebpLLImage(input_data.size(),
                           (uint8_t*)input_data.data(),
                           &xsize,
                           &ysize,
                           &argb_image)) {
      fprintf(stderr, "Error while decoding input image.\n");
      return 1;
    }
  }
  if (!skip_png) {
    WritePng(out_path, xsize, ysize, argb_image);
  }
  free(argb_image);
  return 0;
}
