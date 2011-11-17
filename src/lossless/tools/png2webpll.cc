// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <set>
#include <string>
#include <iterator>
#include <fstream>
#include <streambuf>
#include "../dec/decode.h"
#include "../enc/webp_bit_stream.h"
#include "png.h"
#include "pngconf.h"

// Take in a png data in a string. Give an easy access to the png
// data through the methods.
//
// See http://www.faqs.org/rfcs/rfc2083.html for description of png values.
class PngDecoder {
 public:
  explicit PngDecoder(const std::string &png, bool premultiplied_alpha);
  virtual ~PngDecoder();

  // Cause subsequent calls to GetPixelArgb and GetPixelsInArgbArray to return
  // pixels whose alpha channel contains transparency information (255 == fully
  // transparent) rather than opacity information (255 == fully opaque). This is
  // useful when dealing with libraries like ImageMagick that have an inverted
  // view of things.
  //
  // If you call this method twice, it will invert back to normal.
  void InvertAlphaChannel();

  int width() const;
  int height() const;

  // Returns true when an error occured in processing.
  bool error() const;

  // This returns the pixel at the given location.
  uint32 GetPixelArgb(int x, int y) const;

  // Returns all pixels starting in an array, indexed by y * x_size + x.
  const uint32 *GetPixelsInArgbArray() const;

 private:
  // Convert the png into the internal argb buffer.
  void ConvertToArgb();

  void ConvertToPremultipliedAlpha();

  int GetPixelComponent(const int x,
                        const int y,
                        const int component) const;
  int components() const;

  // Pimpl pattern to hide the hideous third_party/png (libpng) interface.
  struct Impl;
  Impl *pimpl_;

  bool error_;
};

// PngDecoder::Impl hides the png state from the users.
struct PngDecoder::Impl {
  Impl() : png_ptr(0), info_ptr(0), end_info(0), row_pointers(0) {}

  void ConvertToArgb();

  png_structp png_ptr;
  png_infop info_ptr;
  png_infop end_info;
  png_bytep *row_pointers;
  std::vector<uint32> pixels_;
};

int PngDecoder::width() const {
  return png_get_image_width(pimpl_->png_ptr, pimpl_->info_ptr);
}

int PngDecoder::height() const {
  return png_get_image_height(pimpl_->png_ptr, pimpl_->info_ptr);
}

int PngDecoder::components() const {
  return png_get_channels(pimpl_->png_ptr, pimpl_->info_ptr);
}

void PngDecoder::InvertAlphaChannel() {
  // Invert the alpha byte of each pixel. We must shift 255 (= 0xff) left 24
  // bits so that it aligns with the alpha byte of the pixel. Note that
  // subtracting the alpha value from 255 is equivalent to XOR'ing with 255
  // since 255 is 11111111 in binary.
  uint32 *pixels = &pimpl_->pixels_[0];
  const int size = pimpl_->pixels_.size();

  for (int i = 0; i < size; ++i) {
    pixels[i] ^= 0xff << 24;
  }
}

inline int PngDecoder::GetPixelComponent(const int x,
                                         const int y,
                                         const int component) const {
  return pimpl_->row_pointers[y][x * components() + component];
}

uint32 PngDecoder::GetPixelArgb(const int x, const int y) const {
  return pimpl_->pixels_[y * width() + x];
}

const uint32 *PngDecoder::GetPixelsInArgbArray() const {
  return &pimpl_->pixels_[0];
}

static inline uint8 MultiplyAlpha(const uint8 a, const uint8 color) {
  const uint32 rounded_mul = (a * color) + 128;
  return static_cast<uint8>((rounded_mul + (rounded_mul >> 8)) >> 8);
}

void PngDecoder::ConvertToPremultipliedAlpha() {
  uint32 *p = &pimpl_->pixels_[0];
  for (int i = 0; i < pimpl_->pixels_.size(); ++i) {
    const uint32 v = p[i];
    const uint8 a = static_cast<uint8>(v >> 24);
    const uint8 r = MultiplyAlpha(a, (v >> 16) & 0xff);
    const uint8 g = MultiplyAlpha(a, (v >> 8) & 0xff);
    const uint8 b = MultiplyAlpha(a, (v >> 0) & 0xff);
    p[i] = (a << 24) | (r << 16) | (g << 8) | b;
  }
}

void PngDecoder::ConvertToArgb() {
  const int xsize = width();
  const int ysize = height();
  pimpl_->pixels_.clear();
  pimpl_->pixels_.resize(xsize * ysize);

  switch (components()) {
    case 1: {
      // Indexcolor or gray-scale.
      png_bytep trans = 0;
      int num_trans = 0;
      png_color_16* trans_color = 0;
      png_get_tRNS(pimpl_->png_ptr, pimpl_->info_ptr,
                   &trans, &num_trans, &trans_color);

      png_colorp palette;
      int num_palette;
      png_color_8_struct palette_with_alpha[256];
      if (png_get_PLTE(pimpl_->png_ptr, pimpl_->info_ptr,
                       &palette, &num_palette)) {
        // We have a PLTE tag, the image is palettized.
        for (int i = 0; i < num_palette; ++i) {
          palette_with_alpha[i].red = palette[i].red;
          palette_with_alpha[i].green = palette[i].green;
          palette_with_alpha[i].blue = palette[i].blue;
          palette_with_alpha[i].alpha = 255;
        }
      } else {
        // We do not have a PLTE tag, the image is grayscale.
        for (int i = 0; i < 256; ++i) {
          palette_with_alpha[i].red = i;
          palette_with_alpha[i].green = i;
          palette_with_alpha[i].blue = i;
          palette_with_alpha[i].alpha = 255;
        }
      }
      for (int i = 0; i < num_trans; ++i) {
        palette_with_alpha[trans[i]].alpha = 0;
      }
      for (int y = 0; y < ysize; ++y) {
        for (int x = 0; x < xsize; ++x) {
          const uint8 i = GetPixelComponent(x, y, 0);
          const uint8 r = palette_with_alpha[i].red;
          const uint8 g = palette_with_alpha[i].green;
          const uint8 b = palette_with_alpha[i].blue;
          const uint8 a = palette_with_alpha[i].alpha;
          const uint32 pixel = (a << 24) | (r << 16) | (g << 8) | b;
          pimpl_->pixels_[y * xsize + x] = pixel;
        }
      }
      return;
    }
    case 2:
      // Grayscale with alpha.
      for (int y = 0; y < ysize; ++y) {
        for (int x = 0; x < xsize; ++x) {
          const uint8 gray = GetPixelComponent(x, y, 0);
          const uint8 a = GetPixelComponent(x, y, 1);
          const uint32 pixel = (a << 24) | (gray * 0x10101);
          pimpl_->pixels_[y * xsize + x] = pixel;
        }
      }
      return;
    case 3: {
      // RGB
      for (int y = 0; y < ysize; ++y) {
        for (int x = 0; x < xsize; ++x) {
          const uint8 r = GetPixelComponent(x, y, 0);
          const uint8 g = GetPixelComponent(x, y, 1);
          const uint8 b = GetPixelComponent(x, y, 2);
          const uint8 a = 255;
          const uint32 pixel = (a << 24) | (r << 16) | (g << 8) | b;
          pimpl_->pixels_[y * xsize + x] = pixel;
        }
      }
      return;
    }
    case 4: {
      // RGBA
      for (int y = 0; y < ysize; ++y) {
        for (int x = 0; x < xsize; ++x) {
          const uint8 r = GetPixelComponent(x, y, 0);
          const uint8 g = GetPixelComponent(x, y, 1);
          const uint8 b = GetPixelComponent(x, y, 2);
          const uint8 a = GetPixelComponent(x, y, 3);
          const uint32 pixel = (a << 24) | (r << 16) | (g << 8) | b;
          pimpl_->pixels_[y * xsize + x] = pixel;
        }
      }
      return;
    }
    default:
      // Should be never executed.
      error_ = true;
      return;
  }
}

struct PngStringReaderClosure {
  explicit PngStringReaderClosure(const std::string& str) : str_(str), pos_(0) {}
  const std::string& str_;
  int pos_;
};

void StrStreamCallback(png_structp png_ptr, png_bytep data, png_size_t size) {
  PngStringReaderClosure* png_closure =
      static_cast<PngStringReaderClosure*>(png_get_io_ptr(png_ptr));
  memcpy(data, &png_closure->str_[png_closure->pos_], size);
  png_closure->pos_ += size;
}

bool PngDecoder::error() const {
  return error_;
}

PngDecoder::PngDecoder(const std::string &png_data_str,
                       const bool premultiplied_alpha)
    : pimpl_(new Impl),
      error_(false) {
  pimpl_->png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
                                           NULL, NULL, NULL);
  VERIFY(pimpl_->png_ptr != NULL);
  pimpl_->info_ptr = png_create_info_struct(pimpl_->png_ptr);
  VERIFY(pimpl_->info_ptr != NULL);
  pimpl_->end_info = png_create_info_struct(pimpl_->png_ptr);

  PngStringReaderClosure png_reader(png_data_str);
  if (setjmp(png_jmpbuf(pimpl_->png_ptr)) == 0) {
    png_set_read_fn(pimpl_->png_ptr, &png_reader, StrStreamCallback);

    // The png_transforms flags are as follows:
    // packing == convert 1,2,4 bit images,
    // strip == 16 -> 8 bits / channel,
    // shift == use sBIT dynamics, and
    // expand == palettes -> rgb, grayscale -> 8 bit images, tRNS -> alpha.
    const unsigned int png_transforms =
        PNG_TRANSFORM_PACKING |
        PNG_TRANSFORM_EXPAND |
        PNG_TRANSFORM_STRIP_16;

    png_read_png(pimpl_->png_ptr, pimpl_->info_ptr, png_transforms, NULL);
    pimpl_->row_pointers = png_get_rows(pimpl_->png_ptr, pimpl_->info_ptr);
    ConvertToArgb();
    if (premultiplied_alpha) {
      ConvertToPremultipliedAlpha();
    }
  } else {
    // Ok we are here because of the setjmp.
    error_ = true;
  }
}

PngDecoder::~PngDecoder() {
  png_destroy_read_struct(&pimpl_->png_ptr,
                          &pimpl_->info_ptr,
                          &pimpl_->end_info);
  delete pimpl_;
}

static int UniquePixels(int n, const uint32 *argb) {
  std::set<uint32> k;
  k.insert(argb[0]);
  for (int i = 1; i < n; ++i) {
    if (argb[i - 1] == argb[i]) {
      continue;
    }
    if (k.find(argb[i]) == k.end()) {
      k.insert(argb[i]);
      if (k.size() >= 65536) {
        // enough...
        break;
      }
    }
  }
  return k.size();
}

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

static void WriteStringToFile(const char *path, const std::string &str) {
  FILE *fp = fopen(path, "wb");
  VERIFY(fp != 0);
  VERIFY(fwrite(&str[0], str.size(), 1, fp) == 1);
  VERIFY(fclose(fp) == 0);
}

static int FindClosestDiscretized(int a, int bits) {
  int best_goodness = 514;
  int best_val = 0;
  for (int i = -1; i < 2; ++i) {
    int val = a + i * (1 << bits);
    if (val < 0) {
      val = 0;
    }
    if (val > 255) {
      val = 255;
    }
    int candidate = ((val >> bits) << bits) | (val >> (8 - bits));
    // Smallest distance but favor i == 0 over i == -1 and i == 1
    // since that keeps the overall intensity more constant in the
    // images.
    int goodness = 2 * abs(a - candidate) + abs(i);
    if (best_goodness > goodness) {
      best_goodness = goodness;
      best_val = candidate;
    }
  }
  return best_val;
}

static uint32 ClosestDiscretizedArgb(uint32 a, int bits) {
  return (FindClosestDiscretized(a >> 24, bits) << 24) |
      (FindClosestDiscretized((a >> 16) & 0xff, bits) << 16) |
      (FindClosestDiscretized((a >> 8) & 0xff, bits) << 8) |
      (FindClosestDiscretized(a & 0xff, bits));
}

static bool IsFar(uint32 a, uint32 b, int limit) {
  for(int k = 0; k < 4; ++k) {
    int delta = int((a >> (k * 8)) & 0xff) - int((b >> (k * 8)) & 0xff);
    if (delta >= limit || delta <= -limit) {
      return true;
    }
  }
  return false;
}

static void NearLossless(int xsize, int ysize, uint32 *argb, int limit_bits) {
  std::vector<uint32> copy(argb, argb + xsize * ysize);
  int limit = 1 << limit_bits;
  for (int y = 0; y < ysize; ++y) {
    for (int x = 0; x < xsize; ++x) {
      int ix = y * xsize + x;
      // Check that all pixels in 4-connected neighborhood are smooth.
      bool smooth_area = true;
      if (x != 0 && IsFar(copy[ix], copy[ix - 1], limit)) {
        smooth_area = false;
      } else if (y != 0 && IsFar(copy[ix], copy[ix - xsize], limit)) {
        smooth_area = false;
      } else if (x != xsize - 1 && IsFar(copy[ix], copy[ix + 1], limit)) {
        smooth_area = false;
      } else if (y != ysize - 1 && IsFar(copy[ix], copy[ix + xsize], limit)) {
        smooth_area = false;
      }
      if (!smooth_area) {
        argb[ix] = ClosestDiscretizedArgb(argb[ix], limit_bits);
      }
    }
  }
}

int main(int argc, char **argv) {
  char *in_path = NULL;
  char *out_path = NULL;
  bool verbose = false;
  int quality = 95;
  bool verify = true;
  bool show_usage = argc < 2;
  int near_lossless = 0;
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i], "-o")) {
      if (i < argc - 1) {
        out_path = argv[i + 1];
        ++i;
      } else {
        show_usage = true;
      }
    } else if (!strcmp(argv[i], "-c")) {
      if (i < argc - 1) {
        quality = atoi(argv[i + 1]);
        ++i;
      } else {
        show_usage = true;
      }
    } else if (!strcmp(argv[i], "-n")) {
      if (i < argc - 1) {
        near_lossless = atoi(argv[i + 1]);
        ++i;
      } else {
        show_usage = true;
      }
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
            "Usage: %s xyzzy.png [-o xyzzy.webpll] "
            "[-c 0..100] [-n 1..4] [-h]\n"
            "  -c makes output denser and slower (0 == fastest)\n"
            "  -n X applies preprocessing that selectively throws "
            "away X least significant bits\n"
            "  -h shows the usage information\n",
            argv[0]);
    exit(1);
  }
  std::string alternate_out_path;
  if (!out_path) {
    alternate_out_path.assign(in_path);
    alternate_out_path.append(".webpll");
    out_path = &alternate_out_path[0];
  }

  std::string to_file;
  std::string input_file_str = ReadFileToString(in_path);
  size_t orig_size = input_file_str.size();
  PngDecoder png_decoder(input_file_str, false);
  assert(!png_decoder.error());
  const int xsize = png_decoder.width();
  const int ysize = png_decoder.height();
  const uint32 *argb_orig = png_decoder.GetPixelsInArgbArray();

  VERIFY(near_lossless >= 0);
  VERIFY(near_lossless < 5);
  std::vector<uint32> near_lossless_argb(argb_orig, argb_orig + xsize * ysize);
  if (near_lossless) {
    for (int k = near_lossless; k != 0; --k) {
      NearLossless(xsize, ysize, &near_lossless_argb[0], k);
    }
  }

  // Modes:
  // 0: normal predicted mode with lots of entropy codes, good for photos
  // 1: non-color-predicted mode
  // 2: no spatial predict, color
  // 3: mode with a single palette

  bool try_with_small_palette = UniquePixels(xsize * ysize, argb_orig) <= 256;

  int histogram_bits = 3;
  if (quality < 50) ++histogram_bits;
  if (quality < 10) ++histogram_bits;
  if (quality < 1) ++histogram_bits;

  for (int mode = 0; mode < 5; mode++) {
    if (mode >= 3 && !try_with_small_palette) {
      continue;
    }
    char *bytes;
    int n_bytes;
    const bool use_small_palette = mode == 3 || mode == 4;
    const bool use_spatial_predict = mode == 0 || mode == 4;
    const bool use_cross_color_transform = mode == 0 || mode == 2;
    // If the colors fit in the palette, there is rarely a benefit for doing
    // a near lossless using the method implemented here.
    const uint32 *argb_to_compress =
        use_small_palette ?
        png_decoder.GetPixelsInArgbArray() :
        &near_lossless_argb[0];
    EncodeWebpLLImage(png_decoder.width(),
                      png_decoder.height(),
                      argb_to_compress,
                      quality,
                      use_small_palette,
                      use_spatial_predict,
                      4,
                      histogram_bits,
                      use_cross_color_transform,
                      mode == 0 ? 4 : 10,
                      false,  // no error detection bits
                      &n_bytes,
                      &bytes);
    std::string to_file_candidate(bytes, bytes + n_bytes);
    free(bytes);
    if (mode == 0 || to_file.size() > to_file_candidate.size()) {
      to_file = to_file_candidate;
    }
    if (verbose) {
      fprintf(stderr,
              "Strategy %d: %d bytes\n", mode, int(to_file_candidate.size()));
    }
    if (verify) {
      int xsize_verify;
      int ysize_verify;
      uint32* argb_image_verify;
      if (!DecodeWebpLLImage(to_file_candidate.size(),
                             (uint8*)to_file_candidate.data(),
                             &xsize_verify,
                             &ysize_verify,
                             &argb_image_verify)) {
        fprintf(stderr, "Error while decoding input image.\n");
        return 1;
      }
      VERIFY(xsize == xsize_verify);
      VERIFY(ysize == ysize_verify);
      for (int i = 0; i < xsize * ysize; ++i) {
        if (argb_image_verify[i] != argb_to_compress[i]) {
          fprintf(stderr,
                  "pixel at differs for --in=%s mode=%d %08x != %08x",
                 in_path,
                 mode,
                 argb_image_verify[i],
                 argb_to_compress[i]);
          abort();
        }
      }
      free(argb_image_verify);
    }
    if (quality == 0) {
      break;
    }
  }
  WriteStringToFile(out_path, to_file);
  size_t new_size = to_file.size();
  fprintf(stderr, "png = %zd, webpll = %zd bytes (%g)\n",
          orig_size, new_size, double(new_size)/orig_size);
  return 0;
}
