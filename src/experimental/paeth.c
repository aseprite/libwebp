// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Spatial prediction using Paeth filter
//
// Author: Urvang (urvang@google.com)

#include "./paeth.h"
#include <assert.h>
#include <stddef.h>
#include <stdlib.h>

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

//------------------------------------------------------------------------------

static WEBP_INLINE int AbsDiff(int a, int b) {
  return (a > b) ? a - b : b - a;
}

static WEBP_INLINE uint8_t PaethPredictor(uint8_t a, uint8_t b, uint8_t c) {
  const int p = a + b - c;  // Base.
  const int pa = AbsDiff(p, a);
  const int pb = AbsDiff(p, b);
  const int pc = AbsDiff(p, c);

  // Return nearest to base of a, b, c.
  return (pa <= pb && pa <= pc) ? a : (pb <= pc) ? b : c;
}

void PaethFilter(const uint8_t* data, int width, int height, int bpp,
                 int stride, uint8_t* filtered_data) {
  int h;

  assert(data != NULL);
  assert(filtered_data != NULL);
  assert(width > 0);
  assert(height > 0);
  assert(bpp > 0);
  assert(stride >= width * bpp);

  // Filter line-by-line.
  for (h = 0; h < height; ++h) {
    int w;
    const uint8_t* scan_line = data + h * stride;
    uint8_t* out = filtered_data + h * stride;

    if (h == 0) {  // Top scan line (special case).
      for (w = 0; w < bpp; ++w) {
        out[w] = scan_line[w];
      }
      for (w = bpp; w < width * bpp; ++w) {
        // Note: PaethPredictor(scan_line[w - bpp], 0, 0) == scan_line[w - bpp].
        out[w] = scan_line[w] - scan_line[w - bpp];
      }
    } else {  // General case.
      const uint8_t* prev_line = scan_line - stride;
      for (w = 0; w < bpp; ++w) {
        // Note: PaethPredictor(0, prev_line[w], 0) == prev_line[w].
        out[w] = scan_line[w] - prev_line[w];
      }
      for (w = bpp; w < width * bpp; ++w) {
        out[w] = scan_line[w] - PaethPredictor(scan_line[w - bpp], prev_line[w],
                                               prev_line[w - bpp]);
      }
    }
  }
}

void PaethReconstruct(const uint8_t* data, int width, int height, int bpp,
                      int stride, uint8_t* recon_data) {
  int h;

  assert(data != NULL);
  assert(recon_data != NULL);
  assert(width > 0);
  assert(height > 0);
  assert(bpp > 0);
  assert(stride >= width * bpp);

  // Reconstruct line-by-line.
  for (h = 0; h < height; ++h) {
    int w;
    const uint8_t* scan_line = data + h * stride;
    uint8_t* out = recon_data + h * stride;

    if (h == 0) {  // Top scan line (special case).
      for (w = 0; w < bpp; ++w) {
        out[w] = scan_line[w];
      }
      for (w = bpp; w < width * bpp; ++w) {
        // Note: PaethPredictor(out[w - bpp], 0, 0) == out[w - bpp].
        out[w] = scan_line[w] + out[w - bpp];
      }
    } else {  // General case.
      const uint8_t* out_prev = out - stride;
      for (w = 0; w < bpp; ++w) {
        // Note: PaethPredictor(0, out_prev[w], 0) == out_prev[w].
        out[w] = scan_line[w] + out_prev[w];
      }
      for (w = bpp; w < width * bpp; ++w) {
        out[w] = scan_line[w] + PaethPredictor(out[w - bpp], out_prev[w],
                                               out_prev[w - bpp]);
      }
    }
  }
}

//------------------------------------------------------------------------------

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif
