// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)
//
// Models a 2d image of histograms.

#include "histogram_image.h"

#include <stdio.h>
#include <string.h>

#include "histogram.h"
#include "backward_references.h"
#include "../common/integral_types.h"

int BuildHistogramImage(int xsize, int ysize,
                        int histobits,
                        int palettebits,
                        const PixOrCopy *backward_refs,
                        int backward_refs_size,
                        Histogram ***image_arg,
                        int *image_size) {
  int histo_xsize = histobits ? (xsize + (1 << histobits) - 1) >> histobits : 1;
  int histo_ysize = histobits ? (ysize + (1 << histobits) - 1) >> histobits : 1;
  int i;
  int x = 0;
  int y = 0;
  Histogram **image;
  *image_arg = NULL;
  *image_size = histo_xsize * histo_ysize;
  image = (Histogram **)calloc(*image_size, sizeof(*image));
  if (image == NULL) {
    return 0;
  }
  for (i = 0; i < *image_size; ++i) {
    image[i] = (Histogram *)malloc(sizeof(*image[i]));
    if (!image[i]) {
      int k;
      for (k = 0; k < *image_size; ++k) {
        free(image[k]);
      }
      free(image);
      return 0;
    }
    Histogram_Init(image[i], palettebits);
  }
  // x and y trace the position in the image.
  for (i = 0; i < backward_refs_size; ++i) {
    const PixOrCopy v = backward_refs[i];
    const int ix =
        histobits ? (y >> histobits) * histo_xsize + (x >> histobits) : 0;
    Histogram_AddSinglePixOrCopy(image[ix], v);
    x += PixOrCopy_Length(&v);
    while (x >= xsize) {
      x -= xsize;
      ++y;
    }
  }
  *image_arg = image;
  return 1;
}

int CombineHistogramImage(Histogram **in,
                          int in_size,
                          int quality,
                          int palettebits,
                          Histogram ***out_arg,
                          int *out_size) {
  int ok = 1;
  int i;
  unsigned int seed = 0;
  int tries_with_no_success = 0;
  int inner_iters = 10 + quality / 2;
  int iter;
  double *bit_costs = (double *)malloc(in_size * sizeof(*bit_costs));
  Histogram **out = (Histogram **)calloc(in_size, sizeof(*out));
  *out_arg = out;
  *out_size = in_size;
  if (bit_costs == NULL || out == NULL) {
    ok = 0;
    goto exit_label;
  }
  // Copy
  for (i = 0; i < in_size; ++i) {
    Histogram *new_histo = (Histogram *)malloc(sizeof(*new_histo));
    if (new_histo == NULL) {
      ok = 0;
      goto exit_label;
    }
    Histogram_Init(new_histo, palettebits);
    *new_histo = *(in[i]);
    out[i] = new_histo;
    bit_costs[i] = Histogram_EstimateBits(out[i]);
  }
  // Collapse similar histograms.
  for (iter = 0; iter < in_size * 3 && *out_size >= 2; ++iter) {
    double best_val = 0;
    int best_ix0 = 0;
    int best_ix1 = 0;
    // Try a few times.
    int k;
    for (k = 0; k < inner_iters; ++k) {
      // Choose two, build a combo out of them.
      double cost_val;
      Histogram *combo;
      int ix0 = rand_r(&seed) % *out_size;
      int ix1;
      int diff = ((k & 7) + 1) % (*out_size - 1);
      if (diff >= 3) {
        diff = rand_r(&seed) % (*out_size - 1);
      }
      ix1 = (ix0 + diff + 1) % *out_size;
      if (ix0 == ix1) {
        continue;
      }
      combo = (Histogram *)malloc(sizeof(*combo));
      if (combo == NULL) {
        ok = 0;
        goto exit_label;
      }
      Histogram_Init(combo, palettebits);
      *combo = *out[ix0];
      Histogram_Add(combo, out[ix1]);
      cost_val =
          Histogram_EstimateBits(combo) - bit_costs[ix0] - bit_costs[ix1];
      if (best_val > cost_val) {
        best_val = cost_val;
        best_ix0 = ix0;
        best_ix1 = ix1;
      }
      free(combo);
    }
    if (best_val < 0.0) {
      Histogram_Add(out[best_ix0], out[best_ix1]);
      bit_costs[best_ix0] =
          best_val + bit_costs[best_ix0] + bit_costs[best_ix1];
      // Erase (*out)[best_ix1]
      free(out[best_ix1]);
      memmove(&out[best_ix1], &out[best_ix1 + 1],
              (*out_size - best_ix1 - 1) * sizeof(out[0]));
      memmove(&bit_costs[best_ix1], &bit_costs[best_ix1 + 1],
              (*out_size - best_ix1 - 1) * sizeof(bit_costs[0]));
      --(*out_size);
      tries_with_no_success = 0;
    }
    if (++tries_with_no_success >= 50) {
      break;
    }
  }
 exit_label:
  free(bit_costs);
  if (!ok) {
    if (out) {
      int i;
      for (i = 0; i < *out_size; ++i) {
        free(out[i]);
      }
      free(out);
    }
  }
  return ok;
}

// What is the bit cost of moving square_histogram from
// cur_symbol to candidate_symbol.
static double HistogramDistance(const Histogram * const square_histogram,
                                int cur_symbol,
                                int candidate_symbol,
                                Histogram **candidate_histograms) {
  double new_bit_cost;
  double previous_bit_cost;
  Histogram modified;
  if (cur_symbol == candidate_symbol) {
    return 0;  // Going nowhere. No savings.
  }
  previous_bit_cost =
      Histogram_EstimateBits(candidate_histograms[candidate_symbol]);
  if (cur_symbol != -1) {
    previous_bit_cost +=
        Histogram_EstimateBits(candidate_histograms[cur_symbol]);
  }

  Histogram_Init(&modified, square_histogram->palette_code_bits_);
  // Compute the bit cost of the histogram where the data moves to.
  modified = *candidate_histograms[candidate_symbol];
  Histogram_Add(&modified, square_histogram);
  new_bit_cost = Histogram_EstimateBits(&modified);

  // Compute the bit cost of the histogram where the data moves away.
  if (cur_symbol != -1) {
    modified = *candidate_histograms[cur_symbol];
    Histogram_Remove(&modified, square_histogram);
    new_bit_cost += Histogram_EstimateBits(&modified);
  }
  return new_bit_cost - previous_bit_cost;
}

void RefineHistogramImage(Histogram **raw,
                          int raw_size,
                          uint32_t *symbols,
                          int out_size,
                          Histogram **out) {
  int i;
  // Find the best 'out' histogram for each of the raw histograms
  for (i = 0; i < raw_size; ++i) {
    int best_out = 0;
    double best_bits = HistogramDistance(raw[i], symbols[i], 0, out);
    int k;
    for (k = 1; k < out_size; ++k) {
      double cur_bits = HistogramDistance(raw[i], symbols[i], k, out);
      if (cur_bits < best_bits) {
        best_bits = cur_bits;
        best_out = k;
      }
    }
    symbols[i] = best_out;
  }

  // Recompute each out based on raw and symbols.
  for (i = 0; i < out_size; ++i) {
    Histogram_Clear(out[i]);
  }
  for (i = 0; i < raw_size; ++i) {
    Histogram_Add(out[symbols[i]], raw[i]);
  }
}
