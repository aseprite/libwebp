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
#include <vector>

#include "histogram.h"
#include "backward_references.h"
#include "../common/integral_types.h"

void BuildHistogramImage(int xsize, int ysize,
                         int histobits,
                         int palettebits,
                         const LiteralOrCopy *backward_refs,
                         int backward_refs_size,
                         Histogram ***image,
                         int *image_size) {
  int histo_xsize = histobits ? (xsize + (1 << histobits) - 1) >> histobits : 1;
  int histo_ysize = histobits ? (ysize + (1 << histobits) - 1) >> histobits : 1;
  int i;
  *image_size = histo_xsize * histo_ysize;
  *image = (Histogram **)malloc(*image_size * sizeof(Histogram *));
  for (i = 0; i < *image_size; ++i) {
    (*image)[i] = new Histogram(palettebits);
  }
  // x and y trace the position in the image.
  int x = 0;
  int y = 0;
  for (i = 0; i < backward_refs_size; ++i) {
    const LiteralOrCopy &v = backward_refs[i];
    const int ix =
        histobits ? (y >> histobits) * histo_xsize + (x >> histobits) : 0;
    (*image)[ix]->AddSingleLiteralOrCopy(v);
    x += v.Length();
    while (x >= xsize) {
      x -= xsize;
      ++y;
    }
  }
}

void CombineHistogramImage(Histogram **in,
                           int in_size,
                           int quality,
                           int palettebits,
                           Histogram ***out,
                           int *out_size) {
  int i;
  unsigned int seed = 0;
  int tries_with_no_success = 0;
  int inner_iters = 10 + quality / 2;
  int iter;
  // Copy
  *out_size = in_size;
  *out = (Histogram **)malloc(in_size * sizeof(Histogram *));
  std::vector<double> bit_costs(in_size);
  for (i = 0; i < in_size; ++i) {
    Histogram *new_histo = new Histogram(palettebits);
    *new_histo = *(in[i]);
    (*out)[i] = new_histo;
    bit_costs[i] = (*out)[i]->EstimateBits();
  }
  // Collapse similar histograms.
  for (iter = 0; iter < in_size * 3 && *out_size >= 2; ++iter) {
    double best_val = 0;
    int best_ix0 = 0;
    int best_ix1 = 0;
    // Try a few times.
    int k;
    for (k = 0; k < inner_iters; ++k) {
      // Choose two.
      int ix0 = rand_r(&seed) % *out_size;
      int diff = ((k & 7) + 1) % (*out_size - 1);
      if (diff >= 3) {
        diff = rand_r(&seed) % (*out_size - 1);
      }
      const int ix1 = (ix0 + diff + 1) % *out_size;
      if (ix0 == ix1) {
        continue;
      }
      Histogram *combo = new Histogram(palettebits);
      *combo = *(*out)[ix0];
      combo->Add(*(*out)[ix1]);
      const double val = combo->EstimateBits() -
          bit_costs[ix0] - bit_costs[ix1];
      if (best_val > val) {
        best_val = val;
        best_ix0 = ix0;
        best_ix1 = ix1;
      }
      delete combo;
    }
    if (best_val < 0.0) {
      (*out)[best_ix0]->Add(*((*out)[best_ix1]));
      bit_costs[best_ix0] =
          best_val + bit_costs[best_ix0] + bit_costs[best_ix1];
      // Erase (*out)[best_ix1]
      delete (*out)[best_ix1];
      memmove(&(*out)[best_ix1], &(*out)[best_ix1 + 1],
              (*out_size - best_ix1 - 1) * sizeof((*out)[0]));
      --(*out_size);
      bit_costs.erase(bit_costs.begin() + best_ix1);
      tries_with_no_success = 0;
    }
    if (++tries_with_no_success >= 50) {
      break;
    }
  }
}

// What is the bit cost of moving square_histogram from
// cur_symbol to candidate_symbol.
double HistogramDistance(const Histogram &square_histogram,
                         int cur_symbol,
                         int candidate_symbol,
                         std::vector<Histogram *> *candidate_histograms) {
  double new_bit_cost;
  double previous_bit_cost;
  if (cur_symbol == candidate_symbol) {
    return 0;  // going nowhere. no savings.
  }
  previous_bit_cost =
      (*candidate_histograms)[candidate_symbol]->EstimateBits();
  if (cur_symbol != -1) {
    previous_bit_cost += (*candidate_histograms)[cur_symbol]->EstimateBits();
  }

  Histogram *tmp = new Histogram(square_histogram.palette_code_bits_);
  // Compute the bit cost of the histogram where the data moves to.
  *tmp = *(*candidate_histograms)[candidate_symbol];
  tmp->Add(square_histogram);
  new_bit_cost = tmp->EstimateBits();

  // Compute the bit cost of the histogram where the data moves away.
  if (cur_symbol != -1) {
    *tmp = *(*candidate_histograms)[cur_symbol];
    tmp->Remove(square_histogram);
    new_bit_cost += tmp->EstimateBits();
  }
  delete tmp;
  return new_bit_cost - previous_bit_cost;
}

void RefineHistogramImage(Histogram **raw,
                          int raw_size,
                          std::vector<uint32_t> *symbols,
                          std::vector<Histogram *> *out) {
  int i;
  symbols->resize(raw_size);

  // Find the best 'out' histogram for each of the raw histograms
  for (i = 0; i < raw_size; ++i) {
    int best_out = 0;
    double best_bits = HistogramDistance(*raw[i], (*symbols)[i], 0, out);
    int k;
    for (k = 1; k < out->size(); ++k) {
      double cur_bits = HistogramDistance(*raw[i], (*symbols)[i], k, out);
      if (cur_bits < best_bits) {
        best_bits = cur_bits;
        best_out = k;
      }
    }
    (*symbols)[i] = best_out;
  }

  // Recompute each out based on raw and symbols.
  for (i = 0; i < out->size(); ++i) {
    (*out)[i]->Clear();
  }
  for (i = 0; i < raw_size; ++i) {
    (*out)[(*symbols)[i]]->Add(*raw[i]);
  }
}
