// Copyright 2011 Google Inc.
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
                         const std::vector<LiteralOrCopy> &backward_refs,
                         std::vector<Histogram *> *image) {
  int histo_xsize = histobits ? (xsize + (1 << histobits) - 1) >> histobits : 1;
  int histo_ysize = histobits ? (ysize + (1 << histobits) - 1) >> histobits : 1;
  image->resize(histo_xsize * histo_ysize);
  for (int i = 0; i < image->size(); ++i) {
    (*image)[i] = new Histogram(palettebits);
  }
  // x and y trace the position in the image.
  int x = 0;
  int y = 0;
  for (int i = 0; i < backward_refs.size(); ++i) {
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

void CombineHistogramImage(const std::vector<Histogram *> &in,
                           int quality,
                           int palettebits,
                           std::vector<Histogram *> *out) {
  // Copy
  out->resize(in.size());
  std::vector<double> bit_costs(in.size());
  for (int i = 0; i < in.size(); ++i) {
    Histogram *new_histo = new Histogram(palettebits);
    *new_histo = *(in[i]);
    (*out)[i] = new_histo;
    bit_costs[i] = (*out)[i]->EstimateBits();
  }
  // Collapse similar histograms.
  unsigned int seed = 0;
  int tries_with_no_success = 0;
  int inner_iters = 10 + quality / 2;
  for (int iter = 0; iter < in.size() * 3 && out->size() >= 2; ++iter) {
    double best_val = 0;
    int best_ix0 = 0;
    int best_ix1 = 0;
    // Try a few times.
    for (int k = 0; k < inner_iters; ++k) {
      // Choose two.
      int ix0 = rand_r(&seed) % out->size();
      int diff = ((k & 7) + 1) % (out->size() - 1);
      if (diff >= 3) {
        diff = rand_r(&seed) % (out->size() - 1);
      }
      const int ix1 = (ix0 + diff + 1) % out->size();
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
      delete *(out->begin() + best_ix1);
      out->erase(out->begin() + best_ix1);
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
  if (cur_symbol == candidate_symbol) {
    return 0;  // going nowhere. no savings.
  }
  double previous_bit_cost =
      (*candidate_histograms)[candidate_symbol]->EstimateBits();
  if (cur_symbol != -1) {
    previous_bit_cost += (*candidate_histograms)[cur_symbol]->EstimateBits();
  }

  Histogram *tmp = new Histogram(square_histogram.palette_code_bits_);
  // Compute the bit cost of the histogram where the data moves to.
  *tmp = *(*candidate_histograms)[candidate_symbol];
  tmp->Add(square_histogram);
  double new_bit_cost = tmp->EstimateBits();

  // Compute the bit cost of the histogram where the data moves away.
  if (cur_symbol != -1) {
    *tmp = *(*candidate_histograms)[cur_symbol];
    tmp->Remove(square_histogram);
    new_bit_cost += tmp->EstimateBits();
  }
  delete tmp;
  return new_bit_cost - previous_bit_cost;
}

void RefineHistogramImage(const std::vector<Histogram *> &raw,
                          std::vector<uint32> *symbols,
                          std::vector<Histogram *> *out) {
  symbols->resize(raw.size());

  // Find the best 'out' histogram for each of the raw histograms
  for (int i = 0; i < raw.size(); ++i) {
    int best_out = 0;
    double best_bits = HistogramDistance(*raw[i], (*symbols)[i], 0, out);
    for (int k = 1; k < out->size(); ++k) {
      double cur_bits = HistogramDistance(*raw[i], (*symbols)[i], k, out);
      if (cur_bits < best_bits) {
        best_bits = cur_bits;
        best_out = k;
      }
    }
    (*symbols)[i] = best_out;
  }

  // Recompute each out based on raw and symbols.
  for (int i = 0; i < out->size(); ++i) {
    (*out)[i]->Clear();
  }
  for (int i = 0; i < raw.size(); ++i) {
    (*out)[(*symbols)[i]]->Add(*raw[i]);
  }
}
