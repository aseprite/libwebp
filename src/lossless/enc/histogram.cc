// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)

#include "histogram.h"

#include <math.h>
#include <stdio.h>

#include "backward_distance.h"
#include "backward_references.h"
#include "../common/integral_types.h"

namespace {

// A lookup table for small values of log(int) to be used in entropy
// computation.
//
// ", ".join(["%.16ff" % x for x in [0.0]+[log(x) for x in range(1, 256)]])
static const double kLogTable[] = {
  0.0000000000000000, 0.0000000000000000, 0.6931471805599453,
  1.0986122886681098, 1.3862943611198906, 1.6094379124341003,
  1.7917594692280550, 1.9459101490553132, 2.0794415416798357,
  2.1972245773362196, 2.3025850929940459, 2.3978952727983707,
  2.4849066497880004, 2.5649493574615367, 2.6390573296152584,
  2.7080502011022101, 2.7725887222397811, 2.8332133440562162,
  2.8903717578961645, 2.9444389791664403, 2.9957322735539909,
  3.0445224377234230, 3.0910424533583161, 3.1354942159291497,
  3.1780538303479458, 3.2188758248682006, 3.2580965380214821,
  3.2958368660043291, 3.3322045101752038, 3.3672958299864741,
  3.4011973816621555, 3.4339872044851463, 3.4657359027997265,
  3.4965075614664802, 3.5263605246161616, 3.5553480614894135,
  3.5835189384561099, 3.6109179126442243, 3.6375861597263857,
  3.6635616461296463, 3.6888794541139363, 3.7135720667043080,
  3.7376696182833684, 3.7612001156935624, 3.7841896339182610,
  3.8066624897703196, 3.8286413964890951, 3.8501476017100584,
  3.8712010109078911, 3.8918202981106265, 3.9120230054281460,
  3.9318256327243257, 3.9512437185814275, 3.9702919135521220,
  3.9889840465642745, 4.0073331852324712, 4.0253516907351496,
  4.0430512678345503, 4.0604430105464191, 4.0775374439057197,
  4.0943445622221004, 4.1108738641733114, 4.1271343850450917,
  4.1431347263915326, 4.1588830833596715, 4.1743872698956368,
  4.1896547420264252, 4.2046926193909657, 4.2195077051761070,
  4.2341065045972597, 4.2484952420493594, 4.2626798770413155,
  4.2766661190160553, 4.2904594411483910, 4.3040650932041702,
  4.3174881135363101, 4.3307333402863311, 4.3438054218536841,
  4.3567088266895917, 4.3694478524670215, 4.3820266346738812,
  4.3944491546724391, 4.4067192472642533, 4.4188406077965983,
  4.4308167988433134, 4.4426512564903167, 4.4543472962535073,
};

// Faster logarithm for small integers, with the property of log(0) == 0.
inline double FastLog(int v) {
  if (v < sizeof(kLogTable) / sizeof(kLogTable[0])) {
    return kLogTable[v];
  }
  return log(v);
}

}  // namespace

void ConvertPopulationCountTableToBitEstimates(
    int n, const int *population_counts,
    double *output) {
  int sum = 0;
  int nonzeros = 0;
  for (int i = 0; i < n; ++i) {
    sum += population_counts[i];
    if (population_counts[i] != 0) {
      ++nonzeros;
    }
  }
  if (nonzeros <= 1) {
    for (int i = 0; i < n; ++i) {
      output[i] = 0;
    }
    return;
  }
  double log2sum = log2(sum);
  for (int i = 0; i < n; ++i) {
    if (population_counts[i] == 0) {
      output[i] = log2sum;
    } else {
      output[i] = log2sum - log2(population_counts[i]);
    }
  }
}

void Histogram::AddSingleLiteralOrCopy(const LiteralOrCopy &v) {
  if (v.IsLiteral()) {
    ++alpha_[v.Literal(3)];
    ++red_[v.Literal(2)];
    ++literal_[v.Literal(1)];
    ++blue_[v.Literal(0)];
  } else if (v.IsPaletteIx()) {
    int literal_ix = 256 + v.PaletteIx();
    ++literal_[literal_ix];
  } else {
    int code, extra_bits_count, extra_bits_value;
    BackwardLength::Encode(v.Length(),
                           &code,
                           &extra_bits_count,
                           &extra_bits_value);
    ++literal_[256 + (1 << palette_code_bits_) + code];
    BackwardDistance::Encode(v.Distance(),
                             &code,
                             &extra_bits_count,
                             &extra_bits_value);
    ++distance_[code];
  }
}

void Histogram::Build(const LiteralOrCopy *literal_and_length,
                      int n_literal_and_length) {
  Clear();
  for (int i = 0; i < n_literal_and_length; ++i) {
    AddSingleLiteralOrCopy(literal_and_length[i]);
  }
}

double BitsEntropy(const int *array, int n) {
  double retval = 0;
  int sum = 0;
  int nonzeros = 0;
  int max_val = 0;
  for (int i = 0; i < n; ++i) {
    if (array[i] != 0) {
      sum += array[i];
      ++nonzeros;
      retval += array[i] * FastLog(array[i]);
      if (max_val < array[i]) {
        max_val = array[i];
      }
    }
  }
  if (nonzeros <= 1) {
    return 0;
  }
  retval -= sum * FastLog(sum);
  retval *= -1.4426950408889634;  // 1.0 / -FastLog(2);

  // Two symbols, they will be 0 and 1 in a Huffman code.
  // Let's mix in a bit of entropy to favor good clustering when
  // distributions of these are combined.
  if (nonzeros == 2) {
    return 0.99 * sum + 0.01 * retval;
  }

  // No matter what the entropy says, we cannot be better than min_limit
  // with Huffman coding. I am mixing a bit of entropy into the
  // min_limit since it produces much better (~0.5 %) compression results
  // perhaps because of better entropy clustering.
  double mix = 0.627;
  if (nonzeros == 3) {
    mix = 0.95;
  } else if (nonzeros == 4) {
    mix = 0.7;
  }
  double min_limit = 2 * sum - max_val;
  min_limit = mix * min_limit + (1.0 - mix) * retval;
  if (retval < min_limit) {
    return min_limit;
  }
  return retval;
}

double Histogram::EstimateBitsBulk() const {
  double retval = BitsEntropy(&literal_[0], NumLiteralOrCopyCodes()) +
      BitsEntropy(&red_[0], 256) +
      BitsEntropy(&blue_[0], 256) +
      BitsEntropy(&alpha_[0], 256) +
      BitsEntropy(&distance_[0], kDistanceCodes);
  // Compute the extra bits cost.
  for (size_t i = 2; i < kLengthCodes - 2; ++i) {
    retval += (i >> 1) * literal_[256 + (1 << palette_code_bits_) + i + 2];
  }
  for (size_t i = 2; i < kDistanceCodes - 2; ++i) {
    retval += (i >> 1) * distance_[i + 2];
  }
  return retval;
}

double Histogram::EstimateBits() const {
  return EstimateBitsHeader() + EstimateBitsBulk();
}

// Returns the cost encode the rle-encoded entropy code.
// The constants in this function are experimental.
double HuffmanCost(const int *population, int length) {
  // Small bias because Huffman code length is typically not stored in
  // full length.
  static const int kHuffmanCodeOfHuffmanCodeSize = kCodeLengthCodes * 3;
  static const double kSmallBias = 9.1;
  double retval = kHuffmanCodeOfHuffmanCodeSize - kSmallBias;
  int streak = 0;
  int i = 0;
  for (; i < length - 1; ++i) {
    ++streak;
    if (population[i] == population[i + 1]) {
      continue;
    }
 last_streak_hack:
    // population[i] points now to the symbol in the streak of same values.
    if (streak > 3) {
      if (population[i] == 0) {
        retval += 1.5625 + 0.234375 * streak;
      } else {
        retval += 2.578125 + 0.703125 * streak;
      }
    } else {
      if (population[i] == 0) {
        retval += 1.796875 * streak;
      } else {
        retval += 3.28125 * streak;
      }
    }
    streak = 0;
  }
  if (i == length - 1) {
    ++streak;
    goto last_streak_hack;
  }
  return retval;
}

double Histogram::EstimateBitsHeader() const {
  return HuffmanCost(&alpha_[0], 256) +
      HuffmanCost(&red_[0], 256) +
      HuffmanCost(&literal_[0], NumLiteralOrCopyCodes()) +
      HuffmanCost(&blue_[0], 256) +
      HuffmanCost(&distance_[0], kDistanceCodes);
}
