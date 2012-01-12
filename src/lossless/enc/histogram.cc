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
static const float kLogTable[] = {
  0.0000000000000000f, 0.0000000000000000f, 0.6931471805599453f,
  1.0986122886681098f, 1.3862943611198906f, 1.6094379124341003f,
  1.7917594692280550f, 1.9459101490553132f, 2.0794415416798357f,
  2.1972245773362196f, 2.3025850929940459f, 2.3978952727983707f,
  2.4849066497880004f, 2.5649493574615367f, 2.6390573296152584f,
  2.7080502011022101f, 2.7725887222397811f, 2.8332133440562162f,
  2.8903717578961645f, 2.9444389791664403f, 2.9957322735539909f,
  3.0445224377234230f, 3.0910424533583161f, 3.1354942159291497f,
  3.1780538303479458f, 3.2188758248682006f, 3.2580965380214821f,
  3.2958368660043291f, 3.3322045101752038f, 3.3672958299864741f,
  3.4011973816621555f, 3.4339872044851463f, 3.4657359027997265f,
  3.4965075614664802f, 3.5263605246161616f, 3.5553480614894135f,
  3.5835189384561099f, 3.6109179126442243f, 3.6375861597263857f,
  3.6635616461296463f, 3.6888794541139363f, 3.7135720667043080f,
  3.7376696182833684f, 3.7612001156935624f, 3.7841896339182610f,
  3.8066624897703196f, 3.8286413964890951f, 3.8501476017100584f,
  3.8712010109078911f, 3.8918202981106265f, 3.9120230054281460f,
  3.9318256327243257f, 3.9512437185814275f, 3.9702919135521220f,
  3.9889840465642745f, 4.0073331852324712f, 4.0253516907351496f,
  4.0430512678345503f, 4.0604430105464191f, 4.0775374439057197f,
  4.0943445622221004f, 4.1108738641733114f, 4.1271343850450917f,
  4.1431347263915326f, 4.1588830833596715f, 4.1743872698956368f,
  4.1896547420264252f, 4.2046926193909657f, 4.2195077051761070f,
  4.2341065045972597f, 4.2484952420493594f, 4.2626798770413155f,
  4.2766661190160553f, 4.2904594411483910f, 4.3040650932041702f,
  4.3174881135363101f, 4.3307333402863311f, 4.3438054218536841f,
  4.3567088266895917f, 4.3694478524670215f, 4.3820266346738812f,
  4.3944491546724391f, 4.4067192472642533f, 4.4188406077965983f,
  4.4308167988433134f, 4.4426512564903167f, 4.4543472962535073f,
  4.4659081186545837f, 4.4773368144782069f, 4.4886363697321396f,
  4.4998096703302650f, 4.5108595065168497f, 4.5217885770490405f,
  4.5325994931532563f, 4.5432947822700038f, 4.5538768916005408f,
  4.5643481914678361f, 4.5747109785033828f, 4.5849674786705723f,
  4.5951198501345898f, 4.6051701859880918f, 4.6151205168412597f,
  4.6249728132842707f, 4.6347289882296359f, 4.6443908991413725f,
  4.6539603501575231f, 4.6634390941120669f, 4.6728288344619058f,
  4.6821312271242199f, 4.6913478822291435f, 4.7004803657924166f,
  4.7095302013123339f, 4.7184988712950942f, 4.7273878187123408f,
  4.7361984483944957f, 4.7449321283632502f, 4.7535901911063645f,
  4.7621739347977563f, 4.7706846244656651f, 4.7791234931115296f,
  4.7874917427820458f, 4.7957905455967413f, 4.8040210447332568f,
  4.8121843553724171f, 4.8202815656050371f, 4.8283137373023015f,
  4.8362819069514780f, 4.8441870864585912f, 4.8520302639196169f,
  4.8598124043616719f, 4.8675344504555822f, 4.8751973232011512f,
  4.8828019225863706f, 4.8903491282217537f, 4.8978397999509111f,
  4.9052747784384296f, 4.9126548857360524f, 4.9199809258281251f,
  4.9272536851572051f, 4.9344739331306915f, 4.9416424226093039f,
  4.9487598903781684f, 4.9558270576012609f, 4.9628446302599070f,
  4.9698132995760007f, 4.9767337424205742f, 4.9836066217083363f,
  4.9904325867787360f, 4.9972122737641147f, 5.0039463059454592f,
  5.0106352940962555f, 5.0172798368149243f, 5.0238805208462765f,
  5.0304379213924353f, 5.0369526024136295f, 5.0434251169192468f,
  5.0498560072495371f, 5.0562458053483077f, 5.0625950330269669f,
  5.0689042022202315f, 5.0751738152338266f, 5.0814043649844631f,
  5.0875963352323836f, 5.0937502008067623f, 5.0998664278241987f,
  5.1059454739005803f, 5.1119877883565437f, 5.1179938124167554f,
  5.1239639794032588f, 5.1298987149230735f, 5.1357984370502621f,
  5.1416635565026603f, 5.1474944768134527f, 5.1532915944977793f,
  5.1590552992145291f, 5.1647859739235145f, 5.1704839950381514f,
  5.1761497325738288f, 5.1817835502920850f, 5.1873858058407549f,
  5.1929568508902104f, 5.1984970312658261f, 5.2040066870767951f,
  5.2094861528414214f, 5.2149357576089859f, 5.2203558250783244f,
  5.2257466737132017f, 5.2311086168545868f, 5.2364419628299492f,
  5.2417470150596426f, 5.2470240721604862f, 5.2522734280466299f,
  5.2574953720277815f, 5.2626901889048856f, 5.2678581590633282f,
  5.2729995585637468f, 5.2781146592305168f, 5.2832037287379885f,
  5.2882670306945352f, 5.2933048247244923f, 5.2983173665480363f,
  5.3033049080590757f, 5.3082676974012051f, 5.3132059790417872f,
  5.3181199938442161f, 5.3230099791384085f, 5.3278761687895813f,
  5.3327187932653688f, 5.3375380797013179f, 5.3423342519648109f,
  5.3471075307174685f, 5.3518581334760666f, 5.3565862746720123f,
  5.3612921657094255f, 5.3659760150218512f, 5.3706380281276624f,
  5.3752784076841653f, 5.3798973535404597f, 5.3844950627890888f,
  5.3890717298165010f, 5.3936275463523620f, 5.3981627015177525f,
  5.4026773818722793f, 5.4071717714601188f, 5.4116460518550396f,
  5.4161004022044201f, 5.4205349992722862f, 5.4249500174814029f,
  5.4293456289544411f, 5.4337220035542400f, 5.4380793089231956f,
  5.4424177105217932f, 5.4467373716663099f, 5.4510384535657002f,
  5.4553211153577017f, 5.4595855141441589f, 5.4638318050256105f,
  5.4680601411351315f, 5.4722706736714750f, 5.4764635519315110f,
  5.4806389233419912f, 5.4847969334906548f, 5.4889377261566867f,
  5.4930614433405482f, 5.4971682252932021f, 5.5012582105447274f,
  5.5053315359323625f, 5.5093883366279774f, 5.5134287461649825f,
  5.5174528964647074f, 5.5214609178622460f, 5.5254529391317835f,
  5.5294290875114234f, 5.5333894887275203f, 5.5373342670185366f,
  5.5412635451584258f,
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
  retval -= sum * FastLog(sum);
  retval *= -1.4426950408889634;  // 1.0 / -FastLog(2);
  double mix = 0.627;
  if (nonzeros < 5) {
    if (nonzeros <= 1) {
      return 0;
    }
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
    if (nonzeros == 3) {
      mix = 0.95;
    } else {
      mix = 0.7;  // nonzeros == 4.
    }
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
