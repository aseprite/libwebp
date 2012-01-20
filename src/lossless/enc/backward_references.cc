// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)

#include <math.h>
#include <stdio.h>

#include <algorithm>
#include <map>
#include <set>

#include "backward_references.h"
#include "histogram.h"
#include "../common/integral_types.h"
#include "../common/pixel_hasher.h"

static const unsigned char plane_to_code_lut[128] = {
 96,   73,  55,  39,  23,  13,   5,  1,  255, 255, 255, 255, 255, 255, 255, 255,
 101,  78,  58,  42,  26,  16,   8,  2,    0,   3,  9,   17,  27,  43,  59,  79,
 102,  86,  62,  46,  32,  20,  10,  6,    4,   7,  11,  21,  33,  47,  63,  87,
 105,  90,  70,  52,  37,  28,  18,  14,  12,  15,  19,  29,  38,  53,  71,  91,
 110,  99,  82,  66,  48,  35,  30,  24,  22,  25,  31,  36,  49,  67,  83, 100,
 115, 108,  94,  76,  64,  50,  44,  40,  34,  41,  45,  51,  65,  77,  95, 109,
 118, 113, 103,  92,  80,  68,  60,  56,  54,  57,  61,  69,  81,  93, 104, 114,
 119, 116, 111, 106,  97,  88,  84,  74,  72,  75,  85,  89,  98, 107, 112, 117,
};

static const int kMinLength = 2;

int DistanceToPlaneCode(int xsize, int ysize, int dist) {
  int yoffset = dist / xsize;
  int xoffset = dist - yoffset * xsize;
  if (xoffset <= 8 && yoffset < 8) {
    return plane_to_code_lut[yoffset * 16 + 8 - xoffset] + 1;
  } else if (xoffset > xsize - 8 && yoffset < 7) {
    return plane_to_code_lut[(yoffset + 1) * 16 + 8 + (xsize - xoffset)] + 1;
  }
  return dist + 120;
}

static inline int FindMatchLength(const uint32* array1,
                                  const uint32* array2,
                                  const int max_limit) {
  int matched = 0;
  while (matched < max_limit && array1[matched] == array2[matched]) {
    ++matched;
  }
  return matched;
}

class HashChain {
 public:
  explicit HashChain(int size) {
    chain_ = (int *)malloc(size * sizeof(*chain_));
    for (int i = 0; i < size; ++i) {
      chain_[i] = -1;
    }
    for (int i = 0; i < kHashSize; ++i) {
      hash_to_first_index[i] = -1;
    }
  }

  ~HashChain() {
    free(chain_);
  }

  void Insert(const uint32* argb, int32 ix) {
    // Insertion of two pixels at a time.
    const uint64 key = GetPixPair(argb);
    const uint64 hash_code = GetHash64(key);
    chain_[ix] = hash_to_first_index[hash_code];
    hash_to_first_index[hash_code] = ix;
  }

  bool FindCopy(int index, int xsize, const uint32 * __restrict argb,
                int maxlen, int * __restrict offset, int * __restrict len) {
    const uint64 next_two_pixels = GetPixPair(&argb[index]);
    const uint64 hash_code = GetHash64(next_two_pixels);
    *len = 0;
    *offset = 0;
    int prev_length = 0;
    int64 best_val = 0;
    int give_up = 0;
    const int min_pos = (index > kWindowSize) ? index - kWindowSize : 0;
    for (int32 pos = hash_to_first_index[hash_code];
         pos >= min_pos;
         pos = chain_[pos]) {
      ++give_up;
      if (give_up >= 101) {
        if (give_up >= 1001 ||
            best_val >= 0xff0000) {
          break;
        }
      }
      if (*len != 0 && argb[pos + *len - 1] != argb[index + *len - 1]) {
        continue;
      }
      int64 length = FindMatchLength(argb + pos, argb + index, maxlen);
      int64 val = 65536 * length;
      int y = (index - pos) / xsize;
      int x = (index - pos) % xsize;
      // Favoring 2d locality here gives savings for certain images.
      if (y < 8) {
        if (x > xsize / 2) {
          x = xsize - x;
        }
        if (x <= 7 && x >= -8) {
          val -= y * y + x * x;
        } else {
          val -= 9 * 9 + 9 * 9;
        }
      } else {
        val -= 9 * 9 + 9 * 9;
      }
      if (best_val < val) {
        prev_length = length;
        best_val = val;
        *len = length;
        *offset = index - pos;
        if (length >= kMaxLength) {
          return true;
        }
        if ((*offset == 1 || *offset == xsize) && *len >= 128) {
          return true;
        }
      }
    }
    return *len >= kMinLength;
  }

 private:
  static inline uint64 GetHash64(uint64 num) {
    num *= kHashMultiplier;
    num >>= 64 - kHashBits;
    return num;
  }

  static inline uint64 GetPixPair(const uint32 *argb) {
    return (static_cast<uint64>(argb[1]) << 32) | argb[0];
  }

  static const int kHashBits = 20;
  static const int kHashSize = (1 << kHashBits);
  static const uint64 kHashMultiplier = 0xc6a4a7935bd1e995ULL;
  static const int kWindowSize = (1 << 20) - 120;  // A window with 1M pixels
                                                   // (4 megabytes) - 120
                                                   // special codes for short
                                                   // distances.

  int32 hash_to_first_index[kHashSize];  // Stores the most recently added
                                         // position with the given hash value.
  int32 *chain_;  // chain_[pos] stores the previous position with the same
                  // hash value.
};

static inline void PushBackCopy(int distance, int length,
                                std::vector<LiteralOrCopy> *stream) {
  while (length >= kMaxLength) {
    stream->push_back(LiteralOrCopy::CreateCopy(1, kMaxLength));
    length -= kMaxLength;
  }
  if (length > 0) {
    stream->push_back(LiteralOrCopy::CreateCopy(1, length));
  }
}

void BackwardReferencesRle(int xsize, int ysize, const uint32 *argb,
                           std::vector<LiteralOrCopy> *stream) {
  const int pix_count = xsize * ysize;
  stream->reserve(pix_count);
  int streak = 0;
  for (int i = 0; i < pix_count; ++i) {
    if (i >= 1 && argb[i] == argb[i - 1]) {
      ++streak;
    } else {
      PushBackCopy(1, streak, stream);
      streak = 0;
      stream->push_back(LiteralOrCopy::CreateLiteral(argb[i]));
    }
  }
  PushBackCopy(1, streak, stream);
}

void BackwardReferencesHashChain(int xsize, int ysize, bool use_palette,
                                 const uint32 *argb, int palette_bits,
                                 std::vector<LiteralOrCopy> *stream) {
  const int pix_count = xsize * ysize;
  PixelHasherLine hashers(xsize, kRowHasherXSubsampling, palette_bits);
  HashChain *hash_chain = new HashChain(pix_count);
  for (int i = 0; i < pix_count; ) {
    // Alternative#1: Code the pixels starting at 'i' using backward reference.
    int offset = 0;
    int len = 0;
    const int x = i % xsize;
    if (i < pix_count - 1) {  // FindCopy(i,..) reads pixels at [i] and [i + 1].
      const int maxlen = std::min(pix_count - i, kMaxLength);
      hash_chain->FindCopy(i, xsize, argb, maxlen, &offset, &len);
    }
    if (len >= kMinLength) {
      hash_chain->Insert(&argb[i], i);
      // Alternative#2: Insert the pixel at 'i' as literal, and code the
      // pixels starting at 'i + 1' using backward reference.
      int offset2 = 0;
      int len2 = 0;
      if (i < pix_count - 2) {  // FindCopy(i+1,..) reads [i + 1] and [i + 2].
        const int maxlen = std::min(pix_count - (i + 1), kMaxLength);
        hash_chain->FindCopy(i + 1, xsize, argb, maxlen, &offset2,
                             &len2);
        if (len2 > len + 1) {
          // Alternative#2 is a better match. So push pixel at 'i' as literal.
          if (use_palette && hashers.Contains(x, argb[i])) {
            const int ix = hashers.GetIndex(argb[i]);
            stream->push_back(LiteralOrCopy::CreatePaletteIx(ix));
          } else {
            stream->push_back(LiteralOrCopy::CreateLiteral(argb[i]));
          }
          hashers.Insert(x, argb[i]);
          i++;  // Backward reference to be done for next pixel.
          len = len2;
          offset = offset2;
        }
      }
      if (len >= kMaxLength) {
        len = kMaxLength - 1;
      }
      stream->push_back(LiteralOrCopy::CreateCopy(offset, len));
      for (int k = 0; k < len; ++k) {
        hashers.Insert((i + k) % xsize, argb[i + k]);
        if (k != 0 && i + k + 1 < pix_count) {
          // Add to the hash_chain (but cannot add the last pixel).
          hash_chain->Insert(&argb[i + k], i + k);
        }
      }
      i += len;
    } else {
      if (use_palette && hashers.Contains(x, argb[i])) {
        // push pixel as a palette pixel
        int ix = hashers.GetIndex(argb[i]);
        stream->push_back(LiteralOrCopy::CreatePaletteIx(ix));
      } else {
        stream->push_back(LiteralOrCopy::CreateLiteral(argb[i]));
      }
      hashers.Insert(x, argb[i]);
      if (i + 1 < pix_count) {
        hash_chain->Insert(&argb[i], i);
      }
      ++i;
    }
  }
  delete hash_chain;
}

class CostModel {
 public:
  void Build(int xsize, int ysize, int recursion_level, int use_palette,
             const uint32 *argb, int palette_bits) {
    palette_bits_ = palette_bits;

    std::vector<LiteralOrCopy> stream;
    if (recursion_level > 0) {
      BackwardReferencesTraceBackwards(xsize, ysize, recursion_level - 1,
                                       use_palette, argb,
                                       palette_bits, &stream);
    } else {
      BackwardReferencesHashChain(xsize, ysize, use_palette, argb,
                                  palette_bits, &stream);
    }
    Histogram histo(palette_bits);
    for (int i = 0; i < stream.size(); ++i) {
      histo.AddSingleLiteralOrCopy(stream[i]);
    }
    ConvertPopulationCountTableToBitEstimates(
        histo.NumLiteralOrCopyCodes(),
        &histo.literal_[0], &literal_[0]);
    ConvertPopulationCountTableToBitEstimates(
        kNumSymbols, &histo.red_[0], &red_[0]);
    ConvertPopulationCountTableToBitEstimates(
        kNumSymbols, &histo.blue_[0], &blue_[0]);
    ConvertPopulationCountTableToBitEstimates(
        kNumSymbols, &histo.alpha_[0], &alpha_[0]);
    ConvertPopulationCountTableToBitEstimates(
        kDistanceCodes, &histo.distance_[0], &distance_[0]);
  }

  double LiteralCost(uint32 v) const {
    return alpha_[v >> 24] +
        red_[(v >> 16) & 0xff] +
        literal_[(v >> 8) & 0xff] +
        blue_[v & 0xff];
  }

  double PaletteCost(uint32 ix) const {
    int literal_ix = kNumSymbols + ix;
    return literal_[literal_ix];
  }

  double LengthCost(uint32 len) const {
    int code, extra_bits_count, extra_bits_value;
    BackwardLength::Encode(len,
                           &code,
                           &extra_bits_count,
                           &extra_bits_value);
    return literal_[kNumSymbols +
           (1 << palette_bits_) + code] +
           extra_bits_count;
  }

  double DistanceCost(uint32 distance) const {
    int code, extra_bits_count, extra_bits_value;
    BackwardDistance::Encode(distance,
                             &code,
                             &extra_bits_count,
                             &extra_bits_value);
    return distance_[code] + extra_bits_count;
  }

 private:
  static const int kNumSymbols = 256;
  double alpha_[kNumSymbols];
  double red_[kNumSymbols];
  double literal_[kNumSymbols + (1 << kPaletteCodeBitsMax) + kLengthCodes];
  double blue_[kNumSymbols];
  double distance_[kDistanceCodes];
  int palette_bits_;
};

void BackwardReferencesHashChainDistanceOnly(int xsize, int ysize,
                                             int recursive_cost_model,
                                             bool use_palette,
                                             const uint32 *argb,
                                             int palette_bits,
                                             std::vector<uint32> *dist_array) {
  const int pix_count = xsize * ysize;
  dist_array->resize(pix_count);
  std::vector<double> cost(pix_count, 1e100);
  CostModel *cost_model = new CostModel;
  cost_model->Build(xsize, ysize, recursive_cost_model, use_palette, argb,
                    palette_bits);

  PixelHasherLine hashers(xsize, kRowHasherXSubsampling, palette_bits);
  HashChain *hash_chain = new HashChain(pix_count);
  // We loop one pixel at a time, but store all currently best points to
  // non-processed locations from this point.
  (*dist_array)[0] = 0;
  for (int i = 0; i < pix_count; ++i) {
    double prev_cost = 0.0;
    if (i > 0) {
      prev_cost = cost[i - 1];
    }
    for (int shortmax = 0; shortmax < 2; ++shortmax) {
      int offset = 0;
      int len = 0;
      if (i < pix_count - 1) {  // FindCopy reads pixels at [i] and [i + 1].
        const int maxlen = std::min(pix_count - i, shortmax ? 2 : kMaxLength);
        hash_chain->FindCopy(i, xsize, argb, maxlen, &offset, &len);
      }
      if (len >= kMinLength) {
        const int code = DistanceToPlaneCode(xsize, ysize, offset);
        const double distance_cost = prev_cost + cost_model->DistanceCost(code);
        for (int k = 1; k < len; ++k) {
          const double cost_val = distance_cost + cost_model->LengthCost(k);
          if (cost[i + k] > cost_val) {
            cost[i + k] = cost_val;
            (*dist_array)[i + k] = k + 1;
          }
        }
        // This if is for speedup only. It roughly doubles the speed, and
        // makes compression worse by .1 %.
        if (len >= 128 && code < 2) {
          // Long copy for short distances, let's skip the middle
          // lookups for better copies.
          // 1) insert the hashes.
          for (int k = 0; k < len; ++k) {
            hashers.Insert((i + k) % xsize, argb[i + k]);
            if (i + k + 1 < pix_count) {
              // Add to the hash_chain (but cannot add the last pixel).
              hash_chain->Insert(&argb[i + k], i + k);
            }
          }
          // 2) jump.
          i += len - 1;  // for loop does ++i, thus -1 here.
          goto next_symbol;
        }
      }
    }
    if (i < pix_count - 1) {
      hash_chain->Insert(&argb[i], i);
    }
    {
      // inserting a literal pixel
      double cost_val = prev_cost;
      const int x = i % xsize;
      double mul0 = 1.0;
      double mul1 = 1.0;
      if (recursive_cost_model == 0) {
        mul0 = 0.68;
        mul1 = 0.82;
      }
      if (use_palette && hashers.Contains(x, argb[i])) {
        int ix = hashers.GetIndex(argb[i]);
        cost_val += cost_model->PaletteCost(ix) * mul0;
      } else {
        cost_val += cost_model->LiteralCost(argb[i]) * mul1;
      }
      if (cost[i] > cost_val) {
        cost[i] = cost_val;
        (*dist_array)[i] = 1;  // only one is inserted.
      }
      hashers.Insert(x, argb[i]);
    }
 next_symbol: ;
  }
  // Last pixel still to do, it can only be a single step if not reached
  // through cheaper means already.
  delete hash_chain;
  delete cost_model;
}

void TraceBackwards(const std::vector<uint32> &dist_array,
                    std::vector<uint32> *chosen_path) {
  for (int i = dist_array.size() - 1; i >= 0; ) {
    int k = dist_array[i];
    assert(k >= 1);
    chosen_path->push_back(k);
    i -= k;
  }
  std::reverse(chosen_path->begin(), chosen_path->end());
}

void BackwardReferencesHashChainFollowChosenPath(
    int xsize,
    int ysize,
    bool use_palette,
    const uint32 *argb,
    int palette_bits,
    const std::vector<uint32> &chosen_path,
    std::vector<LiteralOrCopy> *stream) {
  const int pix_count = xsize * ysize;
  PixelHasherLine hashers(xsize, kRowHasherXSubsampling, palette_bits);
  HashChain *hash_chain = new HashChain(pix_count);
  int i = 0;
  for (int ix = 0; ix < chosen_path.size(); ++ix) {
    int offset = 0;
    int len = 0;
    int maxlen = chosen_path[ix];
    if (maxlen != 1) {
      hash_chain->FindCopy(i, xsize, argb, maxlen, &offset, &len);
      assert(len == maxlen);
      stream->push_back(LiteralOrCopy::CreateCopy(offset, len));
      for (int k = 0; k < len; ++k) {
        hashers.Insert((i + k) % xsize, argb[i + k]);
        if (i + k + 1 < pix_count) {
          // Add to the hash_chain (but cannot add the last pixel).
          hash_chain->Insert(&argb[i + k], i + k);
        }
      }
      i += len;
    } else {
      if (use_palette && hashers.Contains(i % xsize, argb[i])) {
        // push pixel as a palette pixel
        int ix = hashers.GetIndex(argb[i]);
        VERIFY(hashers.Lookup(i % xsize, ix) == argb[i]);
        stream->push_back(LiteralOrCopy::CreatePaletteIx(ix));
      } else {
        stream->push_back(LiteralOrCopy::CreateLiteral(argb[i]));
      }
      hashers.Insert(i % xsize, argb[i]);
      if (i + 1 < pix_count) {
        hash_chain->Insert(&argb[i], i);
      }
      ++i;
    }
  }
  delete hash_chain;
}


void BackwardReferencesTraceBackwards(int xsize, int ysize,
                                      int recursive_cost_model,
                                      bool use_palette,
                                      const uint32 *argb,
                                      int palette_bits,
                                      std::vector<LiteralOrCopy> *stream) {
  std::vector<uint32> dist_array;
  BackwardReferencesHashChainDistanceOnly(xsize, ysize, recursive_cost_model,
                                          use_palette, argb, palette_bits,
                                          &dist_array);
  std::vector<uint32> chosen_path;
  TraceBackwards(dist_array, &chosen_path);
  BackwardReferencesHashChainFollowChosenPath(xsize, ysize, use_palette,
                                              argb, palette_bits, chosen_path,
                                              stream);
}

void BackwardReferences2DLocality(int xsize, int ysize, int data_size,
                                  LiteralOrCopy *data) {
  for (int i = 0; i < data_size; ++i) {
    if (data[i].IsCopy()) {
      int dist = data[i].argb_or_offset;
      int transformed_dist = DistanceToPlaneCode(xsize, ysize, dist);
      data[i].argb_or_offset = transformed_dist;
    }
  }
}

bool VerifyBackwardReferences(const uint32* argb, int xsize, int ysize,
                              int palette_bits,
                              const std::vector<LiteralOrCopy>& v) {
  PixelHasherLine hashers(xsize, kRowHasherXSubsampling, palette_bits);
  int num_pixels = 0;
  for (int i = 0; i < v.size(); ++i) {
    if (v[i].IsLiteral()) {
      if (argb[num_pixels] != v[i].Argb()) {
        printf("i %d, pixel %d, original: 0x%08x, literal: 0x%08x\n",
               i, num_pixels, argb[num_pixels], v[i].Argb());
        return false;
      }
      hashers.Insert(num_pixels % xsize, argb[num_pixels]);
      ++num_pixels;
    } else if (v[i].IsPaletteIx()) {
      uint32 palette_entry =
          hashers.Lookup(num_pixels % xsize, v[i].PaletteIx());
      if (argb[num_pixels] != palette_entry) {
        printf("i %d, pixel %d, original: 0x%08x, palette_ix: %d, "
               "palette_entry: 0x%08x\n",
               i, num_pixels, argb[num_pixels], v[i].PaletteIx(),
               palette_entry);
        return false;
      }
      hashers.Insert(num_pixels % xsize, argb[num_pixels]);
      ++num_pixels;
    } else if (v[i].IsCopy()) {
      if (v[i].Distance() == 0) {
        printf("Bw reference with zero distance.\n");
        return false;
      }
      for (int k = 0; k < v[i].len; ++k) {
        if (argb[num_pixels] != argb[num_pixels - v[i].Distance()]) {
          printf("i %d, pixel %d, original: 0x%08x, copied: 0x%08x, dist: %d\n",
                 i, num_pixels, argb[num_pixels],
                 argb[num_pixels - v[i].Distance()], v[i].Distance());
          return false;
        }
        hashers.Insert(num_pixels % xsize, argb[num_pixels]);
        ++num_pixels;
      }
    }
  }
  {
    const int pix_count = xsize * ysize;
    if (num_pixels != pix_count) {
      printf("upsala: %d != %d\n", num_pixels, pix_count);
      return false;
    }
  }
  return true;
}

static void ComputePaletteHistogram(const uint32 *argb, int xsize, int ysize,
                                    const std::vector<LiteralOrCopy> &stream,
                                    int palette_bits, Histogram *histo) {
  PixelHasherLine hashers(xsize, kRowHasherXSubsampling, palette_bits);
  int pixel_index = 0;
  for (int i = 0; i < stream.size(); ++i) {
    const LiteralOrCopy &v = stream[i];
    if (v.IsLiteral()) {
      const int x = pixel_index % xsize;
      if (palette_bits != 0 && hashers.Contains(x, argb[pixel_index])) {
        // push pixel as a palette pixel
        const int ix = hashers.GetIndex(argb[pixel_index]);
        histo->AddSingleLiteralOrCopy(LiteralOrCopy::CreatePaletteIx(ix));
      } else {
        histo->AddSingleLiteralOrCopy(v);
      }
    } else {
      histo->AddSingleLiteralOrCopy(v);
    }
    for (int k = 0; k < v.Length(); ++k) {
      hashers.Insert(pixel_index % xsize, argb[pixel_index]);
      ++pixel_index;
    }
  }
  VERIFY(pixel_index == xsize * ysize);
}

// Returns how many bits are to be used for a palette.
int CalculateEstimateForPaletteSize(const uint32 *argb, int xsize, int ysize) {
  int best_palette_bits = -1;
  double lowest_entropy = 1e99;
  std::vector<LiteralOrCopy> stream;
  BackwardReferencesHashChain(xsize, ysize, 0, argb, 0, &stream);
  for (int palette_bits = 0; palette_bits < 12; ++palette_bits) {
    Histogram histo(palette_bits);
    ComputePaletteHistogram(argb, xsize, ysize, stream, palette_bits, &histo);
    double kMakeLargePaletteSlightlyLessFavorable = 4.0;
    double cur_entropy = histo.EstimateBits() +
        kMakeLargePaletteSlightlyLessFavorable * palette_bits;
    if (palette_bits == 0 || cur_entropy < lowest_entropy) {
      best_palette_bits = palette_bits;
      lowest_entropy = cur_entropy;
    }
  }
  return best_palette_bits;
}
