// Copyright 2012 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: Jyrki Alakuijala (jyrki@google.com)
//

#ifdef USE_LOSSLESS_ENCODER

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "./backward_references.h"
#include "./histogram.h"
#include "../utils/color_cache.h"

#define VALUES_IN_BYTE 256

static const uint8_t plane_to_code_lut[128] = {
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

static int DistanceToPlaneCode(int xsize, int dist) {
  const int yoffset = dist / xsize;
  const int xoffset = dist - yoffset * xsize;
  if (xoffset <= 8 && yoffset < 8) {
    return plane_to_code_lut[yoffset * 16 + 8 - xoffset] + 1;
  } else if (xoffset > xsize - 8 && yoffset < 7) {
    return plane_to_code_lut[(yoffset + 1) * 16 + 8 + (xsize - xoffset)] + 1;
  }
  return dist + 120;
}

static WEBP_INLINE int FindMatchLength(const uint32_t* const array1,
                                       const uint32_t* const array2,
                                       const int max_limit) {
  int match_len = 0;
  while (match_len < max_limit && array1[match_len] == array2[match_len]) {
    ++match_len;
  }
  return match_len;
}

#define HASH_BITS 18
#define HASH_SIZE (1 << HASH_BITS)
static const uint64_t kHashMultiplier = 0xc6a4a7935bd1e995ULL;
static const int kWindowSize = (1 << 20) - 120;  // A window with 1M pixels
                                                 // (4 megabytes) - 120
                                                 // special codes for short
                                                 // distances.

static WEBP_INLINE uint64_t GetPixPairHash64(const uint32_t* const argb) {
  uint64_t key = ((uint64_t)(argb[1]) << 32) | argb[0];
  key *= kHashMultiplier;
  key >>= 64 - HASH_BITS;
  return key;
}

typedef struct {
  // Stores the most recently added position with the given hash value.
  int32_t hash_to_first_index_[HASH_SIZE];
  // chain_[pos] stores the previous position with the same hash value
  // for every pixel in the image.
  int32_t* chain_;
} VP8LHashChain;

static int VP8LHashChainInit(VP8LHashChain* const p, int size) {
  int i;
  p->chain_ = (int*)malloc(size * sizeof(*p->chain_));
  if (p->chain_ == NULL) {
    return 0;
  }
  for (i = 0; i < size; ++i) {
    p->chain_[i] = -1;
  }
  for (i = 0; i < HASH_SIZE; ++i) {
    p->hash_to_first_index_[i] = -1;
  }
  return 1;
}

static void VP8LHashChainClear(VP8LHashChain* const p) {
  if (p != NULL) {
    free(p->chain_);
  }
}

static void VP8LHashChainInsert(VP8LHashChain* const p,
                                const uint32_t* const argb, int32_t pos) {
  // Insertion of two pixels at a time.
  const uint64_t hash_code = GetPixPairHash64(argb);
  p->chain_[pos] = p->hash_to_first_index_[hash_code];
  p->hash_to_first_index_[hash_code] = pos;
}

static int HashChainFindCopy(
    const VP8LHashChain* const p, int quality, int index, int xsize,
    const uint32_t* const argb, int maxlen, int* const distance_ptr,
    int* const length_ptr) {
  const uint64_t hash_code = GetPixPairHash64(&argb[index]);
  int prev_length = 0;
  int64_t best_val = 0;
  const int iter_min_mult = (quality < 50) ? 2 : (quality <= 75) ? 4 : 8;
  const int iter_min = -quality * iter_min_mult;
  int iter_cnt = 10 + (quality >> 1);
  const int min_pos = (index > kWindowSize) ? index - kWindowSize : 0;
  int32_t pos;
  int64_t val;
  int best_length = 0;
  int best_distance = 0;
  assert(xsize > 0);
  for (pos = p->hash_to_first_index_[hash_code];
       pos >= min_pos;
       pos = p->chain_[pos]) {
    int curr_length;
    if (iter_cnt < 0) {
      if (iter_cnt < iter_min || best_val >= 0xff0000) {
        break;
      }
    }
    --iter_cnt;
    if (best_length != 0 &&
        argb[pos + best_length - 1] != argb[index + best_length - 1]) {
      continue;
    }
    curr_length = FindMatchLength(argb + pos, argb + index, maxlen);
    if (curr_length < prev_length) {
      continue;
    }
    val = 65536 * curr_length;
    // Favoring 2d locality here gives savings for certain images.
    if (index - pos < 9 * xsize) {
      const int y = (index - pos) / xsize;
      int x = (index - pos) % xsize;
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
      prev_length = curr_length;
      best_val = val;
      best_length = curr_length;
      best_distance = index - pos;
      if (curr_length >= kMaxLength) {
        break;
      }
      if ((best_distance == 1 || best_distance == xsize) &&
          best_length >= 128) {
        break;
      }
    }
  }
  *distance_ptr = best_distance;
  *length_ptr = best_length;
  return best_length >= kMinLength;
}

static WEBP_INLINE void PushBackCopy(VP8LBackwardRefs* const refs, int length) {
  while (length >= kMaxLength) {
    refs->refs[refs->size++] = PixOrCopyCreateCopy(1, kMaxLength);
    length -= kMaxLength;
  }
  if (length > 0) {
    refs->refs[refs->size++] = PixOrCopyCreateCopy(1, length);
  }
}

static void BackwardReferencesRle(
    int xsize, int ysize, const uint32_t* const argb,
    VP8LBackwardRefs* const refs) {
  const int pix_count = xsize * ysize;
  int match_len = 0;
  int i;
  refs->size = 0;
  for (i = 0; i < pix_count; ++i) {
    if (i >= 1 && argb[i] == argb[i - 1]) {
      ++match_len;
    } else {
      PushBackCopy(refs, match_len);
      match_len = 0;
      refs->refs[refs->size++] = PixOrCopyCreateLiteral(argb[i]);
    }
  }
  PushBackCopy(refs, match_len);
}

static int BackwardReferencesHashChain(
    int xsize, int ysize, const uint32_t* const argb,
    int cache_bits, int quality, VP8LBackwardRefs* const refs) {
  int i;
  int ok = 0;
  int cc_init = 0;
  const int use_color_cache = (cache_bits > 0);
  const int pix_count = xsize * ysize;
  VP8LHashChain* hash_chain = (VP8LHashChain*)malloc(sizeof(*hash_chain));
  VP8LColorCache hashers;

  if (hash_chain == NULL) return 0;
  if (!(cc_init = VP8LColorCacheInit(&hashers, cache_bits)) ||
      !VP8LHashChainInit(hash_chain, pix_count)) {
    goto Error;
  }

  refs->size = 0;
  for (i = 0; i < pix_count; ) {
    // Alternative#1: Code the pixels starting at 'i' using backward reference.
    int offset = 0;
    int len = 0;
    if (i < pix_count - 1) {  // FindCopy(i,..) reads pixels at [i] and [i + 1].
      int maxlen = pix_count - i;
      if (maxlen > kMaxLength) {
        maxlen = kMaxLength;
      }
      HashChainFindCopy(hash_chain, quality, i, xsize, argb, maxlen,
                        &offset, &len);
    }
    if (len >= kMinLength) {
      // Alternative#2: Insert the pixel at 'i' as literal, and code the
      // pixels starting at 'i + 1' using backward reference.
      int offset2 = 0;
      int len2 = 0;
      int k;
      VP8LHashChainInsert(hash_chain, &argb[i], i);
      if (i < pix_count - 2) {  // FindCopy(i+1,..) reads [i + 1] and [i + 2].
        int maxlen = pix_count - (i + 1);
        if (maxlen > kMaxLength) {
          maxlen = kMaxLength;
        }
        HashChainFindCopy(hash_chain, quality,
                          i + 1, xsize, argb, maxlen, &offset2, &len2);
        if (len2 > len + 1) {
          // Alternative#2 is a better match. So push pixel at 'i' as literal.
          if (use_color_cache && VP8LColorCacheContains(&hashers, argb[i])) {
            const int ix = VP8LColorCacheGetIndex(&hashers, argb[i]);
            refs->refs[refs->size] = PixOrCopyCreateCacheIdx(ix);
          } else {
            refs->refs[refs->size] = PixOrCopyCreateLiteral(argb[i]);
          }
          ++refs->size;
          VP8LColorCacheInsert(&hashers, argb[i]);
          i++;  // Backward reference to be done for next pixel.
          len = len2;
          offset = offset2;
        }
      }
      if (len >= kMaxLength) {
        len = kMaxLength - 1;
      }
      refs->refs[refs->size++] = PixOrCopyCreateCopy(offset, len);
      for (k = 0; k < len; ++k) {
        VP8LColorCacheInsert(&hashers, argb[i + k]);
        if (k != 0 && i + k + 1 < pix_count) {
          // Add to the hash_chain (but cannot add the last pixel).
          VP8LHashChainInsert(hash_chain, &argb[i + k], i + k);
        }
      }
      i += len;
    } else {
      if (use_color_cache && VP8LColorCacheContains(&hashers, argb[i])) {
        // push pixel as a PixOrCopyCreateCacheIdx pixel
        int ix = VP8LColorCacheGetIndex(&hashers, argb[i]);
        refs->refs[refs->size] = PixOrCopyCreateCacheIdx(ix);
      } else {
        refs->refs[refs->size] = PixOrCopyCreateLiteral(argb[i]);
      }
      ++refs->size;
      VP8LColorCacheInsert(&hashers, argb[i]);
      if (i + 1 < pix_count) {
        VP8LHashChainInsert(hash_chain, &argb[i], i);
      }
      ++i;
    }
  }
  ok = 1;
Error:
  if (cc_init) VP8LColorCacheClear(&hashers);
  VP8LHashChainClear(hash_chain);
  free(hash_chain);
  return ok;
}

// -----------------------------------------------------------------------------

typedef struct {
  double alpha_[VALUES_IN_BYTE];
  double red_[VALUES_IN_BYTE];
  double literal_[PIX_OR_COPY_CODES_MAX];
  double blue_[VALUES_IN_BYTE];
  double distance_[DISTANCE_CODES_MAX];
  int cache_bits_;
} CostModel;

static int BackwardReferencesTraceBackwards(
    int xsize, int ysize, int recursive_cost_model,
    const uint32_t* const argb, int cache_bits, VP8LBackwardRefs* const refs);

static int CostModelBuild(CostModel* const p, int xsize, int ysize,
                          int recursion_level, const uint32_t* const argb,
                          int cache_bits) {
  int ok = 0;
  VP8LHistogram histo;
  VP8LBackwardRefs refs;

  if (!VP8LBackwardRefsAlloc(&refs, xsize * ysize)) goto Error;

  p->cache_bits_ = cache_bits;
  if (recursion_level > 0) {
    if (!BackwardReferencesTraceBackwards(xsize, ysize, recursion_level - 1,
                                          argb, cache_bits, &refs)) {
      goto Error;
    }
  } else {
    const int quality = 100;
    if (!BackwardReferencesHashChain(xsize, ysize, argb, cache_bits, quality,
                                     &refs)) {
      goto Error;
    }
  }
  VP8LHistogramCreate(&histo, &refs, cache_bits);
  VP8LConvertPopulationCountTableToBitEstimates(
      VP8LHistogramNumCodes(&histo),
      &histo.literal_[0], &p->literal_[0]);
  VP8LConvertPopulationCountTableToBitEstimates(
      VALUES_IN_BYTE, &histo.red_[0], &p->red_[0]);
  VP8LConvertPopulationCountTableToBitEstimates(
      VALUES_IN_BYTE, &histo.blue_[0], &p->blue_[0]);
  VP8LConvertPopulationCountTableToBitEstimates(
      VALUES_IN_BYTE, &histo.alpha_[0], &p->alpha_[0]);
  VP8LConvertPopulationCountTableToBitEstimates(
      DISTANCE_CODES_MAX, &histo.distance_[0], &p->distance_[0]);
  ok = 1;

 Error:
  VP8LClearBackwardRefs(&refs);
  return ok;
}

static WEBP_INLINE double GetLiteralCost(const CostModel* const p, uint32_t v) {
  return p->alpha_[v >> 24] +
      p->red_[(v >> 16) & 0xff] +
      p->literal_[(v >> 8) & 0xff] +
      p->blue_[v & 0xff];
}

static WEBP_INLINE double GetCacheCost(const CostModel* const p, uint32_t idx) {
  const int literal_idx = VALUES_IN_BYTE + kLengthCodes + idx;
  return p->literal_[literal_idx];
}

static WEBP_INLINE double GetLengthCost(const CostModel* const p,
                                        uint32_t length) {
  int code, extra_bits_count, extra_bits_value;
  PrefixEncode(length, &code, &extra_bits_count, &extra_bits_value);
  return p->literal_[VALUES_IN_BYTE + code] + extra_bits_count;
}

static WEBP_INLINE double GetDistanceCost(const CostModel* const p,
                                          uint32_t distance) {
  int code, extra_bits_count, extra_bits_value;
  PrefixEncode(distance, &code, &extra_bits_count, &extra_bits_value);
  return p->distance_[code] + extra_bits_count;
}

static int BackwardReferencesHashChainDistanceOnly(
    int xsize, int ysize, int recursive_cost_model, const uint32_t* const argb,
    int cache_bits, uint32_t* const dist_array) {
  int i;
  int ok = 0;
  int cc_init = 0;
  const int quality = 100;
  const int pix_count = xsize * ysize;
  const int use_color_cache = (cache_bits > 0);
  double* cost = (double*)malloc(pix_count * sizeof(*cost));
  CostModel* cost_model = (CostModel*)malloc(sizeof(*cost_model));
  VP8LHashChain* hash_chain = (VP8LHashChain*)malloc(sizeof(*hash_chain));
  VP8LColorCache hashers;

  if (cost == NULL ||
      cost_model == NULL ||
      hash_chain == NULL ||
      !(cc_init = VP8LColorCacheInit(&hashers, cache_bits)) ||
      !VP8LHashChainInit(hash_chain, pix_count)) {
    goto Error;
  }
  CostModelBuild(cost_model, xsize, ysize, recursive_cost_model, argb,
                 cache_bits);
  for (i = 0; i < pix_count; ++i) {
    cost[i] = 1e100;
  }
  // We loop one pixel at a time, but store all currently best points to
  // non-processed locations from this point.
  dist_array[0] = 0;
  for (i = 0; i < pix_count; ++i) {
    double prev_cost = 0.0;
    int shortmax;
    if (i > 0) {
      prev_cost = cost[i - 1];
    }
    for (shortmax = 0; shortmax < 2; ++shortmax) {
      int offset = 0;
      int len = 0;
      if (i < pix_count - 1) {  // FindCopy reads pixels at [i] and [i + 1].
        int maxlen = shortmax ? 2 : kMaxLength;
        if (maxlen > pix_count - i) {
          maxlen = pix_count - i;
        }
        HashChainFindCopy(hash_chain, quality, i, xsize, argb, maxlen,
                          &offset, &len);
      }
      if (len >= kMinLength) {
        const int code = DistanceToPlaneCode(xsize, offset);
        const double distance_cost =
            prev_cost + GetDistanceCost(cost_model, code);
        int k;
        for (k = 1; k < len; ++k) {
          const double cost_val =
              distance_cost + GetLengthCost(cost_model, k);
          if (cost[i + k] > cost_val) {
            cost[i + k] = cost_val;
            dist_array[i + k] = k + 1;
          }
        }
        // This if is for speedup only. It roughly doubles the speed, and
        // makes compression worse by .1 %.
        if (len >= 128 && code < 2) {
          // Long copy for short distances, let's skip the middle
          // lookups for better copies.
          // 1) insert the hashes.
          for (k = 0; k < len; ++k) {
            VP8LColorCacheInsert(&hashers, argb[i + k]);
            if (i + k + 1 < pix_count) {
              // Add to the hash_chain (but cannot add the last pixel).
              VP8LHashChainInsert(hash_chain, &argb[i + k], i + k);
            }
          }
          // 2) jump.
          i += len - 1;  // for loop does ++i, thus -1 here.
          goto next_symbol;
        }
      }
    }
    if (i < pix_count - 1) {
      VP8LHashChainInsert(hash_chain, &argb[i], i);
    }
    {
      // inserting a literal pixel
      double cost_val = prev_cost;
      double mul0 = 1.0;
      double mul1 = 1.0;
      if (recursive_cost_model == 0) {
        mul0 = 0.68;
        mul1 = 0.82;
      }
      if (use_color_cache && VP8LColorCacheContains(&hashers, argb[i])) {
        int ix = VP8LColorCacheGetIndex(&hashers, argb[i]);
        cost_val += GetCacheCost(cost_model, ix) * mul0;
      } else {
        cost_val += GetLiteralCost(cost_model, argb[i]) * mul1;
      }
      if (cost[i] > cost_val) {
        cost[i] = cost_val;
        dist_array[i] = 1;  // only one is inserted.
      }
      VP8LColorCacheInsert(&hashers, argb[i]);
    }
 next_symbol: ;
  }
  // Last pixel still to do, it can only be a single step if not reached
  // through cheaper means already.
  ok = 1;
Error:
  if (cc_init) VP8LColorCacheClear(&hashers);
  VP8LHashChainClear(hash_chain);
  free(hash_chain);
  free(cost_model);
  free(cost);
  return ok;
}

static void TraceBackwards(
    const uint32_t* const dist_array, int dist_array_size,
    uint32_t** const chosen_path, int* const chosen_path_size) {
  int i;
  // Count how many.
  int count = 0;
  for (i = dist_array_size - 1; i >= 0; ) {
    int k = dist_array[i];
    assert(k >= 1);
    ++count;
    i -= k;
  }
  // Allocate.
  *chosen_path_size = count;
  *chosen_path = (uint32_t*)malloc(count * sizeof(*chosen_path));
  // Write in reverse order.
  for (i = dist_array_size - 1; i >= 0; ) {
    int k = dist_array[i];
    assert(k >= 1);
    (*chosen_path)[--count] = k;
    i -= k;
  }
}

static int BackwardReferencesHashChainFollowChosenPath(
    int xsize, int ysize, const uint32_t* const argb, int cache_bits,
    const uint32_t* const chosen_path, int chosen_path_size,
    VP8LBackwardRefs* const refs) {
  const int quality = 100;
  const int pix_count = xsize * ysize;
  const int use_color_cache = (cache_bits > 0);
  int size = 0;
  int i = 0;
  int k;
  int ix;
  int ok = 0;
  int cc_init = 0;
  VP8LHashChain* hash_chain = (VP8LHashChain*)malloc(sizeof(*hash_chain));
  VP8LColorCache hashers;

  if (hash_chain == NULL ||
      !(cc_init = VP8LColorCacheInit(&hashers, cache_bits)) ||
      !VP8LHashChainInit(hash_chain, pix_count)) {
    goto Error;
  }

  refs->size = 0;
  for (ix = 0; ix < chosen_path_size; ++ix, ++size) {
    int offset = 0;
    int len = 0;
    int maxlen = chosen_path[ix];
    if (maxlen != 1) {
      HashChainFindCopy(hash_chain, quality,
                        i, xsize, argb, maxlen, &offset, &len);
      assert(len == maxlen);
      refs->refs[size] = PixOrCopyCreateCopy(offset, len);
      for (k = 0; k < len; ++k) {
        VP8LColorCacheInsert(&hashers, argb[i + k]);
        if (i + k + 1 < pix_count) {
          // Add to the hash_chain (but cannot add the last pixel).
          VP8LHashChainInsert(hash_chain, &argb[i + k], i + k);
        }
      }
      i += len;
    } else {
      if (use_color_cache && VP8LColorCacheContains(&hashers, argb[i])) {
        // push pixel as a color cache index
        const int idx = VP8LColorCacheGetIndex(&hashers, argb[i]);
        refs->refs[size] = PixOrCopyCreateCacheIdx(idx);
      } else {
        refs->refs[size] = PixOrCopyCreateLiteral(argb[i]);
      }
      VP8LColorCacheInsert(&hashers, argb[i]);
      if (i + 1 < pix_count) {
        VP8LHashChainInsert(hash_chain, &argb[i], i);
      }
      ++i;
    }
  }
  assert(size < refs->max_size);
  refs->size = size;
  ok = 1;
Error:
  if (cc_init) VP8LColorCacheClear(&hashers);
  VP8LHashChainClear(hash_chain);
  free(hash_chain);
  return ok;
}

// Returns 1 on success.
static int BackwardReferencesTraceBackwards(
    int xsize, int ysize, int recursive_cost_model, const uint32_t* const argb,
    int cache_bits, VP8LBackwardRefs* const refs) {
  int ok = 0;
  const int dist_array_size = xsize * ysize;
  uint32_t* chosen_path = NULL;
  int chosen_path_size = 0;
  uint32_t* const dist_array =
      (uint32_t*)malloc(dist_array_size * sizeof(*dist_array));
  if (dist_array == NULL) {
    goto Error;
  }
  if (!BackwardReferencesHashChainDistanceOnly(
      xsize, ysize, recursive_cost_model, argb, cache_bits, dist_array)) {
    free(dist_array);
    goto Error;
  }
  TraceBackwards(dist_array, dist_array_size, &chosen_path, &chosen_path_size);
  free(dist_array);
  if (!BackwardReferencesHashChainFollowChosenPath(
      xsize, ysize, argb, cache_bits, chosen_path, chosen_path_size, refs)) {
    goto Error;
  }
  ok = 1;
 Error:
  free(chosen_path);
  return ok;
}

static void BackwardReferences2DLocality(int xsize,
                                         VP8LBackwardRefs* const refs) {
  int i;
  for (i = 0; i < refs->size; ++i) {
    if (PixOrCopyIsCopy(&refs->refs[i])) {
      const int dist = refs->refs[i].argb_or_distance;
      const int transformed_dist = DistanceToPlaneCode(xsize, dist);
      refs->refs[i].argb_or_distance = transformed_dist;
    }
  }
}

int VP8LGetBackwardReferences(int width, int height,
                              const uint32_t* const argb,
                              int quality, int cache_bits, int use_2d_locality,
                              VP8LBackwardRefs* const best) {
  int ok = 0;
  int lz77_is_useful;
  VP8LBackwardRefs refs_rle, refs_lz77;
  const int num_pix = width * height;

  VP8LBackwardRefsAlloc(&refs_rle, num_pix);
  VP8LBackwardRefsAlloc(&refs_lz77, num_pix);
  VP8LInitBackwardRefs(best);
  if (refs_rle.refs == NULL || refs_lz77.refs == NULL) {
 Error1:
    VP8LClearBackwardRefs(&refs_rle);
    VP8LClearBackwardRefs(&refs_lz77);
    goto End;
  }

  if (!BackwardReferencesHashChain(width, height, argb, cache_bits, quality,
                                   &refs_lz77)) {
    goto End;
  }
  // Backward Reference using RLE only.
  BackwardReferencesRle(width, height, argb, &refs_rle);

  {
    double bit_cost_lz77, bit_cost_rle;
    VP8LHistogram* const histo = (VP8LHistogram*)malloc(sizeof(*histo));
    if (histo == NULL) goto Error1;
    // Evaluate lz77 coding
    VP8LHistogramCreate(histo, &refs_lz77, cache_bits);
    bit_cost_lz77 = VP8LHistogramEstimateBits(histo);
    // Evaluate RLE coding
    VP8LHistogramCreate(histo, &refs_rle, cache_bits);
    bit_cost_rle = VP8LHistogramEstimateBits(histo);
    // Decide if LZ77 is useful.
    lz77_is_useful = (bit_cost_lz77 < bit_cost_rle);
    free(histo);
  }

  // Choose appropriate backward reference.
  if (lz77_is_useful) {
    // TraceBackwards is costly. Run it for higher qualities.
    const int try_lz77_trace_backwards = (quality >= 75);
    *best = refs_lz77;   // default guess: lz77 is better
    VP8LClearBackwardRefs(&refs_rle);
    if (try_lz77_trace_backwards) {
      const int recursion_level = (num_pix < 320 * 200) ? 1 : 0;
      VP8LBackwardRefs refs_trace;
      if (!VP8LBackwardRefsAlloc(&refs_trace, num_pix)) {
        goto End;
      }
      if (BackwardReferencesTraceBackwards(
          width, height, recursion_level, argb, cache_bits, &refs_trace)) {
        VP8LClearBackwardRefs(&refs_lz77);
        *best = refs_trace;
      }
    }
  } else {
    VP8LClearBackwardRefs(&refs_lz77);
    *best = refs_rle;
  }

  if (use_2d_locality) {  // Use backward reference with 2D locality.
    BackwardReferences2DLocality(width, best);
  }
  ok = 1;

 End:
  if (!ok) {
    VP8LClearBackwardRefs(best);
  }
  return ok;
}

// Returns 1 on success.
static int ComputeCacheHistogram(
    const uint32_t* const argb, int xsize, int ysize,
    const VP8LBackwardRefs* const refs, int cache_bits,
    VP8LHistogram* const histo) {
  int pixel_index = 0;
  int i;
  uint32_t k;
  VP8LColorCache hashers;
  if (!VP8LColorCacheInit(&hashers, cache_bits)) {
    return 0;
  }
  for (i = 0; i < refs->size; ++i) {
    const PixOrCopy* const v = &refs->refs[i];
    if (PixOrCopyIsLiteral(v)) {
      if (cache_bits != 0 &&
          VP8LColorCacheContains(&hashers, argb[pixel_index])) {
        // push pixel as a cache index
        const int ix = VP8LColorCacheGetIndex(&hashers, argb[pixel_index]);
        const PixOrCopy token = PixOrCopyCreateCacheIdx(ix);
        VP8LHistogramAddSinglePixOrCopy(histo, &token);
      } else {
        VP8LHistogramAddSinglePixOrCopy(histo, v);
      }
    } else {
      VP8LHistogramAddSinglePixOrCopy(histo, v);
    }
    for (k = 0; k < PixOrCopyLength(v); ++k) {
      VP8LColorCacheInsert(&hashers, argb[pixel_index]);
      ++pixel_index;
    }
  }
  assert(pixel_index == xsize * ysize);
  (void)xsize;  // xsize is not used in non-debug compilations otherwise.
  (void)ysize;  // ysize is not used in non-debug compilations otherwise.
  VP8LColorCacheClear(&hashers);
  return 1;
}

// Returns how many bits are to be used for a color cache.
int VP8LCalculateEstimateForCacheSize(
    const uint32_t* const argb, int xsize, int ysize,
    int* const best_cache_bits) {
  int ok = 0;
  int cache_bits;
  double lowest_entropy = 1e99;
  VP8LBackwardRefs refs;
  static const double kSmallPenaltyForLargeCache = 4.0;
  static const int quality = 30;
  if (!VP8LBackwardRefsAlloc(&refs, xsize * ysize) ||
      !BackwardReferencesHashChain(xsize, ysize, argb, 0, quality, &refs)) {
    goto Error;
  }
  for (cache_bits = 0; cache_bits <= kColorCacheBitsMax; ++cache_bits) {
    double cur_entropy;
    VP8LHistogram histo;
    VP8LHistogramInit(&histo, cache_bits);
    ComputeCacheHistogram(argb, xsize, ysize, &refs, cache_bits, &histo);
    cur_entropy = VP8LHistogramEstimateBits(&histo) +
        kSmallPenaltyForLargeCache * cache_bits;
    if (cache_bits == 0 || cur_entropy < lowest_entropy) {
      *best_cache_bits = cache_bits;
      lowest_entropy = cur_entropy;
    }
  }
  ok = 1;
 Error:
  VP8LClearBackwardRefs(&refs);
  return ok;
}

#endif
