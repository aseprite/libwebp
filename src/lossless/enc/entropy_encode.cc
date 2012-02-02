// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)
//
// Flate like entropy encoding (Huffman) for webp lossless.

#include "entropy_encode.h"

#include <stdlib.h>
#include <string.h>

#include "../common/integral_types.h"

namespace {

struct HuffmanTree {
  int total_count_;
  int value_;
  int pool_index_left_;
  int pool_index_right_;
};

// Sort the root nodes, most popular first.
int CompHuffmanTree(const void *vp0, const void *vp1) {
  const HuffmanTree *v0 = (const HuffmanTree *)vp0;
  const HuffmanTree *v1 = (const HuffmanTree *)vp1;
  if (v0->total_count_ > v1->total_count_) {
    return -1;
  } else if (v0->total_count_ < v1->total_count_) {
    return 1;
  } else {
    if (v0->value_ < v1->value_) {
      return -1;
    }
    if (v0->value_ > v1->value_) {
      return 1;
    }
    return 0;
  }
}

void SetDepth(const HuffmanTree &p,
              HuffmanTree *pool,
              uint8 *depth,
              const int level) {
  if (p.pool_index_left_ >= 0) {
    SetDepth(pool[p.pool_index_left_], pool, depth, level + 1);
    SetDepth(pool[p.pool_index_right_], pool, depth, level + 1);
  } else {
    depth[p.value_] = level;
  }
}

}  // namespace

// This function will create a Huffman tree.
//
// The catch here is that the tree cannot be arbitrarily deep.
// Deflate specifies a maximum depth of 15 bits for "code trees"
// and 7 bits for "code length code trees."
//
// count_limit is the value that is to be faked as the minimum value
// and this minimum value is raised until the tree matches the
// maximum length requirement.
//
// This algorithm is not of excellent performance for very long data blocks,
// especially when population counts are longer than 2**tree_limit, but
// we are not planning to use this with extremely long blocks.
//
// See http://en.wikipedia.org/wiki/Huffman_coding
void CreateHuffmanTree(const int* const histogram, int histogram_size,
                       int tree_depth_limit,
                       uint8* const bit_depths) {
  // For block sizes with less than 64k symbols we never need to do a
  // second iteration of this loop.
  // If we actually start running inside this loop a lot, we would perhaps
  // be better off with the Katajainen algorithm.
  for (int count_limit = 1; ; count_limit *= 2) {
    int tree_size = 0;
    for (int i = 0; i < histogram_size; ++i) {
      if (histogram[i]) {
        ++tree_size;
      }
    }
    // 3 * tree_size is enough to cover all the nodes representing a
    // population and all the inserted nodes combining two existing nodes.
    // The tree pool needs 2 * (tree_size - 1) entities, and the
    // tree needs exactly tree_size entities.
    HuffmanTree *tree =
        (HuffmanTree *)malloc(3 * tree_size * sizeof(HuffmanTree));
    {
      int j = 0;
      for (int i = 0; i < histogram_size; ++i) {
        if (histogram[i]) {
          const int count =
              (histogram[i] < count_limit) ? count_limit : histogram[i];
          tree[j].total_count_ = count;
          tree[j].value_ = i;
          tree[j].pool_index_left_ = -1;
          tree[j].pool_index_right_ = -1;
          ++j;
        }
      }
    }
    qsort((void *)tree, tree_size, sizeof(HuffmanTree), CompHuffmanTree);
    HuffmanTree *tree_pool = tree + tree_size;
    int tree_pool_size = 0;
    if (tree_size >= 2) {
      while (tree_size >= 2) {  // Finish when we have only one root.
        tree_pool[tree_pool_size] = tree[tree_size - 1];
        ++tree_pool_size;
        tree_pool[tree_pool_size] = tree[tree_size - 2];
        ++tree_pool_size;
        int count =
            tree_pool[tree_pool_size - 1].total_count_ +
            tree_pool[tree_pool_size - 2].total_count_;
        tree_size -= 2;
        {
          int k = 0;
          // Search for the insertion point.
          for (k = 0; k < tree_size; ++k) {
            if (tree[k].total_count_ <= count) {
              break;
            }
          }
          memmove(tree + (k + 1),
                  tree + k,
                  (tree_size - k) * sizeof(HuffmanTree));
          tree[k].total_count_ = count;
          tree[k].value_ = -1;

          tree[k].pool_index_left_ = tree_pool_size - 1;
          tree[k].pool_index_right_ = tree_pool_size - 2;
          tree_size = tree_size + 1;
        }
      }
      SetDepth(tree[0], tree_pool, bit_depths, 0);
    } else {
      if (tree_size == 1) {
        // Only one element.
        bit_depths[tree[0].value_] = 1;
      }
    }
    free(tree);
    // We need to pack the Huffman tree in tree_depth_limit bits.
    // If this was not successful, add fake entities to the lowest values
    // and retry.
    {
      int max_depth = bit_depths[0];
      for (int j = 1; j < histogram_size; ++j) {
        if (max_depth < bit_depths[j]) {
          max_depth = bit_depths[j];
        }
      }
      if (max_depth <= tree_depth_limit) {
        break;
      }
    }
  }
}

void WriteHuffmanTreeRepetitions(
    const int value,
    const int prev_value,
    int repetitions,
    int *num_symbols,
    uint8 *tree,
    uint8 *extra_bits_data) {
  if (value != prev_value) {
    tree[*num_symbols] = value;
    extra_bits_data[*num_symbols] = 0;
    ++(*num_symbols);
    --repetitions;
  }
  while (repetitions >= 1) {
    if (repetitions < 3) {
      for (int i = 0; i < repetitions; ++i) {
        tree[*num_symbols] = value;
        extra_bits_data[*num_symbols] = 0;
        ++(*num_symbols);
      }
      return;
    } else if (repetitions < 7) {
      // 3 to 6 left
      tree[*num_symbols] = 16;
      extra_bits_data[*num_symbols] = repetitions - 3;
      ++(*num_symbols);
      return;
    } else {
      tree[*num_symbols] = 16;
      extra_bits_data[*num_symbols] = 3;
      ++(*num_symbols);
      repetitions -= 6;
    }
  }
}

void WriteHuffmanTreeRepetitionsZeros(
    const int value,
    int repetitions,
    int *num_symbols,
    uint8 *tree,
    uint8 *extra_bits_data) {
  while (repetitions >= 1) {
    if (repetitions < 3) {
      for (int i = 0; i < repetitions; ++i) {
        tree[*num_symbols] = value;
        extra_bits_data[*num_symbols] = 0;
        ++(*num_symbols);
      }
      return;
    } else if (repetitions < 11) {
      tree[*num_symbols] = 17;
      extra_bits_data[*num_symbols] = repetitions - 3;
      ++(*num_symbols);
      return;
    } else if (repetitions < 139) {
      tree[*num_symbols] = 18;
      extra_bits_data[*num_symbols] = repetitions - 11;
      ++(*num_symbols);
      return;
    } else {
      tree[*num_symbols] = 18;
      extra_bits_data[*num_symbols] = 0x7f;  // 138 repeated 0s
      ++(*num_symbols);
      repetitions -= 138;
    }
  }
}

void CreateCompressedHuffmanTree(const uint8 *depth,
                                 int depth_size,
                                 int *num_symbols,
                                 uint8 *tree,
                                 uint8 *extra_bits_data) {
  int prev_value = 8;  // 8 is the initial value for rle.
  for (uint32 i = 0; i < depth_size;) {
    const int value = depth[i];
    int reps = 1;
    for (uint32 k = i + 1; k < depth_size && depth[k] == value; ++k) {
      ++reps;
    }
    if (value == 0) {
      WriteHuffmanTreeRepetitionsZeros(value, reps,
                                       num_symbols,
                                       tree, extra_bits_data);
    } else {
      WriteHuffmanTreeRepetitions(value, prev_value, reps,
                                  num_symbols,
                                  tree, extra_bits_data);
      prev_value = value;
    }
    i += reps;
  }
}

namespace {

uint32 ReverseBits(int num_bits, uint32 bits) {
  uint32 retval = 0;
  for (int i = 0; i < num_bits; ++i) {
    retval <<= 1;
    retval |= bits & 1;
    bits >>= 1;
  }
  return retval;
}

}  // namespace

void ConvertBitDepthsToSymbols(const uint8 *depth, int len,
                               uint16 *bits) {
  // This function is based on RFC 1951.
  //
  // In deflate, all bit depths are [1..15]
  // 0 bit depth means that the symbol does not exist.

  const int kMaxBits = 16;  // 0..15 are values for bits
  uint32 bl_count[kMaxBits] = { 0 };
  {
    for (uint32 i = 0; i < len; ++i) {
      ++bl_count[depth[i]];
    }
    bl_count[0] = 0;
  }
  uint32 next_code[kMaxBits];
  next_code[0] = 0;
  {
    int code = 0;
    for (int bits = 1; bits < kMaxBits; ++bits) {
      code = (code + bl_count[bits - 1]) << 1;
      next_code[bits] = code;
    }
  }
  for (uint32 i = 0; i < len; ++i) {
    if (depth[i]) {
      bits[i] = ReverseBits(depth[i], next_code[depth[i]]++);
    }
  }
}
