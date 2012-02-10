// Copyright 2012 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// A utility for reading Canonical Huffman Codes.
//
// Author: urvang@google.com (Urvang Joshi)

#include <assert.h>
#include <stdlib.h>
#include "./huffman.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#define NON_EXISTENT_SYMBOL (-1)

void HuffmanTreeNodeInit(HuffmanTreeNode* const node) {
  if (node != NULL) {
    node->symbol_ = NON_EXISTENT_SYMBOL;
    node->child_[0] = NULL;
    node->child_[1] = NULL;
  }
}

HuffmanTreeNode* HuffmanTreeNodeNew(void) {
  HuffmanTreeNode* node = (HuffmanTreeNode*)malloc(sizeof(*node));
  if (node != NULL) HuffmanTreeNodeInit(node);
  return node;
}

void HuffmanTreeRelease(HuffmanTreeNode* const root) {
  if (root == NULL) return;
  HuffmanTreeRelease(root->child_[0]);
  HuffmanTreeRelease(root->child_[1]);
  free(root);
}

// This internal version does NOT check if node is NULL for performance reasons.
static WEBP_INLINE int IsLeafUnsafe(const HuffmanTreeNode* const node) {
  assert(node != NULL);
  return (node->child_[0] == NULL);  // Implies that node->child_[1] == NULL.
}

int HuffmanTreeIsFull(const HuffmanTreeNode* const root) {
  if (root == NULL) {
    return 0;
  } else if (IsLeafUnsafe(root)) {
    return root->symbol_ != NON_EXISTENT_SYMBOL;
  } else {
    return HuffmanTreeIsFull(root->child_[0]) &&
           HuffmanTreeIsFull(root->child_[1]);
  }
}

int HuffmanCodeLengthsToCodes(const uint32_t* const code_lengths,
                              size_t code_lengths_size,
                              int* const huff_codes) {
  size_t symbol;
  uint32_t max_code_length = 0;
  uint32_t code_len;
  int* code_length_hist = NULL;
  int curr_code;
  int* next_codes = NULL;

  if (code_lengths == NULL || code_lengths_size == 0 || huff_codes == NULL) {
    return 0;
  }

  // Calculate max code length.
  for (symbol = 0; symbol < code_lengths_size; ++symbol) {
    if (code_lengths[symbol] > max_code_length) {
      max_code_length = code_lengths[symbol];
    }
  }

  // Memory allocations.
  code_length_hist = (int*)calloc(max_code_length + 1,
                                  sizeof(*code_length_hist));
  if (code_length_hist == NULL) return 0;
  next_codes = (int*)malloc((max_code_length + 1) * sizeof(*next_codes));
  if (next_codes == NULL) {
    free(code_length_hist);
    return 0;
  }

  // Calculate code length histogram.
  if (code_length_hist == NULL) return 0;
  for(symbol = 0; symbol < code_lengths_size; ++symbol) {
    ++code_length_hist[code_lengths[symbol]];
  }
  code_length_hist[0] = 0;

  // Calculate the initial values of 'next_codes' for each code length.
  // next_codes[code_len] denotes the code to be assigned to the next symbol
  // of code length 'code_len'.
  curr_code = 0;
  next_codes[0] = -1;  // Unused, as code length = 0 implies code doesn't exist.
  for (code_len = 1; code_len <= max_code_length; ++code_len) {
    curr_code = (curr_code + code_length_hist[code_len - 1]) << 1;
    next_codes[code_len] = curr_code;
  }
  free(code_length_hist);

  // Get symbols.
  for (symbol = 0; symbol < code_lengths_size; ++symbol) {
    if (code_lengths[symbol] > 0) {
      huff_codes[symbol] = next_codes[code_lengths[symbol]]++;
    } else {
      huff_codes[symbol] = NON_EXISTENT_SYMBOL;
    }
  }
  free(next_codes);
  return 1;
}

int HuffmanTreeAddSymbol(HuffmanTreeNode* const root, int symbol,
                         uint32_t code_length, int code) {
  if (root == NULL || symbol < 0) return 0;

  if (code_length == 0) {
    // Verify we are at a leaf, so that prefix-tree property is not violated.
    if (!IsLeafUnsafe(root)) return 0;
    // Add symbol in this node.
    root->symbol_ = symbol;
    return 1;
  } else {
    if (root->symbol_ != NON_EXISTENT_SYMBOL) {
      // Violates prefix-tree property.
      return 0;
    }
    if (IsLeafUnsafe(root)) {
      // Allocate children.
      root->child_[0] = HuffmanTreeNodeNew();
      if (root->child_[0] == NULL) return 0;
      root->child_[1] = HuffmanTreeNodeNew();
      if (root->child_[1] == NULL) {
        free(root->child_[0]);
        root->child_[0] = NULL;
        return 0;
      }
    }
    {
      // Add symbol in the appropriate subtree.
      const int child_bit = (code >> (code_length - 1)) & 1;
      return HuffmanTreeAddSymbol(root->child_[child_bit], symbol,
                                  code_length - 1, code);
    }
  }
}

int HuffmanTreeBuild(HuffmanTreeNode* const root,
                     const uint32_t* const code_lengths,
                     size_t code_lengths_size) {
  size_t symbol;
  size_t num_symbols = 0;
  size_t root_symbol = 0;

  if (root == NULL || code_lengths == NULL || code_lengths_size == 0) return 0;
  if (!IsLeafUnsafe(root)) return 0;

  // Find out number of symbols & the root symbol.
  for (symbol = 0; symbol < code_lengths_size; ++symbol) {
    if (code_lengths[symbol] > 0) {
      // Note: code length = 0 indicates non-existent symbol.
      ++num_symbols;
      root_symbol = symbol;
    }
  }

  // Build tree.
  if (num_symbols < 2) {  // Trivial case.
    return HuffmanTreeAddSymbol(root, root_symbol, 0, 0);
  } else {  // Normal case.
    int ok = 0;
    // Get Huffman codes from the code lengths.
    int* const codes = (int*)malloc(code_lengths_size * sizeof(*codes));
    if (codes == NULL) goto End;
    if (!HuffmanCodeLengthsToCodes(code_lengths, code_lengths_size, codes)) {
      goto End;
    }
    // Add the Huffman codes to tree.
    for (symbol = 0; symbol < code_lengths_size; ++symbol) {
      if (codes[symbol] != NON_EXISTENT_SYMBOL) {
        if (!HuffmanTreeAddSymbol(root, symbol, code_lengths[symbol],
                                  codes[symbol])) {
          goto End;
        }
      }
    }
    ok = 1;
 End:
    free(codes);
    if (!ok) {  // Reset children.
      HuffmanTreeRelease(root->child_[0]);
      HuffmanTreeRelease(root->child_[1]);
      HuffmanTreeNodeInit(root);
    }
    return ok;
  }
}

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif
