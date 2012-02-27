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

static WEBP_INLINE void HuffmanTreeNodeInit(HuffmanTreeNode* const node) {
  assert(node != NULL);
  node->symbol_ = NON_EXISTENT_SYMBOL;
  node->child_[0] = NULL;
  node->child_[1] = NULL;
}

static WEBP_INLINE HuffmanTreeNode* GetNextNode(HuffmanTree* const tree) {
  if (tree->next_node_idx_ == tree->nodes_count_) {
    return NULL;
  } else {
    HuffmanTreeNode* node = &tree->nodes_[tree->next_node_idx_];
    HuffmanTreeNodeInit(node);
    ++tree->next_node_idx_;
    return node;
  }
}

// This internal version does NOT check if node is NULL for performance reasons.
static WEBP_INLINE int IsLeafUnsafe(const HuffmanTreeNode* const node) {
  assert(node != NULL);
  return (node->child_[0] == NULL);  // Implies that node->child_[1] == NULL.
}

int HuffmanTreeInit(HuffmanTree* const tree, size_t num_leaves) {
  if (tree == NULL || num_leaves == 0) return 0;
  // We allocate maximum possible nodes in the tree at once.
  // Note that a Huffman tree is a full binary tree; and in a full binary tree
  // with L leaves, the total number of nodes N = 2 * L - 1.
  tree->nodes_count_ = 2 * num_leaves - 1;
  tree->nodes_ =
      (HuffmanTreeNode*)malloc(tree->nodes_count_ * sizeof(tree->nodes_[0]));
  if (tree->nodes_ == NULL) return 0;
  HuffmanTreeNodeInit(&tree->nodes_[0]);  // Initialize root.
  tree->next_node_idx_ = 1;
  return 1;
}

HuffmanTree* HuffmanTreeNew(size_t num_leaves) {
  HuffmanTree* tree;
  if (num_leaves == 0) return NULL;
  tree = (HuffmanTree*)malloc(sizeof(*tree));
  if (tree == NULL) return NULL;
  if (!HuffmanTreeInit(tree, num_leaves)) {
    free(tree);
    return NULL;
  }
  return tree;
}

void HuffmanTreeRelease(HuffmanTree* const tree) {
  if (tree != NULL) {
    free(tree->nodes_);
    tree->nodes_ = NULL;
    tree->nodes_count_ = 0;
  }
}

static int TreeIsFullInternal(const HuffmanTree* const tree) {
  return (tree->next_node_idx_ == tree->nodes_count_);
}

int HuffmanTreeIsFull(const HuffmanTree* const tree) {
  return (tree != NULL && TreeIsFullInternal(tree));
}

int HuffmanCodeLengthsToCodes(const uint32_t* const code_lengths,
                              size_t code_lengths_size, int* const huff_codes) {
  size_t symbol;
  uint32_t code_len;
  int code_length_hist[MAX_ALLOWED_CODE_LENGTH + 1] = { 0 };
  int curr_code;
  int next_codes[MAX_ALLOWED_CODE_LENGTH + 1] = { 0 };
  uint32_t max_code_length = 0;

  assert(code_lengths != NULL);
  assert(code_lengths_size > 0);
  assert(huff_codes != NULL);

  // Calculate max code length.
  for (symbol = 0; symbol < code_lengths_size; ++symbol) {
    if (code_lengths[symbol] > max_code_length) {
      max_code_length = code_lengths[symbol];
    }
  }
  if (max_code_length > MAX_ALLOWED_CODE_LENGTH) return 0;

  // Calculate code length histogram.
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

  // Get symbols.
  for (symbol = 0; symbol < code_lengths_size; ++symbol) {
    if (code_lengths[symbol] > 0) {
      huff_codes[symbol] = next_codes[code_lengths[symbol]]++;
    } else {
      huff_codes[symbol] = NON_EXISTENT_SYMBOL;
    }
  }
  return 1;
}

static int HuffmanTreeAddSymbolInternal(HuffmanTree* const tree,
                                        HuffmanTreeNode* const node,
                                        int symbol, uint32_t code_length,
                                        int code) {
  if (code_length == 0) {
    // Verify we are at a leaf, so that prefix-tree property is not violated.
    if (!IsLeafUnsafe(node)) return 0;
    // Add symbol in this node.
    node->symbol_ = symbol;
    return 1;
  } else {
    if (node->symbol_ != NON_EXISTENT_SYMBOL) {
      // Violates prefix-tree property.
      return 0;
    }
    if (IsLeafUnsafe(node)) {
      // Allocate children.
      node->child_[0] = GetNextNode(tree);
      node->child_[1] = GetNextNode(tree);
      if (node->child_[0] == NULL || node->child_[1] == NULL) return 0;
    }
    {
      // Add symbol in the appropriate subtree.
      const int child_bit = (code >> (code_length - 1)) & 1;
      return HuffmanTreeAddSymbolInternal(tree, node->child_[child_bit],
                                          symbol, code_length - 1, code);
    }
  }
}

int HuffmanTreeAddSymbol(HuffmanTree* const tree, int symbol,
                         uint32_t code_length, int code) {
  if (tree == NULL || tree->nodes_ == NULL || symbol < 0) return 0;
  return HuffmanTreeAddSymbolInternal(tree, &tree->nodes_[0], symbol,
                                      code_length, code);
}

int HuffmanTreeBuild(HuffmanTree* const tree,
                     const uint32_t* const code_lengths,
                     size_t code_lengths_size) {
  size_t symbol;
  size_t num_symbols = 0;
  size_t root_symbol = 0;

  if (tree == NULL || code_lengths == NULL || code_lengths_size == 0) return 0;

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
    if (!HuffmanTreeInit(tree, 1)) return 0;
    return HuffmanTreeAddSymbolInternal(tree, &tree->nodes_[0],
                                        root_symbol, 0, 0);
  } else {  // Normal case.
    int ok = 0;

    // Get Huffman codes from the code lengths.
    int* const codes = (int*)malloc(code_lengths_size * sizeof(*codes));
    if (codes == NULL) goto End;

    tree->nodes_ = NULL;
    if (!HuffmanCodeLengthsToCodes(code_lengths, code_lengths_size, codes)) {
      goto End;
    }

    // Initialize the HuffmanTree based on max_code_length.
    if (!HuffmanTreeInit(tree, num_symbols)) return 0;

    // Add the Huffman codes to tree_lcl.
    for (symbol = 0; symbol < code_lengths_size; ++symbol) {
      if (codes[symbol] != NON_EXISTENT_SYMBOL) {
        if (!HuffmanTreeAddSymbolInternal(tree, &tree->nodes_[0], symbol,
                                          code_lengths[symbol],
                                          codes[symbol])) {
          goto End;
        }
      }
    }
    ok = 1;
 End:
    free(codes);
    ok = ok && TreeIsFullInternal(tree);
    if (!ok) HuffmanTreeRelease(tree);
    return ok;
  }
}

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif
