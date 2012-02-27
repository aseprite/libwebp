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

// There are two ways of using this API:
//
// (1) Build a tree given a list of code lengths in symbol order.
// Code Example:
// uint32_t code_lengths[] = // list of code lengths
// size_t code_lengths_size = // size of code length array
// HuffmanTree tree;
// HuffmanTreeBuild(&tree, code_lengths, code_lengths_size);
// // Use 'tree' to decode some data.
// HuffmanTreeRelease(&tree);
//
// (2) Build a tree by adding one symbol at a time.
// Code Example:
// int symbols[] = // List of symbols
// uint32_t code_lengths[] = // Corresponding list of code lengths
// int codes[] = // Corresponding list of codes
// HuffmanTree tree;
// HuffmanTreeInit(&tree, num_symbols);
// for (int i = 0; i < num_symbols; ++i) {
//   HuffmanTreeAddSymbol(tree, symbols[i], code_lengths[i], codes[i]);
// }
// HuffmanTreeIsFull(tree);
// // Use 'tree' to decode some data.
// HuffmanTreeRelease(&tree);

#ifndef WEBP_UTILS_HUFFMAN_H_
#define WEBP_UTILS_HUFFMAN_H_

#include <stddef.h>
#include "../webp/types.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#define MAX_ALLOWED_CODE_LENGTH 15

// A node of a Huffman tree.
typedef struct HuffmanTreeNode HuffmanTreeNode;
struct HuffmanTreeNode {
  int symbol_;
  HuffmanTreeNode* child_[2];
};

// Huffman Tree.
typedef struct HuffmanTree HuffmanTree;
struct HuffmanTree {
  HuffmanTreeNode* nodes_;
  size_t nodes_count_;
  size_t next_node_idx_;
};

// Returns true if the given node is a leaf of the Huffman tree.
static WEBP_INLINE int HuffmanTreeNodeIsLeaf(
    const HuffmanTreeNode* const node) {
  return (node != NULL &&
          node->child_[0] == NULL);  // Implies that node->child_[1] == NULL.
}

// Initializes a Huffman tree given the required num_leaves.
// Returns true on success.
int HuffmanTreeInit(HuffmanTree* const tree, size_t num_leaves);

// Returns a newly allocated Huffman tree on success; and NULL on error.
HuffmanTree* HuffmanTreeNew(size_t num_leaves);

// Releases the nodes of the Huffman tree.
// Note: It does NOT free 'tree' itself.
void HuffmanTreeRelease(HuffmanTree* const tree);

// Converts Huffman code lengths to corresponding Huffman codes.
// 'huff_codes' should be pre-allocated.
// Returns true on success.
int HuffmanCodeLengthsToCodes(const uint32_t* const code_lengths,
                              size_t code_lengths_size, int* const huff_codes);

// Adds the given symbol 'symbol' with Huffman code 'code' and code length
// 'code_length' to the Huffman tree.
// Note: Tree must be pre-initialized with HuffmanTreeInit().
// Returns true on success.
int HuffmanTreeAddSymbol(HuffmanTree* const tree, int symbol,
                         uint32_t code_length, int code);

// Builds a Huffman tree given the 'code_lengths' list (having size
// 'code_lengths_size').
// Returns true on success.
int HuffmanTreeBuild(HuffmanTree* const tree,
                     const uint32_t* const code_lengths,
                     size_t code_lengths_size);

// Returns 1 if the Huffman tree is completely filled.
int HuffmanTreeIsFull(const HuffmanTree* const tree);

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif  // WEBP_UTILS_HUFFMAN_H_
