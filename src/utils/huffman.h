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
// (1) Build a tree with explicitly given list of symbols, codes & code lengths.
//
// Code Example:
// int symbols[] = // List of symbols
// uint32_t code_lengths[] = // Corresponding list of code lengths
// int codes[] = // Corresponding list of codes
// size_t num_symbols = // size of symbols array
// HuffmanTree tree;
// HuffmanTreeBuild(code_lengths, codes, symbols, num_symbols, &tree);
// // Use 'tree' to decode some data.
// HuffmanTreeRelease(&tree);
//
// (2) Build a tree given a list of code lengths in symbol order. Here the
// symbols & codes are implicitly calculated from code lengths.
//
// Code Example:
// uint32_t code_lengths[] = // list of code lengths
// size_t num_symbols = // size of code length array
// HuffmanTree tree;
// HuffmanTreeBuild(code_lengths, NULL, NULL, num_symbols, &tree);
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

// Releases the nodes of the Huffman tree.
// Note: It does NOT free 'tree' itself.
void HuffmanTreeRelease(HuffmanTree* const tree);

// Converts Huffman code lengths to corresponding Huffman codes.
// 'huff_codes' should be pre-allocated.
// Returns true on success.
int HuffmanCodeLengthsToCodes(const uint32_t* const code_lengths,
                              size_t code_lengths_size, int* const huff_codes);

// Builds a Huffman tree given 'num_symbols' and the lists 'code_lengths',
// 'codes' & 'symbols'.
// 'symbols' & 'codes' (both together) can be passed NULL, in which case they
// are calculated implicitly assuming that 'code_lengths' is in symbol order.
// Returns false in one of the following cases:
//   - If tree or code_lengths is NULL or if num_symbols == 0
//   - If exactly one of codes or code_lengths is NULL.
//   - If the given input results in an invalid tree.
// Otherwise returns true.
int HuffmanTreeBuild(const uint32_t* const code_lengths, const int* const codes,
                     const int* const symbols, size_t num_symbols,
                     HuffmanTree* const tree);

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif  // WEBP_UTILS_HUFFMAN_H_
