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

#ifndef WEBP_UTILS_HUFFMAN_H_
#define WEBP_UTILS_HUFFMAN_H_

#include <stddef.h>
#include "../webp/types.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

// A node of a Huffman tree.
// TODO(urvang): For improved performance, store huffman tree in a flat array
// instead of allocating every node dynamically.
typedef struct HuffmanTreeNode HuffmanTreeNode;
struct HuffmanTreeNode {
  int symbol_;
  HuffmanTreeNode* child_[2];
};

// Initializes a Huffman tree node with default values.
void HuffmanTreeNodeInit(HuffmanTreeNode* const node);

// Returns a newly allocated Huffman tree node on success; and NULL on error.
HuffmanTreeNode* HuffmanTreeNodeNew(void);

// Frees up a complete Huffman tree rooted at 'root'.
void HuffmanTreeRelease(HuffmanTreeNode* const root);

// Returns 1 if the given node is a leaf of the Huffman tree.
static WEBP_INLINE int HuffmanTreeNodeIsLeaf(
    const HuffmanTreeNode* const node) {
  return (node != NULL &&
          node->child_[0] == NULL);  // Implies that node->child_[1] == NULL.
}

// Returns 1 if the Huffman tree rooted at 'root' is completely filled.
int HuffmanTreeIsFull(const HuffmanTreeNode* const root);

// Converts Huffman code lengths to corresponding Huffman codes.
// 'huff_codes' should be pre-allocated.
// Returns 1 on success.
int HuffmanCodeLengthsToCodes(const uint32_t* const code_lengths,
                              size_t code_lengths_size,
                              int* const huff_codes);

// Adds the given symbol 'symbol' with Huffman code 'code' and code length
// 'code_length' to the Huffman tree rooted at 'root'.
// Returns 1 on success.
int HuffmanTreeAddSymbol(HuffmanTreeNode* const root, int symbol,
                         uint32_t code_length, int code);

// Builds a Huffman tree given the 'code_lengths' list (having size
// 'code_lengths_size').
// Returns 1 on success.
int HuffmanTreeBuild(HuffmanTreeNode* const root,
                     const uint32_t* const code_lengths,
                     size_t code_lengths_size);

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif  // WEBP_UTILS_HUFFMAN_H_
