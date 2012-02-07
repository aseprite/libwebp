// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)
//
// Flate like entropy encoding (Huffman) for webp lossless

#ifndef WEBP_ENTROPY_ENCODE_H_
#define WEBP_ENTROPY_ENCODE_H_

#include <stdint.h>

// This function will create a Huffman tree.
//
// The (data,length) contains the population counts.
// The tree_limit is the maximum bit depth of the Huffman codes.
//
// The depth contains the tree, i.e., how many bits are used for
// the symbol.
//
// See http://en.wikipedia.org/wiki/Huffman_coding
void CreateHuffmanTree(const int *data,
                       const int length,
                       const int tree_limit,
                       uint8_t *depth);

// Write a huffman tree from bit depths into the deflate representation
// of a Huffman tree. In deflate, the generated Huffman tree is to be
// compressed once more using a Huffman tree.
void CreateCompressedHuffmanTree(const uint8_t *depth, int len,
                                 int *num_symbols,
                                 uint8_t *tree,
                                 uint8_t *extra_bits_data);

// Get the actual bit values for a tree of bit depths.
void ConvertBitDepthsToSymbols(const uint8_t *depth, int len, uint16_t *bits);

#endif  // WEBP_ENTROPY_ENCODE_H_
