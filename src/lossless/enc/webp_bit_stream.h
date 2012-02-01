// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)
//
// Functions to convert the related data structures into the
// actual bit stream.
//
// These functions do bit addressing into a byte array.

#ifndef WEBP_WEBP_BIT_STREAM_H_
#define WEBP_WEBP_BIT_STREAM_H_

#include <string>
#include <vector>

#include "./backward_references.h"
#include "../common/integral_types.h"

struct BitWriter;

// The Huffman trees of flate are coded with Huffman trees.
// Here, we store the Huffman tree to code the real Huffman trees.
void StoreHuffmanTreeOfHuffmanTreeToBitMask(
    const uint8 *code_length_bitdepth, BitWriter* const bw);

// Store a Huffman tree.
void StoreHuffmanTreeToBitMask(
    const std::vector <uint8> &huffman_tree,
    const std::vector <uint8> &huffman_tree_extra_bits,
    const int num_symbols,
    const uint8 *code_length_bitdepth,
    const std::vector <uint16> &code_length_bitdepth_symbols,
    BitWriter* const bw);

// Store the deflated data with appropriate Huffman codes.
void StoreLiteralsAndBackwardReferencesToBitMask(
    const LiteralOrCopy *literal_and_length,
    const int n_literal_and_length,
    const uint8 *combined_bitdepth,
    const std::vector <uint16> &combined_bitdepth_symbols,
    const uint8 *backward_bitdepth,
    const std::vector <uint16> &backward_bitdepth_symbols,
    BitWriter* const bw);

// Store the Huffman code with given bit lengths.
void StoreHuffmanCode(const std::vector<uint8> bit_lengths,
                      bool is_color_code, BitWriter* const bw);

void StoreImageToBitMask(
    const int xsize,
    const int ysize,
    const int histobits,
    const LiteralOrCopy *literal,
    const int n_literal_and_length,
    const std::vector<uint32> &histogram_symbol,
    const std::vector< std::vector<uint8> > &bitdepth,
    const std::vector< std::vector<uint16> > &bit_symbols,
    BitWriter* const bw);

#endif  // WEBP_WEBP_BIT_STREAM_H_
