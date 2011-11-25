// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: szabadka@google.com (Zoltan Szabadka)

#ifndef WEBP_DEC_HUFFMAN_H_
#define WEBP_DEC_HUFFMAN_H_

#include <stdio.h>
#include <algorithm>
#include <vector>

class HuffmanTreeNode {
 public:
  HuffmanTreeNode() : symbol_(-1) {
    child_[0] = child_[1] = NULL;
  }
  ~HuffmanTreeNode() {
    if (child_[0]) delete child_[0];
    if (child_[1]) delete child_[1];
  }

  int symbol() const { return symbol_; }

  HuffmanTreeNode* child(int i) const { return child_[i & 1]; }

  bool IsLeaf() const {
    return child_[0] == NULL && child_[1] == NULL;
  }

  bool IsFull() const {
    if (IsLeaf()) return symbol_ != -1;
    return child_[0]->IsFull() && child_[1]->IsFull();
  }

  bool AddSymbol(int symbol, int length, unsigned int code) {
    if (symbol < 0) {
      printf("Refusing to add invalid symbol %d\n", symbol);
      return false;
    }
    if (length == 0) {
      symbol_ = symbol;
      if (!IsLeaf()) {
        printf("Insering code %x would violate prefix-tree property.\n", code);
        return false;
      }
      return true;
    }
    if (symbol_ != -1) {
      printf("Insering code %x would violate prefix-tree property.\n", code);
      return false;
    }
    if (IsLeaf()) {
      child_[0] = new HuffmanTreeNode();
      child_[1] = new HuffmanTreeNode();
    }
    int next_bit = (code >> (length - 1)) & 1;
    return child_[next_bit]->AddSymbol(symbol, length - 1, code);
  }

  bool BuildTree(const std::vector<int>& code_lengths) {
    if (!IsLeaf()) {
      printf("Attempting to build tree from non-leaf node.\n");
      return false;
    }
    int max_code_length = 0;
    int num_symbols = 0;
    int root_symbol = 0;
    for (int i = 0; i < code_lengths.size(); ++i) {
      max_code_length = std::max(max_code_length, code_lengths[i]);
      if (code_lengths[i] != 0) {
        ++num_symbols;
        root_symbol = i;
      }
    }
    if (num_symbols < 2) {
      return AddSymbol(root_symbol, 0, 0);
    }
    std::vector<int> bl_count(max_code_length + 1);
    for (int i = 0; i < code_lengths.size(); ++i) {
      ++bl_count[code_lengths[i]];
    }
    unsigned int code = 0;
    bl_count[0] = 0;
    std::vector<unsigned int> next_code(max_code_length + 1);
    for (int l = 1; l <= max_code_length; ++l) {
      code = (code + bl_count[l - 1]) << 1;
      next_code[l] = code;
    }
    for (int i = 0; i < code_lengths.size(); ++i) {
      int length = code_lengths[i];
      if (length != 0) {
        if (!AddSymbol(i, length, next_code[length])) return false;
        ++next_code[length];
      }
    }
    if (!IsFull()) {
      printf("Huffman tree is not full after building.\n");
      return false;
    }
    return true;
  }

 private:
  int symbol_;
  HuffmanTreeNode* child_[2];
};
#endif  // WEBP_DEC_HUFFMAN_H_
