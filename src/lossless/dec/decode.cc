// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: szabadka@google.com (Zoltan Szabadka)

#include "decode.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "./bit_reader.h"
#include "./constants.h"
#include "./huffman.h"
#include "./image_transform.h"
#include "../common/color_cache.h"

#define TAG_SIZE 4
#define CHUNK_HEADER_SIZE 8
#define RIFF_HEADER_SIZE 12
#define HEADER_SIZE      (RIFF_HEADER_SIZE + CHUNK_HEADER_SIZE)
#define SIGNATURE_SIZE   1

#define VERIFY(a) \
  if (!(a)) { \
    fprintf(stderr, "Failed at %s:%d\n", __FILE__, __LINE__); \
    perror("Unrecoverable error"); \
    abort(); \
  };

static const int kMaxImageTransforms = 100;

static int DecodeImageInternal(int xsize, int ysize,
                               int read_transforms,
                               BitReader* const br,
                               uint32_t** const argb_image);

// -----------------------------------------------------------------------------
// COMMON FUNCTIONS

static int MetaSize(int size, int bits) {
  return (size + (1 << bits) - 1) >> bits;
}

static void ReadTransform(int* xsize, int* ysize,
                          ImageTransform* transform,
                          BitReader* br) {
  transform->type = (ImageTransformType)ReadBits(br, 2);
  transform->xsize = *xsize;
  transform->ysize = *ysize;
  switch (transform->type) {
    case PREDICTOR_TRANSFORM:
    case CROSS_COLOR_TRANSFORM:
      {
        transform->data = malloc(sizeof(PerTileTransformData));
        PerTileTransformData* data = (PerTileTransformData*)transform->data;
        data->bits = ReadBits(br, 4);
        DecodeImageInternal(MetaSize(*xsize, data->bits),
                            MetaSize(*ysize, data->bits),
                            0, br, &data->image);
      }
      break;
    case SUBTRACT_GREEN:
      // There is no associated data for this transform.
      transform->data = NULL;
      break;
    case COLOR_INDEXING_TRANSFORM:
      {
        int ncolors = ReadBits(br, 8) + 1;
        transform->data = malloc(sizeof(PerTileTransformData));
        PerTileTransformData* data = (PerTileTransformData*)transform->data;
        if (ncolors <= 16) {
          data->bits = 1;
          if (ncolors <= 4) {
            data->bits = 2;
          }
          if (ncolors <= 2) {
            data->bits = 3;
          }
        } else {
          data->bits = 0;
        }
        *xsize = MetaSize(*xsize, data->bits);
        DecodeImageInternal(ncolors, 1, 0, br, &data->image);
        for (int i = 1; i < ncolors; ++i) {
          data->image[i] = Add(data->image[i], data->image[i - 1]);
        }
      }
      break;
  }
}

// -----------------------------------------------------------------------------
// HUFFMAN CODING

static inline uint32_t GetCopyDistance(uint32_t distance_symbol,
                                       BitReader* br) {
  if (distance_symbol < 4) {
    return distance_symbol + 1;
  }
  uint32_t extra_bits = (distance_symbol - 2) >> 1;
  uint32_t offset = (2 + (distance_symbol & 1)) << extra_bits;
  return offset + ReadBits(br, extra_bits) + 1;
}

static inline uint32_t GetCopyLength(uint32_t length_symbol,
                                     BitReader* br) {
  // Length and distance prefixes are encoded the same way.
  return GetCopyDistance(length_symbol, br);
}

static int PlaneCodeToDistance(int xsize, int ysize, int plane_code) {
  if (plane_code > 120) {
    return plane_code - 120;
  }
  const int dist_code = code_to_plane_lut[plane_code - 1];
  const int yoffset = dist_code >> 4;
  const int xoffset = 8 - (dist_code & 0xf);
  return yoffset * xsize + xoffset;
}

static int AlphabetSize(int tree_type, int palette_size) {
  if (tree_type == 0) return 256 + palette_size + kNumLengthSymbols;
  if (tree_type == 1) return 256;
  if (tree_type == 2) return kNumDistanceSymbols;
  return 0;
}

// Five Huffman codes are used at each meta code:
// 1. green + length prefix codes + palette codes,
// 2. alpha,
// 3. red,
// 4. blue, and,
// 5. distance prefix codes.
static const int kHuffmanCodesPerMetaCode = 5;

// For the above 5 indices, there are three tree-types.
// 0 is for green + length prefix codes + palette codes,
// 1 is for alpha, red and blue,
// 2 is for distance prefix codes.
static inline int HuffmanCodeIndexToTreeType(int v) {
  static const char ret[5] = { 0, 1, 1, 1, 2 };
  return ret[v];
}

int GetMetaIndex(int huffman_xsize, int bits, const uint32_t* image,
                 int x, int y) {
  y >>= bits;
  x >>= bits;
  if (bits == 0) {
    return 0;
  }
  return image[y * huffman_xsize + x];
}

static const int kCodeLengthCodes = 19;
static const uint8_t kCodeLengthCodeOrder[kCodeLengthCodes] = {
  17, 18, 0, 1, 2, 3, 4, 5, 16, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
};
static const int kCodeLengthLiterals = 16;
static const int kCodeLengthRepeatCode = 16;
static const int kCodeLengthExtraBits[3] = { 2, 3, 7 };
static const int kCodeLengthRepeatOffsets[3] = { 3, 3, 11 };

// Returns an integer code decoded with the Huffman code.
//
// FillBitWindow(br) needs to be called at minimum every second call
// to ReadSymbol.
static inline int ReadSymbol(const HuffmanTreeNode& root, BitReader* br) {
  const HuffmanTreeNode* node = &root;
  while (!node->IsLeaf()) {
    node = node->child(ReadOneBitUnsafe(br));
  }
  assert(node);
  assert(node->symbol() != -1);
  return node->symbol();
}

static int ReadCodeLengthTree(BitReader* const br,
                              HuffmanTreeNode* root) {
  int code_lengths[kCodeLengthCodes] = { 0 };
  const int num_codes = 4 + ReadBits(br, 4);
  int i;
  if (num_codes > kCodeLengthCodes) return 0;
  for (i = 0; i < num_codes; ++i) {
    code_lengths[kCodeLengthCodeOrder[i]] = ReadBits(br, 3);
  }
  return root->BuildTree(code_lengths, kCodeLengthCodes);
}

static int ReadHuffmanCodeLengths(const HuffmanTreeNode& decoder_root,
                                  BitReader* br,
                                  int code_lengths_size,
                                  int* code_lengths) {
  bool use_length = ReadBits(br, 1);
  int max_length = 0;
  if (use_length) {
    const int length_nbits = (ReadBits(br, 3) + 1) * 2;
    max_length = ReadBits(br, length_nbits) + 2;
  }
  int previous = 8;
  int num_symbols = 0;
  for (int i = 0; i < code_lengths_size; ++i) {
    if (use_length && ++num_symbols > max_length) break;
    FillBitWindow(br);
    int code_length = ReadSymbol(decoder_root, br);
    VERIFY(code_length < kCodeLengthCodes);
    if (code_length < kCodeLengthLiterals) {
      code_lengths[i] = code_length;
      if (code_length != 0) previous = code_length;
    } else {
      int repeat_ix = code_length - kCodeLengthLiterals;
      int extra_bits = kCodeLengthExtraBits[repeat_ix];
      int repeat_offset = kCodeLengthRepeatOffsets[repeat_ix];
      int repeat = ReadBits(br, extra_bits) + repeat_offset;
      bool use_previous = code_length == kCodeLengthRepeatCode;
      for (int k = 0; k < repeat; ++k) {
        code_lengths[i + k] = use_previous ? previous : 0;
      }
      i += repeat - 1;
    }
  }
  return 1;
}

static int ReadHuffmanCode(int alphabet_size,
                           BitReader* const br,
                           HuffmanTreeNode* const root) {
  int ok = 1;
  const int simple_code = ReadBits(br, 1);
  if (simple_code) {
    const int num_symbols = ReadBits(br, 1) + 1;
    const int short_first_code = ReadBits(br, 1);
    // The first code is either 1 bit or 8 bit code.
    const int first_code = ReadBits(br, 1 + 7 * short_first_code);
    ok = root->AddSymbol(first_code, num_symbols - 1, 0);
    // The second code, when present, is always 8 bit long.
    if (ok && num_symbols == 2) {
      ok = root->AddSymbol(ReadBits(br, 8), num_symbols - 1, 1);
    }
    VERIFY(!ok || root->IsFull());
  } else {
    HuffmanTreeNode decoder_root;
    int* code_lengths = (int *)malloc(alphabet_size * sizeof(int));
    memset(code_lengths, 0, alphabet_size * sizeof(code_lengths[0]));
    ok = ReadCodeLengthTree(br, &decoder_root) &&
        ReadHuffmanCodeLengths(decoder_root, br, alphabet_size, code_lengths) &&
        root->BuildTree(code_lengths, alphabet_size);
    free(code_lengths);
    if (!ok) {
      printf("error #2\n");
    }
  }
  return ok;
}

static int DecodeImageInternal(const int original_xsize,
                               const int original_ysize,
                               const int read_transforms,
                               BitReader* const br,
                               uint32_t** const argb_image) {
  int xsize = original_xsize;
  int ysize = original_ysize;
  int ok = 1;

  ImageTransform* transforms =
      (ImageTransform*)malloc(kMaxImageTransforms * sizeof(ImageTransform));
  int num_transforms = 0;
  if (read_transforms) {
    while (ReadBits(br, 1)) {
      ReadTransform(&xsize, &ysize, &transforms[num_transforms], br);
      ++num_transforms;
      VERIFY(num_transforms < kMaxImageTransforms);
    }
  }

  const bool use_palette = ReadBits(br, 1);
  const int palette_code_bits = use_palette ? ReadBits(br, 4) : 0;
  const int palette_size = use_palette ? 1 << palette_code_bits : 0;
  VP8LColorCache* hashers = use_palette ? new VP8LColorCache : NULL;
  if (hashers) {
    VP8LColorCacheInit(hashers, palette_code_bits);
  }

  bool use_meta_codes = ReadBits(br, 1);
  int huffman_bits = 0;
  uint32_t* huffman_image = NULL;
  int num_huffman_trees = kHuffmanCodesPerMetaCode;
  if (use_meta_codes) {
    huffman_bits = ReadBits(br, 4);
    const int huffman_xsize = MetaSize(xsize, huffman_bits);
    const int huffman_ysize = MetaSize(ysize, huffman_bits);
    DecodeImageInternal(huffman_xsize, huffman_ysize, 0, br, &huffman_image);
    const int huffman_pixels = huffman_xsize * huffman_ysize;
    int num_meta_codes = 0;
    for (int i = 0; i < huffman_pixels; ++i) {
      // The actual value is stored in red (high bits) and green (low bits).
      const int index = (huffman_image[i] >> 8) & 0xffff;
      // Strip alpha (in bits [24..17]).
      huffman_image[i] = index;
      if (num_meta_codes <= index) num_meta_codes = index + 1;
    }
    num_huffman_trees = num_meta_codes * kHuffmanCodesPerMetaCode;
  }


  HuffmanTreeNode *htrees = (HuffmanTreeNode *)
      malloc(num_huffman_trees * sizeof(HuffmanTreeNode));
  for (int i = 0; i < num_huffman_trees; ++i) {
    htrees[i].child_[0] = NULL;
    htrees[i].child_[1] = NULL;
    htrees[i].symbol_ = -1;
  }
  for (int i = 0; ok && i < num_huffman_trees; ++i) {
    int type = HuffmanCodeIndexToTreeType(i % 5);
    int alphabet_size = AlphabetSize(type, palette_size);
    ok = ReadHuffmanCode(alphabet_size, br, &htrees[i]);
  }

  uint32_t *image = (uint32_t*)malloc(xsize * ysize * sizeof(uint32_t));
  int x = 0;
  int y = 0;
  int red = 0;
  int blue = 0;
  int alpha = 0xff000000;
  int meta_ix = -1;
  const HuffmanTreeNode *huff_green = 0;
  const HuffmanTreeNode *huff_red = 0;
  const HuffmanTreeNode *huff_blue = 0;
  const HuffmanTreeNode *huff_alpha = 0;
  const HuffmanTreeNode *huff_dist = 0;
  const int huffman_mask = huffman_bits == 0 ? ~0 : (1 << huffman_bits) - 1;
  const int huffman_xsize = (xsize + (1 << huffman_bits) - 1) >> huffman_bits;
  for (int pos = 0; pos < xsize * ysize; ) {
    FillBitWindow(br);
    if ((x & huffman_mask) == 0) {
      // Only update the huffman code when moving from one block to the
      // next.
      const int meta_index = kHuffmanCodesPerMetaCode *
          GetMetaIndex(huffman_xsize, huffman_bits, huffman_image, x, y);
      if (meta_ix != meta_index) {
        meta_ix = meta_index;
        huff_green = &htrees[meta_ix];
        huff_red = &htrees[meta_ix + 1];
        huff_blue = &htrees[meta_ix + 2];
        huff_alpha = &htrees[meta_ix + 3];
        huff_dist = &htrees[meta_ix + 4];
      }
    }
    int green = ReadSymbol(*huff_green, br);
    // Literal
    if (green < 256) {
      red = ReadSymbol(*huff_red, br) << 16;
      FillBitWindow(br);
      blue = ReadSymbol(*huff_blue, br);
      alpha = ReadSymbol(*huff_alpha, br) << 24;

      uint32_t argb = alpha + red + (green << 8) + blue;
      image[pos] = argb;
      if (hashers) VP8LColorCacheInsert(hashers, argb);
      ++x;
      if (x >= xsize) {
        x = 0;
        ++y;
      }
      ++pos;
      continue;
    }

    // Backward reference
    if (green < 256 + kNumLengthSymbols) {
      const int length_symbol = green - 256;
      const int length = GetCopyLength(length_symbol, br);
      // We have read the length code prefix + extra bits for the length,
      // so reading the next 15 bits can exhaust the bit window.
      // We must fill the window before the next read.
      FillBitWindow(br);
      const int dist_symbol = ReadSymbol(*huff_dist, br);
      uint32_t dist = GetCopyDistance(dist_symbol, br);
      dist = PlaneCodeToDistance(xsize, ysize, dist);
      VERIFY(dist <= pos);
      VERIFY(pos + length <= xsize * ysize);
      if (!hashers) {
        for (int i = 0; i < length; ++i) {
          image[pos] = image[pos - dist];
          ++pos;
        }
        x += length;
        while (x >= xsize) {
          x -= xsize;
          ++y;
        }
      } else {
        for (int i = 0; i < length; ++i) {
          image[pos] = image[pos - dist];
          VP8LColorCacheInsert(hashers, image[pos]);
          ++pos;
          ++x;
          if (x >= xsize) {
            x = 0;
            ++y;
          }
        }
      }
      if (pos == xsize * ysize) {
        break;
      }
      const int meta_index = kHuffmanCodesPerMetaCode *
          GetMetaIndex(huffman_xsize, huffman_bits, huffman_image, x, y);
      if (meta_ix != meta_index) {
        meta_ix = meta_index;
        huff_green = &htrees[meta_ix];
        huff_red = &htrees[meta_ix + 1];
        huff_blue = &htrees[meta_ix + 2];
        huff_alpha = &htrees[meta_ix + 3];
        huff_dist = &htrees[meta_ix + 4];
      }
      continue;
    }

    // Palette
    {
      int palette_symbol = green - 256 - kNumLengthSymbols;
      const uint32_t argb = VP8LColorCacheLookup(hashers, palette_symbol);
      image[pos] = argb;
      VP8LColorCacheInsert(hashers, argb);
      ++x;
      if (x >= xsize) {
        x = 0;
        ++y;
      }
      ++pos;
      continue;
    }

    printf("Error: Could not interpret GREEN symbol %d\n", green);
    abort();
  }

  if (hashers) {
    VP8LColorCacheDelete(hashers);
    delete hashers;
  }
  free(huffman_image);
  {
    for (int i = 0; i < num_huffman_trees; ++i) {
      delete htrees[i].child_[0];
      delete htrees[i].child_[1];
    }
    free(htrees);
  }

  *argb_image = image;
  for (int i = num_transforms - 1; i >= 0; --i) {
    ApplyInverseImageTransform(transforms[i], argb_image);
    FreeImageTransformData(transforms[i]);
  }
  free(transforms);
  return ok;
}

int DecodeWebpLLImage(size_t encoded_image_size,
                      const uint8_t* const encoded_image,
                      int* xsize,
                      int* ysize,
                      uint32_t** argb_image) {
  if (encoded_image_size < HEADER_SIZE + SIGNATURE_SIZE) return false;
  const uint8_t* webpll_data = encoded_image + HEADER_SIZE;
  const size_t webpll_size = encoded_image_size - HEADER_SIZE;

  BitReader br;
  InitBitReader(&br, webpll_data, webpll_size);
  const uint8_t sig = ReadBits(&br, 8 * SIGNATURE_SIZE);
  if (sig != 0x64) return false;

  *xsize = ReadBits(&br, 14) + 1;
  *ysize = ReadBits(&br, 14) + 1;
  DecodeImageInternal(*xsize, *ysize, 1, &br, argb_image);
  return true;
}
