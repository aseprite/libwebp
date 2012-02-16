// Copyright 2012 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// main entry for the decoder
//
// Author: Vikas Arora (vikaas.arora@gmail.com)
//         jyrki@google.com (Jyrki Alakuijala)

#include <stdlib.h>

#include "./vp8li.h"
#include "../utils/bit_reader.h"
#include "../utils/color_cache.h"
#include "../utils/huffman.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#define LOSSLESS_MAGIC_BYTE 0x64

static const size_t kHeaderBytes = 5;
static const uint32_t kImageSizeBits = 14;

static const int kCodeLengthLiterals = 16;
static const int kCodeLengthRepeatCode = 16;
static const int kCodeLengthExtraBits[3] = { 2, 3, 7 };
static const int kCodeLengthRepeatOffsets[3] = { 3, 3, 11 };

#define NUM_LENGTH_CODES    24
#define NUM_DISTANCE_CODES  40
#define NUM_CODES_PER_BYTE 256
// -----------------------------------------------------------------------------
//  Five Huffman codes are used at each meta code:
//  1. green + length prefix codes + color cache codes,
//  2. alpha,
//  3. red,
//  4. blue, and,
//  5. distance prefix codes.
typedef enum {
  GREEN = 0,
  RED   = 1,
  BLUE  = 2,
  ALPHA = 3,
  DIST  = 4,
  HUFFMAN_CODES_PER_META_CODE = 5
} HuffIndex;

static const uint16_t kAlphabetSize[HUFFMAN_CODES_PER_META_CODE] = {
  NUM_CODES_PER_BYTE + NUM_LENGTH_CODES,
  NUM_CODES_PER_BYTE, NUM_CODES_PER_BYTE, NUM_CODES_PER_BYTE,
  NUM_DISTANCE_CODES};


#define NUM_CODE_LENGTH_CODES       19
static const uint8_t kCodeLengthCodeOrder[NUM_CODE_LENGTH_CODES] = {
  17, 18, 0, 1, 2, 3, 4, 5, 16, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
};

#define CODE_TO_PLANE_CODES        120
static const uint8_t code_to_plane_lut[CODE_TO_PLANE_CODES] = {
   0x18, 0x07, 0x17, 0x19, 0x28, 0x06, 0x27, 0x29, 0x16, 0x1a,
   0x26, 0x2a, 0x38, 0x05, 0x37, 0x39, 0x15, 0x1b, 0x36, 0x3a,
   0x25, 0x2b, 0x48, 0x04, 0x47, 0x49, 0x14, 0x1c, 0x35, 0x3b,
   0x46, 0x4a, 0x24, 0x2c, 0x58, 0x45, 0x4b, 0x34, 0x3c, 0x03,
   0x57, 0x59, 0x13, 0x1d, 0x56, 0x5a, 0x23, 0x2d, 0x44, 0x4c,
   0x55, 0x5b, 0x33, 0x3d, 0x68, 0x02, 0x67, 0x69, 0x12, 0x1e,
   0x66, 0x6a, 0x22, 0x2e, 0x54, 0x5c, 0x43, 0x4d, 0x65, 0x6b,
   0x32, 0x3e, 0x78, 0x01, 0x77, 0x79, 0x53, 0x5d, 0x11, 0x1f,
   0x64, 0x6c, 0x42, 0x4e, 0x76, 0x7a, 0x21, 0x2f, 0x75, 0x7b,
   0x31, 0x3f, 0x63, 0x6d, 0x52, 0x5e, 0x00, 0x74, 0x7c, 0x41,
   0x4f, 0x10, 0x20, 0x62, 0x6e, 0x30, 0x73, 0x7d, 0x51, 0x5f,
   0x40, 0x72, 0x7e, 0x61, 0x6f, 0x50, 0x71, 0x7f, 0x60, 0x70,
};


//------------------------------------------------------------------------------

static int DecodeImageStream(uint32_t xsize, uint32_t ysize,
                             VP8LDecoder* const dec,
                             argb_t** const decoded_data);

static int ReadImageSize(BitReader* const br,
                         uint32_t* const width, uint32_t* const height) {
  const int signature = VP8LReadBits(br, 8);
  if (signature != LOSSLESS_MAGIC_BYTE) return 0;
  *width = VP8LReadBits(br, kImageSizeBits) + 1;
  *height = VP8LReadBits(br, kImageSizeBits) + 1;
  return 1;
}

int VP8LGetInfo(const uint8_t* data, uint32_t data_size,
                int* width, int* height) {
  if (data_size >= kHeaderBytes) {
    uint32_t w, h;
    BitReader br;
    VP8LInitBitReader(&br, data, kHeaderBytes);
    if (ReadImageSize(&br, &w, &h)) {
      *width = w;
      *height = h;
      return 1;
    } else {
      return 0;  // Could not locate signature.
    }
  } else {
    return 0;         // not enough data
  }
}

static uint32_t SubSampleSize(uint32_t size, uint32_t sampling_bits) {
  return (size + (1 << sampling_bits) - 1) >> sampling_bits;
}

static inline uint32_t GetCopyDistance(uint32_t distance_symbol,
                                       BitReader* const br) {
  uint32_t extra_bits, offset;
  if (distance_symbol < 4) {
    return distance_symbol + 1;
  }
  extra_bits = (distance_symbol - 2) >> 1;
  offset = (2 + (distance_symbol & 1)) << extra_bits;
  return offset + VP8LReadBits(br, extra_bits) + 1;
}

static inline uint32_t GetCopyLength(uint32_t length_symbol,
                                     BitReader* const br) {
  // Length and distance prefixes are encoded the same way.
  return GetCopyDistance(length_symbol, br);
}

static int PlaneCodeToDistance(int xsize, uint32_t plane_code) {
  int dist_code, yoffset, xoffset;
  if (plane_code > CODE_TO_PLANE_CODES) {
    return plane_code - CODE_TO_PLANE_CODES;
  }
  dist_code = code_to_plane_lut[plane_code - 1];
  yoffset = dist_code >> 4;
  xoffset = 8 - (dist_code & 0xf);
  return yoffset * xsize + xoffset;
}

// Decodes the next Huffman code from bit-stream.
// FillBitWindow(br) needs to be called at minimum every second call
// to ReadSymbol.
static inline int ReadSymbol(const HuffmanTreeNode* root, BitReader* const br) {
  const HuffmanTreeNode* node = root;
  while (!HuffmanTreeNodeIsLeaf(node)) {
    node = node->child_[VP8LReadOneBitUnsafe(br)];
  }
  assert(node);
  assert(node->symbol_ >= 0);

  return node->symbol_;
}

static int ReadHuffmanCodeLengths(
    VP8LDecoder* const dec, const uint32_t* const code_length_code_lengths,
    uint32_t num_codes, uint32_t num_symbols, uint32_t** const code_lengths) {
  int ok = 1;
  BitReader* const br = dec->br_;
  uint32_t max_length = 0;
  uint32_t code_idx, sym_cnt;
  uint32_t* code_lengths_lcl = *code_lengths;
  int code_len, prev_code_len;
  const int use_length = VP8LReadBits(br, 1);
  HuffmanTreeNode root;

  HuffmanTreeNodeInit(&root);
  if (!HuffmanTreeBuild(&root, code_length_code_lengths, num_codes)) return 0;

  if (use_length) {
    const int length_nbits = (VP8LReadBits(br, 3) + 1) * 2;
    max_length = VP8LReadBits(br, length_nbits) + 2;
  }
  sym_cnt = 0;
  prev_code_len = 8;
  for (code_idx = 0; code_idx < num_symbols; ++code_idx) {
    if (use_length && ++sym_cnt > max_length) break;
    VP8LFillBitWindow(br);
    code_len = ReadSymbol(&root, br);
    if (code_len < kCodeLengthLiterals) {
      code_lengths_lcl[code_idx] = code_len;
      if (code_len != 0) prev_code_len = code_len;
    } else {
      const int repeat_ix = code_len - kCodeLengthLiterals;
      const int extra_bits = kCodeLengthExtraBits[repeat_ix];
      const int repeat_offset = kCodeLengthRepeatOffsets[repeat_ix];
      const int rep_cnt = VP8LReadBits(br, extra_bits) + repeat_offset;
      const int use_prev = (code_len == kCodeLengthRepeatCode);
      int rep_iter;
      assert(code_idx + rep_cnt <= num_symbols);
      for (rep_iter = 0; rep_iter < rep_cnt; ++rep_iter) {
        code_lengths_lcl[code_idx + rep_iter] = use_prev ? prev_code_len : 0;
      }
      code_idx += (rep_cnt - 1);
    }
  }

  return ok;
}

static int ReadHuffmanCode(int alphabet_size, VP8LDecoder* const dec,
                           HuffmanTreeNode* const root) {
  int ok = 1;
  BitReader* const br = dec->br_;
  const int simple_code = VP8LReadBits(br, 1);
  HuffmanTreeNodeInit(root);

  if (simple_code) {
    // For simple_code case, build the tree by calling 'AddSymbol' with simple
    // codes. Call to Huffman's BuildTree will assign symbol codes as per it's
    // symbol code logic and will not match the expected symbol codes intended
    // for this case.
    const int num_symbols = VP8LReadBits(br, 1) + 1;
    const int nbits = VP8LReadBits(br, 3);
    if (nbits == 0) {
      ok = HuffmanTreeAddSymbol(root, 0, 0, 0);
    } else {
      int sym_idx;
      const int num_bits = (nbits - 1) * 2 + 4;
      for (sym_idx = 0; ok && sym_idx < num_symbols; ++sym_idx) {
        ok = HuffmanTreeAddSymbol(root, VP8LReadBits(br, num_bits),
                                  num_symbols - 1, sym_idx);
      }
    }
  } else {
    uint32_t code_idx;
    uint32_t code_length_code_lengths[NUM_CODE_LENGTH_CODES] = { 0 };
    const uint32_t num_codes = VP8LReadBits(br, 4) + 4;
    uint32_t* code_lengths = NULL;

    if (num_codes > NUM_CODE_LENGTH_CODES) return 0;

    code_lengths = (uint32_t*)calloc(alphabet_size, sizeof(code_lengths[0]));
    if (code_lengths == NULL) {
      return 0;
    }

    for (code_idx = 0; code_idx < num_codes; ++code_idx) {
      code_length_code_lengths[kCodeLengthCodeOrder[code_idx]] =
          VP8LReadBits(br, 3);
    }
    ok = ReadHuffmanCodeLengths(dec, code_length_code_lengths,
                                NUM_CODE_LENGTH_CODES,
                                alphabet_size, &code_lengths);
    ok = HuffmanTreeBuild(root, code_lengths, alphabet_size);

    free(code_lengths);
  }


  return ok;
}

static int ReadHuffmanCodes(
    int xsize, int ysize, VP8LDecoder* const dec,
    int* const color_cache_bits, int* const color_cache_x_subsample_bits,
    uint32_t** const huffman_image, uint32_t* const huffman_subsample_bits,
    uint32_t** const meta_codes, uint32_t* const meta_code_size,
    HuffmanTreeNode** htrees, uint32_t* const num_huffman_trees) {
  int ok = 1;
  int use_color_cache, color_cache_size;
  uint32_t tree_idx;

  uint32_t* huffman_image_lcl = NULL;
  HuffmanTreeNode* htrees_lcl = NULL;
  uint32_t* meta_codes_lcl = NULL;

  BitReader* const br = dec->br_;
  const int use_meta_huff_codes = VP8LReadBits(br, 1);

  *num_huffman_trees = HUFFMAN_CODES_PER_META_CODE;
  if (use_meta_huff_codes) {
    uint32_t hpc, mc, hc, mci;
    uint32_t meta_codes_nbits, num_meta_codes, nbits;
    uint32_t huffman_xsize, huffman_ysize, huffman_pixs;

    *huffman_subsample_bits = VP8LReadBits(br, 4);
    huffman_xsize = SubSampleSize(xsize, *huffman_subsample_bits);
    huffman_ysize = SubSampleSize(ysize, *huffman_subsample_bits);
    huffman_pixs = huffman_xsize * huffman_ysize;
    if (!DecodeImageStream(huffman_xsize, huffman_ysize, dec,
                           &huffman_image_lcl)) return 0;
    for (hpc = 0; hpc < huffman_pixs; ++hpc) {
      // The huffman data is stored in red and green bytes.
      huffman_image_lcl[hpc] >>= 8;
      huffman_image_lcl[hpc] &= 0xffff;
    }
    meta_codes_nbits = VP8LReadBits(br, 4);
    num_meta_codes = VP8LReadBits(br, meta_codes_nbits) + 2;
    nbits = VP8LReadBits(br, 4);
    *meta_code_size = num_meta_codes * HUFFMAN_CODES_PER_META_CODE;
    meta_codes_lcl = (uint32_t*)calloc(
        *meta_code_size, sizeof(meta_codes_lcl[0]));
    if (meta_codes_lcl == NULL) return 0;

    mci = 0;
    for (mc = 0; mc < num_meta_codes; ++mc) {
      for (hc = 0; hc < HUFFMAN_CODES_PER_META_CODE; ++hc) {
        const uint32_t tree_index = VP8LReadBits(br, nbits);
        meta_codes_lcl[mci] = tree_index;
        if (*num_huffman_trees < tree_index + 1) {
          *num_huffman_trees = tree_index + 1;
        }
        ++mci;
      }
    }
  } else {
    int hc;
    *meta_code_size = HUFFMAN_CODES_PER_META_CODE;
    meta_codes_lcl = (uint32_t*)calloc(
        *meta_code_size, sizeof(meta_codes_lcl[0]));
    if (meta_codes_lcl == NULL) return 0;

    for (hc = 0; hc < HUFFMAN_CODES_PER_META_CODE; ++hc) {
      meta_codes_lcl[hc] = hc;
    }
  }

  use_color_cache = VP8LReadBits(br, 1);
  if (use_color_cache) {
    *color_cache_x_subsample_bits = VP8LReadBits(br, 4);
    *color_cache_bits = VP8LReadBits(br, 4);
    color_cache_size = 1 << *color_cache_bits;
  } else {
    *color_cache_x_subsample_bits = 0;
    *color_cache_bits = 0;
    color_cache_size = 0;
  }

  htrees_lcl = (HuffmanTreeNode *)calloc(
      *num_huffman_trees, sizeof(htrees_lcl[0]));
  if (htrees_lcl == NULL) {
    free(meta_codes_lcl);
    return 0;
  }

  for (tree_idx = 0; ok && tree_idx < *num_huffman_trees; ++tree_idx) {
    const int tree_type = tree_idx % HUFFMAN_CODES_PER_META_CODE;
    int alphabet_size = kAlphabetSize[tree_type];
    if (tree_type == 0) {
      alphabet_size += color_cache_size;
    }
    ok = ReadHuffmanCode(alphabet_size, dec, &htrees_lcl[tree_idx]);
  }

  *huffman_image = huffman_image_lcl;
  *meta_codes = meta_codes_lcl;
  *htrees = htrees_lcl;

  return ok;
}

static inline int GetMetaIndex(
    const uint32_t* const image, uint32_t xsize, uint32_t bits, int x, int y) {
  if (bits == 0) return 0;
  return image[xsize * (y >> bits) + (x >> bits)];
}

typedef HuffmanTreeNode* HuffmanTreeNodeArray[HUFFMAN_CODES_PER_META_CODE];

static void UpdateHuffmanSet(
    const uint32_t* const huffman_image, const uint32_t* const meta_codes,
    HuffmanTreeNode* const htrees, uint32_t huffman_xsize,
    uint32_t huffman_subsample_bits, uint32_t x, uint32_t y,
    int* const orig_meta_ix, HuffmanTreeNodeArray* const huffs) {
  const int meta_index = HUFFMAN_CODES_PER_META_CODE *
      GetMetaIndex(huffman_image, huffman_xsize, huffman_subsample_bits, x, y);

  if (*orig_meta_ix != meta_index) {
    HuffmanTreeNode** const huffs_lcl = *huffs;
    huffs_lcl[GREEN] = &htrees[meta_codes[meta_index + GREEN]];
    huffs_lcl[RED] = &htrees[meta_codes[meta_index + RED]];
    huffs_lcl[BLUE] = &htrees[meta_codes[meta_index + BLUE]];
    huffs_lcl[ALPHA] = &htrees[meta_codes[meta_index + ALPHA]];
    huffs_lcl[DIST] = &htrees[meta_codes[meta_index + DIST]];
    *orig_meta_ix = meta_index;
  }
}

static int DecodePixels(
    VP8LDecoder* const dec, uint32_t xsize, uint32_t ysize,
    int color_cache_bits, int color_cache_x_subsample_bits,
    const uint32_t* const huffman_image, uint32_t huffman_subsample_bits,
    const uint32_t* const meta_codes, HuffmanTreeNode* htrees,
    uint32_t** const decoded_data) {
  int ok = 1;
  int red, green, blue;
  int alpha = 0xff000000;
  int meta_ix = -1;
  uint32_t pos;
  uint32_t x = 0;
  uint32_t y = 0;
  uint32_t* data = NULL;
  BitReader* const br = dec->br_;
  // Collection of HUFFMAN_CODES_PER_META_CODE huffman codes.
  HuffmanTreeNode* huffs[HUFFMAN_CODES_PER_META_CODE] = { NULL };
  const int use_color_cache = (color_cache_bits > 0);
  const int color_cache_size = use_color_cache ? 1 << color_cache_bits : 0;
  // Values in range [NUM_CODES_PER_BYTE .. color_cache_limit[ are
  // color cache codes.
  const int color_cache_limit = NUM_CODES_PER_BYTE + color_cache_size;
  const int huffman_mask = (huffman_subsample_bits == 0) ?
      ~0 : (1 << huffman_subsample_bits) - 1;
  const uint32_t huffman_xsize = SubSampleSize(xsize, huffman_subsample_bits);
  VP8LColorCache* color_cache = NULL;
  if (use_color_cache) {
    color_cache = (VP8LColorCache*)malloc(sizeof(*color_cache));
    VP8LColorCacheInit(color_cache, xsize, color_cache_x_subsample_bits,
                       color_cache_bits);
  }

  data = (uint32_t*)calloc(xsize * ysize,sizeof(uint32_t));
  ok = (data != NULL);
  if (!ok) goto End;

  for (pos = 0; pos < xsize * ysize; ) {
    VP8LFillBitWindow(br);

    // Only update the huffman code when moving from one block to the next.
    if ((x & huffman_mask) == 0) {
      UpdateHuffmanSet(huffman_image, meta_codes, htrees, huffman_xsize,
                       huffman_subsample_bits, x, y, &meta_ix, &huffs);
    }

    green = ReadSymbol(huffs[GREEN], br);
    if (green < NUM_CODES_PER_BYTE) {
      // Literal.
      // Decode and save this pixel.
      red = ReadSymbol(huffs[RED], br);
      VP8LFillBitWindow(br);
      blue = ReadSymbol(huffs[BLUE], br);
      alpha = ReadSymbol(huffs[ALPHA], br);

      data[pos] = (alpha << 24) + (red << 16) + (green << 8) + blue;
      if (color_cache) VP8LColorCacheInsert(color_cache, x, data[pos]);

      // Update pos, x & y.
      ++pos; ++x;
      if (x == xsize) {
        ++y; x = 0;
      }
    } else if (green < color_cache_limit) {
      // Color cache.
      // Decode and save this pixel.
      const int color_cache_key = green - NUM_CODES_PER_BYTE;
      argb_t argb;
      ok = VP8LColorCacheLookup(color_cache, x, color_cache_key, &argb);
      if (!ok) goto End;
      data[pos] = argb;
      VP8LColorCacheInsert(color_cache, x, argb);

      // Update pos, x & y.
      ++pos; ++x;
      if (x == xsize) {
        ++y; x = 0;
      }
    } else if (green - color_cache_limit < NUM_LENGTH_CODES) {
      // Backward reference
      int dist_symbol;
      uint32_t i, dist_code, dist;
      int length_sym = green - color_cache_limit;
      const uint32_t length = GetCopyLength(length_sym, br);
      // Here, we have read the length code prefix + extra bits for the length,
      // so reading the next 15 bits can exhaust the bit window.
      // We must fill the window before the next read.
      VP8LFillBitWindow(br);
      dist_symbol = ReadSymbol(huffs[DIST], br);
      dist_code = GetCopyDistance(dist_symbol, br);
      dist = PlaneCodeToDistance(xsize, dist_code);
      assert(dist <= pos);
      assert(pos + length <= xsize * ysize);

      // Fill data for specified (backward-ref) length and update pos, x & y.
      if (color_cache) {
        for (i = 0; i < length; ++i) {
          data[pos] = data[pos - dist];
          VP8LColorCacheInsert(color_cache, x, data[pos]);
          ++pos; ++x;
          if (x == xsize) {
            ++y; x = 0;
          }
        }
      } else {
        for (i = 0; i < length; ++i) {
          data[pos] = data[pos - dist];
          ++pos;
        }
        x += length;
        while(x >= xsize) {
          x -= xsize;
          ++y;
        }
      }
      if (pos == xsize * ysize) {
        break;
      }

      UpdateHuffmanSet(huffman_image, meta_codes, htrees, huffman_xsize,
                       huffman_subsample_bits, x, y, &meta_ix, &huffs);

      continue;
    } else {
      // Code flow should not come here.
      assert(0);
    }
  }
 End:
  if (!ok) {
    free(data);
  } else {
    *decoded_data = data;
  }
  if (color_cache) {
    VP8LColorCacheRelease(color_cache);
    free(color_cache);
  }
  return ok;
}

static inline uint32_t Average2(uint32_t a0, uint32_t a1) {
  return (((a0 ^ a1) & 0xfefefefeL) >> 1) + (a0 & a1);
}

static inline uint32_t Average3(uint32_t a0, uint32_t a1, uint32_t a2) {
  return Average2(Average2(a0, a2), a1);
}

static inline uint32_t Average4(uint32_t a0, uint32_t a1,
                                uint32_t a2, uint32_t a3) {
  return Average2(Average2(a0, a1), Average2(a2, a3));
}

static uint32_t Add(uint32_t a, uint32_t b) {
  // This computes the sum of each component with mod 256.
  const uint32_t alpha_and_green = (a & 0xff00ff00) + (b & 0xff00ff00);
  const uint32_t red_and_blue = (a & 0x00ff00ff) + (b & 0x00ff00ff);
  return (alpha_and_green & 0xff00ff00) | (red_and_blue & 0x00ff00ff);
}

static inline uint32_t Clip255(uint32_t a) {
  if (a < NUM_CODES_PER_BYTE) {
    return a;
  }
  // return 0, when a is a negative integer.
  // return 255, when a is positive.
  return ~a >> 24;
}

static inline int AddSubtractComponentFull(int a, int b, int c) {
  return Clip255(a + b - c);
}

static uint32_t ClampedAddSubtractFull(uint32_t c0, uint32_t c1, uint32_t c2) {
  const int a = AddSubtractComponentFull(c0 >> 24, c1 >> 24, c2 >> 24);
  const int r = AddSubtractComponentFull((c0 >> 16) & 0xff,
                                         (c1 >> 16) & 0xff,
                                         (c2 >> 16) & 0xff);
  const int g = AddSubtractComponentFull((c0 >> 8) & 0xff,
                                         (c1 >> 8) & 0xff,
                                         (c2 >> 8) & 0xff);
  const int b = AddSubtractComponentFull(c0 & 0xff, c1 & 0xff, c2 & 0xff);
  return (a << 24) | (r << 16) | (g << 8) | b;
}

static inline int AddSubtractComponentHalf(int a, int b) {
  return Clip255(a + (a - b) / 2);
}

static uint32_t ClampedAddSubtractHalf(uint32_t c0, uint32_t c1, uint32_t c2) {
  const uint32_t ave = Average2(c0, c1);
  const int a = AddSubtractComponentHalf(ave >> 24, c2 >> 24);
  const int r = AddSubtractComponentHalf((ave >> 16) & 0xff, (c2 >> 16) & 0xff);
  const int g = AddSubtractComponentHalf((ave >> 8) & 0xff, (c2 >> 8) & 0xff);
  const int b = AddSubtractComponentHalf((ave >> 0) & 0xff, (c2 >> 0) & 0xff);
  return (a << 24) | (r << 16) | (g << 8) | b;
}

static uint32_t Select(uint32_t a, uint32_t b, uint32_t c) {
  const int p0 = (int)(a >> 24) + (int)(b >> 24) - (int)(c >> 24);
  const int p1 = (int)((a >> 16) & 0xff) + (int)((b >> 16) & 0xff) -
      (int)((c >> 16) & 0xff);
  const int p2 = (int)((a >> 8) & 0xff) + (int)((b >> 8) & 0xff) -
      (int)((c >> 8) & 0xff);
  const int p3 = (int)(a & 0xff) + (int)(b & 0xff) - (int)(c & 0xff);
  const int pa = abs(p0 - (a >> 24)) +
      abs(p1 - ((a >> 16) & 0xff)) +
      abs(p2 - ((a >> 8) & 0xff)) +
      abs(p3 - (a & 0xff));
  const int pb = abs(p0 - (b >> 24)) +
      abs(p1 - ((b >> 16) & 0xff)) +
      abs(p2 - ((b >> 8) & 0xff)) +
      abs(p3 - (b & 0xff));
  if (pa <= pb) {
    return a;
  } else {
    return b;
  }
}

#define ARGB_BLACK 0xff000000

static argb_t PredictValue(uint32_t pred_mode, int xsize,
                           const argb_t* const argb) {
  switch(pred_mode) {
    case 0: return ARGB_BLACK;
    case 1: return argb[-1];
    case 2: return argb[-xsize];
    case 3: return argb[-xsize + 1];
    case 4: return argb[-xsize - 1];
    case 5: return Average3(argb[-1], argb[-xsize], argb[-xsize + 1]);
    case 6: return Average2(argb[-1], argb[-xsize - 1]);
    case 7: return Average2(argb[-1], argb[-xsize]);
    case 8: return Average2(argb[-xsize - 1], argb[-xsize]);
    case 9: return Average2(argb[-xsize], argb[-xsize + 1]);
    case 10: return Average4(argb[-1], argb[-xsize - 1],
                             argb[-xsize], argb[-xsize + 1]);
    case 11: return Select(argb[-xsize], argb[-1], argb[-xsize - 1]);
    case 12:
      return ClampedAddSubtractFull(argb[-1], argb[-xsize], argb[-xsize - 1]);
    case 13:
      return ClampedAddSubtractHalf(argb[-1], argb[-xsize], argb[-xsize - 1]);
  }
  return 0;
}

static void PredictorInverseTransform(const VP8LTransform* const transform,
                                      argb_t* const decoded_data) {
  size_t row, col;
  uint32_t image_ix = 0;
  const uint32_t tile_mask = (1 << transform->bits_) - 1;
  const uint32_t tile_xsize =
      (transform->xsize_ + tile_mask) >> transform->bits_;

  uint32_t pred_mode = 0;
  argb_t pred = 0;

  // First Row follows the L (mode=1) mode.
  decoded_data[0] = Add(decoded_data[0], ARGB_BLACK);
  for (col = 1; col < transform->xsize_; ++col) {
    decoded_data[col] = Add(decoded_data[col], decoded_data[col - 1]);
  }
  image_ix += transform->xsize_;

  for (row = 1; row < transform->ysize_; ++row) {
    const uint32_t tile_y = row >> transform->bits_;
    uint32_t tile_x = 0;
    for (col = 0; col < transform->xsize_; ++col, ++image_ix) {
      // Pick the appropriate predictor mode (at start of every tile).
      if (!(col & tile_mask)) {
        const int tile_ix = tile_y * tile_xsize + tile_x;
        pred_mode = (transform->data_[tile_ix] >> 8) & 0xff;
        ++tile_x;
      }
      // First col follows the T (mode=2) mode.
      pred = (col == 0) ? decoded_data[image_ix - transform->xsize_] :
          PredictValue(pred_mode, transform->xsize_, decoded_data + image_ix);
      decoded_data[image_ix] = Add(decoded_data[image_ix], pred);
    }
  }
}

static void AddGreenToBlueAndRed(const VP8LTransform* const transform,
                                 argb_t* const decoded_data) {
  int i;
  int num_pixs = transform->xsize_ * transform->ysize_;
  for (i = 0; i < num_pixs; ++i) {
    const argb_t argb = decoded_data[i];
    argb_t green = (argb >> 8) & 0xff;
    argb_t red_blue = argb & 0x00ff00ff;
    red_blue += ((green << 16) + green);
    red_blue &= 0x00ff00ff;
    decoded_data[i] = (argb & 0xff00ff00) + red_blue;
  }
}

static inline uint32_t ColorTransformDelta(signed char color_pred,
                                           signed char color) {
  return (uint32_t)((int)(color_pred) * color) >> 5;
}

static argb_t TransformColor(uint32_t color_pred, argb_t argb) {
  const uint32_t green_to_red = color_pred & 0xff;
  const uint32_t green_to_blue = (color_pred >> 8) & 0xff;
  const uint32_t red_to_blue = (color_pred >> 16) & 0xff;

  const uint32_t green = argb >> 8;
  const uint32_t red = argb >> 16;
  uint32_t new_red = red;
  uint32_t new_blue = argb;
  new_red += ColorTransformDelta(green_to_red, green);
  new_red &= 0xff;
  new_blue += ColorTransformDelta(green_to_blue, green);
  new_blue += ColorTransformDelta(red_to_blue, new_red);
  new_blue &= 0xff;
  return (argb & 0xff00ff00) | (new_red << 16) | (new_blue);
}

static void ColorSpaceInverseTransform(const VP8LTransform* const transform,
                                       argb_t* const decoded_data) {
  size_t row, col;
  uint32_t image_ix = 0;
  const uint32_t tile_mask = (1 << transform->bits_) - 1;
  const uint32_t tile_xsize =
      (transform->xsize_ + tile_mask) >> transform->bits_;
  argb_t color_pred = 0;

  for (row = 0; row < transform->ysize_; ++row) {
    const uint32_t tile_y = row >> transform->bits_;
    uint32_t tile_x = 0;
    for (col = 0; col < transform->xsize_; ++col, ++image_ix) {
      // Pick the appropriate color predictor mode (at start of every tile).
      if (!(col & tile_mask)) {
        const int tile_ix = tile_y * tile_xsize + tile_x;
        color_pred = transform->data_[tile_ix];
        ++tile_x;
      }
      decoded_data[image_ix] = TransformColor(color_pred,
                                              decoded_data[image_ix]);
    }
  }
}

static void ColorIndexingInverseTransform(const VP8LTransform* const transform,
                                          argb_t* const decoded_data) {
  int i;
  const int num_pixs = transform->xsize_ * transform->ysize_;
  for (i = 0; i < num_pixs; ++i) {
    decoded_data[i] = transform->data_[(decoded_data[i] >> 8) & 0xff];
  }
}

static int PixelBundleInverseTransform(const VP8LTransform* const transform,
                                       argb_t** const decoded_data) {
  uint32_t row, col, tile_x;
  const uint32_t bit_depth = 8 >> transform->bits_;
  const uint32_t num_cols = 1 << transform->bits_;
  const uint32_t bit_mask = num_cols - 1;
  const uint32_t xs = SubSampleSize(transform->xsize_, transform->bits_);

  uint32_t* tmp = (uint32_t*)calloc(
      transform->xsize_ * transform->ysize_, sizeof(*tmp));
  if (tmp == NULL) return 0;

  for (row = 0; row < transform->ysize_; ++row) {
    for (tile_x = 0; tile_x < xs; ++tile_x) {
      const argb_t* const rowp = (*decoded_data) + row * xs;
      uint32_t tile_code = (rowp[tile_x] >> 8) & 0xff;
      for (col = 0; col < num_cols; ++col) {
        const uint32_t x_all = tile_x * num_cols + col;
        if (x_all < transform->xsize_) {
          const uint32_t ix = row * transform->xsize_ + x_all;
          const uint32_t green = tile_code & bit_mask;
          tmp[ix] = 0xff000000 | (green << 8);
          tile_code >>= bit_depth;
        }
      }
    }
  }
  free(*decoded_data);
  *decoded_data = tmp;

  return 1;
}

static int ApplyInverseTransforms(VP8LDecoder* const dec, int start_idx,
                                  argb_t** const decoded_data) {
  int ok = 1;
  int n = dec->next_transform_;
  assert(start_idx >= 0);
  while(ok && n > start_idx) {
    const VP8LTransform* const transform = &(dec->transforms_[--n]);
    switch (transform->type_) {
      case SUBTRACT_GREEN:
        AddGreenToBlueAndRed(transform, *decoded_data);
        break;
      case PREDICTOR_TRANSFORM:
        PredictorInverseTransform(transform, *decoded_data);
        break;
      case CROSS_COLOR_TRANSFORM:
        ColorSpaceInverseTransform(transform, *decoded_data);
        break;
      case COLOR_INDEXING_TRANSFORM:
        ColorIndexingInverseTransform(transform, *decoded_data);
        break;
      case PIXEL_BUNDLE_TRANSFORM:
        ok = PixelBundleInverseTransform(transform, decoded_data);
        break;
      default:
        ok = 0;
        break;
    }
  }
  dec->next_transform_ = n;

  return ok;
}

static int ReadTransform(int* const xsize, int* const ysize,
                         VP8LDecoder* const dec) {
  int ok = 1;
  BitReader* const br = dec->br_;
  VP8LTransform* transform = &(dec->transforms_[dec->next_transform_]);
  const ImageTransformType type = (ImageTransformType)VP8LReadBits(br, 3);

  transform->type_ = type;
  transform->xsize_ = *xsize;
  transform->ysize_ = *ysize;
  transform->data_ = NULL;

  switch (type) {
    case PREDICTOR_TRANSFORM:
    case CROSS_COLOR_TRANSFORM:
      transform->bits_ = VP8LReadBits(br, 4);
      ok = DecodeImageStream(SubSampleSize(transform->xsize_, transform->bits_),
                             SubSampleSize(transform->ysize_, transform->bits_),
                             dec, &transform->data_);
      break;
    case COLOR_INDEXING_TRANSFORM:
      {
        const int num_colors = VP8LReadBits(br, 8) + 1;
        ok = DecodeImageStream(num_colors, 1, dec, &transform->data_);
        if (ok) {
          int i;
          for (i = 1; i < num_colors; ++i) {
            transform->data_[i] = Add(transform->data_[i] ,
                                      transform->data_[i - 1]);
          }
        }
      }
      break;
    case PIXEL_BUNDLE_TRANSFORM:
      transform->bits_ = VP8LReadBits(br, 2);
      transform->data_ = NULL;
      *xsize = SubSampleSize(transform->xsize_, transform->bits_);
      break;
    case SUBTRACT_GREEN:
      break;
    default:
      ok = 0;
  }

  if (ok) ++dec->next_transform_;

  return ok;
}

static int DecodeImageStream(uint32_t xsize, uint32_t ysize,
                             VP8LDecoder* const dec,
                             argb_t** const decoded_data) {
  int transform_xsize = xsize;
  int transform_ysize = ysize;
  int ok = 1;
  uint32_t tree_idx;
  uint32_t* huffman_image = NULL;
  uint32_t* meta_codes = NULL;
  uint32_t huffman_subsample_bits = 0;
  uint32_t num_meta_codes = 0;
  uint32_t num_huffman_trees = 0;
  int color_cache_bits = 0;
  int color_cache_x_subsample_bits = 0;

  HuffmanTreeNode* htrees = NULL;
  BitReader* const br = dec->br_;
  int transform_start_idx = dec->next_transform_;

  // Step#1: Read the transforms.
  while(ok && VP8LReadBits(br, 1)) {
    ok = ReadTransform(&transform_xsize, &transform_ysize, dec);
  }

  // Step#2: Read the Huffman codes.
  ok = ok && ReadHuffmanCodes(transform_xsize, transform_ysize, dec,
                              &color_cache_bits, &color_cache_x_subsample_bits,
                              &huffman_image, &huffman_subsample_bits,
                              &meta_codes, &num_meta_codes,
                              &htrees, &num_huffman_trees);

  // Step#3: Use the Huffman trees to decode the LZ77 encoded data.
  ok = ok && DecodePixels(dec, xsize, ysize,
                          color_cache_bits, color_cache_x_subsample_bits,
                          huffman_image, huffman_subsample_bits, meta_codes,
                          htrees, decoded_data);

  free(huffman_image);
  free(meta_codes);
  for (tree_idx = 0; tree_idx < num_huffman_trees; ++tree_idx) {
    HuffmanTreeRelease(&htrees[tree_idx]);
  }
  free(htrees);

  // Step#4: Apply transforms on the decoded data.
  ok = ok && ApplyInverseTransforms(dec, transform_start_idx, decoded_data);

  return ok;
}

int VP8LDecodeImage(VP8LDecoder* const dec, VP8Io* const io, uint32_t offset) {
  uint32_t width, height;
  argb_t* decoded_data = NULL;
  BitReader br;
  assert(dec);

  if (offset > io->data_size) return 0;

  VP8LInitBitReader(&br, io->data + offset, io->data_size - offset);
  if (!ReadImageSize(&br, &width, &height)) return 0;
  dec->br_ = &br;
  dec->next_transform_ = 0;
  if (!DecodeImageStream(width, height, dec, &decoded_data)) return 0;
  dec->argb_ = decoded_data;

  return 1;
}
//------------------------------------------------------------------------------

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif
