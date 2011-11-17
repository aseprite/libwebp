// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: szabadka@google.com (Zoltan Szabadka)

#include "decode.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <map>
#include <string>

#include "bit_stream.h"
#include "constants.h"
#include "huffman.h"
#include "image_transform.h"
#include "../common/integral_types.h"
#include "../common/pixel_hasher.h"

static const int kMaxImageTransforms = 100;

void DecodeImageInternal(const int xsize,
                         const int ysize,
                         BitStream* stream,
                         uint32** argb_image);

namespace {

int MetaSize(int size, int bits) {
  return (size + (1 << bits) - 1) >> bits;
}

int AlphabetSize(int type, int num_green, int palette_size) {
  if (type == 0) return num_green + palette_size + kNumLengthSymbols;
  if (type == 1) return 256;
  if (type == 2) return kNumDistanceSymbols;
  return 0;
}

int GetMetaIndex(int xsize, int bits,
                 const uint32* image,
                 int x, int y) {
  if (bits == 0) return 0;
  int xs = (xsize + (1 << bits) - 1) >> bits;
  return (image[(y >> bits) * xs +(x >> bits)] >> 8) & 0xffff;
}

void ReadMetaCodes(const int num_rba,
                   BitStream* stream,
                   std::vector<int>* meta_codes,
                   std::map<int, int>* tree_types,
                   int* num_trees) {
  int meta_codes_nbits = stream->Read(4);
  int num_meta_codes = stream->Read(meta_codes_nbits) + 2;
  int nbits = stream->Read(4);
  *num_trees = 0;
  for (int i = 0; i < num_meta_codes; ++i) {
    for (int k = 0; k < num_rba + 2; ++k) {
      int tree_index = stream->Read(nbits);
      meta_codes->push_back(tree_index);
      (*tree_types)[tree_index] = (k == 0) ? 0 : (k == num_rba + 1) ? 2 : 1;
      *num_trees = std::max(tree_index + 1, *num_trees);
    }
  }
}

int GetCopyDistance(int distance_symbol, BitStream* stream) {
  if (distance_symbol < 4) {
    return distance_symbol + 1;
  }
  int extra_bits = (distance_symbol - 2) / 2;
  int offset =
      (1 << (extra_bits + 1)) + (distance_symbol % 2) * (1 << extra_bits) + 1;
  return offset + stream->Read(extra_bits);
}

int PlaneCodeToDistance(int xsize, int ysize, int plane_code) {
  if (plane_code > 120) {
    return plane_code - 120;
  }
  int dist_code = code_to_plane_lut[plane_code - 1];
  int yoffset = dist_code / 16;
  int xoffset = dist_code % 16;
  if (xoffset <= 8) {
    xoffset = 8 - xoffset;
  } else {
    xoffset = xsize - (xoffset - 8);
    yoffset -= 1;
  }
  return yoffset * xsize + xoffset;
}

int GetCopyLength(int length_symbol, BitStream* stream) {
  return GetCopyDistance(length_symbol, stream);
}

std::string CodeLengthDebugString(const std::vector<int>& code_lengths) {
  std::string out = "Code Lengths: ";
  for (int i = 0; i < code_lengths.size(); ++i) {
    if (code_lengths[i] > 0) {
      char str[32];
      sprintf(str, " %d:%d", i, code_lengths[i]);
      out.append(std::string(&str[0]));
    }
  }
  return out;
}

int ReadSymbol(const HuffmanTreeNode& root,
               BitStream* stream) {
  const HuffmanTreeNode* node = &root;
  while (!node->IsLeaf()) node = node->child(stream->Read(1));
  VERIFY(node);
  VERIFY(node->symbol() != -1);
  return node->symbol();
}

void ReadCodeLengthTree(BitStream* stream,
                        HuffmanTreeNode* root) {
  std::vector<int> code_lengths(kCodeLengthCodes);
  int hclen = stream->Read(4);
  for (int i = 0; i < hclen + 4; ++i) {
    code_lengths[kCodeLengthCodeOrder[i]] = stream->Read(3);
  }
  if(!root->BuildTree(code_lengths)) {
    printf("error: %s\n", CodeLengthDebugString(code_lengths).c_str());
    abort();
  }
}

void ReadHuffmanCodeLengths(const HuffmanTreeNode& decoder_root,
                            BitStream* stream,
                            std::vector<int>* code_lengths) {
  bool use_length = stream->Read(1);
  int max_length = 0;
  if (use_length) {
    int length_nbits = (stream->Read(3) + 1) * 2;
    max_length = stream->Read(length_nbits) + 2;
  }
  int previous = 8;
  int num_symbols = 0;
  for (int i = 0; i < code_lengths->size(); ++i) {
    if (use_length && ++num_symbols > max_length) break;
    int code_length = ReadSymbol(decoder_root, stream);
    VERIFY(code_length < kCodeLengthCodes);
    if (code_length < kCodeLengthLiterals) {
      (*code_lengths)[i] = code_length;
      if (code_length != 0) previous = code_length;
    } else {
      int repeat_ix = code_length - kCodeLengthLiterals;
      int extra_bits = kCodeLengthExtraBits[repeat_ix];
      int repeat_offset = kCodeLengthRepeatOffsets[repeat_ix];
      int repeat = stream->Read(extra_bits) + repeat_offset;
      bool use_previous = code_length == kCodeLengthRepeatCode;
      for (int k = 0; k < repeat; ++k) {
        (*code_lengths)[i + k] = use_previous ? previous : 0;
      }
      i += repeat - 1;
    }
  }
}

void ReadHuffmanCode(const int alphabet_size,
                     BitStream* stream,
                     HuffmanTreeNode* root) {
  bool simple_code = stream->Read(1);
  if (simple_code) {
    int num_symbols = stream->Read(1) + 1;
    int nbits = stream->Read(3);
    if (nbits == 0) {
      root->AddSymbol(0, 0, 0);
    } else {
      int num_bits = (nbits - 1) * 2 + 4;
      VERIFY(num_bits <= 12);
      for (int i = 0; i < num_symbols; ++i) {
        root->AddSymbol(stream->Read(num_bits), num_symbols - 1, i);
      }
    }
    VERIFY(root->IsFull());
  } else {
    HuffmanTreeNode decoder_root;
    ReadCodeLengthTree(stream, &decoder_root);
    std::vector<int> code_lengths(alphabet_size);
    ReadHuffmanCodeLengths(decoder_root, stream, &code_lengths);
    if (!root->BuildTree(code_lengths)) {
      printf("error #2: %s\n", CodeLengthDebugString(code_lengths).c_str());
      abort();
    }
  }
}

void ReadTransform(int* xsize, int* ysize,
                   ImageTransform* transform,
                   BitStream* stream) {
  transform->type = (ImageTransformType)stream->Read(3);
  transform->xsize = *xsize;
  transform->ysize = *ysize;
  switch (transform->type) {
    case PREDICTOR_TRANSFORM:
    case CROSS_COLOR_TRANSFORM:
      {
        transform->data = malloc(sizeof(PerTileTransformData));
        PerTileTransformData* data = (PerTileTransformData*)transform->data;
        data->bits = stream->Read(4);
        DecodeImageInternal(MetaSize(*xsize, data->bits),
                            MetaSize(*ysize, data->bits),
                            stream, &data->image);
      }
      break;
    case SUBTRACT_GREEN:
      // There is no associated data for this transform.
      transform->data = NULL;
      break;
    case COMPONENT_SUBSAMPLING_TRANSFORM:
      {
        transform->data = malloc(4 * sizeof(int));
        int* bits = (int*)transform->data;
        for (int i = 0; i < 4; ++i) {
          bits[i] = stream->Read(3);
        }
      }
      break;
    case COLOR_INDEXING_TRANSFORM:
      {
        int ncolors = stream->Read(8) + 1;
        transform->data = malloc(sizeof(PerTileTransformData));
        PerTileTransformData* data = (PerTileTransformData*)transform->data;
        DecodeImageInternal(ncolors, 1, stream, &data->image);
        for (int i = 1; i < ncolors; ++i) {
          data->image[i] = Add(data->image[i], data->image[i - 1]);
        }
      }
      break;
    case PIXEL_BUNDLE_TRANSFORM:
      {
        transform->data = malloc(sizeof(PixelBundleTransformData));
        PixelBundleTransformData* data =
            (PixelBundleTransformData*)transform->data;
        data->xbits = stream->Read(3);
        data->ybits = stream->Read(3);
        data->bit_depth = stream->Read(3) + 1;
        *xsize = MetaSize(*xsize, data->xbits);
        *ysize = MetaSize(*ysize, data->ybits);
      }
      break;
    case IMPLICIT_ALPHA_TRANSFORM:
      {
        transform->data = malloc(sizeof(uint32));
        uint32* argb = (uint32*)transform->data;
        *argb = 0;
        for (int k = 0; k < 4; ++k) {
          *argb |= stream->Read(8) << (8 * k);
        }
      }
      break;
  }
}

}  // namespace

void DecodeImageInternal(const int original_xsize,
                         const int original_ysize,
                         BitStream* stream,
                         uint32** argb_image) {
  int xsize = original_xsize;
  int ysize = original_ysize;

  ImageTransform* transforms =
      (ImageTransform*)malloc(kMaxImageTransforms * sizeof(ImageTransform));
  int num_transforms = 0;
  while (stream->Read(1)) {
    ReadTransform(&xsize, &ysize, &transforms[num_transforms], stream);
    ++num_transforms;
    VERIFY(num_transforms < kMaxImageTransforms);
  }

  bool error_detection_bits = stream->Read(1);
  if (error_detection_bits) {
    VERIFY(kMagicByteForErrorDetection == stream->Read(8));
  }

  int num_rba = stream->Read(2);

  bool use_meta_codes = stream->Read(1);
  int huffman_bits = 0;
  uint32* huffman_image;
  std::vector<int> meta_codes;
  std::map<int, int> tree_types;
  int num_huffman_trees = num_rba + 2;
  if (use_meta_codes) {
    huffman_bits = stream->Read(4);
    DecodeImageInternal(MetaSize(xsize, huffman_bits),
                        MetaSize(ysize, huffman_bits),
                        stream, &huffman_image);
    ReadMetaCodes(num_rba, stream, &meta_codes, &tree_types,
                  &num_huffman_trees);
    if (error_detection_bits) {
      VERIFY(kMagicByteForErrorDetection == stream->Read(8));
    }
  } else{
    for (int k = 0; k < num_rba + 2; ++k) {
      meta_codes.push_back(k);
      tree_types[k] = (k == 0) ? 0 : (k == num_rba + 1) ? 2 : 1;
    }
  }

  bool use_palette = stream->Read(1);
  int palette_x_bits = use_palette ? stream->Read(4) : 0;
  int palette_code_bits = use_palette ? stream->Read(4) : 0;
  int palette_size = use_palette ? 1 << palette_code_bits : 0;
  PixelHasherLine* palette = use_palette ?
      new PixelHasherLine(xsize, palette_x_bits, palette_code_bits) : NULL;

  int green_bit_depth = stream->Read(3) + 1;
  int num_green = 1 << green_bit_depth;

  std::vector<HuffmanTreeNode> htrees(num_huffman_trees);
  for (int i = 0; i < htrees.size(); ++i) {
    int type = tree_types[i];
    int alphabet_size = AlphabetSize(type, num_green, palette_size);
    ReadHuffmanCode(alphabet_size, stream, &htrees[i]);
  }

  if (error_detection_bits) {
    VERIFY(kMagicByteForErrorDetection == stream->Read(8));
  }

  *argb_image = (uint32*)malloc(xsize * ysize * sizeof(uint32));
  bool log = false;
  if (log) printf("\nDecoding actual bit stream.\n");
  for (int pos = 0; pos < xsize * ysize; ) {
    int x = pos % xsize;
    int y = pos / xsize;
    int meta_index = GetMetaIndex(xsize, huffman_bits, huffman_image, x, y);
    int meta_ix = meta_index * (num_rba + 2);

    if (log) printf("x = %d y = %d meta_index = %d  ", x, y, meta_index);

    int green = ReadSymbol(htrees[meta_codes[meta_ix]],
                           stream);

    // Literal
    if (green < num_green) {
      int red = (num_rba > 0) ?
          ReadSymbol(htrees[meta_codes[meta_ix + 1]], stream) : 0;
      int blue = (num_rba > 1) ?
          ReadSymbol(htrees[meta_codes[meta_ix + 2]], stream) : 0;
      int alpha = (num_rba > 2) ?
          ReadSymbol(htrees[meta_codes[meta_ix + 3]], stream) : 255;

      if (log) printf("green=0x%02x red=0x%02x blue=0x%02x alpha=0x%02x\n",
                      green, red, blue, alpha);

      uint32 argb = (alpha << 24) + (red << 16) + (green << 8) + blue;
      (*argb_image)[pos] = argb;
      if (palette) palette->Insert(x, argb);
      ++pos;
      continue;
    }

    // Palette
    int palette_symbol = green - num_green;
    if (palette_symbol < palette_size) {
      if (log) printf("palette_symbol = %d\n", palette_symbol);
      VERIFY(palette);
      uint32 argb = palette->Lookup(pos % xsize, palette_symbol);
      (*argb_image)[pos] = argb;
      palette->Insert(x, argb);
      ++pos;
      continue;
    }

    // Backward reference
    int length_symbol = palette_symbol - palette_size;
    if (length_symbol < kNumLengthSymbols) {
      if (log) printf("length_symbol = %d ", length_symbol);
      int length = GetCopyLength(length_symbol, stream);
      if (log) printf("length = %d ", length);
      int dist_symbol = ReadSymbol(htrees[meta_codes[meta_ix + num_rba + 1]],
                                   stream);
      if (log) printf("dist_symbol = %d ", dist_symbol);
      int dist = GetCopyDistance(dist_symbol, stream);
      if (log) printf("plane dist = %d ", dist);
      dist = PlaneCodeToDistance(xsize, ysize, dist);
      if (log) printf("dist = %d\n", dist);
      VERIFY(dist > 0);
      VERIFY(dist <= pos);
      VERIFY(pos + length <= xsize * ysize);
      for (int i = 0; i < length; ++i) {
        (*argb_image)[pos] = (*argb_image)[pos - dist];
        if (palette) palette->Insert(pos % xsize, (*argb_image)[pos]);
        ++pos;
      }
      continue;
    }
    printf("Error: Could not interpret GREEN symbol %d\n", green);
    abort();
  }
  if (use_meta_codes) {
    free(huffman_image);
  }
  if (palette) {
    delete palette;
  }

  if (error_detection_bits) {
    VERIFY(kMagicByteForErrorDetection == stream->Read(8));
  }

  for (int i = num_transforms - 1; i >= 0; --i) {
    ApplyInverseImageTransform(transforms[i], argb_image);
    FreeImageTransformData(transforms[i]);
  }
  free(transforms);
}

bool DecodeWebpLLImage(int encoded_image_size,
                       const uint8* const encoded_image,
                       int* xsize,
                       int* ysize,
                       uint32** argb_image) {
  BitStream stream(encoded_image_size, encoded_image);
  int size_bits = (stream.Read(3) + 1) * 4;
  *xsize = stream.Read(size_bits);
  *ysize = stream.Read(size_bits);
  DecodeImageInternal(*xsize, *ysize, &stream, argb_image);
  return true;
}
