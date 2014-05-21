#if defined(WEBP_USE_AVX2)
#include <stdlib.h>  // for abs()
#include <immintrin.h>

//------------------------------------------------------------------------------
// Quantization
//

#define QFIX2 0
static WEBP_INLINE int DoQuantizeBlock(int16_t in[16], int16_t out[16],
                                       int shift,
                                       const uint16_t* const sharpen,
                                       const VP8Matrix* const mtx) {
  const __m256i max_coeff_2047 = _mm256_set1_epi16(MAX_LEVEL);
  const __m256i zero = _mm256_setzero_si256();
  __m256i coeff;
  __m256i out_reg;
  __m128i packed_out;

  __m256i in_reg = _mm256_load_si256((__m256i*)&in[0]);
  const __m256i iq_reg = _mm256_load_si256((__m256i*)&mtx->iq_[0]);
  const __m256i q_reg = _mm256_load_si256((__m256i*)&mtx->q_[0]);

  // extract sign(in)  (0x0000 if positive, 0xffff if negative)
  const __m256i sign_reg = _mm256_cmpgt_epi16(zero, in_reg);

  // coeff = abs(in)
  coeff = _mm256_abs_epi16(in_reg);

  // coeff = abs(in) + sharpen
  if (sharpen != NULL) {
    const __m256i sharpen_reg = _mm256_loadu_si256((__m256i*)&sharpen[0]);
    coeff = _mm256_add_epi16(coeff0, sharpen_reg);
  }

  // out = (coeff * iQ + B) >> (QFIX + QFIX2 - shift)
  {
    // doing calculations with 32b precision (QFIX=17)
    // out = (coeff * iQ)
    const __m256i coeff_iQREGH = _mm256_mulhi_epu16(coeff, iq_reg);
    const __m256i coeff_iQREGL = _mm256_mullo_epi16(coeff, iq_reg);

    __m256i out_00 = _mm256_unpacklo_epi16(coeff_iQREGL, coeff_iQREGH);
    __m256i out_08 = _mm256_unpackhi_epi16(coeff_iQREGL, coeff_iQERGH);

    // out = (coeff * iQ + B)
    const __m256i bias_00 = _mm256_loadu_si256((__m256i*)&mtx->bias_[0]);
    const __m256i bias_08 = _mm256_loadu_si256((__m256i*)&mtx->bias_[8]);

    out_00 = _mm256_add_epi32(out_00, bias_00);
    out_08 = _mm256_add_epi32(out_08, bias_08);

    // out = QUANTDIV(coeff, iQ, B, QFIX + QFIX2 - shift)
    out_00 = _mm256_srai_epi32(out_00, QFIX + QFIX2 - shift);
    out_08 = _mm256_srai_epi32(out_08, QFIX + QFIX2 - shift);

    // pack result as 16b
    out_reg = _mm256_packs_epi32(out_00, out_08);

    out_reg = _mm256_permute4x64_epi64(out_reg, 216);
    // if (coeff > 2047) coeff = 2047
    out_reg = _mm256_min_epi16(out_reg, max_coeff_2047);
  }

  // get sign back (if (sign[j]) out_n = -out_n)
  out_reg = _mm256_xor_si256(out_reg, sign_reg);
  out_reg = _mm256_sub_epi16(out_reg, sign_reg);

  // in = out * Q
  in_reg = _mm256_mullo_epi16(out_reg, q_reg);

  _mm256_storeu_si256((__m256i*)&in[0], in_reg);

  // zigzag the output before storing it.
  //
  // The zigzag pattern can almost be reproduced with a small sequence of
  // shuffles. After it, we only need to swap the 7th (ending up in third
  // position instead of twelfth) and 8th values.
  {
    __m256i outZ0, outZ8;
    __m128i outZ0_128b, outZ8_128b;
    outZ0 = _mm256_shufflehi_epi16(out_reg,  _MM_SHUFFLE(2, 1, 3, 0));
    outZ0 = _mm256_shuffle_epi32  (outZ0, _MM_SHUFFLE(3, 1, 2, 0));
    outZ0 = _mm256_shufflehi_epi16(outZ0, _MM_SHUFFLE(3, 1, 0, 2));
    outZ8 = _mm256_shufflelo_epi16(out_reg,  _MM_SHUFFLE(3, 0, 2, 1));
    outZ8 = _mm256_shuffle_epi32  (outZ8, _MM_SHUFFLE(3, 1, 2, 0));
    outZ8 = _mm256_shufflelo_epi16(outZ8, _MM_SHUFFLE(1, 3, 2, 0));
    outZ0_128b = _mm256_castsi256_si128(outZ0);
    _mm_storeu_si128((__m256i*)&out[0], outZ0_128b);
    outZ8_128b = _mm256_extractf128_si256(outZ8, 1);
    _mm_storeu_si128((__m128i*)&out[8], outZ8_128b);
    packed_out = _mm256_packs_epi16(outZ0_128b, outZ8_128b);
  }
  {
    const int16_t outZ_12 = out[12];
    const int16_t outZ_3 = out[3];
    out[3] = outZ_12;
    out[12] = outZ_3;
  }

  // detect if all 'out' values are zeroes or not
  return (_mm_movemask_epi8(_mm_cmpeq_epi8(packed_out, zero)) != 0xffff);
}

static int QuantizeBlock(int16_t in[16], int16_t out[16],
                         const VP8Matrix* const mtx) {
  return DoQuantizeBlock(in, out, 0, &mtx->sharpen_[0], mtx);
}

static int QuantizeBlockWHT(int16_t in[16], int16_t out[16],
                            const VP8Matrix* const mtx) {
  return DoQuantizeBlock(in, out, 0, &mtx->sharpen_[0], mtx);
}

#endif   // WEBP_USE_AVX2

//------------------------------------------------------------------------------
// Entry point

extern void VP8EncDspInitAVX2(void);

void VP8EncDspInitAVX2(void) {
#if defined(WEBP_USE_AVX2)
  VP8EncQuantizeBlock = QuantizeBlock;
  VP8EncQuantizeBlockWHT = QuantizeBlockWHT;
#endif   // WEBP_USE_AVX2
}


