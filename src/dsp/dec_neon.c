// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// ARM NEON version of dsp functions and loop filtering.
//
// Author: somnath@google.com (Somnath Banerjee)

#if defined(__GNUC__) && defined(__ARM_NEON__)

#include "../dec/vp8i.h"
#include <arm_neon.h>

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

//-----------------------------------------------------------------------------
// Simple In-loop filtering (Paragraph 15.2)

static WEBP_INLINE uint8x16_t needs_filter_neon(uint8x16_t p1, uint8x16_t p0,
                                                uint8x16_t q0, uint8x16_t q1,
                                                int thresh) {
  uint8x16_t p0_q0,
             p1_q1,
             sum,
             dup_thresh,
             mask;

  dup_thresh = vdupq_n_u8(thresh);

  p0_q0 = vabdq_u8(p0, q0);
  p1_q1 = vabdq_u8(p1, q1);

  p0_q0 = vqaddq_u8(p0_q0, p0_q0); /* abs(p0-q0) * 2 */
  p1_q1 = vshrq_n_u8(p1_q1, 1); /* abs(p1-q1) / 2 */

  sum = vqaddq_u8(p0_q0, p1_q1);

  mask = vcgeq_u8(dup_thresh, sum);

  return mask;
}

static WEBP_INLINE int8x16_t get_base_delta(int8x16_t ps1, int8x16_t ps0,
                                            int8x16_t qs0, int8x16_t qs1) {
  int8x16_t qs0_ps0,
            ps1_qs1,
            sum;

  qs0_ps0 = vqsubq_s8(qs0, ps0);
  ps1_qs1 = vqsubq_s8(ps1, qs1);

  sum = vqaddq_s8(ps1_qs1, qs0_ps0);
  sum = vqaddq_s8(sum, qs0_ps0);
  sum = vqaddq_s8(sum, qs0_ps0);

  return sum;
}

static WEBP_INLINE int8x16_t do_simple_filter_ps0(int8x16_t ps0,
                                                  int8x16_t filter) {
  int8x16_t const3 = vdupq_n_s8(3);

  filter = vqaddq_s8(filter, const3);
  filter = vshrq_n_s8(filter, 3);

  ps0 = vqaddq_s8(ps0, filter);

  return ps0;
}

static WEBP_INLINE int8x16_t do_simple_filter_qs0(int8x16_t qs0,
                                                  int8x16_t filter) {
  /* "rounding" = 1 << (shift-1)
   * 1 << 2 == 4
   */
  filter = vrshrq_n_s8(filter, 3); //(+4) >> 3

  qs0 = vqsubq_s8(qs0, filter);

  return qs0;
}

static void SimpleVFilter16NEON(uint8_t* p, int stride, int thresh) {
  uint8x16_t p1, p0, q0, q1,
             mask, sign_bit;
  int8x16_t ps1, ps0, qs0, qs1,
            ps0_qs0, ps1_qs1,
            filter;

  p1 = vld1q_u8(p - 2 * stride);
  p0 = vld1q_u8(p - stride);
  q0 = vld1q_u8(p);
  q1 = vld1q_u8(p + stride);

  mask = needs_filter_neon(p1, p0, q0, q1, thresh);

  sign_bit = vdupq_n_u8(0x80);
  ps1 = vreinterpretq_s8_u8(veorq_u8(p1, sign_bit));
  ps0 = vreinterpretq_s8_u8(veorq_u8(p0, sign_bit));
  qs0 = vreinterpretq_s8_u8(veorq_u8(q0, sign_bit));
  qs1 = vreinterpretq_s8_u8(veorq_u8(q1, sign_bit));

  filter = get_base_delta(ps1, ps0, qs0, qs1);

  filter = vandq_s8(filter, (int8x16_t)mask);

  qs0 = do_simple_filter_qs0(qs0, filter);
  q0 = veorq_u8((uint8x16_t)qs0, sign_bit);

  ps0 = do_simple_filter_ps0(ps0, filter);
  p0 = veorq_u8(vreinterpretq_u8_s8(ps0), sign_bit);

  vst1q_u8(p - stride, p0);
  vst1q_u8(p, q0);
}

static void SimpleHFilter16NEON(uint8_t* p, int stride, int thresh) {
  uint8x8x4_t top, bottom;
  uint8x16_t p1, p0, q1, q0,
             mask, sign_bit;
  int8x16_t ps1, ps0, qs0, qs1,
            ps0_qs0, ps1_qs1,
            filter;
  uint8x8x2_t op0_oq0_top, op0_oq0_bottom;

  p -= 2;

  top = vld4_lane_u8(p+0*stride, top, 0);
  top = vld4_lane_u8(p+1*stride, top, 1);
  top = vld4_lane_u8(p+2*stride, top, 2);
  top = vld4_lane_u8(p+3*stride, top, 3);
  top = vld4_lane_u8(p+4*stride, top, 4);
  top = vld4_lane_u8(p+5*stride, top, 5);
  top = vld4_lane_u8(p+6*stride, top, 6);
  top = vld4_lane_u8(p+7*stride, top, 7);
  bottom = vld4_lane_u8(p+8*stride, bottom, 0);
  bottom = vld4_lane_u8(p+9*stride, bottom, 1);
  bottom = vld4_lane_u8(p+10*stride, bottom, 2);
  bottom = vld4_lane_u8(p+11*stride, bottom, 3);
  bottom = vld4_lane_u8(p+12*stride, bottom, 4);
  bottom = vld4_lane_u8(p+13*stride, bottom, 5);
  bottom = vld4_lane_u8(p+14*stride, bottom, 6);
  bottom = vld4_lane_u8(p+15*stride, bottom, 7);

  p1 = vcombine_u8(top.val[0], bottom.val[0]);
  p0 = vcombine_u8(top.val[1], bottom.val[1]);
  q0 = vcombine_u8(top.val[2], bottom.val[2]);
  q1 = vcombine_u8(top.val[3], bottom.val[3]);

  mask = needs_filter_neon(p1, p0, q0, q1, thresh);

  sign_bit = vdupq_n_u8(0x80);
  ps1 = vreinterpretq_s8_u8(veorq_u8(p1, sign_bit));
  ps0 = vreinterpretq_s8_u8(veorq_u8(p0, sign_bit));
  qs0 = vreinterpretq_s8_u8(veorq_u8(q0, sign_bit));
  qs1 = vreinterpretq_s8_u8(veorq_u8(q1, sign_bit));

  filter = get_base_delta(ps1, ps0, qs0, qs1);

  filter = vandq_s8(filter, (int8x16_t)mask);

  qs0 = do_simple_filter_qs0(qs0, filter);
  q0 = veorq_u8((uint8x16_t)qs0, sign_bit);

  ps0 = do_simple_filter_ps0(ps0, filter);
  p0 = veorq_u8(vreinterpretq_u8_s8(ps0), sign_bit);

  op0_oq0_top.val[0] = vget_low_u8(p0);
  op0_oq0_top.val[1] = vget_low_u8(q0);
  op0_oq0_bottom.val[0] = vget_high_u8(p0);
  op0_oq0_bottom.val[1] = vget_high_u8(q0);

  p += 1;

  vst2_lane_u8(p+0*stride, op0_oq0_top, 0);
  vst2_lane_u8(p+1*stride, op0_oq0_top, 1);
  vst2_lane_u8(p+2*stride, op0_oq0_top, 2);
  vst2_lane_u8(p+3*stride, op0_oq0_top, 3);
  vst2_lane_u8(p+4*stride, op0_oq0_top, 4);
  vst2_lane_u8(p+5*stride, op0_oq0_top, 5);
  vst2_lane_u8(p+6*stride, op0_oq0_top, 6);
  vst2_lane_u8(p+7*stride, op0_oq0_top, 7);
  vst2_lane_u8(p+8*stride, op0_oq0_bottom, 0);
  vst2_lane_u8(p+9*stride, op0_oq0_bottom, 1);
  vst2_lane_u8(p+10*stride, op0_oq0_bottom, 2);
  vst2_lane_u8(p+11*stride, op0_oq0_bottom, 3);
  vst2_lane_u8(p+12*stride, op0_oq0_bottom, 4);
  vst2_lane_u8(p+13*stride, op0_oq0_bottom, 5);
  vst2_lane_u8(p+14*stride, op0_oq0_bottom, 6);
  vst2_lane_u8(p+15*stride, op0_oq0_bottom, 7);
}

static void SimpleVFilter16iNEON(uint8_t* p, int stride, int thresh) {
  int k;
  for (k = 3; k > 0; --k) {
    p += 4 * stride;
    SimpleVFilter16NEON(p, stride, thresh);
  }
}

static void SimpleHFilter16iNEON(uint8_t* p, int stride, int thresh) {
  int k;
  for (k = 3; k > 0; --k) {
    p += 4;
    SimpleHFilter16NEON(p, stride, thresh);
  }
}

extern void VP8DspInitNEON(void);

void VP8DspInitNEON(void) {
  VP8SimpleVFilter16 = SimpleVFilter16NEON;
  VP8SimpleHFilter16 = SimpleHFilter16NEON;
  VP8SimpleVFilter16i = SimpleVFilter16iNEON;
  VP8SimpleHFilter16i = SimpleHFilter16iNEON;
}

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif   // __GNUC__ && __ARM_NEON__
