#include <assert.h>
#include <math.h>
#include <stdlib.h>  // for abs()

#include "quant.h"
#include "../dsp/dsp.h"

#if defined(WEBP_USE_AVX2)

//------------------------------------------------------------------------------
// Performs trellis-optimized quantization for the current PickBestIntra4 
// which perform 2 iterations at a time.

static int TrellisQuantizeBlock4x2(const VP8Encoder* const enc,
                                int16_t* in, int16_t* out,
                                int ctx0, int coeff_type,
                                const VP8Matrix* const mtx,
                                int lambda) {
  const ProbaArray* const probas = enc->proba_.coeffs_[coeff_type];
  const CostArray* const costs = enc->proba_.level_cost_[coeff_type];
  const int first = (coeff_type == 0) ? 1 : 0;
  Node nodes[16][NUM_NODES];
  ScoreState score_states[2][NUM_NODES];
  ScoreState* ss_cur = &SCORE_STATE(0, MIN_DELTA);
  ScoreState* ss_prev = &SCORE_STATE(1, MIN_DELTA);
  int best_path[3] = {-1, -1, -1};   // store best-last/best-level/best-previous
  score_t best_score;
  int n, m, p, last;

  {
    score_t cost;
    int err;
    const int thresh = mtx->q_[1] * mtx->q_[1] / 4;
    const int last_proba = probas[VP8EncBands[first]][ctx0][0];

    // compute the position of the last interesting coefficient
    last = first - 1;
    for (n = 15 ; n >= first; --n) {
      const int j = kZigzag[n];
      if (j > 7)
          err = in[(j+8)] * in[(j+8)];
      else
          err = in[j] * in[j];
      if (err > thresh) {
        last = n;
        break;
      }
    
    }
    // we don't need to go inspect up to n = 16 coeffs. We can just go up
    // to last + 1 (inclusive) without losing much.
    if (last < 15) ++last;

    // compute 'skip' score. This is the max score one can do.
    cost = VP8BitCost(0, last_proba);
    best_score = RDScoreTrellis(lambda, cost, 0);

    // initialize source node.
    for (m = -MIN_DELTA; m <= MAX_DELTA; ++m) {
      const score_t rate = (ctx0 == 0) ? VP8BitCost(1, last_proba) : 0;
      ss_cur[m].score = RDScoreTrellis(lambda, rate, 0);
      ss_cur[m].costs = costs[VP8EncBands[first]][ctx0];
    }
  }

  // traverse trellis.
  for (n = first; n <= last; ++n) {
    int16_t inref;
    int sign, level0;
    uint32_t coeff0;
    const int j = kZigzag[n];
    const uint32_t Q  = mtx->q_[j];
    const uint32_t iQ = mtx->iq_[j];
    const uint32_t B = BIAS(0x00);     // neutral bias
    // because the traverse trellis is done through PickBestIntra4 that process
    // 2 iterations at a time the data in the array "in" organized like that:
    // iter1,iter2,iter1,iter2.
    // when we need to access the 8th element from the first iteration
    // we need to shift by 8 element (all second iteration)
    if (j > 7)
          inref = in[j+8];
    else
          inref = in[j];
    sign = (inref < 0);
    coeff0 = (sign ? -inref : inref) + mtx->sharpen_[j];
    level0 = QUANTDIV(coeff0, iQ, B);
    if (level0 > MAX_LEVEL) level0 = MAX_LEVEL;

    {   // Swap current and previous score states
      ScoreState* const tmp = ss_cur;
      ss_cur = ss_prev;
      ss_prev = tmp;
    }

    // test all alternate level values around level0.
    for (m = -MIN_DELTA; m <= MAX_DELTA; ++m) {
      Node* const cur = &NODE(n, m);
      int level = level0 + m;
      const int ctx = (level > 2) ? 2 : level;
      const int band = VP8EncBands[n + 1];
      score_t base_score, last_pos_score;
      score_t best_cur_score = MAX_COST;
      int best_prev = 0;   // default, in case

      ss_cur[m].score = MAX_COST;
      ss_cur[m].costs = costs[band][ctx];
      if (level > MAX_LEVEL || level < 0) {   // node is dead?
        continue;
      }

      // Compute extra rate cost if last coeff's position is < 15
      {
        const score_t last_pos_cost =
            (n < 15) ? VP8BitCost(0, probas[band][ctx][0]) : 0;
        last_pos_score = RDScoreTrellis(lambda, last_pos_cost, 0);
      }

      {
        // Compute delta_error = how much coding this level will
        // subtract to max_error as distortion.
        // Here, distortion = sum of (|coeff_i| - level_i * Q_i)^2
        const int new_error = coeff0 - level * Q;
        const int delta_error =
            kWeightTrellis[j] * (new_error * new_error - coeff0 * coeff0);
        base_score = RDScoreTrellis(lambda, 0, delta_error);
      }

      // Inspect all possible non-dead predecessors. Retain only the best one.
      for (p = -MIN_DELTA; p <= MAX_DELTA; ++p) {
        // Dead nodes (with ss_prev[p].score >= MAX_COST) are automatically
        // eliminated since their score can't be better than the current best.
        const score_t cost = VP8LevelCost(ss_prev[p].costs, level);
        // Examine node assuming it's a non-terminal one.
        const score_t score =
            base_score + ss_prev[p].score + RDScoreTrellis(lambda, cost, 0);
        if (score < best_cur_score) {
          best_cur_score = score;
          best_prev = p;
        }
      }
      // Store best finding in current node.
      cur->sign = sign;
      cur->level = level;
      cur->prev = best_prev;
      ss_cur[m].score = best_cur_score;

      // Now, record best terminal node (and thus best entry in the graph).
      if (level != 0) {
        const score_t score = best_cur_score + last_pos_score;
        if (score < best_score) {
          best_score = score;
          best_path[0] = n;                     // best eob position
          best_path[1] = m;                     // best node index
          best_path[2] = best_prev;             // best predecessor
        }
      }
    }
  }

  // Fresh start
  memset(in ,0 ,8 * sizeof(*in));
  memset(in + 16, 0, 8 * sizeof(*in));
  
  memset(out, 0, 16 * sizeof(*out));
  if (best_path[0] == -1) {
    return 0;   // skip!
  }

  {
    // Unwind the best path.
    // Note: best-prev on terminal node is not necessarily equal to the
    // best_prev for non-terminal. So we patch best_path[2] in.
    int16_t inref;
    int nz = 0;
    int best_node = best_path[1];
    n = best_path[0];
    NODE(n, best_node).prev = best_path[2];   // force best-prev for terminal

    for (; n >= first; --n) {
      const Node* const node = &NODE(n, best_node);
      const int j = kZigzag[n];
      out[n] = node->sign ? -node->level : node->level;
      nz |= node->level;
      inref = out[n] * mtx->q_[j];
      if (j > 7)
         in[j+8] = inref;
      else
         in[j] = inref;
      best_node = node->prev;
    }
    return (nz != 0);
  }
}


static void ReconstructIntra4x2(VP8EncIterator* const it,
                             int16_t levels[32],
                             const uint8_t* const src,
                             uint8_t* yuv_out,
                             int mode,
                             VP8ModeScore* rd_tmp1,
                             VP8ModeScore* rd_tmp2) {
  const VP8Encoder* const enc = it->enc_;
  const uint8_t* const ref = it->yuv_p_ + VP8I4ModeOffsets[mode];
  const VP8SegmentInfo* const dqm = &enc->dqm_[it->mb_->segment_];
  WEBP_ALIGNED_DECL(32, int16_t, tmp[32]);
  
  VP8FTransform4x2(src, ref, tmp);
  if (DO_TRELLIS_I4 && it->do_trellis_) {
    const int x = it->i4_ & 3, y = it->i4_ >> 2;
    const int ctx = it->top_nz_[x] + it->left_nz_[y];
    rd_tmp1->nz = (TrellisQuantizeBlock4x2(enc, tmp, levels, ctx, 3, &dqm->y1_,
                              dqm->lambda_trellis_i4_))<< it->i4_;
    
    rd_tmp2->nz = (TrellisQuantizeBlock4x2(enc, tmp+8, levels+16, ctx, 3, &dqm->y1_,
                              dqm->lambda_trellis_i4_))<< it->i4_;
  } else {
    VP8EncQuantizeBlock4x2(tmp, levels, &dqm->y1_, &rd_tmp1->nz, &rd_tmp2->nz);
    rd_tmp1->nz = rd_tmp1->nz << it->i4_;
    rd_tmp2->nz = rd_tmp2->nz << it->i4_;
  }
  VP8ITransform4x2(ref, tmp, yuv_out, 0);
}

//------------------------------------------------------------------------------

static int PickBestIntra4(VP8EncIterator* const it, VP8ModeScore* const rd) {
  const VP8Encoder* const enc = it->enc_;
  const VP8SegmentInfo* const dqm = &enc->dqm_[it->mb_->segment_];
  const int lambda = dqm->lambda_i4_;
  const int tlambda = dqm->tlambda_;
  const uint8_t* const src0 = it->yuv_in_ + Y_OFF;
  uint8_t* const best_blocks = it->yuv_out2_ + Y_OFF;
  int total_header_bits = 0;
  
  VP8ModeScore rd_best;
  if (enc->max_i4_header_bits_ == 0) {
    return 0;
  }
  InitScore(&rd_best);
  rd_best.H = 211;  // '211' is the value of VP8BitCost(0, 145)
  SetRDScore(dqm->lambda_mode_, &rd_best);
  VP8IteratorStartI4(it);
  do {
    WEBP_ALIGNED_DECL(32, uint8_t, tmp_reg[128]);
    WEBP_ALIGNED_DECL(32, uint8_t, bestblock_reg[128]);
    uint8_t* tmp_reg_ptr = tmp_reg;
    uint8_t* bestblock_reg_ptr = bestblock_reg;
     int best_lane = 0;
    const int kNumBlocks = 1;
    VP8ModeScore rd_i4;
    int mode;
    int best_mode = -1;
    const uint8_t* const src = src0 + VP8Scan[it->i4_];
    const uint16_t* const mode_costs = GetCostModeI4(it, rd->modes_i4);

    InitScore(&rd_i4);
    VP8MakeIntra4Preds(it);
    // perform 2 iterations at a time. 
    for (mode = 0; mode < NUM_BMODES; mode+=2) {
      VP8ModeScore rd_tmp1, rd_tmp2;
      int16_t tmp_levels[32];
      int chosen_mode;
      int chosen_lane;

      // Reconstruct
      ReconstructIntra4x2(it, tmp_levels, src, tmp_reg_ptr, mode, &rd_tmp1, &rd_tmp2);
      // Compute RD-score
      VP8SSE4x4x2(src, tmp_reg_ptr,&rd_tmp1.D,&rd_tmp2.D);
     
      if (tlambda)
      {
             VP8TDisto4x4x2(src, tmp_reg_ptr, kWeightY, &rd_tmp1.SD, &rd_tmp2.SD);
             rd_tmp1.SD = MULT_8B(tlambda, rd_tmp1.SD);
             rd_tmp2.SD = MULT_8B(tlambda, rd_tmp2.SD);
      }
      else
      {
            rd_tmp1.SD = 0;
            rd_tmp2.SD = 0;
      }
      
      rd_tmp1.H = mode_costs[mode];
      rd_tmp2.H = mode_costs[mode + 1];
      
      rd_tmp1.R = VP8GetCostLuma4(it, tmp_levels);
      rd_tmp2.R = VP8GetCostLuma4(it, tmp_levels + 16);
      
      if (mode > 0 && IsFlat(tmp_levels, kNumBlocks, FLATNESS_LIMIT_I4)) {
        rd_tmp1.R += FLATNESS_PENALTY * kNumBlocks;
      }
      if (IsFlat(tmp_levels + 16, kNumBlocks, FLATNESS_LIMIT_I4))
      {
        rd_tmp2.R += FLATNESS_PENALTY * kNumBlocks;
      }
      
      SetRDScore(lambda, &rd_tmp1);
      SetRDScore(lambda, &rd_tmp2);
      
      // check best score between first and second iteration.
      if (rd_tmp2.score < rd_tmp1.score)
      {
          CopyScore(&rd_tmp1, &rd_tmp2);
          chosen_lane = 16;
          chosen_mode = mode + 1;
      } 
      else
      {
         chosen_lane = 0;
         chosen_mode = mode;
      }
      if (best_mode < 0 || rd_tmp1.score < rd_i4.score) {
        CopyScore(&rd_i4, &rd_tmp1);
        best_mode = chosen_mode;
        best_lane = chosen_lane;
        SwapPtr(&tmp_reg_ptr, &bestblock_reg_ptr);
        memcpy(rd_best.y_ac_levels[it->i4_], tmp_levels + best_lane, 256);
      }
    }
    SetRDScore(dqm->lambda_mode_, &rd_i4);
    AddScore(&rd_best, &rd_i4);
    
    if (rd_best.score >= rd->score) {
      return 0;
    }
    total_header_bits += (int)rd_i4.H;   // <- equal to mode_costs[best_mode];
    if (total_header_bits > enc->max_i4_header_bits_) {
      return 0;
    }
   
    // Copy selected samples if not in the right place already.
    VP8Copy4x4x2((bestblock_reg_ptr+best_lane), best_blocks + VP8Scan[it->i4_]);
    rd->modes_i4[it->i4_] = best_mode;
    it->top_nz_[it->i4_ & 3] = it->left_nz_[it->i4_ >> 2] = (rd_i4.nz ? 1 : 0);
  } while (VP8IteratorRotateI4(it, best_blocks));

  // finalize state
  CopyScore(rd, &rd_best);
  VP8SetIntra4Mode(it, rd->modes_i4);
  SwapOut(it);
  memcpy(rd->y_ac_levels, rd_best.y_ac_levels, sizeof(rd->y_ac_levels));
  return 1;   // select intra4x4 over intra16x16
}


#endif  // WEBP_USE_AVX2

extern void VP8EncQuantInitAVX2(void);

void VP8EncQuantInitAVX2(void) {
#if defined(WEBP_USE_AVX2)
 VP8EncPickBestIntra4 = PickBestIntra4;
#endif
}


