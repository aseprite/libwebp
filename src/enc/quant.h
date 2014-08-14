#include "./vp8enci.h"
#include "./cost.h"

// If a coefficient was quantized to a value Q (using a neutral bias),
// we test all alternate possibilities between [Q-MIN_DELTA, Q+MAX_DELTA]
// We don't test negative values though.
#define DO_TRELLIS_I4  1
#define MIN_DELTA 0   // how much lower level to try
#define MAX_DELTA 1   // how much higher
#define NUM_NODES (MIN_DELTA + 1 + MAX_DELTA)
#define NODE(n, l) (nodes[(n)][(l) + MIN_DELTA])
#define SCORE_STATE(n, l) (score_states[n][(l) + MIN_DELTA])

#define FLATNESS_LIMIT_I4  3       // I4 mode
#define FLATNESS_PENALTY   140     // roughly ~1bit per block
#define MULT_8B(a, b) (((a) * (b) + 128) >> 8)
//------------------------------------------------------------------------------


// Performs trellis-optimized quantization.

// Trellis node
typedef struct {
  int8_t prev;            // best previous node
  int8_t sign;            // sign of coeff_i
  int16_t level;          // level
} Node;

// Score state
typedef struct {
  score_t score;          // partial RD score
  const uint16_t* costs;  // shortcut to cost tables
} ScoreState;

static const uint8_t kZigzag[16] = {
  0, 1, 4, 8, 5, 2, 3, 6, 9, 12, 13, 10, 7, 11, 14, 15
};

//------------------------------------------------------------------------------
// Distortion measurement
static const uint16_t kWeightY[16] = {
  38, 32, 20, 9, 32, 28, 17, 7, 20, 17, 10, 4, 9, 7, 4, 2
};

static const uint16_t kWeightTrellis[16] = {
#if USE_TDISTO == 0
  16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16
#else
  30, 27, 19, 11,
  27, 24, 17, 10,
  19, 17, 12,  8,
  11, 10,  8,  6
#endif
};

// Init/Copy the common fields in score.
inline static void InitScore(struct VP8ModeScore* const rd) {
  rd->D  = 0;
  rd->SD = 0;
  rd->R  = 0;
  rd->H  = 0;
  rd->nz = 0;
  rd->score = MAX_COST;
}

inline static void CopyScore(VP8ModeScore* const dst, const VP8ModeScore* const src) {
  dst->D  = src->D;
  dst->SD = src->SD;
  dst->R  = src->R;
  dst->H  = src->H;
  dst->nz = src->nz;      // note that nz is not accumulated, but just copied.
  dst->score = src->score;
}

inline static void AddScore(VP8ModeScore* const dst, const VP8ModeScore* const src) {
  dst->D  += src->D;
  dst->SD += src->SD;
  dst->R  += src->R;
  dst->H  += src->H;
  dst->nz |= src->nz;     // here, new nz bits are accumulated.
  dst->score += src->score;
}

//------------------------------------------------------------------------------
// Performs trellis-optimized quantization.

static WEBP_INLINE void SetRDScore(int lambda, VP8ModeScore* const rd) {
  // TODO: incorporate the "* 256" in the tables?
  rd->score = (rd->R + rd->H) * lambda + 256 * (rd->D + rd->SD);
}

static WEBP_INLINE score_t RDScoreTrellis(int lambda, score_t rate,
                                          score_t distortion) {
  return rate * lambda + 256 * distortion;
}

inline static void SwapPtr(uint8_t** a, uint8_t** b) {
  uint8_t* const tmp = *a;
  *a = *b;
  *b = tmp;
}

inline static void SwapOut(VP8EncIterator* const it) {
  SwapPtr(&it->yuv_out_, &it->yuv_out2_);
}

inline static score_t IsFlat(const int16_t* levels, int num_blocks, score_t thresh) {
  score_t score = 0;
  while (num_blocks-- > 0) {      // TODO(skal): refine positional scoring?
    int i;
    for (i = 1; i < 16; ++i) {    // omit DC, we're only interested in AC
      score += (levels[i] != 0);
      if (score > thresh) return 0;
    }
    levels += 16;
  }
  return 1;
}
//------------------------------------------------------------------------------

// return the cost array corresponding to the surrounding prediction modes.
inline static const uint16_t* GetCostModeI4(VP8EncIterator* const it,
                                     const uint8_t modes[16]) {
  const int preds_w = it->enc_->preds_w_;
  const int x = (it->i4_ & 3), y = it->i4_ >> 2;
  const int left = (x == 0) ? it->preds_[y * preds_w - 1] : modes[it->i4_ - 1];
  const int top = (y == 0) ? it->preds_[-preds_w + x] : modes[it->i4_ - 4];
  return VP8FixedCostsI4[top][left];
}

