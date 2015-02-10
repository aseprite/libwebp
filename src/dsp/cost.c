// Copyright 2014 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
// Author: Skal (pascal.massimino@gmail.com)

#include "./dsp.h"
#include "../enc/cost.h"

//------------------------------------------------------------------------------
// Mode costs

static int GetResidualCost(int ctx0, const VP8Residual* const res) {
  int n = res->first;
  // should be prob[VP8EncBands[n]], but it's equivalent for n=0 or 1
  const int p0 = res->prob[n][ctx0][0];
  const uint16_t* t = res->cost[n][ctx0];
  // bit_cost(1, p0) is already incorporated in t[] tables, but only if ctx != 0
  // (as required by the syntax). For ctx0 == 0, we need to add it here or it'll
  // be missing during the loop.
  int cost = (ctx0 == 0) ? VP8BitCost(1, p0) : 0;

  if (res->last < 0) {
    return VP8BitCost(0, p0);
  }
  for (; n < res->last; ++n) {
    const int v = abs(res->coeffs[n]);
    const int b = VP8EncBands[n + 1];
    const int ctx = (v >= 2) ? 2 : v;
    cost += VP8LevelCost(t, v);
    t = res->cost[b][ctx];
  }
  // Last coefficient is always non-zero
  {
    const int v = abs(res->coeffs[n]);
    assert(v != 0);
    cost += VP8LevelCost(t, v);
    if (n < 15) {
      const int b = VP8EncBands[n + 1];
      const int ctx = (v == 1) ? 1 : 2;
      const int last_p0 = res->prob[b][ctx][0];
      cost += VP8BitCost(0, last_p0);
    }
  }
  return cost;
}

static void SetResidualCoeffs(const int16_t* const coeffs,
                              VP8Residual* const res) {
  int n;
  res->last = -1;
  assert(res->first == 0 || coeffs[0] == 0);
  for (n = 15; n >= 0; --n) {
    if (coeffs[n]) {
      res->last = n;
      break;
    }
  }
  res->coeffs = coeffs;
}

static int VariableLevelCost(int level, const uint8_t probas[NUM_PROBAS]) {
  int pattern = VP8LevelCodes[level - 1][0];
  int bits = VP8LevelCodes[level - 1][1];
  int cost = 0;
  int i;
  for (i = 2; pattern; ++i) {
    if (pattern & 1) {
      cost += VP8BitCost(bits & 1, probas[i]);
    }
    bits >>= 1;
    pattern >>= 1;
  }
  return cost;
}

//------------------------------------------------------------------------------
// Pre-calc level costs once for all

static void CalculateLevelCosts(VP8Proba* const proba) {
  int ctype, band, ctx;

  if (!proba->dirty_) return;  // nothing to do.

  for (ctype = 0; ctype < NUM_TYPES; ++ctype) {
    for (band = 0; band < NUM_BANDS; ++band) {
      for (ctx = 0; ctx < NUM_CTX; ++ctx) {
        const uint8_t* const p = proba->coeffs_[ctype][band][ctx];
        uint16_t* const table = proba->level_cost_[ctype][band][ctx];
        const int cost0 = (ctx > 0) ? VP8BitCost(1, p[0]) : 0;
        const int cost_base = VP8BitCost(1, p[1]) + cost0;
        int v;
        table[0] = VP8BitCost(0, p[1]) + cost0;
        for (v = 1; v <= MAX_VARIABLE_LEVEL; ++v) {
          table[v] = cost_base + VariableLevelCost(v, p);
        }
        // Starting at level 67 and up, the variable part of the cost is
        // actually constant.
      }
    }
  }
  proba->dirty_ = 0;
}
//------------------------------------------------------------------------------
// init function

VP8GetResidualCostFunc VP8GetResidualCost;
VP8SetResidualCoeffsFunc VP8SetResidualCoeffs;
VP8CalculateLevelCostsFunc VP8CalculateLevelCosts;

extern void VP8EncDspCostInitMIPS32(void);
extern void VP8EncDspCostInitMIPSdspR2(void);

#if defined(WEBP_USE_SSE2)
extern void VP8SetResidualCoeffsSSE2(const int16_t* const coeffs,
                                     VP8Residual* const res);
#endif  // WEBP_USE_SSE2

static volatile VP8CPUInfo cost_last_cpuinfo_used =
    (VP8CPUInfo)&cost_last_cpuinfo_used;

WEBP_TSAN_IGNORE_FUNCTION void VP8EncDspCostInit(void) {
  if (cost_last_cpuinfo_used == VP8GetCPUInfo) return;

  VP8GetResidualCost = GetResidualCost;
  VP8SetResidualCoeffs = SetResidualCoeffs;
  VP8CalculateLevelCosts = CalculateLevelCosts;

  // If defined, use CPUInfo() to overwrite some pointers with faster versions.
  if (VP8GetCPUInfo != NULL) {
#if defined(WEBP_USE_MIPS32)
    if (VP8GetCPUInfo(kMIPS32)) {
      VP8EncDspCostInitMIPS32();
    }
#endif
#if defined(WEBP_USE_MIPS_DSP_R2)
    if (VP8GetCPUInfo(kMIPSdspR2)) {
      VP8EncDspCostInitMIPSdspR2();
    }
#endif
#if defined(WEBP_USE_SSE2)
    if (VP8GetCPUInfo(kSSE2)) {
      VP8SetResidualCoeffs = VP8SetResidualCoeffsSSE2;
    }
#endif
  }

  cost_last_cpuinfo_used = VP8GetCPUInfo;
}

//------------------------------------------------------------------------------
