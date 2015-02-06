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

#ifndef WEBP_DSP_COST_H_
#define WEBP_DSP_COST_H_

#include "../enc/cost.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*VP8SetResidualCoeffsFunc)(const int16_t* const coeffs,
                                         VP8Residual* const res);
extern VP8SetResidualCoeffsFunc VP8SetResidualCoeffs;

// Cost calculation function.
typedef int (*VP8GetResidualCostFunc)(int ctx0, const VP8Residual* const res);
extern VP8GetResidualCostFunc VP8GetResidualCost;

void VP8EncDspCostInit(void);  // must be called first

//------------------------------------------------------------------------------

#ifdef __cplusplus
}    // extern "C"
#endif

#endif  /* WEBP_DSP_COST_H_ */
