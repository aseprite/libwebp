// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)
//
// Contains a lookup table for computing distance value prefix codes.

#include "backward_distance.h"

// This lut returns the distance value prefix code.
//
// Indexing to this lut is as follows:
// index = 2 * (hb + 1) + shb,
// where hb is the index of the highest bit set (-1 for 0),
// and shb is 0 if the second highest bit of the distance is not set
// and 1 when it is set.
//
// This is addressed by at max 65: 2 * 31 (the position of the highest bit)
// + 1 (the value of the second highest bit) + 2 for the
// distance_code_lut_offset
const BackwardDistance::DistanceCode
    BackwardDistance::distance_code_lut_[66] = {
  {0, 0},    {0, 0},  // -1 for Floor
  {0, 1},    {0, 1},
  {0, 2},    {0, 3},
  {1, 4},    {1, 5},
  {2, 6},    {2, 7},
  {3, 8},    {3, 9},
  {4, 10},   {4, 11},
  {5, 12},   {5, 13},
  {6, 14},   {6, 15},
  {7, 16},   {7, 17},
  {8, 18},   {8, 19},
  {9, 20},   {9, 21},
  {10, 22},  {10, 23},
  {11, 24},  {11, 25},
  {12, 26},  {12, 27},
  {13, 28},  {13, 29},
  {14, 30},  {14, 31},
  {15, 32},  {15, 33},
  {16, 34},  {16, 35},
  {17, 36},  {17, 37},
  {18, 38},  {18, 39},
  {19, 40},  {19, 41},
  {20, 42},  {20, 43},
  {21, 44},  {21, 45},
  {22, 46},  {22, 47},
  {23, 48},  {23, 49},
  {24, 50},  {24, 51},
  {25, 52},  {25, 53},
  {26, 54},  {26, 55},
  {27, 56},  {27, 57},
  {28, 58},  {28, 59},
  {29, 60},  {29, 61},
  {30, 62},  {30, 63},
};

const BackwardDistance::DistanceCode *
    BackwardDistance::distance_code_lut_offset_ = &distance_code_lut_[2];
