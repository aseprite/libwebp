// Copyright 2012 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
//  Utility functions used by the example programs.
//

#ifndef WEBP_EXAMPLES_EXAMPLE_UTIL_H_
#define WEBP_EXAMPLES_EXAMPLE_UTIL_H_

#include "webp/types.h"

// Allocates storage for entire file 'file_name' and returns contents and size
// in 'data' and 'data_size'. Returns 1 on success, 0 otherwise. '*data' should
// be deleted using free().
int ExUtilReadFile(const char* const file_name,
                   const uint8_t** data, size_t* data_size);

#endif  // WEBP_EXAMPLES_EXAMPLE_UTIL_H_
