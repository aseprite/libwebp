// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// Author: jyrki@google.com (Jyrki Alakuijala)

#ifndef WEBP_INTEGRAL_TYPES_H_
#define WEBP_INTEGRAL_TYPES_H_

#include <stdlib.h>
#include <stdio.h>

#define VERIFY(a) \
  if (!(a)) { \
    fprintf(stderr, "Failed at %s:%d\n", __FILE__, __LINE__); \
    perror("Unrecoverable error"); \
    abort(); \
  };

#endif  // WEBP_INTEGRAL_TYPES_H_
