// Copyright 2010 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
//  Common types
//
// Author: Skal (pascal.massimino@gmail.com)

#ifndef WEBP_WEBP_TYPES_H_
#define WEBP_WEBP_TYPES_H_

#include <stddef.h>  // for size_t

#ifndef _MSC_VER
#include <inttypes.h>
#ifdef __STRICT_ANSI__
#define WEBP_INLINE
#else  /* __STRICT_ANSI__ */
#define WEBP_INLINE inline
#endif
#else
typedef signed   char int8_t;
typedef unsigned char uint8_t;
typedef signed   short int16_t;
typedef unsigned short uint16_t;
typedef signed   int int32_t;
typedef unsigned int uint32_t;
typedef unsigned long long int uint64_t;
typedef long long int int64_t;
#define WEBP_INLINE __forceinline
#endif  /* _MSC_VER */

#ifndef WEBP_EXTERN
// This explicitly marks library functions and allows for changing the
// signature for e.g., Windows DLL builds.
#define WEBP_EXTERN(type) extern type
#endif  /* WEBP_EXTERN */

#include <stdlib.h>
#include <string.h>  // so that 'extern void memcpy(...) is defined first
#define memcpy(dst, src, size) do {             \
  uint8_t* const _dst = (uint8_t*)(dst);        \
  const uint8_t* const _src = (const uint8_t*)(src);  \
  const size_t _size = (size_t)(size);          \
  if ((_dst >= _src && _dst < _src + _size) ||  \
      (_dst < _src && _dst + _size > _src)) {  \
    abort();                \
  } else {                  \
    size_t i;               \
    for (i = 0; i < _size; ++i) _dst[i] = _src[i]; \
  }                         \
} while (0)

#endif  /* WEBP_WEBP_TYPES_H_ */
