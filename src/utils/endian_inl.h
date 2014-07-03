// Copyright 2014 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
// Endian related functions.

#ifndef WEBP_UTILS_ENDIAN_INL_H_
#define WEBP_UTILS_ENDIAN_INL_H_

#include "../webp/types.h"

// some endian fix (e.g.: mips-gcc doesn't define __BIG_ENDIAN__)
#if !defined(__BIG_ENDIAN__) && defined(__BYTE_ORDER__) && \
    (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
#define __BIG_ENDIAN__
#endif

//  endian-specific htoleXX() definition
// TODO(skal): move this to config.h
#if defined(_WIN32)
#if !defined(_M_PPC)
#define htole32(x) (x)
#define htole16(x) (x)
#else     // PPC is BIG_ENDIAN
#include <stdlib.h>
#define htole32(x) (_byteswap_ulong((unsigned long)(x)))
#define htole16(x) (_byteswap_ushort((unsigned short)(x)))
#endif    // _M_PPC
#elif defined(__OpenBSD__) || defined(__NetBSD__) || defined(__FreeBSD__) || \
      defined(__DragonFly__)
#include <sys/endian.h>
#elif defined(__APPLE__)
#include <libkern/OSByteOrder.h>
#define htole32 OSSwapHostToLittleInt32
#define htole16 OSSwapHostToLittleInt16
#elif defined(__native_client__) && !defined(__GLIBC__)
// NaCl without glibc is assumed to be little-endian
#define htole32(x) (x)
#define htole16(x) (x)
#elif defined(__QNX__)
#include <net/netbyte.h>
#else     // pretty much all linux and/or glibc
#include <endian.h>
#endif

// gcc 4.3 has builtin functions for swap32/swap64
// TODO(jzern): this should have a corrsponding autoconf check.
#if defined(__GNUC__) && \
           (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 3))
#define HAVE_BUILTIN_BSWAP
#endif

static WEBP_INLINE uint32_t BSwap32(uint32_t x) {
#if defined(HAVE_BUILTIN_BSWAP)
  return __builtin_bswap32(x);
#elif defined(__i386__) || defined(__x86_64__)
  uint32_t swapped_bytes;
  __asm__ volatile("bswap %0" : "=r"(swapped_bytes) : "0"(x));
  return swapped_bytes;
#elif defined(_MSC_VER)
  return (uint32_t)_byteswap_ulong(x);
#else
  return (x >> 24) | ((x >> 8) & 0xff00) | ((x << 8) & 0xff0000) | (x << 24);
#endif  // HAVE_BUILTIN_BSWAP
}

#endif  // WEBP_UTILS_ENDIAN_INL_H_
