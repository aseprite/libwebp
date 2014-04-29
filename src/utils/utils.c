// Copyright 2012 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
// Misc. common utility functions
//
// Author: Skal (pascal.massimino@gmail.com)

#include <stdlib.h>
#include "./utils.h"

// If defined, will print extra info like total memory used, number of
// alloc/free etc. For debugging/tuning purpose only (it's slow).
// #define PRINT_MEM_INFO
#define PRINT_MEM_TRAFFIC   // print fine traffic details

//------------------------------------------------------------------------------
// Checked memory allocation

#if defined(PRINT_MEM_INFO)

#include <stdio.h>
#include <stdlib.h>  // for abort()

static uint64_t num_malloc_calls = 0;
static uint64_t num_calloc_calls = 0;
static uint64_t num_free_calls = 0;

typedef struct MemBlock MemBlock;
struct MemBlock {
  void* ptr_;
  size_t size_;
  MemBlock* next_;
};

static MemBlock* all_blocks = NULL;
static size_t total_mem = 0;
static size_t high_water_mark = 0;

static int exit_registered = 0;
static void PrintMemInfo() {
  fprintf(stderr, "\nMEMORY INFO:\n");
  fprintf(stderr, "num calls to: malloc = %lld\n", num_malloc_calls);
  fprintf(stderr, "              calloc = %lld\n", num_calloc_calls);
  fprintf(stderr, "              free   = %lld\n", num_free_calls);
  fprintf(stderr, "total_mem: %ld\n", total_mem);
  fprintf(stderr, "high-water mark: %ld\n", high_water_mark);
  while (all_blocks != NULL) {
    MemBlock* b = all_blocks;
    all_blocks = b->next_;
    free(b);
  }
}

static void Increment(uint64_t* const v) {
  if (!exit_registered) {
    atexit(PrintMemInfo);
    exit_registered = 1;
  }
  ++*v;
}

static void AddMem(void* ptr, size_t size) {
  if (ptr != NULL) {
    MemBlock* b = (MemBlock*)malloc(sizeof(*b));
    if (b == NULL) abort();
    b->next_ = all_blocks;
    all_blocks = b;
    b->ptr_ = ptr;
    b->size_ = size;
    total_mem += size;
#if defined(PRINT_MEM_TRAFFIC)
    fprintf(stderr, "Mem: %ld (+%ld)\n", total_mem, size);
#endif
    if (total_mem > high_water_mark) high_water_mark = total_mem;
  }
}

static void SubMem(void* ptr) {
  if (ptr != NULL) {
    MemBlock** b = &all_blocks;
    // Inefficient search, but that's just for debugging.
    while (*b != NULL && (*b)->ptr_ != ptr) b = &(*b)->next_;
    if (*b == NULL) {
      fprintf(stderr, "Invalid pointer free! (%p)\n", ptr);
      abort();
    }
    MemBlock* block = *b;
    *b = block->next_;
    total_mem -= block->size_;
#if defined(PRINT_MEM_TRAFFIC)
    fprintf(stderr, "Mem: %ld (-%ld)\n", total_mem, block->size_);
#endif
    free(block);
  }
}

#else
#define Increment(v) do {} while(0)
#define AddMem(p, s) do {} while(0)
#define SubMem(p)    do {} while(0)
#endif

// Returns 0 in case of overflow of nmemb * size.
static int CheckSizeArgumentsOverflow(uint64_t nmemb, size_t size) {
  const uint64_t total_size = nmemb * size;
  if (nmemb == 0) return 1;
  if ((uint64_t)size > WEBP_MAX_ALLOCABLE_MEMORY / nmemb) return 0;
  if (total_size != (size_t)total_size) return 0;
  return 1;
}

void* WebPSafeMalloc(uint64_t nmemb, size_t size) {
  void* ptr;
  Increment(&num_malloc_calls);
  if (!CheckSizeArgumentsOverflow(nmemb, size)) return NULL;
  assert(nmemb * size > 0);
  ptr = malloc((size_t)(nmemb * size));
  AddMem(ptr, nmemb * size);
  return ptr;
}

void* WebPSafeCalloc(uint64_t nmemb, size_t size) {
  void* ptr;
  Increment(&num_calloc_calls);
  if (!CheckSizeArgumentsOverflow(nmemb, size)) return NULL;
  assert(nmemb * size > 0);
  ptr = calloc((size_t)nmemb, size);
  AddMem(ptr, nmemb * size);
  return ptr;
}

void WebPSafeFree(void* const ptr) {
  Increment(&num_free_calls);
  SubMem(ptr);
  free(ptr);
}

//------------------------------------------------------------------------------
