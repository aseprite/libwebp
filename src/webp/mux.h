// Copyright 2011 Google Inc.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
//  RIFF container manipulation for WEBP images.
//
// Author: Urvang (urvang@google.com)

#ifndef WEBP_WEBP_MUX_H_
#define WEBP_WEBP_MUX_H_

#include <string.h>

#include "webp/types.h"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

// Error codes
typedef enum {
  WEBP_MUX_OK                 =  1,
  WEBP_MUX_ERROR              =  0,
  WEBP_MUX_NOT_FOUND          = -1,
  WEBP_MUX_INVALID_ARGUMENT   = -2,
  WEBP_MUX_INVALID_PARAMETER  = -3,
  WEBP_MUX_BAD_DATA           = -4,
  WEBP_MUX_MEMORY_ERROR       = -5
} WebPMuxError;

// Flag values for different features used in vp8x chunk.
typedef enum {
  TILE_FLAG       = 0x00000001,
  ANIMATION_FLAG  = 0x00000002,
  ICCP_FLAG       = 0x00000004,
  META_FLAG       = 0x00000008,
} FeatureFlags;

typedef struct WebPMux WebPMux;   // main opaque object.

//-----------------------------------------------------------------------------
// Life of a Mux object

// Creates an empty mux object.
// copy_data - value 1 indicates data WILL copied in subsequent calls, and
//             value 0 indicates data will NOT be copied.
// Returns:
//   A pointer to the newly created empty mux object.
WebPMux* WebPMuxNew(int copy_data);

// Deletes the mux object.
void WebPMuxDelete(WebPMux* const mux);

//-----------------------------------------------------------------------------
// Writing

// Creates a mux object from raw data given in WebP RIFF format.
// WebPMuxNew() should be called before calling this function.
// Returns:
//   WEBP_MUX_INVALID_ARGUMENT - if mux is NULL.
//   WEBP_MUX_BAD_DATA - if data cannot be read as a mux object.
//   WEBP_MUX_MEMORY_ERROR - on memory allocation error.
//   WEBP_MUX_OK - on success.
WebPMuxError WebPMuxCreate(WebPMux* const mux, const uint8_t* data,
                           uint32_t size);

// Sets the XMP metadata in the mux object. Any existing metadata chunk(s) will
// be removed.
// Returns:
//   WEBP_MUX_MEMORY_ERROR - on memory allocation error.
//   WEBP_MUX_OK - on success.
WebPMuxError WebPMuxSetMetadata(WebPMux* const mux, const uint8_t* data,
                                uint32_t size);

// Sets the color profile in the mux object. Any existing color profile chunk(s)
// will be removed.
// Returns:
//   WEBP_MUX_MEMORY_ERROR - on memory allocation error
//   WEBP_MUX_OK - on success
WebPMuxError WebPMuxSetColorProfile(WebPMux* const mux, const uint8_t* data,
                                    uint32_t size);

// Sets the animation loop count in the mux object. Any existing loop count
// value(s) will be removed.
// Returns:
//   WEBP_MUX_MEMORY_ERROR - on memory allocation error.
//   WEBP_MUX_OK - on success.
WebPMuxError WebPMuxSetLoopCount(WebPMux* const mux, uint32_t loop_count);

// Adds an animation frame to the mux object.
// nth=0 has a special meaning - last position.
// duration is in milliseconds.
// Returns:
//   WEBP_MUX_NOT_FOUND - If we have less than (nth-1) frames before adding.
//   WEBP_MUX_MEMORY_ERROR - on memory allocation error.
//   WEBP_MUX_OK - on success.
WebPMuxError WebPMuxAddFrame(WebPMux* const mux, uint32_t nth,
                             const uint8_t* data, uint32_t size,
                             uint32_t x_offset, uint32_t y_offset,
                             uint32_t duration);

// Adds a tile to the mux object.
// nth=0 has a special meaning - last position.
// Returns:
//   WEBP_MUX_NOT_FOUND - If we have less than (nth-1) tiles before adding.
//   WEBP_MUX_MEMORY_ERROR - on memory allocation error.
//   WEBP_MUX_OK - on success.
WebPMuxError WebPMuxAddTile(WebPMux* const mux, uint32_t nth,
                            const uint8_t* data, uint32_t size,
                            uint32_t x_offset, uint32_t y_offset);

// Adds a chunk with given tag at nth position.
// nth=0 has a special meaning - last position.
// Returns:
//   WEBP_MUX_INVALID_ARGUMENT - if mux is NULL.
//   WEBP_MUX_INVALID_PARAMETER - if tag is invalid.
//   WEBP_MUX_MEMORY_ERROR - on memory allocation error.
//   WEBP_MUX_OK - on success.
WebPMuxError WebPMuxAddNamedData(WebPMux* const mux, uint32_t nth,
                                 const char* const tag, const uint8_t* data,
                                 uint32_t size);

// Deletes the XMP metadata in the mux object.
// Returns:
//   WEBP_MUX_NOT_FOUND - If mux does not contain metadata.
//   WEBP_MUX_OK - on success.
WebPMuxError WebPMuxDeleteMetadata(WebPMux* const mux);

// Deletes the color profile in the mux object.
// Returns:
//   WEBP_MUX_NOT_FOUND - If mux does not contain color profile.
//   WEBP_MUX_OK - on success.
WebPMuxError WebPMuxDeleteColorProfile(WebPMux* const mux);

// Deletes an animation frame from the mux object.
// nth=0 has a special meaning - last position.
// Returns:
//   WEBP_MUX_NOT_FOUND - If there are less than nth frames in the mux object
//                        before deletion.
//   WEBP_MUX_OK - on success.
WebPMuxError WebPMuxDeleteFrame(WebPMux* const mux, uint32_t nth);

// Deletes a tile from the mux object.
// nth=0 has a special meaning - last position
// Returns:
//   WEBP_MUX_NOT_FOUND - If there are less than nth tiles in the mux object
//                        before deletion.
//   WEBP_MUX_OK - on success.
WebPMuxError WebPMuxDeleteTile(WebPMux* const mux, uint32_t nth);

// Assembles all chunks in WebP RIFF format and and returns in output_data.
// This function also validates the mux object.
// This function does NOT free the mux object. The caller function should free
// the mux object after output_data has been consumed.
// Returns:
//   WEBP_MUX_BAD_DATA - if mux object is invalid.
//   WEBP_MUX_INVALID_ARGUMENT - if either mux, output_data or output_size is
//                               NULL.
//   WEBP_MUX_MEMORY_ERROR - on memory allocation error.
//   WEBP_MUX_OK - on success
WebPMuxError WebPMuxAssemble(WebPMux* const mux, uint8_t** output_data,
                             uint32_t* output_size);

//-----------------------------------------------------------------------------
// Reading

// Gets the feature flags from the mux object.
// Use enum 'FeatureFlags' to test for individual features.
// Returns:
//   WEBP_MUX_NOT_FOUND - if VP8X chunk is not present in mux object.
//   WEBP_MUX_BAD_DATA - if VP8X chunk in mux is invalid.
//   WEBP_MUX_OK - on success.
WebPMuxError WebPMuxGetFeatures(const WebPMux* mux, uint32_t* flags);

// Gets a reference to the XMP metadata in the mux object.
// The caller should not free the returned data.
// Returns:
//   WEBP_MUX_NOT_FOUND - if metadata is not present in mux object.
//   WEBP_MUX_OK - on success.
WebPMuxError WebPMuxGetMetadata(const WebPMux* mux, const uint8_t** data,
                                uint32_t* size);

// Gets a reference to the color profile in the mux object.
// The caller should not free the returned data.
// Returns:
//   WEBP_MUX_NOT_FOUND - if color profile is not present in mux object.
//   WEBP_MUX_OK - on success.
WebPMuxError WebPMuxGetColorProfile(const WebPMux* mux, const uint8_t** data,
                                    uint32_t* size);

// Gets the animation loop count from the mux object.
// Returns:
//   WEBP_MUX_NOT_FOUND - if loop chunk is not present in mux object.
//   WEBP_MUX_OK - on success.
WebPMuxError WebPMuxGetLoopCount(const WebPMux* mux, uint32_t* loop_count);

// Gets a reference to the nth animation frame from the mux object.
// nth=0 has a special meaning - last position.
// duration is in milliseconds.
// Returns:
//   WEBP_MUX_NOT_FOUND - if there are less than nth frames in the mux object.
//   WEBP_MUX_BAD_DATA - if nth frame chunk in mux is invalid.
//   WEBP_MUX_OK - on success.
WebPMuxError WebPMuxGetFrame(const WebPMux* mux, uint32_t nth,
                             const uint8_t** data, uint32_t* size,
                             uint32_t* x_offset, uint32_t* y_offset,
                             uint32_t* duration);

// Gets a reference to the nth tile from the mux object.
// The caller should not free the returned data.
// nth=0 has a special meaning - last position.
// Returns:
//   WEBP_MUX_NOT_FOUND - if there are less than nth tiles in the mux object.
//   WEBP_MUX_BAD_DATA - if nth tile chunk in mux is invalid.
//   WEBP_MUX_OK - on success.
WebPMuxError WebPMuxGetTile(const WebPMux* mux, uint32_t nth,
                            const uint8_t** data, uint32_t* size,
                            uint32_t* x_offset, uint32_t* y_offset);

// Gets number of chunks having tag value tag in the mux object.
// Returns:
//   WEBP_MUX_INVALID_ARGUMENT - if mux is NULL.
//   WEBP_MUX_OK - on success.
WebPMuxError WebPMuxNumNamedElements(const WebPMux* mux, const char* tag,
                                     int* num_elements);

// Gets a reference to the nth chunk having tag value tag in the mux object.
// The caller should not free the returned data.
// nth=0 has a special meaning - last position.
// Returns:
//   WEBP_MUX_NOT_FOUND - If there are less than nth named elements in the mux
//                        object.
//   WEBP_MUX_OK - on success.
WebPMuxError WebPMuxGetNamedData(const WebPMux* mux, const char* tag,
                                 uint32_t nth, const uint8_t** data,
                                 uint32_t* size);

//-----------------------------------------------------------------------------

#if defined(__cplusplus) || defined(c_plusplus)
}    // extern "C"
#endif

#endif  /* WEBP_WEBP_MUX_H_ */
