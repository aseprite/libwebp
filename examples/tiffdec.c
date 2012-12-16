// Copyright 2012 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
// TIFF decode.

#include "./tiffdec.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>

#ifdef WEBP_HAVE_TIFF
#include <tiffio.h>

#include "webp/encode.h"
#include "./metadata.h"

static const struct {
  ttag_t tag;
  size_t storage_offset;
} kTIFFMetadataMap[] = {
  { TIFFTAG_ICCPROFILE, METADATA_OFFSET(iccp) },
  { TIFFTAG_XMLPACKET,  METADATA_OFFSET(xmp) },
  { 0, 0 },
};

static int ExtractMetadataFromTIFF(TIFF* const tif, Metadata* const metadata) {
  int i;
  uint32 exif_ifd_offset;

  for (i = 0; kTIFFMetadataMap[i].tag != 0; ++i) {
    MetadataPayload* const payload =
        (MetadataPayload*)((uint8_t*)metadata +
                           kTIFFMetadataMap[i].storage_offset);
    void* tag_data;
    uint32 tag_data_len;

    if (TIFFGetField(tif, kTIFFMetadataMap[i].tag, &tag_data_len, &tag_data) &&
        !MetadataCopy((const char*)tag_data, tag_data_len, payload)) {
      return 0;
    }
  }

  // TODO(jzern): To extract the raw EXIF directory some parsing of it would be
  // necessary to determine the overall size. In addition, value offsets in
  // individual directory entries may need to be updated as, depending on the
  // type, they are file based.
  // Exif 2.2 Section 4.6.2 Tag Structure
  // TIFF Revision 6.0 Part 1 Section 2 TIFF Structure #Image File Directory
  if (TIFFGetField(tif, TIFFTAG_EXIFIFD, &exif_ifd_offset)) {
    fprintf(stderr, "Warning: EXIF extraction from TIFF is unsupported.\n");
  }
  return 1;
}

int ReadTIFF(const char* const filename,
             WebPPicture* const pic, int keep_alpha,
             Metadata* const metadata) {
  TIFF* const tif = TIFFOpen(filename, "r");
  uint32 width, height;
  uint32* raster;
  int ok = 0;
  int dircount = 1;

  if (tif == NULL) {
    fprintf(stderr, "Error! Cannot open TIFF file '%s'\n", filename);
    return 0;
  }

  while (TIFFReadDirectory(tif)) ++dircount;

  if (dircount > 1) {
    fprintf(stderr, "Warning: multi-directory TIFF files are not supported.\n"
                    "Only the first will be used, %d will be ignored.\n",
                    dircount - 1);
  }

  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
  raster = (uint32*)_TIFFmalloc(width * height * sizeof(*raster));
  if (raster != NULL) {
    if (TIFFReadRGBAImageOriented(tif, width, height, raster,
                                  ORIENTATION_TOPLEFT, 1)) {
      const int stride = width * sizeof(*raster);
      pic->width = width;
      pic->height = height;
      // TIFF data is ABGR
#ifdef __BIG_ENDIAN__
      TIFFSwabArrayOfLong(raster, width * height);
#endif
      ok = keep_alpha
         ? WebPPictureImportRGBA(pic, (const uint8_t*)raster, stride)
         : WebPPictureImportRGBX(pic, (const uint8_t*)raster, stride);
    }
    _TIFFfree(raster);
  } else {
    fprintf(stderr, "Error allocating TIFF RGBA memory!\n");
  }

  if (ok) {
    if (keep_alpha == 2) WebPCleanupTransparentArea(pic);
    if (metadata != NULL) {
      ok = ExtractMetadataFromTIFF(tif, metadata);
      if (!ok) fprintf(stderr, "Error extracting TIFF metadata!\n");
    }
  }

  TIFFClose(tif);
  return ok;
}
#else  // !WEBP_HAVE_TIFF
int ReadTIFF(const char* const filename,
             struct WebPPicture* const pic, int keep_alpha,
             struct Metadata* const metadata) {
  (void)filename;
  (void)pic;
  (void)keep_alpha;
  (void)metadata;
  fprintf(stderr, "TIFF support not compiled. Please install the libtiff "
          "development package before building.\n");
  return 0;
}
#endif  // WEBP_HAVE_TIFF

// -----------------------------------------------------------------------------
