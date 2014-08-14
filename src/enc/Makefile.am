noinst_LTLIBRARIES = libwebpencode.la  libwebpencode_avx2.la

libwebpencode_avx2_la_SOURCES =
libwebpencode_avx2_la_SOURCES += quant_avx2.c
libwebpencode_avx2_la_CPPFLAGS = $(libwebpencode_la_CPPFLAGS)
libwebpencode_avx2_la_CFLAGS = $(AM_CFLAGS) $(AVX2_FLAGS)

libwebpencode_la_SOURCES =
libwebpencode_la_SOURCES += alpha.c
libwebpencode_la_SOURCES += analysis.c
libwebpencode_la_SOURCES += backward_references.c
libwebpencode_la_SOURCES += config.c
libwebpencode_la_SOURCES += cost.c
libwebpencode_la_SOURCES += cost.h
libwebpencode_la_SOURCES += filter.c
libwebpencode_la_SOURCES += frame.c
libwebpencode_la_SOURCES += histogram.c
libwebpencode_la_SOURCES += iterator.c
libwebpencode_la_SOURCES += near_lossless.c
libwebpencode_la_SOURCES += picture.c
libwebpencode_la_SOURCES += picture_csp.c
libwebpencode_la_SOURCES += picture_psnr.c
libwebpencode_la_SOURCES += picture_rescale.c
libwebpencode_la_SOURCES += picture_tools.c
libwebpencode_la_SOURCES += quant.c
libwebpencode_la_SOURCES += syntax.c
libwebpencode_la_SOURCES += token.c
libwebpencode_la_SOURCES += tree.c
libwebpencode_la_SOURCES += vp8enci.h
libwebpencode_la_SOURCES += vp8l.c
libwebpencode_la_SOURCES += webpenc.c


libwebpencodeinclude_HEADERS =
libwebpencodeinclude_HEADERS += quant.h
libwebpencodeinclude_HEADERS += ../webp/encode.h
libwebpencodeinclude_HEADERS += ../webp/types.h
noinst_HEADERS =
noinst_HEADERS += ../webp/format_constants.h

libwebpencode_la_LDFLAGS = -lm
libwebpencode_la_CPPFLAGS = $(USE_EXPERIMENTAL_CODE)
libwebpencodeincludedir = $(includedir)/webp

libwebpencode_la_CPPFLAGS = $(USE_EXPERIMENTAL_CODE) $(USE_SWAP_16BIT_CSP)
libwebpencode_la_LDFLAGS = -lm
libwebpencode_la_LIBADD = libwebpencode_avx2.la
