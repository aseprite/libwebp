AM_CPPFLAGS = -I$(top_srcdir)/src

bin_PROGRAMS = dwebp cwebp webpmux png2webp webp2png

dwebp_SOURCES = dwebp.c stopwatch.h
dwebp_CPPFLAGS = $(AM_CPPFLAGS) $(PNG_INCLUDES) $(JPEG_INCLUDES) $(USE_EXPERIMENTAL_CODE)
dwebp_LDADD = ../src/libwebp.la $(PNG_LIBS) $(JPEG_LIBS)

cwebp_SOURCES = cwebp.c stopwatch.h
cwebp_CPPFLAGS = $(AM_CPPFLAGS) $(PNG_INCLUDES) $(JPEG_INCLUDES) $(USE_EXPERIMENTAL_CODE)
cwebp_LDADD = ../src/libwebp.la $(PNG_LIBS) $(JPEG_LIBS)

webpmux_SOURCES = webpmux.c
webpmux_CPPFLAGS = $(AM_CPPFLAGS) $(USE_EXPERIMENTAL_CODE)
webpmux_LDADD = ../src/mux/libwebpmux.la ../src/libwebp.la

png2webp_SOURCES = png2webp.c
png2webp_CPPFLAGS = $(AM_CPPFLAGS) $(PNG_INCLUDES)
png2webp_LDADD = ../src/libwebp.la $(PNG_LIBS) -lm

webp2png_SOURCES = webp2png.c
webp2png_CPPFLAGS = $(AM_CPPFLAGS) $(PNG_INCLUDES)
webp2png_LDADD = ../src/libwebp.la $(PNG_LIBS) -lm

