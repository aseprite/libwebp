AM_CPPFLAGS = -I$(top_srcdir)/src

noinst_PROGRAMS = png2webpll webpll2png

png2webpll_SOURCES = png2webpll.cc
png2webpll_CPPFLAGS = $(AM_CPPFLAGS) $(PNG_INCLUDES) $(USE_EXPERIMENTAL_CODE)
png2webpll_LDADD = ../libwebpll.la $(PNG_LIBS)

webpll2png_SOURCES = webpll2png.cc
webpll2png_CPPFLAGS = $(AM_CPPFLAGS) $(PNG_INCLUDES) $(USE_EXPERIMENTAL_CODE)
webpll2png_LDADD = ../libwebpll.la $(PNG_LIBS)
