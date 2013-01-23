SUBDIRS = dec dsp utils

if !DECODE_ONLY
  SUBDIRS += enc
endif
if WANT_MUX
  SUBDIRS += mux
endif
if WANT_DEMUX
  SUBDIRS += demux
endif

AM_CPPFLAGS = -I$(top_srcdir)/src
lib_LTLIBRARIES = libwebp.la

libwebp_la_SOURCES =
libwebpinclude_HEADERS =
libwebpinclude_HEADERS += webp/decode.h
if !DECODE_ONLY
libwebpinclude_HEADERS += webp/encode.h
endif
libwebpinclude_HEADERS += webp/types.h
noinst_HEADERS =
noinst_HEADERS += webp/format_constants.h

libwebp_la_LIBADD =
libwebp_la_LIBADD += dec/libwebpdecode.la
libwebp_la_LIBADD += dsp/libwebpdsp.la
if !DECODE_ONLY
libwebp_la_LIBADD += enc/libwebpencode.la
endif
libwebp_la_LIBADD += utils/libwebputils.la

# Use '-no-undefined' to declare that libwebp does not depend on any libraries
# other than the ones listed on the command line, i.e., after linking, it will
# not have unresolved symbols. Some platforms (Windows among them) require all
# symbols in shared libraries to be resolved at library creation.
libwebp_la_LDFLAGS = -no-undefined -version-info 4:1:0
libwebpincludedir = $(includedir)/webp

pkgconfig_DATA = libwebp.pc
${pkgconfig_DATA}: ${top_builddir}/config.status
