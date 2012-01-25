AM_CPPFLAGS = -I$(top_srcdir)/src

libwebplldecode_la_SOURCES = bit_reader.h constants.h \
                             decode.cc huffman.h image_transform.h
libwebplldecode_la_CPPFLAGS = $(USE_EXPERIMENTAL_CODE)

noinst_HEADERS = decode.h
noinst_LTLIBRARIES = libwebplldecode.la
