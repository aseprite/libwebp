AM_CPPFLAGS = -I$(top_srcdir)/src

libwebpllencode_la_SOURCES = backward_references.c \
                             bit_writer.h entropy_encode.c entropy_encode.h \
                             histogram.c histogram.h \
                             histogram_image.c histogram_image.h \
                             predictor.c predictor.h \
                             webp_bit_stream.cc
libwebpllencode_la_CPPFLAGS = $(USE_EXPERIMENTAL_CODE)

noinst_HEADERS = backward_references.h
noinst_LTLIBRARIES = libwebpllencode.la
