AM_CPPFLAGS = -I$(top_srcdir)/src

libwebpllencode_la_SOURCES = backward_references.cc \
                             bit_writer.h entropy_encode.cc entropy_encode.h \
                             histogram.cc histogram.h \
                             histogram_image.cc histogram_image.h \
                             predictor.cc predictor.h \
                             webp_bit_stream.cc
libwebpllencode_la_CPPFLAGS = $(USE_EXPERIMENTAL_CODE)

noinst_HEADERS = backward_distance.h backward_references.h webp_bit_stream.h
noinst_LTLIBRARIES = libwebpllencode.la
