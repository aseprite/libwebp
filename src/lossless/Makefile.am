SUBDIRS = common dec enc

AM_CPPFLAGS = -I$(top_srcdir)/src

libwebpll_la_SOURCES =
libwebpll_la_LIBADD = common/libwebpllcommon.la \
                      dec/libwebplldecode.la \
                      enc/libwebpllencode.la

noinst_LTLIBRARIES = libwebpll.la
