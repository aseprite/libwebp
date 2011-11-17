AM_CPPFLAGS = -I$(top_srcdir)/src

libwebpexperimental_la_SOURCES = alpha.h alpha.c \
                                 paeth.h paeth.c \
                                 quant_levels.c
libwebpexperimental_la_CPPFLAGS = $(USE_EXPERIMENTAL_CODE)
libwebpexperimental_la_LDFLAGS = -version-info 0:0:0
libwebpexperimentalinclude_HEADERS = ../webp/types.h
libwebpexperimentalincludedir = $(includedir)/webp

noinst_LTLIBRARIES = libwebpexperimental.la
