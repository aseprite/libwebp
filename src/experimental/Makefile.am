AM_CPPFLAGS = -I$(top_srcdir)/src

libwebpexperimental_la_SOURCES = paeth.h paeth.c \
                                 alpha.h alpha.c
libwebpexperimental_la_CPPFLAGS = $(USE_EXPERIMENTAL_CODE)
libwebpexperimental_la_LDFLAGS = -version-info 0:0:0
libwebpexperimentalinclude_HEADERS = ../webp/types.h
libwebpexperimentalincludedir = $(includedir)/webp

noinst_LTLIBRARIES = libwebpexperimental.la
