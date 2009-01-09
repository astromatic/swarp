# Program Makefile for SWarp
# Copyright (C) 2002-2006 Emmanuel Bertin.
SUBDIRS			= fits wcs
bin_PROGRAMS		= swarp
swarp_SOURCES		= back.c coadd.c data.c field.c fitswcs.c header.c \
			  interpolate.c main.c makeit.c prefs.c projapprox.c \
			  resample.c threads.c weight.c xml.c \
			  back.h coadd.h data.h define.h field.h fitswcs.h \
			  globals.h header.h interpolate.h key.h preflist.h \
			  prefs.h projapprox.h resample.h threads.h types.h \
			  wcscelsys.h weight.h xml.h
swarp_LDADD		= $(top_builddir)/src/fits/libfits.a $(top_builddir)/src/wcs/libwcs_c.a
DATE=`date +"%Y-%m-%d"`
