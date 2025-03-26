#
#				Makefile.am
#
# src Makefile. Process this file with automake to generate a Makefile
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	SWarp
#
#	Copyright:		(C) 2002-2023 IAP/CNRS/SorbonneU
#
#	SWarp is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
# 	(at your option) any later version.
#	SWarp is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#	You should have received a copy of the GNU General Public License
#	along with SWarp.  If not, see <http://www.gnu.org/licenses/>.
#
#	Last modified:		26/04/2023
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBDIRS			= fits wcs
bin_PROGRAMS		= swarp
swarp_SOURCES		= back.c coadd.c data.c dgeo.c field.c fitswcs.c \
			  header.c interpolate.c main.c makeit.c misc.c \
			  prefs.c projapprox.c resample.c threads.c weight.c \
			  xml.c \
			  back.h coadd.h data.h define.h dgeo.h field.h \
			  fitswcs.h globals.h header.h interpolate.h key.h \
			  misc.h preflist.h prefs.h projapprox.h resample.h \
			  threads.h types.h wcscelsys.h weight.h xml.h
swarp_LDADD		= $(srcdir)/fits/libfits.a $(srcdir)/wcs/libwcs_c.a
DATE=`date +"%Y-%m-%d"`

