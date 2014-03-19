/*
*				preflist.h
*
* Configuration keyword definitions.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SWarp
*
*	Copyright:		(C) 2000-2014 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	SWarp is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*	SWarp is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with SWarp. If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		10/03/2014
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "key.h"

#ifndef _FIELD_H_
#include "field.h"
#endif
#ifndef _INTERPOLATE_H_
#include "interpolate.h"
#endif
#ifndef _PREFS_H_
#include "prefs.h"
#endif

#ifndef _XML_H_
#include "xml.h"
#endif

#ifdef	USE_THREADS
#define	THREADS_PREFMAX	THREADS_NMAX
#else
#define	THREADS_PREFMAX	65535
#endif

int idummy;

pkeystruct key[] =
 {
  {"BACK_DEFAULT", P_FLOATLIST, prefs.back_default, 1,7, -BIG, BIG,
   {""}, 1, MAXINFIELD, &prefs.nback_default},
  {"BACK_FILTTHRESH", P_FLOAT, &prefs.back_fthresh, 0,0, -BIG, BIG},
  {"BACK_SIZE", P_INTLIST, prefs.back_size, 1,2000000000, 0.0,0.0,
   {""}, 1, MAXINFIELD, &prefs.nback_size},
  {"BACK_FILTERSIZE", P_INTLIST, prefs.back_fsize, 1,7, 0.0,0.0,
   {""}, 1, MAXINFIELD, &prefs.nback_fsize},
  {"BACK_TYPE", P_KEYLIST, prefs.back_type, 0,0, 0.0,0.0,
   {"AUTO", "MANUAL", ""},
   1, MAXINFIELD, &prefs.nback_type},
  {"BLANK_BADPIXELS", P_BOOL, &prefs.blank_flag},
  {"CELESTIAL_TYPE", P_KEY, &prefs.celsys_type, 0,0, 0.0,0.0,
   {"NATIVE", "PIXEL", "EQUATORIAL", "GALACTIC", "ECLIPTIC", "SUPERGALACTIC",
	""}},
  {"CENTER", P_STRINGLIST, prefs.image_center, 0,0, 0.0, 0.0,
   {""}, 0, INTERP_MAXDIM, &prefs.nimage_center},
  {"CENTER_TYPE", P_KEYLIST, prefs.center_type, 0,0, 0.0,0.0,
   {"MANUAL", "ALL", "MOST",""},
   1, INTERP_MAXDIM, &prefs.ncenter_type},
  {"CLIP_AMPFRAC", P_FLOAT, &prefs.clip_ampfrac, 0, 0, 0.0, BIG},
  {"CLIP_LOGNAME", P_STRING, prefs.clip_logname},
  {"CLIP_SIGMA",   P_FLOAT, &prefs.clip_sigma,   0, 0, 0.0, BIG},
  {"CLIP_WRITELOG", P_BOOL, &prefs.clip_logflag},
  {"COMBINE", P_BOOL, &prefs.combine_flag},
  {"COMBINE_BUFSIZE", P_INT, &prefs.coaddbuf_size, 1, 16384*1024},
  {"COMBINE_TYPE", P_KEY, &prefs.coadd_type, 0,0, 0.0,0.0,
   {"MEDIAN", "AVERAGE", "MIN", "MAX", "WEIGHTED", "CLIPPED",
	"CHI_OLD", "CHI-MODE", "CHI-MEAN", "SUM",
	"WEIGHTED_WEIGHT", "MEDIAN_WEIGHT",
	"AND", "NAND", "OR", "NOR", ""}},
  {"COPY_KEYWORDS", P_STRINGLIST, prefs.copy_keywords, 0,0, 0.0, 0.0,
   {""}, 0, 1024, &prefs.ncopy_keywords},
  {"DELETE_TMPFILES", P_BOOL, &prefs.removetmp_flag},
  {"FSCALASTRO_TYPE", P_KEY, &prefs.fscalastro_type, 0,0, 0.0,0.0,
   {"NONE", "FIXED", "VARIABLE", ""}},
  {"FSCALE_DEFAULT", P_FLOATLIST, prefs.fscale_default, 0,0, -BIG, BIG,
   {""}, 1, MAXINFIELD, &prefs.nfscale_default},
  {"FSCALE_KEYWORD", P_STRING, prefs.fscale_keyword},
  {"GAIN_DEFAULT", P_FLOATLIST, prefs.gain_default, 0,0, 0.0, BIG,
   {""}, 1, MAXINFIELD, &prefs.ngain_default},
  {"GAIN_KEYWORD", P_STRING, prefs.gain_keyword},
  {"SATLEV_KEYWORD", P_STRING, prefs.sat_keyword},
  {"HEADER_ONLY", P_BOOL, &prefs.headeronly_flag},
  {"HEADER_SUFFIX", P_STRING, prefs.head_suffix},
  {"IMAGEOUT_NAME", P_STRING, prefs.outfield_name},
  {"IMAGE_SIZE", P_INTLIST, prefs.image_size, 0, 2000000000, 0.0, 0.0,
   {""}, 1, INTERP_MAXDIM, &prefs.nimage_size},
  {"INTERPOLATE", P_BOOLLIST, prefs.interp_flag, 0,0, 0.0,0.0,
   {""}, 1, MAXINFIELD, &prefs.ninterp_flag},
  {"MEM_MAX", P_INT, &prefs.mem_max, 1, 1000000000},
  {"NNODES", P_INT, &prefs.nnodes, 1, 65535},
  {"NOPENFILES_MAX", P_INT, &prefs.nopenfiles_max, 0, 1000000000},
  {"NTHREADS", P_INT, &prefs.nthreads, 0, THREADS_PREFMAX},
  {"NODE_INDEX", P_INT, &prefs.node_index, -1, 65534},
  {"OVERSAMPLING", P_INTLIST, prefs.oversamp, 0, 2000000000, 0.0,0.0,
   {""}, 1, INTERP_MAXDIM, &prefs.noversamp},
  {"PIXELSCALE_TYPE", P_KEYLIST, prefs.pixscale_type, 0,0, 0.0,0.0,
   {"MANUAL", "MIN", "MAX", "MEDIAN", "FIT", ""},
   1, INTERP_MAXDIM, &prefs.npixscale_type},
  {"PIXEL_SCALE", P_FLOATLIST, prefs.pixscale, 0,0, 0.0, BIG,
   {""}, 1, INTERP_MAXDIM, &prefs.npixscale},
  {"PROJECTION_ERR", P_FLOATLIST, prefs.proj_err, 0,0, 0.0, 1.0,
   {""}, 1, MAXINFIELD, &prefs.nproj_err},
  {"PROJECTION_TYPE", P_STRING, prefs.projection_name},
  {"RESAMPLE", P_BOOL, &prefs.resample_flag},
  {"RESAMPLE_DIR", P_STRING, prefs.resampdir_name},
  {"RESAMPLE_SUFFIX", P_STRING, prefs.resamp_suffix},
  {"RESAMPLING_TYPE", P_KEYLIST, prefs.resamp_type, 0,0, 0.0,0.0,
   {"FLAGS", "NEAREST", "BILINEAR", "LANCZOS2", "LANCZOS3", "LANCZOS4", ""},
   1, INTERP_MAXDIM, &prefs.nresamp_type},
  {"RESCALE_WEIGHTS", P_BOOLLIST, prefs.wscale_flag, 0,0, 0.0,0.0,
   {""}, 1, MAXINFIELD, &prefs.nwscale_flag},
  {"SATLEV_DEFAULT", P_FLOATLIST, prefs.sat_default, 0,0, -BIG, BIG,
   {""}, 1, MAXINFIELD, &prefs.nsat_default},
  {"SUBTRACT_BACK", P_BOOLLIST, prefs.subback_flag, 0,0, 0.0,0.0,
   {""}, 1, MAXINFIELD, &prefs.nsubback_flag},
  {"VERBOSE_TYPE", P_KEY, &prefs.verbose_type, 0,0, 0.0,0.0,
   {"QUIET", "LOG", "NORMAL", "FULL", ""}},
  {"VMEM_DIR", P_STRING, prefs.swapdir_name},
  {"VMEM_MAX", P_INT, &prefs.vmem_max, 1, 1000000000},
  {"WEIGHT_IMAGE", P_STRINGLIST, prefs.inwfield_name, 0,0, 0.0,0.0,
   {""}, 0, MAXINFIELD, &prefs.ninwfield},
  {"WEIGHTOUT_NAME", P_STRING, prefs.outwfield_name},
  {"WEIGHT_SUFFIX", P_STRING, prefs.weight_suffix},
  {"WEIGHT_THRESH", P_FLOATLIST, prefs.weight_thresh, 0,0, 0.0, BIG,
   {""}, 0, MAXINFIELD, &prefs.nweight_thresh},
  {"WEIGHT_TYPE", P_KEYLIST, prefs.weight_type, 0,0, 0.0,0.0,
   {"NONE", "BACKGROUND", "MAP_RMS", "MAP_VARIANCE", "MAP_WEIGHT",""},
   1, MAXINFIELD, &prefs.nweight_type},
  {"WRITE_FILEINFO", P_BOOL, &prefs.writefileinfo_flag},
  {"WRITE_XML", P_BOOL, &prefs.xml_flag},
  {"XML_NAME", P_STRING, prefs.xml_name},
  {"XSL_URL", P_STRING, prefs.xsl_name},
  {""}
 };

char			keylist[sizeof(key)/sizeof(pkeystruct)][32];
const char		notokstr[] = {" \t=,;\n\r\""};

char *default_prefs[] =
 {
"# Default configuration file for " BANNER " " MYVERSION,
"# EB " DATE,
"#",
"#----------------------------------- Output -----------------------------------",
"IMAGEOUT_NAME          coadd.fits      # Output filename",
"WEIGHTOUT_NAME       coadd.weight.fits # Output weight-map filename",
" ",
"HEADER_ONLY            N               # Only a header as an output file (Y/N)?",
"HEADER_SUFFIX          .head           # Filename extension for additional headers",
" ",
"#------------------------------- Input Weights --------------------------------",
" ",
"WEIGHT_TYPE            NONE            # BACKGROUND,MAP_RMS,MAP_VARIANCE",
"                                       # or MAP_WEIGHT",
"*RESCALE_WEIGHTS        Y               # Rescale input weights/variances (Y/N)?",
"WEIGHT_SUFFIX          .weight.fits    # Suffix to use for weight-maps",
"WEIGHT_IMAGE                           # Weightmap filename if suffix not used",
"                                       # (all or for each weight-map)",
"*WEIGHT_THRESH                         # Bad pixel weight-threshold",
" ",
"#------------------------------- Co-addition ----------------------------------",
" ",
"COMBINE                Y               # Combine resampled images (Y/N)?",
"COMBINE_TYPE           MEDIAN          # MEDIAN,AVERAGE,MIN,MAX,WEIGHTED,CLIPPED",
"                                       # CHI-OLD,CHI-MODE,CHI-MEAN,SUM,",
"                                       # WEIGHTED_WEIGHT,MEDIAN_WEIGHT,",
"                                       # AND,NAND,OR or NOR",
"*CLIP_AMPFRAC           0.3             # Fraction of flux variation allowed",
"*                                       # with clipping",
"*CLIP_SIGMA             4.0             # RMS error multiple variation allowed",
"*                                       # with clipping",
"*CLIP_WRITELOG          N               # Write output file with coordinates of",
"*                                       # clipped pixels (Y/N) ",
"*CLIP_LOGNAME           clipped.log     # Name of output file with coordinates",
"*                                       # of clipped pixels",
"*BLANK_BADPIXELS        N               # Set to 0 pixels having a weight of 0",
" ",
"#-------------------------------- Astrometry ----------------------------------",
" ",
"CELESTIAL_TYPE         NATIVE          # NATIVE, PIXEL, EQUATORIAL,",
"                                       # GALACTIC,ECLIPTIC, or SUPERGALACTIC",
"PROJECTION_TYPE        TAN             # Any WCS projection code or NONE",
"PROJECTION_ERR         0.001           # Maximum projection error (in output",
"                                       # pixels), or 0 for no approximation",
"CENTER_TYPE            ALL             # MANUAL, ALL or MOST",
"CENTER         00:00:00.0, +00:00:00.0 # Coordinates of the image center",
"PIXELSCALE_TYPE        MEDIAN          # MANUAL,FIT,MIN,MAX or MEDIAN",
"PIXEL_SCALE            0.0             # Pixel scale",
"IMAGE_SIZE             0               # Image size (0 = AUTOMATIC)",
" ",
"#-------------------------------- Resampling ----------------------------------",
" ",
"RESAMPLE               Y               # Resample input images (Y/N)?",
"RESAMPLE_DIR           .               # Directory path for resampled images",
"RESAMPLE_SUFFIX        .resamp.fits    # filename extension for resampled images",
" ",
"RESAMPLING_TYPE        LANCZOS3        # NEAREST,BILINEAR,LANCZOS2,LANCZOS3",
"                                       # LANCZOS4 (1 per axis) or FLAGS",
"OVERSAMPLING           0               # Oversampling in each dimension",
"                                       # (0 = automatic)",
"INTERPOLATE            N               # Interpolate bad input pixels (Y/N)?",
"                                       # (all or for each image)",
" ",
"FSCALASTRO_TYPE        FIXED           # NONE,FIXED, or VARIABLE",
"FSCALE_KEYWORD         FLXSCALE        # FITS keyword for the multiplicative",
"                                       # factor applied to each input image",
"FSCALE_DEFAULT         1.0             # Default FSCALE value if not in header",
" ",
"GAIN_KEYWORD           GAIN            # FITS keyword for effect. gain (e-/ADU)",
"GAIN_DEFAULT           0.0             # Default gain if no FITS keyword found",
"*                                       # 0 = infinity (all or for each image)",
"*SATLEV_KEYWORD         SATURATE        # FITS keyword for saturation level (ADU)",
"*SATLEV_DEFAULT         50000.0         # Default saturation if no FITS keyword",
" ",
"#--------------------------- Background subtraction ---------------------------",
" ",
"SUBTRACT_BACK          Y               # Subtraction sky background (Y/N)?",
"                                       # (all or for each image)",
" ",
"BACK_TYPE              AUTO            # AUTO or MANUAL",
"                                       # (all or for each image)",
"BACK_DEFAULT           0.0             # Default background value in MANUAL",
"                                       # (all or for each image)",
"BACK_SIZE              128             # Background mesh size (pixels)",
"                                       # (all or for each image)",
"BACK_FILTERSIZE        3               # Background map filter range (meshes)",
"                                       # (all or for each image)",
"*BACK_FILTTHRESH        0.0             # Threshold above which the background-",
"*                                       # map filter operates",
" ",
"#------------------------------ Memory management -----------------------------",
" ",
"VMEM_DIR               .               # Directory path for swap files",
"VMEM_MAX               2047            # Maximum amount of virtual memory (MB)",
"MEM_MAX                256             # Maximum amount of usable RAM (MB)",
"COMBINE_BUFSIZE        256             # RAM dedicated to co-addition(MB)",
" ",
"#------------------------------ Miscellaneous ---------------------------------",
" ",
"DELETE_TMPFILES        Y               # Delete temporary resampled FITS files",
"                                       # (Y/N)?",
"COPY_KEYWORDS          OBJECT          # List of FITS keywords to propagate",
"                                       # from the input to the output headers",
"WRITE_FILEINFO         N               # Write information about each input",
"                                       # file in the output image header?",
"WRITE_XML              Y               # Write XML file (Y/N)?",
"XML_NAME               swarp.xml       # Filename for XML output",
"*XSL_URL                " XSL_URL,
"*                                       # Filename for XSL style-sheet",
"VERBOSE_TYPE           NORMAL          # QUIET,LOG,NORMAL, or FULL",
"*NNODES                 1               # Number of nodes (for clusters)",
"*NODE_INDEX             0               # Node index (for clusters)",
" ",
#ifdef USE_THREADS
"NTHREADS               0               # Number of simultaneous threads for",
"                                       # the SMP version of " BANNER,
"                                       # 0 = automatic",
#else
"NTHREADS               1               # 1 single thread",
#endif
"*NOPENFILES_MAX         512             # Maximum number of files opened by "
					BANNER,
""
 };

