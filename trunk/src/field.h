/*
*				field.h
*
* Include file for field.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SWarp
*
*	Copyright:		(C) 2000-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		30/01/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FITSCAT_H_
#include "fits/fitscat.h"
#endif

#ifndef _FIELD_H_
#define _FIELD_H_

#define		MAXINFIELD	500000	/* Maximum number of input files */

/*------------------------------ field flags --------------------------------*/
#define		DETECT_FIELD	0x01	/* Detection */
#define		MEASURE_FIELD	0x02	/* Measurement */
#define		FLAG_FIELD	0x04	/* Flagging */
#define		RMS_FIELD	0x08	/* Weighting with std deviations */
#define		VAR_FIELD	0x10	/* Weighting with variances */
#define		WEIGHT_FIELD	0x20	/* Weighting with weights */
#define		BACKRMS_FIELD	0x40	/* Weighting from a backrms matrix */

#define		FIELD_READ	0x100	/* Field available for read access */
#define		FIELD_WRITE	0x200	/* Field available for read access */

/*--------------------------------- typedefs --------------------------------*/
typedef enum	{BACK_RELATIVE, BACK_ABSOLUTE}	backenum;

typedef struct field
  {
  char		filename[MAXCHAR];	/* image filename */
  char		*rfilename;		/* pointer to the reduced image name */
  char		hfilename[MAXCHAR];	/* header filename */
  int		headflag;		/* header found? */		
  char		ident[80];		/* field identifier (read from FITS) */
  catstruct	*cat;			/* cat structure */
  tabstruct	*tab;			/* tab structure */
/* ---- main image parameters */
  int		fieldno;		/* pos of parent ima in command line */
  int		frameno;		/* pos in Multi-extension FITS file */
  int		version;		/* filename version */
  int		width, height;		/* x,y size of the field */
  size_t	npix;			/* total number of pixels */
  int		bitpix;			/* Bits per pixel */
  double	gain;			/* conversion factor e-/ADU */
  double	fgain;			/* flux-scaled gain */
  double	saturation;		/* saturation limit in ADU */
  double	fsaturation;		/* flux scale saturation */
  double	exptime;		/* exposure time (s) */
  double	fscale;			/* relative photometric scale */
  double	fascale;		/* relative phot. scale from astrom. */
  double	ngamma;			/* normalized photo gamma */
/* ---- background parameters */
  float		*back;			/* ptr to the background map in mem */
  float		*dback;			/* ptr to the background deriv. map */
  float		*sigma;			/* ptr to the sigma map */
  float		*dsigma;		/* Ptr to the sigma deriv. map */
  int		backw, backh;		/* x,y size of a bkgnd mesh */
  int		nbackp;			/* total nb of pixels per bkgnd mesh */
  int		nbackx, nbacky;		/* x,y number of bkgnd meshes */
  int		nback;			/* total number of bkgnd meshes */
  int		nbackfx, nbackfy;	/* x,y size of bkgnd filtering mask */
  double       	backdefault;		/* default background value */
  double       	backmean;		/* median bkgnd value in image */
  double       	fbackmean;		/* flux-scaled median bkgnd  */
  double       	backsig;		/* median bkgnd rms in image */
  double       	fbacksig;		/* flux-scaled bkgnd rms */
  double	sigfac;			/* scaling RMS factor (for WEIGHTs) */
  PIXTYPE	*pix;			/* pixel data */
  FLAGTYPE	*ipix;			/* flag data */
  PIXTYPE	*backline;		/* current interpolated bkgnd line */
  backenum     	back_type;		/* background type */
/* ---- astrometric parameters */
  struct wcs	*wcs;			/* astrometric data */
  int		flags;			/* flags defining the field type */
  int		cflags;			/* flags defining the conversion */
  double	*rawmin;		/* Starting pixel for coaddition */
  double	*rawmax;		/* Ending pixel for coaddition */
/* ---- image interpolation */
  PIXTYPE	weight_thresh;		/* weight threshold */
  PIXTYPE	var_thresh;		/* variance threshold */
  struct field	*reffield;	       	/* pointer to a reference field */
/* ---- time */
  char		sdate_end[12];		/* SWarp end date */
  char		stime_end[12];		/* SWarp end time */
  int		time_diff;		/* Execution time */
  }	fieldstruct;

/*------------------------------- functions ---------------------------------*/

extern fieldstruct	*inherit_field(char *filename, fieldstruct *reffield,
					int fflags),
			*init_field(fieldstruct **infield, int ninput,
				char *filename),
			*load_field(catstruct *cat, int frameno, int fieldno);

extern void		end_field(fieldstruct *field),
			printinfo_field(fieldstruct *field,
					fieldstruct *wfield),
			scale_field(fieldstruct *field, fieldstruct *reffield,
					int scaleflag);
#endif
