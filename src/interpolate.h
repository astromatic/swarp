/*
*				interpolate.h
*
* Include file for interpolate.c.
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
*	Last modified:		03/02/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FIELD_H_
#include "field.h"
#endif

#ifndef	_INTERPOLATE_H_
#define	_INTERPOLATE_H_

/*-------------------------------- macros -----------------------------------*/
/*------------------------------- constants ---------------------------------*/

#define	INTERP_MAXDIM		10	/* Max. number of image dimensions */
#define	INTERP_MAXKERNELWIDTH	 8	/* Max. range of kernel (pixels) */

/*--------------------------------- typedefs --------------------------------*/
typedef enum {INTERP_FLAGS, INTERP_NEARESTNEIGHBOUR, INTERP_BILINEAR,
		INTERP_LANCZOS2, INTERP_LANCZOS3, INTERP_LANCZOS4}
			interpenum;

/*-------------------------- structure definitions --------------------------*/
typedef struct ikernel
  {
  interpenum	interptype[INTERP_MAXDIM]; /* Interpolation type along axis */
  int		width[INTERP_MAXDIM];  	/* Interpol. kernel size along axis */
  int		nlines;			/* Number of kernel lines */
  PIXTYPE	*buffer;		/* Data processing buffer */
  PIXTYPE	*wbuffer;		/* Weight processing buffer */
  }	ikernelstruct;
/*----------------------- miscellaneous variables ---------------------------*/
/*-------------------------------- protos -----------------------------------*/

extern int	interpolate_ipix(fieldstruct *field, fieldstruct *wfield,
			double *pos, FLAGTYPE *outipix, FLAGTYPE *woutpix),
		interpolate_pix(fieldstruct *field, fieldstruct *wfield,
			ikernelstruct *ikernel, double *pos,
			PIXTYPE *outipix, PIXTYPE *woutpix);

extern ikernelstruct	*init_ikernel(interpenum *interptype, int naxis);

extern void	free_ikernel(ikernelstruct *ikernel);

#endif
