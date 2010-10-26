/*
*				resample.h
*
* Include file for resample.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SWarp
*
*	Copyright:		(C) 2000-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		26/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FITSCAT_H_
#include "fits/fitscat.h"
#endif

#ifndef _FIELD_H_
#include "field.h"
#endif

#ifndef _INTERPOLATE_H_
#include "interpolate.h"
#endif

#ifndef	_RESAMPLE_H_
#define	_RESAMPLE_H_

/*-------------------------------- macros -----------------------------------*/
/*------------------------------- constants ---------------------------------*/

#define	INTERP_MAXDIM		10	/* Max. number of image dimensions */
#define	INTERP_MAXKERNELWIDTH	 8	/* Max. range of kernel (pixels) */

/*--------------------------------- typedefs --------------------------------*/
/*-------------------------- structure definitions --------------------------*/

/*----------------------- miscellaneous variables ---------------------------*/

/*-------------------------------- protos -----------------------------------*/

#ifdef USE_THREADS
extern void	cancel_resample_threads(void);
#endif

extern void	resample_field(fieldstruct **pinfield, fieldstruct **pinwfield,
			fieldstruct *outfield, fieldstruct *outwfield,
			interpenum *interptype);
#endif
