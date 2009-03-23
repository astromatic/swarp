/*
 				resample.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	Swarp
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for resample.c
*
*	Last modify:	14/04/2003
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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
