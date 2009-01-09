/*
 				interpolate.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	Swarp
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for interpolate.c
*
*	Last modify:	27/04/2003
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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
typedef enum {INTERP_NEARESTNEIGHBOUR, INTERP_BILINEAR, INTERP_LANCZOS2,
		INTERP_LANCZOS3, INTERP_LANCZOS4}	interpenum;

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

extern int	interpolate_pix(fieldstruct *field, fieldstruct *wfield,
			ikernelstruct *ikernel, double *pos,
			PIXTYPE *pixout, PIXTYPE *wpixout);

extern ikernelstruct	*init_ikernel(interpenum *interptype, int naxis);

extern void	free_ikernel(ikernelstruct *ikernel);

#endif
