/*
 				coadd.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SWarp
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for coadd.c
*
*	Last modify:	05/08/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _FITSCAT_H_
#include "fits/fitscat.h"
#endif

#ifndef _FIELD_H_
#include "field.h"
#endif

#ifndef	_COADD_H_
#define	_COADD_H_

/*-------------------------------- macros -----------------------------------*/
/*------------------------------- constants ---------------------------------*/
#define	COADDFLAG_OPEN		0x01
#define	COADDFLAG_FINISHED	0x02

/*--------------------------------- typedefs --------------------------------*/
typedef enum {COADD_MEDIAN, COADD_AVERAGE, COADD_MIN, COADD_MAX,
		COADD_WEIGHTED, COADD_CHI2, COADD_SUM}
			coaddenum;	/* Coaddition type */

/*-------------------------- structure definitions --------------------------*/
typedef struct coaddact
  {
  int	line;
  int	fieldno;
  enum {COADDACT_OPEN, COADDACT_CLOSE, COADDACT_LOAD}	com;
  }	coaddactstruct;

/*----------------------- miscellaneous variables ---------------------------*/

/*-------------------------------- protos -----------------------------------*/

extern PIXTYPE	hmedian(PIXTYPE *arr, int n);
extern int	coadd_fields(fieldstruct **infield, fieldstruct **inwfield,
			int ninput, fieldstruct *outfield,
			fieldstruct *outwfield,
			coaddenum coaddtype, PIXTYPE wthresh),
		max_clique(unsigned int *array, int nnode, int **max);
extern void	coadd_movedata(PIXTYPE *linebuf, PIXTYPE *multibuf,
			unsigned int *multinbuf, int npix, int step),
		coadd_movewdata(PIXTYPE *linebuf, PIXTYPE *multiwbuf,
			unsigned int *multinbuf, int npix, int step);
#endif
