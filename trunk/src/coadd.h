/*
*				coadd.h
*
* Include file for coadd.c.
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
		COADD_WEIGHTED, COADD_CLIPPED,
		COADD_CHI_OLD, COADD_CHI_MODE, COADD_CHI_MEAN,
		COADD_SUM, COADD_WEIGHTED_WEIGHT, COADD_MEDIAN_WEIGHT,
		COADD_AND, COADD_NAND, COADD_OR, COADD_NOR}
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

extern int	coadd_fields(fieldstruct **infield, fieldstruct **inwfield,
			int ninput, fieldstruct *outfield,
			fieldstruct *outwfield,
			coaddenum coaddtype, PIXTYPE wthresh),
		max_clique(unsigned int *array, int nnode, int **max);
extern void	coadd_movedata(PIXTYPE *linebuf, PIXTYPE *multibuf,
			unsigned int *multiobuf, unsigned int *multinbuf,
			int npix, int step, int oid),
		coadd_moveidata(FLAGTYPE *lineibuf, FLAGTYPE *multiibuf,
			unsigned int *multinbuf, int npix, int step),
		coadd_movewdata(PIXTYPE *linebuf, PIXTYPE *multiwbuf,
			unsigned int *multinbuf, int npix, int step),
		coadd_movewidata(FLAGTYPE *lineibuf, FLAGTYPE *multiwibuf,
			unsigned int *multinbuf, int npix, int step);
#endif
