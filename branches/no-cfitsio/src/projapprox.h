/*
*				projapprox.h
*
* Include file for projapprox.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SWarp
*
*	Copyright:		(C) 2003-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		19/07/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

#ifndef	_PROJAPPROX_H_
#define	_PROJAPPROX_H_

/*-------------------------------- macros -----------------------------------*/
/*------------------------------- constants ---------------------------------*/

#define	PROJAPP_NGRIDPOINTS    	16	/* Min number of control nodes/dim */
#define	PROJAPP_MINSTEP		16	/* Min dist between nodes (pixels) */
#define	PROJAPP_CHECKOVERSAMP	4	/* node oversampling for checking */

/*--------------------------------- typedefs --------------------------------*/
/*-------------------------- structure definitions --------------------------*/

typedef struct projapp
  {
  int		naxis;		/* Number of dimensions */
  int		lng,lat;	/* Indices for longitude and latitude */
  int		npoints[NAXIS];	/* Number of nodes per dimension */
  int		npointstot;	/* Total number of nodes */
  double	step[NAXIS];	/* Step (in pixels) between each node */
  double	*projpos[NAXIS];/* coordinate value at each node */
  double	*dprojpos2x[NAXIS];/* second derivative along x at each node */
  double	*dprojpos2y[NAXIS];/* second derivative along y at each node */
  double	*projarea;	/* coordinate value at each node */
  double	*dprojarea2x;	/* second derivative along x at each node */
  double	*dprojarea2y;	/* second derivative along y at each node */
  }		projappstruct;


/*----------------------- miscellaneous variables ---------------------------*/
projappstruct	*projapp_init(wcsstruct *wcsin, wcsstruct *wcsout,
			double projmaxerr, int areaflag);

void		projapp_dmap(projappstruct *projapp),
		projapp_end(projappstruct *projapp),
		projapp_line(projappstruct *projapp, double *startposin,
			double step, int npos, double *posout, double *areaout);

/*-------------------------------- protos -----------------------------------*/

#endif
