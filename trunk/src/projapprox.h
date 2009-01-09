/*
 				projapprox.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	Swarp
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for projapprox.c
*
*	Last modify:	03/01/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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
			double projmaxerr, int areaflag, double meanscale);

void		projapp_dmap(projappstruct *projapp),
		projapp_end(projappstruct *projapp),
		projapp_line(projappstruct *projapp, double *startposin,
			double step, int npos, double *posout, double *areaout);

/*-------------------------------- protos -----------------------------------*/

#endif
