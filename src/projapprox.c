/*
*				projapprox.c
*
* Approximate astrometric re-projections (to speed up warping).
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

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#ifdef HAVE_MATHIMF_H
#include <mathimf.h>
#else
#include <math.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "types.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "projapprox.h"

/****** projapp_init *********************************************************
PROTO	projappstruct *projapp_init(wcsstruct *wcsin, wcsstruct *wcsout,
			double projmaxerr, int areaflag)
PURPOSE	Prepare all the necessary data for approximating a reprojection.
INPUT	Input pointer to wcs structure pointer,
	Output wcs structure pointer,
	Maximum reprojection error (pixels) allowed,
	Pixel area flag (triggers mapping of relative pixel areas if !=0).
OUTPUT	Pointer to an allocated projappstruct structure, or NULL if
	approximation failed 
NOTES	Currently limited to 2D (returns NULL otherwise).
AUTHOR	E. Bertin (IAP)
VERSION	19/07/2013
 ***/
projappstruct	*projapp_init(wcsstruct *wcsin, wcsstruct *wcsout,
			double projmaxerr, int areaflag)
  {
   projappstruct	*projapp;
   double		*projpos[NAXIS],
			rawposout[NAXIS],rawpos[NAXIS],rawposmin[NAXIS],
			wcspos[NAXIS],
			stepc[NAXIS],
			*step, *projline,*projlinet,*projappline,*projapplinet,
			*projarea,
			maxerror, cerror, defstep, stepcx, worldc;
   int			linecount[NAXIS], stepcount[NAXIS], npointsc[NAXIS],
			*npoints,
			d,i,j, naxis,naxisnmax, ngridpoints,npointstot,
			npointsmax, npointscx, stepcountx, nlinesc, swapflag;

/* The present version only works in 2D */
  if (wcsin->naxis != 2 || wcsout->naxis != 2)
    return (projappstruct *)NULL;

  QCALLOC(projapp, projappstruct, 1);

  naxis = projapp->naxis = wcsout->naxis;

/* Handle both angular and non-angular coordinates */
  projapp->lng = wcsin->lng;
  projapp->lat = wcsin->lat;

/* Find the largest dimension in pixels */
  naxisnmax = 0;
  for (d=0; d<naxis; d++)
    if (wcsout->naxisn[d] > naxisnmax)
      naxisnmax = wcsout->naxisn[d];

/* Check if lng and lat are swapped between in and out wcs (vicious idea!) */
  swapflag = (((wcsin->lng != wcsout->lng) || (wcsin->lat != wcsout->lat))
	&& (wcsin->lng != wcsin->lat) && (wcsout->lng != wcsout->lat));

/* Loop until the error is small enough */
  ngridpoints = PROJAPP_NGRIDPOINTS;
  npoints = projapp->npoints;
  step = projapp->step;
  for (maxerror = BIG; maxerror > projmaxerr; ngridpoints *= 2)
    {
/*-- Start all over again */
    for (d=0; d<naxis; d++)
       if (projapp->projpos[d])
         {
         free(projapp->projpos[d]);
         projapp->projpos[d] = NULL;
         }
    defstep = naxisnmax /(ngridpoints-1.0);
/*-- Adapt the suggested step to each dimension */
    npointstot = npointsmax = 1;
    for (d=0; d<naxis; d++)
      {
      npoints[d] = (int)(wcsout->naxisn[d] / defstep + 1.0);
      if (npoints[d]<2)
        npoints[d] = 2;
      
      step[d] = wcsout->naxisn[d]/(npoints[d]-1.0);
      npointstot *= npoints[d];
      npointsmax *= wcsout->naxisn[d] / PROJAPP_MINSTEP;
      }

/*-- if too many grid points necessary, forget it! */
    if (npointstot > npointsmax)
      {
      if (ngridpoints > PROJAPP_NGRIDPOINTS)
        warning("Astrometric approximation too inaccurate for ",
		"this re-projection");
      projapp_end(projapp);
      return (projappstruct *)NULL;
      }
    projapp->npointstot = npointstot;
/*-- Allocate memory and set the starting grid coordinate */
    for (d=0; d<naxis; d++)
      {
      QMALLOC(projapp->projpos[d], double, npointstot);
      projpos[d] = projapp->projpos[d];
      rawposout[d] = rawposmin[d] = 0.5;
      linecount[d] = stepcount[d] = 0;
      }

    if (areaflag)
      {
      QMALLOC(projapp->projarea, double, npointstot);
      projarea = projapp->projarea;
      }
    else
      projarea = NULL;		/* to avoid gcc -Wall warnings */

/*-- Fill the arrays */
    for (i=npointstot; i--;)
      {
      if (raw_to_wcs(wcsout, rawposout, wcspos) != RETURN_OK)
/*---- There is a coordinate "outside the sky" */
        {
        warning("Astrometric approximation impossible for ",
		"this re-projection");
        projapp_end(projapp);
        return (projappstruct *)NULL;
        }
      if (swapflag)
        {
        worldc = wcspos[wcsout->lat];
        wcspos[wcsout->lat] = wcspos[wcsin->lat];
        wcspos[wcsin->lat] = worldc;
        }
      if (wcs_to_raw(wcsin, wcspos, rawpos) != RETURN_OK)
/*---- There is a coordinate "outside the sky" */
        {
        warning("Astrometric approximation impossible for ",
		"this re-projection");
        projapp_end(projapp);
        return (projappstruct *)NULL;
        }

      for (d=0; d<naxis; d++)
        *(projpos[d]++) = rawpos[d];
      if (areaflag)
        *(projarea++) = wcs_scale(wcsout,rawposout) / wcs_scale(wcsin,rawpos);
      for (d=0; d<naxis; d++)
        {
        rawposout[d] = rawposmin[d] + (++stepcount[d])*step[d];
        if (++linecount[d]<npoints[d])
          break;
        else
          {
          linecount[d] = stepcount[d] = 0; /* No need to initialize it to 0! */
          rawposout[d] = rawposmin[d];
          }
        }
      }

/*-- Pre-compute secondary derivatives */
    projapp_dmap(projapp);

/*-- Now check the astrometric errors at n times the grid resolution */
    nlinesc = 1;
    for (d=0; d<naxis; d++)
      {
      rawposout[d] = rawposmin[d] = 0.5;
      stepc[d] = step[d]/PROJAPP_CHECKOVERSAMP;
      npointsc[d] = (int)(wcsout->naxisn[d] / stepc[d] + 1.0);
      if (d)
        nlinesc *= npointsc[d];
      linecount[d] = stepcount[d] = 0;
      }
    npointscx = npointsc[0];
    stepcx = stepc[0];
    QMALLOC(projline, double, naxis*npointscx);
    QMALLOC(projappline, double, naxis*npointscx);
    maxerror = 0.0;
    for (i=nlinesc; i--;)
      {
      rawposout[0] = rawposmin[0];
/*---- The approximation */
      projapp_line(projapp, rawposout, stepcx, npointscx, projappline, NULL);
/*---- The exact computation */
      projlinet = projline;
      stepcountx = 0;
      for (j=npointscx; j--; projlinet += naxis)
        {
        raw_to_wcs(wcsout, rawposout, wcspos);
        if (swapflag)
          {
          worldc = wcspos[wcsout->lat];
          wcspos[wcsout->lat] = wcspos[wcsin->lat];
          wcspos[wcsin->lat] = worldc;
          }
        wcs_to_raw(wcsin, wcspos, projlinet);
        rawposout[0] = rawposmin[0] + (++stepcountx)*stepcx;
        }
/*---- The comparison */
      projlinet = projline;
      projapplinet = projappline;
      for (j=npointscx; j--;)
        for (d=0; d<naxis; d++)
          if ((cerror=fabs(*(projapplinet++) - *(projlinet++))) > maxerror)
            maxerror = cerror;
      for (d=1; d<naxis; d++)
        {
        rawposout[d] = rawposmin[d] + (++stepcount[d])*stepc[d];
        if (++linecount[d]<npointsc[d])
          break;
        else
	  {
          linecount[d] = stepcount[d] = 0; /* No need to initialize it to 0! */
          rawposout[d] = rawposmin[d];
          }
        }
      }
    free(projline);
    free(projappline);
    }

  return projapp;
  }


/****** projapp_end ***********************************************************
PROTO	void projapp_end(projappstruct *projapp)
PURPOSE	"End" a reprojection approximation structure.
INPUT	Input pointer to projapp structure pointer,
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	03/01/2008
 ***/
void	projapp_end(projappstruct *projapp)
  {
   int	d;

  for (d=0; d<projapp->naxis; d++)
    {
    if (projapp->projpos[d])
      {
      free(projapp->projpos[d]);
      projapp->projpos[d] = NULL;
      }
    if (projapp->dprojpos2x[d])
      {
      free(projapp->dprojpos2x[d]);
      projapp->dprojpos2x[d] = NULL;
      }
    if (projapp->dprojpos2y[d])
      {
      free(projapp->dprojpos2y[d]);
      projapp->dprojpos2y[d] = NULL;
      }
    }

  free(projapp->projarea);
  free(projapp->dprojarea2x);
  free(projapp->dprojarea2y);

  free(projapp);

  return;
  }


/****** projapp_dmap *********************************************************
PROTO	void projapp_dmap(projappstruct *projapp)
PURPOSE	Prepare bicubic-splines by computing arrays of 2nd derivatives.
INPUT	Input projappstruct pointer.
OUTPUT	-.
NOTES	Currently limited to 2D vectors.
AUTHOR	E. Bertin (IAP)
VERSION	03/01/2008
 ***/
void	projapp_dmap(projappstruct *projapp)
  {
   double	*u, *map,*dmap, *mapt,*dmapt,
		temp;
   int		x,y,d, naxis, nbx,nbxm1,nby,nbym1;

  naxis = projapp->naxis;
  nbx = projapp->npoints[0];
  nbxm1 = nbx - 1;
  nby = projapp->npoints[1];
  nbym1 = nby - 1;
/*-- Computation of 2nd derivatives along x */
  QMALLOC(u, double, nbxm1);	/* temporary array */

  for (d=0; d<naxis; d++)
    {
    QMALLOC(projapp->dprojpos2x[d], double, projapp->npointstot);
    map = projapp->projpos[d];
    dmap = projapp->dprojpos2x[d];
    for (y=nby; y--;)
      {
      mapt = map+1;
      dmapt = dmap;
      map += nbx;
      dmap += nbx;
      if (nbx>1)
        {
        *dmapt = *u = 0.0;	/* "natural" lower boundary condition */
        for (x=nbxm1; --x; mapt++)
          {
          temp = -1.0/(*dmapt+4.0);
          *(++dmapt) = temp;
          temp *= *(u++) - 6.0*(*(mapt+1)+*(mapt-1)-2.0**mapt);
          *u = temp;
          }
        *(++dmapt) = 0.0;	/* "natural" upper boundary condition */
        for (x=nbx-2; x--;)
          {
          temp = *(dmapt--);
          *dmapt = (*dmapt*temp+*(u--))/6.0;
          }
        }
      else
        *dmapt = 0.0;
      }
    }

  if (projapp->projarea)
    {
/*-- The same for relative pixel areas */
    QMALLOC(projapp->dprojarea2x, double, projapp->npointstot);
    map = projapp->projarea;
    dmap = projapp->dprojarea2x;
    for (y=nby; y--;)
      {
      mapt = map+1;
      dmapt = dmap;
      map += nbx;
      dmap += nbx;
      if (nbx>1)
        {
        *dmapt = *u = 0.0;	/* "natural" lower boundary condition */
        for (x=nbxm1; --x; mapt++)
          {
          temp = -1.0/(*dmapt+4.0);
          *(++dmapt) = temp;
          temp *= *(u++) - 6.0*(*(mapt+1)+*(mapt-1)-2.0**mapt);
          *u = temp;
          }
        *(++dmapt) = 0.0;	/* "natural" upper boundary condition */
        for (x=nbx-2; x--;)
          {
          temp = *(dmapt--);
          *dmapt = (*dmapt*temp+*(u--))/6.0;
          }
        }
      else
        *dmapt = 0.0;
      }
    }

  free(u);

/* Derivatives along "y" */
  QMALLOC(u, double, nbym1);	/* temporary array */

  for (d=0; d<naxis; d++)
    {
    QMALLOC(projapp->dprojpos2y[d], double, projapp->npointstot);
    map = projapp->projpos[d];
    dmap = projapp->dprojpos2y[d];
    for (x=0; x<nbx; x++)
      {
      mapt = map++;
      dmapt = dmap++;
      if (nby>1)
        {
        *dmapt = *u = 0.0;	/* "natural" lower boundary condition */
        mapt += nbx;
        for (y=nbym1; --y; mapt+=nbx)
          {
          temp = -1.0/(*dmapt+4.0);
          *(dmapt += nbx) = temp;
          temp *= *(u++) - 6.0*(*(mapt+nbx)+*(mapt-nbx)-2.0**mapt);
          *u = temp;
          }
        *(dmapt+=nbx) = 0.0;	/* "natural" upper boundary condition */
        for (y=nby-2; y--;)
          {
          temp = *dmapt;
          dmapt -= nbx;
          *dmapt = (*dmapt*temp+*(u--))/6.0;
          }
        }
      else
        *dmapt = 0.0;
      }
    }

  if (projapp->projarea)
    {
/*-- The same for relative pixel areas */
    QMALLOC(projapp->dprojarea2y, double, projapp->npointstot);
    map = projapp->projarea;
    dmap = projapp->dprojarea2y;
    for (x=0; x<nbx; x++)
      {
      mapt = map++;
      dmapt = dmap++;
      if (nby>1)
        {
        *dmapt = *u = 0.0;	/* "natural" lower boundary condition */
        mapt += nbx;
        for (y=nbym1; --y; mapt+=nbx)
          {
          temp = -1.0/(*dmapt+4.0);
          *(dmapt += nbx) = temp;
          temp *= *(u++) - 6.0*(*(mapt+nbx)+*(mapt-nbx)-2.0**mapt);
          *u = temp;
          }
        *(dmapt+=nbx) = 0.0;	/* "natural" upper boundary condition */
        for (y=nby-2; y--;)
          {
          temp = *dmapt;
          dmapt -= nbx;
          *dmapt = (*dmapt*temp+*(u--))/6.0;
          }
        }
      else
        *dmapt = 0.0;
      }
    }

  free(u);

  return;
  }


/****** projapp_line *********************************************************
PROTO	void projapp_line(projappstruct *projapp, double *startposin,
		double step, int npos, double *posout, double *areaout)
PURPOSE	Approximate several reprojections on the same line
INPUT	Input projappstruct pointer,
	ptr to input coordinate vector array (of size naxis*nvectors),
	step along "x" (NAXIS1) dimension,
	number of input vectors.
	ptr to output coordinate vector array (where the results are written).
	ptr to output pixel area vector array (where the results are written).
OUTPUT	-.
NOTES	Pixel area are computed only if areaout != NULL.
	Currently limited to 2D vectors.
AUTHOR	E. Bertin (IAP)
VERSION	03/01/2008
 ***/
void	projapp_line(projappstruct *projapp, double *startposin, double step,
		int npos, double *posout, double *areaout)
  {
   double	*node,*nodep, *anode, *blo,*bhi,*dblo,*dbhi, *posoutt,*areaoutt,
		xstep, dx,ddx,cdx, dy,dy3,cdy,cdy3;
   int		d,i,j, xl,yl, x, ax,ax0,dax, nbx,nbxm1,nby, ylstep, nxnodes,
		naxis;

  naxis = projapp->naxis;
  nbx = projapp->npoints[0];
  nbxm1 = nbx - 1;
  nby = projapp->npoints[1];
  node = anode = NULL;		/* To avoid gcc -Wall warnings */

/*-- Prepare interpolation along x */
  if (nbx>1)
    {
    xstep = 1.0/projapp->step[0];
/*-- Reduced start x coordinate */
    dx = (startposin[0]-0.5)*xstep;
    dx -= (xl = (int)dx);
    xstep *= step;	/* input step in reduced units */
    if (xl<0)
      {
      xl = 0;
      dx -= 1.0;
      }
    else if (xl>=nbx-1)
      {
      xl = nbx<2 ? 0 : nbx-2;
      dx += 1.0;
      }
    }
  else
    {
    xl = 0;
    dx = 0.0;
    xstep = 1.0;
    }

/* Prepare interpolation along y */
  if (nby > 1)
    {
/*-- Reduced start y coordinate */
    dy = (startposin[1]-0.5)/projapp->step[1];
    dy -= (yl = (int)dy);
    if (yl<0)
      {
      yl = 0;
      dy -= 1.0;
      }
    else if (yl>=nby-1)
      {
      yl = nby<2 ? 0 : nby-2;
      dy += 1.0;
      }
    cdy = 1.0 - dy;
    dy3 = (dy*dy*dy-dy);
    cdy3 = (cdy*cdy*cdy-cdy);
    ylstep = nbx*yl;
    nxnodes = (int)(dx+npos*xstep) + 1;
    if (nxnodes > (nbx-xl))
      nxnodes = nbx - xl;
    QMALLOC(node, double, nxnodes);	/* Interpolated map */
    if (areaout)
      QMALLOC(anode, double, nxnodes);	/* Interpolated map */
    }
  else
    {
    ylstep = 0.0;			/* To avoid gcc -Wall warnings */
    dy = dy3 = cdy = cdy3 = nxnodes = 0.0;  /* To avoid gcc -Wall warnings */
    }

  for (d=0; d<naxis; d++)
    {
    if (nby > 1)
      {
/*---- Interpolation along y for each node */
      blo = projapp->projpos[d] + ylstep + xl;
      bhi = blo + nbx;
      dblo = projapp->dprojpos2y[d] + ylstep + xl;
      dbhi = dblo + nbx;
      nodep = node;
      for (x=nxnodes; x--;)
        *(nodep++) = cdy**(blo++) + dy**(bhi++) + cdy3**(dblo++)+dy3**(dbhi++);
      }
    else
      node =  projapp->projpos[d] + xl;

    if (nbx > 1)
      {
/*---- Interpolation along x */
      blo = node;
      bhi = blo + 1;
      dblo = projapp->dprojpos2x[d] + xl;
      dbhi = dblo + 1;
      ax0 = xl;
      posoutt = posout+d;
      for (i=0,j=npos; j--; i++)
        {
        ax = (int)(ddx = xl + dx + i*xstep);
        if ((int)(dax = ax-ax0) && ax<nbxm1)
          {
          blo+=dax;
          bhi+=dax;
          dblo+=dax;
          dbhi+=dax;
          ax0 = ax;
          }
        ddx -= (double)ax0;
        cdx = 1.0 - ddx;
        *posoutt = (cdx*(*blo+(cdx*cdx-1)**dblo)+ddx*(*bhi+(ddx*ddx-1)**dbhi));
        posoutt += naxis;
        }
      }
    else
      {
      posoutt = posout+d;
      for (j=npos; j--;)
        *(posout++) = *node;
      }
    }

  if (areaout)
    {
    if (nby > 1)
      {
/*---- Interpolation along y for each node */
      blo = projapp->projarea + ylstep + xl;
      bhi = blo + nbx;
      dblo = projapp->dprojarea2y + ylstep + xl;
      dbhi = dblo + nbx;
      nodep = anode;
      for (x=nxnodes; x--;)
        *(nodep++) = cdy**(blo++) + dy**(bhi++) + cdy3**(dblo++)+dy3**(dbhi++);
      }
    else
      node =  projapp->projarea + xl;

    if (nbx > 1)
      {
/*---- Interpolation along x */
      blo = anode;
      bhi = blo + 1;
      dblo = projapp->dprojarea2x + xl;
      dbhi = dblo + 1;
      ax0 = xl;
      areaoutt = areaout;
      for (i=0,j=npos; j--; i++)
        {
        ax = (int)(ddx = xl + dx + i*xstep);
        if ((int)(dax = ax-ax0) && ax<nbxm1)
          {
          blo+=dax;
          bhi+=dax;
          dblo+=dax;
          dbhi+=dax;
          ax0 = ax;
          }
        ddx -= (double)ax0;
        cdx = 1.0 - ddx;
        *(areaoutt++) = (cdx*(*blo+(cdx*cdx-1)**dblo)
			+ ddx*(*bhi+(ddx*ddx-1)**dbhi));
        }
      }
    else
      {
      areaoutt = areaout;
      for (j=npos; j--;)
        *(areaout++) = *anode;
      }
    }

  if (nby>1)
    {
    if (areaout)
      free(anode);
    free(node);
    }

  return;
  }

