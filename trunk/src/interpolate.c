 /*
 				interpolate.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SWarp
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Function related to vignet manipulations.
*
*	Last modify:	07/01/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#ifdef HAVE_MATHIMF_H
#include <mathimf.h>
#else
#define _GNU_SOURCE
#include <math.h>
#endif
#include	<stdio.h>
#include	<stdlib.h>

#include	"define.h"
#include	"types.h"
#include	"globals.h"
#include	"fits/fitscat.h"
#include	"interpolate.h"

static void	make_kernel(double pos, double *kernel, interpenum interptype);

int		interp_kernwidth[5]={1,2,4,6,8};

/****** interpolate_pix *******************************************************
PROTO	int interpolate_pix(fieldstruct *field, fieldstruct *wfield,
		ikernelstruct ikernel, double *pos,
		PIXTYPE *pixout, PIXTYPE *wpixout)
PURPOSE	Interpolate pixel data through sinc interpolation.
INPUT	Field structure pointer,
	Weight field structure pointer,
	Interpolation kernel structure pointer,
	Position vector,
	Pointer to the output pixel,
	Pointer to the output weight.
OUTPUT	RETURN_OK if pixel falls within the input frame, RETURN_ERROR
	otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	31/01/2006
 ***/
int	interpolate_pix(fieldstruct *field, fieldstruct *wfield,
		ikernelstruct *ikernel, double *pos, PIXTYPE *outpix,
		PIXTYPE *woutpix)

  {
   PIXTYPE		*pixin,*pixout,
			pixval;
   double		dpos[INTERP_MAXDIM],
			kernel_vector[INTERP_MAXKERNELWIDTH];
   double		*kvector,
			max, val;
   long			step[INTERP_MAXDIM];
   long			start, fac;
   int			linecount[INTERP_MAXDIM];
   int			*naxisn,
			i,j,n, ival, naxis, nlines, kwidth,width, badpixflag;

  naxis = field->tab->naxis;
  naxisn = field->tab->naxisn;
  width = field->width;
  start = 0;
  fac = 1;
  for (n=0; n<naxis; n++)
    {
    val = *(pos++);
    width = *(naxisn++);
/*-- Get the integer part of the current coordinate or nearest neighbour */
    ival = (ikernel->interptype[n]==INTERP_NEARESTNEIGHBOUR)?
					(int)(val-0.50001):(int)val;
/*-- Store the fractional part of the current coordinate */
    dpos[n] = val - ival;
/*-- Check if interpolation start/end exceed image boundary... */
    kwidth = ikernel->width[n];
    ival-=kwidth/2;
    if (ival<0 || ival+kwidth<=0 || ival+kwidth>width)
      {
      *outpix = 0.0;
      if (woutpix)
        *woutpix = BIG;
      return RETURN_ERROR;
      }
/*-- Update starting pointer */
    start += ival*fac;
/*-- Update step between interpolated regions */
    step[n] = fac*(width-kwidth);
    linecount[n] = 0.0;
    fac *= width;
    }

/* Update Interpolation kernel vectors */
  make_kernel(*dpos, kernel_vector, ikernel->interptype[0]);
  kwidth = ikernel->width[0];
  nlines = ikernel->nlines;
/* First step: interpolate along NAXIS1 from the data themselves */
  badpixflag = 0;
  pixin = field->pix+start;
  pixout = ikernel->buffer;
  for (j=nlines; j--;)
    {
    val = 0.0;
    kvector = kernel_vector;
    for (i=kwidth; i--;)
      if ((pixval = *(pixin++))>-BIG)
        val += *(kvector++)*pixval;
      else
        {
        badpixflag = 1;
        kvector++;
        }
    *(pixout++) = val;
    for (n=1; n<naxis; n++)
      {
      pixin+=step[n-1];
      if (++linecount[n]<ikernel->width[n])
        break;
      else
        linecount[n] = 0;	/* No need to initialize it to 0! */
      }
    }

/* Now the weight (variance, in fact) map */
  if (wfield)
    {
    pixin = wfield->pix+start;
    pixout = ikernel->wbuffer;
    for (j=nlines; j--;)
      {
      max = 0.0;
      for (i=kwidth; i--;)
        if ((val = *(pixin++))>max)
          max = val;
      *(pixout++) = max;
      for (n=1; n<naxis; n++)
        {
        pixin+=step[n-1];
        if (++linecount[n]<ikernel->width[n])
          break;
        else
          linecount[n] = 0;
        }
      }
    }

/* Second step: interpolate along other axes from the interpolation buffer */
  for (n=1; n<naxis; n++)
    {
    make_kernel(dpos[n], kernel_vector, ikernel->interptype[n]);
    kwidth = ikernel->width[n];
    pixout = pixin = ikernel->buffer;
    for (j = (nlines/=kwidth); j--;)
      {
      val = 0.0;
      kvector = kernel_vector;
      for (i=kwidth; i--;)
        val += *(kvector++)**(pixin++);
      *(pixout++) = val;
     }
/*-- Now the weight (variance, in fact) map */
    if (wfield)
      {
      pixout = pixin = ikernel->wbuffer;
      for (j = nlines; j--;)
        {
        max = 0.0;
        for (i=kwidth; i--;)
          if ((val = *(pixin++))>max)
            max = val;
        *(pixout++) = max;
        }
      }
    }

/* Finally, fill the output pointer(s) */
  *outpix = ikernel->buffer[0];
  if (woutpix)
    *woutpix = wfield? ikernel->wbuffer[0]
			: (badpixflag? BIG :
				(PIXTYPE)(field->backsig*field->backsig));

  return RETURN_OK;
  }


/****** make_kernel **********************************************************
PROTO	void make_kernel(double pos, double *kernel, interpenum interptype)
PURPOSE	Conpute interpolation-kernel data
INPUT	Position,
	Pointer to the output kernel data,
	Interpolation method.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	07/01/2008
 ***/
void	make_kernel(double pos, double *kernel, interpenum interptype)
  {
   double	x, val, sinx1,sinx2,sinx3,cosx1;

  if (interptype == INTERP_NEARESTNEIGHBOUR)
    *kernel = 1;
  else if (interptype == INTERP_BILINEAR)
    {
    *(kernel++) = 1.0-pos;
    *kernel = pos;
    }
  else if (interptype == INTERP_LANCZOS2)
    {
    if (pos<1e-5 && pos>-1e5)
      {
      *(kernel++) = 0.0;
      *(kernel++) = 1.0;
      *(kernel++) = 0.0;
      *kernel = 0.0;
      }
    else
      {
      x = -PI/2.0*(pos+1.0);
#ifdef HAVE_SINCOS
      sincos(x, &sinx1, &cosx1);
#else
      sinx1 = sin(x);
      cosx1 = cos(x);
#endif
      val = (*(kernel++) = sinx1/(x*x));
      x += PI/2.0;
      val += (*(kernel++) = -cosx1/(x*x));
      x += PI/2.0;
      val += (*(kernel++) = -sinx1/(x*x));
      x += PI/2.0;
      val += (*kernel = cosx1/(x*x));
      val = 1.0/val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *kernel *= val;
      }
    }
  else if (interptype == INTERP_LANCZOS3)
    {
    if (pos<1e-5 && pos>-1e5)
      {
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *(kernel++) = 1.0;
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *kernel = 0.0;
      }
    else
      {
      x = -PI/3.0*(pos+2.0);
#ifdef HAVE_SINCOS
      sincos(x, &sinx1, &cosx1);
#else
      sinx1 = sin(x);
      cosx1 = cos(x);
#endif
      val = (*(kernel++) = sinx1/(x*x));
      x += PI/3.0;
      val += (*(kernel++) = (sinx2=-0.5*sinx1-0.866025403785*cosx1)
				/ (x*x));
      x += PI/3.0;
      val += (*(kernel++) = (sinx3=-0.5*sinx1+0.866025403785*cosx1)
				/(x*x));
      x += PI/3.0;
      val += (*(kernel++) = sinx1/(x*x));
      x += PI/3.0;
      val += (*(kernel++) = sinx2/(x*x));
      x += PI/3.0;
      val += (*kernel = sinx3/(x*x));
      val = 1.0/val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *kernel *= val;
      }
    }
  else if (interptype == INTERP_LANCZOS4)
    {
    if (pos<1e-5 && pos>-1e5)
      {
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *(kernel++) = 1.0;
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *kernel = 0.0;
      }
    else
      {
      x = -PI/4.0*(pos+3.0);
#ifdef HAVE_SINCOS
      sincos(x, &sinx1, &cosx1);
#else
      sinx1 = sin(x);
      cosx1 = cos(x);
#endif
      val = (*(kernel++) = sinx1/(x*x));
      x += PI/4.0;
      val +=(*(kernel++) = -(sinx2=0.707106781186*(sinx1+cosx1))
				/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = cosx1/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = -(sinx3=0.707106781186*(cosx1-sinx1))/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = -sinx1/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = sinx2/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = -cosx1/(x*x));
      x += PI/4.0;
      val += (*kernel = sinx3/(x*x));
      val = 1.0/val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *kernel *= val;
      }
    }
  else
    error(EXIT_FAILURE, "*Internal Error*: Unknown interpolation type in ",
		"make_kernel()");

  return;
  }


/****** init_ikernel **********************************************************
PROTO	ikernelstruct	*init_ikernel(interpenum *interptype, int naxis)
PURPOSE	Prepare interpolation operations.
INPUT	Interpolation type for each axis,
	Number of axes.
OUTPUT	Pointer to the newly created ikernel structure.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	16/02/2002
 ***/
ikernelstruct	*init_ikernel(interpenum *interptype, int naxis)
  {
   ikernelstruct	*ikernel;
   int			n;

  QMALLOC(ikernel, ikernelstruct, 1);
  ikernel->nlines = 1;
  for (n=0; n<naxis; n++)
    {
    ikernel->nlines*=(ikernel->width[n]=interp_kernwidth[(int)interptype[n]]);
    ikernel->interptype[n] = interptype[n];
    }
  ikernel->nlines /= ikernel->width[0];
  QMALLOC(ikernel->buffer, PIXTYPE, ikernel->nlines);
  QMALLOC(ikernel->wbuffer,PIXTYPE, ikernel->nlines);

  return ikernel;
  }


/****** free_ikernel **********************************************************
PROTO	void free_kernel(ikernelstruct *ikernel)
PURPOSE	Free buffers allocated for interpolation operations.
INPUT	Interpolation kernel.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	16/03/2002
 ***/
void	free_ikernel(ikernelstruct *ikernel)
  {
  free(ikernel->buffer);
  free(ikernel->wbuffer);
  free(ikernel);

  return;
  }

