/*
*				back.c
*
* Functions dealing with the image background.
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

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#ifdef HAVE_MATHIMF_H
#include <mathimf.h>
#else
#include <math.h>
#endif
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"fits/fitscat.h"
#include	"back.h"
#include	"coadd.h"
#include	"field.h"
#include	"misc.h"
#include	"prefs.h"
#include	"weight.h"

/******************************** make_back **********************************/
/*
Background maps are established from the images themselves; thus we need to
make at least one first pass through the data.
*/
void	make_back(fieldstruct *field, fieldstruct *wfield, int wscale_flag)

  {
   backstruct	*backmesh,*wbackmesh, *bm,*wbm;
   tabstruct	*tab, *wtab;
   PIXTYPE	*buf,*wbuf, *buft,*wbuft;
   size_t	bufsize, bufsize2,
		bufshift, size,meshsize,jumpsize;
   off_t	fcurpos,wfcurpos, wfcurpos2,fcurpos2;
   int		i,j,k,m,n, step, nlines,
		w,bw, bh, nx,ny,nb,
		lflag, nr;
   float	*ratio,*ratiop, *weight, *sigma,
		sratio;

/* If the weight-map is not an external one, no stats are needed for it */
  if (wfield && (wfield->flags&BACKRMS_FIELD))
    wfield= NULL;
  tab = field->tab;
  if (wfield)
    wtab = wfield->tab;
  else
    wtab = NULL;	/* to avoid gcc -Wall warnings */
  w = field->width;
  bw = field->backw;
  bh = field->backh;
  nx = field->nbackx;
  ny = field->nbacky;
  nb = field->nback;

  NFPRINTF(OUTPUT, "Setting up background maps ...");

/* Decide if it is worth displaying progress each 16 lines */

  lflag = (field->width*field->backh >= (size_t)65536);

/* Save current positions in files */

  wfcurpos = wfcurpos2 = 0;	/* to avoid gcc -Wall warnings */
  QFSEEK(tab->cat->file, tab->bodypos, SEEK_SET, field->filename);
  QFTELL(fcurpos, tab->cat->file, field->filename);
  if (wfield)
    {
    QFSEEK(wtab->cat->file, wtab->bodypos, SEEK_SET, wfield->filename);
    QFTELL(wfcurpos, wtab->cat->file, wfield->filename);
    }

/* Allocate a correct amount of memory to store pixels */

  bufsize = (size_t)w*bh;
  meshsize = bufsize;
  nlines = 0;
  if (bufsize > (size_t)BACK_BUFSIZE)
    {
    nlines = BACK_BUFSIZE/w;
    step = (field->backh-1)/nlines+1;
    bufsize = (size_t)(nlines = field->backh/step)*w;
    bufshift = (step/2)*(size_t)w;
    jumpsize = (step-1)*(size_t)w;
    }
  else
    bufshift = jumpsize = 0;	/* to avoid gcc -Wall warnings */

/* Allocate some memory */
  QMALLOC(backmesh, backstruct, nx);		/* background information */
  QMALLOC(buf, PIXTYPE, bufsize);		/* pixel buffer */
  free(field->back);
  QMALLOC(field->back, float, nb);		/* background map */
  free(field->backline);
  QMALLOC(field->backline, PIXTYPE, w);		/* current background line */
  free(field->sigma);
  QMALLOC(field->sigma, float, nb);		/* sigma map */
  if (wfield)
    {
    QMALLOC(wbackmesh, backstruct, nx);		/* background information */
    QMALLOC(wbuf, PIXTYPE, bufsize);		/* pixel buffer */
    free(wfield->back);
    QMALLOC(wfield->back, float, nb);		/* background map */
    free(wfield->sigma);
    QMALLOC(wfield->sigma, float, nb);		/* sigma map */
    wfield->sigfac = 1.0;
    set_weightconv(wfield);			/* Prepare conversion */
    }
  else
    {
    wbackmesh = NULL;
    wbuf = NULL;
    }

/* Loop over the data packets */
  for (j=0; j<ny; j++)
    {
    if (lflag && j)
      NPRINTF(OUTPUT,
	"\33[1M> Setting up background map at line:%7d / %-7d\n\33[1A",
	      j*bh, field->height);
    if (!nlines)
      {
/*---- The image is small enough so that we can make exhaustive stats */
      if (j == ny-1 && field->npix%bufsize)
        bufsize = field->npix%bufsize;
      read_body(tab, buf, bufsize);
      if (wfield)
        {
        read_body(wtab, wbuf, bufsize);
        weight_to_var(wbuf, bufsize);
        }
/*---- Build the histograms */
      backstat(backmesh, wbackmesh, buf, wbuf, bufsize,nx, w, bw,
	wfield?wfield->var_thresh:0.0);
      bm = backmesh;
      for (m=nx; m--; bm++)
        if (bm->mean <= -BIG)
          bm->histo=NULL;
        else
          QCALLOC(bm->histo, int, bm->nlevels);
      if (wfield)
        {
        wbm = wbackmesh;
        for (m=nx; m--; wbm++)
          if (wbm->mean <= -BIG)
            wbm->histo=NULL;
          else
            QCALLOC(wbm->histo, int, wbm->nlevels);
        }
      backhisto(backmesh, wbackmesh, buf, wbuf, bufsize,nx, w, bw,
	wfield?wfield->var_thresh:0.0);
      }
    else
      {
/*---- Image size too big, we have to skip a few data !*/
      QFTELL(fcurpos2, tab->cat->file, field->filename);
      if (wfield)
        QFTELL(wfcurpos2, wtab->cat->file, wfield->filename);
      if (j == ny-1 && (n=field->height%field->backh))
        {
        meshsize = n*(size_t)w;
        nlines = BACK_BUFSIZE/w;
        step = (n-1)/nlines+1;
        bufsize = (nlines = n/step)*(size_t)w;
        bufshift = (step/2)*(size_t)w;
        jumpsize = (step-1)*(size_t)w;
        free(buf);
        QMALLOC(buf, PIXTYPE, bufsize);		/* pixel buffer */
        if (wfield)
          {
          free(wbuf);
          QMALLOC(wbuf, PIXTYPE, bufsize);	/* pixel buffer */
          }
        }

/*---- Read and skip, read and skip, etc... */
      QFSEEK(tab->cat->file, bufshift*tab->bytepix, SEEK_CUR, field->filename);
      buft = buf;
      for (i=nlines; i--; buft += w)
        {
        read_body(tab, buft, w);
        if (i)
          QFSEEK(tab->cat->file, jumpsize*tab->bytepix, SEEK_CUR,
		field->filename);
        }

      if (wfield)
        {
/*------ Read and skip, read and skip, etc... now on the weight-map */
        QFSEEK(wtab->cat->file,bufshift*wtab->bytepix, SEEK_CUR,
		wfield->filename);
        wbuft = wbuf;
        for (i=nlines; i--; wbuft += w)
          {
          read_body(wtab, wbuft, w);
          weight_to_var(wbuft, w);
          if (i)
            QFSEEK(wtab->cat->file, jumpsize*wtab->bytepix, SEEK_CUR,
		wfield->filename);
          }
        }
      backstat(backmesh, wbackmesh, buf, wbuf, bufsize, nx, w, bw,
	wfield?wfield->var_thresh:0.0);
      QFSEEK(tab->cat->file, fcurpos2, SEEK_SET, field->filename);
      bm = backmesh;
      for (m=nx; m--; bm++)
        if (bm->mean <= -BIG)
          bm->histo=NULL;
        else
          QCALLOC(bm->histo, int, bm->nlevels);
      if (wfield)
        {
        QFSEEK(wtab->cat->file, wfcurpos2, SEEK_SET, wfield->filename);
        wbm = wbackmesh;
        for (m=nx; m--; wbm++)
          if (wbm->mean <= -BIG)
            wbm->histo=NULL;
          else
            QCALLOC(wbm->histo, int, wbm->nlevels);
        }
/*---- Build (progressively this time) the histograms */
      for(size=meshsize, bufsize2=bufsize; size>0; size -= bufsize2)
        {
        if (bufsize2>size)
          bufsize2 = size;
        read_body(tab, buf, bufsize2);
        if (wfield)
          {
          read_body(wtab, wbuf, bufsize2);
          weight_to_var(wbuf, bufsize2);
          }
        backhisto(backmesh, wbackmesh, buf, wbuf, bufsize2, nx, w, bw,
		wfield?wfield->var_thresh:0.0);
        }
      }

/*-- Compute background statistics from the histograms */
    bm = backmesh;
    for (m=0; m<nx; m++, bm++)
      {
      k = m+nx*j;
      backguess(bm, field->back+k, field->sigma+k);
      free(bm->histo);
      }
    if (wfield)
      {
      wbm = wbackmesh;
      for (m=0; m<nx; m++, wbm++)
        {
        k = m+nx*j;
        backguess(wbm, wfield->back+k, wfield->sigma+k);
        free(wbm->histo);
        }
      }
    }

/* Free memory */
  free(buf);
  free(backmesh);
  if (wfield)
    {
    free(wbackmesh);
    free(wbuf);
    }

/* Go back to the original position */
  QFSEEK(field->tab->cat->file, fcurpos, SEEK_SET, field->filename);
  if (wfield)
    QFSEEK(wfield->tab->cat->file, wfcurpos, SEEK_SET, wfield->filename);

/* Median-filter and check suitability of the background map */
  NFPRINTF(OUTPUT, "Filtering background map(s) ...");
  filter_back(field);
  if (wfield)
    filter_back(wfield);

/* Compute normalization for variance- or weight-maps*/
  if (wfield && wscale_flag && wfield->flags&(VAR_FIELD|WEIGHT_FIELD))
    {      
    nr = 0;
    QMALLOC(ratio, float, wfield->nback);
    ratiop = ratio;
    weight = wfield->back;
    sigma = field->sigma;
    for (i=wfield->nback; i--; sigma++)
      if ((sratio=*(weight++)) > 0.0
		&& (sratio = *sigma/sqrt(sratio)) > 0.0)
        {
        *(ratiop++) = sratio;
        nr++;
        }
    wfield->sigfac = (double)fqmedian(ratio, nr);
    for (i=0; i<nr && ratio[i]<=0.0; i++);
    if (i<nr)
      wfield->sigfac = (double)fqmedian(ratio+i, nr-i);
    else
      {
      warning("Null or negative global weighting factor:","defaulted to 1");
      wfield->sigfac = 1.0;
      } 

/*-- Update variance threshold */
    wfield->var_thresh *= (PIXTYPE)(wfield->sigfac*wfield->sigfac);

    free(ratio);
    }

/* Compute 2nd derivatives along the y-direction */
  NFPRINTF(OUTPUT, "Computing backgound d-map ...");
  free(field->dback);
  field->dback = make_backspline(field, field->back);
  NFPRINTF(OUTPUT, "Computing backgound-noise d-map ...");
  free(field->dsigma);
  field->dsigma = make_backspline(field, field->sigma);
/* If asked for, force the backmean parameter to the supplied value */
  if (field->back_type == BACK_ABSOLUTE)
    field->backmean = field->backdefault;

  field->fbackmean = field->backmean*field->fascale;
  field->fbacksig = field->backsig*field->fascale;

  return;
  }


/******************************** backstat **********************************/
/*
Compute robust statistical estimators in a row of meshes.
*/
void	backstat(backstruct *backmesh, backstruct *wbackmesh,
		PIXTYPE *buf, PIXTYPE *wbuf, size_t bufsize,
			int n, int w, int bw, PIXTYPE wthresh)

  {
   backstruct	*bm, *wbm;
   double	pix,wpix, sig, mean,wmean, sigma,wsigma, step;
   PIXTYPE	*buft,*wbuft, lcut,wlcut, hcut,whcut;
   int		m,h,x,y, npix,wnpix, offset, lastbite;

  h = bufsize/w;
  bm = backmesh;
  wbm = wbackmesh;
  offset = w - bw;
  step = sqrt(2/PI)*QUANTIF_NSIGMA/QUANTIF_AMIN;
  wmean = wsigma = wlcut = whcut = 0.0;	/* to avoid gcc -Wall warnings */
  for (m = n; m--; bm++,buf+=bw)
    {
    if (!m && (lastbite=w%bw))
      {
      bw = lastbite;
      offset = w-bw;
      }
    mean = sigma = 0.0;
    buft=buf;
/*-- We separate the weighted case at this level to avoid penalty in CPU */
    npix = 0;
    if (wbackmesh)
      {
      wmean = wsigma = 0.0;
      wbuft = wbuf;
      for (y=h; y--; buft+=offset,wbuft+=offset)
        for (x=bw; x--;)
          {
          pix = *(buft++);
          if ((wpix = *(wbuft++)) < wthresh && pix > -BIG)
            {
            wmean += wpix;
            wsigma += wpix*wpix;
            mean += pix;
            sigma += pix*pix;
            npix++;
            }
	  }
      }
    else
      {
      for (y=h; y--; buft+=offset)
        for (x=bw; x--;)
          if ((pix = *(buft++)) > -BIG)
	    {
            mean += pix;
            sigma += pix*pix;
            npix++;
            }
      }

/*-- If not enough valid pixels, discard this mesh */
    if ((float)npix < (float)(bw*h*BACK_MINGOODFRAC))
      {
      bm->mean = bm->sigma = -BIG;
      if (wbackmesh)
        {
        wbm->mean = wbm->sigma = -BIG;
        wbm++;
        wbuf += bw;
        }
      continue;
      }
    if (wbackmesh)
      {
      wmean /= (double)npix;
      wsigma = (sig = wsigma/npix - wmean*wmean)>0.0? sqrt(sig):0.0;
      wlcut = wbm->lcut = (PIXTYPE)(wmean - 2.0*wsigma);
      whcut = wbm->hcut = (PIXTYPE)(wmean + 2.0*wsigma);
      }
    mean /= (double)npix;
    sigma = (sig = sigma/npix - mean*mean)>0.0? sqrt(sig):0.0;
    lcut = bm->lcut = (PIXTYPE)(mean - 2.0*sigma);
    hcut = bm->hcut = (PIXTYPE)(mean + 2.0*sigma);
    mean = sigma = 0.0;
    npix = wnpix = 0;
    buft = buf;
    if (wbackmesh)
      {
      wmean = wsigma = 0.0;
      wbuft=wbuf;
      for (y=h; y--; buft+=offset, wbuft+=offset)
        for (x=bw; x--;)
          {
          pix = *(buft++);
          if ((wpix = *(wbuft++))<wthresh && pix<=hcut && pix>=lcut)
            {
            mean += pix;
            sigma += pix*pix;
            npix++;
            if (wpix<=whcut && wpix>=wlcut)
              {
              wmean += wpix;
              wsigma += wpix*wpix;
              wnpix++;
              }
            }
          }
      }
    else
      for (y=h; y--; buft+=offset)
        for (x=bw; x--;)
          {
          pix = *(buft++);
          if (pix<=hcut && pix>=lcut)
            {
            mean += pix;
            sigma += pix*pix;
            npix++;
            }
          }

    bm->npix = npix;
    mean /= (double)npix;
    sig = sigma/npix - mean*mean;
    sigma = sig>0.0 ? sqrt(sig):0.0;
    bm->mean = mean;
    bm->sigma = sigma;
    if ((bm->nlevels = (int)(step*npix+1)) > QUANTIF_NMAXLEVELS)
      bm->nlevels = QUANTIF_NMAXLEVELS;
    bm->qscale = sigma>0.0? 2*QUANTIF_NSIGMA*sigma/bm->nlevels : 1.0;
    bm->qzero = mean - QUANTIF_NSIGMA*sigma;
    if (wbackmesh)
      {
      wbm->npix = wnpix;
      wmean /= (double)wnpix;
      sig = wsigma/wnpix - wmean*wmean;
      wsigma = sig>0.0 ? sqrt(sig):0.0;
      wbm->mean = wmean;
      wbm->sigma = wsigma;
      if ((wbm->nlevels = (int)(step*wnpix+1)) > QUANTIF_NMAXLEVELS)
        wbm->nlevels = QUANTIF_NMAXLEVELS;
      wbm->qscale = wsigma>0.0? 2*QUANTIF_NSIGMA*wsigma/wbm->nlevels : 1.0;
      wbm->qzero = wmean - QUANTIF_NSIGMA*wsigma;
      wbm++;
      wbuf += bw;
      }
    }

  return;
  }


/******************************** backhisto *********************************/
/*
Compute robust statistical estimators in a row of meshes.
*/
void	backhisto(backstruct *backmesh, backstruct *wbackmesh,
		PIXTYPE *buf, PIXTYPE *wbuf, size_t bufsize,
			int n, int w, int bw, PIXTYPE wthresh)
  {
   backstruct	*bm,*wbm;
   PIXTYPE	*buft,*wbuft,
		pix;
   float	qscale,wqscale, cste,wcste, wpix;
   int		*histo,*whisto;
   int		h,m,x,y, nlevels,wnlevels, lastbite, offset, bin;

  h = bufsize/w;
  bm = backmesh;
  wbm = wbackmesh;
  offset = w - bw;
  for (m=0; m++<n; bm++ , buf+=bw)
    {
    if (m==n && (lastbite=w%bw))
      {
      bw = lastbite;
      offset = w-bw;
      }
/*-- Skip bad meshes */
    if (bm->mean <= -BIG)
      {
      if (wbackmesh)
        {
        wbm++;
        wbuf += bw;
        }
      continue;
      }
    nlevels = bm->nlevels;
    histo = bm->histo;
    qscale = bm->qscale;
    cste = 0.499999 - bm->qzero/qscale;
    buft = buf;
    if (wbackmesh)
      {
      wnlevels = wbm->nlevels;
      whisto = wbm->histo;
      wqscale = wbm->qscale;
      wcste = 0.499999 - wbm->qzero/wqscale;
      wbuft = wbuf;
      for (y=h; y--; buft+=offset, wbuft+=offset)
        for (x=bw; x--;)
          {
          bin = (int)((pix=*(buft++))/qscale + cste);
          if ((wpix = *(wbuft++))<wthresh && pix>-BIG && bin<nlevels && bin>=0)
            {
            (*(histo+bin))++;
            bin = (int)(wpix/wqscale + wcste);
            if (bin>=0 && bin<wnlevels)
              (*(whisto+bin))++;
            }
          }
      wbm++;
      wbuf += bw;
      }
    else
      for (y=h; y--; buft += offset)
        for (x=bw; x--;)
          {
          bin = (int)((pix=*(buft++))/qscale + cste);
          if (bin>=0 && bin<nlevels && pix>-BIG)
            (*(histo+bin))++;
          }
    }

  return;
  }

/******************************* backguess **********************************/
/*
Estimate the background from a histogram;
*/
float	backguess(backstruct *bkg, float *mean, float *sigma)

#define	EPS	(1e-4)	/* a small number */

  {
   int		*histo, *hilow, *hihigh, *histot;
   unsigned int lowsum, highsum, sum;
   double	ftemp, mea, sig, sig1, med, dpix;
   int		i, n, lcut,hcut, nlevelsm1, pix;

/* Leave here if the mesh is already classified as `bad' */
  if (bkg->mean<=-BIG)
    {
    *mean = *sigma = -BIG;
    return -BIG;
    }

  histo = bkg->histo;
  hcut = nlevelsm1 = bkg->nlevels-1;
  lcut = 0;

  sig = 10.0*nlevelsm1;
  sig1 = 1.0;
  mea = med = bkg->mean;
  for (n=100; n-- && (sig>=0.1) && (fabs(sig/sig1-1.0)>EPS);)
    {
    sig1 = sig;
    sum = mea = sig = 0.0;
    lowsum = highsum = 0;
    histot = hilow = histo+lcut;
    hihigh = histo+hcut;
    for (i=lcut; i<=hcut; i++)
      {
      if (lowsum<highsum)
        lowsum += *(hilow++);
      else
        highsum +=  *(hihigh--);
      sum += (pix = *(histot++));
      mea += (dpix = (double)pix*i);
      sig += dpix*i;
      }

    med = hihigh>=histo?
	((hihigh-histo)+0.5+((double)highsum-lowsum)/(2.0*(*hilow>*hihigh?
                                                *hilow:*hihigh)))
       : 0.0;

    if (sum)
      {
      mea /= (double)sum;
      sig = sig/sum - mea*mea;
      }
    sig = sig>0.0?sqrt(sig):0.0;
    lcut = (ftemp=med-3.0*sig)>0.0 ?(int)(ftemp>0.0?ftemp+0.5:ftemp-0.5):0;
    hcut = (ftemp=med+3.0*sig)<nlevelsm1 ?(int)(ftemp>0.0?ftemp+0.5:ftemp-0.5)
								: nlevelsm1;
    }

  *mean = fabs(sig)>0.0? (fabs(bkg->sigma/(sig*bkg->qscale)-1) < 0.0 ?
			    bkg->qzero+mea*bkg->qscale
			    :(fabs((mea-med)/sig)< 0.3 ?
			      bkg->qzero+(2.5*med-1.5*mea)*bkg->qscale
			     :bkg->qzero+med*bkg->qscale))
                       :bkg->qzero+mea*bkg->qscale;

  *sigma = sig*bkg->qscale;


  return *mean;
  }


/******************************* filter_back *********************************/
/*
Median filtering of the background map to remove the contribution from bright
sources.
*/
void	filter_back(fieldstruct *field)

  {
   float	*back,*sigma, *back2,*sigma2, *bmask,*smask, *sigmat,
		d2,d2min, fthresh, med, val,sval;
   int		i,j,px,py, np, nx,ny, npx,npx2, npy,npy2, dpx,dpy, x,y, nmin;

  fthresh = prefs.back_fthresh;
  nx = field->nbackx;
  ny = field->nbacky;
  np = field->nback;
  npx = field->nbackfx/2;
  npy = field->nbackfy/2;
  npy *= nx;

  QMALLOC(bmask, float, (2*npx+1)*(2*npy+1));
  QMALLOC(smask, float, (2*npx+1)*(2*npy+1));
  QMALLOC(back2, float, np);
  QMALLOC(sigma2, float, np);

  back = field->back;
  sigma = field->sigma;
  val = sval = 0.0;			/* to avoid gcc -Wall warnings */

/* Look for `bad' meshes and interpolate them if necessary */
  for (i=0,py=0; py<ny; py++)
    for (px=0; px<nx; px++,i++)
      if ((back2[i]=back[i])<=-BIG)
        {
/*------ Seek the closest valid mesh */
        d2min = BIG;
        nmin = 0;
        for (j=0,y=0; y<ny; y++)
          for (x=0; x<nx; x++,j++)
            if (back[j]>-BIG)
              {
              d2 = (float)(x-px)*(x-px)+(y-py)*(y-py);
              if (d2<d2min)
                {
                val = back[j];
                sval = sigma[j];
                nmin = 1;
                d2min = d2;
                }
              else if (d2==d2min)
                {
                val += back[j];
                sval += sigma[j];
                nmin++;
                }
              }
        back2[i] = nmin? val/nmin: 0.0;
        sigma[i] = nmin? sval/nmin: 1.0;
        }
  memcpy(back, back2, (size_t)np*sizeof(float));

/* Do the actual filtering */
  for (py=0; py<np; py+=nx)
    {
    npy2 = np - py - nx;
    if (npy2>npy)
      npy2 = npy;
    if (npy2>py)
      npy2 = py;
    for (px=0; px<nx; px++)
      {
      npx2 = nx - px - 1;
      if (npx2>npx)
        npx2 = npx;
      if (npx2>px)
        npx2 = px;
      i=0;
      for (dpy = -npy2; dpy<=npy2; dpy+=nx)
        {
        y = py+dpy;
        for (dpx = -npx2; dpx <= npx2; dpx++)
          {
          x = px+dpx;
          bmask[i] = back[x+y];
          smask[i++] = sigma[x+y];
          }
        }
      if (fabs((med=fqmedian(bmask, i))-back[px+py])>=fthresh)
        {
        back2[px+py] = med;
        sigma2[px+py] = fqmedian(smask, i);
        }
      else
        {
        back2[px+py] = back[px+py];
        sigma2[px+py] = sigma[px+py];
        }
      }
    }

  free(bmask);
  free(smask);
  memcpy(back, back2, np*sizeof(float));
  field->backmean = (double)fqmedian(back2, np);
  free(back2);
  memcpy(sigma, sigma2, np*sizeof(float));
  field->backsig = (double)fqmedian(sigma2, np);

  if (field->backsig<=0.0)
    {
    sigmat = sigma2+np;
    for (i=np; i-- && *(--sigmat)>0.0;);
    if (i>=0 && i<(np-1))
      field->backsig = fqmedian(sigmat+1, np-1-i);
    else
      {
      if (field->flags&(DETECT_FIELD|MEASURE_FIELD))
        warning("Image contains mainly constant data; ",
		"I'll try to cope with that...");
      field->backsig = 1.0;
      }
    }

  free(sigma2);


  return;
  }


/******************************* make_backspline *****************************/
/*
Pre-compute 2nd derivatives along the y direction at background nodes.
*/
float *make_backspline(fieldstruct *field, float *map)

  {
   int		x,y, nbx,nby,nbym1;
   float	*dmap,*dmapt,*mapt, *u, temp;

  nbx = field->nbackx;
  nby = field->nbacky;
  nbym1 = nby - 1;
  QMALLOC(dmap, float, field->nback);
  for (x=0; x<nbx; x++)
    {
    mapt = map+x;
    dmapt = dmap+x;
    if (nby>1)
      {
      QMALLOC(u, float, nbym1);	/* temporary array */
      *dmapt = *u = 0.0;	/* "natural" lower boundary condition */
      mapt += nbx;
      for (y=1; y<nbym1; y++, mapt+=nbx)
        {
        temp = -1/(*dmapt+4);
        *(dmapt += nbx) = temp;
        temp *= *(u++) - 6*(*(mapt+nbx)+*(mapt-nbx)-2**mapt);
        *u = temp;
        }
      *(dmapt+=nbx) = 0.0;	/* "natural" upper boundary condition */
      for (y=nby-2; y--;)
        {
        temp = *dmapt;
        dmapt -= nbx;
        *dmapt = (*dmapt*temp+*(u--))/6.0;
        }
      free(u);
      }
    else
      *dmapt = 0.0;
    }

  return dmap;
  }


/******************************** backline **********************************/
/*
Interpolate background at line y (bicubic spline interpolation between
background map vertices).
*/
void	backline(fieldstruct *field, int y, PIXTYPE *line)

  {
   int		i,j,x,yl, nbx,nbxm1,nby, nx,width, ystep, changepoint;
   float	dx,dx0,dy,dy3, cdx,cdy,cdy3, temp, xstep,
		*node,*nodep,*dnode, *blo,*bhi,*dblo,*dbhi, *u;
   PIXTYPE	bval;

  width = field->width;

  if (field->back_type==BACK_ABSOLUTE)
    {
/*-- In absolute background mode, just subtract a cste */
    bval = (PIXTYPE)field->backmean;
    for (i=width; i--;)
      *(line++) = bval;
    return;
    }

  nbx = field->nbackx;
  nbxm1 = nbx - 1;
  nby = field->nbacky;
  if (nby > 1)
    {
    dy = (float)y/field->backh - 0.5;	
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
/*-- Interpolation along y for each node */
    cdy = 1 - dy;
    dy3 = (dy*dy*dy-dy);
    cdy3 = (cdy*cdy*cdy-cdy);
    ystep = nbx*yl;
    blo = field->back + ystep;
    bhi = blo + nbx;
    dblo = field->dback + ystep;
    dbhi = dblo + nbx;
    QMALLOC(node, float, nbx);	/* Interpolated background */
    nodep = node;
    for (x=nbx; x--;)
      *(nodep++) = cdy**(blo++) + dy**(bhi++) + cdy3**(dblo++) + dy3**(dbhi++);

/*-- Computation of 2nd derivatives along x */
    QMALLOC(dnode, float, nbx);	/* 2nd derivative along x */
    if (nbx>1)
      {
      QMALLOC(u, float, nbxm1);	/* temporary array */
      *dnode = *u = 0.0;	/* "natural" lower boundary condition */
      nodep = node+1;
      for (x=nbxm1; --x; nodep++)
        {
        temp = -1/(*(dnode++)+4);
        *dnode = temp;
        temp *= *(u++) - 6*(*(nodep+1)+*(nodep-1)-2**nodep);
        *u = temp;
        }
      *(++dnode) = 0.0;	/* "natural" upper boundary condition */
      for (x=nbx-2; x--;)
        {
        temp = *(dnode--);
        *dnode = (*dnode*temp+*(u--))/6.0;
        }
      free(u);
      dnode--;
      }
    }
  else
    {
/*-- No interpolation and no new 2nd derivatives needed along y */
    node = field->back;
    dnode = field->dback;
    }

/*-- Interpolation along x */
  if (nbx>1)
    {
    nx = field->backw;
    xstep = 1.0/nx;
    changepoint = nx/2;
    dx  = (xstep - 1)/2;	/* dx of the first pixel in the row */
    dx0 = ((nx+1)%2)*xstep/2;	/* dx of the 1st pixel right to a bkgnd node */
    blo = node;
    bhi = node + 1;
    dblo = dnode;
    dbhi = dnode + 1;
    for (x=i=0,j=width; j--; i++, dx += xstep)
      {
      if (i==changepoint && x>0 && x<nbxm1)
        {
        blo++;
        bhi++;
        dblo++;
        dbhi++;
        dx = dx0;
        }
      cdx = 1 - dx;
      *(line++) = (PIXTYPE)(cdx*(*blo+(cdx*cdx-1)**dblo)
			+ dx*(*bhi+(dx*dx-1)**dbhi));
      if (i==nx)
        {
        x++;
        i = 0;
        }
      }
    }
  else
    for (j=width; j--;)
      *(line++) = (PIXTYPE)*node;

  if (nby>1)
    {
    free(node);
    free(dnode);
    }

  return;
  }


/******************************* backrmsline ********************************
PROTO   void backrmsline(fieldstruct *field, int y, PIXTYPE *line)
PURPOSE Bicubic-spline interpolation of the background noise along the current
        scanline (y).
INPUT   Measurement or detection field pointer,
        Current line position. 
        Where to put the data. 
OUTPUT  -.
NOTES   Most of the code is a copy of subbackline(), for optimization reasons.
AUTHOR  E. Bertin (IAP & Leiden & ESO)
VERSION 02/02/98
 ***/
void	backrmsline(fieldstruct *field, int y, PIXTYPE *line)

  {
   int		i,j,x,yl, nbx,nbxm1,nby, nx,width, ystep, changepoint;
   float	dx,dx0,dy,dy3, cdx,cdy,cdy3, temp, xstep,
		*node,*nodep,*dnode, *blo,*bhi,*dblo,*dbhi, *u;

  nbx = field->nbackx;
  nbxm1 = nbx - 1;
  nby = field->nbacky;
  if (nby > 1)
    {
    dy = (float)y/field->backh - 0.5;
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
/*-- Interpolation along y for each node */
    cdy = 1 - dy;
    dy3 = (dy*dy*dy-dy);
    cdy3 = (cdy*cdy*cdy-cdy);
    ystep = nbx*yl;
    blo = field->sigma + ystep;
    bhi = blo + nbx;
    dblo = field->dsigma + ystep;
    dbhi = dblo + nbx;
    QMALLOC(node, float, nbx);	/* Interpolated background */
    nodep = node;
    for (x=nbx; x--;)
      *(nodep++) = cdy**(blo++) + dy**(bhi++) + cdy3**(dblo++) + dy3**(dbhi++);

/*-- Computation of 2nd derivatives along x */
    QMALLOC(dnode, float, nbx);	/* 2nd derivative along x */
    if (nbx>1)
      {
      QMALLOC(u, float, nbxm1);	/* temporary array */
      *dnode = *u = 0.0;	/* "natural" lower boundary condition */
      nodep = node+1;
      for (x=nbxm1; --x; nodep++)
        {
        temp = -1/(*(dnode++)+4);
        *dnode = temp;
        temp *= *(u++) - 6*(*(nodep+1)+*(nodep-1)-2**nodep);
        *u = temp;
        }
      *(++dnode) = 0.0;	/* "natural" upper boundary condition */
      for (x=nbx-2; x--;)
        {
        temp = *(dnode--);
        *dnode = (*dnode*temp+*(u--))/6.0;
        }
      free(u);
      dnode--;
      }
    }
  else
    {
/*-- No interpolation and no new 2nd derivatives needed along y */
    node = field->sigma;
    dnode = field->dsigma;
    }

/*-- Interpolation along x */
  width = field->width;
  if (nbx>1)
    {
    nx = field->backw;
    xstep = 1.0/nx;
    changepoint = nx/2;
    dx  = (xstep - 1)/2;	/* dx of the first pixel in the row */
    dx0 = ((nx+1)%2)*xstep/2;	/* dx of the 1st pixel right to a bkgnd node */
    blo = node;
    bhi = node + 1;
    dblo = dnode;
    dbhi = dnode + 1;
    for (x=i=0,j=width; j--; i++, dx += xstep)
      {
      if (i==changepoint && x>0 && x<nbxm1)
        {
        blo++;
        bhi++;
        dblo++;
        dbhi++;
        dx = dx0;
        }
      cdx = 1 - dx;
      *(line++) = (PIXTYPE)(cdx*(*blo+(cdx*cdx-1)**dblo)
			+ dx*(*bhi+(dx*dx-1)**dbhi));
      if (i==nx)
        {
        x++;
        i = 0;
        }
      }
    }
  else
    for (j=width; j--;)
      *(line++) = (PIXTYPE)*node;

  if (nby>1)
    {
    free(node);
    free(dnode);
    }

  return;
  }


/********************************* end_back **********************************/
/*
Terminate background procedures (mainly freeing memory).
*/
void	end_back(fieldstruct *field)

  {
  free(field->back);
  free(field->dback);
  free(field->sigma);
  free(field->dsigma);
  free(field->backline);

  return;
  }

