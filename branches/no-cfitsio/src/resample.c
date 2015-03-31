/*
*				resample.c
*
* Manage high-level resampling.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SWarp
*
*	Copyright:		(C) 2000-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
#include "data.h"
#include "field.h"
#include "header.h"
#include "interpolate.h"
#include "prefs.h"
#include "projapprox.h"
#include "resample.h"
#ifdef USE_THREADS
#include "threads.h"
#endif
#include "weight.h"
#include "wcs/wcs.h"

#ifdef USE_THREADS
 pthread_t		*thread;
 pthread_mutex_t	linemutex;
 pthread_cond_t		*linecond;
 tabstruct		*ftab, *wtab;
 int			*queue,*proc,*writeflag,
			absline,procline,writeline;
#else
 int			nproc;
#endif

 fieldstruct		*infield,*inwfield, *field, *wfield;
 ikernelstruct		**ikernel;
 wcsstruct		**wcsinp, **wcsoutp;
 projappstruct		*projapp;
 double			rawmin[NAXIS], rawmax[NAXIS],rawpos0[NAXIS],
			stepover[NAXIS],
			**rawposp, **rawbuf, **oversampbuf,**oversampwbuf,
			**rawbufarea,
			ascale;
 PIXTYPE		**outbuf,**outwbuf;
 FLAGTYPE		**outibuf,**outwibuf, **oversampibuf,**oversampwibuf;
 int			**oversampnbuf, *oversamp,
			noversamp, oversampflag, width, height, naxis, nlines,
			approxflag, dispstep, iflag;
 char			padbuf[FBSIZE];

/*------------------------------ function -----------------------------------*/
#ifdef USE_THREADS
static int		pthread_nextline(int l);
static void		*pthread_warp_lines(void *arg);
#endif
static void		warp_line(int p);


/****** resample_field *******************************************************
PROTO	void resample_field(fieldstruct **pinfield, fieldstruct **pinwfield,
			fieldstruct *outfield, fieldstruct *outwfield,
			interpenum *interptype)
PURPOSE	Resample an image.
INPUT	Input pointer to field structure pointer,
	Input pointer weight-map field structure pointer,
	Output total field structure pointer,
	Output total weight-map field structure pointer,
	Interpolation type.
	Output
OUTPUT	-.
NOTES	The structure pointers pointed by pinfield and and pinwfield are
	updated and point to the resampled fields on output.
AUTHOR	E. Bertin (IAP)
VERSION	19/07/2013
 ***/
void	resample_field(fieldstruct **pinfield, fieldstruct **pinwfield,
		fieldstruct *outfield, fieldstruct *outwfield,
		interpenum *interptype)
  {
#ifdef USE_THREADS
   static pthread_attr_t	pthread_attr;
   int				p;
#else
   int				y;
#endif
   wcsstruct		wcsloc,
			*wcs;
   static char		filename[MAXCHAR],filename2[MAXCHAR],
			resampext1[MAXCHAR], resampext2[MAXCHAR];
   char			*pstr;
   double		ascale1, projerr;
   int			d,size, l;

  infield = *pinfield;
  inwfield = *pinwfield;

/* Create new file name */
  strcpy(filename2, infield->rfilename);
  if ((pstr=strrchr(filename2, '.')))
    *pstr = '\0';
  if (infield->version>1)
    sprintf(pstr, "_v%d", infield->version);
  if (infield->frameno)
    sprintf(filename, "%s/%s.%04d%s", prefs.resampdir_name, filename2,
				infield->frameno, prefs.resamp_suffix);
  else
    sprintf(filename, "%s/%s%s", prefs.resampdir_name, filename2,
				prefs.resamp_suffix);

/* We make a copy of the output field */
  field = inherit_field(filename, outfield, FIELD_WRITE);
  field->backmean = infield->backmean;
  field->fbackmean = infield->fbackmean;
  field->backsig = infield->backsig;
  field->fbacksig = infield->fbacksig;
  field->gain = infield->gain;
  field->fgain = infield->fgain;
  field->saturation = infield->saturation;
  field->fsaturation = infield->fsaturation;
  field->exptime = infield->exptime;
  field->fieldno = infield->fieldno;
  field->fascale = infield->fascale;
  field->fscale = infield->fscale;
  field->headflag = infield->headflag;
  strcpy(field->ident, infield->ident);
  iflag = (outfield->bitpix>0);

/* Now modify some characteristics of the output file */
/* We use a small "dirty" trick to frame the output */
  wcsloc = *(outfield->wcs);	/* Copy the CONTENT only */
  wcs = infield->wcs;
  wcsloc.obsdate = wcs->obsdate;
  wcsloc.epoch = wcs->epoch;
  naxis = wcsloc.naxis;
  for (d=0; d<naxis;d++)
    {
    wcsloc.outmin[d] = wcs->outmin[d] < 1? 1:wcs->outmin[d];
    wcsloc.outmax[d] = wcs->outmax[d] > wcsloc.naxisn[d]?
                        (double)wcsloc.naxisn[d] : wcs->outmax[d];
    wcsloc.crpix[d] -= (wcsloc.outmin[d]-1);
    wcsloc.naxisn[d] = (wcsloc.outmax[d]-wcsloc.outmin[d]+1);
    if (wcsloc.naxisn[d]<=0)
/*-- Image out of output bounds: create a dummy one with only 1 pixel */
      {
      wcsloc.naxisn[d]=1;
      wcsloc.outmax[d] = wcsloc.outmin[d];
      }
    }

  wcs = field->wcs = copy_wcs(&wcsloc);
  write_wcs(field->tab, wcs);

/* Update field characteristics */
  width = field->width = field->tab->naxisn[0];
  field->height = 1;
  for (d=1; d<naxis; d++)
    field->height *= field->tab->naxisn[d];
  height = field->height;
  field->npix = field->width*field->height;

/* Add relevant information to FITS headers */
  writefitsinfo_field(field, infield);

/* Write image header */
  if (open_cat(field->cat, WRITE_ONLY) != RETURN_OK)
      error(EXIT_FAILURE, "*Error*: cannot open for writing ", filename);
  if (prefs.removetmp_flag && prefs.combine_flag)
    add_cleanupfilename(filename);
  QFTELL(field->tab->headpos, field->cat->file, filename);
  QFWRITE(field->tab->headbuf, field->tab->headnblock*FBSIZE,
	field->cat->file, filename);
  QFTELL(field->tab->bodypos, field->cat->file, filename);

/* Get resampled suffix extension if available */
  strcpy(resampext1, prefs.resamp_suffix);
  *resampext2 = '\0';
  if ((pstr=strrchr(resampext1, '.')))
    {
    strcpy(resampext2, pstr);
    *pstr = '\0';
    }

/* Now go on with output weight-map */
  if (infield->frameno)
    sprintf(filename, "%s/%s.%04d%s.weight%s", prefs.resampdir_name, filename2,
					infield->frameno,resampext1,
					resampext2);
  else
    sprintf(filename, "%s/%s%s.weight%s", prefs.resampdir_name, filename2,
					resampext1, resampext2);

  wfield = init_weight(filename, field);
  if (inwfield)
    set_weightconv(inwfield);
/* Add relevant information to FITS headers */
  writefitsinfo_field(wfield, inwfield? inwfield : infield);
  wfield->cflags = 0;
  if (inwfield)
    {
    wfield->sigfac = inwfield->sigfac;
    wfield->weight_thresh = inwfield->weight_thresh;
    }

/* Write image header */
  if (open_cat(wfield->cat, WRITE_ONLY) != RETURN_OK)
    error(EXIT_FAILURE, "*Error*: cannot open for writing ", filename);
  if (prefs.removetmp_flag && prefs.combine_flag)
    add_cleanupfilename(filename);
  QFTELL(wfield->tab->headpos, wfield->cat->file, filename);
  QFWRITE(wfield->tab->headbuf, wfield->tab->headnblock*FBSIZE,
	wfield->cat->file, filename);
  QFTELL(wfield->tab->bodypos, wfield->cat->file, filename);

/* Prepare oversampling stuff */
  ascale = 1.0;
  noversamp = 1;
  oversamp = prefs.oversamp;
  for (d=0; d<naxis; d++)
    {
    ascale1 = wcs->wcsscale[d]/infield->wcs->wcsscale[d];
    ascale *= ascale1;
    if (!oversamp[d])
/*-- Automatic oversampling mode */
      oversamp[d] = ascale1>1.0 ? (int)(ascale1+0.5) : 1;
    noversamp *= oversamp[d];
    }
  oversampflag = (noversamp>1);
  if (ascale < 1.0)
    ascale = 1.0;

/* Turn approximation on or off */
  approxflag = ((projerr = prefs.proj_err[infield->fieldno] > 0.0)
    && (projapp = projapp_init(infield->wcs, field->wcs, projerr,
	prefs.fscalastro_type==FSCALASTRO_VARIABLE)));

#ifdef USE_THREADS
/* Set up multi-threading stuff */
/* Number of active threads */
  nproc = prefs.nthreads;
/* Twice more lines than threads to deal with desynchronization issues */
  nlines = 2*nproc;
  QMALLOC(linecond, pthread_cond_t, nlines);
  for (l=0; l<nlines; l++)
    {
    QPTHREAD_COND_INIT(&linecond[l], NULL);
    }
  QPTHREAD_MUTEX_INIT(&linemutex, NULL);
  QPTHREAD_ATTR_INIT(&pthread_attr);
  QPTHREAD_ATTR_SETDETACHSTATE(&pthread_attr, PTHREAD_CREATE_JOINABLE);
#else
  nproc = nlines = 1;
#endif

/* Initialize the astrometric vector */
  for (d=0; d<naxis;d++)
    {
    stepover[d]=1.0/oversamp[d];
    rawmin[d] = rawpos0[d] = 0.5 + 0.5*stepover[d];
    rawmax[d] = (double)wcs->naxisn[d];
    }

/*-- Allocate memory for buffer pointers */
  QMALLOC(rawposp, double *, nlines);
  if (iflag)
    {
    QMALLOC(outibuf, FLAGTYPE *, nlines);
    QMALLOC(outwibuf, FLAGTYPE *, nlines);
    }
  else
    {
    QMALLOC(outbuf, PIXTYPE *, nlines)
    QMALLOC(outwbuf, PIXTYPE *, nlines);
    }
  QMALLOC(rawbuf, double *, nlines);
  QCALLOC(rawbufarea, double *, nlines);
  QMALLOC(ikernel, ikernelstruct *, nlines);
  QMALLOC(wcsinp, wcsstruct *, nlines);
  QMALLOC(wcsoutp, wcsstruct *, nlines);
  if (oversampflag)
    {
    if (iflag)
      {
      QMALLOC(oversampibuf, FLAGTYPE *, nlines)
      QMALLOC(oversampwibuf, FLAGTYPE *, nlines);
      }
    else
      {
      QMALLOC(oversampbuf, double *, nlines)
      QMALLOC(oversampwbuf, double *, nlines);
      }
    QMALLOC(oversampnbuf, int *, nlines);
    }

  for (l=0; l<nlines; l++)
    {
/*-- Allocate memory for the output buffers that contain the final data in */
/*-- internal format (PIXTYPE) */
    if (iflag)
      {
      QMALLOC(outibuf[l], FLAGTYPE, width)
      QMALLOC(outwibuf[l], FLAGTYPE, width);
      }
    else
      {
      QMALLOC(outbuf[l], PIXTYPE, width)
      QMALLOC(outwbuf[l], PIXTYPE, width);
      }
    if (oversampflag)
/*---- Provide memory space for integrating the content of oversampled pixels*/
      {
      if (iflag)
        {
        QMALLOC(oversampibuf[l], FLAGTYPE, width)
        QMALLOC(oversampwibuf[l], FLAGTYPE, width);
        }
      else
        {
        QMALLOC(oversampbuf[l], double, width)
        QMALLOC(oversampwbuf[l], double, width);
        }
      QMALLOC(oversampnbuf[l], int, width);
      }
/*-- Provide memory space for the current astrometric line */
    QMALLOC(rawbuf[l], double, naxis*width);
    if (prefs.fscalastro_type==FSCALASTRO_VARIABLE)
      QMALLOC(rawbufarea[l], double, width)
    QMALLOC(rawposp[l], double, naxis);
/*-- Initialize interpolation kernel */
    ikernel[l] = init_ikernel(interptype, naxis);
/*-- Make copies of the WCS structure (the WCS library is not reentrant) */
    wcsinp[l] = copy_wcs(infield->wcs);
    wcsoutp[l] = copy_wcs(field->wcs);
    }

/* Compute reasonable line display-step */
  dispstep = (int)(nproc*50000.0/noversamp/width);
  if (!dispstep)
    dispstep = 1;

/* Start threads! */
#ifdef USE_THREADS
  QCALLOC(writeflag, int, nlines);
  QCALLOC(queue, int, nlines);
  QMALLOC(proc, int, nproc);
  QMALLOC(thread, pthread_t, nproc);
  ftab = field->tab;
  wtab = wfield->tab;
  writeline = absline = procline = 0;
  for (p=0; p<nproc; p++)
    {
    proc[p] = p;
    QPTHREAD_CREATE(&thread[p], &pthread_attr, &pthread_warp_lines, &proc[p]);
    }
#else
/* The old single-threaded way */
  for (y=0; y<height; y++)
    {
    for (d=1; d<naxis; d++)
        rawposp[0][d] = rawpos0[d];
/*---- Update coordinate vector */
    for (d=1; d<naxis; d++)
      if ((rawpos0[d]+=1.0)<= rawmax[d])
        break;
      else
        rawpos0[d] = rawmin[d];
    if (!(y%dispstep))
      NPRINTF(OUTPUT, "\33[1M> Resampling line:%7d / %-7d\n\33[1A",
	y, height);
    warp_line(0);
    if (iflag)
      {
      write_ibody(field->tab, outibuf[0], width);
      write_ibody(wfield->tab, outwibuf[0], width);
      }
    else
      {
      write_body(field->tab, outbuf[0], width);
      write_body(wfield->tab, outwbuf[0], width);
      }
    }
#endif

#ifdef USE_THREADS
  for (p=0; p<nproc; p++)
    QPTHREAD_JOIN(thread[p], NULL);
  QPTHREAD_MUTEX_DESTROY(&linemutex);
  for (l=0; l<nlines; l++)
    {
    QPTHREAD_COND_DESTROY(&linecond[l]);
    }
  QPTHREAD_ATTR_DESTROY(&pthread_attr);
  free(linecond);
  free(writeflag);
  free(queue);
  free(proc);
  free(thread);
#endif

/* FITS padding*/
  size = PADEXTRA(field->tab->tabsize);
  if (size)
    QFWRITE(padbuf, (size_t)size, field->cat->file, field->filename);
  size = PADEXTRA(wfield->tab->tabsize);
  if (size)
    QFWRITE(padbuf, (size_t)size, wfield->cat->file, wfield->filename);

/* Close files */
  close_cat(field->cat);
  close_cat(wfield->cat);

/* Return back new field pointers */
  *pinfield = field;
  *pinwfield = wfield;

/* Free memory */
  for (l=0; l<nlines; l++)
    {
    free(rawposp[l]);
    free(iflag? (void *)outibuf[l] : (void *)outbuf[l]);
    free(iflag? (void *)outwibuf[l]: (void *)outwbuf[l]);
    if (oversampflag)
      {
      free(iflag? (void *)oversampibuf[l] : (void *)oversampbuf[l]);
      free(iflag? (void *)oversampwibuf[l]: (void *)oversampwbuf[l]);
      free(oversampnbuf[l]);
      }
    free(rawbuf[l]);
    free(rawbufarea[l]);
/*-- Free interpolation kernel */
    free_ikernel(ikernel[l]);
    end_wcs(wcsinp[l]);
    end_wcs(wcsoutp[l]);
    }
  free(rawposp);
  free(iflag? (void *)outibuf : (void *)outbuf);
  free(iflag? (void *)outwibuf: (void *)outwbuf);
  if (oversampflag)
    {
    free(iflag? (void *)oversampibuf : (void *)oversampbuf);
    free(iflag? (void *)oversampwibuf: (void *)oversampwbuf);
    free(oversampnbuf);
    }
  free(rawbuf);
  free(rawbufarea);
  free(ikernel);
  free(wcsinp);
  free(wcsoutp);

  if (approxflag)
    projapp_end(projapp);

  end_field(infield);
  if (inwfield)
    end_field(inwfield);

  return;
  }

#ifdef USE_THREADS


/****** pthread_warp_lines ****************************************************
PROTO	void *pthread_warp_lines(void *arg)
PURPOSE	thread that takes care of resampling image lines
INPUT	Pointer to the thread number.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	20/10/2003
 ***/
void	*pthread_warp_lines(void *arg)
  {
   int	l;

  l = -1;
  while ((l=pthread_nextline(l))!= -1)
    warp_line(l);
  
  pthread_exit(NULL);

  return (void *)NULL;
  }


/****** pthread_nextline ******************************************************
PROTO	int pthread_nextline(int l)
PURPOSE	Return the next available line to be resampled.
INPUT	-.
OUTPUT	Next available line index.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	03/02/2012
 ***/
int	pthread_nextline(int l)
  {
   double	rawpos[NAXIS];
   int		d, q, y;

  QPTHREAD_MUTEX_LOCK(&linemutex);
/* The newly processed line is ready to be written to disk */
  if (l>=0)
    writeflag[l] = 2;
/* If we just finished the "right" line, write it to disk! */
  if (l == writeline)
    {
    while (writeflag[writeline]==2)
      {
      if (iflag)
        {
        write_ibody(ftab, outibuf[writeline], width);
        write_ibody(wtab, outwibuf[writeline], width);
        }
      else
        {
        write_body(ftab, outbuf[writeline], width);
        write_body(wtab, outwbuf[writeline], width);
        }
      writeflag[writeline] = 0;
      QPTHREAD_COND_BROADCAST(&linecond[writeline]);
      writeline = (writeline+1)%nlines;
      }
    }
/* If no more line to process, return a "-1" (meaning exit thread) */
  if ((y=absline++) >= height)
    l=-1;
  else
    {
    l = procline;
    procline = (procline+1)%nlines;
/*------ Update coordinate vector */
    for (d=1; d<naxis; d++)
      rawpos[d] = rawpos0[d];
    for (d=1; d<naxis; d++)
      if ((rawpos0[d]+=1.0)<= rawmax[d])
        break;
      else
        rawpos0[d] = rawmin[d];
/*-- If the next available buffer has not been flushed yet, wait */
    q=++queue[l];
    while (writeflag[l] || --q)
      QPTHREAD_COND_WAIT(&linecond[l], &linemutex);
/*-- Set content */
    if (queue[l])
      queue[l]--;
    for (d=1; d<naxis; d++)
      rawposp[l][d] = rawpos[d];
    writeflag[l] = 1;
    }
  QPTHREAD_MUTEX_UNLOCK(&linemutex);
  if (!(y%dispstep))
    NPRINTF(OUTPUT, "\33[1M> Resampling line:%7d / %-7d\n\33[1A",
	y, height);

  return l;
  }


/****** cancel_resample_threads ***********************************************
PROTO   void clean_resample_threads(void)
PURPOSE Cancel remaining active threads
INPUT   -.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 25/04/2002
 ***/
void    cancel_resample_threads(void)
  {
   int  p;

  for (p=0; p<nproc; p++)
    QPTHREAD_CANCEL(thread[p]);

  return;
  }

#endif

/****** warp_line *************************************************************
PROTO	void warp_line(int p)
PURPOSE	Resample an image line.
INPUT	Thread number.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	05/02/2012
 ***/
void	warp_line(int p)
  {
   wcsstruct		*wcsin,*wcsout;
   double		rawposover[NAXIS], wcspos[NAXIS],
			*rawpos, *rawbufc, *oversampt,*oversampwt, *rawbufareac,
			area, worldc;
   PIXTYPE		*out, *outw,
			pix,pixw;
   FLAGTYPE		*outi,*outwi, *oversampit,*oversampwit,
			ipix,ipixw;
   int			nstepover[NAXIS],stepcount[NAXIS],
			*oversampnt,
			d, o, x, ninput, swapflag;

  if (iflag)
    {
    outi = outibuf[p];
    outwi = outwibuf[p];
    }
  else
    {
    out = outbuf[p];
    outw = outwbuf[p];
    }
  rawpos = rawposp[p];
  wcsin = wcsinp[p];
  wcsout = wcsoutp[p];
  area = infield->fascale;
/* Check if lng and lat are swapped between in and out wcs (vicious idea!) */
  swapflag = (((wcsin->lng != wcsout->lng) || (wcsin->lat != wcsout->lat))
	&& (wcsin->lng != wcsin->lat) && (wcsout->lng != wcsout->lat));

/* Compute all alpha's and delta's for the current line */
  rawpos[0] = rawmin[0];
  if (oversampflag)
    {
/*-- Oversampling is on */
    for (d=0; d<naxis; d++)
      {
      stepcount[d] = nstepover[d] = 0;
      rawposover[d] = rawpos[d];
      }
    if (iflag)
      {
      memset(oversampibuf[p], 0, sizeof(FLAGTYPE)*width);
      memset(oversampwibuf[p], 0, sizeof(FLAGTYPE)*width);
      }
    else
      {
      memset(oversampbuf[p], 0, sizeof(double)*width);
      memset(oversampwbuf[p], 0, sizeof(double)*width);
      }
    memset(oversampnbuf[p], 0, sizeof(int)*width);
    for (o=noversamp; o--; )
      {
      rawbufc = rawbuf[p];
      rawbufareac = rawbufarea[p];
      if (approxflag)
/*------ With approximation */
        projapp_line(projapp, rawposover, 1.0, width, rawbufc, rawbufareac);
      else
        {
/*------ Without approximation */
        for (x=width; x--; *rawposover+=1.0,  rawbufc+=naxis)
          {
          raw_to_wcs(wcsout, rawposover, wcspos);
          if (*wcspos == WCS_NOCOORD)
            *rawbufc = WCS_NOCOORD;
          else
            {
            if (swapflag)
              {
              worldc = wcspos[wcsout->lat];
              wcspos[wcsout->lat] = wcspos[wcsin->lat];
              wcspos[wcsin->lat] = worldc;
              }
            wcs_to_raw(wcsin, wcspos, rawbufc);
            if (rawbufareac)
              *rawbufareac = wcs_scale(wcsout, rawposover)
			/ wcs_scale(wcsin, rawbufc);
            }
          if (rawbufareac)
            rawbufareac++;
          }
	}
/*---- Resample the current line */
      if (iflag)
        {
        oversampit = oversampibuf[p];
        oversampwit = oversampwibuf[p];
        oversampnt = oversampnbuf[p];
        rawbufc = rawbuf[p];
        rawbufareac = rawbufarea[p];
        for (x=width; x--; rawbufc+=naxis)
          {
          if (rawbufareac)
            area = *(rawbufareac++);
          if (*rawbufc != WCS_NOCOORD
		&& interpolate_ipix(infield, inwfield, rawbufc,&ipix,&ipixw)
			== RETURN_OK)
            {
            *(oversampit++) |= ipix;
            *(oversampwit++) |= 1;
            (*(oversampnt++))++;
            }
          else
            {
            oversampit++;
            oversampwit++;
            oversampnt++;
            }
          }
        }
      else
        {
        oversampt = oversampbuf[p];
        oversampwt = oversampwbuf[p];
        oversampnt = oversampnbuf[p];
        rawbufc = rawbuf[p];
        rawbufareac = rawbufarea[p];
        for (x=width; x--; rawbufc+=naxis)
          {
          if (rawbufareac)
            area = *(rawbufareac++);
          if (*rawbufc != WCS_NOCOORD
		&& (interpolate_pix(infield, inwfield, ikernel[p],rawbufc,
			&pix,&pixw),pixw<BIG))
            {
            *(oversampt++) += area * (double)pix;
            *(oversampwt++) += (double)pixw * area*area;
            (*(oversampnt++))++;
            }
          else
            {
            oversampt++;
            oversampwt++;
            oversampnt++;
            }
          }
        }
/*---- Update coordinate vector */
      for (d=0; d<naxis; d++)
        {
        rawposover[d] = rawpos[d] + (++stepcount[d])*stepover[d];
        if (++nstepover[d]<oversamp[d])
          break;
        else
          {
          stepcount[d] = nstepover[d] = 0; /* No need to initialize it to 0! */
          rawposover[d] = rawpos[d];
          }
        }
      }
/*-- Now transfer to the output line */
    if (iflag)
      {
      oversampit = oversampibuf[p];
      oversampwit = oversampwibuf[p];
      oversampnt = oversampnbuf[p];
      for (x=width; x--; oversampit++, oversampwit++)
        {
        if ((ninput = *(oversampnt++)))
          {
          *(outi++) = *oversampit;
/*--------- Convert variances to weight */
          *(outwi++) = *oversampwit;
          }
        else
          *(outwi++) = *(outi++) = 0;
        }
      }
    else
      {
      oversampt = oversampbuf[p];
      oversampwt = oversampwbuf[p];
      oversampnt = oversampnbuf[p];
      for (x=width; x--; oversampt++, oversampwt++)
        {
        if ((ninput = *(oversampnt++)))
          {
          *(out++) = (PIXTYPE)(*oversampt/ninput);
/*--------- Convert variances to weight */
          *(outw++) = (PIXTYPE)((ninput / *oversampwt)*ascale);
          }
        else
          *(out++) = *(outw++) = 0.0;
        }
      }
    }
  else
    {
/*-- No oversampling */
    rawbufc = rawbuf[p];
    rawbufareac = rawbufarea[p];
    if (approxflag)
/*---- With approximation */
      projapp_line(projapp, rawpos, 1.0, width, rawbufc, rawbufareac);
    else
/*---- Without approximation */
      for (x=width; x--; (*rawpos)+=1.0,  rawbufc+=naxis)
        {
        raw_to_wcs(wcsout, rawpos, wcspos);
        if (*wcspos == WCS_NOCOORD)
          *rawbufc = WCS_NOCOORD;
        else
          {
          if (swapflag)
            {
            worldc = wcspos[wcsout->lat];
            wcspos[wcsout->lat] = wcspos[wcsin->lat];
            wcspos[wcsin->lat] = worldc;
            }
          wcs_to_raw(wcsin, wcspos, rawbufc);
          if (rawbufareac)
            *rawbufareac = wcs_scale(wcsout, rawpos)
			/ wcs_scale(wcsin, rawbufc);
          }
        if (rawbufareac)
          rawbufareac++;
        }
/*-- Resample the line */
    rawbufc = rawbuf[p];
    rawbufareac = rawbufarea[p];
    if (iflag)
      for (x=width; x--; rawbufc+=naxis)
        {
        if (rawbufareac)
          area = *(rawbufareac++);
        if (*rawbufc != WCS_NOCOORD)
          interpolate_ipix(infield, inwfield, rawbufc,outi++,outwi++);
        else
          *(outwi++) = *(outi++) = 0;
        }
    else
      for (x=width; x--; rawbufc+=naxis)
        {
        if (rawbufareac)
          area = *(rawbufareac++);
        if (*rawbufc != WCS_NOCOORD)
          {
          interpolate_pix(infield, inwfield, ikernel[p], rawbufc,out,outw);
          *(out++) *= area;
/*------- Convert variance to weight */
          *outw = (*outw < BIG) ? 1.0/(*outw*area*area) : 0.0;
          outw++;
          }
        else
          *(out++) = *(outw++) = 0.0;
        }
    }

  return;
  }

