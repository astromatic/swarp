/*
                                  coadd.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*       Part of:        Swarp
*
*       Author:         E.BERTIN (IAP)
*
*       Contents:       Coaddition routines
*
*       Last modify:    16/07/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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
#include "coadd.h"
#include "data.h"
#include "field.h"
#include "header.h"
#include "interpolate.h"
#include "prefs.h"
#ifdef USE_THREADS
#include "threads.h"
#endif
#include "weight.h"
#include "wcs/wcs.h"

#define	ARRAY_BIT(x,y)		(array[(y/(8*sizeof(unsigned int)))*nnode+x]&(1 \
					<<(y%(8*sizeof(unsigned int)))))
#define	SET_ARRAY_BIT(x,y)	(array[(y/(8*sizeof(unsigned int)))*ninput+x] \
					|= (1<<(y%(8*sizeof(unsigned int)))))

 coaddenum	coadd_type;
 PIXTYPE	*multibuf,*multiwbuf, *outbuf,*outwbuf,
		coadd_wthresh;
 unsigned int	*multinbuf;
 int		coadd_nomax, coadd_width;
 char		padbuf[FBSIZE];

#ifdef USE_THREADS
 pthread_t		*thread,
			movthread;
 pthread_mutex_t	coaddmutex;
 threads_gate_t		*pthread_startgate, *pthread_stopgate,
			*pthread_startgate2, *pthread_stopgate2;
 PIXTYPE		*pthread_linebuf, *pthread_multibuf;
 unsigned int		*pthread_multinbuf;
 int			pthread_bufline, pthread_nbuflines, pthread_npix,
			pthread_step, pthread_endflag, pthread_wdataflag;
#endif

/*------------------------------ function -----------------------------------*/
 static int	coadd_line(int l);

 static PIXTYPE	fast_median(PIXTYPE *arr, int n);
 static int	coadd_load(fieldstruct *field, fieldstruct *wfield,
			PIXTYPE *multibuf, PIXTYPE *multiwbuf,
			unsigned int *multinbuf,
			int *rawpos, int *rawmin, int *rawmax,
			int nlines, int outwidth, int multinmax);
static void	max_clique_recur(unsigned int *array, int nnode, int *old,
			int ne, int ce, int **compsub, int *ncompsub,
			int **best, int *nbest);

#ifdef USE_THREADS
 static void	*pthread_coadd_lines(void *arg),
		*pthread_move_lines(void *arg);
#endif

/******* coadd_fields *********************************************************
PROTO	int coadd_fields(fieldstruct **infield, fieldstruct **inwfield,
			int ninput,
			fieldstruct *outfield, fieldstruct *outwfield,
			coaddenum coaddtype, PIXTYPE wthresh)
PURPOSE	Coadd images.
INPUT	Input field ptr array,
	Input weight field ptr array,
	number of input fields,
	Output field ptr,
	Output weight field ptr,
	Coaddition type.
OUTPUT	RETURN_OK if no error, or RETURN_ERROR in case of non-fatal error(s).
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 16/07/2007
 ***/
int coadd_fields(fieldstruct **infield, fieldstruct **inwfield,	int ninput,
			fieldstruct *outfield, fieldstruct *outwfield,
			coaddenum coaddtype, PIXTYPE wthresh)

  {
#ifdef USE_THREADS
   static pthread_attr_t	pthread_attr;
   int				*proc,
				p;
#endif
   wcsstruct		*wcs;
   PIXTYPE		*emptybuf,*outline,
			*pix,*wpix;
   double		exptime, w,w1,w2, mw, satlev;
   unsigned int		*cflag,*array,
			d, n, n1,n2, flag;
   int			bufmin[NAXIS], bufmax[NAXIS], bufpos[NAXIS],
			rawmax[NAXIS], rawpos[NAXIS], rawpos2[NAXIS],
			min1[NAXIS],max1[NAXIS],
			*ybegbufline,*yendbufline, *maxclique,
			y, y2,dy, ybuf,ybufmax,
			outwidth,multiwidth, width, height,min,max,
			naxis, nlines, nlinesmax,
			nbuflines,nbuflines2,nbuflinesmax, size, omax,omax2,
			offbeg, offend, fieldno;

  coadd_type = coaddtype;
  coadd_wthresh = wthresh;
  naxis = outfield->tab->naxis;
/* The output width is the useful length of the NAXIS1 axis */
  width = outfield->width;
/* The output ``height'' is the product of all other axis lengths */
  height = outfield->height;

/* Find global extrema of mapped image */
  for (d=0; d<naxis; d++)
    {
    min = 1<<30;
    max = 0;
    for (n=0; n<ninput; n++)
      {
      wcs = infield[n]->wcs;
/*---- Update input frame limits in output frame */
      wcs->outmin[d] = (int)floor(outfield->wcs->crpix[d] - wcs->crpix[d] + 1.0);
      if (wcs->outmin[d] < min)
        min = (int)wcs->outmin[d];
      wcs->outmax[d] = wcs->outmin[d] + wcs->naxisn[d] - 1;
      if (wcs->outmax[d] > max)
        max = (int)wcs->outmax[d];
      }
    rawpos[d] = 1;
    rawmax[d] = outfield->wcs->naxisn[d];
    bufmin[d] = min>0 ? (min<rawmax[d] ? min : rawmax[d]) : 1;
    bufmax[d] = max>0 ? (max<rawmax[d] ? max : rawmax[d]) : 1;
    }

  QMALLOC(ybegbufline, int, ninput);
  QMALLOC(yendbufline, int, ninput);
  for (n=0; n<ninput; n++)
    {
    wcs = infield[n]->wcs;
    nbuflines = nbuflines2 = 0;
    nlines = 1;
    for (d=1; d<naxis; d++)
      {
      if (d>1)
        nlines *= (bufmax[d-1] - bufmin[d-1]) + 1;
      nbuflines +=  (wcs->outmin[d]-bufmin[d])*nlines;
      nbuflines2 += (wcs->outmax[d]-bufmin[d])*nlines;
      }
    ybegbufline[n] = nbuflines;
    yendbufline[n] = nbuflines2;
    }

/* Empty pixels at the beginning of each line */
  offbeg = bufmin[0]-1;

/* Empty pixels at the end of each line */
  offend = outfield->wcs->naxisn[0] - bufmax[0];

  coadd_width = outwidth = outfield->wcs->naxisn[0] - offend - offbeg;

/* Build a graph of overlaps */
  QCALLOC(array, unsigned int, ninput*(1+(ninput-1)/(8*sizeof(unsigned int))));
  for (n1=0; n1<ninput; n1++)
    {
    SET_ARRAY_BIT(n1,n1);
    wcs = infield[n1]->wcs;
    for (d=0; d<naxis; d++)
      {
      min1[d] = wcs->outmin[d];
      max1[d] = wcs->outmax[d];
      }
    for (n2=n1+1; n2<ninput; n2++)
      {
      wcs = infield[n2]->wcs;
      flag = 1;
      for (d=0; d<naxis; d++)
        if (wcs->outmin[d] > max1[d] || wcs->outmax[d] < min1[d])
	  {
          flag = 0;
          break;
          }
      if (flag)
        {
        SET_ARRAY_BIT(n1,n2);
	SET_ARRAY_BIT(n2,n1);
        }
      }
    }

/* Find the densest overlap, that is the maximal clique (NP-complete pb!) */
  omax = max_clique(array, ninput, &maxclique);
  free(array);

/* Compute maximum gain and maximum total exposure time */
  exptime = w1 = w2 = mw = 0.0;
  fieldno = -1;
  omax2 = 0;
  for (n=0; n<omax; n++)
    {
/*-- Keep only one extension per file */
    n1 = maxclique[n];
    if (infield[n1]->fieldno == fieldno)
      continue;
    omax2++;
    fieldno = infield[n1]->fieldno;
    exptime += infield[n1]->exptime;
    w = (coaddtype == COADD_WEIGHTED && infield[n1]->backsig > 0.0)?
	1.0/(infield[n1]->backsig*infield[n1]->backsig) : 1.0;
    w1 += w;
    if (infield[n1]->gain > 0.0)
      {
      w2 += w*w/infield[n1]->gain;
      mw += infield[n1]->gain;
      }
    }
  free(maxclique);

  outfield->fieldno = omax2;	/* This is not what is meant for but hey */
  outfield->exptime = exptime;

/* Approximation to the final equivalent gain */
  outfield->gain = 0.0;
  if (coaddtype == COADD_WEIGHTED || coaddtype == COADD_AVERAGE)
    {
    if (w2 > 0.0)
      outfield->gain = w1*w1/w2;
    }
  else if (coaddtype == COADD_MEDIAN)
    {
    if (w2 > 0.0)
      outfield->gain = w1*w1/w2;
    if (omax2 > 2)
      outfield->gain /= PI;
    }
  else if (coaddtype == COADD_SUM)
    {
    if (w2 > 0.0)
      outfield->gain = w1*w1/w2/omax2;
    }
  else
    outfield->gain = mw/omax2;

/* Compute output saturation level (the minimum on all saturation) */
  satlev = BIG;
  for (n=0; n<ninput; n++)
    {
    infield[n]->saturation = (prefs.subback_flag[n])? (infield[n]->saturation
                     - infield[n]->backmean) : infield[n]->saturation;
    if (infield[n]->saturation < satlev)
      satlev = infield[n]->saturation;
    }
  outfield->saturation = satlev;

/* Add relevant information to output FITS headers */
  writefitsinfo_outfield(outfield, *infield);
  writefitsinfo_outfield(outwfield, inwfield? *inwfield : *infield);

  *gstr = '\0';
  NPRINTF(OUTPUT, "-------------- Co-adding frames            \n");
  NPRINTF(OUTPUT, "Maximum overlap density: %d frame%s\n",
	omax2, omax2>1? "s" : "");

  coadd_nomax = omax;
  multiwidth = outwidth*omax;
  nbuflinesmax = (prefs.coaddbuf_size*1024*1024)
	/ (2*(multiwidth+outwidth)*sizeof(PIXTYPE));
  if (nbuflinesmax < 1)
    nbuflinesmax = 1;
#ifdef USE_THREADS
/* Number of active threads */
  nproc = prefs.nthreads;
/* Make sure that the number of multibuffer lines is a multiple of the */
/* number of processes */
  nbuflinesmax += (nproc - (nbuflinesmax%nproc))%nproc;
#endif
/* Allocate memory for the "multi-buffers" storing "packed" pixels from all */
/* images for the current line(s), prior to co-addition */
  QMALLOC(multibuf, PIXTYPE, nbuflinesmax*multiwidth);
  QMALLOC(multiwbuf, PIXTYPE, nbuflinesmax*multiwidth);
  QMALLOC(multinbuf, unsigned int, nbuflinesmax*outwidth);
/* Allocate memory for the output buffers that contain "empty data" */
  QCALLOC(emptybuf, PIXTYPE, width);
  QCALLOC(outline, PIXTYPE, width);
/* Allocate memory for the output buffers that contain the final data in */
/* internal format (PIXTYPE) */
  QMALLOC(outbuf, PIXTYPE, nbuflinesmax*outwidth);
  QMALLOC(outwbuf, PIXTYPE, nbuflinesmax*outwidth);
  QCALLOC(cflag, unsigned int, ninput);

/* Open output file and save header */
  if (open_cat(outfield->cat, WRITE_ONLY) != RETURN_OK)
      error(EXIT_FAILURE, "*Error*: cannot open for writing ",
		outfield->filename);
  QFWRITE(outfield->tab->headbuf, outfield->tab->headnblock*FBSIZE,
	outfield->cat->file, outfield->filename);

/* Open output weight file and save header */
  outwfield->sigfac = (double)1.0;	/* A possible scaling among others */
  if (open_cat(outwfield->cat, WRITE_ONLY) != RETURN_OK)
    error(EXIT_FAILURE, "*Error*: cannot open for writing ",
		outwfield->filename);
  QFWRITE(outwfield->tab->headbuf, outwfield->tab->headnblock*FBSIZE,
	outwfield->cat->file, outwfield->filename);

/* Only for WRITING the weights */
  set_weightconv(outwfield);

/* Start all threads! */
#ifdef USE_THREADS
/* Set up multi-threading stuff */
  QPTHREAD_MUTEX_INIT(&coaddmutex, NULL);
  QPTHREAD_ATTR_INIT(&pthread_attr);
  QPTHREAD_ATTR_SETDETACHSTATE(&pthread_attr, PTHREAD_CREATE_JOINABLE);
  pthread_startgate = threads_gate_init(nproc+1, NULL);
  pthread_stopgate = threads_gate_init(nproc+1, NULL);
  pthread_startgate2 = threads_gate_init(2, NULL);
  pthread_stopgate2 = threads_gate_init(2, NULL);
  QMALLOC(proc, int, nproc);
  QMALLOC(thread, pthread_t, nproc);
  pthread_bufline = pthread_nbuflines = nbuflinesmax;
  pthread_endflag = 0;
/* Start the co-addition threads */
  for (p=0; p<nproc; p++)
    {
    proc[p] = p;
    QPTHREAD_CREATE(&thread[p], &pthread_attr, &pthread_coadd_lines, &proc[p]);
    }
/* Start the data mover thread */
  QPTHREAD_CREATE(&movthread, &pthread_attr, &pthread_move_lines, &p);
#endif

/* Loop over output ``lines'': this can be over more than 1 (Y) dimension */
  for (y=ybufmax=0; y<height; y+=nlines)
    {
    NPRINTF(OUTPUT, "\33[1M> Preparing line:%7d / %-7d\n\33[1A", y+1, height);
/*-- Skip empty lines */
    for (d=naxis; --d;)
      rawpos2[d] = rawpos[d];
    nlinesmax = height-y;
    for (nbuflines=nlines=0; nbuflines<nbuflinesmax && nlines<nlinesmax;
		nlines++)
      {
      for (d=naxis; --d;)
        if (rawpos2[d]<bufmin[d] || rawpos2[d]>bufmax[d])
          break;
      if (d<=0 && !(nbuflines++))
/*------ bufpos[] stores the coordinates of the next non-empty line */
        for (d=naxis; --d;)
            bufpos[d] = rawpos2[d];
      for (d=1; d<naxis; d++)
        if ((++rawpos2[d])<=rawmax[d])
          break;
        else
          rawpos2[d] = 1;
      }
    ybuf = ybufmax;
    ybufmax = ybuf+nbuflines;
    if (multiwidth)
/*---- Initialize output data line */
      memset(multinbuf, 0, outwidth*nbuflines*sizeof(unsigned int));
/*-- Examine the batch of input images for the current output image section */
    NPRINTF(OUTPUT, "\33[1M> Reading   line:%7d / %-7d\n\33[1A", y+1,height);
    for (n=0; n<ninput; n++)
      if (!(cflag[n] & COADDFLAG_FINISHED))
        {
        wcs = infield[n]->wcs;
/*------ Focus on images that begin before the current buffer ends */
        if (ybegbufline[n] < ybufmax)
	  {
/*-------- Discard images entirely below the current lower limit */
          if (yendbufline[n] < ybuf)
            {
            cflag[n] |= COADDFLAG_FINISHED;
            continue;
            }
/*-------- Open images if needed */
          if (!(cflag[n] & COADDFLAG_OPEN))
	    {
            if (open_cat(infield[n]->cat, READ_ONLY) != RETURN_OK)
              error(EXIT_FAILURE, "*Error*: cannot open for reading ",
	      infield[n]->filename);
            if (inwfield[n])
	      {
              if (open_cat(inwfield[n]->cat, READ_ONLY) != RETURN_OK)
			error(EXIT_FAILURE,"*Error*: cannot open for reading ",
				inwfield[n]->filename);
 	      }
            cflag[n] |= COADDFLAG_OPEN;
            dy = ybegbufline[n] - ybuf;
            if (dy<0)
              dy = 0;
	    }
          else
            dy = 0;	/* We just keep reading the file */
          nbuflines2 = yendbufline[n] - ybuf + 1;
          if (nbuflines2 > nbuflines)
            nbuflines2 = nbuflines;
          nbuflines2 -= dy;
/*-------- Refill the buffers with new data */
          if (coadd_load(infield[n], inwfield[n],
			multibuf+dy*multiwidth, multiwbuf+dy*multiwidth,
			multinbuf+dy*outwidth,
			bufpos, bufmin, bufmax, nbuflines2, outwidth, omax)
		!= RETURN_OK)
/*-------- End of the image, we can close the file */
            {
            close_cat(infield[n]->cat);
            if (inwfield[n])
               close_cat(inwfield[n]->cat);
            cflag[n] ^= COADDFLAG_OPEN;
            cflag[n] |= COADDFLAG_FINISHED;
	    }
	  }
	}

    NPRINTF(OUTPUT, "\33[1M> Co-adding line:%7d / %-7d\n\33[1A", y+1,height);
/*-- Now perform the coaddition itself */
#ifdef USE_THREADS
    QPTHREAD_MUTEX_LOCK(&coaddmutex);
    pthread_bufline = 0;
    pthread_nbuflines = nbuflines;
    QPTHREAD_MUTEX_UNLOCK(&coaddmutex);
    threads_gate_sync(pthread_startgate);
/* ( Slave threads process the current buffer data here ) */
    threads_gate_sync(pthread_stopgate);
#else
    for (y2=0; y2<nbuflines; y2++)
      coadd_line(y2);
#endif
    NPRINTF(OUTPUT, "\33[1M> Writing   line:%7d / %-7d\n\33[1A", y+1,height);
/*-- Write the image buffer lines */
    pix = outbuf;
    for (d=naxis; --d;)
      rawpos2[d] = rawpos[d];
    for (y2=nlines; y2--;)
      {
/*---- Skip empty lines */
      for (d=naxis; --d;)
        if (rawpos2[d]<bufmin[d] || rawpos2[d]>bufmax[d])
          break;
      if (d>0)
        write_body(outfield->tab, emptybuf, width);
      else
        {
        memcpy(outline+offbeg, pix, outwidth*sizeof(PIXTYPE));
        write_body(outfield->tab, outline, width);
        pix += outwidth;
	}
/*---- Update coordinate vector */
      for (d=1; d<naxis; d++)
        if ((++rawpos2[d])<=rawmax[d])
          break;
        else
          rawpos2[d] = 1;
      }

    wpix = outwbuf;
/*-- The weight buffer lines */
    for (y2=nlines; y2--;)
      {
/*---- Skip empty lines */
      for (d=naxis; --d;)
        if (rawpos[d]<bufmin[d] || rawpos[d]>bufmax[d])
          break;
      if (d>0)
        write_body(outwfield->tab, emptybuf, width);
      else
        {
        var_to_weight(wpix, outwidth);
        memcpy(outline+offbeg, wpix, outwidth*sizeof(PIXTYPE));
        write_body(outwfield->tab, outline, width);
        wpix += outwidth;
	}
/*---- Update coordinate vector */
      for (d=1; d<naxis; d++)
        if ((++rawpos[d])<=rawmax[d])
          break;
        else
          rawpos[d] = 1;
      }
    }


/* FITS padding*/
  size = PADEXTRA(outfield->tab->tabsize);
  if (size)
    QFWRITE(padbuf, (size_t)size, outfield->cat->file, outfield->filename);
  size = PADEXTRA(outwfield->tab->tabsize);
  if (size)
    QFWRITE(padbuf, (size_t)size, outwfield->cat->file, outwfield->filename);

/* Close files */
  close_cat(outfield->cat);
  close_cat(outwfield->cat);
  for (n = 0; n<ninput; n++)
    {
    close_cat(infield[n]->cat);
    if (inwfield[n])
      close_cat(inwfield[n]->cat);
    }

#ifdef USE_THREADS
  pthread_endflag = 1;
/* (Re-)activate existing threads... */
  threads_gate_sync(pthread_startgate);
  threads_gate_sync(pthread_startgate2);
/* ... and shutdown all threads */
  for (p=0; p<nproc; p++)
    QPTHREAD_JOIN(thread[p], NULL);
  QPTHREAD_JOIN(movthread, NULL);
  threads_gate_end(pthread_startgate);
  threads_gate_end(pthread_stopgate);
  threads_gate_end(pthread_startgate2);
  threads_gate_end(pthread_stopgate2);
  QPTHREAD_MUTEX_DESTROY(&coaddmutex);
  QPTHREAD_ATTR_DESTROY(&pthread_attr);
  free(proc);
  free(thread);
#endif

/* Free Buffers */
  free(cflag);
  free(ybegbufline);
  free(yendbufline);
  free(emptybuf);
  free(outline);
  if (outwidth)
    {
    free(multibuf);
    free(multiwbuf);
    free(multinbuf);
    free(outbuf);
    free(outwbuf);
    }

  return RETURN_OK;
  }


/******* max_clique *********************************************************
PROTO	int max_clique(int *array, int nnode, int nnodes, int **max)
PURPOSE	Return the number of nodes of the maximal clique of a subgraph. 
INPUT	Input graph array,
	Number of graph nodes,
	Pointer to the list of nodes in the maximal clique found (output).
OUTPUT	Number of nodes in the maximal clique found.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 05/08/2006
 ***/
int	max_clique(unsigned int *array, int nnode, int **max)
  {
   int	*all, *compsub,
	c, ncompsub, nmax;

  QMALLOC(all, int, nnode);
  QMALLOC(compsub, int, nnode)
  ncompsub = 0;
  for (c=0; c<nnode; c++)
    all[c] = c;
  *max = NULL;
  nmax = 0;

/* Run the recursive part */
  max_clique_recur(array, nnode, all, 0, nnode, &compsub, &ncompsub, max, &nmax);

  free(compsub);
  free(all);

  return nmax;
  }


/******* max_clique_recur ****************************************************
PROTO	void max_clique_recur(int *array, int nnode, int *old, int ne,int ce,
		int **compsub, int *ncompsub, int **best, int *nbest)
PURPOSE	Recursive part of the Bron & Kerbosch 1973 algorithm
	(Comm. of the ACM 16, 575), with graph access modified to optimize
	memory usage.
INPUT	Input graph array,
	Number of graph nodes,
	Working list of nodes,
	Start node,
	End node,
	Pointer to another working list of nodes,
	Number of working nodes,
	Pointer to the largest list of clique nodes,
	Number of clique nodes.
OUTPUT	-.
NOTES   Recursive function.
AUTHOR  E. Bertin (IAP)
VERSION 05/08/2006
 ***/
static void	max_clique_recur(unsigned int *array, int nnode, int *old,
		int ne, int ce, int **compsub, int *ncompsub, int **best,
		int *nbest)
  {
   int	*new,
	i,j,p,s, count, fixp, minnod, newce, newne, nod, pos, sel;

  QMALLOC(new, int, nnode);
  minnod = ce;
  nod = 0;
  pos = fixp = s = 0;	/* To avoid gcc -Wall warnings */
/* Determine each counter value and look for minimum */
  for (i=0; i<ce && minnod != 0; i++)
    {
    p = old[i];
    count = 0;
/*-- Count disconnections */
    for (j=ne; j<ce && count < minnod; j++)
      {
      if (!ARRAY_BIT(p,old[j]))
        {
        count++;
/*------ Save position of potential candidate */
        pos = j;
        }
      }
/*-- Test new minimum */
    if (count < minnod)
      {
      fixp = p;
      minnod = count;
      if (i<ne)
        s = pos;
      else
        {
        s = i;
/*------ Preincr */
        nod = 1;
        }
      }
    }
/* If fixed point initially chosen from candidates then number of */
/* disconnections will be preincreased by one */
/*-- Backtrackcycle */
  for (nod+=minnod; nod>=1; nod--)
    {
/*-- Interchange */
    p = old[s];
    old[s] = old[ne];
    sel = old[ne] = p;

/*-- Fill new set "not" */
    newne = 0;
    for (i=0; i<ne; i++)
      if (ARRAY_BIT(sel,old[i]))
        new[newne++] = old[i];

/*-- Fill new set "cand" */
    newce = newne;
    for (i=ne+1; i<ce; i++)
      if (ARRAY_BIT(sel,old[i]))
        new[newce++] = old[i];

/*-- Add to compsub */
    (*compsub)[(*ncompsub)++] = sel;

    if (!newce)
      {
      if (*nbest < *ncompsub)
/*------ Found a max clique */
        {
        *nbest = *ncompsub;
        if (*best)
          free(*best);
        QMEMCPY(*compsub, *best, int, *nbest);
        }
      }
    else if (newne < newce)
      max_clique_recur(array, nnode, new, newne, newce, compsub, ncompsub,
		best, nbest);
/*-- Remove from "compsub" */
    (*ncompsub)--;
/*-- Add to "not" */
    ne++;
    if (nod > 1)
/*---- Select a candidate disconnected to the fixed point */
      for (s=ne; ARRAY_BIT(fixp,old[s]); s++);
    }
  free(new);

  return;
  }


#ifdef USE_THREADS

/****** pthread_coadd_lines ****************************************************
PROTO	void *pthread_coadd_lines(void *arg)
PURPOSE	thread that takes care of coadding image "lines"
INPUT	Pointer to the thread number.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	30/07/2006
 ***/
void	*pthread_coadd_lines(void *arg)
  {
   int	bufline, proc;

  bufline = -1;
  proc = *((int *)arg);
  threads_gate_sync(pthread_startgate);
  while (!pthread_endflag)
    {
    QPTHREAD_MUTEX_LOCK(&coaddmutex);
    if (pthread_bufline<pthread_nbuflines)
      {
      bufline = pthread_bufline++;
      QPTHREAD_MUTEX_UNLOCK(&coaddmutex);
      coadd_line(bufline);
      }
    else
      {
      QPTHREAD_MUTEX_UNLOCK(&coaddmutex);
/*---- Wait for the input buffer to be updated */
      threads_gate_sync(pthread_stopgate);
/* ( Master thread process loads and saves new data here ) */
      threads_gate_sync(pthread_startgate);
      }
    }

  pthread_exit(NULL);

  return (void *)NULL;
  }

#endif

/******* coadd_line **********************************************************
PROTO	int coadd_line(int l)
PURPOSE	Coadd a line of pixels.
INPUT	Current line number.
OUTPUT	RETURN_OK if no error, or RETURN_ERROR in case of non-fatal error(s).
NOTES   Requires many global variables (for multithreading).
AUTHOR  E. Bertin (IAP)
VERSION 27/07/2006
 ***/
int coadd_line(int l)

  {
   PIXTYPE		*inpix,*inwpix,*inpixt,*inwpixt,
			*pixstack, *pixt, *outpix,*outwpix;
   double		val, val2, wval,wval2,wval3;
   unsigned int		*inn;
   int			i,x, ninput2, blankflag;

  blankflag = prefs.blank_flag;
  inpix = multibuf+l*coadd_width*coadd_nomax;
  inwpix = multiwbuf+l*coadd_width*coadd_nomax;
  inn = multinbuf+l*coadd_width;
  outpix = outbuf+l*coadd_width;
  outwpix = outwbuf+l*coadd_width;
  switch(coadd_type)
    {
    case COADD_WEIGHTED:
      for (x=coadd_width; x--; inpix+=coadd_nomax, inwpix+=coadd_nomax)
        {
        ninput2 = 0;
        wval = val = val2 = 0.0;
        inpixt = inpix;
        inwpixt = inwpix;
        for (i=*(inn++); i--;)
          {
          val2 = *(inpixt++);
          wval2 = *(inwpixt++);
          if (wval2<coadd_wthresh)
            {
            ninput2++;
            wval += (wval2=1.0/wval2);
            val += val2*wval2;
            }
          }
        if (ninput2)
          {
          *(outwpix++) = (wval= 1.0/wval);
          *(outpix++) = val*wval;
          }
        else
          {
          *(outpix++) = blankflag? 0.0 : val2;
          *(outwpix++) = BIG;
          }
        }
      break;
    case COADD_MEDIAN:
      QMALLOC(pixstack, PIXTYPE, coadd_nomax);
      for (x=coadd_width; x--; inpix+=coadd_nomax, inwpix+=coadd_nomax)
        {
        ninput2 = 0;
        wval = wval3 = val2 = 0.0;
        inpixt = inpix;
        inwpixt = inwpix;
        pixt = pixstack;
        for (i=*(inn++); i--;)
          {
          val2 = *(inpixt++);
          wval2 = *(inwpixt++);        
          if (wval2 < coadd_wthresh)
            {
            *(pixt++) = val2;
            wval += 1.0/sqrt(wval2);
            wval3 += wval2;
            ninput2++;
            }      
          }
        if (ninput2)
          {
          *(outpix++) = fast_median(pixstack, ninput2);
/*-------- We assume Gaussian input noise */
          *(outwpix++) = ninput2>2?
		(PI*ninput2*ninput2/(2*wval*wval*(ninput2+((ninput2&1)?(PI/2-1)
			:(PI-2)))))
		  : wval3/(ninput2*ninput2);
          }
        else
          {
          *(outpix++) = blankflag? 0.0 : val2;
          *(outwpix++) = BIG;
          }
        }
      free(pixstack);
      break;
    case COADD_AVERAGE:
      for (x=coadd_width; x--; inpix+=coadd_nomax, inwpix+=coadd_nomax)
        {
        ninput2 = 0;
        wval = val = val2 = 0.0;
        inpixt = inpix;
        inwpixt = inwpix;
        for (i=*(inn++); i--;)
          {
          val2 = *(inpixt++);
          wval2 = *(inwpixt++);
          if (wval2<coadd_wthresh)
            {
            ninput2++;
            val += val2;
            wval += wval2;
            }
          }
        if (ninput2)
          {
          *(outpix++) = val/ninput2;
          *(outwpix++) = wval/(ninput2*ninput2);
          }
        else
          {
          *(outpix++) = blankflag? 0.0 : val2;
          *(outwpix++) = BIG;
          }
        }
      break;
    case COADD_MIN:
      for (x=coadd_width; x--; inpix+=coadd_nomax, inwpix+=coadd_nomax)
        {
        ninput2 = 0;
        val = BIG;
        wval = val2 = 0.0;
        inpixt = inpix;
        inwpixt = inwpix;
        for (i=*(inn++); i--;)
          {
          val2 = *(inpixt++);
          wval2 = *(inwpixt++);
          if (wval2<coadd_wthresh && val2<val)
            {
            ninput2++;
            val = val2;
            }
          }
        if (ninput2)
          {
          *(outpix++) = val;
          *(outwpix++) = 1.0;
          }
        else
          {
          *(outpix++) = blankflag? 0.0 : val2;
          *(outwpix++) = BIG;
          }
        }
      break;
    case COADD_MAX:
      for (x=coadd_width; x--; inpix+=coadd_nomax, inwpix+=coadd_nomax)
        {
        ninput2 = 0;
        val = -BIG;
        wval = val2 = 0.0;
        inpixt = inpix;
        inwpixt = inwpix;
        for (i=*(inn++); i--;)
          {
          val2 = *(inpixt++);
          wval2 = *(inwpixt++);
          if (wval2<coadd_wthresh && val2>val)
            {
            ninput2++;
            val = val2;
            }
          }
        if (ninput2)
          {
          *(outpix++) = val;
          *(outwpix++) = 1.0;
          }
        else
          {
          *(outpix++) = blankflag? 0.0 : val2;
          *(outwpix++) = BIG;
          }
        }
      break;
    case COADD_CHI2:
      for (x=coadd_width; x--; inpix+=coadd_nomax, inwpix+=coadd_nomax)
        {
        ninput2 = 0;
        wval = val = val2 = 0.0;
        inpixt = inpix;
        inwpixt = inwpix;
        for (i=*(inn++); i--;)
          {
          val2 = *(inpixt++);
          wval2 = *(inwpixt++);
          if (wval2<coadd_wthresh)
            {
            ninput2++;
            val += val2*val2/wval2;
            }
          }
        if (ninput2)
          {
          *(outpix++) = sqrt(val/ninput2);
          *(outwpix++) = 1.0;
          }
        else
          {
          *(outpix++) = blankflag? 0.0 : sqrt(val2*val2);
          *(outwpix++) = BIG;
          }
        }
      break;
    case COADD_SUM:
      for (x=coadd_width; x--; inpix+=coadd_nomax, inwpix+=coadd_nomax)
        {
        ninput2 = 0;
        wval = val = val2 = 0.0;
        inpixt = inpix;
        inwpixt = inwpix;
        for (i=*(inn++); i--;)
          {
          val2 = *(inpixt++);
          wval2 = *(inwpixt++);
          if (wval2<coadd_wthresh)
            {
            ninput2++;
            val += val2;
            wval += wval2;
            }
          }
        if (ninput2)
          {
          *(outpix++) = val;
          *(outwpix++) = wval;
          }
        else
          {
          *(outpix++) = blankflag? 0.0 : val2;
          *(outwpix++) = BIG;
          }
        }
      break;
    default:
      error(EXIT_FAILURE, "*Internal Error*: Unknown Combine option in ",
			"coadd_line()");
    }

  return RETURN_OK;
  }


/******* fast_median **********************************************************
PROTO	PIXTYPE fast_median(PIXTYPE *arr, int n)
PURPOSE	Fast median from an input array, optimized version based on the
	select() routine (Numerical Recipes, 2nd ed. Section 8.5 and
	http://www.eso.org/~ndevilla/median/). If n is even, then the result
	is the average of the 2 "central" values.
INPUT	Input pixel array ptr,
	number of input pixels,
OUTPUT	Value of the median.
NOTES	n must be >0. Warning: changes the order of data (but does not sort
	them)!
AUTHOR	E. Bertin (IAP)
VERSION	27/06/2003
 ***/
#define MEDIAN_SWAP(a,b) \
	{val=a; a=b; b=val;}

PIXTYPE fast_median(PIXTYPE *arr, int n) 
  {
   PIXTYPE	*alow, *ahigh, *amedian, *amiddle, *all, *ahh,
		val, valmax, valmax2;
   int		i, nless;

  if (n==1)
    return *arr;
  else if (n==2)
      return 0.5*(arr[0]+arr[1]);
  else if (n==3)
    {
    if (arr[0]>arr[1])
      MEDIAN_SWAP(*arr, *(arr+1));
    if (arr[1]>arr[2])
      MEDIAN_SWAP(*(arr+1), *(arr+2));
    if (arr[0]>arr[1])
      MEDIAN_SWAP(*arr, *(arr+1));
    return arr[1];
    }
  alow = arr;
  ahigh = arr + n - 1;
  amedian = arr + n/2;
  while (ahigh > (all=alow + 1))
    {
/*-- Find median of low, middle and high items; swap into position low */
    amiddle = alow + (ahigh-alow)/2;
    if (*amiddle > *ahigh)
      MEDIAN_SWAP(*amiddle, *ahigh);
    if (*alow > *ahigh)
      MEDIAN_SWAP(*alow, *ahigh);
    if (*amiddle > *alow)
      MEDIAN_SWAP(*amiddle, *alow);

/*-- Swap low item (now in position middle) into position (low+1) */
    MEDIAN_SWAP(*amiddle, *all);

/*-- Nibble from each end towards middle, swapping items when stuck */
    ahh = ahigh;
    for (;;)
      {
      while (*alow > *(++all));
      while (*(--ahh) > *alow);

      if (ahh < all)
        break;

      MEDIAN_SWAP(*all, *ahh);
      }

/*-- Swap middle item (in position low) back into correct position */
    MEDIAN_SWAP(*alow, *ahh) ;

/*-- Re-set active partition */
    if (ahh <= amedian)
      alow = all;
    if (ahh >= amedian)
      ahigh = ahh - 1;
    }

/* One or two elements left */
  if (ahigh == all && *alow > *ahigh)
    MEDIAN_SWAP(*alow, *ahigh);

  if (n&1)
/*-- Odd case */
    return *amedian;
  else
    {
/*-- Even case */
    valmax2 = *amedian;
    valmax = -BIG;
    alow = arr;
    nless = 0;
    for (i=n/2;i--;)
      if ((val=*(alow++))<valmax2)
        {
        nless++;
        if (val > valmax)
          valmax = val;
        }
    return nless<n/2? *amedian : (*amedian+valmax)/2.0;
    }

  }

#undef MEDIAN_SWAP


/******************************** hmedian ***********************************/
/*
Median using Heapsort algorithm (for float arrays) (based on Num.Rec algo.).
Warning: changes the order of data!
*/
PIXTYPE   hmedian(PIXTYPE *ra, int n)

  {
   int          l, j, ir, i;
   float        rra;


  if (n<2)
    return *ra;
  ra--;
  for (l = ((ir=n)>>1)+1;;)
    {
    if (l>1)
      rra = ra[--l];
    else
      {
      rra = ra[ir];
      ra[ir] = ra[1];
      if (--ir == 1)
        {
        ra[1] = rra;
        return n&1? ra[n/2+1] : (ra[n/2]+ra[n/2+1])/2.0;
        }
      }
    for (j = (i=l)<<1; j <= ir;)
      {
      if (j < ir && ra[j] < ra[j+1])
        ++j;
      if (rra < ra[j])
        {
        ra[i] = ra[j];
        j += (i=j);
        }
      else
        j = ir + 1;
      }
    ra[i] = rra;
    }

/* (the 'return' is inside the loop!!) */
  }


/******* coadd_load **********************************************************
PROTO	int coadd_load(fieldstruct *field, fieldstruct *wfield,
			PIXTYPE *multibuf, PIXTYPE *multiwbuf,
			unsigned int *multinbuf,
			int *rawpos, int *bufmin, int *bufmax,
			int nlines, int outwidth, int multinmax)
PURPOSE	Load images and weights to coadd in the current buffer.
INPUT	Input field ptr array,
	Input weight field ptr array,
	Science image buffer,
	Weight image buffer,
	Number-of-pixels buffer,
	Array of current coordinates in output frame,
	Array of minimum coordinates in output frame,
	Array of maximum coordinates in output frame,
	Total number of lines (can be higher than the number of buffers lines),
	Pixel buffer line width,
	Number of overlapping images.
OUTPUT	RETURN_ERROR in case no more data are worth reading,
	RETURN_OK otherwise.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 02/08/2006
 ***/
int	coadd_load(fieldstruct *field, fieldstruct *wfield,
			PIXTYPE *multibuf, PIXTYPE *multiwbuf,
			unsigned int *multinbuf,
			int *rawpos, int *bufmin, int *bufmax,
			int nbuflines, int outwidth, int multinmax)
  {
    wcsstruct		*wcs;
    OFF_T		offset, pixcount;
    PIXTYPE		*linebuf, *line;
    unsigned int	*multinbuf2,
			d, y, nflag, naxis;
    int			rawpos2[NAXIS],
			ival, inoffset, inbeg, muloffset, width;
#ifdef USE_THREADS
    unsigned int	threadstep;
#endif

  wcs = field->wcs;
  naxis = wcs->naxis;
  inoffset = 0;
  inbeg = wcs->outmin[0] - bufmin[0];
  if (inbeg<0)
    {
    inoffset = -inbeg;
    inbeg = 0;
    }
  muloffset = inbeg*multinmax;
  width = wcs->outmax[0] - bufmin[0] + 1;
  if (width > (outwidth - inbeg))
    width = outwidth - inbeg;
  if (width<=0)
    return RETURN_ERROR;	/* the image is not included in the output */
  else if (width > field->width)
    width = field->width;

/* First, the data themselves */
#ifdef USE_THREADS
  pthread_wdataflag = 0;
  pthread_npix = width;
  pthread_step = multinmax;
  threadstep = 0;
  QMALLOC(linebuf, PIXTYPE, 2*field->width);
#else
  QMALLOC(linebuf, PIXTYPE, field->width);
  line = linebuf;
#endif
  for (d=naxis; --d;)
    rawpos2[d] = rawpos[d];
  multinbuf2 = multinbuf;
  nflag = 1;
  for (y=nbuflines; y;)
    {
/*-- Check that the present line is contained in the field */
    for (d=naxis; --d;)
      if (rawpos2[d]<wcs->outmin[d] || rawpos2[d]>wcs->outmax[d])
        break;
    if (d>0)
      nflag = 1;
    else
      {
      if (nflag)
/*------ well now we are back within the frame; we might need to fseek() */
        {
        nflag = 0;
        offset = 0;
        pixcount = 1;
        for (d=1; d<naxis; d++)
          {
          pixcount *= wcs->naxisn[d-1];
          ival = rawpos2[d] - wcs->outmin[d];
          if (ival > 0)
            offset += (OFF_T)ival*pixcount;
          }
        QFSEEK(field->cat->file,
		field->tab->bodypos+offset*field->tab->bytepix,
		SEEK_SET, field->filename);
        }
#ifdef USE_THREADS
      line = linebuf+(threadstep&1)*field->width;
      read_body(field->tab, line, field->width);
      if (threadstep++)
        threads_gate_sync(pthread_stopgate2);
      pthread_linebuf = line+inoffset;
      pthread_multibuf = multibuf+muloffset;
      pthread_multinbuf = multinbuf2+inbeg;
      threads_gate_sync(pthread_startgate2);
#else
      read_body(field->tab, line, field->width);
      coadd_movedata(line+inoffset,
		multibuf+muloffset, multinbuf2+inbeg, width, multinmax);
#endif
      multibuf += outwidth*multinmax;
      multinbuf2 += outwidth;
      y--;
      }

/*---- Update coordinate vector */
    for (d=1; d<naxis; d++)
      if ((++rawpos2[d])<=bufmax[d])
        break;
      else
        rawpos2[d] = bufmin[d];
    }

/* Now the weights */
  for (d=naxis; --d;)
    rawpos2[d] = rawpos[d];
  multinbuf2 = multinbuf;
  nflag = 1;
  for (y=nbuflines; y;)
    {
/*-- Check that the present line is contained in the field */
    for (d=naxis; --d;)
      if (rawpos2[d]<wcs->outmin[d] || rawpos2[d]>wcs->outmax[d])
        break;
    if (d>0)
      nflag = 1;
    else
      {
#ifdef USE_THREADS
      line = linebuf+(threadstep&1)*field->width;
#endif
      if (wfield)
        {
        if (nflag)
/*-------- well now we are back withing the frame; we might need to fseek() */
          {
          nflag = 0;
          offset = 0;
          pixcount = 1;
          for (d=1; d<naxis; d++)
            {
            pixcount *= wcs->naxisn[d-1];
            ival = rawpos2[d] - wcs->outmin[d];
            if (ival > 0)
              offset += ival*pixcount;
            }
          QFSEEK(wfield->cat->file,
		wfield->tab->bodypos+offset*wfield->tab->bytepix,
		SEEK_SET, wfield->filename);
          }
        read_body(wfield->tab, line, field->width);
        }
#ifdef USE_THREADS
      if (threadstep++)
        threads_gate_sync(pthread_stopgate2);
      pthread_wdataflag = 1;
      pthread_linebuf = (wfield? (line+inoffset) : NULL);
      pthread_multibuf = multiwbuf+muloffset;
      pthread_multinbuf = multinbuf2+inbeg;
      threads_gate_sync(pthread_startgate2);
#else
      coadd_movewdata(wfield? (line+inoffset) : NULL,
		multiwbuf+muloffset, multinbuf2+inbeg, width, multinmax);
#endif
      multiwbuf += outwidth*multinmax;
      multinbuf2 += outwidth;
      y--;
      }
/*-- Update coordinate vector */
    for (d=1; d<naxis; d++)
      if ((++rawpos2[d])<=bufmax[d])
        break;
      else
        rawpos2[d] = bufmin[d];
    }

#ifdef USE_THREADS
  if (threadstep++)
    threads_gate_sync(pthread_stopgate2);
#endif

  free(linebuf);

/*  Check whether some data remain to be read; return RETURN_ERROR otherwise */
  for (d=naxis; --d;)
    if (rawpos2[d]>wcs->outmax[d])
      return RETURN_ERROR;
    else if (rawpos2[d]<wcs->outmax[d])
      break;

  return RETURN_OK;
  }


#ifdef USE_THREADS

/****** pthread_move_lines ***************************************************
PROTO	void *pthread_move_lines(void *arg)
PURPOSE	thread that takes care of moving image "lines" from the input to the
	co-addition buffer
INPUT	Pointer to the thread number.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	21/10/2003
 ***/
void	*pthread_move_lines(void *arg)
  {
   int	proc;

  proc = *((int *)arg);
  threads_gate_sync(pthread_startgate2);
  while (!pthread_endflag)
    {
    if (pthread_wdataflag)
      coadd_movewdata(pthread_linebuf, pthread_multibuf, pthread_multinbuf,
		pthread_npix, pthread_step);
    else
      coadd_movedata(pthread_linebuf, pthread_multibuf, pthread_multinbuf,
		pthread_npix, pthread_step);
    threads_gate_sync(pthread_stopgate2);
/*-- ( Master thread loads new data here ) */
    threads_gate_sync(pthread_startgate2);
    }

  pthread_exit(NULL);

  return (void *)NULL;
  }

#endif

/******* coadd_movedata ******************************************************
PROTO	void coadd_movedata(PIXTYPE *linebuf, PIXTYPE *multibuf,
			int *multinbuf, int npix, int step)
PURPOSE	Move data from the input load buffer to the co-addition buffer.
INPUT	Input Buffer,
	co-addition buffer,
	number-of-inputs buffer,
	number of pixels to process,
	step in pixels of the co-addition buffer.
OUTPUT	-.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 03/07/2003
 ***/
void	coadd_movedata(PIXTYPE *linebuf, PIXTYPE *multibuf,
			unsigned int *multinbuf, int npix, int step)
  {
   PIXTYPE	*pix, *multi;
   unsigned int	*multin,
		x;

  pix = linebuf;
  multi = multibuf;
  multin = multinbuf;
  for (x=npix; x--;  multi += step)
    *(multi+*(multin++)) = *(pix++);

  return;
  }


/******* coadd_movewdata *****************************************************
PROTO	void coadd_movewdata(PIXTYPE *linebuf, PIXTYPE *multiwbuf,
			int *multinbuf, int npix, int step)
PURPOSE	Move weight data from the input load buffer to the co-addition buffer.
INPUT	Input Buffer,
	co-addition weight buffer,
	number-of-inputs buffer,
	number of pixels to process,
	step in pixels of the co-addition buffer.
OUTPUT	-.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 01/08/2006
 ***/
void	coadd_movewdata(PIXTYPE *linebuf, PIXTYPE *multiwbuf,
			unsigned int *multinbuf, int npix, int step)
  {
   PIXTYPE	*wpix, *multiw,
		val;
   unsigned int	*multin,
		x, n;

  multiw = multiwbuf;
  multin = multinbuf;
  if (linebuf)
    {
/*-- Resampled weight-map present (always the case when RESAMPLE is Y) */
     wpix = linebuf;
     for (x=npix; x--;  multiw += step)
       {
       n = *multin;
       if ((val=*(wpix++)) > 0.0)
         {
         *(multiw+n) = 1.0/val;
         (*(multin++))++;
         }
       else if (!n)
         {
         *(multiw+n) = BIG;
         (*(multin++))++;
         }
       else
         multin++;
       }
     }
  else
/*-- Resampled weight-map not present: weight is 1 within frame limits */
    for (x=npix; x--; multiw += step)
      {
      *(multiw+*multin) =  1.0;
      (*(multin++))++;
      }
  return;
  }

