/*
*				makeit.c
*
* Main program.
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
#include <time.h>

#include "define.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "back.h"
#include "data.h"
#include "field.h"
#include "header.h"
#include "misc.h"
#include "prefs.h"
#include "resample.h"
#include "xml.h"

#define	NFIELD	128	/* Increment in the number of fields */

static int	selectext(char *filename);
time_t		thetime, thetime2;

/********************************** makeit ***********************************/
void	makeit(void)
  {
   fieldstruct		**infield,**inwfield, *outfield,*outwfield;
   catstruct		*cat, *wcat;
   tabstruct		*tab;
   struct tm		*tm;
   double		dtime, dtimef;
   char			*rfilename;
   int		       	*next;
   int			i,j,k,l, ninfield, ntinfield,ntinfield2,
			nfield,	jima,jweight, version;

/* Install error logging */
  error_installfunc(write_error);

/* Processing start date and time */
  thetime = time(NULL);
  tm = localtime(&thetime);
  dtime = counter_seconds();
  sprintf(prefs.sdate_start,"%04d-%02d-%02d",
        tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
  sprintf(prefs.stime_start,"%02d:%02d:%02d",
        tm->tm_hour, tm->tm_min, tm->tm_sec);

  NFPRINTF(OUTPUT, "");
  QPRINTF(OUTPUT,
        "----- %s %s started on %s at %s with %d thread%s\n\n",
                BANNER,
                MYVERSION,
                prefs.sdate_start,
                prefs.stime_start,
                prefs.nthreads,
                prefs.nthreads>1? "s":"");

/* Install the signal-catching routines for temporary file cleanup */
#ifdef USE_THREADS
  install_cleanup(cancel_resample_threads);
#else
  install_cleanup(NULL);
#endif
/* Load input images */
  ninfield = prefs.ninfield;
/* First check input files and count FITS extensions when available */
  ntinfield = k = 0;
  nfield = NFIELD;
  QCALLOC(next, int, ninfield);
  QMALLOC(infield, fieldstruct *, nfield);
  QMALLOC(inwfield, fieldstruct *, nfield);
  NFPRINTF(OUTPUT, "Examining input data ...")
  for (i=0; i<ninfield; i++)
    {
/*-- Test if the filename contains a bracket indicating a particular extension*/
    jima = selectext(prefs.infield_name[i]);
    if (!(cat=read_cat(prefs.infield_name[i])))
      {
      sprintf(gstr, "*Error*: %s not found", prefs.infield_name[i]);
      error(EXIT_FAILURE, gstr,"");
      }
    if (jima >= cat->ntab)
        error(EXIT_FAILURE, "Not enough valid FITS image extensions in ",
		prefs.infield_name[i]);
/*-- Examine all extensions */
    wcat = NULL;
    jweight= RETURN_ERROR;		/* to avoid gcc -Wall warnings */
    if (prefs.weight_type[i] && prefs.weight_type[i] != WEIGHT_FROMBACK)
      {
      jweight = selectext(prefs.inwfield_name[i]);
      if (!(wcat=read_cat(prefs.inwfield_name[i])))
        {
        sprintf(gstr, "*Error*: %s not found", prefs.inwfield_name[i]);
        error(EXIT_FAILURE, gstr,"");
        }
      if (jweight >= wcat->ntab)
        error(EXIT_FAILURE, "Not enough valid FITS image extensions in ",
		prefs.inwfield_name[i]);
      }
    tab=cat->tab;
    for (j=0; j<cat->ntab; j++,tab=tab->nexttab)
      {
      if ((jima>=0 && j!=jima)
	|| (jima<0 && (!tab->naxis || (tab->tfields && tab->bitpix==8))))
        continue;
      if (k >= nfield)
        {
        nfield += NFIELD;
        QREALLOC(infield, fieldstruct *, nfield);
        QREALLOC(inwfield,fieldstruct *, nfield);
        }
      infield[k] = load_field(cat, j, i);
      inwfield[k] = load_weight(wcat, infield[k], jweight<0? j:jweight, i,
				prefs.weight_type[i]);
      next[i]++;
      k++;
      }
/*-- Put version to reduced filenames (to avoid duplicated resamps later) */
    if (k)
      {
      version = 1;
      rfilename = infield[k-1]->rfilename;
      for (l=0; l<ntinfield; l++)
        {
/*------ Check only the 1st valid image extension matching current filename */
        if ((infield[l]->frameno==0 || infield[l]->frameno==1)
		&& !strcmp(infield[l]->rfilename, rfilename))
          version++;
        }
      if (version>1)
        {
        ntinfield2 = ntinfield + next[i];
        for (l=ntinfield; l<ntinfield2; l++)
          infield[l]->version = version;
        }
      }
    ntinfield += next[i];
    if (!next[i])
      warning("No suitable data found in ", cat->filename);
    free_cat(&cat, 1);
    if (wcat)
      free_cat(&wcat, 1);
    }

/* Initialize the XML stack */
  if (prefs.xml_flag)
    init_xml(ntinfield+1);

/* Create output image (but nothing written to disk yet) */
  outwfield = NULL;
  NFPRINTF(OUTPUT, "Creating NEW output image ...")
  outfield = init_field(infield, ntinfield, prefs.outfield_name);
  NFPRINTF(OUTPUT, "Creating NEW weight-map ...")
  outwfield = init_weight(prefs.outwfield_name, outfield);
  NFPRINTF(OUTPUT, "")
  QPRINTF(OUTPUT, "------- Output File %s:\n", outfield->rfilename);
  printinfo_field(outfield, outwfield);
  QPRINTF(OUTPUT, "\n");

/* The first field in the XML stack is the output field */
  if (prefs.xml_flag)
    update_xml(outfield, outwfield);

/* HEADER_ONLY option: write the output FITS header and exit */
  if (prefs.headeronly_flag)
    {
/*-- Open output file and save header */
    if (open_cat(outfield->cat, WRITE_ONLY) != RETURN_OK)
      error(EXIT_FAILURE, "*Error*: cannot open for writing ",
		outfield->filename);
/*-- Add relevant information to FITS header */
    if (prefs.node_index==0)
      {
      writefitsinfo_outfield(outfield, *infield);
      QFWRITE(outfield->tab->headbuf, outfield->tab->headnblock*FBSIZE,
	outfield->cat->file, outfield->filename);
      }
    goto the_end;
    }

/* Read and transform the data */
  NFPRINTF(OUTPUT, "Loading input data ...")
  k = 0;
  for (i=0; i<ninfield; i++)
    {
/*-- Processing start date and time */
    for (j=0; j<next[i]; j++, k++)
      {
      dtimef = counter_seconds();
/*---- Display some info */
      if (!j)
        {
        NFPRINTF(OUTPUT, "")
        QPRINTF(OUTPUT, "-------------- File %s:\n", infield[k]->rfilename);
        }
/*---- Compute projected limits and scaling in output frame */
      if (prefs.resample_flag)
        {
        frame_wcs(infield[k]->wcs, outfield->wcs);
        scale_field(infield[k],outfield,prefs.fscalastro_type!=FSCALASTRO_NONE);
        }

      printinfo_field(infield[k], inwfield[k]);

      if (prefs.resample_flag)
        {
/*------ Open input files */
        if (open_cat(infield[k]->cat, READ_ONLY) != RETURN_OK)
          error(EXIT_FAILURE, "*Error*: Cannot re-open ", infield[k]->filename);
        if (inwfield[k])
          {
          if (open_cat(inwfield[k]->cat, READ_ONLY) != RETURN_OK)
            error(EXIT_FAILURE, "*Error*: Cannot re-open ",
			inwfield[k]->filename);
          }
/*------ Pre-compute the background map */
        if (prefs.outfield_bitpix<0)
          {
          FPRINTF(OUTPUT, "\n");
          make_back(infield[k], inwfield[k], prefs.wscale_flag[i]);
          FPRINTF(OUTPUT, "\n");
          }
        }
      if (inwfield[k])
        sprintf(gstr, "   Weight scale: %.7g", inwfield[k]->sigfac);
      else
        *gstr = '\0';
      QPRINTF(OUTPUT, "    Background: %.7g   RMS: %.7g%s\n\n",
		infield[k]->backmean, infield[k]->backsig, gstr);

      if (prefs.resample_flag)
        {
/*------ Read (and convert) the weight data */
        if (inwfield[k])
          {
          sprintf(gstr, "Reading %s ...", inwfield[k]->filename);
          NFPRINTF(OUTPUT, gstr)
          read_weight(inwfield[k]);
          }
/*------ Read (and convert) the data */
        sprintf(gstr, "Reading %s", infield[k]->filename);
        NFPRINTF(OUTPUT, gstr)
        read_data(infield[k], inwfield[k], prefs.outfield_bitpix);
/*------ Resample the data (no need to close catalogs) */
        sprintf(gstr, "Resampling %s ...", infield[k]->filename);
        NFPRINTF(OUTPUT, gstr)
        resample_field(&infield[k], &inwfield[k], outfield, outwfield,
		prefs.resamp_type);
        }
      thetime2 = time(NULL);
      tm = localtime(&thetime2);
      sprintf(infield[k]->sdate_end,"%04d-%02d-%02d",
		tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
      sprintf(infield[k]->stime_end,"%02d:%02d:%02d",
		tm->tm_hour, tm->tm_min, tm->tm_sec);
      infield[k]->time_diff = counter_seconds() - dtimef;
      if (prefs.xml_flag)
        update_xml(infield[k], inwfield[k]);
      }
    }

  if (!prefs.combine_flag)
    goto the_end;

/* Apply flux scaling to input images */
  for (k=0; k<ntinfield; k++)
    {
    if (prefs.coadd_type==COADD_WEIGHTED_WEIGHT
	|| prefs.coadd_type==COADD_MEDIAN_WEIGHT)
      {
      infield[k]->cat->tab->bscale /= (infield[k]->fscale*infield[k]->fscale);
      infield[k]->fbackmean /= (infield[k]->fscale*infield[k]->fscale);
      infield[k]->fbacksig /= (infield[k]->fscale*infield[k]->fscale);
      infield[k]->fgain *= (infield[k]->fscale*infield[k]->fscale);
      infield[k]->fsaturation /= (infield[k]->fscale*infield[k]->fscale);
      }
    else
      {
      infield[k]->cat->tab->bscale *= infield[k]->fscale;
      infield[k]->fbackmean *= infield[k]->fscale;
      infield[k]->fbacksig *= infield[k]->fscale;
      infield[k]->fgain /= infield[k]->fscale;
      infield[k]->fsaturation *= infield[k]->fscale;
      }
    if (inwfield[k])
      inwfield[k]->cat->tab->bscale /= (infield[k]->fscale*infield[k]->fscale);
    }

/* Go! */
  coadd_fields(infield, inwfield, ntinfield, outfield, outwfield,
		prefs.coadd_type, BIG);

the_end:
/* Update the output field meta-data */
  if (prefs.xml_flag)
    update_xml(outfield, outwfield);

/* Close files and free memory */
  NFPRINTF(OUTPUT, "Closing files ...")
  for (k=0; k<ntinfield; k++)
    {
    end_field(infield[k]);
    if (inwfield[k])
      end_field(inwfield[k]);
    }
  free(next);
  free(infield);
  free(inwfield);

  end_field(outfield);
  end_field(outwfield);
  cleanup_files();

/* Processing end date and time */
  thetime2 = time(NULL);
  tm = localtime(&thetime2);
  sprintf(prefs.sdate_end,"%04d-%02d-%02d",
        tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
  sprintf(prefs.stime_end,"%02d:%02d:%02d",
        tm->tm_hour, tm->tm_min, tm->tm_sec);
  prefs.time_diff = counter_seconds() - dtime;

/* Write XML */
  if (prefs.xml_flag)
    {
    write_xml(prefs.xml_name);
    end_xml();
    }

  return;
  }


/****** selectext ************************************************************
PROTO 	int selectext(char *filename)
PURPOSE	Return the user-selected extension number [%d] from the file name.
INPUT	Filename character string.
OUTPUT	Extension number, or RETURN_ERROR if nos extension specified.
NOTES	The bracket and its extension number are removed from the filename if
	found.
AUTHOR	E. Bertin (IAP)
VERSION	09/10/2007
 ***/
static int	selectext(char *filename)
  {
   char	*bracl,*bracr;
   int	next;

  if (filename && (bracl=strrchr(filename, '[')))
    {
    *bracl = '\0';
    if ((bracr=strrchr(bracl+1, ']')))
      *bracr = '\0';
    next = strtol(bracl+1, NULL, 0);
    return next;
    }

  return RETURN_ERROR;
  }


/****** write_error ********************************************************
PROTO	void	write_error(char *msg1, char *msg2)
PURPOSE	Manage files in case of a catched error
INPUT	a character string,
	another character string
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	26/07/2006
 ***/
void    write_error(char *msg1, char *msg2)
  {
   char		error[MAXCHAR];

  sprintf(error, "%s%s", msg1,msg2);
  if (prefs.xml_flag)
    write_xmlerror(prefs.xml_name, error);
  end_xml();

  return;
  }


