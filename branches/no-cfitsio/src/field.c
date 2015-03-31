/*
*				field.c
*
* Manage fields (input images).
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
*	Last modified:		06/06/2013
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
#include "globals.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "back.h"
#include "coadd.h"
#include "data.h"
#include "field.h"
#include "header.h"
#include "key.h"
#include "misc.h"
#include "prefs.h"
#include "wcs/wcs.h"

/****** load_field ************************************************************
PROTO	fieldstruct *load_field(catstruct *cat, int frameno, int fieldno)
PURPOSE	Initialize a field structure (in read mode)
INPUT	Cat structure,
	FITS extension number in file (0=primary)
	Field number in coaddition (for various config settings)
	Field flags,
OUTPUT	The new field pointer if OK, NULL otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	25/08/2010
 ***/
fieldstruct	*load_field(catstruct *cat, int frameno, int fieldno)

  {
   tabstruct	*tab, *intab;
   fieldstruct	*field;
   char		*pstr;
   int		i;

/* First allocate memory for the new field (and nullify pointers) */
  QCALLOC(field,fieldstruct, 1);
  field->flags = FIELD_READ;
/* Set conversion flags */
  field->cflags = (prefs.interp_flag[fieldno]?CONVERT_INTERP:0)
                |(prefs.subback_flag[fieldno]?CONVERT_BACKSUB:0);
  field->frameno = frameno;
  field->fieldno = fieldno;
  field->cat = new_cat(1);
  strcpy(field->cat->filename, cat->filename);
  strcpy (field->filename, cat->filename);
  intab = cat->tab;
  for (i=frameno; i--;)
    intab = intab->nexttab;
  copy_tab_fromptr(intab, field->cat, 0);
  tab = field->tab = field->cat->tab;
  tab->cat = field->cat;

/* A short, "relative" version of the filename */
  if (!(field->rfilename = strrchr(field->filename, '/')))
    field->rfilename = field->filename;
  else
    field->rfilename++;

/* Create a file name with a "header" extension */
  strcpy(field->hfilename, field->filename);
  if (!(pstr = strrchr(field->hfilename, '.')))
    pstr = field->hfilename+strlen(field->hfilename);
  sprintf(pstr, "%s", prefs.head_suffix);

  sprintf(gstr, "Looking for %s ...", field->rfilename);
  NFPRINTF(OUTPUT, gstr);

/* Insert additional header informations from the "header" file */
  field->headflag = !read_aschead(field->hfilename, frameno, tab);

  if (tab->naxis<1)
    error(EXIT_FAILURE, "*Error*: Zero-dimensional table in ",field->filename);

/* Force data to be at least 2D */
  if (tab->naxis<2)
    {
    tab->naxis = 2;
    QREALLOC(tab->naxisn, int, 2);
    tab->naxisn[1] = 1;
    }

/* Some defaults */
  field->fascale = 1.0;
  field->fscale = prefs.fscale_default[fieldno];
  field->gain = prefs.gain_default[fieldno];
  field->saturation = prefs.sat_default[fieldno];

/* Force input celestial system to "PIXEL" if requested by user */
  if (prefs.celsys_type == CELSYS_PIXEL)
    for (i=0; i<tab->naxis; i++)
      {
      sprintf(gstr, "CTYPE%-3d", i+1);
      fitswrite(tab->headbuf, gstr, "PIXEL", H_STRING, T_STRING);
      }

/* Read WCS information in FITS header */
  field->wcs = read_wcs(tab);

/* Read additional field-related information in FITS header */
  readfitsinfo_field(field, tab);

/* Initialize effective (flux-scaled) levels */
  field->fgain = field->gain;
  field->fsaturation = field->saturation;
  field->fbackmean = field->backmean;
  field->fbacksig = field->backsig;

/* Set field width and field height (the latter can be "virtual") */
  field->width = tab->naxisn[0];
  field->height = 1;
  for (i=1; i<tab->naxis; i++)
    field->height *= tab->naxisn[i];
  field->npix = field->width*field->height;

/*-- Background */
  field->backw = prefs.back_size[fieldno]<field->width ?
					prefs.back_size[fieldno]
					: field->width;
  field->backh = prefs.back_size[fieldno]<field->height ?
					prefs.back_size[fieldno]
					: field->height;
  field->nbackp = field->backw * field->backh;
  if ((field->nbackx = (field->width-1)/field->backw + 1) < 1)
    field->nbackx = 1;
  if ((field->nbacky = (field->height-1)/field->backh + 1) < 1)
    field->nbacky = 1;
  field->nback = field->nbackx * field->nbacky;
  field->nbackfx = field->nbackx>1 ? prefs.back_fsize[fieldno] : 1;
  field->nbackfy = field->nbacky>1 ? prefs.back_fsize[fieldno] : 1;
/* Set the back_type flag if absolute background is selected */
  field->back_type = prefs.back_type[fieldno];
  field->backdefault = prefs.back_default[fieldno];

/* Check flux scale */
  if (field->fscale == 0.0)
    {
    warning(field->filename, " has flux scale = 0: I will take 1 instead");
    field->fscale = 1.0;
    }

  return field;
  }


/****** inherit_field ********************************************************
PROTO	fieldstruct *inherit_field(char *filename, fieldstruct *reffield,
				int flags)
PURPOSE	Make a copy of a field structure.
INPUT	Reference field pointer.
	flags.
OUTPUT	The new field pointer if OK, NULL otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	24/01/2012
 ***/
fieldstruct	*inherit_field(char *filename, fieldstruct *reffield,
				int flags)
  {
   fieldstruct	*field;

/* First allocate memory for the new field (and nullify pointers) */
  QCALLOC(field, fieldstruct, 1);

/* Copy what is important and reset the remaining */
  *field = *reffield;
  field->cat = new_cat(1);
  inherit_cat(reffield->cat, field->cat);
  field->tab = field->cat->tab;
  field->tab->cat = field->cat;
  field->flags = flags;
/* We don't need this */
  field->back = NULL;
  field->dback = NULL;
  field->sigma = NULL;
  field->dsigma = NULL;
  field->ipix = NULL;
  field->pix = NULL;
  field->backline = NULL;
  field->wcs = NULL;
  field->rawmin = NULL;
  field->rawmax = NULL;
  field->reffield =reffield;

  strcpy(field->filename, filename);
  strcpy(field->cat->filename, filename);

  return field;
  }


/****** end_field ************************************************************
PROTO	void end_field(fieldstruct *field)
PURPOSE	Free a field structure.
INPUT	field structure pointer.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	16/04/2000
 ***/
void	end_field(fieldstruct *field)

  {
/* Check first that a tab structure is present */
  if (field->tab)
    {
/*-- Terminate astrometry */
    end_wcs(field->wcs);
/*-- End memory mapping (if allocated) */
    free_body(field->tab);

/*-- Close cat if not already done */
    if (field->cat)
      free_cat(&field->cat, 1);
    field->tab = NULL;
    }

  end_back(field);
  field->pix = NULL;
  field->ipix = NULL;
  free(field);

  return;
  }


/****** printinfo_field ******************************************************
PROTO	void printinfo_field(fieldstruct *field)
PURPOSE	Print info about a field
INPUT	Pointer to the field.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	08/03/2007
 ***/
void	printinfo_field(fieldstruct *field, fieldstruct *wfield)

  {
   wcsstruct		*wcs;
   char			stra[16], strd[16];
   static double	pixpos[NAXIS],
			wcspos[NAXIS];
   double		pixscale;
   int			i;
  
/* Information about the file */
  if (field->frameno)
      sprintf(gstr, "Extension #%d:", field->frameno);
    else
      *gstr ='\0';
  QPRINTF(OUTPUT, "  %s  \"%.20s\"  %s  %s  %dx%d  %d bits (%s)\n",
	gstr, *field->ident? field->ident: "no ident",
	wfield? "WEIGHTED" : "unweighted",	
	field->headflag? "EXT. HEADER" : "no ext. header",	
	field->width, field->height, field->tab->bytepix*8,
	field->tab->bitpix>0?
		(field->tab->compress_type!=COMPRESS_NONE ?
			"compressed":"integers") : "floats");

/* Astrometry */
  wcs = field->wcs;

/* Find field center */
  for (i=0; i<wcs->naxis; i++)
    pixpos[i] = (wcs->naxisn[i]+1.0)/2.0;
  raw_to_wcs(wcs, pixpos, wcspos);
  if (wcs->lat != wcs->lng)
    {
    pixscale = wcs->wcsscale[wcs->lng]*DEG;
    QPRINTF(OUTPUT, "    Center: %s %s   %.3g'x%.3g'  Scale: %.4g ''/pixel\n",
	degtosexal(wcspos[wcs->lng], stra),
	degtosexde(wcspos[wcs->lat], strd),
	wcs->naxisn[wcs->lng]*pixscale/ARCMIN,
	wcs->naxisn[wcs->lat]*pixscale/ARCMIN,
	pixscale/ARCSEC);
    }
  else if (wcs->naxis >= 2)
    {
    QPRINTF(OUTPUT,
	"    Center: %.3g,%.3g   %.3gx%.3g  Scale: %.4gx%.4g /pixel\n",
	wcspos[0],
	wcspos[1],
	wcs->naxisn[0]*wcs->wcsscale[0],
	wcs->naxisn[1]*wcs->wcsscale[1],
	wcs->wcsscale[0],
	wcs->wcsscale[1]);
    }
  else
    QPRINTF(OUTPUT, "    Center: %.3g   %.3g  Scale: %.3g /pixel\n",
	wcspos[0],
	wcs->naxisn[0]*wcs->wcsscale[0],
	wcs->wcsscale[0]);

    for (i=0; i<wcs->naxis; i++)
    if (i==wcs->lat || i==wcs->lng)
      {
      pixscale = wcs->wcsscale[i]*DEG;
      break;
      }

/* Photometry */
  QPRINTF(OUTPUT, "    Gain: %.3g e-/ADU   Flux scaling (astrom/photom): "
	"%.4g X / %.4g X\n",
	field->gain,
	field->fascale,
	field->fscale);

  return;
  }


/******* init_field ********************************************************
PROTO	fieldstruct *init_field(fieldstruct **infield, int ninput,
				char *filename)

PURPOSE	Automatically set appropriate output field parameters according to the
	prefs and a set of input fields.
INPUT	Input field ptr array,
	number of input fields,
	Filename.
OUTPUT	Pointer to the new output field.
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION 06/06/2013
 ***/
fieldstruct *init_field(fieldstruct **infield, int ninput, char *filename)
  {
   extern char		celsysname[][2][8];

   fieldstruct		*field;
   tabstruct		*tab;
   wcsstruct		*wcs;
   double		pixscale[NAXIS],
			val, epoch, obsdate;
   float		*scale;
   char			*ctype[NAXIS],
			*pstr;
   int			*axis,
			i,j,n,npstr, lng,lat, naxis, countmin0,countmax0,
			countmin,countmax;

/* First allocate memory for the new field */
  QCALLOC(field,fieldstruct, 1);
  field->flags = FIELD_WRITE;

  strcpy (field->filename, filename);
/* A short, "relative" version of the filename */
  if (!(field->rfilename = strrchr(field->filename, '/')))
    field->rfilename = field->filename;
  else
    field->rfilename++;
/* Create a file name with a "header" extension */
  strcpy(field->hfilename, filename);
  if (!(pstr = strrchr(field->hfilename, '.')))
    pstr = field->hfilename+strlen(field->hfilename);
  sprintf(pstr, "%s", prefs.head_suffix);

  field->cat = new_cat(1);
  init_cat(field->cat);
  strcpy(field->cat->filename, field->filename);
  field->tab = tab = field->cat->tab;
  tab->cat = field->cat;

/*-------------------------- Set the astrometry -----------------------------*/
/* Check that angular axes are the same and adopt them for output */

  lng = lat = -1;
  naxis = 0;				/* to avoid gcc -Wall warnings */
  for (i=0; i<ninput; i++)
    {
    if (lng == -1)
      lng = infield[i]->wcs->lng;
    if (lat == -1)
      lat = infield[i]->wcs->lat;
    if (!i)
      naxis = infield[i]->wcs->naxis;
    else if (infield[i]->wcs->naxis != naxis)
      error(EXIT_FAILURE, "*Error*: Mismatched number of axes in ",
				infield[i]->filename);
    }

/* Make sure that longitudes are before latitudes */
/*
  if (lng>lat)
    {
    i = lng;
    lng = lat;
    lat = i;
    }
*/
  QMALLOC(axis, int, ninput*naxis);

  for (i=0; i<ninput; i++)
    {
    for (n=0; n<naxis; n++)
      axis[i*naxis+n] = n;
    if (lng != -1 && infield[i]->wcs->lng != -1)
      axis[i*naxis+lng] = infield[i]->wcs->lng;
    if (lat != -1 && infield[i]->wcs->lat != -1)
      axis[i*naxis+lat] = infield[i]->wcs->lat;
    }


  QMALLOC(scale, float, ninput);
  QMALLOC(tab->naxisn, int, naxis);
/* Produce a floating point output */
  field->bitpix = tab->bitpix = prefs.outfield_bitpix;

  if (prefs.resample_flag)
    {
     double	wcsmin[NAXIS],wcsmax[NAXIS],
		wcsmin1,wcsmax1,wcsmin2,wcsmax2;

    tab->naxis = naxis;

/*-- Create a new WCS structure */
/*-- Copy the types */
    for (n=0; n<naxis; n++)
      {
      QMALLOC(ctype[n], char, 16);
      strncpy(ctype[n], infield[0]->wcs->ctype[axis[0+n]], 8);
      }

/*---- Change the Celestial system if needed */
    if (lng!=lat)
      {
      if (prefs.celsys_type == CELSYS_PIXEL)
        {
        strcpy(ctype[lng], "PIXEL");
        strcpy(ctype[lat], "PIXEL");
        }
      else if (prefs.celsys_type != CELSYS_NATIVE)
        {
        n = prefs.celsys_type - 2;	/* Drop "NATIVE" and "PIXEL" options */
        strncpy(ctype[lng], celsysname[n][0], 4);
        strncpy(ctype[lat], celsysname[n][1], 4);
        }
      }

    field->wcs = wcs = create_wcs(ctype, NULL, NULL, NULL, NULL, naxis);

/*-- Copy the units */
    for (n=0; n<naxis; n++)
      {
      if (n==lng)
        strncpy(wcs->cunit[lng], infield[0]->wcs->cunit[axis[0+lng]], 8);
      else if (n==lat)
        strncpy(wcs->cunit[lat], infield[0]->wcs->cunit[axis[0+lat]], 8);
      free(ctype[n]);
      }

    for (n=0; n<naxis; n++)
      {
      wcsmin[n] =  wcsmax[n] = 0.0;	/* to avoid gcc -Wall warnings */
      switch (prefs.center_type[n])
        {
        case CENTER_MOST:
/*-------- Find (in a dumb way) the densest part of the coadd */
          countmin0 = countmax0 = 0;
          for (j=0; j<ninput; j++)
            {
            wcsmin1 = infield[j]->wcs->wcsmin[axis[j*naxis+n]];
            wcsmax1 = infield[j]->wcs->wcsmax[axis[j*naxis+n]];
            countmin = countmax = 0;
            for (i=0; i<ninput; i++)
              {
              wcsmin2 = infield[i]->wcs->wcsmin[axis[i*naxis+n]];
              wcsmax2 = infield[i]->wcs->wcsmax[axis[i*naxis+n]];
/*------------ Test for the lower limit */
              if ((n==lng && fcmp_0_p360(wcsmin1, wcsmin2)) || wcsmin1>wcsmin2)
                countmin++;
              if ((n==lng && fcmp_0_p360(wcsmax2, wcsmin1)) || wcsmax2>wcsmin1)
                countmin++;
/*------------ Test for the upper limit */
              if ((n==lng && fcmp_0_p360(wcsmax2, wcsmax1)) || wcsmax2>wcsmax1)
                countmax++;
              if ((n==lng && fcmp_0_p360(wcsmax1, wcsmin2)) || wcsmax1>wcsmin2)
                countmax++;
              }

/*---------- Update the "most intersecting" limits */
            if (countmin>countmin0)
              {
              wcsmin[n] = wcsmin1;
              countmin0 = countmin;
              }
            if (countmax>countmax0)
              {
              wcsmax[n] = wcsmax1;
              countmax0 = countmax;
              }
            }
          wcs->crval[n] = (wcsmin[n]+wcsmax[n])/2.0;
          if (n==lng)
            wcs->crval[n] = wcsmax[n]>wcsmin[n]? fmod_0_p360(wcs->crval[n])
					: fmod_0_p360(wcs->crval[n] - 180.0);
          break;

        case CENTER_ALL:
          wcsmin[n] = infield[0]->wcs->wcsmin[axis[n]];
          wcsmax[n] = infield[0]->wcs->wcsmax[axis[n]];
          for (i=1; i<ninput; i++)
            {
            wcsmin2 = infield[i]->wcs->wcsmin[axis[i*naxis+n]];
            wcsmax2 = infield[i]->wcs->wcsmax[axis[i*naxis+n]];
/*---------- Test for the lower limit */
            if ((n==lng && fcmp_0_p360(wcsmin[n],wcsmin2)) || wcsmin[n]>wcsmin2)
              wcsmin[n] = wcsmin2;
/*---------- Test for the upper limit */
            if ((n==lng && fcmp_0_p360(wcsmax2,wcsmax[n])) || wcsmax2>wcsmax[n])
              wcsmax[n] = wcsmax2;
            }
          wcs->crval[n] = (wcsmin[n]+wcsmax[n])/2.0;
          if (n==lng)
            wcs->crval[n] = wcsmax[n]>wcsmin[n]? fmod_0_p360(wcs->crval[n])
					: fmod_0_p360(wcs->crval[n] - 180.0);
          break;

        case CENTER_MANUAL:
          wcsmin[n] = infield[0]->wcs->wcsmin[axis[n]];
          wcsmax[n] = infield[0]->wcs->wcsmax[axis[n]];
          for (i=1; i<ninput; i++)
            {
            wcsmin2 = infield[i]->wcs->wcsmin[axis[i*naxis+n]];
            wcsmax2 = infield[i]->wcs->wcsmax[axis[i*naxis+n]];
/*---------- Test for the lower limit */
            if ((n==lng && fcmp_0_p360(wcsmin[n],wcsmin2)) || wcsmin[n]>wcsmin2)
              wcsmin[n] = wcsmin2;
/*---------- Test for the upper limit */
            if ((n==lng && fcmp_0_p360(wcsmax2,wcsmax[n])) || wcsmax2>wcsmax[n])
              wcsmax[n] = wcsmax2;
            }
/*-------- Handled swapped ra, dec axes (e.g. SDSS) */
          npstr = n;
          if (lng>lat)
            {
            if (n==lng)
              npstr = lat;
            else if (n==lat)
              npstr = lng;
            }
          pstr = prefs.image_center[npstr];
  	  wcs->crval[n] = strchr(pstr, ':') ?
			  (n==lng?sextodegal(pstr):sextodegde(pstr))
			: atof(pstr);
         break;

        default:
          error(EXIT_FAILURE,
	      "*Internal Error*: Unknown area type in ", "init_field()");
          break;
        }

      pixscale[n] = 0.0;
      switch (prefs.pixscale_type[n])
        {
        case PIXSCALE_MIN:
          pixscale[n] = BIG;
          for (i=0; i<ninput; i++)
            if (infield[i]->wcs->wcsscale[axis[i*naxis+n]] < pixscale[n])
              pixscale[n] = infield[i]->wcs->wcsscale[axis[i*naxis+n]];
          break;
        case PIXSCALE_MAX:
          pixscale[n] = -BIG;
          for (i=0; i<ninput; i++)
            if (infield[i]->wcs->wcsscale[axis[i*naxis+n]] > pixscale[n])
              pixscale[n] = infield[i]->wcs->wcsscale[axis[i*naxis+n]];
          break;
        case PIXSCALE_MEDIAN:
          for (i=0; i<ninput; i++)
            scale[i] = (float)infield[i]->wcs->wcsscale[axis[i*naxis+n]];
          pixscale[n] = (double)fqmedian(scale, ninput);
          break;
        case PIXSCALE_FIT:
          if (prefs.center_type[n] == CENTER_MANUAL)
            error(EXIT_FAILURE,
	      "*Error*: cannot choose an output pixel scale if ",
		"CENTER_TYPE is in MANUAL mode");
          if (!prefs.image_size[n])
            error(EXIT_FAILURE,
	      "*Error*: cannot choose an output pixel scale if ",
		"IMAGE_SIZE is not provided");

          pixscale[n] = (wcsmax[n] - wcsmin[n]) / (double)prefs.image_size[n];
          break;

        case PIXSCALE_MANUAL:
          pixscale[n] = prefs.pixscale[n];
          if ((n==lng || n==lat) && lat!=lng)
            pixscale[n] *= (ARCSEC/DEG);
          break;

        default:
          error(EXIT_FAILURE,
	      "*Internal Error*: Unknown pixel scale type in ","init_field()");
          break;
        }

      if (!pixscale[n])
        error(EXIT_FAILURE, "*Error*: null output pixel size ", "");

/*---- Following the tradition in astronomy, CDELTx decreases with longitude */
      wcs->cd[n*(naxis+1)] = wcs->cdelt[n] = (n==lng && lng != lat)?
		-pixscale[n]
		:pixscale[n];
      }

/*-- Poorman's boundary computation */
    if (wcs->celsysconvflag)
      {
      if ((prefs.center_type[lng] != CENTER_MANUAL
	|| prefs.center_type[lat] != CENTER_MANUAL))
        eq_to_celsys(wcs,  wcs->crval);
      eq_to_celsys(wcs, wcsmin);
      eq_to_celsys(wcs, wcsmax);
      for (n=0; n<naxis; n++)
        if (wcsmax[n]<wcsmin[n])
          {
          val = wcsmin[n];
          wcsmin[n] = wcsmax[n];
          wcsmax[n] = val;
          }
      }

    for (n=0; n<naxis; n++)
      {
/*---- Compute image size. It is necessarily > 0 */
      if (prefs.image_size[n])
        {
        tab->naxisn[n] = wcs->naxisn[n] = prefs.image_size[n];
        wcs->crpix[n] = (wcs->naxisn[n]+1.0)/2.0;
        }
      else
        {
        if (prefs.center_type[n] == CENTER_MANUAL)
          {
/*-------- This is a special case where CRPIX is not NAXISN/2 */
          
          val = wcs->cdelt[n] >=0 ? wcs->crval[n] - wcsmin[n] :
		wcsmax[n] - wcs->crval[n];
          if ((n==lng || n==lat) && lat!=lng)
            {
            if (val<-180.0)
              val += 360.0;
            wcs->crpix[n] = (int)(fmod(val, 180.0)/pixscale[n]) + 1.0;
/*-------- Add a 5% margin in field size */
            tab->naxisn[n] = wcs->naxisn[n]
		= (int)(fmod(wcsmax[n]-wcsmin[n]+360.0, 360.0)*1.05/pixscale[n])
		  + 1;
            }
          else
            {
            wcs->crpix[n] = (int)(val/pixscale[n]) + 1.0;
            tab->naxisn[n] = wcs->naxisn[n]
		= (int)(fabs(wcsmax[n]-wcsmin[n])*1.05/pixscale[n]) + 1;
	    }
          }
        else
          {
/*-------- Add a 5% margin in field size */
          if (lat!=lng && (n==lng || n==lat))
            tab->naxisn[n] = wcs->naxisn[n]
		= (int)(fmod_0_p360(wcsmax[n]-wcsmin[n])*1.05/pixscale[n])
		  + 1;
          else
            tab->naxisn[n] = wcs->naxisn[n]
		= (int)(fabs(wcsmax[n]-wcsmin[n])*1.05/pixscale[n]) + 1;
          wcs->crpix[n] = (wcs->naxisn[n]+1.0)/2.0;
          }
        }
      }

/*-- The special case of longitude */

    if (lat!=lng)
      {
      if (strcmp(prefs.projection_name, "NONE"))
        {
        strcpy(wcs->ctype[lng]+5, prefs.projection_name);
        strcpy(wcs->ctype[lat]+5, prefs.projection_name);
        }
      else
        {
        strcpy(wcs->ctype[lng], "");
        strcpy(wcs->ctype[lat], "");
        }
      if (!prefs.image_size[lng])
        {
        if (prefs.center_type[lng] == CENTER_MANUAL)
          {
          val = cos((wcs->crval[lat] + wcs->cdelt[lat]
			* ((wcs->naxisn[lat]+1.0)/2.0-wcs->crpix[lat]))*DEG);
          tab->naxisn[lng] = wcs->naxisn[lng]
		= (int)((wcs->naxisn[lng]-1)*val)+1;
          wcs->crpix[lng] *= val;
          wcs->crpix[lng] = (int)(wcs->crpix[lng]+0.4999);
          }
        else
          {
          tab->naxisn[lng] = wcs->naxisn[lng]
		= (int)((wcs->naxisn[lng]-1)*cos(wcs->crval[lat]*DEG))+1;
          wcs->crpix[lng] = (wcs->naxisn[lng]+1.0)/2.0;
          }
        }
      else if (prefs.pixscale_type[lng] == PIXSCALE_FIT)
        wcs->cd[lng*(naxis+1)] = (wcs->cdelt[lng]*=cos(wcs->crval[lat]*DEG));
/*---- Make pixel scales equal in alpha and delta */
      if (prefs.pixscale[lng] == prefs.pixscale[lat])
        {
        if ((val = fabs(wcs->cdelt[lng]/wcs->cdelt[lat])) < 1.0)
          wcs->cd[lng*(naxis+1)] = (wcs->cdelt[lng] /= val);
        else
          wcs->cd[lat*(naxis+1)] = (wcs->cdelt[lat] *= val);
        }

/*---- No negative longitude!! */
      wcs->crval[lng] = fmod(wcs->crval[lng]+360.0, 360.0);
      }

/*-- Default equinox, RA-DEC-sys, longpole and latpole */
    wcs->equinox = 2000.0;
    wcs->radecsys = RDSYS_ICRS;
    wcs->longpole = 999.0;
    wcs->latpole = 999.0;
    }
  else
    {
/*--- No resampling: everything is in the input headers */
/*-----------------------------------------------*/
/*--- Take the first image header as a reference */
     double	rawcenter[NAXIS], wcscenter[NAXIS],
		crpixmin,crpixmax,crpixmin1,crpixmax1,crpixmin2,crpixmax2;
     int	naxisnmax, center_flag;

    wcs = field->wcs = copy_wcs(infield[0]->wcs);
    center_flag = 0;
    for (n=0; n<naxis; n++)
      { 
      switch (prefs.center_type[n])
        {
        case CENTER_MOST:
/*-------- Find (in a dumb way) the densest part of the coadd */
          countmin0 = countmax0 = 0;
          crpixmax = crpixmin = 0.0;	/* to avoid gcc -Wall warnings */
          for (j=0; j<ninput; j++)
            {
            crpixmax1 = infield[j]->wcs->crpix[n];
            crpixmin1 = infield[j]->wcs->crpix[n]-infield[j]->wcs->naxisn[n]+1;
            countmin = countmax = 0;
            for (i=0; i<ninput; i++)
              {
              crpixmax2 = infield[i]->wcs->crpix[n];
              crpixmin2 = infield[i]->wcs->crpix[n]
				-infield[i]->wcs->naxisn[n]+1;
/*------------ Test for the lower limit */
              if (crpixmax2>=crpixmax1)
                countmax++;
              if (crpixmin2<crpixmax1)
                countmax++;
/*------------ Test for the upper limit */
              if (crpixmin2<=crpixmin1)
                countmin++;
              if (crpixmax2>crpixmin1)
                countmin++;
              }

/*---------- Update the "most intersecting" limits */
            if (countmin>countmin0)
              {
              crpixmin = crpixmin1;
              countmin0 = countmin;
              }
            if (countmax>countmax0)
              {
              crpixmax = crpixmax1;
              countmax0 = countmax;
              }
            }
          wcs->crpix[n] = crpixmax;
          wcs->naxisn[n] = (int)(crpixmax-crpixmin + 1.01);
          if (wcs->naxisn[n]<1)
            wcs->naxisn[n] = 1;
          tab->naxisn[0] = wcs->naxisn[n];
          break;

        case CENTER_ALL:
          crpixmax = -BIG;
          crpixmin = BIG;
          for (i=0; i<ninput; i++)
            {
            crpixmax2 = infield[i]->wcs->crpix[n];
            crpixmin2 = crpixmax2 - infield[i]->wcs->naxisn[n]+1;
/*---------- Test for the lower limit */
            if (crpixmax2>crpixmax)
              crpixmax = crpixmax2;
/*---------- Test for the upper limit */
            if (crpixmin2<crpixmin)
              crpixmin = crpixmin2;
            }
          wcs->crpix[n] = crpixmax;
          naxisnmax = (int) (crpixmax - crpixmin + 1.01);
          tab->naxisn[n] = wcs->naxisn[n] = (naxisnmax>1)? naxisnmax : 1;
          break;

        case CENTER_MANUAL:
          center_flag = 1;
          break;

        default:
          error(EXIT_FAILURE,
	      "*Internal Error*: Unknown area type in ", "init_field()");
        break;
        }
/*---- Manual Image size */
      if (prefs.image_size[n])
        {
        if (!center_flag)
          wcs->crpix[n] -= (wcs->naxisn[n] - prefs.image_size[n])/2;
        tab->naxisn[n] = wcs->naxisn[n] = prefs.image_size[n];
        }
      }
    if (center_flag)
      {
      for (n=0; n<naxis; n++)
        {
        npstr = n;
        if (lng>lat)
          {
          if (n==lng)
            npstr = lat;
          else if (n==lat)
            npstr = lng;
          }
        pstr = prefs.image_center[npstr];
        wcscenter[n] = strchr(pstr, ':') ?
			  (n==lng?sextodegal(pstr):sextodegde(pstr))
			: atof(pstr);
	}
      wcs_to_raw(wcs, wcscenter, rawcenter);
      for (n=0; n<naxis; n++)
        wcs->crpix[n] += (int)(wcs->naxisn[n]/2 - rawcenter[n] +0.9999);
      }
    }

  free(axis);

/* Compute mean epoch */
  epoch = 0.0;
  obsdate = BIG;
  for (j=0; j<ninput; j++)
    {
    epoch += infield[j]->wcs->epoch;
    if (infield[j]->wcs->obsdate<obsdate)
      obsdate = infield[j]->wcs->obsdate;
    }
  field->wcs->epoch = epoch / ninput;
  field->wcs->obsdate = obsdate;
  update_head(tab);
  write_wcs(tab, wcs);
/* Insert additional header informations from the "header" file */
  if (read_aschead(field->hfilename, 0, tab))
    {
/*-- No external header: update WCS internal structures only */
    if (!wcs->wcsprm)
      QCALLOC(wcs->wcsprm, struct wcsprm, 1);
/*-- Test if the WCS is recognized and a celestial pair is found */
    wcsset(wcs->naxis,(const char(*)[9])wcs->ctype, wcs->wcsprm);

/*-- Initialize other WCS structures */
    init_wcs(wcs);
/*-- Find the range of coordinates */
    range_wcs(wcs);
    }
  else
    {
    warning("FITS header data read in ", field->hfilename);
/*-- Drop the current WCS structure and update from what's in the header */
    end_wcs(wcs);
    field->wcs = wcs = read_wcs(tab);
    field->headflag = 1;
    }

  field->width = tab->naxisn[0];
  field->height = 1;
  for (i=1; i<tab->naxis; i++)
    field->height *= tab->naxisn[i];
  field->npix = field->width*field->height;

/* Default flux scale and gain and saturation */
  field->fscale = field->fascale = 1.0;
  field->gain = 0.0;
  field->saturation = 0.0;

  free(scale);

  return field;
  }


/******* scale_field *********************************************************
PROTO	void scale_field(fieldstruct *field, fieldstruct *reffield,
			int scaleflag)
PURPOSE	Compute the flux-scaling factor for each input field.
INPUT	Field ptr,
	Reference field ptr,
	Scaling flag.
OUTPUT	-.
NOTES   If scaleflag is set, pixel values are scaled in such a way that fluxes
	are conserved (but not surface brightness), otherwise, the scaling
	only applies to the output gain.
AUTHOR  E. Bertin (IAP)
VERSION 25/08/2010
 ***/
void	scale_field(fieldstruct *field, fieldstruct *reffield, int scaleflag)
  {
   wcsstruct		*wcs;
   static double	raw[NAXIS], wcspos[NAXIS], wcspos2[NAXIS];
   double		*wcsscale,
			inscale,outscale;
   int			i, naxis, lng,lat, lngin,latin;

  wcs = reffield->wcs;
  naxis = wcs->naxis;
  lngin = lng = wcs->lng;
  latin = lat = wcs->lat;
  if (lng != -1 && field->wcs->lng != -1)
    lngin = field->wcs->lng;
  if (lat != -1 && field->wcs->lat != -1)
    latin = field->wcs->lat;
/* ``Go'' to the position where scale has been computed on input image */
  for (i=0; i<naxis; i++)
    wcspos[i] = field->wcs->wcsscalepos[i];
  wcspos[lng] = field->wcs->wcsscalepos[lngin];
  wcspos[lat] = field->wcs->wcsscalepos[latin];
  wcs_to_raw(wcs, wcspos, raw);
  wcsscale = field->wcs->wcsscale;
/* Compute scaling factors for input and output images */
  inscale = outscale = 1.0;
  for (i=0; i<naxis; i++)
    {
    if ((i==lng || i==lat) && lng!=lat)
      outscale *= sqrt(wcs_scale(wcs, raw));
    else
      {
      raw[i] += 1.0;
      raw_to_wcs(wcs, raw, wcspos2);
      outscale *= fabs(wcspos2[i] - wcspos[i]);
      raw[i] -= 1.0;
      }
    inscale *= wcsscale[i];
    }

  if (inscale!=0.0 && outscale!=0.0)
    {
    if (scaleflag)
      {
      field->fascale = (outscale/inscale);
      field->fsaturation = field->saturation*field->fascale;
      }
    else
      field->fgain = field->gain*(inscale/outscale);
    }

  return;
  }


