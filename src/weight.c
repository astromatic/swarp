/*
*				weight.c
*
* Manage weight-maps.
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
*	Last modified:		30/04/2013
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
#include	"back.h"
#include	"field.h"
#include	"fits/fitscat.h"
#include	"prefs.h"
#include	"weight.h"

fieldstruct	*weight_reffield;
PIXTYPE		weight_fac, weight_thresh;
long		weight_pixcount;
int		weight_type, weight_width, weight_y;

/******* load_weight *********************************************************
PROTO	fieldstruct load_weight(catstruct *cat, fieldstruct *reffield,
			int frameno, weightenum wtype)
PURPOSE	Load a weight-map field
INPUT	Cat structurep pointer,
	Reference field pointer,
	FITS extension no,
	Weight type.
OUTPUT	RETURN_OK if no error, or RETURN_ERROR in case of non-fatal error(s).
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 05/10/2010
 ***/
fieldstruct	*load_weight(catstruct *cat, fieldstruct *reffield,
			int frameno, int fieldno, weightenum wtype)

  {
   fieldstruct	*wfield;
   tabstruct	*tab, *intab;
   int		i, wflags;

  wflags = 0;	/* to avoid gcc -Wall warnings */
  switch(wtype)
    {
    case WEIGHT_NONE:
      return NULL;
    case WEIGHT_FROMBACK:
      wflags = BACKRMS_FIELD;
      break;
    case WEIGHT_FROMRMSMAP:
      wflags = RMS_FIELD;
      break;

    case WEIGHT_FROMVARMAP:
      wflags = VAR_FIELD;
      break;
    case WEIGHT_FROMWEIGHTMAP:
      wflags = WEIGHT_FIELD;
      break;

    default:
      error(EXIT_FAILURE,
	"*Internal Error*: Unknown weight-map type in ", "load_weight()");
      break;
    }

  if (wtype == WEIGHT_FROMBACK)
/*-- In this special case, one needs to CREATE a new weight-map FITS file */
    wfield = inherit_field(reffield->filename, reffield, FIELD_READ | wflags);
  else
    {
/*-- First allocate memory for the new field (and nullify pointers) */
    QCALLOC(wfield,fieldstruct, 1);
    wfield->flags = FIELD_READ|wflags;
    wfield->fieldno = fieldno;
    wfield->frameno = frameno;
    strcpy (wfield->filename, cat->filename);
/*-- A short, "relative" version of the filename */
    if (!(wfield->rfilename = strrchr(wfield->filename, '/')))
      wfield->rfilename = wfield->filename;
    else
      wfield->rfilename++;

    sprintf(gstr, "Looking for %s ...", wfield->rfilename);
    NFPRINTF(OUTPUT, gstr);

/*-- Check the image exists and read important info (image size, etc...) */
    if (frameno >= cat->ntab)
      error(EXIT_FAILURE, "*Internal Error*: FITS extension unavailable in ",
		wfield->filename);
    wfield->cat = new_cat(1);
    strcpy(wfield->cat->filename, wfield->filename);
    intab = cat->tab;
    for (i=frameno; i--;)
      intab = intab->nexttab;
    copy_tab_fromptr(intab, wfield->cat, 0);
    tab = wfield->tab = wfield->cat->tab;
    tab->cat = wfield->cat;

    if (tab->naxis<1)
      error(EXIT_FAILURE, "*Error*: Zero-dimensional table in ",
		wfield->filename);

/*-- Force data to be at least 2D */
    if (tab->naxis<2)
      {
      tab->naxis = 2;
      QREALLOC(tab->naxisn, int, 2);
      tab->naxisn[1] = 1;
      }

/*-- Set field width and field height (the latter can be "virtual") */
    wfield->width = tab->naxisn[0];
    wfield->height = 1;
    for (i=1; i<tab->naxis; i++)
      wfield->height *= tab->naxisn[i];
    wfield->npix = wfield->width*wfield->height;

    if ((wfield->width!=reffield->width)||(wfield->height!=reffield->height))
      error(EXIT_FAILURE,
		"*Error*: image and weight map have different sizes","");
    }

/*-- Background */
  wfield->backw = reffield->backw;
  wfield->backh = reffield->backh;
  wfield->nbackp = reffield->nbackp;
  wfield->nbackx = reffield->nbackx;
  wfield->nbacky = reffield->nbacky;
  wfield->nback = reffield->nback;
  wfield->nbackfx = reffield->nbackfx;
  wfield->nbackfy = reffield->nbackfy;
/* Set the back_type flag if absolute background is selected */
  wfield->back_type = reffield->back_type;
  wfield->backdefault = reffield->backdefault;

/* Default normalization factor (will be changed if necessary later) */
  wfield->sigfac = 1.0;
  set_weightconv(wfield);
  wfield->weight_thresh = prefs.weight_thresh[fieldno];
  wfield->var_thresh = wfield->weight_thresh;
  weight_to_var(&wfield->var_thresh, 1);

  return wfield;
  }


/******* init_weight *********************************************************
PROTO	fieldstruct init_weight(char *filename, fieldstruct *reffield,
			weightenum wtype)
PURPOSE	Create a weight-map field
INPUT	Weight-map filename
	Reference field pointer,
	Weight type,
OUTPUT	RETURN_OK if no error, or RETURN_ERROR in case of non-fatal error(s).
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 05/10/2010
 ***/
fieldstruct	*init_weight(char *filename, fieldstruct *reffield)

  {
   fieldstruct	*wfield;

  wfield = inherit_field(filename, reffield, WEIGHT_FIELD|reffield->flags);

/* Default normalization factor */
  wfield->sigfac = 1.0;
  set_weightconv(wfield);
  wfield->var_thresh = BIG;

  return wfield;
  }


/******* read_weight *********************************************************
PROTO	void read_weight(fieldstruct *wfield)
PURPOSE	Read weights and store them in internal format (calibrated variance).
INPUT	Input weight field ptr,
OUTPUT	-.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 03/02/2012
 ***/
void	read_weight(fieldstruct *wfield)
  {
  set_weightconv(wfield);
  if (wfield->bitpix>0)
    wfield->ipix = alloc_ibody(wfield->tab, NULL);
  else
    wfield->pix = alloc_body(wfield->tab, weight_to_var);

  return;
  }


/******* set_weightconv ******************************************************
PROTO	void set_weightconv(fieldstruct *wfield)
PURPOSE	Set current weight conversion factor and flags.
INPUT	Weight field ptr.
OUTPUT	-.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 30/04/2013
 ***/
void	set_weightconv(fieldstruct *wfield)
  {
  weight_reffield = (wfield->flags & BACKRMS_FIELD) ? wfield->reffield:NULL;
  weight_pixcount = 0;
  weight_y = 0;
  weight_width = wfield->width;
  weight_type = wfield->flags&(BACKRMS_FIELD|RMS_FIELD|VAR_FIELD|WEIGHT_FIELD);
  weight_fac = (PIXTYPE)(wfield->sigfac*wfield->sigfac);
  weight_thresh = wfield->weight_thresh;

  return;
  }


/******* weight_to_var *******************************************************
PROTO	void weight_to_var(PIXTYPE *data, int npix)
PURPOSE	Turn input weights into a variance map.
INPUT	Input data ptr,
	Number of pixels.
OUTPUT	-.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 30/04/2013
 ***/
void	weight_to_var(PIXTYPE *data, int npix)

  {
   static PIXTYPE	*weight_backdata;
   int			i, n;

  switch(weight_type)
    {
    case BACKRMS_FIELD:
      n = weight_pixcount;
      if (npix==1)
        {
        *data *= *data;
        break;
        }
      for (i=npix; i--;)
        {
/*------ New line */
        if (!(n++%weight_width))
          {
          backrmsline(weight_reffield, weight_y++, weight_reffield->backline);
          weight_backdata = weight_reffield->backline;
          }
        *(data++) = *(weight_backdata++);
        }
      weight_pixcount = n;
      break;
    case RMS_FIELD:
      for (i=npix; i--; data++)
        *data *= *data<weight_thresh? *data : BIG;
      break;
    case VAR_FIELD:
      for (i=npix; i--; data++)
        *data *= *data<weight_thresh? weight_fac : BIG;
      break;
    case WEIGHT_FIELD:
      for (i=npix; i--; data++)
        *data = *data>weight_thresh? weight_fac/(*data) : BIG;
      break;
    default:
      error(EXIT_FAILURE,
	"*Internal Error*: Unknown weight-map type in ", "weight_to_var()");
      break;
    }

  return;
  }


/******* var_to_weight *******************************************************
PROTO	void var_to_weight(PIXTYPE *data, int npix)
PURPOSE	Turn variance into a weight map.
INPUT	Input data ptr,
	Number of pixels.
OUTPUT	-.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 05/10/2010
 ***/
void	var_to_weight(PIXTYPE *data, int npix)

  {
   int			i;

  switch(weight_type)
    {
    case RMS_FIELD:
      for (i=npix; i--; data++)
        if (*data<BIG)
          *data = sqrt(*data);
      break;
    case VAR_FIELD:
      break;
    case WEIGHT_FIELD:
     for (i=npix; i--; data++)
      if (*data <BIG)
        *data = weight_fac/(*data);
      else
        *data = 0.0;
      break;
    default:
      error(EXIT_FAILURE,
	"*Internal Error*: Unknown weight-map type in ", "var_to_weight()");
      break;
    }

  return;
  }

