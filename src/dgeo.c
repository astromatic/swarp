/**
* @file         dgeo.c
* @brief        Manage differential geometry maps (to correct pixel positions)
* @date         12/02/2015
* @copyright
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*       This file part of:      SWarp
*
*       Copyright:              (C) 2020 IAP/CNRS/SorbonneU
*
*       License:                GNU General Public License
*
*       SWarp is free software: you can redistribute it and/or modify
*       it under the terms of the GNU General Public License as published by
*       the Free Software Foundation, either version 3 of the License, or
*       (at your option) any later version.
*       SWarp is distributed in the hope that it will be useful,
*       but WITHOUT ANY WARRANTY; without even the implied warranty of
*       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*       GNU General Public License for more details.
*       You should have received a copy of the GNU General Public License
*       along with SWarp. If not, see <http://www.gnu.org/licenses/>.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"dgeo.h"
#include	"field.h"
#include	"fits/fitscat.h"
#include	"prefs.h"

/******* load_dgeo ***********************************************************
PROTO	fieldstruct load_dgeo(catstruct *dcat, fieldstruct *reffield,
			int frameno, int fieldno, weightenum dgeotype)
PURPOSE	Load a differential geometry map field
INPUT	Cat structure pointer,
	Reference field pointer,
	FITS extension no,
	dgeo type.
OUTPUT	RETURN_OK if no error, or RETURN_ERROR in case of non-fatal error(s).
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 04/11/2020
 ***/
fieldstruct *load_dgeo(catstruct *cat, fieldstruct *reffield,
			int frameno, int fieldno, dgeoenum dgeotype)
  {
   fieldstruct	*dgeofield;
   tabstruct	*tab, *intab;
   int		i, dgeoflags;

  dgeoflags = 0;	/* to avoid gcc -Wall warnings */
  switch(dgeotype)
    {
    case DGEO_NONE:
      return NULL;
    case DGEO_PIXEL:
      dgeoflags = DGEO_FIELD;
      break;

    default:
      error(EXIT_FAILURE,
	"*Internal Error*: Unknown differential geometry map type in ",
		"load_dgeo()");
      break;
    }

/* First allocate memory for the new field (and nullify pointers) */
  QCALLOC(dgeofield, fieldstruct, 1);
  dgeofield->flags = FIELD_READ | dgeoflags;
  dgeofield->fieldno = fieldno;
  dgeofield->frameno = frameno;
  strcpy (dgeofield->filename, cat->filename);
/* A short, "relative" version of the filename */
  if (!(dgeofield->rfilename = strrchr(dgeofield->filename, '/')))
      dgeofield->rfilename = dgeofield->filename;
    else
      dgeofield->rfilename++;

  sprintf(gstr, "Looking for %s ...", dgeofield->rfilename);
  NFPRINTF(OUTPUT, gstr);

/* Check that the image exists and read important info (image size, etc...) */
  if (frameno >= cat->ntab)
    error(EXIT_FAILURE, "*Internal Error*: FITS extension unavailable in ",
		dgeofield->filename);
  dgeofield->cat = new_cat(1);
  strcpy(dgeofield->cat->filename, dgeofield->filename);
  intab = cat->tab;
  for (i=frameno; i--;)
    intab = intab->nexttab;
  copy_tab_fromptr(intab, dgeofield->cat, 0);
  tab = dgeofield->tab = dgeofield->cat->tab;
  tab->cat = dgeofield->cat;

  if (tab->naxis != 3)
    error(EXIT_FAILURE, "*Error*: Differential geometry map should be 3D in ",
		dgeofield->filename);

  if (tab->naxisn[2] != 2)
    error(EXIT_FAILURE, "*Error*: dgeo map should have 2 components per pixel ",
		dgeofield->filename);
/* Set field width and field height (the latter can be "virtual") */
  dgeofield->width = tab->naxisn[0];
  dgeofield->height = tab->naxisn[1];
  dgeofield->npix = dgeofield->width * dgeofield->height * 2;

  if ((dgeofield->width != reffield->width)
    	|| (dgeofield->height != reffield->height))
    error(EXIT_FAILURE,
		"*Error*: image and dgeo map have different sizes","");

  return dgeofield;
  }


/******* read_dgeo *********************************************************
PROTO	void read_dgeo(fieldstruct *dgeofield)
PURPOSE	Read differential geometry map.
INPUT	Input dgeo field ptr.
OUTPUT	-.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 04/11/2020
 ***/
void	read_dgeo(fieldstruct *dgeofield)
  {
  dgeofield->pix = alloc_body(dgeofield->tab, NULL);

  return;
  }


