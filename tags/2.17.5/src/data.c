 /*
 				data.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SWarp
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Reading and conversion of data.
*
*	Last modify:	09/11/2003
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
#include "globals.h"
#include "fits/fitscat.h"
#include "back.h"
#include "data.h"
#include "field.h"

/* Conversion routine globals */
fieldstruct	*convert_field;

PIXTYPE		*convert_pix, *convert_varpix, *convert_backupbuf,
		*convert_backup,
		convert_varthresh;

long		convert_pixcount;

int		*convert_ytimeoutbuf, *convert_ytimeout,
		convert_xtimeout,convert_xtimeout0,convert_ytimeout0,
		convert_width, convert_y,
		convert_interpflag, convert_backsubflag;

/******* read_data **********************************************************
PROTO	void read_data(fieldstruct *field, fieldstruct *wfield)
PURPOSE	Read data and store them in internal format (interpolated,
	background-subtracted, and dynamic-compressed).
INPUT	Input field ptr,
	Input weight field ptr.
OUTPUT	-.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 14/06/2000
 ***/
void	read_data(fieldstruct *field, fieldstruct *wfield)
  {

/* Prepare static (global) variables */
  convert_width = field->width;
  convert_pixcount = 0;
  convert_y = 0;
  convert_field = field;

/* Prepare interpolation */
  if (wfield && (field->cflags & CONVERT_INTERP))
    {
    convert_interpflag = 1;
    convert_varpix = wfield->pix;
    convert_varthresh = wfield->weight_thresh;
    QCALLOC(convert_ytimeoutbuf, int, convert_width);
    QMALLOC(convert_backupbuf, PIXTYPE, convert_width);
    }
  else
    convert_interpflag = 0;

/* Prepare background subtraction */
  if (field->cflags & CONVERT_BACKSUB)
    convert_backsubflag = 1;
  else
    convert_backsubflag = 0;

  field->pix = alloc_body(field->tab, convert_data);
  if (!field->pix)
    error(EXIT_FAILURE, "*Error*: Not enough memory ",
	  "(either in RAM or mapped to disk)");

/* Free allocated buffers */
  if (convert_interpflag)
    {
    free(convert_ytimeoutbuf);
    free(convert_backupbuf);
    }

  return;
  }


/******* convert_data ********************************************************
PROTO	void convert_data(PIXTYPE *pix, int npix)
PURPOSE	Read data and store them in internal format (interpolated,
	background-subtracted, and dynamic-compressed).
INPUT	Input field ptr,
	Number of pixels.
OUTPUT	-.
NOTES   This routine is intended to be repeatedly called from alloc_body()
	over the SAME IMAGE (no interleaving). It makes use of many static
	variables which are kept between calls.
AUTHOR  E. Bertin (IAP)
VERSION 09/11/2003
 ***/
void	convert_data(PIXTYPE *pix, int npix)
  {
   static PIXTYPE	*convert_backdata;
   PIXTYPE		*data, *vardata;
   int			i,n;

/* Data interpolation using the weight map */
  if (convert_interpflag)
    {
    data = pix;
    n = convert_pixcount;
    vardata = convert_varpix;
    for (i=npix; i--;)
      {
/*---- New line */
      if (!(n++%convert_width))
        {
        convert_xtimeout = 0;	/* As if first pixel is already interpolated */
        convert_ytimeout = convert_ytimeoutbuf;
        convert_backup = convert_backupbuf;
        }
/*---- Check if interpolation is needed */
      if (*(vardata++)>=convert_varthresh || *data < -BIG)
        {
/*------ Check if the previous pixel was already interpolated */
        if (!convert_xtimeout)
          {
          if (*convert_ytimeout)
            {
            (*convert_ytimeout)--;
            *data = *convert_backup;
            }
          }
        else
          {
          convert_xtimeout--;
          *data = *(data-1);
          }
        }
      else
        {
        convert_xtimeout = convert_xtimeout0;
        *convert_ytimeout = convert_ytimeout0;
        }
      *(convert_backup++) = *(data++);
      }
    }

/* Background-subtraction */
  if (convert_backsubflag)
    {
    data = pix;
    n = convert_pixcount;
    for (i=npix; i--; data++)
      {
/*---- New line */
      if (!(n++%convert_width))
        {
        backline(convert_field, convert_y++, convert_field->backline);
        convert_backdata = convert_field->backline;
        }
      if (*data>-BIG)
        *data -= *convert_backdata;
      convert_backdata++;
      }
    }

  convert_pixcount += npix;

  return;
  }


