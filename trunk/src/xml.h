/*
*				xml.h
*
* Include file for xml.c.
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

#ifndef _FITSCAT_H_
#include "fits/fitscat.h"
#endif

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

/*----------------------------- Internal constants --------------------------*/
#ifndef XSL_URL
#define	XSL_URL	"."
#endif

/*--------------------------------- typedefs --------------------------------*/
typedef struct
  {
  int		fieldno;
  int		headflag;				/* external header? */
  int		extension;
  char 		ext_date[16],ext_time[16];		/* date and time */
  float		ext_elapsed;				/* processing time */
  char		ident[80];				/* identifiant */
  float		exptime;				/* exposure time */
  float		backmean;				/* mean background */
  float		backsig;				/* mean back stddev */
  float		sigfac;					/* mean weight scaling*/
  float		weight_thresh;				/* weight threshold */
  float		gain;					/* gain (e-/ADU) */
  float		saturation;				/* saturation (ADU) */
  float		fscale;					/* photometric scaling*/
  float		fascale;				/* astrometric scaling*/
  int		naxis;					/* number of axes */
  float		pixscale;				/* pixel scale (deg2) */
  int		celsys;					/* celestial system */
  double	centerpos[NAXIS];			/* center coordinates */
  double	equinox;				/* equinox of coords */
  double	epoch;					/* epoch of coords */
  double	obsdate;				/* observation date */
  }	xmlstruct;

/*------------------------------- functions ---------------------------------*/

extern int		end_xml(void),
			init_xml(int next),
			update_xml(fieldstruct *field, fieldstruct *wfield),
			write_xml(char *filename),
			write_xml_header(FILE *file),
			write_xml_meta(FILE *file, char *error);

extern void		write_xmlerror(char *filename, char *error);
