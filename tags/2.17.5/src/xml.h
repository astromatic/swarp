 /*
 				xml.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SWarp
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	XML logging.
*
*	Last modify:	25/06/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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
/*
#define	XSL_URL	"file:///home/bertin/sources/sex/xsl/sex.xsl"
*/
#ifndef XSL_URL
#define	XSL_URL	"."
#endif
/* Alternate XSLT file at TERAPIX: */
/* will not work with recent browsers because of security limitations */
/*
#define	XSL_URL_ALT	"http://terapix.iap.fr/cplt/xsl/sex.xsl"
*/
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
  }	xmlstruct;

/*------------------------------- functions ---------------------------------*/

extern int		end_xml(void),
			init_xml(int next),
			update_xml(fieldstruct *field, fieldstruct *wfield),
			write_xml(char *filename),
			write_xml_header(FILE *file),
			write_xml_meta(FILE *file, char *error);

extern void		write_xmlerror(char *filename, char *error);
