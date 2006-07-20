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
*	Last modify:	20/07/2006
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
  int		frameno;
  int		extension;
  char 		ext_date[16],ext_time[16];		/* date and time */
  float		ext_elapsed;				/* processing time */
  char		imagename[MAXCHAR];			/* image filename*/
  char		weightname[MAXCHAR];			/* weight filename */
  char		ident[80];				/* identifiant */
  float		backmean;				/* mean background */
  float		backsig;				/* mean back stddev */
  float		sigfac;					/* mean weight scaling*/
  float		gain;					/* gain (e-/ADU) */
  float		fscale;					/* photometric scaling*/
  float		fascale;				/* astrometric scaling*/
  float		pixscale;				/* pixel scale (deg2) */
  double	centerpos[NAXIS]			/* center coordinates */
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
