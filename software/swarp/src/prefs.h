 /*
 				prefs.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SWarp
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for prefs.c.
*
*	Last modify:	03/01/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _FIELD_H_
#include "field.h"
#endif

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

#ifndef _COADD_H_
#include "coadd.h"
#endif

#ifndef _INTERPOLATE_H_
#include "interpolate.h"
#endif

#ifndef _WEIGHT_H_
#include "weight.h"
#endif

#ifndef _PREFS_H_
#define _PREFS_H_

/*----------------------------- Internal constants --------------------------*/

#define         MAXCHARL	16384	/* max. nb of chars in a string list */
#define		MAXLIST		256	/* max. nb of list members */
#define		MAXLISTSIZE	2000000 /* max size of list */

/*--------------------------------- typedefs --------------------------------*/
/*------------------------------- preferences -------------------------------*/
typedef struct
  {
  char		**command_line;		/* Command line */
  int		ncommand_line;		/* nb of params */
  char		prefs_name[MAXCHAR];	/* prefs filename*/
  char		*(infield_name[MAXINFIELD]);/* Filename(s) of input images */
  int		ninfield;		/* Number of input images */
  char		*(inwfield_name[MAXINFIELD]);/* Filename(s) of input weights */
  int		ninwfield;		/* Number of input weight-maps */
  char		head_suffix[MAXCHAR];	/* Generic suffix for FITS headers */
  char		resamp_suffix[MAXCHAR];	/* Generic suffix for resampled FITS */
  char		weight_suffix[MAXCHAR];	/* Generic suffix for input weights */
  char		outfield_name[MAXCHAR];	/* Output image filename */
  char		outwfield_name[MAXCHAR];/* Output weight-map filename */
  weightenum	weight_type[MAXINFIELD];/* Weight type */
  int		nweight_type;		/* nb of params */
  double	weight_thresh[MAXINFIELD];/* Weight threshold */
  int		nweight_thresh;		/* nb of params */
  int		interp_flag[MAXINFIELD];/* Interpolation flag */
  int		ninterp_flag;		/* nb of params */
  int		subback_flag[MAXINFIELD];/* Background-subtraction flag */
  int		nsubback_flag;		/* nb of params */
  int		oversamp[INTERP_MAXDIM];/* Oversampling (per pixel/ per dim) */
  int		noversamp;		/* nb of params */
  backenum	back_type[MAXINFIELD];	/* Background subtraction type */
  int		nback_type;		/* nb of params */
  int		back_size[MAXINFIELD];	/* Background mesh size */
  int		nback_size;		/* nb of params */
  int		back_fsize[MAXINFIELD];	/* Background filter size */
  int		nback_fsize;		/* nb of params */
  double	back_default[MAXINFIELD];/* Default background in MANUAL */
  int		nback_default;		/* nb of params */
  double	back_fthresh;		/* Background filter threshold */
  char		gain_keyword[MAXCHAR];	/* FITS keyword for gain */
  double	gain_default[MAXINFIELD];/* Default gain (e-/ADU) */
  int		ngain_default;		/* nb of params */
  char		sat_keyword[MAXCHAR];	/* FITS keyword for saturation */
  double	sat_default[MAXINFIELD];/* Default saturation (ADU) */
  int		nsat_default;		/* nb of params */
  char		fscale_keyword[MAXCHAR];/* FITS keyword for flux scale */
  double	fscale_default[MAXINFIELD];/* Default flux scale */
  int		nfscale_default;		/* nb of params */
  enum {FSCALASTRO_NONE, FSCALASTRO_FIXED, FSCALASTRO_VARIABLE}
		fscalastro_type;	/* Astrometric flux-scaling type */
  interpenum	resamp_type[INTERP_MAXDIM];/* Image resampling method */
  int		nresamp_type;		/* nb of params */
  coaddenum	coadd_type;		/* Coaddition type */
  int		blank_flag;		/* Blank pixels with a weight of 0? */
/* Output image coordinates */
  char		projection_name[MAXCHAR];/* Projection WCS code */
  celsysenum	celsys_type;		/* Celestial system type */
  enum {CENTER_MANUAL, CENTER_ALL, CENTER_MOST}
		center_type[INTERP_MAXDIM];/* Centering type */
  int		ncenter_type;		/* nb of params */
  char		*(image_center[INTERP_MAXDIM]);/* Center coordinates */
  int		nimage_center;		/* nb of params */
  enum {PIXSCALE_MANUAL,PIXSCALE_MIN,PIXSCALE_MAX,PIXSCALE_MEDIAN,PIXSCALE_FIT}
		pixscale_type[INTERP_MAXDIM];/* Pixel scale type */
  int		npixscale_type;		/* nb of params */
  double	pixscale[INTERP_MAXDIM];/* Pixel scales */
  int		npixscale;		/* nb of params */
  int		image_size[INTERP_MAXDIM];/* Pixel scales */
  int		nimage_size;		/* nb of params */
  double	proj_err[MAXINFIELD];	/* max astrom approximation error */
  int		nproj_err;		/* nb of params */

/* Temporary files */
  int		removetmp_flag;		/* Remove temporary FITS files ? */
/* Virtual memory handling */
  int		mem_max;		/* Max amount of allocatable RAM */ 
  int		vmem_max;		/* Max amount of allocatable VMEM */ 
  char		swapdir_name[MAXCHAR];	/* Name of virtual mem directory */

  char		resampdir_name[MAXCHAR];/* Name of resampling directory */
  int		coaddbuf_size;		/* Amount of RAM for coadd buffer */
/* Multithreading */
  int		nthreads;		/* Number of active threads */
/* Misc */
  int		combine_flag;		/* Write coadded image? */
  int		headeronly_flag;	/* Restrict output to a header? */
  int		resample_flag;		/* Resample input images? */
  int		writefileinfo_flag;	/* Write info for each input file ? */
  char		*(copy_keywords[1024]);	/* FITS keywords to be propagated */
  int		ncopy_keywords;		/* nb of params */
  int		nnodes;			/* Number of nodes (for clusters) */  
  int		node_index;		/* Node index (for multiprocessing) */ 
  enum {QUIET, NORM, WARN, FULL}	verbose_type;	/* display type */
  int		xml_flag;		/* Write XML file? */
  char		xml_name[MAXCHAR];	/* XML file name */
  char		xsl_name[MAXCHAR];	/* XSL file name (or URL) */
  char		sdate_start[12];	/* SWarp start date */
  char		stime_start[12];	/* SWarp start time */
  char		sdate_end[12];		/* SWarp end date */
  char		stime_end[12];		/* SWarp end time */
  double	time_diff;		/* Execution time */
  }	prefstruct;

prefstruct	prefs;

/*-------------------------------- protos -----------------------------------*/
extern char	*list_to_str(char *listname);

extern int	cistrcmp(char *cs, char *ct, int mode);

extern void	dumpprefs(int state),
		readprefs(char *filename,char **argkey,char **argval,int narg),
		useprefs(void);

#endif

