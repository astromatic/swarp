 /*
 				back.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SWarp
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	functions dealing with background computation.
*
*	Last modify:	10/12/2004
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _FIELD_H_
#include "field.h"
#endif

/*----------------------------- Internal constants --------------------------*/
#define	BACK_BUFSIZE		(8*1024*1024)	/* bkgnd buffer */
#define	BACK_MINGOODFRAC	0.5		/* min frac with good weights*/
#define	QUANTIF_NSIGMA		5		/* histogram limits */
#define	QUANTIF_NMAXLEVELS	4096		/* max nb of quantif. levels */
#define	QUANTIF_AMIN		4		/* min nb of "mode pixels" */

/* NOTES:
One must have:		BACK_BUFSIZE >= MAXPICSIZE
			0 < QUANTIF_NSIGMA <= 10
			QUANTIF_AMIN > 0
*/

/*------------------------------- structures --------------------------------*/
/* Background info */
typedef struct structback
  {
  float		mode, mean, sigma;	/* Background mode, mean and sigma */
  int		*histo;			/* Pointer to a histogram */
  int		nlevels;		/* Nb of histogram bins */
  float		qzero, qscale;		/* Position of histogram */
  float		lcut, hcut;		/* Histogram cuts */
  int		npix;			/* Number of pixels involved */
  }	backstruct;


/*------------------------------- functions ---------------------------------*/
extern void	backhisto(backstruct *, backstruct *, PIXTYPE *, PIXTYPE *,
			size_t, int, int, int, PIXTYPE),
		backline(fieldstruct *, int, PIXTYPE *),
		backstat(backstruct *, backstruct *, PIXTYPE *, PIXTYPE *,
			size_t, int, int, int, PIXTYPE),
		backrmsline(fieldstruct *, int, PIXTYPE *),
		end_back(fieldstruct *),
		filter_back(fieldstruct *),
		make_back(fieldstruct *, fieldstruct *);

extern float	backguess(backstruct *, float *, float *),
		*make_backspline(fieldstruct *, float *);


