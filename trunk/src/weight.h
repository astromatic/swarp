/*
 				weight.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SWarp
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for weight.c.
*
*	Last modify:	21/11/2003
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _WEIGHT_H_
#define _WEIGHT_H_

typedef	enum {WEIGHT_NONE, WEIGHT_FROMBACK, WEIGHT_FROMRMSMAP,
		WEIGHT_FROMVARMAP, WEIGHT_FROMWEIGHTMAP}
		weightenum;		/* WEIGHT_IMAGE type */

/*---------------------------------- protos --------------------------------*/

extern fieldstruct	*init_weight(char *filename, fieldstruct *reffield),
			*load_weight(catstruct *cat, fieldstruct *reffield,
				int frameno, int fieldno, weightenum wtype);

extern void		read_weight(fieldstruct *wfield),
			set_weightconv(fieldstruct *wfield),
			var_to_weight(PIXTYPE *data, int npix),
			weight_to_var(PIXTYPE *data, int npix);

#endif

