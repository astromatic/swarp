 /*
 				data.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SWarp
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	include for data.c.
*
*	Last modify:	16/04/2000
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _FIELD_H_
#include "field.h"
#endif

/*--------------------------- conversion flags ------------------------------*/
#define		CONVERT_INTERP	0x01	/* Interpolation */
#define		CONVERT_BACKSUB	0x02	/* Background subtraction */
#define		CONVERT_DYNCOMP	0x04	/* Dynamic compression */


/*------------------------------- functions ---------------------------------*/

extern void		convert_data(PIXTYPE *pix, int npix),
			read_data(fieldstruct *field, fieldstruct *wfield);

