 /*
 				header.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SWarp
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for header.c.
*
*	Last modify:	21/04/2003
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _FITSCAT_H_
#include "fits/fitscat.h"
#endif

#ifndef _FIELD_H_
#include "field.h"
#endif

/*------------------------------- functions ---------------------------------*/
extern int		read_aschead(char *filename, int frameno,
					tabstruct *tab);

extern void		readfitsinfo_field(fieldstruct *field, tabstruct *tab),
			writefitsinfo_field(fieldstruct *field,
					fieldstruct *infield),
			writefitsinfo_outfield(fieldstruct *field,
					fieldstruct *infield);
