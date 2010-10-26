/*
*				weight.h
*
* Include file for weight.c.
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

