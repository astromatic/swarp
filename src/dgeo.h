/**
* @file         dgeo.h
* @brief        Include file for dgeo.c.
* @date         04/11/2020
* @copyright
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*       This file part of:      SWarp
*
*       Copyright:              (C) 2020 IAP/CNRS/SorbonneU
*
*       License:                GNU General Public License
*
*       SExtractor is free software: you can redistribute it and/or modify
*       it under the terms of the GNU General Public License as published by
*       the Free Software Foundation, either version 3 of the License, or
*       (at your option) any later version.
*       SExtractor is distributed in the hope that it will be useful,
*       but WITHOUT ANY WARRANTY; without even the implied warranty of
*       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*       GNU General Public License for more details.
*       You should have received a copy of the GNU General Public License
*       along with SExtractor. If not, see <http://www.gnu.org/licenses/>.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FITSCAT_H_
#include "fits/fitscat.h"
#endif

#ifndef _FIELD_H_
#include "field.h"
#endif

#ifndef _DGEO_H_
#define _DGEO_H_

typedef	enum {DGEO_NONE, DGEO_PIXEL}	dgeoenum; /* Diff. Geo map type */

//----------------------------- Internal constants ----------------------------
//---------------------------------- Protos -----------------------------------
extern fieldstruct	*load_dgeo(catstruct *cat, fieldstruct *reffield,
				int frameno, int fieldno, dgeoenum dgeotype);

extern void		read_dgeo(fieldstruct *dgeofield);
#endif

