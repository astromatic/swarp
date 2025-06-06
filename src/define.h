/*
*				define.h
*
* Global definitions.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SWarp
*
*	Copyright:		(C) 2002-2021 IAP/CNRS/SorbonneU
*	          		(C)	2021-2023 CFHT/CNRS
*	          		(C) 2023-2025 CEA/AIM/UParisSaclay
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
*	Last modified:		04/04/2025
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/* Check if we are using a configure script here */
#ifndef HAVE_CONFIG_H
#define		VERSION		"1.x"
#define		DATE		"2003-x-x"
#define		THREADS_NMAX	1024		/* max. number of threads */
#endif


/*------------------------ what, who, when and where ------------------------*/

#define		BANNER		"SWarp"
#define		MYVERSION	VERSION
#define		COPYRIGHT	"2002-2025 IAP/CNRS/SorbonneU/CFHT/CEA/AIM/UParisSaclay"
#define		DISCLAIMER	BANNER " comes with ABSOLUTELY NO WARRANTY\n" \
		"You may redistribute copies of " BANNER "\n" \
		"under the terms of the GNU General Public License."
#define		AUTHORS		"Emmanuel BERTIN <bertin@universite-paris-saclay.fr>"
#define		WEBSITE		"https://astromatic.net/software/swarp/"
#define		INSTITUTE	"AIM https://irfu.cea.fr/dap/"

/*----------------------------- Physical constants --------------------------*/

#ifndef PI
#define PI      	3.1415926535898
#endif

/*----------------------------- Internal constants --------------------------*/

#define		BIG		1e+30		/* a huge number */
#define		OUTPUT		stderr		/* where all msgs are sent */
#define		MAXCHAR		512		/* max. number of characters */

/*--------------------- in case of missing constants ------------------------*/

#ifndef         SEEK_SET
#define         SEEK_SET        0
#endif
#ifndef         SEEK_CUR
#define         SEEK_CUR        1
#endif

#ifndef EXIT_SUCCESS
#define 	EXIT_SUCCESS	0
#endif
#ifndef EXIT_FAILURE
#define		EXIT_FAILURE	-1
#endif

/*---------------------------- return messages ------------------------------*/

#define		RETURN_OK		0
#define		RETURN_ERROR		(-1)
#define		RETURN_FATAL_ERROR	(-2)

/*------------------------------- Other Macros ------------------------------*/

#ifdef _LARGEFILE_SOURCE
#define	FSEEKO	fseeko
#define	FTELLO	ftello
#else
#define	FSEEKO	fseek
#define	FTELLO	ftell
#endif

#define	DEXP(x)	exp(2.30258509299*(x))	/* 10^x */

#define QFREAD(ptr, size, afile, fname) \
		if (fread(ptr, (size_t)(size), (size_t)1, afile)!=1) \
		  error(EXIT_FAILURE, "*Error* while reading ", fname)

#define QFWRITE(ptr, size, afile, fname) \
		if (fwrite(ptr, (size_t)(size), (size_t)1, afile)!=1) \
		  error(EXIT_FAILURE, "*Error* while writing ", fname)

#define	QFSEEK(afile, offset, pos, fname) \
		if (FSEEKO(afile, (offset), pos)) \
		  error(EXIT_FAILURE,"*Error*: file positioning failed in ", \
			fname)

#define	QFTELL(afile, pos, fname) \
		if ((pos=FTELLO(afile))==-1) \
		  error(EXIT_FAILURE,"*Error*: file position unknown in ", \
			fname)

#define	QCALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)calloc((size_t)(nel),sizeof(typ)))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%zd elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define	QMALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%zd elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define	QMALLOC16(ptr, typ, nel) \
		{if (posix_memalign((void **)&ptr, 16, (size_t)(nel)*sizeof(typ))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%zd elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define	QREALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)realloc(ptr, (size_t)(nel)*sizeof(typ))))\
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%zd elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define QMEMCPY(ptrin, ptrout, typ, nel) \
		{if (ptrin) \
                  {if (!(ptrout = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		     { \
		     sprintf(gstr, #ptrout " (" #nel "=%zd elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		     error(EXIT_FAILURE,"Could not allocate memory for ",gstr);\
                     }; \
                   memcpy(ptrout, ptrin, (size_t)(nel)*sizeof(typ)); \
                   }; \
                 }

#define	RINT(x)	(int)(floor(x+0.5))

#define	NPRINTF		if (prefs.verbose_type == NORM) fprintf

#define	NFPRINTF(w,x)	{if (prefs.verbose_type == NORM) \
				fprintf(w, "\33[1M> %s\n\33[1A",x); \
			else if (prefs.verbose_type == FULL) \
				fprintf(w, "%s\n", x);}

#define	FPRINTF		if (prefs.verbose_type == FULL)	fprintf

#define QPRINTF		if (prefs.verbose_type != QUIET)	fprintf
