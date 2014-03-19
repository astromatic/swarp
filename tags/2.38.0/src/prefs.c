/*
*				prefs.c
*
* Run-time configuration.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SWarp
*
*	Copyright:		(C) 2000-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		22/03/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include	<ctype.h>

#ifdef HAVE_MATHIMF_H
#include <mathimf.h>
#else
#include <math.h>
#endif

#ifdef HAVE_GETRLIMIT
#include 	<sys/time.h>
#include 	<sys/resource.h>
#endif

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<unistd.h>
#if defined(USE_THREADS) \
&& (defined(__APPLE__) || defined(FREEBSD) || defined(NETBSD))	/* BSD, Apple */
 #include	<sys/types.h>
 #include	<sys/sysctl.h>
#elif defined(USE_THREADS) && defined(HAVE_MPCTL)		/* HP/UX */
 #include	<sys/mpctl.h>
#endif
#include	"define.h"
#include	"globals.h"
#include	"fits/fitscat.h"
#include	"fitswcs.h"
#include 	"prefs.h"
#include	"preflist.h"

/********************************* dumpprefs ********************************/
/*
Print the default preference parameters.
*/
void    dumpprefs(int state)
  {
   char **dp;

  dp = default_prefs;
  while (**dp)
    if (**dp != '*')
      printf("%s\n",*(dp++));
    else if (state)
      printf("%s\n",*(dp++)+1);
    else
      dp++;
  return;
  }


/********************************* readprefs ********************************/
/*
Read a configuration file in ``standard'' format (see the SExtractor
documentation)
*/
void    readprefs(char *filename, char **argkey, char **argval, int narg)

  {
   FILE		*infile;
   char		*cp, str[MAXCHARL], *keyword, *value, **dp;
   int		i, ival, nkey, warn, argi, flagc, flagd, flage, flagz;
   double	dval;
#ifdef	HAVE_GETENV
   static char	value2[MAXCHARL],envname[MAXCHAR];
   char		*dolpos, *listbuf;
#endif


  if ((infile = fopen(filename,"r")) == NULL)
    {
    flage = 1;
    warning(filename, " not found, using internal defaults");
    }
  else
    flage = 0;

/*Build the keyword-list from pkeystruct-array */

  for (i=0; key[i].name[0]; i++)
    strcpy(keylist[i], key[i].name);
  keylist[i][0] = '\0';


/*Scan the configuration file*/

  listbuf = NULL;
  argi=0;
  flagc = 0;
  flagd = 1;
  dp = default_prefs;
  for (warn=0;;)
    {
    if (flagd)
      {
      if (**dp)
        {
        if (**dp=='*')
          strcpy(str, *(dp++)+1);
        else
          strcpy(str, *(dp++));
        }
      else
        flagd = 0;
      }
    if (!flagc && !flagd)
      if (flage || !fgets(str, MAXCHARL, infile))
        flagc=1;

    if (flagc)
      {
      if (argi<narg)
        {
        sprintf(str, "%s %s", argkey[argi], argval[argi]);
        argi++;
        }
      else
        break;
      }

    keyword = strtok(str, notokstr);
    if (keyword && keyword[0]!=0 && keyword[0]!=(char)'#')
      {
     if (warn>=10)
        error(EXIT_FAILURE, "*Error*: No valid keyword found in ", filename);
      nkey = findkeys(keyword, keylist, FIND_STRICT);
      if (nkey!=RETURN_ERROR)
        {
        value = strtok((char *)NULL, notokstr);
#ifdef	HAVE_GETENV
/*------ Expansion of environment variables (preceded by '$') */
        if (value && (dolpos=strchr(value, '$')))
          {
           int	nc;
           char	*valuet,*value2t, *envval;

          value2t = value2;
          valuet = value;
          while (dolpos)
            {
            while (valuet<dolpos)
              *(value2t++) = *(valuet++);	/* verbatim copy before '$' */
            if (*(++valuet) == (char)'{')
              valuet++;
            strncpy(envname, valuet, nc=strcspn(valuet,"}/:\"\'\\"));
            *(envname+nc) = (char)'\0';
            if (*(valuet+=nc) == (char)'}')
              valuet++;
            if (!(envval=getenv(envname)))
              error(EXIT_FAILURE, "Environment variable not found: ",
				envname);
            while(*envval)			/* Copy the ENV content */
              *(value2t++) = *(envval++);
            while(*valuet && *valuet!=(char)'$')/* Continue verbatim copy */
              *(value2t++) = *(valuet++);
            if (*valuet)
              dolpos = valuet;
            else
              {
              dolpos = NULL;
              *value2t = (char)'\0';
              }
	    }

          value = strtok(value2, notokstr);
          }
#endif
        switch(key[nkey].type)
          {
          case P_FLOAT:
            if (!value || value[0]==(char)'#')
              error(EXIT_FAILURE, keyword," keyword has no value!");
            if (*value=='@')
              value = listbuf = list_to_str(value+1);
            dval = atof(value);
            if (dval>=key[nkey].dmin && dval<=key[nkey].dmax)
              *(double *)(key[nkey].ptr) = dval;
            else
              error(EXIT_FAILURE, keyword," keyword out of range");
            break;

          case P_INT:
            if (!value || value[0]==(char)'#')
              error(EXIT_FAILURE, keyword," keyword has no value!");
            if (*value=='@')
              value = listbuf = list_to_str(value+1);
            ival = strtol(value, NULL, 0);
            if (ival>=key[nkey].imin && ival<=key[nkey].imax)
              *(int *)(key[nkey].ptr) = ival;
            else
              error(EXIT_FAILURE, keyword, " keyword out of range");
            break;

          case P_STRING:
            if (!value || value[0]==(char)'#')
              error(EXIT_FAILURE, keyword," string is empty!");
            if (*value=='@')
              value = listbuf = list_to_str(value+1);
            strcpy((char *)key[nkey].ptr, value);
            break;

          case P_BOOL:
            if (!value || value[0]==(char)'#')
              error(EXIT_FAILURE, keyword," keyword has no value!");
            if (*value=='@')
              value = listbuf = list_to_str(value+1);
            if ((cp = strchr("yYnN", (int)value[0])))
              *(int *)(key[nkey].ptr) = (tolower((int)*cp)=='y')?1:0;
            else
              error(EXIT_FAILURE, keyword, " value must be Y or N");
            break;

          case P_KEY:
            if (!value || value[0]==(char)'#')
              error(EXIT_FAILURE, keyword," keyword has no value!");
            if (*value=='@')
              value = listbuf = list_to_str(value+1);
            if ((ival = findkeys(value, key[nkey].keylist,FIND_STRICT))
			!= RETURN_ERROR)
              *(int *)(key[nkey].ptr) = ival;
            else
              error(EXIT_FAILURE, keyword, " set to an unknown keyword");
            break;

          case P_BOOLLIST:
            if (value && *value=='@')
              value = strtok(listbuf = list_to_str(value+1), notokstr);
            for (i=0; i<MAXLIST&&value&&value[0]!=(char)'#'; i++)
              {
              if (i>=key[nkey].nlistmax)
                error(EXIT_FAILURE, keyword, " has too many members");
              if ((cp = strchr("yYnN", (int)value[0])))
                ((int *)(key[nkey].ptr))[i] = (tolower((int)*cp)=='y')?1:0;
              else
                error(EXIT_FAILURE, keyword, " value must be Y or N");
              value = strtok((char *)NULL, notokstr);
              }
            if (i<key[nkey].nlistmin)
              error(EXIT_FAILURE, keyword, " list has not enough members");
            *(key[nkey].nlistptr) = i;
            break;

          case P_INTLIST:
            if (value && *value=='@')
              value = strtok(listbuf = list_to_str(value+1), notokstr);
            for (i=0; i<MAXLIST&&value&&value[0]!=(char)'#'; i++)
              {
              if (i>=key[nkey].nlistmax)
                error(EXIT_FAILURE, keyword, " has too many members");
              ival = strtol(value, NULL, 0);
              if (ival>=key[nkey].imin && ival<=key[nkey].imax)
                ((int *)key[nkey].ptr)[i] = ival;
              else
                error(EXIT_FAILURE, keyword, " keyword out of range");
              value = strtok((char *)NULL, notokstr);
              }
            if (i<key[nkey].nlistmin)
              error(EXIT_FAILURE, keyword, " list has not enough members");
            *(key[nkey].nlistptr) = i;
            break;

          case P_FLOATLIST:
            if (value && *value=='@')
              value = strtok(listbuf = list_to_str(value+1), notokstr);
            for (i=0; i<MAXLIST&&value&&value[0]!=(char)'#'; i++)
              {
              if (i>=key[nkey].nlistmax)
                error(EXIT_FAILURE, keyword, " has too many members");
              dval = atof(value);
              if (dval>=key[nkey].dmin && dval<=key[nkey].dmax)
                ((double *)key[nkey].ptr)[i] = dval;
              else
                error(EXIT_FAILURE, keyword, " keyword out of range");
              value = strtok((char *)NULL, notokstr);
              }
            if (i<key[nkey].nlistmin)
              error(EXIT_FAILURE, keyword, " list has not enough members");
            *(key[nkey].nlistptr) = i;
            break;

          case P_KEYLIST:
            if (value && *value=='@')
              value = strtok(listbuf = list_to_str(value+1), notokstr);
            for (i=0; i<MAXLIST && value && value[0]!=(char)'#'; i++)
              {
              if (i>=key[nkey].nlistmax)
                error(EXIT_FAILURE, keyword, " has too many members");
	      if ((ival = findkeys(value, key[nkey].keylist, FIND_STRICT))
			!= RETURN_ERROR)
                ((int *)(key[nkey].ptr))[i] = ival;
              else
                error(EXIT_FAILURE, keyword, " set to an unknown keyword");
              value = strtok((char *)NULL, notokstr);
              }
            if (i<key[nkey].nlistmin)
              error(EXIT_FAILURE, keyword, " list has not enough members");
            *(key[nkey].nlistptr) = i;
            break;

          case P_STRINGLIST:
            if (value && *value=='@')
              value = strtok(listbuf = list_to_str(value+1), notokstr);
            if (!value || value[0]==(char)'#')
              {
              value = "";
              flagz = 1;
              }
            else
              flagz = 0;
            for (i=0; i<MAXLIST && value && value[0]!=(char)'#'; i++)
              {
              if (i>=key[nkey].nlistmax)
                error(EXIT_FAILURE, keyword, " has too many members");
              free(((char **)key[nkey].ptr)[i]);
              QMALLOC(((char **)key[nkey].ptr)[i], char, MAXCHAR);
              strcpy(((char **)key[nkey].ptr)[i], value);
              value = strtok((char *)NULL, notokstr);
              }
            if (i<key[nkey].nlistmin)
              error(EXIT_FAILURE, keyword, " list has not enough members");
            *(key[nkey].nlistptr) = flagz?0:i;
            break;

          default:
            error(EXIT_FAILURE, "*Internal ERROR*: Type Unknown",
				" in readprefs()");
            break;
          }
        if (listbuf)
          {
          free(listbuf);
          listbuf = NULL;
          }
        key[nkey].flag = 1;
        }
      else
        {
        warning(keyword, " keyword unknown");
        warn++;
        }
      }
    }

  for (i=0; key[i].name[0]; i++)
    if (!key[i].flag)
      error(EXIT_FAILURE, key[i].name, " configuration keyword missing");
  if (!flage)
    fclose(infile);

  return;
  }


/********************************* findkeys **********************************/
/*
find an item within a list of keywords.
*/
int	findkeys(char *str, char keyw[][32], int mode)

  {
  int i;

  for (i=0; keyw[i][0]; i++)
    if (!cistrcmp(str, keyw[i], mode))
      return i;

  return RETURN_ERROR;
  }


/******************************* cistrcmp ***********************************/
/*
case-insensitive strcmp.
*/
int	cistrcmp(char *cs, char *ct, int mode)

  {
   int	i, diff;

  if (mode)
    {
    for (i=0; cs[i]&&ct[i]; i++)
      if ((diff=tolower((int)cs[i])-tolower((int)ct[i])))
        return diff;
    }
  else
    {
    for (i=0; cs[i]||ct[i]; i++)
      if ((diff=tolower((int)cs[i])-tolower((int)ct[i])))
        return diff;
    }

  return 0;
  }


/****** list_to_str **********************************************************
PROTO	char	*list_to_str(char *listname)
PURPOSE	Read the content of a file and convert it to a long string.
INPUT	File name.
OUTPUT	Pointer to an allocated string, or NULL if something went wrong.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	06/02/2008
 ***/
char	*list_to_str(char *listname)
  {
   FILE	*fp;
   char		liststr[MAXCHAR],
		*listbuf, *str;
   int		l, bufpos, bufsize;

  if (!(fp=fopen(listname,"r")))
    error(EXIT_FAILURE, "*Error*: File not found: ", listname);
  bufsize = 8*MAXCHAR;
  QMALLOC(listbuf, char, bufsize);
  for (bufpos=0; fgets(liststr,MAXCHAR,fp);)
    for (str=NULL; (str=strtok(str? NULL: liststr, "\n\r\t "));)
      {
      if (bufpos>MAXLISTSIZE)
        error(EXIT_FAILURE, "*Error*: Too many parameters in ", listname);
      l = strlen(str)+1;
      if (bufpos+l > bufsize)
        {
        bufsize += 8*MAXCHAR;
        QREALLOC(listbuf, char, bufsize);
        }
      if (bufpos)
        listbuf[bufpos-1] = ' ';
      strcpy(listbuf+bufpos, str);
      bufpos += l;
      }
  fclose(fp);

  return listbuf;
  }


/********************************* useprefs *********************************/
/*
Update various structures according to the prefs.
*/
void	useprefs(void)

  {
   unsigned short	ashort=1;
   char			*pstr;
   int			i,j, weight_flag;
#ifdef USE_THREADS
   int			nproc;
#endif


/* Test if byteswapping will be needed */
  bswapflag = *((char *)&ashort);

/* Multithreading */
#ifdef USE_THREADS
  if (!prefs.nthreads)
    {
/*-- Get the number of processors for parallel builds */
/*-- See, e.g. http://ndevilla.free.fr/threads */
    nproc = -1;
#if defined(_SC_NPROCESSORS_ONLN)		/* AIX, Solaris, Linux */
    nproc = (int)sysconf(_SC_NPROCESSORS_ONLN);
#elif defined(_SC_NPROCESSORS_CONF)
    nproc = (int)sysconf(_SC_NPROCESSORS_CONF);
#elif defined(__APPLE__) || defined(FREEBSD) || defined(NETBSD)	/* BSD, Apple */
    {
     int	mib[2];
     size_t	len;

     mib[0] = CTL_HW;
     mib[1] = HW_NCPU;
     len = sizeof(nproc);
     sysctl(mib, 2, &nproc, &len, NULL, 0);
     }
#elif defined (_SC_NPROC_ONLN)			/* SGI IRIX */
    nproc = sysconf(_SC_NPROC_ONLN);
#elif defined(HAVE_MPCTL)			/* HP/UX */
    nproc =  mpctl(MPC_GETNUMSPUS_SYS, 0, 0);
#endif

    if (nproc>0)
      prefs.nthreads = nproc;
    else
      {
      prefs.nthreads = 2;
      warning("Cannot find the number of CPUs on this system:",
		"NTHREADS defaulted to 2");
      }
    }
#else
  if (prefs.nthreads != 1)
    {
    prefs.nthreads = 1;
    warning("NTHREADS != 1 ignored: ",
	"this build of " BANNER " is single-threaded");
    }
#endif

#ifdef HAVE_GETRLIMIT
/* Check that the maximum number of open files is < what the system allows */
  if (prefs.combine_flag)
    {
     struct rlimit	rlim;
    if (prefs.nopenfiles_max && getrlimit(RLIMIT_NOFILE, &rlim) != -1
	&& rlim.rlim_cur< RLIM_INFINITY
	&& prefs.nopenfiles_max > (unsigned int)rlim.rlim_cur)
      {
      sprintf(gstr, "%d", (unsigned int)rlim.rlim_cur);
      warning("NOPENFILES_MAX larger than the system limit currently set at ",
		gstr);
      }
    }
#endif

/* Trigger integer mode if RESAMPLE is Y and RESAMPLING_TYPE is set to FLAGS */
  prefs.outfield_bitpix = BP_FLOAT;
  for (i=0; i<prefs.nresamp_type; i++)
    if ((prefs.resample_flag) && prefs.resamp_type[i] == INTERP_FLAGS)
      {
      prefs.outfield_bitpix = BP_LONG;
      for (j=0; j<prefs.nresamp_type; j++)
        prefs.resamp_type[j] = INTERP_FLAGS;
      break;
      }

/* Trigger integer mode if COMBINE is Y and boolean arithmetics are selected */
  if (prefs.combine_flag)
    {
    if (prefs.coadd_type == COADD_AND
	|| prefs.coadd_type == COADD_NAND
	|| prefs.coadd_type == COADD_OR
	|| prefs.coadd_type == COADD_NOR)
      {
      prefs.outfield_bitpix = BP_LONG;
      for (j=0; j<prefs.nresamp_type; j++)
        prefs.resamp_type[j] = INTERP_FLAGS;
      }
    else if (prefs.outfield_bitpix == BP_LONG)
      {
      for (j=0; j<prefs.nresamp_type; j++)
        prefs.coadd_type = COADD_OR;
      warning("COMBINE_TYPE incompatible with RESAMPLING_TYPE FLAGS: ",
	"Forcing to OR");
      }
    }   

/* Force weight flag */
  weight_flag = 0;
  for (i=0; i<prefs.nweight_type; i++)
    weight_flag |= (prefs.weight_type[i]!=WEIGHT_NONE);

/* If Weights are needed... */
  if (weight_flag)
    {
/*-- Use the WEIGHT_SUFFIX to identify the weight-maps */
    if (!prefs.ninwfield)
      {
      for (i=0; i<prefs.ninfield; i++)
        {
        QMALLOC(prefs.inwfield_name[i], char, MAXCHAR);
/*------ Create a file name with a new extension */
        strcpy(prefs.inwfield_name[i], prefs.infield_name[i]);
        if (!(pstr = strrchr(prefs.inwfield_name[i], '.')))
          pstr = prefs.inwfield_name[i]+strlen(prefs.inwfield_name[i]);
        sprintf(pstr, "%s", prefs.weight_suffix);
        }
      prefs.ninwfield = prefs.ninfield;
      }

    if (prefs.nweight_type > prefs.ninwfield)
      {
      for (i=0; i<prefs.nweight_type; i++)
        {
        if (prefs.weight_type[i] == WEIGHT_NONE
		|| prefs.weight_type[i] == WEIGHT_FROMBACK)
          {
/*-------- If the background map is internal, shift the next filenames ...*/
          for (j=prefs.ninwfield; j>i; j--)
            prefs.inwfield_name[j] = prefs.inwfield_name[j-1];
/*-------- ... and replace the current one with a dummy one */
          QMALLOC(prefs.inwfield_name[i], char, MAXCHAR);
          sprintf(prefs.inwfield_name[i], "INTERNAL");
          prefs.ninwfield++;
          }
        }      
/*---- Now check that we haven't gone too far!! */
      if (prefs.ninwfield > prefs.nweight_type)
        error(EXIT_FAILURE, "*Error*: the number of WEIGHT_TYPEs and ",
		"weight-maps do not match");
      }

    if (prefs.ninwfield != prefs.ninfield)
      {
/*---- Weight-maps given through the WEIGHT_IMAGE keyword */
      if (prefs.ninwfield == 1)
	{
	  warning("Several input images and a single weight-map found: ",
		"applying the same weight-map to all images");
	  prefs.ninwfield = prefs.ninfield;
	  for (i=1; i<prefs.ninwfield; i++)
	    {
	      QMALLOC(prefs.inwfield_name[i], char, MAXCHAR);
	      strcpy(prefs.inwfield_name[i],prefs.inwfield_name[0]);
	    }
	}
      else
        error(EXIT_FAILURE, "*Error*: the number of input images and ",
		"weight-maps do not match");
      }
/*-- Weight types */
    for (i=prefs.nweight_type; i<prefs.ninwfield; i++)
      prefs.weight_type[i] = prefs.weight_type[prefs.nweight_type-1];
    prefs.nweight_type = prefs.ninwfield;

/*-- Weight rescaling flag */
    for (i=prefs.nwscale_flag; i<prefs.ninfield; i++)
      prefs.wscale_flag[i] = prefs.wscale_flag[prefs.nwscale_flag-1];
    prefs.nwscale_flag = prefs.ninfield;
    }

/* If flags and types are less than images, copy the last one */
/* Background subtraction flag */
  for (i=prefs.nsubback_flag; i<prefs.ninfield; i++)
    prefs.subback_flag[i] = prefs.subback_flag[prefs.nsubback_flag-1];
  prefs.nsubback_flag = prefs.ninfield;
/* Interpolation flag */
  for (i=prefs.ninterp_flag; i<prefs.ninfield; i++)
    prefs.interp_flag[i] = prefs.interp_flag[prefs.ninterp_flag-1];
  prefs.ninterp_flag = prefs.ninfield;
/* Background subtraction type */
  for (i=prefs.nback_type; i<prefs.ninfield; i++)
    prefs.back_type[i] = prefs.back_type[prefs.nback_type-1];
  prefs.nback_type = prefs.ninfield;
/* Gain default value */
  for (i=prefs.ngain_default; i<prefs.ninfield; i++)
    prefs.gain_default[i] = prefs.gain_default[prefs.ngain_default-1];
  prefs.ngain_default = prefs.ninfield;
/* Saturation default value */
  for (i=prefs.nsat_default; i<prefs.ninfield; i++)
    prefs.sat_default[i] = prefs.sat_default[prefs.nsat_default-1];
  prefs.nsat_default = prefs.ninfield;
/* Flux scale default value */
  for (i=prefs.nfscale_default; i<prefs.ninfield; i++)
    prefs.fscale_default[i] = prefs.fscale_default[prefs.nfscale_default-1];
  prefs.nfscale_default = prefs.ninfield;
/* Background mesh sizes */
  for (i=prefs.nback_size; i<prefs.ninfield; i++)
    prefs.back_size[i] = prefs.back_size[prefs.nback_size-1];
  prefs.nback_size = prefs.ninfield;
/* Background filter sizes */
  for (i=prefs.nback_fsize; i<prefs.ninfield; i++)
    prefs.back_fsize[i] = prefs.back_fsize[prefs.nback_fsize-1];
  prefs.nback_fsize = prefs.ninfield;
/* Background values */
  for (i=prefs.nback_default; i<prefs.ninfield; i++)
    prefs.back_default[i] = prefs.back_default[prefs.nback_default-1];
  prefs.nback_default = prefs.ninfield;

/* Default weight-threshold values */
  if (!prefs.nweight_thresh)
    {
    prefs.weight_thresh[0] = (prefs.weight_type[0]==WEIGHT_FROMWEIGHTMAP)?
		      0.0:BIG;
    prefs.nweight_thresh = 1;
    }
  else
    for (i=prefs.nweight_thresh; i<prefs.ninfield; i++)
      prefs.weight_thresh[i] = prefs.weight_thresh[prefs.nweight_thresh-1];
  prefs.nweight_thresh = prefs.ninfield;

/* Projection approximation error */
  for (i=prefs.nproj_err; i<prefs.ninfield; i++)
    prefs.proj_err[i] = prefs.proj_err[prefs.nproj_err-1];
  prefs.nproj_err = prefs.ninfield;

/* Interpolation method */
  for (i=prefs.nresamp_type; i<INTERP_MAXDIM; i++)
    prefs.resamp_type[i] = prefs.resamp_type[prefs.nresamp_type-1];
  prefs.nresamp_type = INTERP_MAXDIM;
/* Interpolation oversampling */
  for (i=prefs.noversamp; i<INTERP_MAXDIM; i++)
    prefs.oversamp[i] = prefs.oversamp[prefs.noversamp-1];
  prefs.noversamp = INTERP_MAXDIM;
/* Centering type */
  for (i=prefs.ncenter_type; i<INTERP_MAXDIM; i++)
    prefs.center_type[i] = prefs.center_type[prefs.ncenter_type-1];
  prefs.ncenter_type = INTERP_MAXDIM;
/* Center */
  for (i=prefs.nimage_center; i<INTERP_MAXDIM; i++)
    {
    QMALLOC(prefs.image_center[i], char, MAXCHAR);
    strcpy(prefs.image_center[i], prefs.image_center[prefs.nimage_center-1]);
    }
  prefs.nimage_center = INTERP_MAXDIM;
/* Pixel scale type */
  for (i=prefs.npixscale_type; i<INTERP_MAXDIM; i++)
    prefs.pixscale_type[i] = prefs.pixscale_type[prefs.npixscale_type-1];
  prefs.npixscale_type = INTERP_MAXDIM;
/* Pixel scale */
  for (i=prefs.npixscale; i<INTERP_MAXDIM; i++)
    prefs.pixscale[i] = prefs.pixscale[prefs.npixscale-1];
  prefs.npixscale = INTERP_MAXDIM;
/* Image size */
  for (i=prefs.nimage_size; i<INTERP_MAXDIM; i++)
    prefs.image_size[i] = prefs.image_size[prefs.nimage_size-1];
  prefs.nimage_size = INTERP_MAXDIM;

/* Check projection type */
  if (wcs_supproj(prefs.projection_name))
    error(EXIT_FAILURE, "*Error*: Unsupported projection type: ",
	prefs.projection_name);

/* Set Virtual memory parameters (for alloc_body()) */
  set_swapdir(prefs.swapdir_name);
  set_maxram((size_t)prefs.mem_max*1024*1024);
  set_maxvram((size_t)prefs.vmem_max*1024*1024);

/* Copy FITS keywords: check and correct lengths */
  for (i=0; i<prefs.ncopy_keywords; i++)
    if ((j=strlen(pstr=prefs.copy_keywords[i]))!=8)
      {
      pstr[8] = '\0';
      for (pstr+=j; j<8; j++)
        *(pstr++) = ' ';
      }

  return;
  }

