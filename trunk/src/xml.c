/*
*				xml.c
*
* Handle XML metada.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SWarp
*
*	Copyright:		(C) 2000-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		07/06/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "define.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "field.h"
#include "fitswcs.h"
#include "key.h"
#include "prefs.h"
#include "xml.h"

extern time_t		thetime,thetime2;	/* from makeit.c */
extern pkeystruct	key[];			/* from preflist.h */
extern char		keylist[][32];		/* from preflist.h */
xmlstruct		xmlref;
xmlstruct		*xmlstack = NULL;
int			nxml=0, nxmlmax=0;

/****** init_xml ************************************************************
PROTO	int	init_xml(int next)
PURPOSE	Initialize a set of meta-data kept in memory before being written to the
	XML file
INPUT	Number of image extensions.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	26/07/2006
 ***/
int	init_xml(int next)
  {
  QMALLOC(xmlstack, xmlstruct, next);
  nxml = 0;
  nxmlmax = next;

  return EXIT_SUCCESS;
  }


/****** end_xml ************************************************************
PROTO	int	end_xml(void)
PURPOSE	Free the set of meta-data kept in memory.
INPUT	-.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	12/07/2006
 ***/
int	end_xml(void)
  {
  free(xmlstack);

  return EXIT_SUCCESS;
  }


/****** update_xml ***********************************************************
PROTO	int	update_xml(fieldstruct *field, fieldstruct *wfield)
PURPOSE	Update a set of meta-data kept in memory before being written to the
	XML file
INPUT	-.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
VERSION	11/04/2011
 ***/
int	update_xml(fieldstruct *field, fieldstruct *wfield)
  {
   xmlstruct	*x;
   double	pixpos[NAXIS], wcspos[NAXIS];
   int		d;

  if (nxml < nxmlmax)
    x = &xmlstack[nxml++];
  else
    x = &xmlstack[0];	/* Extra calls update the meta-data of output frame */
  x->fieldno = field->fieldno;
  x->extension = field->frameno;
  strcpy(x->ext_date, field->sdate_end);
  strcpy(x->ext_time, field->stime_end);
  x->ext_elapsed = field->time_diff;
  strcpy(x->ident, field->ident);
  x->exptime = field->exptime;
  x->backmean = field->backmean;
  x->backsig = field->backsig;
  x->sigfac = wfield? wfield->sigfac : 1.0;
  x->weight_thresh = wfield? wfield->weight_thresh : 0.0;
  x->gain = field->gain;
  x->saturation = field->saturation;
  x->fscale = field->fscale;
  x->fascale = field->fascale;
  x->pixscale = field->wcs->pixscale*DEG/ARCSEC;
  x->equinox = field->wcs->equinox;
  x->epoch = field->wcs->epoch;
  x->obsdate = field->wcs->obsdate;
  x->naxis = field->wcs->naxis;
  x->celsys = (int)(field->wcs->celsysconvflag? field->wcs->celsys : -1);
  x->headflag = field->headflag;
  for (d=0; d<field->wcs->naxis; d++)
    pixpos[d] = (field->wcs->naxisn[d]+1.0)/2.0;
  raw_to_wcs(field->wcs, pixpos, wcspos);
  for (d=0; d<field->wcs->naxis; d++)
    x->centerpos[d] = wcspos[d];

  return EXIT_SUCCESS;
  }


/****** write_xml ************************************************************
PROTO	int	write_xml(char *filename)
PURPOSE	Save meta-data to an XML file/stream.
INPUT	XML file name.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	26/07/2006
 ***/
int	write_xml(char *filename)
  {
   FILE		*file;
   int			pipe_flag;

  pipe_flag = 0;
  if (!strcmp(prefs.xml_name, "STDOUT"))
    {
    file = stdout;
    pipe_flag = 1;
    }
  else if (!(file = fopen(prefs.xml_name, "w")))
    return RETURN_ERROR;

  write_xml_header(file);
  write_xml_meta(file, (char *)NULL);

  fprintf(file, "</RESOURCE>\n");
  fprintf(file, "</VOTABLE>\n");

  if (!pipe_flag)
    fclose(file);

  return RETURN_OK;
  }


/****** write_xml_header ******************************************************
PROTO	int	write_xml_header(FILE *file)
PURPOSE	Save an XML-VOtable header to an XML file/stream
INPUT	file or stream pointer.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
VERSION	26/07/2006
 ***/
int	write_xml_header(FILE *file)
  {
   char		sysname[16],
		*filename, *rfilename;

/* A short, "relative" version of the filename */
  filename = prefs.outfield_name;
  if (!(rfilename = strrchr(filename, '/')))
    rfilename = filename;
  else
    rfilename++;

  fprintf(file, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  fprintf(file, "<?xml-stylesheet type=\"text/xsl\" href=\"%s\"?>\n",
	prefs.xsl_name);
  fprintf(file, "<VOTABLE "
	"xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
	"xsi:noNamespaceSchemaLocation="
	"\"http://www.ivoa.net/xml/VOTable/v1.1\">\n");
  fprintf(file, "<DESCRIPTION>produced by %s</DESCRIPTION>\n", BANNER);
  fprintf(file, "<!-- VOTable description at "
	"http://www.ivoa.net/Documents/latest/VOT.html -->\n");
  fprintf(file, "<RESOURCE ID=\"%s\" name=\"%s\">\n", BANNER, rfilename);
  fprintf(file, " <DESCRIPTION>Data related to %s"
	"</DESCRIPTION>\n", BANNER);
  fprintf(file, " <INFO name=\"QUERY_STATUS\" value=\"OK\" />\n");
  switch(prefs.celsys_type)
    {
    case CELSYS_PIXEL:
      sprintf(sysname, "xy");
      break;
    case CELSYS_GALACTIC:
      sprintf(sysname, "galactic");
      break;
    case CELSYS_ECLIPTIC:
      sprintf(sysname, "ecl_FK5");
      break;
    case CELSYS_SUPERGALACTIC:
      sprintf(sysname, "supergalactic");
      break;
    case CELSYS_EQUATORIAL:
    case CELSYS_NATIVE:
    default:
      sprintf(sysname, "ICRS");
      break;
    }

  fprintf(file, " <COOSYS ID=\"J2000\" equinox=\"J2000\""
	" epoch=\"J%.10g\" system=\"%s\"/>\n",
	nxml? (xmlstack[0].epoch? xmlstack[0].epoch: 2000.0) : 2000.0, sysname);

  return RETURN_OK;
  }


/****** write_xml_meta ********************************************************
PROTO	int	write_xml_meta(FILE *file, char *error)
PURPOSE	Save meta-data to an XML-VOTable file or stream
INPUT	Pointer to the output file (or stream),
	Pointer to an error msg (or NULL).
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	07/06/2011
 ***/
int	write_xml_meta(FILE *file, char *error)
  {
   xmlstruct		*x;
   struct tm		*tm;
   char			sysname[16],
			*pspath,*psuser, *pshost, *str;
   int			d,f,n, naxis;

/* Processing date and time if msg error present */
  if (error)
    {
    thetime2 = time(NULL);
    tm = localtime(&thetime2);
    sprintf(prefs.sdate_end,"%04d-%02d-%02d",
        tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
    sprintf(prefs.stime_end,"%02d:%02d:%02d",
        tm->tm_hour, tm->tm_min, tm->tm_sec);
    prefs.time_diff = difftime(thetime2, thetime);
    }

/* Username */
  psuser = pspath = pshost = NULL;
#ifdef HAVE_GETENV
  if (!(psuser=getenv("USERNAME")))	/* Cygwin,... */
    psuser = getenv("LOGNAME");		/* Linux,... */
  pspath = getenv("PWD");
  pshost = getenv("HOSTNAME");
#endif

  naxis = (nxml? xmlstack[0].naxis : 2);
  fprintf(file, " <RESOURCE ID=\"MetaData\" name=\"MetaData\">\n");
  fprintf(file, "  <DESCRIPTION>%s meta-data</DESCRIPTION>\n", BANNER);
  fprintf(file, "  <INFO name=\"QUERY_STATUS\" value=\"OK\" />\n");
  fprintf(file, "  <PARAM name=\"Software\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.title;meta.software\" value=\"%s\"/>\n",
	BANNER);
  fprintf(file, "  <PARAM name=\"Version\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.version;meta.software\" value=\"%s\"/>\n",
	MYVERSION);
  fprintf(file, "  <PARAM name=\"Soft_URL\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.ref.url;meta.software\" value=\"%s\"/>\n",
	WEBSITE);
  fprintf(file, "  <PARAM name=\"Soft_Auth\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.bib.author;meta.software\" value=\"%s\"/>\n",
	"Emmanuel Bertin");
  fprintf(file, "  <PARAM name=\"Soft_Ref\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.bib.bibcode;meta.software\" value=\"%s\"/>\n",
	"2002ASPC..281..228B");
  fprintf(file, "  <PARAM name=\"NThreads\" datatype=\"int\""
	" ucd=\"meta.number;meta.software\" value=\"%d\"/>\n",
    	prefs.nthreads);
  fprintf(file, "  <PARAM name=\"Date\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"time.event.end;meta.software\" value=\"%s\"/>\n",
	prefs.sdate_end);
  fprintf(file, "  <PARAM name=\"Time\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"time.event.end;meta.software\" value=\"%s\"/>\n",
	prefs.stime_end);
  fprintf(file, "  <PARAM name=\"Duration\" datatype=\"float\""
	" ucd=\"time.event;meta.software\" value=\"%.2f\" unit=\"s\"/>\n",
	prefs.time_diff);

  fprintf(file, "  <PARAM name=\"User\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.curation\" value=\"%s\"/>\n",
	psuser);
  fprintf(file, "  <PARAM name=\"Host\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.curation\" value=\"%s\"/>\n",
	pshost);
  fprintf(file, "  <PARAM name=\"Path\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset\" value=\"%s\"/>\n",
	pspath);

  fprintf(file,
	"  <PARAM name=\"Image_Name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"obs.image;meta.fits\" value=\"%s\"/>\n", prefs.outfield_name);

  fprintf(file,
	"  <PARAM name=\"Weight_Name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"obs.image;meta.fits\" value=\"%s\"/>\n", prefs.outwfield_name);

  if (error)
    {
    fprintf(file, "\n  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	"!!!!!!!!!!!!!!!!!!!! -->\n");
    fprintf(file, "  <!-- !!!!!!!!!!!!!!!!!!!!!! an Error occured"
	" !!!!!!!!!!!!!!!!!!!!! -->\n");
    fprintf(file, "  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	"!!!!!!!!!!!!!!!!!!!! -->\n");
    fprintf(file,"  <PARAM name=\"Error_Msg\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta\" value=\"%s\"/>\n", error);
    fprintf(file, "  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	"!!!!!!!!!!!!!!!!!!!! -->\n");
    fprintf(file, "  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	"!!!!!!!!!!!!!!!!!!!! -->\n\n");
    }
  else
    {
    fprintf(file, "  <PARAM name=\"NOverlap_Max\" datatype=\"int\""
	" ucd=\"meta.number;stat.max;obs.image\" value=\"%d\"/>\n",
	xmlstack[0].fieldno);
    fprintf(file, "  <PARAM name=\"ExpTime_Max\" datatype=\"float\""
	" ucd=\"time.expo;stat.max;obs.image\" value=\"%g\" unit=\"s\"/>\n",
	xmlstack[0].exptime);
    fprintf(file, "  <PARAM name=\"Gain_Max\" datatype=\"float\""
	" ucd=\"instr.calib;stat.max;obs.image\" value=\"%g\""
	" unit=\"photon/ADU\"/>\n",
	xmlstack[0].gain);
    }

/* Meta-data for each extension */
  fprintf(file, "  <TABLE ID=\"Input_Image_Data\" name=\"Input_Image_Data\">\n");
  fprintf(file, "   <DESCRIPTION>Data gathered by %s for every FITS"
	" input image</DESCRIPTION>\n", BANNER);
  fprintf(file, "   <!-- NFrames may be 0"
	" if an error occurred early in the processing -->\n");
  fprintf(file, "   <PARAM name=\"NFrames\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n",
	nxmlmax>0? nxmlmax-1 : 0);
  fprintf(file, "   <!-- CurrFrame may differ from NFrames"
	" if an error occurred -->\n");
  fprintf(file, "   <PARAM name=\"CurrFrame\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n",
	nxml>0? nxml-1 : 0);
  fprintf(file, "   <FIELD name=\"Frame_Index\" datatype=\"int\""
        " ucd=\"meta.record\"/>\n");
  fprintf(file, "   <FIELD name=\"Image_Name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"obs.image;meta.fits\"/>\n");
  fprintf(file, "   <FIELD name=\"Weight_Name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"obs.image;meta.fits\"/>\n");
  fprintf(file, "   <FIELD name=\"External_Header\" datatype=\"boolean\""
	" ucd=\"meta.code;obs.param\"/>\n");
  fprintf(file, "   <FIELD name=\"Image_Ident\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"meta.id;obs\"/>\n");
  fprintf(file, "   <FIELD name=\"Extension\" datatype=\"int\""
        " ucd=\"meta.record\"/>\n");
  fprintf(file, "   <FIELD name=\"Date\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.record;time.event.end\"/>\n");
  fprintf(file, "   <FIELD name=\"Time\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.record;time.event.end\"/>\n");
  fprintf(file, "   <FIELD name=\"Duration\" datatype=\"float\""
	" ucd=\"meta.record;time.event.end\" unit=\"s\"/>\n");
  fprintf(file, "   <FIELD name=\"Background_Mean\" datatype=\"float\""
	" ucd=\"instr.skyLevel;obs.image;stat.median\" unit=\"adu\"/>\n");
  fprintf(file, "   <FIELD name=\"Background_StDev\" datatype=\"float\""
	" ucd=\"stat.stdev;obs.image;stat.median\" unit=\"adu\"/>\n");
  fprintf(file, "   <FIELD name=\"Subtract_Back\" datatype=\"boolean\""
	" ucd=\"meta.code;obs.param\"/>\n");
  fprintf(file, "   <FIELD name=\"Back_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;obs.param\"/>\n");
  fprintf(file, "   <FIELD name=\"Back_Size\" datatype=\"int\""
	" ucd=\"obs.param\" unit=\"pix\"/>\n");
  fprintf(file, "   <FIELD name=\"Back_FilterSize\" datatype=\"int\""
	" ucd=\"obs.param\" unit=\"pix\"/>\n");
  fprintf(file, "   <FIELD name=\"Back_Default\" datatype=\"float\""
	" ucd=\"obs.param\" unit=\"adu\"/>\n");
  fprintf(file, "   <FIELD name=\"Weight_Type\" datatype=\"boolean\""
	" arraysize=\"*\" ucd=\"stat.weight;meta.code\"/>\n");
  fprintf(file, "   <FIELD name=\"Rescale_Weights\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"stat.weight;meta.code\"/>\n");
  fprintf(file, "   <FIELD name=\"Weight_Thresh\" datatype=\"float\""
	" ucd=\"instr.sensitivity;obs.param\" unit=\"adu\"/>\n");
  fprintf(file, "   <FIELD name=\"Weight_Scaling\" datatype=\"float\""
	" ucd=\"arith.factor;obs.image;stat.median\"/>\n");
  fprintf(file, "   <FIELD name=\"Interpolate\" datatype=\"boolean\""
	" ucd=\"meta.code;obs.param\"/>\n");
  fprintf(file, "   <FIELD name=\"Gain\" datatype=\"float\""
	" ucd=\"instr.calib;obs.image\" unit=\"photon/adu\"/>\n");
  fprintf(file, "   <FIELD name=\"Saturation\" datatype=\"float\""
	" ucd=\"instr.calib;obs.image\" unit=\"adu\"/>\n");
  fprintf(file, "   <FIELD name=\"ExpTime\" datatype=\"float\""
	" ucd=\"time.expo;obs.image\" unit=\"s\"/>\n");
  fprintf(file, "   <FIELD name=\"Photometric_Flux_Scaling\" datatype=\"float\""
	" ucd=\"phot.calib;obs.image\"/>\n");
  fprintf(file, "   <FIELD name=\"Astrometric_Flux_Scaling\" datatype=\"float\""
	" ucd=\"phot.calib;obs.image\"/>\n");
  fprintf(file, "   <FIELD name=\"Field_Coordinates\" datatype=\"double\""
	" arraysize=\"%d\" ucd=\"pos.eq;obs.image\" unit=\"%s\"/>\n",
	naxis, nxml? (xmlstack[0].celsys >=0? "deg":"pix") : "deg");
  fprintf(file, "   <FIELD name=\"Pixel_Scale\" datatype=\"float\""
	" ucd=\"instr.pixel;obs.image;stat.mean\" unit=\"arcsec\"/>\n");
  fprintf(file, "   <FIELD name=\"ObsDate\" datatype=\"double\""
	" ucd=\"time.start;obs\" unit=\"yr\"/>\n");
  fprintf(file, "   <FIELD name=\"Equinox\" datatype=\"double\""
	" ucd=\"time.equinox;obs\" unit=\"yr\"/>\n");
  fprintf(file, "   <FIELD name=\"Epoch\" datatype=\"double\""
	" ucd=\"time.epoch;obs\" unit=\"yr\"/>\n");
  fprintf(file, "   <FIELD name=\"COOSYS\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"meta.code;pos\"/>\n");
  fprintf(file, "   <DATA><TABLEDATA>\n");
  for (n=1; n<nxml; n++)
    {
    x = &xmlstack[n];
    f = x->fieldno;
    fprintf(file, "    <TR>\n"
	"     <TD>%d</TD><TD>%s</TD><TD>%s</TD><TD>%c</TD><TD>%s</TD>\n"
	"     <TD>%d</TD><TD>%s</TD><TD>%s</TD><TD>%.2f</TD>\n"
	"     <TD>%g</TD><TD>%g</TD><TD>%c</TD><TD>%s</TD>"
	"<TD>%d</TD><TD>%d</TD><TD>%g</TD>\n"
	"     <TD>%s</TD><TD>%c</TD><TD>%g</TD><TD>%g</TD><TD>%c</TD>\n"
	"     <TD>%g</TD><TD>%g</TD><TD>%g</TD><TD>%g</TD>\n"
	"     <TD>%g</TD><TD>%.10g",
	n,
	prefs.infield_name[f],
	(prefs.inwfield_name[f] && *prefs.inwfield_name[f])?
		prefs.inwfield_name[f] : "(null)",
        x->headflag? 'T' : 'F',
	(x->ident && *(x->ident)) ? x->ident : "(null)",
	x->extension,
	x->ext_date,
	x->ext_time,
	x->ext_elapsed,
	x->backmean,
	x->backsig,
        prefs.subback_flag[f]? 'T' : 'F',
    	key[findkeys("BACK_TYPE", keylist,
		FIND_STRICT)].keylist[prefs.back_type[f]],
        prefs.back_size[f],
        prefs.back_fsize[f],
        prefs.back_default[f],
    	key[findkeys("WEIGHT_TYPE", keylist,
		FIND_STRICT)].keylist[prefs.weight_type[f]],
        prefs.wscale_flag[f]? 'T' : 'F',
	x->weight_thresh,
	x->sigfac,
        prefs.interp_flag[f]? 'T' : 'F',
	x->gain,
	x->saturation,
	x->exptime,
	x->fscale,
	x->fascale,
	x->centerpos[0]);
    switch(x->celsys)
      {
      case CELSYS_PIXEL:
        sprintf(sysname, "xy");
        break;
      case CELSYS_GALACTIC:
        sprintf(sysname, "galactic");
        break;
      case CELSYS_ECLIPTIC:
        sprintf(sysname, "ecl_FK5");
        break;
      case CELSYS_SUPERGALACTIC:
        sprintf(sysname, "supergalactic");
        break;
      case CELSYS_EQUATORIAL:
      case CELSYS_NATIVE:
      default:
        sprintf(sysname, "ICRS");
        break;
      }
    for (d=1; d<naxis; d++)
      fprintf(file, " %.10g", x->centerpos[d]);
    fprintf(file,
	"</TD><TD>%g</TD><TD>%.10g</TD><TD>%.10g</TD><TD>%.10g</TD>\n"
	"     <TD>%s</TD>\n    </TR>\n",
	x->pixscale, x->obsdate, x->equinox, x->epoch, sysname);
    }
  fprintf(file, "   </TABLEDATA></DATA>\n");
  fprintf(file, "  </TABLE>\n");

/* Warnings */
  fprintf(file, "  <TABLE ID=\"Warnings\" name=\"Warnings\">\n");
  fprintf(file,
	"   <DESCRIPTION>%s warnings (limited to the last %d)</DESCRIPTION>\n",
	BANNER, WARNING_NMAX);
  fprintf(file, "   <FIELD name=\"Date\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta;time.event.end\"/>\n");
  fprintf(file, "   <FIELD name=\"Time\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta;time.event.end\"/>\n");
  fprintf(file, "   <FIELD name=\"Msg\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta\"/>\n");
  fprintf(file, "   <DATA><TABLEDATA>\n");
  for (str = warning_history(); *str; str = warning_history())
    fprintf(file, "    <TR><TD>%10.10s</TD><TD>%8.8s</TD><TD>%s</TD></TR>\n",
	str, str+11, str+22);
  fprintf(file, "   </TABLEDATA></DATA>\n");
  fprintf(file, "  </TABLE>\n");

/* Configuration file */
  fprintf(file, "  <RESOURCE ID=\"Config\" name=\"Config\">\n");
  fprintf(file, "   <DESCRIPTION>%s configuration</DESCRIPTION>\n", BANNER);
  fprintf(file,
	"   <PARAM name=\"Command_Line\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"obs.param\" value=\"%s",
	prefs.command_line[0]);
  for (n=1; n<prefs.ncommand_line; n++)
    fprintf(file, " %s", prefs.command_line[n]);
  fprintf(file, "\"/>\n");
  fprintf(file,
	"   <PARAM name=\"Prefs_Name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"obs.param;meta.file\" value=\"%s\"/>\n",
	prefs.prefs_name);

  if (!error)
    {
    fprintf(file,
	"   <PARAM name=\"ImageOut_Name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset;meta.file\" value=\"%s\"/>\n",
	prefs.outfield_name);
    fprintf(file,
	"   <PARAM name=\"WeightOut_Name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset;meta.file\" value=\"%s\"/>\n",
	prefs.outwfield_name);

    fprintf(file,
	"   <PARAM name=\"Header_Only\" datatype=\"boolean\""
	" ucd=\"meta.code;obs.param\" value=\"%c\"/>\n",
    	prefs.headeronly_flag? 'T':'F');
    fprintf(file,
	"   <PARAM name=\"Header_Suffix\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset;meta.file\" value=\"%s\"/>\n",
	prefs.head_suffix);

    fprintf(file,
	"   <PARAM name=\"Weight_Suffix\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset;meta.file\" value=\"%s\"/>\n",
	prefs.weight_suffix);
    fprintf(file,
	"   <PARAM name=\"Combine\" datatype=\"boolean\""
	" ucd=\"meta.code\" value=\"%c\"/>\n",
    	prefs.combine_flag? 'T':'F');
    fprintf(file,
	"   <PARAM name=\"Combine_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;obs.param\" value=\"%s\"/>\n",
    	key[findkeys("COMBINE_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.coadd_type]);
    fprintf(file,
	"   <PARAM name=\"Blank_BadPixels\" datatype=\"boolean\""
	" ucd=\"meta.code;obs.param\" value=\"%c\"/>\n",
    	prefs.blank_flag? 'T':'F');
    fprintf(file,
	"   <PARAM name=\"Celestial_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta;pos;obs.param\" value=\"%s\"/>\n",
    	key[findkeys("CELESTIAL_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.celsys_type]);
    fprintf(file,
	"   <PARAM name=\"Projection_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta;pos;obs.param\" value=\"%s\"/>\n",
    	prefs.projection_name);

    fprintf(file, "   <PARAM name=\"Projection_Err\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"pos.angDistance;stat.max;obs.param\""
	" value=\"%g",
	prefs.proj_err[1] != prefs.proj_err[0]? 2:1, prefs.proj_err[0]);
    if (prefs.proj_err[1] != prefs.proj_err[0])
      fprintf(file, " %g", prefs.proj_err[1]);
    fprintf(file, "\" unit=\"pix\"/>\n");

    fprintf(file,
	"   <PARAM name=\"Center_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;pos;obs.param\" value=\"%s",
	key[findkeys("CENTER_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.center_type[0]]);
    if (prefs.center_type[1] != prefs.center_type[0])
      fprintf(file, ",%s",
    	key[findkeys("CENTER_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.center_type[1]]);
    fprintf(file, "\"/>\n");

    fprintf(file,
	"   <PARAM name=\"Center\" datatype=\"double\" arraysize=\"2\""
	" ucd=\"pos;obs.param\" value=\"%15.10g",
	xmlstack[0].centerpos[0]);
    for (d=1; d<naxis; d++)
      fprintf(file, " %15.10g", xmlstack[0].centerpos[d]);
    fprintf(file, "\" unit=\"%s\"/>\n",
	nxml? (xmlstack[0].celsys >=0? "deg":"pix") : "deg");

    fprintf(file,
	"   <PARAM name=\"PixelScale_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;pos;obs.param\" value=\"%s",
	key[findkeys("PIXELSCALE_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.pixscale_type[0]]);
    if (prefs.pixscale_type[1] != prefs.pixscale_type[0])
      fprintf(file, ",%s",
    	key[findkeys("PIXELSCALE_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.pixscale_type[1]]);
    fprintf(file, "\"/>\n");

    fprintf(file, "   <PARAM name=\"Pixel_Scale\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"instr.scale;instr.pixel;obs.param\""
	" value=\"%g",
	prefs.pixscale[1] != prefs.pixscale[0]? 2:1, prefs.pixscale[0]);
    if (prefs.pixscale[1] != prefs.pixscale[0])
      fprintf(file, " %g", prefs.pixscale[1]);
    fprintf(file, "\"/>\n");

    fprintf(file, "   <PARAM name=\"Image_Size\" datatype=\"int\""
	" arraysize=\"%d\" ucd=\"instr.scale;instr.pixel;obs.param\""
	" value=\"%d",
	prefs.image_size[1] != prefs.image_size[0]? 2:1, prefs.image_size[0]);
    if (prefs.image_size[1] != prefs.image_size[0])
      fprintf(file, " %d", prefs.image_size[1]);
    fprintf(file, "\"/>\n");

    fprintf(file,
	"   <PARAM name=\"Resample\" datatype=\"boolean\""
	" ucd=\"meta.code\" value=\"%c\"/>\n",
    	prefs.resample_flag? 'T':'F');
    fprintf(file,
	"   <PARAM name=\"Resample_Dir\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset\" value=\"%s\"/>\n",
	prefs.resampdir_name);
    fprintf(file,
	"   <PARAM name=\"Resample_Suffix\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset;meta.file\" value=\"%s\"/>\n",
	prefs.resamp_suffix);

    fprintf(file,
	"   <PARAM name=\"Resampling_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;pos;obs.param\" value=\"%s",
	key[findkeys("RESAMPLING_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.resamp_type[0]]);
    if (prefs.resamp_type[1] != prefs.resamp_type[0])
      fprintf(file, ",%s",
    	key[findkeys("RESAMPLING_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.resamp_type[1]]);
    fprintf(file, "\"/>\n");

    fprintf(file, "   <PARAM name=\"Oversampling\" datatype=\"int\""
	" arraysize=\"%d\" ucd=\"arith.factor;instr.pixel;obs.param\""
	" value=\"%d",
	prefs.oversamp[1] != prefs.oversamp[0]? 2:1, prefs.oversamp[0]);
    if (prefs.oversamp[1] != prefs.oversamp[0])
      fprintf(file, " %d", prefs.oversamp[1]);
    fprintf(file, "\"/>\n");

    fprintf(file,
	"   <PARAM name=\"FScalAstro_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"phot.calib;obs.param\" value=\"%s\"/>\n",
    	key[findkeys("FSCALASTRO_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.fscalastro_type]);
    fprintf(file,
	"   <PARAM name=\"FScale_Keyword\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;phot.calib;obs.param\" value=\"%s\"/>\n",
	prefs.fscale_keyword);

    fprintf(file, "   <PARAM name=\"FScale_Default\" datatype=\"float\""
	" ucd=\"phot.calib;obs.param\" value=\"%g\"/>\n",
	prefs.fscale_default[0]);

    fprintf(file,
	"   <PARAM name=\"Gain_Keyword\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;phot.calib;obs.param\" value=\"%s\"/>\n",
	prefs.gain_keyword);

    fprintf(file, "   <PARAM name=\"Gain_Default\" datatype=\"float\""
	" ucd=\"phot.calib;obs.param\" value=\"%g\"/>\n",
	prefs.gain_default[0]);

    fprintf(file,
	"   <PARAM name=\"SatLev_Keyword\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;instr.saturation;obs.param\" value=\"%s\"/>\n",
	prefs.sat_keyword);

    fprintf(file, "   <PARAM name=\"SatLev_Default\" datatype=\"float\""
	" ucd=\"instr.saturation;obs.param\" value=\"%g\"/>\n",
	prefs.sat_default[0]);

    fprintf(file,
	"   <PARAM name=\"Subtract_Back\" datatype=\"boolean\""
	" ucd=\"meta.code\" value=\"%c\"/>\n",
    	prefs.subback_flag[0]? 'T':'F');
    fprintf(file,
	"   <PARAM name=\"Back_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;obs.param\" value=\"%s\"/>\n",
    	key[findkeys("FSCALASTRO_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.back_type[0]]);
    fprintf(file, "   <PARAM name=\"Back_Size\" datatype=\"int\""
	" ucd=\"obs.param\" value=\"%d\"/>\n",
	prefs.back_size[0]);
    fprintf(file, "   <PARAM name=\"Back_FilterSize\" datatype=\"int\""
	" ucd=\"obs.param\" value=\"%d\"/>\n",
	prefs.back_fsize[0]);
    fprintf(file, "   <PARAM name=\"Back_Default\" datatype=\"float\""
	" ucd=\"obs.param\" value=\"%g\" unit=\"adu\"/>\n",
	prefs.back_default[0]);
    fprintf(file, "   <PARAM name=\"Back_FiltThresh\" datatype=\"float\""
	" ucd=\"phot.count;arith.ratio;obs.param\" value=\"%g\"/>\n",
    	prefs.back_fthresh);

    fprintf(file,
	"   <PARAM name=\"VMem_Dir\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta\" value=\"%s\"/>\n",
	prefs.swapdir_name);
    fprintf(file,
	"   <PARAM name=\"VMem_Max\" datatype=\"int\""
	" ucd=\"meta.number;stat.max\" value=\"%d\" unit=\"Mbyte\"/>\n",
	prefs.vmem_max);
    fprintf(file,
	"   <PARAM name=\"Mem_Max\" datatype=\"int\""
	" ucd=\"meta.number;stat.max\" value=\"%d\" unit=\"Mbyte\"/>\n",
	prefs.mem_max);
    fprintf(file,
	"   <PARAM name=\"Combine_BufSize\" datatype=\"int\""
	" ucd=\"meta.number;stat.max\" value=\"%d\" unit=\"Mbyte\"/>\n",
	prefs.coaddbuf_size);

    fprintf(file,
	"   <PARAM name=\"Delete_TmpFiles\" datatype=\"boolean\""
	" ucd=\"meta.code\" value=\"%c\"/>\n",
    	prefs.resample_flag? 'T':'F');

    if (prefs.ncopy_keywords)
      {
      fprintf(file,
	"   <PARAM name=\"Copy_Keywords\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta;meta.fits\" value=\"%s",
    	prefs.copy_keywords[0]);
      for (n=1; n<prefs.ncopy_keywords; n++)
        fprintf(file, ",%s", prefs.copy_keywords[n]);
      fprintf(file, "\"/>\n");
      }

    fprintf(file,
	"   <PARAM name=\"Write_FileInfo\" datatype=\"boolean\""
	" ucd=\"meta.code\" value=\"%c\"/>\n",
    	prefs.writefileinfo_flag? 'T':'F');

    fprintf(file,
	"   <PARAM name=\"Verbose_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code\" value=\"%s\"/>\n",
    	key[findkeys("VERBOSE_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.verbose_type]);
    }

  fprintf(file, "  </RESOURCE>\n");
  fprintf(file, " </RESOURCE>\n");

  return RETURN_OK;
  }


/****** write_xmlerror ******************************************************
PROTO	int	write_xmlerror(char *error)
PURPOSE	Save meta-data to a simplified XML file in case of a catched error
INPUT	a character string.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	26/07/2006
 ***/
void	write_xmlerror(char *filename, char *error)
  {
   FILE			*file;
   int			pipe_flag;

  pipe_flag = 0;
  if (!strcmp(filename, "STDOUT"))
    {
    file = stdout;
    pipe_flag = 1;
    }
  else if (!(file = fopen(filename, "w")))
    return;

  write_xml_header(file);
  write_xml_meta(file, error);

  fprintf(file, "</RESOURCE>\n");
  fprintf(file, "</VOTABLE>\n");

  if (!pipe_flag)
    fclose(file);

  return;
  }


