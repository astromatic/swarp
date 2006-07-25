 /*
				xml.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SWarp
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	XML logging.
*
*	Last modify:	24/07/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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

extern time_t		thetimet,thetimet2;	/* from makeit.c */
extern pkeystruct	key[];			/* from preflist.h */
extern char		keylist[][32];		/* from preflist.h */
xmlstruct		xmlref;
xmlstruct		*xmlstack = NULL;
int			nxml=0, nxmlmax=0;

/****** init_xml ************************************************************
PROTO	int	init_xml(void)
PURPOSE	Initialize a set of meta-data kept in memory before being written to the
	XML file
INPUT	Number of image extensions.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	03/07/2006
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
VERSION	24/07/2006
 ***/
int	update_xml(fieldstruct *field, fieldstruct *wfield)
  {
   xmlstruct	*x;
   double	pixpos[NAXIS], wcspos[NAXIS];
   int		d;

  if (nxml >= nxmlmax)
    error(EXIT_FAILURE,"*Internal Error*: too many extensions in XML stack","");
  x = &xmlstack[nxml++];
  x->fieldno = field->fieldno;
  strcpy(x->imagename, field->filename); 
  strcpy(x->weightname, wfield? wfield->filename : "(null)"); 
  x->extension = field->frameno;
  strcpy(x->ext_date, field->sdate_end);
  strcpy(x->ext_time, field->stime_end);
  x->ext_elapsed = field->time_diff;
  strcpy(x->ident, field->ident); 
  x->backmean = field->backmean;
  x->backsig = field->backsig;
  x->sigfac = field->sigfac;
  x->gain = field->gain;
  x->fscale = field->fscale;
  x->fascale = field->fascale;
  x->pixscale = field->pixscale;
  x->epoch = field->epoch;
  x->celsys = (int)(field->wcs->celsysconvflag? field->wcs->wcscelsys : -1);
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
VERSION	14/07/2006
 ***/
int	write_xml(char *filename)
  {
   FILE		*file;

  if (!(file = fopen(prefs.xml_name, "w")))
    return RETURN_ERROR;

  write_xml_header(file);
  write_vo_fields(file);

  fprintf(file, "   <DATA>\n");
  if (prefs.cat_type == FITS_LDAC || prefs.cat_type == FITS_TPX
	|| prefs.cat_type == FITS_10)
    fprintf(file,
	"   <FITS extnum=\"%d\"><STREAM href=\"%s%s\" /> </FITS>",
	prefs.cat_type == FITS_10? 1:2,
	prefs.cat_name[0] == '/'? "file://" : "file:",
	prefs.cat_name);
  fprintf(file, "   </DATA>\n");
  fprintf(file, "  </TABLE>\n");

  write_xml_meta(file, (char *)NULL);

  fprintf(file, "</RESOURCE>\n");
  fprintf(file, "</VOTABLE>\n");

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
VERSION	14/07/2006
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
  fprintf(file, " <DESCRIPTION>Catalog of sources extracted with %s"
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
	" epoch=\"J2000\" system=\"%s\"/>\n", sysname);
  fprintf(file, " <TABLE ID=\"Source_List\" name=\"%s/out\">\n", rfilename);
  fprintf(file,
	"  <DESCRIPTION>Table of sources detected in image</DESCRIPTION>\n");
  fprintf(file,
	"  <!-- Now comes the definition of each %s parameter -->\n", BANNER);

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
VERSION	25/07/2006
 ***/
int	write_xml_meta(FILE *file, char *error)
  {
   char			*pspath,*psuser, *pshost, *str;
   struct tm		*tm;
   int			d,n, naxis;

/* Processing date and time if msg error present */
  if (error)
    {
    thetimet2 = time(NULL);
    tm = localtime(&thetimet2);
    sprintf(prefs.sdate_end,"%04d-%02d-%02d",
        tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
    sprintf(prefs.stime_end,"%02d:%02d:%02d",
        tm->tm_hour, tm->tm_min, tm->tm_sec);
    prefs.time_diff = difftime(thetimet2, thetimet);
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
	" ucd=\"time.event;meta.software\" value=\"%.0f\" unit=\"s\"/>\n",
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

/* Meta-data for each extension */
  fprintf(file, "  <TABLE ID=\"Input_Image_Data\" name=\"Input_Image_Data\">\n");
  fprintf(file, "   <DESCRIPTION>Data gathered by %s for every FITS"
	" input image</DESCRIPTION>\n", BANNER);
  fprintf(file, "   <!-- NFrames may be 0"
	" if an error occurred early in the processing -->\n");
  fprintf(file, "   <PARAM name=\"NFrames\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n",
	nxmlmax);
  fprintf(file, "   <!-- CurrFrame may differ from NFrames"
	" if an error occurred -->\n");
  fprintf(file, "   <PARAM name=\"CurrFrame\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n",
	nxml);
  fprintf(file, "   <FIELD name=\"Frame_Index\" datatype=\"int\""
        " ucd=\"meta.record\"/>\n");
  fprintf(file, "   <FIELD name=\"Image_Name\" datatype=\"*\""
	" ucd=\"obs.image;meta.fits\"/>\n");
  fprintf(file, "   <FIELD name=\"Weight_Name\" datatype=\"*\""
	" ucd=\"obs.image;meta.fits\"/>\n");
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
  fprintf(file, "   <FIELD name=\"Weight_Type\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"stat.weight;meta.code\"/>\n");
  fprintf(file, "   <FIELD name=\"Weight_Thresh\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"instr.sensitivity;obs.param\" unit=\"adu\"/>\n");
  fprintf(file, "   <FIELD name=\"Weight_Scaling\" datatype=\"float\""
	" ucd=\"arith.factor;obs.image;stat.median\"/>\n");
  fprintf(file, "   <FIELD name=\"Gain\" datatype=\"float\""
	" ucd=\"instr.calib;obs.image\" unit=\"photon/adu\"/>\n");
  fprintf(file, "   <FIELD name=\"Photometric_Flux_Scaling\" datatype=\"float\""
	" ucd=\"phot.calib;obs.image\"/>\n");
  fprintf(file, "   <FIELD name=\"Astrometric_Flux_Scaling\" datatype=\"float\""
	" ucd=\"phot.calib;obs.image\"/>\n");
  fprintf(file, "   <FIELD name=\"Field_Coordinates\" datatype=\"double\""
	" arraysize=\"%d\" ucd=\"phot.eq;obs.image\" unit=\"%s\"/>\n",
	naxis, nxml? (xmlstack[0].celsys >=0? "deg":"pix") : "deg");
  fprintf(file, "   <FIELD name=\"Pixel_Scale\" datatype=\"float\""
	" ucd=\"instr.scale;obs.image;stat.mean\" unit=\"arcsec\"/>\n");
  fprintf(file, "   <FIELD name=\"Epoch\" datatype=\"double\""
	" ucd=\"time.epoch;obs\" unit=\"yr\"/>\n",
	prefs.nimage_name);
  fprintf(file, "   <DATA><TABLEDATA>\n");
  for (n=0; n<nxml; n++)
    {
    x = &xmlstack[n];
    f = x->frameno;
    fprintf(file, "    <TR>\n"
	"     <TD>%d</TD><TD>%s</TD><TD>%s</TD><TD>%s</TD>\n"
	"      <TD>%d</TD><TD>%s</TD><TD>%s</TD><TD>%.0f</TD>\n"
	"      <TD>%g</TD><TD>%g</TD><TD>%c</TD><TD>%s</TD>"
	"<TD>%d</TD><TD>%d</TD><TD>%g</TD>\n"
	"      <TD>%s</TD><TD>%g</TD><TD>%c</TD>\n"
	"      <TD>%g</TD><TD>%g</TD><TD>%g</TD>\n      ",
	n+1,
	x->image_name,
	x->weight_name,
	x->ident,
	x->extension,
	x->ext_date,
	x->ext_time,
	x->ext_elapsed,
	x->backmean,
	x->backsig,
        prefs.subback_flag[f]? 'T' : 'F',
    	key[findkeys("BACK_TYPE",keylist,
		FIND_STRICT)].keylist[prefs.back_type[f]]),
        prefs.back_size[f],
        prefs.back_fsize[f];
        prefs.back_default[f];
    	key[findkeys("WEIGHT_TYPE", keylist,
		FIND_STRICT)].keylist[prefs.weight_type[f]),
	prefs.weight_thresh[f],
	x->sigfac,
        prefs.interp_type[f],
	x->gain,
	x->fscale,
	x->fascale);
    for (d=0; d<naxis; d++)
      fprintf(file, "<TD>%15.10g</TD>", x->centerpos[d]);
    fprintf(file, "<TD>%g</TD><TD>%g</TD>\n",
	x->pixscale, x->epoch);
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
	prefs.header_suffix);

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
			FIND_STRICT)].keylist[prefs.center_type[0]]););
    if (prefs.center_type[1] != prefs.center_type[0])
      fprintf(file, ",%s",
    	key[findkeys("CENTER_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.center_type[1]]);
    fprintf(file, "\"/>\n");

    fprintf(file,
	"   <PARAM name=\"Center\" datatype=\"double\" arraysize=\"2\""
	" ucd=\"pos;obs.param\" value=\"%g %g\" unit=\"%s\"/>\n",
	strchr(prefs.image_center[0], ':') ?
		sextodegal(prefs.image_center[0]) : atof(prefs.image_center[0]),
	strchr(prefs.image_center[1], ':') ?
		sextodegde(prefs.image_center[1]) : atof(prefs.image_center[1]));

    fprintf(file,
	"   <PARAM name=\"PixelScale_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;pos;obs.param\" value=\"%s",
	key[findkeys("PIXELSCALE_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.pixscale_type[0]]););
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
	" value=\"%g",
	prefs.image_size[1] != prefs.image_size[0]? 2:1, prefs.image_size[0]);
    if (prefs.image_size[1] != prefs.image_size[0])
      fprintf(file, " %g", prefs.image_size[1]);
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
			FIND_STRICT)].keylist[prefs.resamp_type[0]]););
    if (prefs.resamp_type[1] != prefs.resamp_type[0])
      fprintf(file, ",%s",
    	key[findkeys("RESAMPLING_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.resamp_type[1]]);
    fprintf(file, "\"/>\n");

    fprintf(file, "   <PARAM name=\"Oversampling\" datatype=\"int\""
	" arraysize=\"%d\" ucd=\"arith.factor;instr.pixel;obs.param\""
	" value=\"%g",
	prefs.oversamp[1] != prefs.oversamp[0]? 2:1, prefs.oversamp[0]);
    if (prefs.oversamp[1] != prefs.oversamp[0])
      fprintf(file, " %g", prefs.oversamp[1]);
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
	"   <PARAM name=\"Subtract_Back\" datatype=\"boolean\""
	" ucd=\"meta.code\" value=\"%c\"/>\n",
    	prefs.subtract_flag[0]? 'T':'F');
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
VERSION	14/07/2006
 ***/
void	write_xmlerror(char *filename, char *error)
  {
   FILE			*file;

  if (!(file = fopen(filename, "w")))
    return;

  write_xml_header(file);

  fprintf(file, " </TABLE>\n");

  write_xml_meta(file, error);

  fprintf(file, "</RESOURCE>\n");
  fprintf(file, "</VOTABLE>\n");

  fclose(file);

  return;
  }


