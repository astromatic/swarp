<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE xsl:stylesheet [
	<!ENTITY nbsp "&#160;">
	<!ENTITY deg "&#176;">
	<!ENTITY amin "&#180;">
	<!ENTITY asec "&#168;">
	]>

<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<!-- *********************** Global XSL template ************************** -->
 <xsl:template match="/">
  <xsl:variable name="date" select="/VOTABLE/RESOURCE/RESOURCE[@name='MetaData']/PARAM[@name='Date']/@value"/>
  <xsl:variable name="time" select="/VOTABLE/RESOURCE/RESOURCE[@name='MetaData']/PARAM[@name='Time']/@value"/>
  <html>

<!-- HTML head  -->

   <head>

<!-- javascript -->

<!--  <script type="text/javascript" language="javascript"> -->

    <script src="http://terapix.iap.fr/cplt/xsl/sorttable.js"/>

    <style type="text/css">
     p.sansserif {font-family: sans-serif}
     body {background-color: white}
     mono {font-family: monospace}
     elen {font-family: monospace; font-size: 100%; font-weight: bold; color: green }
     elep {font-family: monospace; font-size: 100%; font-weight: bold; color: red }
     el {font-family: monospace; font-size: 100%; color: black}
     a {text-decoration: none}
     table.sortable a.sortheader
      {
      background-color:#FFEECC;
      color: black;
      font-weight: bold;
      font-size: 80%;
      text-decoration: none;
      display: button;
      }
     table.sortable span.sortarrow
      {
      color: black;
      font-weight: bold;
      text-decoration: none;
      }
     table.sortable a.sortheader.sub
      {
      vertical-align: sub;
      }
     </style>

     <title>
      Processing summary on <xsl:value-of select="$date"/> at <xsl:value-of select="$time"/>
     </title>
    </head>

<!-- HTML body -->

    <BODY>
     <TABLE BORDER="0" CELLPADDING="0" CELLSPACING="0" WIDTH="100%">
      <TR>
       <TD ALIGN="LEFT">
        <TABLE BORDER="0">
         <TR>
          <TD ALIGN="CENTER">
           <IMG SRC="http://terapix.iap.fr/cplt/xsl/terapixLogo.png" ALT="Terapix"/>
          </TD>
          <TD ALIGN="CENTER">
           <IMG SRC="http://terapix.iap.fr/cplt/xsl/terapixTitle.png" ALT="Logo"/>
          </TD>
          <TD ALIGN="CENTER">
           <FONT color="#669933">
            <B> Processing summary</B>
           </FONT>
          </TD>
          <TD ALIGN="CENTER">
           <IMG SRC="http://terapix.iap.fr/cplt/xsl/terapixPicture.gif" ALT="Terapix banner"/>
          </TD>
         </TR>
        </TABLE>
       </TD>
      </TR>
      <TR>
       <TD>
        <TABLE BORDER="0" WIDTH="100%" BGCOLOR="#000000">
         <TR>
          <TH BGCOLOR="#000000" ALIGN="LEFT"><FONT SIZE="-1" COLOR="#FFFFFF"> Home > Tools > Data reduction</FONT></TH>
         </TR>
        </TABLE>
       </TD>
      </TR>
     </TABLE>
    <xsl:call-template name="VOTable"/>
   </BODY>
  </html>
 </xsl:template>
<!-- **************** Generic XSL template for VOTables ****************** -->
 <xsl:template name="VOTable">
  <xsl:for-each select="/VOTABLE">
   <xsl:call-template name="Resource"/>
  </xsl:for-each>
 </xsl:template>
<!-- *************** Generic XSL template for Resources ****************** -->
 <xsl:template name="Resource">
  <xsl:for-each select="RESOURCE">
   <xsl:choose>
    <xsl:when test="@ID='SWarp'">
     <xsl:call-template name="swarp"/>
    </xsl:when>
   </xsl:choose>
  </xsl:for-each>
 </xsl:template>
<!-- ********************** XSL template for SWarp *********************** -->
 <xsl:template name="swarp">
  <xsl:for-each select="RESOURCE[@ID='MetaData']">
   <xsl:call-template name="RunInfo"/>
   <xsl:for-each select="TABLE[@ID='Input_Image_Data']">
    <xsl:call-template name="Input_Image_Data"/>
   </xsl:for-each>
   <xsl:for-each select="RESOURCE[@ID='Config']">
    <xsl:call-template name="Config"/>
   </xsl:for-each>
   <xsl:for-each select="TABLE[@ID='Warnings']">
    <xsl:call-template name="Warnings"/>
   </xsl:for-each>
  </xsl:for-each>
 </xsl:template>
<!-- ************* Generic XSL RunInfo template for MetaData ************* -->
 <xsl:template name="RunInfo">
  <p>
<!-- Software name, version, date, time and number of threads -->
   <a>
    <xsl:attribute name="href">
     <xsl:value-of select="PARAM[@name='Soft_URL']/@value"/>
    </xsl:attribute>
    <b>
     <xsl:value-of select="PARAM[@name='Software']/@value"/>&nbsp;<xsl:value-of select="PARAM[@name='Version']/@value"/>
    </b>
   </a>
   started on
   <b><xsl:value-of select="PARAM[@name='Date']/@value"/></b>
   at
   <b><xsl:value-of select="PARAM[@name='Time']/@value"/></b>
   with
   <b><xsl:value-of select="PARAM[@name='NThreads']/@value"/></b>
   thread<xsl:if test="PARAM[@name='NThreads']/@value &gt; 1">s</xsl:if>

<!-- Run time -->
   <xsl:variable name="duration" select="PARAM[@name='Duration']/@value"/>
   (run time:
    <b>
     <xsl:choose> 
      <xsl:when test="$duration &gt; 3600.0">
       <xsl:value-of
	select='concat(string(floor($duration div 3600)),
	" h ", format-number(floor(($duration div 60) mod 60.0), "00"),
	" min")'/>
      </xsl:when>
      <xsl:otherwise>
       <xsl:choose>
        <xsl:when test="$duration &gt; 60.0">
         <xsl:value-of
	  select='concat(format-number(floor($duration div 60),"##"),
	  " min ", format-number(floor($duration mod 60.0), "00")," s")'/>
        </xsl:when>
        <xsl:otherwise>
         <xsl:value-of select='concat(string($duration), " s")'/>
        </xsl:otherwise>
       </xsl:choose>
      </xsl:otherwise>
     </xsl:choose>
    </b>)
    <br />
   by user <b><xsl:value-of select="PARAM[@name='User']/@value"/></b>
   from <b><xsl:value-of select="PARAM[@name='Host']/@value"/></b>
   in <b><mono><xsl:value-of select="PARAM[@name='Path']/@value"/></mono></b>
  </p>
  <p>
   <b style="color: red"><xsl:if test="PARAM[@name='Error_Msg']/@value &gt; 0">
    An Error occured!!! </xsl:if>
   <xsl:value-of select="PARAM[@name='Error_Msg']/@value"/></b>
  </p>
  <p>
  <sans-serif><i>click to expand or hide tables</i></sans-serif>
  </p>
 </xsl:template>
<!-- ********************** XSL template for Input_Image_Data ************** -->
  <xsl:template name="Input_Image_Data">
   <xsl:variable name="index" select="count(FIELD[@name='Frame_Index']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="name" select="count(FIELD[@name='Image_Name']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="weight" select="count(FIELD[@name='Weight_Name']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="headflag" select="count(FIELD[@name='External_Header']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ident" select="count(FIELD[@name='Image_Ident']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="next" select="count(FIELD[@name='Extension']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="date" select="count(FIELD[@name='Date']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="time" select="count(FIELD[@name='Time']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="duration" select="count(FIELD[@name='Duration']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="backmean" select="count(FIELD[@name='Background_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="backdev" select="count(FIELD[@name='Background_StDev']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="subback" select="count(FIELD[@name='Subtract_Back']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="backtype" select="count(FIELD[@name='Back_Type']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="backsize" select="count(FIELD[@name='Back_Size']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="backfsize" select="count(FIELD[@name='Back_FilterSize']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="backdef" select="count(FIELD[@name='Back_Default']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="weighttype" select="count(FIELD[@name='Weight_Type']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="weighttres" select="count(FIELD[@name='Weight_Thresh']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="weightscale" select="count(FIELD[@name='Weight_Scaling']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="interp" select="count(FIELD[@name='Interpolate']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="gain" select="count(FIELD[@name='Gain']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="saturation" select="count(FIELD[@name='Saturation']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="exptime" select="count(FIELD[@name='ExpTime']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="photfscale" select="count(FIELD[@name='Photometric_Flux_Scaling']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astrfscale" select="count(FIELD[@name='Astrometric_Flux_Scaling']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fieldcoord" select="count(FIELD[@name='Field_Coordinates']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pixscale" select="count(FIELD[@name='Pixel_Scale']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="equinox" select="count(FIELD[@name='Equinox']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="epoch" select="count(FIELD[@name='Epoch']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="coosys" select="count(FIELD[@name='COOSYS']/preceding-sibling::FIELD)+1"/>
   <p>
    <BUTTON type="button" style="background:#CCEECC; font-family: sans-serif; font-weight: bold;" onclick="showhideTable('SWarp')">
     Summary Table on Input Files
    </BUTTON>
    <TABLE class="sortable" id="SWarp" BORDER="2" style="display: none">
     <TR>
      <TH BGCOLOR="#FFEECC">Index</TH>
      <TH BGCOLOR="#FFEECC">Image Name</TH>
      <TH BGCOLOR="#FFEECC">Weight Name</TH>
      <TH BGCOLOR="#FFEECC">External Header</TH>
      <TH BGCOLOR="#FFEECC">Identifier</TH>
      <TH BGCOLOR="#FFEECC">Next</TH>
      <TH BGCOLOR="#FFEECC">Date</TH>
      <TH BGCOLOR="#FFEECC">Time</TH>
      <TH BGCOLOR="#FFEECC">Duration</TH>
      <TH BGCOLOR="#FFEECC">Background Mean</TH>
      <TH BGCOLOR="#FFEECC">Background Deviation</TH>
      <TH BGCOLOR="#FFEECC">Subtract Backgrond</TH>
      <TH BGCOLOR="#FFEECC">Background Type</TH>
      <TH BGCOLOR="#FFEECC">BACK_SIZE</TH>
      <TH BGCOLOR="#FFEECC">BACK_FILTERSIZE</TH>
      <TH BGCOLOR="#FFEECC">Background Default</TH>
      <TH BGCOLOR="#FFEECC">Weight Type</TH>
      <TH BGCOLOR="#FFEECC">Weight Threshold</TH>
      <TH BGCOLOR="#FFEECC">Weight Scaling</TH>
      <TH BGCOLOR="#FFEECC">Interpolate</TH>
      <TH BGCOLOR="#FFEECC">Gain</TH>
      <TH BGCOLOR="#FFEECC">Saturation</TH>
      <TH BGCOLOR="#FFEECC">Exposure Time</TH>
      <TH BGCOLOR="#FFEECC">Photo Flux Scale</TH>
      <TH BGCOLOR="#FFEECC">Astro Flux Scale</TH>
      <TH BGCOLOR="#FFEECC">Field Coordinates</TH>
      <TH BGCOLOR="#FFEECC">Pixel Scale</TH>
      <TH BGCOLOR="#FFEECC">Equinox</TH>
      <TH BGCOLOR="#FFEECC">Epoch</TH>
      <TH BGCOLOR="#FFEECC">Coordinate System</TH>
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$index]"/></el>
        </td>
        <td  BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$name]"/></el>
        </td>
        <td  BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$weight]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <xsl:choose>
          <xsl:when test="TD[$headflag] = 'T'">
           <elen>H</elen>
          </xsl:when>
         </xsl:choose>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$ident]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$next]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$date]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$time]"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <xsl:choose> 
          <xsl:when test="$duration &gt; 3600.0">
           <el><xsl:value-of
	    select='concat(format-number(floor(TD[$duration] div 3600),"##00"),
	    " h ", format-number(floor((TD[$duration] div 60) mod 60.0),"##00"),
	    " min")'/></el>
          </xsl:when>
          <xsl:otherwise>
           <xsl:choose>
            <xsl:when test="$duration &gt; 60.0">
             <el><xsl:value-of
	      select='concat(format-number(floor(TD[$duration] div 60),"##00"),
	      " min ", format-number(floor(TD[$duration] mod 60.0),"##00")," s")'/></el>
            </xsl:when>
            <xsl:otherwise>
             <el><xsl:value-of select='concat(format-number(TD[$duration],"##00"), " s")'/></el>
            </xsl:otherwise>
           </xsl:choose>
          </xsl:otherwise>
         </xsl:choose>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$backmean],'##0000.00')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
          <el><xsl:value-of select="format-number(TD[$backdev],'##00.00')"/></el>
       </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <xsl:choose>
          <xsl:when test="TD[$subback] = 'T'">
           <elen>Y</elen>
          </xsl:when>
          <xsl:otherwise>
           <elen>N</elen>
          </xsl:otherwise>
         </xsl:choose>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$backtype]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$backsize]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$backfsize]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$backdef]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$weighttype]"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$weighttres]"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$weightscale]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <xsl:choose>
          <xsl:when test="TD[$interp] = 'T'">
           <elen>Y</elen>
          </xsl:when>
          <xsl:otherwise>
           <elen>N</elen>
          </xsl:otherwise>
         </xsl:choose>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$gain],'##0.000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$saturation],'#########')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select='concat(TD[$exptime]," s")'/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$photfscale]"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$astrfscale]"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$fieldcoord]"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$pixscale],'##0.000')"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$equinox]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$epoch]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$coosys]"/></el>
        </td>
       </tr>
      </xsl:for-each>
     </xsl:for-each>
    </TABLE>
   </p>
 </xsl:template>

<!-- ********************** XSL template for Config File ********************** -->
  <xsl:template name="Config">
   <p>
    <BUTTON type="button" style="background:#CCEECC; font-family: sans-serif; font-weight: bold;" onclick="showhideTable('config')">
     Configuration File: <xsl:value-of select="PARAM[@name='Prefs_Name']/@value"/>
    </BUTTON>
    <TABLE id="config" class="sortable" style="display: none">
     <TR>
      <TH BGCOLOR="#FFEECC">Config Parameter</TH>
      <TH BGCOLOR="#FFEECC">Value</TH>
     </TR>
     <xsl:for-each select="PARAM[position()>2]">
      <tr BGCOLOR="#EEEEEE">
       <td><el><xsl:value-of select="@name"/></el></td>
       <td><el><xsl:value-of select="@value"/></el></td>
      </tr>
     </xsl:for-each>
    </TABLE>
   </p>
   <p>
    <BUTTON type="button" style="background:#CCEECC; font-family: monospace; font-weight: bold: font-size: 80%;" onclick="showhideTable('commandline')">
     Command Line
    </BUTTON>
    <TABLE id="commandline" style="display: none">
     <TR>
      <TD BGCOLOR="#FFEECC" style="font-size: 80%;"><el>Command Line: <xsl:value-of select="PARAM[@name='Command_Line']/@value"/></el></TD>
     </TR>
    </TABLE>
   </p>
  </xsl:template>

<!-- ********************** XSL template for Warnings ********************** -->
  <xsl:template name="Warnings">
   <xsl:variable name="date" select="count(FIELD[@name='Date']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="time" select="count(FIELD[@name='Time']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="msg" select="count(FIELD[@name='Msg']/preceding-sibling::FIELD)+1"/>
   <p>
    <BUTTON type="button" style="background:#CCEECC; font-family: monospace; font-weight: bold: font-size: 80%;" onclick="showhideTable('warnings')">
     Warnings (limited to the last 100)
    </BUTTON>
    <TABLE id="warnings" style="display: none">
     <TR style="font-size: 80%;">
      <TH BGCOLOR="#FFEECC">Date</TH>
      <TH BGCOLOR="#FFEECC">Time</TH>
      <TH BGCOLOR="#FFEECC">Message</TH>
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
        <td  BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$date]"/></el>
        </td>
        <td BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$time]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$msg]"/></el>
        </td>
       </tr>
      </xsl:for-each>
     </xsl:for-each>
    </TABLE>
   </p>
 </xsl:template>

 <xsl:template name="Rest">
</xsl:template>

</xsl:stylesheet>
