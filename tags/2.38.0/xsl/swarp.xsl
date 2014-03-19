<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE xsl:stylesheet [
	<!ENTITY nbsp "&#160;">
	<!ENTITY deg "&#176;">
	<!ENTITY amin "&#180;">
	<!ENTITY asec "&#168;">
        <!ENTITY darr "&#8595;">
	]>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<!-- 
#				swarp.xsl
#
# Global XSL template
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	SWarp
#
#	Copyright:		(C) 2005-2011 Emmanuel Bertin - IAP/CNRS/UPMC
#
#	License:		GNU General Public License
#
#	SWarp is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
# 	(at your option) any later version.
#	SWarp is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#	You should have received a copy of the GNU General Public License
#	along with SWarp. If not, see <http://www.gnu.org/licenses/>.
#
#	Last modified:		07/06/2011
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -->

 <xsl:template match="/">
  <xsl:variable name="date" select="/VOTABLE/RESOURCE/RESOURCE[@name='MetaData']/PARAM[@name='Date']/@value"/>
  <xsl:variable name="time" select="/VOTABLE/RESOURCE/RESOURCE[@name='MetaData']/PARAM[@name='Time']/@value"/>
  <HTML>
   <HEAD>
    <link rel="shortcut icon" type="image/x-icon" href="http://astromatic.net/xsl/favicon.ico" />
    <script type="text/javascript" src="http://astromatic.net/xsl/sorttable.js"/>

    <style type="text/css">
     p {
      font-family: sans-serif;
      }
     p.italic {font-style: italic}
     body {
      margin: 10px;
      background-color: #e0e0e0;
      background-image: url("http://astromatic.net/xsl/body_bg.jpg");
      background-repeat: repeat-x;
      background-position: top;
      min-width:662px;
      }
     mono {font-family: monospace}
     elen {
      font-family: monospace;
      font-weight: bold;
      color: green
      }
     elep {
      font-family: monospace;
      font-weight: bold;
      color: red
      }
     el {
      font-family: monospace;
      font-size: 100%;
      color: black;
      }
     elm {
      font-family: monospace;
      font-size: 67%;
      white-space: nowrap;
      }
     a {text-decoration: none; font-style: bold; color: #476674}
     a:hover {text-decoration: underline;}
     #header {
      padding: 5px;
      min-width: 662px;
      background-image: url("http://astromatic.net/xsl/astromaticleft.png");
      background-repeat: repeat-x;
      background-position: left top;
      text-align: left;
      font-size: 1.2em;
      margin: 0 0 30px 0;
      color:#d3e7f0;
      font-weight: bold;
      }
     th {
      background-color:#d3e7f0;
      border-top: 1px solid white;
      border-left: 1px solid white;
      border-right: 1px solid #476674;
      border-bottom: 1px solid #476674;
      -moz-border-radius: 3px;
      -khtml-border-radius: 3px;
      -webkit-border-radius: 3px;
      border-radius: 3px;
      padding: 2px;
      line-height: 12px;
      }
     td {
      background-color:#f2f4f4;
      padding-left: 2px;
      padding-right: 2px;
      }
     table.sortable {
      border-top: 1px solid #476674;
      border-left: 1px solid #476674;
      border-right: 1px solid white;
      border-bottom: 1px solid white;
      -moz-border-radius: 3px;
      -khtml-border-radius: 3px;
      -webkit-border-radius: 3px;
      border-radius: 3px;
      }
     table.sortable a.sortheader {
      background-color:#d3e7f0;
      font-weight: bold;
      font-size: 80%;
      text-decoration: none;
      display: button;
      }

     table.sortable span.sortarrow {
      color: black;
      font-weight: bold;
      text-decoration: blink;
      }
     table.sortable a.sortheader.sub {vertical-align: sub}
     </style>

     <title>
      Processing summary on <xsl:value-of select="$date"/> at <xsl:value-of select="$time"/>
     </title>
    </HEAD>
    <BODY>
     <div id="header">
      <a href="/"><img style="vertical-align: middle; border:0px" src="http://astromatic.net/xsl/astromatic.png" title="Astromatic home" alt="Astromatic.net" /></a>  Processing summary
     </div>
     <xsl:call-template name="VOTable"/>
   </BODY>
  </HTML>
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
         <xsl:value-of select='concat(format-number($duration,"##0.0"), " s")'/>
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
   <xsl:variable name="rescaleweights" select="count(FIELD[@name='Rescale_Weights']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="weightscale" select="count(FIELD[@name='Weight_Scaling']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="weighttres" select="count(FIELD[@name='Weight_Thresh']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="interp" select="count(FIELD[@name='Interpolate']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="gain" select="count(FIELD[@name='Gain']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="saturation" select="count(FIELD[@name='Saturation']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="exptime" select="count(FIELD[@name='ExpTime']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="photfscale" select="count(FIELD[@name='Photometric_Flux_Scaling']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="astrfscale" select="count(FIELD[@name='Astrometric_Flux_Scaling']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fieldcoord" select="count(FIELD[@name='Field_Coordinates']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pixscale" select="count(FIELD[@name='Pixel_Scale']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="obsdate" select="count(FIELD[@name='ObsDate']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="equinox" select="count(FIELD[@name='Equinox']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="coosys" select="count(FIELD[@name='COOSYS']/preceding-sibling::FIELD)+1"/>
   <p>
    <BUTTON type="button" title="click to expand" onclick="showhideTable('SWarp')">
     Summary Table on Input Files&nbsp;&darr;
    </BUTTON>
    <TABLE class="sortable" id="SWarp" style="display: none">
     <TR>
      <TH>Index</TH>
      <TH>Image Name</TH>
      <TH>Weight Name</TH>
      <TH>External Header</TH>
      <TH>Identifier</TH>
      <TH>Next</TH>
      <TH>Date</TH>
      <TH>Time</TH>
      <TH>Duration</TH>
      <TH>Background Mean</TH>
      <TH>Background Deviation</TH>
      <TH>Subtract Backgrond</TH>
      <TH>Background Type</TH>
      <TH>BACK_SIZE</TH>
      <TH>BACK_FILTERSIZE</TH>
      <TH>Background Default</TH>
      <TH>Weight Type</TH>
      <TH>Rescale Weights</TH>
      <TH>Weight Scaling</TH>
      <TH>Weight Threshold</TH>
      <TH>Interpolate</TH>
      <TH>Gain</TH>
      <TH>Saturation</TH>
      <TH>Exposure Time</TH>
      <TH>Photo Flux Scale</TH>
      <TH>Astro Flux Scale</TH>
      <TH>Field Coordinates</TH>
      <TH>Pixel Scale</TH>
      <TH>Observation Date</TH>
      <TH>Equinox</TH>
      <TH>Coordinate System</TH>
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
        <td align="center">
         <el><xsl:value-of select="TD[$index]"/></el>
        </td>
        <td >
         <el><xsl:value-of select="TD[$name]"/></el>
        </td>
        <td >
         <el><xsl:value-of select="TD[$weight]"/></el>
        </td>
        <td align="center">
         <xsl:choose>
          <xsl:when test="TD[$headflag] = 'T'">
           <elen>H</elen>
          </xsl:when>
         </xsl:choose>
        </td>
        <td align="center">
         <el><xsl:value-of select="TD[$ident]"/></el>
        </td>
        <td align="center">
         <el><xsl:value-of select="TD[$next]"/></el>
        </td>
        <td align="center">
         <el><xsl:value-of select="TD[$date]"/></el>
        </td>
        <td align="center">
         <el><xsl:value-of select="TD[$time]"/></el>
        </td>
        <td align="right">
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
             <el><xsl:value-of select='concat(format-number(TD[$duration],"##0.0"), " s")'/></el>
            </xsl:otherwise>
           </xsl:choose>
          </xsl:otherwise>
         </xsl:choose>
        </td>
        <td align="right">
         <el><xsl:value-of select="format-number(TD[$backmean],'##0000.00')"/></el>
        </td>
        <td align="right">
          <el><xsl:value-of select="format-number(TD[$backdev],'##00.00')"/></el>
       </td>
        <td align="center">
         <xsl:choose>
          <xsl:when test="TD[$subback] = 'T'">
           <elen>Y</elen>
          </xsl:when>
          <xsl:otherwise>
           <elen>N</elen>
          </xsl:otherwise>
         </xsl:choose>
        </td>
        <td align="center">
         <el><xsl:value-of select="TD[$backtype]"/></el>
        </td>
        <td align="center">
         <el><xsl:value-of select="TD[$backsize]"/></el>
        </td>
        <td align="center">
         <el><xsl:value-of select="TD[$backfsize]"/></el>
        </td>
        <td align="center">
         <el><xsl:value-of select="TD[$backdef]"/></el>
        </td>
        <td align="center">
         <el><xsl:value-of select="TD[$weighttype]"/></el>
        </td>
        <td align="center">
         <xsl:choose>
          <xsl:when test="TD[$rescaleweights] = 'T'">
           <elen>Y</elen>
          </xsl:when>
          <xsl:otherwise>
           <elen>N</elen>
          </xsl:otherwise>
         </xsl:choose>
        </td>
        <td align="right">
         <el><xsl:value-of select="TD[$weightscale]"/></el>
        </td>
        <td align="right">
         <el><xsl:value-of select="TD[$weighttres]"/></el>
        </td>
        <td align="center">
         <xsl:choose>
          <xsl:when test="TD[$interp] = 'T'">
           <elen>Y</elen>
          </xsl:when>
          <xsl:otherwise>
           <elen>N</elen>
          </xsl:otherwise>
         </xsl:choose>
        </td>
        <td align="right">
         <el><xsl:value-of select="format-number(TD[$gain],'##0.000')"/></el>
        </td>
        <td align="right">
         <el><xsl:value-of select="format-number(TD[$saturation],'#########')"/></el>
        </td>
        <td align="right">
         <el><xsl:value-of select='concat(TD[$exptime]," s")'/></el>
        </td>
        <td align="right">
         <el><xsl:value-of select="TD[$photfscale]"/></el>
        </td>
        <td align="right">
         <el><xsl:value-of select="TD[$astrfscale]"/></el>
        </td>
        <td align="right">
         <el><xsl:value-of select="TD[$fieldcoord]"/></el>
        </td>
        <td align="right">
         <el><xsl:value-of select="format-number(TD[$pixscale],'##0.000')"/></el>
        </td>
        <td align="center">
         <el><xsl:value-of select="TD[$obsdate]"/></el>
        </td>
        <td align="center">
         <el><xsl:value-of select="TD[$equinox]"/></el>
        </td>
        <td align="center">
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
    <BUTTON type="button" title="click to expand" onclick="showhideTable('config')">
     Configuration File:
     <B><xsl:value-of select="PARAM[@name='Prefs_Name']/@value"/></B>
     &darr;
    </BUTTON>
    <TABLE id="config" class="sortable" style="display: none">
     <TR>
      <TH>Config Parameter</TH>
      <TH>Value</TH>
     </TR>
     <xsl:for-each select="PARAM[position()>2]">
      <tr>
       <td><el><xsl:value-of select="@name"/></el></td>
       <td><el><xsl:value-of select="@value"/></el></td>
      </tr>
     </xsl:for-each>
    </TABLE>
   </p>
   <p>
    <BUTTON type="button" title="click to expand" onclick="showhideTable('commandline')">
     Command Line&nbsp;&darr;
    </BUTTON>
    <TABLE id="commandline" style="display: none">
     <TR>
      <TD style="font-size: 80%;"><el><xsl:value-of select="PARAM[@name='Command_Line']/@value"/></el></TD>
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
    <BUTTON type="button" title="click to expand" onclick="showhideTable('warnings')">
     Warnings (limited to the last 100)&nbsp;&darr;
    </BUTTON>
    <TABLE id="warnings" class="sortable" style="display: none">
     <TR style="font-size: 80%;">
      <TH>Date</TH>
      <TH>Time</TH>
      <TH>Message</TH>
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
        <td >
         <el><xsl:value-of select="TD[$date]"/></el>
        </td>
        <td>
         <el><xsl:value-of select="TD[$time]"/></el>
        </td>
        <td align="center">
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
