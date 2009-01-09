/*
 				wcscelsys.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	LDACTools+
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for fitswcs.c
*
*	Last modify:	21/07/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*-------------------------------- constants --------------------------------*/

/* Equatorial coordinates of origin and pole and rotation sign of equatorial,*/
/* galactic, ecliptic and supergalactic reference frames, from Allen Astron. */
/* Quantities, 4th ed. */

char	celsysname[][2][8] = {  {"RA--", "DEC-"},
				{"GLON", "GLAT"},
				{"ELON", "ELAT"},
				{"SLON", "SLAT"},
				{""}};
double	celsysorig[][2] = {	{0.0, 0.0},
				{266.40499625, -28.93617242},
				{0.0, 0.0},
				{42.29235, 59.52315}},
	celsyspole[][2] = {	{0.0, 90.0},
				{192.85948123, 27.12825120},
				{273.85261111, 66.99111111},
				{283.7514, 15.70480}},
/* Note: the code to handle the rotation sign is not yet implemented!!! */
	celsyssign[]	= {	 1.0,
				 1.0,
				 1.0,
				 1.0};

