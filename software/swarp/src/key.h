 /*
 				key.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SWarp
*
*	Author:		E.BERTIN (IAP, Leiden & ESO)
*
*	Contents:	Keyword structure.
*
*	Last modify:	29/02/2000
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*--------------------------------- constants -------------------------------*/

#define         FIND_STRICT     0
#define         FIND_NOSTRICT   1

/*--------------------------- structure definitions -------------------------*/
/* Preference keyword */
typedef struct
  {
  char		name[16];
  enum  {P_FLOAT, P_INT, P_STRING, P_BOOL, P_KEY, P_INTLIST, P_FLOATLIST,
	P_BOOLLIST, P_KEYLIST, P_STRINGLIST} type;
  void		*ptr;			/* Pointer to the keyword value */
  int		imin, imax;		/* Range for int's */
  float		dmin, dmax;		/* Range for floats */
  char		keylist[16][16];	/* List of keywords */
  int           nlistmin;		/* Minimum number of list members */
  int           nlistmax; 		/* Maximum number of list members */
  int		*nlistptr;		/* Ptr to store the nb of read params*/
  int		flag;
  }	pkeystruct;

/*---------------------------------- protos --------------------------------*/

int	findkeys(char *str, char key[][16], int mode);

