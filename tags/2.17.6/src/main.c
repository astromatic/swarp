 /*
 				main.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SWarp
*
*	Author:		E.Bertin & C.Marmo (IAP)
*
*	Contents:	Parsing of the command line.
*
*	Last modify:	08/01/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include	<ctype.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#ifdef HAVE_MPI
#include	<mpi.h>
#endif

#include	"define.h"
#include	"globals.h"
#include	"fits/fitscat.h"
#include	"prefs.h"

#define		SYNTAX \
BANNER " image1 [image2 ...][@image_list1 [@image_list2 ...]]\n" \
"\t\t[-c <configuration_file>][-<keyword> <value>]\n" \
"> to dump a default configuration file: " BANNER " -d \n" \
"> to dump a default extended configuration file: " BANNER " -dd \n"

extern const char	notokstr[];

/********************************** main ************************************/
int	main(int argc, char *argv[])

  {
   char		**argkey, **argval,
		*str, *listbuf;
   int		a, narg, nim, ntok, opt, opt2, bufpos,bufsize;

#ifdef HAVE_MPI
  MPI_Init (&argc,&argv);
#endif

  if (argc<2)
    {
    fprintf(OUTPUT, "\n		%s  version %s (%s)\n", BANNER,MYVERSION,DATE);
    fprintf(OUTPUT, "\nby %s\n", COPYRIGHT);
    fprintf(OUTPUT, "visit %s\n", WEBSITE);
    error(EXIT_SUCCESS, "SYNTAX: ", SYNTAX);
    }
  QMALLOC(argkey, char *, argc);
  QMALLOC(argval, char *, argc);

/*default parameters */
  prefs.command_line = argv;
  prefs.ncommand_line = argc;
  prefs.ninfield = 1;
  prefs.infield_name[0] = "image";
  strcpy(prefs.prefs_name, "default.swarp");
  narg = nim = 0;
  listbuf = (char *)NULL;
  bufpos = 0;
  bufsize = MAXCHAR*1000;

  for (a=1; a<argc; a++)
    {
    if (*(argv[a]) == '-')
      {
      opt = (int)argv[a][1];
      if (strlen(argv[a])<4 || opt == '-')
        {
        opt2 = (int)tolower((int)argv[a][2]);
        if (opt == '-')
          {
          opt = opt2;
          opt2 = (int)tolower((int)argv[a][3]);
          }
        switch(opt)
          {
          case 'c':
            if (a<(argc-1))
              strcpy(prefs.prefs_name, argv[++a]);
            break;
          case 'd':
            dumpprefs(opt2=='d' ? 1 : 0);
            exit(EXIT_SUCCESS);
            break;
          case 'v':
            printf("%s version %s (%s)\n", BANNER,MYVERSION,DATE);
            exit(EXIT_SUCCESS);
            break;
          case 'h':
          default:
            error(EXIT_SUCCESS,"SYNTAX: ", SYNTAX);
          }
        }
      else
        {
        argkey[narg] = &argv[a][1];
        argval[narg++] = argv[++a];
        }       
      }
    else
      {
/*---- The input image filename(s) */
      for(; (a<argc) && (*argv[a]!='-'); a++)
        {
        str = (*argv[a] == '@'? listbuf=list_to_str(argv[a]+1) : argv[a]);
        for (ntok=0; (str=strtok(ntok?NULL:str, notokstr)); nim++,ntok++)
          if (nim<MAXINFIELD)
            prefs.infield_name[nim] = str;
          else
            error(EXIT_FAILURE, "*Error*: Too many input images: ", str);
        }
      prefs.ninfield = nim;
      a--;
      }
    }
  prefs.ninfield = nim;
  readprefs(prefs.prefs_name, argkey, argval, narg);
  useprefs();

  free(argkey);
  free(argval);

#ifdef HAVE_MPI
  MPI_Comm_size(MPI_COMM_WORLD,&prefs.nnodes);
  MPI_Comm_rank(MPI_COMM_WORLD,&prefs.node_index);
#endif

  makeit();

  free(listbuf);

  NFPRINTF(OUTPUT, "");
  NPRINTF(OUTPUT, "> All done (in %.0f s)\n", prefs.time_diff);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  exit(EXIT_SUCCESS);
  }

