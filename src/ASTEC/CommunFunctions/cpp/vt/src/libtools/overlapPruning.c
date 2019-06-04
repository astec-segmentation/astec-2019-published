/*************************************************************************
 * overlapPruning.c -
 *
 * $Id: overlapPruning.c,v 1.0 2015/04/14 14:10:51 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2015/04/14
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>
#include <vt_overlapTools.h>
#define TABLENGTH 2000


typedef struct local_par {
  vt_names names;
  double epsilon;
} local_par;



/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );

static int _verbose_ = 0;



static char *usage = "[image-in] [image-out] [-e %lf]\n\
\t [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'file-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -e %lf : threshold value for the pruning (default = 0)\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];



int main( int argc, char *argv[] )
{

  local_par par;
  vt_image *imAB;



  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );

  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );

  /*--- lecture des images d'entree ---*/
  imAB = _VT_Inrimage( par.names.in );
  if ( imAB == (vt_image*)NULL )
    VT_ErrorParse("unable to read input image\n", 0);

  /*--- operations eventuelles sur les images d'entree ---*/
  if ( par.names.inv == 1 )
      VT_InverseImage( imAB );

  if ( par.names.swap == 1 )
      VT_SwapImage( imAB );

  /* analyse */
  if(VT_OverlapPruning(imAB, par.epsilon) != 1) {
      VT_FreeImage(imAB);
      VT_ErrorParse("unexpected error during the overlap pruning operation\n", 0);
  }

  /* imwrite */
  VT_WriteInrimageWithName(imAB, par.names.out);

  /*--- liberations memoire ---*/

  VT_FreeImage(imAB);

  return( 0 );
}






static void VT_Parse( int argc,
                      char *argv[],
                      local_par *par )
{
  int i, nb;
  char text[STRINGLENGTH];
  int status;

  if ( VT_CopyName( program, argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);
  if ( argc == 1 ) VT_ErrorParse("\n", 0 );

  /*--- lecture des parametres ---*/
  i = 1; nb = 0;
  while ( i < argc ) {
    if ( argv[i][0] == '-' ) {
      if ( argv[i][1] == '\0' ) {
        if ( nb == 0 ) {
          /*--- standart input ---*/
          strcpy( par->names.in, "<" );
          nb += 1;
        }
      }
      /*--- arguments generaux ---*/
      else if ( strcmp ( argv[i], "-help" ) == 0 ) {
        VT_ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 ) {
        _VT_VERBOSE_ = 1;
        _verbose_++;
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 ) {
        _VT_DEBUG_ = 1;
      }

      /*--- traitement eventuel de l'image d'entree ---*/
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
        par->names.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
        par->names.swap = 1;
      }


      /* Parametres de calcul */

      /*else if ( strcmp ( argv[i], "-probability" ) == 0 ) {
        par->measure = PROBABILITY;
      }
      else if ( strcmp ( argv[i], "-dice" ) == 0 ) {
        par->measure = DICE;
      }
      else if ( strcmp ( argv[i], "-jaccard" ) == 0 ) {
        par->measure = JACCARD;
      }*/

      else if ( strcmp ( argv[i], "-e") == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -e...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->epsilon));
        if ( status <= 0 ) VT_ErrorParse( "parsing -e...\n", 0 );
      }

      /*--- option inconnue ---*/
      else {
        sprintf(text,"unknown option %s\n",argv[i]);
        VT_ErrorParse(text, 0);
      }
    }
    /*--- saisie des noms d'images ---*/
    else if ( argv[i][0] != 0 ) {
      if ( nb == 0 ) {
        strncpy( par->names.in, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else if ( nb == 1 ) {
        strncpy( par->names.out, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else
        VT_ErrorParse("too much file names when parsing\n", 0 );
    }
    i += 1;
  }

  /*--- s'il n'y a pas assez de noms ... ---*/
  if (nb == 0) {
      strcpy( par->names.in,  "<" );  /* standart input */
      strcpy( par->names.out, ">" );  /* standart output */
  }
  if (nb == 1) {
    strcpy( par->names.out, ">" );  /* standart output */
  }

}




static void VT_ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"Erreur : %s",str);
  exit( 1 );
}




static void VT_InitParam( local_par *par )
{
  VT_Names( &(par->names) );
  par->epsilon = 0.0;
}

