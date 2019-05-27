/*************************************************************************
 * segmentationOverlapping.c -
 *
 * $Id: segmentationOverlapping.c,v 1.0 2015/04/08 18:22:51 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2015/04/08
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>
#include <vt_overlapTools.h>

#define TABLENGTH 2000



typedef struct local_par {
  vt_names names;
  double r;
  int flag_real;
  int bckgrdA;
  int bckgrdB;
  measureType measure;
} local_par;



/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );

static int _verbose_ = 0;



static char *usage = "[image-in] [image-ext] [image-out]\n\
\t [-r | -rr %lf] [-rv %d] [-bckgrdA|-bA %d] [-bckgrdB|-bB %d] [-probability|-dice|-jaccard] [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-out' est absent, on prendra stdout\n\
\t -r|-rr %lf : erosion parameter in real coordinates (default is 0)\n\
\t -rv %d : erosion parameter in voxel coordinates (default is 0)\n\
\t -bckgrdA|-bA : background parameter for image in (default is 0)\n\
\t -bckgrdB|-bB : background parameter for image ext (default is 0)\n\
\t -probability|-dice|-jaccard : method used for the connexity function\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];



int main( int argc, char *argv[] )
{

  local_par par;
  vt_image *imA;
  vt_image *imB;
  vt_image imAB;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );

  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );

  /*--- lecture des images d'entree ---*/
  imA = _VT_Inrimage( par.names.in );
  if ( imA == (vt_image*)NULL )
    VT_ErrorParse("unable to read input image A\n", 0);
  imB = _VT_Inrimage( par.names.ext );
  if ( imB == (vt_image*)NULL ) {
    VT_FreeImage(imA);
    VT_ErrorParse("unable to read input image B\n", 0);
  }

  /*--- operations eventuelles sur les images d'entree ---*/
  if ( par.names.inv == 1 )  {
      VT_InverseImage( imA );
      VT_InverseImage( imB );
  }
  if ( par.names.swap == 1 ) {
      VT_SwapImage( imA );
      VT_SwapImage( imB );
  }


  /* analysis */
  if (VT_SegmentationOverlapping(imA, imB, &imAB, par.measure, par.r, par.flag_real, par.bckgrdA, par.bckgrdB) != 1)
  {
      VT_FreeImage(imA);
      VT_FreeImage(imB);
      VT_ErrorParse("unable to do the comparison between image-in and image-ext\n", 0);
  }

  /* writing */

  VT_WriteInrimageWithName(&imAB, par.names.out);


  /*--- liberations memoire ---*/

  VT_FreeImage(imA);
  VT_FreeImage(imB);
  VT_FreeImage(&imAB);

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

      else if ( strcmp ( argv[i], "-bckgrdA" ) == 0 || strcmp ( argv[i], "-bA" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -bckgrdA...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->bckgrdA ));
        if ( status <= 0 ) VT_ErrorParse( "parsing -bckgrdA...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-bckgrdB" ) == 0 || strcmp ( argv[i], "-bB" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -bckgrdB...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->bckgrdB ));
        if ( status <= 0 ) VT_ErrorParse( "parsing -bckgrdB...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-r" ) == 0 || strcmp ( argv[i], "-rr" ) == 0 ) {
        i += 1;
        par->flag_real = 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -r...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->r));
        if ( status <= 0 ) VT_ErrorParse( "parsing -r...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-rv" ) == 0 ){
        i += 1;
        par->flag_real = 0;
        if ( i >= argc)    VT_ErrorParse( "parsing -rv...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->r));
        if ( status <= 0 ) VT_ErrorParse( "parsing -rv...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-probability" ) == 0 ) {
        par->measure = PROBABILITY;
      }
      else if ( strcmp ( argv[i], "-dice" ) == 0 ) {
        par->measure = DICE;
      }
      else if ( strcmp ( argv[i], "-jaccard" ) == 0 ) {
        par->measure = JACCARD;
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
        strncpy( par->names.ext, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else if ( nb == 2 ) {
        strncpy( par->names.out, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else
        VT_ErrorParse("too much file names when parsing\n", 0 );
    }
    i += 1;
  }

  /*--- s'il n'y a pas assez de noms ... ---*/
  if (nb == 0 || nb == 1) VT_ErrorParse("not enought file names when parsing\n", 0 );
  if (nb == 2) {
      strcpy( par->names.out,  ">" );  /* standart output */
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
  par->r=0;
  par->bckgrdA=0;
  par->bckgrdB=0;
  par->flag_real=1;
  par->measure=PROBABILITY;
}

