/*************************************************************************
 * overlapAnalysis.c -
 *
 * $Id: overlapAnalysis.c,v 1.0 2015/04/14 11:01:51 gael Exp $
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
  outputFormat format;
  analysisMethod method;
  double lt;
  double ht;
} local_par;



/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );

static int _verbose_ = 0;



static char *usage = "[image-in] [file-out]\n\
\t [-light|-complete|-matrix] [-max|-nonnull] [-inv] [-swap] [-v] [-D] [-help]";
/*
 * \t [-light|-complete|-matrix] [-max|-nonnull] [-[lt|sb] %lf] [-[ht|sh] %lf] [-inv] [-swap] [-v] [-D] [-help]";
 */

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'file-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -light|-complete|-matrix : output format (light-complete : txt file ; matrix : image file)\n\
\t -light : each line of the output file gives the number of labels of A - the number of labels of B per connected components\n\
\t -complete : each line of the output file gives the list of labels of A - the list of labels of B per connected components\n\
\t -matrix : same size as the input image, where the non-null values are identifying the connected labels\n\
\t -max|-nonnull : analysis method for the graph construction and connected components extraction (default : nonnull)\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

/*
 * \t -[lt|sb] %lf: skips overlapping values under this low threshold parameter\n\t -[ht|sh] %lf: correspondences are built only if the overlapping value is higher than this high threshold parameter\n
 */

static char program[STRINGLENGTH];



int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *imAB;
  vt_image imL;
  int sizeA, sizeB;
  unsigned short int LAB;


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

  if (VT_OverlapAnalysis(imAB, &imL, par.method, par.lt, par.ht, &LAB, &sizeA, &sizeB) != 1) {
      VT_FreeImage(imAB);
      VT_ErrorParse("unable to analyse the overlap image\n", 0);
  }



  /*--- Ecriture resultat ---*/
  if (par.format == MATRIX)
    VT_WriteInrimageWithName(&imL, par.names.out);
  else {
    char *output=NULL;
    /* output will contain the character array that has to be written */
    FILE *fichier=NULL;
    VT_WriteOverlapAnalysis(&imL, &output, LAB, sizeA, sizeB, par.format);
    if (strcmp(par.names.out, ">") == 0)
      fichier = stdout;
    else
      fichier = fopen(par.names.out, "w");
    if (fichier == NULL){
        VT_FreeImage(imAB);
        VT_FreeImage(&imL);
        if (output)
        {
          free(output);
          output = NULL;
        }
        VT_ErrorParse("unable to convert into text file the overlap analysis\n", 0);
    }
    fprintf(fichier, "%s", output);
    fclose(fichier);
    if (output)
    {
      free(output);
      output = NULL;
    }
  }


  /*--- liberations memoire ---*/

  VT_FreeImage(imAB);
  VT_FreeImage(&imL);

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

      else if ( strcmp ( argv[i], "-light" ) == 0 ) {
        par->format = LIGHT;
      }
      else if ( strcmp ( argv[i], "-complete" ) == 0 ) {
        par->format = COMPLETE;
      }
      else if ( strcmp ( argv[i], "-matrix" ) == 0 ) {
        par->format = MATRIX;
      }

      else if ( strcmp ( argv[i], "-max" ) == 0 ) {
        par->method = MAX;
      }
      else if ( strcmp ( argv[i], "-nonnull" ) == 0 ) {
        par->method = NONNULL;
      }
      else if ( strcmp ( argv[i], "-sh" ) == 0 || strcmp ( argv[i], "-ht" ) == 0) {
        /* par->method = THRESHOLD;
         */
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -ht...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->ht));
        if ( status <= 0 ) VT_ErrorParse( "parsing -ht...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-sb" ) == 0 || strcmp ( argv[i], "-lt" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -lt...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->lt));
        if ( status <= 0 ) VT_ErrorParse( "parsing -lt...\n", 0 );
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
  par->format = LIGHT;
  par->method = MAX;
  par->lt = -1.0;
  par->ht = -1.0;
}

