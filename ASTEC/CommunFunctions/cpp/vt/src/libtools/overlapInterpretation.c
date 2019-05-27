/*************************************************************************
 * overlapInterpretation.c -
 *
 * $Id: overlapPruning.c,v 1.0 2015/05/06 17:20:51 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2015/05/06
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_overlapTools.h>
/*
 * #include <vt_common.h>
 */

#define TABLENGTH 2000

typedef enum {
    VOXEL,
    CELL
} quantificationMode;


typedef struct local_par {
  vt_names fileNames;
  vt_names imageNames;
  vt_names imageOutNames;
  quantificationMode mode;
  int flagIn;
  int flagExt;
  int flagInOut;
  int flagExtOut;
} local_par;



/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );

static int _verbose_ = 0;



static char *usage = "[file-in] [file-out] [-voxel|-cell]\n\
\t [-image-in %s] [-image-ext %s] [-image-in-out %s] [-image-ext-out %s] [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t si 'file-in' est '-', on prendra stdin\n\
\t si 'file-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -pixel : computes the statistics with respect to pixel (need for -image-in and -image-ext options)\n\
\t -cell : computes the statistics with respect to cells\n\
\t -image-in %s : segmentation image in\n\
\t -image-ext %s : segmentation image ext\n\
\t -image-in-out %s : segmentation image in output name (re-labelled cells with respect to their correspondance)\n\
\t -image-ext-out %s : segmentation image ext output name (re-labelled cells with respect to their correspondance)\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];



int main( int argc, char *argv[] )
{

  local_par par;
  char *output=NULL;
  FILE *fp;
  FILE *fileout = NULL;

  vt_image *imA=NULL, *imB=NULL;
  vt_image imAout, imBout;


  int nLabelsA, nLabelsB;
  int bijectionA=0;
  int bijectionB=0;
  int divergentA=0;
  int divergentB=0;
  int Afond=0;
  int Bfond=0;
  int sous_segA[ORDERMAX];
  int sous_segB[ORDERMAX];
  int sur_segA[ORDERMAX];
  int sur_segB[ORDERMAX];


  long long vox_bijectionA=0;
  long long vox_bijectionB=0;
  long long vox_divergentA=0;
  long long vox_divergentB=0;
  long long vox_Afond=0;
  long long vox_Bfond=0;
  long long vox_sous_segA[ORDERMAX];
  long long vox_sous_segB[ORDERMAX];
  long long vox_sur_segA[ORDERMAX];
  long long vox_sur_segB[ORDERMAX];

  long long *volume_bijA=NULL;
  long long *volume_bijB=NULL;

  double *proportion_recouvrementA=NULL;
  double *proportion_recouvrementB=NULL;

  long long volume_totalA=0;
  long long volume_totalB=0;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );

  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );

  /*--- lecture des images d'entree ---*/
  if (par.flagInOut && !par.flagIn) VT_ErrorParse("option with output image needs the option with input image\n", 0);
  if (par.flagExtOut && !par.flagExt) VT_ErrorParse("option with output image needs the option with input image\n", 0);

  if (par.flagIn)
  {
      imA = _VT_Inrimage( par.imageNames.in );
      if ( imA == (vt_image*)NULL )
        VT_ErrorParse("unable to read input image in\n", 0);
  }
  if (par.flagExt)
  {
      imB = _VT_Inrimage( par.imageNames.ext );
      if ( imB == (vt_image*)NULL ) {
        if (par.flagIn) {
            VT_FreeImage(imA);
        }
        VT_ErrorParse("unable to read input image ext\n", 0);
      }
  }



  /*--- operations eventuelles sur les images d'entree ---*/

  /* tests de compatibilite d'images */
  /*if (imAB->type != FLOAT )
  {
      VT_FreeImage(imAB);
      VT_ErrorParse("input image expected type is FLOAT\n", 0);
  }
  array = (float ***) imAB->array;

  sizeA=imAB->dim.x;
  sizeB=imAB->dim.y;
  nSlices=imAB->dim.z;
  */


  /* analyse */

  if (strcmp(par.fileNames.in, "<") != 0 )
      fp = fopen( par.fileNames.in, "r" );
  else
      fp = stdin;
  if ( fp == NULL ) {
    if (par.flagIn) VT_FreeImage(imA);
    if (par.flagExt) VT_FreeImage(imB);
    VT_ErrorParse("unable to read input file\n", 0);
  }


  if (VT_OverlapInterpretation(fp, imA, imB, par.flagIn, par.flagExt, &imAout, &imBout, par.flagInOut, par.flagExtOut, &nLabelsA, &nLabelsB, &Afond, &Bfond,
                               &volume_totalA, &volume_totalB, &bijectionA, &bijectionB, sous_segA, sous_segB, sur_segA, sur_segB,
                               &divergentA, &divergentB,
                               &vox_Afond, &vox_Bfond, &vox_bijectionA, &vox_bijectionB, vox_sous_segA, vox_sous_segB, vox_sur_segA, vox_sur_segB,
                               &vox_divergentA, &vox_divergentB,
                               &volume_bijA, &volume_bijB, &proportion_recouvrementA, &proportion_recouvrementB) != 1) {
      if (par.flagIn) VT_FreeImage(imA);
      if (par.flagExt) VT_FreeImage(imB);
      VT_ErrorParse("error during the overlap interpretation", 0);
  }


  /* imwrite */

  if(par.flagInOut) {
    VT_WriteInrimageWithName(&imAout, par.imageOutNames.in);
    VT_FreeImage(&imAout);
  }

  if(par.flagExtOut) {
    VT_WriteInrimageWithName(&imBout, par.imageOutNames.ext);
    VT_FreeImage(&imBout);
  }

  /* file out */

  if (strcmp(par.fileNames.out, ">") == 0)
    fileout = stdout;
  else
    fileout = fopen( par.fileNames.out, "w" );
  if (fileout == NULL) {
      if (par.flagIn) VT_FreeImage(imA);
      if (par.flagExt) VT_FreeImage(imB);

      if(volume_bijA)
        free(volume_bijA);
      if(volume_bijB)
        free(volume_bijB);
      if(proportion_recouvrementA)
        free(proportion_recouvrementA);
      if(proportion_recouvrementB)
        free(proportion_recouvrementB);
      VT_ErrorParse("Unable to write the output file\n", 0);
  }


  VT_WriteOverlapInterpretationCell(&output, nLabelsA, nLabelsB, bijectionA, bijectionB, sous_segA, sous_segB, sur_segA, sur_segB, Afond, Bfond, divergentA, divergentB);
  if(output) {
      fprintf(fileout, "%s", output);
      free(output);
      output=NULL;
  }

  if (par.flagIn && par.flagExt) {
    VT_WriteOverlapInterpretationVoxel(&output, vox_bijectionA, vox_bijectionB, vox_sous_segA, vox_sous_segB, vox_sur_segA, vox_sur_segB, vox_Afond, vox_Bfond, vox_divergentA, vox_divergentB);
    if(output) {
      fprintf(fileout, "%s", output);
      free(output);
      output=NULL;
    }
  }

  if (par.flagIn && par.flagExt) {
      VT_WriteOverlapInterpretationVoxel2(&output, volume_bijA, volume_bijB, proportion_recouvrementA, proportion_recouvrementB, bijectionA, bijectionB, volume_totalA, volume_totalB);
      if (output) {
          fprintf(fileout, "%s", output);
          free(output);
          output=NULL;
      }
  }


  if (fileout) {
      fclose(fileout); fileout = NULL;
  }

  if (strcmp(par.fileNames.in, ">") != 0 && fp)
  {
      fclose(fp);
      fp = NULL;
  }

  /*--- liberations memoire ---*/

  if (par.flagIn) VT_FreeImage(imA);
  if (par.flagExt) VT_FreeImage(imB);

  if(volume_bijA)
    free(volume_bijA);
  if(volume_bijB)
    free(volume_bijB);
  if(proportion_recouvrementA)
    free(proportion_recouvrementA);
  if(proportion_recouvrementB)
    free(proportion_recouvrementB);

  return( 0 );
}






static void VT_Parse( int argc,
                      char *argv[],
                      local_par *par )
{
  int i, nb;
  char text[STRINGLENGTH];

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
          strcpy( par->fileNames.in, "<" );
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
        par->fileNames.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
        par->fileNames.swap = 1;
      }


      /* Parametres de calcul */

      else if ( strcmp ( argv[i], "-voxel" ) == 0 ) {
        par->mode = VOXEL;
      }
      else if ( strcmp ( argv[i], "-cell" ) == 0 ) {
        par->mode= CELL;
      }

      else if ( strcmp ( argv[i], "-image-in") == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -image-in...\n", 0 );
        strncpy( par->imageNames.in, argv[i], STRINGLENGTH );
        par->flagIn = 1;
      }
      else if ( strcmp ( argv[i], "-image-ext") == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -image-ext...\n", 0 );
        strncpy( par->imageNames.ext, argv[i], STRINGLENGTH );
        par->flagExt = 1;
      }

      else if ( strcmp ( argv[i], "-image-in-out") == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -image-in-out...\n", 0 );
        strncpy( par->imageOutNames.in, argv[i], STRINGLENGTH );
        par->flagInOut = 1;
      }
      else if ( strcmp ( argv[i], "-image-ext-out") == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -image-ext-out...\n", 0 );
        strncpy( par->imageOutNames.ext, argv[i], STRINGLENGTH );
        par->flagExtOut = 1;
      }

      /*--- option inconnue ---*/
      else {
        sprintf(text,"unknown option %s\n",argv[i]);
        VT_ErrorParse(text, 0);
      }
    }
    /*--- saisie des noms de fichiers ---*/
    else if ( argv[i][0] != 0 ) {
      if ( nb == 0 ) {
        strncpy( par->fileNames.in, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else if ( nb == 1 ) {
        strncpy( par->fileNames.out, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else
        VT_ErrorParse("too much file names when parsing\n", 0 );
    }
    i += 1;
  }

  /*--- s'il n'y a pas assez de noms ... ---*/
  if (nb == 0) {
      strcpy( par->fileNames.in,  "<" );  /* standart input */
      strcpy( par->fileNames.out, ">" );  /* standart output */
  }
  if (nb == 1) {
    strcpy( par->fileNames.out, ">" );  /* standart output */
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
  VT_Names( &(par->fileNames) );
  VT_Names( &(par->imageNames) );
  VT_Names( &(par->imageOutNames) );
  par->mode=CELL;
  par->flagIn = 0;
  par->flagExt = 0;
  par->flagInOut = 0;
  par->flagExtOut = 0;
}


