/*************************************************************************
 * pointSetDistances.c -
 *
 * $Id: pointSetDistances.c,v 1.0 2015/12/11 09:45:51 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2015/12/11
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>
/*#include <vt_planeRegistration.h>*/
#include <vt_tree.h>

typedef struct local_par {
  vt_names names; /*image In, image Ext*/
  vt_names namesOut; /*image Out, file Out*/
  int flagFileOut;
  /*int flagImOut;
  //int flagImIn;
  //int flagImExt;
  */
} local_par;


/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );

static int _verbose_ = 0;


static char *usage = "[image-in] [image-ext] [image-out] [-file-out %s]\n\
\t [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si '-image-out' est '-', on prendra stdout\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];

int main( int argc, char *argv[] )
{

  local_par par;
  vt_image *imageIn;
  vt_image *imageExt;
  vt_image imResiduals;
  vt_pointList p_in;
  vt_pointList p_ext;

  double sizx=0,sizy=0,sizz=0;
  int dimx=0,dimy=0,dimz=0;

  double *residuals=NULL;

  int i;

  FILE* fichier = NULL;

  /*--- initialisation des parametres ---*/
  MT_InitParam( &par );

  /*--- lecture des parametres ---*/
  MT_Parse( argc, argv, &par );

  /*--- lecture de l'image d'entree in ---*/
  imageIn = _VT_Inrimage( par.names.in );
  if ( imageIn == (vt_image*)NULL )
    MT_ErrorParse("unable to read input image in\n", 0);
  sizx=imageIn->siz.x;
  sizy=imageIn->siz.y;
  sizz=imageIn->siz.z;
  dimx=imageIn->dim.x;
  dimy=imageIn->dim.y;
  dimz=imageIn->dim.z;

  /*fprintf(stdout, "%d %d %d \n%f %f %f\n", dimx, dimy, dimz,
  //        sizx, sizy, sizz);
  */

  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.names.inv == 1 )  VT_InverseImage( imageIn );
  if ( par.names.swap == 1 ) VT_SwapImage( imageIn );

  if (VT_AllocPointListWithImage(&p_in, imageIn) != 1) {
      VT_FreeImage( imageIn );
      VT_Free( (void**)&imageIn);
      MT_ErrorParse("unable to allocate pointset image for image-in\n", 0);
      return(0);
  }

  if(VT_ExtractPointListWithImage(*imageIn, &p_in) != 1)
  {
      VT_FreeImage( imageIn );
      VT_Free( (void**)&imageIn);
      VT_FreePointList(&p_in);
      MT_ErrorParse("unable to compute pointset image for image-in\n", 0);
      return(0);
  }

  /*--- liberations memoires partielles ---*/
  VT_FreeImage(imageIn);

  /*--- lecture de l'image d'entree ext ---*/
  imageExt = _VT_Inrimage( par.names.ext );
  if ( imageExt == (vt_image*)NULL )
    MT_ErrorParse("unable to read input image ext\n", 0);

  if ((int)imageExt->dim.x != dimx || (int)imageExt->dim.y != dimy || (int)imageExt->dim.z != dimz ||
          imageExt->siz.x != sizx || imageExt->siz.y != sizy || imageExt->siz.z != sizz )
  {
      fprintf(stderr, "Warning: input images do not have the same dimensions / resolutions...\n");
  }

  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.names.inv == 1 )  VT_InverseImage( imageExt );
  if ( par.names.swap == 1 ) VT_SwapImage( imageExt );

  if (VT_AllocPointListWithImage(&p_ext, imageExt) != 1) {
      VT_FreeImage( imageExt );
      VT_Free( (void**)&imageExt);
      VT_FreePointList(&p_in);
      MT_ErrorParse("unable to allocate pointset image for image-ext\n", 0);
      return(0);
  }

  if(VT_ExtractPointListWithImage(*imageExt, &p_ext) != 1)
  {
      VT_FreeImage( imageExt );
      VT_Free( (void**)&imageExt);
      VT_FreePointList(&p_in);
      VT_FreePointList(&p_ext);
      MT_ErrorParse("unable to compute pointset image for image-ext\n", 0);
      return(0);
  }

  /*--- liberations memoires partielles ---*/
  VT_FreeImage(imageExt);

  /*--- Construct KD-Tree for points ext ---*/

  tree_node *node_ext;
  vt_kdtree(&node_ext, p_ext, 0);

  /*--- calculs ---*/

  if (VT_ResidualsVector(p_in, node_ext, &residuals) != 1)
  {
    if(p_in.list)
      VT_FreePointList(&p_in);
    if (p_ext.list)
      VT_FreePointList(&p_ext);
    clearTree(&node_ext);
    MT_ErrorParse("Error while computing residuals\n", 0);
  }

  /* Image of residuals */

  /*--- init et allocation image de sortie ---*/

  VT_InitImage( &imResiduals, par.namesOut.out, dimx, dimy,
                dimz, FLOAT );

  if ( VT_AllocImage( &(imResiduals) ) != 1 ) {
    if(p_in.list)
      VT_FreePointList(&p_in);
    if (p_ext.list)
      VT_FreePointList(&p_ext);
    clearTree(&node_ext);
    MT_ErrorParse("problem while allocating imResiduals\n", 0 );
  }
  imResiduals.siz.x = sizx;
  imResiduals.siz.y = sizy;
  imResiduals.siz.z = sizz;

  fprintf(stdout, "%d %d %d \n%f %f %f\n", (int)imResiduals.dim.x, (int)imResiduals.dim.y, (int)imResiduals.dim.z,
          imResiduals.siz.x, imResiduals.siz.y, imResiduals.siz.z);

  VT_SetResidualImage(&imResiduals, p_in, residuals);


  VT_WriteInrimage(&imResiduals);
  VT_FreeImage(&imResiduals);

  /* output text file */
  fichier = NULL;
  if(par.namesOut.ext[0] != '\0') /* residuals */
  {
    int *ind;

    /* residuals vector sorting */

    if (VT_TriRapideVector(residuals, p_in.n, &ind) != 1)
    {
      free(residuals); residuals = NULL;
      if(p_in.list)
        VT_FreePointList(&p_in);
      if (p_ext.list)
        VT_FreePointList(&p_ext);
      clearTree(&node_ext);
      MT_ErrorParse("Error while sorting residuals\n", 0);
    }

    /* sorted residuals writting in output txt file */

      fichier = fopen(par.namesOut.ext, "w");
      if( fichier == NULL)
      {
        free(residuals); residuals = NULL;
        if(p_in.list)
          VT_FreePointList(&p_in);
        if (p_ext.list)
          VT_FreePointList(&p_ext);
        clearTree(&node_ext);
        free(ind); ind=NULL;
        MT_ErrorParse("Erreur pour l'ecriture des residus dans un fichier\n", 0);
      }
      for (i = 0 ; i<p_in.n ; i++)
        fprintf(fichier, "%f\n", residuals[ind[i]]);
      fclose(fichier);
      free(ind); ind = NULL;
  }


  /* free memory*/

  free(residuals); residuals = NULL;
  if(p_in.list)
    VT_FreePointList(&p_in);
  if (p_ext.list)
    VT_FreePointList(&p_ext);
  clearTree(&node_ext);

  return( 0 );

}



static void MT_Parse( int argc,
                      char *argv[],
                      local_par *par )
{
  int i, nb;
  char text[STRINGLENGTH];
  /*int status;*/

  if ( VT_CopyName( program, argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);
  if ( argc == 1 ) MT_ErrorParse("\n", 0 );

  /*--- lecture des parametres ---*/
  i = 1;
  nb = 0;
  while ( i < argc ) {
    if ( argv[i][0] == '-' ) {

      /*--- arguments generaux ---*/
      if ( strcmp ( argv[i], "-help" ) == 0 ) {
        MT_ErrorParse("\n", 1);
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

      else if ( strcmp ( argv[i], "-residuals-out" ) == 0 || strcmp ( argv[i], "-residuals" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -residuals...\n",0);
        sprintf(par->namesOut.ext, "%s", argv[i]);
        par->flagFileOut=1;
      }



      /*--- option inconnue ---*/
      else {
        sprintf(text,"unknown option %s\n",argv[i]);
        MT_ErrorParse(text, 0);
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
        strncpy( par->namesOut.out, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else {
        MT_ErrorParse("too much file names when parsing\n", 0 );
      }
    }
    i += 1;
  }

  /*--- s'il n'y a pas assez de noms ... ---*/
  if (nb == 0) {
      MT_ErrorParse("not enough file names when parsing\n", 0 );
  }
  if (nb == 1)
    strcpy( par->names.ext, ">" );  /* standart ext */


}




static void MT_ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"Erreur : %s",str);
  exit( 1 );
}








static void MT_InitParam( local_par *par )
{
  VT_Names( &(par->names) );
  VT_Names( &(par->namesOut) );
  par->flagFileOut=0;
}

