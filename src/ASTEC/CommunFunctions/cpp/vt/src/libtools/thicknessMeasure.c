/*************************************************************************
 * thicknessMeasure.c -
 *
 * $Id: thicknessMeasure.c,v 1.0 2016/01/12 10:42:51 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2016/01/12
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>
#include <vt_tree.h>

typedef struct local_par {
  vt_names names; /* image In*/
  vt_names namesOut; /* image Out, files Out*/
  int flagImIn;
  int flagNN;
  int label1;
  int label2;
} local_par;


/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );

static int _verbose_ = 0;


static char *usage = "[image-in|file-in] [-label-1|-l1 %d] [-label-2|-l2 %d]\n\
\t [-nn|-NN %s] [-residuals-1-2|-r12 %s] [-residuals-2-1|-r21 %s]\n\
\t [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' | 'file-in' est '-', on prendra stdin\n\
\t si '-image-out' est '-', on prendra stdout\n\
\t -inv : inverse 'image-in'\n\
\t -label-1 %d : label from which thickness is measured (default = 1)\n\
\t -label-2 %d : label into which thickness is measured (default = 2)\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];

int main( int argc, char *argv[] )
{

  local_par par;
  vt_image *imageIn;
  vt_image imNearestNeighborsImage;
  vt_pointList *p_in_list = NULL;

  double sizx=0,sizy=0,sizz=0;
  int dimx=0,dimy=0,dimz=0;

  double *residuals12=NULL;
  double *residuals21=NULL;

  int i;
  int nClasses;

  FILE* fichier = NULL;

  /*--- initialisation des parametres ---*/
  MT_InitParam( &par );

  /*--- lecture des parametres ---*/
  MT_Parse( argc, argv, &par );

  if (par.flagImIn)
  {
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

      /*--- operations eventuelles sur l'image d'entree ---*/
      if ( par.names.inv == 1 )  VT_InverseImage( imageIn );
      if ( par.names.swap == 1 ) VT_SwapImage( imageIn );

      if (VT_AllocPointListsWithImage(&p_in_list, &nClasses, imageIn) != 1) {
          VT_FreeImage( imageIn );
          MT_ErrorParse("unable to allocate point lists image for image-in\n", 0);
          return(0);
      }

      if (nClasses != 3)
      {
          VT_FreeImage( imageIn );
          for (i = 0 ; i < nClasses ; i++) {
              VT_FreePointList(&(p_in_list[i]));
          }
          free(p_in_list);
          p_in_list = NULL;
          MT_ErrorParse("unexpected number of extracted classes (waiting for image of labels 0, 1 and 2)\n", 0);
          return(0);
      }

      if(VT_ExtractPointListsWithImage(*imageIn, &p_in_list, nClasses) != 1)
      {
          VT_FreeImage( imageIn );
          for (i = 0 ; i < nClasses ; i++) {
              VT_FreePointList(&(p_in_list[i]));
          }
          free(p_in_list);
          p_in_list = NULL;
          MT_ErrorParse("unable to compute point lists with image-in\n", 0);
          return(0);
      }
      /*fprintf(stdout, "Point lists extracted\n");*/

      /*--- liberations memoires partielles ---*/
      VT_FreeImage(imageIn);
  }
  else
  {
      /*--- allocation des listes de points in ---*/
      if ( VT_AllocPointListsWithFileName(&p_in_list, &nClasses, par.names.in) != 1) {
          MT_ErrorParse("unable to allocate pointset image for file-in\n", 0);
          return(0);
      }

      if (nClasses != 3)
      {
          for (i = 0 ; i < nClasses ; i++) {
              VT_FreePointList(&(p_in_list[i]));
          }
          free(p_in_list);
          p_in_list = NULL;
          MT_ErrorParse("unexpected number of extracted classes (waiting for image of labels 0, 1 and 2)\n", 0);
          return(0);
      }

      /*--- extraction des listes de points in ---*/
      if(VT_ExtractPointListsWithFileName(par.names.in, &p_in_list, nClasses) != 1)
      {
        for (i = 0 ; i < nClasses ; i++) {
            VT_FreePointList(&(p_in_list[i]));
        }
        free(p_in_list);
        p_in_list = NULL;
        MT_ErrorParse("unable to compute pointset image for file-in\n", 0);
        return(0);
      }

  }


  /*--- Construct KD-Trees for points in ---*/

  tree_node **nodes_in;
  vt_kdtrees(&nodes_in, p_in_list, nClasses);

  /*--- calculs ---*/

  int i1=par.label1, i2=par.label2;
  vt_point *p, q;
  int j;
  double d;

  /* L1 vers L2 */
  vt_pointList p_sub_2;
  if (VT_AllocPointList(&p_sub_2, p_in_list[i1].n) != 1)
  {
      for (i = 0 ; i< nClasses ; i++)
        clearTree(&(nodes_in[i]));
      free(nodes_in);
      nodes_in=NULL;
      for (i = 0 ; i < nClasses ; i++) {
          VT_FreePointList(&(p_in_list[i]));
      }
      free(p_in_list);
      p_in_list = NULL;
      MT_ErrorParse("Unable to allocate the pointList structure... \n", 0);
  }
  for (j = 0 ; j< p_in_list[i1].n ; j++)
  {
      p = &(p_in_list[i1].list[j]);
      VT_NNS(*p, nodes_in[i2], &q, &d);
      VT_SetPointListIndex(&p_sub_2, j, q);
  }

  /* L2 vers L1 */
  vt_pointList p_sub_1;
  if (VT_AllocPointList(&p_sub_1, p_in_list[i2].n) != 1)
  {
      for (i = 0 ; i< nClasses ; i++)
        clearTree(&(nodes_in[i]));
      free(nodes_in);
      nodes_in=NULL;
      for (i = 0 ; i < nClasses ; i++) {
          VT_FreePointList(&(p_in_list[i]));
      }
      free(p_in_list);
      p_in_list = NULL;
      VT_FreePointList(&p_sub_2);
      MT_ErrorParse("Unable to allocate the pointList structure... \n", 0);
  }

  for (j = 0 ; j< p_in_list[i2].n ; j++)
  {
      p = &(p_in_list[i2].list[j]);
      VT_NNS(*p, nodes_in[i1], &q, &d);
      VT_SetPointListIndex(&p_sub_1, j, q);
  }

  /* Nearest Neighbors image initialisation + allocation */
  if (par.flagNN)
  {
    unsigned char ***theNNarray;
    VT_InitImage(&imNearestNeighborsImage,par.namesOut.in,dimx,dimy,dimz,(int)UCHAR);
    if (VT_AllocImage(&imNearestNeighborsImage) != 1)
    {
        for (i = 0 ; i< nClasses ; i++)
          clearTree(&(nodes_in[i]));
        free(nodes_in);
        nodes_in=NULL;
        for (i = 0 ; i < nClasses ; i++) {
            VT_FreePointList(&(p_in_list[i]));
        }
        free(p_in_list);
        p_in_list = NULL;
        VT_FreePointList(&(p_sub_1));
        VT_FreePointList(&(p_sub_2));
        MT_ErrorParse("Unable to allocate an output image, exiting...\n",0);
    }
    imNearestNeighborsImage.siz.x=sizx;
    imNearestNeighborsImage.siz.y=sizy;
    imNearestNeighborsImage.siz.z=sizz;
    theNNarray=(unsigned char ***)imNearestNeighborsImage.array;

    /* Image setting */

    for (i=0; i<p_sub_2.n; i++)
    {
        int indx,indy,indz;
        indx=(int)((p_sub_2.list[i].x)/sizx);
        indy=(int)((p_sub_2.list[i].y)/sizy);
        indz=(int)((p_sub_2.list[i].z)/sizz);
        theNNarray[indz][indy][indx]=(unsigned char) i2;
    }

    for (i=0; i<p_sub_1.n; i++)
    {
        int indx,indy,indz;
        indx=(int)((p_sub_1.list[i].x)/sizx);
        indy=(int)((p_sub_1.list[i].y)/sizy);
        indz=(int)((p_sub_1.list[i].z)/sizz);
        theNNarray[indz][indy][indx]=(unsigned char) i1;
    }

    /* Image writing */

    VT_WriteInrimage(&imNearestNeighborsImage);
    VT_FreeImage(&imNearestNeighborsImage);
  }

  /*--- Residual vectors ---*/
  /*--- label1 -> label 2 ---*/
  if (VT_ResidualsVector(p_in_list[i1], nodes_in[i2], &residuals12) != 1)
  {
    for (i = 0 ; i< nClasses ; i++)
      clearTree(&(nodes_in[i]));
    free(nodes_in);
    nodes_in=NULL;
    for (i = 0 ; i < nClasses ; i++) {
        VT_FreePointList(&(p_in_list[i]));
    }
    free(p_in_list);
    p_in_list = NULL;
    VT_FreePointList(&(p_sub_1));
    VT_FreePointList(&(p_sub_2));
    MT_ErrorParse("Error while computing residuals from label 1 to label 2\n", 0);
  }
  /* output text file */
  fichier = NULL;
  if(par.namesOut.ext[0] != '\0') /* residuals12*/
  {
    int *ind;
    /* residuals vector sorting */
    if (VT_TriRapideVector(residuals12, p_in_list[i1].n, &ind) != 1)
    {
      free(residuals12); residuals12 = NULL;
      for (i = 0 ; i< nClasses ; i++)
        clearTree(&(nodes_in[i]));
      free(nodes_in);
      nodes_in=NULL;
      for (i = 0 ; i < nClasses ; i++) {
          VT_FreePointList(&(p_in_list[i]));
      }
      free(p_in_list);
      p_in_list = NULL;
      VT_FreePointList(&(p_sub_1));
      VT_FreePointList(&(p_sub_2));
      MT_ErrorParse("Error while sorting residuals from label 1 to label 2\n", 0);
    }
    fichier = fopen(par.namesOut.ext, "w");
    if( fichier == NULL)
    {
        free(ind);  ind = NULL;
        free(residuals12); residuals12 = NULL;
        for (i = 0 ; i< nClasses ; i++)
          clearTree(&(nodes_in[i]));
        free(nodes_in);
        nodes_in=NULL;
        for (i = 0 ; i < nClasses ; i++) {
            VT_FreePointList(&(p_in_list[i]));
        }
        free(p_in_list);
        p_in_list = NULL;
        VT_FreePointList(&(p_sub_1));
        VT_FreePointList(&(p_sub_2));
        MT_ErrorParse("Erreur pour l'ecriture des residus dans un fichier\n", 0);
    }
    for (i = 0 ; i<p_in_list[i1].n ; i++)
        fprintf(fichier, "%f\n", residuals12[ind[i]]);
    fclose(fichier);
    free(ind); ind = NULL;
  }

  /*--- label2 -> label 1 ---*/
  if (VT_ResidualsVector(p_in_list[i2], nodes_in[i1], &residuals21) != 1)
  {
    free(residuals12); residuals12 = NULL;
    for (i = 0 ; i< nClasses ; i++)
      clearTree(&(nodes_in[i]));
    free(nodes_in);
    nodes_in=NULL;
    for (i = 0 ; i < nClasses ; i++) {
        VT_FreePointList(&(p_in_list[i]));
    }
    free(p_in_list);
    p_in_list = NULL;
    VT_FreePointList(&(p_sub_1));
    VT_FreePointList(&(p_sub_2));
    MT_ErrorParse("Error while computing residuals from label 1 to label 2\n", 0);
  }

  /* output text file */
  fichier = NULL;
  if(par.namesOut.out[0] != '\0') /* residuals21*/
  {
    int *ind;
    /* residuals vector sorting */
    if (VT_TriRapideVector(residuals21, p_in_list[i2].n, &ind) != 1)
    {
      free(residuals12); residuals12 = NULL;
      free(residuals21); residuals21 = NULL;
      for (i = 0 ; i< nClasses ; i++)
        clearTree(&(nodes_in[i]));
      free(nodes_in);
      nodes_in=NULL;
      for (i = 0 ; i < nClasses ; i++) {
          VT_FreePointList(&(p_in_list[i]));
      }
      free(p_in_list);
      p_in_list = NULL;
      VT_FreePointList(&(p_sub_1));
      VT_FreePointList(&(p_sub_2));
      MT_ErrorParse("Error while sorting residuals from label 2 to label 1\n", 0);
    }
    fichier = fopen(par.namesOut.out, "w");
    if( fichier == NULL)
    {
        free(ind);  ind = NULL;
        free(residuals12); residuals12 = NULL;
        free(residuals21); residuals21 = NULL;
        for (i = 0 ; i< nClasses ; i++)
          clearTree(&(nodes_in[i]));
        free(nodes_in);
        nodes_in=NULL;
        for (i = 0 ; i < nClasses ; i++) {
            VT_FreePointList(&(p_in_list[i]));
        }
        free(p_in_list);
        p_in_list = NULL;
        VT_FreePointList(&(p_sub_1));
        VT_FreePointList(&(p_sub_2));
        MT_ErrorParse("Erreur pour l'ecriture des residus dans un fichier\n", 0);
    }
    for (i = 0 ; i<p_in_list[i2].n ; i++)
        fprintf(fichier, "%f\n", residuals21[ind[i]]);
    fclose(fichier);
    free(ind); ind = NULL;
  }


  /* free memory*/

  for (i = 0 ; i< nClasses ; i++)
    clearTree(&(nodes_in[i]));
  free(nodes_in);
  nodes_in=NULL;
  for (i = 0 ; i < nClasses ; i++) {
      VT_FreePointList(&(p_in_list[i]));
  }
  free(p_in_list);
  p_in_list = NULL;
  VT_FreePointList(&(p_sub_1));
  VT_FreePointList(&(p_sub_2));
  free(residuals12); residuals12=NULL;
  free(residuals21); residuals21=NULL;

  return( 0 );

}



static void MT_Parse( int argc,
                      char *argv[],
                      local_par *par )
{
  int i, nb;
  char text[STRINGLENGTH];
  int status;

  if ( VT_CopyName( program, argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);
  if ( argc == 1 ) MT_ErrorParse("\n", 0 );

  /*--- lecture des parametres ---*/
  i = 1;
  nb = 0;
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

      else if ( strcmp ( argv[i], "-in-points" ) == 0 || strcmp ( argv[i], "-flo-points" ) == 0 ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-flo-points", 0 );
        strncpy( par->names.in, argv[i], STRINGLENGTH );
        par->flagImIn=0;
        nb++;
      }

      else if ( strcmp ( argv[i], "-nn" ) == 0 || strcmp ( argv[i], "-NN" ) == 0 ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-NN", 0 );
        strncpy( par->namesOut.in, argv[i], STRINGLENGTH );
        par->flagNN=1;
        nb++;
      }

      else if ( strcmp ( argv[i], "-residuals-1-2" ) == 0 || strcmp ( argv[i], "-r12" ) == 0 ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-residuals-1-2", 0 );
        strncpy( par->namesOut.ext, argv[i], STRINGLENGTH );
        nb++;
      }

      else if ( strcmp ( argv[i], "-residuals-2-1" ) == 0 || strcmp ( argv[i], "-r21" ) == 0 ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-residuals-2-1", 0 );
        strncpy( par->namesOut.out, argv[i], STRINGLENGTH );
        nb++;
      }

      /* Parametres de calcul */


      else if ( strcmp ( argv[i], "-label-1" ) == 0 || strcmp ( argv[i], "-l1" ) == 0 ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-label-1", 0 );
        status = sscanf( argv[i], "%d", &(par->label1) );
        if ( status <= 0 ) MT_ErrorParse( "-label-1", 0 );
      }
      else if ( strcmp ( argv[i], "-label-2" ) == 0 || strcmp ( argv[i], "-l2" ) == 0 ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-label-2", 0 );
        status = sscanf( argv[i], "%d", &(par->label2) );
        if ( status <= 0 ) MT_ErrorParse( "-label-2", 0 );
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
          strncpy( par->names.out, argv[i], STRINGLENGTH );
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
      strcpy( par->names.in, ">" );  /* standart input */
      strcpy( par->names.out, "<" );  /* standart output */
  }
  if (nb == 1)
    strcpy( par->names.out, "<" );  /* standart output */

  char *point;
  if((point = strrchr(par->names.in,'.')) != NULL ) {
     if(strcmp(point,".txt") == 0) {
         par->flagImIn=0;
     }
  }

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
  par->flagImIn=1;
  par->flagNN=0;
  par->label1=1;
  par->label2=2;
}

