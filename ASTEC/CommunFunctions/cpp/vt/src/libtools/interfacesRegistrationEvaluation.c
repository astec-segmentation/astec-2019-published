/*************************************************************************
 * interfacesRegistrationEvaluation.c -
 *
 * $Id: interfacesRegistrationEvaluation.c,v 1.0 2016/01/15 11:34:51 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2016/01/15
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>
#include <vt_tree.h>
#include <bal-transformation.h>

typedef struct local_par {
  vt_names names; /* image In, image Ext, trsf */
  vt_names namesOut; /* image Out, files Out */
  int flagImIn;
  int flagImExt;
  int label1;
  int label2;
  int flagInitTrsf;
} local_par;


/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );

static int _verbose_ = 0;


static char *usage = "[image-in] [image-ext] [-trsf %s] [-label-1|-l1 %d] [-label-2|-l2 %d]\n\
\t \n\
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

  vt_pointList *p_in_list = NULL;
  /*vt_image *imageExt;
  vt_pointList *p_ext_list = NULL;
  vt_pointList *p_ext_list_trsf = NULL;

  double T_ext_in[16];
  */
  double sizx=0,sizy=0,sizz=0;
  int dimx=0,dimy=0,dimz=0;

  /*double *residuals12=NULL;
  //double *residuals21=NULL;
  */
  int i;
  int nClasses;

  /*FILE* fichier = NULL;*/

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

  /*
  if (par.flagImExt)
  {
      // lecture de l'image d'entree ext
      imageExt = _VT_Inrimage( par.names.ext );
      if ( imageExt == (vt_image*)NULL )
        MT_ErrorParse("unable to read input image ext\n", 0);

      // operations eventuelles sur l'image d'entree
      if ( par.names.inv == 1 )  VT_InverseImage( imageIn );
      if ( par.names.swap == 1 ) VT_SwapImage( imageIn );

      if (VT_AllocPointListsWithImage(&p_ext_list, &nClasses, imageExt) != 1) {
          VT_FreeImage( imageExt );
          MT_ErrorParse("unable to allocate point lists image for image-ext\n", 0);
          return(0);
      }

      if (nClasses != 3)
      {
          VT_FreeImage( imageExt );
          for (i = 0 ; i < 3 ; i++) {
              VT_FreePointList(&(p_in_list[i]));
          }
          free(p_in_list);
          p_in_list = NULL;
          for (i = 0 ; i < nClasses ; i++) {
              VT_FreePointList(&(p_ext_list[i]));
          }
          free(p_ext_list);
          p_ext_list = NULL;
          MT_ErrorParse("unexpected number of extracted classes (waiting for image of labels 0, 1 and 2)\n", 0);
          return(0);
      }

      if(VT_ExtractPointListsWithImage(*imageExt, &p_ext_list, nClasses) != 1)
      {
          VT_FreeImage( imageExt );
          for (i = 0 ; i < 3 ; i++) {
              VT_FreePointList(&(p_in_list[i]));
          }
          free(p_in_list);
          p_in_list = NULL;
          for (i = 0 ; i < nClasses ; i++) {
              VT_FreePointList(&(p_ext_list[i]));
          }
          free(p_ext_list);
          p_ext_list = NULL;
          MT_ErrorParse("unable to compute point lists with image-in\n", 0);
          return(0);
      }

      // liberations memoires partielles
      VT_FreeImage(imageExt);
  }
  else
  {
      // allocation des listes de points ext
      if ( VT_AllocPointListsWithFileName(&p_ext_list, &nClasses, par.names.ext) != 1) {
          MT_ErrorParse("unable to allocate pointset image for file-ext\n", 0);
          return(0);
      }

      if (nClasses != 3)
      {
          for (i = 0 ; i < 3 ; i++) {
              VT_FreePointList(&(p_in_list[i]));
          }
          free(p_in_list);
          p_in_list = NULL;
          for (i = 0 ; i < nClasses ; i++) {
              VT_FreePointList(&(p_ext_list[i]));
          }
          free(p_ext_list);
          p_ext_list = NULL;
          MT_ErrorParse("unexpected number of extracted classes (waiting for image of labels 0, 1 and 2)\n", 0);
          return(0);
      }

      // extraction des listes de points ext
      if(VT_ExtractPointListsWithFileName(par.names.ext, &p_ext_list, nClasses) != 1)
      {
        for (i = 0 ; i < 3 ; i++) {
            VT_FreePointList(&(p_in_list[i]));
        }
        free(p_in_list);
        p_in_list = NULL;
        for (i = 0 ; i < nClasses ; i++) {
            VT_FreePointList(&(p_ext_list[i]));
        }
        free(p_ext_list);
        p_ext_list = NULL;
        MT_ErrorParse("unable to compute pointset image for file-in\n", 0);
        return(0);
      }

  }
  bal_transformation theTrsf;
  if(par.flagInitTrsf)
  {
      BAL_InitTransformation( &theTrsf );
      // Lecture trsf initiale
      if ( BAL_ReadTransformation( &theTrsf, par.names.out ) != 1 ) {
        for (i = 0 ; i < 3 ; i++) {
            VT_FreePointList(&(p_in_list[i]));
        }
        free(p_in_list);
        p_in_list = NULL;
        for (i = 0 ; i < nClasses ; i++) {
            VT_FreePointList(&(p_ext_list[i]));
        }
        free(p_ext_list);
        p_ext_list = NULL;
        fprintf( stderr, "%s: unable to read 'real' transformation '%s'\n", program, par.names.out );
        return (-1);
      }
      double *m = theTrsf.mat.m;
      for(i=0; i<16; i++)
        T_ext_in[i] = m[i];
  }
  else
  {
      for (i=0; i<16;i++)
        T_ext_in[i]=0.0;
      for(i=0;i<4;i++)
        T_ext_in[i*5]=1.0;
  }
    */

  /*--- Points ext transformation ---*/
  /*if (VT_AllocPointListsFromPointLists(&p_ext_list_trsf, nClasses, p_ext_list) != 1)
  {
      for (i = 0 ; i < 3 ; i++) {
          VT_FreePointList(&(p_in_list[i]));
      }
      free(p_in_list);
      p_in_list = NULL;
      for (i = 0 ; i < nClasses ; i++) {
          VT_FreePointList(&(p_ext_list[i]));
      }
      free(p_ext_list);
      p_ext_list = NULL;
      fprintf( stderr, "%s: unable to read 'real' transformation '%s'\n", program, par.names.out );
      return (-1);
  }
  for (i=0; i<nClasses ; i++)
    VT_TransformPointList(T_ext_in, p_ext_list[i], &(p_ext_list_trsf[i]));
  */

  /*--- Construct KD-Trees for points in ---*/

  tree_node **nodes_in;
  vt_kdtrees(&nodes_in, p_in_list, nClasses);
  /*tree_node **nodes_ext_trsf;
  //vt_kdtrees(&nodes_ext_trsf, p_ext_list_trsf, nClasses);
  */
  /*--- calculs ---*/

  int i1=par.label1, i2=par.label2;
  vt_point *p, q;
  int j;
  double d;

  /*---- Image IN ----*/
  /* L1 vers L2 */
  vt_coupleList c_in_12;
  if (VT_AllocCoupleList(&c_in_12, p_in_list[i1].n) != 1)
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
      /*for (i = 0 ; i< nClasses ; i++)
        clearTree(&(nodes_ext_trsf[i]));
      free(nodes_ext_trsf);
      nodes_ext_trsf=NULL;
      for (i = 0 ; i < nClasses ; i++) {
          VT_FreePointList(&(p_ext_list[i]));
      }
      free(p_ext_list);
      p_ext_list = NULL;
      for (i = 0 ; i < nClasses ; i++) {
          VT_FreePointList(&(p_ext_list_trsf[i]));
      }
      free(p_ext_list_trsf);
      p_ext_list_trsf = NULL;
      */
      MT_ErrorParse("Unable to allocate the coupleList structure... \n", 0);
  }
  for (j = 0 ; j< p_in_list[i1].n ; j++)
  {
      p = &(p_in_list[i1].list[j]);
      VT_NNS(*p, nodes_in[i2], &q, &d);
      VT_CopyPoint(&(c_in_12.list[j].start), *p);
      VT_CopyPoint(&(c_in_12.list[j].end), q);

  }


  /* L2 vers L1 */
  vt_coupleList c_in_21;
  if (VT_AllocCoupleList(&c_in_21, p_in_list[i2].n) != 1)
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
      /*for (i = 0 ; i< nClasses ; i++)
        clearTree(&(nodes_ext_trsf[i]));
      free(nodes_ext_trsf);
      nodes_ext_trsf=NULL;
      for (i = 0 ; i < nClasses ; i++) {
          VT_FreePointList(&(p_ext_list[i]));
      }
      free(p_ext_list);
      p_ext_list = NULL;
      for (i = 0 ; i < nClasses ; i++) {
          VT_FreePointList(&(p_ext_list_trsf[i]));
      }
      free(p_ext_list_trsf);
      p_ext_list_trsf = NULL;
      */
      VT_FreeCoupleList(&c_in_12);

      MT_ErrorParse("Unable to allocate the coupleList structure... \n", 0);
  }
  for (j = 0 ; j< p_in_list[i2].n ; j++)
  {
      p = &(p_in_list[i2].list[j]);
      VT_NNS(*p, nodes_in[i1], &q, &d);
      VT_CopyPoint(&(c_in_21.list[j].start), *p);
      VT_CopyPoint(&(c_in_21.list[j].end), q);
  }

  /* Write Couple Lists */
  if ( par.namesOut.in[0] != '\0' )
      VT_WriteCoupleList( &c_in_12, par.namesOut.in );
  if ( par.namesOut.ext[0] != '\0' )
      VT_WriteCoupleList( &c_in_21, par.namesOut.ext );

  /* TESTS */
  vt_coupleList c;
  /*VT_AllocCoupleListFromFileName(&c, par.namesOut.in);
  //VT_ExtractCoupleListFromFileName(&c, par.namesOut.in);
  //if ( par.namesOut.out != '\0' )
  //    VT_WriteCoupleList( &c, par.namesOut.out );
  //VT_FreeCoupleList(&c);

  //vt_pointList ptest;
  //fprintf(stdout, "oucou ; c_in_12.n = %d\n", c_in_12.n);
  //VT_CoupleListToMedianList(c_in_12, &ptest);
  //for (i=0;i<10;i++)
  //{
  //    printPoint(ptest.list[i]);
  //}
  //VT_FreePointList(&ptest);
  */
  VT_ExtractBijectiveCoupleList(&c, c_in_12, c_in_21, 0);
  if ( par.namesOut.out[0] != '\0' )
      VT_WriteCoupleList( &c, par.namesOut.out );

  vt_image imOut;
  VT_InitImage(&imOut, "fooImg.hdr", dimx, dimy, dimz, (int)UCHAR);
  imOut.siz.x=sizx;
  imOut.siz.y=sizy;
  imOut.siz.z=sizz;
  if (VT_AllocImage(&imOut) != 1)
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
      VT_FreeCoupleList(&c_in_12);
      VT_FreeCoupleList(&c_in_21);
      VT_FreeCoupleList(&c);
      fprintf(stderr, "Unable to allocate output image\n");
      return(-1);
  }
  if ( VT_SetImageFromCoupleList(&imOut, c) != 1)
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
      VT_FreeCoupleList(&c_in_12);
      VT_FreeCoupleList(&c_in_21);
      VT_FreeCoupleList(&c);
      VT_FreeImage(&imOut);
      fprintf(stderr, "Unable to set output image from couple list\n");
      return(-1);
  }
  VT_WriteInrimage(&imOut);
  VT_FreeImage(&imOut);
  VT_WriteCoupleListDistanceHistogram(c,"fooHist.txt");
  VT_FreeCoupleList(&c);





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
  VT_FreeCoupleList(&c_in_12);
  VT_FreeCoupleList(&c_in_21);
  /*free(residuals12); residuals12=NULL;
  //free(residuals21); residuals21=NULL;
  */

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
        if ( nb == 1 ) {
          /*--- standart input ---*/
          strcpy( par->names.ext, "<" );
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
      else if ( strcmp ( argv[i], "-in-points" ) == 0 || strcmp ( argv[i], "-flo-points" ) == 0 ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-ext-points", 0 );
        strncpy( par->names.ext, argv[i], STRINGLENGTH );
        par->flagImExt=0;
        nb++;
      }

      /*else if ( strcmp ( argv[i], "-nn" ) == 0 || strcmp ( argv[i], "-NN" ) == 0 ) {
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
      }*/

      /* Parametres de calcul */

      else if ( strcmp ( argv[i], "-trsf" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -trsf...\n",0);
        sprintf(par->names.out, "%s", argv[i]);
        par->flagInitTrsf=1;
      }

      else if ( strcmp ( argv[i], "-tmp1" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -tmp1...\n",0);
        sprintf(par->namesOut.in, "%s", argv[i]);
      }

      else if ( strcmp ( argv[i], "-tmp2" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -tmp2...\n",0);
        sprintf(par->namesOut.ext, "%s", argv[i]);
      }

      else if ( strcmp ( argv[i], "-tmp3" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -tmp3...\n",0);
        sprintf(par->namesOut.out, "%s", argv[i]);
      }

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
          strncpy( par->names.ext, argv[i], STRINGLENGTH );
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
      MT_ErrorParse("not enought file names when parsing\n", 0 );
  }
  if (nb == 1) {
      strcpy( par->names.ext, ">" );  /* standart input */     
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
  par->flagImExt=1;
  par->label1=1;
  par->label2=2;
  par->flagInitTrsf=0;
}

