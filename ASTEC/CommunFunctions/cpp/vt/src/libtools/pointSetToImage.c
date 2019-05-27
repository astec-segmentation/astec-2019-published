/*************************************************************************
 * pointSetToImage.c -
 *
 * $Id: pointSetToImage.c,v 1.0 2016/02/17 11:31:51 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2016/02/17
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>
#include <vt_tree.h>
#include <float.h>

typedef struct local_par {
  vt_names names; /* image In / fileIn, subset file out */
  int flagImIn;
  int flagImTemplate;
  int dimx;
  int dimy;
  int dimz;
} local_par;


/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );

static int _verbose_ = 0;


static char *usage = "[image-in | file-in] [image-out] [-template | -t %s | -dim %d %d %d]\n\
\t [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si '-image-out' est '-', on prendra stdout\n\
\t -template %s : template image for output\n\
\t -dim %d %d %d : output image dimensions (default: 400x400x400) if no template image\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];

int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *imageIn;
  vt_image *imageTemplate;
  vt_image imres;
  vt_pointList *p_in_list = NULL;

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

      /*--- operations eventuelles sur l'image d'entree ---*/
      if ( par.names.inv == 1 )  VT_InverseImage( imageIn );
      if ( par.names.swap == 1 ) VT_SwapImage( imageIn );

      if (VT_AllocPointListsWithImage(&p_in_list, &nClasses, imageIn) != 1) {
          VT_FreeImage( imageIn );
          MT_ErrorParse("unable to allocate point lists image for image-in\n", 0);
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

  /* Init result image */
  if(par.flagImTemplate)
  {
      /*--- lecture de l'image template ---*/
      imageTemplate = _VT_Inrimage( par.names.ext );
      if ( imageTemplate == (vt_image*)NULL )
        MT_ErrorParse("unable to read input image ext\n", 0);
      VT_InitFromImage(&imres, imageTemplate, par.names.out, (int)UCHAR);
      if (VT_AllocImage(&imres) != 1)
      {
          for (i = 0 ; i < nClasses ; i++) {
              VT_FreePointList(&(p_in_list[i]));
          }
          free(p_in_list);
          p_in_list = NULL;
          VT_FreeImage(imageTemplate);
          VT_FreeImage(&imres);
          MT_ErrorParse("unable to allocate output image\n", 0);
      }
      imres.siz.x=imageTemplate->siz.x;
      imres.siz.y=imageTemplate->siz.y;
      imres.siz.z=imageTemplate->siz.z;
      VT_FreeImage(imageTemplate);
  }
  else
  {
      double INFINITY=DBL_MAX;
      double xmin=INFINITY,ymin=INFINITY,zmin=INFINITY;
      double xmax=-INFINITY,ymax=-INFINITY,zmax=-INFINITY;
      double sizx,sizy,sizz;
      /* set {x,y,z}{min,max} */
      for (i=0 ; i<nClasses ; i++)
      {
          vt_pointList *p = &(p_in_list[i]);
          vt_point *P;
          int j;
          for (j=0;j<p->n;j++)
          {
              P=&(p->list[j]);
              xmin=(xmin<P->x) ? xmin : P->x;
              ymin=(ymin<P->y) ? ymin : P->y;
              zmin=(zmin<P->z) ? zmin : P->z;

              xmax=(xmax>P->x) ? xmax : P->x;
              ymax=(ymax>P->y) ? ymax : P->y;
              zmax=(zmax>P->z) ? zmax : P->z;
          }

      }
      if (xmin==INFINITY || ymin==INFINITY || zmin==INFINITY ||
              xmax==-INFINITY || ymax==-INFINITY || zmax==-INFINITY )
      {
          for (i = 0 ; i < nClasses ; i++) {
              VT_FreePointList(&(p_in_list[i]));
          }
          free(p_in_list);
          p_in_list = NULL;
          MT_ErrorParse("unexpected extraction result of extremal values\n", 0);
      }
      /* trsf point-lists */
      for (i=0 ; i<nClasses ; i++)
      {
          vt_pointList *p = &(p_in_list[i]);
          vt_point *P;
          int j;
          for (j=0;j<p->n;j++)
          {
              P=&(p->list[j]);
              P->x-=xmin;
              P->y-=ymin;
              P->z-=zmin;
              ymin=(ymin<P->y) ? ymin : P->y;
              zmin=(zmin<P->z) ? zmin : P->z;

              xmax=(xmax>P->x) ? xmax : P->x;
              ymax=(ymax>P->y) ? ymax : P->y;
              zmax=(zmax>P->z) ? zmax : P->z;
          }

      }
      sizx=(xmax-xmin)/(par.dimx-1);
      sizy=(ymax-ymin)/(par.dimy-1);
      sizz=(zmax-zmin)/(par.dimz-1);

      VT_InitImage(&imres, par.names.out, par.dimx, par.dimy, par.dimz, (int)UCHAR);
      if (VT_AllocImage(&imres) != 1)
      {
          for (i = 0 ; i < nClasses ; i++) {
              VT_FreePointList(&(p_in_list[i]));
          }
          free(p_in_list);
          p_in_list = NULL;
          VT_FreeImage(&imres);
          MT_ErrorParse("unable to allocate output image\n", 0);
      }
      imres.siz.x=sizx;
      imres.siz.y=sizy;
      imres.siz.z=sizz;

  }

  for (i=0 ; i<nClasses ; i++)
      VT_SetImageFromPointList(&imres, p_in_list[i]);

  /*--- Write result ---*/
  VT_WriteInrimage(&imres);

  /* free memory*/
  VT_FreeImage(&imres);
  for (i = 0 ; i < nClasses ; i++) {
      VT_FreePointList(&(p_in_list[i]));
  }
  free(p_in_list);
  p_in_list = NULL;

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
          /*--- standart output ---*/
          strcpy( par->names.out, ">" );
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


      /* Parametres de calcul */

      else if ( strcmp ( argv[i], "-template" ) == 0 || strcmp ( argv[i], "-t" ) == 0 ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-template", 0 );
        strncpy( par->names.ext, argv[i], STRINGLENGTH );
        par->flagImTemplate=1;

      }

      else if ( strcmp ( argv[i], "-dim" ) == 0  ) {
          i ++;
          if ( i >= argc)    MT_ErrorParse( "-dim", 0 );
          status = sscanf( argv[i], "%d", &(par->dimx) );
          if ( status <= 0 ) MT_ErrorParse( "-dim", 0 );
          i ++;
          if ( i >= argc)    MT_ErrorParse( "-dim", 0 );
          status = sscanf( argv[i], "%d", &(par->dimy) );
          if ( status <= 0 ) MT_ErrorParse( "-dim", 0 );
          i ++;
          if ( i >= argc)    MT_ErrorParse( "-dim", 0 );
          status = sscanf( argv[i], "%d", &(par->dimz) );
          if ( status <= 0 ) MT_ErrorParse( "-dim", 0 );
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
      MT_ErrorParse("not enought file names when parsing\n", 0 );
  }
  if (nb == 1) {
      strcpy( par->names.out, ">" );  /* standart output */
  }


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
  par->flagImIn=1;
  par->flagImTemplate=0;
  par->dimx=300;
  par->dimy=300;
  par->dimz=300;
}

