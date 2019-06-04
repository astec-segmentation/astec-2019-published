/*************************************************************************
 * membrane.c -
 *
 * $Id: membrane.c,v 1.0 2013/08/12 10:34:51 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2013/08/12
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>

typedef struct local_par {
  vt_names names;
  int oblic;
} local_par;

typedef struct {
  double x;
  double y;
} point_t;

typedef point_t* point_ptr_t;

/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );
/*
 * static double ccw(point_t* p1, point_t* p2, point_t* p3);
 * static void convex_hull(point_t* points, ssize_t npoints, point_ptr_t** out_hull, ssize_t* out_hullsize);
 */

static int _verbose_ = 0;





static char *usage = "[image-in] [image-out] \n\
\t [-inv] [-swap] [-v] [-D] [-help]";
/* [-oblic] */

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];







int main( int argc, char *argv[] )
{

  local_par par;
  vt_image *image;
  vt_image imres;
  /*  int flag_3D = 1; */
  char name[DOUBLESTRINGLENGTH];
  sprintf( name, "%s", par.names.out );


  /*--- initialisation des parametres ---*/
  MT_InitParam( &par );

  /*--- lecture des parametres ---*/
  MT_Parse( argc, argv, &par );

  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL )
    MT_ErrorParse("unable to read input image\n", 0);


  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.names.inv == 1 )  VT_InverseImage( image );
  if ( par.names.swap == 1 ) VT_SwapImage( image );


  /* if ( image->dim.z == 1 )
   *  flag_3D = 0;
   */



  /* Alloc image sortie */

  int t = (int) UCHAR;
  sprintf( name, "%s", par.names.out );
  VT_InitImage( &imres, name, image->dim.x, image->dim.y,
              image->dim.z, t );
  if ( VT_AllocImage( &(imres) ) != 1 ) {
    MT_ErrorParse("problem while allocating imres\n", 0 );
    VT_FreeImage( image );
    return( -1 );
  }

  imres.siz.x=image->siz.x;
  imres.siz.y=image->siz.y;
  imres.siz.z=image->siz.z;

  /*--- calculs ---*/

  unsigned short int ***imgarray=NULL;

  unsigned char ***imgarrayu8=NULL;
  unsigned char ***resarray;

  switch (image->type) {
    case SCHAR:
    case UCHAR:
      imgarrayu8= (unsigned char ***)image->array;
      break;
    case SSHORT:
    case USHORT:
      imgarray= (unsigned short int ***)image->array;
      break;
    case TYPE_UNKNOWN:
    default:
      VT_Error("image type unknown or not supported for this program",program);
      return( 0 );
  }

  resarray= (unsigned char ***)imres.array;

  int z,y,x, Z,Y,X;
  int val;

  for (z=0 ; z<(int)image->dim.z ; z++)
  for (y=0 ; y<(int)image->dim.y ; y++)
  for (x=0 ; x<(int)image->dim.x ; x++) {
      resarray[z][y][x]=(unsigned char) 0;
  }

  if(par.oblic == 0 )
  {
    if (_verbose_)
      fprintf(stdout, "Petit Z...\n");
    /* Petit Z */
    val=0;
    for (z=0 ; z<(int)image->dim.z ; z++)
    {
      for (y=0 ; y<(int)image->dim.y ; y++) {
        for (x=0 ; x<(int)image->dim.x ; x++) {
          switch(image->type) {
          case SCHAR:
          case UCHAR:
              val = (int) imgarrayu8[z][y][x];
              break;
          case SSHORT:
          case USHORT:
              val = (int) imgarray[z][y][x];
              break;
          default:
              VT_Error("image type unknown or not supported for this program",program);
              return( 0 );
          }
          if(val!=0) break;
        }
        if(val!=0) break;
      }
      if(val!=0) break;
    }

    for (Z=0 ; Z<z ; Z++)
    for (y=0 ; y<(int)image->dim.y ; y++)
    for (x=0 ; x<(int)image->dim.x ; x++) {
      resarray[Z][y][x]=(unsigned char) 255;
    }

    if (_verbose_)
      fprintf(stdout, "Grand Z...\n");
    /* Grand Z */
    val=0;
    for (z=image->dim.z-1 ; z>=0 ; z--)
    {
      for (y=0 ; y<(int)image->dim.y ; y++) {
        for (x=0 ; x<(int)image->dim.x ; x++) {
          switch(image->type) {
          case SCHAR:
          case UCHAR:
              val = (int) imgarrayu8[z][y][x];
              break;
          case SSHORT:
          case USHORT:
              val = (int) imgarray[z][y][x];
              break;
          default:
              VT_Error("image type unknown or not supported for this program",program);
              return( 0 );
          }
          if(val!=0) break;
        }
        if(val!=0) break;
      }
      if(val!=0) break;
    }

    for (Z=z+1 ; Z<(int)image->dim.z ; Z++)
    for (y=0 ; y<(int)image->dim.y ; y++)
    for (x=0 ; x<(int)image->dim.x ; x++) {
      resarray[Z][y][x]=(unsigned char) 255;
    }

    if (_verbose_)
      fprintf(stdout, "Petit Y...\n");
    /* Petit Y sur chaque z */
    for (z=0 ; z<(int)image->dim.z ; z++)
    {
      val=0;
      for (y=0 ; y<(int)image->dim.y ; y++) {
        for (x=0 ; x<(int)image->dim.x ; x++) {
          switch(image->type) {
          case SCHAR:
          case UCHAR:
              val = (int) imgarrayu8[z][y][x];
              break;
          case SSHORT:
          case USHORT:
              val = (int) imgarray[z][y][x];
              break;
          default:
              VT_Error("image type unknown or not supported for this program",program);
              return( 0 );
          }
          if(val!=0) break;
        }
        if(val!=0) break;
      }
      for (Y=0 ; Y<y ; Y++)
      for (x=0 ; x<(int)image->dim.x ; x++) {
        resarray[z][Y][x]=(unsigned char) 255;
      }
    }

    if (_verbose_)
      fprintf(stdout, "Grand Y...\n");
    /* Grand Y sur chaque z */
    for (z=0 ; z<(int)image->dim.z ; z++)
    {
      val=0;
      for (y=image->dim.y-1 ; y>=0 ; y--) {
        for (x=0 ; x<(int)image->dim.x ; x++) {
          switch(image->type) {
          case SCHAR:
          case UCHAR:
              val = (int) imgarrayu8[z][y][x];
              break;
          case SSHORT:
          case USHORT:
              val = (int) imgarray[z][y][x];
              break;
          default:
              VT_Error("image type unknown or not supported for this program",program);
              return( 0 );
          }
          if(val!=0) break;
        }
        if(val!=0) break;
      }
      for (Y=y+1 ; Y<(int)image->dim.y ; Y++)
      for (x=0 ; x<(int)image->dim.x ; x++) {
        resarray[z][Y][x]=(unsigned char) 255;
      }
    }

    if (_verbose_)
      fprintf(stdout, "Petit X...\n");
    /* Petit X sur chaque z */
    for (z=0 ; z<(int)image->dim.z ; z++)
    {
      val=0;
      for (x=0 ; x<(int)image->dim.x ; x++) {
        for (y=0 ; y<(int)image->dim.y ; y++) {
          switch(image->type) {
          case SCHAR:
          case UCHAR:
              val = (int) imgarrayu8[z][y][x];
              break;
          case SSHORT:
          case USHORT:
              val = (int) imgarray[z][y][x];
              break;
          default:
              VT_Error("image type unknown or not supported for this program",program);
              return( 0 );
          }
          if(val!=0) break;
        }
        if(val!=0) break;
      }
      for (X=0 ; X<x ; X++)
      for (y=0 ; y<(int)image->dim.y ; y++) {
        resarray[z][y][X]=(unsigned char) 255;
      }
    }

    if (_verbose_)
      fprintf(stdout, "Grand X...\n");
    /* Grand X sur chaque z */
    for (z=0 ; z<(int)image->dim.z ; z++)
    {
      val=0;
      for (x=image->dim.x-1 ; x>=0 ; x--) {
        for (y=0 ; y<(int)image->dim.y ; y++) {
          switch(image->type) {
          case SCHAR:
          case UCHAR:
              val = (int) imgarrayu8[z][y][x];
              break;
          case SSHORT:
          case USHORT:
              val = (int) imgarray[z][y][x];
              break;
          default:
              VT_Error("image type unknown or not supported for this program",program);
              return( 0 );
          }
          if(val!=0) break;
        }
        if(val!=0) break;
      }
      for (X=x+1 ; X<(int)image->dim.x ; X++)
      for (y=0 ; y<(int)image->dim.y ; y++) {
        resarray[z][y][X]=(unsigned char) 255;
      }
    }
  }
  else
  {
      VT_Error("this functionality has not been implemented yet\n",program);
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_FreeImage( &imres );
      return( 0 );
  }




  VT_WriteInrimage( &(imres) );

  /*--- liberations memoires ---*/
  VT_FreeImage( image );
  VT_Free( (void**)&image );
  VT_FreeImage( &imres );

  return( 0 ); /* returns the number of deleted labels */
}








static void MT_Parse( int argc,
                      char *argv[],
                      local_par *par )
{
  int i, nb;
  char text[STRINGLENGTH];

  if ( VT_CopyName( program, argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);
  if ( argc == 1 ) MT_ErrorParse("\n", 0 );

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
      else if ( strcmp ( argv[i], "-oblic" ) == 0 ) {
        par->oblic = 1;
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
      else
        MT_ErrorParse("too much file names when parsing\n", 0 );
    }
    i += 1;
  }

  /*--- s'il n'y a pas assez de noms ... ---*/
  if (nb == 0) {
    strcpy( par->names.in,  "<" );  /* standart input */
    strcpy( par->names.out, ">" );  /* standart output */
  }
  if (nb == 1)
    strcpy( par->names.out, ">" );  /* standart output */


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
  par->oblic = 0;
}



/* Three points are a counter-clockwise turn if ccw > 0, clockwise if
 * ccw < 0, and collinear if ccw = 0 because ccw is a determinant that
 * gives the signed area of the triangle formed by p1, p2 and p3.
 *
static double ccw(point_t* p1, point_t* p2, point_t* p3)
{
  return (p2->x - p1->x)*(p3->y - p1->y) - (p2->y - p1->y)*(p3->x - p1->x);
}

 * Returns a list of points on the convex hull in counter-clockwise order.
 * Note: the last point in the returned list is the same as the first one.
 *
static void convex_hull(point_t* points, ssize_t npoints, point_ptr_t** out_hull, ssize_t* out_hullsize)
{
  point_ptr_t* hull;
  ssize_t i, t, k = 0;

  hull = *out_hull;

  * lower hull *
  for (i = 0; i < npoints; ++i) {
    while (k >= 2 && ccw(hull[k-2], hull[k-1], &points[i]) <= 0) --k;
    hull[k++] = &points[i];
  }

  * upper hull *
  for (i = npoints-2, t = k+1; i >= 0; --i) {
    while (k >= t && ccw(hull[k-2], hull[k-1], &points[i]) <= 0) --k;
    hull[k++] = &points[i];
  }

  *out_hull = hull;
  *out_hullsize = k;
}
*/
