/*************************************************************************
 * reechSlices.c -
 *
 * $Id: reechSlices.c,v 1.1 2000/10/20 11:05:31 greg Exp $
 *
 * Copyright (c) INRIA 2000
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * 
 *
 * ADDITIONS, CHANGES
 *
 */

#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>

#include <vt_common.h>
#include <vt_tube3D.h>
#include <vt_link.h>
#include <vt_tube2Dmatlab.h>
#include <vt_histo.h>

typedef struct {
  double pt[3];
  double theta;
  double phi;
  double frame[9];
  double mat[16];
} typePoint;
    

typedef enum {
  _TRILIN_,
  _CSPLINE_
} enumTypeReech;



typedef struct local_par {
  vt_names names;

  char genericAngleName[STRINGLENGTH]; 
  char intensityImage[STRINGLENGTH]; 

  double voxelSize;
  int dim;

  int lowThreshold;
  int highThreshold;

  enumTypeReech typeReech;

  float surfThreshold;

} local_par;




/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[image-in] [image-out] [-angles %s]\n\
\t [-image|-i %s] [-lt|-sb %d] [-ht|-sh %d] [-cc %d] [-dim %d] [-vs %lf]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n\
 $Revision: 1.1 $ $Date: 2000/10/20 11:05:31 $ $Author: greg $\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image;
  char name[DOUBLESTRINGLENGTH];

  int c, i, j;
  int nbPts = 0;
  int firstPoint[3];
  int* theIntPts = (int*)NULL;
  vt_3m m;

  typePoint* thePts = (typePoint* )NULL;

  FILE *f, *fopen();
  int fd;
  float *t;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );

  if ( par.intensityImage[0] != '\0' ) {
    image =_VT_Inrimage( par.intensityImage );
    if ( image == (vt_image*)NULL ) {
      free( thePts );
      VT_ErrorParse("unable to read intensity input image\n", 0);
    }
  }

  i = strlen( par.names.out )-1;
  for ( ; i >= 0 && par.names.out[i] != '/' ; i-- ) 
    sprintf( name, "%s.raw", par.names.out );
  fd = creat( name, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
  sprintf( name, "%s.m", par.names.out );
  f = fopen( name, "w" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );
  fprintf( f, "%%\n" );
  fprintf( f, "%%" );
  for ( j=0; j<argc; j++ ) fprintf( f, " %s", argv[j] );
  fprintf( f, "\n" );
  fprintf( f, "%%\n" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );
  fprintf( f, "echo off\n" );
  fprintf( f, "fid = fopen('%s.raw', 'r' );\n", &(par.names.out[i+1]) );
  fprintf( f, "figure;\n" );
  fprintf( f, "hold on;\n" );

  if ( par.intensityImage[0] != '\0' ) {

    VT_3DDrawImage( image, fd, f );
    VT_FreeImage( image );
    VT_Free( (void**)&image );
  }



  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);


  (void)VT_3m( image, (vt_image*)NULL, &m );


  if ( par.lowThreshold < (int)(m.min+0.5) ) 
    par.lowThreshold = (int)(m.min+0.5);
  if ( par.lowThreshold < 1 )
    par.lowThreshold = 1;
  
   if ( par.highThreshold > (int)(m.min+0.5) ) 
    par.highThreshold = (int)(m.min+0.5);
  if ( par.highThreshold < 1 )
    par.highThreshold = 1;
 

  for ( c=par.lowThreshold; c<=par.highThreshold; c++ ) {
  
    fprintf( stderr, "processing %2d/%d\n", c, (int)(m.max+0.5) );

    if ( VT_ExtractFirstPoint( image, firstPoint, c, c ) != 1 ) continue;
    
    theIntPts = VT_ExtractLinkedCurve( image, firstPoint, &nbPts, c, c );

    if ( theIntPts == (int*)NULL  ) continue;
  
    t = (float*)malloc( nbPts*sizeof(float) );

    for ( i=0; i<nbPts; i++ ) t[i] = image->siz.x * theIntPts[3*i+0];
    if ( write( fd, t, nbPts*sizeof(float) ) == -1 ) {
      fprintf( stderr, "%s: error when writing\n", program );
    }
    fprintf( f, "XCENTER%d = fread( fid, %d, 'float%lu' );\n", c, nbPts, 8*sizeof( float ) );
    for ( i=0; i<nbPts; i++ ) t[i] = image->siz.y * theIntPts[3*i+1];
    if ( write( fd, t, nbPts*sizeof(float) ) == -1 ) {
      fprintf( stderr, "%s: error when writing\n", program );
    }
    fprintf( f, "YCENTER%d = fread( fid, %d, 'float%lu' );\n", c, nbPts, 8*sizeof( float ) );
    for ( i=0; i<nbPts; i++ ) t[i] = image->siz.z * theIntPts[3*i+2];
    if ( write( fd, t, nbPts*sizeof(float) ) == -1 ) {
      fprintf( stderr, "%s: error when writing\n", program );
    }
    fprintf( f, "ZCENTER%d = fread( fid, %d, 'float%lu' );\n", c, nbPts, 8*sizeof( float ) );
    /*
      y     yellow        .     point              -     solid
      m     magenta       o     circle             :     dotted
      c     cyan          x     x-mark             -.    dashdot 
      r     red           +     plus               --    dashed   
      g     green         *     star
      b     blue          s     square
      w     white         d     diamond
      k     black   
    */
    switch ( c%6 ) {
    case 1 :
      fprintf( f, "plot3( XCENTER%d, YCENTER%d, ZCENTER%d, 'r-', 'LineWidth', 3 );\n", c, c, c );
      break;
    case 2 :
      fprintf( f, "plot3( XCENTER%d, YCENTER%d, ZCENTER%d, 'b-', 'LineWidth', 3 );\n", c, c, c );
      break;
    case 3 :
      fprintf( f, "plot3( XCENTER%d, YCENTER%d, ZCENTER%d, 'g-', 'LineWidth', 3 );\n", c, c, c );
      break;
    case 4 :
      fprintf( f, "plot3( XCENTER%d, YCENTER%d, ZCENTER%d, 'y-', 'LineWidth', 3 );\n", c, c, c );
      break;
    case 5 :
      fprintf( f, "plot3( XCENTER%d, YCENTER%d, ZCENTER%d, 'c-', 'LineWidth', 3 );\n", c, c, c );
      break;
    case 0 :
      fprintf( f, "plot3( XCENTER%d, YCENTER%d, ZCENTER%d, 'm-', 'LineWidth', 3 );\n", c, c, c );
      break;
    }

    free(t);
    free( theIntPts );
    
  }

  fprintf( f, "view(3);\n" );
  fprintf( f, "grid;\n" );
  fprintf( f, "axis equal;\n" );
  fprintf( f, "axis vis3d;\n" );
  fprintf( f, "hold off;\n" );
  fprintf( f, "fclose( fid );\n" );
  close( fd );
  fclose( f );


  VT_FreeImage( image );
  VT_Free( (void**)&image );






  return( 1 );
}















static void VT_Parse( int argc, 
		      char *argv[], 
		      local_par *par )
{
  int i, nb, status;
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






      else if ( strcmp ( argv[i], "-vs" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -vs...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->voxelSize) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -vs...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-dim" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -dim...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->dim) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -dim...\n", 0 );
      }




      else if ( strcmp ( argv[i], "-sb" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -sb...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->lowThreshold) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -sb...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-lt" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -lt...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->lowThreshold) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -lt...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-ht" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -ht...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->highThreshold) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -ht...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-sh" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -sh...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->highThreshold) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -sh...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-surf" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -surf...\n", 0 );
	status = sscanf( argv[i],"%f",&(par->surfThreshold) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -surf...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-cc" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -cc...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->highThreshold) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -cc...\n", 0 );
	par->lowThreshold = par->highThreshold;
      }


      else if ( strcmp ( argv[i], "-trilin" ) == 0 ) {
	par->typeReech = _TRILIN_;
      }
      else if ( strcmp ( argv[i], "-cspline" ) == 0 ) {
	par->typeReech = _CSPLINE_;
      }



      /*---  ---*/
      else if ( strcmp ( argv[i], "-angles" ) == 0 || 
		strcmp ( argv[i], "-a" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -angles...\n", 0 );
	strncpy( par->genericAngleName, argv[i], STRINGLENGTH );  
      }
      else if ( strcmp ( argv[i], "-image" ) == 0 || 
		strcmp ( argv[i], "-i" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -image...\n", 0 );
	strncpy( par->intensityImage, argv[i], STRINGLENGTH );  
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
  if (nb == 1)
    strcpy( par->names.out, ">" );  /* standart output */
  
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

  par->genericAngleName[0] = '\0';
  par->intensityImage[0] = '\0';

  par->voxelSize = 0.2;
  par->dim = 33;

  par->lowThreshold  = 1;
  par->highThreshold = 255;

  par->typeReech = _TRILIN_;

  par->surfThreshold = 50;

}
