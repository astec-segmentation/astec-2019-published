/*************************************************************************
 * reechSlices.c -
 *
 * $Id: reechSlices.c,v 1.4 2002/09/05 17:15:06 greg Exp $
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
\t [-cspline | -trilin]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n\
 $Revision: 1.4 $ $Date: 2002/09/05 17:15:06 $ $Author: greg $\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image, imres, imgrd, imlap;
  char name[DOUBLESTRINGLENGTH];

  int i;
  int nbPts = 0;
  int firstPoint[3], theDim[3];
  int* theIntPts = (int*)NULL;

  typePoint* thePts = (typePoint* )NULL;

  typeCSplineCoefficients *theCSplineCoeff = NULL;
  typeCSplinesInSlice auxInSlice;



  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);





  fprintf( stderr, " construction de la liste de points \n" );





  /* on cherche une chaine 
   */

  if ( VT_ExtractFirstPoint( image, firstPoint, 
			     par.lowThreshold, par.highThreshold ) != 1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );    
    VT_ErrorParse("unable to extract first point\n", 0);
  }

  theIntPts = VT_ExtractLinkedCurve( image, firstPoint, &nbPts, 
			     par.lowThreshold, par.highThreshold );

  if ( theIntPts == (int*)NULL  ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );    
    VT_ErrorParse("unable to extract linked curve\n", 0 );
  }
  
  VT_FreeImage( image );
  VT_Free( (void**)&image );







  if ( 0 ) {
    for ( i=0; i<nbPts; i++ ) {
      fprintf( stderr, "%4d/%d   (%3d %3d %3d)  \n", i, nbPts, 
	       theIntPts[3*i], theIntPts[3*i+1], theIntPts[3*i+2] );
    }
  }





  


  /* on remplit la structure des points 
   */
  thePts = (typePoint* )malloc( nbPts*sizeof( typePoint ) );
  for ( i=0; i<nbPts; i++ ) {
    thePts[i].pt[0] = theIntPts[3*i+0];
    thePts[i].pt[1] = theIntPts[3*i+1];
    thePts[i].pt[2] = theIntPts[3*i+2];
  }



  sprintf( name, "%s.theta.inr", par.genericAngleName );
  image = _VT_Inrimage( name );
  if ( image == (vt_image*)NULL ) {
    sprintf( name, "%s.theta.inr", par.genericAngleName );
    image = _VT_Inrimage( name );
    if ( image == (vt_image*)NULL ) {
      free( theIntPts );
      free( thePts );
      VT_ErrorParse("unable to read theta input image\n", 0);
    }
  }
  switch( image->type ) {
  default :
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    free( theIntPts );
    free( thePts );
    VT_ErrorParse("unable to deal with theta input image type \n", 0);
    return( 1 );
  case FLOAT :
    {
      float ***theBuf = (float***)image->array;
      for ( i=0; i<nbPts; i++ ) {
	thePts[i].theta = theBuf[ theIntPts[3*i+2] ][ theIntPts[3*i+1] ][ theIntPts[3*i+0] ];
      }
    }
    break;
  }
  
  VT_FreeImage( image );
  VT_Free( (void**)&image );



  sprintf( name, "%s.phi.inr", par.genericAngleName );
  image = _VT_Inrimage( name );
  if ( image == (vt_image*)NULL ) {
    free( theIntPts );
    free( thePts );
    VT_ErrorParse("unable to read phi input image\n", 0);
  }
  switch( image->type ) {
  default :
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    free( theIntPts );
    free( thePts );
    VT_ErrorParse("unable to deal with phi input image type \n", 0);
    return( 1 );
  case FLOAT :
    {
      float ***theBuf = (float***)image->array;
      for ( i=0; i<nbPts; i++ ) {
	thePts[i].phi = theBuf[ theIntPts[3*i+2] ][ theIntPts[3*i+1] ][ theIntPts[3*i+0] ];
      }
    }
    break;
  }
  
  VT_FreeImage( image );
  VT_Free( (void**)&image );

  


  
  free( theIntPts );
  theIntPts = (int*)NULL;





  

  image =_VT_Inrimage( par.intensityImage );
  if ( image == (vt_image*)NULL ) {
    free( thePts );
    VT_ErrorParse("unable to read intensity input image\n", 0);
  }

  


  /*--- initialisation de l'image resultat ---*/
  sprintf( name, "%s.inr", par.names.out );
  VT_Image( &imres );
  VT_InitImage( &imres, name, par.dim, par.dim, nbPts, image->type );
  imres.siz.x = imres.siz.y = imres.siz.z = par.voxelSize;
  
  
  if ( VT_AllocImage( &imres ) != 1 ) {
    free( thePts );
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }


  if ( par.typeReech == _CSPLINE_ ) {
    theDim[0] = image->dim.x;
    theDim[1] = image->dim.y;
    theDim[2] = image->dim.z;
    theCSplineCoeff = ComputeCSplineCoefficients( image->buf, image->type, theDim );
    if ( theCSplineCoeff == NULL ) {
      free( thePts );
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_FreeImage( &imres );
      VT_ErrorParse("unable to allocate spline coeff\n", 0 );
    }

    if ( VT_AllocCSplinesInSlice( &auxInSlice, &imres ) != 1 ) {
      free( thePts );
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_FreeImage( &imres );
      FreeTypeCSplineCoefficients( &theCSplineCoeff );
      VT_ErrorParse("unable to allocate aux\n", 0 );
    }

    sprintf( name, "%s.grad.inr", par.names.out );
    VT_Image( &imgrd );
    VT_InitImage( &imgrd, name, par.dim, par.dim, nbPts, FLOAT );
    imgrd.siz.x = imgrd.siz.y = imgrd.siz.z = par.voxelSize;

    sprintf( name, "%s.lapl.inr", par.names.out );
    VT_Image( &imlap );
    VT_InitImage( &imlap, name, par.dim, par.dim, nbPts, FLOAT );
    imlap.siz.x = imlap.siz.y = imlap.siz.z = par.voxelSize;

    if ( VT_AllocImage( &imgrd ) != 1 ) {
      free( thePts );
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_FreeImage( &imres );
      FreeTypeCSplineCoefficients( &theCSplineCoeff );
      VT_FreeCSplinesInSlice( &auxInSlice );
      VT_ErrorParse("unable to allocate output image\n", 0);
    }

    if ( VT_AllocImage( &imlap ) != 1 ) {
      free( thePts );
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_FreeImage( &imres );
      FreeTypeCSplineCoefficients( &theCSplineCoeff );
      VT_FreeCSplinesInSlice( &auxInSlice );
      VT_FreeImage( &imgrd );
      VT_ErrorParse("unable to allocate output image\n", 0);
    }




  }





  /**************************************************
  {
    int x, y, z;
    int nb = 0;
    u8 ***theBuf = (u8***)image->array;
    float ***thetaBuf = (float***)imTheta->array;
    float ***phiBuf = (float***)imPhi->array;
    double pt[3];

    for ( z=0; z<image->dim.z; z++ )
    for ( y=0; y<image->dim.y; y++ )
    for ( x=0; x<image->dim.x; x++ ) {

      if ( theBuf[z][y][x] == 0 ) continue;

      fprintf( stderr, "%4d/%d   (%3d %3d %3d)  \n", nb, nbPts, x, y, z );

      pt[0] = x;
      pt[1] = y;
      pt[2] = z;
      VT_ReechTriLinSlice( &imres, nb, imInten, pt, (double)thetaBuf[z][y][x],
		     (double)phiBuf[z][y][x],  
		     (double*)NULL, (double*)NULL, (double*)NULL );
      nb++;
    }

  }
  **************************************************/


  for ( i=0; i<nbPts; i++ ) {
    if ( 0 ) {
      fprintf( stderr, "%4d/%d   (%f %f %f)  \n", i, nbPts, 
	       thePts[i].pt[0], thePts[i].pt[1], thePts[i].pt[2]  );
    }
    /*
    VT_ReechTriLinSlice( &imres, i, image, thePts[i].pt, thePts[i].theta,
		   thePts[i].phi, 
		   (double*)NULL, (double*)NULL, (double*)NULL );
    */

    switch ( par.typeReech ) {
    default :
    case _TRILIN_ :

      if ( i == 0 ) {
	VT_ReechTriLinSlice( &imres, i, image, thePts[i].pt, thePts[i].theta,
			     thePts[i].phi, thePts[i].frame, (double*)NULL,
			     thePts[i].mat );
      } else {
	VT_ReechTriLinSlice( &imres, i, image, thePts[i].pt, thePts[i].theta,
			     thePts[i].phi, thePts[i].frame, thePts[i-1].frame, 
			     thePts[i].mat );
      }
      break;

    case _CSPLINE_ :
      if ( i == 0 ) {
	VT_ReechCSplineSlice( &imres, &imgrd, &imlap,
			      i, theCSplineCoeff, &auxInSlice, image, 
			      thePts[i].pt, thePts[i].theta,
			      thePts[i].phi, thePts[i].frame, (double*)NULL,
			      thePts[i].mat );
      } else {
	VT_ReechCSplineSlice( &imres, &imgrd, &imlap,
			      i, theCSplineCoeff, &auxInSlice, image, 
			      thePts[i].pt, thePts[i].theta,
			      thePts[i].phi, thePts[i].frame, thePts[i-1].frame, 
			      thePts[i].mat );
      }
      break;

    }
    if ( 0 ) {
      fprintf( stderr, "   vecteur 1   ( %f %f %f )\n", 
	       thePts[i].frame[0], thePts[i].frame[3], thePts[i].frame[6] );
      fprintf( stderr, "   vecteur 2   ( %f %f %f )\n", 
	       thePts[i].frame[1], thePts[i].frame[4], thePts[i].frame[7] );
      fprintf( stderr, "   vecteur dir ( %f %f %f )\n", 
	       thePts[i].frame[2], thePts[i].frame[5], thePts[i].frame[8] );
    }

    if ( 0 ) {
      fprintf( stderr, "%d/%d (%f %f) -> (%f %f %f)\n",
	       i, nbPts, (imres.dim.x-1)/2.0, (imres.dim.y-1)/2.0, 
	       thePts[i].mat[ 0] *(imres.dim.x-1)/2.0 + thePts[i].mat[ 1] *(imres.dim.y-1)/2.0 + thePts[i].mat[ 3],
	       thePts[i].mat[ 4] *(imres.dim.x-1)/2.0 + thePts[i].mat[ 5] *(imres.dim.y-1)/2.0 + thePts[i].mat[ 7],
	       thePts[i].mat[ 8] *(imres.dim.x-1)/2.0 + thePts[i].mat[ 9] *(imres.dim.y-1)/2.0 + thePts[i].mat[11] );
    }

  }

















  {
    FILE *f, *fopen();
    int i, j, fd;
    float *vx, *vy, *vz;
    int x, y;
    
    vx = (float*)malloc( 3 * imres.dim.x * imres.dim.y * sizeof( float ) );
    vy = vz = vx;
    vy += imres.dim.x * imres.dim.y;
    vz += 2*imres.dim.x * imres.dim.y;

    i = strlen( par.names.out )-1;
    for ( ; i >= 0 && par.names.out[i] != '/' ; i-- ) 
      ;

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
    VT_3DDrawImage( image, fd, f );

    for ( i=0; i<nbPts; i++ ) {
      
      /* on ajoute +1 car les indices commencent a 1 en matlab
       */
      for ( j=0, y=0; y<(int)imres.dim.y; y++ )
  for ( x=0; x<(int)imres.dim.x; x++, j++ ) {
	  vx[j] = thePts[i].mat[ 0] * x + thePts[i].mat[ 1] * y + thePts[i].mat[ 3] + 1;
	  vy[j] = thePts[i].mat[ 4] * x + thePts[i].mat[ 5] * y + thePts[i].mat[ 7] + 1;
	  vz[j] = thePts[i].mat[ 8] * x + thePts[i].mat[ 9] * y + thePts[i].mat[11] + 1;
	}
      
      
      if ( write( fd, vx, 3 * imres.dim.x * imres.dim.y * sizeof( float ) ) == -1 ) {
	fprintf( stderr, "%s: error when writing\n", program );
      }
      fprintf( f, "VX%d = fread( fid, [%lu,%lu], 'float%lu' );\n", i,
	       imres.dim.x, imres.dim.y, 8*sizeof( float ) );
      fprintf( f, "VY%d = fread( fid, [%lu,%lu], 'float%lu' );\n",  i,
	       imres.dim.x, imres.dim.y, 8*sizeof( float ) );
      fprintf( f, "VZ%d = fread( fid, [%lu,%lu], 'float%lu' );\n",  i,
	       imres.dim.x ,imres.dim.y, 8*sizeof( float ) );
      fprintf( f, "%% H%d = slice( double(IMAGE), VX%d, VY%d, VZ%d );\n", i, i, i, i );
      fprintf( f, "%% set( H%d, 'FaceColor', 'interp', 'EdgeColor', 'none' );\n", i );
      fprintf( f, "HC%d = contourslice( double(IMAGE), VX%d, VY%d, VZ%d, [%f %f] );\n", 
	       i, i, i, i, par.surfThreshold, par.surfThreshold );
      fprintf( f, "set( HC%d, 'EdgeColor', 'blue', 'LineWidth', 3 );\n", i );    
    }

    fprintf( f, "\n\n" );
    


    /* on ajoute +1 car les indices commencent a 1 en matlab
     */
    {
      float *t = NULL;
     
      t = (float*)malloc( nbPts*sizeof(float) );
      
      for ( i=0; i<nbPts; i++ ) t[i] = thePts[i].pt[0] + 1;
      if ( write( fd, t, nbPts*sizeof(float) ) == -1 ) {
	fprintf( stderr, "%s: error when writing\n", program );
      }
      fprintf( f, "XCENTER = fread( fid, %d, 'float%lu' );\n", nbPts, 8*sizeof( float ) );
      
      for ( i=0; i<nbPts; i++ ) t[i] = thePts[i].pt[1] + 1;
      if ( write( fd, t, nbPts*sizeof(float) ) == -1 ) {
	fprintf( stderr, "%s: error when writing\n", program );
      }
      fprintf( f, "YCENTER = fread( fid, %d, 'float%lu' );\n", nbPts, 8*sizeof( float ) );
      
      for ( i=0; i<nbPts; i++ ) t[i] = thePts[i].pt[2] + 1;
      if ( write( fd, t, nbPts*sizeof(float) ) == -1 ) {
	fprintf( stderr, "%s: error when writing\n", program );
      }
      fprintf( f, "ZCENTER = fread( fid, %d, 'float%lu' );\n", nbPts, 8*sizeof( float ) );
      
      fprintf( f, "plot3( XCENTER, YCENTER, ZCENTER, 'r-', 'LineWidth', 3 );\n" );
      
      free(t);
    }

    fprintf( f, "hold off;\n" );
    fprintf( f, "fclose( fid );\n" );
    close( fd );
    fclose( f );
  }

















  if ( par.typeReech == _CSPLINE_ ) {
    FreeTypeCSplineCoefficients( &theCSplineCoeff );
    VT_FreeCSplinesInSlice( &auxInSlice );
    VT_WriteInrimage( &imgrd );
    VT_FreeImage( &imgrd );
    VT_WriteInrimage( &imlap );
    VT_FreeImage( &imlap );
  }


  
  VT_FreeImage( image );
  VT_Free( (void**)&image );
  free( thePts );



  /*--- ecriture de l'image resultat ---*/
  if ( VT_WriteInrimage( &imres ) == -1 ) {
    VT_FreeImage( &imres );
    VT_ErrorParse("unable to write output image\n", 0);
  }
  
  /*--- liberations memoires ---*/
  VT_FreeImage( &imres );
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
