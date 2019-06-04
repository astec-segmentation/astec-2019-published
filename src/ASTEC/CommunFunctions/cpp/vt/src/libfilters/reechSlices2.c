/*************************************************************************
 * reechSlices.c -
 *
 * $Id: reechSlices2.c,v 1.3 2002/12/09 16:35:25 greg Exp $
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

#include <vt_common.h>
#include <vt_tube3D.h>
#include <vt_link.h>
#include <vt_tube2Dmatlab.h>
#include <transfo.h>

typedef struct {

  double pt_mm[3];
  double pt_voxel[3];

  /* repere
     (#2 #5 #8) est le vecteur "directeur" (ie dans la direction du vaisseau)
     (#0 #3 #6) et (#1 #4 #7) sont les deux vecteurs orthogonaux
  */
  double frame_mm[9];

  double transfo_slice_mm[16];

} typePoint;
    

typedef enum {
  _TRILIN_,
  _CSPLINE_
} enumTypeReech;


typedef enum {
  _DEFAULT_,
  _PRISME_
} enumTypePropagation;



typedef struct local_par {
  vt_names names;

  char genericAngleName[STRINGLENGTH]; 
  char intensityImage[STRINGLENGTH]; 

  double voxelSize;
  int dim;

  int lowThreshold;
  int highThreshold;

  enumTypeReech typeReech;
  enumTypePropagation typePropagation;
  double tolerance_vectorial_product;

  float surfThreshold;

} local_par;




/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[image-in] [image-out] [-angles %s]\n\
\t [-image|-i %s] [-lt|-sb %d] [-ht|-sh %d] [-cc %d] [-dim %d] [-vs %lf]\n\
\t [-cspline | -trilin] [-tol %lf] [-propag prisme | default] \n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n\
 $Revision: 1.3 $ $Date: 2002/12/09 16:35:25 $ $Author: greg $\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image, imres;
  char name[DOUBLESTRINGLENGTH];

  int i;
  int nbPts = 0;
  int firstPoint[3], theDim[3];
  int* theIntPts = (int*)NULL;

  typePoint* thePts = (typePoint* )NULL;
  double *theta = NULL, *phi = NULL;
  double curve_voxel_size[3];
  double angle_voxel_size[3];
  double norm, vdir[3];
  double mat[16];


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
  curve_voxel_size[0] = image->siz.x;
  curve_voxel_size[1] = image->siz.y;
  curve_voxel_size[2] = image->siz.z;
  theIntPts = VT_ExtractLinkedCurve( image, firstPoint, &nbPts, 
			     par.lowThreshold, par.highThreshold );

  if ( theIntPts == (int*)NULL  || nbPts <= 0 ) {
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
     les points sont en mm
   */
  thePts = (typePoint* )malloc( nbPts*sizeof( typePoint ) );
  for ( i=0; i<nbPts; i++ ) {
    thePts[i].pt_mm[0] = theIntPts[3*i+0] * curve_voxel_size[0];
    thePts[i].pt_mm[1] = theIntPts[3*i+1] * curve_voxel_size[1];
    thePts[i].pt_mm[2] = theIntPts[3*i+2] * curve_voxel_size[2];
  }
  theta = (double*)malloc( 2 * nbPts*sizeof( double ) );
  phi = theta;
  phi += nbPts;



  sprintf( name, "%s.theta.inr", par.genericAngleName );
  image = _VT_Inrimage( name );
  if ( image == (vt_image*)NULL ) {
    sprintf( name, "%s.theta.inr.gz", par.genericAngleName );
    image = _VT_Inrimage( name );
    if ( image == (vt_image*)NULL ) {
      free( theIntPts );
      free( thePts );
      VT_ErrorParse("unable to read theta input image\n", 0);
    }
  }
  angle_voxel_size[0] = image->siz.x;
  angle_voxel_size[1] = image->siz.y;
  angle_voxel_size[2] = image->siz.z;
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
	theta[i] = theBuf[ theIntPts[3*i+2] ][ theIntPts[3*i+1] ][ theIntPts[3*i+0] ];
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
	phi[i] = theBuf[ theIntPts[3*i+2] ][ theIntPts[3*i+1] ][ theIntPts[3*i+0] ];
      }
    }
    break;
  }
  VT_FreeImage( image );
  VT_Free( (void**)&image );

  


  
  free( theIntPts );
  theIntPts = (int*)NULL;
  





  fprintf( stderr, " calcul des reperes (en mm)\n" );

  /* calcul des vecteurs directeurs
   */
  for ( i=0; i<nbPts; i++ ) {
    SphericalAnglesToUnitVector( theta[i], phi[i], vdir );
    vdir[0] *= angle_voxel_size[0];
    vdir[1] *= angle_voxel_size[1];
    vdir[2] *= angle_voxel_size[2];
    norm = sqrt( vdir[0]*vdir[0] + vdir[1]*vdir[1] + vdir[2]*vdir[2] );
    vdir[0] /= norm;
    vdir[1] /= norm;
    vdir[2] /= norm;
    thePts[i].frame_mm[0] = thePts[i].frame_mm[3] = thePts[i].frame_mm[6] = 0.0;
    thePts[i].frame_mm[1] = thePts[i].frame_mm[4] = thePts[i].frame_mm[7] = 0.0;
    thePts[i].frame_mm[2] = vdir[0];
    thePts[i].frame_mm[5] = vdir[1];
    thePts[i].frame_mm[8] = vdir[2];
  }
  /* calcul des reperes
     on a:
     - les points en mm            -> thePts[i].pt_mm[]
     - le vecteur directeur en mm  -> 
     - la derniere colonne de la matrice de transformation -> thePts[i].frame_mm[]
   */
  VT_ComputeRealFrame( thePts[0].frame_mm, NULL );


  switch( par.typePropagation ) {
  default :
  case _DEFAULT_ :
    for ( i=1; i<nbPts; i++ ) {
      VT_ComputeRealFrame( thePts[i].frame_mm, 
			   thePts[i-1].frame_mm );
    }
    break;
  case _PRISME_ :
    for ( i=1; i<nbPts; i++ ) {
      VT_ComputeRealFrameWithPivot( thePts[i].frame_mm, 
				    thePts[i-1].frame_mm,
				    par.tolerance_vectorial_product );
    }
    break;
  }



  /* calcul des transformations
     pour reechantillonner chacune des coupes
     du plan (en pixels) vers l'espace reel
     le centre du plan doit etre le point de la courbe

     M_3D_reel(id en mm) = MAT * M_2D_pixel
   */
  for ( i=0; i<nbPts; i++ ) {
    thePts[i].transfo_slice_mm[ 0] = par.voxelSize * thePts[i].frame_mm[0];
    thePts[i].transfo_slice_mm[ 1] = par.voxelSize * thePts[i].frame_mm[1];
    thePts[i].transfo_slice_mm[ 2] = par.voxelSize * thePts[i].frame_mm[2];

    thePts[i].transfo_slice_mm[ 4] = par.voxelSize * thePts[i].frame_mm[3];
    thePts[i].transfo_slice_mm[ 5] = par.voxelSize * thePts[i].frame_mm[4];
    thePts[i].transfo_slice_mm[ 6] = par.voxelSize * thePts[i].frame_mm[5];

    thePts[i].transfo_slice_mm[ 8] = par.voxelSize * thePts[i].frame_mm[6];
    thePts[i].transfo_slice_mm[ 9] = par.voxelSize * thePts[i].frame_mm[7];
    thePts[i].transfo_slice_mm[10] = par.voxelSize * thePts[i].frame_mm[8];

    thePts[i].transfo_slice_mm[12] = 0.0;
    thePts[i].transfo_slice_mm[13] = 0.0;
    thePts[i].transfo_slice_mm[14] = 0.0;
    thePts[i].transfo_slice_mm[15] = 1.0;

    thePts[i].transfo_slice_mm[ 3] = thePts[i].pt_mm[0] - 
      ( thePts[i].transfo_slice_mm[ 0] * (par.dim-1)/2.0 
      + thePts[i].transfo_slice_mm[ 1] * (par.dim-1)/2.0 );
    thePts[i].transfo_slice_mm[ 7] = thePts[i].pt_mm[1] - 
      ( thePts[i].transfo_slice_mm[ 4] * (par.dim-1)/2.0 
      + thePts[i].transfo_slice_mm[ 5] * (par.dim-1)/2.0 );
    thePts[i].transfo_slice_mm[11] = thePts[i].pt_mm[2] - 
      ( thePts[i].transfo_slice_mm[ 8] * (par.dim-1)/2.0 
      + thePts[i].transfo_slice_mm[ 9] * (par.dim-1)/2.0 );

  }

  {
    FILE *ftrsf;
    int j;
      
    sprintf( name, "%s.transfos", par.names.out );
    ftrsf = fopen( name, "w" );
    
    fprintf( ftrsf, "#" );
    for ( i=0; i<argc; i++ ) fprintf( ftrsf, " %s", argv[i] );
    fprintf( ftrsf, "\n" );
    fprintf( ftrsf, "\n" );

    for ( i=0; i<nbPts; i++ ) {
      for (j=0; j<16; j++) 
	fprintf( ftrsf, " %g", thePts[i].transfo_slice_mm[j] );
      fprintf( ftrsf, "\n" );
    }

    fclose( ftrsf );
  }


  fprintf( stderr, " reechantillonnage\n" );



  image =_VT_Inrimage( par.intensityImage );
  if ( image == (vt_image*)NULL ) {
    free( thePts );
    VT_ErrorParse("unable to read intensity input image\n", 0);
  }

  


  /*--- initialisation de l'image resultat ---*/
  sprintf( name, "%s.inr", par.names.out );
  VT_Image( &imres );
  VT_InitImage( &imres, name, par.dim, par.dim, nbPts, image->type );
  imres.siz.x = imres.siz.y = par.voxelSize;
  imres.siz.z = par.voxelSize;
  
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

  }





  for ( i=0; i<nbPts; i++ ) {

    fprintf( stderr, " processing slice #%3d/%d\r",
	     i+1, nbPts);

    mat[ 0] = thePts[i].transfo_slice_mm[ 0] / image->siz.x;
    mat[ 1] = thePts[i].transfo_slice_mm[ 1] / image->siz.x;
    mat[ 2] = thePts[i].transfo_slice_mm[ 2] / image->siz.x;
    mat[ 3] = thePts[i].transfo_slice_mm[ 3] / image->siz.x;
    mat[ 4] = thePts[i].transfo_slice_mm[ 4] / image->siz.y;
    mat[ 5] = thePts[i].transfo_slice_mm[ 5] / image->siz.y;
    mat[ 6] = thePts[i].transfo_slice_mm[ 6] / image->siz.y;
    mat[ 7] = thePts[i].transfo_slice_mm[ 7] / image->siz.y;
    mat[ 8] = thePts[i].transfo_slice_mm[ 8] / image->siz.z;
    mat[ 9] = thePts[i].transfo_slice_mm[ 9] / image->siz.z;
    mat[10] = thePts[i].transfo_slice_mm[10] / image->siz.z;
    mat[11] = thePts[i].transfo_slice_mm[11] / image->siz.z;

    /*
    VT_ReechTriLinSlice( &imres, i, image, thePts[i].pt, thePts[i].theta,
		   thePts[i].phi, 
		   (double*)NULL, (double*)NULL, (double*)NULL );
    */

    switch ( par.typeReech ) {
    default :
    case _TRILIN_ :

      VT_ReechTriLinSlice2( &imres, i, image, mat );
      break;

    case _CSPLINE_ :
      VT_ReechCSplineSlice2( &imres, i, theCSplineCoeff, &auxInSlice, mat );
      break;

    }
    
  }

















  if ( par.typeReech == _CSPLINE_ ) {
    FreeTypeCSplineCoefficients( &theCSplineCoeff );
    VT_FreeCSplinesInSlice( &auxInSlice );
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


      else if ( strcmp ( argv[i], "-tol" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -tol...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->tolerance_vectorial_product) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -tol...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-propagation" ) == 0 || 
		strcmp ( argv[i], "-propag" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -propagation...\n", 0 );
	if ( strcmp ( argv[i], "default" ) == 0 ||
	     strcmp ( argv[i], "defaut" ) == 0 ) {
	  par->typePropagation = _DEFAULT_;
	}
	else if ( strcmp ( argv[i], "prisme" ) == 0 ) {
	  par->typePropagation = _PRISME_;
	}
	else {
	  sprintf(text,"unknown propagation option '%s'\n",argv[i]);
	  VT_ErrorParse(text, 0);
	}
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
  par->typePropagation = _DEFAULT_;
  par->tolerance_vectorial_product = 0.05;

  par->surfThreshold = 50;

}
