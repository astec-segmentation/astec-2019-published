/*************************************************************************
 * minimum.c -
 *
 * $Id: test-contours-3D.c,v 1.7 2006/05/16 09:33:34 greg Exp $
 *
 * Copyright (c) INRIA 1999
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
#include <vt_isocontours.h>

#include <vt_tube2Dmatlab.h>
#include <vt_contoursMatlab.h>
#include <vt_tubeutils.h>



typedef enum {
  _2D_,
  _3D_
} prismeTypeOutput;



typedef struct local_par {
  vt_names names;

  int min_slice;
  int max_slice;
  double threshold;

  char trsfs[STRINGLENGTH];

  int draw_frames;

  char cnts[STRINGLENGTH];
  double tolerance_vectorial_product;
  prismeTypeOutput prisme_output;
  
} local_par;




/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[image-in] [image-out]\n\
\t [-slice %d | [-min-slice %d & -max-slice %d]] [-th %lf]\n\
\t [-trsfs %s] [-matlab %s] [-draw-frames]\n\
\t [-cnts %s] [-tol %f] [-2D] [-3D]\n\
\t [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
 si 'image-in' est '-', on prendra stdin\n\
 si 'image-out' est absent, on prendra stdout\n\
 si les deux sont absents, on prendra stdin et stdout\n\
 -slice %d:     process a single slice\n\
 -min-slice %d: first slice to be processed\n\
 -max-slice %d: last slice to be processed\n\
 -th %lf:       threshold for the iso-contours\n\
 -trsfs %s:     file of 3D transformations such that\n\
                M_3D_real(id in mm) = MAT * M_2D_slice(id in pixel)\n\
\n\
 -matlab %s:    generated matlab files (.m and .raw)\n\
 -draw-frames:  draw frames (in matlab output, versatile)\n\
\n\
 -cnts %s: Prisme compatible output\n\
 -2D:      generates '%s.next.cnt', '%s.prev.cnt', and '%s.3D.cnt'\n\
     '%s.next.cnt' contains 2D coordinates of contours #i in a local frame\n\
     compatible with a local frame of contours #(i+1) in\n\
     '%s.prev.cnt'. '%s.3D.cnt' contains 3D real coordinates of contours\n\
-tol:      error (to determine if two successive cross-sections intersect\n\
     (for '-2D' mode in Prisme compatible output)\n\
-3D:       generates '%s.2D.cnt' and '%s.2D.trsf'. '%s.2D.cnt' contains 2D\n\
     coordinates of contours in their native frame, while '%s.2D.trsf'\n\
     contains the transformation towards the real frame. Ie for contour #i\n\
     M_3D_reel(id en mm) = MAT(i) * M_2D_slice(id en pixel)\n\
 -v : mode verbose\n\
 -D : mode debug\n\
\n\
 $Revision: 1.7 $ $Date: 2006/05/16 09:33:34 $ $Author: greg $\n";

static char program[STRINGLENGTH];




typedef struct {
  int z;
  int imat;
  double *mat;
  int icnt;
  typeContour2D *with_prev;
  typeContour2D *theCnt;
  typeContour2D *with_next;
  double sliceCenter[2];
  double sliceCorner[2];
} typeSingleContour;



void _DrawSliceBorders( typeSingleContour *s, FILE *f );
void _DrawSliceFrame( typeSingleContour *s, FILE *f );
int _ComputeJointFrames( typeSingleContour *s1,
			 typeSingleContour *s2,
			 double *c1,
			 double *x1,
			 double *y1,
			 double *c2,
			 double *x2,
			 double *y2,
			 double tolerance );
void _DrawJointFrames( typeSingleContour *s1,
		       typeSingleContour *s2,
		       FILE *f, double tol );
int _ReCompute2DContours( typeSingleContour *s1,
			  typeSingleContour *s2,
			  double tol );
void _PrintContourInCnt( typeContour2D *c, double *mat,
			 int n,
			 double z,
			 int print_z, FILE *f );
int InverseMat4x4( double *matrice, double *inv );





int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image, imres;

  typeStructure structure;
  typeSlice *slice;

  int n, ncontours = 0;
  int z;
  int j, i;
  int next;
  
  double identity[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  
  double **mat = NULL;
  int ntrsfs = 0;

  int nsinglecontours;
  typeSingleContour *theSingleContour = NULL;

  char matlab_type[20];

  int m, n_step = 4;
  int n_realcontours;
  double *xr, *yr, *zr;

  double *csurf;
  double *ci, cmax;
  int cn;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);
  
  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.names.inv == 1 )  VT_InverseImage( image );
  if ( par.names.swap == 1 ) VT_SwapImage( image );



  
  /*--- initialisation de l'image resultat ---*/
  VT_Image( &imres );
  VT_InitFromImage( &imres, image, par.names.out, UCHAR );

  if ( VT_AllocImage( &imres ) != 1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }

  VT_ComputeNumberOfContoursPerCell( image, &imres, par.threshold );

  /*--- ecriture de l'image resultat ---*/
  if ( par.names.out[0] != '\0' && par.names.out[0] != '>' ) {
    if ( VT_WriteInrimage( &imres ) == -1 ) {
      VT_FreeImage( image );
      VT_FreeImage( &imres );
      VT_Free( (void**)&image );
      VT_ErrorParse("unable to write output image\n", 0);
    }
  }

  VT_FreeImage( &imres );


  if ( par.min_slice < 0 ) par.min_slice = 0;
  if ( par.min_slice >= (int)image->dim.z ) par.min_slice = image->dim.z - 1;
  if ( par.max_slice < 0 ) par.max_slice = 0;
  if ( par.max_slice >= (int)image->dim.z ) par.max_slice = image->dim.z - 1;



  ncontours = 0;
  initStructure( &structure );

  for ( z = par.min_slice; z <= par.max_slice; z++ ) {

    slice = VT_Compute2DIsoContours( image, z, par.threshold );
    if ( slice == NULL ) {
      printf( "bof (z=%d)\n", z );
    }
    
    if ( 0 ) printf( "found %d contours in slice %d\n", slice->n, z  );

    ncontours += slice->n;

    if ( addSliceToStructure( slice, &structure ) != 1 ) {
      VT_ErrorParse("error when adding slice to structure\n", 0);
    }

  }
  printf( "found %d contours\n", ncontours );
  


  /* les derniers contours peuvent etre vides, 
     ils ont ete initialises mais pas remplis
  */
  /*
  for (n=0, j=0; j<listOfContours.n; j++ ) {
    if ( listOfContours.theCts[j].n <= 0 ) {
      printf( "contour #%3d is empty\n", j );
      n++;
    }
  }
  if ( n > 0 )
    printf( "found %d empty contours\n", n );
  */


  if ( 0 ) {
    printStructure( &structure, 2 );
  }





  /* transformations
     pour reechantillonner chacune des coupes
     du plan (en pixels) vers l'espace reel
     le centre du plan doit etre le point de la courbe

     M_3D_reel(id en mm) = MAT * M_2D_slice(id en pixel)
   */
  if ( par.trsfs[0] != '\0' )
    mat = _ReadTransformationsFile( par.trsfs, &ntrsfs );




  /*
  for (j=0; j<listOfContours.n; j++ ) {
    printf( "%3d #=%d z=%f\n",j, listOfContours.theCts[j].n, 
  }
  */
  


  if ( par.max_slice >= par.min_slice ) {

    nsinglecontours = par.max_slice - par.min_slice + 1;
    theSingleContour = (typeSingleContour *)malloc( nsinglecontours
						    * sizeof(typeSingleContour) ); 
    if ( theSingleContour == NULL ) {
      VT_ErrorParse("unable to allocate single contours\n", 0);
    }
    
    for ( n = 0; n < nsinglecontours; n++ ) {
      theSingleContour[n].z = par.min_slice + n;
      theSingleContour[n].imat = -1;
      theSingleContour[n].mat = NULL;
      theSingleContour[n].icnt = -1;
      theSingleContour[n].theCnt = NULL;
      theSingleContour[n].with_prev = NULL;
      theSingleContour[n].with_next = NULL;
    }
    
    
    /* get the contours
       - identify the matrix
       - identify one contour corresponding to the slice
         - it should be closed
	 - it should contain the centre of the slice
	 - it is the smallest
     */
    for ( n = 0; n < nsinglecontours; n++ ) {

      if ( mat != NULL ) {
	theSingleContour[n].imat = theSingleContour[n].z;
	theSingleContour[n].mat  = mat[ theSingleContour[n].z ];
      }

      for (j=0; j<structure.n; j++ ) {

	if ( structure.theSlices[j]->z < theSingleContour[n].z - 0.1 ||
	     structure.theSlices[j]->z > theSingleContour[n].z + 0.1 )
	  continue;

	if ( structure.theSlices[j]->n == 0 ) continue;
	
	for ( i=0; i<structure.theSlices[j]->n; i ++ ) {

	  if ( structure.theSlices[j]->theContours[i]->n == 0 )
	    continue;
	  if ( structure.theSlices[j]->theContours[i]->topology == _OPEN_ )
	    continue;
	  if ( _PointInContour( structure.theSlices[j]->theContours[i],
				(image->dim.x-1)/2.0, (image->dim.y-1)/2.0 ) != 1 )
	    continue;
	
	  if ( theSingleContour[n].theCnt == NULL ) {
	    theSingleContour[n].icnt = j;
	    theSingleContour[n].theCnt = structure.theSlices[j]->theContours[i];
	  }
	  else {
	    fprintf( stdout, "slice %d has already one valid contour (i=%d)\n",
		     theSingleContour[n].z, theSingleContour[n].icnt );
	    if ( theSingleContour[n].theCnt->n 
		 > structure.theSlices[j]->theContours[i]->n ) {
	      fprintf( stdout, "   switch for smaller contour #%d -> #%d\n", 
		       theSingleContour[n].icnt, j );
	      theSingleContour[n].icnt = j;
	      theSingleContour[n].theCnt = structure.theSlices[j]->theContours[i];
	    }
	    else
	      fprintf( stdout, "   not accepting additional contour #%d\n", j );
	  }

	}
      }

    }





    /* print in 3D for prisme
     */

    if ( par.cnts[0] != '\0' && par.prisme_output == _3D_ ) {

      {
        char name[DOUBLESTRINGLENGTH];
	int s = 0;

	FILE *f2d, *ftrsf, *fopen();
	sprintf( name, "%s.2D.cnt", par.cnts );
	f2d = fopen( name, "w" );
	sprintf( name, "%s.2D.trsf", par.cnts );
	ftrsf = fopen( name, "w" );

	for ( n = 0; n < nsinglecontours; n++ ) {
	  if ( theSingleContour[n].theCnt == NULL ) continue;
	  s++;
	}
	fprintf( f2d, "S %d\n", s );

	for ( n = 0; n < nsinglecontours; n++ ) {
	  if ( theSingleContour[n].theCnt == NULL ) continue;
	  _PrintContourInCnt( theSingleContour[n].theCnt, 
			      NULL,
			      theSingleContour[n].theCnt->n,
			      (double)theSingleContour[n].z,
			      0, f2d );
	  for (i=0; i<16; i++ ) {
	    fprintf( ftrsf, "%f", theSingleContour[n].mat[i] );
	    if ( i<15 ) fprintf( ftrsf, " " );
	  }
	  fprintf( ftrsf, "\n" );
	}
	
	fclose( f2d );
	fclose( ftrsf);
      }


    }



    
    xr = (double*)malloc( 3*nsinglecontours*sizeof( double ) );
    yr = xr; yr += nsinglecontours;
    zr = yr; zr += nsinglecontours;
    n_realcontours = 0;

    csurf  = (double*)malloc( nsinglecontours*sizeof( double ) );
    ci     = (double*)malloc( nsinglecontours*sizeof( double ) );
    cn = 0;

    /* put contours and matrices
       in millimeters
       measures section area
    */
    for ( n = 0; n < nsinglecontours; n++ ) {
      
      if ( theSingleContour[n].theCnt != NULL ) {

	for ( j=0; j<theSingleContour[n].theCnt->n; j++ ) {
	  theSingleContour[n].theCnt->thePts[j].x *= image->siz.x;
	  theSingleContour[n].theCnt->thePts[j].y *= image->siz.y;
	}
	theSingleContour[n].z                   *= image->siz.z;
	
	if ( theSingleContour[n].mat != NULL ) {
	  theSingleContour[n].mat[ 0] /= image->siz.x;
	  theSingleContour[n].mat[ 4] /= image->siz.x;
	  theSingleContour[n].mat[ 8] /= image->siz.x;

	  theSingleContour[n].mat[ 1] /= image->siz.y;
	  theSingleContour[n].mat[ 5] /= image->siz.y;
	  theSingleContour[n].mat[ 9] /= image->siz.y;

	  theSingleContour[n].mat[ 2] /= image->siz.z;
	  theSingleContour[n].mat[ 6] /= image->siz.z;
	  theSingleContour[n].mat[10] /= image->siz.z;
	}

	theSingleContour[n].sliceCenter[0] = image->siz.x 
	  * ( image->dim.x - 1 ) / 2.0;
	theSingleContour[n].sliceCenter[1] = image->siz.y 
	  * ( image->dim.y - 1 ) / 2.0;
	theSingleContour[n].sliceCorner[0] = image->siz.x * ( - 0.5 );
	theSingleContour[n].sliceCorner[1] = image->siz.y * ( - 0.5 );
	
	xr[ n_realcontours ] = 
	  theSingleContour[n].mat[ 0] * theSingleContour[n].sliceCenter[0] +
	  theSingleContour[n].mat[ 1] * theSingleContour[n].sliceCenter[0] +
	  theSingleContour[n].mat[ 3]; 
	yr[ n_realcontours ] = 
	  theSingleContour[n].mat[ 4] * theSingleContour[n].sliceCenter[0] +
	  theSingleContour[n].mat[ 5] * theSingleContour[n].sliceCenter[0] +
	  theSingleContour[n].mat[ 7]; 
	zr[ n_realcontours ] = 
	  theSingleContour[n].mat[ 8] * theSingleContour[n].sliceCenter[0] +
	  theSingleContour[n].mat[ 9] * theSingleContour[n].sliceCenter[0] +
	  theSingleContour[n].mat[11]; 
	n_realcontours++;

	csurf[ cn ] = surfaceContour2D( theSingleContour[n].theCnt,
					theSingleContour[n].sliceCenter[0],
					theSingleContour[n].sliceCenter[1] );
	ci[ cn ] = n;
	cn ++;


      }
      
    }


    /* transform in 2D for prisme
     */

    if ( par.cnts[0] != '\0' && par.prisme_output == _2D_ ) {

      for ( n = 0; n < nsinglecontours; n++ ) {
	if ( theSingleContour[n].theCnt == NULL ) continue;
	
	next = n ;
	for ( j = n+1; n == next && j < nsinglecontours; j++ ) {
	  if ( theSingleContour[j].theCnt != NULL ) next = j;
	}
	if ( next == n ) continue;
	
	/* on cree les contours
	 */
	_ReCompute2DContours( &(theSingleContour[n]), &(theSingleContour[next]), 
			      par.tolerance_vectorial_product );
      }

      {
        char name[DOUBLESTRINGLENGTH];
	int s = 0;

	FILE *f3d, *fprev, *fnext, *fopen();
	sprintf( name, "%s.3D.cnt", par.cnts );
	f3d = fopen( name, "w" );
	sprintf( name, "%s.prev.cnt", par.cnts );
	fprev = fopen( name, "w" );
	sprintf( name, "%s.next.cnt", par.cnts ); 
	fnext = fopen( name, "w" );

	for ( n = 0; n < nsinglecontours; n++ ) {
	  if ( theSingleContour[n].theCnt == NULL ) continue;
	  s++;
	}
	fprintf( f3d, "S %d\n", s );
	fprintf( fprev, "S %d\n", s );
	fprintf( fnext, "S %d\n", s );

	for ( n = 0; n < nsinglecontours; n++ ) {
	  if ( theSingleContour[n].theCnt == NULL ) continue;
	  _PrintContourInCnt( theSingleContour[n].theCnt, 
			      theSingleContour[n].mat,
			      theSingleContour[n].theCnt->n,
			      (double)theSingleContour[n].z,
			      1, f3d );
	  _PrintContourInCnt( theSingleContour[n].with_prev, 
			      NULL,
			      theSingleContour[n].theCnt->n,
			      (double)theSingleContour[n].z,
			      0, fprev );
	  _PrintContourInCnt( theSingleContour[n].with_next, 
			      NULL,
			      theSingleContour[n].theCnt->n,
			      (double)theSingleContour[n].z,
			      0, fnext );
	}
	
	fclose( f3d );
	fclose( fprev );
	fclose( fnext );
      }


    }




    if ( par.names.ext[0] != '\0' ) {
      char name[DOUBLESTRINGLENGTH];
      FILE *f, *fopen();
      int startname, fd;
      
      startname = strlen( par.names.ext )-1;
      for ( ; startname >= 0 && par.names.ext[startname] != '/' ; startname-- ) 
	;
      
      sprintf( name, "%s.raw", par.names.ext );
      fd = creat( name, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
      
      sprintf( name, "%s.m", par.names.ext );
      f = fopen( name, "w" );
      
      fprintf( f, "\n\n\n" );
      fprintf( f, "%%\n" );
      fprintf( f, "%%" );
      for ( j=0; j<argc; j++ ) fprintf( f, " %s", argv[j] );
      fprintf( f, "\n" );
      fprintf( f, "%%\n" );
      fprintf( f, "\n\n\n" );
      
      fprintf( f, "echo off\n" );
      fprintf( f, "fid = fopen('%s.raw', 'r' );\n", &(par.names.ext[startname+1]) );
      
      fprintf( f, "\n\n\n" );
      if ( write( fd, xr, n_realcontours * sizeof( double ) ) == -1 ) {
	fprintf( stderr, "%s: error when writing\n", program );
      }
      if ( write( fd, yr, n_realcontours * sizeof( double ) ) == -1 ) {
	fprintf( stderr, "%s: error when writing\n", program );
      }
      if ( write( fd, zr, n_realcontours * sizeof( double ) ) == -1 ) {
	fprintf( stderr, "%s: error when writing\n", program );
      }
      fprintf( f, "XCENTER = fread( fid, %d, 'float%lu' );\n", 
	       n_realcontours, 8*sizeof(double) );
      fprintf( f, "YCENTER = fread( fid, %d, 'float%lu' );\n", 
	       n_realcontours, 8*sizeof(double) );
      fprintf( f, "ZCENTER = fread( fid, %d, 'float%lu' );\n", 
	       n_realcontours, 8*sizeof(double) );

      fprintf( f, "\n\n\n" );
      if ( write( fd, ci, cn * sizeof( double ) ) == -1 ) {
	fprintf( stderr, "%s: error when writing\n", program );
      }
      if ( write( fd, csurf, cn * sizeof( double ) ) == -1 ) {
	fprintf( stderr, "%s: error when writing\n", program );
      }
      cmax = csurf[0];
      for ( j=1; j<cn; j++ ) if ( cmax < csurf[j] ) cmax = csurf[j];
      fprintf( f, "CI = fread( fid, %d, 'float%lu' );\n", 
	      cn, 8*sizeof(double) );
      fprintf( f, "CSURF = fread( fid, %d, 'float%lu' );\n", 
	      cn, 8*sizeof(double) );
      fprintf( f, "\n" );
      fprintf( f, "%% figure;\n" );
      fprintf( f, "%% hold on;\n" );
      fprintf( f, "%% plot( CI, CSURF, 'b-' );\n" );
      fprintf( f, "%% hold off;\n" );
      fprintf( f, "\n" );
      fprintf( f, "MAKEDIAM=1;\n" );
      fprintf( f, "MAKECONTOURS=1;\n" );
      fprintf( f, "dfig = figure;\n" );
      fprintf( f, "hold on;\n" );
      fprintf( f, "plot( CI, 2 * sqrt( CSURF / (3.1416) ), 'b-', 'Linewidth', 3 );\n" );
      fprintf( f, "axis( [%f %f 0 %d] );\n", ci[0], ci[cn-1], (int)(2.0*sqrt(cmax/3.1415)+1) );
      fprintf( f, "grid;\n" );
      fprintf( f, "set(gca,'Linewidth', 3, 'FontSize', 16, 'FontWeight', 'bold');\n" );
      fprintf( f, "hold off;\n" );
      fprintf( f, "if MAKEDIAM ~= 0\n" );
      fprintf( f, "   print(dfig, '-djpeg80','%s-diam' );\n", par.names.ext );
      fprintf( f, "end\n" );
	
      fprintf( f, "\n\n\n" );
      fprintf( f, "MAKEMOV1=0;\n" );
      fprintf( f, "MAKEPPM1=0;\n" );
      fprintf( f, "MAKEJPG1=0;\n" );
      fprintf( f, "\n\n\n" );
      fprintf( f, "figure;\n" );
      fprintf( f, "hold on;\n" );
      fprintf( f, "view(3);\n" );
      fprintf( f, "grid;\n" );
      fprintf( f, "plot3( XCENTER, YCENTER, ZCENTER, 'r', 'Linewidth', 3 );\n" );
      fprintf( f, "axis equal;\n" );
      fprintf( f, "axis vis3d;\n" );
      fprintf( f, "hold off;\n" );
      fprintf( f, "\n\n\n" );

      fprintf( f, "fig = figure;\n" );
      fprintf( f, "hold on;\n" );
      if ( mat != NULL ) {
	fprintf( f, "view(3);\n" );
	fprintf( f, "%% axis equal;\n" );
	fprintf( f, "%% axis vis3d;\n" );     
	fprintf( f, "grid;\n" );
      }   
      fprintf( f, "plot3( XCENTER, YCENTER, ZCENTER, 'r', 'Linewidth', 3 );\n" );
      fprintf( f, "if MAKEMOV1 ~= 0\n" );
      fprintf( f, "   MOV1(1) = getframe(fig);\n" );
      fprintf( f, "end\n" );
      fprintf( f, "\n" );
      
      if ( mat == NULL && par.min_slice == par.max_slice )
	VT_2DDrawImage( image, par.min_slice, fd, f );


      for ( m=1, n = 0; n < nsinglecontours; m++, n++ ) {

	if ( theSingleContour[n].theCnt != NULL ) {
	
	  next = n ;
	  for ( j = n+1; n == next && j < nsinglecontours; j++ ) {
	    if ( theSingleContour[j].theCnt != NULL ) next = j;
	  }
	  
	  printf( "z = %3d contour #%3d (n=%d)\n",
		  theSingleContour[n].z, theSingleContour[n].icnt,
		  theSingleContour[n].theCnt->n );
	  
	  
	  if ( theSingleContour[n].mat == NULL ) {
	    if ( par.min_slice == par.max_slice )
	      MAT_DrawContour2D( theSingleContour[n].theCnt, fd, f, theSingleContour[n].icnt );
	    else 
	      MAT_DrawContour2Din3D( theSingleContour[n].theCnt, 
				     (double)0.0, identity,
				     fd, f, theSingleContour[n].icnt, NULL );
	  }
	  else {
	    MAT_DrawContour2Din3D( theSingleContour[n].theCnt,
				   (double)0.0, theSingleContour[n].mat,
				   fd, f, theSingleContour[n].icnt, NULL );
	    
	    if ( par.cnts[0] != '\0' ) {
	      if ( par.draw_frames ||
		   (par.max_slice - par.min_slice >= 1 
		    && par.max_slice - par.min_slice <= 3) ) {
		_DrawSliceFrame( &(theSingleContour[n]), f );
		_DrawSliceBorders( &(theSingleContour[n]), f );
		if ( next != n ) {
		  _DrawJointFrames( &(theSingleContour[n]), &(theSingleContour[next]), 
				    f, par.tolerance_vectorial_product );
		}
	      }
	    } /* if ( par.cnts[0] != '\0' ) */
	    
	  }

	}
	fprintf( f, "if MAKEPPM1 ~= 0\n" );
	fprintf( f, "   print(fig, '-dppm','vessel-contours-%d' );\n", m+1 );
	fprintf( f, "end\n" );
	fprintf( f, "if MAKEJPG1 ~= 0\n" );
	fprintf( f, "   print(fig, '-djpeg80','vessel-contours-%d' );\n", m+1 );
	fprintf( f, "end\n" );
	fprintf( f, "if MAKEMOV1 ~= 0\n" );
	fprintf( f, "   MOV1(%d) = getframe(fig);\n", m+1 );
	fprintf( f, "end\n" );
	fprintf( f, "\n" );

      }

      fprintf( f, "\n" );
      fprintf( f, "axis equal;\n" );
      fprintf( f, "axis vis3d;\n" );
      fprintf( f, "set(gca,'Linewidth', 3, 'FontSize', 16, 'FontWeight', 'bold');\n" );
      fprintf( f, "hold off;\n" );
      fprintf( f, "if MAKECONTOURS ~= 0\n" );
      fprintf( f, "   print(fig, '-djpeg80','%s-cont' );\n", par.names.ext );
      fprintf( f, "end\n" );
      fprintf( f, "fclose( fid);\n" );

      fprintf( f, "\n\n\n" );
      fprintf( f, "if MAKEMOV1 ~= 0\n" );
      fprintf( f, "   movie2avi( MOV1, 'vessel-contours.avi' );\n" );
      fprintf( f, "   figure;\n" );
      fprintf( f, "   movie(MOV1);\n" );
      fprintf( f, "\n\n\n" );
      fprintf( f, "\n" );
      fprintf( f, "end\n" );

      fprintf( f, "\n\n\n" );
      
      if ( 0 ) {

	fprintf( f, "fim = fopen('%s', 'r' );\n", par.names.in );
	fprintf( f, "HEADER = fread( fim, 256, 'uint8' );\n" );
	fprintf( f, "\n\n\n" );
	
	switch( image->type ) {
	default :
	  break;
	case UCHAR :  sprintf( matlab_type, "uint8" );   break;
	case USHORT : sprintf( matlab_type, "uint16" );  break;
	case SSHORT : sprintf( matlab_type, "int16" );   break;
	case FLOAT :  sprintf( matlab_type, "float32" ); break;
	case DOUBLE : sprintf( matlab_type, "float64" ); break;
	}
	
	for ( n=0; n<par.min_slice; n++ )
	  fprintf( f, "FOO = fread( fim, [%lu %lu], '%s' );", image->dim.x, image->dim.y, matlab_type );

	fprintf( f, "XSLICE = ones(%lu,1) * [1:%lu] * %f;\n", image->dim.y, image->dim.x, image->siz.x ); 
	fprintf( f, "YSLICE = [1:%lu]' * ones(1,%lu) * %f;\n", image->dim.y, image->dim.x, image->siz.y  ); 
	fprintf( f, "ZSLICE = ones(%lu,1) * ones(1,%lu);\n", image->dim.y, image->dim.x ); 
	fprintf( f, "\n" );
	
	for ( n = 0; n < nsinglecontours; n++ ) {
	  fprintf( f, "\n" );
	  fprintf( f, "XSLICE%d = XSLICE;\n", n );
	  fprintf( f, "YSLICE%d = YSLICE;\n", n );
	  fprintf( f, "ZSLICE%d = ZSLICE;\n", n );
	  fprintf( f, "for i=1:%lu\n", image->dim.x );
	  fprintf( f, "for j=1:%lu\n", image->dim.y );
	  fprintf( f, "  XSLICE%d(j,i) = %f * XSLICE(j,i) + %f * YSLICE(j,i) + %f;\n", 
		   n, theSingleContour[n].mat[0], theSingleContour[n].mat[1], theSingleContour[n].mat[3] );
	  fprintf( f, "  YSLICE%d(j,i) = %f * XSLICE(j,i) + %f * YSLICE(j,i) + %f;\n", 
		   n, theSingleContour[n].mat[4], theSingleContour[n].mat[5], theSingleContour[n].mat[7] );
	  fprintf( f, "  ZSLICE%d(j,i) = %f * XSLICE(j,i) + %f * YSLICE(j,i) + %f;\n", 
		   n, theSingleContour[n].mat[8], theSingleContour[n].mat[9], theSingleContour[n].mat[11] );
	  fprintf( f, "end;\n" );
	  fprintf( f, "end;\n" );
	  
	  fprintf( f, "CSLICE%d = fread( fim, [%lu %lu], '%s' );\n", 
		   n, image->dim.x, image->dim.y, matlab_type );
	}

	fprintf( f, "\n" );
	fprintf( f, "fclose( fim);\n" );

	fprintf( f, "\n\n\n" );
	fprintf( f, "MAKEMOV2=0;\n" );
	fprintf( f, "MAKEPPM2=0;\n" );
	fprintf( f, "MAKEJPG2=0;\n" );
	fprintf( f, "\n\n\n" );
	fprintf( f, "fig2 = figure;\n" );
	fprintf( f, "hold on;\n" );
	fprintf( f, "colormap( 'gray' );\n" );
	if ( mat != NULL ) {
	  fprintf( f, "view(3);\n" );
	  fprintf( f, "%% axis equal;\n" );
	  fprintf( f, "%% axis vis3d;\n" );     
	  fprintf( f, "grid;\n" );
	}
	fprintf( f, "plot3( XCENTER, YCENTER, ZCENTER, 'r', 'Linewidth', 3 );\n" );
	fprintf( f, "if MAKEMOV2 ~= 0\n" );
	fprintf( f, "   MOV2(1) = getframe(fig2);\n"  );
	fprintf( f, "end\n" );
	fprintf( f, "\n" );

	for ( m=1, n = 0; n < nsinglecontours; m++, n+=n_step ) {
	  if ( 1 || theSingleContour[n].theCnt != NULL ) {
	    fprintf( f, "surf( XSLICE%d, YSLICE%d, ZSLICE%d, CSLICE%d' );\n", n, n, n, n );
	    fprintf( f, "shading interp;\n" );
	  }
	fprintf( f, "if MAKEPPM2 ~= 0\n" );
	fprintf( f, "   print(fig2, '-dppm','vessel-slices-%d' );\n", m+1 );
	fprintf( f, "end\n" );
	fprintf( f, "if MAKEJPG2 ~= 0\n" );
	fprintf( f, "   print(fig2, '-djpeg80','vessel-slices-%d' );\n", m+1 );
	fprintf( f, "end\n" );
	fprintf( f, "if MAKEMOV2 ~= 0\n" );
	fprintf( f, "   MOV2(%d) = getframe(fig2);\n", m+1 );
	fprintf( f, "end\n" );
	fprintf( f, "\n" );
	}

	fprintf( f, "\n" );
	fprintf( f, "hold off;\n" );
	fprintf( f, "\n\n\n" );
	fprintf( f, "\n\n\n" );
	fprintf( f, "if MAKEMOV2 ~= 0\n" );
	fprintf( f, "   movie2avi( MOV2, 'vessel-slices.avi' );\n" );
	fprintf( f, "   figure;\n" );
	fprintf( f, "   movie(MOV2);\n" );
	fprintf( f, "\n\n\n" );
	fprintf( f, "\n" );
	fprintf( f, "end\n" );

      }

    }

  }

  /*--- liberations memoires ---*/
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


      else if ( strcmp ( argv[i], "-th" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -th...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->threshold) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -th...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-slice" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -slice...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->min_slice) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -slice...\n", 0 );
	par->max_slice = par->min_slice;
      }
      else if ( strcmp ( argv[i], "-min-slice" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -min-slice...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->min_slice) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -min-slice...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-max-slice" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -max-slice...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->max_slice) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -max-slice...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-matlab" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -matlab...\n", 0 );
	strncpy( par->names.ext, argv[i], STRINGLENGTH ); 
      }
      else if ( strcmp ( argv[i], "-draw-frames" ) == 0 ) {
	par->draw_frames = 1;
      }
      else if ( strcmp ( argv[i], "-trsfs" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -trsfs...\n", 0 );
	strncpy( par->trsfs, argv[i], STRINGLENGTH ); 
      }

      else if ( strcmp ( argv[i], "-cnts" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -cnts...\n", 0 );
	strncpy( par->cnts, argv[i], STRINGLENGTH ); 
      }
      else if ( strcmp ( argv[i], "-2D" ) == 0 ) {
	par->prisme_output = _2D_;
      }
      else if ( strcmp ( argv[i], "-3D" ) == 0 ) {
	par->prisme_output = _3D_;
      }
      
      else if ( strcmp ( argv[i], "-tol" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -tol...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->tolerance_vectorial_product) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -tol...\n", 0 );
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

  par->min_slice = 0;
  par->max_slice = 10000;
  par->threshold = 100;

  par->trsfs[0] = '\0';

  par->draw_frames = 0;

  par->cnts[0] = '\0';
  par->tolerance_vectorial_product = 0.05;
  par->prisme_output = _2D_;
}




/***********************************************************************
           b     blue          .     point              -     solid
           g     green         o     circle             :     dotted
           r     red           x     x-mark             -.    dashdot 
           c     cyan          +     plus               --    dashed   
           m     magenta       *     star
           y     yellow        s     square
           k     black         d     diamond
                               v     triangle (down)
                               ^     triangle (up)
                               <     triangle (left)
                               >     triangle (right)
                               p     pentagram
                               h     hexagram
**********************************************************************/

void _MatVect( double *res, double *mat, double *vect )
{
  res[0] = mat[ 0] *vect[0] + mat[ 1] * vect[1] + mat[ 2] * vect[2] + mat [ 3];
  res[1] = mat[ 4] *vect[0] + mat[ 5] * vect[1] + mat[ 6] * vect[2] + mat [ 7];
  res[2] = mat[ 8] *vect[0] + mat[ 9] * vect[1] + mat[10] * vect[2] + mat [11];
}


void _DrawSliceBorders( typeSingleContour *s, FILE *f )
{
  double pt[3];
  double pt00[3], pt01[3], pt11[3], pt10[3];
  
  pt[0] = s->sliceCorner[0];
  pt[1] = s->sliceCorner[1];
  pt[2] = 0;
  
  _MatVect( pt00, s->mat, pt );

  pt[0] +=  2 * s->sliceCenter[0];

  _MatVect( pt01, s->mat, pt );

  pt[1] +=  2 * s->sliceCenter[1];
  
  _MatVect( pt11, s->mat, pt );
  
  pt[0] -=  2 * s->sliceCenter[0];
 
  _MatVect( pt10, s->mat, pt );

  fprintf( f, "plot3( [%f %f %f %f %f], [%f %f %f %f %f], [%f %f %f %f %f], 'k-' );\n",
	   pt00[0], pt01[0], pt11[0], pt10[0], pt00[0],
	   pt00[1], pt01[1], pt11[1], pt10[1], pt00[1],
	   pt00[2], pt01[2], pt11[2], pt10[2], pt00[2] );
  fprintf( f, "\n" );

}
  
void _DrawSliceFrame( typeSingleContour *s, FILE *f )
{
  double pt[3];
  double ptc[3];
  double l = 1;

  pt[0] = s->sliceCenter[0];
  pt[1] = s->sliceCenter[1];
  pt[2] = 0;

  _MatVect( ptc, s->mat, pt );

  fprintf( f, "plot3( [%f %f], [%f %f], [%f %f], 'r-', 'Linewidth', 3 );\n",
	   ptc[0], ptc[0] + l * s->mat[2],
	   ptc[1], ptc[1] + l * s->mat[6],
	   ptc[2], ptc[2] + l * s->mat[10] );
  fprintf( f, "plot3( [%f %f], [%f %f], [%f %f], 'b-', 'Linewidth', 2 );\n",
	   ptc[0], ptc[0] + l * s->mat[0],
	   ptc[1], ptc[1] + l * s->mat[4],
	   ptc[2], ptc[2] + l * s->mat[8] );
  fprintf( f, "plot3( [%f %f], [%f %f], [%f %f], 'b--', 'Linewidth', 2 );\n",
	   ptc[0], ptc[0] + l * s->mat[1],
	   ptc[1], ptc[1] + l * s->mat[5],
	   ptc[2], ptc[2] + l * s->mat[9] );
}


int _ComputeJointFrames( typeSingleContour *s1, /* prev */
			 typeSingleContour *s2, /* next */
			 double *c1,
			 double *x1,
			 double *y1,
			 double *c2,
			 double *x2,
			 double *y2,
			 double tolerance )
{
  double *mat1 = s1->mat;
  double *mat2 = s2->mat;
  double n;
  double rc1[3];
  double rc2[3];
  double inv[16];
  double pt[3];
  
  double a, dp;
  double v[3];

  double tolerance_scalar_product = 0.001;

  /* y' = z1 vectoriel z2
   */
  y1[0] = mat1[ 6] * mat2[10] - mat1[10] * mat2[ 6];
  y1[1] = mat1[10] * mat2[ 2] - mat1[ 2] * mat2[10];
  y1[2] = mat1[ 2] * mat2[ 6] - mat1[ 6] * mat2[ 2];

  n = sqrt( y1[0]*y1[0] + y1[1]*y1[1] + y1[2]*y1[2] );

  /* plans paralleles ou presque,
     les deux bases doivent etre quasiment identiques
     par construction
   */
  if ( n < tolerance ) {


    x1[0] = mat1[0];
    x1[1] = mat1[4];
    x1[2] = mat1[8];

    y1[0] = mat1[1];
    y1[1] = mat1[5];
    y1[2] = mat1[9];
    
    x2[0] = mat2[0];
    x2[1] = mat2[4];
    x2[2] = mat2[8];

    y2[0] = mat2[1];
    y2[1] = mat2[5];
    y2[2] = mat2[9];
    
    /* on projette c1 sur le second plan
     */
    
    dp = mat1[2] * mat2[2] + mat1[6] * mat2[6] + mat1[10] * mat2[10];
    if ( fabs( dp ) < tolerance_scalar_product ) {
      fprintf( stderr, "produit scalaire trop faible ( ?? ) = %g\n", dp );
      return( -1 );
    }

    c1[0] = c1[1] = c1[2] = 0.0;
    _MatVect( rc1, mat1, c1 );
    _MatVect( rc2, mat2, s2->sliceCenter );
  
    v[0] = rc2[0] - rc1[0];
    v[1] = rc2[1] - rc1[1];
    v[2] = rc2[2] - rc1[2];
    
    a = (v[0] * mat2[2] + v[1] * mat2[6] + v[2] * mat2[10]) / dp;
    
    pt[0] = rc1[0] + a * mat1[2];
    pt[1] = rc1[1] + a * mat1[6];
    pt[2] = rc1[2] + a * mat1[10];
    
    InverseMat4x4( mat2, inv );
    _MatVect( c2, inv, pt );
    
    return( 0 );
  }

  
  /* normalisation de y
   */
  y1[0] /= n;
  y1[1] /= n;
  y1[2] /= n;

  y2[0] = y1[0];
  y2[1] = y1[1];
  y2[2] = y1[2];

  /* x'i = y' vectoriel zi
   */  
  x1[0] = y1[1] * mat1[10] - y1[2] * mat1[ 6];
  x1[1] = y1[2] * mat1[ 2] - y1[0] * mat1[10];
  x1[2] = y1[0] * mat1[ 6] - y1[1] * mat1[ 2];

  n = sqrt( x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2] );

  x1[0] /= n;
  x1[1] /= n;
  x1[2] /= n;

  x2[0] = y2[1] * mat2[10] - y2[2] * mat2[ 6];
  x2[1] = y2[2] * mat2[ 2] - y2[0] * mat2[10];
  x2[2] = y2[0] * mat2[ 6] - y2[1] * mat2[ 2];

  n = sqrt( x2[0]*x2[0] + x2[1]*x2[1] + x2[2]*x2[2] );

  x2[0] /= n;
  x2[1] /= n;
  x2[2] /= n;

  /* on cherche un point de la droite commune
     par exemple un point du plan 1
     M = c1 + a x'1
     l'equation du plan 2 est  Mc2 . z2 = 0
     donc a = c1c2 . z2 / x'1 . z2
  */

  dp = x1[0] * mat2[2] + x1[1] * mat2[6] + x1[2] * mat2[10];
  if ( fabs( dp ) < tolerance_scalar_product ) {
    fprintf( stderr, "produit scalaire trop faible = %g\n", dp );
    return( -1 );
  }
  
  _MatVect( rc1, mat1, s1->sliceCenter );
  _MatVect( rc2, mat2, s2->sliceCenter );
  
  v[0] = rc2[0] - rc1[0];
  v[1] = rc2[1] - rc1[1];
  v[2] = rc2[2] - rc1[2];

  a = (v[0] * mat2[2] + v[1] * mat2[6] + v[2] * mat2[10]) / dp;

  pt[0] = rc1[0] + a * x1[0];
  pt[1] = rc1[1] + a * x1[1];
  pt[2] = rc1[2] + a * x1[2];

  InverseMat4x4( mat1, inv );
  _MatVect( c1, inv, pt );
  c1[2] = 0;
  
  InverseMat4x4( mat2, inv );
  _MatVect( c2, inv, pt );
  c2[2] = 0;
  
  return( 1 );
}



void _DrawJointFrames( typeSingleContour *s1,
		       typeSingleContour *s2,
		       FILE *f, double tol )
{
  double c1[3] = {0.0,0.0,0.0};
  double x1[3] = {0.0,0.0,0.0};
  double y1[3] = {0.0,0.0,0.0};
  double c2[3] = {0.0,0.0,0.0};
  double x2[3] = {0.0,0.0,0.0};
  double y2[3] = {0.0,0.0,0.0};

  double pt[3];
  double l = 1;
  double c = 2;

  _ComputeJointFrames( s1, s2, c1, x1, y1, c2, x2, y2, tol );
  
  _MatVect( pt, s1->mat, c1 );

  fprintf( f, "plot3( [%f %f], [%f %f], [%f %f], 'r-', 'Linewidth', 3 );\n",
	   pt[0], pt[0] + l * s1->mat[2],
	   pt[1], pt[1] + l * s1->mat[6],
	   pt[2], pt[2] + l * s1->mat[10] );

  fprintf( f, "plot3( [%f %f], [%f %f], [%f %f], 'g-', 'Linewidth', 3 );\n",
	   pt[0], pt[0] + l * x1[0],
	   pt[1], pt[1] + l * x1[1],
	   pt[2], pt[2] + l * x1[2] );
  
  fprintf( f, "plot3( [%f %f], [%f %f], [%f %f], 'b--', 'Linewidth', 3 );\n",
	   pt[0], pt[0] + c * l * y1[0],
	   pt[1], pt[1] + c * l * y1[1],
	   pt[2], pt[2] + c * l * y1[2] );
  
  _MatVect( pt, s2->mat, c2 );

  fprintf( f, "plot3( [%f %f], [%f %f], [%f %f], 'r-', 'Linewidth', 3 );\n",
	   pt[0], pt[0] + l * s2->mat[2],
	   pt[1], pt[1] + l * s2->mat[6],
	   pt[2], pt[2] + l * s2->mat[10] );
  
  fprintf( f, "plot3( [%f %f], [%f %f], [%f %f], 'g-', 'Linewidth', 3 );\n",
	   pt[0], pt[0] + l * x2[0],
	   pt[1], pt[1] + l * x2[1],
	   pt[2], pt[2] + l * x2[2] );
  
  fprintf( f, "plot3( [%f %f], [%f %f], [%f %f], 'c-', 'Linewidth', 3 );\n",
	   pt[0], pt[0] + c * l * y2[0],
	   pt[1], pt[1] + c * l * y2[1],
	   pt[2], pt[2] + c * l * y2[2] );
  

}


void _ReComputeOne2DContour( typeSingleContour *s,
			     typeContour2D *old,
			     typeContour2D *new,
			     double *c,
			     double *x,
			     double *y )
{
  double mat[6];
  int n;

  mat[0] = x[0]*s->mat[0] + x[1]*s->mat[4] + x[2]*s->mat[8];
  mat[1] = y[0]*s->mat[0] + y[1]*s->mat[4] + y[2]*s->mat[8];
  mat[2] = c[0] - mat[0] * c[0] - mat[1] * c[1];

  mat[3] = x[0]*s->mat[1] + x[1]*s->mat[5] + x[2]*s->mat[9];
  mat[4] = y[0]*s->mat[1] + y[1]*s->mat[5] + y[2]*s->mat[9];
  mat[5] = c[1] - mat[3] * c[0] - mat[4] * c[1];
  
  for ( n=0; n<old->n; n++ ) {
    new->thePts[n].x = mat[0]*old->thePts[n].x + mat[1]*old->thePts[n].y + mat[2];
    new->thePts[n].y = mat[3]*old->thePts[n].x + mat[4]*old->thePts[n].y + mat[5];
  }
}


int _ReCompute2DContours( typeSingleContour *s1, /* prev */
			  typeSingleContour *s2, /* next */
			  double tol )
{
  char *proc = "_ReCompute2DContours";
  double c1[3] = {0.0,0.0,0.0};
  double x1[3] = {0.0,0.0,0.0};
  double y1[3] = {0.0,0.0,0.0};
  double c2[3] = {0.0,0.0,0.0};
  double x2[3] = {0.0,0.0,0.0};
  double y2[3] = {0.0,0.0,0.0};

  s1->with_next = copyContour2D( s1->theCnt );
  s2->with_prev = copyContour2D( s2->theCnt );

  if ( s1->with_next == NULL || s2->with_prev == NULL ) {
    fprintf( stderr, "allocation failed\n" );
    return( 0 );
  }
  
  switch( _ComputeJointFrames( s1, s2, c1, x1, y1, c2, x2, y2, tol ) ) {
  default :
    fprintf( stderr, "%s: case not handled in switch ??? \n", proc );
    return( 0 );
  case -1 :
    fprintf( stderr, "%s: this case should not occur ... \n", proc );
    return( 0 );
  case 0 :
    fprintf( stderr, "%s: parallel slices\n", proc );
    break;
  case 1 :
    break;
  }
  return( 1 );
  _ReComputeOne2DContour( s1, s1->theCnt, s1->with_next, c1, x1, y1 );
  _ReComputeOne2DContour( s2, s2->theCnt, s2->with_prev, c2, x2, y2 );
  return( 1 );
}



void _PrintContourInCnt( typeContour2D *c,
			 double *mat,
			 int nvertices,
			 double z,
			 int print_z, /* print z coordinates of points */
			 FILE *f )
{
  int i;
  fprintf( f, "v %d z %g\n", nvertices, z );
  fprintf( f, "{\n" );

  if ( c == NULL ) {
    for (i=0; i<nvertices; i++ ) {
      fprintf( f, "0 0" );
      if ( print_z ) fprintf( f, " 0" );
      fprintf( f, "\n" );
    }
  }
  else {

    if ( mat == NULL ) {
      for (i=0; i<c->n; i++ ) {
	fprintf( f, "%g %g", c->thePts[i].x, c->thePts[i].y );
	if ( print_z ) fprintf( f, " %g", z );
	fprintf( f, "\n" );
      }
    }
    else {
      for (i=0; i<c->n; i++ ) {
	fprintf( f, "%g ",
		 mat[0]*c->thePts[i].x + mat[1]*c->thePts[i].y + mat[ 3] );
	fprintf( f, "%g",
		 mat[4]*c->thePts[i].x + mat[5]*c->thePts[i].y + mat[ 7] );
	if ( print_z ) 
	  fprintf( f, " %g",
		   mat[8]*c->thePts[i].x + mat[9]*c->thePts[i].y + mat[11] );
	fprintf( f, "\n" );
      }
    }

  }
  fprintf( f, "}\n" );
}













#define TINY 1e-12
int InverseMat4x4( double *matrice, double *inv )
{
  register int i, j, k;
  int kmax, rang = 4;
  register double c, max;
  double mat [16];
  
  for (i=0; i<16; i++ ) {
    mat[i] = matrice[i] ;
    inv[i] = 0.0;
  }
  inv[0] = inv[5] = inv[10] = inv[15] = 1.0;
  
  for ( j=0; j<4; j++ ) {
    if ( (mat[j*4+j] > (-TINY)) && (mat[j*4+j] < TINY) ) {
      /* recherche du plus grand element non nul sur la colonne j */
      kmax = j;
      max = 0.0;
      for (k=j+1; k<4; k++ ) {
	c = ( mat[k*4+j] > 0.0 ) ? mat[k*4+j] : (-mat[k*4+j]) ;
	if ( (c > TINY) && (c > max) ) { max = c; kmax = k; }
      }
      if ( kmax == j ) {
	/* la ligne est nulle */
	rang --;
      } else {
	/* sinon, on additionne */
	for ( i=0; i<4; i++ ) {
	  mat[j*4+i] += mat[kmax*4+i];
	  inv[j*4+i] += inv[kmax*4+i];
	}
      }
    }
    if ( (mat[j*4+j] < (-TINY)) || (mat[j*4+j] > TINY) ) {
      /* les autres lignes */
      for (k=0; k<4; k++) {
	if ( k != j ) {
	  c = mat[k*4 + j] / mat[j*4 + j];
	  for ( i=0; i<4; i++ ) {
	    mat[k*4 + i] -= c * mat[j*4 + i];
	    inv[k*4 + i] -= c * inv[j*4 + i];
	  }
	}
      }
      /* la ligne */
      c = mat[j*4 + j];
      for ( i=0; i<4; i++ ) {
	mat[j*4 + i] /= c;
	inv[j*4 + i] /= c;
      }
    }
  }

  return( rang );
}
