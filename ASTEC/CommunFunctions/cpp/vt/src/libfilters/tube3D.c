/*************************************************************************
 * tube3D.c -
 *
 * $Id: tube3D.c,v 1.8 2001/06/08 15:10:50 greg Exp $
 *
 * Copyright (c) INRIA 2000
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Mon Sep 11 18:20:24 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 */

#include <tube3D.h>


static int _verbose_ = 0;


void _VerboseInTube3D()
{
  if ( _verbose_ <= 0 ) _verbose_ = 1;
  else _verbose_ ++;
}
void _NoVerboseInTube3D()
{
  _verbose_ = 0;
}


/* number of added points at each border
   of a line for the recursive filtering
*/
static int border_length[3] = { 20, 20, 20 };







/* pre-computation for the response computation
   c2 <- cos( i * 2 * pi / (| radius * 2 * pi | + 1)
   c3 <- sin( i * 2 * pi / (| radius * 2 * pi | + 1)
   the two other fields are used during the response 
   computation
   f  <- 1 if the point belongs to the image (0 else)
   v  <- computed response value
 */
typedef struct {
  double c2;
  double c3;
  double v;
  char   f;
} typeCirclePoint;








static void compute_3D_line_response_in_one_slice( float *theXX, 
						   float *theYY, 
						   float *theZZ,
						   float *theXY, 
						   float *theXZ, 
						   float *theYZ, 
						   float *theX,  
						   float *theY,  
						   float *theZ,
						   typeCirclePoint *thePts,  
						   int nbPts, double radius,
						   int dimx, int dimy, 
						   int dimz, int slice,
						   float *theRep, 
						   float *theTheta, 
						   float *thePhi );

static void unit_vector_to_spherical_angles( const double *v, 
					     double *theta,
					     double *phi );
static void spherical_angles_to_units_vectors( const double theta,
					       const double phi,
					       double *v1, 
					       double *v2,
					       double *v3 );




















int combine_3D_line_response_at_several_scales( float *resResponse,
						float *resTheta,
						float *resPhi,
						float **bufResponse,
						float **bufTheta,
						float **bufPhi,
						int *theDim,
						int nbScales )
{
  char *proc = "combine_3D_line_response_at_several_scales";
  int n, i, v;

  /* test of args
   */
  if ( bufResponse == NULL || 
       bufTheta == NULL || bufPhi == NULL || resResponse == NULL || 
       resTheta == NULL || resPhi == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: NULL buffers in args\n", proc );
    return( -1 );
  }
  
  if ( nbScales <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: no specified scales\n", proc );
    return( -1 );
  }

  if ( theDim[2] == 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to deal with 2D images\n", proc );
    return( -1 );
  }

  v = theDim[0]*theDim[1]*theDim[2];

  for ( n = 0; n < nbScales; n++ ) {
    
    if ( bufResponse[n] == NULL ||  bufTheta[n] == NULL || bufPhi[n] == NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: NULL buffers in args for scale #%2d\n", proc, n );
      return( -1 );
    }

    if ( n == 0 ) {
      (void)memcpy( resResponse, bufResponse[n], v*sizeof( float ) );
      (void)memcpy( resTheta,    bufTheta[n],    v*sizeof( float ) );
      (void)memcpy( resPhi,      bufPhi[n],    v*sizeof( float ) );
    } else {
      for ( i=0; i<v; i++ ) {
	if ( resResponse[i] >= bufResponse[n][i] ) continue;
	resResponse[i] = bufResponse[n][i];
	resTheta[i]    = bufTheta[n][i];
	resPhi[i]      = bufPhi[n][i];
      }
    }

  }

  return( 1 );
}












int compute_3D_line_response_at_several_scales( void *theBuf,
						bufferType typeBuf,
						float **bufResponse,
						float **bufTheta,
						float **bufPhi,
						int *theDim,
						double *scale,
						int nbScales )
{
  char *proc = "compute_3D_line_response_at_several_scales";
  int n;

  /* test of args
   */
  if ( theBuf == NULL || bufResponse == NULL || 
       bufTheta == NULL || bufPhi == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: NULL buffers in args\n", proc );
    return( -1 );
  }
  
  if ( nbScales <= 0 ||  scale == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: no specified scales\n", proc );
    return( -1 );
  }

  if ( theDim[2] == 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to deal with 2D images\n", proc );
    return( -1 );
  }

  for ( n = 0; n < nbScales; n++ ) {
    
    if ( bufResponse[n] == NULL ||  bufTheta[n] == NULL || bufPhi[n] == NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: NULL buffers in args for scale #%2d\n", proc, n );
      return( -1 );
    }
    if ( compute_3D_line_response_at_one_scale( theBuf, typeBuf, 
						bufResponse[n], bufTheta[n], 
						bufPhi[n], theDim, scale[n] ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: computation failed at scale #%2d = %f\n", proc, n, scale[n] );
      
    }
    return( -1 );
  }

  return( 1 );
}






/* calcul incremental d'une reponse multi-echelle

   
 */




int compute_3D_line_response_at_multiple_scales( void *theBuf,
						 bufferType typeBuf,
						 float *bufResponse,
						 float *bufTheta,
						 float *bufPhi,
						 int *theDim,
						 double *scale,
						 int nbScales )
{
  char *proc = "compute_3D_line_response_at_multiple_scales";
  int n, i, v;

  float *allocatedBuf = NULL;
  float *tmpResponse;
  float *tmpTheta;
  float *tmpPhi;

  /* test of args
   */
  if ( theBuf == NULL || bufResponse == NULL || 
       bufTheta == NULL || bufPhi == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: NULL buffers in args\n", proc );
    return( -1 );
  }
  
  if ( nbScales <= 0 ||  scale == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: no specified scales\n", proc );
    return( -1 );
  }

  if ( theDim[2] == 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to deal with 2D images\n", proc );
    return( -1 );
  }

  /* one scale 
   */
  if ( compute_3D_line_response_at_one_scale( theBuf, typeBuf, 
					      bufResponse, bufTheta, bufPhi, theDim, scale[0] ) != 1 ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: computation failed at scale %f\n", proc, scale[0] );
    }
    return( -1 );
  }

  if ( nbScales == 1 ) return( 1 );
    

  /* several scales
   */
  
  allocatedBuf = (float*)malloc( 3 * theDim[0]*theDim[1]*theDim[2] * sizeof(float) );
  if ( allocatedBuf == NULL ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: auxiliary buffer allocation failed\n", proc );
    }
    return( -1 );
  }

  tmpResponse = allocatedBuf;  
  tmpTheta    = tmpResponse;   tmpTheta += theDim[0]*theDim[1]*theDim[2];
  tmpPhi      = tmpResponse;   tmpPhi   += theDim[0]*theDim[1]*theDim[2];




  /* on choisit une approche incrementale pour construire 
     la reponse multi-echelle
  */
  v = theDim[0]*theDim[1]*theDim[2];
  for ( n = 1; n < nbScales; n++ ) {
    
    if ( compute_3D_line_response_at_one_scale( theBuf, typeBuf, 
					       tmpResponse, tmpTheta, tmpPhi, theDim, scale[n] ) != 1 ) {
      if ( _verbose_ ) {
	fprintf( stderr, "%s: computation failed at scale #%2d = %f\n", proc, n, scale[n] );
      }
      free( allocatedBuf );
      return( -1 );
    }
    
    /* maxima computation
     */
    for (i=0; i<v; i++ ) {
      if ( bufResponse[i] >= tmpResponse[i] ) continue;
      bufResponse[i] = tmpResponse[i];
      bufTheta[i]    = tmpTheta[i];
      bufPhi[i]      = tmpPhi[i];
    }

  }




  free( allocatedBuf );
  return( 1 );
}















/* calcul de la reponse a une echelle

   le vecteur 'direction de la structure tubulaire'
   est donne par 2 angles (theta, phi) et vaut

   ( x )   ( cos theta sin phi )
   ( y ) = ( sin theta sin phi )
   ( z )   ( cos phi )

   
 */
int compute_3D_line_response_at_one_scale( void *theBuf,
					   bufferType typeBuf,
					   float *bufResponse,
					   float *bufTheta,
					   float *bufPhi,
					   int *theDim,
					   double scale )
{
  char *proc = "compute_3D_line_response_at_one_scale";

  float *allocatedBuf = NULL;
  bufferType typeAllocatedBuf = FLOAT;

  float *buf_X;
  float *buf_Y;
  float *buf_Z;
  float *tmp_Z0;
  float *tmp_Z1;
  float *buf_ZZ;

  float *buf_XX;
  float *buf_YY;
  float *buf_XY;
  float *buf_XZ;
  float *buf_YZ;

  int dimx = theDim[0];
  int dimy = theDim[1];
  int dimz = theDim[2];
  int dimxy = dimx * dimy;

  float theCoeffs[3];
  filterType theFilter = GAUSSIAN_DERICHE;
  derivativeOrder theDerivatives[3] = { NODERIVATIVE, NODERIVATIVE, NODERIVATIVE };
  int dim[3];

  double theta = 1.0;
  int nb_circle_points = 0;
  typeCirclePoint *circle_points = NULL;
  int p;

  int slice;

  float *tmpResponse;



  /* test of args
   */
  if ( theBuf == NULL || bufResponse == NULL || bufTheta == NULL || bufPhi == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: NULL buffers in args\n", proc );
    return( -1 );
  }
  
  if ( scale <= 0.0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: can not deal with negative scale (%f)\n", proc, scale );
    return( -1 );
  }

  if ( theDim[2] == 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to deal with 2D images\n", proc );
    return( -1 );
  }


  /* pour faire les calculs, on a besoin 
     en 3D :
     - des 3 images 3D gradients Ix, Iy, Iz
     - de  2 images 3D intermediaires contenant
       * lissage selon Z      => Ixx, Iyy, Ixy
       * derivee 1ere selon Z => Ixz, Iyz
     - d'une image 3D Izz
     - de 5 images 2D => Ixx, Iyy, Ixy, Ixz, Iyz
  */
  allocatedBuf = (float*)malloc( (6 * dimz + 5 ) * dimxy * sizeof(float) );
  if ( allocatedBuf == NULL ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: auxiliary buffer allocation failed\n", proc );
    }
    return( -1 );
  }
  buf_X = allocatedBuf;
  buf_Y  = buf_X;   buf_Y  +=   dimz*dimxy;
  buf_Z  = buf_X;   buf_Z  += 2*dimz*dimxy;
  tmp_Z0 = buf_X;   tmp_Z0 += 3*dimz*dimxy;
  tmp_Z1 = buf_X;   tmp_Z1 += 4*dimz*dimxy;
  buf_ZZ = buf_X;   buf_ZZ += 5*dimz*dimxy;
  buf_XX = buf_X;   buf_XX += 6*dimz*dimxy;
  buf_YY = buf_X;   buf_YY += 6*dimz*dimxy +   dimxy;
  buf_XY = buf_X;   buf_XY += 6*dimz*dimxy + 2*dimxy;
  buf_XZ = buf_X;   buf_XZ += 6*dimz*dimxy + 3*dimxy;
  buf_YZ = buf_X;   buf_YZ += 6*dimz*dimxy + 4*dimxy;







  /* 3D filtering
   */
  theCoeffs[0] = theCoeffs[1] = theCoeffs[2] = scale;

  theDerivatives[0] = NODERIVATIVE;
  theDerivatives[1] = NODERIVATIVE;
  theDerivatives[2] = DERIVATIVE_0;
  
  if( _verbose_ >= 2 ) fprintf( stderr, " 3D filtering (.,.,0)\n" );
  
  if ( RecursiveFilterOnBuffer( theBuf, typeBuf,
				tmp_Z0, typeAllocatedBuf,
				theDim, border_length,
				theDerivatives, theCoeffs, theFilter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (.,.,0)\n", proc );
    free( allocatedBuf );
    return( -1 );
  }



  theDerivatives[0] = DERIVATIVE_0;
  theDerivatives[1] = DERIVATIVE_1;
  theDerivatives[2] = NODERIVATIVE;
  
  if( _verbose_ >= 2 ) fprintf( stderr, " 3D filtering (0,1,.)\n" );
  
  if ( RecursiveFilterOnBuffer( tmp_Z0, typeAllocatedBuf,
				buf_Y,  typeAllocatedBuf,
				theDim, border_length,
				theDerivatives, theCoeffs, theFilter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (0,1,0)\n", proc );
    free( allocatedBuf );
    return( -1 );
  }



  theDerivatives[0] = DERIVATIVE_1;
  theDerivatives[1] = DERIVATIVE_0;
  theDerivatives[2] = NODERIVATIVE;

  if( _verbose_ >= 2 ) fprintf( stderr, " 3D filtering (1,0,.)\n" );

  if ( RecursiveFilterOnBuffer( tmp_Z0, typeAllocatedBuf,
				buf_X,  typeAllocatedBuf,
				theDim, border_length,
				theDerivatives, theCoeffs, theFilter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (1,0,0)\n", proc );
    free( allocatedBuf );
    return( -1 );
  }



  theDerivatives[0] = NODERIVATIVE;
  theDerivatives[1] = NODERIVATIVE;
  theDerivatives[2] = DERIVATIVE_1;

  if( _verbose_ >= 2 ) fprintf( stderr, " 3D filtering (.,.,1)\n" );

  if ( RecursiveFilterOnBuffer( theBuf, typeBuf,
				tmp_Z1, typeAllocatedBuf,
				theDim, border_length,
				theDerivatives, theCoeffs, theFilter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (.,.,1)\n", proc );
    free( allocatedBuf );
    return( -1 );
  }



  theDerivatives[0] = DERIVATIVE_0;
  theDerivatives[1] = DERIVATIVE_0;
  theDerivatives[2] = NODERIVATIVE;

  if( _verbose_ >= 2 ) fprintf( stderr, " 3D filtering (0,0,.)\n" );

  if ( RecursiveFilterOnBuffer( tmp_Z1, typeAllocatedBuf,
				buf_Z,  typeAllocatedBuf,
				theDim, border_length,
				theDerivatives, theCoeffs, theFilter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (0,0,1)\n", proc );
    free( allocatedBuf );
    return( -1 );
  }



  theDerivatives[0] = DERIVATIVE_0;
  theDerivatives[1] = DERIVATIVE_0;
  theDerivatives[2] = DERIVATIVE_2;

  if( _verbose_ >= 2 ) fprintf( stderr, " 3D filtering (0,0,2)\n" );

  if ( RecursiveFilterOnBuffer( theBuf,  typeBuf,
				buf_ZZ,  typeAllocatedBuf,
				theDim, border_length,
				theDerivatives, theCoeffs, theFilter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (0,0,2)\n", proc );
    free( allocatedBuf );
    return( -1 );
  }
  


  /* for an efficient computation of the response
     pre-computation of the circle's points 
  */
     
  nb_circle_points = (int)(2.0 * 3.1415926536 * theta * scale) + 1;
  if ( nb_circle_points < 4 )     nb_circle_points = 4;
  if ( nb_circle_points % 4 > 0 ) nb_circle_points += 4 - nb_circle_points % 4;
  
  circle_points = (typeCirclePoint *)malloc( nb_circle_points*sizeof(typeCirclePoint) );
  if ( circle_points == NULL ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: circle's point allocation failed\n", proc );
    }
    free( allocatedBuf );
    return( -1 );
  }

  circle_points[0].c2 = 1.0;
  circle_points[0].c3 = 0.0;
  for ( p = 1; p < nb_circle_points; p ++ ) {
    circle_points[p].c2 = cos( p * 2.0 * 3.1415926536 / (double)nb_circle_points );
    circle_points[p].c3 = sin( p * 2.0 * 3.1415926536 / (double)nb_circle_points );
  }
  
  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "%s: build %d points with radius = %f\n", proc, nb_circle_points, theta * scale );
  }





  
  /* slice by slice computation
     of the 3D response
  */
  theDerivatives[2] = NODERIVATIVE;
  dim[0] = theDim[0];
  dim[1] = theDim[1];
  dim[2] = 1;

  for ( slice = 0; slice < theDim[2]; slice ++ ) {

    if( _verbose_ >= 2 ) fprintf( stderr, " processing slice #%3d with scale %9.6f    \r", slice, scale );

    /* 2D filtering
     */
    theDerivatives[0] = DERIVATIVE_2;
    theDerivatives[1] = DERIVATIVE_0;
    
    if ( RecursiveFilterOnBuffer( (void*)(&(tmp_Z0[slice*dimxy])), typeAllocatedBuf,
				  buf_XX,  typeAllocatedBuf,
				  dim, border_length,
				  theDerivatives, theCoeffs, theFilter ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: error in computating derivatives (2,0,0) in slice %d\n", 
		 proc, slice );
      free( circle_points );
      free( allocatedBuf );
      return( -1 );
    }
    
  
    theDerivatives[0] = DERIVATIVE_0;
    theDerivatives[1] = DERIVATIVE_2;
    
    if ( RecursiveFilterOnBuffer( (void*)(&(tmp_Z0[slice*dimxy])), typeAllocatedBuf,
				  buf_YY,  typeAllocatedBuf,
				  dim, border_length,
				  theDerivatives, theCoeffs, theFilter ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: error in computating derivatives (0,2,0) in slice %d\n", 
		 proc, slice );
      free( circle_points );
      free( allocatedBuf );
      return( -1 );
    }
    
    
    theDerivatives[0] = DERIVATIVE_1;
    theDerivatives[1] = DERIVATIVE_1;
    
    if ( RecursiveFilterOnBuffer( (void*)(&(tmp_Z0[slice*dimxy])), typeAllocatedBuf,
				  buf_XY,  typeAllocatedBuf,
				  dim, border_length,
				  theDerivatives, theCoeffs, theFilter ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: error in computating derivatives (1,1,0) in slice %d\n", 
		 proc, slice );
      free( circle_points );
      free( allocatedBuf );
      return( -1 );
    }
    
    
    theDerivatives[0] = DERIVATIVE_1;
    theDerivatives[1] = DERIVATIVE_0;
    
    if ( RecursiveFilterOnBuffer( (void*)(&(tmp_Z1[slice*dimxy])), typeAllocatedBuf,
				  buf_XZ,  typeAllocatedBuf,
				  dim, border_length,
				  theDerivatives, theCoeffs, theFilter ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: error in computating derivatives (1,0,1) in slice %d\n", 
		 proc, slice );
      free( circle_points );
      free( allocatedBuf );
      return( -1 );
    }
    
    
    theDerivatives[0] = DERIVATIVE_0;
    theDerivatives[1] = DERIVATIVE_1;
    
    if ( RecursiveFilterOnBuffer( (void*)(&(tmp_Z1[slice*dimxy])), typeAllocatedBuf,
				  buf_YZ,  typeAllocatedBuf,
				  dim, border_length,
				  theDerivatives, theCoeffs, theFilter ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: error in computating derivatives (0,1,1) in slice %d\n", 
		 proc, slice );
      free( circle_points );
      free( allocatedBuf );
      return( -1 );
    }
    




    /* Response computation
     */
    compute_3D_line_response_in_one_slice( buf_XX, buf_YY, &(buf_ZZ[slice*dimxy]), 
					   buf_XY, buf_XZ, buf_YZ,
					   buf_X,  buf_Y,  buf_Z,
					   circle_points, nb_circle_points, theta*scale,
					   dimx, dimy, dimz, slice,
					   &(bufResponse[slice*dimxy]), 
					   &(bufTheta[slice*dimxy]), 
					   &(bufPhi[slice*dimxy]) );


    /* for multi-scale purpose and comparison between 
       scale, the response has to be multiplicated
       by the scale.
       More precisely, the nth derivatives have to be 
       multiplicated by sigma^n
    */

    tmpResponse = &(bufResponse[slice*dimxy]);
    for ( p=0; p<dimxy; p++ ) tmpResponse[p] *= scale;

    


  }
  if( _verbose_ >= 2 ) fprintf( stderr, "\n" );








  free( circle_points );
  free( allocatedBuf );
  return( 1 );
}




























/* response computation
   all the derivatives being computed

   2nd derivatives (args #1 -> #6) are 2D buffers
   1st derivatives (args #7 -> #9) are 3D buffers
   output buffers (args #17 -> args #19)
   
*/
static void compute_3D_line_response_in_one_slice( float *theXX, float *theYY, 
						   float *theZZ,
						   float *theXY, float *theXZ, 
						   float *theYZ, 
						   float *theX,  float *theY,  
						   float *theZ,
						   typeCirclePoint *thePts,  
						   int nbPts, double radius,
						   int dimx, int dimy, int dimz, 
						   int slice,
						   float *theRep, float *theTheta,
						   float *thePhi )
{
  int x, y, i;

  double hessien[9];
  double valprop[3];
  double vecprop[9];
  double v[3];


  int n;
  double vx, vy, vz;
  double gx, gy, gz;
  double rx, ry, rz;
  double dx, dy, dz;
  int ix, iy, iz;

  double coef[8];
  int j, indx[8];


  double dxdy, dxdz, dydz, dxdydz;
  double v4, v5, v6;

  int nbPosPts;
  int nbValPts;
  double sumValPts;
  double sumPosPts;

  double theta, phi;

  
  for ( i = 0, y = 0; y < dimy; y ++ )
  for ( x = 0; x < dimx; x ++, i ++ ) {
    


    /* calcul des valeurs et vecteurs propres
     */
    hessien[0] = theXX[i];
    hessien[1] = hessien[3] = theXY[i];
    hessien[2] = hessien[6] = theXZ[i];
    hessien[4] = theYY[i];
    hessien[5] = hessien[7] = theYZ[i];
    hessien[8] = theZZ[i];


    
    theRep[i] = theTheta[i] = thePhi[i] = 0;



    if ( _ComputeEigensOfSymetricSquareMatrix( hessien, valprop, vecprop, 3 ) != 1 ) {
      continue;
    }
    if ( _SortEigensInAbsIncreasingOrder( valprop, vecprop, 3 ) != 1 ) {
      continue;
    }




    
    /* le point appartient a un tube blanc sur fond
       noir si les deux plus grandes sont negatives
       et beaucoup plus grandes que la plus petites
       
       On a fabs(valprop[0]) <= fabs(valprop[1]) <= fabs(valprop[2])
       
       On veut :
       1.        valprop[1] < 0         &&   valprop[2] < 0
       2.        fabs(valprop[0]) * 2   <=   fabs(valprop[1])

    */
    if ( valprop[1] >= 0 || valprop[2] >= 0 ) {
      continue;
    }
    if ( fabs(valprop[0]) * 2 > fabs(valprop[1]) )  {
      continue;
    }






    /* OK, ici on regarde le point 
       le vecteur de la direction du vaisseau est celui
       associe a la plus petite valeur propre
       ie vecprop[0,3,6]
       les deux autres sont vecprop[1,4,7] et vecprop[2,5,8]
    */

    nbPosPts = nbValPts = 0;
    sumValPts = sumPosPts = 0;

    for ( n=0; n<nbPts; n++ ) {

      thePts[n].f = 0;
      thePts[n].v = 0;

      vx = thePts[n].c2 * vecprop[1] + thePts[n].c3 * vecprop[2];
      rx = x + radius * vx;
      ix = (int)rx;
      if ( rx <= 0.0 || ix >= dimx-1 ) {
	continue;
      }

      vy = thePts[n].c2 * vecprop[4] + thePts[n].c3 * vecprop[5];
      ry = y + radius * vy;
      iy = (int)ry;
      if ( ry <= 0.0 || iy >= dimy-1 ) {
	continue;
      }

      vz = thePts[n].c2 * vecprop[7] + thePts[n].c3 * vecprop[8];
      rz = slice + radius * vz;
      iz = (int)rz;
      if ( rz <= 0.0 || iz >= dimz-1 ) {
	continue;
      }



      dx = rx - ix;
      dy = ry - iy;
      dz = rz - iz;

      dxdy = dx*dy;
      dxdz = dx*dz;
      dydz = dy*dz;
      dxdydz = dxdy*dz;

      v6 = dxdz-dxdydz;
      v5 = dxdy-dxdydz;
      v4 = dx-dxdy-v6;

      coef[0] = dxdydz;              indx[0] = (iz+1)*dimx*dimy + (iy+1)*dimx + ix+1;
      coef[1] = (dydz-dxdydz);       indx[1] = indx[0] - 1;
      coef[2] = v6;                  indx[2] = indx[0] - dimx;
      coef[3] = (dz-dydz-v6);        indx[3] = indx[2] - 1;
      coef[4] = v5;                  indx[4] = indx[0] - dimx*dimy;
      coef[5] = (dy-dydz-v5);        indx[5] = indx[4] - 1;
      coef[6] = v4;                  indx[6] = indx[4] - dimx;
      coef[7] = (1-dy-dz+dydz-v4);   indx[7] = indx[6] - 1;

      /*
      gx = 0;
      gx += dxdydz        * theX[iz+1][iy+1][ix+1];
      gx += (dydz-dxdydz) * theX[iz+1][iy+1][ix  ];
      gx += v6            * theX[iz+1][iy  ][ix+1];
      gx += (dz-dydz-v6)  * theX[iz+1][iy  ][ix  ];
      gx += v5            * theX[iz  ][iy+1][ix+1];
      gx += (dy-dydz-v5)  * theX[iz  ][iy+1][ix  ];
      gx += v4            * theX[iz  ][iy  ][ix+1];
      gx += (1-dy-dz+dydz-v4) * theX[iz  ][iy][ix];

      gy = 0;
      gy += dxdydz        * theY[iz+1][iy+1][ix+1];
      gy += (dydz-dxdydz) * theY[iz+1][iy+1][ix  ];
      gy += v6            * theY[iz+1][iy  ][ix+1];
      gy += (dz-dydz-v6)  * theY[iz+1][iy  ][ix  ];
      gy += v5            * theY[iz  ][iy+1][ix+1];
      gy += (dy-dydz-v5)  * theY[iz  ][iy+1][ix  ];
      gy += v4            * theY[iz  ][iy  ][ix+1];
      gy += (1-dy-dz+dydz-v4) * theY[iz  ][iy][ix];

      gz = 0;
      gz += dxdydz        * theZ[iz+1][iy+1][ix+1];
      gz += (dydz-dxdydz) * theZ[iz+1][iy+1][ix  ];
      gz += v6            * theZ[iz+1][iy  ][ix+1];
      gz += (dz-dydz-v6)  * theZ[iz+1][iy  ][ix  ];
      gz += v5            * theZ[iz  ][iy+1][ix+1];
      gz += (dy-dydz-v5)  * theZ[iz  ][iy+1][ix  ];
      gz += v4            * theZ[iz  ][iy  ][ix+1];
      gz += (1-dy-dz+dydz-v4) * theZ[iz  ][iy][ix];
      */

      gx = gy = gz = 0;
      for ( j=0; j<8; j++ ) {
	gx += coef[j] * theX[ indx[j] ];
	gy += coef[j] * theY[ indx[j] ];
	gz += coef[j] * theZ[ indx[j] ];
      }

      thePts[n].f = 1;
      thePts[n].v = - (gx * vx + gy * vy + gz * vz);

      nbValPts ++;
      sumValPts += thePts[n].v;

      if ( thePts[n].v > 0 ) {
	nbPosPts ++;
	sumPosPts += thePts[n].v;
      }
			       
    }

    
    


    /* on a la liste des points
       on pourrait calculer la reponse que l'on veut
    */
    
    if ( sumValPts <= 0.0 ) continue;

    v[0] = vecprop[0];
    v[1] = vecprop[3];
    v[2] = vecprop[6];
    
    /* comme reponse on prend la moyenne
     */

    theRep[i] = sumValPts / (double)nbValPts;
    unit_vector_to_spherical_angles( v, &theta, &phi );
    theTheta[i] = theta;
    thePhi[i]   = phi;

  }


}






int compute_3D_line_response_extrema( float *theExtrema, 
				      float *theResponse,
				      float *theTheta,
				      float *thePhi,
				      int *theDim )
{
  char *proc = "compute_3D_line_response_extrema";

  int dimx = theDim[0];
  int dimy = theDim[1];
  int dimz = theDim[2];
  int i, x, y, z;

  double v1[3], v2[3], v3[3];
  
  double rep, r;

  double theta, phi;

  double vx=0.0, vy=0.0, vz=0.0;
  
  double c = sqrt(2.0)/2.0;

  double rx, ry, rz;
  double dx, dy, dz;
  int ix, iy, iz;
  int index;
  int dir, is_a_extrema;
  double dxdy, dxdz, dydz, dxdydz;
  double v4, v5, v6;



  /* test of args
   */
  if ( theExtrema == NULL || theResponse == NULL || 
       theTheta == NULL || thePhi == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: NULL buffers in args\n", proc );
    return( -1 );
  }
  if ( theExtrema == theResponse ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: extrema and response buffers must be different\n", proc );
    return( -1 );
  }

  for ( i=0, z=0; z<dimz; z++ ) 
  for ( y=0; y<dimy; y++ ) 
  for ( x=0; x<dimx; x++, i++ ) {

    theta = (double)theTheta[i];
    phi   = (double)thePhi[i];
    
    theExtrema[i] = 0.0;
    if ( theResponse[i] <= 0.0 ) continue;

    rep = theResponse[i];

    /* v1 est le vecteur directeur du vaisseau
       v2 et v3 sont deux vecteurs orthogonaux
    */
    spherical_angles_to_units_vectors( theta, phi, v1, v2, v3 );
    
    for ( is_a_extrema=1, dir=0; dir<8 && is_a_extrema; dir++ ) {
      
      switch( dir ) {
      default :
      case 0 :
	vx = v2[0];    vy = v2[1];    vz = v2[2];    break;
      case 1 :
	vx = -v2[0];   vy = -v2[1];   vz = -v2[2];   break;
      case 2 :
	vx = v3[0];    vy = v3[1];    vz = v3[2];    break;
      case 3 :
	vx = -v3[0];   vy = -v3[1];   vz = -v3[2];   break;
      case 4 :
	vx =  c * v2[0] + c * v3[0];
	vy =  c * v2[1] + c * v3[1];
	vz =  c * v2[2] + c * v3[2];
	break;
      case 5 :
	vx = -c * v2[0] + c * v3[0];
	vy = -c * v2[1] + c * v3[1];
	vz = -c * v2[2] + c * v3[2];
	break;
       case 6 :
	vx =  c * v2[0] - c * v3[0];
	vy =  c * v2[1] - c * v3[1];
	vz =  c * v2[2] - c * v3[2];
	break;
      case 7 :
	vx = -c * v2[0] - c * v3[0];
	vy = -c * v2[1] - c * v3[1];
	vz = -c * v2[2] - c * v3[2];
	break;
      }

      rx = x + vx;
      ix = (int)rx;
      if ( rx <= 0.0 || ix >= dimx-1 ) {
	is_a_extrema = 0;
	continue;
      }
      ry = y + vy;
      iy = (int)ry;
      if ( ry <= 0.0 || iy >= dimy-1 ) {
	is_a_extrema = 0;
	continue;
      }
      rz = z + vz;
      iz = (int)rz;
      if ( rz <= 0.0 || iz >= dimz-1 ) {
	is_a_extrema = 0;
	continue;
      }


      dx = rx - ix;
      dy = ry - iy;
      dz = rz - iz;

      dxdy = dx*dy;
      dxdz = dx*dz;
      dydz = dy*dz;
      dxdydz = dxdy*dz;

      v6 = dxdz-dxdydz;
      v5 = dxdy-dxdydz;
      v4 = dx-dxdy-v6;

      index = (iz+1)*dimx*dimy + (iy+1)*dimx + (ix+1);
      r = 0;
      /*
      r += dxdydz        * theResponse[iz+1][iy+1][ix+1];
      r += (dydz-dxdydz) * theResponse[iz+1][iy+1][ix  ];
      r += v6            * theResponse[iz+1][iy  ][ix+1];
      r += (dz-dydz-v6)  * theResponse[iz+1][iy  ][ix  ];
      r += v5            * theResponse[iz  ][iy+1][ix+1];
      r += (dy-dydz-v5)  * theResponse[iz  ][iy+1][ix  ];
      r += v4            * theResponse[iz  ][iy  ][ix+1];
      r += (1-dy-dz+dydz-v4) * theResponse[iz  ][iy][ix];
      */
      r += dxdydz        * theResponse[index];
      r += (dydz-dxdydz) * theResponse[index-1];
      r += v6            * theResponse[index-dimx];
      r += (dz-dydz-v6)  * theResponse[index-dimx-1];
      index -= dimx*dimy;
      r += v5            * theResponse[index];
      r += (dy-dydz-v5)  * theResponse[index-1];
      r += v4            * theResponse[index-dimx];
      r += (1-dy-dz+dydz-v4) * theResponse[index-dimx-1];
      
      if ( r >= rep ) {
	is_a_extrema = 0;
	continue;
      }
      
    }
    
    if ( is_a_extrema ) {
      theExtrema[i] = rep;
    }

  }

  return( 1 );
}
























/* Transformation du vecteur unitaire (x,y,z) en 2 angles (theta,phi)

   ( x )   ( cos theta sin phi )
   ( y ) = ( sin theta sin phi )
   ( z )   ( cos phi )
*/
static void unit_vector_to_spherical_angles( const double *v, 
				  double *theta,
				  double *phi )
{
  double n;
  
  /* from % `man acos`
     
    ...
    The acos() and acosf() functions compute the principal value of the arc
    cosine of x in the interval [0,pi] radians. The value of x must be in the
    domain [-1,1].
    ...
    The atan2() and atan2f() functions compute the principal value of the arc
    tangent of y/x, in the interval [-pi,pi] radians.  The sign of atan2() and
    atan2f() is determined by the sign of y.  The value of atan2(y,x) is com-
    puted as follows where f is the number of fraction bits associated with the
    data type.
    _____________________________________________________
    Value of Input Arguments    Angle Returned
    _____________________________________________________
    x = 0 or y/x > 2**(f+1)     pi/2 * (sign y)
    x > 0 and y/x <= 2**(f+1)   atan(y/x)
    x < 0 and y/x <= 2**(f+1)   pi * (sign y) + atan(y/x)
    _____________________________________________________
    ...
  */

  n = sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
  *phi   = acos( v[2]/n );
  *theta = atan2( v[1]/n, v[0]/n );
}



/* Transformation de 2 angles (theta,phi) en repere orthonorme

   le vecteur lie au 2 angles est 
   ( x )   ( cos theta sin phi )
   ( y ) = ( sin theta sin phi )
   ( z )   ( cos phi )

   un vecteur orthonorme est
   (  sin theta )
   ( -cos theta )
   ( 0 )
   
   et le troisieme vecteur (par produit vectoriel) est alors
   ( cos theta cos phi )
   ( sin theta cos phi )
   ( -sin phi )

   Le repere est donc direct.

*/

static void spherical_angles_to_units_vectors( const double theta,
				    const double phi,
				    double *v1, 
				    double *v2,
				    double *v3 )
{
  double sp, cp, ct, st;

  sp = sin( phi );
  cp = cos( phi );
  st = sin( theta );
  ct = cos( theta );
  
  v1[0] = sp * ct;
  v1[1] = sp * st;
  v1[2] = cp;
  
  v2[0] =  st;
  v2[1] = -ct;
  v2[2] = 0.0;
  
  v3[0] =  cp * ct;
  v3[1] =  cp * st;
  v3[2] = -sp;
}





































typedef enum {
  _WILLBEDELETED   = 220,
  _DELETABLE       = 230,
  _UNDELETABLE     = 255
} status_type_point_to_be_thinned;


typedef struct {
  int x;
  int y;
  int z;
  int i;
  status_type_point_to_be_thinned status;
  int is_inside;
} type_point_to_be_thinned;

/* construction d'une liste de points pour l'amincissement
 */

static type_point_to_be_thinned *build_list_of_points_to_be_thinned( unsigned char *input_buffer,
								     unsigned char *output_buffer,
								     int *theDim,
								     int *nb_of_points )
{
  char *proc = "build_list_of_points_to_be_thinned";
  int dimx = theDim[0];
  int dimy = theDim[1];
  int dimz = theDim[2];
  int x, y, z, i;
  int iz, iy;
  int nb_allocated_points = 0;
  int nb_points = 0;
  int nb_stack_of_points = dimx * dimy;
  type_point_to_be_thinned *theList = NULL, *tmp = NULL;


  *nb_of_points = -1;

  for ( i=0, z=0; z<dimz; z++ ) {
    iz =  ( z == 0 || z == dimz-1 ) ? 0 : 1;
    for ( y=0; y<dimy; y++ ) {
      iy = ( iz == 0 || y == 0 || y == dimy-1 ) ? 0 : 1;
      for ( x=0; x<dimx; x++, i++ ) {

	/* background
	 */
	if ( input_buffer[i] == 0 ) {
	  output_buffer[i] = 0;
	  continue;
	}

	/* the point is to be added to the list
	 */
	output_buffer[i] = _DELETABLE;

	if ( nb_points >= nb_allocated_points ) {
	  tmp = (type_point_to_be_thinned*)malloc( (nb_allocated_points+nb_stack_of_points)*sizeof(type_point_to_be_thinned) );
	  if ( tmp == NULL ) {
	    if ( theList != NULL ) free( theList );
	    if ( _verbose_ ) {
	      fprintf( stderr, "%s: allocation #%d failed\n", 
		       proc, 1+nb_allocated_points/nb_stack_of_points );
	    }
	    return( NULL );
	  }
	  if ( theList != NULL ) {
	    (void)memcpy( tmp, theList, (nb_allocated_points)*sizeof(type_point_to_be_thinned) );
	    free( theList );
	  }
	  theList = tmp;
	  tmp     = NULL;
	  nb_allocated_points += nb_stack_of_points;
	}

	/* adding the point
	 */
	theList[ nb_points ].x = x;
	theList[ nb_points ].y = y;
	theList[ nb_points ].z = z;
	theList[ nb_points ].i = i;
	theList[ nb_points ].status = _DELETABLE;
	theList[ nb_points ].is_inside = ( iy == 0 || x == 0 || x == dimx-1 ) ? 0 : 1;
	nb_points ++;
      }
    }
  }

  *nb_of_points = nb_points;
  return( theList );
}



















/* amincissement d'une image binaire contenant des objets filaires 
   "epais"
*/

int thin_3D_thick_lines( unsigned char *input_buffer,
			 unsigned char *output_buffer,
			 int *theDim )
{
  char *proc = "thin_3D_thick_lines";

  type_point_to_be_thinned *theList = NULL;
  int nb_of_points = 0;
  int dimx = theDim[0];
  int dimy = theDim[1];
  int dimz = theDim[2];

  /*
   - Le voisinage est numerote comme suit
        0  1  2  -   9 10 11  -  18 19 20
        3  4  5  -  12 13 14  -  21 22 23
        6  7  8  -  15 16 17  -  24 25 26
   - dans une direction d'amicissement, il y a une 
     6- ou une 18-epaisseur si 
     #0 == 0 (c'est un point du fond)
     il existe #i != 0, i=1...5 
  */
  int neighborhood[27];
  int thin_along_up_dir[6] = { 22,  4,  1,  3,  5,  7 };
  int thin_along_bo_dir[6] = {  4, 22, 19, 21, 23, 25 };
  int thin_along_no_dir[6] = { 16, 10,  1,  9, 11, 19 };
  int thin_along_so_dir[6] = { 10, 16,  7, 15, 17, 25 };
  int thin_along_ea_dir[6] = { 14, 12,  3,  9, 15, 21 };
  int thin_along_we_dir[6] = { 12, 14,  5, 11, 17, 23 };
  int *thin_dir;
  int i_dir = 0;

  int offset[27];
  int outside_value = 0;

  int n, nb_points = 0;
  int nb_deleted_points = 0;
  int nb_deleted_points_in_one_dir;
  int nb_neighbors;
  
  int i, j, k;
  


  /* initialisation des offsets
   */
  offset[12] = -1;
  offset[13] =  0;
  offset[14] =  1;
  for ( i=-1; i<=1; i++ ) {
    offset[10+i] = offset[13+i] - dimx;
    offset[16+i] = offset[13+i] + dimx;
  }
  for ( i=-4; i<=4; i++ ) {
    offset[ 4+i] = offset[13+i] - dimx*dimy;
    offset[22+i] = offset[13+i] + dimx*dimy;
    
  }



  /* construction de la liste de points
   */
  theList = build_list_of_points_to_be_thinned( input_buffer, output_buffer, theDim, &nb_of_points );
  if ( theList == NULL ) {
    if ( nb_of_points == 0 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: the input buffer is empty\n", proc );
    } else {
      if ( _verbose_ )
	fprintf( stderr, "%s: error while building the list of points\n", proc );
    }
    return( -1 );
  }



  /* main loop 
   */
  nb_points = nb_of_points;

  do {
    
    nb_deleted_points = 0;

    /* lock for thicknesses in each direction
     */
    for ( i_dir = 0; i_dir < 6 && nb_points > 0; i_dir ++ ) {
      
      switch( i_dir ) {
      default :
      case 0 :   thin_dir = thin_along_up_dir;   break;
      case 1 :   thin_dir = thin_along_bo_dir;   break;
      case 2 :   thin_dir = thin_along_no_dir;   break;
      case 3 :   thin_dir = thin_along_so_dir;   break;
      case 4 :   thin_dir = thin_along_ea_dir;   break;
      case 5 :   thin_dir = thin_along_we_dir;   break;
      }

      for ( n=0; n<nb_points; n++ ) {


	if ( theList[n].status != _DELETABLE )
	  continue;

	/* get the neighborhood
	 */
	if ( theList[n].is_inside ) {
	  for ( i=0; i<27; i++ )
	    neighborhood[ i ] = output_buffer[ theList[n].i + offset[i] ];
	} else {
	  for ( k = -1; k < 2; k++ ) 
	  for ( j = -1; j < 2; j++ ) 
	  for ( i = -1; i < 2; i++ ) {
	    if ( theList[n].x+i < 0 || theList[n].x+i >= dimx ||
		 theList[n].y+j < 0 || theList[n].y+j >= dimy ||
		 theList[n].z+k < 0 || theList[n].y+k >= dimz ) {
	      neighborhood[13 + k*9 + j*3 +i] = outside_value;
	    } else {
	      neighborhood[13 + k*9 + j*3 +i] = output_buffer[theList[n].i + k*dimx*dimy + j*dimx + i];
	    }
	  }
	}


	/* look for directional thickness
	   we want 
	   - neighborhood[ thin_dir[0] ]  = 0 : background
	   - neighborhood[ thin_dir[i] ] != 0 for at least one i in 1..5
	*/
	if ( neighborhood[ thin_dir[0] ] != 0 ) continue;
	if ( neighborhood[ thin_dir[1] ] == 0 && 
	     neighborhood[ thin_dir[2] ] == 0 && 
	     neighborhood[ thin_dir[3] ] == 0 && 
	     neighborhood[ thin_dir[4] ] == 0 && 
	     neighborhood[ thin_dir[5] ] == 0 ) continue;


	/* if the point is not simple
	   process the next point
	*/
	if ( IsA3DPointSimple( neighborhood ) != 1 ) {
	  continue;
	}

	/* count the points in the neighborhood
	   do not count the central (ie current) point
	 */
	for ( nb_neighbors=0, i=0; i<27; i++ )
	  if ( neighborhood[ i ] )  nb_neighbors ++;
	nb_neighbors --;

	/* in the neighborhood, some points are already to be deleted
	 */
	for ( i=0; i<27; i++ )
	  if ( neighborhood[ i ] == _WILLBEDELETED )
	    neighborhood[ i ] = 0;

	/* test again
	 */
	if ( IsA3DPointSimple( neighborhood ) != 1 ) {
	  continue;
	}

	/* ok, the point can be deleted
	   just check the end condition
	   if nb_neighbors == 1, it is an end point
	*/
	if ( nb_neighbors == 1 ) {
	  theList[n].status           = _UNDELETABLE;
	  output_buffer[theList[n].i] = _UNDELETABLE;
	  continue;
	}

	theList[n].status           = _WILLBEDELETED;
	output_buffer[theList[n].i] = _WILLBEDELETED;
      }

      /* re-process the list
       */
      nb_deleted_points_in_one_dir = 0;
      for ( i=0, n=0; n<nb_points; n++ ) {
	switch( theList[n].status ) {
	case _DELETABLE :
	  /* the point is still to be processed 
	   */
	  theList[ i++ ] = theList[ n ];
	  break;
	case _WILLBEDELETED :
	   /* the point is to be removed from the image
	   */
	  output_buffer[theList[n].i] = 0;
	  nb_deleted_points_in_one_dir ++;
	case _UNDELETABLE :
	  /* the point is to be removed from the list
	     => do nothing
	   */
	  break;
	}
      }

      /* update variables
       */
      nb_points = i;
      nb_deleted_points += nb_deleted_points_in_one_dir;


    } /* loop of directions */
    
    /* re-process the list
     */
    for ( n=0; n<nb_points; n++ ) {
      output_buffer[theList[n].i] = _UNDELETABLE;
    }

  } while ( nb_deleted_points > 0 && nb_points > 0 ); /* end of main loop */







  free( theList );
  return( 1 );
}





