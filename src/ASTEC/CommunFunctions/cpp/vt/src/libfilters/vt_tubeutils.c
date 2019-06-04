/*************************************************************************
 * vt_tubeutils.c -
 *
 * $Id: vt_tubeutils.c,v 1.5 2002/09/05 17:15:06 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Tue Jul 11 15:46:00 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <vt_tubeutils.h>


double ** _ReadTransformationsFile( char *name, int *nb )
{
  char *proc = "_ReadTransformationsFile";

  FILE *f;
  int readIsNotFinished = 1;
  int l = 512;
  char *str = NULL;
  int i, status;
  
  double mat[16];
  int stack = 50;
  int r=0, a=0;
  double *amat = NULL;
  double *tmat = NULL;
  
  char *tmp;
  double **pmat = NULL;
  
  


  str = (char*)malloc( l );
  f = fopen( name, "r" );
  if ( f == NULL ) {
    fprintf( stderr, "%s: unable to open '%s'\n", proc, name );
    return( NULL );
  }


  do {

    if ( fgets( str, l, f ) == NULL ) {
      readIsNotFinished = 0;
      continue;
    }


    for ( i = 0; i<l 
	    && (str[i] == ' ' || str[i] == '\t'); i++ )
      ;
    if ( i == l || str[i] == '#' || str[i] == '\n' ) {
      continue;
    }
    

    if ( str[i] != '-' && (str[i] < '0' || str[i] > '9') ) continue;

    status = sscanf( &(str[i]), "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		     &mat[ 0], &mat[ 1], &mat[ 2], &mat[ 3],
		     &mat[ 4], &mat[ 5], &mat[ 6], &mat[ 7],
		     &mat[ 8], &mat[ 9], &mat[10], &mat[11],
		     &mat[12], &mat[13], &mat[14], &mat[15] );
    if ( status != 16 ) continue;
    
    if ( r <= a ) {
      tmat = (double*)malloc( (a+stack)*16*sizeof(double) );
      if ( a > 0 ) {
	memcpy( tmat, amat, a*16*sizeof(double) );
	free( amat );
      }
      amat = tmat;
      a = a+stack;
    }
    memcpy( &(amat[r*16]), mat, 16*sizeof(double) );
    r ++;
    
  } while ( readIsNotFinished == 1 );

  fclose( f );
  free( str );

  fprintf( stderr, "%s: read %d transformations in '%s'\n",
	   proc, r, name );


  tmp = (char*)malloc( r * (sizeof(double*) + 16*sizeof(double)) );
  pmat = (double**)tmp;
  tmp += r * sizeof(double*);
  memcpy( tmp, amat, r*16*sizeof(double) );
  for ( i=0; i<r; i++ )
    pmat[i] = &( ((double*)tmp)[i*16] );
  free( amat );
  
  *nb = r;

  return( pmat );

}




float _GetInterpolated2DValue( vt_image *theIm,
			       double x, double y, int z )
{
  char *proc = "_GetInterpolated2DValue";
  int ix, iy;
  double leftdx, leftdy;
  double righdx, righdy;
  double val = 0.0;

  if ( x < 0.0 || x >= theIm->dim.x-1 ) return( 0.0 );
  if ( y < 0.0 || y >= theIm->dim.y-1 ) return( 0.0 );

  ix = (int)x;
  iy = (int)y;
  leftdx = x - ix;
  leftdy = y - iy;
  righdx = 1.0 - leftdx;
  righdy = 1.0 - leftdy;
  
#define INTERPOLATE2DVALUE( TYPE ) {         \
  TYPE ***theBuf = (TYPE ***)theIm->array; \
  val  = righdx*righdy * (double)theBuf[z][iy][ix];     \
  val += leftdx*righdy * (double)theBuf[z][iy][ix+1];   \
  val += righdx*leftdy * (double)theBuf[z][iy+1][ix];   \
  val += leftdx*leftdy * (double)theBuf[z][iy+1][ix+1]; \
  }

  switch ( theIm->type ) {
  default :
    fprintf( stderr, "%s: such type not handled in switch\n", proc );
    exit( 2 );
  case UCHAR :
    INTERPOLATE2DVALUE( u8 );
    break;
  case FLOAT :
    INTERPOLATE2DVALUE( r32 );
    break;
  }
  return( val );
}

















float _GetInterpolated3DValue( vt_image *theIm,
			       double x, double y, double z )
{
  char *proc = "_GetInterpolated3DValue";
  int ix, iy, iz;
  double leftdx, leftdy, leftdz;
  double righdx, righdy, righdz;
  double val = 0.0;

  if ( x < 0.0 || x >= theIm->dim.x-1 ) return( 0.0 );
  if ( y < 0.0 || y >= theIm->dim.y-1 ) return( 0.0 );
  if ( z < 0.0 || z >= theIm->dim.z-1 ) return( 0.0 );

  ix = (int)x;
  iy = (int)y;
  iz = (int)z;
  leftdx = x - ix;
  leftdy = y - iy;
  leftdz = z - iz;
  righdx = 1.0 - leftdx;
  righdy = 1.0 - leftdy;
  righdz = 1.0 - leftdz;
  
#define INTERPOLATE3DVALUE( TYPE ) {         \
  TYPE ***theBuf = (TYPE ***)theIm->array; \
  val  = righdx*righdy*righdz * (double)theBuf[iz][iy][ix];     \
  val += leftdx*righdy*righdz * (double)theBuf[iz][iy][ix+1];   \
  val += righdx*leftdy*righdz * (double)theBuf[iz][iy+1][ix];   \
  val += leftdx*leftdy*righdz * (double)theBuf[iz][iy+1][ix+1]; \
  val += righdx*righdy*leftdz * (double)theBuf[iz+1][iy][ix];     \
  val += leftdx*righdy*leftdz * (double)theBuf[iz+1][iy][ix+1];   \
  val += righdx*leftdy*leftdz * (double)theBuf[iz+1][iy+1][ix];   \
  val += leftdx*leftdy*leftdz * (double)theBuf[iz+1][iy+1][ix+1]; \
  }

  switch ( theIm->type ) {
  default :
    fprintf( stderr, "%s: such type not handled in switch\n", proc );
    exit( 2 );
  case UCHAR :
    INTERPOLATE3DVALUE( u8 );
    break;
  case USHORT :
    INTERPOLATE3DVALUE( u16 );
    break;
  case FLOAT :
    INTERPOLATE3DVALUE( r32 );
    break;
  }
  return( val );
}

















typeConvolutionPoint *_Build2DConvolutionMask( int r,
					       double sigma,
					       int *nb )
{
  typeConvolutionPoint *tmp = (typeConvolutionPoint *)NULL;
  int size = (2*r+1) * (2*r+1);
  int i, x, y;
  int d2;
  double sum = 0.0;
  
  *nb = 0;
  if ( r <= 0 ) return( tmp );

  tmp = (typeConvolutionPoint *)malloc( size * sizeof(typeConvolutionPoint) );

  i=0;
  for ( y = -r; y <= r ; y ++ )
  for ( x = -r; x <= r ; x ++ ) {
    d2 = x*x + y*y;
    if ( d2 > r*r ) continue;
    tmp[i].dx = x;
    tmp[i].dy = y;
    tmp[i].dz = 0;
    tmp[i].c = exp( - d2 / ( 2.0 * sigma * sigma ) ) ;
    sum += tmp[i].c;
    i ++;
  }
  
  for ( x=0; x<i; x++ ) 
    tmp[x].c /= sum;
  
  *nb = i;
  return( tmp );
}
