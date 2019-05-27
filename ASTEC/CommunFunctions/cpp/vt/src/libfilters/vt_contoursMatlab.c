
#include <stdlib.h>
#include <unistd.h>

#include <vt_contoursMatlab.h>


void MAT_DrawContour2D( typeContour2D *c,
			int rawDataFileDesc,
			FILE *commandFile,
			int d )
{
  char *proc = "MAT_DrawContour2D";
  int i;
  double *t;

  t = (double*)malloc( c->n * sizeof(double) );

  /* on ajoute 1 car l'origine est a (1,1) en matlab
   */

  for ( i=0; i<c->n; i++ ) t[i] = c->thePts[i].x + 1;
  if ( write( rawDataFileDesc, t, c->n * sizeof(double) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }
  for ( i=0; i<c->n; i++ ) t[i] = c->thePts[i].y + 1;
  if ( write( rawDataFileDesc, t, c->n * sizeof(double) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }

  free( t );

  fprintf( commandFile, "XC%d = fread( fid, %d, 'float%lu' );\n", 
	   d, c->n, 8*sizeof(double) );
  fprintf( commandFile, "YC%d = fread( fid, %d, 'float%lu' );\n", 
	   d, c->n, 8*sizeof(double) );
  fprintf( commandFile, "plot( XC%d, YC%d );\n", d, d );
  if ( c->topology == _CLOSED_ )
    fprintf( commandFile, "plot( [XC%d(%d) XC%d(1)], [YC%d(%d) YC%d(1)] );\n",
	     d, c->n, d, d, c->n, d );
}




void MAT_DrawContour2Din3D( typeContour2D *c,
			    double z,
			    double *mat,
			    int rawDataFileDesc,
			    FILE *commandFile,
			    int d, char *options )
{
  char * proc = "MAT_DrawContour2Din3D";
  int i;
  double *t;

  t = (double*)malloc( 3 * c->n * sizeof(double) );

  if ( mat != NULL ) {
    for ( i=0; i<c->n; i++ ) {
      t[       i] = mat[0]*c->thePts[i].x + mat[1]*c->thePts[i].y + mat[ 2]*z + mat[ 3];
      t[  c->n+i] = mat[4]*c->thePts[i].x + mat[5]*c->thePts[i].y + mat[ 6]*z + mat[ 7];
      t[2*c->n+i] = mat[8]*c->thePts[i].x + mat[9]*c->thePts[i].y + mat[10]*z + mat[11];
    }
  }

  else {

    for ( i=0; i<c->n; i++ ) {
      t[       i] = c->thePts[i].x;
      t[  c->n+i] = c->thePts[i].y;
      t[2*c->n+i] = z;
     }

    /* on ajoute 1 car l'origine est a (1,1) en matlab
     */
    for ( i=0; i<c->n; i++ ) {
      t[       i] += 1;
      t[  c->n+i] += 1;
      /* t[2*c->n+i] += 1; */
    }

  }


  if ( write( rawDataFileDesc, t, 3 * c->n * sizeof(double) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }
  
  free( t );

  fprintf( commandFile, "XC%d = fread( fid, %d, 'float%lu' );\n", 
	   d, c->n, 8*sizeof(double) );
  fprintf( commandFile, "YC%d = fread( fid, %d, 'float%lu' );\n", 
	   d, c->n, 8*sizeof(double) );
  fprintf( commandFile, "ZC%d = fread( fid, %d, 'float%lu' );\n", 
	   d, c->n, 8*sizeof(double) );
  if ( options == NULL )
    fprintf( commandFile, "plot3( XC%d, YC%d, ZC%d );\n", d, d, d );
  else 
    fprintf( commandFile, "plot3( XC%d, YC%d, ZC%d, %s );\n", d, d, d, options );
  if ( c->topology == _CLOSED_ ) {
    if ( options == NULL )
      fprintf( commandFile, "plot3( [XC%d(%d) XC%d(1)], [YC%d(%d) YC%d(1)], [ZC%d(%d) ZC%d(1)] );\n",
	       d, c->n, d, d, c->n, d, d, c->n, d );
    else
      fprintf( commandFile, "plot3( [XC%d(%d) XC%d(1)], [YC%d(%d) YC%d(1)], [ZC%d(%d) ZC%d(1)], %s );\n",
	       d, c->n, d, d, c->n, d, d, c->n, d, options );
  }
}



void MAT_DrawContour3D( typeContour3D *c,
			double *mat,
			int rawDataFileDesc,
			FILE *commandFile,
			int d )
{
  char *proc = "MAT_DrawContour3D";
  int i;
  double *t;

  t = (double*)malloc( 3 * c->n * sizeof(double) );

  if ( mat != NULL ) {
    for ( i=0; i<c->n; i++ ) {
      t[       i] = mat[0]*c->thePts[i].x + mat[1]*c->thePts[i].y + mat[ 2]*c->thePts[i].z + mat[ 3];
      t[  c->n+i] = mat[4]*c->thePts[i].x + mat[5]*c->thePts[i].y + mat[ 6]*c->thePts[i].z + mat[ 7];
      t[2*c->n+i] = mat[8]*c->thePts[i].x + mat[9]*c->thePts[i].y + mat[10]*c->thePts[i].z + mat[11];
    }
  }

  else {

    for ( i=0; i<c->n; i++ ) {
      t[       i] = c->thePts[i].x;
      t[  c->n+i] = c->thePts[i].y;
      t[2*c->n+i] = c->thePts[i].z;
     }

    /* on ajoute 1 car l'origine est a (1,1) en matlab
     */
    for ( i=0; i<c->n; i++ ) {
      t[       i] += 1;
      t[  c->n+i] += 1;
      /* t[2*c->n+i] += 1; */
    }

  }


  if ( write( rawDataFileDesc, t, 3 * c->n * sizeof(double) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }
  
  free( t );

  fprintf( commandFile, "XC%d = fread( fid, %d, 'float%lu' );\n", 
	   d, c->n, 8*sizeof(double) );
  fprintf( commandFile, "YC%d = fread( fid, %d, 'float%lu' );\n", 
	   d, c->n, 8*sizeof(double) );
  fprintf( commandFile, "ZC%d = fread( fid, %d, 'float%lu' );\n", 
	   d, c->n, 8*sizeof(double) );
  fprintf( commandFile, "plot3( XC%d, YC%d, ZC%d );\n", d, d, d );
  if ( c->topology == _CLOSED_ )
    fprintf( commandFile, "plot3( [XC%d(%d) XC%d(1)], [YC%d(%d) YC%d(1)], [ZC%d(%d) ZC%d(1)] );\n",
	     d, c->n, d, d, c->n, d, d, c->n, d );
}
