#include <vt_common.h>
#include <vt_isocontours.h>

#include <vt_tube2Dmatlab.h>
#include <vt_contoursMatlab.h>
#include <vt_tubeutils.h>

int main( int argc, char *argv[] )
{
  vt_image image;
  double xc, yc;
  int i,j;
  float ***buf;
  
  double r, s;
  typeSlice *slice;

  FILE *f, *fopen();
  int fd;

  VT_Image( &image );
  VT_InitImage( &image, "test_radius.inr", 512, 512, 1, FLOAT );
  VT_AllocImage( &image );

  xc = 0.5 + (image.dim.x -1)/ 2;
  yc = 0.5 + (image.dim.y -1)/ 2;

  printf( "\n xc =%f yc =%f \n\n", xc, yc );

  buf = (float***)image.array;
  for (i=0; i<image.dim.x; i++)
  for (j=0; j<image.dim.y; j++)
    buf[0][j][i] = sqrt( (j-yc)*(j-yc) + (i-xc)*(i-xc) );

  fd = creat( "test_radius.raw", S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
  f = fopen( "test_radius.m", "w" );
  fprintf( f, "echo off\n" );
  fprintf( f, "fid = fopen('test_radius.raw', 'r' );\n" );

  fprintf( f, "\n\n\n" );
  fprintf( f, "t=0:0.1:2*pi+0.1;\n" );
  fprintf( f, "figure;\n" );
  fprintf( f, "hold on;\n" );
  
  for ( r=1.0; r < 10.0; r+= 1 ) {
    slice = VT_Compute2DIsoContours( &image, 0, r );

    MAT_DrawContour2D( slice->theContours[0], fd, f, (int)(r+0.5) );
    fprintf( f, "plot( 1+%f+%f*cos(t), 1+%f+%f*sin(t), 'r' );\n",
	     xc, r, yc, r );

    s = surfaceContour2D( slice->theContours[0], 0, 0 );
    printf( "r=%f surface=%f radius=%f\n",
	    r, s, sqrt( s/3.1416 ) );
    s = surfaceContour2D( slice->theContours[0], xc, yc );
    printf( "   r=%f surface=%f radius=%f\n",
	    r, s, sqrt( s/3.1416 ) );

    freeSlice( slice );
    free( slice );
  }

  fprintf( f, "\n\n\n" );
  fprintf( f, "fclose(fid);\n" );

  
  VT_WriteInrimage( &image );
}
