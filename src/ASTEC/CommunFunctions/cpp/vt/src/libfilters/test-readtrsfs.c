#include <stdio.h>
#include <vt_tubeutils.h>

int main( int argc, char *argv[] )
{
  int nb = 0;
  int i, j;
  double **mat;

  mat = _ReadTransformationsFile( argv[1], &nb );
  for ( i=0; i<nb; i++ ) {
    for (j=0; j<16;j++ )
      printf( " %g", mat[i][j] );
    printf("\n" );
  }
}
