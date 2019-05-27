#include <stdio.h>
#include <vt_contours.h>

int main( int argc, char *argv[] )
{
  typeStructure structure;

  initStructure( &structure );
  readStructure( &structure, argv[1] );

  writeStructure( &structure, argv[2] );
}
