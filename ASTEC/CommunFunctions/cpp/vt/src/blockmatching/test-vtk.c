/*************************************************************************
 * test-copy.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2012, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mon Nov 19 17:45:00 CET 2012
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#include <stdio.h>
#include <stdlib.h>

#include <bal-vtk.h>



int main(int argc, char *argv[])
{
  bal_vtkDataSet mesh;
  union {
    unsigned char uc[2];
    unsigned short us;
  } twobytes;
  twobytes.us = 255;

  if ( twobytes.uc[1] == 0 ) {
    fprintf( stderr, "LITTLE\n" );
  }
  else {
    fprintf( stderr, "BIG\n" );
  }

  if ( argc <= 2 ) {
    fprintf( stderr, "must provide a file name and a destination\n" );
    exit( 1 );
  }

  BAL_InitVtkDataSet( &mesh );

  if ( BAL_ReadVtkDataSet( &mesh, argv[1] ) != 1 ) {
    BAL_FreeVtkDataSet( &mesh );
    fprintf( stderr, "error when reading\n" );
    return( -1 );
  }

  BAL_FprintfVtkDataSet( stderr, &mesh );

  if ( BAL_WriteVtkDataSet( &mesh, argv[2], _VTK_ASCII_DATA_ ) != 1 ) {
    BAL_FprintfVtkDataSet( stderr, &mesh );
    BAL_FreeVtkDataSet( &mesh );
    fprintf( stderr, "error when write\n" );
    return( -1 );
  }



  BAL_FreeVtkDataSet( &mesh );

  return( 0 );
}
