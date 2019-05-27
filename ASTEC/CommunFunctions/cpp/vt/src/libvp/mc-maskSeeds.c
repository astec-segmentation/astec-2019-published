/*************************************************************************
 * mc-maskSeeds.c - template for executable creation
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mer 25 jul 2018 18:20:30 CEST
 *
 * ADDITIONS, CHANGES
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include <api-maskSeeds.h>

#include <vt_common.h>





static int _verbose_ = 1;





/* static function definitions
 */

static char *_Array2Str( int argc, char *argv[] );
static char *_BaseName( char *p );
static double _GetTime();
static double _GetClock();






int main( int argc, char *argv[] )
{
  lineCmdParamMaskSeeds par;
  vt_image *seed_image;
  vt_image *cell_image;
  char *lineoptions;


  double time_init = _GetTime();
  double time_exit;
  double clock_init = _GetClock();
  double clock_exit;


  /* parameter initialization
   */
  API_InitParam_maskSeeds( &par );



  /* parameter parsing
   */
  if ( argc <= 1 )
      API_ErrorParse_maskSeeds( _BaseName( argv[0] ), (char*)NULL, 0 );
  API_ParseParam_maskSeeds( 1, argc, argv, &par );
  


  /* input image reading
   */
  seed_image = _VT_Inrimage( par.input_seed_name );
  if ( seed_image == (vt_image*)NULL ) {
      API_ErrorParse_maskSeeds( _BaseName( argv[0] ), "unable to read input seed image ...\n", 0 );
  }

  cell_image = _VT_Inrimage( par.input_cell_name );
  if ( cell_image == (vt_image*)NULL ) {
      VT_FreeImage( seed_image );
      VT_Free( (void**)&seed_image );
      API_ErrorParse_maskSeeds( _BaseName( argv[0] ), "unable to read input cell image ...\n", 0 );
  }


  /* API call
   */

  lineoptions = _Array2Str( argc, argv );
  if ( lineoptions == (char*)NULL ) {
      VT_FreeImage( cell_image );
      VT_Free( (void**)&cell_image );
      VT_FreeImage( seed_image );
      VT_Free( (void**)&seed_image );
      API_ErrorParse_maskSeeds( _BaseName( argv[0] ), "unable to translate command line options ...\n", 0 );
  }

  if ( API_maskSeeds( seed_image, cell_image, seed_image, lineoptions, (char*)NULL ) != 1 ) {
      free( lineoptions );
      VT_FreeImage( cell_image );
      VT_Free( (void**)&cell_image );
      VT_FreeImage( seed_image );
      VT_Free( (void**)&seed_image );
      API_ErrorParse_maskSeeds( _BaseName( argv[0] ), "some error occurs during processing ...\n", 0 );
  }



  /* memory freeing
   */
  free( lineoptions );
  VT_FreeImage( cell_image );
  VT_Free( (void**)&cell_image );



  /* output image writing
   */

  if ( VT_WriteInrimageWithName( seed_image, par.output_name ) == -1 ) {
      VT_FreeImage( seed_image );
      VT_Free( (void**)&seed_image );
    API_ErrorParse_maskSeeds( _BaseName( argv[0] ), "unable to write output image ...\n", 0 );
  }



  /* memory freeing
   */
  VT_FreeImage( seed_image );
  VT_Free( (void**)&seed_image );




  time_exit = _GetTime();
  clock_exit = _GetClock();

  if ( par.print_time ) { 
    fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName( argv[0] ), time_exit - time_init );
    fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
    fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );
  }


  return( 0 );
}





/************************************************************
 *
 * static functions
 *
 ************************************************************/



static char *_Array2Str( int argc, char *argv[] )
{
  char *proc = "_Array2Str";
  int i, l;
  char *s, *t;

  if ( argc <= 1 || argv == (char**)NULL ) {
    if ( _verbose_ >= 2 )
      fprintf( stderr, "%s: no options in argv[]\n", proc );
    return( (char*)NULL );
  }

  /* there are argc-1 strings
   * compute the sum of string lengths from 1 to argc-1
   * + number of interval between successive strings (argc-2)
   * + 1 to add a trailing '\0'
   */
  for ( l=argc-1, i=1; i<argc; i++ ) {
    l += strlen( argv[i] );
  }

  s = (char*)malloc( l * sizeof( char ) );
  if ( s == (char*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( (char*)NULL );
  }

  for ( t=s, i=1; i<argc; i++ ) {
    (void)strncpy( t, argv[i], strlen( argv[i] ) );
    t += strlen( argv[i] );
    if ( i < argc-1 ) {
      *t = ' ';
      t++;
    }
    else {
      *t = '\0';
    }
  }

  return( s );
}



static char *_BaseName( char *p )
{
  int l;
  if ( p == (char*)NULL ) return( (char*)NULL );
  l = strlen( p ) - 1;
  while ( l >= 0 && p[l] != '/' ) l--;
  if ( l < 0 ) l = 0;
  if ( p[l] == '/' ) l++;
  return( &(p[l]) );
}



static double _GetTime()
{
  struct timeval tv;
  gettimeofday(&tv, (void *)0);
  return ( (double) tv.tv_sec + tv.tv_usec*1e-6 );
}



static double _GetClock() 
{
  return ( (double) clock() / (double)CLOCKS_PER_SEC );
}
