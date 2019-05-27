/*************************************************************************
 * membrane.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2016, all rights reserved
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * Ven 2 dec 2016 13:36:43 CEST
 *
 * ADDITIONS, CHANGES
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include <vtmalloc.h>

#include <api-membrane.h>

#include <vt_common.h>





static int _verbose_ = 1;





/* static function definitions
 */

static char *_Array2Str( int argc, char *argv[] );
static char *_BaseName( char *p );
/*static double _GetTime();
static double _GetClock();
*/





int main( int argc, char *argv[] )
{
  lineCmdParamMembrane par;
  vt_image *image;
  vt_image *mask = NULL;
  char *lineoptions;

  vt_3Dimres imsRes;
  int flag_3D = 1;
/*
  clock_t start, stop;
  double elapsed;

  start = clock();

  double time_init = _GetTime();
  double time_exit;
  double clock_init = _GetClock();
  double clock_exit;
*/

  clock_t start, stop;
  double elapsed;

  start = clock();

  /* parameter initialization
   */
  API_InitParam_membrane( &par );



  /* parameter parsing
   */
  if ( argc <= 1 )
      API_ErrorParse_membrane( _BaseName( argv[0] ), (char*)NULL, 0 );
  API_ParseParam_membrane( 1, argc, argv, &par );
  


  /* input image reading
   */
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) {
      API_ErrorParse_membrane( _BaseName( argv[0] ), "unable to read input image ...\n", 0 );
  }

  if ( par.names.inv == 1 )  VT_InverseImage( image );
  if ( par.names.swap == 1 ) VT_SwapImage( image );


  if ( par.dimension == 2 || image->dim.z == 1 )
    flag_3D = 0;

  if (par.flagMask == 1)
  {
      mask = _VT_Inrimage( par.names.ext );
      if ( mask == (vt_image*)NULL || mask->dim.x != image->dim.x || mask->dim.y != image->dim.y ||mask->dim.z != image->dim.z ) {
          VT_FreeImage( image);
          API_ErrorParse_membrane( _BaseName( argv[0] ), "unable to read mask image ...\n", 0 );
      }
  }


  /* output image creation
   */

  if ( VT_Alloc3DImres( &imsRes, image, par.names.out, flag_3D ) != 1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    if (par.flagMask == 1) {
      VT_FreeImage( mask );
      VT_Free( (void**)&mask );
    }
    API_ErrorParse_membrane(_BaseName( argv[0] ), "unable to allocate response images ...\n", 0 );
  }


  /* API call response function
   */

  lineoptions = _Array2Str( argc, argv );
  if ( lineoptions == (char*)NULL )
      API_ErrorParse_membrane( _BaseName( argv[0] ), "unable to translate command line options ...\n", 0 );

  if ( API_membrane( image, &imsRes, mask, lineoptions, (char*)NULL ) != 1 ) {
      vtfree( lineoptions );
      VT_Free3DImres( &imsRes );
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      if (par.flagMask == 1) {
        VT_FreeImage( mask );
        VT_Free( (void**)&mask );
      }
      API_ErrorParse_membrane( _BaseName( argv[0] ), "some error occurs during processing ...\n", -1 );
  }


  /* output image writing
   */

  if ( par.writeImages )
      VT_Write3DImres( &imsRes, flag_3D );
  else
      VT_Write3DImresAngles( &imsRes, flag_3D );

  /* API call extrema
   */

  if ( API_extrema( &imsRes, mask, par.flagMask, lineoptions, (char*)NULL ) != 1 ) {
      vtfree( lineoptions );
      VT_Free3DImres( &imsRes );
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      if (par.flagMask == 1) {
        VT_FreeImage( mask );
        VT_Free( (void**)&mask );
      }
      API_ErrorParse_membrane( _BaseName( argv[0] ), "some error occurs during processing ...\n", -1 );
  }

  vtfree( lineoptions );

  /* output image writing
   */

  if ( VT_WriteInrimage( &(imsRes.imTheta) ) == -1 ) {
    VT_Free3DImres( &imsRes );
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    if (par.flagMask == 1) {
      VT_FreeImage( mask );
      VT_Free( (void**)&mask );
    }
    API_ErrorParse_membrane( _BaseName( argv[0] ), "unable to write output image ...\n", -1 );
  }
  
  /* memory freeing
   */
  VT_Free3DImres( &imsRes );
  VT_FreeImage( image );
  VT_Free( (void**)&image );
  if (par.flagMask == 1) {
    VT_FreeImage( mask );
    VT_Free( (void**)&mask );
  }



  if ( par.trace_allocations ) {
    fprintfVtMallocTrace( stderr );
    clearVtMalloc();
  }

/*
  time_exit = _GetTime();
  clock_exit = _GetClock();

  if ( par.print_time ) { 
    fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName( argv[0] ), time_exit - time_init );
    fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
    fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );
  }
*/

  stop = clock();
  elapsed = (double)(stop-start)/CLOCKS_PER_SEC;

  if ( par.print_time ) {
    fprintf(stdout, "Elapsed time : \t%.1fn\n", elapsed);
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

  s = (char*)vtmalloc( l * sizeof( char ), "s", proc );
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


/*
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
*/
