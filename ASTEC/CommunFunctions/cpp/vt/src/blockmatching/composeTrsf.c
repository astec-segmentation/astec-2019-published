/*************************************************************************
 * composeTrsf.c - template for executable creation
 *
 * $Id$
 *
 * Copyright (c) INRIA 2015, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Jeu  2 jul 2015 13:53:42 CEST
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

#include <api-composeTrsf.h>







static int _verbose_ = 1;





/* static function definitions
 */

static char *_Array2Str( int argc, char *argv[] );
static char *_BaseName( char *p );
static double _GetTime();
static double _GetClock();






int main( int argc, char *argv[] )
{
  lineCmdParamComposeTrsf par;
  char *lineoptions;
  int n;


  double time_init = _GetTime();
  double time_exit;
  double clock_init = _GetClock();
  double clock_exit;


  /* parameter initialization
   */
  API_InitParam_composeTrsf( &par );



  /* parameter parsing
   */
  if ( argc <= 1 )
      API_ErrorParse_composeTrsf( _BaseName( argv[0] ), (char*)NULL, 0 );
  API_ParseParam_composeTrsf( 1, argc, argv, &par );
  

  /* names building
   */
  n = par.input_names.n_data;
  if ( _verbose_ >= 2 && n > 0 )
      fprintf( stderr, "%s: there are %d transformations names from the command line\n",
               _BaseName( argv[0] ), n );
  if ( par.input_list != (char*)NULL && par.input_list[0] != '\0' ) {
      if ( buildStringListFromFile( par.input_list, &(par.input_names) ) != 1 ) {
          API_FreeParam_composeTrsf( &par );
          API_ErrorParse_composeTrsf( _BaseName( argv[0] ), "unable to build transformation names from list ...\n", 0 );
      }
      if ( _verbose_ >= 2 && par.input_names.n_data > n )
          fprintf( stderr, "%s: there are %d transformation names from the list file\n",
                         _BaseName( argv[0] ), par.input_names.n_data - n );
      n = par.input_names.n_data;
  }
  if ( par.input_format != (char*)NULL && par.input_format[0] != '\0' ) {
      if ( buildStringListFromFormat( par.input_format, par.firstindex, par.lastindex, &(par.input_names) ) != 1 ) {
          API_FreeParam_composeTrsf( &par );
          API_ErrorParse_composeTrsf( _BaseName( argv[0] ), "unable to build transformation names from format ...\n", 0 );
      }
      if ( _verbose_ >= 2 && par.input_names.n_data > n )
          fprintf( stderr, "%s: there are %d transformation names from format\n",
                         _BaseName( argv[0] ), par.input_names.n_data - n );
      n = par.input_names.n_data;
  }

  if ( _verbose_ >= 2 )
      printStringList( stderr, &(par.input_names), "Input transformations" );


  /* API call
   */

  lineoptions = _Array2Str( argc, argv );
  if ( lineoptions == (char*)NULL ) {
      API_FreeParam_composeTrsf( &par );
      API_ErrorParse_composeTrsf( _BaseName( argv[0] ), "unable to translate command line options ...\n", 0 );
  }

  if ( API_INTERMEDIARY_composeTrsf( &(par.input_names), par.output_name,
                                     par.template_name,
                                    lineoptions, (char*)NULL ) != 1 ) {
      vtfree( lineoptions );
      API_FreeParam_composeTrsf( &par );
      API_ErrorParse_composeTrsf( _BaseName( argv[0] ), "some error occurs during processing ...\n", -1 );
  }
  vtfree( lineoptions );



  if ( par.trace_allocations ) {
    fprintfVtMallocTrace( stderr );
    clearVtMalloc();
  }

  time_exit = _GetTime();
  clock_exit = _GetClock();

  if ( par.print_time ) { 
    fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName( argv[0] ), time_exit - time_init );
    fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
    fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );
  }

  API_FreeParam_composeTrsf( &par );


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
