/*************************************************************************
 * test-2DcellProperties.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Jeu 25 oct 2018 17:12:22 CEST
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <bal-2DcellPropertiesTest.h>



static char *program = NULL;

static char *usage = "[-lmin %d] [-lmax %d]\n\
        [-amin %d] [-amax %d] [-adelta %d]\n\
        [-shape disk|square]\n\
        [-surface 4-edges 4-neighbors]";

static char *detail = "";



static int _verbose_ = 1;


static void _ErrorParse( char *str, int flag );
static char *_BaseName( char *p );
static void _WriteCmdLine( int argc, char *argv[], char *name );





/***************************************************
 *
 *
 *
 ***************************************************/





int main(int argc, char *argv[])
{
  int i, n;

  int lMin = 10;
  int lMax = 10;

  int angleMin = 0;
  int angleMax = 0;
  int angleDelta = 5;

  enumShape shape = _SQUARE_;
  enumSurfaceEstimation surface = _OUTER_4NEIGHBORS_;


  /***************************************************
   *
   * parsing parameters
   *
   ***************************************************/
  program = argv[0];

  /* no arguments
   */
  if ( argc == 1 ) _ErrorParse( NULL, 0 );

  for ( i=1; i<argc; i++ ) {
    if ( ( strcmp ( argv[i], "-help") == 0 )
        || ( strcmp ( argv[i], "-h") == 0 && argv[i][2] == '\0' )
        || ( strcmp ( argv[i], "--help") == 0 )
        || ( strcmp ( argv[i], "--h") == 0 && argv[i][3] == '\0' ) ) {
      _ErrorParse( NULL, 1 );
    }

    else if ( strcmp ( argv[i], "-lmin" ) == 0 ) {
      i ++;
      if ( i >= argc) _ErrorParse( "parsing -lmin", 0 );
      n = sscanf( argv[i], "%d", &lMin );
      if ( n != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "... -lmin argument: '%s' ?!\n", argv[i] );
        _ErrorParse( "parsing -lmin", 0 );
      }
    }
    else if ( strcmp ( argv[i], "-lmax" ) == 0 ) {
      i ++;
      if ( i >= argc) _ErrorParse( "parsing -lmax", 0 );
      n = sscanf( argv[i], "%d", &lMax );
      if ( n != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "... -lmax argument: '%s' ?!\n", argv[i] );
        _ErrorParse( "parsing -lmax", 0 );
      }
    }

    else if ( strcmp ( argv[i], "-amin" ) == 0 ) {
      i ++;
      if ( i >= argc) _ErrorParse( "parsing -amin", 0 );
      n = sscanf( argv[i], "%d", &angleMin );
      if ( n != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "... -amin argument: '%s' ?!\n", argv[i] );
        _ErrorParse( "parsing -amin", 0 );
      }
    }
    else if ( strcmp ( argv[i], "-amax" ) == 0 ) {
      i ++;
      if ( i >= argc) _ErrorParse( "parsing -amax", 0 );
      n = sscanf( argv[i], "%d", &angleMax );
      if ( n != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "... -amax argument: '%s' ?!\n", argv[i] );
        _ErrorParse( "parsing -amax", 0 );
      }
    }
    else if ( strcmp ( argv[i], "-adelta" ) == 0 ) {
      i ++;
      if ( i >= argc) _ErrorParse( "parsing -adelta", 0 );
      n = sscanf( argv[i], "%d", &angleDelta );
      if ( n != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "... -adelta argument: '%s' ?!\n", argv[i] );
        _ErrorParse( "parsing -adelta", 0 );
      }
    }

    else if ( strcmp ( argv[i], "-shape" ) == 0 ) {
      i ++;
      if ( i >= argc) _ErrorParse( "parsing -shape", 0 );
      if ( strcmp ( argv[i], "square" ) == 0 ) {
          shape = _SQUARE_;
      }
      else if ( strcmp ( argv[i], "disk" ) == 0 ) {
          shape = _DISK_;
      }
      else {
          if ( _verbose_ )
            fprintf( stderr, "... -shape argument: '%s' ?!\n", argv[i] );
          _ErrorParse( "parsing -shape", 0 );
      }
    }
    else if ( strcmp ( argv[i], "-surface" ) == 0 ) {
      i ++;
      if ( i >= argc) _ErrorParse( "parsing -surface", 0 );
      if ( strcmp ( argv[i], "4-edges" ) == 0 ) {
          surface = _4EDGES_;
      }
      else if ( strcmp ( argv[i], "4-neighbors" ) == 0 ) {
          surface = _OUTER_4NEIGHBORS_;
      }
      else {
          if ( _verbose_ )
            fprintf( stderr, "... -surface argument: '%s' ?!\n", argv[i] );
          _ErrorParse( "parsing -surface", 0 );
      }
    }


    else {
      if ( _verbose_ )
        fprintf( stderr, "... unknown option/argument: '%s'\n", argv[i] );
      _ErrorParse( "parsing arguments", 0 );
    }
  }

  /***************************************************
   *
   * parsing done
   *
   ***************************************************/
  _WriteCmdLine( argc, argv, (char*)NULL );


  if ( BAL_2DcellPropertiesShapeTest( lMin, lMax, angleMin, angleMax, angleDelta, shape, surface ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "... some error occurs\n" );
      exit( -1 );
  }


  return( 0 );
}





/***************************************************
 *
 *
 *
 ***************************************************/





static void _ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage: %s %s\n",_BaseName(program), usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s\n",detail);
  if ( str != NULL ) (void)fprintf(stderr,"Error: %s\n",str);
  exit( 1 );
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


static void _WriteCmdLine( int argc, char *argv[], char *name )
{
  char filename[STRINGLENGTH];
  int i;
  FILE *fcmdl;

  if ( name == (char*)NULL || name[0] == '\0' )
      sprintf( filename,"%s-%d-cmdline.log", _BaseName(argv[0]), getpid() );
  else
      sprintf( filename,"%s", name );

  fcmdl = fopen( filename, "w" );
  if ( fcmdl == NULL ) {
    if ( _verbose_ ) {
        fprintf( stderr, "%s: unable to open '%s' for writing\n", _BaseName( argv[0] ), filename );
    }
    return;
  }

  for ( i=0; i<argc; i++ ) {
    fprintf( fcmdl, "%s", argv[i] );
    if ( i < argc-1 )
      fprintf( fcmdl, " " );
  }
  fprintf( fcmdl, "\n" );

  fclose( fcmdl );
}
