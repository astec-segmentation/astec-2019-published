/*************************************************************************
 * test-3DcellProperties.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Lun 26 nov 2018 16:59:46 CET
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

#include <vtmalloc.h>

#include <bal-lineartrsf-tools.h>

#include <bal-3DcellPropertiesTest.h>



static char *program = NULL;

static char *usage = "[-ntests %d] [-sigma %f] [-L %d]\n\
        [-shape disk|cube|square|sphere]\n\
        [-s|-surface|-surface-estimation 6n|6-neighbors|lindblad|windreich]\n\
        [-image %s -minsigma %f -maxsigma %f]\n\
        [-seed %ld] [-scilab %s]";

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

  int L = 10;
  int nTests = 1;
  float sigma = 0.0;
  float maxSigma = 0.0;
  long int randomSeed = 0;
  char *basename = (char*)NULL;
  char *imagename = (char*)NULL;

  bal_image theIm;

  enumShape shape = _SQUARE_;
  enumSurfaceEstimation surface = _OUTER_6NEIGHBORS_;

  typeCellSequence cellProperties;


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

    else if ( strcmp ( argv[i], "-ntests" ) == 0 ) {
      i ++;
      if ( i >= argc) _ErrorParse( "parsing -ntests", 0 );
      n = sscanf( argv[i], "%d", &nTests );
      if ( n != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "... -ntests argument: '%s' ?!\n", argv[i] );
        _ErrorParse( "parsing -ntests", 0 );
      }
    }

    else if ( strcmp ( argv[i], "-sigma" ) == 0
              || strcmp ( argv[i], "-minsigma" ) == 0 ) {
      i ++;
      if ( i >= argc) _ErrorParse( "parsing -minsigma", 0 );
      n = sscanf( argv[i], "%f", &sigma );
      if ( n != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "... -minsigma argument: '%s' ?!\n", argv[i] );
        _ErrorParse( "parsing -minsigma", 0 );
      }
    }

    else if ( strcmp ( argv[i], "-maxsigma" ) == 0 ) {
      i ++;
      if ( i >= argc) _ErrorParse( "parsing -maxsigma", 0 );
      n = sscanf( argv[i], "%f", &maxSigma );
      if ( n != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "... -maxsigma argument: '%s' ?!\n", argv[i] );
        _ErrorParse( "parsing -maxsigma", 0 );
      }
    }

    else if ( (strcmp ( argv[i], "-L" ) == 0 && argv[i][2] == '\0')
              || (strcmp ( argv[i], "-l" ) == 0 && argv[i][2] == '\0') ) {
      i ++;
      if ( i >= argc) _ErrorParse( "parsing -L", 0 );
      n = sscanf( argv[i], "%d", &L );
      if ( n != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "... -L argument: '%s' ?!\n", argv[i] );
        _ErrorParse( "parsing -L", 0 );
      }
    }


    else if ( strcmp ( argv[i], "-shape" ) == 0 ) {
      i ++;
      if ( i >= argc) _ErrorParse( "parsing -shape", 0 );
      if ( strcmp ( argv[i], "square" ) == 0
           || strcmp ( argv[i], "cube" ) == 0 ) {
          shape = _SQUARE_;
      }
      else if ( strcmp ( argv[i], "disk" ) == 0
                || strcmp ( argv[i], "sphere" ) == 0 ) {
          shape = _DISK_;
      }
      else {
          if ( _verbose_ )
            fprintf( stderr, "... -shape argument: '%s' ?!\n", argv[i] );
          _ErrorParse( "parsing -shape", 0 );
      }
    }

    else if ( (strcmp ( argv[i], "-s" ) == 0 && argv[i][2] == '\0')
              || (strcmp ( argv[i], "-surface" ) == 0)
              || (strcmp ( argv[i], "-surface-estimation" ) == 0) ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -surface-estimation ...\n", 0 );
        if ( (strcmp ( argv[i], "6n" ) == 0 && argv[i][2] == '\0')
             || strcmp ( argv[i], "6-neighbors" ) == 0  ) {
            surface = _OUTER_6NEIGHBORS_;
        }
        else if ( strcmp ( argv[i], "lindblad" ) == 0 ) {
            surface = _LINDBLAD_;
        }
        else if ( strcmp ( argv[i], "windreich" ) == 0 ) {
            surface = _WINDREICH_;
        }
        else {
            fprintf( stderr, "unknown surface estimation type: '%s'\n", argv[i] );
            _ErrorParse( "parsing -surface-estimation ...", 0 );
        }
    }

    else if ( strcmp ( argv[i], "-image" ) == 0 && argv[i][6] == '\0' ) {
      i ++;
      if ( i >= argc) _ErrorParse( "parsing -image", 0 );
      imagename = argv[i];
    }

    else if ( strcmp ( argv[i], "-scilab" ) == 0 && argv[i][7] == '\0' ) {
      i ++;
      if ( i >= argc) _ErrorParse( "parsing -scilab", 0 );
      basename = argv[i];
    }

    else if ( strcmp ( argv[i], "-seed" ) == 0 && argv[i][5] == '\0' ) {
      i ++;
      if ( i >= argc) _ErrorParse( "parsing -seed", 0 );
      n = sscanf( argv[i], "%ld", &randomSeed );
      if ( n != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "... -seed argument: '%s' ?!\n", argv[i] );
        _ErrorParse( "parsing -seed", 0 );
      }
      _SetRandomSeed( randomSeed );
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
  if ( 0 ) _WriteCmdLine( argc, argv, (char*)NULL );

  if ( 0 ) fprintf( stderr, "random seed = %ld\n", _GetRandomSeed() );

  BAL_InitCellSequence( &cellProperties );
  if ( BAL_AllocCellSequence( &cellProperties, nTests ) != 1 ) {
    if ( _verbose_)
      fprintf( stderr, "... allocation of property list failed\n" );
    exit( -1 );
  }

  if ( imagename == (char*)NULL ) {

      if ( BAL_3DcellPropertiesShapeTest( &cellProperties, L, sigma, shape, surface, basename ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "... some error occurs\n" );
          exit( -1 );
      }

  }
  else {
      if ( BAL_ReadImage( &theIm, imagename, 0 ) != 1 ) {
        if ( _verbose_ )
            fprintf( stderr, "... can not read '%s'\n", imagename );
          exit( -1 );
      }
      if ( BAL_3DcellPropertiesImageTest( &cellProperties, &theIm, sigma, maxSigma, surface, basename ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "... some error occurs\n" );
          exit( -1 );
      }

  }

  BAL_FreeCellSequence( &cellProperties );


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
