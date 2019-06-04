/*************************************************************************
 * vt_maskSeeds.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Jeu 26 jul 2018 10:21:24 CEST
 *
 * ADDITIONS, CHANGES
 *
 *
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vtmalloc.h>

#include <vt_copy.h>

#include <vt_maskSeeds.h>


static int _verbose_ = 1;


/*************************************************************************
 *
 *
 *
 *************************************************************************/

typedef struct _seed {
  int first_cell_id;
  int toberemoved;
} _seed;


int VT_SelectSeedsinCells( vt_image *imseeds, vt_image *imcells, vt_image *imres )
{
  char *proc = "VT_SelectSeedsinCells";
  int nmaxseeds = 0;
  _seed *seedList = (_seed*)NULL;
  int n, r;
  size_t i, v;

  switch( imseeds->type ) {
  default :
      if ( _verbose_ )
          fprintf( stderr, "%s: such seed image type not handled yet\n", proc );
      return( -1 );
  case UCHAR :
      nmaxseeds = 256;
      break;
  case SSHORT :
      nmaxseeds = 32768;
      break;
  case USHORT :
      nmaxseeds = 65536;
      break;
  }

  seedList = (_seed*)vtmalloc( nmaxseeds * sizeof(_seed), "seedList", proc );
  if ( seedList == (_seed*)NULL ) {
      if ( _verbose_ )
          fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
  }

  for ( n=0; n<nmaxseeds; n++ ) {
      seedList[n].first_cell_id = -1;
      seedList[n].toberemoved = 0;
  }

  v = imseeds->dim.x * imseeds->dim.y * imseeds->dim.z;


#define _CELLCASE( TYPE ) {                                            \
  TYPE *theCells = (TYPE*)imcells->buf;                                \
  for ( i=0; i<v; i++ ) {                                              \
    if ( theSeeds[i] == 0 ) continue;                                  \
    if ( seedList[ theSeeds[i] ].toberemoved == 1 ) continue;          \
    if ( seedList[ theSeeds[i] ].first_cell_id == -1 ) {               \
      seedList[ theSeeds[i] ].first_cell_id = theCells[i];             \
    }                                                                  \
    else if ( seedList[ theSeeds[i] ].first_cell_id == theCells[i] ) { \
      continue;                                                        \
    }                                                                  \
    else {                                                             \
      seedList[ theSeeds[i] ].toberemoved = 1;                         \
    }                                                                  \
  }                                                                    \
}

#define _SEEDCASE( TYPE ) {             \
  TYPE *theSeeds = (TYPE*)imseeds->buf; \
  switch( imcells->type ) {             \
  default :                             \
    vtfree( seedList );                 \
    if ( _verbose_ )                    \
        fprintf( stderr, "%s: such cell image type not handled yet\n", proc ); \
    return( -1 );                       \
  case UCHAR :                          \
    _CELLCASE( u8 );                    \
    break;                              \
  case SSHORT :                         \
    _CELLCASE( s16 );                   \
    break;                              \
  case USHORT :                         \
    _CELLCASE( u16 );                   \
    break;                              \
  }                                     \
}

  switch( imseeds->type ) {
  default :
    vtfree( seedList );
    if ( _verbose_ )
        fprintf( stderr, "%s: such seed image type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    _SEEDCASE( u8 );
    break;
  case SSHORT :
    _SEEDCASE( s16 );
    break;
  case USHORT :
    _SEEDCASE( u16 );
    break;
  }



  for ( r=0, n=0; n<nmaxseeds; n++ ) {
      if ( seedList[n].toberemoved == 1 ) r++;
  }

  if ( _verbose_ ) {
    if ( r == 0 ) {
      if ( _verbose_)
        fprintf( stderr, "%s: there is no seed to be removed\n", proc );
    }
    else {
      if ( _verbose_)
        fprintf( stderr, "%s: there are %d seeds to be removed:", proc, r );
      for ( r=0, n=0; n<nmaxseeds; n++ ) {
        if ( seedList[n].toberemoved == 1 ) {
          if ( r > 0 ) fprintf( stderr, "," );
          fprintf( stderr, " %d", n );
          r++;
        }
      }
    }
  }



  if ( VT_CopyImage( imseeds, imres ) != 1 ) {
      vtfree( seedList );
      if ( _verbose_ )
          fprintf( stderr, "%s: error when copying seed image\n", proc );
      return( -1 );
  }

#define _RESCASE( TYPE ) {                          \
  TYPE *resBuf = (TYPE*)imres->buf;                 \
  for ( i=0; i<v; i++ ) {                           \
      if ( resBuf[i] == 0 ) continue;               \
      if ( seedList[ resBuf[i] ].toberemoved == 1 ) \
          resBuf[i] = 0;                            \
  }                                                 \
}

  if ( r > 0 ) {
    switch( imres->type ) {
    default :
      vtfree( seedList );
      if ( _verbose_ )
          fprintf( stderr, "%s: such result image type not handled yet\n", proc );
      return( -1 );
    case UCHAR :
      _RESCASE( u8 );
      break;
    case SSHORT :
      _RESCASE( s16 );
      break;
    case USHORT :
      _RESCASE( u16 );
      break;
    }
  }


  vtfree( seedList );


  return( 1 );
}
