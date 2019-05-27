/*************************************************************************
 * bal_2DcellPropertiesTest.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Jeu 25 oct 2018 11:38:04 CEST
 *
 * ADDITIONS, CHANGES
 *
 */


/* random(), srandom(): Feature Test Macro Requirements for glibc
 * _SVID_SOURCE || _BSD_SOURCE || _XOPEN_SOURCE >= 500 || _XOPEN_SOURCE && _XOPEN_SOURCE_EXTENDED
 *
 * compilation with [gcc (GCC) 5.3.1 20151207 (Red Hat 5.3.1-2)] yields
 * "_BSD_SOURCE and _SVID_SOURCE are deprecated, use _DEFAULT_SOURCE"
 */
#define _XOPEN_SOURCE
#define _XOPEN_SOURCE_EXTENDED


#include <math.h>

#include <vtmalloc.h>

#include <bal-transformation-tools.h>

#include <bal-2DcellPropertiesTest.h>


static int _verbose_ = 1;


/************************************************************
 *
 *
 *
 ************************************************************/



int _count4Edges( bal_image *theIm, int lFore, int lBack )
{
  char *proc = "_count4Edges";
  size_t x, y, z;
  int n = 0;

  if ( 0 ) fprintf( stderr, "%s\n", proc );

#define _COUNT4EDGES( TYPE ) {   \
  u8 ***theBuf = (u8***)theIm->array;  \
  for ( z=0; z<theIm->nplanes; z++ )   \
  for ( y=0; y<theIm->nrows; y++ )     \
  for ( x=0; x<theIm->ncols; x++ ) {   \
    if ( x < theIm->ncols-1 ) {        \
      if ((theBuf[z][y][x] == lFore && theBuf[z][y][x+1] == lBack) \
              || (theBuf[z][y][x] == lBack && theBuf[z][y][x+1] == lFore)) \
          n ++;                        \
    }                                  \
    if ( y < theIm->nrows-1 ) {        \
      if ((theBuf[z][y][x] == lFore && theBuf[z][y+1][x] == lBack) \
              || (theBuf[z][y][x] == lBack && theBuf[z][y+1][x] == lFore)) \
          n ++;                        \
    }                                  \
    if ( z < theIm->nplanes-1 ) {      \
      if ((theBuf[z][y][x] == lFore && theBuf[z+1][y][x] == lBack) \
              || (theBuf[z][y][x] == lBack && theBuf[z+1][y][x] == lFore)) \
          n ++;                        \
    }                                  \
  }                                    \
}

  switch ( theIm->type ) {
  default :
    fprintf( stderr, "%s: type not handled in switch\n", proc );
    return( -1 );
  case UCHAR :
    _COUNT4EDGES( u8 );
    break;
  case USHORT :
    _COUNT4EDGES( u16 );
    break;
  }

  return( n );
}





int _count4Neighbors( bal_image *theIm, int lFore, int lBack )
{
  char *proc = "_count4Neighbors";
  size_t x, y, z;
  int n = 0;

  if ( 0 ) fprintf( stderr, "%s\n", proc );

#define _COUNT4NEIGHBORS( TYPE ) {   \
  u8 ***theBuf = (u8***)theIm->array;  \
  for ( z=0; z<theIm->nplanes; z++ )   \
  for ( y=0; y<theIm->nrows; y++ )     \
  for ( x=0; x<theIm->ncols; x++ ) {   \
    if ( theBuf[z][y][x] == lFore ) continue; \
    if ( theBuf[z][y][x] != lBack ) continue; \
    if ( x > 0 && theBuf[z][y][x-1] == lFore ) n++; \
    else if ( x < theIm->ncols-1 && theBuf[z][y][x+1] == lFore ) n++; \
    else if ( y > 0 && theBuf[z][y-1][x] == lFore ) n++; \
    else if ( y < theIm->nrows-1 && theBuf[z][y+1][x] == lFore ) n++; \
    else if ( z > 0 && theBuf[z-1][y][x] == lFore ) n++; \
    else if ( z < theIm->nplanes-1 && theBuf[z][y+1][x] == lFore ) n++; \
  }                                    \
}

  switch ( theIm->type ) {
  default :
    fprintf( stderr, "%s: type not handled in switch\n", proc );
    return( -1 );
  case UCHAR :
    _COUNT4NEIGHBORS( u8 );
    break;
  case USHORT :
    _COUNT4NEIGHBORS( u16 );
    break;
  }

  return( n );
}



/************************************************************
 *
 *
 *
 ************************************************************/




int BAL_2DcellPropertiesShapeTest( int dlMin, int dlMax, int angleMin, int angleMax, int angleDelta,
                            enumShape shape, enumSurfaceEstimation surface )
{
  char *proc = "BAL_2DcellPropertiesShapeTest";

  bal_image theIm;
  char imageName[256];
  float *voxelSize = (float*)NULL;
  u8 ***theBuf = (u8 ***)NULL;

  bal_transformation theTrsf;
  char trsfName[256];
  double *m;
  double center[2];

  int x, y, l, dl;
  float angle, dx, dy, d, dmax=10;

  int s1, s2;

  voxelSize = (float*)vtmalloc( (dlMax-dlMin+1)*sizeof(float), "voxelSize", proc );
  if ( voxelSize == (float*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }

  for ( dl=dlMin; dl<=dlMax; dl++ ) {

    l = 2*dl+1;

    for ( angle=angleMin; angle<=angleMax; angle += angleDelta ) {

      switch( shape ) {
      default :
        if ( _verbose_ )
          fprintf( stderr, "%s: such shape not handled yet\n", proc );
        vtfree( voxelSize );
        return( -1 );
      case _SQUARE_ :
        sprintf( imageName, "square_%06d_%03.0f.mha", l, angle );
        break;
      case _DISK_ :
        sprintf( imageName, "disk_%06d_%03.0f.mha", l, angle );
        break;
      }
      sprintf( trsfName, "trsf_%06d_%03.0f.trsf", l, angle );

      if ( BAL_AllocFullImage( &theIm, imageName, 3*l, 3*l, 1, 1, 1.0/(float)l, 1.0/(float)l, 1.0/(float)l, UCHAR ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: image allocation error\n", proc );
          vtfree( voxelSize );
          return( -1 );
      }

      switch( shape ) {
      default :
          if ( _verbose_ )
            fprintf( stderr, "%s: such shape not handled yet\n", proc );
          BAL_FreeImage( &theIm );
          vtfree( voxelSize );
          return( -1 );
      case _SQUARE_ :
          theBuf = (u8 ***)theIm.array;
          for ( y=0; y<3*l; y++ )
          for ( x=0; x<3*l; x++ ) {
            if ( y<l || 2*l<=y || x<l || 2*l<=x )
                theBuf[0][y][x] = 0;
            else
                theBuf[0][y][x] = 255;
          }
          break;
      case _DISK_ :
          theBuf = (u8 ***)theIm.array;
          center[0] = (3*l-1)/2.0;
          center[1] = (3*l-1)/2.0;
          dmax = 10;
          for ( y=0; y<3*l; y++ )
          for ( x=0; x<3*l; x++ ) {
            dx = (float)x - center[0];
            dy = (float)y - center[1];
            d = sqrt( dx*dx + dy*dy ) * 1.0/(float)l;
            if ( d > 0.5 ) {
                theBuf[0][y][x] = 0;
                if ( dmax > d ) dmax = d;
            }
            else
                theBuf[0][y][x] = 255;
          }
          break;
      }

      if ( angle > 0.0 ) {
        BAL_InitTransformation( &theTrsf );
        if ( BAL_AllocTransformation( &theTrsf, RIGID_2D, &theIm ) != 1 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: transformation allocation error\n", proc );
            BAL_FreeImage( &theIm );
            vtfree( voxelSize );
            return( -1 );
        }

        theTrsf.mat.m[0] = cos( -(angle*3.1415926535/180.0) );
        theTrsf.mat.m[1] = -sin( -(angle*3.1415926535/180.0) );
        theTrsf.mat.m[4] = sin( -(angle*3.1415926535/180.0) );
        theTrsf.mat.m[5] = cos( -(angle*3.1415926535/180.0) );

        m = theIm.to_real.m;
        center[0] = m[ 0] * (3*l-1)/2.0 + m[ 1] * (3*l-1)/2.0;
        center[1] = m[ 4] * (3*l-1)/2.0 + m[ 5] * (3*l-1)/2.0;
        theTrsf.mat.m[3] = center[0] - theTrsf.mat.m[0] * center[0] - theTrsf.mat.m[1] * center[1];
        theTrsf.mat.m[7] = center[1] - theTrsf.mat.m[4] * center[0] - theTrsf.mat.m[5] * center[1];
        /*
        theTrsf.mat.m[3] += ( (double)random()/(double)(RAND_MAX) - 0.5) * 1.0/(float)l;
        theTrsf.mat.m[7] += ( (double)random()/(double)(RAND_MAX) - 0.5) * 1.0/(float)l;
        */

        if ( BAL_ResampleImage( &theIm, &theIm, &theTrsf, NEAREST ) != 1 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: resampling error\n", proc );
            BAL_FreeTransformation( &theTrsf );
            BAL_FreeImage( &theIm );
            vtfree( voxelSize );
            return( -1 );
        }

        if ( 1 ) {
          if ( BAL_WriteTransformation( &theTrsf, trsfName ) != 1 ) {
              if ( _verbose_ )
                fprintf( stderr, "%s: error when writing transformation '%s'\n", proc, trsfName );
              BAL_FreeTransformation( &theTrsf );
              BAL_FreeImage( &theIm );
              vtfree( voxelSize );
              return( -1 );
          }
        }

        BAL_FreeTransformation( &theTrsf );
      }

      if ( 1 ) {
        if ( BAL_WriteImage( &theIm, imageName ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: error when writing image '%s'\n", proc, imageName );
          BAL_FreeImage( &theIm );
          vtfree( voxelSize );
          return( -1 );
        }
      }

      switch( surface ) {
      default :
      case _4EDGES_ :
          s1 = _count4Edges(&theIm, 255, 0);
          s2 = _count4Edges(&theIm, 0, 255);
          break;
      case _OUTER_4NEIGHBORS_ :
          s1 = _count4Neighbors(&theIm, 255, 0);
          s2 = _count4Neighbors(&theIm, 0, 255);
          break;
      }

      switch( shape ) {
      default :
      case _SQUARE_ :
        fprintf( stderr, "voxel=%6.4f, angle=%03.0f", 1.0/(float)l, angle );
        fprintf( stderr, ", l(%d,%d)=%d [%f]", 255, 0, s1, s1*1.0/(float)l );
        fprintf( stderr, ", l(%d,%d)=%d [%f]", 0, 255, s2, s2*1.0/(float)l );
        fprintf( stderr, "\n" );
        break;
      case _DISK_ :
        fprintf( stderr, "voxel=%6.4f, angle=%03.0f", 1.0/(float)l, angle );
        fprintf( stderr, ", l(%d,%d)=%d [%f]", 255, 0, s1, s1*1.0/(float)l );
        fprintf( stderr, ", l(%d,%d)=%d [%f]", 0, 255, s2, s2*1.0/(float)l );
        fprintf( stderr, ", rnext=%8.6f", dmax );
        fprintf( stderr, "\n" );
        break;
      }

      BAL_FreeImage( &theIm );
    }
  }

  vtfree( voxelSize );
  return( 1 );
}




