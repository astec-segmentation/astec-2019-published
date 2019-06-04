/*************************************************************************
 * vt_bound.c - extraction de parametres de box sur des parties numerotees
 *
 * $Id: vt_bound.c,v 1.13 2013/08/06 17:50:52 gael Exp $
 *
 * DESCRIPTION:
 *
 *
 * Copyright Gregoire Malandain, INRIA
 *
 *
 * ADDITIONS, CHANGES:
 *
 *
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#include <vt_bound.h>


static int _verbose_ = 0;

void _VerboseInBound()
{
  _verbose_ = 1;
}
void _NoVerboseInBound()
{
  _verbose_ = 0;
}










int _CreateArrayOfParametersFromImage( vt_image *theIm,
                                       int slice,
                                       typeBoxParameters **thePar )
{
  char *proc="_CreateArrayOfParametersFromImage";
  int n;
  int s=slice;
  typeBoxParameters *theComp = (typeBoxParameters *)NULL;


  /* 2D ou 3D ?
   */
  if ( theIm->dim.z == 1 ) s = 0;
  if ( (theIm->dim.z > 1) && (s < 0 || s >= (int)theIm->dim.z) ) s = -1;


  /* c'est pas optimal. mais ca ira bien pour l'instant
     on alloue n+1 composantes pour adresser chaque
     element par son intensite
   */
  n = _MaxValueInImage( theIm, s );
  if ( n == 0 ) {
    if ( _VT_VERBOSE_ || _VT_DEBUG_ ) {
      fprintf( stderr, "%s: null image.\n", proc );
    }
    return( 0 );
  }

  theComp = (typeBoxParameters *)malloc( (n+1)*sizeof(typeBoxParameters) );
  if ( theComp == (typeBoxParameters *)NULL ) {
    if ( _VT_VERBOSE_ || _VT_DEBUG_ ) {
      fprintf( stderr, "%s: can not allocate components list.\n", proc );
    }
    return( -1 );
  }

  _InitArrayOfParametersFromImage( theIm, theComp, n );
  _FillArrayOfParametersFromImage( theIm, s, theComp );



  *thePar = theComp;
  return( n );
}



















void _InitArrayOfParametersFromImage( vt_image *theIm,
                                     typeBoxParameters *thePar,
                                     int n )
{
  int i;
  for ( i=0; i<=n; i++ ) {
    thePar[i].ptmin[0] = theIm->dim.x-1;
    thePar[i].ptmin[1] = theIm->dim.y-1;
    thePar[i].ptmin[2] = theIm->dim.z-1;
    thePar[i].ptmax[0] = 0;
    thePar[i].ptmax[1] = 0;
    thePar[i].ptmax[2] = 0;
    thePar[i].volume = 0;
  }
}


void _FillArrayOfParametersFromImage( vt_image *theIm,
                                      int slice,
                                      typeBoxParameters *thePar )
{
  char *proc = "_FillArrayOfParametersFromImage";
  int s=slice;
  int x, y, z, color;

  /* 2D ou 3D ?
   */
  if ( theIm->dim.z == 1 ) s = 0;
  if ( (theIm->dim.z > 1) && (s < 0 || s >= (int)theIm->dim.z) ) s = -1;

  switch ( theIm->type ) {
  case USHORT :
    {
      u16 *** theBuf = (u16 ***)theIm->array;

      /* cas 2D
       */
      if ( theIm->dim.z == 1 || s >= 0 ) {
        for ( y=0; y<(int)theIm->dim.y; y++ )
        for ( x=0; x<(int)theIm->dim.x; x++ ) {

          if ( theBuf[s][y][x] == 0 ) continue;

          color = theBuf[s][y][x];

          thePar[ color ].volume ++;
          if ( x < thePar[ color ].ptmin[0] ) thePar[ color ].ptmin[0] = x;
          if ( y < thePar[ color ].ptmin[1] ) thePar[ color ].ptmin[1] = y;
          if ( x > thePar[ color ].ptmax[0] ) thePar[ color ].ptmax[0] = x;
          if ( y > thePar[ color ].ptmax[1] ) thePar[ color ].ptmax[1] = y;
          thePar[ color ].ptmin[2] = thePar[ color ].ptmax[2] = s;
        }
        break;
      }

      /* cas 3D
       */
      for ( z=0; z<(int)theIm->dim.z; z++ )
      for ( y=0; y<(int)theIm->dim.y; y++ )
      for ( x=0; x<(int)theIm->dim.x; x++ ) {

        if ( theBuf[z][y][x] == 0 ) continue;

        color = theBuf[z][y][x];

        thePar[ color ].volume ++;
        if ( x < thePar[ color ].ptmin[0] ) thePar[ color ].ptmin[0] = x;
        if ( y < thePar[ color ].ptmin[1] ) thePar[ color ].ptmin[1] = y;
        if ( z < thePar[ color ].ptmin[2] ) thePar[ color ].ptmin[2] = z;
        if ( x > thePar[ color ].ptmax[0] ) thePar[ color ].ptmax[0] = x;
        if ( y > thePar[ color ].ptmax[1] ) thePar[ color ].ptmax[1] = y;
        if ( z > thePar[ color ].ptmax[2] ) thePar[ color ].ptmax[2] = z;

      }

    }
    break;
  default :
    if ( _VT_VERBOSE_ || _VT_DEBUG_ )
      fprintf( stderr, "%s: such image type is not handled yet.\n", proc );
    return;
  }

}










int _MaxValueInImage( vt_image *theIm,
                      int slice )
{
  char *proc = "_MaxValueInImage";
  int v, offset=0;
  int volume = theIm->dim.x * theIm->dim.y * theIm->dim.z;
  int max=0, s = slice;


  /* 2D ou 3D ?
   */
  if ( theIm->dim.z == 1 ) s = 0;
  if ( (theIm->dim.z > 1) && (s < 0 || s >= (int)theIm->dim.z) ) s = -1;
  /* cas 2D
   */
  if ( theIm->dim.z == 1 || s >= 0 ) {
    offset = s * theIm->dim.x * theIm->dim.y;
    volume = theIm->dim.x * theIm->dim.y;
  }

  switch ( theIm->type ) {
  case UCHAR :
    {
      u8 *theBuf = (u8*)theIm->buf;
      max = theBuf[offset];
      for ( v=1; v<volume; v++ )
        if ( theBuf[offset+v] > max ) max = theBuf[offset+v];
    }
    break;
  case USHORT :
    {
      u16 *theBuf = (u16*)theIm->buf;
      max = theBuf[offset];
      for ( v=1; v<volume; v++ )
        if ( theBuf[offset+v] > max ) max = theBuf[offset+v];
    }
    break;
  default :
    if ( _VT_VERBOSE_ || _VT_DEBUG_ )
      fprintf( stderr, "%s: such image type is not handled yet.\n", proc );
    return( 0 );
  }

  return( max );
}








