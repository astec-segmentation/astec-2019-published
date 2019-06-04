/*************************************************************************
 * vt_link.c -
 *
 * $Id: vt_link.c,v 1.4 2006/05/16 09:33:34 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Wed Oct 11 09:22:02 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <vt_link.h>

static int _PTS_ARRAY_SIZE_ = 10000;
static int _verbose_ = 1;

int *VT_ExtractLinkedCurve( vt_image *theIm,
			    int *firstPt,
			    int *nbPts,
			    int sb,
			    int sh )
{
  char *proc = "VT_ExtractLinkedCurve";
  int *thePts = (int*)NULL;
  
  int dimx = theIm->dim.x;
  int dimy = theIm->dim.y;
  int dimz = theIm->dim.z;
  int nb, nv, v[26][3];
  int x, y, z;
  int i, j, k;

  *nbPts = 0;
  
  thePts = (int*)malloc( _PTS_ARRAY_SIZE_*3*sizeof(int) );
  if ( thePts == (int*)NULL ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate list of points\n", proc );
    return( (int*)NULL );
  }



  nb = 0;
  
  x = thePts[0] = firstPt[0];
  y = thePts[1] = firstPt[1];
  z = thePts[2] = firstPt[2];

  nb ++;



  switch ( theIm->type ) {
  default :
    if ( _verbose_ ) 
      fprintf( stderr, "%s: image type not handled yet\n", proc );
    free( thePts );
    return( (int*)NULL );
  case UCHAR :
    {
      u8 ***theBuf = (u8 ***)theIm->array;

      do {
	/* on cherche les voisins de (x,y,z)
	 */
	nv = 0;
	for ( k = -1; k <= 1; k ++ ) {
	  if ( z+k < 0 || z+k >= dimz ) continue;
	  for ( j = -1; j <= 1; j ++ ) {
	    if ( y+j < 0 || y+j >= dimy ) continue;
	    for ( i = -1; i <= 1; i ++ ) {
	      if ( x+i < 0 || x+i >= dimx ) continue;
	      if ( i == 0 && j == 0 && k == 0 ) continue;
	      if ( theBuf[z+k][y+j][x+i] < sb ) continue;
	      if ( theBuf[z+k][y+j][x+i] > sh ) continue;
	      v[nv][0] = x+i;
	      v[nv][1] = y+j;
	      v[nv][2] = z+k;
	      nv ++;
	    }
	  }
	}
	
	/* on analyse les voisins de (x, y, z)
	 */
	if ( nv <= 0 || nv > 2 ) {
	  
	  if ( _verbose_ ) 
	    fprintf( stderr, "%s: Warning, %d neighbors were found\n", proc, nv );
	  
	} else  if ( nv == 1 ) {
	  
	  if ( nb == 1 ) {
	    x = thePts[3*nb+0] = v[0][0];
	    y = thePts[3*nb+1] = v[0][1];
	    z = thePts[3*nb+2] = v[0][2];
	    nb ++;
	  } else 
	    nv = 0;
	  
	} else  if ( nv == 2 ) {
	  
	  if ( thePts[3*(nb-2)+0] == v[0][0] &&
	       thePts[3*(nb-2)+1] == v[0][1] &&
	       thePts[3*(nb-2)+2] == v[0][2] ) {
	    x = thePts[3*nb+0] = v[1][0];
	    y = thePts[3*nb+1] = v[1][1];
	    z = thePts[3*nb+2] = v[1][2];
	  } else {
	    x = thePts[3*nb+0] = v[0][0];
	    y = thePts[3*nb+1] = v[0][1];
	    z = thePts[3*nb+2] = v[0][2];
	  }
	  nb ++;
	}
	/* 1. debut de liste :
	   on entre avec nb = 1
	   on trouve     nv = 1
	   on ajoute le voisin => nb = 2
	   2.  nv = 2
	*/
      } while ( (nv == 2) || ( nv == 1 && nb == 2 ) );
    }
    break;

  case USHORT :
    {
      u16 ***theBuf = (u16 ***)theIm->array;

      do {
	/* on cherche les voisins de (x,y,z)
	 */
	nv = 0;
	for ( k = -1; k <= 1; k ++ ) {
	  if ( z+k < 0 || z+k >= dimz ) continue;
	  for ( j = -1; j <= 1; j ++ ) {
	    if ( y+j < 0 || y+j >= dimy ) continue;
	    for ( i = -1; i <= 1; i ++ ) {
	      if ( x+i < 0 || x+i >= dimx ) continue;
	      if ( i == 0 && j == 0 && k == 0 ) continue;
	      if ( theBuf[z+k][y+j][x+i] < sb ) continue;
	      if ( theBuf[z+k][y+j][x+i] > sh ) continue;
	      v[nv][0] = x+i;
	      v[nv][1] = y+j;
	      v[nv][2] = z+k;
	      nv ++;
	    }
	  }
	}
	
	/* on analyse les voisins de (x, y, z)
	 */
	if ( nv <= 0 || nv > 2 ) {
	  
	  if ( _verbose_ ) 
	    fprintf( stderr, "%s: Warning, %d neighbors were found\n", proc, nv );
	  
	} else  if ( nv == 1 ) {
	  
	  if ( nb == 1 ) {
	    x = thePts[3*nb+0] = v[0][0];
	    y = thePts[3*nb+1] = v[0][1];
	    z = thePts[3*nb+2] = v[0][2];
	    nb ++;
	  } else 
	    nv = 0;
	  
	} else  if ( nv == 2 ) {
	  
	  if ( thePts[3*(nb-2)+0] == v[0][0] &&
	       thePts[3*(nb-2)+1] == v[0][1] &&
	       thePts[3*(nb-2)+2] == v[0][2] ) {
	    x = thePts[3*nb+0] = v[1][0];
	    y = thePts[3*nb+1] = v[1][1];
	    z = thePts[3*nb+2] = v[1][2];
	  } else {
	    x = thePts[3*nb+0] = v[0][0];
	    y = thePts[3*nb+1] = v[0][1];
	    z = thePts[3*nb+2] = v[0][2];
	  }
	  nb ++;
	}
	/* 1. debut de liste :
	   on entre avec nb = 1
	   on trouve     nv = 1
	   on ajoute le voisin => nb = 2
	   2.  nv = 2
	*/
      } while ( (nv == 2) || ( nv == 1 && nb == 2 ) );
    }
    break;

  }
  
  if ( nv > 2 ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: linking failed for image %s\n\t stop at a junction (nv=%d)\n",
	       proc, theIm->name, nv );
    }
  }


  *nbPts = nb;
  return( thePts );
}









int VT_ExtractFirstPoint( vt_image *theIm,
			  int *firstPt,
			  int sb, int sh )
{
  char *proc = "VT_ExtractFirstPoint";
  int dimx = theIm->dim.x;
  int dimy = theIm->dim.y;
  int dimz = theIm->dim.z;
  int nv;
  int x, y, z;
  int i, j, k;

  firstPt[0] = firstPt[1] = firstPt[2] = 0;

  switch ( theIm->type ) {
  default :
    if ( _verbose_ ) 
      fprintf( stderr, "%s: image type not handled yet\n", proc );
    return( 0 );
  case UCHAR :
    {
      u8 ***theBuf = (u8 ***)theIm->array;

      for ( z=0; z<dimz; z++ )
      for ( y=0; y<dimy; y++ )
      for ( x=0; x<dimx; x++ ) {
    
	if ( theBuf[z][y][x] < sb ) continue;
	if ( theBuf[z][y][x] > sh ) continue;
	
	/* on cherche les voisins de (x,y,z)
	 */
	nv = 0;
	for ( k = -1; k <= 1; k ++ ) {
	  if ( z+k < 0 || z+k >= dimz ) continue;
	  for ( j = -1; j <= 1; j ++ ) {
	    if ( y+j < 0 || y+j >= dimy ) continue;
	    for ( i = -1; i <= 1; i ++ ) {
	      if ( x+i < 0 || x+i >= dimx ) continue;
	      if ( i == 0 && j == 0 && k == 0 ) continue;
	      if ( theBuf[z+k][y+j][x+i] < sb ) continue;
	      if ( theBuf[z+k][y+j][x+i] > sh ) continue;
	      nv ++;
	    }
	  }
	}

	if ( nv == 1 ) {
	  firstPt[0] = x;
	  firstPt[1] = y;
	  firstPt[2] = z;
	  return( 1 );
	}
      }
    }
    break;

  case USHORT :
    {
      u16 ***theBuf = (u16 ***)theIm->array;

      for ( z=0; z<dimz; z++ )
      for ( y=0; y<dimy; y++ )
      for ( x=0; x<dimx; x++ ) {
    
	if ( theBuf[z][y][x] < sb ) continue;
	if ( theBuf[z][y][x] > sh ) continue;
	
	/* on cherche les voisins de (x,y,z)
	 */
	nv = 0;
	for ( k = -1; k <= 1; k ++ ) {
	  if ( z+k < 0 || z+k >= dimz ) continue;
	  for ( j = -1; j <= 1; j ++ ) {
	    if ( y+j < 0 || y+j >= dimy ) continue;
	    for ( i = -1; i <= 1; i ++ ) {
	      if ( x+i < 0 || x+i >= dimx ) continue;
	      if ( i == 0 && j == 0 && k == 0 ) continue;
	      if ( theBuf[z+k][y+j][x+i] < sb ) continue;
	      if ( theBuf[z+k][y+j][x+i] > sh ) continue;
	      nv ++;
	    }
	  }
	}

	if ( nv == 1 ) {
	  firstPt[0] = x;
	  firstPt[1] = y;
	  firstPt[2] = z;
	  return( 1 );
	}
      }
    }
    break;

  }
  return( 0 );
}
