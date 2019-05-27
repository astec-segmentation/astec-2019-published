/*************************************************************************
 * ccparameters.c - extraction de parametres sur des parties numerotees
 *
 * $Id: ccparameters.c,v 1.2 2001/04/03 10:27:21 greg Exp $
 *
 * DESCRIPTION: 
 *
 *
 *
 * AUTHOR:
 * Gregoire Malandain
 *
 * CREATION DATE:
 * Sun Feb 11 10:45:53 MET 2001
 *
 * Copyright Gregoire Malandain, INRIA
 *
 *
 * ADDITIONS, CHANGES:
 *
 *
 */

#include <stdlib.h>
extern long random(void);

#include <stdio.h>
#include <math.h>
#include <string.h>


#include <ccparameters.h>




static int _verbose_ = 0;

void _VerboseInCcParameters()
{
  _verbose_ = 1;
}
void _NoVerboseInCcParameters()
{
  _verbose_ = 0;
}







typedef struct {
  int ptmin[3];
  int ptmax[3];
} typeExtendedParameter;




typedef struct {
  double x;
  double y;
  double z;
} typePoint;







typedef struct {

  int n;
  int nbAllocs;

  typePoint *point; /* calcul     du diametre max: point "normal" */
  typePoint *projs; /* estimation du diametre med: point "projetes" 
		       dans la direction du diametre max 
		    */
  double *abscisse; /* estimation du diametre min */

} typeListe;




typedef enum {
  _ESTIMATE_,
  _COMPUTE_,
} enumMaxFeret;

static enumMaxFeret typeMaxFeret = _ESTIMATE_;
static double _EPSILON_FOR_NORMS_ = 0.0001;






static double _ComputeMaxSquareDistanceInList( typePoint *theList,
					    int n,
					    typePoint *first,
					    typePoint *last );
static double _EstimateMaxSquareDistanceInList( typePoint *theList,
					     int n,
					     typePoint *first,
					     typePoint *last  );
static double _ComputeMaxDifferenceInList( double *theList,
					   int n );




static void _ComputeFeretDiametersIn3DList( typeListe *theList,
					    typeParameter *thePar );



static void _FillListWithBorderPoints( unsigned short int *theBuf,
				       const int *theDim,
				       typeListe *theList,
				       int color,
				       int *minCorner,
				       int *maxCorner );


#ifdef _UNUSED_
static void _PrintTypeListe( typeListe *theListe );
#endif
static void _InitTypeListe( typeListe *theListe );
static void _FreeTypeListe( typeListe *theListe );
static int _AllocTypeListe( typeListe *theListe, int n );
static void _ApplyVoxelDimensionsToList( typeListe *theListe,
					 const double *sizes );






typeParameter *ComputeParameterFromLabels( unsigned short int *theBuf,
					   const int *theDim,
					   const double *voxelSizes,
					   const int nbLabels)
{
  typeParameter *thePar = NULL;
  typeExtendedParameter *theExt = NULL;
  typeListe theList;
  int x, y, z, i, n;

  if ( nbLabels <= 0 ) return( NULL );

  /* first allocations 
   */
  thePar = (typeParameter *)malloc( (nbLabels+1)*sizeof(typeParameter) );
  if ( thePar == NULL ) return( NULL );
  theExt = (typeExtendedParameter *)malloc( (nbLabels+1)*sizeof(typeExtendedParameter) );
  if ( theExt == NULL ) {
    free( thePar );
    return( NULL );
  }

  /* initialisations 
   */
  for ( n=0; n<=nbLabels; n++ ) {

    thePar[n].volume = 0;
    thePar[n].maxDiameter = 0.0;
    thePar[n].medDiameter = 0.0;
    thePar[n].minDiameter = 0.0;
    
    theExt[n].ptmin[0] = theDim[0]-1;
    theExt[n].ptmin[1] = theDim[1]-1;
    theExt[n].ptmin[2] = theDim[2]-1;
    theExt[n].ptmax[0] = 0;
    theExt[n].ptmax[1] = 0;
    theExt[n].ptmax[2] = 0;
  }

  /* first parameters
   */
  for ( i=0, z=0; z<theDim[2]; z++ ) 
  for ( y=0; y<theDim[1]; y++ ) 
  for ( x=0; x<theDim[0]; x++, i++ ) {
  
    if ( theBuf[i] == 0 )       continue;
    if ( theBuf[i] > nbLabels ) continue;
    
    n = theBuf[i];
    
    thePar[n].volume ++;
    if ( x < theExt[n].ptmin[0] ) theExt[n].ptmin[0] = x;
    if ( y < theExt[n].ptmin[1] ) theExt[n].ptmin[1] = y;
    if ( z < theExt[n].ptmin[2] ) theExt[n].ptmin[2] = z;
    if ( x > theExt[n].ptmax[0] ) theExt[n].ptmax[0] = x;
    if ( y > theExt[n].ptmax[1] ) theExt[n].ptmax[1] = y;
    if ( z > theExt[n].ptmax[2] ) theExt[n].ptmax[2] = z;
    
  }



  /* maximum volume and last allocation 
   */
  i = thePar[1].volume;
  for ( n=2; n<=nbLabels; n++ )
    if ( i < thePar[n].volume ) i = thePar[n].volume;
  _InitTypeListe( &theList );
  if ( _AllocTypeListe( &theList, i ) != 1 ) {
    free( theExt );
    free( thePar );
    return( NULL );
  }

  /* more parameters
   */
  for ( n=1; n<=nbLabels; n++ ) {
    _FillListWithBorderPoints( theBuf, theDim, &theList, n, 
			       theExt[n].ptmin, theExt[n].ptmax );
    _ApplyVoxelDimensionsToList( &theList, voxelSizes );
    _ComputeFeretDiametersIn3DList( &theList, &thePar[n] );
  }


  _FreeTypeListe( &theList );
  free( theExt );
  return( thePar );
}




















void _SetFeretDiameterCalculusToComputation()
{
  typeMaxFeret = _COMPUTE_;
}
void _SetFeretDiameterCalculusToEstimation()
{
  typeMaxFeret = _ESTIMATE_;
}



















static double _ComputeMaxSquareDistanceInList( typePoint *theList,
					       int n,
					       typePoint *first,
					       typePoint *last )
{
  int i, j;
  double max, d;

  if ( n == 1 ) {
    *first = *last = theList[0];
    return( (double)0.0 );
  }

  max = d = 0.0;
  for ( i=0; i<n-1; i++ ) 
  for ( j=i+1; j<n; j++ ) {
    d = (theList[i].z - theList[j].z) * (theList[i].z - theList[j].z) 
      + (theList[i].y - theList[j].y) * (theList[i].y - theList[j].y) 
      + (theList[i].x - theList[j].x) * (theList[i].x - theList[j].x);
    if ( d > max ) {
      max = d;
      *first = theList[i];
      *last  = theList[j];
    }
  }

  return( max );
}




















































static double _EstimateMaxSquareDistanceInList( typePoint *theList,
						int nbpts,
						typePoint *first,
						typePoint *last  )
{
  int i, i2, j, n;
  typePoint p1, p2;
  typePoint pf, pl;
  double max, old, m, d;
  int endIsReached = 0;
  int nb=nbpts;
  
  
  if ( nb == 1 ) {
    *first = *last = theList[0];
    return( (double)0.0 );
  } 

  max = 0;
  

  /* on tire un point au hasard dans les n premiers
     on cherche le plus eloigne
     on echange les points
     on itere jusqu'a stabilite
  */

  i2 = (int)( (nb - 1) * ((double)random() / (double)2147483647.0) + 0.5);

  do {

    m = 0;

    do {

      p1 = theList[i2];  
      theList[i2] = theList[nb-1]; theList[nb-1] = p1;
      nb --;
      
      old = m;
      for ( i=0; i<nb; i++ ) {
	d = (theList[i].z - p1.z) * (theList[i].z - p1.z) 
	  + (theList[i].y - p1.y) * (theList[i].y - p1.y) 
	  + (theList[i].x - p1.x) * (theList[i].x - p1.x);
	if ( d > m ) {
	  m = d;
	  p2 = theList[i];
	  i2 = i;
	  pf = p1;
	  pl = p2;
	}
      }

    } while ( m > old );
    
  
    if ( m > max ) {
      /* on a un maximum
       */
      max    = m;
      *first = pf;
      *last  = pl;
      

      /* on cherche le point le plus eloigne du cercle
       */
      i2 = -1;
      m = 0.0;
      for ( i=0; i<nb; i++ ) {
	d = (theList[i].z - pf.z) * (theList[i].z - pl.z) 
	  + (theList[i].y - pf.y) * (theList[i].y - pl.y) 
	  + (theList[i].x - pf.x) * (theList[i].x - pl.x);
	if ( m < d ) {
	  i2 = i;
	  m  = d;
	}
      }

      if ( i2 < 0 ) return( max );
      
    } else {

      endIsReached = 1;
    }

  } while( endIsReached == 0 );
  

  /* on met les points en dehors du cercle
     devant la liste
  */
  n = -1;
  for ( i=0; i<nb; i++ ) {
    d = (theList[i].z - pf.z) * (theList[i].z - pl.z) 
      + (theList[i].y - pf.y) * (theList[i].y - pl.y) 
      + (theList[i].x - pf.x) * (theList[i].x - pl.x);
    if ( d > 0 ) {
      n ++;
      p1 = theList[n]; theList[n] = theList[i]; theList[i] = p1;
    }
  }
  
  if (n <0) return( max );


  

  /* on finit par une recherche exhaustive
   */
  for ( i=0; i<=n; i++ ) 
  for ( j=i+1; j<nb; j++ ) {
    d = (theList[i].z - theList[j].z) * (theList[i].z - theList[j].z) 
      + (theList[i].y - theList[j].y) * (theList[i].y - theList[j].y) 
      + (theList[i].x - theList[j].x) * (theList[i].x - theList[j].x);
    if ( d > max ) {
      max = d;
      *first = theList[i];
      *last  = theList[j];
    }
  }

  return( max );
}























static double _ComputeMaxDifferenceInList( double *theList,
					   int n )
{
  int i;
  double max, min;

  min = max = theList[0];
  for ( i=1; i<n; i++ ) {
    if ( theList[i] < min ) min = theList[i];
    else if ( theList[i] > max ) max = theList[i];
  }
  return( max - min );
}

























static void _ComputeFeretDiametersIn3DList( typeListe *theList,
					    typeParameter *thePar )
{
  typePoint pt1, pt2;
  double maxDirection[3];
  double medDirection[3];
  double minDirection[3];
  double n;
  int i;
  double ps;



  switch ( typeMaxFeret ) {
  default :
  case _COMPUTE_ :
    (void)_ComputeMaxSquareDistanceInList( theList->point, theList->n,
					   &pt1, &pt2 );
    break;
  case _ESTIMATE_ :
    (void)_EstimateMaxSquareDistanceInList( theList->point, theList->n,
					    &pt1, &pt2 );
    break;
  }


  /* directions
   */
  maxDirection[0] = (double)(pt2.x - pt1.x);
  maxDirection[1] = (double)(pt2.y - pt1.y);
  maxDirection[2] = (double)(pt2.z - pt1.z);

  n = maxDirection[0]*maxDirection[0] + maxDirection[1]*maxDirection[1] 
    + maxDirection[2]*maxDirection[2];

  if ( _verbose_ ) {
    printf("MAX = (%9.5g %9.5g %9.5g) - (%9.5g %9.5g %9.5g) = %g\n",
	   pt1.x, pt1.y, pt1.z, pt2.x, pt2.y, pt2.z, n );
  }

  /* c'est un point isole
   */
  if ( n <= _EPSILON_FOR_NORMS_ ) {
    thePar->maxDiameter = 0.0;
    maxDirection[0] = maxDirection[1] = maxDirection[2] = 0.0;
    thePar->medDiameter = 0.0;
    medDirection[0] = medDirection[1] = medDirection[2] = 0.0;
    thePar->minDiameter = 0.0;
    minDirection[0] = minDirection[1] = minDirection[2] = 0.0;
    return;
  }

  n = sqrt( n );
  maxDirection[0] /= n;
  maxDirection[1] /= n;
  maxDirection[2] /= n;

  thePar->maxDiameter = n;

  if ( _verbose_ ) {
    printf( "   DIR = (%9.5g %9.5g %9.5g)\n", 
	    maxDirection[0], maxDirection[1], maxDirection[2] );
  }
  
  
  /* on enleve la projection du point sur la premiere direction
     ainsi tous les points seront coplanaires
  */
  for ( i=0; i<theList->n; i++ ) {
    ps = theList->point[i].x * maxDirection[0] +
      theList->point[i].y * maxDirection[1] +
      theList->point[i].z * maxDirection[2];
    theList->projs[i].x = theList->point[i].x - ps * maxDirection[0];
    theList->projs[i].y = theList->point[i].y - ps * maxDirection[1];
    theList->projs[i].z = theList->point[i].z - ps * maxDirection[2];
  }


  switch ( typeMaxFeret ) {
  default :
  case _COMPUTE_ :
    (void)_ComputeMaxSquareDistanceInList( theList->projs, theList->n,
					   &pt1, &pt2 );
    break;
  case _ESTIMATE_ :
    (void)_EstimateMaxSquareDistanceInList( theList->projs, theList->n,
					    &pt1, &pt2 );
    break;
  }

  medDirection[0] = (double)(pt2.x - pt1.x);
  medDirection[1] = (double)(pt2.y - pt1.y);
  medDirection[2] = (double)(pt2.z - pt1.z);

  n = medDirection[0]*medDirection[0] + medDirection[1]*medDirection[1] 
    + medDirection[2]*medDirection[2];

  if ( _verbose_ ) {
    printf("MED = (%9.5g %9.5g %9.5g) - (%9.5g %9.5g %9.5g) = %g\n",
	   pt1.x, pt1.y, pt1.z, pt2.x, pt2.y, pt2.z, n );
  }


  /* c'est un segment de droite
   */
  if ( n <= _EPSILON_FOR_NORMS_ ) {
    thePar->medDiameter = 0.0;
    medDirection[0] = medDirection[1] = medDirection[2] = 0.0;
    thePar->minDiameter = 0.0;
    minDirection[0] = minDirection[1] = minDirection[2] = 0.0;
    return;
  }

  n = sqrt( n );
  medDirection[0] /= n;
  medDirection[1] /= n;
  medDirection[2] /= n;

  thePar->medDiameter = n;


  if ( _verbose_ ) {
    printf( "   DIR = (%9.5g %9.5g %9.5g)\n", 
	    medDirection[0], medDirection[1], medDirection[2] );
  }


  /* 3eme direction
     c'est le produit vectoriel des deux autres
   */
  minDirection[0] = maxDirection[1] * medDirection[2]
                          - maxDirection[2] * medDirection[1];
  minDirection[1] = maxDirection[2] * medDirection[0]
                          - maxDirection[0] * medDirection[2];
  minDirection[2] = maxDirection[0] * medDirection[1]
                          - maxDirection[1] * medDirection[0];

  n = minDirection[0]*minDirection[0] + minDirection[1]*minDirection[1] 
    + minDirection[2]*minDirection[2];
  
  if ( n <= _EPSILON_FOR_NORMS_ ) {
    thePar->minDiameter = 0.0;
    minDirection[0] = minDirection[1] = minDirection[2] = 0.0;
    return;
  }

  n = sqrt( n );

  minDirection[0] /= n;
  minDirection[1] /= n;
  minDirection[2] /= n;

  if ( _verbose_ ) {
    printf( "MINDIR = (%9.5g %9.5g %9.5g)\n", 
	    minDirection[0], minDirection[1], minDirection[2] );
  }


  /* abscisses curvilignes
   */
  for ( i=0; i<theList->n; i++ )
    theList->abscisse[i] = theList->point[i].x * minDirection[0] 
      + theList->point[i].y * minDirection[1] 
      + theList->point[i].z * minDirection[2];
  
  thePar->minDiameter = _ComputeMaxDifferenceInList( theList->abscisse, 
    theList->n );

}








/* Pas de test en cas de depassement dans
   la liste
*/

static void _FillListWithBorderPoints( unsigned short int *theBuf,
				       const int *theDim,
				       typeListe *theList,
				       int color,
				       int *minCorner,
				       int *maxCorner )
{
  int i, x, y, z;

  int dz = theDim[0]*theDim[1];
  int dy = theDim[0];
  
  theList->n = 0;



  /* do we have to do tests ?
   */
  if ( minCorner[0] >= 1 && maxCorner[0] < theDim[0]-1 &&
       minCorner[1] >= 1 && maxCorner[1] < theDim[1]-1 &&
       minCorner[2] >= 1 && maxCorner[2] < theDim[2]-1 ) {
    
    for ( z=minCorner[2]; z<=maxCorner[2] && z<theDim[2]; z ++ ) {
      for ( y=minCorner[1]; y<=maxCorner[1] && y<theDim[1]; y ++ ) {
	i = z * theDim[0]*theDim[1] + y * theDim[0] + minCorner[0];
	for ( x=minCorner[0]; x<=maxCorner[0] && x<theDim[0]; x ++, i++ ) {
	  
	  if ( theBuf[i] != color ) continue;
	  if ( theBuf[i-dz] != color || theBuf[i+dz] != color || 
	       theBuf[i-dy] != color || theBuf[i+dy] != color ||
	       theBuf[i- 1] != color || theBuf[i+ 1] != color ) {
	    theList->point[theList->n].x = x;
	    theList->point[theList->n].y = y;
	    theList->point[theList->n].z = z;
	    theList->n ++;
	  }

	}
      }
    }
  } 

  /* here we have do to tests
   */
  else {

    for ( z=minCorner[2]; z<=maxCorner[2] && z<theDim[2]; z ++ ) {
      for ( y=minCorner[1]; y<=maxCorner[1] && y<theDim[1]; y ++ ) {
	i = z * theDim[0]*theDim[1] + y * theDim[0] + minCorner[0];
	for ( x=minCorner[0]; x<=maxCorner[0] && x<theDim[0]; x ++, i++ ) {
	  
	  if ( theBuf[i] != color ) continue;
	  if ( x-1 >= 0 ) {
	    if ( theBuf[i- 1] != color ) {
	      theList->point[theList->n].x = x;
	      theList->point[theList->n].y = y;
	      theList->point[theList->n].z = z;
	      theList->n ++;
	      continue;
	    }
	  }
	  if ( x+1 < theDim[0] ) {
	    if ( theBuf[i+ 1] != color ) {
	      theList->point[theList->n].x = x;
	      theList->point[theList->n].y = y;
	      theList->point[theList->n].z = z;
	      theList->n ++;
	      continue;
	    }
	  }
	  if ( y-1 >= 0 ) {
	    if ( theBuf[i-dy] != color ) {
	      theList->point[theList->n].x = x;
	      theList->point[theList->n].y = y;
	      theList->point[theList->n].z = z;
	      theList->n ++;
	      continue;
	    }
	  }
	  if ( y+1 < theDim[1] ) {
	    if ( theBuf[i+dy] != color ) {
	      theList->point[theList->n].x = x;
	      theList->point[theList->n].y = y;
	      theList->point[theList->n].z = z;
	      theList->n ++;
	      continue;
	    }
	  }
	  if ( z-1 >= 0 ) {
	    if ( theBuf[i-dz] != color ) {
	      theList->point[theList->n].x = x;
	      theList->point[theList->n].y = y;
	      theList->point[theList->n].z = z;
	      theList->n ++;
	      continue;
	    }
	  }
	  if ( z+1 < theDim[2] ) {
	    if ( theBuf[i+dz] != color ) {
	      theList->point[theList->n].x = x;
	      theList->point[theList->n].y = y;
	      theList->point[theList->n].z = z;
	      theList->n ++;
	      continue;
	    }
	  }

	}
      }
    }
  }

  return;
  
}















#ifdef _UNUSED_
static void _PrintTypeListe( typeListe *theListe )
{
  int i;
  for ( i=0; i<theListe->n; i++ ) {
    printf("%5d : %5.2f %5.2f %5.2f\n", i, theListe->point[i].x, 
	   theListe->point[i].y, theListe->point[i].z );
  }
}
#endif






static void _InitTypeListe( typeListe *theListe )
{
  theListe->n = 0;
  theListe->nbAllocs = 0;
  theListe->point = (typePoint *)NULL;
  theListe->projs = (typePoint *)NULL;
  theListe->abscisse = (double*)NULL;
}








static void _FreeTypeListe( typeListe *theListe )
{
  if ( theListe->nbAllocs > 0 ) {
    if ( theListe->point != (typePoint *)NULL ) 
      free( theListe->point );
    if ( theListe->projs != (typePoint *)NULL ) 
      free( theListe->projs );
    if ( theListe->abscisse != (double*)NULL )
      free( theListe->abscisse );
  }
  theListe->n = 0;
  theListe->nbAllocs = 0;
  theListe->point = (typePoint *)NULL;
  theListe->projs = (typePoint *)NULL;
  theListe->abscisse = (double*)NULL;
}







static int _AllocTypeListe( typeListe *theListe, int n )
{
  char *proc = "_AllocTypeListe";
  typePoint *p = (typePoint *)NULL;
  double *d = (double*)NULL;

  
  if ( n <= theListe->nbAllocs ) return( 1 );
  if ( n <= 0 )  return( 1 );


  p = (typePoint *)malloc( n * sizeof( typePoint ) );
  if ( p == (typePoint *)NULL ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: can not allocate 1st new points list (size=%d).\n", 
	       proc, n );
    }
    return( 0 );
  }
  if ( theListe->nbAllocs > 0 && theListe->point != (typePoint *)NULL ) 
    free( theListe->point );
  
  theListe->point = p;


  
  p = (typePoint *)malloc( n * sizeof( typePoint ) );
  if ( p == (typePoint *)NULL ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: can not allocate 2nd new points list (size=%d).\n", 
	       proc, n );
    }
    return( 0 );
  }
  if ( theListe->nbAllocs > 0 && theListe->projs != (typePoint *)NULL ) 
    free( theListe->projs );
  
  theListe->projs = p;



  d = (double*)malloc( n * sizeof( double ) );
  if ( d == (double*)NULL ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: can not allocate 3rd new points list (size=%d).\n", 
	       proc, n );
    }
    return( 0 );
  }
  if ( theListe->nbAllocs > 0 && theListe->abscisse != (double *)NULL ) 
    free( theListe->abscisse );
  
  theListe->abscisse = d;

  theListe->nbAllocs = n;
  return( 1 );
}






static void _ApplyVoxelDimensionsToList( typeListe *theList,
					 const double *sizes )
{
  int i;
  
  if ( sizes == NULL ) return;
  for ( i=0; i<theList->n; i++ ) {
    theList->point[i].x *= sizes[0];
    theList->point[i].y *= sizes[1];
    theList->point[i].z *= sizes[2];
  }
  
}
