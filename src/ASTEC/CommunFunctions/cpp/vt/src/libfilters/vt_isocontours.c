#include <math.h>

#include <vt_isocontours.h>
#include <vt_extract.h>


/* les directions d'entrees sorties
   in a cube  (x,y)   - (x+1,y)
                |           |
              (x,y+1) - (x+1,y+1)
   XINF -> edge = (x,y)   - (x+1,y)
   XSUP -> edge = (x,y+1) - (x+1,y+1)
   YINF -> edge = (x,y)   - (x,y+1)
   YSUP -> edge = (x+1,y) - (x+1,y+1)
 */
#define XINF 0
#define XSUP 1
#define YINF 2
#define YSUP 3



/*
  1 - 2 
  |   |
  8 - 4
 
  0000 1000 0100 1100
  0010 1010 0110 1110
  0001 1001 0101 1101
  0011 1011 0111 1111
*/
  
static int nbcontours[16] = { 0, 1, 1, 1,
			      1, 2, 1, 1,
			      1, 1, 2, 1,
			      1, 1, 1, 0 };



typeListOfContours3D *VT_Compute2DIsoContoursInXYSlice( vt_image *theIm, 
					     int z, /* slice to be processed */
					     double threshold )
{
  char *proc = "VT_Compute2DIsoContoursInXYSlice";
  typeSlice *theSlice;
  typeListOfContours3D *theCnt3D = NULL;
  typeContour3D *c;
  int i, n;

  theSlice = VT_Compute2DIsoContours( theIm, z, threshold );
  if ( theSlice == NULL ) {
    return( NULL ) ;
  }

  theCnt3D = (typeListOfContours3D*)malloc( sizeof(typeListOfContours3D) );
  if ( theCnt3D == NULL ) {
    fprintf( stderr, "%s: allocation failed\n", proc );
    freeSlice( theSlice );
    free( theSlice );
    return( NULL );
  }
  initListOfContours3D( theCnt3D );

  for ( i=0; i<theSlice->n; i++ ) {
    c = (typeContour3D *)malloc( sizeof( typeContour3D ) );
    initContour3D( c );
    allocContour3D( c, theSlice->theContours[i]->n );
    for ( n=0; n<theSlice->theContours[i]->n; n++ ) {
      c->thePts[n].x = theSlice->theContours[i]->thePts[n].x;
      c->thePts[n].y = theSlice->theContours[i]->thePts[n].y;
      c->thePts[n].z = z;
    }
    c->n = theSlice->theContours[i]->n;
    c->topology = theSlice->theContours[i]->topology;
    addContour3DToListOfContours3D( c, theCnt3D );
  }

  freeSlice( theSlice );
  free( theSlice );
  return( theCnt3D );
}


typeListOfContours3D *VT_Compute2DIsoContoursInXZSlice( vt_image *theIm, 
					     int y, /* slice to be processed */
					     double threshold )
{
  char *proc = "VT_Compute2DIsoContoursInXZSlice";
  vt_image extIm;
  vt_ipt extCorner;
  vt_image tmpIm;
  typeListOfContours3D *theCnt3D = NULL;
  typeSlice *theSlice;
  typeContour3D *c;
  int i, n;

  extCorner.x = 0;
  extCorner.y = y;
  extCorner.z = 0;

  VT_InitImage( &extIm, NULL, theIm->dim.x, 1, theIm->dim.z, theIm->type );
  if ( VT_AllocImage( &extIm ) != 1 ) {
    fprintf( stderr, "%s: allocation of auxiliary image #1 failed\n", proc );
    return( NULL );
  }

  VT_InitImage( &tmpIm, NULL, theIm->dim.x, theIm->dim.z, 1, theIm->type );
  if ( VT_AllocImage( &tmpIm ) != 1 ) {
    fprintf( stderr, "%s: allocation of auxiliary image #2 failed\n", proc );
    VT_FreeImage( &extIm );
    return( NULL );
  }

  if ( VT_Extract( &extIm, theIm, &extCorner ) != 1 ) {
    fprintf( stderr, "%s: extraction of auxiliary image #1 failed\n", proc );
    VT_FreeImage( &tmpIm );
    VT_FreeImage( &extIm );
    return( 0 );
  }

  memcpy( tmpIm.buf, extIm.buf, VT_SizeImage( &tmpIm ) );
  VT_FreeImage( &extIm );

  theSlice = VT_Compute2DIsoContours( &tmpIm, 0, threshold );
  VT_FreeImage( &tmpIm );

  if ( theSlice == NULL ) {
    return( NULL ) ;
  }

  theCnt3D = (typeListOfContours3D*)malloc( sizeof(typeListOfContours3D) );
  if ( theCnt3D == NULL ) {
    fprintf( stderr, "%s: allocation failed\n", proc );
    freeSlice( theSlice );
    free( theSlice );
    return( NULL );
  }
  initListOfContours3D( theCnt3D );

  for ( i=0; i<theSlice->n; i++ ) {
    c = (typeContour3D *)malloc( sizeof( typeContour3D ) );
    initContour3D( c );
    allocContour3D( c, theSlice->theContours[i]->n );
    for ( n=0; n<theSlice->theContours[i]->n; n++ ) {
      c->thePts[n].x = theSlice->theContours[i]->thePts[n].x;
      c->thePts[n].y = y;
      c->thePts[n].z = theSlice->theContours[i]->thePts[n].y;
    }
    c->n = theSlice->theContours[i]->n;
    c->topology = theSlice->theContours[i]->topology;
    addContour3DToListOfContours3D( c, theCnt3D );
  }

  freeSlice( theSlice );
  free( theSlice );
  return( theCnt3D );
}


typeListOfContours3D *VT_Compute2DIsoContoursInYZSlice( vt_image *theIm, 
					     int x, /* slice to be processed */
					     double threshold )
{
  char *proc = "VT_Compute2DIsoContoursInXZSlice";
  vt_image extIm;
  vt_ipt extCorner;
  vt_image tmpIm;
  typeListOfContours3D *theCnt3D = NULL;
  typeSlice *theSlice;
  typeContour3D *c;
  int i, n;

  extCorner.x = x;
  extCorner.y = 0;
  extCorner.z = 0;

  VT_InitImage( &extIm, NULL, 1, theIm->dim.y, theIm->dim.z, theIm->type );
  if ( VT_AllocImage( &extIm ) != 1 ) {
    fprintf( stderr, "%s: allocation of auxiliary image #1 failed\n", proc );
    return( NULL );
  }

  VT_InitImage( &tmpIm, "tmp.inr", theIm->dim.y, theIm->dim.z, 1, theIm->type );
  if ( VT_AllocImage( &tmpIm ) != 1 ) {
    fprintf( stderr, "%s: allocation of auxiliary image #2 failed\n", proc );
    VT_FreeImage( &extIm );
    return( NULL );
  }

  if ( VT_Extract( &extIm, theIm, &extCorner ) != 1 ) {
    fprintf( stderr, "%s: extraction of auxiliary image #1 failed\n", proc );
    VT_FreeImage( &tmpIm );
    VT_FreeImage( &extIm );
    return( 0 );
  }

  memcpy( tmpIm.buf, extIm.buf, VT_SizeImage( &tmpIm ) );
  VT_FreeImage( &extIm );

  if ( 0 ) VT_WriteInrimage( &tmpIm );

  theSlice = VT_Compute2DIsoContours( &tmpIm, 0, threshold );
  VT_FreeImage( &tmpIm );

  if ( theSlice == NULL ) {
    return( NULL ) ;
  }

  theCnt3D = (typeListOfContours3D*)malloc( sizeof(typeListOfContours3D) );
  if ( theCnt3D == NULL ) {
    fprintf( stderr, "%s: allocation failed\n", proc );
    freeSlice( theSlice );
    free( theSlice );
    return( NULL );
  }
  initListOfContours3D( theCnt3D );

  for ( i=0; i<theSlice->n; i++ ) {
    c = (typeContour3D *)malloc( sizeof( typeContour3D ) );
    initContour3D( c );
    allocContour3D( c, theSlice->theContours[i]->n );
    for ( n=0; n<theSlice->theContours[i]->n; n++ ) {
      c->thePts[n].x = x;
      c->thePts[n].y = theSlice->theContours[i]->thePts[n].x;
      c->thePts[n].z = theSlice->theContours[i]->thePts[n].y;
    }
    c->n = theSlice->theContours[i]->n;
    c->topology = theSlice->theContours[i]->topology;
    addContour3DToListOfContours3D( c, theCnt3D );
  }

  freeSlice( theSlice );
  free( theSlice );
  return( theCnt3D );
}


typeSlice *VT_Compute2DIsoContours( vt_image *theIm, 
				    int z, /* slice to be processed */
				    double threshold )
{
  int x, y;
  vt_image immarks;
  u8 *** theMarks;
  char *proc = "VT_Compute2DIsoContours";

  typeSlice *slice = NULL;
  
  VT_InitFromImage( &immarks, theIm, "marks.inr", UCHAR );
  immarks.dim.z = 1;
  if ( VT_AllocImage( &immarks ) != 1 ) {
    return( NULL );
  }
  memset( immarks.buf, 3, theIm->dim.x * theIm->dim.y );
  theMarks = (u8***)immarks.array;
  
  /* the x edge from (x,y) to (x+1,y) is encoded as (& 1) in [x,y]
     the y edge from (x,y) to (x,y+1) is encoded as (& 2) in [x,y]
  */

  slice = (typeSlice *)malloc( sizeof(typeSlice) );
  if ( slice == NULL ) {
    VT_FreeImage( &immarks );
    return( NULL );
  }
  initSlice( slice );

  switch( theIm->type ) {

  default :  
    fprintf( stderr, "%s: type not handled in switch\n", proc );
    free( slice );
    VT_FreeImage( &immarks );
    return( NULL );

  case UCHAR :
    {
      u8***theBuf = (u8***)theIm->array;
      /* look for a first point
       */
      for ( y=0; y<(int)theIm->dim.y - 1; y++ ) {
	
  for ( x=0; x<(int)theIm->dim.x - 1; x++ ) {
	  if ( theBuf[z][y][x] < threshold ) {
	    if ( theBuf[z][y][x+1] >= threshold 
		 && theMarks[0][y][x] & (unsigned char)1 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, XINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	    }
	    if ( theBuf[z][y+1][x] >= threshold 
		 && theMarks[0][y][x] & (unsigned char)2 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, YINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	    }
	  }
	  else {
	    if ( theBuf[z][y][x+1] < threshold 
		 && theMarks[0][y][x] & (unsigned char)1 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, XINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	    }
	    if ( theBuf[z][y+1][x] < threshold 
		 && theMarks[0][y][x] & (unsigned char)2 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, YINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	    }
	  }
	}

	/* dernier point de la ligne
	 */
	x = theIm->dim.x - 1;
	if ( theBuf[z][y][x] < threshold ) {
	  if ( theBuf[z][y+1][x] >= threshold 
	       && theMarks[0][y][x] & (unsigned char)2 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, YINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	  }
	}
	else {
	  if ( theBuf[z][y+1][x] < threshold 
	       && theMarks[0][y][x] & (unsigned char)2 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, YINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	  }
	}

      } /* fin du traitement des lignes */

      /* arete horizontale pour la derniere ligne
       */
      y = theIm->dim.y - 1;
      for ( x=0; x<(int)theIm->dim.x - 1; x++ ) {
	if ( theBuf[z][y][x] < threshold ) {
	  if ( theBuf[z][y][x+1] >= threshold 
	       && theMarks[0][y][x] & (unsigned char)1 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, XINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	  }
	}
	else {
	  if ( theBuf[z][y][x+1] < threshold 
	       && theMarks[0][y][x] & (unsigned char)1 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, XINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	  }
	}
      }
    }
    /* end case UCHAR : */
    break;

  case USHORT :
    {
      u16***theBuf = (u16***)theIm->array;
      /* look for a first point
       */
      for ( y=0; y<(int)theIm->dim.y - 1; y++ ) {
	
  for ( x=0; x<(int)theIm->dim.x - 1; x++ ) {
	  if ( theBuf[z][y][x] < threshold ) {
	    if ( theBuf[z][y][x+1] >= threshold 
		 && theMarks[0][y][x] & (unsigned char)1 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, XINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	    }
	    if ( theBuf[z][y+1][x] >= threshold 
		 && theMarks[0][y][x] & (unsigned char)2 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, YINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	    }
	  }
	  else {
	    if ( theBuf[z][y][x+1] < threshold 
		 && theMarks[0][y][x] & (unsigned char)1 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, XINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	    }
	    if ( theBuf[z][y+1][x] < threshold 
		 && theMarks[0][y][x] & (unsigned char)2 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, YINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	    }
	  }
	}

	/* dernier point de la ligne
	 */
	x = theIm->dim.x - 1;
	if ( theBuf[z][y][x] < threshold ) {
	  if ( theBuf[z][y+1][x] >= threshold 
	       && theMarks[0][y][x] & (unsigned char)2 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, YINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	  }
	}
	else {
	  if ( theBuf[z][y+1][x] < threshold 
	       && theMarks[0][y][x] & (unsigned char)2 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, YINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	  }
	}

      } /* fin du traitement des lignes */

      /* arete horizontale pour la derniere ligne
       */
      y = theIm->dim.y - 1;
      for ( x=0; x<(int)theIm->dim.x - 1; x++ ) {
	if ( theBuf[z][y][x] < threshold ) {
	  if ( theBuf[z][y][x+1] >= threshold 
	       && theMarks[0][y][x] & (unsigned char)1 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, XINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	  }
	}
	else {
	  if ( theBuf[z][y][x+1] < threshold 
	       && theMarks[0][y][x] & (unsigned char)1 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, XINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	  }
	}
      }
    }
    /* end case USHORT : */
    break;

  case SSHORT :
    {
      s16***theBuf = (s16***)theIm->array;
      /* look for a first point
       */
      for ( y=0; y<(int)theIm->dim.y - 1; y++ ) {
	
  for ( x=0; x<(int)theIm->dim.x - 1; x++ ) {
	  if ( theBuf[z][y][x] < threshold ) {
	    if ( theBuf[z][y][x+1] >= threshold 
		 && theMarks[0][y][x] & (unsigned char)1 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, XINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	    }
	    if ( theBuf[z][y+1][x] >= threshold 
		 && theMarks[0][y][x] & (unsigned char)2 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, YINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	    }
	  }
	  else {
	    if ( theBuf[z][y][x+1] < threshold 
		 && theMarks[0][y][x] & (unsigned char)1 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, XINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	    }
	    if ( theBuf[z][y+1][x] < threshold 
		 && theMarks[0][y][x] & (unsigned char)2 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, YINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	    }
	  }
	}

	/* dernier point de la ligne
	 */
	x = theIm->dim.x - 1;
	if ( theBuf[z][y][x] < threshold ) {
	  if ( theBuf[z][y+1][x] >= threshold 
	       && theMarks[0][y][x] & (unsigned char)2 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, YINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	  }
	}
	else {
	  if ( theBuf[z][y+1][x] < threshold 
	       && theMarks[0][y][x] & (unsigned char)2 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, YINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	  }
	}

      } /* fin du traitement des lignes */

      /* arete horizontale pour la derniere ligne
       */
      y = theIm->dim.y - 1;
      for ( x=0; x<(int)theIm->dim.x - 1; x++ ) {
	if ( theBuf[z][y][x] < threshold ) {
	  if ( theBuf[z][y][x+1] >= threshold 
	       && theMarks[0][y][x] & (unsigned char)1 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, XINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	  }
	}
	else {
	  if ( theBuf[z][y][x+1] < threshold 
	       && theMarks[0][y][x] & (unsigned char)1 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, XINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	  }
	}
      }
    }
    /* end case SSHORT : */
    break;

   case FLOAT :
    {
      float***theBuf = (float***)theIm->array;
      /* look for a first point
       */
      for ( y=0; y<(int)theIm->dim.y - 1; y++ ) {
	
  for ( x=0; x<(int)theIm->dim.x - 1; x++ ) {
	  if ( theBuf[z][y][x] < threshold ) {
	    if ( theBuf[z][y][x+1] >= threshold 
		 && theMarks[0][y][x] & (unsigned char)1 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, XINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	    }
	    if ( theBuf[z][y+1][x] >= threshold 
		 && theMarks[0][y][x] & (unsigned char)2 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, YINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	    }
	  }
	  else {
	    if ( theBuf[z][y][x+1] < threshold 
		 && theMarks[0][y][x] & (unsigned char)1 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, XINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	    }
	    if ( theBuf[z][y+1][x] < threshold 
		 && theMarks[0][y][x] & (unsigned char)2 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, YINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	    }
	  }
	}

	/* dernier point de la ligne
	 */
	x = theIm->dim.x - 1;
	if ( theBuf[z][y][x] < threshold ) {
	  if ( theBuf[z][y+1][x] >= threshold 
	       && theMarks[0][y][x] & (unsigned char)2 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, YINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	  }
	}
	else {
	  if ( theBuf[z][y+1][x] < threshold 
	       && theMarks[0][y][x] & (unsigned char)2 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, YINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	  }
	}

      } /* fin du traitement des lignes */

      /* arete horizontale pour la derniere ligne
       */
      y = theIm->dim.y - 1;
      for ( x=0; x<(int)theIm->dim.x - 1; x++ ) {
	if ( theBuf[z][y][x] < threshold ) {
	  if ( theBuf[z][y][x+1] >= threshold 
	       && theMarks[0][y][x] & (unsigned char)1 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, XINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	  }
	}
	else {
	  if ( theBuf[z][y][x+1] < threshold 
	       && theMarks[0][y][x] & (unsigned char)1 ) {
	      if ( VT_ExtractContour2DAndAddToSlice( theIm, &immarks, x, y, z,
						     threshold, XINF, slice ) != 1 ) {
		free( slice );
		VT_FreeImage( &immarks );
		return( NULL );
	      }
	  }
	}
      }
    }
    /* end case FLOAT : */
    break;
 }


  if ( 0 ) VT_WriteInrimage( &immarks );
  VT_FreeImage( &immarks );
  
  slice->z = z;

  return( slice );
}



int VT_ExtractContour2DAndAddToSlice( vt_image *theIm, 
				      vt_image *theMarks,
				      int startx, int starty, int startz,
				      double threshold, 
				      int EDGE,
				      typeSlice *slice )
{
  typeContour2D *c;

  c = VT_ExtractContour2D( theIm, theMarks, startx, starty, startz, threshold, EDGE );
  if ( c == NULL ) return( 0 );
  if( addContour2DToSlice( c, slice ) != 1 ) {
    freeContour2D( c );
    free( c );
    return( 0 );
  }
  return( 1 );
}



typeContour2D *VT_ExtractContour2D( vt_image *theIm, 
				    vt_image *theMarks,
				    int startx, int starty, int startz,
				    double threshold, 
				    int EDGE )
{
  char *proc = "VT_ExtractContour";
  typeContour2D *resCont, *theCont, theCont1, theCont2;
  int ix = startx;
  int iy = starty;
  const int iz = startz;
  double x, y;
  u8 *** bufMarks = (u8***)theMarks->array;

  int nparcours, nmaxparcours;
  int finparcours;
  int CEDGE, NEDGE;
  double middle, cube[2][2];
  int configuration;
  int i, n;

  initContour2D( &theCont1 );
  initContour2D( &theCont2 );

  resCont = (typeContour2D*)malloc( sizeof( typeContour2D ) );
  if ( resCont == NULL ) return( NULL );
  initContour2D( resCont );

  
  /* compute first point 
   */
  x = ix;
  y = iy;
  
  switch ( theIm->type ) {

  default : 
    free( resCont );
    return( NULL );

  case UCHAR : 
    {
      u8***theBuf = (u8***)theIm->array;
      switch ( EDGE ) {
      case XINF :
	x = ix + ( (double)threshold - (double)theBuf[iz][iy][ix] ) /
	  ( (double)theBuf[iz][iy][ix+1] - (double)theBuf[iz][iy][ix] );
	break;
      case YINF :
	y = iy + ( (double)threshold - (double)theBuf[iz][iy][ix] ) /
	  ( (double)theBuf[iz][iy+1][ix] - (double)theBuf[iz][iy][ix] );
	break;
      }
    }
    break;
  case USHORT : 
    {
      u16***theBuf = (u16***)theIm->array;
      switch ( EDGE ) {
      case XINF :
	x = ix + ( (double)threshold - (double)theBuf[iz][iy][ix] ) /
	  ( (double)theBuf[iz][iy][ix+1] - (double)theBuf[iz][iy][ix] );
	break;
      case YINF :
	y = iy + ( (double)threshold - (double)theBuf[iz][iy][ix] ) /
	  ( (double)theBuf[iz][iy+1][ix] - (double)theBuf[iz][iy][ix] );
	break;
      }
    }
    break;
  case SSHORT : 
    {
      s16***theBuf = (s16***)theIm->array;
      switch ( EDGE ) {
      case XINF :
	x = ix + ( (double)threshold - (double)theBuf[iz][iy][ix] ) /
	  ( (double)theBuf[iz][iy][ix+1] - (double)theBuf[iz][iy][ix] );
	break;
      case YINF :
	y = iy + ( (double)threshold - (double)theBuf[iz][iy][ix] ) /
	  ( (double)theBuf[iz][iy+1][ix] - (double)theBuf[iz][iy][ix] );
	break;
      }
    }
    break;
   case FLOAT : 
    {
      float***theBuf = (float***)theIm->array;
      switch ( EDGE ) {
      case XINF :
	x = ix + ( (double)threshold - (double)theBuf[iz][iy][ix] ) /
	  ( (double)theBuf[iz][iy][ix+1] - (double)theBuf[iz][iy][ix] );
	break;
      case YINF :
	y = iy + ( (double)threshold - (double)theBuf[iz][iy][ix] ) /
	  ( (double)theBuf[iz][iy+1][ix] - (double)theBuf[iz][iy][ix] );
	break;
      }
    }
    break;
  }

  if ( addPointToContour2D( x, y, &theCont1 ) != 1 ) {
    free( resCont );
    return( NULL );
  }

  switch ( EDGE ) {
  case XINF : bufMarks[0][iy][ix] &= ~(unsigned char)1; break;
  case YINF : bufMarks[0][iy][ix] &= ~(unsigned char)2; break;
  }

  if ( 0 ) {
    printf( " point #1 : " );
    printf( "%4d %4d %4d", ix, iy, iz );
    switch( EDGE ) {
    case XINF : printf( " <-> %4d %4d %4d", ix+1, iy, iz ); break;
    case YINF : printf( " <-> %4d %4d %4d", ix, iy+1, iz ); break;
    }
    printf( " = %f %f %d\n", x ,y, iz );
  }


  /*
    les departs d'aretes horizontales ne peuvent avoir lieu que sur
    la premiere ligne, sinon c'est forcement vertical
  */
  switch( EDGE ) {
  case XINF :
    if ( starty != 0 ) {
      fprintf( stderr, "probleme in %s\n", proc );
    }
  }
  

  /* pour les aretes horizontales (a priori seulement la premiere ligne)
     un seul parcours suffit
     pour les autres (les aretes verticales), il faut prevoir 2 parcours
     pour les courbes ouvertes
  */
  
  /* nombre max de pacours 
   */
  switch( EDGE ) {
  default :
  case XINF : nmaxparcours = 1; break;
  case YINF :
    if ( startx > 0 && startx < (int)theIm->dim.x-1 ) nmaxparcours = 2;
    else nmaxparcours = 1; 
    break;
  }

  for ( nparcours = 0; nparcours < nmaxparcours; nparcours ++ ) {
    

    switch( nparcours ) {
    default : theCont = NULL; break;
    case 0  : theCont = &theCont1; break;
    case 1  : theCont = &theCont2; break;
    }
    
    /* cube to be processed
       CEDGE is the edge from which the first point has be extracted
       (ix,iy,iz) is the left upper corner

       les directions d'entrees sorties
       in a cube  (x,y)   - (x+1,y)
                    |           |
                  (x,y+1) - (x+1,y+1)
       XINF -> edge = (x,y)   - (x+1,y)
       XSUP -> edge = (x,y+1) - (x+1,y+1)
       YINF -> edge = (x,y)   - (x,y+1)
       YSUP -> edge = (x+1,y) - (x+1,y+1)
    */

    CEDGE = EDGE;
    switch( EDGE ) {
    default :
    case XINF : 
      break;
    case  YINF :
      if ( (startx == (int)theIm->dim.x-1) || nparcours == 1 ) {
	ix = startx-1;
	iy = starty;
	CEDGE = YSUP;
      }
    }
    
    /*
    printf( " parcours #%d from (%d %d %d)\n ", nparcours+1, ix, iy, iz );
    */
    
    finparcours = 0;
    do {

      /* extract configuration
       */
      switch ( theIm->type ) {
      default : free( resCont ); return( NULL );
      case UCHAR :
	{
	  u8***theBuf = (u8***)theIm->array;
	  cube[0][0] = (double)theBuf[iz][iy  ][ix];
	  cube[0][1] = (double)theBuf[iz][iy  ][ix+1];
	  cube[1][0] = (double)theBuf[iz][iy+1][ix];
	  cube[1][1] = (double)theBuf[iz][iy+1][ix+1];
	}
	break;
      case USHORT :
	{
	  u16***theBuf = (u16***)theIm->array;
	  cube[0][0] = (double)theBuf[iz][iy  ][ix];
	  cube[0][1] = (double)theBuf[iz][iy  ][ix+1];
	  cube[1][0] = (double)theBuf[iz][iy+1][ix];
	  cube[1][1] = (double)theBuf[iz][iy+1][ix+1];
	}
	break;
      case SSHORT :
	{
	  s16***theBuf = (s16***)theIm->array;
	  cube[0][0] = (double)theBuf[iz][iy  ][ix];
	  cube[0][1] = (double)theBuf[iz][iy  ][ix+1];
	  cube[1][0] = (double)theBuf[iz][iy+1][ix];
	  cube[1][1] = (double)theBuf[iz][iy+1][ix+1];
	}
	break;
       case FLOAT :
	{
	  float***theBuf = (float***)theIm->array;
	  cube[0][0] = (double)theBuf[iz][iy  ][ix];
	  cube[0][1] = (double)theBuf[iz][iy  ][ix+1];
	  cube[1][0] = (double)theBuf[iz][iy+1][ix];
	  cube[1][1] = (double)theBuf[iz][iy+1][ix+1];
	}
	break;
      }
      configuration = 0;
      if ( cube[0][0] >= threshold )  configuration |= 1;
      if ( cube[0][1] >= threshold )  configuration |= 2;
      if ( cube[1][1] >= threshold )  configuration |= 4;
      if ( cube[1][0] >= threshold )  configuration |= 8;

      NEDGE = -1;

      /*
      if ( 1 ) {
	fprintf( stdout, "#%d/%d configuration = %2d", 
		 nparcours+1, nmaxparcours, configuration );
	fprintf( stdout, " (%3d %3d %3d)", ix, iy, iz );
	fprintf( stdout, " arete courante =" );
	switch( CEDGE ) {
	default : break;
	case XINF : fprintf( stdout, "XINF" ); break;
	case XSUP : fprintf( stdout, "XSUP" ); break;
	case YINF : fprintf( stdout, "YINF" ); break;
	case YSUP : fprintf( stdout, "YSUP" ); break;
	}
	fprintf( stdout, "\n" );
	fprintf( stdout, "      %9f", cube[0][0] );
	if ( CEDGE == XINF ) fprintf( stdout, " === " );
	else                 fprintf( stdout, " --- " );
	fprintf( stdout, "%9f\n", cube[0][1] );
	if ( CEDGE == YINF )    fprintf( stdout, "          ||            |\n" );
	else if (CEDGE == YSUP) fprintf( stdout, "          |             ||\n" );
        else                    fprintf( stdout, "          |             |\n" );
	if ( CEDGE == YINF )    fprintf( stdout, "          ||            |\n" );
	else if (CEDGE == YSUP) fprintf( stdout, "          |             ||\n" );
        else                    fprintf( stdout, "          |             |\n" );
	fprintf( stdout, "      %9f", cube[1][0] );
	if ( CEDGE == XSUP ) fprintf( stdout, " === " );
	else                 fprintf( stdout, " --- " );
	fprintf( stdout, "%9f\n", cube[1][1] );
      }
      */


      switch ( configuration ) {

      case 0 :
      case 15 :
	fprintf( stderr, "%s: no contours in this cube ?\n", proc );
	fprintf( stdout, "      config = %d, threshold = %f\n",
		 configuration, threshold );
	if ( 0 ) {
	  fprintf( stdout, "      %9f", cube[0][0] );
	  if ( CEDGE == XINF ) fprintf( stdout, " === " );
	  else                 fprintf( stdout, " --- " );
	  fprintf( stdout, "%9f\n", cube[0][1] );
	  if ( CEDGE == YINF ) fprintf( stdout, "          ||            |\n" );
	  else                 fprintf( stdout, "          |             ||\n" );
	  if ( CEDGE == YINF ) fprintf( stdout, "          ||            |\n" );
	  else                 fprintf( stdout, "          |             ||\n" );
	  fprintf( stdout, "      %9f", cube[1][0] );
	  if ( CEDGE == XSUP ) fprintf( stdout, " === " );
	  else                 fprintf( stdout, " --- " );
	  fprintf( stdout, "%9f\n", cube[1][1] );
	}
	free( resCont );
	return( NULL );

      case 1 :
      case 14 :
	switch ( CEDGE ) {
	case XSUP :
	case YSUP :
	  fprintf( stderr, "%s: bad configuration, case 1/14\n", proc );
	  free( resCont );
	  return( NULL );
	case XINF :
	  NEDGE = YINF; break;
	case YINF :
	  NEDGE = XINF; break;
	}
	break;

      case 2 :
      case 13 :
	switch ( CEDGE ) {
	case XSUP :
	case YINF :
	  fprintf( stderr, "%s: bad configuration, case 1/13\n", proc );
	  free( resCont );
	  return( NULL );
	case XINF :
	  NEDGE = YSUP; break;
	case YSUP :
	  NEDGE = XINF; break;
	}
	break;

      case 3 :
      case 12 :
	switch ( CEDGE ) {
	case XSUP :
	case XINF :
	  fprintf( stderr, "%s: bad configuration, case 3/12\n", proc );
	  free( resCont );
	  return( NULL );
	case YINF :
	  NEDGE = YSUP; break;
	case YSUP :
	  NEDGE = YINF; break;
	}
	break;

      case 4 :
      case 11 :
	switch ( CEDGE ) {
	case XINF :
	case YINF :
	  fprintf( stderr, "%s: bad configuration, case 4/11\n", proc );
	  free( resCont );
	  return( NULL );
	case XSUP :
	  NEDGE = YSUP; break;
	case YSUP :
	  NEDGE = XSUP; break;
	}
	break;

      case 6 :
      case 9 :
	switch ( CEDGE ) {
	case YSUP :
	case YINF :
	  fprintf( stderr, "%s: bad configuration, case 6/9\n", proc );
	  free( resCont );
	  return( NULL );
	case XSUP :
	  NEDGE = XINF; break;
	case XINF :
	  NEDGE = XSUP; break;
	}
	break;
	
      case 7 :
      case 8 :
	switch ( CEDGE ) {
	case YSUP :
	case XINF :
	  fprintf( stderr, "%s: bad configuration, case 6/9\n", proc );
	  free( resCont );
	  return( NULL );
	case XSUP :
	  NEDGE = YINF; break;
	case YINF :
	  NEDGE = XSUP; break;
	}
	break;
	

      case 5 :
      case 10 :
	middle = ( cube[0][0] + cube[0][1] + cube[1][0] + cube[1][1] ) / 4.0;
	switch ( CEDGE ) {
	case XINF :
	  if ( cube[0][0] >= threshold ) {
	    NEDGE = ( middle >= threshold ) ? YSUP : YINF;
	  } 
	  else {
	    NEDGE = ( middle >= threshold ) ? YINF : YSUP;
	  }
	  break;
	case YINF :
	  if ( cube[0][0] >= threshold ) {
	    NEDGE = ( middle >= threshold ) ? XSUP : XINF;
	  } 
	  else {
	    NEDGE = ( middle >= threshold ) ? XINF : XSUP;
	  }
	  break;
	case XSUP :
	  if ( cube[0][0] >= threshold ) {
	    NEDGE = ( middle >= threshold ) ? YINF : YSUP;
	  } 
	  else {
	    NEDGE = ( middle >= threshold ) ? YSUP : YINF;
	  }
	  break;
	case YSUP :
	  if ( cube[0][0] >= threshold ) {
	    NEDGE = ( middle >= threshold ) ? XINF : XSUP;
	  } 
	  else {
	    NEDGE = ( middle >= threshold ) ? XSUP : XINF;
	  }
	  break;
	}

      }

      /*
      if ( starty == 154 && iz == 60 ) {
	fprintf( stdout, "                                     " );
	fprintf( stdout, " arete suivante =" );
	switch( NEDGE ) {
	default : break;
	case XINF : fprintf( stdout, "XINF" ); break;
	case XSUP : fprintf( stdout, "XSUP" ); break;
	case YINF : fprintf( stdout, "YINF" ); break;
	case YSUP : fprintf( stdout, "YSUP" ); break;
	}
	fprintf( stdout, "\n" );
      }
      */
      /* on connait l'arete a traiter
	 est-elle deja traitee ? => contour ferme
      */

      switch ( NEDGE ) {
      case XINF :
	if ( !(bufMarks[0][iy][ix] & (unsigned char)1) ) finparcours = 1;
	break;
      case XSUP :
	if ( !(bufMarks[0][iy+1][ix] & (unsigned char)1) ) finparcours = 1;
	break;
      case YINF :
	if ( !(bufMarks[0][iy][ix] & (unsigned char)2) ) finparcours = 1;
	break;
      case YSUP :
	if ( !(bufMarks[0][iy][ix+1] & (unsigned char)2) ) finparcours = 1;
	break;
      }

      if ( finparcours == 1 ) {
	if ( nparcours == 0 ) theCont->topology = _CLOSED_;
	nmaxparcours = 1;
	continue;
      }


      /* on doit alors ajouter le point 
	 - calcul de sa position
	 - marquage de l'arete
       */
      switch ( NEDGE ) {
      case XINF :
	x = ix + ( threshold - cube[0][0] ) / ( cube[0][1] - cube[0][0] );
	y = iy;
	bufMarks[0][iy][ix] &= ~(unsigned char)1;
	break;
      case XSUP :
	x = ix + ( threshold - cube[1][0] ) / ( cube[1][1] - cube[1][0] );
	y = iy + 1;
	bufMarks[0][iy+1][ix] &= ~(unsigned char)1;
	break;
      case YINF :
	x = ix;
	y = iy + ( threshold - cube[0][0] ) / ( cube[1][0] - cube[0][0] );
	bufMarks[0][iy][ix] &= ~(unsigned char)2;
	break;
      case YSUP :
	x = ix + 1;
	y = iy + ( threshold - cube[0][1] ) / ( cube[1][1] - cube[0][1] );
	bufMarks[0][iy][ix+1] &= ~(unsigned char)2;
	break;
      }
      
      if ( addPointToContour2D( x, y, theCont ) != 1 ) {
	free( resCont );
	return( NULL );
      }

      /*
      if ( 0 ) {
	printf( "   point : " );
	switch( NEDGE ) {
	case XINF : 
	  printf( "%4d %4d %4d", ix, iy, iz );
	  printf( " <-> %4d %4d %4d", ix+1, iy, iz ); break;
	case YINF : 
	  printf( "%4d %4d %4d", ix, iy, iz );
	  printf( " <-> %4d %4d %4d", ix, iy+1, iz ); break;
	case XSUP : 
	  printf( "%4d %4d %4d", ix, iy+1, iz );
	  printf( " <-> %4d %4d %4d", ix+1, iy+1, iz ); break;
	case YSUP : 
	  printf( "%4d %4d %4d", ix+1, iy, iz );
	  printf( " <-> %4d %4d %4d", ix+1, iy+1, iz ); break;
	}
	printf( " = %f %f %d\n",x ,y, iz );
      }
      */

      /* condition geometrique d'arret 
       */
      switch ( NEDGE ) {
      case YINF :
	if ( ix == 0 ) finparcours = 1;
	break;
      case YSUP :
  if ( ix+1 == (int)theIm->dim.x-1 ) finparcours = 1;
	break;
      case XINF :
	if ( iy == 0 ) finparcours = 1;
	break;
      case XSUP :
  if ( iy+1 == (int)theIm->dim.y-1 ) finparcours = 1;
	break;
      }
      if ( finparcours == 1 ) {
	theCont->topology = _OPEN_;
	continue;
      }


      /* passage a la configuration suivante
       */
      switch ( NEDGE ) {
      case XINF :   iy --;   CEDGE = XSUP;   break;
      case XSUP :   iy ++;   CEDGE = XINF;   break;
      case YINF :   ix --;   CEDGE = YSUP;   break;
      case YSUP :   ix ++;   CEDGE = YINF;   break;
     }
      
      
    } while ( finparcours == 0 );

  } /* fin de for ( nparcours = 0; nparcours < nmaxparcours; nparcours ++ ) */


  if ( 0 ) {
    printf( " liste #1 = %d    liste #2 = %d\n", theCont1.n, theCont2.n );
  }

  /* contour ouvert 
     => fusion des deux listes
  */
  if ( theCont2.n > 0 ) {
    if ( allocContour2D( resCont, theCont1.n+theCont2.n ) != 1 ) {
      free( resCont );
      return( NULL );
    }
    resCont->topology = _OPEN_;
    
    for ( i=0, n=theCont2.n-1; n>=0; i++, n-- ) {
      resCont->thePts[i].x = theCont2.thePts[n].x;
      resCont->thePts[i].y = theCont2.thePts[n].y;
    }
    for ( n=0; n<theCont1.n; i++, n++ ) {
      resCont->thePts[i].x = theCont1.thePts[n].x;
      resCont->thePts[i].y = theCont1.thePts[n].y;
    }
    resCont->n = theCont1.n+theCont2.n;

    freeContour2D( &theCont1 );
    freeContour2D( &theCont2 );

  }
  else {

    if ( allocContour2D( resCont, theCont1.n ) != 1 ) {
      free( resCont );
      return( NULL );
    }
    for ( i=0, n=0; n<theCont1.n; i++, n++ ) {
      resCont->thePts[i].x = theCont1.thePts[n].x;
      resCont->thePts[i].y = theCont1.thePts[n].y;
    }
    resCont->n = theCont1.n;
    resCont->topology = theCont1.topology;
    
    freeContour2D( &theCont1 );

  }

  return( resCont );
}










/* compute if a point is inside a closed contour
   or not
*/

int _PointInContour( typeContour2D *c,
		     double x, double y ) 
{
  typePoint2D m, r, w, u;
  typePoint2D p1, p2;
  double nr, nw;
  double n, prod, d, lambda;
  int p, nb_left, nb_right;
  

  /* on definit une droite partant d'une point,
     on prend un vecteur quelconque
  */
  m.x = x;
  m.y = y;
  r.x = 1;
  r.y = 0.1;
  nr = sqrt( r.x * r.x + r.y * r.y );

  nb_left = nb_right = 0;

  for ( p=0; p<c->n; p++ ) {
    p1.x = c->thePts[p].x;
    p1.y = c->thePts[p].y;
    p2.x = ( p == c->n-1 ) ? c->thePts[0].x : c->thePts[p+1].x;
    p2.y = ( p == c->n-1 ) ? c->thePts[0].y : c->thePts[p+1].y;

    w.x = p2.x - p1.x;
    w.y = p2.y - p1.y;
    nw = sqrt( w.x * w.x + w.y * w.y );

    d = ( r.x * w.y - r.y * w.x );
    if ( fabs ( d / (nr * nw) ) < 0.00001 ) {
      /*
       * les vecteurs r et w sont colineaires
       */
      if ( fabs( (p1.x - m.x) * r.y - (p1.y - m.y) * r.x ) / nr < 0.00001 ) {
	fprintf( stdout, "Warning (1): I should do something ..\n" );
      }
      continue;
    }
    
    /*
     * lambda / m + lambda * r appartient a la droite p1-p2
     * u = m + lambda * r - p1
     */

    lambda = ( (m.y - p1.y) * w.x - (m.x - p1.x) * w.y ) / d;
    if ( lambda == 0 ) {
      fprintf( stderr, "Warning(2): I should do something ..\n" );
      continue;
    }
    /*
      lambda > 0 : on est a gauche
      lambda < 0 :  "  "  " droite
    */
    
    u.x = m.x + lambda * r.x - p1.x;
    u.y = m.y + lambda * r.y - p1.y;
    
    /* 
     * u in [p1 p2]
     * => u.w > 0 && |u| < |w|
     */
    prod = u.x * w.x + u.y * w.y;
    if ( prod == 0.0 ) {
      fprintf( stderr, "Warning(3): I should do something ..\n" );
      continue;
    }
    if ( prod > 0.0 ) {
      n = sqrt( u.x * u.x + u.y * u.y );
      if ( n == nw )  {
	fprintf( stderr, "Warning(4): I should do something ..\n" );
	continue;
      }
      if ( n < nw ) {
	if ( lambda < 0 ) nb_right ++;
	if ( lambda > 0 ) nb_left  ++;
      }
    }
  }

  if ( nb_right % 2 == 1 && nb_left % 2 == 1 ) return( 1 );
  return( 0 );
}

























static int _ContoursPerCell( vt_image *theIm,
			     int x, int y, int z,
			     double threshold )
{
  int ith;
  int n;

  ith = (threshold>0) ? (int)(threshold+0.5) : (int)(threshold-0.5);
  n = 0;
  
  switch ( theIm->type ) {
  default :
    return( -1 );
  case UCHAR :
    {
      u8 ***theBuf = (u8***)theIm->array;
      if ( theBuf[z][y][x] >= ith )     n |= 1;
      if ( theBuf[z][y][x+1] >= ith )   n |= 2;
      if ( theBuf[z][y+1][x+1] >= ith ) n |= 4;
      if ( theBuf[z][y+1][x] >= ith )   n |= 8;
    }
    break;
  case USHORT :
    {
      u16 ***theBuf = (u16***)theIm->array;
      if ( theBuf[z][y][x] >= ith )     n |= 1;
      if ( theBuf[z][y][x+1] >= ith )   n |= 2;
      if ( theBuf[z][y+1][x+1] >= ith ) n |= 4;
      if ( theBuf[z][y+1][x] >= ith )   n |= 8;
    }
    break;
  case SSHORT :
    {
      s16 ***theBuf = (s16***)theIm->array;
      if ( theBuf[z][y][x] >= ith )     n |= 1;
      if ( theBuf[z][y][x+1] >= ith )   n |= 2;
      if ( theBuf[z][y+1][x+1] >= ith ) n |= 4;
      if ( theBuf[z][y+1][x] >= ith )   n |= 8;
    }
    break;
  case FLOAT :
    {
      r32 ***theBuf = (r32***)theIm->array;
      if ( theBuf[z][y][x] >= threshold )     n |= 1;
      if ( theBuf[z][y][x+1] >= threshold )   n |= 2;
      if ( theBuf[z][y+1][x+1] >= threshold ) n |= 4;
      if ( theBuf[z][y+1][x] >= threshold )   n |= 8;
    }
    break;
 }
  switch( nbcontours[n]) {
  case 0 :
  case 1 : return( nbcontours[n] );
  case 2 : return( 11 );
  }
  return( 111 );
}

int VT_ComputeNumberOfContoursPerCell( vt_image *theIm,
				       vt_image *resIm,
				       double threshold )
{
  u8 ***resBuf = (u8***)resIm->array;
  int x,y,z;

  if ( resIm->type != UCHAR ) return( 0 );
  
  for ( z=0; z<(int)resIm->dim.z; z++ )
  for ( y=0; y<(int)resIm->dim.y; y++ )
  for ( x=0; x<(int)resIm->dim.x; x++ )
    resBuf[z][y][x] = 0;

  for ( z=0; z<(int)resIm->dim.z; z++ ) {
    for ( y=0; y<(int)resIm->dim.y-1; y++ )
    for ( x=0; x<(int)resIm->dim.x-1; x++ ) {
      resBuf[z][y][x] = _ContoursPerCell( theIm, x, y, z, threshold );
    }
  }

  return( 1 );
}


