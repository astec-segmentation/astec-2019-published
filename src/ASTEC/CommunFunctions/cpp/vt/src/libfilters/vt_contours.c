
#include <stdio.h>
#include <stdlib.h>
/* #include <malloc.h> */
#include <string.h>
#include <math.h>
#include <unistd.h>

#include <vt_contours.h>

static int _verbose_ = 1;


void initContour2D( typeContour2D *c )
{
  c->topology = _OPEN_;
  c->n = 0;
  c->nalloc = 0;
  c->thePts = NULL;
}

int allocContour2D( typeContour2D *c, int n )
{
  char *proc = "allocContour2D";
  typePoint2D *thePts = NULL;
  if ( n <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad number of points (%d)\n", proc, n );
    return( -1 );
  }

  thePts = (typePoint2D *)malloc( n *sizeof(typePoint2D) );
  if ( thePts == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed (n=%d)\n", proc, n );
    return( -1 );
  }

  if ( c->thePts != NULL && c->nalloc > 0 ) {
    if ( c->n > 0 )
      memcpy( thePts, c->thePts, c->n * sizeof( typePoint2D ) );
    free( c->thePts );
  }

  c->thePts = thePts;
  c->nalloc = n;

  return( 1 );
}

void freeContour2D( typeContour2D *c )
{
  if ( c == NULL ) return;
  c->n = 0;
  if ( c->thePts != NULL ) free( c->thePts );
  initContour2D( c );
}

void printContour2D( typeContour2D *c, int p )
{
  int n;
  if ( c == NULL ) return;
  if( 0 ) fprintf( stdout, "\n" );
  fprintf( stdout, " points %d (alloc=%d)\n", c->n, c->nalloc );
  if ( p > 0 ) {
    for ( n=0; n<c->n; n++ ) 
      fprintf( stdout, "  #%3d: %f %f\n", n, c->thePts[n].x, c->thePts[n].y );
    fprintf( stdout, "\n" );
  }
}

int addPointToContour2D( double x, double y, typeContour2D *c )
{
  int nalloc;
  typePoint2D *thePts;
  
  if ( c == NULL ) return( -1 );

  if ( c->n == c->nalloc ) {

   nalloc = c->nalloc;
   nalloc += _POINTS_;
   thePts = (typePoint2D *)malloc( nalloc * sizeof(typePoint2D) );
   if ( thePts == NULL ) return( -1 );
   
   if ( c->n > 0 ) {
     memcpy( thePts, c->thePts, c->n*sizeof(typePoint2D) );
   }

   c->nalloc = nalloc;
   free( c->thePts );

   c->thePts = thePts;

  }

  c->thePts[c->n].x = x;
  c->thePts[c->n].y = y;
  c->n ++;

  return( 1 );
}

typeContour2D *copyContour2D( typeContour2D *c )
{
  typeContour2D *newc = NULL;

  if ( c == NULL ) return( NULL );

  newc = (typeContour2D *)malloc( sizeof( typeContour2D ) );
  if ( newc == NULL ) return( NULL );

  initContour2D( newc );

  *newc = *c;
  newc->thePts = NULL;
  newc->thePts = (typePoint2D *)malloc( c->n * sizeof(typePoint2D) );
  if ( newc->thePts == NULL ) {
    free( newc );
    return( NULL );
  }
  newc->nalloc = newc->n = c->n;
  memcpy( newc->thePts, c->thePts, c->n * sizeof(typePoint2D) );

  newc->topology = c->topology;
  return( newc );
}



double surfaceContour2D( typeContour2D *c, double x, double y )
{
  double s = 0;
  double ds;
  int n;

  if ( c == NULL || 
       c->topology == _OPEN_ ||
       c->n == 0 ||
       c->thePts == NULL )
    return( s );

  for ( n=1; n<c->n; n++ ) {
    ds = c->thePts[n-1].x * c->thePts[n].y 
      +  c->thePts[n].x   * y
      +  x                * c->thePts[n-1].y 
      -  x                * c->thePts[n].y 
      -  c->thePts[n].x   * c->thePts[n-1].y 
      -  c->thePts[n-1].x * y;
    /*
    printf( " (%f %f) - (%f %f) - (%f %f) = %f\n",
	    c->thePts[n-1].x, c->thePts[n-1].y,
	    c->thePts[n].x, c->thePts[n].y,
	    x, y, ds );
    */
    s += ds;
  }
  ds = c->thePts[c->n-1].x * c->thePts[0].y 
    +  c->thePts[0].x    * y
    +  x                 * c->thePts[c->n-1].y 
    -  x                 * c->thePts[0].y 
    -  c->thePts[0].x    * c->thePts[c->n-1].y 
    -  c->thePts[c->n-1].x * y;
  /*
  printf( " (%f %f) - (%f %f) - (%f %f) = %f\n",
	  c->thePts[c->n-1].x, c->thePts[c->n-1].y,
	  c->thePts[0].x, c->thePts[0].y,
	  x, y, ds );
  */
 
  s += ds;
  s /= 2.0;
  if ( s < 0 ) return( -s );
  return( s );
}




void initContour3D( typeContour3D *c )
{
  c->topology = _CLOSED_;
  c->n = 0;
  c->nalloc = 0;
  c->thePts = NULL;
}

int allocContour3D( typeContour3D *c, int n )
{
  char *proc = "allocContour3D";
  typePoint3D *thePts = NULL;
  if ( n <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad number of points (%d)\n", proc, n );
    return( -1 );
  }

  thePts = (typePoint3D *)malloc( n *sizeof(typePoint3D) );
  if ( thePts == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed (n=%d)\n", proc, n );
    return( -1 );
  }

  if ( c->thePts != NULL && c->nalloc > 0 ) {
    if ( c->n > 0 )
      memcpy( thePts, c->thePts, c->n * sizeof( typePoint3D ) );
    free( c->thePts );
  }

  c->thePts = thePts;
  c->nalloc = n;

  return( 1 );
}

void freeContour3D( typeContour3D *c )
{
  if ( c == NULL ) return;
  c->n = 0;
  if ( c->thePts != NULL ) free( c->thePts );
  initContour3D( c );
}



void printContour3D( typeContour3D *c, int p )
{
  int n;
  if ( c == NULL ) return;
  fprintf( stdout, "\n" );
  fprintf( stdout, "points %d (alloc=%d)\n", c->n, c->nalloc );
  if ( p > 0 ) {
    for ( n=0; n<c->n; n++ ) 
      fprintf( stdout, "#%3d: %f %f %f\n", p, c->thePts[p].x, 
	       c->thePts[p].y, c->thePts[p].z );
  }
}

int addPointToContour3D( double x, double y, double z,
			 typeContour3D *c )
{
  int nalloc;
  typePoint3D *thePts;
  
  if ( c == NULL ) return( -1 );

  if ( c->n == c->nalloc ) {

   nalloc = c->nalloc;
   nalloc += _POINTS_;
   thePts = (typePoint3D *)malloc( nalloc * sizeof(typePoint3D) );
   if ( thePts == NULL ) return( -1 );
   
   if ( c->n > 0 ) {
     memcpy( thePts, c->thePts, c->n*sizeof(typePoint3D) );
   }

   c->nalloc = nalloc;
   free( c->thePts );

   c->thePts = thePts;

  }

  c->thePts[c->n].x = x;
  c->thePts[c->n].y = y;
  c->thePts[c->n].y = z;
  c->n ++;

  return( 1 );
}

typeContour3D *copyContour3D( typeContour3D *c )
{
  typeContour3D *newc = NULL;

  if ( c == NULL ) return( NULL );

  newc = (typeContour3D *)malloc( sizeof( typeContour3D ) );
  if ( newc == NULL ) return( NULL );

  initContour3D( newc );

  *newc = *c;
  newc->thePts = NULL;
  newc->thePts = (typePoint3D *)malloc( c->n * sizeof(typePoint3D) );
  if ( newc->thePts == NULL ) {
    free( newc );
    return( NULL );
  }
  newc->nalloc = newc->n = c->n;
  memcpy( newc->thePts, c->thePts, c->n * sizeof(typePoint3D) );

  newc->topology = c->topology;
  return( newc );
}


int matlab_writeContour3D( FILE *mfile, int rfile, typeContour3D *c,
			   double *size,
			   char *desc, char *options )
{
  char *proc = "matlab_writeContour3D";
  double *t = NULL;
  int i;

  if ( c == NULL || c->n <= 0 ) return( 0 );
  t = (double*)malloc( c->n * sizeof( double ) );
  if ( t == NULL ) return( -1 );
  
  for ( i=0; i<c->n; i++ ) t[i] = c->thePts[i].x;
  if ( write( rfile, t, c->n * sizeof( double ) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }
  for ( i=0; i<c->n; i++ ) t[i] = c->thePts[i].y;
  if ( write( rfile, t, c->n * sizeof( double ) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }
  for ( i=0; i<c->n; i++ ) t[i] = c->thePts[i].z;
  if ( write( rfile, t, c->n * sizeof( double ) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }

  if ( desc != NULL ) {
    fprintf( mfile, "X%s = fread( fid, %d, 'float%lu' );\n", desc, c->n,  8*sizeof(double) );
    fprintf( mfile, "Y%s = fread( fid, %d, 'float%lu' );\n", desc, c->n,  8*sizeof(double) );
    fprintf( mfile, "Z%s = fread( fid, %d, 'float%lu' );\n", desc, c->n,  8*sizeof(double) );
    if ( size != NULL ) 
      fprintf( mfile, "plot3( %f * X%s, %f * Y%s, %f * Z%s", 
	       size[0], desc, size[1], desc, size[2], desc );
    else
      fprintf( mfile, "plot3( X%s, Y%s, Z%s", desc, desc, desc );
  }
  else {
    fprintf( mfile, "X = fread( fid, %d, 'float%lu' );\n", c->n,  8*sizeof(double) );
    fprintf( mfile, "Y = fread( fid, %d, 'float%lu' );\n", c->n,  8*sizeof(double) );
    fprintf( mfile, "Z = fread( fid, %d, 'float%lu' );\n", c->n,  8*sizeof(double) );
    if ( size != NULL ) 
      fprintf( mfile, "plot3( %f * X, %f * Y, %f * Z", size[0], size[1], size[2] );
    else
      fprintf( mfile, "plot3( X, Y, Z" );
  }
  if ( options != NULL ) fprintf( mfile, ", %s", options );
  fprintf( mfile, " );\n" );
  
  fprintf( mfile, "\n" );
  
  return( 1 );
}
















void initSlice( typeSlice *s )
{
  s->z = 0;
  s->mat = NULL;
  
  s->n = 0;
  s->nalloc = 0;
  s->theContours = NULL;
}  

int addContour2DToSlice( typeContour2D *c,
			 typeSlice *s )
{
  char *proc = "addContour2DToSlice";
  int cindex;
  int n, nalloc;
  typeContour2D **theContours = NULL;

  cindex = s->n;

  if ( cindex >= s->nalloc ) {

    nalloc = s->nalloc + _CONTOURS2D_;
    theContours = (typeContour2D **)malloc( nalloc * sizeof(typeContour2D*) );
    if ( theContours == NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate theContours\n", proc );
      return( -1 );
    }
    for ( n=0; n<nalloc; n++ ) theContours[n] = NULL;

    if ( s->nalloc > 0 && s->theContours != NULL ) {
      memcpy( theContours, s->theContours, s->nalloc * sizeof(typeContour2D*) );
      free( s->theContours );
    }

    s->theContours = theContours;
    s->nalloc = nalloc;
    
  }

  s->theContours[ cindex ] = c;
  s->n ++;

  return( 1 );
}

void freeSlice( typeSlice *s )
{
  int n;
  if ( s->mat != NULL ) free( s->mat );
  s->mat = NULL;

  if ( s->theContours == NULL ) {
    s->n = 0;
    s->nalloc = 0;
    return;
  }

  for ( n=0; n<s->n; n++ ) {
    freeContour2D( s->theContours[n] );
    free( s->theContours[n] );
  }
  free( s->theContours );
  initSlice( s );
}  

void printSlice( typeSlice *s, int p )
{
  int n;
  if ( s == NULL ) return;
  if( 0 ) fprintf( stdout, "\n" );
  if ( p > 0 && s->theContours != NULL) {
    for ( n=0; n<s->n; n++ ) {
      printContour2D( s->theContours[n], p-1 );
    }
  }
}


















void initStructure( typeStructure *s )
{
  s->n = 0;
  s->nalloc = 0;
  s->theSlices = NULL;
}

int addSliceToStructure( typeSlice *s,
			 typeStructure *structure )
{
  char * proc = "addContour2DToStructure";
  int sindex;
  int n, nalloc;
  typeSlice **theSlices = NULL;

  sindex = structure->n;

  if ( sindex >= structure->nalloc ) {
    
    nalloc = structure->nalloc + _SLICES_;
    theSlices = (typeSlice **)malloc( nalloc * sizeof(typeSlice*) );
    if ( theSlices == NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate slices\n", proc );
      return( -1 );
    }
    for ( n=0; n<nalloc; n++ ) theSlices[n] = NULL;
    
    if ( structure->nalloc > 0 && structure->theSlices != NULL ) {
      memcpy( theSlices, structure->theSlices, 
	      structure->nalloc * sizeof(typeSlice*) );
      free( structure->theSlices );
    }
    
    structure->theSlices = theSlices;
    structure->nalloc = nalloc;
    
  }
    
  structure->theSlices[ sindex ] = s;
  structure->n ++;
 
  return( 1 );
}

int addContour2DToStructure( typeContour2D *c,
			     typeStructure *structure, double z )
{
  char * proc = "addContour2DToStructure";
  int s, sindex = -1;
  double z_error = 0.0000001;
  typeSlice *slice = NULL;
  
  for ( s=0; s<structure->n && sindex == -1; s++ ) {
    if ( fabs( z - structure->theSlices[s]->z ) < z_error ) sindex = s;
  }

  if ( sindex == -1 ) {

    slice = (typeSlice *)malloc( sizeof(typeSlice) );
    if ( slice == NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate slice\n", proc );
      return( -1 );
    }
    initSlice( slice );
    slice->z = z;
    
    if ( addSliceToStructure( slice, structure ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to add slice to structure \n", proc );
      free( slice );
      return( -1 );
    }
    sindex = structure->n - 1;
  }

  return ( addContour2DToSlice( c, structure->theSlices[ sindex ] ) );
}

void freeStructure( typeStructure *s )
{
  int n;
  if ( s->theSlices == NULL ) {
    s->n = 0;
    s->nalloc = 0;
    return;
  }

  for ( n=0; n<s->n; n++ ) {
    freeSlice( s->theSlices[n] );
    free( s->theSlices[n] );
  }
  free( s->theSlices );
  initStructure( s );
}

void printStructure( typeStructure *s, int p )
{
  int n;
  if ( s == NULL ) return;
  fprintf( stdout, "\n" );
  fprintf( stdout, "slices %d (alloc=%d)\n", s->n, s->nalloc );
  if ( p > 0 && s->theSlices != NULL ) {
    for ( n=0; n<s->n; n++ ) {
      fprintf( stdout, "slice #%3d (z=%f)\n", n, s->theSlices[n]->z );
      printSlice( s->theSlices[n], p-1 );
    }
  }
  fprintf( stdout, "\n" );  
}









void initListOfContours3D( typeListOfContours3D *l )
{
  l->n = 0;
  l->nalloc = 0;
  l->theContours = NULL;
}

int addContour3DToListOfContours3D( typeContour3D *c, 
				    typeListOfContours3D *l )
{
  char *proc = "addContour3DToListOfContours3D";
  int cindex;
  int n, nalloc;
  typeContour3D **theContours = NULL;

  cindex = l->n;

  if ( cindex >= l->nalloc ) {

    nalloc = l->nalloc + _CONTOURS3D_;
    theContours = (typeContour3D **)malloc( nalloc * sizeof(typeContour3D*) );
    if ( theContours == NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate theContours\n", proc );
      return( -1 );
    }
    for ( n=0; n<nalloc; n++ ) theContours[n] = NULL;

    if ( l->nalloc > 0 && l->theContours != NULL ) {
      memcpy( theContours, l->theContours, l->nalloc * sizeof(typeContour3D*) );
      free( l->theContours );
    }

    l->theContours = theContours;
    l->nalloc = nalloc;
    
  }

  l->theContours[ cindex ] = c;
  l->n ++;

  return( 1 );
}

typeListOfContours3D* copyListOfContours3D( typeListOfContours3D *l )
{
  typeListOfContours3D *r = NULL;
  typeContour3D *c;
  int n;

  if ( l == NULL || l->n <= 0 ) return( NULL );

  r = (typeListOfContours3D *)malloc( sizeof( typeListOfContours3D ) );
  if ( r == NULL ) {
    return( NULL );
  }

  initListOfContours3D( r );

  r->theContours = (typeContour3D **)malloc( l->n * sizeof(typeContour3D*) );
  if ( r->theContours == NULL ) {
    free( r );
    return( NULL );
  }

  r->n = r->nalloc = l->n;
  for ( n=0; n<l->n; n ++ ) r->theContours[n] = NULL;
  
  for ( n=0; n<l->n; n ++ ) {
    c = copyContour3D( l->theContours[n] );
    if ( c == NULL ) {
      freeListOfContours3D( r );
      free( r );
      return( NULL );
    }
    r->theContours[n] = c;
  }
  
  return( r );
}

void freeListOfContours3D( typeListOfContours3D *l )
{
  int n;
  
  if ( l->theContours == NULL ) {
    l->n = 0;
    l->nalloc = 0;
    return;
  }

  for ( n=0; n<l->n; n++ ) {
    freeContour3D( l->theContours[n] );
    free( l->theContours[n] );
  }
  free( l->theContours );
  initListOfContours3D( l );
}

int matlab_writeListOfContours3D( FILE *mfile, int rfile, typeListOfContours3D *l,
				  double *size,
				  char *desc, char *options )
{
  int i;
  char d[512];
  if ( l == NULL || l->n <=0 ) return( 0 );
  for (i=0; i<l->n; i++ ) {
    if ( desc == NULL ) sprintf( d, "%d", i );
    else sprintf( d, "%s%d", desc, i );
    matlab_writeContour3D( mfile, rfile, l->theContours[i], size, d, options );
  }

  return( 1 );
}






/*********************************************************************************
 ********************************************************************************/



int cnt_readStructure( typeStructure *structure, FILE *f )
{
  char *proc = "cnt_readStructure";
  char str[256];
  int s, slices = 0;
  int v, vertices;
  double z;
  int c, p=0;  
  typePoint2D *pts = NULL;
  typeContour2D *cnt = NULL;
  
  if ( fgets( str, 256, f ) == NULL ) {
    fprintf( stderr, "%s: error when reading\n", proc );
  }
  if ( sscanf( str, "S %d", &slices ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to get slices number\n", proc );
    return( -1 );
  }

  for ( s=0; s<slices; s++ ) {
    
    if ( fgets( str, 256, f ) == NULL ) {
      fprintf( stderr, "%s: error when reading\n", proc );
    }

    if ( sscanf( str, "v %d z %lf", &vertices, &z ) != 2 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to get vertices number and/or Z\n", proc );
      return( -1 );
    }
    
    if ( vertices > 0 ) {
      pts = (typePoint2D *)malloc( vertices * sizeof( typePoint2D ) );
      if ( pts == NULL ) {
	if ( _verbose_ )
	  fprintf( stderr, "%s: unable to allocate auxiliary points array\n", proc );
	return( -1 );
      }
    }

    c = v = 0;
    do {

      if ( fgets( str, 256, f ) == NULL ) {
	fprintf( stderr, "%s: error when reading\n", proc );
      }

      /* new contour
       */
      if ( str[0] == '{' ) {

	cnt = (typeContour2D *)malloc( sizeof(typeContour2D) );
	if ( cnt == NULL ) {
	  if ( _verbose_ )
	    fprintf( stderr, "%s: unable to allocate new contour\n", proc );	
	  if ( pts != NULL ) free( pts );
	  return( -1 );
	}
	initContour2D( cnt );
	p = 0;
      }
      
      /* end of contour 
       */
      else if ( str[0] == '}' ) {

	if ( allocContour2D( cnt, p ) != 1 ) {
	  if ( _verbose_ )
	    fprintf( stderr, "%s: unable to allocate points array\n", proc );
	  if ( pts != NULL ) free( pts );
	  if ( cnt != NULL ) free( cnt );
	  return( -1 );
	}
	memcpy( cnt->thePts, pts, p*sizeof( typePoint2D ) );
	cnt->n = p;

	if ( addContour2DToStructure( cnt, structure, z ) != 1 ) {
	  if ( _verbose_ )
	    fprintf( stderr, "%s: unable to add contour to structure\n", proc );
	  if ( pts != NULL ) free( pts );
	  freeContour2D( cnt );
	  if ( cnt != NULL ) free( cnt );
	  return( -1 );
	}
	cnt = NULL;
	c++;
      }

      /* point of contour
       */
      else {
	if ( sscanf( str, "%lf %lf", &pts[p].x,  &pts[p].y ) != 2 ) {
	  if ( _verbose_ )
	    fprintf( stderr, "%s: unable to read point #p (cont #%d, slice #%d)\n", 
		     proc, c, s );
	  if ( pts != NULL ) free( pts );
	  if ( cnt != NULL ) free( cnt );
	  return( -1 );
	}
	p++;
	v++;
      }
      
    } while ( v<vertices || cnt != NULL );



    if ( pts != NULL ) free( pts );
    pts = NULL;

  }
  return( 1 );
}



int act2D_readStructure( typeStructure *s __attribute__ ((unused)),
                         FILE *f __attribute__ ((unused)) )
{
  char *proc = "act2D_readStructure";
  if ( _verbose_ )
    fprintf( stderr, "%s: not implemented yet\n", proc );
  return( 0 );
}



int act3D_readStructure( typeStructure *s __attribute__ ((unused)),
                         FILE *f __attribute__ ((unused)) )
{
  char *proc = "act3D_readStructure";
  if ( _verbose_ )
    fprintf( stderr, "%s: not implemented yet\n", proc );
  return( 0 );
}





int cnt_writeStructure( typeStructure *structure, FILE *f )
{
  int s, v, c, p;
  
  fprintf( f, "S %d\n", structure->n );

  for ( s=0; s<structure->n; s++ ) {

    for ( v=0, c=0; c<structure->theSlices[s]->n; c++ )
      v += structure->theSlices[s]->theContours[c]->n;
    
    fprintf( f, "v %d z %g\n", v, structure->theSlices[s]->z );

    for ( c=0; c<structure->theSlices[s]->n; c++ ) {
      fprintf( f, "{\n" );
      for ( p=0; p<structure->theSlices[s]->theContours[c]->n; p++ )
	fprintf( f, "%g %g\n", 
		 structure->theSlices[s]->theContours[c]->thePts[p].x,
		 structure->theSlices[s]->theContours[c]->thePts[p].y );
      fprintf( f, "}\n" );
    }

  }
  return( 1 );
}



int act2D_writeStructure( typeStructure *s __attribute__ ((unused)),
                          FILE *f __attribute__ ((unused)) )
{
  char *proc = "act2D_writeStructure";
  if ( _verbose_ )
    fprintf( stderr, "%s: not implemented yet\n", proc );
  return( 0 );
}



int act3D_writeStructure( typeStructure *structure, FILE *f )
{
  int s, n, c, p;

  fprintf(f, "# ActiveContour3D File version 1\n");
  fprintf(f, "CONTENT\n");
  fprintf(f, "Number of meshes = 1\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "Big endian\n");
  fprintf(f, "END\n");
  fprintf(f, "ActiveContour3D dummy\n");
  fprintf(f, "Main\n");
  fprintf(f, "END\n");
  
  for ( n=0, s=0; s<structure->n; s++ )
    n += structure->theSlices[s]->n;
  fprintf(f, "%d\n", n );

  for ( s=0; s<structure->n; s++ ) {

    for ( c=0; c<structure->theSlices[s]->n; c++ ) {

      switch( structure->theSlices[s]->theContours[c]->topology ) {
      case _CLOSED_ : fprintf(f, "1 " ); break;
      default :       fprintf(f, "0 " );
      }
      fprintf(f, "%d\n", structure->theSlices[s]->theContours[c]->n );
      for ( p=0; p<structure->theSlices[s]->theContours[c]->n; p++ )
	fprintf( f, "%g %g %g\n",
		 structure->theSlices[s]->theContours[c]->thePts[p].x,
		 structure->theSlices[s]->theContours[c]->thePts[p].y,
		 structure->theSlices[s]->z );
    }
  }
  
  return( 1 );
}



int matlab_writeStructure( typeStructure *s __attribute__ ((unused)),
                           char *filename __attribute__ ((unused)) )
{
  char *proc = "matlab_writeStructure";
  if ( _verbose_ )
    fprintf( stderr, "%s: not implemented yet\n", proc );
  return( 0 );
}



int readStructure( typeStructure *s,
		   char *filename )
{
  char *proc = "readStructure";
  int (*_read_structure_)( typeStructure *, FILE * ) = NULL;
  FILE *file = NULL;
  int r;

  if ( filename == NULL ) return( -1 );
  if ( strlen( filename ) <= 4 ) return( -1 );

  if ( strncmp( filename+strlen(filename)-3, "cnt", 3 ) == 0 ) {
    _read_structure_ = & cnt_readStructure;
  }
  else if ( strncmp( filename+strlen(filename)-5, "act2D", 5 ) == 0 ) {
    _read_structure_ = & act2D_readStructure;
  }
  else if ( strncmp( filename+strlen(filename)-5, "act3D", 5 ) == 0 ) {
    _read_structure_ = & act3D_readStructure;
  }
  else {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unknown file format for '%s'\n", proc, filename );
    return( -1 );
  }
  
  file = fopen( filename, "r" );
  if ( file == NULL ) {
    fprintf( stderr, "%s: unable to open '%s'\n", proc, filename );
    return( -1 );
  }
  r = (*_read_structure_)( s, file );
  fclose( file );

  return( r );
}

int writeStructure( typeStructure *s,
		   char *filename )
{
  char *proc = "writeStructure";
  int (*_write_structure_)( typeStructure *, FILE * ) = NULL;
  FILE *file = NULL;
  int r;

  if ( filename == NULL ) return( -1 );
  if ( strlen( filename ) <= 4 ) return( -1 );

  if ( strncmp( filename+strlen(filename)-3, "cnt", 3 ) == 0 ) {
    _write_structure_ = & cnt_writeStructure;
  }
  else if ( strncmp( filename+strlen(filename)-5, "act2D", 5 ) == 0 ) {
    _write_structure_ = & act2D_writeStructure;
  }
  else if ( strncmp( filename+strlen(filename)-5, "act3D", 5 ) == 0 ) {
    _write_structure_ = & act3D_writeStructure;
  }
  else if ( strncmp( filename+strlen(filename)-2, ".m", 2 ) == 0 ) {
    return( matlab_writeStructure( s, filename ) );
  }
  else {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unknown file format for '%s'\n", proc, filename );
    return( -1 );
  }
  
  file = fopen( filename, "w" );
  if ( file == NULL ) {
    fprintf( stderr, "%s: unable to open '%s'\n", proc, filename );
    return( -1 );
  }
  r = (*_write_structure_)( s, file );
  fclose( file );

  return( r );
}
