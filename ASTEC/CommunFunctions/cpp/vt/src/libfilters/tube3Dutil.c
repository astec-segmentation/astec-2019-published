/*************************************************************************
 * tube3Dutil.c -
 *
 * $Id: tube3Dutil.c,v 1.2 2001/06/08 15:10:52 greg Exp $
 *
 * Copyright (c) INRIA 2000
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Mon Sep 11 18:20:24 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 */

#include <tube3Dutil.h>

#define _TRACE_ 0

static int _verbose_ = 0;


void _VerboseInTube3Dutil()
{
  if ( _verbose_ <= 0 ) _verbose_ = 1;
  else _verbose_ ++;
}
void _NoVerboseInTube3Dutil()
{
  _verbose_ = 0;
}










typedef enum {
  UNKNOWN,
  COMPONENT,
  JUNCTION
} enumType;

#define MAXNEIGHBORS 8

typedef struct {

  int label;
  int size;

  int xmin;
  int xmax;
  int ymin;
  int ymax;
  int zmin;
  int zmax;

  enumType type;

  int nbNeighbors;
  int labelsNeighbors[MAXNEIGHBORS];

} typeCC;











/* destruction des branches trop courtes
   
   l'image binaire contient des courbes (3D ou 2D)
   deja aminicies. On peut donc caracteriser
   en comptant le nombre de voisins :
   1 : point simple (fin de courbe)
   2 : point de courbe simple 
   2 ou + : point de jonction

*/

int remove_small_simple_curves( unsigned char *input_buffer,
				unsigned char *output_buffer,
				int *theDim,
				int minimal_length )
{
  char *proc = "remove_small_simple_curves";

  unsigned short int *tmp_buffer=NULL, *cc_buffer = NULL;
  int dimx = theDim[0];
  int dimy = theDim[1];
  int dimz = theDim[2];
  int dimxy = theDim[0]*theDim[1];
  int x, y, z, i, j, k, l, m;
  int iz, iy;
  int offset[27];
  int nb_neighbors;

  int nb_curves, nb_remaining_curves;
  int nb_junctions, nb_remaining_junctions;

  typeCC *theCC = (typeCC *)NULL;
  int n, e, b, jct;
  int neighbors[27];
  int min, max;
  int label;
  


  /* initialisation des offsets
   */
  offset[12] = -1;
  offset[13] =  0;
  offset[14] =  1;
  for ( i=-1; i<=1; i++ ) {
    offset[10+i] = offset[13+i] - dimx;
    offset[16+i] = offset[13+i] + dimx;
  }
  for ( i=-4; i<=4; i++ ) {
    offset[ 4+i] = offset[13+i] - dimx*dimy;
    offset[22+i] = offset[13+i] + dimx*dimy;
    
  }

  
  /* allocations des buffers auxiliaires
   */
  tmp_buffer = (unsigned short int *)calloc( dimz*dimxy, sizeof( unsigned short int ) );
  if ( tmp_buffer == NULL ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate first auxiliary buffer\n", proc );
    return( -1 );
  }
  cc_buffer = (unsigned short int *)calloc( dimz*dimxy, sizeof( unsigned short int ) );
  if ( cc_buffer == NULL ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate second auxiliary buffer\n", proc );
    free( tmp_buffer );
    return( -1 );
  }
  

  /* initialisations des buffers auxiliaires
     pour le calcul des composantes connexes
  */

  for ( i=0, z=0; z<dimz; z++ ) {
    iz =  ( z == 0 || z == dimz-1 ) ? 0 : 1;
    for ( y=0; y<dimy; y++ ) {
      iy = ( iz == 0 || y == 0 || y == dimy-1 ) ? 0 : 1;
      for ( x=0; x<dimx; x++, i++ ) {

	/* background
	 */
	if ( input_buffer[i] == 0 ) {
	  continue;
	}

	/* count neighbors,
	   the central point should not be counted
	   then remove 1 at the end
	*/
	nb_neighbors = 0;
	/* are we inside or outside ?
	 */
	if ( iy == 0 || x == 0 || x == dimx-1 ) {
	  for ( l = -1; l < 2; l++ ) 
	  for ( k = -1; k < 2; k++ ) 
	  for ( j = -1; j < 2; j++ ) {
	    if ( x+j >= 0 &&  x+j < dimx &&
		 y+k >= 0 &&  y+k < dimy &&
		 z+l >= 0 &&  z+l < dimz )
	      if ( input_buffer[ (z+l)*dimxy + (y+k)*dimy + (x+j) ] ) nb_neighbors ++;
	  }
	} else {
	  for ( j=0; j<27; j++ )
	    if ( input_buffer[ i + offset[j] ] ) nb_neighbors ++;
	}
	nb_neighbors --;

	/* admissible values for nb_neighbors are in [0..8]
	   nb_neighbors = 0 : isolated point 
	   nb_neighbors = 1 : end point
	   nb_neighbors = 2 : simple curve point
	   nb_neighbors > 2 : curves junction point
	*/
	if ( nb_neighbors >= 0 && nb_neighbors <= 2 )
	  cc_buffer[i] = 255;
	else
	  tmp_buffer[i] = 255;
	
      }
    }
  }


  /* count connected components for junctions
     junctions will be labeled from 1 to nb_junctions
   */
  nb_junctions = CountConnectedComponentsWithAllParams( tmp_buffer, USHORT,
						     tmp_buffer, USHORT,
						     theDim, 1.0,
						     26, /* connectivity */
						     1, /* minimal size of component */
						     0, /* all components */
						     0  /* no binary output */ );

  if ( nb_junctions == -1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when counting junctions\n", proc );
    free( cc_buffer );
    free( tmp_buffer );
    return( -1 );
  }

  /* no junctions
     => there isn't any "small" spurious branches
     it's done
  */
  if ( nb_junctions == 0 ) {
    free( cc_buffer );
    free( tmp_buffer );
    if ( input_buffer != output_buffer ) {
      (void)memcpy( output_buffer, input_buffer, dimx*dimy*dimz );
    }
    return( 1 );
  }
  

  
  /* count connected components for simple curves
     curves will be labeled from 1 to nb_curves
     curves are then relabeled in decreasing order
     with respect to size
   */
  nb_curves = CountConnectedComponentsWithAllParams( cc_buffer, USHORT,
						     cc_buffer, USHORT,
						     theDim, 1.0,
						     26, /* connectivity */
						     1, /* minimal size of component */
						     0, /* all components */
						     0  /* no binary output */ );
  if ( nb_curves == -1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when counting simple curves\n", proc );
    free( cc_buffer );
    free( tmp_buffer );
    return( -1 );
  }
  if ( nb_curves == 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: no simple curves found in image\n", proc );
    free( cc_buffer );
    free( tmp_buffer );
    return( -1 );
  }

  if ( RelabelConnectedComponentsByDecreasingSize( cc_buffer, USHORT,  theDim ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to relabel simple curves\n", proc );
    free( cc_buffer );
    free( tmp_buffer );
    return( -1 );
  }
  







  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "%s: found %d curves and %d junctions\n", proc, nb_curves, nb_junctions );
  }

  /* translate labels and copy junctions into
     the simple curve image
     it will save space
  */
  if ( nb_curves + nb_junctions > 65535 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: too many  curves (%d) and junctions (%d) in image\n",
	       proc, nb_curves, nb_junctions );
    free( cc_buffer );
    free( tmp_buffer );
    return( -1 );
  }

  for ( i=0; i<dimz*dimxy; i++ ) {
    if ( tmp_buffer[i] == 0 ) continue;
    cc_buffer[i] = tmp_buffer[i] + nb_curves;
  }

  free( tmp_buffer );
  




  /* here we've got
     from 1           to nb_curves              : simple curves
     from 1+nb_curves to nb_junctions+nb_curves : junctions
  */




  /* array of components
   */

  theCC = (typeCC*)malloc( (nb_junctions+nb_curves+1)*sizeof(typeCC) );
  if ( theCC == (typeCC *)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, " %s: can not allocate auxiliary array.\n", proc );
    free( cc_buffer );
    return( -1 );
  }
  


  for ( i=0; i<=nb_junctions+nb_curves; i++ ) {

    theCC[i].label = i;
    theCC[i].size = 0;

    theCC[i].xmin = dimx - 1;
    theCC[i].xmax = 0;
    theCC[i].ymin = dimy - 1;
    theCC[i].ymax = 0;
    theCC[i].zmin = dimz - 1;
    theCC[i].zmax = 0;

    if ( i >= 1 && i <= nb_curves ) 
      theCC[i].type = COMPONENT;
    else if ( i >= nb_curves+1 && i <= nb_curves+nb_junctions ) 
      theCC[i].type = JUNCTION;
    else
      theCC[i].type = UNKNOWN;
    
    theCC[i].nbNeighbors = 0;
    for ( j=0; j<MAXNEIGHBORS; j++ )
      theCC[i].labelsNeighbors[j] = 0;

  }





  
  for ( i=0, z=0; z<dimz; z++ ) {
    iz =  ( z == 0 || z == dimz-1 ) ? 0 : 1;
    for ( y=0; y<dimy; y++ ) {
      iy = ( iz == 0 || y == 0 || y == dimy-1 ) ? 0 : 1;
      for ( x=0; x<dimx; x++, i++ ) {

	label = cc_buffer[i];
	if ( label == 0 ) continue;
	
	/* increment size
	   update bounding box
	*/
	theCC[ label ].size ++;
	
	if ( theCC[ label ].xmin > x ) theCC[ label ].xmin = x;
	if ( theCC[ label ].xmax < x ) theCC[ label ].xmax = x;
	if ( theCC[ label ].ymin > y ) theCC[ label ].ymin = y;
	if ( theCC[ label ].ymax < y ) theCC[ label ].ymax = y;
	if ( theCC[ label ].zmin > z ) theCC[ label ].zmin = z;
	if ( theCC[ label ].zmax < z ) theCC[ label ].zmax = z;


	/* looking for neighbor component
	 */
	if ( iy == 0 || x == 0 || x == dimx-1 ) {
	  for ( n = 0, l = -1; l < 2; l++ ) 
	  for ( k = -1; k < 2; k++ ) 
	  for ( j = -1; j < 2; j++, n++ ) {
	    neighbors[n] = 0;
	    if ( x+j >= 0 &&  x+j < dimx &&
		 y+k >= 0 &&  y+k < dimy &&
		 z+l >= 0 &&  z+l < dimz ) {
	      neighbors[n] = cc_buffer[ (z+l)*dimxy + (y+k)*dimy + (x+j) ];
	    }
	  }
	} else {
	  for ( j=0; j<27; j++ ) {
	    neighbors[j] = cc_buffer[ i + offset[j] ];
	  }
	}

	for ( j=0; j<27; j++ ) {
	  if ( neighbors[j] == 0 || neighbors[j] == label ) continue;
	  /* les 3 lignes ci-dessous seraient a supprimer
	     pour que tout se passe bien en presence de boucles.
	  */
	  /*
	    for ( e=0, n=0; n < theCC[ label ].nbNeighbors; n++ )
	    if ( theCC[ label ].labelsNeighbors[n] == neighbors[j] ) e = 1;
	    if ( e == 1 ) continue;
	  */
	  theCC[ label ].labelsNeighbors[ theCC[ label ].nbNeighbors ] = neighbors[j];
	  theCC[ label ].nbNeighbors ++;
	}

	
      }
    }
  }


  if ( _TRACE_ ) {
    for ( i=1; i<=nb_junctions+nb_curves; i++ ) {
      if ( i <= nb_curves )
	printf( " COMPONENT #%-5d: ", i );
      else
	printf( " JUNCTION  #%-5d: ", i );
      printf( "size=%3d, ", theCC[ i ].size );
      printf( "(%d) = ", theCC[ i ].nbNeighbors );
      printf( "{ " );
      for ( n=0; n < theCC[ i ].nbNeighbors; n++ )
	printf( "%5d ", theCC[ i ].labelsNeighbors[n] );
      printf( "}\n" );
    }
    printf( "\n" );
  }






  
  /* processing the simple curves
   */

  nb_remaining_curves    = nb_curves;
  nb_remaining_junctions = nb_junctions;

  for ( b = nb_curves; (b >= 1) && (theCC[ b ].size <= minimal_length); b-- ) {
  
    /* on ne peut traiter que celles qui sont independantes
     */
    if ( theCC[ b ].label != b ) continue;

    /* on ne peut traiter que celles qui sont terminales
     */
    if ( theCC[ b ].nbNeighbors != 1 ) continue;


    /* here we've got a component that can be deleted
       it is connected to one single junction, #jct
     */
    jct = theCC[ b ].labelsNeighbors[0];


    if ( _TRACE_ ) {
      printf( " deleting component #%-5d ... processing junction #%-5d", b, jct );
    }


    /* delete the component
     */

    theCC[ b ].label = 0;
    for ( iz=theCC[ b ].zmin*dimxy,   z=theCC[ b ].zmin; z<=theCC[ b ].zmax; z++, iz+=dimxy )
    for ( iy=iz+theCC[ b ].ymin*dimx, y=theCC[ b ].ymin; y<=theCC[ b ].ymax; y++, iy+=dimx  )
    for ( i=iy+theCC[ b ].xmin,       x=theCC[ b ].xmin; x<=theCC[ b ].xmax; x++, i++ ) {
      if ( cc_buffer[i] == b ) cc_buffer[i] = 0;
    }

    theCC[ b ].size  = 0;
    theCC[ b ].label = 0;
    theCC[ b ].nbNeighbors = 0;
    nb_remaining_curves --;


    /* now we have to process the junction
       first of all, some points of the junction can be deleted
    */

    if ( _TRACE_ ) {
      printf( " (%d->", theCC[ jct ].size );
    }

    do {

      m = 0;
      
      for ( iz=theCC[jct].zmin*dimxy,   z=theCC[jct].zmin; z<=theCC[jct].zmax; z++, iz+=dimxy )
      for ( iy=iz+theCC[jct].ymin*dimx, y=theCC[jct].ymin; y<=theCC[jct].ymax; y++, iy+=dimx  )
      for ( i=iy+theCC[jct].xmin,       x=theCC[jct].xmin; x<=theCC[jct].xmax; x++, i++ ) {
	
	if ( cc_buffer[i] != jct ) continue;
	
	for ( n = 0, l = -1; l < 2; l++ ) 
	for ( k = -1; k < 2; k++ ) 
	for ( j = -1; j < 2; j++, n++ ) {
	  neighbors[n] = 0;
	  if ( x+j >= 0 &&  x+j < dimx &&
	       y+k >= 0 &&  y+k < dimy &&
	       z+l >= 0 &&  z+l < dimz )
	    neighbors[n] = cc_buffer[ i+offset[n] ];
	}
	
	if ( IsA3DPointSimple( neighbors ) == 1 ) {
	  cc_buffer[ i ] = 0;
	  m++;
	}
	
      }
      
      theCC[ jct ].size -= m;

    } while ( m > 0 );

    if ( _TRACE_ ) {
      printf( "%d)\n", theCC[ jct ].size );
    }



    /* second, we count the remaining neighbors 
       of the junction
    */
    for ( e=0, n=0; n<theCC[ jct ].nbNeighbors; n++ ) {
      if ( theCC[ jct ].labelsNeighbors[n] != b ) {
	theCC[ jct ].labelsNeighbors[e] = theCC[ jct ].labelsNeighbors[n];
	e++;
      }
    }
    if ( e == theCC[ jct ].nbNeighbors ) {
      fprintf( stderr, "%s: FATAL ERROR, component #%d is not a neighbor of junction #%d\n",
	       proc, b, jct );
    } 
    theCC[ jct ].nbNeighbors = e;



    /* several case :
       - first case: it is still a junction
       - second case: it is no more a junction
         * first exception: a branche connected to a circle
           the circle could be deleted (depending on how
	   we have extract neighbors, which is not completely clear)
	   and then the junction could have disappeared
	   It should have one neighbor in this case
    */

    if ( theCC[ jct ].nbNeighbors > 2 ) continue;
    

    nb_remaining_junctions --;


    if ( theCC[ jct ].nbNeighbors == 1 ) {
      if ( theCC[ jct ].size != 0 ) {
	fprintf( stderr, "%s: FATAL ERROR, junction #%d is not completely deleted\n",
		 proc, jct );
	
	min = theCC[ theCC[ jct ].labelsNeighbors[0] ].label;
	theCC[ min ].size += theCC[ jct ].size;
	theCC[ jct ].label = min;
	theCC[ jct ].type = COMPONENT;
	theCC[ jct ].size = 0;
	
      } else {
	theCC[ jct ].label = 0;
      }
      
      theCC[ jct ].nbNeighbors = 0;
      
      continue;
    }
    




    /* it is no more a junction
       it has two neighbors
       - in general it has 2 differents neighbors
       - in case of loop (a circle + a branche)
         it may have 2 identical neighbors (the branche has been deleted)

       do some relabeling
       the junction and its neighbors are fused and get the same label
       two neighbors: we first get the two labels (#min and #max)
                      process the second exception (2 identical neighbors)
    */
    



    if ( theCC[ theCC[ jct ].labelsNeighbors[0] ].label 
	 > theCC[ theCC[ jct ].labelsNeighbors[1] ].label ) {
      max = theCC[ theCC[ jct ].labelsNeighbors[0] ].label;
      min = theCC[ theCC[ jct ].labelsNeighbors[1] ].label;
    } else {
      min = theCC[ theCC[ jct ].labelsNeighbors[0] ].label;
      max = theCC[ theCC[ jct ].labelsNeighbors[1] ].label;
    }
    
    

    if ( _TRACE_ ) {
      printf( "    merging junction #%-5d and components #%-5d and #%-5d\n",
	      jct, min, max );
    }


    /* add the junction to component #min
     */

    theCC[ min ].size += theCC[ jct ].size;
    if ( theCC[ min ].xmin > theCC[ jct ].xmin ) theCC[ min ].xmin = theCC[ jct ].xmin;
    if ( theCC[ min ].xmax < theCC[ jct ].xmax ) theCC[ min ].xmax = theCC[ jct ].xmax;
    if ( theCC[ min ].ymin > theCC[ jct ].ymin ) theCC[ min ].ymin = theCC[ jct ].ymin;
    if ( theCC[ min ].ymax < theCC[ jct ].ymax ) theCC[ min ].ymax = theCC[ jct ].ymax;
    if ( theCC[ min ].zmin > theCC[ jct ].zmin ) theCC[ min ].zmin = theCC[ jct ].zmin;
    if ( theCC[ min ].zmax < theCC[ jct ].zmax ) theCC[ min ].zmax = theCC[ jct ].zmax;
    theCC[ jct ].label = min;
    theCC[ jct ].type = COMPONENT;
    theCC[ jct ].size = 0;
    theCC[ jct ].nbNeighbors = 0;

    
    if ( min == max ) {
      continue;
    }

    


    /* we have to merge component #max with #min
       1. manage the equivalence table
       2. manage remaining neighbors of #max and #min
          if #max has two neighbors (sat #jct and #e), 
	     the #jct neighbor of #min is changed into #e
	     and the #max neighbor of #e is changed into #min
          if #max has one neighbor, #jct is just deleted
             from the neighbors of #min
     */
    nb_remaining_curves --;

    
    for ( i=max; i<=nb_curves+nb_junctions; i++ ) {
      if ( theCC[i].label == max ) theCC[i].label = min;
    }
    
    theCC[ min ].size += theCC[ max ].size;
    if ( theCC[ min ].xmin > theCC[ max ].xmin ) theCC[ min ].xmin = theCC[ max ].xmin;
    if ( theCC[ min ].xmax < theCC[ max ].xmax ) theCC[ min ].xmax = theCC[ max ].xmax;
    if ( theCC[ min ].ymin > theCC[ max ].ymin ) theCC[ min ].ymin = theCC[ max ].ymin;
    if ( theCC[ min ].ymax < theCC[ max ].ymax ) theCC[ min ].ymax = theCC[ max ].ymax;
    if ( theCC[ min ].zmin > theCC[ max ].zmin ) theCC[ min ].zmin = theCC[ max ].zmin;
    if ( theCC[ min ].zmax < theCC[ max ].zmax ) theCC[ min ].zmax = theCC[ max ].zmax;
    /* theCC[ max ].label = min; */
    /* theCC[ max ].type = COMPONENT; */
    theCC[ max ].size = 0;
    /* theCC[ max ].nbNeighbors = 0; */

    

    
    if ( theCC[ max ].nbNeighbors == 2 ) {

      if ( theCC[ max ].labelsNeighbors[0] == jct ) e = theCC[ max ].labelsNeighbors[1];
      else                                          e = theCC[ max ].labelsNeighbors[0];
      if ( theCC[ min ].labelsNeighbors[0] == jct ) theCC[ min ].labelsNeighbors[0] = e;
      else                                          theCC[ min ].labelsNeighbors[1] = e;
      for ( n=0; n<theCC[ e ].nbNeighbors; n++ )
	if ( theCC[ e ].labelsNeighbors[n] == max ) theCC[ e ].labelsNeighbors[n] = min;

    } else {
      
      for ( e=0, n=0; n<theCC[ min ].nbNeighbors; n++ ) {
	if ( theCC[ min ].labelsNeighbors[n] != jct ) {
	  theCC[ min ].labelsNeighbors[e] = theCC[ min ].labelsNeighbors[n];
	  e++;
	}
      }
      if ( e == theCC[ min ].nbNeighbors ) {
	fprintf( stderr, "%s: FATAL ERROR, junction #%d is not a neighbor of component #%d\n",
		 proc, jct, min );
      }
      theCC[ min ].nbNeighbors = e;

      /* should be #min examined?
       */
      if ( min >= b && theCC[ min ].size <= minimal_length ) {
	b = min + 1;
      }
      
    }

    theCC[ max ].nbNeighbors = 0;

  }


  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "%s: %d curves and %d junctions are remaining\n", proc, 
	     nb_remaining_curves, nb_remaining_junctions );
  }



  if ( _TRACE_ ) {
    printf( "\n" );
    printf( " remaining curves = %d\n", nb_remaining_curves );
    printf( "\n" );
    for ( i=1; i<=nb_junctions+nb_curves; i++ ) {
      if ( theCC[ i ].size == 0 ) continue;
      if ( i <= nb_curves )
	printf( " COMPONENT #%-5d: ", i );
      else
	printf( " JUNCTION  #%-5d: ", i );
      printf( "size=%3d, ", theCC[ i ].size );
      printf( "(%d) = ", theCC[ i ].nbNeighbors );
      printf( "{ " );
      for ( n=0; n < theCC[ i ].nbNeighbors; n++ )
	printf( "%5d ", theCC[ i ].labelsNeighbors[n] );
      printf( "}" );
      printf( "\n" );
    }
    printf( "\n" );
  }




  
  /* relabel ...
   */
  for ( j=0, i=1; i<=nb_junctions+nb_curves; i++ ) {
    if ( theCC[i].label == i )
      theCC[i].label = ++j;
    else 
      theCC[i].label = theCC[ theCC[i].label ].label;
  }

  for ( j=0, i=1; i<=nb_junctions+nb_curves; i++ ) {
    for ( n=0; n < theCC[ i ].nbNeighbors; n++ )
      theCC[ i ].labelsNeighbors[n] = theCC[ theCC[ i ].labelsNeighbors[n] ].label;
  }

  if ( _TRACE_ ) {
    for ( i=0; i<dimz*dimxy; i++ ) {
      if ( cc_buffer[i] == 0 ) {
	output_buffer[i] = 0;
	continue;
      }
      output_buffer[i] = theCC[ cc_buffer[i] ].label;
    }
  } else {
    for ( i=0; i<dimz*dimxy; i++ ) {
      if ( cc_buffer[i] == 0 ) {
	output_buffer[i] = 0;
	continue;
      }
      output_buffer[i] = 255;
    }
  }


  if ( _TRACE_ ) {
    printf( "\n" );
    printf( " remaining curves = %d\n", nb_remaining_curves );
    printf( "\n" );
    for ( i=1; i<=nb_junctions+nb_curves; i++ ) {
      if ( theCC[ i ].size == 0 ) continue;
      if ( i <= nb_curves )
	printf( " COMPONENT #%-5d: ", i );
      else
	printf( " JUNCTION  #%-5d: ", i );
      printf( "size=%3d, ", theCC[ i ].size );
      printf( "label=%3d, ", theCC[ i ].label );
      printf( "(%d) = ", theCC[ i ].nbNeighbors );
      printf( "{ " );
      for ( n=0; n < theCC[ i ].nbNeighbors; n++ )
	printf( "%5d ", theCC[ i ].labelsNeighbors[n] );
      printf( "}" );
      printf( "\n" );
      if ( theCC[ i ].size > 0 ) {
	printf( "                   " );
	printf( "x=[%3d-%3d] y=[%3d-%3d] z=[%3d-%3d]",
		theCC[ i ].xmin, theCC[ i ].xmax, 
		theCC[ i ].ymin, theCC[ i ].ymax, 
		theCC[ i ].zmin, theCC[ i ].zmax );
      }
      printf( "\n" );
    }
    printf( "\n" );
  }




  free( theCC );
  free( cc_buffer );
  return( 1 );

}















/* destruction des composantes contenant des jonctions
   ou des composantes trop courtes
   
   l'image binaire contient des courbes (3D ou 2D)
   deja aminicies. On peut donc caracteriser
   en comptant le nombre de voisins :
   1 : point simple (fin de courbe)
   2 : point de courbe simple 
   2 ou + : point de jonction

*/

int remove_non_simple_components( unsigned char *input_buffer,
				  unsigned char *output_buffer,
				  int *theDim,
				  int minimal_length )
{
  char *proc = "remove_non_simple_components";
  unsigned short int *cc_buffer = NULL;
  int dimx = theDim[0];
  int dimy = theDim[1];
  int dimz = theDim[2];
  int dimxy = theDim[0]*theDim[1];
  int x, y, z, i, j, k, l;
  int iy, iz;

  int offset[27];
  int nb_neighbors;
  int nb_components;

  int *cc_flag = NULL;
  int *cc_extr = NULL;
  
  /* initialisation des offsets
   */
  offset[12] = -1;
  offset[13] =  0;
  offset[14] =  1;
  for ( i=-1; i<=1; i++ ) {
    offset[10+i] = offset[13+i] - dimx;
    offset[16+i] = offset[13+i] + dimx;
  }
  for ( i=-4; i<=4; i++ ) {
    offset[ 4+i] = offset[13+i] - dimx*dimy;
    offset[22+i] = offset[13+i] + dimx*dimy;
    
  }

  
  /* allocations des buffers auxiliaires
   */
  cc_buffer = (unsigned short int *)calloc( dimz*dimxy, sizeof( unsigned short int ) );
  if ( cc_buffer == NULL ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate auxiliary buffer\n", proc );
    return( -1 );
  }


  /* count connected components
   */
  nb_components = CountConnectedComponentsWithAllParams( input_buffer, UCHAR,
							 cc_buffer, USHORT,
							 theDim, 1.0,
							 26, /* connectivity */
							 minimal_length, 
							 0, /* all components */
							 0  /* no binary output */ );

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "%s: found %d connected components\n", proc, nb_components );
  }


  /* allocations and initializations of auxiliary arrays
   */
  cc_extr = cc_flag = (int *)calloc( 2*(nb_components+1), sizeof( int ) );
  if ( cc_flag == NULL ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate flag's auxiliary buffer\n", proc );
    free( cc_buffer );
    return( -1 );
  }
  cc_extr += (nb_components+1);
  

  
  for ( i=0, z=0; z<dimz; z++ ) {
    iz =  ( z == 0 || z == dimz-1 ) ? 0 : 1;
    for ( y=0; y<dimy; y++ ) {
      iy = ( iz == 0 || y == 0 || y == dimy-1 ) ? 0 : 1;
      for ( x=0; x<dimx; x++, i++ ) {
	
	if ( cc_buffer[i] == 0 ) continue;
	
	nb_neighbors = 0;
	
	if ( iy == 0 || x == 0 || x == dimx-1 ) {
	  for ( l = -1; l < 2; l++ ) 
	  for ( k = -1; k < 2; k++ ) 
	  for ( j = -1; j < 2; j++ ) {
	    if ( x+j >= 0 &&  x+j < dimx &&
		 y+k >= 0 &&  y+k < dimy &&
		 z+l >= 0 &&  z+l < dimz )
	      if ( cc_buffer[ (z+l)*dimxy + (y+k)*dimy + (x+j) ] ) nb_neighbors ++;
	  }
	} else {
	  for ( j=0; j<27; j++ )
	    if ( cc_buffer[ i + offset[j] ] ) nb_neighbors ++;
	}
	nb_neighbors --;

	/* admissible values for nb_neighbors are in [0..8]
	   nb_neighbors = 0 : isolated point 
	   nb_neighbors = 1 : end point
	   nb_neighbors = 2 : simple curve point
	   nb_neighbors > 2 : curves junction point
	*/

	if      ( nb_neighbors > 2 )  cc_flag[ cc_buffer[i] ] = 1;
	else if ( nb_neighbors == 1 ) cc_extr[ cc_buffer[i] ] ++;
      }
    }
  }

  cc_flag[0] = 1;
  
  for ( i=1; i<=nb_components; i++ ) {
    if ( cc_extr[i] != 2 ) cc_flag[i] = 1;
    if ( cc_flag[i] == 1 ) nb_components--;
  }
  
  fprintf( stderr, "%s: found %d valid connected components\n", proc, nb_components );


  for ( i=0; i<dimz*dimy*dimx; i++ ) {
    output_buffer[i] = ( cc_flag[ cc_buffer[i] ] == 1 ) ? 0 : 255;
  }
  

  free( cc_flag );
  free( cc_buffer );
  return( 1 );
}









/* extract parameters from isolated simple curves
 */


typeFiberParameter * compute_fibers_parameters( unsigned char *input_buffer,
						int *theDim,
						int minimal_length,
						int *nb_fibers )
{
  char *proc = "compute_fibers_parameters";
  unsigned short int *cc_buffer = NULL;
  int dimx = theDim[0];
  int dimy = theDim[1];
  int dimz = theDim[2];
  int dimxy = theDim[0]*theDim[1];
  int x, y, z, i, j, k, l;
  int iy, iz;

  int offset[27];
  int nb_neighbors;
  int nb_components;

  int *cc_flag = NULL;
  typeFiberParameter *tmp_param = NULL, *cc_param = NULL;


  *nb_fibers = 0;

  
  /* initialisation des offsets
   */
  offset[12] = -1;
  offset[13] =  0;
  offset[14] =  1;
  for ( i=-1; i<=1; i++ ) {
    offset[10+i] = offset[13+i] - dimx;
    offset[16+i] = offset[13+i] + dimx;
  }
  for ( i=-4; i<=4; i++ ) {
    offset[ 4+i] = offset[13+i] - dimx*dimy;
    offset[22+i] = offset[13+i] + dimx*dimy;
    
  }

  
  /* allocations des buffers auxiliaires
   */
  cc_buffer = (unsigned short int *)calloc( dimz*dimxy, sizeof( unsigned short int ) );
  if ( cc_buffer == NULL ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate auxiliary buffer\n", proc );
    return( NULL );
  }


  /* count connected components
   */
  nb_components = CountConnectedComponentsWithAllParams( input_buffer, UCHAR,
							 cc_buffer, USHORT,
							 theDim, 1.0,
							 26, /* connectivity */
							 minimal_length, 
							 0, /* all components */
							 0  /* no binary output */ );

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "%s: found %d connected components\n", proc, nb_components );
  }


  /* allocations and initializations of auxiliary arrays
   */
  cc_flag = (int *)calloc( (nb_components+1), sizeof( int ) );
  if ( cc_flag == NULL ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate flag's auxiliary buffer\n", proc );
    free( cc_buffer );
    return( NULL );
  }

  tmp_param = (typeFiberParameter *)malloc( (nb_components+1)*sizeof( typeFiberParameter ) );
  if ( tmp_param == NULL ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate parameter's auxiliary buffer\n", proc );
    free( cc_flag );
    free( cc_buffer );
    return( NULL );
  }

  for ( i=0; i<= nb_components; i++ ) {
    tmp_param[i].size = 0;
    tmp_param[i].extremity_1[0] = -1;
    tmp_param[i].extremity_1[1] = -1;
    tmp_param[i].extremity_1[2] = -1;
    tmp_param[i].extremity_2[0] = -1;
    tmp_param[i].extremity_2[1] = -1;
    tmp_param[i].extremity_2[2] = -1;
  }


  
  for ( i=0, z=0; z<dimz; z++ ) {
    iz =  ( z == 0 || z == dimz-1 ) ? 0 : 1;
    for ( y=0; y<dimy; y++ ) {
      iy = ( iz == 0 || y == 0 || y == dimy-1 ) ? 0 : 1;
      for ( x=0; x<dimx; x++, i++ ) {
	
	if ( cc_buffer[i] == 0 ) continue;
	
	nb_neighbors = 0;
	
	if ( iy == 0 || x == 0 || x == dimx-1 ) {
	  for ( l = -1; l < 2; l++ ) 
	  for ( k = -1; k < 2; k++ ) 
	  for ( j = -1; j < 2; j++ ) {
	    if ( x+j >= 0 &&  x+j < dimx &&
		 y+k >= 0 &&  y+k < dimy &&
		 z+l >= 0 &&  z+l < dimz )
	      if ( cc_buffer[ (z+l)*dimxy + (y+k)*dimy + (x+j) ] ) nb_neighbors ++;
	  }
	} else {
	  for ( j=0; j<27; j++ )
	    if ( cc_buffer[ i + offset[j] ] ) nb_neighbors ++;
	}
	nb_neighbors --;

	/* admissible values for nb_neighbors are in [0..8]
	   nb_neighbors = 0 : isolated point 
	   nb_neighbors = 1 : end point
	   nb_neighbors = 2 : simple curve point
	   nb_neighbors > 2 : curves junction point
	*/

	tmp_param[ cc_buffer[i] ].size ++;
	
	if ( nb_neighbors == 2 || nb_neighbors == 0 ) continue;
	
	if ( nb_neighbors > 2 ) {
	  cc_flag[ cc_buffer[i] ] = 1;
	  continue;
	}

	if ( tmp_param[ cc_buffer[i] ].extremity_1[0] < 0 ) {
	  tmp_param[ cc_buffer[i] ].extremity_1[0] = x;
	  tmp_param[ cc_buffer[i] ].extremity_1[1] = y;
	  tmp_param[ cc_buffer[i] ].extremity_1[2] = z;
	} else if ( tmp_param[ cc_buffer[i] ].extremity_2[0] < 0 ) {
	  tmp_param[ cc_buffer[i] ].extremity_2[0] = x;
	  tmp_param[ cc_buffer[i] ].extremity_2[1] = y;
	  tmp_param[ cc_buffer[i] ].extremity_2[2] = z;
	} else {
	  cc_flag[ cc_buffer[i] ] = 1;
	}
	

      }
    }
  }

  free( cc_buffer );


  *nb_fibers = nb_components;
  for ( i=1; i<=nb_components; i++ ) {
    if ( tmp_param[ i ].extremity_1[0] < 0 || tmp_param[ i ].extremity_2[0] < 0 )
      cc_flag[i] = 1;
    if ( cc_flag[i] == 1 ) (*nb_fibers)--;
  }
  fprintf( stderr, "%s: found %d valid connected components\n", proc, (*nb_fibers) );

  if ( (*nb_fibers) == nb_components ) {
    free( cc_flag );
    return( tmp_param );
  }


  
  cc_param = (typeFiberParameter *)malloc( ((*nb_fibers)+1)*sizeof( typeFiberParameter ) );
  if ( cc_param == NULL ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate second parameter's auxiliary buffer\n", proc );
    free( tmp_param );
    free( cc_flag );
    return( NULL );
  }

  cc_param[0] = tmp_param[0];
  for ( j=1, i=1; i<=nb_components; i++ ) {
    if ( cc_flag[i] == 1 ) continue;
    cc_param[j++] = tmp_param[i];
  }

  free( tmp_param );
  free( cc_flag );
  return( cc_param );
 
}
