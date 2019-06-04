/*************************************************************************
 * bal-cellProperties.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Mer 24 oct 2018 16:20:12 CEST
 *
 * ADDITIONS, CHANGES
 *
 */

/* cbrt(): Feature Test Macro Requirements for glibc (see feature_test_macros(7)):
 * _ISOC99_SOURCE || _POSIX_C_SOURCE >= 200112L
 *    || _XOPEN_SOURCE >= 500
 *    || // Since glibc 2.19: // _DEFAULT_SOURCE
 *    || // Glibc versions <= 2.19: // _BSD_SOURCE || _SVID_SOURCE
 */
#define _XOPEN_SOURCE 500

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <eigens.h>
#include <vtmalloc.h>

#include <surface-lindblad.h>
#include <surface-windreich.h>

#include <bal-cellProperties.h>



static int _verbose_ = 1;
static int _debug_ = 0;


void BAL_SetVerboseInBalCellProperties( int v )
{
  _verbose_ = v;
}

void BAL_IncrementVerboseInBalCellProperties(  )
{
  _verbose_ ++;
}

void BAL_DecrementVerboseInBalCellProperties(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}



void BAL_SetDebugInBalCellProperties( int d )
{
  _debug_ = d;
}

void BAL_IncrementDebugInBalCellProperties(  )
{
  _debug_ ++;
}

void BAL_DecrementDebugInBalCellProperties(  )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}





/************************************************************
 *
 *
 *
 ************************************************************/



void BAL_InitUpdateList( typeUpdateList *l )
{
  l->data = (int*)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;
}



void BAL_FreeUpdateList( typeUpdateList *l )
{
  if ( l->data != (int*)NULL ) vtfree( l->data );
  BAL_InitUpdateList( l );
}


static int _updates_to_be_allocated_ = 10;


int BAL_AddUpdateToList( typeUpdateList *l, int u )
{
  char *proc = "BAL_AddUpdateToList";
  int i, n;
  int *data = (int*)NULL;
  int s = l->n_allocated_data;

  /* is the update already existing?
   */
  for ( n=0; n<l->n_data; n++ ) {
    if ( l->data[n] == u ) {
      return( 1 );
    }
  }

  /* create update
   */
  if ( l->n_data == l->n_allocated_data ) {
    s += _updates_to_be_allocated_;
    data  = (int*)vtmalloc( s*sizeof(int), "data", proc );
    if ( data == (int*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    }
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(int) );
      vtfree( l->data );
    }
    for ( i=l->n_data; i<l->n_allocated_data; i++ )
        data[i] = 0;
    l->n_allocated_data = s;
    l->data = data;
  }

  l->data[l->n_data] = u;
  l->n_data ++;

  return( 1 );
}



void BAL_FprintfUpdateList( FILE *f, typeUpdateList *l )
{
  int i;

  fprintf( f, "[" );
  for ( i=0; i<l->n_data; i++ ) {
    fprintf( f, "%d", l->data[i] );
    if ( i < l->n_data-1 )
      fprintf( f, ", " );
  }
  fprintf( f, "]" );
  fprintf( f, "\n" );
}



int _isUpdateInList( typeUpdateList *l, int u )
{
  int i;
  if ( l == (typeUpdateList*)NULL) return( 0 );
  if ( l->n_data <= 0 ) return( 0 );
  for ( i=0; i<l->n_data; i++ ) {
    if ( l->data[i] == u ) return( 1 );
  }
  return( 0 );
}





/************************************************************
 *
 *
 *
 ************************************************************/





static char *_Property( enumProperty property )
{
  switch( property ) {
  default : return( "Unknown" );
  case _NONE_ : return( "_NONE_" );
  case _LINEAGE_ : return( "_LINEAGE_" );
  case _FORWARD_NEIGHBORS_ : return( "_FORWARD_NEIGHBORS_" );
  case _BACKWARD_NEIGHBORS_ : return( "_BACKWARD_NEIGHBORS_" );
  case _H_MIN_ : return( "_H_MIN_" );
  case _VOLUME_ : return( "_VOLUME_" );
  case _SURFACE_ : return( "_SURFACE_" );
  case _SIGMA_ : return( "_SIGMA_" );
  case _LABEL_IN_TIME_ : return( "_LABEL_IN_TIME_" );
  case _BARYCENTER_ : return( "_BARYCENTER_" );
  case _PRINCIPAL_VALUE_ : return( "_PRINCIPAL_VALUE_" );
  case _PRINCIPAL_VECTOR_ : return( "_PRINCIPAL_VECTOR_" );
  case _FATE_ : return( "_FATE_" );
  case _FATE2_ : return( "_FATE2_" );
  case _FATE3_ : return( "_FATE3_" );
  case _FATE4_ : return( "_FATE4_" );
  case _ALL_CELLS_ : return( "_ALL_CELLS_" );
  case _NAME_ : return( "_NAME_" );
  case _HISTORY_ : return( "_HISTORY_" );
  case _CONTACT_SURFACE_ : return( "_CONTACT_SURFACE_" );
  }
  return( "this should not occur" );
}



static char *_PropertyName( enumProperty property )
{
  switch( property ) {
  default : return( "Unknown" );
  case _NONE_ : return( "_NONE_" );
  case _LINEAGE_ : return( "cell lineage" );
  case _FORWARD_NEIGHBORS_ : return( "cell forward neighbors" );
  case _BACKWARD_NEIGHBORS_ : return( "cell backward neighbors" );
  case _H_MIN_ : return( "cell h-min values" );
  case _VOLUME_ : return( "cell volume" );
  case _SURFACE_ : return( "cell surface" );
  case _SIGMA_ : return( "cell sigma values" );
  case _LABEL_IN_TIME_ : return( "cell label in time" );
  case _BARYCENTER_ : return( "cell barycenter" );
  case _PRINCIPAL_VALUE_ : return( "cell principal values" );
  case _PRINCIPAL_VECTOR_ : return( "cell principal vectors" );
  case _FATE_ : return( "cell fate" );
  case _FATE2_ : return( "cell fate #2" );
  case _FATE3_ : return( "cell fate #3" );
  case _FATE4_ : return( "cell fate #4" );
  case _ALL_CELLS_ : return( "all cell labels" );
  case _NAME_ : return( "cell name" );
  case _HISTORY_ : return( "cell history" );
  case _CONTACT_SURFACE_ : return( "cell contact surfaces" );
  }
  return( "this should not occur" );
}



void BAL_InitPropertyList( typePropertyList *o )
{
  int i;
  o->n_allocated_data = _MAX_PROPERTIES_;
  o->n_data = 0;
  for ( i=0; i<_MAX_PROPERTIES_; i++ )
      o->data[i] = _NONE_;
}



void BAL_DefaultPropertyList( typePropertyList *o )
{
  BAL_InitPropertyList( o );
  o->data[o->n_data++] = _VOLUME_;
  o->data[o->n_data++] = _SURFACE_;
  o->data[o->n_data++] = _COMPACTNESS_;
  o->data[o->n_data++] = _BARYCENTER_;
  o->data[o->n_data++] = _PRINCIPAL_VALUE_;
  o->data[o->n_data++] = _PRINCIPAL_VECTOR_;
  o->data[o->n_data++] = _CONTACT_SURFACE_;
  o->data[o->n_data++] = _ALL_CELLS_;
  o->data[o->n_data++] = _LINEAGE_;
}



static int _addPropertyToList( typePropertyList *l, enumProperty p )
{
  char *proc = "_addPropertyToList";
  int i;

  for ( i=0; i<l->n_data; i++ ) {
    if ( l->data[i] == p ) return( 1 );
  }
  if ( l->n_data < l->n_allocated_data ) {
    l->data[ l->n_data ] = p;
    l->n_data ++;
    return( 1 );
  }
  if ( _verbose_ )
    fprintf( stderr, "%s: unable to add property to list\n", proc );
  return( -1 );
}



void _removePropertyFromList( typePropertyList *l, enumProperty p )
{
  int i, j;
  for ( i=0, j=0; i<l->n_data; i++ ) {
    if ( l->data[ i ] != p ) {
      l->data[ j++ ] = l->data[ i ];
    }
  }
  l->n_data = j;
}



void BAL_FprintfPropertyList( FILE *f, typePropertyList *o, char *str )
{
  int i, j, l = 5;

  if ( str != (char*)NULL ) fprintf( f, "%s", str );
  else fprintf( f, "* Property list" );
  fprintf( f, "\n" );

  if ( 0 ) fprintf( f, "   - n_allocated_data = %d\n", o->n_allocated_data );
  if ( 0 ) fprintf( f, "   - n_data = %d\n", o->n_data );

  if ( str != (char*)NULL ) l = strlen(str);
  for ( i=0; i<o->n_data; i++ ) {
    if ( o->data[i] == _NONE_ ) continue;
    for (j=0; j<l; j++ ) fprintf( f, " " );
    fprintf( f, "data[%2d] = '%s'\n", i,  _PropertyName(o->data[i]) );
  }
}



int _isPropertyInList( typePropertyList *o, enumProperty p )
{
  int i;
  if ( o == (typePropertyList*)NULL) return( 1 );
  if ( o->n_data <= 0 ) return( 1 );
  for ( i=0; i<o->n_data; i++ ) {
    if ( o->data[i] == p ) return( 1 );
  }
  return( 0 );
}



void BAL_GetWritePropertyList( typePropertyList *writeList, typePropertyList *readList, typePropertyList *outputList )
{
  char *proc = "BAL_GetWritePropertyList";
  int i;

  BAL_InitPropertyList( writeList );

  if ( readList == (typePropertyList *)NULL || readList->n_data <= 0 ) {
    for ( i=0; i<outputList->n_data; i++ )
      writeList->data[i] = outputList->data[i];
    writeList->n_data = outputList->n_data;
  }
  else if ( outputList == (typePropertyList *)NULL || outputList->n_data <= 0 ) {
    for ( i=0; i<readList->n_data; i++ )
      writeList->data[i] = readList->data[i];
    writeList->n_data = readList->n_data;
  }
  else {
    for ( i=0; i<outputList->n_data; i++ ) {
      if ( _isPropertyInList( readList, outputList->data[i] ) == 1 )
          writeList->data[ writeList->n_data++ ] = outputList->data[i];
      else {
        if ( _verbose_ ) {
          fprintf( stderr, "%s: property '%s' was not in read properties, discard it\n", proc, _PropertyName(outputList->data[i]) );
        }
      }
    }
  }
}



/************************************************************
 *
 *
 *
 ************************************************************/



static void _cellLocalId( int id, int *cell, int *time )
{
  *time = id / 10000;
  *cell = id % 10000;
}



static int _cellGlobalId( int cell, int time )
{
    return( time * 10000 + cell );
}



/************************************************************
 *
 *
 *
 ************************************************************/



static void _initSpatialNeighbor( typeSpatialNeighbor *n )
{
  n->label = -1;
  n->surface = 0;
}


static void _initSpatialNeighborList( typeSpatialNeighborList *l )
{
  l->data = (typeSpatialNeighbor*)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;
}



static void _freeSpatialNeighborList( typeSpatialNeighborList *l )
{
  if ( l->data != (typeSpatialNeighbor*)NULL ) vtfree( l->data );
  _initSpatialNeighborList( l );
}



static int _compareTemporalNeighbor ( const void * a, const void * b )
{
    typeTemporalNeighbor *da = (typeTemporalNeighbor*)a;
    typeTemporalNeighbor *db = (typeTemporalNeighbor*)b;

    if ( da->npts > db->npts ) return( -1 );
    if ( da->npts < db->npts ) return( 1 );
    return( 0 );
}



static void _initTemporalNeighbor( typeTemporalNeighbor *d )
{
  d->label = -1;
  d->time = -1;
  d->globalId = -1;
  d->npts = 0;
}



static void _initTemporalNeighborList( typeTemporalNeighborList *l )
{
  l->data = (typeTemporalNeighbor*)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;
}



static void _freeTemporalNeighborList( typeTemporalNeighborList *l )
{
  if ( l->data != (typeTemporalNeighbor*)NULL ) vtfree( l->data );
  _initTemporalNeighborList( l );
}



static void _initCell( typeCell *c )
{
  c->ptmin[0] = c->ptmin[1] = c->ptmin[2] = 0;
  c->ptmax[0] = c->ptmax[1] = c->ptmax[2] = 0;

  c->npts = 0;
  c->surface = 0.0;
  c->compactness = 0.0;

  c->voxel_barycenter[0] = 0.0;
  c->voxel_barycenter[1] = 0.0;
  c->voxel_barycenter[2] = 0.0;

  c->voxel_covariance_matrix[0] = 0.0;
  c->voxel_covariance_matrix[1] = 0.0;
  c->voxel_covariance_matrix[2] = 0.0;
  c->voxel_covariance_matrix[3] = 0.0;
  c->voxel_covariance_matrix[4] = 0.0;
  c->voxel_covariance_matrix[5] = 0.0;
  c->voxel_covariance_matrix[6] = 0.0;
  c->voxel_covariance_matrix[7] = 0.0;
  c->voxel_covariance_matrix[8] = 0.0;

  c->voxel_eigenvalues[0] = 0.0;
  c->voxel_eigenvalues[1] = 0.0;
  c->voxel_eigenvalues[2] = 0.0;

  c->voxel_eigenvectors[0] = 0.0;
  c->voxel_eigenvectors[1] = 0.0;
  c->voxel_eigenvectors[2] = 0.0;
  c->voxel_eigenvectors[3] = 0.0;
  c->voxel_eigenvectors[4] = 0.0;
  c->voxel_eigenvectors[5] = 0.0;
  c->voxel_eigenvectors[6] = 0.0;
  c->voxel_eigenvectors[7] = 0.0;
  c->voxel_eigenvectors[8] = 0.0;

  _initSpatialNeighborList( &(c->neighbors) );

  _initTemporalNeighborList( &(c->forward) );
  _initTemporalNeighborList( &(c->backward) );
  _initTemporalNeighborList( &(c->lineage) );
  _initTemporalNeighborList( &(c->history) );

  c->fate = (char*)NULL;
  c->fate2 = (char*)NULL;
  c->fate3 = (char*)NULL;
  c->fate4 = (char*)NULL;

  c->name = (char*)NULL;
}



static void _freeCell( typeCell *c )
{
  _freeSpatialNeighborList( &(c->neighbors) );

  _freeTemporalNeighborList( &(c->forward) );
  _freeTemporalNeighborList( &(c->backward) );
  _freeTemporalNeighborList( &(c->lineage) );
  _freeTemporalNeighborList( &(c->history) );

  if ( c->fate != (char*)NULL ) vtfree( c->fate );
  if ( c->fate2 != (char*)NULL ) vtfree( c->fate2 );
  if ( c->fate3 != (char*)NULL ) vtfree( c->fate3 );
  if ( c->fate4 != (char*)NULL ) vtfree( c->fate4 );

  if ( c->name != (char*)NULL ) vtfree( c->name );

  _initCell( c );
}



static void _initCellPartial( typeCell *c, typePropertyList *theProperties )
{

  c->ptmin[0] = c->ptmin[1] = c->ptmin[2] = 0;
  c->ptmax[0] = c->ptmax[1] = c->ptmax[2] = 0;

  if ( _isPropertyInList( theProperties, _VOLUME_ ) ) {
    c->npts = 0;
  }

  if ( _isPropertyInList( theProperties, _SURFACE_ ) ) {
    c->surface = 0;
  }

  if ( _isPropertyInList( theProperties, _COMPACTNESS_ ) ) {
    c->compactness = 0;
  }


  if ( _isPropertyInList( theProperties, _BARYCENTER_ ) ) {
    c->npts = 0;

    c->voxel_barycenter[0] = 0.0;
    c->voxel_barycenter[1] = 0.0;
    c->voxel_barycenter[2] = 0.0;
  }

  if ( _isPropertyInList( theProperties, _PRINCIPAL_VALUE_ )
       || _isPropertyInList( theProperties, _PRINCIPAL_VECTOR_ ) ) {
    c->npts = 0;

    c->voxel_barycenter[0] = 0.0;
    c->voxel_barycenter[1] = 0.0;
    c->voxel_barycenter[2] = 0.0;

    c->voxel_covariance_matrix[0] = 0.0;
    c->voxel_covariance_matrix[1] = 0.0;
    c->voxel_covariance_matrix[2] = 0.0;
    c->voxel_covariance_matrix[3] = 0.0;
    c->voxel_covariance_matrix[4] = 0.0;
    c->voxel_covariance_matrix[5] = 0.0;
    c->voxel_covariance_matrix[6] = 0.0;
    c->voxel_covariance_matrix[7] = 0.0;
    c->voxel_covariance_matrix[8] = 0.0;

    c->voxel_eigenvalues[0] = 0.0;
    c->voxel_eigenvalues[1] = 0.0;
    c->voxel_eigenvalues[2] = 0.0;

    c->voxel_eigenvectors[0] = 0.0;
    c->voxel_eigenvectors[1] = 0.0;
    c->voxel_eigenvectors[2] = 0.0;
    c->voxel_eigenvectors[3] = 0.0;
    c->voxel_eigenvectors[4] = 0.0;
    c->voxel_eigenvectors[5] = 0.0;
    c->voxel_eigenvectors[6] = 0.0;
    c->voxel_eigenvectors[7] = 0.0;
    c->voxel_eigenvectors[8] = 0.0;
 }

  if ( _isPropertyInList( theProperties, _CONTACT_SURFACE_ ) ) {
    _freeSpatialNeighborList( &(c->neighbors) );
  }

  if ( _isPropertyInList( theProperties, _LINEAGE_ ) ) {
    _freeTemporalNeighborList( &(c->forward) );
    _freeTemporalNeighborList( &(c->backward) );
    _freeTemporalNeighborList( &(c->lineage) );
  }
  if ( _isPropertyInList( theProperties, _HISTORY_ ) ) {
    _freeTemporalNeighborList( &(c->history) );
  }
}





static int _spatial_neighbors_to_be_allocated_ = 10;


static int _addSpatialNeighborToCell( typeCell *c, int neighbor, float contrib )
{
  char *proc = "_addSpatialNeighborToCell";
  int i, n;
  typeSpatialNeighborList *l = &(c->neighbors);
  typeSpatialNeighbor *data = (typeSpatialNeighbor*)NULL;
  int s = l->n_allocated_data;

  /* is the neighbor already existing?
   */
  if ( 0 ) fprintf( stderr, "test %2d vs %2d/%2d neighbors\n", neighbor, l->n_data, l->n_allocated_data );
  for ( n=0; n<l->n_data; n++ ) {
    if ( l->data[n].label == neighbor ) {
      l->data[n].surface += contrib;
      return( 1 );
    }
  }

  /* create neighbor
   */
  if ( l->n_data == l->n_allocated_data ) {
    s += _spatial_neighbors_to_be_allocated_;
    data  = (typeSpatialNeighbor*)vtmalloc( s*sizeof(typeSpatialNeighbor), "data", proc );
    if ( data == (typeSpatialNeighbor*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    }
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(typeSpatialNeighbor) );
      vtfree( l->data );
    }
    for ( i=l->n_data; i<l->n_allocated_data; i++ )
        _initSpatialNeighbor( &(data[i]) );
    l->n_allocated_data = s;
    l->data = data;
  }

  l->data[l->n_data].label = neighbor;
  l->data[l->n_data].surface = contrib;
  l->n_data ++;

  return( 1 );
}





static int _temporal_neighbors_to_be_allocated_ = 10;


static int _addTemporalNeighborToList( typeTemporalNeighborList *l, int cellGlobalId, int contrib )
{
  char *proc = "_addTemporalNeighborToCell";
  int i, n;
  typeTemporalNeighbor *data = (typeTemporalNeighbor*)NULL;
  int s = l->n_allocated_data;
  int cellLabel, cellTime;

  _cellLocalId( cellGlobalId, &cellLabel, &cellTime );

  /* is the idcell already existing?
   */
  for ( n=0; n<l->n_data; n++ ) {
    if ( l->data[n].globalId == cellGlobalId ) {
      l->data[n].npts += contrib;
      return( 1 );
    }
  }

  /* create neighbor
   */
  if ( l->n_data == l->n_allocated_data ) {
    s += _temporal_neighbors_to_be_allocated_;
    data  = (typeTemporalNeighbor*)vtmalloc( s*sizeof(typeTemporalNeighbor), "data", proc );
    if ( data == (typeTemporalNeighbor*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    }
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(typeTemporalNeighbor) );
      vtfree( l->data );
    }
    for ( i=l->n_data; i<l->n_allocated_data; i++ )
        _initTemporalNeighbor( &(data[i]) );
    l->n_allocated_data = s;
    l->data = data;
  }

  l->data[l->n_data].label = cellLabel;
  l->data[l->n_data].time = cellTime;
  l->data[l->n_data].globalId = cellGlobalId;
  l->data[l->n_data].npts = contrib;
  l->n_data ++;

  return( 1 );
}



static int _addTemporalNeighborToCellLineage( typeCell *c, int cellGlobalId, int contrib )
{
  char *proc = "_addTemporalNeighborToCellLineage";
  if ( _addTemporalNeighborToList( &(c->lineage), cellGlobalId, contrib ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to add lineage %d to cell\n", proc, cellGlobalId );
    return( -1 );
  }
  return( 1 );
}



static int _addTemporalNeighborToCellHistory( typeCell *c, int cellGlobalId, int contrib )
{
  char *proc = "_addTemporalNeighborToCellHistory";

  if ( _addTemporalNeighborToList( &(c->lineage), cellGlobalId, contrib ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to add history %d to cell\n", proc, cellGlobalId );
    return( -1 );
  }
  return( 1 );
}



static int _addBackwardNeighborToCell( typeCell *c, int cellGlobalId, int contrib )
{
  char *proc = "_addBackwardNeighborToCell";
  if ( _addTemporalNeighborToList( &(c->backward), cellGlobalId, contrib ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to add mother %d to cell\n", proc, cellGlobalId );
    return( -1 );
  }
  return( 1 );
}





static int _addForwardNeighborToCell( typeCell *c, int cellGlobalId, int contrib )
{
  char *proc = "_addForwardNeighborToCell";
  if ( _addTemporalNeighborToList( &(c->forward), cellGlobalId, contrib ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to add daughter %d to cell\n", proc, cellGlobalId );
    return( -1 );
  }
  return( 1 );
}





/************************************************************
 *
 * image-based properties
 *
 ************************************************************/





void BAL_InitCellImage( typeCellImage *l )
{
  l->data = (typeCell*)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;
  l->previous_acquisition_time = 0;
  l->acquisition_time = 0;
  l->next_acquisition_time = 0;
}



void BAL_FreeCellImage( typeCellImage *l )
{
  int i;

  for ( i=0; i<l->n_data; i++ ) {
      _freeCell( &(l->data[i]) );
  }

  if ( l->data != (typeCell*)NULL ) vtfree( l->data );
  BAL_InitCellImage( l );
}



int BAL_AllocCellImage( typeCellImage *l , int n )
{
  char *proc = "BAL_AllocCellImage";
  int i;

  l->data = (typeCell*)vtmalloc( n * sizeof(typeCell), "l->data", proc );

  if ( l->data == (typeCell*)NULL ) {
      if ( _verbose_ )
          fprintf( stderr, "%s: unable to allocate data\n", proc );
    return( -1 );
  }
  l->n_data = n;
  l->n_allocated_data = n;

  for ( i=0; i<l->n_data; i++ )
      _initCell( &(l->data[i]) );

  return( 1 );
}



int BAL_InitCellImageFromImage( typeCellImage *theList,  bal_image *theIm )
{
  char *proc = "BAL_InitCellImageFromImage";
  size_t i, v;
  int ncells;

  /* number of cells
   */
  v = theIm->ncols * theIm->nrows * theIm->nplanes;

#define _NUMBER_OF_CELLS( TYPE ) {                \
  TYPE *theBuf = (TYPE*)theIm->data;              \
  for ( i=0; i<v; i++ ) {                         \
    if ( ncells < theBuf[i] ) ncells = theBuf[i]; \
  }                                               \
}

  ncells = 0;
  switch ( theIm->type ) {
  default :
    fprintf( stderr, "%s: type not handled in switch\n", proc );
    return( -1 );
  case UCHAR :
    _NUMBER_OF_CELLS( u8 );
    break;
  case SSHORT :
    _NUMBER_OF_CELLS( s16 );
    break;
  case USHORT :
    _NUMBER_OF_CELLS( u16 );
    break;
  }

  if ( _debug_ >= 2 )
    fprintf( stderr, "%s: found %d cells\n", proc, ncells );

  if ( ncells <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: NULL image\n", proc );
    return( 0 );
  }

  if ( theList->n_data > 0 ) {
    if ( theList->n_data < ncells+1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: already allocated list with too few cells\n", proc );
      return( -1 );
    }
    else {
      /* no allocation to be done
       * assume this is all right
       */
      return( ncells );
    }
  }

  /* allocate list of cells
   */
  if ( BAL_AllocCellImage( theList, ncells+1 ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }

  return( ncells );
}





static int _cells_to_be_allocated_ = 50;


static typeCell* _getCellFromCellImage( typeCellImage *l, int cell )
{
  char *proc = "_getCellFromCellImage";
  int i;
  typeCell *data = (typeCell*)NULL;
  int s = l->n_allocated_data;

  if ( cell < l->n_data )
      return( &(l->data[cell]) );

  if ( cell < l->n_allocated_data ) {
    l->n_data = cell+1;
    return( &(l->data[cell]) );
  }

  /* create cell
   */
  while ( cell >= s )
      s += _cells_to_be_allocated_;

  data  = (typeCell*)vtmalloc( s*sizeof(typeCell), "data", proc );
  if ( data == (typeCell*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( (typeCell*)NULL );
  }
  if ( l->n_allocated_data > 0 ) {
    (void)memcpy( data, l->data, l->n_allocated_data*sizeof(typeCell) );
    vtfree( l->data );
  }
  for ( i=l->n_data; i<s; i++ ) {
      _initCell( &(data[i]) );
  }
  l->n_allocated_data = s;
  l->n_data = cell+1;
  l->data = data;

  return( &(l->data[cell]) );
}





/************************************************************
 *
 * sequence-based properties
 *
 ************************************************************/





void BAL_InitCellSequence( typeCellSequence *l )
{
  l->data = (typeCellImage*)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;
}



void BAL_FreeCellSequence( typeCellSequence *l )
{
  int i;

  for ( i=0; i<l->n_data; i++ ) {
      BAL_FreeCellImage( &(l->data[i]) );
  }

  if ( l->data != (typeCellImage*)NULL ) vtfree( l->data );
  BAL_InitCellSequence( l );
}



int BAL_AllocCellSequence( typeCellSequence *l , int n )
{
  char *proc = "BAL_AllocCellSequence";
  int i;

  l->data = (typeCellImage*)vtmalloc( n * sizeof(typeCellImage), "l->data", proc );

  if ( l->data == (typeCellImage*)NULL ) {
      if ( _verbose_ )
          fprintf( stderr, "%s: unable to allocate data\n", proc );
    return( -1 );
  }
  l->n_data = n;
  l->n_allocated_data = n;

  for ( i=0; i<l->n_data; i++ )
      BAL_InitCellImage( &(l->data[i]) );

  return( 1 );
}





static int _cell_images_to_be_allocated_ = 10;

static typeCellImage* _getCellImageFromSequence( typeCellSequence *l, int acquisition )
{
  char *proc = "_getCellImageFromSequence";
  int i;
  typeCellImage *data = (typeCellImage*)NULL;
  int s = l->n_allocated_data;

  if ( acquisition < l->n_data )
      return( &(l->data[acquisition]) );

  if ( acquisition < l->n_allocated_data ) {
    l->n_data = acquisition+1;
    return( &(l->data[acquisition]) );
  }

  /* create cell image
   */
  while ( acquisition >= s ) {
      s += _cell_images_to_be_allocated_;
  }

  data  = (typeCellImage*)vtmalloc( s*sizeof(typeCellImage), "data", proc );
  if ( data == (typeCellImage*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( (typeCellImage*)NULL );
  }
  if ( l->n_allocated_data > 0 ) {
    (void)memcpy( data, l->data, l->n_allocated_data*sizeof(typeCellImage) );
    vtfree( l->data );
  }
  for ( i=l->n_data; i<s; i++ ) {
      BAL_InitCellImage( &(data[i]) );
      data[i].previous_acquisition_time = i-1;
      data[i].acquisition_time = i;
      data[i].next_acquisition_time = i+1;
  }
  l->n_allocated_data = s;
  l->n_data = acquisition+1;
  l->data = data;

  return( &(l->data[acquisition]) );
}





int BAL_InitCellSequencePartial( typeCellSequence *l,
                                 typePropertyList *theProperties, typeUpdateList *theIndexes )
{
  char *proc = "BAL_InitCellSequencePartial";
  typePropertyList lineageProperty;
  typePropertyList lineageMotherProperty;
  int i, j;
  typeCellImage *cellImage;

  BAL_InitPropertyList( &lineageProperty );
  BAL_InitPropertyList( &lineageMotherProperty );
  if ( _isPropertyInList( theProperties, _LINEAGE_ ) ) {
    lineageProperty.data[0] = _LINEAGE_;
    lineageProperty.n_data = 1;
  }

  for ( i=0; i<theIndexes->n_data; i++ ) {

    if ( theIndexes->data[i] < 0 || theIndexes->data[i] >= l->n_data ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: index %d out of sequence bounds, skip it\n", proc, theIndexes->data[i] );
      continue;
    }

    /* image #i
     */
    cellImage = &(l->data[ theIndexes->data[i] ]);
    if ( cellImage->n_data > 0 ) {
      if ( _verbose_ >= 2 )
        fprintf( stderr, "%s: partial initialisation of image #%d properties\n",
                 proc, theIndexes->data[i] );
      for ( j=0; j<cellImage->n_data; j++ ) {
        _initCellPartial( &(cellImage->data[j]), theProperties );
      }
    }

    /* image #i-1
     */
    if ( lineageProperty.n_data > 0 ) {
      if ( theIndexes->data[i] > 0 && _isUpdateInList( theIndexes, theIndexes->data[i]-1 ) == 0 ) {
        cellImage = &(l->data[ theIndexes->data[i]-1 ]);
        if ( cellImage->n_data > 0 ) {
          if ( _verbose_ >= 2 )
            fprintf( stderr, "%s: initialisation of image #%d forward lineage\n",
                     proc, theIndexes->data[i]-1 );
          for ( j=0; j<cellImage->n_data; j++ )
              _initCellPartial( &(cellImage->data[j]), &lineageProperty );
        }
      }
    }

  }

  return( 1 );
}





/************************************************************
 *
 * some writing facilities (for debug purposes)
 *
 ************************************************************/





void _fprintfCell( FILE *f, typeCell *c )
{
  int n;
  if ( 0 ) {
    fprintf( f, "box=[(%3d, %3d, %3d),(%3d, %3d, %3d)] ", c->ptmin[0],
            c->ptmin[1], c->ptmin[2], c->ptmax[0], c->ptmax[1], c->ptmax[2] );
  }
  fprintf( f, "npts=%7d ", c->npts );
  fprintf( f, "surface=%8.2f ", c->surface );
  if ( 0 ) {
      fprintf( f, "center=(%5.1f, %5.1f, %5.1f) ", c->voxel_barycenter[0],
          c->voxel_barycenter[1], c->voxel_barycenter[2] );
  }
  fprintf( f, "\n" );
  if ( 0 ) {
      fprintf( f, "covariance = %8.2f %8.2f %8.2f\n", c->voxel_covariance_matrix[0],
              c->voxel_covariance_matrix[1], c->voxel_covariance_matrix[2] );
      fprintf( f, "             %8.2f %8.2f %8.2f\n", c->voxel_covariance_matrix[3],
              c->voxel_covariance_matrix[4], c->voxel_covariance_matrix[5] );
      fprintf( f, "             %8.2f %8.2f %8.2f\n", c->voxel_covariance_matrix[6],
              c->voxel_covariance_matrix[7], c->voxel_covariance_matrix[8] );
      fprintf( f, "eigenvalues = (%8.2f, %8.2f, %8.2f) ", c->voxel_eigenvalues[0],
              c->voxel_eigenvalues[1], c->voxel_eigenvalues[2] );
      fprintf( f, "\n" );
      fprintf( f, "eigenvectors = %8.2f %8.2f %8.2f\n", c->voxel_eigenvectors[0],
              c->voxel_eigenvectors[1], c->voxel_eigenvectors[2] );
      fprintf( f, "               %8.2f %8.2f %8.2f\n", c->voxel_eigenvectors[3],
              c->voxel_eigenvectors[4], c->voxel_eigenvectors[5] );
      fprintf( f, "               %8.2f %8.2f %8.2f\n", c->voxel_eigenvectors[6],
              c->voxel_eigenvectors[7], c->voxel_eigenvectors[8] );
  }
  if ( 0 ) {
      if ( c->neighbors.n_data > 0 ) {
          fprintf( f, "  " );
          for ( n=0; n<c->neighbors.n_data; n++ ) {
            fprintf( f, "neigh[%3d]=%6f", c->neighbors.data[n].label, c->neighbors.data[n].surface );
            if ( n < c->neighbors.n_data - 1 )
              fprintf( f, ", " );
          }
          fprintf( f, " \n" );
        }
  }
  if ( 1 ) {
        if ( c->backward.n_data > 0 ) {
          fprintf( f, "  " );
          for ( n=0 ; n<c->backward.n_data; n++ ) {
            fprintf( f, "back[%3d]=%6d", c->backward.data[n].label, c->backward.data[n].npts );
            if ( n < c->backward.n_data - 1 )
              fprintf( f, ", " );
          }
          fprintf( f, " \n" );
        }
  }
  if ( 1 ) {
      if ( c->forward.n_data > 0 ) {
        fprintf( f, "  " );
        for ( n=0 ; n<c->forward.n_data; n++ ) {
          fprintf( f, "forw[%3d]=%6d", c->forward.data[n].label, c->forward.data[n].npts );
          if ( n < c->forward.n_data - 1 )
            fprintf( f, ", " );
        }
        fprintf( f, " \n" );
      }
  }
  if ( 1 ) {
        if ( c->lineage.n_data > 0 ) {
          fprintf( f, "  " );
          for ( n=0 ; n<c->lineage.n_data; n++ ) {
            fprintf( f, "lineage[%3d]=%6d", c->lineage.data[n].label, c->lineage.data[n].npts );
            if ( n < c->lineage.n_data - 1 )
              fprintf( f, ", " );
          }
          fprintf( f, " \n" );
        }
  }

}



void BAL_FprintfCellImage( FILE *f, typeCellImage *l )
{
  int i;
  for ( i=0; i<l->n_data; i++ ) {
    if ( l->data[i].npts > 0 || l->data[i].surface > 0 ) {
      fprintf( f, "--------------------------------------------------\n" );
      fprintf( f, "#%3d cell %3d ", l->acquisition_time, i );
      _fprintfCell( f, &(l->data[i]) );
      if ( 0 ) fprintf( f, "\n" );
    }
  }
  fprintf( f, "--------------------------------------------------\n" );
}



void BAL_FprintfCellSequence( FILE *f, typeCellSequence *l )
{
  int i;
  for ( i=0; i<l->n_data; i++ ) {
    if ( l->data[i].n_data <= 0 ) continue;
    fprintf( f, "==== image %3d ===================================\n", l->data[i].acquisition_time );
    BAL_FprintfCellImage( f, &(l->data[i]) );
    fprintf( f, "==================================================\n" );
  }
}





/************************************************************
 *
 *
 *
 ************************************************************/





static typeCell* _getCellFromSequence( typeCellSequence *l, int id )
{
  char *proc = "_getCellFromSequence";
  int cellid, acquisition_time;
  typeCellImage *cellImage;
  typeCell *cell;

  _cellLocalId( id, &cellid, &acquisition_time );

  cellImage = _getCellImageFromSequence( l, acquisition_time );
  if ( cellImage == (typeCellImage*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to get cell image #%d from sequence\n",
               proc, acquisition_time );
    return( (typeCell*)NULL );
  }

  cell = _getCellFromCellImage( cellImage, cellid );
  if ( cell == (typeCell*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to get cell #%d from cell image\n",
               proc, cellid );
    return( (typeCell*)NULL );
  }

  return( cell );
}





/************************************************************
 *
 *
 *
 ************************************************************/





void BAL_FprintfEnumSurfaceEstimation( FILE *f, enumSurfaceEstimation s )
{
    switch( s ) {
    default :
        fprintf( f, "unknown surface estimation type" );
        break;
    case _4EDGES_ :
        fprintf( f, "_4EDGES_" );
        break;
    case _OUTER_4NEIGHBORS_ :
        fprintf( f, "_OUTER_4NEIGHBORS_" );
        break;
    case _OUTER_6NEIGHBORS_ :
        fprintf( f, "_OUTER_6NEIGHBORS_" );
        break;
    case _LINDBLAD_ :
        fprintf( f, "_LINDBLAD_" );
        break;
    case _WINDREICH_ :
        fprintf( f, "_WINDREICH_" );
        break;
    }
}



static int _compareint ( const void * a, const void * b )
{
    int *da = (int*)a;
    int *db = (int*)b;

    if ( (*da) > (*db) ) return( -1 );
    if ( (*da) < (*db) ) return( 1 );
    return( 0 );
}



static int _Outer6NeighborsSurfaceEstimation( bal_image *theIm, typeCellImage *theList )
{
  char *proc = "_Outer6NeighborsSurfaceEstimation";
  size_t x, y, z;
  int label, n, neighbors[6], nneighbors;

  /* compute cell properties
   */

  /* the contact surface of cell #cell with cell #c is defined as the
   * number of voxel of cell #c that are 6-neighbors of cell #cell
   *
   * this is an outer surface contact.
   *
   * this is computed by examinating each point, see whether it is 6-neighors,
   * and, if so, count it as contact point for its neighbors
   *
   * if we want to count background as contact surface for cells,
   * we have to consider background points
   * so do not skip them
   *
   * if there are neighboring labels, this point count in the contact surface
   * recall: historical way of counting contact surfaces of cell c is to dilate
   * cell c with a 6-structuring element and to count other labels
   *
   * Hence, if [z][y][x] is a 6-neighbor of an other label, it counts as 1 more contact point,
   *
   * If there are more than 2 neighbors, they are sorted to make sure that they are not processed twice
   *
   */

#define _OUTSIDE6NEIGHBORS( TYPE ) {                                                    \
  TYPE ***theBuf = (TYPE***)theIm->array;                                               \
  for ( z=0; z<theIm->nplanes; z++ ) {                                                  \
    /* fprintf( stderr, "Z=%lu\n", z); */                                               \
    for ( y=0; y<theIm->nrows; y++ )                                                    \
    for ( x=0; x<theIm->ncols; x++ ) {                                                  \
      label = theBuf[z][y][x];                                                          \
      /* if ( label == 0 || label == 1 ) continue;                                      \
       */                                                                               \
      nneighbors = 0;                                                                   \
      if ( x > 0 && theBuf[z][y][x-1] != label ) neighbors[ nneighbors++ ] = theBuf[z][y][x-1];                \
      if ( x < theIm->ncols-1 && theBuf[z][y][x+1] != label ) neighbors[ nneighbors++ ] = theBuf[z][y][x+1];   \
      if ( y > 0 && theBuf[z][y-1][x] != label ) neighbors[ nneighbors++ ] = theBuf[z][y-1][x];                \
      if ( y < theIm->nrows-1 && theBuf[z][y+1][x] != label ) neighbors[ nneighbors++ ] = theBuf[z][y+1][x];   \
      if ( z > 0 && theBuf[z-1][y][x] != label ) neighbors[ nneighbors++ ] = theBuf[z-1][y][x];                \
      if ( z < theIm->nplanes-1 && theBuf[z+1][y][x] != label ) neighbors[ nneighbors++ ] = theBuf[z+1][y][x]; \
      if ( nneighbors > 0 ) {                                                           \
        if ( nneighbors > 2 )                                                           \
          qsort( neighbors, nneighbors, sizeof(int), &_compareint );                    \
        for ( n=0; n<nneighbors; n++ ) {                                                \
            /* we had the point (x,y,z) as contact point to its neighbors,              \
             * so we can skip background labels                                         \
             */                                                                         \
            if ( neighbors[n] == 0 || neighbors[n] == 1 ) continue;                     \
            if ( n>=1 && neighbors[n] == neighbors[n-1] ) continue;                     \
            /* add 1 to surface contact of cell #neighbors[n] with cell #label          \
             */                                                                         \
            if( 0 ) fprintf( stderr, "[%d,%d,%d]=%d has neighbor %d\n", (int)x, (int)y, (int)z, label, neighbors[n] ); \
            if ( _addSpatialNeighborToCell( &(theList->data[ neighbors[n] ]), label, 1.0 ) != 1 ) { \
              if ( _verbose_ )                                                          \
                fprintf( stderr, "%s: unable to add neighbor %d to cell %d\n", proc, label, neighbors[n] ); \
              return( -1 );                                                             \
            }                                                                           \
        }                                                                               \
      }                                                                                 \
    }                                                                                   \
  }                                                                                     \
}

  switch ( theIm->type ) {
  default :
    fprintf( stderr, "%s: type not handled in switch\n", proc );
    return( -1 );
  case UCHAR :
    _OUTSIDE6NEIGHBORS( u8 );
    break;
  case USHORT :
    _OUTSIDE6NEIGHBORS( u16 );
    break;
  }

  return( 1 );
}




static int _WindreichSurfaceEstimation( bal_image *theIm, typeCellImage *theList )
{
  char *proc = "_WindreichSurfaceEstimation";
  size_t x, y, z;
  int label;
  int i, j, k, n;
  int neighborhood[27];
  int neighbors[6], nneighbors;
  float ds;

  /* Voxel-based surface area estimation: from theory to practice,
   * G. Windreich, N. Kiryati, G. Lohmann,
   * Pattern Recognition, 36(11):2531-2541, (2003)
   */

#define _WINDREICHCONFIGURATION( TYPE ) {                           \
  TYPE ***theBuf = (TYPE***)theIm->array;                           \
  for ( z=0; z<theIm->nplanes; z++ )                                \
  for ( y=0; y<theIm->nrows; y++ )                                  \
  for ( x=0; x<theIm->ncols; x++ ) {                                \
    label = theBuf[z][y][x];                                        \
    if ( label == 0 || label == 1 ) continue;                       \
    /* 26-neighborhood                                              \
     */                                                             \
    for ( n=0, k=-1; k<=1; k++ )                                    \
    for ( j=-1; j<=1; j++ )                                         \
    for ( i=-1; i<=1; i++, n++ ) {                                  \
        if ( 0 <= (int)x+i && (int)x+i < (int)theIm->ncols          \
             && 0 <= (int)y+j && (int)y+j < (int)theIm->nrows       \
             && 0 <= (int)z+k && (int)z+k < (int)theIm->nplanes )   \
            neighborhood[n] = theBuf[(int)z+k][(int)y+j][(int)x+i]; \
        else                                                        \
            neighborhood[n] = 0;                                    \
    }                                                               \
    /* 6-neighbors                                                  \
     */                                                             \
    nneighbors = 0;                                                 \
    if ( neighborhood[ 4] != label ) neighbors[ nneighbors++ ] = neighborhood[ 4]; \
    if ( neighborhood[10] != label ) neighbors[ nneighbors++ ] = neighborhood[10]; \
    if ( neighborhood[12] != label ) neighbors[ nneighbors++ ] = neighborhood[12]; \
    if ( neighborhood[14] != label ) neighbors[ nneighbors++ ] = neighborhood[14]; \
    if ( neighborhood[16] != label ) neighbors[ nneighbors++ ] = neighborhood[16]; \
    if ( neighborhood[22] != label ) neighbors[ nneighbors++ ] = neighborhood[22]; \
    if ( nneighbors == 0 ) continue;                                \
    /* surface voxel                                                \
     */                                                             \
    ds = _surface_windreich( neighborhood );                        \
    if ( ds < 0.0 ) {                                               \
      if ( _verbose_ )                                              \
        fprintf( stderr, "%s: weird surface contribution\n", proc ); \
      return( -1 );                                                 \
    }                                                               \
    if ( nneighbors == 1 ) {                                        \
      if ( _addSpatialNeighborToCell( &(theList->data[ label ]), neighbors[0], ds ) != 1 ) { \
        if ( _verbose_ )                                            \
          fprintf( stderr, "%s: unable to add neighbor %d to cell %d\n", proc, neighbors[0], label ); \
        return( -1 );                                               \
      }                                                             \
      continue;                                                     \
    }                                                               \
    if ( nneighbors > 2 )                                           \
      qsort( neighbors, nneighbors, sizeof(int), &_compareint );    \
    for ( i=0; i<nneighbors; i++ ) {                                \
      if ( i>=1 && neighbors[i] == neighbors[i-1] ) continue;       \
      for (n=0, j=i; j<nneighbors; j++ )                            \
          if ( neighbors[j] == neighbors[i] ) n++;                  \
      if ( _addSpatialNeighborToCell( &(theList->data[ label ]), neighbors[i], ds * (float)n / (float)nneighbors ) != 1 ) { \
          if ( _verbose_ )                                          \
            fprintf( stderr, "%s: unable to add neighbor %d to cell %d\n", proc, neighbors[i], label ); \
          return( -1 );                                             \
      }                                                             \
    }                                                               \
  }                                                                 \
}

  switch ( theIm->type ) {
  default :
    fprintf( stderr, "%s: type not handled in switch\n", proc );
    return( -1 );
  case UCHAR :
    _WINDREICHCONFIGURATION( u8 );
    break;
  case USHORT :
    _WINDREICHCONFIGURATION( u16 );
    break;
  }

  return( 1 );
}





/* configuration
 * 0 1  -  4 5
 * 2 3  -  6 7
 */

float _lindbladConfigurationEstimation( int *binary )
{
  char *proc = "_lindbladConfigurationEstimation";
  int i, n;

  for ( n=0, i=0; i<8; i++ )
    if ( binary[i] > 0 ) n++;

  if ( n == 0 ) {
    return( 0.0 );
  }
  else if ( n == 1 ) {
    return( 0.636 );
  }
  else if ( n == 2 ) {
    if ( (binary[0] > 0 && binary[1] > 0)
         || (binary[0] > 0 && binary[2] > 0)
         || (binary[0] > 0 && binary[4] > 0)
         || (binary[1] > 0 && binary[3] > 0)
         || (binary[1] > 0 && binary[5] > 0)
         || (binary[2] > 0 && binary[3] > 0)
         || (binary[2] > 0 && binary[6] > 0)
         || (binary[3] > 0 && binary[7] > 0)
         || (binary[4] > 0 && binary[5] > 0)
         || (binary[4] > 0 && binary[6] > 0)
         || (binary[5] > 0 && binary[7] > 0)
         || (binary[6] > 0 && binary[7] > 0) )
        return( 0.669 );
    if ( (binary[0] > 0 && binary[3] > 0)
         || (binary[0] > 0 && binary[5] > 0)
         || (binary[0] > 0 && binary[6] > 0)
         || (binary[1] > 0 && binary[2] > 0)
         || (binary[1] > 0 && binary[4] > 0)
         || (binary[1] > 0 && binary[7] > 0)
         || (binary[2] > 0 && binary[4] > 0)
         || (binary[2] > 0 && binary[7] > 0)
         || (binary[3] > 0 && binary[5] > 0)
         || (binary[3] > 0 && binary[6] > 0)
         || (binary[4] > 0 && binary[7] > 0)
         || (binary[5] > 0 && binary[6] > 0) )
        return( 1.272 );
    if ( (binary[0] > 0 && binary[7] > 0)
         || (binary[1] > 0 && binary[6] > 0)
         || (binary[2] > 0 && binary[5] > 0)
         || (binary[3] > 0 && binary[4] > 0) )
        return( 1.272 );
    if ( _verbose_ )
      fprintf( stderr, "%s: weird, this case (n==2) is not handled\n", proc );
    return( 0.0 );
  }
  else {
      return( _surface_lindblad(binary) );
  }

  if ( _verbose_ )
    fprintf( stderr, "%s: weird, this case (n>2) should not occur\n", proc );
  return( 0.0 );

}



typedef struct _typeNeighbor {
    int label;
    int n;
} _typeNeighbor;


int _lindbladConfiguration( int *configuration,
                         _typeNeighbor *neighbor, int nlabel,
                         typeCellImage *theList )
{
  char *proc = "_lindbladConfiguration";
  int i, l, m;
  int binary[8];
  float ds;

  if ( nlabel <= 1 ) return( 1 );

  /* 2 neighbors case
   */
  if ( nlabel == 2 ) {
    m = ( neighbor[0].n <= neighbor[1].n ) ? neighbor[0].label : neighbor[1].label;
    for ( i=0; i<8; i++ )
      binary[i] = (configuration[i] ==  m) ? 1 : 0;
    ds = _lindbladConfigurationEstimation( binary );
    if ( _addSpatialNeighborToCell( &(theList->data[ neighbor[0].label ]),
                             neighbor[1].label, ds ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to add neighbor %d to cell %d\n",
                   proc, neighbor[1].label, neighbor[0].label );
        return( -1 );
    }
    if ( _addSpatialNeighborToCell( &(theList->data[ neighbor[1].label ]),
                             neighbor[0].label, ds ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to add neighbor %d to cell %d\n",
                   proc, neighbor[0].label, neighbor[1].label );
        return( -1 );
    }
    return( 1 );
  }

  /* >2 neighbors case
   */
  for ( l=0; l<nlabel; l++ ) {
    for ( i=0; i<8; i++ )
      binary[i] = (configuration[i] ==  neighbor[l].label) ? 1 : 0;
    ds = _lindbladConfigurationEstimation( binary );
    for ( i=0; i<nlabel; i++ ) {
      if ( neighbor[i].label == neighbor[l].label ) continue;
      if ( _addSpatialNeighborToCell( &(theList->data[ neighbor[l].label ]),
                               neighbor[i].label,
                               ds*neighbor[i].n/(8.0-neighbor[l].n) ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to add neighbor %d to cell %d\n",
                     proc, neighbor[i].label, neighbor[l].label );
          return( -1 );
      }
    }
  }

  return( 1 );
}



static int _LindbladSurfaceEstimation( bal_image *theIm, typeCellImage *theList )
{
  char *proc = "_LindbladSurfaceEstimation";
  size_t x, y, z;
  int i, nlabel;
  int configuration[8], label[8];
  _typeNeighbor neighbor[8];

  /* Surface area estimation of digitized 3D objects using weighted local configurations,
   * Joakim Lindblad,
   * Image and Vision Computing, 23(2):111-122, 2005
   */

#define _LINDBLADCONFIGURATION( TYPE ) {                   \
  TYPE ***theBuf = (TYPE***)theIm->array;                  \
  for ( z=0; z<theIm->nplanes-1; z++ )                     \
  for ( y=0; y<theIm->nrows-1; y++ )                       \
  for ( x=0; x<theIm->ncols-1; x++ ) {                     \
    configuration[0] = label[0] = theBuf[z  ][y  ][x  ];   \
    configuration[1] = label[1] = theBuf[z  ][y  ][x+1];   \
    configuration[2] = label[2] = theBuf[z  ][y+1][x  ];   \
    configuration[3] = label[3] = theBuf[z  ][y+1][x+1];   \
    configuration[4] = label[4] = theBuf[z+1][y  ][x  ];   \
    configuration[5] = label[5] = theBuf[z+1][y  ][x+1];   \
    configuration[6] = label[6] = theBuf[z+1][y+1][x  ];   \
    configuration[7] = label[7] = theBuf[z+1][y+1][x+1];   \
    qsort( label, 8, sizeof(int), &_compareint );          \
    if ( label[0] == label[7] ) continue;                  \
    for ( i=0; i<8; i++ ) {                                \
      neighbor[i].label = -1;                              \
      neighbor[i].n = 0;                                   \
    }                                                      \
    for ( nlabel=-1, i=0; i<8; i++ ) {                     \
      if ( i==0 || (i>=1 && label[i] != label[i-1]) ) {    \
        nlabel++;                                          \
        neighbor[nlabel].label = label[i];                 \
        neighbor[nlabel].n = 1;                            \
      }                                                    \
      else {                                               \
        neighbor[nlabel].n ++;                             \
      }                                                    \
    }                                                      \
    nlabel ++;                                             \
    /* there are nlabel different labels                   \
     */                                                    \
    if ( _lindbladConfiguration( configuration, neighbor, nlabel, theList ) != 1 ) { \
      if ( _verbose_ )                                     \
        fprintf( stderr, "%s: error when processing point (%lu,%lu,%lu)\n", proc, x, y, z ); \
      return( -1 );                                        \
    }                                                      \
  }                                                        \
}

  switch ( theIm->type ) {
  default :
    fprintf( stderr, "%s: type not handled in switch\n", proc );
    return( -1 );
  case UCHAR :
    _LINDBLADCONFIGURATION( u8 );
    break;
  case USHORT :
    _LINDBLADCONFIGURATION( u16 );
    break;
  }

  return( 1 );
}





/************************************************************
 *
 *
 *
 ************************************************************/



static int _CellImageVolume( bal_image *theIm, typeCellImage *theList )
{
  char *proc = "_CellImageVolume";
  size_t x, y, z;
  int label;
  typeCell *cell = (typeCell*)NULL;
  int ncells;


  ncells = BAL_InitCellImageFromImage( theList, theIm );
  if ( ncells <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: list initialisation failed\n", proc );
    return( -1 );
  }

  /* this initialisation is mandatory for update purposes
   */
  for ( label = 0; label <= ncells; label ++ ) {
    cell = &(theList->data[ label ]);
    cell->npts = 0;
  }

#define _CELLVOLUME( TYPE ) {                  \
  TYPE ***theBuf = (TYPE***)theIm->array;      \
  for ( z=0; z<theIm->nplanes; z++ ) {         \
    for ( y=0; y<theIm->nrows; y++ )           \
    for ( x=0; x<theIm->ncols; x++ ) {         \
      label = theBuf[z][y][x];                 \
      /* if we want to count background as contact surface for cells,\
       * we have to consider background points     \
       * so do not skip them                       \
       * if ( label == 0 || label == 1 ) continue; \
       */                                      \
      cell = &(theList->data[ label ]);        \
      cell->npts ++;                           \
    }                                          \
  }                                            \
}


  /* compute cell basic volume
   *
   */
  switch ( theIm->type ) {
  default :
    fprintf( stderr, "%s: type not handled in switch\n", proc );
    return( -1 );
  case UCHAR :
    _CELLVOLUME( u8 );
    break;
  case USHORT :
    _CELLVOLUME( u16 );
    break;
  }
  return( 1 );
}





int BAL_CellImageProperties( bal_image *theIm, typeCellImage *theList,
                             typePropertyList *propertyList,
                             enumSurfaceEstimation surfaceType )
{
  char *proc = "BAL_CellImageProperties";
  size_t x, y, z, i;
  int label;
  typeCell *cell = (typeCell*)NULL;
  int ncells;


  ncells = BAL_InitCellImageFromImage( theList, theIm );
  if ( ncells <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: list initialisation failed\n", proc );
    return( -1 );
  }

  /* this initialisation is mandatory for update purposes
   */
  for ( label = 0; label <= ncells; label ++ ) {
    cell = &(theList->data[ label ]);
    cell->npts = 0;
    for ( i=0; i<3; i++ )
        cell->voxel_barycenter[i] = 0;
    for ( i=0; i<9; i++ )
        cell->voxel_covariance_matrix[i] = 0;
  }

#define _FIRSTCELLPROPERTIES( TYPE ) {         \
  TYPE ***theBuf = (TYPE***)theIm->array;      \
  for ( z=0; z<theIm->nplanes; z++ ) {         \
    for ( y=0; y<theIm->nrows; y++ )           \
    for ( x=0; x<theIm->ncols; x++ ) {         \
      label = theBuf[z][y][x];                 \
      /* if we want to count background as contact surface for cells,\
       * we have to consider background points     \
       * so do not skip them                       \
       * if ( label == 0 || label == 1 ) continue; \
       */                                      \
      cell = &(theList->data[ label ]);        \
      if ( cell->npts == 0 ) {                 \
        cell->ptmin[0] = cell->ptmax[0] = x;   \
        cell->ptmin[1] = cell->ptmax[1] = y;   \
        cell->ptmin[2] = cell->ptmax[2] = z;   \
      }                                        \
      else {                                   \
        if ( cell->ptmin[0] > (int)x ) cell->ptmin[0] = x; \
        if ( cell->ptmin[1] > (int)y ) cell->ptmin[1] = y; \
        if ( cell->ptmin[2] > (int)z ) cell->ptmin[2] = z; \
        if ( cell->ptmax[0] < (int)x ) cell->ptmax[0] = x; \
        if ( cell->ptmax[1] < (int)y ) cell->ptmax[1] = y; \
        if ( cell->ptmax[2] < (int)z ) cell->ptmax[2] = z; \
      }                                        \
      cell->npts ++;                           \
                                               \
      cell->voxel_barycenter[0] += x;          \
      cell->voxel_barycenter[1] += y;          \
      cell->voxel_barycenter[2] += z;          \
                                               \
      cell->voxel_covariance_matrix[0] += x*x; \
      cell->voxel_covariance_matrix[1] += x*y; \
      cell->voxel_covariance_matrix[2] += x*z; \
      cell->voxel_covariance_matrix[3] += y*x; \
      cell->voxel_covariance_matrix[4] += y*y; \
      cell->voxel_covariance_matrix[5] += y*z; \
      cell->voxel_covariance_matrix[6] += z*x; \
      cell->voxel_covariance_matrix[7] += z*y; \
      cell->voxel_covariance_matrix[8] += z*z; \
    }                                          \
  }                                            \
}


  /* compute cell basic properties
   *
   */
  switch ( theIm->type ) {
  default :
    fprintf( stderr, "%s: type not handled in switch\n", proc );
    return( -1 );
  case UCHAR :
    _FIRSTCELLPROPERTIES( u8 );
    break;
  case USHORT :
    _FIRSTCELLPROPERTIES( u16 );
    break;
  }


  /* compute inertia parameters
   * recompute barycenter
   */
  for ( label = 0; label <= ncells; label ++ ) {
    cell = &(theList->data[ label ]);
    if ( cell->npts == 0 ) continue;
    for ( i=0; i<3; i++ )
        cell->voxel_barycenter[i] /= (double)cell->npts;
  }


  if ( _isPropertyInList( propertyList, _PRINCIPAL_VALUE_ ) == 1
       || _isPropertyInList( propertyList, _PRINCIPAL_VECTOR_ ) == 1 ) {

    for ( label = 0; label <= ncells; label ++ ) {
      cell = &(theList->data[ label ]);
      if ( cell->npts == 0 ) continue;

      for ( i=0; i<9; i++ )
          cell->voxel_covariance_matrix[i] /= (double)cell->npts;

      /* center covariance matrix
       */
      cell->voxel_covariance_matrix[0] -= cell->voxel_barycenter[0]*cell->voxel_barycenter[0];
      cell->voxel_covariance_matrix[1] -= cell->voxel_barycenter[0]*cell->voxel_barycenter[1];
      cell->voxel_covariance_matrix[2] -= cell->voxel_barycenter[0]*cell->voxel_barycenter[2];
      cell->voxel_covariance_matrix[3] -= cell->voxel_barycenter[1]*cell->voxel_barycenter[0];
      cell->voxel_covariance_matrix[4] -= cell->voxel_barycenter[1]*cell->voxel_barycenter[1];
      cell->voxel_covariance_matrix[5] -= cell->voxel_barycenter[1]*cell->voxel_barycenter[2];
      cell->voxel_covariance_matrix[6] -= cell->voxel_barycenter[2]*cell->voxel_barycenter[0];
      cell->voxel_covariance_matrix[7] -= cell->voxel_barycenter[2]*cell->voxel_barycenter[1];
      cell->voxel_covariance_matrix[8] -= cell->voxel_barycenter[2]*cell->voxel_barycenter[2];

      if ( _ComputeEigensOfSymetricSquareMatrix( cell->voxel_covariance_matrix,
                                                 cell->voxel_eigenvalues,
                                                 cell->voxel_eigenvectors, 3 ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: error when computing covariance eigenvectors of cell %d\n", proc, label );
      }
      if ( _SortEigensInAbsDecreasingOrder( cell->voxel_eigenvalues,
                                            cell->voxel_eigenvectors, 3 ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: error when sorting covariance eigenvalues of cell %d\n", proc, label );
      }
    }
  }


  /* contact neighbors
   */
  if ( _isPropertyInList( propertyList, _SURFACE_ ) == 1
       || _isPropertyInList( propertyList, _COMPACTNESS_ ) == 1
       || _isPropertyInList( propertyList, _CONTACT_SURFACE_ ) == 1 ) {
    switch( surfaceType ) {
    default :
        if ( _verbose_ )
          fprintf( stderr, "%s: contact surface estimation type not implemented yet\n", proc );
        return( -1 );
    case _OUTER_6NEIGHBORS_ :
        if ( _Outer6NeighborsSurfaceEstimation( theIm, theList ) != 1 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: contact surface estimation failed\n", proc );
            return( -1 );
        }
        break;
    case _LINDBLAD_ :
        if ( _LindbladSurfaceEstimation( theIm, theList ) != 1 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: contact surface estimation failed\n", proc );
            return( -1 );
        }
        break;
    case _WINDREICH_ :
        if ( _WindreichSurfaceEstimation( theIm, theList ) != 1 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: contact surface estimation failed\n", proc );
            return( -1 );
        }
        break;
    }
    for ( label = 0; label <= ncells; label ++ ) {
      cell = &(theList->data[ label ]);
      cell->surface = 0.0;
      for ( i=0; i<(size_t)cell->neighbors.n_data; i++ ) {
        cell->surface += cell->neighbors.data[i].surface;
      }

      if ( cell->surface > 0.0 )
          cell->compactness = cbrt( (double)cell->npts ) / sqrt( cell->surface );
    }
  }

  return( 1 );
}





/************************************************************
 *
 *
 *
 ************************************************************/





int BAL_CellImageForwardIntersection( bal_image *theMother, bal_image *theDaughter, typeCellImage *theList )
{
  char *proc = "BAL_CellImageForwardIntersection";
  size_t x, y, z;
  int i;

  if ( theMother->type != theDaughter->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: images have different types\n", proc );
    return( -1 );
  }

  if ( theMother->nplanes != theDaughter->nplanes
       || theMother->nrows != theDaughter->nrows
       || theMother->ncols != theDaughter->ncols
       || theMother->vdim != theDaughter->vdim ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: images have different dimensions\n", proc );
    return( -1 );
  }

  /* initialize cell list, but also their volume
   * (tests could be done on volumes)
   */
  if ( _CellImageVolume( theMother, theList ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: list initialisation failed\n", proc );
    return( -1 );
  }

  /* re-init forward lineages
   */
  for ( i=0; i<theList->n_data; i++ ) {
      _freeTemporalNeighborList( &(theList->data[i].forward) );
  }


  /* keep background labels in intersection
   */
#define _IMAGELINEAGE( TYPE ) {                 \
  TYPE ***motBuf = (TYPE***)theMother->array;   \
  TYPE ***dauBuf = (TYPE***)theDaughter->array; \
  for ( z=0; z<theMother->nplanes; z++ )        \
  for ( y=0; y<theMother->nrows; y++ )          \
  for ( x=0; x<theMother->ncols; x++ ) {        \
    if ( 0 && (dauBuf[z][y][x] == 0 || dauBuf[z][y][x] == 1) ) continue; \
    if ( 0 && (motBuf[z][y][x] == 0 || motBuf[z][y][x] == 1) ) continue; \
    if ( _addForwardNeighborToCell( &(theList->data[motBuf[z][y][x]]),   \
                             _cellGlobalId(dauBuf[z][y][x], theList->next_acquisition_time), 1 ) != 1 ) { \
      if ( _verbose_ )                                            \
          fprintf( stderr, "%s: unable to add daughter %d to cell %d\n", \
                   proc, dauBuf[z][y][x], motBuf[z][y][x] );      \
        return( -1 );                           \
    }                                           \
  }                                             \
}

  switch ( theMother->type ) {
  default :
    fprintf( stderr, "%s: type not handled in switch\n", proc );
    return( -1 );
  case UCHAR :
    _IMAGELINEAGE( u16 );
    break;
  case USHORT :
    _IMAGELINEAGE( u16 );
    break;
  }

  /* sort forward neighbors
   */
  for ( i=0; i<theList->n_data; i++ ) {
    if ( theList->data[i].forward.n_data > 0 )
      qsort( theList->data[i].forward.data, theList->data[i].forward.n_data, sizeof(typeTemporalNeighbor), &_compareTemporalNeighbor );
  }

  return( 1 );
}





int BAL_CellImageBackwardIntersection( typeCellImage *theMother, typeCellImage *theDaughter )
{
  char *proc = "BAL_CellImageBackwardIntersection";
  int c, d;
  typeCell *cell;
  typeTemporalNeighbor *forward, *firstBackward;

  /* re-init backward lineages
   */
  for ( c=0; c<theDaughter->n_data; c++ ) {
      _freeTemporalNeighborList( &(theDaughter->data[c].backward) );
  }

  for ( c=0; c<theMother->n_data; c++ ) {
    if ( theMother->data[c].npts <= 0 ) continue;
    for ( d=0; d<theMother->data[c].forward.n_data; d++ ) {
        forward = &(theMother->data[c].forward.data[d]);
        cell = &(theDaughter->data[ forward->label ]);
        if ( _addBackwardNeighborToCell( cell, _cellGlobalId(c, theMother->acquisition_time), forward->npts ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to add mother %d to cell %d\n", proc, c, forward->label );
          return( -1 );
        }
    }
  }

  /* sort backward neighbors
   */
  for ( c=0; c<theMother->n_data; c++ ) {
      _freeTemporalNeighborList( &(theMother->data[c].lineage) );
  }

  for ( c=0; c<theDaughter->n_data; c++ ) {
    if ( theDaughter->data[c].backward.n_data > 0 ) {
      qsort( theDaughter->data[c].backward.data, theDaughter->data[c].backward.n_data, sizeof(typeTemporalNeighbor), &_compareTemporalNeighbor );
      firstBackward = &(theDaughter->data[c].backward.data[0]);
      cell = &(theMother->data[ firstBackward->label ]);
      if ( _addTemporalNeighborToCellLineage( cell, _cellGlobalId(c, theDaughter->acquisition_time), firstBackward->npts ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to add lineage %d to cell %d\n", proc, c, firstBackward->label );
        return( -1 );
      }
    }
  }

  return( 1 );
}





int BAL_CellSequenceBackwardIntersection( typeCellSequence *l )
{
  char *proc = "BAL_CellSequenceBackwardIntersection";
  int i;

  for ( i=0; i<l->n_data-1; i++ ) {
    if ( l->data[i].n_data <= 0 ) continue;

    /* retro-propagate cell intersections
     */
    if ( BAL_CellImageBackwardIntersection( &(l->data[i]), &(l->data[i+1]) ) != 1 ) {
      if ( _verbose_ )
          fprintf( stderr, "%s: error when extracting backward filiation\n", proc );
        return( -1 );
    }
  }

  return( 1 );
}





/************************************************************
 *
 *
 *
 ************************************************************/





typedef struct _diagCell {
  int idImage;
  int idCell;

  int cell_volume;

  int idMother;
  int mother_volume;
  int mother_intersection;

  float intersectionMotherRatio;
  float intersectionCellRatio;

  int ndaughters;
} _diagCell;



typedef struct _diagCellList {
  _diagCell *data;
  int n_data;
} _diagCellList;





static void _initDiagCell( _diagCell *c )
{
  c->idImage = 0;
  c->idCell = 0;

  c->cell_volume = 0;

  c->idMother = 0;
  c->mother_volume = 0;
  c->mother_intersection = 0;

  c->intersectionMotherRatio = 0.0;
  c->intersectionCellRatio = 0.0;

  c->ndaughters = 0;
}





static void _initDiagCellList( _diagCellList *l )
{
  l->data = (_diagCell*)NULL;
  l->n_data = 0;
}





static void _freeDiagCellList( _diagCellList *l )
{
  if ( l->data != (_diagCell*)NULL ) vtfree( l->data );
  _initDiagCellList( l );
}





static int _allocDiagCellList( _diagCellList *l, int n )
{
  char *proc = "_allocDiagCellList";
  int i;

  l->data = (_diagCell*)vtmalloc( n * sizeof(_diagCell), "l->data", proc );

  if ( l->data == (_diagCell*)NULL ) {
      if ( _verbose_ )
          fprintf( stderr, "%s: unable to allocate data\n", proc );
    return( -1 );
  }
  l->n_data = n;

  for ( i=0; i<l->n_data; i++ )
      _initDiagCell( &(l->data[i]) );

  return( 1 );
}





static void _fillDiagCellListFromCellImage( _diagCellList *l,
                                            typeCellImage *previousCellImage,
                                            typeCellImage *cellImage,
                                            int firstIndex )
{
  int i;
  for ( i=0; i<cellImage->n_data; i++ ) {
    l->data[firstIndex+i].idImage = cellImage->acquisition_time;
    l->data[firstIndex+i].idCell = i;
    l->data[firstIndex+i].cell_volume = cellImage->data[i].npts;

    if ( cellImage->data[i].backward.n_data > 0 ) {
        l->data[firstIndex+i].idMother = cellImage->data[i].backward.data[0].label;
        l->data[firstIndex+i].mother_intersection = cellImage->data[i].backward.data[0].npts;
        if ( previousCellImage != (typeCellImage*)NULL ) {
            l->data[firstIndex+i].mother_volume = previousCellImage->data[ l->data[firstIndex+i].idMother ].npts;
        }
        l->data[firstIndex+i].intersectionCellRatio = (float)l->data[firstIndex+i].mother_intersection / (float)l->data[firstIndex+i].cell_volume;
        l->data[firstIndex+i].intersectionMotherRatio = (float)l->data[firstIndex+i].mother_intersection / (float)l->data[firstIndex+i].mother_volume;
    }
    else {
      l->data[firstIndex+i].idMother = -1;
    }

    l->data[firstIndex+i].ndaughters = cellImage->data[i].lineage.n_data;
  }
}





static int _allocDiagCellListFromCellImage( _diagCellList *l, typeCellImage *cellImage )
{
  char *proc = "_allocDiagCellListFromCellImage";

  if ( _allocDiagCellList( l, cellImage->n_data ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( -1 );
  }
  _fillDiagCellListFromCellImage( l, (typeCellImage*)NULL, cellImage, 0 );

  return( 1 );
}




static int _allocDiagCellListFromCellSequence( _diagCellList *l, typeCellSequence *cellSequence )
{
  char *proc = "_allocDiagCellListFromCellSequence";
  int n, i;

  for ( n=0, i=0; i<cellSequence->n_data; i++ ) {
      n += cellSequence->data[i].n_data;
  }

  if ( _allocDiagCellList( l, n ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( -1 );
  }

  for ( n=0, i=0; i<cellSequence->n_data; i++ ) {
    if ( n == 0 ) {
        _fillDiagCellListFromCellImage( l, (typeCellImage*)NULL, &(cellSequence->data[i]), n );
    }
    else {
        _fillDiagCellListFromCellImage( l, &(cellSequence->data[i-1]), &(cellSequence->data[i]), n );
    }
    n += cellSequence->data[i].n_data;
  }

  return( 1 );
}




static int _diagCompareVolumes ( const void * a, const void * b )
{
    _diagCell *da = (_diagCell*)a;
    _diagCell *db = (_diagCell*)b;

    /* sort volumes in increasing order
     */
    if ( da->cell_volume < db->cell_volume ) return( -1 );
    if ( da->cell_volume > db->cell_volume ) return( 1 );
    return( 0 );
}





static int _diagCompareIMR ( const void * a, const void * b )
{
    _diagCell *da = (_diagCell*)a;
    _diagCell *db = (_diagCell*)b;

    /* sort volumes in increasing order
     */
    if ( da->intersectionMotherRatio < db->intersectionMotherRatio ) return( -1 );
    if ( da->intersectionMotherRatio > db->intersectionMotherRatio ) return( 1 );
    return( 0 );
}



static int _diagCompareICR ( const void * a, const void * b )
{
    _diagCell *da = (_diagCell*)a;
    _diagCell *db = (_diagCell*)b;

    /* sort volumes in increasing order
     */
    if ( da->intersectionCellRatio < db->intersectionCellRatio ) return( -1 );
    if ( da->intersectionCellRatio > db->intersectionCellRatio ) return( 1 );
    return( 0 );
}





static int _cellImageVolumeDiagnosis( FILE *f, typeCellImage *cellImage, int nmin, int nmax )
{
  char *proc = "_cellImageVolumeDiagnosis";
  _diagCellList diagList;
  int i, j;

  _initDiagCellList( &diagList );
  if ( cellImage->n_data <= 0 ) return( 1 );

  if ( _allocDiagCellListFromCellImage( &diagList, cellImage ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }

  qsort( diagList.data, diagList.n_data, sizeof(_diagCell), &_diagCompareVolumes );

  fprintf( f, "  = image #%3d\n", cellImage->acquisition_time );
  for ( i=0; diagList.data[i].cell_volume <= 0; i++ ) {
      ;
  }
  for ( j=0; j<nmin; j++, i++ ) {
    fprintf( f, "    cell #%3d: volume %8d\n", diagList.data[i].idCell, diagList.data[i].cell_volume );
  }
  fprintf( f, "    ...\n" );
  for ( i=diagList.n_data-nmax; i<diagList.n_data; i++ ) {
    fprintf( f, "    cell #%3d: volume %8d\n", diagList.data[i].idCell, diagList.data[i].cell_volume );
  }

  _freeDiagCellList( &diagList );
  return( 1 );
}





static int _cellSequenceVolumeDiagnosis( FILE *f, typeCellSequence *cellProperties, int nmin, int nmax )
{
  char *proc = "_cellSequenceVolumeDiagnosis";
  int i, j;
  _diagCellList diagList;

  if ( 0 ) {
    fprintf( f, "=== image-based diagnosis\n" );
    for ( i=0; i<cellProperties->n_data; i++ ) {
      if ( cellProperties->data[i].n_data <= 0 ) continue;
      if ( _cellImageVolumeDiagnosis( f, &(cellProperties->data[i]), 2, 2 ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: volume checking failed for entry #%d (acquisition=%d)\n", proc, i, cellProperties->data[i].acquisition_time );
        return( 1 );
      }
    }
  }

  _initDiagCellList( &diagList );
  if ( _allocDiagCellListFromCellSequence( &diagList, cellProperties ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }

  qsort( diagList.data, diagList.n_data, sizeof(_diagCell), &_diagCompareVolumes );

  fprintf( f, "=== sequence-based diagnosis\n" );
  for ( i=0; diagList.data[i].cell_volume <= 0; i++ ) {
      ;
  }
  for ( j=0; j<nmin && i<diagList.n_data; i++ ) {
    if ( diagList.data[i].idCell == 0 || diagList.data[i].idCell == 1 )
        continue;
    fprintf( f, "    image #%3d, cell #%3d: volume %8d\n", diagList.data[i].idImage, diagList.data[i].idCell, diagList.data[i].cell_volume );
    j++;
  }

  fprintf( f, "    ...\n" );

  for ( j=1, i=diagList.n_data; j<nmax && i>=0; i-- ) {
      if ( diagList.data[i].idCell == 0 || diagList.data[i].idCell == 1 )
          continue;
      j++;
  }
  for ( ; i<diagList.n_data; i++ ) {
      if ( diagList.data[i].idCell == 0 || diagList.data[i].idCell == 1 )
          continue;
      fprintf( f, "    image #%3d, cell #%3d: volume %8d\n", diagList.data[i].idImage, diagList.data[i].idCell, diagList.data[i].cell_volume );
  }

  _freeDiagCellList( &diagList );
  return( 1 );
}





static int _cellSequenceLineageDiagnosis( FILE *f, typeCellSequence *cellProperties )
{
  char *proc = "_cellSequenceLineageDiagnosis";
  _diagCellList diagList;
  int i, j, n ,firstImage;
  typeCell *cell;

  _initDiagCellList( &diagList );
  if ( _allocDiagCellListFromCellSequence( &diagList, cellProperties ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }

  for ( i=0; cellProperties->data[i].n_data==0; i++ )
      ;
  firstImage = i;


  fprintf( f, "=== appearing cells === \n" );
  for ( n=0, i=0; i<diagList.n_data; i++ ) {
    if ( diagList.data[i].idImage == firstImage ) continue;
    if ( diagList.data[i].idCell == 0 || diagList.data[i].idCell == 1 ) continue;
    if ( diagList.data[i].cell_volume == 0 ) continue;
    if ( diagList.data[i].idMother == -1 ) continue;
    if ( diagList.data[i].idMother == 0 || diagList.data[i].idMother == 1 )
        fprintf( f, "    image #%3d, cell #%3d\n", diagList.data[i].idImage, diagList.data[i].idCell );
  }

  fprintf( f, "=== weird daughters number (<1 || >2) === \n" );
  for ( i=0; i<diagList.n_data; i++ ) {
    if ( diagList.data[i].idImage == cellProperties->n_data-1 ) continue;
    if ( diagList.data[i].idCell == 0 || diagList.data[i].idCell == 1 ) continue;
    if ( diagList.data[i].cell_volume == 0 ) continue;
    if ( diagList.data[i].ndaughters != 1 && diagList.data[i].ndaughters != 2 ) {
        fprintf( f, "    image #%3d, cell #%3d: %2d daughters:", diagList.data[i].idImage, diagList.data[i].idCell,
                 diagList.data[i].ndaughters );
        cell = &(cellProperties->data[diagList.data[i].idImage].data[diagList.data[i].idCell]);
        for ( j=0; j<cell->lineage.n_data; j++ ) {
          fprintf( f, " %3d (vol=%d)", cell->lineage.data[j].label, cell->lineage.data[j].npts );
          if ( j<cell->lineage.n_data-1 ) fprintf( f, ",");
        }
        fprintf( f, "\n" );
    }
  }

  qsort( diagList.data, diagList.n_data, sizeof(_diagCell), &_diagCompareIMR );

  fprintf( f, "=== (cell intersection mother) / mother === \n" );
  for ( i=0, n=0; i<diagList.n_data && n<20; i++ ) {
    if ( diagList.data[i].idImage == firstImage ) continue;
    if ( diagList.data[i].idMother == -1 ) continue;
    if ( diagList.data[i].idCell == 0 || diagList.data[i].idCell == 1 ) continue;
    if ( diagList.data[i].idMother == 0 || diagList.data[i].idMother == 1 ) continue;
    fprintf( f, "    image #%3d, cell #%3d: mother #%3d, vol_intersection=%8d, vol_mother=%8d, ratio=%8.6f\n",
             diagList.data[i].idImage, diagList.data[i].idCell, diagList.data[i].idMother,
             diagList.data[i].mother_intersection, diagList.data[i].mother_volume,
             diagList.data[i].intersectionMotherRatio );
    n++;
  }
  fprintf( f, "    ...\n" );

  qsort( diagList.data, diagList.n_data, sizeof(_diagCell), &_diagCompareICR );

  fprintf( f, "=== (cell intersection mother) / cell === \n" );
  for ( i=0, n=0; i<diagList.n_data && n<20; i++ ) {
    if ( diagList.data[i].idImage == firstImage ) continue;
    if ( diagList.data[i].idMother == -1 ) continue;
    if ( diagList.data[i].idCell == 0 || diagList.data[i].idCell == 1 ) continue;
    if ( diagList.data[i].idMother == 0 || diagList.data[i].idMother == 1 ) continue;
    fprintf( f, "    image #%3d, cell #%3d: mother #%3d, vol_intersection=%8d, vol_cell=%8d, ratio=%8.6f\n",
             diagList.data[i].idImage, diagList.data[i].idCell, diagList.data[i].idMother,
             diagList.data[i].mother_intersection, diagList.data[i].cell_volume,
             diagList.data[i].intersectionCellRatio );
    n++;
  }
  fprintf( f, "    ...\n" );


  return( 1 );
}





int BAL_CellSequenceDiagnosis( FILE *f, typeCellSequence *cellProperties, typePropertyList *properties )
{
  char *proc = "BAL_CellSequenceDiagnosis";
  int i;

  for ( i=0; i<properties->n_data; i++ ) {
    switch( properties->data[i] ) {
    default :
      break;
    case _VOLUME_ :
      fprintf( f, "============= VOLUME CHECKING =============\n" );
      if ( _cellSequenceVolumeDiagnosis( f, cellProperties, 20, 10 ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: error when checking cell volumes\n", proc );
        return( -1 );
      }
      fprintf( f, "===========================================\n" );
      break;
    case _LINEAGE_ :
        fprintf( f, "================= LINEAGE =================\n" );
        if ( _cellSequenceLineageDiagnosis( f, cellProperties ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: error when checking lineage\n", proc );
          return( -1 );
        }
        fprintf( f, "===========================================\n" );
      break;
    }
  }

  return( 1 );
}





/************************************************************
 *
 * XML writing facilities
 *
 ************************************************************/





static void _XML_MetaHeader( FILE *f, int d )
{
    int n;

    for (n=0; n<d; n++) fprintf( f, "  " );
    fprintf(f, "<data>\n"); return;
}



static void _XML_MetaFooter( FILE *f, int d )
{
    int n;

    for (n=0; n<d; n++) fprintf( f, "  " );
    fprintf(f, "</data>\n"); return;
}



static void _XML_Header( FILE *f, enumProperty property, int spaces, int cr )
{
  char *proc = "_XML_Header";
  int n;

  for (n=0; n<spaces; n++) fprintf( f, "  " );

  switch( property ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: property '%s' not handled yet\n", proc, _Property(property) );
    return;
  case _LINEAGE_ :               fprintf(f, "<cell_lineage>"); break;
  case _FORWARD_NEIGHBORS_ :     fprintf(f, "<cell_forward_neighbors>"); break;
  case _BACKWARD_NEIGHBORS_ :    fprintf(f, "<cell_backward_neighbors>"); break;
  case _H_MIN_ :                 fprintf(f, "<cell_h_min>"); break;
  case _VOLUME_ :                fprintf(f, "<cell_volume>"); break;
  case _SURFACE_ :               fprintf(f, "<cell_surface>"); break;
  case _COMPACTNESS_ :           fprintf(f, "<cell_compactness>"); break;
  case _SIGMA_ :                 fprintf(f, "<cell_sigma>"); break;
  case _LABEL_IN_TIME_ :         fprintf(f, "<cell_labels_in_time>"); break;
  case _BARYCENTER_ :            fprintf(f, "<cell_barycenter>"); break;
  case _PRINCIPAL_VALUE_ :       fprintf(f, "<cell_principal_values>"); break;
  case _PRINCIPAL_VECTOR_ :      fprintf(f, "<cell_principal_vectors>"); break;
  case _FATE_ :                  fprintf(f, "<cell_fate>"); break;
  case _FATE2_ :                 fprintf(f, "<cell_fate_2>"); break;
  case _FATE3_ :                 fprintf(f, "<cell_fate_3>"); break;
  case _FATE4_ :                 fprintf(f, "<cell_fate_4>"); break;
  case _ALL_CELLS_ :             fprintf(f, "<all_cells>"); break;
  case _NAME_ :                  fprintf(f, "<cell_name>"); break;
  case _HISTORY_ :               fprintf(f, "<cell_history>"); break;
  case _CONTACT_SURFACE_ :       fprintf(f, "<cell_contact_surface>"); break;
  }

  if ( cr > 0 ) fprintf(f, "\n");
}



static void _XML_Footer( FILE *f, enumProperty property, int d )
{
  char *proc = "_XML_Footer";
  int n;

  for (n=0; n<d; n++) fprintf( f, "  " );

  switch( property ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: property '%s' not handled yet\n", proc, _Property(property) );
    return;
  case _LINEAGE_ :               fprintf(f, "</cell_lineage>\n"); return;
  case _FORWARD_NEIGHBORS_ :     fprintf(f, "</cell_forward_neighbors>\n"); return;
  case _BACKWARD_NEIGHBORS_ :    fprintf(f, "</cell_backward_neighbors>\n"); return;
  case _H_MIN_ :                 fprintf(f, "</cell_h_min>\n"); return;
  case _VOLUME_ :                fprintf(f, "</cell_volume>\n"); return;
  case _SURFACE_ :               fprintf(f, "</cell_surface>\n"); return;
  case _COMPACTNESS_ :           fprintf(f, "</cell_compactness>\n"); return;
  case _SIGMA_ :                 fprintf(f, "</cell_sigma>\n"); return;
  case _LABEL_IN_TIME_ :         fprintf(f, "</cell_labels_in_time>\n"); return;
  case _BARYCENTER_ :            fprintf(f, "</cell_barycenter>\n"); return;
  case _PRINCIPAL_VALUE_ :       fprintf(f, "</cell_principal_values>\n"); return;
  case _PRINCIPAL_VECTOR_ :      fprintf(f, "</cell_principal_vectors>\n"); return;
  case _FATE_ :                  fprintf(f, "</cell_fate>\n"); return;
  case _FATE2_ :                 fprintf(f, "</cell_fate_2>\n"); return;
  case _FATE3_ :                 fprintf(f, "</cell_fate_3>\n"); return;
  case _FATE4_ :                 fprintf(f, "</cell_fate_4>\n"); return;
  case _ALL_CELLS_ :             fprintf(f, "</all_cells>\n"); return;
  case _NAME_ :                  fprintf(f, "</cell_name>\n"); return;
  case _HISTORY_ :               fprintf(f, "</cell_history>\n"); return;
  case _CONTACT_SURFACE_ :       fprintf(f, "</cell_contact_surface>\n"); return;
  }

  if ( _verbose_ )
    fprintf( stderr, "%s: this should not be reached, property was '%s'\n", proc, _Property(property) );
  return;
}



static void _XML_WriteCellImageProperty( FILE *f, typeCellImage *l, enumProperty property, int d )
{
  char *proc = "_XML_WriteCellImageProperty";
  int i, j, n;

  /* skip cell #0 and #1
   */
  for ( i=2; i<l->n_data; i++ ) {

    /* skip empty cells
     */
    switch( property ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: property '%s' not handled yet (1)\n", proc, _Property(property) );
      break;
    case _LINEAGE_ :
        if ( l->data[i].lineage.n_data <= 0 ) continue;
        break;
    case _FORWARD_NEIGHBORS_ :
        if ( l->data[i].forward.n_data <= 0 ) continue;
        break;
    case _BACKWARD_NEIGHBORS_ :
        if ( l->data[i].backward.n_data <= 0 ) continue;
        break;
    case _VOLUME_ :
        if ( l->data[i].npts <= 0 ) continue;
        break;
    case _SURFACE_ :
        if ( l->data[i].surface <= 0.0 ) continue;
        break;
    case _COMPACTNESS_ :
        if ( l->data[i].compactness <= 0.0 ) continue;
        break;
    case _BARYCENTER_ :
    case _PRINCIPAL_VALUE_ :
    case _PRINCIPAL_VECTOR_ :
        if ( l->data[i].npts <= 0 ) continue;
        break;
    case _FATE_ :
        if ( l->data[i].fate == (char*)NULL ) continue;
        break;
    case _FATE2_ :
        if ( l->data[i].fate2 == (char*)NULL ) continue;
        break;
    case _FATE3_ :
        if ( l->data[i].fate3 == (char*)NULL ) continue;
        break;
    case _FATE4_ :
        if ( l->data[i].fate4 == (char*)NULL ) continue;
        break;
    case _NAME_ :
        if ( l->data[i].name == (char*)NULL ) continue;
        break;
    case _HISTORY_ :
        if ( l->data[i].history.n_data <= 0 ) continue;
        break;
    case _CONTACT_SURFACE_ :
        if ( l->data[i].neighbors.n_data <= 0 ) continue;
        break;
    }

    switch( property ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: property '%s' not handled yet (2)\n", proc, _Property(property) );
      break;
    case _LINEAGE_ :
    case _HISTORY_ :
      break;
    case _FORWARD_NEIGHBORS_ :
    case _BACKWARD_NEIGHBORS_ :
    case _VOLUME_ :
    case _SURFACE_ :
    case _COMPACTNESS_ :
    case _BARYCENTER_ :
    case _PRINCIPAL_VALUE_ :
    case _PRINCIPAL_VECTOR_ :
    case _FATE_ :
    case _FATE2_ :
    case _FATE3_ :
    case _FATE4_ :
    case _NAME_ :
    case _CONTACT_SURFACE_ :
      if ( (property == _FORWARD_NEIGHBORS_ && l->data[i].forward.n_data <= 0)
           || (property == _BACKWARD_NEIGHBORS_ && l->data[i].backward.n_data <= 0) )
        break;
      for (n=0; n<d; n++) fprintf( f, "  " );
      fprintf( f, "<cell cell-id=\"%d\">", _cellGlobalId(i, l->acquisition_time) );
      break;
    }


    switch( property ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: property '%s' not handled yet (3)\n", proc, _Property(property) );
      break;
    case _LINEAGE_ :
      if ( l->data[i].lineage.n_data > 0 ) {
        for (n=0; n<d; n++) fprintf( f, "  " );
        fprintf( f, "<cell cell-id=\"%d\">", _cellGlobalId(i, l->acquisition_time) );
        fprintf( f, "[" );
        for ( j=0; j<l->data[i].lineage.n_data; j++ ) {
          fprintf( f, "%d", _cellGlobalId(l->data[i].lineage.data[j].label,
                                    l->next_acquisition_time) );
          if ( j<l->data[i].lineage.n_data-1 )
            fprintf( f, ", " );
        }
        fprintf( f, "]" );
        fprintf( f, "</cell>\n" );
      }
      break;
    case _FORWARD_NEIGHBORS_ :
      if ( l->data[i].forward.n_data > 0 ) {
        fprintf( f, "\n" );
        for ( j=0; j<l->data[i].forward.n_data; j++ ) {
          for (n=0; n<d+1; n++) fprintf( f, "  " );
          fprintf( f, "<cell cell-id=\"%d\">", _cellGlobalId(l->data[i].forward.data[j].label,
                                                       l->next_acquisition_time) );
          fprintf( f, "%d", l->data[i].forward.data[j].npts );
          fprintf( f, "</cell>\n" );
        }
        for (n=0; n<d; n++) fprintf( f, "  " );
      }
      break;
    case _BACKWARD_NEIGHBORS_ :
      if ( l->data[i].backward.n_data > 0 ) {
        fprintf( f, "\n" );
        for ( j=0; j<l->data[i].backward.n_data; j++ ) {
          for (n=0; n<d+1; n++) fprintf( f, "  " );
          fprintf( f, "<cell cell-id=\"%d\">", _cellGlobalId(l->data[i].backward.data[j].label,
                                                       l->previous_acquisition_time) );
          fprintf( f, "%d", l->data[i].backward.data[j].npts );
          fprintf( f, "</cell>\n" );
        }
        for (n=0; n<d; n++) fprintf( f, "  " );
      }
      break;
    case _H_MIN_ :
      if ( _verbose_ )
        fprintf( stderr, "%s: property '%s' not handled yet\n", proc, _Property(property) );
      break;
    case _VOLUME_ :
      fprintf( f, "%d", l->data[i].npts );
      break;
    case _SURFACE_ :
      fprintf( f, "%f", l->data[i].surface );
      break;
    case _COMPACTNESS_ :
      fprintf( f, "%f", l->data[i].compactness );
      break;
    case _SIGMA_ :
    case _LABEL_IN_TIME_ :
      if ( _verbose_ )
        fprintf( stderr, "%s: property '%s' not handled yet\n", proc, _Property(property) );
      break;
    case _BARYCENTER_ :
      fprintf( f, "[%f, %f, %f]", l->data[i].voxel_barycenter[0], l->data[i].voxel_barycenter[1], l->data[i].voxel_barycenter[2] );
      break;
    case _PRINCIPAL_VALUE_ :
      fprintf( f, "[%f, %f, %f]", l->data[i].voxel_eigenvalues[0], l->data[i].voxel_eigenvalues[1], l->data[i].voxel_eigenvalues[2] );
      break;
    case _PRINCIPAL_VECTOR_ :
        fprintf( f, "[[%f, %f, %f], [%f, %f, %f], [%f, %f, %f]]",
                 l->data[i].voxel_eigenvectors[0], l->data[i].voxel_eigenvectors[3], l->data[i].voxel_eigenvectors[6],
                 l->data[i].voxel_eigenvectors[1], l->data[i].voxel_eigenvectors[4], l->data[i].voxel_eigenvectors[7],
                 l->data[i].voxel_eigenvectors[2], l->data[i].voxel_eigenvectors[5], l->data[i].voxel_eigenvectors[8] );
        break;
    case _FATE_ :
        fprintf( f, "'%s'", l->data[i].fate );
        break;
    case _FATE2_ :
        fprintf( f, "'%s'", l->data[i].fate2 );
        break;
    case _FATE3_ :
        fprintf( f, "'%s'", l->data[i].fate3 );
        break;
    case _FATE4_ :
        fprintf( f, "'%s'", l->data[i].fate4 );
        break;
    case _ALL_CELLS_ :
        if ( _verbose_ )
          fprintf( stderr, "%s: property '%s' not handled yet\n", proc, _Property(property) );
        break;
    case _NAME_ :
        fprintf( f, "'%s'", l->data[i].name );
        break;
    case _HISTORY_ :
      if ( l->data[i].history.n_data > 0 ) {
        for (n=0; n<d; n++) fprintf( f, "  " );
        fprintf( f, "<cell cell-id=\"%d\">", _cellGlobalId(i, l->acquisition_time) );
        fprintf( f, "[" );
        for ( j=0; j<l->data[i].lineage.n_data; j++ ) {
          fprintf( f, "%d", l->data[i].history.data[j].globalId );
          if ( j<l->data[i].history.n_data-1 )
            fprintf( f, ", " );
        }
        fprintf( f, "]" );
        fprintf( f, "</cell>\n" );
      }
      break;
    case _CONTACT_SURFACE_ :
      if ( l->data[i].neighbors.n_data > 0 ) {
        fprintf( f, "\n" );
        for ( j=0; j<l->data[i].neighbors.n_data; j++ ) {
          for (n=0; n<d+1; n++) fprintf( f, "  " );
          fprintf( f, "<cell cell-id=\"%d\">", _cellGlobalId(l->data[i].neighbors.data[j].label,
                                                       l->acquisition_time) );
          fprintf( f, "%f", l->data[i].neighbors.data[j].surface );
          fprintf( f, "</cell>\n" );
        }
        for (n=0; n<d; n++) fprintf( f, "  " );
      }
      break;
    }

    switch( property ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: property '%s' not handled yet (4)\n", proc, _Property(property) );
      break;
    case _LINEAGE_ :
    case _HISTORY_ :
      break;
    case _FORWARD_NEIGHBORS_ :
    case _BACKWARD_NEIGHBORS_ :
    case _VOLUME_ :
    case _SURFACE_ :
    case _COMPACTNESS_ :
    case _BARYCENTER_ :
    case _PRINCIPAL_VALUE_ :
    case _PRINCIPAL_VECTOR_ :
    case _FATE_ :
    case _FATE2_ :
    case _FATE3_ :
    case _FATE4_ :
    case _NAME_ :
    case _CONTACT_SURFACE_ :
      if ( (property == _FORWARD_NEIGHBORS_ && l->data[i].forward.n_data <= 0)
           || (property == _BACKWARD_NEIGHBORS_ && l->data[i].backward.n_data <= 0) )
        break;
      fprintf( f, "</cell>\n" );
      break;
    }


  }
}


/* print header, content, and footer
 * for the property
 */
extern void _XML_FprintfCellImageProperty( FILE *f, typeCellImage *l, enumProperty property, int d )
{
  int n, j, k=0;

  switch( property ) {
  default :
      break;
  case _ALL_CELLS_ :
      _XML_Header( f, property, d, 0 );
      fprintf( f, "[" );
      for ( j=2; j<l->n_data; j++ ) {
        /* skip empty cells
         */
        if ( l->data[j].npts <= 0 ) continue;
        /* extra prints
         */
        if ( k>0 ) fprintf( f, ", " );
        if ( k>0 && k%10 == 0 ) {
          fprintf( f, "\n" );
          for (n=0; n<d+1; n++) fprintf( f, "  " );
          fprintf( f, " " );
        }
        /* print cell
         */
        fprintf( f, "%d", _cellGlobalId(j, l->acquisition_time) );
        k++;
      }
      fprintf( f, "]" );
      _XML_Footer( f, property, 0 );
      return;
  }

  _XML_Header( f, property, d, 1 );
  _XML_WriteCellImageProperty( f, l, property,  d+1 );
  _XML_Footer( f, property, d );
}



void XML_FprintfCellImageProperty( FILE *f, typeCellImage *l, typePropertyList *output )
{
  int n;

  /* one single feature
   */
  if ( output != (typePropertyList*)NULL && output->n_data == 1 ) {
    _XML_FprintfCellImageProperty( f, l, output->data[0], 0 );
  }
  else {
    _XML_MetaHeader( f, 0 );
    /* all features
     */
    if ( output == (typePropertyList*)NULL || output->n_data == 0 ) {
      _XML_FprintfCellImageProperty( f, l, _VOLUME_, 1 );
      _XML_FprintfCellImageProperty( f, l, _SURFACE_, 1 );
      _XML_FprintfCellImageProperty( f, l, _COMPACTNESS_, 1 );
      _XML_FprintfCellImageProperty( f, l, _BARYCENTER_, 1 );
      _XML_FprintfCellImageProperty( f, l, _PRINCIPAL_VALUE_, 1 );
      _XML_FprintfCellImageProperty( f, l, _PRINCIPAL_VECTOR_, 1 );
      _XML_FprintfCellImageProperty( f, l, _CONTACT_SURFACE_, 1 );
      _XML_FprintfCellImageProperty( f, l, _LINEAGE_, 1 );
      _XML_FprintfCellImageProperty( f, l, _ALL_CELLS_, 1 );
    }
    /* selection of features
     */
    else {
      for ( n=0; n<output->n_data; n++ )
        _XML_FprintfCellImageProperty( f, l, output->data[n], 1 );
    }
    _XML_MetaFooter( f, 0 );
  }
}





void _XML_FprintfCellSequenceProperty( FILE *f, typeCellSequence *l, enumProperty property, int d )
{
  int i;
  int n, j, k=0;

  switch( property ) {
  default :
      break;
  case _ALL_CELLS_ :
      _XML_Header( f, property, d, 0 );
      fprintf( f, "[" );
      for (i=0; i<l->n_data; i++) {
          for ( j=2; j<l->data[i].n_data; j++ ) {
            /* skip empty cells
             */
            if ( l->data[i].data[j].npts <= 0 ) continue;
            /* extra prints
             */
            if ( k>0 ) fprintf( f, ", " );
            if ( k>0 && k%10 == 0 ) {
              fprintf( f, "\n" );
              for (n=0; n<d+1; n++) fprintf( f, "  " );
              fprintf( f, " " );
            }
            /* print cell
             */
            fprintf( f, "%d", _cellGlobalId(j, l->data[i].acquisition_time) );
            k++;
          }
      }
      fprintf( f, "]" );
      _XML_Footer( f, property, 0 );
      return;
  }

  _XML_Header( f, property, d, 1 );
  for (i=0; i<l->n_data; i++)
    _XML_WriteCellImageProperty( f, &(l->data[i]), property,  d+1 );
  _XML_Footer( f, property, d );
}





void XML_FprintfCellSequenceProperty( FILE *f, typeCellSequence *l, typePropertyList *output )
{
  int n;
  typePropertyList localProp;
  typePropertyList *ptrProp = (typePropertyList*)NULL;

  /* one single feature
   */
  if ( output != (typePropertyList*)NULL && output->n_data == 1 ) {
    _XML_FprintfCellSequenceProperty( f, l, output->data[0], 0 );
    return;
  }

  /* multiple features
   * encapsulation of all properties
   */
  BAL_InitPropertyList( &localProp );
  if ( output == (typePropertyList*)NULL || output->n_data == 0 ) {
    BAL_DefaultPropertyList( &localProp );
    ptrProp = &localProp;
  }
  else {
    ptrProp = output;
  }

  _XML_MetaHeader( f, 0 );

  for ( n=0; n<ptrProp->n_data; n++ )
    _XML_FprintfCellSequenceProperty( f, l, ptrProp->data[n], 1 );

  _XML_MetaFooter( f, 0 );

}





/************************************************************
 *
 * XML reading facilities
 *
 ************************************************************/







static int _hasString( char *p, char *s )
{
   char *ptr = p;

   while ( *ptr != '\0' && *ptr != '\n' && *ptr != *s )
       ptr++;

   if ( *ptr == '\0' || *ptr == '\n' )
       return( 0 );

   if ( strlen( ptr ) < strlen( s ) )
       return( 0 );

   if ( strncmp( ptr, s, strlen(s) ) == 0 )
       return( 1 );

   return( 0 );
}



char *_getStrValue( char *str, int *length )
{
  int s, l;

  *length = -1;
  if ( str == (char*)NULL ) return( (char*)NULL );

  for ( s=0; s<(int)strlen(str) && str[s] != '\''; s++ )
      ;
  s++;
  if ( s >= (int)strlen(str) ) return( (char*)NULL );

  for ( l=s+1; l<(int)strlen(str) && str[l] != '\''; l++ )
      ;
  l--;
  if ( l >= (int)strlen(str) ) return( (char*)NULL );

  *length = l-s+1;
  return( &(str[s]) );
}



void _comma2point( char *thestr, char *resstr )
{
  char *s1 = thestr;
  char *s2 = resstr;

  while ( *s1 != '\0' ) {
    if ( *s1 != ',' ) {
      *s2++ = *s1++;
      continue;
    }
    if ( s1 == thestr ) {
      *s2++ = *s1++;
      continue;
    }
    if ( '0' <= s1[-1] && s1[-1] <= '9' && '0' <= s1[1] && s1[1] <= '9') {
      *s2++ = '.';
      s1++;
    }
    else {
      *s2++ = *s1++;
    }
  }
}



static int _XML_ReadCellSequenceProperty( FILE *f, typeCellSequence *l, typePropertyList *o,
                                          enumProperty property, char *endstring )
{
  char *proc = "_XML_ReadCellSequenceProperty";
  char str[2048], str2[2048], *ptr;

  int cellid, cellid2;
  int cellnb, acquisition_time;
  typeCell *cell;
  typeCell *daughterCell, *historyCell;

  int nbread;
  int i;
  int volume;
  double dvect[9];
  float fval;
  int ival;
  char *ptrstrval;
  char strval[256];



  /* add property to list
   */
  if ( _addPropertyToList( o, property ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to add property to list\n", proc );
    return( -1 );
  }


  while ( fgets( str, 2048, f ) != NULL ) {

    /* is reading done
     */
    if ( _hasString( str, endstring ) == 1 )
      return( 1 );

    ptr = str;
    while ( *ptr != '\0' && *ptr != '\n' && *ptr != '<' )
        ptr++;
    if ( *ptr == '\0' || *ptr == '\n' ) {
        if ( _verbose_ ) {
          fprintf( stderr, "%s: weird string\n", proc );
          fprintf( stderr, "\t %s", str );
        }
        return( -1 );
    }

    /* not done? read property
     */
    switch( property ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such property not handled yet\n", proc );
      return( -1 );

    case _LINEAGE_ :
      if ( sscanf( ptr, "<cell cell-id=\"%d\">", &cellid ) != 1 ) {
        if ( _verbose_ ) {
          fprintf( stderr, "%s: unable to translate string\n", proc );
          fprintf( stderr, "\t %s", ptr );
        }
        return( -1 );
      }
      cell = _getCellFromSequence( l, cellid );
      if ( cell == (typeCell*)NULL ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to get cell when reading lineage\n", proc );
        return( -1 );
      }

      while ( *ptr != '\0' && *ptr != '\n' && *ptr != '[' )
          ptr++;
      if ( *ptr == '\0' || *ptr == '\n' ) {
          if ( _verbose_ ) {
            fprintf( stderr, "%s: weird string (lineage #1)\n", proc );
            fprintf( stderr, "\t %s", str );
          }
          return( -1 );
      }
      ptr++;

      do {
        if ( sscanf( ptr, "%d", &cellid2 ) != 1 ) {
          if ( _verbose_ ) {
            fprintf( stderr, "%s: unable to translate string\n", proc );
            fprintf( stderr, "\t %s", ptr );
          }
          return( -1 );
        }

        /* create the daughter cell
         * else it won't be if only the lineage is read
         */
        daughterCell = _getCellFromSequence( l, cellid2 );
        if ( daughterCell == (typeCell*)NULL ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to get daughter cell when reading lineage\n", proc );
          return( -1 );
        }

        if ( _addTemporalNeighborToCellLineage( cell, cellid2, 1 ) != 1 ) {
          if ( _verbose_ ) {
            fprintf( stderr, "%s: unable to add daughter #%d to cell\n", proc, cellid2 );
          }
          return( -1 );
        }

        while ( '0' <= *ptr && *ptr <= '9' ) ptr++;

        if ( *ptr == ',' ) {
          ptr ++;
          while ( *ptr != '\0' && *ptr != '\n' && (*ptr < '0' || '9' < *ptr) )
              ptr++;
          if ( *ptr == '\0' || *ptr == '\n' ) {
              if ( _verbose_ ) {
                fprintf( stderr, "%s: weird string (lineage #2)\n", proc );
                fprintf( stderr, "\t %s", str );
              }
              return( -1 );
          }
        }
        else if ( *ptr == ']' ) {
            break;
        }
        else {
          if ( _verbose_ ) {
            fprintf( stderr, "%s: weird string (lineage #3)\n", proc );
            fprintf( stderr, "\t str='%s' ptr='%s'", str, ptr );
          }
          return( -1 );
        }
      } while ( 1 );
      break;

    case _FORWARD_NEIGHBORS_ :
      if ( sscanf( ptr, "<cell cell-id=\"%d\">", &cellid ) != 1 ) {
        if ( _verbose_ ) {
          fprintf( stderr, "%s: unable to translate string\n", proc );
          fprintf( stderr, "\t %s", ptr );
        }
        return( -1 );
      }
      cell = _getCellFromSequence( l, cellid );
      if ( cell == (typeCell*)NULL ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to get cell when reading forward neighbors\n", proc );
        return( -1 );
      }

      do {

        if ( fgets( str, 512, f ) == NULL ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unexpected end of file while reading forward neighbors\n", proc );
          return( -1 );
        }

        ptr = str;
        while ( *ptr != '\0' && *ptr != '\n' && *ptr != '<' )
            ptr++;
        if ( *ptr == '\0' || *ptr == '\n' ) {
            if ( _verbose_ ) {
              fprintf( stderr, "%s: weird string\n", proc );
              fprintf( stderr, "\t %s", str );
            }
            return( -1 );
        }

        if ( _hasString( str, "</cell>" ) == 1 ) break;

        if ( sscanf( ptr, "<cell cell-id=\"%d\">%d</cell>", &cellid2, &ival ) != 2 ) {
          if ( _verbose_ ) {
            fprintf( stderr, "%s: unable to translate string\n", proc );
            fprintf( stderr, "\t %s", ptr );
          }
          return( -1 );
        }
        if ( _addForwardNeighborToCell( cell, cellid2, ival ) != 1 ) {
          if ( _verbose_ ) {
            fprintf( stderr, "%s: unable to add forward neighbor #%d to cell\n", proc, cellid2 );
          }
          return( -1 );
        }
      } while ( 1 );
      break;

    case _BACKWARD_NEIGHBORS_ :
      if ( sscanf( ptr, "<cell cell-id=\"%d\">", &cellid ) != 1 ) {
        if ( _verbose_ ) {
          fprintf( stderr, "%s: unable to translate string\n", proc );
          fprintf( stderr, "\t %s", ptr );
        }
        return( -1 );
      }
      cell = _getCellFromSequence( l, cellid );
      if ( cell == (typeCell*)NULL ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to get cell when reading backward neighbors\n", proc );
        return( -1 );
      }

      do {

        if ( fgets( str, 512, f ) == NULL ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unexpected end of file while reading backward neighbors\n", proc );
          return( -1 );
        }

        ptr = str;
        while ( *ptr != '\0' && *ptr != '\n' && *ptr != '<' )
            ptr++;
        if ( *ptr == '\0' || *ptr == '\n' ) {
            if ( _verbose_ ) {
              fprintf( stderr, "%s: weird string\n", proc );
              fprintf( stderr, "\t %s", str );
            }
            return( -1 );
        }

        if ( _hasString( str, "</cell>" ) == 1 ) break;

        if ( sscanf( ptr, "<cell cell-id=\"%d\">%d</cell>", &cellid2, &ival ) != 2 ) {
          if ( _verbose_ ) {
            fprintf( stderr, "%s: unable to translate string\n", proc );
            fprintf( stderr, "\t %s", ptr );
          }
          return( -1 );
        }
        if ( _addBackwardNeighborToCell( cell, cellid2, ival ) != 1 ) {
          if ( _verbose_ ) {
            fprintf( stderr, "%s: unable to add backward neighbor #%d to cell\n", proc, cellid2 );
          }
          return( -1 );
        }
      } while ( 1 );
      break;

    case _H_MIN_ :
      if ( _verbose_ )
        fprintf( stderr, "%s: such property not handled yet\n", proc );
      return( -1 );

    case _VOLUME_ :
      if ( sscanf( ptr, "<cell cell-id=\"%d\">%d</cell>", &cellid, &volume ) != 2 ) {
        if ( _verbose_ ) {
          fprintf( stderr, "%s: unable to translate string\n", proc );
          fprintf( stderr, "\t %s", ptr );
        }
        return( -1 );
      }
      cell = _getCellFromSequence( l, cellid );
      if ( cell == (typeCell*)NULL ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to get cell when reading volumes\n", proc );
        return( -1 );
      }
      cell->npts = volume;
      break;

    case _SURFACE_ :
      _comma2point( ptr, str2 );
      if ( sscanf( str2, "<cell cell-id=\"%d\">%f</cell>", &cellid, &fval ) != 2 ) {
        if ( _verbose_ ) {
          fprintf( stderr, "%s: unable to translate string\n", proc );
          fprintf( stderr, "\t %s", str2 );
        }
        return( -1 );
      }
      cell = _getCellFromSequence( l, cellid );
      if ( cell == (typeCell*)NULL ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to get cell when reading surfaces\n", proc );
        return( -1 );
      }
      cell->surface = fval;
      break;

    case _COMPACTNESS_ :
      _comma2point( ptr, str2 );
      if ( sscanf( str2, "<cell cell-id=\"%d\">%f</cell>", &cellid, &fval ) != 2 ) {
        if ( _verbose_ ) {
          fprintf( stderr, "%s: unable to translate string\n", proc );
          fprintf( stderr, "\t %s", str2 );
        }
        return( -1 );
      }
      cell = _getCellFromSequence( l, cellid );
      if ( cell == (typeCell*)NULL ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to get cell when reading compactness\n", proc );
        return( -1 );
      }
      cell->compactness = fval;
      break;

    case _SIGMA_ :
    case _LABEL_IN_TIME_ :
        if ( _verbose_ )
          fprintf( stderr, "%s: such property not handled yet\n", proc );
        return( -1 );

    case _BARYCENTER_ :
      _comma2point( ptr, str2 );
      if ( sscanf( str2, "<cell cell-id=\"%d\">[%lf, %lf, %lf]</cell>",
                   &cellid, &(dvect[0]), &(dvect[1]), &(dvect[2]) ) != 4 ) {
        if ( _verbose_ ) {
          fprintf( stderr, "%s: unable to translate string\n", proc );
          fprintf( stderr, "\t %s", str2 );
        }
        return( -1 );
      }
      cell = _getCellFromSequence( l, cellid );
      if ( cell == (typeCell*)NULL ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to get cell when reading barycenters\n", proc );
        return( -1 );
      }
      for ( i=0; i<3; i++ ) cell->voxel_barycenter[i] = dvect[i];
      break;

    case _PRINCIPAL_VALUE_ :
      _comma2point( ptr, str2 );
      if ( sscanf( str2, "<cell cell-id=\"%d\">[%lf, %lf, %lf]</cell>",
                   &cellid, &(dvect[0]), &(dvect[1]), &(dvect[2]) ) != 4 ) {
        if ( _verbose_ ) {
          fprintf( stderr, "%s: unable to translate string\n", proc );
          fprintf( stderr, "\t %s", str2 );
        }
        return( -1 );
      }
      cell = _getCellFromSequence( l, cellid );
      if ( cell == (typeCell*)NULL ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to get cell when reading principal values\n", proc );
        return( -1 );
      }
      for ( i=0; i<3; i++ ) cell->voxel_eigenvalues[i] = dvect[i];
      break;

    case _PRINCIPAL_VECTOR_ :
      _comma2point( ptr, str2 );
      nbread = sscanf( str2, "<cell cell-id=\"%d\">[[%lf, %lf, %lf], [%lf, %lf, %lf], [%lf, %lf, %lf]]</cell>",
                       &cellid, &(dvect[0]), &(dvect[3]), &(dvect[6]),
                       &(dvect[1]), &(dvect[4]), &(dvect[7]),
                       &(dvect[2]), &(dvect[5]), &(dvect[8]) );
      if ( nbread != 10 ) {
        if ( _verbose_ ) {
          fprintf( stderr, "%s: unable to translate string\n", proc );
          fprintf( stderr, "\t %s", str2 );
        }
        return( -1 );
      }
      cell = _getCellFromSequence( l, cellid );
      if ( cell == (typeCell*)NULL ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to get cell when reading principal vectors\n", proc );
        return( -1 );
      }
      for ( i=0; i<9; i++ ) cell->voxel_eigenvectors[i] = dvect[i];
      break;

    case _FATE_ :
        if ( sscanf( ptr, "<cell cell-id=\"%d\">'%s'</cell>", &cellid, strval ) != 2 ) {
          if ( _verbose_ ) {
            fprintf( stderr, "%s: unable to translate string\n", proc );
            fprintf( stderr, "\t %s", ptr );
          }
          return( -1 );
        }
        cell = _getCellFromSequence( l, cellid );
        if ( cell == (typeCell*)NULL ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to get cell when reading volumes\n", proc );
          return( -1 );
        }
        ptrstrval = _getStrValue( str, &ival );
        if ( ptrstrval == (char*)NULL ) {
            if ( _verbose_ ) {
              fprintf( stderr, "%s: empty fate ?!\n", proc );
              fprintf( stderr, "\t %s", ptr );
            }
            return( -1 );
        }
        if ( cell->fate != (char*)NULL ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: fate already exists ?!\n", proc );
            return( -1 );
        }
        cell->fate = (char*)vtmalloc( ival+1, "cell->fate", proc );
        if ( cell->fate == (char*)NULL ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: allocation failed\n", proc );
            return( -1 );
        }
        (void)strncpy( cell->fate, ptrstrval, ival );
        cell->fate[ival] = '\0';
        break;

    case _FATE2_ :
        if ( sscanf( ptr, "<cell cell-id=\"%d\">'%s'</cell>", &cellid, strval ) != 2 ) {
          if ( _verbose_ ) {
            fprintf( stderr, "%s: unable to translate string\n", proc );
            fprintf( stderr, "\t %s", ptr );
          }
          return( -1 );
        }
        cell = _getCellFromSequence( l, cellid );
        if ( cell == (typeCell*)NULL ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to get cell when reading volumes\n", proc );
          return( -1 );
        }
        ptrstrval = _getStrValue( str, &ival );
        if ( ptrstrval == (char*)NULL ) {
            if ( _verbose_ ) {
              fprintf( stderr, "%s: empty fate ?!\n", proc );
              fprintf( stderr, "\t %s", ptr );
            }
            return( -1 );
        }
        if ( cell->fate2 != (char*)NULL ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: fate already exists ?!\n", proc );
            return( -1 );
        }
        cell->fate2 = (char*)vtmalloc( ival+1, "cell->fate2", proc );
        if ( cell->fate2 == (char*)NULL ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: allocation failed\n", proc );
            return( -1 );
        }
        (void)strncpy( cell->fate2, ptrstrval, ival );
        cell->fate2[ival] = '\0';
        break;

    case _FATE3_ :
        if ( sscanf( ptr, "<cell cell-id=\"%d\">'%s'</cell>", &cellid, strval ) != 2 ) {
          if ( _verbose_ ) {
            fprintf( stderr, "%s: unable to translate string\n", proc );
            fprintf( stderr, "\t %s", ptr );
          }
          return( -1 );
        }
        cell = _getCellFromSequence( l, cellid );
        if ( cell == (typeCell*)NULL ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to get cell when reading volumes\n", proc );
          return( -1 );
        }
        ptrstrval = _getStrValue( str, &ival );
        if ( ptrstrval == (char*)NULL ) {
            if ( _verbose_ ) {
              fprintf( stderr, "%s: empty fate ?!\n", proc );
              fprintf( stderr, "\t %s", ptr );
            }
            return( -1 );
        }
        if ( cell->fate3 != (char*)NULL ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: fate already exists ?!\n", proc );
            return( -1 );
        }
        cell->fate3 = (char*)vtmalloc( ival+1, "cell->fate3", proc );
        if ( cell->fate3 == (char*)NULL ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: allocation failed\n", proc );
            return( -1 );
        }
        (void)strncpy( cell->fate3, ptrstrval, ival );
        cell->fate3[ival] = '\0';
        break;

    case _FATE4_ :
        if ( sscanf( ptr, "<cell cell-id=\"%d\">'%s'</cell>", &cellid, strval ) != 2 ) {
          if ( _verbose_ ) {
            fprintf( stderr, "%s: unable to translate string\n", proc );
            fprintf( stderr, "\t %s", ptr );
          }
          return( -1 );
        }
        cell = _getCellFromSequence( l, cellid );
        if ( cell == (typeCell*)NULL ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to get cell when reading volumes\n", proc );
          return( -1 );
        }
        ptrstrval = _getStrValue( str, &ival );
        if ( ptrstrval == (char*)NULL ) {
            if ( _verbose_ ) {
              fprintf( stderr, "%s: empty fate ?!\n", proc );
              fprintf( stderr, "\t %s", ptr );
            }
            return( -1 );
        }
        if ( cell->fate4 != (char*)NULL ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: fate already exists ?!\n", proc );
            return( -1 );
        }
        cell->fate4 = (char*)vtmalloc( ival+1, "cell->fate4", proc );
        if ( cell->fate4 == (char*)NULL ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: allocation failed\n", proc );
            return( -1 );
        }
        (void)strncpy( cell->fate4, ptrstrval, ival );
        cell->fate4[ival] = '\0';
        break;

    case _ALL_CELLS_ :
        if ( _verbose_ )
          fprintf( stderr, "%s: such property not handled yet\n", proc );
        return( -1 );

    case _NAME_ :
        if ( sscanf( ptr, "<cell cell-id=\"%d\">'%s'</cell>", &cellid, strval ) != 2 ) {
          if ( _verbose_ ) {
            fprintf( stderr, "%s: unable to translate string\n", proc );
            fprintf( stderr, "\t %s", ptr );
          }
          return( -1 );
        }
        cell = _getCellFromSequence( l, cellid );
        if ( cell == (typeCell*)NULL ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to get cell when reading volumes\n", proc );
          return( -1 );
        }
        ptrstrval = _getStrValue( str, &ival );
        if ( ptrstrval == (char*)NULL ) {
            if ( _verbose_ ) {
              fprintf( stderr, "%s: empty name ?!\n", proc );
              fprintf( stderr, "\t %s", ptr );
            }
            return( -1 );
        }
        if ( cell->name != (char*)NULL ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: name already exists ?!\n", proc );
            return( -1 );
        }
        cell->name = (char*)vtmalloc( ival+1, "cell->name", proc );
        if ( cell->name == (char*)NULL ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: allocation failed\n", proc );
            return( -1 );
        }
        (void)strncpy( cell->name, ptrstrval, ival );
        cell->name[ival] = '\0';
        break;

    case _HISTORY_ :
        if ( sscanf( ptr, "<cell cell-id=\"%d\">", &cellid ) != 1 ) {
          if ( _verbose_ ) {
            fprintf( stderr, "%s: unable to translate string\n", proc );
            fprintf( stderr, "\t %s", ptr );
          }
          return( -1 );
        }
        cell = _getCellFromSequence( l, cellid );
        if ( cell == (typeCell*)NULL ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to get cell when reading lineage\n", proc );
          return( -1 );
        }

        while ( *ptr != '\0' && *ptr != '\n' && *ptr != '[' )
            ptr++;
        if ( *ptr == '\0' || *ptr == '\n' ) {
            if ( _verbose_ ) {
              fprintf( stderr, "%s: weird string (history #1)\n", proc );
              fprintf( stderr, "\t %s", str );
            }
            return( -1 );
        }
        ptr++;

        do {
          if ( sscanf( ptr, "%d", &cellid2 ) != 1 ) {
            if ( _verbose_ ) {
              fprintf( stderr, "%s: unable to translate string\n", proc );
              fprintf( stderr, "\t %s", ptr );
            }
            return( -1 );
          }

          /* create the daughter cell
           * else it won't be if only the lineage is read
           */
          historyCell = _getCellFromSequence( l, cellid2 );
          if ( historyCell == (typeCell*)NULL ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to get daughter cell when reading lineage\n", proc );
            return( -1 );
          }

          if ( _addTemporalNeighborToCellHistory( cell, cellid2, 1 ) != 1 ) {
            if ( _verbose_ ) {
              fprintf( stderr, "%s: unable to add history #%d to cell\n", proc, cellid2 );
            }
            return( -1 );
          }

          while ( '0' <= *ptr && *ptr <= '9' ) ptr++;

          if ( *ptr == ',' ) {
            ptr ++;
            while ( *ptr != '\0' && *ptr != '\n' && (*ptr < '0' || '9' < *ptr) )
                ptr++;
            if ( *ptr == '\0' || *ptr == '\n' ) {
                if ( _verbose_ ) {
                  fprintf( stderr, "%s: weird string (history #2)\n", proc );
                  fprintf( stderr, "\t %s", str );
                }
                return( -1 );
            }
          }
          else if ( *ptr == ']' ) {
              break;
          }
          else {
            if ( _verbose_ ) {
              fprintf( stderr, "%s: weird string (history #3)\n", proc );
              fprintf( stderr, "\t str='%s' ptr='%s'", str, ptr );
            }
            return( -1 );
          }
        } while ( 1 );
        break;

    case _CONTACT_SURFACE_ :
      if ( sscanf( ptr, "<cell cell-id=\"%d\">", &cellid ) != 1 ) {
        if ( _verbose_ ) {
          fprintf( stderr, "%s: unable to translate string\n", proc );
          fprintf( stderr, "\t %s", ptr );
        }
        return( -1 );
      }
      cell = _getCellFromSequence( l, cellid );
      if ( cell == (typeCell*)NULL ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to get cell when reading contact surfaces\n", proc );
        return( -1 );
      }

      do {

        if ( fgets( str, 512, f ) == NULL ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unexpected end of file while reading contact surfaces\n", proc );
          return( -1 );
        }

        ptr = str;
        while ( *ptr != '\0' && *ptr != '\n' && *ptr != '<' )
            ptr++;
        if ( *ptr == '\0' || *ptr == '\n' ) {
            if ( _verbose_ ) {
              fprintf( stderr, "%s: weird string\n", proc );
              fprintf( stderr, "\t %s", str );
            }
            return( -1 );
        }

        if ( _hasString( str, "</cell>" ) == 1 ) break;

        _comma2point( ptr, str2 );
        if ( sscanf( str2, "<cell cell-id=\"%d\">%f</cell>", &cellid2, &fval ) != 2 ) {
          if ( _verbose_ ) {
            fprintf( stderr, "%s: unable to translate string\n", proc );
            fprintf( stderr, "\t %s", str2 );
          }
          return( -1 );
        }
        _cellLocalId( cellid2, &cellnb, &acquisition_time );
        if ( _addSpatialNeighborToCell( cell, cellnb, fval ) != 1 ) {
          if ( _verbose_ ) {
            fprintf( stderr, "%s: unable to add neighbor #%d to cell\n", proc, cellnb );
          }
          return( -1 );
        }
      } while ( 1 );
      break;

    }
  }

  return( 1 );
}





int XML_FscanfCellSequenceProperty( FILE *f, typeCellSequence *l, typePropertyList *o )
{
  char *proc = "XML_FscanfCellSequenceProperty";
  char str[2048];
  int ret;
  int line = 0;


  while ( fgets( str, 2048, f ) != NULL ) {

    /* skip comment (to be done)
     */
    if ( _hasString( str, "</data>" ) == 1 ) {
        return( 0 );
    }
    else if ( _hasString( str, "<data>" ) == 1 ) {
      /* loop to read single properties
       */
      do {
        ret = XML_FscanfCellSequenceProperty( f, l, o );
      } while ( ret == 1 );
      if ( ret == 0 ) break;
      if ( ret <= 0 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: some error occurs when parsing file\n", proc );
        return( -1 );
      }
    }
    else if ( _hasString( str, "<lineage_tree>" ) == 1 ) {
      if ( _debug_ >= 2 ) fprintf( stderr, "%s: reading <lineage_tree>\n", proc );
      if ( _XML_ReadCellSequenceProperty( f, l, o, _LINEAGE_, "</lineage_tree>" ) != 1 ) {
        BAL_FreeCellSequence( l );
        if ( _verbose_ )
          fprintf( stderr, "%s: ... error when reading lineage tree\n", proc );
        return( -1 );
      }
    }
    else if ( _hasString( str, "<cell_lineage>" ) == 1 ) {
      if ( _debug_ >= 2 ) fprintf( stderr, "%s: reading <cell_lineage>\n", proc );
      if ( _XML_ReadCellSequenceProperty( f, l, o, _LINEAGE_, "</cell_lineage>" ) != 1 ) {
        BAL_FreeCellSequence( l );
        if ( _verbose_ )
          fprintf( stderr, "%s: ... error when reading lineage tree\n", proc );
        return( -1 );
      }
    }
    else if ( _hasString( str, "<cell_forward_neighbors>" ) == 1 ) {
      if ( _debug_ >= 2 ) fprintf( stderr, "%s: reading <cell_forward_neighbors>\n", proc );
      if ( _XML_ReadCellSequenceProperty( f, l, o, _FORWARD_NEIGHBORS_, "</cell_forward_neighbors>" ) != 1 ) {
        BAL_FreeCellSequence( l );
        if ( _verbose_ )
          fprintf( stderr, "%s: ... error when reading forward neighbors\n", proc );
        return( -1 );
      }
    }
    else if ( _hasString( str, "<cell_backward_neighbors>" ) == 1 ) {
      if ( _debug_ >= 2 ) fprintf( stderr, "%s: reading <cell_backward_neighbors>\n", proc );
      if ( _XML_ReadCellSequenceProperty( f, l, o, _BACKWARD_NEIGHBORS_, "</cell_backward_neighbors>" ) != 1 ) {
        BAL_FreeCellSequence( l );
        if ( _verbose_ )
          fprintf( stderr, "%s: ... error when reading backward neighbors\n", proc );
        return( -1 );
      }
    }
    else if ( _hasString( str, "<cell_h_min>" ) == 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: ... <cell_h_min> not handled yet\n", proc );
      return( -1 );
    }
    else if ( _hasString( str, "<cell_volume>" ) == 1 ) {
      if ( _debug_ >= 2 ) fprintf( stderr, "%s: reading <cell_volume>\n", proc );
      if ( _XML_ReadCellSequenceProperty( f, l, o, _VOLUME_, "</cell_volume>" ) != 1 ) {
        BAL_FreeCellSequence( l );
        if ( _verbose_ )
          fprintf( stderr, "%s: ... error when reading volumes\n", proc );
        return( -1 );
      }
    }
    else if ( _hasString( str, "<cell_surface>" ) == 1 ) {
      if ( _debug_ >= 2 ) fprintf( stderr, "%s: reading <cell_surface>\n", proc );
      if ( _XML_ReadCellSequenceProperty( f, l, o, _SURFACE_, "</cell_surface>" ) != 1 ) {
        BAL_FreeCellSequence( l );
        if ( _verbose_ )
          fprintf( stderr, "%s: ... error when reading surfaces\n", proc );
        return( -1 );
      }
    }
    else if ( _hasString( str, "<cell_compactness>" ) == 1 ) {
      if ( _debug_ >= 2 ) fprintf( stderr, "%s: reading <cell_compactness>\n", proc );
      if ( _XML_ReadCellSequenceProperty( f, l, o, _COMPACTNESS_, "</cell_compactness>" ) != 1 ) {
        BAL_FreeCellSequence( l );
        if ( _verbose_ )
          fprintf( stderr, "%s: ... error when reading compactness\n", proc );
        return( -1 );
      }
    }
    else if ( _hasString( str, "<cell_sigma>" ) == 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: ... <cell_sigma> not handled yet\n", proc );
        return( -1 );
    }
    else if ( _hasString( str, "<cell_labels_in_time>" ) == 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: ... <cell_labels_in_time> not handled yet\n", proc );
      return( -1 );
    }
    else if ( _hasString( str, "<cell_barycenter>" ) == 1 ) {
      if ( _debug_ >= 2 ) fprintf( stderr, "%s: reading <cell_barycenter>\n", proc );
      if ( _XML_ReadCellSequenceProperty( f, l, o, _BARYCENTER_, "</cell_barycenter>" ) != 1 ) {
        BAL_FreeCellSequence( l );
        if ( _verbose_ )
          fprintf( stderr, "%s: ... error when reading barycenters\n", proc );
        return( -1 );
      }
    }
    else if ( _hasString( str, "<cell_principal_values>" ) == 1 ) {
      if ( _debug_ >= 2 ) fprintf( stderr, "%s: reading <cell_principal_values>\n", proc );
      if ( _XML_ReadCellSequenceProperty( f, l, o, _PRINCIPAL_VALUE_, "</cell_principal_values>" ) != 1 ) {
        BAL_FreeCellSequence( l );
        if ( _verbose_ )
          fprintf( stderr, "%s: ... error when reading principal values\n", proc );
        return( -1 );
      }
    }
    else if ( _hasString( str, "<cell_principal_vectors>" ) == 1 ) {
      if ( _debug_ >= 2 ) fprintf( stderr, "%s: reading <cell_principal_vectors>\n", proc );
      if ( _XML_ReadCellSequenceProperty( f, l, o, _PRINCIPAL_VECTOR_, "</cell_principal_vectors>" ) != 1 ) {
        BAL_FreeCellSequence( l );
        if ( _verbose_ )
          fprintf( stderr, "%s: ... error when reading principal vectors\n", proc );
        return( -1 );
      }
    }
    else if ( _hasString( str, "<cell_fate>" ) == 1 ) {
        if ( _debug_ >= 2 ) fprintf( stderr, "%s: reading <cell_fate>\n", proc );
        if ( _XML_ReadCellSequenceProperty( f, l, o, _FATE_, "</cell_fate>" ) != 1 ) {
          BAL_FreeCellSequence( l );
          if ( _verbose_ )
            fprintf( stderr, "%s: ... error when reading cell fate\n", proc );
          return( -1 );
        }
    }
    else if ( _hasString( str, "<cell_fate_2>" ) == 1 ) {
        if ( _debug_ >= 2 ) fprintf( stderr, "%s: reading <cell_fate_2>\n", proc );
        if ( _XML_ReadCellSequenceProperty( f, l, o, _FATE2_, "</cell_fate_2>" ) != 1 ) {
          BAL_FreeCellSequence( l );
          if ( _verbose_ )
            fprintf( stderr, "%s: ... error when reading cell fate #2\n", proc );
          return( -1 );
        }
    }
    else if ( _hasString( str, "<cell_fate_3>" ) == 1 ) {
        if ( _debug_ >= 2 ) fprintf( stderr, "%s: reading <cell_fate_3>\n", proc );
        if ( _XML_ReadCellSequenceProperty( f, l, o, _FATE3_, "</cell_fate_3>" ) != 1 ) {
          BAL_FreeCellSequence( l );
          if ( _verbose_ )
            fprintf( stderr, "%s: ... error when reading cell fate #3\n", proc );
          return( -1 );
        }
    }
    else if ( _hasString( str, "<cell_fate_4>" ) == 1 ) {
        if ( _debug_ >= 2 ) fprintf( stderr, "%s: reading <cell_fate_4>\n", proc );
        if ( _XML_ReadCellSequenceProperty( f, l, o, _FATE4_, "</cell_fate_4>" ) != 1 ) {
          BAL_FreeCellSequence( l );
          if ( _verbose_ )
            fprintf( stderr, "%s: ... error when reading cell fate #4\n", proc );
          return( -1 );
        }
    }
    else if ( _hasString( str, "<all_cells>" ) == 1 ) {
      if ( _addPropertyToList( o, _ALL_CELLS_ ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to add property to list\n", proc );
        return( -1 );
      }
      do {
        if ( fgets( str, 512, f ) == NULL ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unexpected end of file while reading all cells\n", proc );
          return( -1 );
        }
        if ( _hasString( str, "</all_cells>" ) == 1 ) break;
      } while ( 1 );
    }
    else if ( _hasString( str, "<cell_name>" ) == 1 ) {
        if ( _debug_ >= 2 ) fprintf( stderr, "%s: reading <cell_name>\n", proc );
        if ( _XML_ReadCellSequenceProperty( f, l, o, _NAME_, "</cell_name>" ) != 1 ) {
          BAL_FreeCellSequence( l );
          if ( _verbose_ )
            fprintf( stderr, "%s: ... error when reading cell name\n", proc );
          return( -1 );
        }
    }
    else if ( _hasString( str, "<cell_history>" ) == 1 ) {
        if ( _debug_ >= 2 ) fprintf( stderr, "%s: reading <cell_history>\n", proc );
        if ( _XML_ReadCellSequenceProperty( f, l, o, _HISTORY_, "</cell_history>" ) != 1 ) {
          BAL_FreeCellSequence( l );
          if ( _verbose_ )
            fprintf( stderr, "%s: ... error when reading cell history\n", proc );
          return( -1 );
        }
    }
    else if ( _hasString( str, "<cell_contact_surface>" ) == 1 ) {
      if ( _debug_ >= 2 ) fprintf( stderr, "%s: reading <cell_contact_surface>\n", proc );
      if ( _XML_ReadCellSequenceProperty( f, l, o, _CONTACT_SURFACE_, "</cell_contact_surface>" ) != 1 ) {
        BAL_FreeCellSequence( l );
        if ( _verbose_ )
          fprintf( stderr, "%s: ... error when reading contact surfaces\n", proc );
        return( -1 );
      }
    }
    else {
        if ( _verbose_ ) {
          fprintf( stderr, "%s: ... unknown property\n", proc );
          fprintf( stderr, "\t reading %s\n", str );
        }
        return( -1 );
    }

    line ++;
  }

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "%s: %d cell images were indexed\n", proc, l->n_data );
  }


  /* do some post-processing
   */
  l->data[ l->n_data-1 ].next_acquisition_time = 0;

  /*
  if ( _isPropertyInList( o, _LINEAGE_ ) == 1 ) {
    BAL_CellSequenceBackwardLineage( l );

    for ( i=0; i<l->n_data; i++ ) {
      if ( l->data[i].n_data <= 0 ) continue;
      for ( c=0; c<l->data[i].n_data; c++ ) {
        l->data[i].data[c].mothers.n_selected_data = l->data[i].data[c].mothers.n_data;
        l->data[i].data[c].daughters.n_selected_data = l->data[i].data[c].daughters.n_data;
      }
    }
  }
  */

  return( 1 );
}











