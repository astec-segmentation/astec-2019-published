/*************************************************************************
 * bal-cellProperties.h -
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

#ifndef _bal_cellproperties_h_
#define _bal_cellproperties_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <bal-image.h>


extern void BAL_SetVerboseInBalCellProperties( int v );
extern void BAL_IncrementVerboseInBalCellProperties(  );
extern void BAL_DecrementVerboseInBalCellProperties(  );


extern void BAL_SetDebugInBalCellProperties( int d );
extern void BAL_IncrementDebugInBalCellProperties(  );
extern void BAL_DecrementDebugInBalCellProperties(  );



typedef struct typeUpdateList {
  int n_allocated_data;
  int n_data;
  int *data;
} typeUpdateList;

extern void BAL_InitUpdateList( typeUpdateList *l );
extern void BAL_FreeUpdateList( typeUpdateList *l );
extern int BAL_AddUpdateToList( typeUpdateList *l, int u );
extern void BAL_FprintfUpdateList( FILE *f, typeUpdateList *l );
extern int _isUpdateInList( typeUpdateList *l, int u );


typedef enum enumProperty {
  _NONE_,
  _LINEAGE_,
  _FORWARD_NEIGHBORS_,
  _BACKWARD_NEIGHBORS_,
  _H_MIN_,
  _VOLUME_,
  _SURFACE_,
  _COMPACTNESS_,
  _SIGMA_,
  _LABEL_IN_TIME_,
  _BARYCENTER_,
  _PRINCIPAL_VALUE_,
  _PRINCIPAL_VECTOR_,
  _FATE_,
  _FATE2_,
  _FATE3_,
  _FATE4_,
  _ALL_CELLS_,
  _NAME_,
  _HISTORY_,
  _CONTACT_SURFACE_
} enumProperty;

#define _MAX_PROPERTIES_ 30

typedef struct typePropertyList {
  int n_allocated_data;
  int n_data;
  enumProperty data[_MAX_PROPERTIES_];
} typePropertyList;

extern void BAL_InitPropertyList( typePropertyList *o );
extern void BAL_DefaultPropertyList( typePropertyList *o );
extern void _removePropertyFromList( typePropertyList *l, enumProperty p );
extern void BAL_FprintfPropertyList( FILE *f, typePropertyList *o, char *str );
extern int _isPropertyInList( typePropertyList *o, enumProperty p );
extern void BAL_GetWritePropertyList( typePropertyList *writeList, typePropertyList *readList, typePropertyList *outputList );





/* to count the contact surface
 */
typedef struct typeSpatialNeighbor {
  int label;
  float surface;
} typeSpatialNeighbor;



typedef struct typeSpatialNeighborList {
  typeSpatialNeighbor *data;
  int n_data;
  int n_allocated_data;
} typeSpatialNeighborList;



/* for the lineage
 */
typedef struct typeTemporalNeighbor {
  int label;
  int time;
  int globalId;
  int npts;
} typeTemporalNeighbor;



typedef struct typeTemporalNeighborList {
  typeTemporalNeighbor *data;
  int n_data;
  int n_allocated_data;
} typeTemporalNeighborList;



/* cell structure
 */
typedef struct typeCell {

  /* bounding box
   */
  int ptmin[3];
  int ptmax[3];

  /* #points
   */
  int npts;
  float surface;
  float compactness;

  /* inertia parameters
   */
  double voxel_barycenter[3];
  double voxel_covariance_matrix[9];

  double voxel_eigenvalues[3];
  double voxel_eigenvectors[9];

  /* neighbor list
   */
  typeSpatialNeighborList neighbors;

  /* lineage lists
   */
  typeTemporalNeighborList forward;
  typeTemporalNeighborList backward;

  typeTemporalNeighborList lineage;
  typeTemporalNeighborList history;

  /* fates
   */
  char *fate;
  char *fate2;
  char *fate3;
  char *fate4;

  /* name
   */
  char *name;

} typeCell;



/* cell list for one image
 */

typedef struct typeCellImage {
  typeCell *data;
  int n_data;
  int n_allocated_data;
  int previous_acquisition_time;
  int acquisition_time;
  int next_acquisition_time;
} typeCellImage;

extern void BAL_InitCellImage( typeCellImage *l );
extern void BAL_FreeCellImage( typeCellImage *l );
extern int BAL_AllocCellImage( typeCellImage *l , int n );
extern int BAL_InitCellImageFromImage( typeCellImage *theList,  bal_image *theIm );



typedef struct typeCellSequence {
  typeCellImage *data;
  int n_data;
  int n_allocated_data;
} typeCellSequence;

extern void BAL_InitCellSequence( typeCellSequence *l );
extern void BAL_FreeCellSequence( typeCellSequence *l );
extern int BAL_AllocCellSequence( typeCellSequence *l , int n );
extern int BAL_InitCellSequencePartial( typeCellSequence *l,
                                 typePropertyList *theProperties, typeUpdateList *theIndexes );


extern void _fprintfCell( FILE *f, typeCell *c );
extern void BAL_FprintfCellImage( FILE *f, typeCellImage *l );
extern void BAL_FprintfCellSequence( FILE *f, typeCellSequence *l );


typedef enum enumSurfaceEstimation {
    _4EDGES_,
    _OUTER_4NEIGHBORS_,
    _OUTER_6NEIGHBORS_,
    _LINDBLAD_,
    _WINDREICH_
} enumSurfaceEstimation;

extern void BAL_FprintfEnumSurfaceEstimation( FILE *f, enumSurfaceEstimation s );



extern int BAL_CellImageProperties( bal_image *theIm, typeCellImage *theList,
                                    typePropertyList *propertyList,
                                    enumSurfaceEstimation surfaceType );


extern int BAL_CellImageForwardIntersection( bal_image *theMother, bal_image *theDaughter, typeCellImage *theList );
extern int BAL_CellImageBackwardIntersection( typeCellImage *theMother, typeCellImage *theDaughter );
extern int BAL_CellSequenceBackwardIntersection( typeCellSequence *l );

extern int BAL_CellSequenceDiagnosis( FILE *f, typeCellSequence *cellProperties,
                                      typePropertyList *properties );



extern void XML_FprintfCellImageProperty( FILE *f, typeCellImage *l, typePropertyList *output );
extern void XML_FprintfCellSequenceProperty( FILE *f, typeCellSequence *l, typePropertyList *output );

extern int XML_FscanfCellSequenceProperty( FILE *f, typeCellSequence *l, typePropertyList *o );


#ifdef __cplusplus
}
#endif

#endif
