/*************************************************************************
 * bal-transformation-list-tools.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mar 21 jan 2014 18:01:34 CET
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#ifndef BAL_TRANSFORMATION_LIST_TOOLS_H
#define BAL_TRANSFORMATION_LIST_TOOLS_H

#ifdef __cplusplus
extern "C" {
#endif





#include <bal-image.h>
#include <bal-transformation.h>



extern void BAL_SetVerboseInBalTransformationListTools( int v );
extern void BAL_IncrementVerboseInBalTransformationListTools(  );
extern void BAL_DecrementVerboseInBalTransformationListTools(  );
extern void BAL_SetDebugInBalTransformationListTools( int v );
extern void BAL_IncrementDebugInBalTransformationListTools(  );
extern void BAL_DecrementDebugInBalTransformationListTools(  );



typedef struct typeBox {
  bal_floatPoint voxel;
  bal_floatPoint min;
  bal_floatPoint max;
  bal_floatPoint fov;
} typeBox;

typedef struct typeBoxList {
  typeBox *data;
  int n_data;
} typeBoxList;


extern void _initBoxList( typeBoxList *l );
extern void _freeBoxList( typeBoxList *l );
extern int _allocBoxList( typeBoxList *l, int n );

extern void _fprintfBoxList( FILE *f, typeBoxList *l, char *desc );


extern int _fillBoxesWithImageName( typeBoxList *theBox, bufferType *type, char *theName,
                                    bal_transformationList *theTrsf, int threshold,
                                    int refTrsf, bal_doublePoint *voxelSize );

extern int _fillBoxesWithImageNames( typeBoxList *theBox, bufferType *type,
                                     stringList *theName,
                                     bal_transformationList *theTrsf, int threshold,
                                     int refTrsf, int refTemplate  );



extern int BAL_ChangeTransformationList( bal_transformationList *theList,
                                         typeBoxList *boxList,
                                         bal_transformationList *resList,
                                         bal_image *resIm,
                                         int refTrsf,
                                         int *margin,
                                         int *extend );




#ifdef __cplusplus
}
#endif

#endif
