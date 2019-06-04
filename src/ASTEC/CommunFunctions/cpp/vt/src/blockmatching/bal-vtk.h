/*************************************************************************
 * bal-vtk.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Lun 22 jan 2018 18:22:04 CET
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#ifndef BAL_VTK_H
#define BAL_VTK_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>

#include <bal-stddef.h>



extern void BAL_SetVerboseInBalVtk( int v );
extern void BAL_IncrementVerboseInBalVtk(  );
extern void BAL_DecrementVerboseInBalVtk(  );

extern void BAL_SetDebugInBalVtk( int v );
extern void BAL_IncrementDebugInBalVtk( );
extern void BAL_DecrementDebugInBalVtk( );



typedef enum enumVtkDataEncoding {
   _VTK_BINARY_DATA_,
   _VTK_ASCII_DATA_
} enumVtkDataEncoding;



typedef enum enumVtkAttributeType {
  _VTK_NO_ATTRIBUTE_,
  _VTK_CELL_ATTRIBUTE_,
  _VTK_POINT_ATTRIBUTE_
} enumVtkAttributeType;



typedef enum enumVtkData{
    _VTK_UNKNOWN_,
    _VTK_POINTS_,
    _VTK_LINES_,
    _VTK_POLYGONS_,
    _VTK_TRIANGLE_STRIPS_,
    _VTK_VERTICES_,
    _VTK_SCALARS_,
    _VTK_COLOR_SCALARS_,
    _VTK_LOOKUP_TABLE_,
    _VTK_VECTORS_,
    _VTK_NORMALS_,
    _VTK_TEXTURE_COORDINATES_,
    _VTK_TENSORS_,
    _VTK_FIELD_
} enumVtkData;



typedef struct bal_vtkData {
  enumVtkAttributeType attribute;
  enumVtkData type;
  int n_primitives;
  int n_data;
  void *data;
  char *name;
  int dataLength;
  int iodone;
} bal_vtkData;



typedef struct bal_vtkDataSet {
  int n_data;
  int n_allocated_data;
  bal_vtkData *data;

  enumUnitTransfo unit;
  typeVoxelSize vx;          /* real voxel size in X dimension */
  typeVoxelSize vy;          /* real voxel size in Y dimension */
  typeVoxelSize vz;          /* real voxel size in Z dimension */

} bal_vtkDataSet;



extern void BAL_InitVtkDataSet( bal_vtkDataSet *l );
extern void BAL_FreeVtkDataSet( bal_vtkDataSet *l );

extern void BAL_FprintfVtkDataSet( FILE *f, bal_vtkDataSet *s );

extern int BAL_ReadVtkDataSet( bal_vtkDataSet *m, char *name );
extern int BAL_WriteVtkDataSet( bal_vtkDataSet *dataSet, char *name, enumVtkDataEncoding dataEncoding );

extern void BAL_VtkDataSetFlipNormals( bal_vtkDataSet *l );
extern int BAL_CopyVtkDataSet( bal_vtkDataSet *theSet, bal_vtkDataSet *resSet );



#ifdef __cplusplus
}
#endif

#endif
