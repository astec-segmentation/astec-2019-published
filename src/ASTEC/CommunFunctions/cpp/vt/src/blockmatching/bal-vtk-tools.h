/*************************************************************************
 * bal-vtk-tools.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Lun 29 jan 2018 10:26:27 CET
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#ifndef BAL_VTK_TOOLS_H
#define BAL_VTK_TOOLS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>

#include <bal-stddef.h>
#include <bal-vtk.h>
#include <bal-transformation.h>
#include <bal-transformation-tools.h>


extern void BAL_SetVerboseInBalVtkTools( int v );
extern void BAL_IncrementVerboseInBalVtkTools(  );
extern void BAL_DecrementVerboseInBalVtkTools(  );

extern void BAL_SetDebugInBalVtkTools( int v );
extern void BAL_IncrementDebugInBalVtkTools(  );
extern void BAL_DecrementDebugInBalVtkTools(  );

extern int BAL_TransformVtk( bal_vtkDataSet *theSet, bal_vtkDataSet *resSet, bal_transformation *theTrsf );



#ifdef __cplusplus
}
#endif

#endif
