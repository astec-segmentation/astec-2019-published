/*************************************************************************
 * bal-transformation-cell.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Jeu 29 nov 2018 15:38:39 CET
 *
 * ADDITIONS, CHANGES
 *
 */

#ifndef _bal_transformation_cell_h_
#define _bal_transformation_cell_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <bal-image.h>
#include <bal-transformation.h>

extern void BAL_SetVerboseInBalTransformationCell( int v );
extern void BAL_IncrementVerboseInBalTransformationCell(  );
extern void BAL_DecrementVerboseInBalTransformationCell(  );

extern int BAL_ResampleCellImage( bal_image *image, bal_image *resImage,
                           bal_image *weightImage,
                           bal_transformation *theTr,
                                  float sigma );

#ifdef __cplusplus
}
#endif

#endif
