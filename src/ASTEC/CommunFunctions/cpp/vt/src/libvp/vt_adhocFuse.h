/*************************************************************************
 * vt_adhocFuse.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Ven  2 fév 2018 16:40:20 CET
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _vt_adhocfuse_h_
#define _vt_adhocfuse_h_

#ifdef __cplusplus
extern "C" {
#endif


#include <vt_image.h>




typedef enum enumExtremaMethod {
  _GLOBAL_,
  _CELL_,
  _CELL_BORDER_,
  _CELL_INTERIOR_,
  _VOXEL_
} enumExtremaMethod;






extern int AdhocFuse_RescaleIntensityImage( vt_image *intImage,
                                            vt_image *minMskImage,
                                            vt_image *maxMskImage,
                                     vt_image *segImage,
                                     vt_image *resIntImage,
                                     vt_image *resMinImage,
                                     vt_image *resMaxImage,
                                     int *resMin,
                                     int *resMax,
                                     float percentileMin,
                                     float percentileMax,
                                     enumExtremaMethod methodMin,
                                     enumExtremaMethod methodMax,
                                     float sigma );

#ifdef __cplusplus
}
#endif

#endif
