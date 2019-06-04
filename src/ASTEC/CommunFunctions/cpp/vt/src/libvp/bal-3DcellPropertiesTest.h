/*************************************************************************
 * bal_3DcellPropertiesTest.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Lun 26 nov 2018 11:42:23 CET
 *
 * ADDITIONS, CHANGES
 *
 */

#ifndef _bal_3dcellpropertiestest_h_
#define _bal_3dcellpropertiestest_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <bal-cellProperties.h>
#include <bal-2DcellPropertiesTest.h>


extern int BAL_3DcellPropertiesShapeTest( typeCellSequence *theSequence,
                                     int dl,
                                     float sigma,
                                     enumShape shape,
                                     enumSurfaceEstimation surfaceEstimationType,
                                     char *basename );

extern int BAL_3DcellPropertiesImageTest( typeCellSequence *theSequence,
                                   bal_image *theIm,
                                   float minSigma, float maxSigma,
                              enumSurfaceEstimation surfaceEstimationType,
                              char *basename );

#ifdef __cplusplus
}
#endif

#endif
