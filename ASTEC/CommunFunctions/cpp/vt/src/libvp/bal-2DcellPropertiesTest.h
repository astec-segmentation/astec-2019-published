/*************************************************************************
 * bal_2DcellPropertiesTest.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Jeu 25 oct 2018 11:38:04 CEST
 *
 * ADDITIONS, CHANGES
 *
 */

#ifndef _bal_2dcellpropertiestest_h_
#define _bal_2dcellpropertiestest_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <bal-image.h>
#include <bal-cellProperties.h>



typedef enum enumShape {
    _SQUARE_,
    _DISK_
} enumShape;



extern int BAL_2DcellPropertiesShapeTest( int dlMin, int dlMax, int angleMin, int angleMax, int angleDelta,
                                   enumShape shape, enumSurfaceEstimation surface );



#ifdef __cplusplus
}
#endif

#endif
