/*************************************************************************
 * vt_maskSeeds.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Jeu 26 jul 2018 10:21:24 CEST
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _vt_maskseeds_h_
#define _vt_maskseeds_h_

#ifdef __cplusplus
extern "C" {
#endif


#include <vt_image.h>


extern int VT_SelectSeedsinCells( vt_image *imseeds, vt_image *imcells, vt_image *res);



#ifdef __cplusplus
}
#endif

#endif
