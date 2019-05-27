/*************************************************************************
 * api-linearCombination.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2019, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Mer  6 mar 2019 11:17:18 CET
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _api_linearcombination_h_
#define _api_linearcombination_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <string-tools.h>
#include <typedefs.h>

extern int API_linearCombination( stringList *weightNames, stringList *imageNames, char *resName, ImageType type );

#ifdef __cplusplus
}
#endif

#endif
