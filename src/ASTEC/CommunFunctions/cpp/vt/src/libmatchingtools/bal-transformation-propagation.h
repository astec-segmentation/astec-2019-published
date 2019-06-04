/*************************************************************************
 * bal-transformation-propagation.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2015, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mer 28 oct 2015 22:19:40 CET
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#ifndef BAL_TRANSFORMATION_PROPAGATION_H
#define BAL_TRANSFORMATION_PROPAGATION_H

#ifdef __cplusplus
extern "C" {
#endif






#include <bal-transformation.h>


extern void BAL_SetVerboseInBalTransformationPropagation( int v );
extern void BAL_IncrementVerboseInBalTransformationPropagation(  );
extern void BAL_DecrementVerboseInBalTransformationPropagation(  );

extern void BAL_SetDebugInBalTransformationPropagation( int d );
extern void BAL_IncrementDebugInBalTransformationPropagation(  );
extern void BAL_DecrementDebugInBalTransformationPropagation(  );

extern int BAL_EstimateTransformationByPropagation( bal_transformationArray *theTrsfs,
                                             bal_transformationList *resTrsfs,
                                             int referenceindex );


#ifdef __cplusplus
}
#endif

#endif
