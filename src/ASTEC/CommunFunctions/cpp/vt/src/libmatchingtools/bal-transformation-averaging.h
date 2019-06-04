/*************************************************************************
 * bal-transformation-averaging.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Dim 16 jui 2013 16:56:47 CEST
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#ifndef BAL_TRANSFORMATION_AVERAGING_H
#define BAL_TRANSFORMATION_AVERAGING_H

#ifdef __cplusplus
extern "C" {
#endif






#include <bal-transformation.h>
#include <bal-estimator.h>


typedef enum {
  _MEAN_,
  _ROBUST_MEAN_
} enumTransformationAveraging;



extern void BAL_SetVerboseInBalTransformationAveraging( int v );
extern void BAL_IncrementVerboseInBalTransformationAveraging(  );
extern void BAL_DecrementVerboseInBalTransformationAveraging(  );

extern void BAL_SetDebugInBalTransformationAveraging( int d );
extern void BAL_IncrementDebugInBalTransformationAveraging(  );
extern void BAL_DecrementDebugInBalTransformationAveraging(  );

extern int BAL_ComputeAverageTransformation( bal_transformationList *listTrsf,
					     bal_transformation *resTrsf,
					     bal_estimator *estimator );


extern int BAL_EstimateTransformationByAveraging( bal_transformationArray *theTrsfs,
					       bal_transformationList *resTrsfs,
					       bal_estimator *estimator,
					       double sigma,
					       int n_iteration,
					       int useinverse );




#ifdef __cplusplus
}
#endif

#endif
