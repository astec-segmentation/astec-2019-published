/*************************************************************************
 * vt_tubeutils.h -
 *
 * $Id: vt_tubeutils.h,v 1.4 2002/09/05 17:15:06 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Tue Jul 11 15:46:40 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _vt_tubeutils_h_
#define _vt_tubeutils_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
  /* #include <malloc.h> */

#include <vt_image.h>
#include <vt_common.h>


extern double ** _ReadTransformationsFile( char *name, int *nb );



extern float _GetInterpolated2DValue( vt_image *theIm,
			       double x, double y, int z );

extern float _GetInterpolated3DValue( vt_image *theIm,
			       double x, double y, double z );






typedef struct {
  int dx, dy, dz;
  double c;
} typeConvolutionPoint;



extern typeConvolutionPoint *_Build2DConvolutionMask( int r,
					       double sigma,
					       int *nb );


#ifdef __cplusplus
}
#endif

#endif /* _vt_tubeutils_h_ */

