/*************************************************************************
 * vt_bound.h - extraction de parametres de box sur des parties numerotees
 *
 * $Id: vt_bound.h,v 1.10 2013/08/06 17:50:54 gael Exp $
 *
 * DESCRIPTION:
 *
 *
 *
 *
 *
 * AUTHOR:
 *
 *
 * CREATION DATE:
 * Aug 6 2013
 *
 * Copyright Gregoire Malandain, INRIA
 *
 *
 * ADDITIONS, CHANGES:
 *
 *
 */

#ifndef _vt_bound_h_
#define _vt_bound_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <vt_common.h>

extern void _VerboseInBound();
extern void _NoVerboseInBound();



typedef struct {

  int ptmin[3];
  int ptmax[3];

  int volume;

} typeBoxParameters;






extern int _CreateArrayOfParametersFromImage( vt_image *theIm,
                                       int slice,
                                       typeBoxParameters **thePar );





extern int _MaxValueInImage( vt_image *theIm,
                      int slice );
extern void _InitArrayOfParametersFromImage( vt_image *theIm,
                                     typeBoxParameters *thePar,
                                     int n );
extern void _FillArrayOfParametersFromImage( vt_image *theIm,
                                      int slice,
                                      typeBoxParameters *thePar );

#ifdef __cplusplus
}
#endif

#endif /* _vt_bound_h_ */
