/*************************************************************************
 * vt_tube2Dmatlab.h -
 *
 * $Id: vt_tube2Dmatlab.h,v 1.7 2002/04/19 15:41:39 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Tue Jul 11 11:39:48 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _vt_tube2Dmatlab_h_
#define _vt_tube2Dmatlab_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
  /* #include <malloc.h> */

#include <vt_image.h>
#include <vt_common.h>

#include <vt_tube2D.h>


extern void VT_2DDrawImage( vt_image *theIm,
			    int slice,
			    int rawDataFileDesc,
			    FILE *commandFile );


extern void VT_3DDrawImage( vt_image *theIm,
			    int rawDataFileDesc,
			    FILE *commandFile );


extern void VT_2DDrawWeightedVectors( vt_image *imWeight,
			       vt_image *imTheta,
			       int rawDataFileDesc,
			       FILE *commandFile );


extern void VT_Old3DDrawWeightedVectors( vt_image *imWeight,
			       vt_image *imTheta,
			       vt_image *imPhi,
			       int rawDataFileDesc,
			       FILE *commandFile );
extern void VT_3DDrawWeightedVectors( vt_image *imWeight,
			       vt_image *imTheta,
			       vt_image *imPhi,
			       int rawDataFileDesc,
			       FILE *commandFile, double threshold );

#ifdef __cplusplus
}
#endif

#endif /* _vt_tube2Dmatlab_h_ */

