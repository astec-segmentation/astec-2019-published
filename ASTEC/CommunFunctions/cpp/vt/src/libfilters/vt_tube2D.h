/*************************************************************************
 * vt_tube2D.h -
 *
 * $Id: vt_tube2D.h,v 1.6 2006/05/16 09:33:34 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Mon Jul 10 16:46:12 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _vt_tube2D_h_
#define _vt_tube2D_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
  /* #include <malloc.h> */

#include <vt_image.h>
#include <vt_common.h>

#include <vt_tubeutils.h>


typedef enum {
  KRISSIAN,
  FRANGI
} methodType;


typedef struct vt_2Dtensor {
  vt_image imxx;
  vt_image imxy;
  vt_image imyy;

  /* valeurs propres
   */
  vt_image imvp1;
  vt_image imvp2;
  /* angles des vecteurs propres
     v = ( cos theta, sin theta )
  */
  vt_image imtheta1;
  vt_image imtheta2;  
} vt_2Dtensor;




typedef struct vt_2Dimages {
  /* filtrage
   */
  vt_image imx;
  vt_image imy;

  vt_2Dtensor hessien;
  /* reponse, extrema
   */
  vt_image imr;
  vt_image ime;

} vt_2Dimages;



typedef enum {
  _BLACK_,
  _WHITE_
} enumStructureColor;


typedef struct {
  vt_image imTheta;
  vt_image imPhi;
  vt_image imRep;
  vt_image imScale;
} vt_3Dimres;






extern int VT_2DTensorVoting( vt_2Dtensor *par,
		       vt_image *imWeight,
		       vt_image *imTheta,
		       vt_image *imMask,
		       double largeSigma,
		       double smallSigma,
		       double multSigma );

 
extern int VT_2DTensorGaussianVoting( vt_2Dtensor *par,
		       vt_image *imRes,
		       vt_image *imTheta,
		       vt_image *imMask,
		       double sigma );

extern int VT_Compute2DMultiScale( vt_image *theIm,
			    vt_3Dimres *imsRes,
			    double scale1,
			    double scale2,
			    int nbscales,
			    enumStructureColor color,
			    methodType mode );

extern void VT_Compute2DExtrema( vt_3Dimres *imsRes,
			  vt_image *imExt );

extern void VT_Compute2DMaskedExtrema( vt_image *imRes, 
			  vt_image *imTheta,
			  vt_image *imMask );

extern void VT_Compute2DResponse( vt_2Dimages *par, enumStructureColor color,
				  double theta, double sigma );

extern void VT_Compute2DResponseFrangi( vt_2Dimages *par, double sigma );

extern void VT_FilterOn2DEigenValues( vt_2Dimages *par, enumStructureColor color );

void VT_Compute2DEigenVectors( vt_2Dtensor *par,
			       vt_image *imMask );

extern int VT_Filter2Dimages( vt_image *im, vt_2Dimages *par, float *theCoeffs );

extern int VT_Filter2Dimages2ndOrder( vt_image *im, vt_2Dimages *par, float *theCoeffs );



extern int  VT_Alloc2Dtensor( vt_2Dtensor *par, vt_image *im, char *genericname );
extern void VT_Free2Dtensor ( vt_2Dtensor *par );
extern void VT_Write2Dtensor( vt_2Dtensor *par );

extern int  VT_Alloc2Dimages( vt_2Dimages *par, vt_image *im, char *genericname );
extern void VT_Free2Dimages ( vt_2Dimages *par );
extern void VT_Write2Dimages( vt_2Dimages *par );

extern int VT_Alloc3DImres( vt_3Dimres *par, vt_image *im, 
			    char *genericname, int alloc_phi );
extern void VT_Write3DImres( vt_3Dimres *par, int write_phi );
extern void VT_Write3DImresAngles( vt_3Dimres *par, int write_phi );
extern void VT_Free3DImres( vt_3Dimres *par );

#ifdef __cplusplus
}
#endif

#endif /* _vt_tube2D_h_ */

