/*************************************************************************
 * vt_tube3D.h -
 *
 * $Id: vt_tube3D.h,v 1.9 2006/05/16 09:33:34 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Mon Sep 11 18:21:59 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _vt_tube3D_h_
#define _vt_tube3D_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>
#include <vt_tube2D.h>
#include <reech4x4.h>
#include <cspline.h>


extern int VT_GetVerboseInVtTube3D( );
extern void VT_SetVerboseInVtTube3D( int v );
extern void VT_IncrementVerboseInVtTube3D(  );
extern void VT_DecrementVerboseInVtTube3D(  );
extern int VT_GetDebugInVtTube3D( );
extern void VT_SetDebugInVtTube3D( int v );
extern void VT_IncrementDebugInVtTube3D(  );
extern void VT_DecrementDebugInVtTube3D(  );





typedef struct {
  double c2;
  double c3;
  double v;
  char   f;
} typeCirclePoint;


extern typeCirclePoint *VT_BuildCircle( double radius, int *n );









typedef struct {

  /* tableaux 2D
   */
  float **theXX;
  float **theYY;
  float **theZZ;
  float **theXY;
  float **theXZ;
  float **theYZ;

  /* tableaux 3D
   */
  float ***theX;
  float ***theY;
  float ***theZ;

  /* resultats 2D
   */
  float **theRep;
  float **theTheta;
  float **thePhi;

} typeResponseInSlice;



extern void VT_Compute3DresponseInSlice( typeResponseInSlice *aux,
					typeCirclePoint *thePts,
					int nbPts,
					 double radius,
					int dimx,
					int dimy,
					int dimz,
					int slice,
					enumStructureColor color,
					methodType mode );











typedef struct {
  /* le gradient
   */
  vt_image imx;
  vt_image imy;
  vt_image imz;
  /* les images filtrees en z
     qui ne sont pas filtrees selon x et y
     1. lissage selon Z
     2. derivee premiere selon Z
   */
  vt_image tmp0;
  vt_image tmp1;
  /* l'image derivee seconde selon z
   */
  vt_image imzz;
} vt_3Dimages;



extern int  VT_Alloc3Dimages( vt_3Dimages *par, vt_image *im, char *genericname );
extern void VT_Write3DImages( vt_3Dimages *par );
extern void VT_Free3DImages ( vt_3Dimages *par );
extern int  VT_Filter3Dimages( vt_image *im, vt_3Dimages *par, float *theCoeffs );
extern int  VT_Filter3Dimages2ndOrder( vt_image *im, vt_3Dimages *par, float *theCoeffs );












extern int VT_Compute3DMultiScale( vt_image *theIm,
			    vt_3Dimres *imsRes,
			    double scale1,
			    double scale2,
			    int nbscales,
				   enumStructureColor color,
				   methodType mode );

extern void VT_Compute3DExtrema( vt_3Dimres *imsRes,
				 vt_image *imExt );








typedef struct {
  /* les images construites
     a partir de l'image 3D lissee selon z
  */
  vt_image imxx;
  vt_image imyy;
  vt_image imxy;
  /* les images construites
     a partir de l'image 3D derivee selon z
  */
  vt_image imxz;
  vt_image imyz;

} vt_2Dimauxs;



extern int VT_Alloc2Dimauxs( vt_2Dimauxs *par, vt_image *im, char *genericname );
extern void VT_Write2Dimauxs( vt_2Dimauxs *par );
extern void VT_Free2Dimauxs( vt_2Dimauxs *par );
/* Calcul de derivees dans le plan
   - calcul de xx
   - calcul de yy
   - calcul de xy
   - calcul de xz
   - calcul de yz
 */
extern int VT_Filter2Dimauxs( vt_3Dimages *ims, vt_2Dimauxs *par, 
			      float *theCoeffs, int slice );









/* calcule le repere (dans le monde reel)
   permettant de reechantillonner une coupe orthogonale
   a partir de la direction du vaisseau (dans le volume image)
   
   La direction du vaisseau doit etre la troisieme direction
   (qui sera associee au Z de la coupe reechantillonnee)
 */


extern void VT_ComputeNextRealFrame( double thetaInImage,
			      double phiInImage,
			      double *ptSize,
			      double *newRealFrame,
			      double *oldRealFrame );

extern void VT_ComputeRealFrame( double *frame,
			  double *previousFrame );

extern void VT_ComputeRealFrameWithPivot( double *frame,
					  double *previousFrame, 
					  double tolerance );

extern int VT_ReechTriLinSlice( vt_image *theSlice, 
				int slice,
				vt_image *theIm,
				double *ptInImage,
			  double thetaInImage,
			  double phiInImage,
			  double *newRealFrame,
			  double *oldRealFrame,
			  double *newMat );

extern int VT_ReechTriLinSlice2( vt_image *theSlice,
			 int slice,
			 vt_image *theIm,
			 double *mat );


typedef struct {
  vt_image imxx;
  vt_image imxy;
  vt_image imxz;
  vt_image imyy;
  vt_image imyz;
  vt_image imzz;
  vt_image imx;
  vt_image imy;
  vt_image imz;
  vt_image d1;
  vt_image d2;
} typeCSplinesInSlice;

extern int VT_AllocCSplinesInSlice( typeCSplinesInSlice *par,
			     vt_image *theIm );


extern void VT_FreeCSplinesInSlice( typeCSplinesInSlice *par );




extern int VT_ReechCSplineSlice( vt_image *theSlice,
			  vt_image *theGrad,
			  vt_image *theLapl,
			  int slice,
			  typeCSplineCoefficients *theCoeff,
			  typeCSplinesInSlice *aux,
			  vt_image *theIm,
			  double *ptInImage,
			  double thetaInImage,
			  double phiInImage,
			  double *newRealFrame,
			  double *oldRealFrame,
			  double *newMat );
extern int VT_ReechCSplineSlice2( vt_image *theSlice,
			  int slice,
			  typeCSplineCoefficients *theCoeff,
			  typeCSplinesInSlice *aux,
			  double *mat );

#ifdef __cplusplus
}
#endif

#endif /* _vt_tube3D_h_ */

