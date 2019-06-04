/*************************************************************************
 * mt_membrane3D.h -
 *
 * $Id: mt_membrane3D.h,v 2.0 2013/10/22 11:10:00 gael Exp $
 *
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 * 
 * CREATION DATE: 
 * 2013/06/17
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _mt_membrane3D_h_
#define _mt_membrane3D_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <mt_membrane2D.h>
#include <vt_tube3D.h>

#define NELEMS(n) (sizeof(n) / sizeof (*n))


extern void MT_SetVerboseInMtMembrane3D( int v );
extern void MT_IncrementVerboseInMtMembrane3D(  );
extern void MT_DecrementVerboseInMtMembrane3D(  );



/*
  import vt_tube3D.h :
  typedef - typeResponseInSlice
          - typeCirclePoint
          - vt_3Dimages
          - vt_2Dimauxs
          - typeCSplinesInSlice
          
  fonction- VT_BuildCircle
          - VT_Compute3DresponseInSlice
          - VT_Filter3Dimages
          - VT_Compute3DMultiScale
          - VT_Compute3DExtrema
          - VT_Filter2Dimauxs
          - VT_ComputeNextRealFrame
          - VT_ComputeRealFrame
          - VT_ComputeRealFrameWithPivot
          - VT_ReechTriLinSlice
          - VT_ReechTriLinSlice2
          - VT_ReechCSplineSlice
          - VT_ReechCSplineSlice2
*/


typedef enum enumMode {
  MODELBASED,
  ACME
} enumMode;


typedef struct {

  vt_image imxx;
  vt_image imyy;
  vt_image imzz;
  vt_image imxy;
  vt_image imxz;
  vt_image imyz;

  vt_image imvp1;
  vt_image imvp2;
  vt_image imvp3;

  vt_image imtheta1;
  vt_image imphi1;
  vt_image imtheta2;
  vt_image imphi2;
  vt_image imtheta3;
  vt_image imphi3;
  
  vt_image iszero;

} vt_3Dtensor;



extern int  VT_Alloc3DtensorFromImage( vt_3Dtensor *par, vt_image *im,
          char *genericname );
extern int  VT_Alloc3Dtensor( vt_3Dtensor *par, char *n, int x, int y, int z,
          int t );
extern void VT_Free3Dtensor ( vt_3Dtensor *par );
extern void VT_Write3Dtensor( vt_3Dtensor *par );
extern int VT_Write3DtensorWithNames( vt_3Dtensor *par,
                                      char *tensorname,
                                      char *eigenvaluename,
                                      char *anglename,
                                      char *binaryname,
                                      char *suffix );
extern int VT_Write3DtensorWithName( vt_3Dtensor *par, char *genericname);





typedef struct mt_angles {
  double **angles;
  int Nangles;
  double **vectors;
  double maxScalarProduct;
} mt_angles;

extern void MT_InitAngles( mt_angles *a );
extern void MT_FreeAngles( mt_angles *a );


extern int MT_Compute3DAngles( mt_angles *angles, int Niter );



extern int MT_Compute3DMultiScale( vt_image *theIm,
                            vt_3Dimres *imsRes,
                vt_image *mask,
                int flagMask,
                            double scale1,
                            double scale2,
                            int nbscales,
                            double zfact,
                                   enumStructureColor color,
                                   enumMode mode,
                                   int hsp);

extern void MT_Compute3DExtrema( vt_3Dimres *imsRes, vt_image *mask, int flagMask,
                                 double zfact,
                                 vt_image *imExt );



extern int MT_SampleBin(vt_image *imageBin, double sample);

extern void MT_SetAdditionalImagesForParallelVoting( int i );

extern int MT_Compute3DTensorVoting( vt_3Dtensor *theTensor,
	  vt_image **imBin,  double scale, double zfact, int Niter, int Nangles, int Nsticks,
	  enumTVmode mode, hessianMode initHessian, char *parName, int writeImages );




extern void MT_ComputeTensorSurfaceExtrema( vt_3Dtensor *imTensor,
    vt_image *imExt, double zfact  );

extern void MT_ComputeTensorBallExtrema( vt_3Dtensor *imTensor,
    vt_image *imExt );

extern void MT_ComputeTensorLineExtrema( vt_3Dtensor *imTensor,
    vt_image *imExt );





extern int MT_Compute3DTestFields( vt_3Dtensor **tfields, int *dimFields,
          mt_angles *angles, double r, double alpha);

                    
extern int MT_Compute3DStickFields( vt_3Dtensor **sfields, int *dimFields, 
          mt_angles *angles, double *theCoeffs, enumTVmode mode);


extern int MT_Compute3DBallFieldFromStick( vt_3Dtensor *bfield,
          vt_3Dtensor *sfields, int *dimFields,
          int Nangles); 

extern int MT_Compute3DPlateFields( vt_3Dtensor **pfields, int *dimFields,
          mt_angles *angles, double *theCoeffs, int Nsticks,
          enumTVmode mode);




#ifdef __cplusplus
}
#endif

#endif /* _mt_membrane3D_h_ */

