/*************************************************************************
 * vt_planeRegistration.h -
 *
 * $Id: vt_planeRegistration.h,v 1.0 2014/07/23 15:15:00 gael Exp $
 *
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2014/09/17
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef VT_PLANEREGISTRATION_H
#define VT_PLANEREGISTRATION_H


#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <transfo.h>
#include <eigens.h>
#include <math.h>
  /* #include <vt_tree.h> */

#include <vt_image.h>

#include <bal-estimator.h>

#define TABLENGTH 2000

typedef struct {
    double x;
    double y;
    double z;
    double weight;
    int label;
} vt_barycentre;

typedef struct {
    vt_barycentre *barycentres;
    int n;
    double vx;
    double vy;
    double vz;
} vt_image_barycentres;

extern int VT_ComputeLabelBarycentres(vt_image labels, vt_image_barycentres* b, int background);

extern int VT_NCell(vt_image img, int background);

extern int VT_FindLabel(vt_image_barycentres b, int val);

  /* extern int VT_InitImageBarycentres(vt_image_barycentres* b, vt_image image); // definie en statique */

extern int VT_TransformImageBarycentres(double *T, vt_image_barycentres b, vt_image_barycentres *b_trsf);

extern int VT_CopyImageBarycentres(vt_image_barycentres b, vt_image_barycentres *b_copy);

extern int VT_SetImageBarycentresWithList(vt_image_barycentres* b, char* filename, int background);

extern int VT_AllocImageBarycentres(vt_image_barycentres* b, int n);

extern int VT_AllocImageBarycentresWithImage(vt_image_barycentres* b, vt_image labels, int background);

extern void VT_FreeImageBarycentres(vt_image_barycentres* b);

extern void VT_PrintImageBarycentres(vt_image_barycentres b);

extern void VT_ImageBarycentresCentralLayersExtraction(vt_image_barycentres* b, double* n, double delta);

extern int VT_Barycentre(vt_image* gradient, double* barycentre);

extern void VT_BaryImgLabels(vt_image_barycentres b, double* G);

extern void VT_RotationMatrix(double *n_new, double* n, double* R);

extern void VT_RotationMatrixWithAngle(double *u, double theta, double* R);

extern void VT_RotationMatrix180(double *n, double *R);

extern void VT_Translation(double *G_new, double *G, double *R, double *t);

extern int VT_ComputePairsOfPoints(vt_image_barycentres b_1, vt_image_barycentres b_2, int **pairs, int *npairs, double *dnorm);

extern int VT_ComputePointTrsf(int **pairs, int n, vt_image_barycentres b_in, vt_image_barycentres b_ext,
                               double *r, double *sd, double *T_tmp, int flag_trsf, bal_estimator estimator);


extern int VT_AllocImagePointsWithImage(vt_image_barycentres* b, vt_image img);

extern int VT_ExtractImagePoints(vt_image img, vt_image_barycentres* b);

extern int VT_ComputeCP(vt_image_barycentres b_in, vt_image_barycentres b_ext, double *T_ext_in_init, double **T_ext_in, double *residual, int flag_affine, bal_estimator estimator) ;

#ifdef __cplusplus
}
#endif

#endif /*  VT_PLANEREGISTRATION_H */
