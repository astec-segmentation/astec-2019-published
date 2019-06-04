/*************************************************************************
 * vt_directionDistribution.h -
 *
 * $Id: vt_directionDistribution.h,v 1.0 2014/07/23 11:00:00 gael Exp $
 *
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2014/07/23
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _vt_directionDistribution_h_
#define _vt_directionDistribution_h_


#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <transfo.h>
#include <eigens.h>

#include <vt_image.h>

#define PI 3.141592654

typedef struct {
  int x;
  int y;
  int z;
  double n[3];
  double coef;
  double weight;
} vt_weighted_vector;


/*
extern int VT_totalLeastSquarePlane(vt_image* imageIn, 
									vt_image* imageTht, 
									vt_image* imagePhi, 
									double L, 
									double sigma,
									double *n, 
									double *E);
*/
  extern int VT_totalLeastSquarePlane(vt_weighted_vector* image, /* nuage de points */
									int nb,
									double L, 
									double sigma,
									double *n, 
									double *E);
extern double gaussianFun(double x, double sigma2x2);

extern int sphereToNormal(void*** array, bufferType t, int* dim, double* n );

extern void VT_FreeWeightedVector(vt_weighted_vector *v);

  extern int VT_AllocateAndBuildWeightedVector(	vt_image* imageIn,	/* nuage de points */
						vt_image* imageTht, 	/* angle spherique theta */
						vt_image* imagePhi, 	/* angle spherique phi */
						double L, 		/* demi-largeur du sous-espace */
						double *n, 		/* parametres du plan initial */
						vt_weighted_vector **image,	/* liste a allouer */
						int *nb	/* taille de la liste */
								);

extern double VT_MeasureSquareError(vt_weighted_vector *v, int nb, double *n, double L, double sigma_p);

extern int VT_AddPlaneToImage(vt_image* img, double *n, void * value);

extern int VT_ExtractMaxima(vt_image* histo, double** normalx, double** normaly, double** normalz, double **values, int *N);


#ifdef __cplusplus
}
#endif

#endif /* _vt_directionDistribution_h_ */
