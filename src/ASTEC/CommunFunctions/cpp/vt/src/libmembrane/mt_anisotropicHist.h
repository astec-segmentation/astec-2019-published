/*************************************************************************
 * mt_anisotropicHist.h -
 *
 * $Id: mt_anisotropicHist.h,v 1.0 2013/10/22 17:07:00 gael Exp $
 *
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2013/10/22
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _mt_anisotropicHist_h_
#define _mt_anisotropicHist_h_


#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include <vt_image.h>

typedef enum {
  _MSF_,
  _MT_,
  _UNKNOWN_
}  vt_typeSeuil;

typedef enum {
  _NO_TYPE_,
  _GAUSS_,
  _RAYLEIGH_,
  _RAYLEIGH_POS_,
  _RAYLEIGH_CENTERED_,
  _Xn_
}  vt_LM_type;

typedef enum {
  _HIST_PROJ_,
  _HIST_PROJ_2_,
  _HIST_N_,
  _HIST_N_2_
} vt_histmode;

typedef enum {
  _NO_BIN_,
  _BIN_PROJ_,
  _BIN_DOT_,
  _BIN_DOT2_,
  _BIN_AND_,
  _BIN_OR_,
  _BIN_ALL_
} vt_binmode;

typedef struct {
  vt_image *imrep;
  vt_image *imtheta;
  vt_image *imphi;
  vt_image *immask;
} vt_Rep_Angles;

typedef struct {
  int ncell;
  int initncell;
  double min;
  int initmin;
  double max;
  int initmax;
  double precision;
  int initprecision;
} vt_histparam;

typedef struct {
  int ncell;
  double lower;
  double upper;
  double precision;
  double total;
} vt_histDescriptor;

typedef struct {
  double *histX;
  double *histY;
  double *histZ;
  vt_histDescriptor descriptorX;
  vt_histDescriptor descriptorY;
  vt_histDescriptor descriptorZ;
} vt_anisotropicHisto;

typedef struct {
    double lmin;
    double lmax;
    int initlmin;
    int initlmax;
    double sigma_rayleigh_centered;
    double coef_rayleigh_centered;
} vt_levenbergConfig;

extern void MT_InitHistParam ( vt_histparam *histparam);

extern void MT_InitLevenParam ( vt_levenbergConfig *par);
extern void MT_PrintLevenParam ( vt_levenbergConfig par);


extern int MT_ComputeHisto3D( vt_Rep_Angles imageIn,
                             vt_anisotropicHisto *histo,
                             vt_histmode mode,
                             vt_histparam inithistparam,
                             int normalize,
                             double lmin, 
                             int initlmin, 
                             double lmax,
                             int initlmax, 
                             int automatic,
                             vt_levenbergConfig *par,
                             int flag_mask,
                             int wi,
                             vt_image *imResAnisotropes);

extern int MT_ComputeRepAnisotropic(vt_Rep_Angles imageIn,
                                    vt_image *imRepAnisotropes);

extern int MT_WriteHisto3D(vt_anisotropicHisto histo,
                            int cptLevenberg,
                            vt_LM_type *LM_types,
                            int nbLeven,
                            double *parLevenOutX,
                            double *parLevenOutY,
                            double *parLevenOutZ,
                            double **thresholds,
                            double *ratioIn,
                            int nratio,
                            double *thX,
                            double *thY,
                            double *thZ,
                            double *MSF,
                            int nMSF,
                            double *MT,
                            int nMT,
                            vt_binmode *binarise,
                            int nbin,
                            double **ratioOut,
                            char *filename,
                            char *prefix);

extern double MT_ComputeHistoN(double *buf,
                            double *bufn,
                            int dim,
                            vt_histmode mode,
                            double *hist,
                            vt_histDescriptor desc);

  /* extern int  MT_InitHistParam ( vt_histparam *histparam); */

extern int  MT_InitDescriptor ( vt_histDescriptor *desc,
                            double l,
                            double u,
                            int n,
                            vt_histparam initparam);

extern int  MT_InitDescriptorWithNCell ( vt_histDescriptor *desc,
                            double l,
                            double u,
                            int ncell);

extern int MT_ThresholdsFromHisto( double **thres,
                             vt_anisotropicHisto histo,
                             double *rates, double nrates);

extern int MT_Binarise(vt_Rep_Angles imIn,
                            vt_image *imResAnisotropes,
                            vt_image *imBin,
                            char *name,
                            double thX,
                            double thY,
                            double thZ,
                            vt_binmode bin,
                            double *r);

extern int MT_BinariseProgressive(vt_Rep_Angles imIn,
                            vt_image *imResAnisotropes,
                            vt_image *imBin,
                            char *name,
                            double thXbeg,
                            double thYbeg,
                            double thZbeg,
                            double thXend,
                            double thYend,
                            double thZend,
                            vt_binmode bin,
                            double *r);

extern void MT_Smooth(double *hist, double *smooth_hist, int n_elt);

extern int  VT_AllocRep_AnglesWithImageRep( vt_Rep_Angles *par);
extern void VT_InitRep_Angles ( vt_Rep_Angles *par );
extern void VT_FreeRep_Angles ( vt_Rep_Angles *par );

extern int  VT_AllocAnisotropicHisto ( vt_anisotropicHisto *par);
extern void VT_FreeAnisotropicHisto ( vt_anisotropicHisto *par );


#ifdef __cplusplus
}
#endif

#endif /* _mt_anisotropicHist_h_ */
