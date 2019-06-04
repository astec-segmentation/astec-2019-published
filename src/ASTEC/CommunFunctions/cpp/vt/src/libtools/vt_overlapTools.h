/*************************************************************************
 * vt_overlapTools.h - outils pour la comparaison de segmentations
 *
 * $Id: vt_overlapTools.h,v 1.1 2015/05/13 11:57:54 gael Exp $
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
 * May 13 2015
 *
 * Copyright Gael Michelin, INRIA
 *
 *
 * ADDITIONS, CHANGES:
 *
 *
 */

#ifndef _vt_overlap_tools_h_
#define _vt_overlap_tools_h_

#ifdef __cplusplus
extern "C" {
#endif


#include <vt_common.h>
#include <vt_neighborhood.h>


#define ORDERMAX 5


typedef enum {
    PROBABILITY,
    DICE,
    JACCARD
} measureType;

typedef enum {
    LIGHT,
    COMPLETE,
    MATRIX
} outputFormat;

typedef enum {
    MAX,
    NONNULL,
    THRESHOLD
} analysisMethod;

/* bbox : inputs x, y, z must be of type int *{x,y,z}[2];
 */
extern int bbox(vt_image *img, int **x, int **y, int **z, long long **volumes, int *nlab);

extern void setLabel(vt_image *imgIn, vt_image *imgOut, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int labelIn, int labelOut);

extern long long volume(vt_image *img, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int label);

extern int recouvrement(vt_image *imgIn, vt_image *imgExt,
                        int xminIn, int xmaxIn, int yminIn, int ymaxIn, int zminIn, int zmaxIn,
                        int xminExt, int xmaxExt, int yminExt, int ymaxExt, int zminExt, int zmaxExt,
                        int labelIn, int labelExt,
                        long long *volumeIn, long long *volumeExt,
                        double *recouvrementIn, double *recouvrementExt);

extern float overlapFun ( int ab, int a, int b, measureType measure );

extern int Neighborhood2Int ( Neighborhood N );

extern int VT_SegmentationOverlapping( vt_image *imA, vt_image *imB, vt_image *imAB, measureType measure, double r, int flag_real, int bckgrdA, int bckgrdB);

extern int VT_OverlapPruning( vt_image *im, double epsilon);

extern int VT_OverlapAnalysis(vt_image *imAB, vt_image *imL, analysisMethod method, double lt, double ht, unsigned short int *LAB, int *sizeA, int *sizeB);

extern void VT_WriteOverlapAnalysis(vt_image *imL, char **output, unsigned short int LAB, int sizeA, int sizeB, outputFormat format);

extern void VT_WriteOverlapInterpretationVoxel2(char **output, long long  *volume_bijA, long long  *volume_bijB, double *proportion_recouvrementA, double *proportion_recouvrementB,
                                          long long  bijectionA, long long  bijectionB, long long  volume_totalA, long long  volume_totalB);

extern void VT_WriteOverlapInterpretationCell(char **output, int nLabelsA, int nLabelsB, int bijectionA, int bijectionB,
                                              int *sous_segA, int *sous_segB, int *sur_segA, int *sur_segB, int Afond, int Bfond, int divergentA, int divergentB);

extern void VT_WriteOverlapInterpretationVoxel(char **output, long long  vox_bijectionA, long long  vox_bijectionB,
                                              long long  *vox_sous_segA, long long  *vox_sous_segB, long long  *vox_sur_segA, long long  *vox_sur_segB, long long  vox_Afond, long long  vox_Bfond, long long  vox_divergentA, long long  vox_divergentB);

extern int VT_OverlapInterpretation(FILE *fp, vt_image *imA, vt_image *imB, int flagIn, int flagExt, vt_image* imAout, vt_image *imBout, int flagInOut, int flagExtOut,
                                    int *nLabelsA, int *nLabelsB, int *Afond, int *Bfond, long long  *volume_totalA, long long  *volume_totalB,
                                    int *bijectionA, int *bijectionB, int *sous_segA, int *sous_segB, int *sur_segA, int *sur_segB,
                                    int *divergentA, int *divergentB,
                                    long long *_vox_Afond, long long *_vox_Bfond,
                                    long long *_vox_bijectionA, long long *_vox_bijectionB, long long *vox_sous_segA, long long *vox_sous_segB, long long *vox_sur_segA, long long *vox_sur_segB,
                                    long long *_vox_divergentA, long long *_vox_divergentB,
                                    long long **volume_bijA, long long **volume_bijB, double **proportion_recouvrementA, double **proportion_recouvrementB);

#ifdef __cplusplus
}
#endif

#endif /* _vt_overlap_tools_h_ */
