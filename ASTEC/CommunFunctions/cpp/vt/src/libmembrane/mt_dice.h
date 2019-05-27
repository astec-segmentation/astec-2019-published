/*************************************************************************
 * mt_dice.h -
 *
 * $Id: mt_dice.h,v 1.6 2013/07/30 11:30:00 gael Exp $
 *
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2013/07/30
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _mt_dice_h_
#define _mt_dice_h_


#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include <vt_image.h>

extern int MT_ComputeTab3D(vt_image imageIn,
                    vt_image imageExt,
                    int ***table,
                    int **inHist,
                    int **extHist,
                    int **indicesIn,
                    int **indicesExt,
                    int *sizeIn,
                    int *sizeExt
                    );


extern int MT_ComputeDice3D( vt_image imageIn,
                             vt_image imageExt,
                             char *filename);

extern int MT_ComputeDice3DTrsf( vt_image imageIn,
                                 vt_image imageExt,
                              double *T_ext_in,
                                 int **indicesInPtr,
                                int **indicesExtPtr,
                                int *nonEmptyInPtr,
                                int *nonEmptyExtPtr,
                                double ***DicePtr);


extern int MT_ComputeDice2D( vt_image imageIn,
                             vt_image imageExt,
                             char *filename );

extern int MT_ComputeConfusionTab3D( vt_image labelIn,
                             vt_image imageExt,
                             char *filename );

extern int MT_ComputeIntersection3D( vt_image imageIn,
                             vt_image imageExt,
                             char *filename);

extern int MT_ComputeSymmetryDice3D( vt_image *imageIn,
                             int x_pos,
                             char *filename);

extern int MT_ComputeSymmetryDice3DBis( vt_image *imageIn,
                             int x_pos,
                             double *r);

extern int MT_ComputeSymmetryDiceOblic3D( vt_image *imageIn,
                             double *n,
                             char *filename);

extern int MT_ComputeSymmetryDiceOblic3DBis( vt_image *imageIn,
                             double *n,
                             double *r);

extern int VT_ComputePairDices(vt_image *imageIn, vt_image *imageExt, double *T_ext_in, int *labelsIn, int *labelsExt, int npairs,  double *dices);


#ifdef __cplusplus
}
#endif

#endif /* _mt_dice_h_ */
