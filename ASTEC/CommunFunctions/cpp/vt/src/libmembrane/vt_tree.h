/*************************************************************************
 * vt_tree.h -
 *
 * $Id: vt_tree.h,v 1.0 2015/07/27 15:15:00 gael Exp $
 *
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2015/07/27
 *
 * ADDITIONS, CHANGES
 *
 *
 */


#ifndef VT_TREE_H
#define VT_TREE_H


#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <vt_common.h>
#include <bal-estimator.h>



typedef struct {
    double x;
    double y;
    double z;
    double weight;
    int value;
    int attribute;
} vt_point;

typedef struct {
    vt_point *list;
    int n;
} vt_pointList;

typedef struct tree_node{
    int axis;
    vt_point median;
    struct tree_node *leftChild;
    struct tree_node *rightChild;
} tree_node;

extern int VT_NPoints(vt_image *img);

extern int VT_AllocPointListsWithImage(vt_pointList **p, int *n, vt_image *image);

extern int VT_AllocPointListWithImage(vt_pointList *p, vt_image *image) ;

extern int VT_AllocPointListWithFileName(vt_pointList *p, char *name);

extern int VT_AllocPointListsWithFileName(vt_pointList **p, int *nClasses, char *name);

extern int VT_AllocPointList(vt_pointList *p, int n);

extern int VT_ExtractPointListWithImage(vt_image img, vt_pointList* p);

extern int VT_ExtractPointListsWithImage(vt_image img, vt_pointList** p_list, int nClasses);

extern int VT_ExtractPointListWithFileName(char *name, vt_pointList* p);

extern int VT_ExtractPointListsWithFileName(char *name, vt_pointList** p, int nClasses);

extern void VT_FreePointList(vt_pointList *p);

extern int VT_TreeComputeCP(vt_pointList p_in, vt_pointList p_ext, tree_node *node_ext, double *T_ext_in_init, double *T_ext_in, double *residual, int flag_affine, bal_estimator estimator) ;

extern int VT_TreesComputeCP(vt_pointList p_in, vt_pointList *p_ext_list, tree_node **nodes_ext, int nClasses, double *T_ext_in_init, double *T_ext_in, double *residual, int flag_affine, bal_estimator estimator) ;

extern void tri_rapide(vt_pointList p, int axis, int premier, int dernier, int *ind);

extern void tri_rapide_vector(double *v, int premier, int dernier, int *ind);

extern int VT_TriRapideVector(double *v, int n, int **Ind);

extern void vt_kdtree(tree_node **node, vt_pointList p, int depth);

extern void vt_kdtrees(tree_node ***nodeList, vt_pointList *p, int n);

extern int beforeAfterSubsets(vt_pointList p, int axis, vt_pointList *pL, vt_pointList *pR, vt_point *median);

extern void printTree(tree_node *node);

extern void clearTree(tree_node **p);

extern int VT_ComputePointListsTrsf(int **pairs, int n, vt_pointList p_1, vt_pointList p_2,
                                    double *r, double *sd, double *T,
                                    int flag_trsf, bal_estimator estimator);

extern int VT_ComputePointListsListTrsf(int **pairs, int n, vt_pointList p_1, vt_pointList *p_list,
                                    double *r, double *sd, double *T,
                                    int flag_trsf, bal_estimator estimator);

extern void VT_PointListBarycentreInitTrsf(vt_pointList p, double *T, double *B);

extern void VT_TreeBarycentre(tree_node *node_ext, double *B_ext, double *weight);

extern int VT_TransformPointList(double *T, vt_pointList p, vt_pointList *p_trsf);

extern int VT_ComputePairsOfPointsWithTree(vt_pointList p_1, tree_node *node_2, int **pairs, int *npairs, double *dnorm);

extern int VT_ComputePairsOfPointsWithTrees(vt_pointList p_1, tree_node **node_list, int **pairs, int *npairs, int nClasses, double *dnorm);

extern int VT_NNS(vt_point q, tree_node *node, vt_point *best, double *d);

extern int VT_ResidualsVector(vt_pointList p_in, tree_node *node_ext, double *T_ext_in, double **residuals);

extern int VT_LabelsResidualsVector(vt_pointList p_in, tree_node **nodes_ext, double *T_ext_in, double **residuals);

extern void VT_SetResidualImage(vt_image *imResiduals, vt_pointList p_in, double *residuals);

#ifdef __cplusplus
}
#endif

#endif /*  VT_TREE_H */
