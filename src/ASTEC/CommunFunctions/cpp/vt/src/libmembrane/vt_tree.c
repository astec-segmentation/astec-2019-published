/*
 * vt_tree.c -
 *
 * $Id: vt_tree.c,v 1.0 2015/07/27 10:52:34 gael Exp $
 *
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2015/07/23
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <vtmalloc.h>

#include <bal-field.h>
#include <bal-field-tools.h>
#include <bal-point.h>
#include <bal-transformation.h>
#include <bal-transformation-tools.h>

#include <vt_tree.h>





int VT_NPoints(vt_image *img)
{
    char *proc="VT_NPoints";

    unsigned char *bufu8=NULL;
    unsigned short int *bufu16=NULL;

    size_t i;
    int n=0;

    switch (img->type) {
    default:
            fprintf(stderr, "%s: such image type not handled yet\n", proc);
        return(-1);
    case UCHAR:
    case SCHAR:
        bufu8=(unsigned char *)img->buf;
        for (i=0 ; i<img->dim.x*img->dim.y*img->dim.z ; i++)
            if(bufu8[i] != (unsigned char) 0)
                n++;
        break;
    case USHORT:
    case SSHORT:
        bufu16=(unsigned short int *)img->buf;
        for (i=0 ; i<img->dim.x*img->dim.y*img->dim.z ; i++)
            if(bufu16[i] != (unsigned char) 0)
                n++;
        break;
    }

    return(n);
}

static int VT_NClasses(int **histo, int *N, vt_image *img)
{
    char *proc="VT_NClasses";

    unsigned char *bufu8=NULL;
    unsigned short int *bufu16=NULL;

    size_t i;
    int n=0;
    int *h = NULL;
    int v = 0;

    switch (img->type) {
    default:
        fprintf(stderr, "%s: such image type not handled yet\n", proc);
        return(-1);
    case UCHAR:
    case SCHAR:
        bufu8=(unsigned char *)img->buf;
        break;
    case USHORT:
    case SSHORT:
        bufu16=(unsigned short int *)img->buf;
        break;
    }
    for (i=0 ; i<img->dim.x*img->dim.y*img->dim.z ; i++) {
        switch (img->type) {
        default:
            fprintf(stderr, "%s: such image type not handled yet\n", proc);
            return(-1);
        case UCHAR:
        case SCHAR:
            v = (int) bufu8[i];
            break;
        case USHORT:
        case SSHORT:
            v = (int) bufu16[i];
            break;
        }
        if (n<v) n = v;
    }
    n++;
    h = vtmalloc(n * sizeof(int), "h", proc );
    for (i=0 ; i<(size_t)n ; i++) h[i] = 0;
    for (i=0 ; i<img->dim.x*img->dim.y*img->dim.z ; i++) {
        switch (img->type) {
        default:
            fprintf(stderr, "%s: such image type not handled yet\n", proc);
            return(-1);
        case UCHAR:
        case SCHAR:
            v = (int) bufu8[i];
            break;
        case USHORT:
        case SSHORT:
            v = (int) bufu16[i];
            break;
        }
        if (v == 0) continue;
        h[v] = h[v] + 1;
    }

    *N = n;
    *histo = h;
    /* fprintf(stdout, "%s: nClasses = %d\n", proc, n); */
    /* for (i=0;i<n;i++)    { */
    /*     fprintf(stdout, "h[%d] = %d\t", i, h[i]); */
    /* } */
    /* fprintf(stdout, "\n"); */
    return(1);
}

static int VT_NClassesFile(int **histo, int *N, char *name)
{
    char *proc="VT_NClassesFile";

    int i;

    int *h = NULL;

    FILE *f;
    char line[512];
    float f1, f2, f3, f4;
    int o;
    int n = 0;


    if ((f = fopen (name, "r")) == NULL) {
      fprintf( stderr, "%s: unable to open '%s' for reading\n", proc, name );
      return( -1 );
    }


    while ( fgets(line, 512, f) != NULL ) {
        o = sscanf( line, "%f %f %f %f\n", &f1, &f2, &f3, &f4);
        if ( o == 4 ) {
            n = ((int)f4>n) ? (int)f4 : n;
        }
    }
    n++;

    fclose( f );

    if ( n<0 ) {
      fprintf( stderr, "%s: negative number of points detected in '%s'\n", proc, name );
      return( -1 );
    }

    h = vtmalloc(n * sizeof(int), "h", proc );

    for (i=0 ; i<n ; i++) h[i] = 0;


    if ((f = fopen (name, "r")) == NULL) {
      fprintf( stderr, "%s: unable to open '%s' for reading\n", proc, name );
      return( -1 );
    }


    while ( fgets(line, 512, f) != NULL ) {
        o = sscanf( line, "%f %f %f %f\n", &f1, &f2, &f3, &f4);
        if ( o == 4 ) {
            if ((int)f4 < 0) {
                fprintf(stderr, "%s: unexpected negative class value in file %s, line is %s\n", proc, name, line);
                vtfree(h);
                fclose( f );
                return(0);
            }
            if ((int)f4 == 0) continue;
            h[(int)f4] = h[(int)f4] + 1;
        }
    }

    fclose( f );

    *N = n;
    *histo = h;
    /* fprintf(stdout, "%s: nClasses = %d\n", proc, n); */
    /* for (i=0;i<n;i++)    { */
    /*     fprintf(stdout, "h[%d] = %d\t", i, h[i]); */
    /* } */
    /* fprintf(stdout, "\n"); */
    return(1);
}


static int VT_InitPointListWithImage(vt_pointList* p, vt_image image)
{
    char *proc="VT_InitPointListWithImage";

    unsigned char *bufu8=NULL;
    unsigned short int *bufu16=NULL;

    size_t i;
    int val;
    int n=0;

    switch (image.type) {
    default:
            fprintf(stderr, "%s: such image type not handled yet\n", proc);
        return(-1);
    case UCHAR:
    case SCHAR:
        bufu8=(unsigned char *)image.buf;
        break;
    case USHORT:
    case SSHORT:
        bufu16=(unsigned short int *)image.buf;
        break;
    }

    for (i=0 ; i<image.dim.x*image.dim.y*image.dim.z ; i++) {
        switch (image.type){
        default:
                fprintf(stderr, "%s: such image type not handled yet\n", proc);
            return(-1);
        case UCHAR:
        case SCHAR:
            val=(int)bufu8[i];
            break;
        case USHORT:
        case SSHORT:
            val=(int)bufu16[i];
            break;
        }
        if (val==0) continue;
        n++;
    }

    if(n != p->n)
    {
            fprintf(stderr, "%s: incompatible vt_pointList field n = %d and vt_image number of points = %d\n", proc, p->n, n);
        return(0);
    }

    for (i=0; i<(size_t)p->n; i++) {
        vt_point *thePoint = &(p->list[i]);
        thePoint->attribute=0;
        thePoint->value=0;
        thePoint->weight=0;
        thePoint->x=0;
        thePoint->y=0;
        thePoint->z=0;
    }

    return (1);
}

int VT_AllocPointList(vt_pointList *p, int n)
{
    char *proc="VT_AllocPointList";
    if (n<0) {
        fprintf(stderr, "%s: negative number of points\n", proc);
        return(0);
    }
    if (n>0) {
      p->list = vtmalloc(n*sizeof(vt_point), "p->list", proc );
      if (p->list == (vt_point*)NULL) {
        fprintf(stderr, "%s: allocation failed\n", proc);
        return(0);
      }
    }
    else {
        p->list = (vt_point *)NULL;
    }
    p->n = n;
    return(1);

}

int VT_AllocPointListsWithImage(vt_pointList **p, int *n, vt_image *image)
{
    char *proc="VT_AllocPointListsWithImage";
    int i, N=0;
    int *histo=NULL;
    if (VT_NClasses(&histo, &N, image) != 1) {
        fprintf(stderr, "%s: unable to compute number of classes in the image\n", proc);
        return(0);
    }

    if(N<0) {
        fprintf(stderr, "%s: negative number of classes\n", proc);
        return(0);
    }
    if(N==0) {
        *p = NULL;
        *n = N;
        fprintf(stderr, "%s: zero classes found\n", proc);
        return(1);
    }
    vt_pointList *P = vtmalloc(N * sizeof(vt_pointList), "P", proc );

    /* fprintf(stdout, "%s : N = %d\n", proc, N); */

    for (i = 0; i < N ; i++) {
      /* fprintf(stdout, "%s : histo[%d] = %d\n", proc, i, histo[i]); */
        VT_AllocPointList(&(P[i]), histo[i]);
    }

    *p = P;
    *n=N;
    vtfree(histo);

    return(1);
}

int VT_AllocPointListWithImage(vt_pointList *p, vt_image *image)
{

    char *proc="VT_AllocPointListWithImage";

    int n=VT_NPoints(image);

    if (n<0) {
        fprintf(stderr, "%s: negative number of points\n", proc);
        return(0);
    }
    if (n==0) {
        p->n = 0;
        p->list = (vt_point*)NULL;
        return(1);
    }
    p->list = vtmalloc(n*sizeof(vt_point), "p->list", proc );
    if (p->list == (vt_point*)NULL) {
        fprintf(stderr, "%s: allocation failed\n", proc);
        return(0);
    }
    p->n = n;
    return(VT_InitPointListWithImage(p, *image));

}

int VT_AllocPointListsWithFileName(vt_pointList **p, int *n, char *name)
{
    char *proc="VT_AllocPointListsWithFileName";

    int i, N=0;
    int *histo=NULL;

    if (VT_NClassesFile(&histo, &N, name) != 1) {
        fprintf(stderr, "%s: unable to compute number of classes in the image\n", proc);
        return(0);
    }

    if(N<0) {
        fprintf(stderr, "%s: negative number of classes\n", proc);
        return(0);
    }
    if(N==0) {
        *p = NULL;
        *n = N;
        fprintf(stderr, "%s: zero classes found\n", proc);
        return(1);
    }

    vt_pointList *P = vtmalloc(N * sizeof(vt_pointList), "P", proc );

    for (i = 0; i < N ; i++) {
        VT_AllocPointList(&(P[i]), histo[i]);
    }

    *p = P;
    *n=N;
    vtfree(histo);

    return(1);

}

int VT_AllocPointListWithFileName(vt_pointList *p, char *name)
{
    char *proc="VT_AllocPointListWithFileName";

    FILE *f;
    char line[512];
    float f1, f2, f3, f4;
    int o;
    int n = 0;
    int i;

    if ((f = fopen (name, "r")) == NULL) {
      fprintf( stderr, "%s: unable to open '%s' for reading\n", proc, name );
      return( -1 );
    }


    while ( fgets(line, 512, f) != NULL ) {
        o = sscanf( line, "%f %f %f %f\n", &f1, &f2, &f3, &f4);
        if ( ( o == 3 || o == 4 ) ) {
            n+=1;
        }
    }

    fclose( f );

    if ( n<0 ) {
      fprintf( stderr, "%s: negative number of points detected in '%s'\n", proc, name );
      return( -1 );
    }

    if (n==0) {
        p->n = 0;
        p->list = (vt_point*)NULL;
        return(1);
    }

    p->list = vtmalloc(n*sizeof(vt_point), "p->list", proc );
    if (p->list == (vt_point*)NULL) {
        fprintf(stderr, "%s: allocation failed\n", proc);
        return(0);
    }
    p->n = n;

    for (i=0; i<p->n; i++) {
        vt_point *thePoint = &(p->list[i]);
        thePoint->attribute=0;
        thePoint->value=0;
        thePoint->weight=0;
        thePoint->x=0;
        thePoint->y=0;
        thePoint->z=0;
    }

    return(1);
}



int VT_ExtractPointListWithImage(vt_image img, vt_pointList* p)
{
    char *proc="VT_ExtractPointListWithImage";

    unsigned char ***arrayu8=NULL;
    unsigned short int ***arrayu16=NULL;

    int i, j, k, val;
    int ind;

    vt_point *thePoint;

    switch (img.type) {
    default:
            fprintf(stderr, "%s: such image type not handled yet\n", proc);
        return(-1);
    case UCHAR:
    case SCHAR:
        arrayu8=(unsigned char ***)img.array;
        break;
    case USHORT:
    case SSHORT:
        arrayu16=(unsigned short int ***)img.array;
        break;
    }

    ind=0;
    for (k=0 ; k<(int)img.dim.z ; k++)
    for (j=0 ; j<(int)img.dim.y ; j++)
    for (i=0 ; i<(int)img.dim.x ; i++) {
        switch (img.type){
        default:
                fprintf(stderr, "%s: such image type not handled yet\n", proc);
            return(-1);
        case UCHAR:
        case SCHAR:
            val=(int)arrayu8[k][j][i];
            break;
        case USHORT:
        case SSHORT:
            val=(int)arrayu16[k][j][i];
            break;
        }
        if(val==0) continue;

        thePoint=&(p->list[ind]);

        thePoint->weight = 1;
        thePoint->attribute = ind;
        thePoint->value = val;
        thePoint->x = i*img.siz.x;
        thePoint->y = j*img.siz.y;
        thePoint->z = k*img.siz.z;
        ind++;
    }

    fprintf(stderr, "%s : ind = %d \t p->n = %d\n", proc, ind, p->n);

    return(1);
}

int VT_ExtractPointListsWithImage(vt_image img, vt_pointList** P_list, int nClasses)
{
    char *proc="VT_ExtractPointListsWithImage";

    unsigned char ***arrayu8=NULL;
    unsigned short int ***arrayu16=NULL;
    vt_pointList* p;
    vt_pointList* p_list = *P_list;
    int i, j, k, val;
    int *ind;

    ind = vtmalloc (nClasses * sizeof(int), "ind", proc );
    if (ind==NULL)
    {
        fprintf(stderr, "%s: unable to allocate a vector\n", proc);
        return(-1);
    }
    for (i = 0 ; i < nClasses ; i++)
        ind[i]=0;

    vt_point *thePoint;

    switch (img.type) {
    default:
            fprintf(stderr, "%s: such image type not handled yet\n", proc);
        return(-1);
    case UCHAR:
    case SCHAR:
        arrayu8=(unsigned char ***)img.array;
        break;
    case USHORT:
    case SSHORT:
        arrayu16=(unsigned short int ***)img.array;
        break;
    }

    /* fprintf(stderr, "%s: test\n", proc); */

    /* ind=0; */
    for (k=0 ; k<(int)img.dim.z ; k++)
    for (j=0 ; j<(int)img.dim.y ; j++)
    for (i=0 ; i<(int)img.dim.x ; i++) {
        switch (img.type){
        default:
                fprintf(stderr, "%s: such image type not handled yet\n", proc);
            return(-1);
        case UCHAR:
        case SCHAR:
            val=(int)arrayu8[k][j][i];
            break;
        case USHORT:
        case SSHORT:
            val=(int)arrayu16[k][j][i];
            break;
        }
        if(val==0) continue;

        /* if(ind[val]<20) { */
        /*     fprintf(stderr, "%s: val = %d\tind[%d] = %d\n", proc, val, val, ind[val]); */
        /* } */

        p = &(p_list[val]);

        /* if(ind[val]<20) { */
        /*     fprintf(stderr, "%s: p_list[%d].n = %d\n", proc, val, p->n ); */
        /* } */
        thePoint=&(p->list[ind[val]]);

        thePoint->weight = 1;
        thePoint->attribute = ind[val];
        thePoint->value = val;
        thePoint->x = i*img.siz.x;
        thePoint->y = j*img.siz.y;
        thePoint->z = k*img.siz.z;
        ind[val] = ind[val]+1;
    }


    fprintf(stderr, "%s : \n", proc);
    for (i=0;i<nClasses;i++) {
        p = &(p_list[i]);
        if (p->n == 0 && ind[i] == 0) continue;
        fprintf(stderr, "     class %d : \t ind = %d \t p->n = %d\n", i, ind[i], p->n);
    }

    vtfree(ind);

    return(1);
}


int VT_ExtractPointListWithFileName(char *name, vt_pointList* p)
{
    char *proc="VT_ExtractPointListWithFileName";
    int ind;
    int o;
    vt_point *thePoint;
    ind=0;

    FILE *f;
    char line[512];
    float f1, f2, f3, f4;


    if ((f = fopen (name, "r")) == NULL) {
      fprintf( stderr, "%s: unable to open '%s' for reading\n", proc, name );
      return( -1 );
    }


    while ( fgets(line, 512, f) != NULL ) {
        o = sscanf( line, "%f %f %f %f\n", &f1, &f2, &f3, &f4);
        if ( o  == 3 || o == 4 ) {
            thePoint=&(p->list[ind]);

            thePoint->weight = 1;
            thePoint->value = (o == 4) ? f4 : 255;
            thePoint->attribute = ind;
            thePoint->x = f1;
            thePoint->y = f2;
            thePoint->z = f3;
            ind++;
        }
    }

    fclose( f );

    fprintf(stderr, "%s : ind = %d \t p->n = %d\n", proc, ind, p->n);

    if (ind != p->n) {
        fprintf(stderr, "%s: abnormal number of points found (%d instead of the %d expected)\n", proc, ind, p->n);
        return(0);
    }

    /*for (ind = 0 ; ind < p->n ; ind++)
    {
        thePoint = &(p->list[ind]);
        fprintf(stdout, "Point[%d].{x,y,z} = { %f, %f, %f }\n", ind, thePoint->x, thePoint->y, thePoint->z);
        fprintf(stdout, "         .value} = %d\n", thePoint->value);
        fprintf(stdout, "         .attribute} = %d\n", thePoint->attribute);
        fprintf(stdout, "         .weight} = %f\n", thePoint->weight);
    }*/

    return(1);
}


int VT_ExtractPointListsWithFileName(char *name, vt_pointList** P_list, int nClasses)
{
    char *proc="VT_ExtractPointListsWithFileName";

    int o;
    vt_point *thePoint;

    FILE *f;
    char line[512];
    float f1, f2, f3, f4;


    vt_pointList* p;
    vt_pointList* p_list = *P_list;
    int i;
    int *ind;

    ind = vtmalloc (nClasses * sizeof(int), "ind", proc );
    if (ind==NULL)
    {
        fprintf(stderr, "%s: unable to allocate a vector\n", proc);
        return(-1);
    }
    for (i = 0 ; i < nClasses ; i++)
        ind[i]=0;

    if ((f = fopen (name, "r")) == NULL) {
      fprintf( stderr, "%s: unable to open '%s' for reading\n", proc, name );
      return( -1 );
    }


    while ( fgets(line, 512, f) != NULL ) {
        o = sscanf( line, "%f %f %f %f\n", &f1, &f2, &f3, &f4);
        if ( o == 4 ) {

            if((int)f4<0) {
                fprintf(stderr, "%s: unexpected negative class value in %s, line %s\n", proc, name, line);
                vtfree(ind);
                fclose( f );
                return(0);
            }
            if((int)f4==0) continue;

            p = &(p_list[(int)f4]);

            thePoint=&(p->list[ind[(int)f4]]);

            thePoint->weight = 1;
            thePoint->value = (o == 4) ? (int)f4 : 255;
            thePoint->attribute = ind[(int)f4];
            thePoint->x = f1;
            thePoint->y = f2;
            thePoint->z = f3;
            ind[(int)f4] = ind[(int)f4]+1;
        }
    }

    fclose( f );

    /* fprintf(stderr, "%s: test\n", proc); */


    fprintf(stderr, "%s : \n", proc);
    for (i=0;i<nClasses;i++) {
        p = &(p_list[i]);
        if (p->n == 0 && ind[i] == 0) continue;
        fprintf(stderr, "     class %d : \t ind = %d \t p->n = %d\n", i, ind[i], p->n);
    }

    vtfree(ind);
    return(1);
}


void VT_FreePointList(vt_pointList *p)
{
    if (p->list)
      vtfree(p->list);
    p=(vt_pointList*) NULL;
}



int VT_TreeComputeCP(vt_pointList p_in, vt_pointList p_ext, tree_node *node_ext, double *T_ext_in_init, double *T_ext_in, double *residual, int flag_affine, bal_estimator estimator)
{


  char *proc="VT_TreeComputeCP";

  int *pairs[2];
  int npairs;

  double sd;
  double r;

  vt_pointList p_in_trsf;

  /* Allocation de la liste des points transformes */

  if (VT_AllocPointList(&p_in_trsf, p_in.n) != 1) {
    fprintf(stderr, "%s: unable to allocate point list p_in_trsf...\n", proc);
    return(0);
  }

  /* Transformation initiale des points de p_in */

  if (VT_TransformPointList(T_ext_in_init, p_in, &p_in_trsf) != 1) {
      VT_FreePointList(&p_in_trsf);
      fprintf(stderr, "%s: unable to apply initial transformation of p_in...\n", proc);
      return(0);
  }

  /* Calcul des appariements de points */

  if (VT_ComputePairsOfPointsWithTree(p_in_trsf, node_ext, pairs, &npairs, &r) != 1)
  {
      VT_FreePointList(&p_in_trsf);
      fprintf(stderr, "%s: unable to compute the pairings...\n", proc);
      return(0);
  }


  /* Calcul de la nouvelle transformation */

  if( VT_ComputePointListsTrsf(pairs, npairs, p_in, p_ext, &r, &sd, T_ext_in, flag_affine, estimator) != 1) {
      VT_FreePointList(&p_in_trsf);
      vtfree(pairs[0]); pairs[0]=NULL;
      vtfree(pairs[1]); pairs[1]=NULL;
      fprintf(stderr, "%s: unable to compute new transformation...\n", proc);
      return(0);
  }



  *residual = r;

  /* Frees */
  VT_FreePointList(&p_in_trsf);
  vtfree(pairs[0]); pairs[0]=NULL;
  vtfree(pairs[1]); pairs[1]=NULL;

  return(1);
}

int VT_TreesComputeCP(vt_pointList p_in, vt_pointList *p_ext_list, tree_node **nodes_ext, int nClasses, double *T_ext_in_init, double *T_ext_in, double *residual, int flag_affine, bal_estimator estimator)
{

    char *proc="VT_TreesComputeCP";


    double sd;
    double r;

    /* int i; */

    int *pairs[2];

    int npairs;

    vt_pointList p_in_trsf;

    /* Allocation de la liste des points transformes */

    if (VT_AllocPointList(&p_in_trsf, p_in.n) != 1) {
      fprintf(stderr, "%s: unable to allocate point list p_in_trsf...\n", proc);
      return(0);
    }

    /* Transformation initiale des points de p_in */

    if (VT_TransformPointList(T_ext_in_init, p_in, &p_in_trsf) != 1) {
        VT_FreePointList(&p_in_trsf);
        fprintf(stderr, "%s: unable to apply initial transformation of p_in...\n", proc);
        return(0);
    }

    /* Calcul des appariements de points */

    if (VT_ComputePairsOfPointsWithTrees(p_in_trsf, nodes_ext, pairs, &npairs, nClasses, &r) != 1)
    {
        VT_FreePointList(&p_in_trsf);
        fprintf(stderr, "%s: unable to compute the pairings...\n", proc);
        return(0);
    }

    /* Copie éventuelle de la transformation initiale */

    /* if ( flag_affine == 2 )
      for (i=0 ; i<16 ; i++)
        T_ext_in_init_cpy[i]=T_ext_in_init[i]; */

    /* Calcul de la nouvelle transformation */

    /* if( (flag_affine != 2 && VT_ComputePointListsListTrsf(pairs, npairs, p_in, p_ext_list, &r, &sd, T_ext_in, flag_affine, estimator) != 1) ||
       (flag_affine == 2 && VT_ComputePointListsListTrsf(pairs, npairs, p_in_trsf, p_ext_list, &r, &sd, T_ext_in, flag_affine, estimator) != 1) ) { */
    if( VT_ComputePointListsListTrsf(pairs, npairs, p_in, p_ext_list, &r, &sd, T_ext_in, flag_affine, estimator) != 1 ) {
        VT_FreePointList(&p_in_trsf);
        fprintf(stderr, "%s: unable to compute new transformation...\n", proc);
        vtfree(pairs[0]); pairs[0]=NULL;
        vtfree(pairs[1]); pairs[1]=NULL;
        return(0);
    }

    /*if ( flag_affine == 2 )
    {
        for (i = 0 ; i < 3 ; i++)
        {
          for (j = 0 ; j < 3 ; j++)
              T_ext_in[j+4*i] = T_ext_in_init_cpy[j+4*i];
          T_ext_in[4+4*i] += T_ext_in_init_cpy[j+4*i];
        }
    }*/

    *residual = r;

    /* Frees */
    VT_FreePointList(&p_in_trsf);
    vtfree(pairs[0]); pairs[0]=NULL;
    vtfree(pairs[1]); pairs[1]=NULL;

    return(1);
}



/*partitionner(tableau T, entier premier, entier dernier, pivot)
echanger T[pivot] et T[dernier]
j := premier
pour i de premier a dernier - 1
    si T[i] <= T[dernier] alors
        echanger T[i] et T[j]
        j := j + 1
echanger T[dernier] et T[j]
renvoyer j
*/
static int partitionner(vt_pointList p, int axis, int premier, int dernier, int pivot, int *ind)
{
   int i, j, tmp;
   vt_point *pi, *pp;
   int echanger;

   tmp=ind[dernier];
   ind[dernier] = ind[pivot];
   ind[pivot]=tmp;
   j = premier;

   pp = &(p.list[ind[dernier]]);

   for (i=premier ; i<dernier ; i++)
   {
       echanger = 0;
       pi = &(p.list[ind[i]]);

       switch(axis)
       {
       default:
       case 0 :
	 /*  X axis */
           if (pi->x <= pp->x)
               echanger = 1;
           break;
       case 1 :
	 /*  Y axis */
           if (pi->y <= pp->y)
               echanger = 1;
           break;
       case 2 :
	 /*  Z axis */
           if (pi->z <= pp->z)
               echanger = 1;
           break;
       }

       if (echanger)
       {
           tmp = ind[j];
           ind[j] = ind[i];
           ind[i] = tmp;
           j++;
       }
   }
   tmp = ind[j];
   ind[j] = ind[dernier];
   ind[dernier] = tmp;


   return(j);
}

/*partitionner(tableau T, entier premier, entier dernier, pivot)
echanger T[pivot] et T[dernier]
j := premier
pour i de premier a dernier - 1
    si T[i] <= T[dernier] alors
        echanger T[i] et T[j]
        j := j + 1
echanger T[dernier] et T[j]
renvoyer j
*/
static int partitionner_vector(double *v, int premier, int dernier, int pivot, int *ind)
{
   int i, j, tmp;
   double vi, vp;

   tmp=ind[dernier];
   ind[dernier] = ind[pivot];
   ind[pivot]=tmp;
   j = premier;

   vp = v[ind[dernier]];

   for (i=premier ; i<dernier ; i++)
   {
       vi = v[ind[i]];

       if (vi<=vp)
       {
           tmp = ind[j];
           ind[j] = ind[i];
           ind[i] = tmp;
           j++;
       }
   }
   tmp = ind[j];
   ind[j] = ind[dernier];
   ind[dernier] = tmp;

   return(j);
}


int choix_pivot( int premier, int dernier )
{
  /* return(dernier); */
    return (premier + (rand () % (dernier-premier+1)));
}

/*tri_rapide(tableau T, entier premier, entier dernier)
    debut
        si premier < dernier alors
            pivot := choix_pivot(T, premier, dernier)
            pivot := partitionner(T, premier, dernier, pivot)
            tri_rapide(T, premier, pivot-1)
            tri_rapide(T, pivot+1, dernier)
        fin si
    fin
*/
void tri_rapide(vt_pointList p, int axis, int premier, int dernier, int *ind)
{
  /*  Rq : le tableau ind doit contenir toutes les valeurs entre 0 et p.n-1, eventuellement permutees */
  /*  tri_rapide retournera le tableau ind avec les indices tels que p.list[ind[:]] soit trie selon l'axe X (axis=0), Y (axis=1) ou Z (axis=2) */
    int pivot;



    if(premier < dernier)
    {
        pivot = choix_pivot(premier, dernier);
        pivot = partitionner(p, axis, premier, dernier, pivot, ind);
        tri_rapide(p, axis, premier, pivot-1, ind);
        tri_rapide(p, axis, pivot+1, dernier, ind);
    }
}

void tri_rapide_vector(double *v, int premier, int dernier, int *ind)
{
  /*  Rq : le tableau ind doit contenir toutes les valeurs entre 0 et n-1, eventuellement permutees */
  /*  tri_rapide_vector retournera le tableau ind avec les indices tels que v[ind[:]] soit trie par ordre croissant */

  /* fprintf(stderr, "   tri_rapide_vector : %d, %d\n", premier, dernier); */

    int pivot;

    if(premier < dernier)
    {
        pivot = choix_pivot(premier, dernier);
        pivot = partitionner_vector(v, premier, dernier, pivot, ind);
        tri_rapide_vector(v, premier, pivot-1, ind);
        tri_rapide_vector(v, pivot+1, dernier, ind);
    }
}

int VT_TriRapideVector(double *v, int n, int **Ind)
{
    char *proc="VT_TriRapideVector";

    /*  VT_TriRapideVector : Ind pointera sur le tableau d'indices ind tels que v[ind[:]] soit trie par ordre croissant */
    int i;
    int *ind = vtmalloc(n * sizeof(int), "ind", proc );

    if (ind==NULL)
    {
        fprintf(stderr, "%s: unable to allocate a vector\n", proc);
        return(0);
    }

    for (i = 0 ; i<n ; i++) {
        ind[i] = i;
    }

    tri_rapide_vector(v, 0, n-1, ind);

    *Ind=ind;
    return(1);
}


int beforeAfterSubsets(vt_pointList p, int axis, vt_pointList *pL, vt_pointList *pR, vt_point *median)
{
    char *proc="beforeAfterSubsets";


    if (p.n == 0)    {
        pL = NULL;
        pR = NULL;
        return(1);
    }

    int m = (int)(p.n/2); /*  m = indice de l'elementt median apres tri */
    int i;

    vt_point *point;
    vt_point *cible;

    int *ind = vtmalloc(p.n * sizeof(int), "ind", proc );


    if (ind == NULL)
    {
        fprintf(stderr, "%s: unable to allocate index table of length %d\n", proc, p.n);
        return(0);
    }


    for (i=0 ; i<p.n ; i++) ind[i]=i; /*  init ind table (faire un vt_triRapide qui engloberait l'initialisation et les appels recursifs?) */
    tri_rapide(p, axis, 0, p.n-1, ind);



    if (VT_AllocPointList(pL, m) != 1) {
        fprintf(stderr, "%s: unable to allocate point list left\n", proc);
        return(0);
    }


    if (VT_AllocPointList(pR, p.n-1-m) != 1) {
        VT_FreePointList(pL);
        fprintf(stderr, "%s: unable to allocate point list right\n", proc);
        return(0);
    }


    point = &(p.list[ind[m]]);
    median->attribute = point->attribute;
    median->weight = point->weight;
    median->x = point->x;
    median->y = point->y;
    median->z = point->z;


    for (i = 0 ; i<m ; i++)
    {
        point = &(p.list[ind[i]]);
        cible = &(pL->list[i]);
        cible->attribute=point->attribute;
        cible->weight=point->weight;
        cible->x = point->x;
        cible->y = point->y;
        cible->z = point->z;
    }


    for (i = 0 ; i<p.n-1-m ; i++)
    {
        point = &(p.list[ind[i+m+1]]);
        cible = &(pR->list[i]);
        cible->attribute=point->attribute;
        cible->weight=point->weight;
        cible->x = point->x;
        cible->y = point->y;
        cible->z = point->z;
    }


    vtfree(ind);

    return(1);
}






void vt_kdtree(tree_node **Node, vt_pointList p, int depth)
{
    char *proc="vt_kdtree";


    tree_node *node;
    vt_point med;
    vt_point *median;
    vt_pointList pL, pR;
    int axis = depth % 3; /*  3-D tree */

    node = vtmalloc(sizeof(tree_node), "node", proc );

    if(p.n == 0) {
        node = NULL;
        *Node = node;
        return;
    }




    if (beforeAfterSubsets(p, axis, &pL, &pR, &med) != 1)
    {
        node = NULL;
        *Node = node;
        fprintf(stderr, "%s: error while computing the kd tree at depth = %d\n", proc, depth);
        return;
    }

    median = &(node->median);

    median->x = med.x;
    median->y = med.y;
    median->z = med.z;

    median->attribute = med.attribute;
    median->weight = med.weight;
    node->axis = axis;

    vt_kdtree(&(node->leftChild), pL, depth+1);

    vt_kdtree(&(node->rightChild), pR, depth+1);


    *Node = node;

    VT_FreePointList(&pL);
    VT_FreePointList(&pR);

}


void vt_kdtrees(tree_node ***Nodes, vt_pointList *pList, int n)
{
    char *proc="vt_kdtrees";
    int i;

    tree_node **nodes = (*Nodes);
    nodes = vtmalloc( n * sizeof(tree_node *), "nodes", proc );
    if (nodes == NULL) {
        fprintf(stderr, "%s : unable to allocate nodees vector\n", proc);
        return;
    }
    for (i=0 ; i<n ; i++) nodes[i]=NULL;
    /* vt_pointList p; */

    int depth = 0;
    /* int axis = depth % 3; // 3-D tree */


    for (i=0 ; i<n ; i++) {

      /* tree_node *node = nodes[i]; */

      vt_kdtree(&(nodes[i]), pList[i], depth);

      /*vt_point med;
      vt_point *median;
      vt_pointList pL, pR;

      p = pList[i];

      Node = nodeList[i];
      node = vtmalloc(sizeof(tree_node), "node", proc );

      if(p.n == 0) {
          node = NULL;
          *Node = node;
          continue;
      }

      if (beforeAfterSubsets(p, axis, &pL, &pR, &med) != 1)
      {
          node = NULL;
          *Node = node;
          fprintf(stderr, "%s: error while computing the kd tree associated to value %d\n", proc, i);
          continue;
      }

      median = &(node->median);

      median->x = med.x;
      median->y = med.y;
      median->z = med.z;

      median->attribute = med.attribute;
      median->weight = med.weight;
      node->axis = axis;

      vt_kdtree(&(node->leftChild), pL, depth+1);

      vt_kdtree(&(node->rightChild), pR, depth+1);

      */
      /* *Node = node; */

      /* VT_FreePointList(&pL); */
      /* VT_FreePointList(&pR); */
    }

    *Nodes = nodes;
    return;
}


void printTree(tree_node *node)
{
    vt_point *pt;

    if(!node) return;

    if(node->leftChild)  printTree(node->leftChild);

    pt = &(node->median);

    printf("median = { %1.1f, %1.1f, %1.1f }, axis = %d, weight = %1.1f, attribute = %d, value = %d\n", pt->x, pt->y, pt->z, node->axis, pt->weight, pt->attribute, pt->value);

    if(node->rightChild) printTree(node->rightChild);
}


void clearTree(tree_node **p)
{

    tree_node *tmpP = *p;

    /* fprintf(stdout, "clearTree: ") */

    if(!tmpP) return;

    if(tmpP->leftChild)  clearTree(&tmpP->leftChild);

    if(tmpP->rightChild) clearTree(&tmpP->rightChild);

    vtfree(tmpP);

    *p = NULL;
}





int VT_ComputePointListsListTrsf(int **pairs, int n, vt_pointList p_1, vt_pointList *p_list,
                                    double *r, double *sd, double *T,
                                    int flag_trsf, bal_estimator estimator)
{
    char *proc="VT_ComputePointListsListTrsf";
    int i;

    vt_point theP;
    vt_pointList p_2;

    /*  Debut des hostilites */



    bal_transformation theResTransformation;

    bal_typeFieldPointList theFloatingPoints;
    bal_typeFieldPointList theReferencePoints;

    enumTypeTransfo transformation_type=RIGID_3D;
    if (flag_trsf==1) transformation_type=AFFINE_3D;
    if (flag_trsf==2) transformation_type=TRANSLATION_3D;
    if (flag_trsf==3) transformation_type=SIMILITUDE_3D;

    enumUnitTransfo points_unit = REAL_UNIT;


    FIELD theField;



    BAL_InitTypeFieldPointList( &theFloatingPoints );
    BAL_InitTypeFieldPointList( &theReferencePoints );
    BAL_InitTransformation( &theResTransformation );


    for (i=0; i<n ; i++)
    {
        if (pairs[1][i]<0) continue;

        bal_typeFieldPoint c;
        theP=p_1.list[pairs[0][i]];
        c.x = theP.x;
        c.y = theP.y;
        c.z = theP.z;
        if ( BAL_AddTypeFieldPointToTypeFieldPointList( &theReferencePoints, &c ) != 1 ) {
              fprintf( stderr, "%s: unable to add measure #%d (=%f %f %f) to list\n",
               proc, i, c.x, c.y, c.z );
            if (i>0){
                BAL_FreeTypeFieldPointList( &theFloatingPoints );
                BAL_FreeTypeFieldPointList( &theReferencePoints );
            }
            return( -1 );
        }

        p_2=p_list[theP.value];
        theP=p_2.list[pairs[1][i]];
        c.x = theP.x;
        c.y = theP.y;
        c.z = theP.z;
        if ( BAL_AddTypeFieldPointToTypeFieldPointList( &theFloatingPoints, &c ) != 1 ) {
              fprintf( stderr, "%s: unable to add measure #%d (=%f %f %f) to list\n",
               proc, i, c.x, c.y, c.z );
            BAL_FreeTypeFieldPointList( &theFloatingPoints );
            if (i>0)
                BAL_FreeTypeFieldPointList( &theReferencePoints );
            return( -1 );
        }
    }

    if ( theFloatingPoints.n_data != theReferencePoints.n_data ) {
      fprintf( stderr, "%s: list of points have different lengths\n", proc );
      fprintf( stderr, "   list 1 contains %d points\n", theFloatingPoints.n_data );
      fprintf( stderr, "   list 2 contains %d points\n", theReferencePoints.n_data );
      if ( theFloatingPoints.n_data > theReferencePoints.n_data ){
        theFloatingPoints.n_data = theReferencePoints.n_data;
        fprintf( stderr, "   uses the first %d of both lists\n", theReferencePoints.n_data );
      }
      else  {
        theReferencePoints.n_data = theFloatingPoints.n_data;
        fprintf( stderr, "   uses the first %d of both lists\n", theFloatingPoints.n_data );
      }
    }

    /*  processing */

    /*  initializating transformation */

    theResTransformation.type = transformation_type;
    switch ( transformation_type ) {
    default :
      fprintf(stderr, "such transformation type non handled yet\n");
      return(0);
      break;

    case TRANSLATION_2D :
    case TRANSLATION_3D :
    case TRANSLATION_SCALING_2D :
    case TRANSLATION_SCALING_3D :
    case RIGID_2D :
    case RIGID_3D :
    case SIMILITUDE_2D :
    case SIMILITUDE_3D :
    case AFFINE_2D :
    case AFFINE_3D :
      if ( BAL_AllocTransformation( &theResTransformation, transformation_type, NULL ) != 1 ) {
        BAL_FreeTypeFieldPointList( &theFloatingPoints );
        BAL_FreeTypeFieldPointList( &theReferencePoints );
        fprintf(stderr, "%s: error when allocating linear transformation\n" , proc);
        return(0);
      }
      break;
    }


    if ( BAL_AllocateField( &theField, theFloatingPoints.n_data ) != 1 ) {
      BAL_FreeTypeFieldPointList( &theFloatingPoints );
      BAL_FreeTypeFieldPointList( &theReferencePoints );
      BAL_FreeTransformation( &theResTransformation );
      fprintf( stderr,"%s: error when allocating displacement field\n", proc);
      return(0);
    }
    /* TODO VERIFIER ICI */
    theField.vx = 1;
    theField.vy = 1;
    theField.vz = 1;


    theFloatingPoints.unit =  points_unit;
    theReferencePoints.unit =  points_unit;


    /* won't work in the points are in voxel units, with different
       voxel sizes
    */
    if ( BAL_ComputePairingFieldFromTypeFieldPointList( &theField,
                           &theFloatingPoints,
                           &theReferencePoints ) != 1 ) {
      BAL_FreeField( &theField );
      BAL_FreeTypeFieldPointList( &theFloatingPoints );
      BAL_FreeTypeFieldPointList( &theReferencePoints );
      BAL_FreeTransformation( &theResTransformation );
      fprintf( stderr,"%s: error when computing displacement field\n", proc);
      return(0);
    }


    /* transformation computation procedure assumes that
       the points are in voxel units (in the reference referential)
    */
    if ( points_unit == REAL_UNIT ) {
      for ( i=0; i<(int)theField.n_computed_pairs; i++ ) {
        theField.data[i].origin.x /= theField.vx;
        theField.data[i].origin.y /= theField.vy;
        theField.data[i].origin.z /= theField.vz;
        theField.data[i].vector.x /= theField.vx;
        theField.data[i].vector.y /= theField.vy;
        theField.data[i].vector.z /= theField.vz;
      }
    }

    /* computes transformation
     */

    if ( BAL_ComputeIncrementalTransformation( &theResTransformation,
                           &theField, &(estimator) ) != 1 ) {
      BAL_FreeField( &theField );
      BAL_FreeTypeFieldPointList( &theFloatingPoints );
      BAL_FreeTypeFieldPointList( &theReferencePoints );
      BAL_FreeTransformation( &theResTransformation );
      fprintf( stderr,"%s: error when computing transformation\n", proc);
    }


    BAL_FreeTypeFieldPointList( &theFloatingPoints );
    BAL_FreeTypeFieldPointList( &theReferencePoints );



    /* computes residuals
     */
    if ( 1 ) {
      if ( BAL_ComputeTransformationResiduals( &theResTransformation,
                           &theField ) != 1 ) {
        BAL_FreeTransformation( &theResTransformation );
        fprintf(stderr, "%s: error when computing residuals\n", proc);
      }
      {
        double sum = 0.0;
        double s1=0.0;

        for ( i=0; i<(int)theField.n_selected_pairs; i++ ){
            sum += sqrt( theField.pointer[i]->error );
            s1+=theField.pointer[i]->error;
        }

        if(0) {
            fprintf( stdout, "# average on %lu points: \t%f\n", theField.n_selected_pairs , sum / (double)theField.n_selected_pairs );
            fprintf( stdout, "----------------------------------------\n");
        }
        *r = sum/ ((double)theField.n_selected_pairs ); /*  residu moyen */
        *sd= sqrt(s1/(double)theField.n_selected_pairs - (sum/(double)theField.n_selected_pairs) * (sum/(double)theField.n_selected_pairs));
      }
    }

    BAL_FreeField( &theField );


    /* writing transformation
     */

      int l,c;
      l=theResTransformation.mat.l;
      c=theResTransformation.mat.c;
      if (l!=4 || c!=4) {
          fprintf( stderr, "%s: unexpected size of transformation matrix\n", proc);
          return(0);
      }


      for (i=0; i<c*l; i++) {
          T[i]=theResTransformation.mat.m[i];
      }

    /* BAL_PrintTransformation(stderr, &theResTransformation, NULL); */
    BAL_FreeTransformation( &theResTransformation );

    /*  Fin des hostilites */

    return(1);
}

int VT_ComputePointListsTrsf(int **pairs, int n, vt_pointList p_1, vt_pointList p_2,
                        double *r, double *sd, double *T,
                        int flag_trsf, bal_estimator estimator)
{
    char *proc="VT_ComputePointListsTrsf";
    int i;

    vt_point theP;


    /*  Debut des hostilites */



    bal_transformation theResTransformation;

    bal_typeFieldPointList theFloatingPoints;
    bal_typeFieldPointList theReferencePoints;

    enumTypeTransfo transformation_type=RIGID_3D;
    if (flag_trsf==1) transformation_type=AFFINE_3D;
    if (flag_trsf==2) transformation_type=TRANSLATION_3D;
    if (flag_trsf==3) transformation_type=SIMILITUDE_3D;

    enumUnitTransfo points_unit = REAL_UNIT;


    FIELD theField;



    BAL_InitTypeFieldPointList( &theFloatingPoints );
    BAL_InitTypeFieldPointList( &theReferencePoints );
    BAL_InitTransformation( &theResTransformation );


    for (i=0; i<n ; i++)
    {
        bal_typeFieldPoint c;
        theP=p_1.list[pairs[0][i]];
        c.x = theP.x;
        c.y = theP.y;
        c.z = theP.z;
        if ( BAL_AddTypeFieldPointToTypeFieldPointList( &theReferencePoints, &c ) != 1 ) {
              fprintf( stderr, "%s: unable to add measure #%d (=%f %f %f) to list\n",
               proc, i, c.x, c.y, c.z );
            if (i>0){
                BAL_FreeTypeFieldPointList( &theFloatingPoints );
                BAL_FreeTypeFieldPointList( &theReferencePoints );
            }
            return( -1 );
        }

        theP=p_2.list[pairs[1][i]];
        c.x = theP.x;
        c.y = theP.y;
        c.z = theP.z;
        if ( BAL_AddTypeFieldPointToTypeFieldPointList( &theFloatingPoints, &c ) != 1 ) {
              fprintf( stderr, "%s: unable to add measure #%d (=%f %f %f) to list\n",
               proc, i, c.x, c.y, c.z );
            BAL_FreeTypeFieldPointList( &theFloatingPoints );
            if (i>0)
                BAL_FreeTypeFieldPointList( &theReferencePoints );
            return( -1 );
        }
    }

    if ( theFloatingPoints.n_data != theReferencePoints.n_data ) {
      fprintf( stderr, "%s: list of points have different lengths\n", proc );
      fprintf( stderr, "   list 1 contains %d points\n", theFloatingPoints.n_data );
      fprintf( stderr, "   list 2 contains %d points\n", theReferencePoints.n_data );
      if ( theFloatingPoints.n_data > theReferencePoints.n_data ){
        theFloatingPoints.n_data = theReferencePoints.n_data;
        fprintf( stderr, "   uses the first %d of both lists\n", theReferencePoints.n_data );
      }
      else  {
        theReferencePoints.n_data = theFloatingPoints.n_data;
        fprintf( stderr, "   uses the first %d of both lists\n", theFloatingPoints.n_data );
      }
    }

    /*  processing */

    /*  initializating transformation */

    theResTransformation.type = transformation_type;
    switch ( transformation_type ) {
    default :
      fprintf(stderr, "such transformation type non handled yet\n");
      return(0);
      break;

    case TRANSLATION_2D :
    case TRANSLATION_3D :
    case TRANSLATION_SCALING_2D :
    case TRANSLATION_SCALING_3D :
    case RIGID_2D :
    case RIGID_3D :
    case SIMILITUDE_2D :
    case SIMILITUDE_3D :
    case AFFINE_2D :
    case AFFINE_3D :
      if ( BAL_AllocTransformation( &theResTransformation, transformation_type, NULL ) != 1 ) {
        BAL_FreeTypeFieldPointList( &theFloatingPoints );
        BAL_FreeTypeFieldPointList( &theReferencePoints );
        fprintf(stderr, "%s: error when allocating linear transformation\n" , proc);
        return(0);
      }
      break;
    }


    if ( BAL_AllocateField( &theField, theFloatingPoints.n_data ) != 1 ) {
      BAL_FreeTypeFieldPointList( &theFloatingPoints );
      BAL_FreeTypeFieldPointList( &theReferencePoints );
      BAL_FreeTransformation( &theResTransformation );
      fprintf( stderr,"%s: error when allocating displacement field\n", proc);
      return(0);
    }

    theField.vx = 1;
    theField.vy = 1;
    theField.vz = 1;


    theFloatingPoints.unit =  points_unit;
    theReferencePoints.unit =  points_unit;


    /* won't work in the points are in voxel units, with different
       voxel sizes
    */
    if ( BAL_ComputePairingFieldFromTypeFieldPointList( &theField,
                           &theFloatingPoints,
                           &theReferencePoints ) != 1 ) {
      BAL_FreeField( &theField );
      BAL_FreeTypeFieldPointList( &theFloatingPoints );
      BAL_FreeTypeFieldPointList( &theReferencePoints );
      BAL_FreeTransformation( &theResTransformation );
      fprintf( stderr,"%s: error when computing displacement field\n", proc);
      return(0);
    }


    /* transformation computation procedure assumes that
       the points are in voxel units (in the reference referential)
    */
    if ( points_unit == REAL_UNIT ) {
      for ( i=0; i<(int)theField.n_computed_pairs; i++ ) {
        theField.data[i].origin.x /= theField.vx;
        theField.data[i].origin.y /= theField.vy;
        theField.data[i].origin.z /= theField.vz;
        theField.data[i].vector.x /= theField.vx;
        theField.data[i].vector.y /= theField.vy;
        theField.data[i].vector.z /= theField.vz;
      }
    }

    /* computes transformation
     */

    if ( BAL_ComputeIncrementalTransformation( &theResTransformation,
                           &theField, &(estimator) ) != 1 ) {
      BAL_FreeField( &theField );
      BAL_FreeTypeFieldPointList( &theFloatingPoints );
      BAL_FreeTypeFieldPointList( &theReferencePoints );
      BAL_FreeTransformation( &theResTransformation );
      fprintf( stderr,"%s: error when computing transformation\n", proc);
    }


    BAL_FreeTypeFieldPointList( &theFloatingPoints );
    BAL_FreeTypeFieldPointList( &theReferencePoints );



    /* computes residuals
     */
    if ( 1 ) {
      if ( BAL_ComputeTransformationResiduals( &theResTransformation,
                           &theField ) != 1 ) {
        BAL_FreeTransformation( &theResTransformation );
        fprintf(stderr, "%s: error when computing residuals\n", proc);
      }
      {
        double sum = 0.0;
        double s1=0.0;

        for ( i=0; i<(int)theField.n_selected_pairs; i++ ){
            sum += sqrt( theField.pointer[i]->error );
            s1+=theField.pointer[i]->error;
        }

        if(0) {
            fprintf( stdout, "# average on %lu points: \t%f\n", theField.n_selected_pairs , sum / (double)theField.n_selected_pairs );
            fprintf( stdout, "----------------------------------------\n");
        }
        *r = sum/ ((double)theField.n_selected_pairs ); /*  residu moyen */
        *sd= sqrt(s1/(double)theField.n_selected_pairs - (sum/(double)theField.n_selected_pairs) * (sum/(double)theField.n_selected_pairs));
      }
    }

    BAL_FreeField( &theField );


    /* writing transformation
     */

      int l,c;
      l=theResTransformation.mat.l;
      c=theResTransformation.mat.c;
      if (l!=4 || c!=4) {
          fprintf( stderr, "%s: unexpected size of transformation matrix\n", proc);
          return(0);
      }


      for (i=0; i<c*l; i++) {
          T[i]=theResTransformation.mat.m[i];
      }

    BAL_FreeTransformation( &theResTransformation );

    /*  Fin des hostilites */

    return(1);
}

void VT_PointListBarycentreInitTrsf(vt_pointList p, double *T, double *B)
{
    vt_point *theP;
    int i;
    double w=0;
    for (i=0; i<p.n; i++)
    {
      theP=&(p.list[i]);
      B[0]+=(theP->x*T[0]+theP->y*T[1]+theP->z*T[2]+T[3])*theP->weight;
      B[1]+=(theP->x*T[4]+theP->y*T[5]+theP->z*T[6]+T[7])*theP->weight;
      B[2]+=(theP->x*T[8]+theP->y*T[9]+theP->z*T[10]+T[11])*theP->weight;
      w+=theP->weight;
    }
    if (w != 0)
      for (i=0 ; i<3 ; i++)
        B[i]/=w;
}

void VT_TreeBarycentreRec(tree_node *node, double *B, double *w)
{
    if (node == NULL) return;

    VT_TreeBarycentreRec(node->leftChild, B, w);
    VT_TreeBarycentreRec(node->rightChild, B, w);

    vt_point *theP=&(node->median);

    B[0]+=theP->x*theP->weight;
    B[1]+=theP->y*theP->weight;
    B[2]+=theP->z*theP->weight;

    *w = *w+theP->weight;
}

void VT_TreeBarycentre(tree_node *node, double *B, double *weight)
{
    double w=0;
    B[0]=0; B[1]=0; B[2]=0;
    VT_TreeBarycentreRec(node, B, &w);
    if (w!=0) {
      B[0] /= w;
      B[1] /= w;
      B[2] /= w;
    }
    *weight=w;
}

int VT_TransformPointList(double *T, vt_pointList p, vt_pointList *p_trsf)
{
    char *proc="VT_TransformPointList";
    int i;
    vt_point *theP, *theTrsf;
    double xval,yval,zval;
    if( p.n != p_trsf->n)
    {
        fprintf(stderr, "%s: conflicting field n: original = %d, trsf = %d\n", proc, p.n, p_trsf->n);
        return(0);
    }
    for (i=0; i<p.n; i++)
    {
        theP=&(p.list[i]);
        xval=theP->x;
        yval=theP->y;
        zval=theP->z;
        theTrsf=&(p_trsf->list[i]);
        theTrsf->x=(xval*T[0]+yval*T[1]+zval*T[2]+T[3]);
        theTrsf->y=(xval*T[4]+yval*T[5]+zval*T[6]+T[7]);
        theTrsf->z=(xval*T[8]+yval*T[9]+zval*T[10]+T[11]);
        theTrsf->weight=theP->weight;
        theTrsf->value=theP->value;
        theTrsf->attribute=theP->attribute;
    }

    return(1);
}

int VT_ComputePairsOfPointsWithTree(vt_pointList p_1, tree_node *node_2, int **p, int *npairs, double *dnorm)
{
    char *proc="VT_ComputePairsOfPointsWithTree";
    int n=p_1.n;
    int i;
    vt_point theP;
    int *pairs[2];

    int N;
    double d =0.0;
    double *D, dsum;


    D = vtmalloc( n*sizeof(double), "D", proc );
    if(D==NULL) {

        fprintf(stderr, "%s: unable to allocate distance vector\n", proc);
        return(0);
    }
    pairs[0] = vtmalloc( n*sizeof(int), "pairs[0]", proc );
    if(pairs[0]==NULL) {
        vtfree(D);
        D=NULL;

        fprintf(stderr, "%s: unable to allocate pairs first vector\n", proc);
        return(0);
    }
    pairs[1] = vtmalloc( n*sizeof(int), "pairs[1]", proc );
    if(pairs[1]==NULL) {
        vtfree(D);
        D=NULL;
        vtfree(pairs[0]);
        pairs[0]=NULL;

        fprintf(stderr, "%s: unable to allocate pairs second vector\n", proc);
        return(0);
    }

    for (i=0 ; i<p_1.n ; i++) {
        theP=p_1.list[i];
        vt_point tmp;
        VT_NNS(theP, node_2, &tmp, &d);
        int j=tmp.attribute;
        if (j<0)
        {
            vtfree(D);
            D=NULL;
            vtfree(pairs[0]);
            pairs[0]=NULL;
            vtfree(pairs[1]);
            pairs[1]=NULL;

            fprintf(stderr, "%s: unexpected result of VT_NNS...\n", proc);
            return(0);
        }
        pairs[0][i]=i;
        pairs[1][i]=j;
        D[i]=d;
    }

    N=n;
    p[0] = vtmalloc(N*sizeof(int), "p[0]", proc );
    if(p[0]==NULL) {
        vtfree(D);
        D=NULL;
        vtfree(pairs[0]);
        pairs[0]=NULL;
        vtfree(pairs[1]);
        pairs[1]=NULL;

        fprintf(stderr, "%s: unable to allocate the p[0] vector...\n", proc);
        return(0);
    }
    p[1] = vtmalloc(N*sizeof(int), "p[1]", proc );
    if(p[1]==NULL) {
        vtfree(D);
        D=NULL;
        vtfree(pairs[0]);
        pairs[0]=NULL;
        vtfree(pairs[1]);
        pairs[1]=NULL;
        vtfree(p[0]);
        p[0]=NULL;

        fprintf(stderr, "%s: unable to allocate the p[1] vector...\n", proc);
        return(0);
    }
    dsum=0;
    for(i=0; i<N; i++) {
        p[0][i]=pairs[0][i];
        p[1][i]=pairs[1][i];
        dsum+=D[i];
    }
    if(N>0)
        dsum/=N;
    else dsum=-1;

    vtfree(D);
    D=NULL;
    vtfree(pairs[0]);
    pairs[0]=NULL;
    vtfree(pairs[1]);
    pairs[1]=NULL;
    /* vtfree(IND); */
    /* IND=NULL; */

    *dnorm=dsum;
    *npairs=N;
    return(1);
}

int VT_ComputePairsOfPointsWithTrees(vt_pointList p_1, tree_node **node_list, int **p, int *npairs, int nClasses, double *dnorm)
{
    char *proc="VT_ComputePairsOfPointsWithTrees";
    int n=p_1.n;
    int i;
    vt_point theP;
    int *pairs[2];

    int N, nApparied=0;
    double d, dsum = 0;


    pairs[0] = vtmalloc( n*sizeof(int), "pairs[0]", proc );
    if(pairs[0]==NULL) {
        fprintf(stderr, "%s: unable to allocate pairs first vector\n", proc);
        return(0);
    }
    pairs[1] = vtmalloc( n*sizeof(int), "pairs[1]", proc );
    if(pairs[1]==NULL) {
        vtfree(pairs[0]);
        pairs[0]=NULL;
        fprintf(stderr, "%s: unable to allocate pairs second vector\n", proc);
        return(0);
    }

    for (i=0 ; i<p_1.n ; i++) {
        theP=p_1.list[i];
        vt_point tmp;
        if (theP.value >= nClasses || VT_NNS(theP, node_list[theP.value], &tmp, &d) <= 0)
        {
            pairs[0][i]=i;
            pairs[1][i]=-1;
            fprintf(stderr, "%s: warning : unexpected result of VT_NNS : class %d may not exist in ext point set\n", proc, theP.value);
        }
        else {
            int j=tmp.attribute;
            if (j<0)
            {
                vtfree(pairs[0]);
                pairs[0]=NULL;
                vtfree(pairs[1]);
                pairs[1]=NULL;

                fprintf(stderr, "%s: unexpected result of VT_NNS \n", proc);
                return(0);
            }
            pairs[0][i]=i;
            pairs[1][i]=j;
            dsum += d;
            nApparied++;
        }
    }

    N=n;
    p[0] = vtmalloc(N*sizeof(int), "p[0]", proc );
    if(p[0]==NULL) {
        vtfree(pairs[0]);
        pairs[0]=NULL;
        vtfree(pairs[1]);
        pairs[1]=NULL;

        fprintf(stderr, "%s: unable to allocate the p[0] vector...\n", proc);
        return(0);
    }
    p[1] = vtmalloc(N*sizeof(int), "p[1]", proc );
    if(p[1]==NULL) {
        vtfree(pairs[0]);
        pairs[0]=NULL;
        vtfree(pairs[1]);
        pairs[1]=NULL;
        vtfree(p[0]);
        p[0]=NULL;

        fprintf(stderr, "%s: unable to allocate the p[1] vector...\n", proc);
        return(0);
    }
    for(i=0; i<N; i++) {
        p[0][i]=pairs[0][i];
        p[1][i]=pairs[1][i];
    }
    if(nApparied>0)
        dsum/=nApparied;
    else dsum=-1;

    vtfree(pairs[0]);
    pairs[0]=NULL;
    vtfree(pairs[1]);
    pairs[1]=NULL;

    *dnorm=dsum;
    *npairs=N;
    return(1);
}

static double squareDist(vt_point p, vt_point q)
{
    return((p.x-q.x)*(p.x-q.x)+(p.y-q.y)*(p.y-q.y)+(p.z-q.z)*(p.z-q.z));
}

static void NNS(vt_point q, tree_node *node, vt_point *best, double *w)
{
    vt_point *p;
    int search_first;


    p = &(node->median);
    search_first = 0;
    double qval, pval;
    switch(node->axis) {
    default:
    case 0:
        qval = q.x;
        pval = p->x;
        break;
    case 1:
        qval = q.y;
        pval = p->y;
        break;
    case 2:
        qval = q.z;
        pval = p->z;
        break;
    }
    if ( qval < pval )
        search_first = 0;
    else
        search_first = 1;

    if (search_first == 0)
    {
        if (node->leftChild) {
          NNS(q, node->leftChild, best, w);
        }
        double _w = squareDist(q, node->median);
        if(_w < (*w))
        {
          best->attribute = p->attribute;
          best->weight = p->weight;
          best->value = p->value;
          best->x = p->x;
          best->y = p->y;
          best->z = p->z;
          *w = _w;
        }
        if (sqrt(*w) > fabs(qval-pval) && node->rightChild)
        {
            NNS(q, node->rightChild,best,w);
        }
    }
    else
    {
        if (node->rightChild) {
          NNS(q, node->rightChild, best, w);
        }
        double _w = squareDist(q, node->median);
        if(_w < (*w))
        {
          best->attribute = p->attribute;
          best->weight = p->weight;
          best->value = p->value;
          best->x = p->x;
          best->y = p->y;
          best->z = p->z;
          *w = _w;
        }
        if (sqrt(*w) > fabs(qval-pval) && node->leftChild)
        {
            NNS(q, node->leftChild,best,w);
        }
    }
    return;
}

int VT_NNS(vt_point q, tree_node *node, vt_point *best, double *d)
{
    double w = 1000000.0;
    if (node==NULL) return(-1);
    NNS(q, node, best, &w);
    *d = sqrt(w);

    return(1);
}



int VT_ResidualsVector(vt_pointList p_in, tree_node *node_ext, double *T_ext_in, double **residuals)
{
    char *proc="VT_ResidualsVector";

    int i;
    vt_point P, Q;
    vt_pointList p_in_trsf;

    double d;
    double *r = NULL;

    r = vtmalloc( p_in.n * sizeof(double), "r", proc );

    if (r==NULL)
    {
        fprintf(stderr, "%s: unable to allocate residuals vector\n", proc);
        return(0);
    }

    if (VT_AllocPointList(&p_in_trsf, p_in.n) != 1) {
      vtfree(r); r=NULL;
      fprintf(stderr, "%s: unable to allocate point list p_in_trsf...\n", proc);
      return(0);
    }


    if (VT_TransformPointList(T_ext_in, p_in, &p_in_trsf) != 1) {
        vtfree(r); r=NULL;
        VT_FreePointList(&p_in_trsf);
        fprintf(stderr, "%s: unable to transform the point list\n", proc);
        return(0);
    }

    for (i=0 ; i<p_in.n ; i++)
    {
        Q = p_in_trsf.list[i];
        if(VT_NNS(Q, node_ext, &P, &d) != 1)
        {
            vtfree(r); r=NULL;
            VT_FreePointList(&p_in_trsf);
            fprintf(stderr, "%s: unable to allocate residuals vector\n", proc);
            return(0);
        }
        r[i]=d;
    }

    *residuals = r;
    VT_FreePointList(&p_in_trsf);
    return(1);
}



int VT_LabelsResidualsVector(vt_pointList p_in, tree_node **nodes_ext, double *T_ext_in, double **residuals)
{
    char *proc="VT_LabelsResidualsVector";

    int i;
    vt_point P, Q;

    double d;
    double *r = NULL;
    vt_pointList p_in_trsf;

    r = vtmalloc( p_in.n * sizeof(double), "r", proc );

    if (r==NULL)
    {
        fprintf(stderr, "%s: unable to allocate residuals vector\n", proc);
        return(0);
    }

    if (VT_AllocPointList(&p_in_trsf, p_in.n) != 1) {
      vtfree(r); r=NULL;
      fprintf(stderr, "%s: unable to allocate point list p_in_trsf...\n", proc);
      return(0);
    }

    if (VT_TransformPointList(T_ext_in, p_in, &p_in_trsf) != 1) {
        vtfree(r); r=NULL;
        VT_FreePointList(&p_in_trsf);
        fprintf(stderr, "%s: unable to transform the point list\n", proc);
        return(0);
    }

    for (i=0 ; i<p_in.n ; i++)
    {
        Q = p_in_trsf.list[i];
        if(VT_NNS(Q, nodes_ext[Q.value], &P, &d) != 1)
        {
            vtfree(r); r=NULL;
            VT_FreePointList(&p_in_trsf);
            fprintf(stderr, "%s: unable to allocate residuals vector\n", proc);
            return(0);
        }
        r[i]=d;
    }

    *residuals = r;
    VT_FreePointList(&p_in_trsf);

    return(1);
}

void VT_SetResidualImage(vt_image *imResiduals, vt_pointList p_in, double *residuals)
{
    int ind;
    int i,j,k;
    double x,y,z;
    vt_point *P;
    float ***array;
    switch(imResiduals->type ) {
    default:
        fprintf(stderr, "unexpected imResiduals parameter type : waiting for FLOAT type\n");
        return;
    case FLOAT:
        array=(float***)imResiduals->array;
    }

    for (ind = 0 ; ind < p_in.n ; ind++)
    {
        P = &(p_in.list[ind]);
        x=P->x;
        y=P->y;
        z=P->z;
        i=(int)(x/imResiduals->siz.x);
        j=(int)(y/imResiduals->siz.y);
        k=(int)(z/imResiduals->siz.z);
        if (i<0 || i>=(int)imResiduals->dim.x ||
            j<0 || j>=(int)imResiduals->dim.y ||
            k<0 || k>=(int)imResiduals->dim.z )
        {
            continue;
        }
        array[k][j][i] = (float)(residuals[ind]);
    }

    return;
}
