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

#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <vt_tree.h>



/*#include <bal-field.h>
//#include <bal-field-tools.h>
//#include <bal-point.h>
//#include <bal-transformation.h>
//#include <bal-transformation-tools.h>
*/
int VT_NPoints(vt_image *img)
{
    char *proc="VT_NPoints";

    unsigned char *bufu8=NULL;
    unsigned short int *bufu16=NULL;

    unsigned int i;
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

    unsigned int i;
    unsigned int n=0;
    int *h = NULL;
    unsigned int v = 0;

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
    h = malloc(n * sizeof(int));
    for (i=0 ; i<n ; i++) h[i] = 0;
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
    /*fprintf(stdout, "%s: nClasses = %d\n", proc, n);
    //for (i=0;i<n;i++)    {
    //    fprintf(stdout, "h[%d] = %d\t", i, h[i]);
    //}
    //fprintf(stdout, "\n");
    */
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

    h = malloc(n * sizeof(int));

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
                free(h);
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
    /*fprintf(stdout, "%s: nClasses = %d\n", proc, n);
    //for (i=0;i<n;i++)    {
    //    fprintf(stdout, "h[%d] = %d\t", i, h[i]);
    //}
    //fprintf(stdout, "\n");
    */
    return(1);
}


static int VT_InitPointListWithImage(vt_pointList* p, vt_image image)
{
    char *proc="VT_InitPointListWithImage";

    unsigned char *bufu8=NULL;
    unsigned short int *bufu16=NULL;

    unsigned int i;
    int val;
    unsigned int n=0;

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

    if(n != (unsigned int)p->n)
    {
            fprintf(stderr, "%s: incompatible vt_pointList field n = %d and vt_image number of points = %d\n", proc, p->n, n);
        return(0);
    }

    for (i=0; i<(unsigned int)p->n; i++) {
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
      p->list =malloc(n*sizeof(vt_point));
      if (p->list == (vt_point*)NULL) {
        fprintf(stderr, "%s: malloc failed\n", proc);
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
    vt_pointList *P = malloc(N * sizeof(vt_pointList));
    if (P == NULL)
    {
        fprintf(stderr, "%s: Unable to allocate point lists\n", proc);
        return(0);
    }

    for (i = 0; i < N ; i++) {
        /*fprintf(stdout, "%s : histo[%d] = %d\n", proc, i, histo[i]);*/
        VT_AllocPointList(&(P[i]), histo[i]);
    }

    *p = P;
    *n=N;
    free(histo);

    return(1);
}

int VT_AllocPointListsWithPointLists(vt_pointList **p, int n, vt_pointList *q)
{
    char *proc="VT_AllocPointListsWithPointLists";
    int i;
    vt_pointList *P = malloc(n * sizeof(vt_pointList));
    if (P == NULL)
    {
        fprintf(stderr, "%s: Unable to allocate point lists\n", proc);
        return(0);
    }
    for (i = 0; i < n ; i++) {
        VT_AllocPointList(&(P[i]), q[i].n);
    }
    *p = P;
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
    p->list =malloc(n*sizeof(vt_point));
    if (p->list == (vt_point*)NULL) {
        fprintf(stderr, "%s: malloc failed\n", proc);
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

    vt_pointList *P = malloc(N * sizeof(vt_pointList));
    if (P == NULL)
    {
        free(histo);
        fprintf(stderr, "%s: Unable to allocate point lists\n", proc);
        return(0);
    }

    for (i = 0; i < N ; i++) {
        VT_AllocPointList(&(P[i]), histo[i]);
    }

    *p = P;
    *n=N;
    free(histo);

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

    p->list =malloc(n*sizeof(vt_point));
    if (p->list == (vt_point*)NULL) {
        fprintf(stderr, "%s: malloc failed\n", proc);
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
        thePoint->x = (i+0.5)*img.siz.x;
        thePoint->y = (j+0.5)*img.siz.y;
        thePoint->z = (k+0.5)*img.siz.z;
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

    ind = malloc (nClasses * sizeof(int));
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

    /*fprintf(stderr, "%s: test\n", proc);

    //ind=0;
    */
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

        /*if(ind[val]<20) {
        //    fprintf(stderr, "%s: val = %d\tind[%d] = %d\n", proc, val, val, ind[val]);
        //}
        */
        p = &(p_list[val]);

        /*if(ind[val]<20) {
        //    fprintf(stderr, "%s: p_list[%d].n = %d\n", proc, val, p->n );
        //}
        */
        thePoint=&(p->list[ind[val]]);

        thePoint->weight = 1;
        thePoint->attribute = ind[val];
        thePoint->value = val;
        thePoint->x = (i+0.5)*img.siz.x;
        thePoint->y = (j+0.5)*img.siz.y;
        thePoint->z = (k+0.5)*img.siz.z;
        ind[val] = ind[val]+1;
    }


    fprintf(stderr, "%s : \n", proc);
    for (i=0;i<nClasses;i++) {
        p = &(p_list[i]);
        if (p->n == 0 && ind[i] == 0) continue;
        fprintf(stderr, "     class %d : \t ind = %d \t p->n = %d\n", i, ind[i], p->n);
    }

    free(ind);

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

int VT_WritePointListWithFileName(vt_pointList *p, char *name)
{
    char *proc="VT_WritePointListWithFileName";
    int i;
    FILE *f;
    vt_point *thePoint;
    if ((f = fopen (name, "w")) == NULL) {
      fprintf( stderr, "%s: unable to open '%s' for writing\n", proc, name );
      return( -1 );
    }

    for (i=0 ; i<p->n ; i++)
    {
        thePoint = &(p->list[i]);
        fprintf(f, "%f %f %f %d", thePoint->x, thePoint->y, thePoint->z, thePoint->value);
    }

    fclose(f);

    return(1);
}

int VT_WritePointListsWithFileName(vt_pointList *plist, int nClasses, char *name)
{
    char *proc="VT_WritePointListsWithFileName";
    int i,j;
    FILE *f;
    vt_point *thePoint;
    vt_pointList *p;
    if ((f = fopen (name, "w")) == NULL) {
      fprintf( stderr, "%s: unable to open '%s' for writing\n", proc, name );
      return( -1 );
    }

    /*fprintf(stdout, "nclasses=%d\n", nClasses);*/
    for (j=0 ; j<nClasses ; j++)
    {
        p=&(plist[j]);
        /*fprintf(stdout, "classe[%d]->n=%d\n", j, p->n);*/
        for (i=0; i<p->n ; i++)
        {
          thePoint=&(p->list[i]);
          fprintf(f, "%f %f %f %d\n", thePoint->x, thePoint->y, thePoint->z, thePoint->value);
        }
    }

    fclose(f);

    return(1);
}


int VT_SetPointListIndex(vt_pointList *p, int ind, vt_point q)
{
    if(ind>=p->n)
    {
        return(0);
    }

    vt_point *thePoint = &(p->list[ind]);
    thePoint->weight = q.weight;
    thePoint->value = q.value;
    thePoint->attribute = ind;
    thePoint->x = q.x;
    thePoint->y = q.y;
    thePoint->z = q.z;

    return(1);
}

void VT_CopyPoint(vt_point *point, vt_point q)
{

    point->weight = q.weight;
    point->value = q.value;
    point->attribute = q.attribute;
    point->x = q.x;
    point->y = q.y;
    point->z = q.z;

    return;
}

void printPoint(vt_point point)
{
    fprintf(stdout, "point ( %f %f %f ) -> { value = %d ; weight = %f ; attribute = %d }\n",
            point.x,point.y,point.z,point.value,point.weight,point.attribute);

    return;
}

double pointsDistance2(vt_point p1, vt_point p2)
{
    return((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z));
}

double pointsDistance(vt_point p1, vt_point p2)
{
    return(sqrt(pointsDistance2(p1, p2)));
}


void VT_MultiplyPointListCoordinates(vt_pointList *p, double coef)
{
    int i;
    vt_point *thePoint;
    for (i=0 ; i<p->n ; i++)
    {
        thePoint = &(p->list[i]);
        thePoint->x *= coef;
        thePoint->y *= coef;
        thePoint->z *= coef;
    }
}

void VT_MultiplyPointListsCoordinates(vt_pointList *plist, int nClasses, double coef)
{
    int i,j;
    vt_point *thePoint;
    for (j=0 ; j<nClasses ; j++) {
      for (i=0 ; i<plist[j].n ; i++)
      {
        thePoint = &(plist[j].list[i]);
        thePoint->x *= coef;
        thePoint->y *= coef;
        thePoint->z *= coef;
      }
    }
}

void VT_ApplyTrsfCoef(vt_pointList p_ref, vt_pointList *p_flo, double coef)
{
    /*char *proc="VT_ApplyTrsfCoef";*/
    tree_node *n_ref;
    vt_kdtree(&n_ref, p_ref, 0);
    vt_point *p,q;
    int j;
    double d;
    /*
    vt_coupleList c;
    if (VT_AllocCoupleList(&c, p_flo->n) != 1)
    {
        fprintf(stderr, "%s: unable to allocate couple list\n", proc);
        return;
    }
    for (j = 0 ; j< p_flo.n ; j++)
    {
        p = &(p_flo->list[j]);
        VT_NNS(*p, n_ref, &q, &d);
        VT_CopyPoint(&(c.list[j].start), *p);
        VT_CopyPoint(&(c.list[j].end), q);
    }*/
    for (j=0 ; j < p_flo->n ; j++)
    {
        p=&(p_flo->list[j]);
        VT_NNS(*p, n_ref, &q, &d);
        p->x=coef*p->x+(1-coef)*q.x;
        p->y=coef*p->y+(1-coef)*q.y;
        p->z=coef*p->z+(1-coef)*q.z;
    }

    /*VT_FreeCoupleList(&c);*/
    clearTree(&n_ref);
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

    ind = malloc (nClasses * sizeof(int));
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
                free(ind);
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

    /*fprintf(stderr, "%s: test\n", proc);*/


    fprintf(stderr, "%s : \n", proc);
    for (i=0;i<nClasses;i++) {
        p = &(p_list[i]);
        if (p->n == 0 && ind[i] == 0) continue;
        fprintf(stderr, "     class %d : \t ind = %d \t p->n = %d\n", i, ind[i], p->n);
    }

    free(ind);
    return(1);
}


void VT_FreePointList(vt_pointList *p)
{
    if (p->list)
      free(p->list);
    p=(vt_pointList*) NULL;
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
           /* X axis */
           if (pi->x <= pp->x)
               echanger = 1;
           break;
       case 1 :
           /* Y axis */
           if (pi->y <= pp->y)
               echanger = 1;
           break;
       case 2 :
           /* Z axis */
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


int choix_pivot( int premier, int dernier)
{
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
    /* Rq : le tableau ind doit contenir toutes les valeurs entre 0 et p.n-1, eventuellement permutees
    // tri_rapide retournera le tableau ind avec les indices tels que p.list[ind[:]] soit trie selon l'axe X (axis=0), Y (axis=1) ou Z (axis=2)
    */
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
    /* Rq : le tableau ind doit contenir toutes les valeurs entre 0 et n-1, eventuellement permutees
    // tri_rapide_vector retournera le tableau ind avec les indices tels que v[ind[:]] soit trie par ordre croissant

    //fprintf(stderr, "   tri_rapide_vector : %d, %d\n", premier, dernier);
    */
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

    /* VT_TriRapideVector : Ind pointera sur le tableau d'indices ind tels que v[ind[:]] soit trie par ordre croissant*/
    int i;
    int *ind = malloc(n * sizeof(int));

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

    int m = (int)(p.n/2); /* m = indice de l'elementt median apres tri */
    int i;

    vt_point *point;
    vt_point *cible;

    int *ind = malloc(p.n * sizeof(int));


    if (ind == NULL)
    {
        fprintf(stderr, "%s: unable to allocate index table of length %d\n", proc, p.n);
        return(0);
    }


    for (i=0 ; i<p.n ; i++) ind[i]=i; /* init ind table (faire un vt_triRapide qui engloberait l'initialisation et les appels recursifs?)*/
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
    median->value = point->value;
    median->x = point->x;
    median->y = point->y;
    median->z = point->z;


    for (i = 0 ; i<m ; i++)
    {
        point = &(p.list[ind[i]]);
        cible = &(pL->list[i]);
        cible->attribute=point->attribute;
        cible->weight=point->weight;
        cible->value=point->value;
        cible->x = point->x;
        cible->y = point->y;
        cible->z = point->z;
    }


    for (i = 0 ; i<p.n-1-m ; i++)
    {
        point = &(p.list[ind[i+m+1]]);
        cible = &(pR->list[i]);
        cible->attribute=point->attribute;
        cible->value=point->value;
        cible->weight=point->weight;
        cible->x = point->x;
        cible->y = point->y;
        cible->z = point->z;
    }


    free(ind);

    return(1);
}






void vt_kdtree(tree_node **Node, vt_pointList p, int depth)
{
    char *proc="vt_kdtree";


    tree_node *node;
    vt_point med;
    vt_point *median;
    vt_pointList pL, pR;
    int axis = depth % 3; /* 3-D tree */

    node = malloc(sizeof(tree_node));

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
    median->value = med.value;
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
    nodes = malloc(n * sizeof(tree_node *));
    if (nodes == NULL) {
        fprintf(stderr, "%s : unable to allocate nodees vector\n", proc);
        return;
    }
    for (i=0 ; i<n ; i++) nodes[i]=NULL;
    /*vt_pointList p; */

    int depth = 0;
    /*int axis = depth % 3; // 3-D tree*/


    for (i=0 ; i<n ; i++) {

      /*tree_node *node = nodes[i];*/

      vt_kdtree(&(nodes[i]), pList[i], depth);

      /*vt_point med;
      vt_point *median;
      vt_pointList pL, pR;

      p = pList[i];

      Node = nodeList[i];
      node = malloc(sizeof(tree_node));

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


      // *Node = node;

      // VT_FreePointList(&pL);
      // VT_FreePointList(&pR);
      */
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

    /*fprintf(stdout, "clearTree: ")*/

    if(!tmpP) return;

    if(tmpP->leftChild)  clearTree(&tmpP->leftChild);

    if(tmpP->rightChild) clearTree(&tmpP->rightChild);

    free(tmpP);

    *p = NULL;
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
    if( p.n != p_trsf->n)
    {
        fprintf(stderr, "%s: conflicting field n: original = %d, trsf = %d\n", proc, p.n, p_trsf->n);
        return(0);
    }
    for (i=0; i<p.n; i++)
    {
        theP=&(p.list[i]);
        theTrsf=&(p_trsf->list[i]);
        theTrsf->x=(theP->x*T[0]+theP->y*T[1]+theP->z*T[2]+T[3]);
        theTrsf->y=(theP->x*T[4]+theP->y*T[5]+theP->z*T[6]+T[7]);
        theTrsf->z=(theP->x*T[8]+theP->y*T[9]+theP->z*T[10]+T[11]);
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
    double d = 0.0;
    double *D, dsum;


    D=malloc(n*sizeof(double));
    if(D==NULL) {

        fprintf(stderr, "%s: unable to allocate distance vector\n", proc);
        return(0);
    }
    pairs[0]=malloc(n*sizeof(int));
    if(pairs[0]==NULL) {
        free(D);
        D=NULL;

        fprintf(stderr, "%s: unable to allocate pairs first vector\n", proc);
        return(0);
    }
    pairs[1]=malloc(n*sizeof(int));
    if(pairs[1]==NULL) {
        free(D);
        D=NULL;
        free(pairs[0]);
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
            free(D);
            D=NULL;
            free(pairs[0]);
            pairs[0]=NULL;
            free(pairs[1]);
            pairs[1]=NULL;

            fprintf(stderr, "%s: unexpected result of VT_NNS...\n", proc);
            return(0);
        }
        pairs[0][i]=i;
        pairs[1][i]=j;
        D[i]=d;
    }

    N=n;
    p[0]=malloc(N*sizeof(int));
    if(p[0]==NULL) {
        free(D);
        D=NULL;
        free(pairs[0]);
        pairs[0]=NULL;
        free(pairs[1]);
        pairs[1]=NULL;

        fprintf(stderr, "%s: unable to allocate the p[0] vector...\n", proc);
        return(0);
    }
    p[1]=malloc(N*sizeof(int));
    if(p[1]==NULL) {
        free(D);
        D=NULL;
        free(pairs[0]);
        pairs[0]=NULL;
        free(pairs[1]);
        pairs[1]=NULL;
        free(p[0]);
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

    free(D);
    D=NULL;
    free(pairs[0]);
    pairs[0]=NULL;
    free(pairs[1]);
    pairs[1]=NULL;
    /*free(IND);
    //IND=NULL;
    */
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


    pairs[0]=malloc(n*sizeof(int));
    if(pairs[0]==NULL) {
        fprintf(stderr, "%s: unable to allocate pairs first vector\n", proc);
        return(0);
    }
    pairs[1]=malloc(n*sizeof(int));
    if(pairs[1]==NULL) {
        free(pairs[0]);
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
                free(pairs[0]);
                pairs[0]=NULL;
                free(pairs[1]);
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
    p[0]=malloc(N*sizeof(int));
    if(p[0]==NULL) {
        free(pairs[0]);
        pairs[0]=NULL;
        free(pairs[1]);
        pairs[1]=NULL;

        fprintf(stderr, "%s: unable to allocate the p[0] vector...\n", proc);
        return(0);
    }
    p[1]=malloc(N*sizeof(int));
    if(p[1]==NULL) {
        free(pairs[0]);
        pairs[0]=NULL;
        free(pairs[1]);
        pairs[1]=NULL;
        free(p[0]);
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

    free(pairs[0]);
    pairs[0]=NULL;
    free(pairs[1]);
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



int VT_ResidualsVector(vt_pointList p_in, tree_node *node_ext, double **residuals)
{
    char *proc="VT_ResidualsVector";

    int i;
    vt_point P, Q;

    double d;
    double *r = NULL;

    r = malloc(p_in.n * sizeof(double));

    if (r==NULL)
    {
        fprintf(stderr, "%s: unable to allocate residuals vector\n", proc);
        return(0);
    }


    for (i=0 ; i<p_in.n ; i++)
    {
        Q = p_in.list[i];
        if(VT_NNS(Q, node_ext, &P, &d) != 1)
        {
            free(r); r=NULL;
            fprintf(stderr, "%s: unable to allocate residuals vector\n", proc);
            return(0);
        }
        r[i]=d;
    }

    *residuals = r;
    return(1);
}



void VT_SetResidualImage(vt_image *imResiduals, vt_pointList p_in, double *residuals)
{
    /*int tmp=0;
    */
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
            /*tmp++;
            */
            continue;
        }
        array[k][j][i] = (float)(residuals[ind]);
    }
    /*fprintf(stdout, "tmp = %d\n", tmp);*/
    return;
}

void VT_SetImageFromPointList(vt_image *img, vt_pointList p)
{
    int ind;
    int i,j,k;
    double x,y,z;
    vt_point *P;
    unsigned char ***array;
    switch(img->type ) {
    default:
        fprintf(stderr, "unexpected image parameter type : waiting for UCHAR type\n");
        return;
    case UCHAR:
        array=(unsigned char***)img->array;
    }

    for (ind = 0 ; ind < p.n ; ind++)
    {
        P = &(p.list[ind]);
        x=P->x;
        y=P->y;
        z=P->z;
        i=(int)(x/img->siz.x);
        j=(int)(y/img->siz.y);
        k=(int)(z/img->siz.z);
        if (i<0 || i>=(int)img->dim.x ||
            j<0 || j>=(int)img->dim.y ||
            k<0 || k>=(int)img->dim.z )
        {
            continue;
        }
        array[k][j][i] = (unsigned char)(P->value);
    }
    return;
}


int VT_AllocCoupleList(vt_coupleList *c, int n)
{
    char *proc="VT_AllocCoupleList";
    c->list = NULL;
    c->list = malloc(n*sizeof(vt_couple));
    if (n>0 && c->list == NULL) {
        fprintf(stderr, "%s: unable to allocate a couple list\n",proc);
        return(0);
    }
    c->n = n;
    return(1);
}


int VT_BuildCoupleListFromPointLists(vt_coupleList *c, vt_pointList starts, vt_pointList ends)
{
    char *proc="VT_BuildCoupleListFromPointLists";
    int i;
    vt_couple *couple;
    if (starts.n != ends.n || c->n != starts.n)
    {
        fprintf(stderr, "%s: input lists must be of same size. Aborting...\n", proc);
        return(0);
    }
    for (i = 0 ; i<c->n ; i++)
    {
        couple=&(c->list[i]);
        VT_CopyPoint(&(couple->start), starts.list[i]);
        VT_CopyPoint(&(couple->end), ends.list[i]);
        couple->start.attribute=i;
        couple->end.attribute=i;
    }
    return(1);
}

int VT_AllocCoupleListFromFileName(vt_coupleList *c, char *name)
{
    char *proc="VT_AllocCoupleListWithFileName";

    FILE *f;
    char line[512];
    float f1, f2, f3, f4, f5, f6;
    int o;
    int n = 0;
    int i;

    if ((f = fopen (name, "r")) == NULL) {
      fprintf( stderr, "%s: unable to open '%s' for reading\n", proc, name );
      return( -1 );
    }


    while ( fgets(line, 512, f) != NULL ) {
        o = sscanf( line, "( %f %f %f ) ( %f %f %f )\n", &f1, &f2, &f3, &f4, &f5, &f6);
        if ( o == 6 ) {
            n+=1;
        }
    }

    fclose( f );

    if ( n<0 ) {
      fprintf( stderr, "%s: negative number of points detected in '%s'\n", proc, name );
      return( -1 );
    }

    if (n==0) {
        c->n = 0;
        c->list = (vt_couple*)NULL;
        return(1);
    }

    c->list =malloc(n*sizeof(vt_couple));
    if (c->list == (vt_couple*)NULL) {
        fprintf(stderr, "%s: malloc failed\n", proc);
        return(0);
    }
    c->n = n;

    for (i=0; i<c->n; i++) {
        vt_point *thePoint = &(c->list[i].start);
        thePoint->attribute=0;
        thePoint->value=0;
        thePoint->weight=0;
        thePoint->x=0;
        thePoint->y=0;
        thePoint->z=0;
        thePoint = &(c->list[i].end);
        thePoint->attribute=0;
        thePoint->value=0;
        thePoint->weight=0;
        thePoint->x=0;
        thePoint->y=0;
        thePoint->z=0;
    }

    return(1);
}


int VT_ExtractCoupleListFromFileName(vt_coupleList *c, char *name)
{
    char *proc="VT_ExtractCoupleListFromFileName";
    int ind;
    int o;
    vt_point *thePoint;
    ind=0;

    FILE *f;
    char line[512];
    float f1, f2, f3, f4, f5, f6;

    if ((f = fopen (name, "r")) == NULL) {
      fprintf( stderr, "%s: unable to open '%s' for reading\n", proc, name );
      return( -1 );
    }


    while ( fgets(line, 512, f) != NULL ) {
        o = sscanf( line, "( %f %f %f ) ( %f %f %f )\n", &f1, &f2, &f3, &f4, &f5, &f6);
        if ( o  == 6 ) {
            thePoint=&(c->list[ind].start);
            thePoint->weight = 1; /* default */
            thePoint->value =  255; /* default */
            thePoint->attribute = ind;
            thePoint->x = f1;
            thePoint->y = f2;
            thePoint->z = f3;

            thePoint=&(c->list[ind].end);
            thePoint->weight = 1; /* default */
            thePoint->value =  255; /* default */
            thePoint->attribute = ind;
            thePoint->x = f4;
            thePoint->y = f5;
            thePoint->z = f6;
            ind++;
        }
    }

    fclose( f );

    fprintf(stderr, "%s : ind = %d \t c->n = %d\n", proc, ind, c->n);

    if (ind != c->n) {
        fprintf(stderr, "%s: abnormal number of points found (%d instead of the %d expected)\n", proc, ind, c->n);
        return(0);
    }

    return(1);
}

void printCouple(vt_couple c)
{
    fprintf(stdout, "start = { %f %f %f } \tend = { %f %f %f }\n", c.start.x, c.start.y, c.start.z,
            c.end.x, c.end.y, c.end.z);
}

int VT_WriteCoupleList(vt_coupleList *c, char *name)
{
    char *proc="VT_WriteCoupleList";
    int i;
    FILE *fichier = fopen(name, "w");
    vt_point *start, *end;
    if( fichier == NULL)
    {
        fprintf(stderr, "%s: Erreur pour l'ecriture du fichier %s\n",proc, name);
        return(0);
    }
    for (i = 0 ; i < c->n ; i++)
    {
        start = &(c->list[i].start);
        end = &(c->list[i].end);
        fprintf(fichier, "( %f %f %f ) ( %f %f %f )\n", start->x, start->y, start->z,
                end->x, end->y, end->z);
    }
    return(1);
}

void VT_FreeCoupleList(vt_coupleList *c)
{
    free(c->list);
    c->list=NULL;
    c->n=0;
}


int VT_CoupleListToStartList(vt_coupleList c, vt_pointList *p)
{
    char *proc="VT_CoupleListToStartList";
    int i;
    vt_point *thePoint;
    p->list=malloc(c.n * sizeof(vt_point));
    if(p->list == NULL)
    {
        fprintf(stderr, "%s: unable to allocate the point list\n", proc);
        return(0);
    }
    p->n=c.n;
    for (i = 0 ; i < c.n ; i++)
    {
        thePoint=&(p->list[i]);
        thePoint->value=c.list[i].start.value;
        thePoint->attribute=i;
        thePoint->weight=c.list[i].start.weight;
        thePoint->x=c.list[i].start.x;
        thePoint->y=c.list[i].start.y;
        thePoint->z=c.list[i].start.z;
    }
    return(1);
}

int VT_CoupleListToEndList(vt_coupleList c, vt_pointList *p)
{
    char *proc="VT_CoupleListToEndList";
    int i;
    vt_point *thePoint;
    p->list=malloc(c.n * sizeof(vt_point));
    if(p->list == NULL)
    {
        fprintf(stderr, "%s: unable to allocate the point list\n", proc);
        return(0);
    }
    p->n=c.n;
    for (i = 0 ; i < c.n ; i++)
    {
        thePoint=&(p->list[i]);
        thePoint->value=c.list[i].end.value;
        thePoint->attribute=i;
        thePoint->weight=c.list[i].end.weight;
        thePoint->x=c.list[i].end.x;
        thePoint->y=c.list[i].end.y;
        thePoint->z=c.list[i].end.z;
    }
    return(1);
}

int VT_CoupleListToMedianList(vt_coupleList c, vt_pointList *p)
{
    char *proc="VT_CoupleListToMedianList";
    int i;
    p->list=malloc(c.n * sizeof(vt_point));
    if(p->list == NULL)
    {
        fprintf(stderr, "%s: unable to allocate the point list\n", proc);
        return(0);
    }
    p->n=c.n;
    vt_point *thePoint;
    for (i = 0 ; i < c.n ; i++)
    {
        thePoint=&(p->list[i]);
        thePoint->value=c.list[i].start.value;
        thePoint->attribute=i;
        thePoint->weight=(c.list[i].start.weight+c.list[i].end.weight)*0.5;
        thePoint->x=(c.list[i].start.x+c.list[i].end.x)*0.5;
        thePoint->y=(c.list[i].start.y+c.list[i].end.y)*0.5;
        thePoint->z=(c.list[i].start.z+c.list[i].end.z)*0.5;
    }
    return(1);
}

double VT_CouplesScalarProduct(vt_couple c1, vt_couple c2)
{
    double x1=c1.end.x-c1.start.x;
    double x2=c2.end.x-c2.start.x;
    double y1=c1.end.y-c1.start.y;
    double y2=c2.end.y-c2.start.y;
    double z1=c1.end.z-c1.start.z;
    double z2=c2.end.z-c2.start.z;
    return(x1*x2+y1*y2+z1*z2);
}

double VT_CoupleNorm2(vt_couple c)
{
    return((c.end.x-c.start.x)*(c.end.x-c.start.x)+(c.end.y-c.start.y)*(c.end.y-c.start.y)+(c.end.z-c.start.z)*(c.end.z-c.start.z));
}

double VT_CoupleNorm(vt_couple c)
{
    return(sqrt(VT_CoupleNorm2(c)));
}

double VT_CouplesAngle(vt_couple c1, vt_couple c2)
{
    double normsProduct=sqrt(VT_CoupleNorm2(c1)*VT_CoupleNorm2(c2));
    double cos=VT_CouplesScalarProduct(c1,c2)/normsProduct;
    return(acos(cos));
}

int VT_ExtractBijectiveCoupleList(vt_coupleList *c, vt_coupleList c1, vt_coupleList c2, double epsilon2)
{
    char *proc="VT_ExtractBijectiveCoupleList";
    int i, j, ind=0;
    vt_point *pt1, pt2, *pt3;
    vt_pointList p2;
    tree_node *n2;
    double d = 0.0;

    if (VT_AllocCoupleList(c, c1.n) != 1)
    {
        fprintf(stderr, "%s: unable to allocate a couple list\n", proc);
        return(0);
    }

    if (VT_CoupleListToStartList(c2, &p2) != 1)
    {
        fprintf(stderr, "%s: unable to extract the start list from a couple list\n", proc);
        VT_FreeCoupleList(c);
        return(0);
    }

    vt_kdtree(&n2, p2, 0);

    for (i=0 ; i<c->n ; i++)
    {
        pt1=&(c1.list[i].start);
        VT_NNS(c1.list[i].end, n2, &pt2, &d);
        /*if (i<10)
        {
            fprintf(stdout, "pt1: \n");
            printPoint(*pt1);
            fprintf(stdout, "pt1.end: \n");
            printPoint(c1.list[i].end);
            fprintf(stdout, "pt2 (d = %f): \n", d);
            printPoint(pt2);
        }*/
        if (d > ((epsilon2<1e-3) ? epsilon2 : 1e-3))
        {
            /*if (i<10)
            //    fprintf(stdout, "FUCK %f\n", (epsilon2<1e-3) ? epsilon2 : 1e-3);
            */
            continue;
        }
        j=pt2.attribute;
        pt3=&(c2.list[j].end);

        /*if (i<10)
        {
            fprintf(stdout, "pt3: \n");
            printPoint(*pt3);
            fprintf(stdout, "\n");
        }*/
        if(pointsDistance2(*pt1, *pt3) > epsilon2)
        {
            continue;
        }
        VT_SetCoupleListIndex(c, ind++, *pt1, pt2);
    }
    c->n = ind;
    fprintf(stdout,"%s: ind = %d\n",proc, ind);
    VT_FreePointList(&p2);
    clearTree(&n2);
    return(1);
}

int VT_SetCoupleListIndex(vt_coupleList *c, int ind, vt_point p, vt_point q)
{
    vt_point *theStart = &(c->list[ind].start);
    vt_point *theEnd = &(c->list[ind].end);

    VT_CopyPoint(theStart, p);
    VT_CopyPoint(theEnd, q);

    theStart->attribute=ind;
    theEnd->attribute=ind;

    return(1);
}

int VT_SetImageFromCoupleList(vt_image *img, vt_coupleList c)
{
    char *proc="VT_SetImageFromCoupleList";
    int i, j, k, ind;
    double x, y, z;
    vt_point *P;
    unsigned char ***arrayU8=(unsigned char ***)img->array;
    unsigned short int ***arrayU16=(unsigned short int ***)img->array;
    float ***array32=(float ***)img->array;

    for (ind = 0 ; ind < c.n ; ind++)
    {
        /* Start */
        P = &(c.list[ind].start);
        x=P->x;
        y=P->y;
        z=P->z;
        i=(int)(x/img->siz.x);
        j=(int)(y/img->siz.y);
        k=(int)(z/img->siz.z);
        if (ind < 0) {
            fprintf(stdout, "ind %d : ", ind);
            printCouple(c.list[ind]);
            fprintf(stdout, "[i, j, k, value] ----> [ %d %d %d %d]\n\n",i,j,k,P->value);
        }
        if (i>=0 && i<(int)img->dim.x &&
            j>=0 && j<(int)img->dim.y &&
            k>=0 && k<(int)img->dim.z )
        {
            switch (img->type) {
            default:
                fprintf(stderr, "%s: unexpected image type, aborting.\n", proc);
                return(-1);
            case UCHAR:
            case SCHAR:
                arrayU8[k][j][i] = (unsigned char)(P->value);
                break;
            case USHORT:
            case SSHORT:
                arrayU16[k][j][i] = (unsigned short int)(P->value);
                break;
            case FLOAT:
                array32[k][j][i] = (float)(P->value);
                break;
            }
        }
        /* End */
        P = &(c.list[ind].end);
        x=P->x;
        y=P->y;
        z=P->z;
        i=(int)(x/img->siz.x);
        j=(int)(y/img->siz.y);
        k=(int)(z/img->siz.z);
        if (i>=0 && i<(int)img->dim.x &&
            j>=0 && j<(int)img->dim.y &&
            k>=0 && k<(int)img->dim.z )
        {
            switch (img->type) {
            default:
                fprintf(stderr, "%s: unexpected image type, aborting.\n", proc);
                return(-1);
            case UCHAR:
            case SCHAR:
                arrayU8[k][j][i] = (unsigned char)(P->value);
                break;
            case USHORT:
            case SSHORT:
                arrayU16[k][j][i] = (unsigned short int)(P->value);
                break;
            case FLOAT:
                array32[k][j][i] = (float)(P->value);
                break;
            }
        }
    }

    return(1);
}

int VT_WriteCoupleListDistanceHistogram(vt_coupleList c, char *name)
{
    char *proc="VT_WriteCoupleListDistanceHistogram";
    int i;
    FILE *f=fopen(name, "w");
    if (f == NULL)
    {
        fprintf(stderr, "%s: unable to create output file '%s'\n", proc, name);
        return(0);
    }
    for (i = 0 ; i < c.n ; i++)
    {
        fprintf(f, "%f\n", VT_CoupleNorm(c.list[i]));
    }
    fclose(f);
    return(1);
}
