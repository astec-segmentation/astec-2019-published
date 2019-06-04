/*
 * vt_planeRegistration.c -
 *
 * $Id: vt_planeRegistration.c,v 1.0 2014/07/23 15:17:34 gael Exp $
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



#include <convert.h>
#include <vtmalloc.h>

#include <vt_common.h>

/* #include <bal-estimator.h> */
#include <bal-field.h>
#include <bal-field-tools.h>
#include <bal-point.h>
#include <bal-transformation.h>
#include <bal-transformation-tools.h>

#include <vt_tree.h> /* required for VT_NPoints() definition */

#include <vt_planeRegistration.h>



static int _verbose_ = 0;


/* static double EPS=1e-16; */



int VT_Barycentre(vt_image* gradient, double* barycentre)
{
  char *proc="VT_Barycentre";
  size_t i,j,k;
  double weight;
  double sum=0;

  unsigned char ***arrayu8=NULL;
  unsigned short int ***arrayu16=NULL;
  float ***array32=NULL;

  switch (gradient->type) {
  default:
      if (_verbose_)
          fprintf(stderr, "%s: such image type not handled yet\n", proc);
      return(0);
  case UCHAR:
  case SCHAR:
      arrayu8 = (unsigned char ***) gradient->array;
      break;
  case USHORT:
  case SSHORT:
      arrayu16 = (unsigned short int ***) gradient->array;
      break;
  case FLOAT:
      array32 = (float ***) gradient->array;
      break;
  }

  barycentre[0]=0;
  barycentre[1]=0;
  barycentre[2]=0;

  for (k=0 ; k<gradient->dim.z ; k++)
  for (j=0 ; j<gradient->dim.y ; j++)
  for (i=0 ; i<gradient->dim.x ; i++) {
    switch (gradient->type) {
    default:
        if (_verbose_)
            fprintf(stderr, "%s: such image type not handled yet\n", proc);
        return(0);
    case UCHAR:
    case SCHAR:
        weight = (double)arrayu8[k][j][i];
        break;
    case USHORT:
    case SSHORT:
        weight = (double)arrayu16[k][j][i];
        break;
    case FLOAT:
        weight = (double)array32[k][j][i];
        break;
    }

    barycentre[0] += weight * i * gradient->siz.x;
    barycentre[1] += weight * j * gradient->siz.y;
    barycentre[2] += weight * k * gradient->siz.z;
    sum += weight;
  }
  if (sum != 0) {
    barycentre[0]/=sum;
    barycentre[1]/=sum;
    barycentre[2]/=sum;
  }


  return(1);
}

void VT_BaryImgLabels(vt_image_barycentres b, double* G)
{
    double w=0;
    int i;
    vt_barycentre theB;

    G[0]=0;G[1]=0;G[2]=0;
    for (i=0; i<b.n; i++)
    {
        theB=b.barycentres[i];
        w+=theB.weight;
        G[0]+=theB.weight*theB.x;
        G[1]+=theB.weight*theB.y;
        G[2]+=theB.weight*theB.z;
    }
    if (w!= 0) {
      G[0]*=b.vx/w;
      G[1]*=b.vy/w;
      G[2]*=b.vz/w;
    }
}


/*  Computes rotation matrix new->old coordinates to align plane normal vector to X axis */
void VT_RotationMatrix(double *n_new, double* n, double* R)
{
    double ux, uy, uz, c, s;
    ux=n[1]*n_new[2]-n[2]*n_new[1];
    uy=-n[0]*n_new[2]+n[2]*n_new[0];
    uz=n[0]*n_new[1]-n[1]*n_new[0]; /*  vecteur u = vecteur de rotation */
    s=sqrt(ux*ux+uy*uy+uz*uz);
    if(s<1e-6)
    {
      R[0]=R[4]=R[8]=1.0;
      R[1]=R[2]=R[3]=R[5]=R[6]=R[7]=0.0;
    }
    else {
      ux/=s;
      uy/=s;
      uz/=s;
      c=-(n[0]*n_new[0]+n[1]*n_new[1]+n[2]*n_new[2])/(sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])*sqrt(n_new[0]*n_new[0]+n_new[1]*n_new[1]+n_new[2]*n_new[2]));
      s=sqrt(1-c*c);
      /*
      R=[                                                       ...
      ux^2+(1-ux^2)*c   ux*uy*(1-c)-uz*s    ux*uz*(1-c)+uy*s  ; ...
      ux*uy*(1-c)+uz*s  uy^2+(1-uy^2)*c     uy*uz*(1-c)-ux*s  ; ...
      ux*uz*(1-c)-uy*s  uy*uz*(1-c)+ux*s    uz^2+(1-uz^2)*c   ; ...
      ];
      */
      R[0]=ux*ux+(1-ux*ux)*c;
      R[1]=ux*uy*(1-c)-uz*s;
      R[2]=ux*uz*(1-c)+uy*s;
      R[3]=ux*uy*(1-c)+uz*s;
      R[4]=uy*uy+(1-uy*uy)*c;
      R[5]=uy*uz*(1-c)-ux*s;
      R[6]=ux*uz*(1-c)-uy*s;
      R[7]=uy*uz*(1-c)+ux*s;
      R[8]=uz*uz+(1-uz*uz)*c;
    }
}

void VT_RotationMatrixWithAngle(double *u, double theta, double* R)
{
    double ux, uy, uz, c, s;
    ux=u[0];
    uy=u[1];
    uz=u[2];    /* u : vecteur de rotation */
    s=sqrt(ux*ux+uy*uy+uz*uz);
    if(s<1e-6)
    {
      R[0]=R[4]=R[8]=1.0;
      R[1]=R[2]=R[3]=R[5]=R[6]=R[7]=0.0;
    }
    else {
      ux/=s;
      uy/=s;
      uz/=s;
      c=cos(theta);
      s=sin(theta);
      /*
      R=[                                                       ...
      ux^2+(1-ux^2)*c   ux*uy*(1-c)-uz*s    ux*uz*(1-c)+uy*s  ; ...
      ux*uy*(1-c)+uz*s  uy^2+(1-uy^2)*c     uy*uz*(1-c)-ux*s  ; ...
      ux*uz*(1-c)-uy*s  uy*uz*(1-c)+ux*s    uz^2+(1-uz^2)*c   ; ...
      ];
      */
      R[0]=ux*ux+(1-ux*ux)*c;
      R[1]=ux*uy*(1-c)-uz*s;
      R[2]=ux*uz*(1-c)+uy*s;
      R[3]=ux*uy*(1-c)+uz*s;
      R[4]=uy*uy+(1-uy*uy)*c;
      R[5]=uy*uz*(1-c)-ux*s;
      R[6]=ux*uz*(1-c)-uy*s;
      R[7]=uy*uz*(1-c)+ux*s;
      R[8]=uz*uz+(1-uz*uz)*c;
    }
}

void VT_RotationMatrix180(double *n, double *R)
{
  double ux, uy, uz,s ;/* , c, s; */
    int imin=0;
    if(fabs(n[0])>fabs(n[1])) imin=1;
    if(fabs(n[imin])>fabs(n[2])) imin=2;
    switch(imin) {
    default:
    case 0:
        ux=0;
        uy=-n[2];
        uz=n[1];
        break;
    case 1:
        ux=-n[2];
        uy=0;
        uz=n[0];
        break;
    case 2:
        ux=n[1];
        uy=-n[0];
        uz=0;
        break;
    }
    s=sqrt(ux*ux+uy*uy+uz*uz);
    if(s<1e-6)
    {
      R[0]=R[4]=R[8]=-1.0;
      R[1]=R[2]=R[3]=R[5]=R[6]=R[7]=0.0;
    }
    else {
      ux/=s;
      uy/=s;
      uz/=s;
      /* c=-1; */
      /* s=0; */
      /*
      R=[                                                       ...
      ux^2+(1-ux^2)*c   ux*uy*(1-c)-uz*s    ux*uz*(1-c)+uy*s  ; ...
      ux*uy*(1-c)+uz*s  uy^2+(1-uy^2)*c     uy*uz*(1-c)-ux*s  ; ...
      ux*uz*(1-c)-uy*s  uy*uz*(1-c)+ux*s    uz^2+(1-uz^2)*c   ; ...
      ];
      */
      R[0]=ux*ux-(1-ux*ux);
      R[1]=ux*uy*2;
      R[2]=ux*uz*2;
      R[3]=ux*uy*2;
      R[4]=uy*uy-(1-uy*uy);
      R[5]=uy*uz*2;
      R[6]=ux*uz*2;
      R[7]=uy*uz*2;
      R[8]=uz*uz-(1-uz*uz);
    }
}

void VT_Translation(double *G_new, double *G, double *R, double *t)
{
    t[0]=G[0]-(R[0]*G_new[0]+R[1]*G_new[1]+R[2]*G_new[2]);
    t[1]=G[1]-(R[3]*G_new[0]+R[4]*G_new[1]+R[5]*G_new[2]);
    t[2]=G[2]-(R[6]*G_new[0]+R[7]*G_new[1]+R[8]*G_new[2]);
}

int VT_FindLabel(vt_image_barycentres b, int val)
{
    int i;
    vt_barycentre theBary;
    for (i=0 ; i<b.n ; i++)
    {
        theBary=b.barycentres[i];
        if (val==theBary.label)
            return(i);
    }
    return (-1);
}





int VT_NCell(vt_image img, int background)
{
    char *proc="VT_NCell";

    unsigned char *bufu8=NULL;
    unsigned short int *bufu16=NULL;

    size_t i;
    int j, val, valold;
    int n=0;
    int *labels;
    int delta=256;

    int add;

    switch (img.type) {
    default:
        if(1 || _verbose_)
            fprintf(stderr, "%s: such image type not handled yet\n", proc);
        return(-1);
    case UCHAR:
    case SCHAR:
        bufu8=(unsigned char *)img.buf;
        break;
    case USHORT:
    case SSHORT:
        bufu16=(unsigned short int *)img.buf;
        break;
    }



    labels = vtmalloc(delta*sizeof(int), "labels", proc );
    if(labels==(int*)NULL)
    {
        if(1 || _verbose_)
            fprintf(stderr, "%s: problem while allocating label vector\n", proc);
        return(-1);
    }


    valold=background;
    for (i=0 ; i<img.dim.x*img.dim.y*img.dim.z ; i++) {
        switch (img.type){
        default:
            if(1 || _verbose_)
                fprintf(stderr, "%s: such image type not handled yet\n", proc);
            vtfree(labels);
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

        if (val==valold || val==background) continue;
        valold=val;

        add=1;
        for (j=0 ; j<n ; j++) {
            if (val!=labels[j]) continue;
            add=0;
            break;
        }
        if (add==0) continue;

        if (n%delta == 0)
        {
            labels= vtrealloc( labels, (delta+n)*sizeof(int), "labels", proc );
            if(labels==(int*)NULL)
            {
                if(1 || _verbose_)
                    fprintf(stderr, "%s: problem while re-allocating label vector\n", proc);
                return(-1);
            }
        }
        labels[n++]=val;
    }

    vtfree(labels);
    return (n);
}

int VT_ComputeLabelBarycentres(vt_image labels, vt_image_barycentres* b, int background)
{
    char *proc="VT_ComputeLabelBarycentres";

    unsigned char ***arrayu8=NULL;
    unsigned short int ***arrayu16=NULL;

    size_t i, j, k;
    int val, valold;
    int ind;

    vt_barycentre *theBary;

    switch (labels.type) {
    default:
        if(_verbose_)
            fprintf(stderr, "%s: such image type not handled yet\n", proc);
        return(-1);
    case UCHAR:
    case SCHAR:
        arrayu8=(unsigned char ***)labels.array;
        break;
    case USHORT:
    case SSHORT:
        arrayu16=(unsigned short int ***)labels.array;
        break;
    }

    valold=background;
    ind=0;
    for (k=0 ; k<labels.dim.z ; k++)
    for (j=0 ; j<labels.dim.y ; j++)
    for (i=0 ; i<labels.dim.x ; i++) {
        switch (labels.type){
        default:
            if(_verbose_)
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
        if(val==background) continue;
        if(val!=valold) {
            ind=VT_FindLabel(*b, val);
            if(ind<0)
            {
                if(_verbose_)
                    fprintf(stderr, "%s: did not find the index of label # %d\n", proc, val);
                return(-1);
            }
        }
        valold=val;
        theBary=&(b->barycentres[ind]);
        if(theBary->label != val) {
            if(_verbose_)
                fprintf(stderr, "%s: unexpected returned index # %d for value %d\n", proc, ind, val);
            return(-1);
        }
        theBary->weight += 1;
        theBary->x += i;
        theBary->y += j;
        theBary->z += k;
    }

    for (ind=0 ; ind<b->n ; ind++) {
        theBary=&(b->barycentres[ind]);
        if (theBary->weight==0) continue;
	/*         theBary->x *= labels.siz.x/theBary->weight; */
	/*         theBary->y *= labels.siz.y/theBary->weight; */
	/*         theBary->z *= labels.siz.z/theBary->weight; */
        theBary->x /= theBary->weight;
        theBary->y /= theBary->weight;
        theBary->z /= theBary->weight;
        /*theBary->weight *= b->vx * b->vy * b->vz;*/
    }

    return(1);
}


int VT_ExtractImagePoints(vt_image img, vt_image_barycentres* b)
{
    char *proc="VT_ExtractImagePoints";

    unsigned char ***arrayu8=NULL;
    unsigned short int ***arrayu16=NULL;

    size_t i, j, k;
    int val;
    int ind;

    vt_barycentre *theBary;

    switch (img.type) {
    default:
        if(_verbose_)
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
    for (k=0 ; k<img.dim.z ; k++)
    for (j=0 ; j<img.dim.y ; j++)
    for (i=0 ; i<img.dim.x ; i++) {
        switch (img.type){
        default:
            if(_verbose_)
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

        theBary=&(b->barycentres[ind++]);

        theBary->weight += 1;
        theBary->x += i;
        theBary->y += j;
        theBary->z += k;
    }

    fprintf(stderr, "%s : ind = %d \t b->n = %d\n", proc, ind, b->n);

    for (ind=0 ; ind<b->n ; ind++) {
        theBary=&(b->barycentres[ind]);
        if (theBary->weight==0) continue;
        theBary->x /= theBary->weight;
        theBary->y /= theBary->weight;
        theBary->z /= theBary->weight;
        theBary->weight *= b->vx * b->vy * b->vz;
    }

    return(1);
}


void VT_PrintImageBarycentres(vt_image_barycentres b)
{
    int i;
    vt_barycentre theBary;
    fprintf(stdout, "Number of labels : %d\n", b.n);
    fprintf(stdout, "Voxelsize : %f %f %f\n", b.vx, b.vy, b.vz);
    if(b.n>0)
    {
        theBary=b.barycentres[0];
        fprintf(stdout, "Gravity centers  : Pt=[ %f  %f  %f ]\t Label=%d\t Weight=%d\n", theBary.x, theBary.y, theBary.z, theBary.label, (int)theBary.weight);
        for (i=1; i<b.n; i++)
        {
            theBary=b.barycentres[i];
            fprintf(stdout, "                    Pt=[ %f  %f  %f ]\t Label=%d\t Weight=%d\n", theBary.x, theBary.y, theBary.z, theBary.label, (int)theBary.weight);
        }
    }
}

static int VT_InitImageBarycentresWithList(vt_image_barycentres *b, char *filename, int background)
{
    char *proc="VT_InitImageBarycentresWithList";

    int n=0;

    int label, out;
    float x,y,z, vol;
    float vx, vy, vz;
    char keyword[512];
    b->vx=1;
    b->vy=1;
    b->vz=1;
    vt_barycentre *theBary;
    FILE *f;
    char line[512];



    if ((f = fopen (filename, "r")) == NULL) {

        fprintf(stderr, "%s: file %s not existing...", proc, filename);
        return (0);
    }

    if ( fgets(line, 512, f) == NULL ) {
        fprintf(stderr,  "%s: an error occured while reading file.\n", proc);
    }
    out=sscanf( line, "%d %f %f %f %f\n", &label, &x, &y, &z, &vol );
    while ( ( out < 4 )  ) {
      if ( sscanf(line, "%s %f %f %f\n", keyword, &vx, &vy, &vz) == 4 && (strcmp(keyword,"voxelsize") == 0 || strcmp(keyword,"Voxelsize") == 0)) {
          b->vx=vx;
          b->vy=vy;
          b->vz=vz;
          if (_verbose_)
            fprintf(stdout, "%s: Barycentres coordinates under voxel size of ( %f %f %f )\n", proc, b->vx, b->vy, b->vz);
      }
      if ( fgets(line, 512, f) == NULL ) {
          fprintf( stderr, "%s: an error occured while reading file %s.\n",proc, filename);
          return(0);
      }
      out=sscanf( line, "%d %f %f %f %f\n", &label, &x, &y, &z, &vol );
    }
    if (label != background) {
        theBary = &(b->barycentres[n++]);
        theBary->label=label;
        theBary->x=x;
        theBary->y=y;
        theBary->z=z;
        if(out==5)
            theBary->weight=vol;
        else
            theBary->weight=0;
    }
    while ( n<b->n && fgets(line, 512, f) != NULL ) {
        out=sscanf( line, "%d %f %f %f %f\n", &label, &x, &y, &z, &vol );
        if ( out == 4 || out == 5 ) {
            if (label != background) {
                theBary = &(b->barycentres[n++]);
                theBary->label=label;
                theBary->x=x;
                theBary->y=y;
                theBary->z=z;
                if (out==5) {
                    theBary->weight=vol;
                }
                else {
                    theBary->weight=0;
                }
            }
        }
    }

    fclose( f );

    if(n != b->n)
    {
        if (1 || _verbose_)
            fprintf(stderr, "%s: incompatible vt_barycentre_image field n = %d and file number of cells = %d\n", proc, b->n, n);
        return(0);
    }
    return (1);
}





static int VT_InitImageBarycentres(vt_image_barycentres* b, vt_image image, int background)
{
    char *proc="VT_InitImageBarycentres";
    unsigned char *bufu8=NULL;
    unsigned short int *bufu16=NULL;

    size_t i;
    int j, val, valold;
    int n=0;
    int *labels;
    int delta=256;

    int add;

    b->vx=image.siz.x;
    b->vy=image.siz.y;
    b->vz=image.siz.z;

    switch (image.type) {
    default:
        if(1 || _verbose_)
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



    labels = vtmalloc(delta*sizeof(int), "labels", proc );
    if(labels==(int*)NULL)
    {
        if(1 || _verbose_)
            fprintf(stderr, "%s: problem while allocating label vector\n", proc);
        return(-1);
    }

    valold=background;
    for (i=0 ; i<image.dim.x*image.dim.y*image.dim.z ; i++) {
        switch (image.type){
        default:
            if(1 || _verbose_)
                fprintf(stderr, "%s: such image type not handled yet\n", proc);
            vtfree(labels);
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

        if (val==valold || val==background) continue;
        valold=val;

        add=1;
        for (j=0 ; j<n ; j++) {
            if (val!=labels[j]) continue;
            add=0;
            break;
        }
        if (add==0) continue;

        if (n%delta == 0)
        {
            labels= vtrealloc(labels, (delta+n)*sizeof(int), "labels", proc );
            if(labels==(int*)NULL)
            {
                if(1 || _verbose_)
                    fprintf(stderr, "%s: problem while re-allocating label vector\n", proc);
                return(-1);
            }
        }
        labels[n++]=val;
    }
    if(n != b->n)
    {
        if (1 || _verbose_)
            fprintf(stderr, "%s: incompatible vt_barycentre_image field n = %d and vt_image number of cells = %d\n", proc, b->n, n);
        vtfree(labels);
        return(0);
    }

    /*  To improve */
    for (i=0; i<(size_t)b->n; i++) {
        vt_barycentre *theBary = &(b->barycentres[i]);
        theBary->label=labels[i];
        theBary->weight=0;
        theBary->x=0;
        theBary->y=0;
        theBary->z=0;
    }

    vtfree(labels);
    return (1);
}

static int VT_InitImagePoints(vt_image_barycentres* b, vt_image image)
{
    char *proc="VT_InitImagePoints";

    unsigned char *bufu8=NULL;
    unsigned short int *bufu16=NULL;

    size_t i;
    int val;
    int n=0;

    b->vx=image.siz.x;
    b->vy=image.siz.y;
    b->vz=image.siz.z;

    switch (image.type) {
    default:
        if(1 || _verbose_)
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
            if(1 || _verbose_)
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

    if(n != b->n)
    {
        if (1 || _verbose_)
            fprintf(stderr, "%s: incompatible vt_barycentre_image field n = %d and vt_image number of points = %d\n", proc, b->n, n);
        return(0);
    }

    for (i=0; i<(size_t)b->n; i++) {
        vt_barycentre *theBary = &(b->barycentres[i]);
        theBary->label=i;
        theBary->weight=0;
        theBary->x=0;
        theBary->y=0;
        theBary->z=0;
    }

    return (1);
}


int VT_TransformImageBarycentres(double *T, vt_image_barycentres b, vt_image_barycentres *b_trsf)
{
    char *proc="VT_TransformImageBarycentres";
    int i;
    vt_barycentre *theB, *theTrsf;
    if( b.n != b_trsf->n)
    {
        if(_verbose_)
            fprintf(stderr, "%s: conflicting field n: original = %d, trsf = %d\n", proc, b.n, b_trsf->n);
        return(0);
    }
    for (i=0; i<b.n; i++)
    {
        theB=&(b.barycentres[i]);
        theTrsf=&(b_trsf->barycentres[i]);
        theTrsf->x=(theB->x*b.vx*T[0]+theB->y*b.vy*T[1]+theB->z*b.vz*T[2]+T[3])/b.vx;
        theTrsf->y=(theB->x*b.vx*T[4]+theB->y*b.vy*T[5]+theB->z*b.vz*T[6]+T[7])/b.vy;
        theTrsf->z=(theB->x*b.vx*T[8]+theB->y*b.vy*T[9]+theB->z*b.vz*T[10]+T[11])/b.vz;
        theTrsf->label=theB->label;
        theTrsf->weight=theB->weight;
    }
    b_trsf->vx=b.vx;
    b_trsf->vy=b.vy;
    b_trsf->vz=b.vz;

    return(1);
}

int VT_CopyImageBarycentres(vt_image_barycentres b, vt_image_barycentres *b_copy)
{
    char *proc="VT_CopyImageBarycentres";
    int i;
    vt_barycentre *theB, *theCopy;
    if( b.n != b_copy->n)
    {
        if(_verbose_)
            fprintf(stderr, "%s: conflicting field n: original = %d, copy = %d\n", proc, b.n, b_copy->n);
        return(0);
    }
    for (i=0; i<b.n; i++)
    {
        theB=&(b.barycentres[i]);
        theCopy=&(b_copy->barycentres[i]);
        theCopy->x=theB->x;
        theCopy->y=theB->y;
        theCopy->z=theB->z;
        theCopy->label=theB->label;
        theCopy->weight=theB->weight;
    }
    b_copy->vx=b.vx;
    b_copy->vy=b.vy;
    b_copy->vz=b.vz;
    return(1);
}





int VT_AllocImageBarycentresWithImage(vt_image_barycentres* b, vt_image labels, int background)
{
    char *proc="VT_AllocImageBarycentresWithImage";
    int n=VT_NCell(labels, background);
    if (n<=0) {
        fprintf(stderr, "%s: negative number of cells\n", proc);
        return(0);
    }
    b->barycentres = vtmalloc( n*sizeof(vt_barycentre), "b->barycentres", proc );
    if (b->barycentres == (vt_barycentre*)NULL) {
        fprintf(stderr, "%s: allocation failed\n", proc);
        return(0);
    }
    b->n = n;
    return(VT_InitImageBarycentres(b, labels, background));
}





int VT_SetImageBarycentresWithList(vt_image_barycentres* b, char *filename, int background)
{
    char *proc="VT_SetImageBarycentresWithList";
    int n=0;
    int out;
    FILE *f;
    if ((f = fopen (filename, "r")) != NULL) {
      /* File */
        char line[512];
        int label;
        float x, y, z, vol;
        while ( fgets(line, 512, f) != NULL ) {
            /*fprintf(stdout, "%s: line %d : %s", filename, n, line);*/
            out=sscanf( line, "%d %f %f %f %f\n", &label, &x, &y, &z, &vol);
            if ( (out == 4 || out == 5 ) ) {
                if (label != background)
                    n+=1;
            }
        }

        if ( n==0 ) {
            fprintf(stderr,  "%s: Found 0 line matching with the format %%d %%f %%f %%f [%%f] in file %s.\n", proc, filename);
            return(0);
        }
        if (n >= TABLENGTH)
        {
            fprintf(stderr,  "%s: Found %d line matching with the format %%d %%f %%f %%f [%%f] in file %s (limit is %d).\n", proc, n, filename, TABLENGTH);
            return(0);
        }
        fclose( f );
    }
    else {
        fprintf(stderr,  "%s: unable to read input label file %s.\n", proc, filename);
        return(0);
    }
    if (n<=0) {
        fprintf(stderr, "%s: negative number of cells\n", proc);
        return(0);
    }
    b->barycentres = vtmalloc( n*sizeof(vt_barycentre), "b->barycentres", proc );
    if (b->barycentres == (vt_barycentre*)NULL) {
        fprintf(stderr, "%s: allocation failed\n", proc);
        return(0);
    }
    b->n = n;
    return(VT_InitImageBarycentresWithList(b, filename, background));
}





int VT_AllocImagePointsWithImage(vt_image_barycentres* b, vt_image img)
{
    char *proc="VT_AllocImagePointsWithImage";
    int n=VT_NPoints(&img);
    if (n<=0) {
        fprintf(stderr, "%s: negative number of points\n", proc);
        return(0);
    }
    b->barycentres = vtmalloc( n*sizeof(vt_barycentre), "b->barycentres", proc );
    if (b->barycentres == (vt_barycentre*)NULL) {
        fprintf(stderr, "%s: allocation failed\n", proc);
        return(0);
    }
    b->n = n;
    return(VT_InitImagePoints(b, img));
}





int VT_AllocImageBarycentres(vt_image_barycentres* b, int n)
{
    char *proc = "VT_AllocImageBarycentres";
    if (n<=0) return(0);
    b->barycentres = vtmalloc( n*sizeof(vt_barycentre), "b->barycentres", proc );
    if (b->barycentres == (vt_barycentre*)NULL)
        return(0);
    b->n = n;
    return(1);
}




void VT_FreeImageBarycentres(vt_image_barycentres* b)
{
    vtfree(b->barycentres);
    b=(vt_image_barycentres*) NULL;
}

void VT_ImageBarycentresCentralLayersExtraction(vt_image_barycentres* b, double* n, double delta)
{
    vt_barycentre *theB, *theTmp;
    double d;
    int N;
    int i;
    N=0;
    for (i=0; i<b->n; i++) {
        theB=&(b->barycentres[i]);
        d=n[0]*theB->x*b->vx+n[1]*theB->y*b->vy+n[2]*theB->z*b->vz+n[3];
        if( fabs(d)<delta) {
            theTmp=&(b->barycentres[N++]);
            theTmp->x=theB->x;
            theTmp->y=theB->y;
            theTmp->z=theB->z;
            theTmp->weight=theB->weight;
            theTmp->label=theB->label;
        }
    }
    b->n=N;
}

static int VT_FindCloser(vt_barycentre bary, double *siz, vt_image_barycentres b, double *d)
{
    vt_barycentre theB;
    double distMin=-1;
    double dist2;
    int iMin=0, i;

    for (i=0 ; i<b.n ; i++)
    {
        theB=b.barycentres[i];
        dist2=(bary.x*siz[0]-theB.x*b.vx)*(bary.x*siz[0]-theB.x*b.vx)+
              (bary.y*siz[1]-theB.y*b.vy)*(bary.y*siz[1]-theB.y*b.vy)+
              (bary.z*siz[2]-theB.z*b.vz)*(bary.z*siz[2]-theB.z*b.vz);
        if (dist2<distMin || distMin<0)
        {
            distMin=dist2;
            iMin=i;
        }
    }
    *d=sqrt(distMin);
    return(iMin);
}


/* VT_ComputePairsOfPoints :
    calcul des appariements de labels selon un critere de plus proche distance
    - liste des indices de labels de b_1 associes aux indices de labels de b_2 les plus proches
    - concatenee avec la liste des indices de labels de b_2 associes aux indices de labels de b_1 les plus proches
*/
int VT_ComputePairsOfPoints(vt_image_barycentres b_1, vt_image_barycentres b_2, int **p, int *npairs, double* dnorm)
{
    char *proc="VT_ComputePairsOfPoints";
    int n=b_1.n+b_2.n;
    int i;
    vt_barycentre theB;
    int *pairs[2];
    int *IND;
    int N;
    double d, *D, dsum;
    double siz[3];

    D = vtmalloc(n*sizeof(double), "D", proc );
    if(D==NULL) {
        if(_verbose_)
            fprintf(stderr, "%s: unable to allocate distance vector\n", proc);
        return(0);
    }
    pairs[0] = vtmalloc(n*sizeof(int), "pairs[0]", proc );
    if(pairs[0]==NULL) {
        vtfree(D);
        D=NULL;
        if(_verbose_)
            fprintf(stderr, "%s: unable to allocate pairs first vector\n", proc);
        return(0);
    }
    pairs[1] = vtmalloc(n*sizeof(int), "pairs[1]", proc );
    if(pairs[1]==NULL) {
        vtfree(D);
        D=NULL;
        vtfree(pairs[0]);
        pairs[0]=NULL;
        if(_verbose_)
            fprintf(stderr, "%s: unable to allocate pairs second vector\n", proc);
        return(0);
    }

    siz[0]=b_1.vx;
    siz[1]=b_1.vy;
    siz[2]=b_1.vz;
    for (i=0 ; i<b_1.n ; i++) {
        theB=b_1.barycentres[i];

        int j=VT_FindCloser(theB, siz, b_2, &d);
        if (j<0)
        {
            vtfree(D);
            D=NULL;
            vtfree(pairs[0]);
            pairs[0]=NULL;
            vtfree(pairs[1]);
            pairs[1]=NULL;
            if(_verbose_)
                fprintf(stderr, "%s: unexpected result of closer label finding...\n", proc);
            return(0);
        }
        pairs[0][i]=i;
        pairs[1][i]=j;
        D[i]=d;
    }

    siz[0]=b_2.vx;
    siz[1]=b_2.vy;
    siz[2]=b_2.vz;
    for (i=0 ; i<b_2.n ; i++) {
        theB=b_2.barycentres[i];

        int j=VT_FindCloser(theB, siz, b_1, &d);
        if (j<0)
        {
            vtfree(D);
            D=NULL;
            vtfree(pairs[0]);
            pairs[0]=NULL;
            vtfree(pairs[1]);
            pairs[1]=NULL;
            if(_verbose_)
                fprintf(stderr, "%s: unexpected result of closer label finding...\n", proc);
            return(0);
        }
        pairs[0][i+b_1.n]=j;
        pairs[1][i+b_1.n]=i;
        D[i+b_1.n]=d;
    }

    IND = vtmalloc( n*sizeof(int), "IND", proc );
    if (IND==NULL)
    {
        vtfree(D);
        D=NULL;
        vtfree(pairs[0]);
        pairs[0]=NULL;
        vtfree(pairs[1]);
        pairs[1]=NULL;
        if(_verbose_)
            fprintf(stderr, "%s: unable to allocate a vector...\n", proc);
        return(0);
    }

    /* MATLAB SCRIPT: */
    /* for i=1:size(P,2), I=i:size(P,2); ind=find(P(1,I)==P(1,i)); if numel(ind)>1, for j=2:numel(ind), if P(2,i)==P(2,I(ind(j))), Q=[Q, P(:,i)]; end, end, end, end */

    N=0;
    for (i=0; i<n-1 ; i++) {
        int k;
        for (k=1 ; k<n-i ; k++) {
            if(1) {
              if(pairs[0][i]!=pairs[0][i+k]) continue;
              if(pairs[1][i]!=pairs[1][i+k]) continue;
            }
            IND[N++]=i;
            break;
        }
    }
    p[0] = vtmalloc( N*sizeof(int), "p[0]", proc );
    if(p[0]==NULL) {
        vtfree(D);
        D=NULL;
        vtfree(pairs[0]);
        pairs[0]=NULL;
        vtfree(pairs[1]);
        pairs[1]=NULL;
        if(_verbose_)
            fprintf(stderr, "%s: unable to allocate the p[0] vector...\n", proc);
        return(0);
    }
    p[1] = vtmalloc( N*sizeof(int), "p[1]", proc );
    if(p[1]==NULL) {
        vtfree(D);
        D=NULL;
        vtfree(pairs[0]);
        pairs[0]=NULL;
        vtfree(pairs[1]);
        pairs[1]=NULL;
        vtfree(p[0]);
        p[0]=NULL;
        if(_verbose_)
            fprintf(stderr, "%s: unable to allocate the p[1] vector...\n", proc);
        return(0);
    }
    dsum=0;
    for(i=0; i<N; i++) {
        p[0][i]=pairs[0][IND[i]];
        p[1][i]=pairs[1][IND[i]];
        dsum+=D[IND[i]];
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
    vtfree(IND);
    IND=NULL;

    *dnorm=dsum;
    *npairs=N;
    return(1);
}






int VT_ComputePointTrsf(int **pairs, int n, vt_image_barycentres b_1, vt_image_barycentres b_2,
                        double *r, double *sd, double *T,
                        int flag_trsf, bal_estimator estimator)
{
    char *proc="VT_ComputePointTrsf";
    int i;

    vt_barycentre theB;


    /*  Debut des hostilites */


    bal_typeFieldPoint c;

    bal_transformation theResTransformation;
    /* bal_image template; */

    bal_typeFieldPointList theFloatingPoints;
    bal_typeFieldPointList theReferencePoints;

    enumTypeTransfo transformation_type=RIGID_3D;
    if (flag_trsf==1) transformation_type=AFFINE_3D;

    enumUnitTransfo points_unit = REAL_UNIT;


    FIELD theField;


    BAL_InitTypeFieldPointList( &theFloatingPoints );
    BAL_InitTypeFieldPointList( &theReferencePoints );
    BAL_InitTransformation( &theResTransformation );

    for (i=0; i<n ; i++)
    {
        theB=b_1.barycentres[pairs[0][i]];
        c.x = theB.x*b_1.vx;
        c.y = theB.y*b_1.vy;
        c.z = theB.z*b_1.vz;
        if ( BAL_AddTypeFieldPointToTypeFieldPointList( &theReferencePoints, &c ) != 1 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to add measure #%d (=%f %f %f) to list\n",
               proc, i, c.x, c.y, c.z );
            if (i>0){
                BAL_FreeTypeFieldPointList( &theFloatingPoints );
                BAL_FreeTypeFieldPointList( &theReferencePoints );
            }
            return( -1 );
        }

        theB=b_2.barycentres[pairs[1][i]];
        c.x = theB.x*b_2.vx;
        c.y = theB.y*b_2.vy;
        c.z = theB.z*b_2.vz;
        if ( BAL_AddTypeFieldPointToTypeFieldPointList( &theFloatingPoints, &c ) != 1 ) {
            if ( _verbose_ )
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
        double s1=0.0;/* , s2=0.0; */

        /* for ( i=0; i<theField.n_computed_pairs; i++ ) { */
	/* fprintf( f, "%f\n", sqrt( theField.data[i].error ) ); */
	/* the residu = sqrt( theField.data[i].error ); */
        /* } */
        /* fprintf( f, "\n" ); */
        for ( i=0; i<(int)theField.n_selected_pairs; i++ ){
            sum += sqrt( theField.pointer[i]->error );
            s1+=theField.pointer[i]->error;
            /* s2+=theField.pointer[i]->error; */
        }

        if(_verbose_) {
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

    if ( 1 ) {
      int l,c;
      /* int j; */
      l=theResTransformation.mat.l;
      c=theResTransformation.mat.c;
      if (l!=4 || c!=4) {
          fprintf( stderr, "%s: unexpected size of transformation matrix\n", proc);
          return(0);
      }
      /*       for (j=0; j<l; j++) */
      for (i=0; i<c*l; i++) {
          T[i]=theResTransformation.mat.m[i];
      }
      /* if ( BAL_WriteTransformation( &theResTransformation, p.result_real_transformation ) != 1 ) { */
      /*   BAL_FreeTransformation( &theResTransformation ); */
      /*   fprintf( stderr, "%s: unable to write real result transformation '%s'\n", */
      /*        program, p.result_real_transformation  ); */
      /*   _ErrorParse( "error when writing transformation", 0 ); */
      /* } */
    }

    /*if ( p.result_voxel_transformation != NULL ) {
      template.vx = b_2.vx;
      template.vy = b_2.vy;
      template.vz = b_2.vz;
      if ( BAL_ChangeTransformationToVoxelUnit( &template, &template,
                            &theResTransformation,
                            &theResTransformation ) != 1 ) {
        BAL_FreeTransformation( &theResTransformation );
        _ErrorParse( "unable to convert 'real' transformation into the 'voxel' world", 0 );
      }

      if ( BAL_WriteTransformation( &theResTransformation, p.result_voxel_transformation ) != 1 ) {
        BAL_FreeTransformation( &theResTransformation );
        fprintf( stderr, "%s: unable to write voxel result transformation '%s'\n",
             program, p.result_voxel_transformation);
        _ErrorParse( "error when writing transformation", 0 );
      }
    }*/

    BAL_FreeTransformation( &theResTransformation );




    /*  Fin des hostilites */

    return(1);
}


int VT_ComputeCP(vt_image_barycentres b_in, vt_image_barycentres b_ext, double *T_ext_in_init, double **T_ext_in, double *residual, int flag_affine, bal_estimator estimator)
{
  char *proc="VT_ComputeCP";

  int *pairs[2];
  int npairs;

  double sd;
  double r;

  vt_image_barycentres b_in_trsf;

  /* Allocation des images des barycentres transformes */

  if (VT_AllocImageBarycentres(&b_in_trsf, b_in.n) != 1) {
    fprintf(stderr, "%s: unable to allocate image of barycentres b_in_trsf...\n", proc);
    return(0);
  }

  /* Transformation initiale des points de b_in */

  if (VT_TransformImageBarycentres(T_ext_in_init, b_in, &b_in_trsf) != 1) {
      VT_FreeImageBarycentres(&b_in_trsf);
      fprintf(stderr, "%s: unable to apply initial transformation of b_in...\n", proc);
      return(0);
  }

  /* Calcul des appariements de points */

  if (VT_ComputePairsOfPoints(b_in_trsf, b_ext, pairs, &npairs, &r) != 1)
  {
      VT_FreeImageBarycentres(&b_in_trsf);
      fprintf(stderr, "%s: unable to compute the pairings...\n", proc);
      return(0);
  }


  /* Calcul de la nouvelle transformation */

  if( VT_ComputePointTrsf(pairs, npairs, b_in, b_ext, &r, &sd, *T_ext_in, flag_affine, estimator) != 1) {
      VT_FreeImageBarycentres(&b_in_trsf);
      vtfree(pairs[0]); pairs[0]=NULL;
      vtfree(pairs[1]); pairs[1]=NULL;
      fprintf(stderr, "%s: unable to compute new transformation...\n", proc);
      return(0);
  }



  *residual = r;

  /* Frees */
  VT_FreeImageBarycentres(&b_in_trsf);



  return(1);
}
