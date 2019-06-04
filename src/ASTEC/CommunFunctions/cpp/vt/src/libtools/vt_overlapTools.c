/*************************************************************************
 * vt_overlapTools.c - outils pour la comparaison de segmentations
 *
 * $Id: vt_overlapTools.c,v 1.1 2015/05/13 13:58:52 gael Exp $
 *
 * DESCRIPTION:
 *
 *
 * Copyright Gael Michelin, INRIA
 *
 *
 * ADDITIONS, CHANGES:
 *
 *
 */

/* getline(): Feature Test Macro Requirements for glibc
 * Since glibc 2.10:
 *     _POSIX_C_SOURCE >= 200809L || _XOPEN_SOURCE >= 700
 * Before glibc 2.10:
 *     _GNU_SOURCE
 */
#define _XOPEN_SOURCE 700

#include <math.h>
#include <stdio.h>
#include <vt_common.h>
#include <vt_overlapTools.h>
#include <vt_morpho.h>
#include <morphotools.h>




/* bbox : inputs x, y, z must be of type int *{x,y,z}[2];
 */
int bbox(vt_image *img, int **x, int **y, int **z, long long **_volumes, int *nlab)
{
  int i,j,k;
  int val;
  int max = 0;
  long long *volumes=NULL;

  unsigned char ***array8=NULL;
  unsigned short int ***array16=NULL;

  switch(img->type) {
  default:
      return (0);
  case UCHAR:
  case SCHAR:
      array8 = (unsigned char ***) img->array;
      for (k=0 ; k<(int)img->dim.z ; k++)
      for (j=0 ; j<(int)img->dim.y ; j++)
      for (i=0 ; i<(int)img->dim.x ; i++)
        if (max< (int) array8[k][j][i]) max = (int) array8[k][j][i];
      break;
  case USHORT:
  case SSHORT:
      array16 = (unsigned short int ***) img->array;
      for (k=0 ; k<(int)img->dim.z ; k++)
      for (j=0 ; j<(int)img->dim.y ; j++)
      for (i=0 ; i<(int)img->dim.x ; i++)
        if (max< (int) array16[k][j][i]) max = (int) array16[k][j][i];
      break;
  }

  volumes = malloc((max+1)*sizeof(long long));
  if(volumes==NULL)
  {
    fprintf(stderr, "unable to allocate volume vector in bbox function\n");
    return(0);
  }
  for (i=0 ; i<=max; i++) volumes[i]=0;
  for (i=0 ; i<2; i++) {
      x[i] = malloc((max+1)*sizeof(int));
      y[i] = malloc((max+1)*sizeof(int));
      z[i] = malloc((max+1)*sizeof(int));
      if (x[i]==NULL || y[i]==NULL || z[i]==NULL)
      {
          for (j=0 ; j<=i; j++) {
              if (x[j]!=NULL) free(x[j]);
              if (y[j]!=NULL) free(y[j]);
              if (z[j]!=NULL) free(z[j]);
          }
          free(volumes);
          fprintf(stderr, "unable to allocate bounding box vectors in bbox function\n");
          return(0);
      }
      for (j = 0; j<=max; j++) {
          x[i][j]=-1;
          y[i][j]=-1;
          z[i][j]=-1;
      }
  }

  for (k=0 ; k<(int)img->dim.z ; k++)
  for (j=0 ; j<(int)img->dim.y ; j++)
  for (i=0 ; i<(int)img->dim.x ; i++)
  {
    if (array8!=NULL)
      val = (int)array8[k][j][i];
    else
      val = (int)array16[k][j][i];

    volumes[val]++;
    if (x[0][val] < 0 || x[0][val] > i) x[0][val] = i;
    if (x[1][val] < 0 || x[1][val] < i) x[1][val] = i;

    if (y[0][val] < 0 || y[0][val] > j) y[0][val] = j;
    if (y[1][val] < 0 || y[1][val] < j) y[1][val] = j;

    if (z[0][val] < 0 || z[0][val] > k) z[0][val] = k;
    if (z[1][val] < 0 || z[1][val] < k) z[1][val] = k;

  }

  *_volumes = volumes;
  *nlab = max+1;
  return (0);
}

void setLabel(vt_image *imgIn, vt_image *imgOut, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int labelIn, int labelOut)
{
  int x, y, z;
  unsigned char ***array8=NULL;
  unsigned short int ***array16=NULL;
  unsigned char ***out8=NULL;
  unsigned short int ***out16=NULL;
  switch(imgIn->type) {
  default:
      return;
  case UCHAR:
  case SCHAR:
      array8 = (unsigned char ***) imgIn->array;
      break;
  case USHORT:
  case SSHORT:
      array16 = (unsigned short int ***) imgIn->array;
      break;
  }
  switch(imgOut->type) {
  default:
      return;
  case UCHAR:
  case SCHAR:
      out8 = (unsigned char ***) imgOut->array;
      break;
  case USHORT:
  case SSHORT:
      out16 = (unsigned short int ***) imgOut->array;
      break;
  }
  for (z=zmin ; z<=zmax ; z++)
  for (y=ymin ; y<=ymax ; y++)
  for (x=xmin ; x<=xmax ; x++) {
    if (array8 != NULL ) {
      if ((int)array8[z][y][x] != labelIn) continue;
    }
    else {
      if ((int)array16[z][y][x] != labelIn) continue;
    }
    if (out8 != NULL)
      out8[z][y][x] = (unsigned char) labelOut;
    else
      out16[z][y][x] = (unsigned short int) labelOut;
  }
}

long long volume(vt_image *img, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int label)
{
  int x, y, z;
  long long vol=0;
  unsigned char ***array8=NULL;
  unsigned short int ***array16=NULL;
  switch(img->type) {
  default:
      return (0);
  case UCHAR:
  case SCHAR:
      array8 = (unsigned char ***) img->array;
      break;
  case USHORT:
  case SSHORT:
      array16 = (unsigned short int ***) img->array;
      break;
  }
  for (z=zmin ; z<=zmax ; z++)
  for (y=ymin ; y<=ymax ; y++)
  for (x=xmin ; x<=xmax ; x++) {
    if (array8 != NULL ) {
      if ((int)array8[z][y][x] != label) continue;
    }
    else {
      if ((int)array16[z][y][x] != label) continue;
    }
    vol++;
  }

  return(vol);
}

int recouvrement(vt_image *imgIn, vt_image *imgExt,
                        int xminIn, int xmaxIn, int yminIn, int ymaxIn, int zminIn, int zmaxIn,
                        int xminExt, int xmaxExt, int yminExt, int ymaxExt, int zminExt, int zmaxExt,
                        int labelIn __attribute__ ((unused)),
                 int labelExt __attribute__ ((unused)),
                        long long *volumeIn, long long *volumeExt,
                        double *recouvrementIn, double *recouvrementExt)
{
    int x, y, z;
    long long volIn=0;
    long long volExt=0;
    long long volInter=0;

    unsigned char ***arrayIn8=NULL;
    unsigned short int ***arrayIn16=NULL;
    unsigned char ***arrayExt8=NULL;
    unsigned short int ***arrayExt16=NULL;

    int xmin = (xminIn<xminExt) ? xminIn : xminExt;
    int ymin = (yminIn<yminExt) ? yminIn : yminExt;
    int zmin = (zminIn<zminExt) ? zminIn : zminExt;
    int xmax = (xmaxIn<xmaxExt) ? xmaxIn : xmaxExt;
    int ymax = (ymaxIn<ymaxExt) ? ymaxIn : ymaxExt;
    int zmax = (zmaxIn<zmaxExt) ? zmaxIn : zmaxExt;

    int boolIn=0, boolExt=0;

    switch(imgIn->type) {
    default:
        return(0);
    case UCHAR:
    case SCHAR:
        arrayIn8 = (unsigned char ***) imgIn->array;
        break;
    case USHORT:
    case SSHORT:
        arrayIn16 = (unsigned short int ***) imgIn->array;
        break;
    }
    switch(imgExt->type) {
    default:
        return(0);
    case UCHAR:
    case SCHAR:
        arrayExt8 = (unsigned char ***) imgExt->array;
        break;
    case USHORT:
    case SSHORT:
        arrayExt16 = (unsigned short int ***) imgExt->array;
        break;
    }

    for (z=zmin ; z<=zmax ; z++)
    for (y=ymin ; y<=ymax ; y++)
    for (x=xmin ; x<=xmax ; x++) {

      if (x>=xminIn && x<=xmaxIn && y>=yminIn && y<=ymaxIn && z>=zminIn && z<=zmaxIn) {
        if ((arrayIn8 != NULL && (int)arrayIn8[z][y][x] == labelIn) || (arrayIn16 != NULL && (int)arrayIn16[z][y][x] == labelIn)) {
          boolIn=1;
          volIn++;
        }
        else boolIn=0;
      }
      else boolIn=0;


      if (x>=xminExt && x<=xmaxExt && y>=yminExt && y<=ymaxExt && z>=zminExt && z<=zmaxExt) {
        if ((arrayExt8 != NULL && (int)arrayExt8[z][y][x] == labelExt) || (arrayExt16 != NULL && (int)arrayExt16[z][y][x] == labelExt)) {
          boolExt=1;
          volExt++;
        }
        else boolExt=0;
      }
      else boolExt=0;

      if (boolIn && boolExt)
          volInter++;
    }

    *volumeIn = volIn;
    *volumeExt= volExt;
    *recouvrementIn = ((double) volInter)/volIn;
    *recouvrementExt = ((double) volInter)/volExt;

    return(1);
}


extern float overlapFun ( int ab, int a, int b, measureType measure )
{
    switch (measure) {
    default:
    case PROBABILITY:
        if (a==0) return (0.0);
        return(((float)ab)/a);
        break;
    case DICE:
        if (a+b==0) return (0.0);
        return((2.0*ab)/(a+b));
        break;
    case JACCARD:
        if (a+b-ab==0) return (0.0);
        return(((float)ab)/(a+b-ab));
        break;
    }
    return(0);
}


extern int Neighborhood2Int ( Neighborhood N )
{
  int connectivity = 26;
  switch ( N ) {
  case N04 :
    connectivity = 4; break;
  case N06 :
    connectivity = 6; break;
  case N08 :
    connectivity = 8; break;
  case N10 :
    connectivity = 10; break;
  case N18 :
    connectivity = 18; break;
  case N26 :
    connectivity = 26; break;
  }
  return( connectivity );
}


extern int VT_SegmentationOverlapping( vt_image *imA, vt_image *imB, vt_image *imAB, measureType measure, double R, int flag_real, int bckgrdA, int bckgrdB)
{
    vt_image imAr;
    vt_image imBr;
    int a, b, ab;
    int nslices;


    unsigned char *bufAu8=NULL, *bufBu8=NULL;
    unsigned char *bufAru8=NULL, *bufBru8=NULL;
    unsigned short int *bufAu16=NULL, *bufBu16=NULL;
    unsigned short int *bufAru16=NULL, *bufBru16=NULL;

    int i;
    int j;
    int x,y,z;
    int aMax=0, bMax=0;

    int r;

    int nRows, nCols;
    int **AB=NULL, **BA=NULL;
    int iA,iAr,jB,jBr;

    char name[256];


    /* tests de compatibilite d'images */


    if (imA->dim.x != imB->dim.x || imA->dim.y != imB->dim.y || imA->dim.z != imB->dim.z )
    {
        fprintf(stderr, "unmatching input image dimensions\n");
        return(0);
    }
    if (imA->siz.x != imB->siz.x || imA->siz.y != imB->siz.y || imA->siz.z != imB->siz.z )
    {
        fprintf(stderr,"unmatching input image voxel sizes, comparison algorithm may return unexpected values\n");
    }

    /* types d'images */
    switch (imA->type) {
    case UCHAR:
    case SCHAR:
        bufAu8=(unsigned char*)imA->buf;
        for (i=0 ; i<(int)(imA->dim.x*imA->dim.y*imA->dim.z) ; i++)
            if (aMax<(int)bufAu8[i]) aMax = (int)bufAu8[i];
        break;
    case USHORT:
    case SSHORT:
        bufAu16=(unsigned short int*)imA->buf;
        for (i=0 ; i<(int)(imA->dim.x*imA->dim.y*imA->dim.z) ; i++)
            if (aMax<(int)bufAu16[i]) aMax = (int)bufAu16[i];
        break;
    default:
        VT_FreeImage(imA);
        VT_FreeImage(imB);
        fprintf(stderr, "unexpected input image A type\n");
        return(0);

    }
    switch (imB->type) {
    case UCHAR:
    case SCHAR:
        bufBu8=(unsigned char*)imB->buf;
        for (i=0 ; i<(int)(imB->dim.x*imB->dim.y*imB->dim.z) ; i++)
            if (bMax<(int)bufBu8[i]) bMax = (int)bufBu8[i];
        break;
    case USHORT:
    case SSHORT:
        bufBu16=(unsigned short int*)imB->buf;
        for (i=0 ; i<(int)(imB->dim.x*imB->dim.y*imB->dim.z) ; i++)
            if (bMax<(int)bufBu16[i]) bMax = (int)bufBu16[i];
        break;
    default:
        VT_FreeImage(imA);
        VT_FreeImage(imB);
        fprintf(stderr, "unexpected input image B type\n");
        return(0);
    }


    /* Alloc & init images erodees */

    if (R<=0) {
        switch (imA->type) {
        case UCHAR:
        case SCHAR:
            for (i=0 ; i<(int)(imA->dim.x*imA->dim.y*imA->dim.z) ; i++)
                if(bufAu8[i] == (unsigned char) bckgrdA)
                    bufAu8[i]=(unsigned char) 0;
            bufAru8=bufAu8;
            break;
        case USHORT:
        case SSHORT:
            for (i=0 ; i<(int)(imA->dim.x*imA->dim.y*imA->dim.z) ; i++)
                if(bufAu16[i] == (unsigned short int) bckgrdA)
                    bufAu16[i]=(unsigned short int) 0;
            bufAru16=bufAu16;
            break;
        default:
            VT_FreeImage(imA);
            VT_FreeImage(imB);
            fprintf(stderr, "unexpected input image A type\n");
            return(0);
        }
        switch (imB->type) {
        case UCHAR:
        case SCHAR:
            for (i=0 ; i<(int)(imA->dim.x*imA->dim.y*imA->dim.z) ; i++)
                if(bufBu8[i] == (unsigned char) bckgrdB)
                    bufBu8[i]=(unsigned char) 0;
            bufBru8=bufBu8;
            break;
        case USHORT:
        case SSHORT:
            for (i=0 ; i<(int)(imA->dim.x*imA->dim.y*imA->dim.z) ; i++)
                if(bufBu16[i] == (unsigned short int ) bckgrdB)
                    bufBu16[i]=(unsigned short int) 0;
            bufBru16=bufBu16;
            break;
        default:
            VT_FreeImage(imA);
            VT_FreeImage(imB);
            fprintf(stderr, "unexpected input image B type\n");
            return(0);
        }
    }
    else {
      int t = (int) imA->type;
      sprintf( name, "%s", "Aerod.inr" );
      VT_InitFromImage(&imAr, imA, name, t);
      if ( VT_AllocImage( &(imAr) ) != 1 ) {
          VT_FreeImage( imA );
          VT_FreeImage( imB );
          fprintf(stderr, "problem while allocating the eroded of A\n");
          return(0);
      }
      t = (int) imB->type;
      sprintf( name, "%s", "Berod.inr" );
      VT_InitFromImage(&imBr, imB, name, t);
      if ( VT_AllocImage( &(imBr) ) != 1 ) {
          VT_FreeImage( imA );
          VT_FreeImage( imB );
          VT_FreeImage( &imAr );
          fprintf(stderr, "problem while allocating the eroded of B\n");
          return(0);
      }
      switch(imA->type){
      default:
      case UCHAR:
      case SCHAR:
          bufAru8 = (unsigned char *)imAr.buf;
          for (i=0 ; i<(int)(imA->dim.x*imA->dim.y*imA->dim.z) ; i++) {
              if(bufAu8[i] == (unsigned char) bckgrdA) bufAu8[i]=(unsigned char) 0;
              bufAru8[i]=bufAu8[i];
          }
          break;
      case USHORT:
      case SSHORT:
          bufAru16 = (unsigned short int *)imAr.buf;
          for (i=0 ; i<(int)(imA->dim.x*imA->dim.y*imA->dim.z) ; i++){
              if(bufAu16[i] == (unsigned short int) bckgrdA) bufAu16[i]=(unsigned short int) 0;
              bufAru16[i]=bufAu16[i];
          }
          break;
      }
      switch(imB->type){
      default:
      case UCHAR:
      case SCHAR:
          bufBru8 = (unsigned char *)imBr.buf;
          for (i=0 ; i<(int)(imA->dim.x*imA->dim.y*imA->dim.z) ; i++){
              if(bufBu8[i] == (unsigned char) bckgrdB) bufBu8[i]=(unsigned char) 0;
              bufBru8[i]=bufBu8[i];
          }
          break;
      case USHORT:
      case SSHORT:
          bufBru16 = (unsigned short int *)imBr.buf;
          for (i=0 ; i<(int)(imA->dim.x*imA->dim.y*imA->dim.z) ; i++){
              if(bufBu16[i] == (unsigned short int) bckgrdB) bufBu16[i]=(unsigned short int) 0;
              bufBru16[i]=bufBu16[i];
          }
          break;
      }



      /* images erosion */

      i=0;
      for (z=0 ; z<(int)imA->dim.z ; z++)
      for (y=0 ; y<(int)imA->dim.y ; y++)
      for (x=0 ; x<(int)imA->dim.x ; x++) {

          if(x<(int)imA->dim.x-1 ) {
              switch (imA->type ) {
              default:
              case UCHAR:
              case SCHAR:
                  if (bufAu8[i] != bufAu8[i+1]) {
                      bufAru8[i]   = (unsigned char) 0;
                      bufAru8[i+1] = (unsigned char) 0;
                  }
                  break;
              case USHORT:
              case SSHORT:
                  if (bufAu16[i] != bufAu16[i+1]) {
                      bufAru16[i]   = (unsigned short int) 0;
                      bufAru16[i+1] = (unsigned short int) 0;
                  }
                  break;
              }
              switch (imB->type ) {
              default:
              case UCHAR:
              case SCHAR:
                  if (bufBu8[i] != bufBu8[i+1]) {
                      bufBru8[i]   = (unsigned char) 0;
                      bufBru8[i+1] = (unsigned char) 0;
                  }
                  break;
              case USHORT:
              case SSHORT:
                  if (bufBu16[i] != bufBu16[i+1]) {
                      bufBru16[i]   = (unsigned short int) 0;
                      bufBru16[i+1] = (unsigned short int) 0;
                  }
                  break;
              }
          }

          if(y<(int)imA->dim.y-1 ) {
              switch (imA->type ) {
              default:
              case UCHAR:
              case SCHAR:
                  if (bufAu8[i] != bufAu8[i+imA->dim.x]) {
                      bufAru8[i]   = (unsigned char) 0;
                      bufAru8[i+imA->dim.x] = (unsigned char) 0;
                  }
                  break;
              case USHORT:
              case SSHORT:
                  if (bufAu16[i] != bufAu16[i+imA->dim.x]) {
                      bufAru16[i]   = (unsigned short int) 0;
                      bufAru16[i+imA->dim.x] = (unsigned short int) 0;
                  }
                  break;
              }
              switch (imB->type ) {
              default:
              case UCHAR:
              case SCHAR:
                  if (bufBu8[i] != bufBu8[i+imA->dim.x]) {
                      bufBru8[i]   = (unsigned char) 0;
                      bufBru8[i+imA->dim.x] = (unsigned char) 0;
                  }
                  break;
              case USHORT:
              case SSHORT:
                  if (bufBu16[i] != bufBu16[i+imA->dim.x]) {
                      bufBru16[i]   = (unsigned short int) 0;
                      bufBru16[i+imA->dim.x] = (unsigned short int) 0;
                  }
                  break;
              }
          }

          if(z<(int)imA->dim.z-1 ) {
              switch (imA->type ) {
              default:
              case UCHAR:
              case SCHAR:
                  if (bufAu8[i] != bufAu8[i+imA->dim.x*imA->dim.y]) {
                      bufAru8[i]   = (unsigned char) 0;
                      bufAru8[i+imA->dim.x*imA->dim.y] = (unsigned char) 0;
                  }
                  break;
              case USHORT:
              case SSHORT:
                  if (bufAu16[i] != bufAu16[i+imA->dim.x*imA->dim.y]) {
                      bufAru16[i]   = (unsigned short int) 0;
                      bufAru16[i+imA->dim.x*imA->dim.y] = (unsigned short int) 0;
                  }
                  break;
              }
              switch (imB->type ) {
              default:
              case UCHAR:
              case SCHAR:
                  if (bufBu8[i] != bufBu8[i+imA->dim.x*imA->dim.y]) {
                      bufBru8[i]   = (unsigned char) 0;
                      bufBru8[i+imA->dim.x*imA->dim.y] = (unsigned char) 0;
                  }
                  break;
              case USHORT:
              case SSHORT:
                  if (bufBu16[i] != bufBu16[i+imA->dim.x*imA->dim.y]) {
                      bufBru16[i]   = (unsigned short int) 0;
                      bufBru16[i+imA->dim.x*imA->dim.y] = (unsigned short int) 0;
                  }
                  break;
              }
          }

          i++;
      }

      if (flag_real == 0)
        r = (int)R - 1;
      else {
        if ( R > 0 && (imA->siz.x != imA->siz.y || imA->siz.x != imA->siz.z ) )
          fprintf(stderr, "warning: erosion is anosotropic since image voxel resolution is anisotropic\n");
        r = (int)(R / imA->siz.x)-1;
      }


      if(r>0)
      {
            int theDim[]={imA->dim.x, imA->dim.y, imA->dim.z};
            typeStructuringElement SE;
            initStructuringElement(&SE);
            SE.nbIterations = 1;
            SE.connectivity = Neighborhood2Int(N06);
            SE.radius = r;

            switch (imA->type) {
            default:
            case UCHAR:
            case SCHAR:
                if (morphologicalErosion(bufAru8, bufAru8, imA->type, theDim, &SE) != 1) {
                    freeStructuringElement( &SE );
                    VT_FreeImage(imA);
                    VT_FreeImage(imB);
                    if (R > 0) {
                      VT_FreeImage(&imAr);
                      VT_FreeImage(&imBr);
                    }
                    fprintf(stderr, "morphological erosion process failed for image A\n");
                }
                break;
            case USHORT:
            case SSHORT:
                if (morphologicalErosion(bufAru16, bufAru16, imA->type, theDim, &SE) != 1) {
                    freeStructuringElement( &SE );
                    VT_FreeImage(imA);
                    VT_FreeImage(imB);
                    if (R > 0) {
                      VT_FreeImage(&imAr);
                      VT_FreeImage(&imBr);
                    }
                    fprintf(stderr, "morphological erosion process failed for image A\n");
                }
                break;
            }

            switch (imB->type) {
            default:
            case UCHAR:
            case SCHAR:
                if (morphologicalErosion(bufBru8, bufBru8, imB->type, theDim, &SE) != 1) {
                    freeStructuringElement( &SE );
                    VT_FreeImage(imA);
                    VT_FreeImage(imB);
                    if (R > 0) {
                      VT_FreeImage(&imAr);
                      VT_FreeImage(&imBr);
                    }
                    fprintf(stderr, "morphological erosion process failed for image B\n");
                }
                break;
            case USHORT:
            case SSHORT:
                if (morphologicalErosion(bufBru16, bufBru16, imA->type, theDim, &SE) != 1) {
                    freeStructuringElement( &SE );
                    VT_FreeImage(imA);
                    VT_FreeImage(imB);
                    if (R > 0) {
                      VT_FreeImage(&imAr);
                      VT_FreeImage(&imBr);
                    }
                    fprintf(stderr, "morphological erosion process failed for image B\n");
                }
                break;
            }
            freeStructuringElement( &SE );
      }
    }



    /*--- Allocation des tableaux d'intersections ---*/

    /* allocating tabulars AB & BA
     */


    nRows=aMax+1;
    nCols=bMax+1;

    AB=malloc(nRows*sizeof(int*));
    if(AB==NULL) {
        VT_FreeImage(imA);
        VT_FreeImage(imB);
        if (R > 0) {
          VT_FreeImage(&imAr);
          VT_FreeImage(&imBr);
        }
        fprintf(stderr, "Unable to allocate the intersections tabular AB\n");
    }
    for (i=0 ; i<nRows ; i++)
    {
        AB[i]=malloc(nCols*sizeof(int));
        if(AB[i]==NULL){
            VT_FreeImage(imA);
            VT_FreeImage(imB);
            if (R > 0) {
              VT_FreeImage(&imAr);
              VT_FreeImage(&imBr);
            }
            for (j=0 ; j<i ; j++)
            {
                free(AB[j]);
                AB[j]=NULL;
            }
            free(AB);
            AB=NULL;
            fprintf(stderr, "Unable to allocate the intersections tabular AB\n");
        }
    }
    if (measure == PROBABILITY) {
      BA=malloc(nCols*sizeof(int*));
      if(BA==NULL) {
        VT_FreeImage(imA);
        VT_FreeImage(imB);
        if (R > 0) {
          VT_FreeImage(&imAr);
          VT_FreeImage(&imBr);
        }
        fprintf(stderr, "Unable to allocate the intersections tabular BA\n");
      }
      for (i=0 ; i<nCols ; i++)
      {
        BA[i]=malloc(nRows*sizeof(int));
        if(BA[i]==NULL){
            VT_FreeImage(imA);
            VT_FreeImage(imB);
            if (R > 0) {
              VT_FreeImage(&imAr);
              VT_FreeImage(&imBr);
            }
            for (j=0 ; j<i ; j++)
            {
                free(BA[j]);
                BA[j]=NULL;
            }
            free(BA);
            BA=NULL;
            fprintf(stderr, "Unable to allocate the intersections tabular BA\n");
        }
      }
    }


    /*--- initialisation des tableaux d'intersection ---*/

    /* initializing tabulars AB & BA
     */

    for (i = 0 ; i<nRows ; i++)
    for (j = 0 ; j<nCols ; j++)
      AB[i][j] = 0;
    if (measure == PROBABILITY)
      for (j = 0 ; j<nCols ; j++)
      for (i = 0 ; i<nRows ; i++)
        BA[j][i] = 0;


    /*--- remplissage des tableaux d'intersection ---*/
    /*
     * Pour la mesure PROBABILITY:
     *
     * AB[i][j] = nb de voxels tq imAerode[vox]==i && imB[vox]==j
     * BA[j][i] = nb de voxels tq imA[vox]==i && imBerod[vox]==j
     * AB[i][0] = histogramme des voxels tq imAerode[vox]==i
     * AB[0][j] = histogramme des voxels tq imB[vox]==j
     * BA[j][0] = histogramme des voxels tq imBerode[vox]==j
     * BA[0][i] = histogramme des voxels tq imA[vox]==i
     *
     * Pour les mesures DICE/JACCARD:
     *
     * AB[i][j] = nb de voxels tq imAerode[vox]==i && imB[vox]==j
     * AB[i][0] = histogramme des voxels tq imAerode[vox]==i
     * AB[0][j] = histogramme des voxels tq imB[vox]==j
     *
     */

    /* filling tabulars AB & BA
     */

    for (x = 0 ; x<(int)(imA->dim.x*imA->dim.y*imA->dim.z) ; x++) {
        switch (imA->type) {
        default:
        case UCHAR:
        case SCHAR:
            iAr=(int)bufAru8[x];
            break;
        case USHORT:
        case SSHORT:
            iAr=(int)bufAru16[x];
            break;
        }
        switch (imB->type) {
        default:
        case UCHAR:
        case SCHAR:
            jB=(int)bufBu8[x];
            break;
        case USHORT:
        case SSHORT:
            jB=(int)bufBu16[x];
            break;
        }

        switch (imA->type) {
        default:
        case UCHAR:
        case SCHAR:
            iA=(int)bufAu8[x];
            break;
        case USHORT:
        case SSHORT:
            iA=(int)bufAu16[x];
            break;
        }
        switch (imB->type) {
        default:
        case UCHAR:
        case SCHAR:
            jBr=(int)bufBru8[x];
            break;
        case USHORT:
        case SSHORT:
            jBr=(int)bufBru16[x];
            break;
        }

        /* if(aMax<iA) aMax=iA;
         * if(bMax<jB) bMax=jB;
         */

        if (measure == PROBABILITY) {
          if (iAr>0)
            AB[iAr][0]+=1;
          if (jB>0)
            AB[0][jB]+=1;
          if (iAr>0 && jB>0)
            AB[iAr][jB]+=1;

          if (iA>0)
            BA[0][iA]+=1;
          if (jBr>0)
            BA[jBr][0]+=1;
          if (iA>0 && jBr>0)
            BA[jBr][iA]+=1;

        }
        else {
          if (iAr>0)
              AB[iAr][0]+=1;
          if (jBr>0)
              AB[0][jBr]+=1;
          if (iAr>0 && jBr>0)
                AB[iAr][jBr]+=1;
        }


    }


    /*--- Matrices de recouvrement ---*/

    /* float *bufAB;
     * allocation of imAB
     */

    if (measure == PROBABILITY) nslices = 2;
    else nslices = 1;

    VT_InitImage(imAB, "tmp", aMax+1, bMax+1, nslices, (int)FLOAT);

    if (VT_AllocImage(imAB)!=1) {
        VT_FreeImage(imA);
        VT_FreeImage(imB);
        if (R > 0) {
          VT_FreeImage(&imAr);
          VT_FreeImage(&imBr);
        }
        for (i=0; i<nRows ; i++)
        {
            free(AB[i]); AB[i]=NULL;
        }
        free(AB); AB=NULL;
        if (measure == PROBABILITY)
          for (i=0; i<nCols ; i++)
          {
            free(BA[i]); BA[i]=NULL;
            free(BA); BA=NULL;
          }
        fprintf(stderr, "unable to allocate the A/B overlapping image\n");
    }



    /* bufAB=(float*)imAB.buf;
     */
    float ***arrayAB=(float***)imAB->array;

    /* init 1st row & 1st col
     */
    for (i=0 ; i<=aMax ; i++) {
      /* bufAB[i*(bMax+1)] = (BA[0][i] == 0) ? 0.0 : 1.0;
       * bufAB[i*(bMax+1)+(aMax+1)*(bMax+1)] = (BA[0][i] == 0) ? 0.0 : 1.0;
       */
      if (nslices > 1) {
        arrayAB[0][0][i] = (BA[0][i] == 0) ? 0.0 : 1.0;
        arrayAB[1][0][i] = (BA[0][i] == 0) ? 0.0 : 1.0;
      }
      else {
        arrayAB[0][0][i] = (AB[i][0] == 0) ? 0.0 : 1.0;
      }
    }
    for (j=0 ; j<=bMax ; j++) {
      arrayAB[0][j][0] = (AB[0][j] == 0) ? 0.0 : 1.0;
      if (nslices > 1)
        arrayAB[1][j][0]  = (AB[0][j] == 0) ? 0.0 : 1.0;
    }


    /* filling imAB[0]
     */
    for (i=1 ; i<=aMax ; i++) {
      a = AB[i][0];
      for (j=1 ; j<=bMax ; j++)  {
        /* recouvrement de i / j
         */
        b = AB[0][j];
        ab = AB[i][j];
        /* if(i==560 && j==724)
         *  fprintf(stdout, "pixel %d, %d : a = %d   b = %d   ab = %d   val = %f   \n", i, j, a, b, ab, overlapFun(ab,a,b, par.measure));
         * bufAB[j+i*(bMax+1)] = overlapFun(ab,a,b, par.measure);
         */
        arrayAB[0][j][i] = overlapFun(ab,a,b, measure);
      }
    }
    /* filling imAB[1]
     */
    if (measure == PROBABILITY) {
        for (j=1 ; j<=bMax ; j++) {
        b = BA[j][0];
        for (i=1 ; i<=aMax ; i++)  {
          /* recouvrement de j / i
           */
          a = BA[0][i];
          ab = BA[j][i];
          /* if(i==560 && j==724) {
           *    fprintf(stdout, "pixel %d, %d : a = %d   b = %d   ab = %d   val = %f   index = %d\n", i, j, a, b, ab, overlapFun(ab,b,a, par.measure), j+i*(bMax+1)+(aMax+1)*(bMax+1));
           * bufAB[j+i*(bMax+1)+(aMax+1)*(bMax+1)] = overlapFun(ab,b,a, par.measure);
           * }
           */
          arrayAB[1][j][i] = overlapFun(ab,b,a, measure);
        }
      }
    }


    /* liberation memoire */
    for (i=0; i<nRows ; i++)
    {
        free(AB[i]); AB[i]=NULL;
    }
    free(AB); AB=NULL;
    if (measure == PROBABILITY)
    for (i=0; i<nCols ; i++)
    {
        free(BA[i]); BA[i]=NULL;
    }
    free(BA); BA=NULL;
    if (R > 0) {
      VT_FreeImage(&imAr);
      VT_FreeImage(&imBr);
    }

    return (1);
}



extern int VT_OverlapPruning(vt_image *im, double epsilon)
{
    float ***array=NULL;

    int i,j,k;
    int sizeA, sizeB;
    int nSlices;

    /* tests de compatibilite d'images */
    if (im->type != FLOAT )
    {
        fprintf(stderr, "input image expected type is FLOAT\n");
        return(0);
    }
    array = (float ***) im->array;

    sizeA=im->dim.x;
    sizeB=im->dim.y;
    nSlices=im->dim.z;

    /* analyse */

    for (k=0 ; k<nSlices ; k++)
    for (j=1 ; j<sizeB ; j++)
    for (i=1 ; i<sizeA ; i++)
        if (array[k][j][i]<epsilon)
            array[k][j][i] = 0.0;

    return (1);
}


extern int VT_OverlapAnalysis(vt_image *imAB, vt_image *imL, analysisMethod method, double lt, double ht, unsigned short int *lab, int *sizea, int *sizeb)
{

    vt_image imG;
    float ***array=NULL;
    unsigned char ***garray=NULL;
    unsigned short int ***larray=NULL;
    int i,j,k=0;
    int sizeA, sizeB;
    int nSlices;
    int *row=NULL, *col=NULL;
    int nrow, ncol;
    int r, c;
    float max;
    unsigned short int LAB;

    /* tests de compatibilite d'images */
    if (imAB->type != FLOAT )
    {
        fprintf(stderr, "input image expected type is FLOAT\n");
    }
    array = (float ***) imAB->array;

    sizeA=imAB->dim.x;
    sizeB=imAB->dim.y;
    nSlices=imAB->dim.z;

    /* Allocations */
    VT_InitImage(&imG, "imG", sizeA, sizeB, 1, (int)UCHAR);
    if (VT_AllocImage(&imG) != 1)
    {
        fprintf(stderr, "unable to allocate the G image\n");
    }
    garray = (unsigned char ***)imG.array;
    row = malloc(sizeA*sizeB*sizeof(int));
    if (row==NULL)
    {
        VT_FreeImage(&imG);
        fprintf(stderr, "unable to allocate the row vector\n");
    }
    col = malloc(sizeA*sizeB*sizeof(int));
    if (col==NULL)
    {
        VT_FreeImage(&imG);
        free(row); row=NULL;
        fprintf(stderr, "unable to allocate the col vector\n");
    }
    VT_InitImage(imL, "imL", sizeA, sizeB, 1, (int)USHORT);
    if (VT_AllocImage(imL) != 1)
    {
        VT_FreeImage(&imG);
        free(row); row=NULL;
        free(col); col=NULL;
        fprintf(stderr, "unable to allocate the L image\n");
    }
    larray = (unsigned short int ***)imL->array;


    /* analyse */

    if(lt < 0) lt = 1.0;
    if(ht < 0) ht = 0.0;

    /* G
     */
    for (j=0 ; j<sizeB ; j++) garray[0][j][0] = (array[0][j][0] == 0.0) ? (unsigned char) 0 : (unsigned char) 1;
    for (i=0 ; i<sizeA ; i++) garray[0][0][i] = (array[0][0][i] == 0.0) ? (unsigned char) 0 : (unsigned char) 1;
    switch (method) {
    default:
    case NONNULL :
      for (j=1 ; j<sizeB ; j++)
      for (i=1 ; i<sizeA ; i++) {
        garray[0][j][i] = (unsigned char)0;
        if (array[0][j][i] > ht ) {
          /*for (k = 1 ; k<sizeB ; k++)
          {
            if (k==j || array[0][k][i] < lt)
              continue;
            break;
          }
          if (k<sizeB)
            continue;
          for (k = 1 ; k<sizeA ; k++)
          {
            if (k==i || array[0][j][k] < lt)
              continue;
            break;
          }
          if (k<sizeA)
            continue;
          */
          garray[0][j][i] = (unsigned char)1;
          continue;
        }
        if ((nSlices>1 && array[1][j][i]) > ht)
        {
          /*for (k = 1 ; k<sizeB ; k++)
          {
            if (k==j || array[1][k][i] < lt)
              continue;
            break;
          }
          if (k<sizeB)
            continue;
          for (k = 1 ; k<sizeA ; k++)
          {
            if (k==i || array[1][j][k] < lt)
              continue;
            break;
          }
          if (k<sizeA)
            continue;
          */
          garray[0][j][i] = (unsigned char)1;
        }
      }
      break;
    case MAX :
      for (i=1 ; i<sizeA ; i++) {
        max = 0;
        for (j=1 ; j<sizeB ; j++) {
          garray[0][j][i] = (unsigned char)0;
          if (array[0][j][i] > max)
          {
              max = array[0][j][i];
              k = j;
          }
        }
        if (max>0)
          garray[0][k][i] = (unsigned char)255;
      }
      for (j=1 ; j<sizeB ; j++) {
        max = 0;
        for (i=1 ; i<sizeA ; i++) {
          if (array[nSlices-1][j][i] > max)
          {
              max = array[nSlices-1][j][i];
              k = i;
          }
        }
        if (max>0)
          garray[0][j][k] = (unsigned char)255;
      }
      break;
    }

    /* L
     */
    LAB = (unsigned short int)0;

    for (i=0 ; i<sizeA ; i++)
     larray[0][0][i] = (unsigned short int) garray[0][0][i] ;

    for (j=1 ; j<sizeB ; j++) {
     larray[0][j][0] = (unsigned short int) garray[0][j][0] ;
     for (i=1 ; i<sizeA ; i++) {
      if (garray[0][j][i] == (unsigned char)0) continue;
      if (larray[0][j][i] != 0) continue;
      /* LAB = LAB+1;
       */
      larray[0][j][i] = ++LAB;
      col[0] = i;
      ncol = 1;
      row[0] = j;
      nrow = 1;

      while (ncol>0) {
        for (c = 0 ; c<ncol ; c++)
        for (k = 1 ; k<sizeB ; k++)
        if (garray[0][k][col[c]] != (unsigned char)0 && larray[0][k][col[c]] == (unsigned short int)0)  {
          larray[0][k][col[c]] = LAB;
          row[nrow++]=k;
        }
        ncol = 0;

        for (r = 0 ; r<nrow ; r++)
        for (k = 1 ; k<sizeA ; k++)
        if (garray[0][row[r]][k] != (unsigned char)0 && larray[0][row[r]][k] == (unsigned short int)0)  {
          larray[0][row[r]][k] = LAB;
          col[ncol++]=k;
        }
        nrow = 0;
      }
     }
    }

    *lab = LAB;
    *sizea = sizeA;
    *sizeb = sizeB;

    /*--- liberations memoire ---*/

    VT_FreeImage(&imG);
    free(row); row=NULL;
    free(col); col=NULL;

    return(1);
}

extern void VT_WriteOverlapAnalysis(vt_image *imL, char **output, unsigned short int LAB, int sizeA, int sizeB, outputFormat format)
{
    char *program = "VT_WriteOverlapAnalysis";
    int i, j, k;
    int *nlabelsA, *nlabelsB, **labelsA, **labelsB;
    unsigned short int ***larray=(unsigned short int ***)imL->array;
    char *out;
    char tmp[20];

    int nchar;
    switch (format) {
    default:
    case COMPLETE:
        nchar = (5+5) * sizeA + (5+5) * sizeB; /* 5 pour le nb de caracteres max des labels, et 5 pour les espaces, tirets et retour chariots...
                                                */
        break;
    case LIGHT:
        nchar = (4+5+5)*LAB; /* 4 : ' ' + '-' + ' ' + '\n'
                              * 5 & 5 : nombre de caracteres max necessaire pour ecrire chacun des nombres...
                              */
        break;
    }


    nlabelsA=malloc(LAB*sizeof(int));
    if (nlabelsA == NULL) return;
    nlabelsB=malloc(LAB*sizeof(int));
    if (nlabelsB == NULL)
    {
        free(nlabelsA);
        fprintf(stderr, "%s: unable to allocate memory\n", program);
        return;
    }
    labelsA=malloc(LAB*sizeof(int*));
    if (labelsA == NULL)
    {
        free(nlabelsA);
        free(nlabelsB);
        fprintf(stderr, "%s: unable to allocate memory\n", program);
        return;
    }
    labelsB=malloc(LAB*sizeof(int*));
    if (labelsB == NULL)
    {
        free(nlabelsA);
        free(nlabelsB);
        free(labelsA);
        fprintf(stderr, "%s: unable to allocate memory\n", program);
        return;
    }

    for (i=0 ; i<LAB ; i++)
    {
        labelsA[i] = malloc(sizeA*sizeof(int));
        if (labelsA[i] == NULL) {
            for (j=0 ; j<i ; j++)
            {
                free(labelsA[j]);
                free(labelsB[j]);
            }
            free(labelsA);
            free(labelsB);
            free(nlabelsA);
            free(nlabelsB);
            fprintf(stderr, "%s: unable to allocate memory\n", program);
            return;
        }
        labelsB[i] = malloc(sizeB*sizeof(int));
        if (labelsB[i] == NULL) {
            for (j=0 ; j<i ; j++)
            {
                free(labelsA[j]);
                free(labelsB[j]);
            }
            free(labelsA[i]);
            free(labelsA);
            free(labelsB);
            free(nlabelsA);
            free(nlabelsB);
            fprintf(stderr, "%s: unable to allocate memory\n", program);
            return;
        }
        nlabelsA[i]= 0;
        nlabelsB[i]= 0;
    }

    out = malloc(nchar*sizeof(char));
    if (out == NULL) {
        for (i=0 ; i<LAB ; i++)
        {
            free(labelsA[i]);
            free(labelsB[i]);
        }
        free(labelsA);
        free(labelsB);
        free(nlabelsA);
        free(nlabelsB);
        fprintf(stderr, "%s: unable to allocate memory\n", program);
        return;
    }
    out[0] = '\0';


    for (j=1 ; j<sizeB ; j++)
    for (i=1 ; i<sizeA ; i++) {
        /*
         */
        int  l = (int)larray[0][j][i];
        if (l == 0) continue;
        /* A
         */
        for (k = 0 ; k<nlabelsA[l-1] ; k++) {
            if(labelsA[l-1][k] == i)
                break;
        }
        if (k == nlabelsA[l-1]) {
            labelsA[l-1][nlabelsA[l-1]] = i;
            nlabelsA[l-1] += 1;
        }
        /* B
         */
        for (k = 0 ; k<nlabelsB[l-1] ; k++) {
            if(labelsB[l-1][k] == j)
                break;
        }
        if (k == nlabelsB[l-1]) {
            labelsB[l-1][nlabelsB[l-1]] = j;
            nlabelsB[l-1] += 1;
        }
    }

    for (i=0 ; i<LAB ; i++) {
      if (format == COMPLETE)
      {
        for (k=0 ; k<nlabelsA[i] ; k++) {
            sprintf(tmp, "%d ", labelsA[i][k]);
            strcat(out, tmp);
        }
        sprintf(tmp, "- ");
        strcat(out, tmp);
        for (k=0 ; k<nlabelsB[i] ; k++) {
            sprintf(tmp, "%d ", labelsB[i][k]);
            strcat(out, tmp);
        }
        sprintf(tmp, "\n");
        strcat(out, tmp);
      }
      if (format == LIGHT)
      {
        sprintf(tmp, "%d - %d\n", nlabelsA[i], nlabelsB[i]);
        strcat(out, tmp);
      }
    }

    if(format == LIGHT) {
        int NA, NB;
        NA = sizeA-1;
        NB = sizeB-1;
        for(i=0; i<LAB; i++) NA -= nlabelsA[i];
        for(i=0; i<LAB; i++) NB -= nlabelsB[i];
        for(i=1; i<sizeA; i++)
          if (larray[0][0][i] == (unsigned short int) 0) NA--;
        for(j=1; j<sizeB; j++)
          if (larray[0][j][0] == (unsigned short int) 0) NB--;
        for (i=0 ; i<NA ; i++) {
            sprintf(tmp, "1 - 0\n");
            strcat(out, tmp);
        }
        for (i=0 ; i<NB ; i++) {
            sprintf(tmp, "0 - 1\n");
            strcat(out, tmp);
        }
    }
    if (format == COMPLETE) {
        for (i=1 ; i<sizeA ; i++) {
            if (larray[0][0][i] == (unsigned short int) 0) continue;
            for (j=1 ; j<sizeB ; j++)
                if (larray[0][j][i] != 0) break;
            if (j != sizeB) continue;
            sprintf(tmp, "%d - \n", i);
            strcat(out, tmp);
        }
        for (j=1 ; j<sizeB ; j++) {
            if (larray[0][j][0] == (unsigned short int) 0) continue;
            for (i=1 ; i<sizeA ; i++)
                if (larray[0][j][i] != 0) break;
            if (i != sizeA) continue;
            sprintf(tmp, " - %d\n", j);
            strcat(out, tmp);
        }
    }
    for (i=0 ; i<LAB ; i++)
    {
        free(labelsA[i]);
        free(labelsB[i]);
    }
    free(labelsA);
    free(labelsB);
    free(nlabelsA);
    free(nlabelsB);

    *output = out;
}


extern void VT_WriteOverlapInterpretationVoxel2(char **output, long long  *volume_bijA, long long  *volume_bijB, double *proportion_recouvrementA, double *proportion_recouvrementB,
                                          long long  bijectionA, long long  bijectionB, long long  volume_totalA, long long  volume_totalB)
{
    char *program="VT_WriteOverlapInterpretationVoxel2";
    double moyenneA=0, devA=0;
    double moyenneB=0, devB=0;
    long long  sum_volume_bijA=0, sum_volume_bijB=0;
    int i;
    char *out;
    char tmp[100];

    out = malloc(100 * sizeof(char) * 5); /* 100 : size max per line ; 5 : number of lines.
                                           */
    if (out == NULL)
    {
        fprintf(stderr, "%s: unable to allocate memory\n", program);
        return;
    }
    out[0] = '\0';

    for (i=0 ; i<bijectionA ; i++) {
        sum_volume_bijA += volume_bijA[i];
        moyenneA += proportion_recouvrementA[i];
    }
    moyenneA /= bijectionA;
    for (i=0 ; i<bijectionA ; i++) {
        devA += (proportion_recouvrementA[i]-moyenneA)*(proportion_recouvrementA[i]-moyenneA);
    }
    devA /= bijectionA;
    devA = sqrt(devA);

    for (i=0 ; i<bijectionB ; i++) {
        sum_volume_bijB += volume_bijB[i];
        moyenneB += proportion_recouvrementB[i];
    }
    moyenneB /= bijectionB;
    for (i=0 ; i<bijectionB ; i++) {
        devB += (proportion_recouvrementB[i]-moyenneB)*(proportion_recouvrementB[i]-moyenneB);
    }
    devB /= bijectionB;
    devB = sqrt(devB);


    sprintf(tmp, "\nStatistiques voxelliques :\n\n");
    strcat(out, tmp);
    sprintf(tmp, "in :\n taux de recouvrement moyen : %f\n ecart-type : %f\n", moyenneA, devA);
    strcat(out, tmp);

    sprintf(tmp, "ext :\n taux de recouvrement moyen : %f\n ecart-type : %f\n", moyenneB, devB);
    strcat(out, tmp);

    sprintf(tmp, "in :\n volume en bijection / volume total : %f\n", ((double)sum_volume_bijA)/volume_totalA);
    strcat(out, tmp);

    sprintf(tmp, "ext :\n volume en bijection / volume total : %f\n", ((double)sum_volume_bijB)/volume_totalB);
    strcat(out, tmp);


    *output = out;
}

extern void VT_WriteOverlapInterpretationCell(char **output, int nLabelsA, int nLabelsB, int bijectionA, int bijectionB,
                                              int *sous_segA, int *sous_segB, int *sur_segA, int *sur_segB, int Afond, int Bfond, int divergentA, int divergentB)
{
    char *program = "VT_WriteOverlapInterpretationCell";
    char *out;
    char tmp[400];


    out = malloc(200 * sizeof(char) * 5); /* 100 : size max per line ; 5 : number of lines.
                                           */
    if (out == NULL)
    {
        fprintf(stderr, "%s: unable to allocate memory\n", program);
        return;
    }
    out[0] = '\0';

    sprintf(tmp, "Statistiques sur les cellules :\n\n");
    strcat(out, tmp);

    sprintf(tmp, "in :\n #cells : %d\n bijections : %f%%\n sous-seg : %f%%\n sur-seg : %f%%\n dans fond : %f%%\n misc. : %f%%\n\n", nLabelsA, ((float)bijectionA)*100/nLabelsA,
            ((float)sous_segA[0]+ sous_segA[1]+ sous_segA[2]+ sous_segA[3]+ sous_segA[4])*100/nLabelsA,
            ((float)sur_segA[0]+ sur_segA[1]+ sur_segA[2]+ sur_segA[3]+ sur_segA[4])*100/nLabelsA,
            ((float)Afond)*100/nLabelsA, ((float)divergentA)*100/nLabelsA);
    strcat(out, tmp);

    sprintf(tmp, "ext :\n #cells : %d\n bijections : %f%%\n sous-seg : %f%%\n sur-seg : %f%%\n dans fond : %f%%\n misc. : %f%%\n\n", nLabelsB, ((float)bijectionB)*100/nLabelsB,
            ((float)sous_segB[0]+ sous_segB[1]+ sous_segB[2]+ sous_segB[3]+ sous_segB[4])*100/nLabelsB,
            ((float)sur_segB[0]+ sur_segB[1]+ sur_segB[2]+ sur_segB[3]+ sur_segB[4])*100/nLabelsB,
            ((float)Bfond)*100/nLabelsB, ((float)divergentB)*100/nLabelsB);
    strcat(out, tmp);
/*
    sprintf(tmp, "in :\n #cells : %d\n bijections : %d\n sous-seg : %d %d %d %d %d\n sur-seg : %d %d %d %d %d\n dans fond : %d\n misc. : %d\n\n", nLabelsA, bijectionA,
            sous_segA[0], sous_segA[1], sous_segA[2], sous_segA[3], sous_segA[4],
            sur_segA[0], sur_segA[1], sur_segA[2], sur_segA[3], sur_segA[4],
            Afond, divergentA);
    strcat(out, tmp);

    sprintf(tmp, "ext :\n #cells : %d\n bijections : %d\n sous-seg : %d %d %d %d %d\n sur-seg : %d %d %d %d %d\n dans fond : %d\n misc. : %d\n\n", nLabelsB, bijectionB,
            sous_segB[0], sous_segB[1], sous_segB[2], sous_segB[3], sous_segB[4],
            sur_segB[0], sur_segB[1], sur_segB[2], sur_segB[3], sur_segB[4],
            Bfond, divergentB);
    strcat(out, tmp);
*/

    *output = out;
}

extern void VT_WriteOverlapInterpretationVoxel(char **output, long long  vox_bijectionA, long long  vox_bijectionB,
                                              long long  *vox_sous_segA, long long  *vox_sous_segB, long long  *vox_sur_segA, long long  *vox_sur_segB, long long  vox_Afond, long long  vox_Bfond, long long  vox_divergentA, long long  vox_divergentB)
{
    char *program = "VT_WriteOverlapInterpretationVoxel";
    char *out;
    char tmp[500];


    out = malloc(200 * sizeof(char) * 5); /* 100 : size max per line ; 5 : number of lines.
                                           */
    if (out == NULL)
    {
        fprintf(stderr, "%s: unable to allocate memory\n", program);
        return;
    }
    out[0] = '\0';

    long long nVoxTotA=0;
    long long nVoxTotB=0;
    int i;

    for (i=0;i<ORDERMAX;i++) {
        nVoxTotA+=vox_sous_segA[i];
        nVoxTotB+=vox_sous_segB[i];
        nVoxTotA+=vox_sur_segA[i];
        nVoxTotB+=vox_sur_segB[i];
    }
    nVoxTotA+=vox_bijectionA;
    nVoxTotB+=vox_bijectionB;
    nVoxTotA+=vox_Afond;
    nVoxTotB+=vox_Bfond;
    nVoxTotA+=vox_divergentA;
    nVoxTotB+=vox_divergentB;

    sprintf(tmp, "Statistiques sur les voxels :\n\n");
    strcat(out, tmp);

    sprintf(tmp, "in :\n bijections : %f%%\n sous-seg : %f%%\n sur-seg : %f%%\n dans fond : %f%%\n misc. : %f%%\n\n", ((float)vox_bijectionA)*100/nVoxTotA,
            ((float)(vox_sous_segA[0]+ vox_sous_segA[1]+ vox_sous_segA[2]+ vox_sous_segA[3]+ vox_sous_segA[4]))*100/nVoxTotA,
            ((float)(vox_sur_segA[0]+ vox_sur_segA[1]+ vox_sur_segA[2]+ vox_sur_segA[3]+ vox_sur_segA[4]))*100/nVoxTotA,
            ((float)vox_Afond)*100/nVoxTotA, ((float)vox_divergentA)*100/nVoxTotA);

    /*fprintf(stdout, "A:\nbijections: %lli / %lli = %f %%\n", vox_bijectionA, nVoxTotA, ((float) vox_bijectionA)/nVoxTotA);
    for (i=0 ; i<ORDERMAX; i++)
        fprintf(stdout, "sous-seg[%d]: %lli / %lli = %f %%\n", i, vox_sous_segA[i], nVoxTotA, ((float) vox_sous_segA[i])/nVoxTotA);
    for (i=0 ; i<ORDERMAX; i++)
        fprintf(stdout, "sur-seg[%d]: %lli / %lli = %f %%\n", i, vox_sur_segA[i], nVoxTotA, ((float) vox_sur_segA[i])/nVoxTotA);
    fprintf(stdout, "dans fond: %lli / %lli = %f %%\n", vox_Afond, nVoxTotA, ((float) vox_Afond)/nVoxTotA);
    fprintf(stdout, "misc.: %lli / %lli = %f %%\n", vox_divergentA, nVoxTotA, ((float) vox_divergentA)/nVoxTotA);

    fprintf(stdout, "\nB:\nbijections: %lli / %lli = %f %%\n", vox_bijectionB, nVoxTotB, ((float) vox_bijectionB)/nVoxTotB);
    for (i=0 ; i<ORDERMAX; i++)
        fprintf(stdout, "sous-seg[%d]: %lli / %lli = %f %%\n", i, vox_sous_segB[i], nVoxTotB, ((float) vox_sous_segB[i])/nVoxTotB);
    for (i=0 ; i<ORDERMAX; i++)
        fprintf(stdout, "sur-seg[%d]: %lli / %lli = %f %%\n", i, vox_sur_segB[i], nVoxTotB, ((float) vox_sur_segB[i])/nVoxTotB);
    fprintf(stdout, "dans fond: %lli / %lli = %f %%\n", vox_Bfond, nVoxTotB, ((float) vox_Bfond)/nVoxTotB);
    fprintf(stdout, "misc.: %lli / %lli = %f %%\n", vox_divergentB, nVoxTotB, ((float) vox_divergentB)/nVoxTotB);
    */

    strcat(out, tmp);

    sprintf(tmp, "ext :\n bijections : %f%%\n sous-seg : %f%%\n sur-seg : %f%%\n dans fond : %f%%\n misc. : %f%%\n\n", ((float)vox_bijectionB)*100/nVoxTotB,
            ((float)(vox_sous_segB[0]+ vox_sous_segB[1]+ vox_sous_segB[2]+ vox_sous_segB[3]+ vox_sous_segB[4]))*100/nVoxTotB,
            ((float)(vox_sur_segB[0]+ vox_sur_segB[1]+ vox_sur_segB[2]+ vox_sur_segB[3]+ vox_sur_segB[4]))*100/nVoxTotB,
            ((float)vox_Bfond)*100/nVoxTotB, ((float)vox_divergentB)*100/nVoxTotB);
    strcat(out, tmp);


    *output = out;
}

extern int VT_OverlapInterpretation(FILE *fp, vt_image *imA, vt_image *imB, int flagIn, int flagExt, vt_image* imAout, vt_image *imBout, int flagInOut, int flagExtOut,
                                    int *_nLabelsA, int *_nLabelsB, int *_Afond, int *_Bfond, long long  *_volume_totalA, long long  *_volume_totalB,
                                    int *_bijectionA, int *_bijectionB, int *sous_segA, int *sous_segB, int *sur_segA, int *sur_segB,
                                    int *_divergentA, int *_divergentB,
                                    long long *_vox_Afond, long long *_vox_Bfond,
                                    long long *_vox_bijectionA, long long *_vox_bijectionB, long long *vox_sous_segA, long long *vox_sous_segB, long long *vox_sur_segA, long long *vox_sur_segB,
                                    long long *_vox_divergentA, long long *_vox_divergentB,
                                    long long **_volume_bijA, long long **_volume_bijB, double **_proportion_recouvrementA, double **_proportion_recouvrementB)
{
    char *proc="VT_OverlapInterpretation";


    int i;

    int bijectionA=0;
    int bijectionB=0;
    int divergentA=0;
    int divergentB=0;
    int Afond=0;
    int Bfond=0;

    for (i=0 ; i<ORDERMAX ; i++) {
        sous_segA[i]=0;
        sous_segB[i]=0;
        sur_segA[i]=0;
        sur_segB[i]=0;
    }

    long long vox_bijectionA=0;
    long long vox_bijectionB=0;
    long long vox_divergentA=0;
    long long vox_divergentB=0;
    long long vox_Afond=0;
    long long vox_Bfond=0;

    for (i=0 ; i<ORDERMAX ; i++) {
        vox_sous_segA[i]=0;
        vox_sous_segB[i]=0;
        vox_sur_segA[i]=0;
        vox_sur_segB[i]=0;
    }

    long long *volume_bijA=NULL;
    long long *volume_bijB=NULL;

    double *proportion_recouvrementA=NULL;
    double *proportion_recouvrementB=NULL;

    long long volume_totalA=0;
    long long volume_totalB=0;

    int nLabelsA=0, nLabelsB=0;

    int *bboxXA[2];
    int *bboxYA[2];
    int *bboxZA[2];
    int *bboxXB[2];
    int *bboxYB[2];
    int *bboxZB[2];
    long long *volumesA;
    long long *volumesB;

    ssize_t read;
    char *line = NULL;
    size_t n = 0;

    if (flagIn)
    {
        bbox(imA, bboxXA,bboxYA,bboxZA,&volumesA, &nLabelsA);
        volume_bijA = malloc(nLabelsA * sizeof(long long));
        if (volume_bijA == NULL) {
            fprintf(stderr, "%s: unable to allocate the array of volumes for image in\n", proc);
            return(0);
        }
        proportion_recouvrementA = malloc(nLabelsA * sizeof(double));
        if (proportion_recouvrementA == NULL) {
            free(volume_bijA); volume_bijA = NULL;
            free(volumesA); volumesA=NULL;
            fprintf(stderr, "%s: unable to allocate the array of overlapping proportions for image in\n", proc);
            return(0);
        }
    }

    if (flagExt)
    {
        bbox(imB, bboxXB,bboxYB,bboxZB,&volumesB, &nLabelsB);
        volume_bijB = malloc(nLabelsB * sizeof(long long));
        if ( volume_bijB == NULL ) {
          if(volume_bijA) free(volume_bijA);
          if(volumesA) free(volumesA);
          if(proportion_recouvrementA) free(proportion_recouvrementA);
          free(volumesB); volumesB=NULL;
          fprintf(stderr, "%s: unable to allocate the array of volumes for image ext\n", proc);
          return(0);
        }
        proportion_recouvrementB = malloc(nLabelsB * sizeof(double));
        if ( proportion_recouvrementB == NULL ) {
          if(volume_bijA) { free(volume_bijA); volume_bijA = NULL; }
          if(volumesA) free(volumesA);
          if(proportion_recouvrementA) { free(proportion_recouvrementA); proportion_recouvrementA = NULL; }
          free(volumesB);
          free(volume_bijB); volume_bijB = NULL;
          fprintf(stderr, "%s: unable to allocate the array of overlapping proportions for image ext\n", proc);
          return(0);
        }
    }



    if (flagInOut)
    {
        VT_InitFromImage(imAout, imA, "imInOut", (int)UCHAR);
        if (VT_AllocImage(imAout) != 1)
        {
            if(volume_bijA) free(volume_bijA);
            if(proportion_recouvrementA) free(proportion_recouvrementA);
            if(volume_bijB) free(volume_bijB);
            if(proportion_recouvrementB) free(proportion_recouvrementB);
            if (volumesA) free(volumesA);
            if (volumesB) free(volumesB);
            fprintf(stderr, "%s: unable to allocate the image-in-out\n", proc);
            return(0);
        }
    }
    if (flagExtOut)
    {
        /* Init
         */
        VT_InitFromImage(imBout, imB, "imExtOut", (int)UCHAR);
        if (VT_AllocImage(imBout) != 1)
        {
            if (flagInOut) VT_FreeImage(imAout);
            if(volume_bijA) free(volume_bijA);
            if(proportion_recouvrementA) free(proportion_recouvrementA);
            if(volume_bijB) free(volume_bijB);
            if(proportion_recouvrementB) free(proportion_recouvrementB);
            if (volumesA) free(volumesA);
            if (volumesB) free(volumesB);
            fprintf(stderr, "%s: unable to allocate the image-ext-out\n", proc);
            return(0);
        }
    }

    int nline = 0;

    while (( read = getline(&line, &n, fp )) != -1) {
        ++nline;

        int flagAB=0; /* 0 : seg A ; 1 : seg B
                       */

        int lab = 0;
        int cptA = 0;
        int cptB = 0;
        int *listeA = malloc(sizeof(int)*read);
        int *listeB = malloc(sizeof(int)*read);

        i=0;
        while(i < read)
        {
            if (line[i]<'0' || line[i]>'9') {
                if(line[i] == '-') {
                    flagAB = 1;
                }
                else {
                    if (lab == 0) { i++; continue; }
                    if (flagAB == 0) listeA[cptA++] = lab;
                    else listeB[cptB++] = lab;
                }
                lab = 0;
            }
            else
            {   /* Construction du label caractere par caractere (badass loop)
                 */
                lab *= 10;
                lab += (int) (line[i] - '0');
            }
            if ( line[i] == '-')
                flagAB = 1;

            i++;
        }



        /* Statistiques */

        int labelA=0;
        int labelB=0;
        if (cptA == 1 && cptB == 1) /* Bijection
                                     */
        {
            bijectionA+=cptA;
            bijectionB+=cptB;


            if (flagIn && flagExt ) {

                  vox_bijectionA+=volumesA[listeA[0]];
                  vox_bijectionB+=volumesB[listeB[0]];

                  long long volA, volB;
                  double recA, recB;
                  recouvrement(imA, imB,
                               bboxXA[0][listeA[0]], bboxXA[1][listeA[0]],
                               bboxYA[0][listeA[0]], bboxYA[1][listeA[0]],
                               bboxZA[0][listeA[0]], bboxZA[1][listeA[0]],
                               bboxXB[0][listeB[0]], bboxXB[1][listeB[0]],
                               bboxYB[0][listeB[0]], bboxYB[1][listeB[0]],
                               bboxZB[0][listeB[0]], bboxZB[1][listeB[0]],
                               listeA[0], listeB[0],
                               &volA, &volB, &recA, &recB );
                  volume_bijA[bijectionA-1] = volA;
                  volume_bijB[bijectionB-1] = volB;
                  proportion_recouvrementA[bijectionA-1] = recA;
                  proportion_recouvrementB[bijectionB-1] = recB;
                  labelA = (int) (100*recA);
                  labelB = (int) (100*recB);
                  volume_totalA+=volA;
                  volume_totalB+=volB;
            }
            else {
                labelA=100;
                labelB=100;
            }
        }


        if (cptA == 1 && cptB > 1) /* Sous seg A/B
                                    */
        {

            int k=(cptB-2>=ORDERMAX) ? ORDERMAX-1 : cptB-2;
            sous_segA[k]+=cptA;
            sur_segB[k] +=cptB;
            /* labelA=100+cptB;
             * labelB=200+cptB;
             */
            labelA=100+((cptB-2>=ORDERMAX) ? ORDERMAX : cptB);
            labelB=200+((cptB-2>=ORDERMAX) ? ORDERMAX : cptB);

            if(flagIn && flagExt) {
              volume_totalA+=volume(imA, bboxXA[0][listeA[0]], bboxXA[1][listeA[0]], bboxYA[0][listeA[0]], bboxYA[1][listeA[0]], bboxZA[0][listeA[0]], bboxZA[1][listeA[0]],
                    listeA[0]);
              vox_sous_segA[k]+=volumesA[listeA[0]];
              int l;
              for (l = 0 ; l<cptB ; l++) {
                vox_sur_segB[k]+=volumesB[listeB[l]];
                volume_totalB+=volume(imB, bboxXB[0][listeB[l]], bboxXB[1][listeB[l]], bboxYB[0][listeB[l]], bboxYB[1][listeB[l]], bboxZB[0][listeB[l]], bboxZB[1][listeB[l]],
                        listeB[l]);
              }
            }

        }



        if (cptA > 1 && cptB == 1) /* Sous seg B/A
                                    */
        {
            int k=(cptA-2>=ORDERMAX) ? ORDERMAX-1 : cptA-2;
            sur_segA[k] +=cptA;
            sous_segB[k]+=cptB;
            labelA=200+((cptA-2>=ORDERMAX) ? ORDERMAX : cptA);
            labelB=100+((cptA-2>=ORDERMAX) ? ORDERMAX : cptA);

            if(flagIn && flagExt) {
              int l;
              for (l = 0 ; l<cptA ; l++) {
                volume_totalA+=volume(imA, bboxXA[0][listeA[l]], bboxXA[1][listeA[l]], bboxYA[0][listeA[l]], bboxYA[1][listeA[l]], bboxZA[0][listeA[l]], bboxZA[1][listeA[l]],
                        listeA[l]);
                vox_sur_segA[k]+=volumesA[listeA[l]];
              }
              vox_sous_segB[k]+=volumesB[listeB[0]];
              volume_totalB+=volume(imB, bboxXB[0][listeB[0]], bboxXB[1][listeB[0]], bboxYB[0][listeB[0]], bboxYB[1][listeB[0]], bboxZB[0][listeB[0]], bboxZB[1][listeB[0]],
                    listeB[0]);
            }
        }
        if (cptA > 1 && cptB > 1) /* Divergence
                                   */
        {
            int k;
            divergentA+=cptA;
            divergentB+=cptB;
            labelA=250;
            labelB=250;
            if(flagIn && flagExt) {
              for (k = 0 ; k<cptA ; k++) {
                  vox_divergentA+=volumesA[listeA[k]];
                  volume_totalA+=volume(imA, bboxXA[0][listeA[k]], bboxXA[1][listeA[k]], bboxYA[0][listeA[k]], bboxYA[1][listeA[k]], bboxZA[0][listeA[k]], bboxZA[1][listeA[k]],
                        listeA[k]);
              }
              for (k = 0 ; k<cptB ; k++){
                vox_divergentB+=volumesB[listeB[k]];
                volume_totalB+=volume(imB, bboxXB[0][listeB[k]], bboxXB[1][listeB[k]], bboxYB[0][listeB[k]], bboxYB[1][listeB[k]], bboxZB[0][listeB[k]], bboxZB[1][listeB[k]],
                        listeB[k]);
              }
            }
        }

        if (cptA == 1 && cptB == 0) /* A dans fond
                                     */
        {
            Afond+=cptA;
            labelA=255;
            if (flagIn && flagExt) {
                volume_totalA+=volume(imA, bboxXA[0][listeA[0]], bboxXA[1][listeA[0]], bboxYA[0][listeA[0]], bboxYA[1][listeA[0]], bboxZA[0][listeA[0]], bboxZA[1][listeA[0]],
                        listeA[0]);
                vox_Afond+=volumesA[listeA[0]];
            }
        }
        if (cptA == 0 && cptB == 1) /* B dans fond
                                     */
        {
            Bfond+=cptB;
            labelB=255;
            if (flagIn && flagExt){
                vox_Bfond+=volumesB[listeB[0]];
                volume_totalB+=volume(imB, bboxXB[0][listeB[0]], bboxXB[1][listeB[0]], bboxYB[0][listeB[0]], bboxYB[1][listeB[0]], bboxZB[0][listeB[0]], bboxZB[1][listeB[0]],
                      listeB[0]);
            }
        }


        if (flagInOut) {
          int k;
          for (k = 0 ; k<cptA ; k++)
            setLabel(imA, imAout, bboxXA[0][listeA[k]], bboxXA[1][listeA[k]], bboxYA[0][listeA[k]], bboxYA[1][listeA[k]], bboxZA[0][listeA[k]], bboxZA[1][listeA[k]],
                     listeA[k], labelA);
        }
        if (flagExtOut) {
          int k;
          for (k = 0 ; k<cptB ; k++)
            setLabel(imB, imBout, bboxXB[0][listeB[k]], bboxXB[1][listeB[k]], bboxYB[0][listeB[k]], bboxYB[1][listeB[k]], bboxZB[0][listeB[k]], bboxZB[1][listeB[k]],
                     listeB[k], labelB);
        }


        /* Liberation memoire */

        free(listeA);
        free(listeB);
    }

    /* free bbox */
    if (flagIn) {
      for (i=0; i<2; i++) {
       free(bboxXA[i]);
       free(bboxYA[i]);
       free(bboxZA[i]);
      }
      free(volumesA);
    }

    if (flagExt){
      for (i=0; i<2; i++) {
        free(bboxXB[i]);
        free(bboxYB[i]);
        free(bboxZB[i]);
      }
      free(volumesB);
    }

    *_nLabelsA = nLabelsA;
    *_nLabelsB = nLabelsB;
    *_Afond = Afond;
    *_Bfond = Bfond;
    *_volume_totalA = volume_totalA;
    *_volume_totalB = volume_totalB;
    *_bijectionA = bijectionA;
    *_bijectionB = bijectionB;
    *_divergentA = divergentA;
    *_divergentB = divergentB;
    *_volume_bijA = volume_bijA;
    *_volume_bijB = volume_bijB;
    *_proportion_recouvrementA = proportion_recouvrementA;
    *_proportion_recouvrementB = proportion_recouvrementB;

    *_vox_bijectionA = vox_bijectionA;
    *_vox_bijectionB = vox_bijectionB;
    *_vox_divergentA = vox_divergentA;
    *_vox_divergentB = vox_divergentB;
    *_vox_Afond = vox_Afond;
    *_vox_Bfond = vox_Bfond;

    return(1);
}
