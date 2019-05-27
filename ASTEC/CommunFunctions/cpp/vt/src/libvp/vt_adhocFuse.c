/*************************************************************************
 * vt_adhocFuse.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Ven  2 fév 2018 16:40:20 CET
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <chamferdistance.h>
#include <chamferdistance-mask.h>
#include <histogram.h>
#include <linearFiltering.h>
#include <vtmalloc.h>

#include <ccparameters.h>

#include <vt_adhocFuse.h>



static int _verbose_ = 1;





/*************************************************************************
 *
 *
 *
 *************************************************************************/



static int _getPercentileFromCumulative( typeHistogram *cumul, float percentile )
{
  char *proc = "_getPercentileFromCumulative";
  s32 *theIndex;
  r32 *theHisto;
  int i;


  if ( cumul->xaxis.typeIndex != SINT ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: such histogram index type not handled yet\n", proc );
    return( -33333 );
  }

  if ( cumul->typeHisto != FLOAT ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: such histogram data type not handled yet\n", proc );
    return( -33333 );
  }

  theIndex = (s32*)cumul->xaxis.index;
  theHisto = (r32*)cumul->data;

  for ( i=0; i<cumul->xaxis.dim; i++ ) {
    if ( i == 0 && percentile < theHisto[i] ) {
        return( theIndex[i] );
    }
    else if ( (0 < i && i < cumul->xaxis.dim-1)
              && theHisto[i] <= percentile && percentile < theHisto[i+1] ) {
      if ( percentile - theHisto[i] <= theHisto[i+1] - percentile )
        return( theIndex[i] );
      else
        return( theIndex[i+1] );
    }
    else if ( i == cumul->xaxis.dim-1 ) {
      return( theIndex[i] );
    }
  }

  return( -33333 );
}







static int _computeCellBasedExtremum( vt_image *intImage,
                                      vt_image *segImage,
                                      vt_image *resImage,
                                      typeHistogram *histo,
                                      typeHistogram *pdf,
                                      typeHistogram *cumul,
                                      float percentile,
                                      enumExtremaMethod method )

{
  char *proc = "_computeCellBasedExtremum";
  typeParameter *theCC = (typeParameter *)NULL;
  int nCC = 0;
  size_t x, y, z;

  s32 min = histo->xaxis.min.val_s32;
  s32 max = histo->xaxis.max.val_s32;
  u32* theHisto = (u32*)histo->data;

  size_t ptmin[3];
  size_t ptmax[3];
  int i, c;

  int value;





  if ( segImage == (vt_image*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: NULL input segmentation image, nothing to do\n", proc );
    return( 1 );
  }

  if ( resImage == (vt_image*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: NULL output image, nothing to do\n", proc );
    return( 1 );
  }

  if ( histo->xaxis.typeIndex != SINT ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: such histogram index type not handled yet\n", proc );
    return( -1 );
  }

  if ( histo->typeHisto != UINT ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: such histogram type not handled yet\n", proc );
    return( -1 );
  }



  /* processing numbered cells
   */
  theCC = ComputeParameterFromLabels( segImage, &nCC );
  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "%s: found %d connect components in '%s'\n",
      proc, nCC, segImage->name );
  }





  switch( segImage->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such segmentation image type not handled yet\n", proc );
    return( -1 );
  case USHORT :
    {
      u16 ***segBuf = (u16***)segImage->array;
      switch( intImage->type ) {
      default :
        if ( _verbose_ )
          fprintf( stderr, "%s: such intensity image type not handled yet\n", proc );
        return( -1 );
      case  USHORT :
        {
          u16 ***intBuf = (u16***)intImage->array;
          u16 ***resBuf = (u16***)resImage->array;

          for ( c=1; c<=nCC; c++ ) {

            if ( theCC[c].volume <= 0 ) continue;

            zeroHistogram( histo );

            /* bounding box
             */
            for ( i=0; i<3; i++ ) {
              ptmin[i] = theCC[c].ptmin[i];
              ptmax[i] = theCC[c].ptmax[i];
            }


            switch ( method ) {
            default :
              if ( _verbose_ )
                fprintf( stderr, "%s: such method not handled yet\n", proc );
              return( -1 );

            case _CELL_ :

                /* whole cell
                 */
                for ( z=ptmin[2]; z<=ptmax[2]; z++ )
                for ( y=ptmin[1]; y<=ptmax[1]; y++ )
                for ( x=ptmin[0]; x<=ptmax[0]; x++ ) {
                  if ( segBuf[z][y][x] == c ) {
                      if ( intBuf[z][y][x] < min ) theHisto[0] ++;
                      else if ( intBuf[z][y][x] > max ) theHisto[histo->xaxis.dim-1] ++;
                      else theHisto[ (int)(((float)(intBuf[z][y][x])-min)/histo->xaxis.binlength + 0.5) ] ++;
                  }
                }
                break;

            case _CELL_BORDER_ :
                /* bounding box augmentee de 1
                 */
                for ( i=0; i<3; i++ ) {
                  if ( ptmin[i] > 0 ) ptmin[i] -= 1;
                  if ( i == 0 && ptmax[i] < segImage->dim.x-1 ) ptmax[i] += 1;
                  if ( i == 1 && ptmax[i] < segImage->dim.y-1 ) ptmax[i] += 1;
                  if ( i == 2 && ptmax[i] < segImage->dim.z-1 ) ptmax[i] += 1;
                }

                /* borders
                 */
                for ( z=ptmin[2]; z<=ptmax[2]; z++ )
                for ( y=ptmin[1]; y<=ptmax[1]; y++ )
                for ( x=ptmin[0]; x<=ptmax[0]; x++ ) {
                  if ( segBuf[z][y][x] == c ) {
                    if ( (x > 0 && segBuf[z][y][x-1] != c)
                         || (x < segImage->dim.x-1 && segBuf[z][y][x+1] != c)
                         || (y > 0 && segBuf[z][y-1][x] != c)
                         || (y < segImage->dim.y-1 && segBuf[z][y+1][x] != c)
                         || (z > 0 && segBuf[z-1][y][x] != c)
                         || (z < segImage->dim.z-1 && segBuf[z+1][y][x] != c) ) {
                      if ( intBuf[z][y][x] < min ) theHisto[0] ++;
                      else if ( intBuf[z][y][x] > max ) theHisto[histo->xaxis.dim-1] ++;
                      else theHisto[ (int)(((float)(intBuf[z][y][x])-min)/histo->xaxis.binlength + 0.5) ] ++;
                    }
                  }
                  else if ( (x > 0 && segBuf[z][y][x-1] == c)
                            || (x < segImage->dim.x-1 && segBuf[z][y][x+1] == c)
                            || (y > 0 && segBuf[z][y-1][x] != c)
                            || (y < segImage->dim.y-1 && segBuf[z][y+1][x] == c)
                            || (z > 0 && segBuf[z-1][y][x] != c)
                            || (z < segImage->dim.z-1 && segBuf[z+1][y][x] == c) ) {
                    if ( intBuf[z][y][x] < min ) theHisto[0] ++;
                    else if ( intBuf[z][y][x] > max ) theHisto[histo->xaxis.dim-1] ++;
                    else theHisto[ (int)(((float)(intBuf[z][y][x])-min)/histo->xaxis.binlength + 0.5) ] ++;
                  }
                }
                break;

            case _CELL_INTERIOR_ :

                /* inner cell -> minimum
                 */
                for ( z=ptmin[2]; z<=ptmax[2]; z++ )
                for ( y=ptmin[1]; y<=ptmax[1]; y++ )
                for ( x=ptmin[0]; x<=ptmax[0]; x++ ) {
                  if ( segBuf[z][y][x] == c ) {
                    if ( (x > 0 && segBuf[z][y][x-1] != c)
                         || (x < segImage->dim.x-1 && segBuf[z][y][x+1] != c)
                         || (y > 0 && segBuf[z][y-1][x] != c)
                         || (y < segImage->dim.y-1 && segBuf[z][y+1][x] != c)
                         || (z > 0 && segBuf[z-1][y][x] != c)
                         || (z < segImage->dim.z-1 && segBuf[z+1][y][x] != c) ) {
                      continue;
                    }
                    if ( intBuf[z][y][x] < min ) theHisto[0] ++;
                    else if ( intBuf[z][y][x] > max ) theHisto[histo->xaxis.dim-1] ++;
                    else theHisto[ (int)(((float)(intBuf[z][y][x])-min)/histo->xaxis.binlength + 0.5) ] ++;
                  }
                }
                break;
            }

            if ( pdfHistogram( pdf, histo ) != 1 ) {
              if ( _verbose_ )
                fprintf( stderr, "%s: error when computing pdf histogram\n", proc );
              return( -1 );
            }
            if ( cumulative1DHistogram( cumul, pdf) != 1 ) {
                if ( _verbose_ )
                  fprintf( stderr, "%s: error when computing cumulative histogram\n", proc );
                return( -1 );
            }

            value = _getPercentileFromCumulative( cumul, percentile );

            for ( z=ptmin[2]; z<=ptmax[2]; z++ )
            for ( y=ptmin[1]; y<=ptmax[1]; y++ )
            for ( x=ptmin[0]; x<=ptmax[0]; x++ ) {
              if ( segBuf[z][y][x] != c ) continue;
              resBuf[z][y][x] = value;
            }

          }
        }
        /* end of intImage->type = USHORT
         */
        break;
      }
    }
    break;
  }

  vtfree( theCC );

  return( 1 );
}





static int _computeVoxelBasedExtremum( vt_image *intImage,
                          vt_image *mskImage,
                          vt_image *resImage )
{
  char *proc = "_computeVoxelBasedExtremum";
  size_t i, size;
  typeChamferMask theMask;
  int chamfercoef_dim = 0;
  int chamfercoef_size = 5;
  int chamfercoef_max = 5;
  double chamfercoef_voxel_size[3];
  u16 *theDist = (u16 *)NULL;
  int theDim[3];


  /* initialization
   */

  theDim[0] = intImage->dim.x;
  theDim[1] = intImage->dim.y;
  theDim[2] = intImage->dim.z;
  size = intImage->dim.x * intImage->dim.y * intImage->dim.z;


#define _VOXELINTTYPE( TYPE ) {        \
  TYPE *intBuf = (TYPE*)intImage->buf; \
  TYPE *resBuf = (TYPE*)resImage->buf; \
  for ( i=0; i<size; i++, mskBuf++, intBuf++, resBuf++ ) { \
    if ( *mskBuf == 0 ) {              \
      *resBuf = 0;                     \
      continue;                        \
    }                                  \
    *resBuf = *intBuf;                 \
  }                                    \
}



#define _VOXELMSKTYPE( TYPE ) {        \
  TYPE *mskBuf = (TYPE*)mskImage->buf; \
  switch( intImage->type ) {           \
  default :                            \
    if ( _verbose_ )                   \
      fprintf( stderr, "%s: such reconstructed image type not handled yet\n", proc ); \
    return( -1 );                      \
  case UCHAR :                         \
    _VOXELINTTYPE( u8 );               \
    break;                             \
  case USHORT :                        \
    _VOXELINTTYPE( u16 );              \
    break;                             \
  }                                    \
}


  switch( mskImage->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such mask image type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    _VOXELMSKTYPE( u8 );
    break;
  case USHORT :
    _VOXELMSKTYPE( u16 );
    break;
  }

  /* chamfer mask computation
   */
  initChamferMask( &theMask );
  chamfercoef_dim = ( intImage->dim.z > 1 ) ? 3 : 2;
  chamfercoef_voxel_size[0] = 1.0;
  chamfercoef_voxel_size[1] = intImage->siz.y / intImage->siz.x;
  chamfercoef_voxel_size[2] = intImage->siz.z / intImage->siz.x;

  if ( buildChamferMask( chamfercoef_voxel_size, chamfercoef_dim, chamfercoef_size,
                         chamfercoef_max, _UNDEFINED_DISTANCE_, &theMask ) != 1 ) {
    freeChamferMask( &theMask );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when computing chamfer mask\n", proc );
    return( -1 );
  }


  /* distance buffer allocation
   */
  theDist = (unsigned short int*)vtmalloc( size * sizeof(u16), "theDist", proc );
  if ( theDist == (unsigned short int*)NULL ) {
    freeChamferMask( &theMask );
    if ( _verbose_ )
      fprintf( stderr, "%s: distance buffer allocation error\n", proc );
    return( -1 );
  }

  /* skiz
   */
  if ( skizWithChamferDistance( resImage->buf, resImage->type, theDist, USHORT,
                                theDim, &theMask, 0 ) != 1 ) {
    vtfree( theDist );
    freeChamferMask( &theMask );
    if ( _verbose_ )
      fprintf( stderr, "%s: error in skiz computation\n", proc );
    return( -1 );
  }


  vtfree( theDist );
  freeChamferMask( &theMask );

  return( 1 );
}










static int _extremaImage( vt_image *intImage,
                            vt_image *mskImage,
                            vt_image *segImage,
                            vt_image *resImage,
                            int *res,
                            float percentile,
                            enumExtremaMethod method )
{
  char *proc = "_extremaImage";
  size_t i, size;
  typeHistogram histo, pdf, cumul;
  int theDim[3];

  size = intImage->dim.x * intImage->dim.y * intImage->dim.z;

  switch( method ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such method not handled yet\n", proc );
    return( -1 );
  case _GLOBAL_ :
      if ( _verbose_ )
        fprintf( stderr, "%s: use '_GLOBAL_' extremum computation method with %f %%\n", proc, percentile );
      break;
  case _CELL_ :
      if ( _verbose_ )
        fprintf( stderr, "%s: use '_CELL_' extremum computation method with %f %%\n", proc, percentile );
      break;
  case _CELL_BORDER_ :
      if ( _verbose_ )
        fprintf( stderr, "%s: use '_CELL_BORDER_' extremum computation method with %f %%\n", proc, percentile );
      break;
  case _CELL_INTERIOR_ :
      if ( _verbose_ )
        fprintf( stderr, "%s: use '_CELL_INTERIOR_' extremum computation method with %f %%\n", proc, percentile );
      break;
  case _VOXEL_ :
      if ( _verbose_ )
        fprintf( stderr, "%s: use '_VOXEL_' extremum computation method with %f %%\n", proc, percentile );
      break;
  }


  /* VOXEL method
   */
  switch( method ) {
  default :
    break;

  case _VOXEL_ :
    if ( _computeVoxelBasedExtremum( intImage, mskImage, resImage ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error during computation\n", proc );
      return( -1 );
    }
    break;
  }



  /* other methods
   */
  initHistogram( &histo );
  initHistogram( &pdf );
  initHistogram( &cumul );

  theDim[0] = intImage->dim.x * intImage->dim.v;
  theDim[1] = intImage->dim.y;
  theDim[2] = intImage->dim.z;



  /* histogram allocation
   */
  if ( mskImage != (vt_image*)NULL && mskImage->buf != (void*)NULL ) {
      if ( alloc1DHistogramFromImage( &histo, intImage->buf, intImage->type,
                                       mskImage->buf, mskImage->type, (double*)NULL, theDim ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: error when allocating/computing histogram\n", proc );
        return( -1 );
      }
  }
  else {
    if ( alloc1DHistogramFromImage( &histo, intImage->buf, intImage->type,
                                     (void*)NULL, TYPE_UNKNOWN, (double*)NULL, theDim ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when allocating/computing histogram\n", proc );
      return( -1 );
    }
  }

  if ( allocHistogramFromHistogram( &pdf, &histo, FLOAT ) != 1 ) {
    freeHistogram( &histo );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when allocating pdf histogram\n", proc );
    return( -1 );
  }

  if ( allocHistogramFromHistogram( &cumul, &pdf, FLOAT ) != 1 ) {
    freeHistogram( &pdf );
    freeHistogram( &histo );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when allocating cumulative histogram\n", proc );
    return( -1 );
  }



  /* GLOBAL method
   */
  switch( method ) {
  default :
    break;

  case _GLOBAL_ :
    if ( pdfHistogram( &pdf, &histo ) != 1 ) {
      freeHistogram( &cumul );
      freeHistogram( &pdf );
      freeHistogram( &histo );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when computing pdf histogram\n", proc );
      return( -1 );
    }
    if ( cumulative1DHistogram( &cumul, &pdf) != 1 ) {
        freeHistogram( &cumul );
        freeHistogram( &pdf );
        freeHistogram( &histo );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when computing cumulative histogram\n", proc );
        return( -1 );
    }

    *res = _getPercentileFromCumulative( &cumul, percentile );
    if ( *res == -33333 ) {
        freeHistogram( &cumul );
        freeHistogram( &pdf );
        freeHistogram( &histo );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when computing value from cumulative histogram\n", proc );
        return( -1 );
    }

    freeHistogram( &cumul );
    freeHistogram( &pdf );
    freeHistogram( &histo );

    if ( resImage != (vt_image*)NULL ) {
      switch ( resImage->type ) {
      default :
        if ( _verbose_ )
          fprintf( stderr, "%s: such image type not handled yet\n", proc );
        return( -1 );
      case USHORT :
        {
          u16 *resBuf = (u16*)resImage->buf;
          for ( i=0; i<size; i++, resBuf++ )
            *resBuf = *res;
        }
        break;
      }
    }

    return( 1 );
  }




  switch( method ) {
  default :
    break;

  case _CELL_ :
  case _CELL_BORDER_ :
  case _CELL_INTERIOR_ :
    if ( segImage == (vt_image*)NULL || segImage->buf == (void*)NULL ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: NULL segmentation image\n", proc );
        return( -1 );
    }
    if ( _computeCellBasedExtremum( intImage, segImage, resImage,
                                    &histo, &pdf, &cumul,
                                 percentile, method ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error during computation\n", proc );
      return( -1 );
    }
    break;
  }

  freeHistogram( &cumul );
  freeHistogram( &pdf );
  freeHistogram( &histo );

  return( 1 );
}











static int _rescaleIntImage( vt_image *theIm, vt_image *resIm,
                             vt_image *minIm, vt_image *maxIm,
                             int defaultMin, int defaultMax )
{
  char *proc = "_rescaleIntImage";
  size_t i, size;
  float a, b;
  int min, max, iy;

  size = theIm->dim.x * theIm->dim.y * theIm->dim.z;

  switch( resIm->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such output image ('%s') type not handled yet\n", proc, resIm->name );
    return( -1 );
  case UCHAR :
    {
      u8 *resBuf = (u8*)resIm->buf;
      int resMin = 0;
      int resMax = 255;
      switch( theIm->type ) {
      default :
        if ( _verbose_ )
          fprintf( stderr, "%s: such output image type not handled yet\n", proc );
        return( -1 );
      case USHORT :
        {
          u16 *theBuf = (u16*)theIm->buf;
          u16 *minBuf = (u16*)NULL;
          u16 *maxBuf = (u16*)NULL;
          if ( minIm != (vt_image*)NULL ) minBuf = (u16*)minIm->buf;
          if ( maxIm != (vt_image*)NULL ) maxBuf = (u16*)maxIm->buf;
          for ( i=0; i<size; i++ ) {
            min = ( minIm != (vt_image*)NULL ) ? (minBuf[i]) : defaultMin;
            max = ( maxIm != (vt_image*)NULL ) ? (maxBuf[i]) : defaultMax;
            a = (float)(resMax - resMin) / (float)(max-min);
            b = (float)resMin - a * (float)min;
            iy = (int)( a * (float)(theBuf[i]) + b + 0.5 );
            if ( iy <= resMin ) {
              resBuf[i] = resMin;
            }
            else if ( resMax <= iy ) {
              resBuf[i] = resMax;
            }
            else {
              resBuf[i] = iy;
            }
          }
        }
        break;
      }
    }
    /* end of case resIm->type == UCHAR
     */
    break;
  }

  return( 1 );
}





int AdhocFuse_RescaleIntensityImage( vt_image *intImage,
                                     vt_image *minMskImage,
                                     vt_image *maxMskImage,
                                     vt_image *segImage,
                                     vt_image *resIntImage,
                                     vt_image *resMinImage,
                                     vt_image *resMaxImage,
                                     int *resMin,
                                     int *resMax,
                                     float percentileMin,
                                     float percentileMax,
                                     enumExtremaMethod methodMin,
                                     enumExtremaMethod methodMax,
                                     float sigma )
{
  char *proc = "AdhocFuse_RescaleIntensityImage";

  vt_image tmpMinImage;
  vt_image *ptrMinImage = (vt_image*)NULL;

  vt_image tmpMaxImage;
  vt_image *ptrMaxImage = (vt_image*)NULL;

  typeFilteringCoefficients filter[3];
  int theDims[3];
  int borderLengths[3] = {0, 0, 0};

  int backgroundValue = 1;
  size_t i, v;



  VT_Image( &tmpMinImage );
  switch( methodMin ) {
  default :
  case _GLOBAL_ :
    ptrMinImage = resMinImage;
    break;
  case _CELL_ :
  case _CELL_BORDER_ :
  case _CELL_INTERIOR_ :
  case _VOXEL_ :
    if ( resMinImage == (vt_image*)NULL || resMinImage->buf == (void*)NULL ) {
      VT_InitFromImage( &tmpMinImage, intImage, (char*)NULL, intImage->type );
      if ( VT_AllocImage( &tmpMinImage ) != 1 ) {
        if ( _verbose_ ) {
          fprintf( stderr, "%s: minimum intensity image allocation failed\n", proc );
        }
        return( -1 );
      }
      ptrMinImage = &tmpMinImage;
    }
    else {
      ptrMinImage = resMinImage;
    }
  }

  VT_Image( &tmpMaxImage );
  switch( methodMax ) {
  default :
  case _GLOBAL_ :
    ptrMaxImage = resMaxImage;
    break;
  case _CELL_ :
  case _CELL_BORDER_ :
  case _CELL_INTERIOR_ :
  case _VOXEL_ :
    if ( resMaxImage == (vt_image*)NULL || resMaxImage->buf == (void*)NULL ) {
      VT_InitFromImage( &tmpMaxImage, intImage, (char*)NULL, intImage->type );
      if ( VT_AllocImage( &tmpMaxImage ) != 1 ) {
        if ( ptrMinImage == &tmpMinImage ) VT_FreeImage( &tmpMinImage );
        if ( _verbose_ ) {
          fprintf( stderr, "%s: maximum intensity image allocation failed\n", proc );
        }
        return( -1 );
      }
      ptrMaxImage = &tmpMaxImage;
    }
    else {
      ptrMaxImage = resMaxImage;
    }
  }


  /* the cell label image may contain 0
   * due to resampling issues
   */
  if ( segImage != (vt_image*)NULL ) {
    v = segImage->dim.x * segImage->dim.y * segImage->dim.z;

    switch( segImage->type ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such segmentation image type not handled yet\n", proc );
      return( -1 );
    case USHORT :
      {
        u16 *segBuf = (u16*)segImage->buf;
        for ( i=0; i<v; i++ ) {
          if ( segBuf[i] == 0 ) segBuf[i] = backgroundValue;
        }
      }
      break;
    }
  }


  if  ( _extremaImage( intImage, minMskImage, segImage, ptrMinImage,
                       resMin, percentileMin, methodMin ) != 1 ) {
    if ( ptrMaxImage == &tmpMaxImage ) VT_FreeImage( &tmpMaxImage );
    if ( ptrMinImage == &tmpMinImage ) VT_FreeImage( &tmpMinImage );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when computing minimum image/value\n", proc );
    return( -1 );
  }

  if  ( _extremaImage( intImage, maxMskImage, segImage, ptrMaxImage,
                       resMax, percentileMax, methodMax ) != 1 ) {
    if ( ptrMaxImage == &tmpMaxImage ) VT_FreeImage( &tmpMaxImage );
    if ( ptrMinImage == &tmpMinImage ) VT_FreeImage( &tmpMinImage );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when computing maximum image/value\n", proc );
    return( -1 );
  }

  if ( sigma > 0 ) {

    initFilteringCoefficients( &(filter[0]) );
    initFilteringCoefficients( &(filter[1]) );
    initFilteringCoefficients( &(filter[2]) );

    filter[0].derivative = DERIVATIVE_0;
    filter[1].derivative = DERIVATIVE_0;
    filter[2].derivative = ( intImage->dim.z == 1 ) ? NODERIVATIVE : DERIVATIVE_0;
    filter[0].coefficient = sigma;
    filter[1].coefficient = sigma;
    filter[2].coefficient = sigma;

    theDims[0] = intImage->dim.x;
    theDims[1] = intImage->dim.y;
    theDims[2] = intImage->dim.z;

    if ( ptrMinImage != (vt_image*)NULL ) {
      if ( separableLinearFiltering( ptrMinImage->buf, ptrMinImage->type,
                                     ptrMinImage->buf, ptrMinImage->type, theDims,
                                     borderLengths, filter ) != 1 ) {
        if ( ptrMaxImage == &tmpMaxImage ) VT_FreeImage( &tmpMaxImage );
        if ( ptrMinImage == &tmpMinImage ) VT_FreeImage( &tmpMinImage );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when smoothing minimum image\n", proc );
        return( -1 );
      }
    }

    if ( ptrMaxImage != (vt_image*)NULL ) {
      if ( separableLinearFiltering( ptrMaxImage->buf, ptrMaxImage->type,
                                     ptrMaxImage->buf, ptrMaxImage->type, theDims,
                                     borderLengths, filter ) != 1 ) {
        if ( ptrMaxImage == &tmpMaxImage ) VT_FreeImage( &tmpMaxImage );
        if ( ptrMinImage == &tmpMinImage ) VT_FreeImage( &tmpMinImage );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when smoothing maximum image\n", proc );
        return( -1 );
      }
    }

  }

  if ( _rescaleIntImage( intImage, resIntImage, ptrMinImage, ptrMaxImage,
                         *resMin, *resMax ) != 1 ) {
    if ( ptrMaxImage == &tmpMaxImage ) VT_FreeImage( &tmpMaxImage );
    if ( ptrMinImage == &tmpMinImage ) VT_FreeImage( &tmpMinImage );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when rescaling intensity image\n", proc );
    return( -1 );
  }

  if ( ptrMaxImage == &tmpMaxImage ) VT_FreeImage( &tmpMaxImage );
  if ( ptrMinImage == &tmpMinImage ) VT_FreeImage( &tmpMinImage );

  return( 1 );
}





















