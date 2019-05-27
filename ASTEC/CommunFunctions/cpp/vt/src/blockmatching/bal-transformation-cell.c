/*************************************************************************
 * bal-transformation-cell.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Jeu 29 nov 2018 15:38:39 CET
 *
 * ADDITIONS, CHANGES
 *
 */





#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vtmalloc.h>

#include <bal-transformation-cell.h>
#include <bal-transformation-inversion.h>
#include <bal-transformation-tools.h>


static int _verbose_ = 1;

void BAL_SetVerboseInBalTransformationCell( int v )
{
  _verbose_ = v;
}

void BAL_IncrementVerboseInBalTransformationCell(  )
{
  _verbose_ ++;
}

void BAL_DecrementVerboseInBalTransformationCell(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}





/*******************************************************************
 *
 *
 *
 *******************************************************************/


typedef struct typeCell {

    /* bounding box
     */
    int ptmin[3];
    int ptmax[3];

    /* #points
     */
    int npts;

} typeCell;



static void _initCell( typeCell *c )
{
  c->ptmin[0] = c->ptmin[1] = c->ptmin[2] = 0;
  c->ptmax[0] = c->ptmax[1] = c->ptmax[2] = 0;

  c->npts = 0;
}


typedef struct typeCellList {
  typeCell *data;
  int n_data;
} typeCellList;



static void _initCellList( typeCellList *l )
{
  l->data = (typeCell*)NULL;
  l->n_data = 0;
}



static void _freeCellList( typeCellList *l )
{
  if ( l->data != (typeCell*)NULL ) vtfree( l->data );
  _initCellList( l );
}



static int _allocCellList( typeCellList *l , int n )
{
  char *proc = "_allocCellList";
  int i;

  l->data = (typeCell*)vtmalloc( n * sizeof(typeCell), "l->data", proc );

  if ( l->data == (typeCell*)NULL ) {
      if ( _verbose_ )
          fprintf( stderr, "%s: unable to allocate data\n", proc );
    return( -1 );
  }
  l->n_data = n;

  for ( i=0; i<l->n_data; i++ )
      _initCell( &(l->data[i]) );

  return( 1 );
}



int _extractCellFeatures( bal_image *theIm, typeCellList *theList )
{
  char *proc = "_extractCellFeatures";
  int ncells = 0;
  size_t i, v, x, y, z;
  typeCell *cell;


  /* extract information
   * about cells
   */

  v = theIm->ncols * theIm->nrows * theIm->nplanes;

#define _NUMBER_OF_CELLS( TYPE ) {                \
  TYPE *theBuf = (TYPE*)theIm->data;              \
  for ( ncells=0, i=0; i<v; i++ ) {               \
    if ( ncells < theBuf[i] ) ncells = theBuf[i]; \
  }                                               \
}

  switch ( theIm->type ) {
  default :
    fprintf( stderr, "%s: type not handled in switch\n", proc );
    return( -1 );
  case UCHAR :
    _NUMBER_OF_CELLS( u8 );
    break;
  case SSHORT :
    _NUMBER_OF_CELLS( s16 );
    break;
  case USHORT :
    _NUMBER_OF_CELLS( u16 );
    break;
  }

  if ( _verbose_ >= 2 )
    fprintf( stderr, "%s: largest cell label is %d\n", proc, ncells );

  if ( ncells <= 0 ) return( 0 );

  /* allocate list of cells
   */
  if ( _allocCellList( theList, ncells+1 ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }

#define _CELLPROPERTIES( TYPE ) {              \
  TYPE ***theBuf = (TYPE***)theIm->array;      \
  for ( z=0; z<theIm->nplanes; z++ ) {         \
    for ( y=0; y<theIm->nrows; y++ )           \
    for ( x=0; x<theIm->ncols; x++ ) {         \
      cell = &(theList->data[ theBuf[z][y][x] ]); \
      if ( cell->npts == 0 ) {                 \
        cell->ptmin[0] = cell->ptmax[0] = x;   \
        cell->ptmin[1] = cell->ptmax[1] = y;   \
        cell->ptmin[2] = cell->ptmax[2] = z;   \
      }                                        \
      else {                                   \
        if ( cell->ptmin[0] > (int)x ) cell->ptmin[0] = x; \
        if ( cell->ptmin[1] > (int)y ) cell->ptmin[1] = y; \
        if ( cell->ptmin[2] > (int)z ) cell->ptmin[2] = z; \
        if ( cell->ptmax[0] < (int)x ) cell->ptmax[0] = x; \
        if ( cell->ptmax[1] < (int)y ) cell->ptmax[1] = y; \
        if ( cell->ptmax[2] < (int)z ) cell->ptmax[2] = z; \
      }                                        \
      cell->npts ++;                           \
    }                                          \
  }                                            \
}

  switch ( theIm->type ) {
  default :
    fprintf( stderr, "%s: type not handled in switch\n", proc );
    return( -1 );
  case UCHAR :
    _CELLPROPERTIES( u8 );
    break;
  case USHORT :
    _CELLPROPERTIES( u16 );
    break;
  }

  return( 1 );
}







/*******************************************************************
 *
 *
 *
 *******************************************************************/

static int _updateCellImageGeometry( bal_image *cellIm, bal_image *theIm, int *ptmin, int margin )
{
  char *proc = "_updateCellImageGeometry";
  bal_transformation theTrsf;
  bal_transformation invTrsf;
  _MATRIX sform_to_voxel;

  /* set the image geometry
   */

  BAL_InitTransformation( &theTrsf );

  if ( BAL_AllocTransformation( &theTrsf, AFFINE_3D, (bal_image*)NULL ) != 1 ) {
    if ( _verbose_)
      fprintf( stderr, "%s: unable to allocate transformation \n", proc);
    return( -1 );
  }
  BAL_SetTransformationToIdentity( &theTrsf );
  theTrsf.transformation_unit = VOXEL_UNIT;

  theTrsf.mat.m[ 3] = ptmin[0]-margin;
  theTrsf.mat.m[ 7] = ptmin[1]-margin;
  theTrsf.mat.m[11] = ptmin[2]-margin;

  BAL_InitTransformation( &invTrsf );
  if ( BAL_AllocTransformation( &invTrsf, AFFINE_3D, (bal_image*)NULL ) != 1 ) {
    BAL_FreeTransformation( &theTrsf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate inverse transformation \n", proc );
    return( -1 );
  }

  if ( BAL_InverseTransformation( &theTrsf, &invTrsf ) != 1 ) {
    BAL_FreeTransformation( &invTrsf );
    BAL_FreeTransformation( &theTrsf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to invert transformation \n", proc );
    return( -1 );
  }

  /* re-computation of Qform matrices
   */

  _mult_mat( &(invTrsf.mat), &(theIm->to_voxel), &(cellIm->to_voxel) );

  if ( InverseMat4x4( cellIm->to_voxel.m, cellIm->to_real.m ) != 4 ) {
    BAL_FreeTransformation( &invTrsf );
    BAL_FreeTransformation( &theTrsf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to invert qform matrix\n", proc );
    return( -1 );
  }

  cellIm->geometry = _BAL_QFORM_GEOMETRY_;
  cellIm->qform_code = theIm->qform_code;
  if ( cellIm->qform_code == 0 )
    cellIm->qform_code = 1;

  /* re-computation of Sform matrices
   */

  _init_mat( &(sform_to_voxel) );
  if ( _alloc_mat( &(sform_to_voxel), 4, 4 ) != 1 ) {
    BAL_FreeTransformation( &invTrsf );
    BAL_FreeTransformation( &theTrsf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate sform 'to_voxel' matrix\n", proc );
    return( -1 );
  }

  if ( InverseMat4x4( theIm->sform_to_real.m, sform_to_voxel.m ) != 4 ) {
    _free_mat( &(sform_to_voxel) );
    BAL_FreeTransformation( &invTrsf );
    BAL_FreeTransformation( &theTrsf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to invert sform real matrix\n", proc );
    return( -1 );
  }

  _mult_mat( &(invTrsf.mat), &sform_to_voxel, &sform_to_voxel );

  if ( InverseMat4x4( sform_to_voxel.m, cellIm->sform_to_real.m ) != 4 ) {
    _free_mat( &(sform_to_voxel) );
    BAL_FreeTransformation( &invTrsf );
    BAL_FreeTransformation( &theTrsf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to invert sform voxel matrix\n", proc );
    return( -1 );
  }

  cellIm->sform_code = theIm->sform_code;
  if ( cellIm->sform_code == 0 )
    cellIm->sform_code = 1;

  _free_mat( &(sform_to_voxel) );
  BAL_FreeTransformation( &invTrsf );
  BAL_FreeTransformation( &theTrsf );

  return( 1 );
}



static int _extractUniqueCellImage( bal_image *weightIm, bal_image *theIm,
                              typeCellList *theList, int c, int margin )
{
  char *proc = "_extractUniqueCellImage";
  typeCell *cell;
  int dimx, dimy, dimz;
  int i, j, k, x, y, z;
  r32 ***weightBuf;

  cell = &(theList->data[c]);
  dimx = cell->ptmax[0] - cell->ptmin[0] + 1 + 2 * margin;
  dimy = cell->ptmax[1] - cell->ptmin[1] + 1 + 2 * margin;
  dimz = cell->ptmax[2] - cell->ptmin[2] + 1 + 2 * margin;

  if ( BAL_AllocFullImage( weightIm, (char*)NULL, dimx, dimy, dimz, theIm->vdim,
                           theIm->vx, theIm->vy, theIm->vz, FLOAT ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation of weight sub-image for cell %d failed\n", proc, c );
    return( -1 );
  }

  weightBuf = (r32***)weightIm->array;

#define _UNIQUECELLEXTRACTION( TYPE ) {        \
  TYPE ***theBuf = (TYPE***)theIm->array;      \
  for ( z=cell->ptmin[2], k=margin; z<=cell->ptmax[2]; z++, k++ ) {   \
    for ( y=cell->ptmin[1], j=margin; y<=cell->ptmax[1]; y++, j++ )   \
    for ( x=cell->ptmin[0], i=margin; x<=cell->ptmax[0]; x++, i++ ) { \
      if ( theBuf[z][y][x] == c )              \
        weightBuf[k][j][i] = 1.0;              \
    }                                          \
  }                                            \
}

  switch ( theIm->type ) {
  default :
    BAL_FreeImage( weightIm );
    fprintf( stderr, "%s: type not handled in switch\n", proc );
    return( -1 );
  case UCHAR :
    _UNIQUECELLEXTRACTION( u8 );
    break;
  case USHORT :
    _UNIQUECELLEXTRACTION( u16 );
    break;
  }

  /* set the image geometry
   */
  if ( _updateCellImageGeometry( weightIm, theIm, cell->ptmin, margin ) != 1 ) {
    BAL_FreeImage( weightIm );
    if ( _verbose_)
      fprintf( stderr, "%s: unable to update weight image geometry \n", proc);
    return( -1 );
  }

  return( 1 );
}










static int _updateUniqueCellImage( bal_image *trsfWeightIm, bal_image *resCellIm, bal_image *resWeightIm,
                             int c, typeCell *theBoundingBox )
{
  char *proc = "_updateUniqueCellImage";
  size_t i, j, k;
  r32 ***trsfWeightBuf;
  r32 ***resWeightBuf;
  size_t i0, j0, k0;

  trsfWeightBuf = (r32***)trsfWeightIm->array;
  resWeightBuf = (r32***)resWeightIm->array;

  i0 = theBoundingBox->ptmin[0];
  j0 = theBoundingBox->ptmin[1];
  k0 = theBoundingBox->ptmin[2];

  if ( k0+trsfWeightIm->nplanes > resCellIm->nplanes || j0+trsfWeightIm->nrows > resCellIm->nrows || i0+trsfWeightIm->ncols > resCellIm->ncols ) {
      fprintf( stderr, "%s: transformed cell image and result image have incompatible dimensions\n", proc );
      return( -1 );
  }
  if ( resWeightIm->nplanes != resCellIm->nplanes || resWeightIm->nrows != resCellIm->nrows || resWeightIm->ncols != resCellIm->ncols ) {
      fprintf( stderr, "%s: maximum weight image and result image have different dimensions\n", proc );
      return( -1 );
  }



#define _UNIQUECELLUPDATE( TYPE ) {                       \
  TYPE ***resCellBuf = (TYPE***)resCellIm->array;   \
  for ( k=0; k<trsfWeightIm->nplanes; k++ )         \
  for ( j=0; j<trsfWeightIm->nrows; j++ )           \
  for ( i=0; i<trsfWeightIm->ncols; i++ ) {         \
    if ( trsfWeightBuf[k][j][i] <= resWeightBuf[k0+k][j0+j][i0+i] || trsfWeightBuf[k][j][i] <= 0.000001 ) \
      continue;                                     \
    resWeightBuf[k0+k][j0+j][i0+i] = trsfWeightBuf[k][j][i]; \
    resCellBuf[k0+k][j0+j][i0+i] = c;                        \
  }                                                 \
}

  switch ( resCellIm->type ) {
  default :
    fprintf( stderr, "%s: type not handled in switch\n", proc );
    return( -1 );
  case UCHAR :
    _UNIQUECELLUPDATE( u8 );
    break;
  case USHORT :
    _UNIQUECELLUPDATE( u16 );
    break;
  }

  return( 1 );
}








static int _allocTrsfImage( bal_image *trsfWeightImage,
                             typeCell *trsfBoundingBox,
                             bal_transformation *theTrsf, bal_image *resImage,
                             bal_image *weightImage )
{
  char *proc = "_allocTrsfImages";
  bal_transformation invTrsf;
  bal_floatPoint tmpPt, minPt, maxPt, thePt[8];
  double *mat;
  int i;
  int dimx, dimy, dimz;

  /* inverse transformation computation
   */
  BAL_InitTransformation( &invTrsf );
  if ( BAL_AllocTransformation( &invTrsf, AFFINE_3D, (bal_image*)NULL ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate inverse transformation \n", proc );
    return( -1 );
  }

  switch ( theTrsf->type ) {
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
    return( -1 );
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
    if ( BAL_AllocTransformation( &invTrsf, theTrsf->type, (bal_image *)NULL ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate inverse transformation (linear case)\n", proc );
      return( -1 );
    }
    break;
  }

  if ( BAL_InverseTransformation( theTrsf, &invTrsf ) != 1 ) {
    BAL_FreeTransformation( &invTrsf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to invert transformation \n", proc );
    return( -1 );
  }

  /* bounding box definition
   */
  thePt[0].x = 0;                      thePt[0].y = 0;                      thePt[0].z = 0;
  thePt[1].x = weightImage->ncols - 1; thePt[1].y = 0;                      thePt[1].z = 0;
  thePt[2].x = 0;                      thePt[2].y = weightImage->nrows - 1; thePt[2].z = 0;
  thePt[3].x = weightImage->ncols - 1; thePt[3].y = weightImage->nrows - 1; thePt[3].z = 0;
  thePt[4].x = 0;                      thePt[4].y = 0;                      thePt[4].z = weightImage->nplanes - 1;
  thePt[5].x = weightImage->ncols - 1; thePt[5].y = 0;                      thePt[5].z = weightImage->nplanes - 1;
  thePt[6].x = 0;                      thePt[6].y = weightImage->nrows - 1; thePt[6].z = weightImage->nplanes - 1;
  thePt[7].x = weightImage->ncols - 1; thePt[7].y = weightImage->nrows - 1; thePt[7].z = weightImage->nplanes - 1;

  /* bounding box transformation
   */
  for ( i=0; i<8; i++ ) {
    /* extracted cell image geometry: voxel to real
     */
    mat = weightImage->to_real.m;
    tmpPt.x = mat[ 0] * thePt[i].x + mat[ 1] * thePt[i].y + mat[ 2] * thePt[i].z + mat[ 3];
    tmpPt.y = mat[ 4] * thePt[i].x + mat[ 5] * thePt[i].y + mat[ 6] * thePt[i].z + mat[ 7];
    tmpPt.z = mat[ 8] * thePt[i].x + mat[ 9] * thePt[i].y + mat[10] * thePt[i].z + mat[11];
    thePt[i].x = tmpPt.x;
    thePt[i].y = tmpPt.y;
    thePt[i].z = tmpPt.z;
    /* from extracted cell image geometry
     * to result image geometry
     */
    if ( BAL_TransformFloatPoint( &(thePt[i]), &(thePt[i]), &invTrsf ) != 1 ) {
      BAL_FreeTransformation( &invTrsf );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute point transformation for corner #%d\n", proc, i );
      return( -1 );
    }
    /* result image geometry: real to voxel
     */
    mat = resImage->to_voxel.m;
    tmpPt.x = mat[ 0] * thePt[i].x + mat[ 1] * thePt[i].y + mat[ 2] * thePt[i].z + mat[ 3];
    tmpPt.y = mat[ 4] * thePt[i].x + mat[ 5] * thePt[i].y + mat[ 6] * thePt[i].z + mat[ 7];
    tmpPt.z = mat[ 8] * thePt[i].x + mat[ 9] * thePt[i].y + mat[10] * thePt[i].z + mat[11];
    thePt[i].x = tmpPt.x;
    thePt[i].y = tmpPt.y;
    thePt[i].z = tmpPt.z;
  }

  BAL_FreeTransformation( &invTrsf );

  /* resulting bounding box
   */
  minPt.x = maxPt.x = thePt[0].x;
  minPt.y = maxPt.y = thePt[0].y;
  minPt.z = maxPt.z = thePt[0].z;
  for ( i=1; i<8; i++ ) {
    if ( minPt.x > thePt[i].x ) minPt.x = thePt[i].x;
    if ( minPt.y > thePt[i].y ) minPt.y = thePt[i].y;
    if ( minPt.z > thePt[i].z ) minPt.z = thePt[i].z;
    if ( maxPt.x < thePt[i].x ) maxPt.x = thePt[i].x;
    if ( maxPt.y < thePt[i].y ) maxPt.y = thePt[i].y;
    if ( maxPt.z < thePt[i].z ) maxPt.z = thePt[i].z;
  }

  _initCell( trsfBoundingBox );

  trsfBoundingBox->ptmin[0] = 0; trsfBoundingBox->ptmax[0] = resImage->ncols-1;
  trsfBoundingBox->ptmin[1] = 0; trsfBoundingBox->ptmax[1] = resImage->nrows-1;
  trsfBoundingBox->ptmin[2] = 0; trsfBoundingBox->ptmax[2] = resImage->nplanes-1;

  if ( 0 <= minPt.x && minPt.x <= resImage->ncols-2 ) {
    trsfBoundingBox->ptmin[0] = (int)(minPt.x + 0.5);
  }
  if ( 0 <= minPt.y && minPt.y <= resImage->nrows-2 ) {
    trsfBoundingBox->ptmin[1] = (int)(minPt.y + 0.5);
  }
  if ( 0 <= minPt.z && minPt.z <= resImage->nplanes-2 ) {
    trsfBoundingBox->ptmin[2] = (int)(minPt.z + 0.5);
  }

  if ( 1 <= maxPt.x && maxPt.x <= resImage->ncols-1 ) {
    trsfBoundingBox->ptmax[0] = (int)(maxPt.x + 0.5);
  }
  if ( 1 <= maxPt.y && maxPt.y <= resImage->nrows-1 ) {
    trsfBoundingBox->ptmax[1] = (int)(maxPt.y + 0.5);
  }
  if ( 1 <= maxPt.z && maxPt.z <= resImage->nplanes-1 ) {
    trsfBoundingBox->ptmax[2] = (int)(maxPt.z + 0.5);
  }

  if ( 0 ) {
    fprintf( stderr, "   ... transformed bounding box = [%f %f %f] x [%f %f %f]\n",
             minPt.x, minPt.y, minPt.z, maxPt.x, maxPt.y, maxPt.z );
    fprintf( stderr, "   ... new voxel bounding box = [%d %d %d] x [%d %d %d]\n",
             trsfBoundingBox->ptmin[0], trsfBoundingBox->ptmin[1], trsfBoundingBox->ptmin[2],
             trsfBoundingBox->ptmax[0], trsfBoundingBox->ptmax[1], trsfBoundingBox->ptmax[2] );
  }

  /* image allocation
   */
  dimx = trsfBoundingBox->ptmax[0] - trsfBoundingBox->ptmin[0] + 1;
  dimy = trsfBoundingBox->ptmax[1] - trsfBoundingBox->ptmin[1] + 1;
  dimz = trsfBoundingBox->ptmax[2] - trsfBoundingBox->ptmin[2] + 1;

  if ( BAL_AllocFullImage( trsfWeightImage, (char*)NULL, dimx, dimy, dimz, resImage->vdim,
                           resImage->vx, resImage->vy, resImage->vz, weightImage->type ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation of result weight sub-image failed\n", proc );
    return( -1 );
  }
  /* set the image geometry
   */
  if ( _updateCellImageGeometry( trsfWeightImage, resImage, trsfBoundingBox->ptmin, 0 ) != 1 ) {
    BAL_FreeImage( trsfWeightImage );
    if ( _verbose_)
      fprintf( stderr, "%s: unable to update result weight image geometry \n", proc);
    return( -1 );
  }

  return( 1 );
}




/*******************************************************************
 *
 * end of static resampling procedures
 *
 *******************************************************************/





int BAL_ResampleCellImage( bal_image *image, bal_image *resImage,
                           bal_image *weightImage,
                           bal_transformation *theTr,
                           float sigma )
{
  char *proc = "BAL_ResampleCellImage";
  typeCellList theList;

  bal_transformation *ptrTr;
  bal_transformation theId;

  bal_image tmpResImage;
  bal_image *ptrResImage = (bal_image*)NULL;

  bal_image *ptrWeightImage = (bal_image*)NULL;
  bal_image globalWeightImage;

  bal_image tmpWeightImage;
  typeCell tmpBoundingBox;

  bal_image extractedWeightImage;

  int cell;
  int margin = 10;

  bal_doublePoint theSigma;



  /* require a positive sigma
   */
  if ( sigma <= 0.0 ) {
      if ( _verbose_ )
          fprintf( stderr, "%s: switch to nearest neighbor resampling\n", proc );
      return( BAL_ResampleImage( image, resImage, theTr, NEAREST ) );
  }

  /* identity transformation, if none
   * could be image to image transformation
   */
  BAL_InitTransformation( &theId );
  if ( theTr == (bal_transformation*)NULL ) {
      if ( BAL_AllocTransformation( &theId, RIGID_3D, (bal_image*)NULL ) != 1 ) {
          if ( _verbose_ )
              fprintf( stderr, "%s: auxiliary transformation allocation failed\n", proc );
          return( -1 );
      }
      if ( BAL_ComputeImageToImageTransformation( resImage, image, &theId ) != 1 ) {
          if ( _verbose_ )
              fprintf( stderr, "%s: unable to compute image to image transformation\n", proc );
          return( -1 );
      }
      ptrTr = &theId;
  }
  else {
      ptrTr = theTr;
  }

  /* auxiliary result image if required
   */
  BAL_InitImage( &tmpResImage, (char*)NULL, 0, 0, 0, 0, TYPE_UNKNOWN );
  if ( resImage == image ) {
    if ( BAL_AllocImageFromImage( &tmpResImage, (char*)NULL, image, image->type ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: auxiliary result image allocation failed\n", proc );
      return( -1 );
    }
    ptrResImage = &tmpResImage;
  }
  else {
      ptrResImage = resImage;
  }

  /* cell properties: bounding boxes, ...
   */
  _initCellList( &theList );
  if ( _extractCellFeatures( image, &theList ) != 1 ) {
      if ( resImage == image ) BAL_FreeImage( &tmpResImage );
      if ( _verbose_ ) {
          fprintf( stderr, "%s: image filled by 0\n", proc );
      }
      return( -1 );
  }

  /* cumulative (in the maximum sense) weight image
   * keep memory of maximum of weights (result image geometry)
   */
  if ( weightImage != (bal_image*)NULL ) {
    ptrWeightImage = weightImage;
  }
  else {
    if ( BAL_AllocImageFromImage( &globalWeightImage, (char*)NULL, resImage, FLOAT ) != 1 ) {
      _freeCellList( &theList );
      if ( resImage == image ) BAL_FreeImage( &tmpResImage );
      if ( _verbose_ )
        fprintf( stderr, "%s: auxiliary weight image allocation failed\n", proc );
      return( -1 );
    }
    ptrWeightImage = &globalWeightImage;
  }

  /* last parameters
   */

  theSigma.x = sigma;
  theSigma.y = sigma;
  theSigma.z = sigma;


  /* loop over cell(s)
   */

  for ( cell=0; cell<theList.n_data; cell++ ) {

    if ( theList.data[cell].npts <= 0 ) continue;

    if ( _verbose_ >= 3 ) {
      fprintf( stderr, "... processing cell #%3d\n", cell );
    }

    /* extract weight and cell sub-image(s)
     */
    BAL_InitImage( &extractedWeightImage, (char*)NULL, 0, 0, 0, 0, TYPE_UNKNOWN );

    if ( _extractUniqueCellImage( &extractedWeightImage, image, &theList, cell, margin ) != 1 ) {
      if ( weightImage == (bal_image*)NULL ) BAL_FreeImage( &globalWeightImage );
      _freeCellList( &theList );
      if ( resImage == image ) BAL_FreeImage( &tmpResImage );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when extracting weight sub-image\n", proc );
      return( -1 );
    }

    if ( _allocTrsfImage( &tmpWeightImage, &tmpBoundingBox, ptrTr, resImage, &extractedWeightImage ) != 1 ) {
      if ( weightImage == (bal_image*)NULL ) BAL_FreeImage( &globalWeightImage );
      _freeCellList( &theList );
      if ( resImage == image ) BAL_FreeImage( &tmpResImage );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when allocating result sub-image\n", proc );
      return( -1 );
    }

    /* smooth weight sub-image
     */
    if ( sigma > 0.0 ) {
        if ( BAL_SmoothImage( &extractedWeightImage, &theSigma ) != 1 ) {
          BAL_FreeImage( &extractedWeightImage );
          BAL_FreeImage( &tmpWeightImage );
          if ( weightImage == (bal_image*)NULL ) BAL_FreeImage( &globalWeightImage );
          _freeCellList( &theList );
            if ( resImage == image ) BAL_FreeImage( &tmpResImage );
            if ( _verbose_ )
              fprintf( stderr, "%s: error when smoothing sub-image of cell %d\n", proc, cell );
            return( -1 );
        }
    }

    /* transform weight sub-image
     */
    if ( BAL_ResampleImage( &extractedWeightImage, &tmpWeightImage, ptrTr, LINEAR ) != 1 ) {
      BAL_FreeImage( &extractedWeightImage );
      BAL_FreeImage( &tmpWeightImage );
      if ( weightImage == (bal_image*)NULL ) BAL_FreeImage( &globalWeightImage );
      _freeCellList( &theList );
      if ( resImage == image ) BAL_FreeImage( &tmpResImage );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when resampling weight sub-image\n", proc );
      return( -1 );
    }

    BAL_FreeImage( &extractedWeightImage );

    if ( _updateUniqueCellImage( &tmpWeightImage, ptrResImage, ptrWeightImage, cell, &tmpBoundingBox ) != 1 ) {
      BAL_FreeImage( &tmpWeightImage );
      if ( weightImage == (bal_image*)NULL ) BAL_FreeImage( &globalWeightImage );
      _freeCellList( &theList );
        if ( resImage == image ) BAL_FreeImage( &tmpResImage );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when updating with cell %d\n", proc, cell );
        return( -1 );
    }

    BAL_FreeImage( &tmpWeightImage );

  }
  /* loop is over
   */

  if ( weightImage == (bal_image*)NULL ) BAL_FreeImage( &globalWeightImage );
  _freeCellList( &theList );

  if ( resImage == image ) {
    if ( BAL_CopyImage( &tmpResImage, resImage ) != 1 ) {
      BAL_FreeImage( &tmpResImage );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when copying result image\n", proc  );
      return( -1 );
    }
    BAL_FreeImage( &tmpResImage );
  }

  return( 1 );
}
