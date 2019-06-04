/*************************************************************************
 * bal-transformation-list-tools.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mar 21 jan 2014 18:00:31 CET
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <chunks.h>
#include <vtmalloc.h>

#include <bal-transformation-compose.h>
#include <bal-transformation-inversion.h>
#include <bal-transformation-tools.h>

#include <bal-transformation-list-tools.h>





static int _verbose_ = 1;
static int _debug_ = 0;
static int _warning_ = 0;





void BAL_SetVerboseInBalTransformationListTools( int v )
{
  _verbose_ = v;
}

void BAL_IncrementVerboseInBalTransformationListTools(  )
{
  _verbose_ ++;
}

void BAL_DecrementVerboseInBalTransformationListTools(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}

void BAL_SetDebugInBalTransformationListTools( int v )
{
  _debug_ = v;
}

void BAL_IncrementDebugInBalTransformationListTools(  )
{
  _debug_ ++;
}

void BAL_DecrementDebugInBalTransformationListTools(  )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}





/**********************************************************************
 *
 * box management
 *
 **********************************************************************/



static void _initBox( typeBox *b )
{
  b->min.x = 0.0;
  b->min.y = 0.0;
  b->min.z = 0.0;

  b->max.x = 0.0;
  b->max.y = 0.0;
  b->max.z = 0.0;

  b->fov.x = 0.0;
  b->fov.y = 0.0;
  b->fov.z = 0.0;
}



void _initBoxList( typeBoxList *l )
{
  l->data = (typeBox*)NULL;
  l->n_data = 0;
}



void _freeBoxList( typeBoxList *l )
{
  if ( l->data != (typeBox*)NULL )
    vtfree( l->data );
  _initBoxList( l );
}



int _allocBoxList( typeBoxList *l, int n )
{
  char *proc = "_allocBoxList";
  int i;

  l->data = (typeBox*)vtmalloc( n * sizeof(typeBox), "l->data", proc );
  if ( l->data == (typeBox*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }

  l->n_data = n;

  for ( i=0; i<n; i++ )
    _initBox( &(l->data[i]) );

  return( 1 );
}





/**********************************************************************
 *
 * box I/O
 *
 **********************************************************************/


static void _fprintfBox( FILE *f, typeBox *b, char *desc )
{
  if ( desc != (char*)NULL )
    fprintf( f, " * box '%s'\n", desc );
  fprintf( stderr, "   - MIN = %f , %f , %f\n", b->min.x, b->min.y, b->min.z );
  fprintf( stderr, "   - MAX = %f , %f , %f\n", b->max.x, b->max.y, b->max.z );
  fprintf( stderr, "   - FOV = %f x %f x %f\n", b->fov.x, b->fov.y, b->fov.z );
}



void _fprintfBoxList( FILE *f, typeBoxList *l, char *desc )
{
  int i;
  char str[16];

  if ( desc != (char*)NULL )
    fprintf( f, "box list '%s'\n", desc );

  for ( i=0; i<l->n_data; i++ ) {
    sprintf( str, "#%d", i );
    _fprintfBox( f, &(l->data[i]), str );
  }
}





/**********************************************************************
 *
 * box utilities
 *
 **********************************************************************/







static void _maxCorners( typeBoxList *theList,
                     bal_floatPoint *mincorner,
                     bal_floatPoint *maxcorner )
{
  char *proc = "_maxCorners";
  int i;
  bal_integerPoint minindex;
  bal_integerPoint maxindex;

  mincorner->x = theList->data[0].min.x;
  mincorner->y = theList->data[0].min.y;
  mincorner->z = theList->data[0].min.z;

  maxcorner->x = theList->data[0].max.x;
  maxcorner->y = theList->data[0].max.y;
  maxcorner->z = theList->data[0].max.z;

  minindex.x = maxindex.x = 0;
  minindex.y = maxindex.y = 0;
  minindex.z = maxindex.z = 0;

  for ( i=0; i<theList->n_data; i++ ) {

    if ( mincorner->x > theList->data[i].min.x ) {
      mincorner->x = theList->data[i].min.x;
      minindex.x = i;
    }
    if ( mincorner->y > theList->data[i].min.y ) {
      mincorner->y = theList->data[i].min.y;
      minindex.y = i;
    }
    if ( mincorner->z > theList->data[i].min.z ) {
      mincorner->z = theList->data[i].min.z;
      minindex.z = i;
    }

    if ( maxcorner->x < theList->data[i].max.x ) {
      maxcorner->x = theList->data[i].max.x;
      maxindex.x = i;
    }
    if ( maxcorner->y < theList->data[i].max.y ) {
      maxcorner->y = theList->data[i].max.y;
      maxindex.y = i;
    }
    if ( maxcorner->z < theList->data[i].max.z ) {
      maxcorner->z = theList->data[i].max.z;
      maxindex.z = i;
    }

  }


  if ( _debug_ >= 2 ) {
    fprintf( stderr, "%s: Computed corners\n", proc );
    fprintf( stderr, "   left corner  = [%9.3f %9.3f %9.3f]\n",
             mincorner->x, mincorner->y, mincorner->z );
    fprintf( stderr, "   right corner = [%9.3f %9.3f %9.3f]\n",
             maxcorner->x, maxcorner->y, maxcorner->z );
  }
  if ( _debug_ >= 3 ) {
    fprintf( stderr, "%s: Computed extremums\n", proc );
    fprintf( stderr, "   left indices  = [%3d %3d %3d]\n",
             minindex.x, minindex.y, minindex.z );
    fprintf( stderr, "   right indices = [%3d %3d %3d]\n",
             maxindex.x, maxindex.y, maxindex.z );
  }
}



static void _updateMaxCorners( bal_floatPoint *pt,
                           bal_floatPoint *mincorner,
                           bal_floatPoint *maxcorner )
{
  if ( mincorner->x > pt->x ) { mincorner->x = pt->x; }
  else if ( maxcorner->x < pt->x ) { maxcorner->x = pt->x; }
  if ( mincorner->y > pt->y ) { mincorner->y = pt->y; }
  else if ( maxcorner->y < pt->y ) { maxcorner->y = pt->y; }
  if ( mincorner->z > pt->z ) { mincorner->z = pt->z; }
  else if ( maxcorner->z < pt->z ) { maxcorner->z = pt->z; }
}





/**********************************************************************
 *
 * box building
 *
 **********************************************************************/


/* the box is the whole field of view
 */

static int _fillBoxWithWholeImage( typeBox *b, bal_image *theIm, bal_transformation *theTrsf )
{
  char *proc = "_fillBoxWithWholeImage";
  bal_floatPoint theCorner[2], thePt, resCorner;
  int x, y, z;

  if ( _warning_ )
    fprintf( stderr, "%s: no QForm compliant yet\n", proc );

  /* corners (voxel coordinates)
   */
  theCorner[0].x = theCorner[0].y = theCorner[0].z = -0.5;
  theCorner[1].x = theIm->ncols - 0.5;
  theCorner[1].y = theIm->nrows - 0.5;
  theCorner[1].z = theIm->nplanes - 0.5;

  /* corners (real coordinates)
   */
  theCorner[0].x *= theIm->vx;
  theCorner[0].y *= theIm->vy;
  theCorner[0].z *= theIm->vz;

  theCorner[1].x *= theIm->vx;
  theCorner[1].y *= theIm->vy;
  theCorner[1].z *= theIm->vz;

  if ( theTrsf == (bal_transformation*)NULL ) {
    b->min = theCorner[0];
    b->max = theCorner[1];
  }
  else {
    for ( z=0; z<2; z++ )
    for ( y=0; y<2; y++ )
    for ( x=0; x<2; x++ ) {
      thePt.x = theCorner[x].x;
      thePt.y = theCorner[y].y;
      thePt.z = theCorner[z].z;
      if ( BAL_TransformFloatPoint( &thePt, &resCorner, theTrsf ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to transform point\n", proc );
        return( -1 );
      }
      if ( x == 0 && y == 0 && z == 0 ) {
        b->min = resCorner;
        b->max = resCorner;
      }
      else {
        _updateMaxCorners( &resCorner, &(b->min), &(b->max) );
      }
    }
  }

  /* field of view
   */
  b->fov.x = b->max.x - b->min.x;
  b->fov.y = b->max.y - b->min.y;
  b->fov.z = b->max.z - b->min.z;

  /* voxel size
   */
  b->voxel.x = theIm->vx;
  b->voxel.y = theIm->vy;
  b->voxel.z = theIm->vz;

  return( 1 );
}





static int _fillBoxWithImage( typeBox *b, bal_image *theIm, bal_transformation *theTrsf, int t )
{
  char *proc = "_fillBoxWithImage";
  bal_floatPoint thePt;
  bal_floatPoint inc[2], theInc, theAdd;
  double *mat;
  size_t i, n=0;
  size_t x, y, z;

  if ( _warning_ )
    fprintf( stderr, "%s: no QForm compliant yet\n", proc );

  /* is the threhold effective ?
   */
  switch( theIm->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    _fillBoxWithWholeImage( b, theIm, theTrsf );
    break;
  case SCHAR :
    if ( t > -128 && t <= 127 ) break;
    _fillBoxWithWholeImage( b, theIm, theTrsf );
    return( 1 );
  case UCHAR :
    if ( t > 0 && t <= 255 ) break;
    _fillBoxWithWholeImage( b, theIm, theTrsf );
    return( 1 );
  case SSHORT :
    if ( t > -32768 && t <= 32767 ) break;
    _fillBoxWithWholeImage( b, theIm, theTrsf );
    return( 1 );
  case USHORT :
    if ( t > 0 && t <= 65535 ) break;
    _fillBoxWithWholeImage( b, theIm, theTrsf );
    return( 1 );
  }




  /* corners (real coordinates of voxel centers)
   */

#define _BOXWITHTHRESHOLDS( TYPE, MIN, MAX ) {   \
  TYPE *theBuf = (TYPE*)theIm->data;             \
  if ( t <= MIN || t > MAX ) {                   \
    _fillBoxWithWholeImage( b, theIm, theTrsf ); \
    return( 1 );                                 \
  }                                              \
  for ( i=0, n=0, z=0; z<theIm->nplanes; z++ )   \
  for ( y=0; y<theIm->nrows; y++ )               \
  for ( x=0; x<theIm->ncols; x++, i++ ) {        \
    if ( theBuf[i] < t ) continue;               \
    thePt.x = (float)x * theIm->vx;              \
    thePt.y = (float)y * theIm->vy;              \
    thePt.z = (float)z * theIm->vz;              \
    if ( theTrsf != (bal_transformation*)NULL ) { \
      if ( BAL_TransformFloatPoint( &thePt, &thePt, theTrsf ) != 1 ) { \
        if ( _verbose_ )                         \
          fprintf( stderr, "%s: unable to transform point\n", proc ); \
         return( -1 );                           \
      }                                          \
    }                                            \
    if ( n == 0 ) {                              \
      b->min = thePt;                            \
      b->max = thePt;                            \
    }                                            \
    else {                                       \
      _updateMaxCorners( &thePt, &(b->min), &(b->max) ); \
    }                                            \
    n ++;                                        \
  }                                              \
}

  switch ( theIm->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    _fillBoxWithWholeImage( b, theIm, theTrsf );
    break;
  case SCHAR :
    _BOXWITHTHRESHOLDS( s8, -128, 127 );
    break;
  case UCHAR :
    _BOXWITHTHRESHOLDS( u8, 0, 255 );
    break;
  case SSHORT :
    _BOXWITHTHRESHOLDS( s16, -32768, 32767 );
    break;
  case USHORT :
    _BOXWITHTHRESHOLDS( u16, 0, 65535 );
    break;
  }

  if ( n == 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: no point above the threshold %d\n", proc, t );
    _fillBoxWithWholeImage( b, theIm, theTrsf );
    return( 1 );
  }

  /* box is defined with voxel centers
   * get a better estimate
   */
  inc[0].x = inc[0].y = inc[0].z = -0.5;
  inc[1].x = inc[1].y = inc[1].z =  0.5;

  inc[0].x *= theIm->vx;
  inc[0].y *= theIm->vy;
  inc[0].z *= theIm->vz;

  inc[1].x *= theIm->vx;
  inc[1].y *= theIm->vy;
  inc[1].z *= theIm->vz;

  for ( z=0; z<2; z++ )
  for ( y=0; y<2; y++ )
  for ( x=0; x<2; x++ ) {
    theInc.x = inc[x].x;
    theInc.y = inc[y].y;
    theInc.z = inc[z].z;
    if ( theTrsf == (bal_transformation*)NULL ) {
      thePt.x = b->min.x + theInc.x;
      thePt.y = b->min.y + theInc.y;
      thePt.z = b->min.z + theInc.z;
      _updateMaxCorners( &thePt, &(b->min), &(b->max) );
      thePt.x = b->max.x + theInc.x;
      thePt.y = b->max.y + theInc.y;
      thePt.z = b->max.z + theInc.z;
      _updateMaxCorners( &thePt, &(b->min), &(b->max) );
    }
    else {
      switch( theTrsf->type ) {
      default :
        if ( _verbose_ )
          fprintf( stderr, "%s: such transformation type not implemented yet\n", proc );
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
        mat = theTrsf->mat.m;
        theAdd.x = mat[ 0] * theInc.x + mat[ 1] * theInc.y + mat[ 2] * theInc.z;
        theAdd.y = mat[ 4] * theInc.x + mat[ 5] * theInc.y + mat[ 6] * theInc.z;
        theAdd.z = mat[ 8] * theInc.x + mat[ 9] * theInc.y + mat[10] * theInc.z;

        thePt.x = b->min.x + theAdd.x;
        thePt.y = b->min.y + theAdd.y;
        thePt.z = b->min.z + theAdd.z;
        _updateMaxCorners( &thePt, &(b->min), &(b->max) );
        thePt.x = b->max.x + theAdd.x;
        thePt.y = b->max.y + theAdd.y;
        thePt.z = b->max.z + theAdd.z;
        _updateMaxCorners( &thePt, &(b->min), &(b->max) );
        break;
      }
    }
  }

  /* field of view
   */
  b->fov.x = b->max.x - b->min.x;
  b->fov.y = b->max.y - b->min.y;
  b->fov.z = b->max.z - b->min.z;

  /* voxel size
   */
  b->voxel.x = theIm->vx;
  b->voxel.y = theIm->vy;
  b->voxel.z = theIm->vz;

  return( 1 );
}










int _fillBoxesWithImageName( typeBoxList *theBox, bufferType *type, char *theName,
                             bal_transformationList *theTrsf, int threshold, int refTrsf,
                             bal_doublePoint *voxelSize )
{
  char *proc = "_fillBoxesWithImageName";
  bal_image theIm;
  int i;
  bal_transformation invTrsf;

  if ( BAL_ReadImage( &theIm, theName, 0 ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to read image '%s'\n", proc, theName );
    return( -1 );
  }

  if ( voxelSize != (bal_doublePoint *)NULL ) {
    if ( voxelSize->x > 0 && voxelSize->y > 0 && voxelSize->z > 0 ) {
      theIm.vx = voxelSize->x;
      theIm.vy = voxelSize->y;
      theIm.vz = voxelSize->z;
    }
  }

  *type = theIm.type;

  BAL_InitTransformation( &invTrsf );
  if ( BAL_AllocTransformation( &invTrsf, theTrsf->data[0].type, (bal_image *)NULL ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate inverse transformation\n", proc );
    return( -1 );
  }


  for ( i=0; i<theTrsf->n_trsfs; i++ ) {

    /* inverse transformation T_{i<-ref}
     * => T_{ref<-i}
     */
    if ( BAL_InverseTransformation(  &(theTrsf->data[i]), &invTrsf ) != 1 ) {
      BAL_FreeTransformation( &invTrsf );
      BAL_FreeImage( &theIm );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to invert transformation #%d\n", proc, i );
      return( -1 );
    }

    /* is there a reference transformation 'r'?
     * if yes, compose T_{r<-ref} o T_{ref<-i} = T_{r<-i}
     */
    if ( 0 <= refTrsf && refTrsf < theTrsf->n_trsfs ) {
      if ( BAL_TransformationComposition( &invTrsf,
                                          &(theTrsf->data[refTrsf]), &invTrsf ) != 1 ) {
        BAL_FreeTransformation( &invTrsf );
        BAL_FreeImage( &theIm );
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to compose with inverse transformation #%d\n", proc, i );
        return( -1 );
      }
    }

    if ( _fillBoxWithImage( &(theBox->data[i]), &theIm, &invTrsf, threshold ) != 1 ) {
      BAL_FreeTransformation( &invTrsf );
      BAL_FreeImage( &theIm );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to fill box #%d\n", proc, i );
      return( -1 );
    }
  }

  BAL_FreeTransformation( &invTrsf );
  BAL_FreeImage( &theIm );
  return( 1 );
}





typedef struct _parallelFillBoxesParam {
    typeBoxList *theBox;
    stringList *theName;
    bal_transformationList *theTrsf;
    int threshold;
    int refTrsf;
    int refTemplate;

    bufferType refType;
} _parallelFillBoxesParam;



static void *_parallelFillBoxesWithImageNames( void *par )
{
  char *proc = "_parallelFillBoxesWithImageNames";
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _parallelFillBoxesParam *p = (_parallelFillBoxesParam*)parameter;

  typeBoxList *theBox = p->theBox;
  stringList *theName = p->theName;
  bal_transformationList *theTrsf = p->theTrsf;

  size_t i;
  bal_image theIm;
  bal_transformation invTrsf;

  BAL_InitTransformation( &invTrsf );
  if ( BAL_AllocTransformation( &invTrsf, theTrsf->data[0].type, (bal_image *)NULL ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate inverse transformation\n", proc );
    chunk->ret = -1;
    return( (void*)NULL );
  }

  for ( i=first; i<=last; i++ ) {

    if ( _debug_ || _verbose_ >= 2 ) {
      fprintf( stderr, "%s: processing '%s'\n", proc, theName->data[i] );
    }

    if ( BAL_ReadImage( &theIm, theName->data[i], 0 ) != 1 ) {
      BAL_FreeTransformation( &invTrsf );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to read image '%s'\n", proc, theName->data[i] );
      chunk->ret = -1;
      return( (void*)NULL );
    }

    /* inverse transformation T_{i<-ref}
     * => T_{ref<-i}
     */
    if ( BAL_InverseTransformation(  &(theTrsf->data[i]), &invTrsf ) != 1 ) {
      BAL_FreeImage( &theIm );
      BAL_FreeTransformation( &invTrsf );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to invert transformation #%lu\n", proc, i );
      chunk->ret = -1;
      return( (void*)NULL );
    }

    /* is there a reference transformation 'r'?
     * if yes, compose T_{r<-ref} o T_{ref<-i} = T_{r<-i}
     */
    if ( 0 <= p->refTrsf && p->refTrsf < theTrsf->n_trsfs ) {
      if ( BAL_TransformationComposition( &invTrsf,
                                          &(theTrsf->data[p->refTrsf]), &invTrsf ) != 1 ) {
        BAL_FreeImage( &theIm );
        BAL_FreeTransformation( &invTrsf );
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to compose with inverse transformation #%lu\n", proc, i );
        chunk->ret = -1;
        return( (void*)NULL );
      }
    }

    if ( _fillBoxWithImage( &(theBox->data[i]), &theIm, &invTrsf, p->threshold ) != 1 ) {
      BAL_FreeImage( &theIm );
      BAL_FreeTransformation( &invTrsf );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to fill box #%lu\n", proc, i );
      chunk->ret = -1;
      return( (void*)NULL );
    }

    if ( (int)i == p->refTemplate ) p->refType = theIm.type;

    BAL_FreeImage( &theIm );
  }

  BAL_FreeTransformation( &invTrsf );

  chunk->ret = 1;
  return( (void*)NULL );
}





int _fillBoxesWithImageNames( typeBoxList *theBox, bufferType *type,
                              stringList *theName,
                              bal_transformationList *theTrsf, int threshold,
                              int refTrsf, int refTemplate )
{
  char *proc = "_fillBoxesWithImageNames";
  typeChunks chunks;
  _parallelFillBoxesParam *aux;
  int n;

  /* initialisations
   * do not forget to set the minimal numbers of elements in one chunk to one
   * (processing one single image per processor is ok!)
   */
  initChunks( &chunks );
  setMinElementsInChunks( 1 );
  if ( buildChunks( &chunks, 0, theTrsf->n_trsfs-1, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return( -1 );
  }

  aux = (_parallelFillBoxesParam*)vtmalloc( chunks.n_allocated_chunks * sizeof(_parallelFillBoxesParam), "aux", proc );
  if ( aux == (_parallelFillBoxesParam*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate auxiliary variables\n", proc );
    freeChunks( &chunks );
    return( -1 );
  }

  for ( n=0; n<chunks.n_allocated_chunks; n++ ) {
    aux[n].theBox = theBox;
    aux[n].theName = theName;
    aux[n].theTrsf = theTrsf;
    aux[n].threshold = threshold;
    aux[n].refTrsf = refTrsf;
    aux[n].refTemplate = refTemplate;
    aux[n].refType = TYPE_UNKNOWN;
    chunks.data[n].parameters = (void*)(&aux[n]);
  }

  if ( processChunks( &_parallelFillBoxesWithImageNames, &chunks, proc ) != 1 ) {
    vtfree( aux );
    freeChunks( &chunks );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to recompute transformations\n", proc );
    return( -1 );
  }

  for ( n=0; n<chunks.n_allocated_chunks; n++ ) {
      if ( aux[n].refType != TYPE_UNKNOWN )
          *type = aux[n].refType;
  }

  vtfree( aux );
  freeChunks( &chunks );

  return( 1 );
}





/**********************************************************************
 *
 *
 *
 **********************************************************************/





/* les transformations de depart sont du type T_{i<-ref}
 *
 * boxList contient les boites englobantes
 * des images d'entree transformees par
 * les T_{r<-ref} o T_{i<-ref}^{-1} = T_{r<-i}
 *
   2. on calcule la translation T_{r<-newref} qui permet de mettre
      toutes les images dans un seul template, avec I_r a une
      frontiere de pixel (translation entiere)
   3. on calcule les transformations resultats
      T_{i<-ref} o T_{r<-ref}^{-1} o T_{r<-newref}
*/


int BAL_ChangeTransformationList( bal_transformationList *theList,
                                  typeBoxList *boxList,
                                  bal_transformationList *resList,
                                  bal_image *resIm,
                                  int refTrsf,
                                  int *margin,
                                  int *extend )
{
  char *proc = "BAL_ChangeTransformationList";

  bal_floatPoint mincorner;
  bal_floatPoint maxcorner;
  bal_floatPoint maxfov;

  bal_integerPoint dim;

  float vx = resIm->vx;
  float vy = resIm->vy;
  float vz = resIm->vz;

  bal_transformation invTrsf;
  bal_transformation tmpTrsf;

  int i;



  /* maximal bounding box
   * - mincorner and maxcorner define the bounding box that
   *   contain all transformed images (in real coordinates)
   * - maxfov contains the maximal fov (in real coordinates)
   */

  _maxCorners( boxList, &mincorner, &maxcorner );



  /* add margins if required
   */

  if ( extend[0] > 0 ) {
    mincorner.x -= margin[0] * vx;
    maxcorner.x += margin[0] * vx;
  }
  if ( extend[1] > 0 ) {
    mincorner.y -= margin[1] * vy;
    maxcorner.y += margin[1] * vy;
  }
  if ( extend[2] > 0 ) {
    mincorner.z -= margin[2] * vz;
    maxcorner.z += margin[2] * vz;
  }

  maxfov.x = maxcorner.x - mincorner.x;
  maxfov.y = maxcorner.y - mincorner.y;
  maxfov.z = maxcorner.z - mincorner.z;



  /* template image dimension computation
   */
  dim.x = (int)( maxfov.x / vx + 0.5 );
  dim.y = (int)( maxfov.y / vy + 0.5 );
  dim.z = (int)( maxfov.z / vz + 0.5 );



  /* template image initialization
   */
  if ( BAL_InitFullImage( resIm, (char *)NULL, dim.x, dim.y, dim.z, 1, vx, vy, vz, resIm->type ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to initialize result template\n", proc );
    return( -1 );
  }



  /* transformations
   */

  BAL_InitTransformation( &tmpTrsf );
  if ( BAL_AllocTransformation( &tmpTrsf, TRANSLATION_3D, (bal_image *)NULL ) != 1 ) {
    BAL_FreeImage( resIm );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate auxiliary transformation\n", proc );
    return( -1 );
  }
  BAL_SetTransformationToIdentity( &tmpTrsf );
  if ( extend[0] > 0 )
      tmpTrsf.mat.m[ 3] = mincorner.x;
  if ( extend[1] > 0 )
      tmpTrsf.mat.m[ 7] = mincorner.y;
  if ( extend[2] > 0 )
      tmpTrsf.mat.m[11] = mincorner.z;


  /* is there a reference transformation ?
   */
  if ( 0 <= refTrsf && refTrsf < theList->n_trsfs ) {

    BAL_InitTransformation( &invTrsf );
    if ( BAL_AllocTransformation( &invTrsf, theList->data[0].type, (bal_image *)NULL ) != 1 ) {
      BAL_FreeImage( resIm );
      BAL_FreeTransformation( &tmpTrsf );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate inverse transformation\n", proc );
      return( -1 );
    }
    if ( BAL_InverseTransformation(  &(theList->data[refTrsf]), &invTrsf ) != 1 ) {
      BAL_FreeImage( resIm );
      BAL_FreeTransformation( &tmpTrsf );
      BAL_FreeTransformation( &invTrsf );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to invert reference transformation #%d\n", proc, refTrsf );
      return( -1 );
    }
    /* compose transformations
       int BAL_TransformationComposition( bal_transformation *t_res,
       bal_transformation *t1,
       bal_transformation *t2 )
       with t_res = t1 o t2
    */
    if ( BAL_TransformationComposition( &tmpTrsf, &invTrsf, &tmpTrsf ) != 1 ) {
      BAL_FreeImage( resIm );
      BAL_FreeTransformation( &tmpTrsf );
      BAL_FreeTransformation( &invTrsf );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compose transformations (preparation)\n", proc );
      return( -1 );
    }

    BAL_FreeTransformation( &invTrsf );

  }


  for ( i=0; i<theList->n_trsfs; i++ ) {
    if ( BAL_TransformationComposition( &(resList->data[i]),
                                        &(theList->data[i]), &tmpTrsf ) != 1 ) {
      BAL_FreeImage( resIm );
      BAL_FreeTransformation( &tmpTrsf );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compose transformations (#%d)\n", proc, i );
      return( -1 );
    }
  }

  BAL_FreeTransformation( &tmpTrsf );

  return( 1 );
}


