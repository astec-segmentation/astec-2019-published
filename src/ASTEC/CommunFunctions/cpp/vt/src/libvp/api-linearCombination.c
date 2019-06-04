/*************************************************************************
 * api-linearCombination.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2019, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Mer  6 mar 2019 11:17:18 CET
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <chunks.h>

#include <vt_copy.h>
#include <vt_image.h>
#include <vt_inrimage.h>

#include <api-linearCombination.h>


static int _verbose_ = 1;

typedef struct {
  void *bufWeight;
  bufferType typeWeight;
  void *bufImage;
  bufferType typeImage;
  r32 *sumWeight;
  r32 *sumImage;
} _addParam;


static void *_addProcedure( void *par )
{
  char *proc = "_addProcedure";
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _addParam *p = (_addParam*)parameter;
  size_t j;

  switch ( p->typeWeight ) {
  default :
    if ( _verbose_ ) {
      fprintf( stderr, "%s: such type not handled for weight image\n", proc );
    }
    chunk->ret = -1;
    return( (void*)NULL );
  case FLOAT :
    {
      r32 *bufWeight = (r32*)p->bufWeight;
      switch( p->typeImage ) {
      default :
        if ( _verbose_ ) {
          fprintf( stderr, "%s: such type not handled for image\n", proc );
        }
        chunk->ret = -1;
        return( (void*)NULL );
      case USHORT :
        {
          u16 *bufImage = (u16*)p->bufImage;
          for ( j=first; j<=last; j++ ) {
            p->sumImage[j] += bufWeight[j] * (r32)bufImage[j];
            p->sumWeight[j] += bufWeight[j];
          }
        }
        break;
      } /* end of switch( theIm->type ) */
    }
    break;
  }
  chunk->ret = 1;
  return( (void*)NULL );
}



static void *_initProcedure( void *par )
{
  char *proc = "_initProcedure";
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _addParam *p = (_addParam*)parameter;
  size_t j;

  switch ( p->typeWeight ) {
  default :
    if ( _verbose_ ) {
      fprintf( stderr, "%s: such type not handled for weight image\n", proc );
    }
    chunk->ret = -1;
    return( (void*)NULL );
  case FLOAT :
    {
      r32 *bufWeight = (r32*)p->bufWeight;
      switch( p->typeImage ) {
      default :
        if ( _verbose_ ) {
          fprintf( stderr, "%s: such type not handled for image\n", proc );
        }
        chunk->ret = -1;
        return( (void*)NULL );
      case USHORT :
        {
          u16 *bufImage = (u16*)p->bufImage;
          for ( j=first; j<=last; j++ ) {
            p->sumImage[j] = bufWeight[j] * (r32)bufImage[j];
            p->sumWeight[j] = bufWeight[j];
          }
        }
        break;
      } /* end of switch( theIm->type ) */
    }
    break;
  }
  chunk->ret = 1;
  return( (void*)NULL );
}


static void *_normProcedure( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _addParam *p = (_addParam*)parameter;
  size_t j;

  for ( j=first; j<=last; j++ ) {
    p->sumImage[j] /= p->sumWeight[j];
  }

  chunk->ret = 1;
  return( (void*)NULL );
}





int API_linearCombination( stringList *weightNames, stringList *imageNames, char *resName, ImageType resType )
{
  char *proc = "API_linearCombination";
  vt_image sumImage;
  vt_image localImage;
  vt_image sumWeight;
  vt_image *theIm;
  vt_image *theWeight;
  int i, j;
  size_t v;
  ImageType localType = TYPE_UNKNOWN;

  _addParam param;
  typeChunks chunks;

  if ( weightNames->n_data <= 0 || (weightNames->n_data != imageNames->n_data) ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: bad numbers of names [weights=%d, images=%d]\n", proc, weightNames->n_data, imageNames->n_data );
    }
    return( -1 );
  }

  VT_Image( &sumImage );
  VT_Image( &sumWeight );
  initChunks( &chunks );

  for ( i=0; i<imageNames->n_data; i++ ) {

    if ( _verbose_ >= 2 )
      fprintf( stderr, "processing couple #%d\n", i );

    if ( _verbose_ >= 3 )
      fprintf( stderr, "   ... reading image '%s' ... ", imageNames->data[i] );

    theIm = _VT_Inrimage( imageNames->data[i] );
    if ( theIm == (vt_image*)NULL ) {
      VT_FreeImage( &sumImage );
      VT_FreeImage( &sumWeight );
      if ( _verbose_ ) {
        fprintf( stderr, "%s: unable to read image '%s'\n", proc, imageNames->data[i] );
      }
      return( -1 );
    }

    if ( _verbose_ >= 3 )
      fprintf( stderr, " done \n");
    if ( _verbose_ >= 3 )
      fprintf( stderr, "   ... reading weight '%s' ... ", weightNames->data[i] );

    theWeight = _VT_Inrimage( weightNames->data[i] );
    if ( theWeight == (vt_image*)NULL ) {
      VT_FreeImage( theIm );
      VT_Free( (void**)&theIm );
      VT_FreeImage( &sumImage );
      VT_FreeImage( &sumWeight );
      if ( _verbose_ ) {
        fprintf( stderr, "%s: unable to read image '%s'\n", proc, weightNames->data[i] );
      }
      return( -1 );
    }

    if ( _verbose_ >= 3 )
      fprintf( stderr, " done \n");

    if ( theIm->dim.x != theWeight->dim.x || theIm->dim.y != theWeight->dim.y || theIm->dim.z != theWeight->dim.z ) {
      VT_FreeImage( theWeight );
      VT_Free( (void**)&theWeight );
      VT_FreeImage( theIm );
      VT_Free( (void**)&theIm );
      VT_FreeImage( &sumImage );
      VT_FreeImage( &sumWeight );
      if ( _verbose_ ) {
        fprintf( stderr, "%s: images '%s' and '%s' have different dimensions\n", proc, imageNames->data[i], weightNames->data[i] );
      }
      return( -1 );
    }



    if ( i == 0 ) {

      localType = theIm->type;
      VT_InitFromImage( &sumImage, theIm, "", FLOAT );
      VT_InitFromImage( &sumWeight, theIm, "", FLOAT );
      v = theIm->dim.x * theIm->dim.y * theIm->dim.z;
      if ( buildChunks( &chunks, 0, v-1, proc ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to compute chunks\n", proc );
        return( -1 );
      }

      if ( VT_AllocImage( &sumImage ) != 1 ) {
        freeChunks( &chunks );
        VT_FreeImage( theWeight );
        VT_Free( (void**)&theWeight );
        VT_FreeImage( theIm );
        VT_Free( (void**)&theIm );
        VT_FreeImage( &sumImage );
        VT_FreeImage( &sumWeight );
        if ( _verbose_ ) {
          fprintf( stderr, "%s: allocation error\n", proc );
        }
        return( -1 );
      }
      if ( VT_AllocImage( &sumWeight ) != 1 ) {
        freeChunks( &chunks );
        VT_FreeImage( theWeight );
        VT_Free( (void**)&theWeight );
        VT_FreeImage( theIm );
        VT_Free( (void**)&theIm );
        VT_FreeImage( &sumImage );
        VT_FreeImage( &sumWeight );
        if ( _verbose_ ) {
          fprintf( stderr, "%s: allocation error\n", proc );
        }
        return( -1 );
      }

      param.bufWeight = theWeight->buf;
      param.typeWeight = theWeight->type;
      param.bufImage = theIm->buf;
      param.typeImage = theIm->type;
      param.sumWeight = (r32*)sumWeight.buf;
      param.sumImage = (r32*)sumImage.buf;

      for ( j=0; j<chunks.n_allocated_chunks; j++ )
        chunks.data[j].parameters = (void*)(&param);

      if ( processChunks( &_initProcedure, &chunks, proc ) != 1 ) {
        freeChunks( &chunks );
        VT_FreeImage( theWeight );
        VT_Free( (void**)&theWeight );
        VT_FreeImage( theIm );
        VT_Free( (void**)&theIm );
        VT_FreeImage( &sumImage );
        VT_FreeImage( &sumWeight );
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to run initialization\n", proc );
        return( -1 );
      }

    } else {

      if ( theIm->dim.x != sumImage.dim.x || theIm->dim.y != sumImage.dim.y || theIm->dim.z != sumImage.dim.z ) {
        freeChunks( &chunks );
        VT_FreeImage( theWeight );
        VT_Free( (void**)&theWeight );
        VT_FreeImage( theIm );
        VT_Free( (void**)&theIm );
        VT_FreeImage( &sumImage );
        VT_FreeImage( &sumWeight );
        if ( _verbose_ ) {
          fprintf( stderr, "%s: images have different dimensions\n", proc );
        }
        return( -1 );
      }

      param.bufWeight = theWeight->buf;
      param.typeWeight = theWeight->type;
      param.bufImage = theIm->buf;
      param.typeImage = theIm->type;
      param.sumWeight = (r32*)sumWeight.buf;
      param.sumImage = (r32*)sumImage.buf;

      for ( j=0; j<chunks.n_allocated_chunks; j++ )
        chunks.data[j].parameters = (void*)(&param);

      if ( processChunks( &_addProcedure, &chunks, proc ) != 1 ) {
        freeChunks( &chunks );
        VT_FreeImage( theWeight );
        VT_Free( (void**)&theWeight );
        VT_FreeImage( theIm );
        VT_Free( (void**)&theIm );
        VT_FreeImage( &sumImage );
        VT_FreeImage( &sumWeight );
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to run initialization\n", proc );
        return( -1 );
      }

    }

    VT_FreeImage( theWeight );
    VT_Free( (void**)&theWeight );
    VT_FreeImage( theIm );
    VT_Free( (void**)&theIm );
  }

  if ( processChunks( &_normProcedure, &chunks, proc ) != 1 ) {
    freeChunks( &chunks );
    VT_FreeImage( &sumImage );
    VT_FreeImage( &sumWeight );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to run initialization\n", proc );
    return( -1 );
  }

  VT_FreeImage( &sumWeight );

  if ( resType != TYPE_UNKNOWN )
    localType = resType;

  if ( sumImage.type == localType ) {
    if ( VT_WriteInrimageWithName( &sumImage, resName ) != 1 ) {
      VT_FreeImage( &sumImage );
      if ( _verbose_ ) {
        fprintf( stderr, "%s: error when writing '%s'\n", proc, resName);
      }
      return( -1 );
    }
  }
  else {
    VT_InitFromImage( &localImage, &sumImage, "", localType );
    if ( VT_AllocImage( &localImage ) != 1 ) {
      VT_FreeImage( &sumImage );
      if ( _verbose_ ) {
        fprintf( stderr, "%s: allocation error\n", proc );
      }
      return( -1 );
    }
    if ( VT_CopyImage( &sumImage, &localImage ) != 1 ) {
      VT_FreeImage( &sumImage );
      if ( _verbose_ ) {
        fprintf( stderr, "%s: error when copying result\n", proc );
      }
      return( -1 );
    }
    if ( VT_WriteInrimageWithName( &localImage, resName ) != 1 ) {
      VT_FreeImage( &sumImage );
      if ( _verbose_ ) {
        fprintf( stderr, "%s: error when writing '%s'\n", proc, resName);
      }
      return( -1 );
    }
    VT_FreeImage( &localImage );
  }


  VT_FreeImage( &sumImage );
  return( 1 );
}
