/*************************************************************************
 * api-cellProperties.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mer 24 oct 2018 16:07:57 CEST
 *
 * ADDITIONS, CHANGES
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <chunks.h>
#include <string-tools.h>
#include <vtmalloc.h>

#include <bal-blockmatching.h>
#include <bal-transformation-compose.h>
#include <bal-transformation-tools.h>

#include <api-cellProperties.h>






static int _verbose_ = 1;
static int _debug_ = 0;


static void _API_ParseParam_cellProperties( char *str, lineCmdParamCellProperties *p );



/************************************************************
 *
 *
 *
 ************************************************************/



void _initUpdateList( typeUpdateList *l )
{
    l->n_allocated_data = 0;
    l->n_data = 0;
    l->data = (int*)NULL;
}



/************************************************************
 *
 * main APIs
 *
 ************************************************************/




static int _writeCellSequenceDiagnosis( char *output_diagnosis_name,
                        typeCellSequence *cellProperties, typePropertyList *properties )
{
  char *proc = "output_diagnosis_name";
  typePropertyList propertyList;
  typePropertyList *propertyPtr = (typePropertyList*)NULL;
  FILE *f = (FILE*)NULL;

  if ( output_diagnosis_name == (char*)NULL || output_diagnosis_name[0] == '\0' ) {
    f = stdout;
  }
  else if ( output_diagnosis_name != (char*)NULL && output_diagnosis_name[0] != '\0' ) {
    f = fopen( output_diagnosis_name, "w" );
    if ( f == (FILE*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when opening '%s'\n", proc, output_diagnosis_name );
      return( -1 );
    }
  }

  if ( f == (FILE*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: diagnosis will be be computed\n", proc );
    return( 1 );
  }

  BAL_InitPropertyList( &propertyList );

  if ( properties == (typePropertyList*)NULL || properties->n_data == 0 ) {
    propertyList.data[propertyList.n_data++] = _VOLUME_;
    propertyList.data[propertyList.n_data++] = _LINEAGE_;
    propertyPtr = &propertyList;
  }
  else {
    propertyPtr = properties;
  }

  if ( BAL_CellSequenceDiagnosis( f, cellProperties, propertyPtr ) != 1 ) {
    if ( f != stdout ) fclose( f );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when computing diagnosis\n", proc );
    return( -1 );
  }

  if ( f != stdout ) fclose( f );

  return( 1 );
}








typedef struct _propertiesParam {
  stringList *segmentationName;
  typePropertyList *propertyList;
  enumSurfaceEstimation surfaceEstimationType;
  typeCellSequence *cellProperties;
} _propertiesParam;



static void * _cellSequenceProperties( void *par )
{
  char *proc = "_cellSequenceProperties";
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _propertiesParam *p = (_propertiesParam*)parameter;
  stringList *segmentationName = p->segmentationName;
  typePropertyList *propertyList = p->propertyList;
  enumSurfaceEstimation surfaceEstimationType = p->surfaceEstimationType;
  typeCellSequence *cellProperties = p->cellProperties;

  size_t i;

  bal_image theSeg;

  BAL_InitImage( &theSeg, NULL, 0, 0, 0, 0, UCHAR );

  for ( i=first; i<=last; i++ ) {

    if ( (int)i >= segmentationName->n_data || (int)i >= cellProperties->n_data )
      continue;

    if ( _verbose_ >= 2 )
      fprintf( stderr, "... properties of image #%3lu\n", i );

    if ( BAL_ReadImage( &theSeg, segmentationName->data[i], 0 ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to read segmentation image '%s'\n", proc, segmentationName->data[i] );
      chunk->ret = -1;
      return( (void*)NULL );
    }

    /* single image properties
     */

    if ( BAL_CellImageProperties( &theSeg, &(cellProperties->data[i]),
                                  propertyList, surfaceEstimationType ) != 1 ) {
        BAL_FreeImage( &theSeg );
        if ( _verbose_ )
          fprintf( stderr, "%s: some error occurs during processing image '%s'\n", proc, segmentationName->data[i]  );
        chunk->ret = -1;
        return( (void*)NULL );
    }

    BAL_FreeImage( &theSeg );

  }

  chunk->ret = 1;
  return( (void*)NULL );
}





typedef struct _lineageParam {
  stringList *fusionName;
  stringList *segmentationName;
  typePropertyList *propertyList;
  typeCellSequence *cellProperties;
  int normalisation;
  bal_blockmatching_pyramidal_param *affineRegistration;
  bal_blockmatching_pyramidal_param *nonlinearRegistration;
} _lineageParam;



static void * _cellSequenceLineage( void *par )
{
  char *proc = "_cellSequenceLineage";
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _lineageParam *p = (_lineageParam*)parameter;
  stringList *fusionName = p->fusionName;
  stringList *segmentationName = p->segmentationName;
  typePropertyList *propertyList = p->propertyList;
  typeCellSequence *cellProperties = p->cellProperties;
  int normalisation = p->normalisation;
  bal_blockmatching_pyramidal_param *affineRegistration = p->affineRegistration;
  bal_blockmatching_pyramidal_param *nonlinearRegistration = p->nonlinearRegistration;

  size_t i;

  bal_image theSeg, nextSeg, resampleSeg;
  bal_image *ptrTheSeg = (bal_image*)NULL;
  bal_image *ptrNextSeg = (bal_image*)NULL;
  bal_image *ptrLineageSeg = (bal_image*)NULL;

  bal_image theFuse, nextFuse;
  bal_image *ptrTheFuse = (bal_image*)NULL;
  bal_image *ptrNextFuse = (bal_image*)NULL;

  bal_image *ptrTmp = (bal_image*)NULL;

  bal_transformation *ptrAffineTrsf = (bal_transformation *)NULL;
  bal_transformation *ptrNonlinearTrsf = (bal_transformation *)NULL;
  bal_transformation *ptrResampleTrsf;
  bal_transformation theComposedTrsf;
  int ptrResampleTrsfIsAllocated = 0;



  BAL_InitImage( &theSeg, NULL, 0, 0, 0, 0, UCHAR );
  BAL_InitImage( &nextSeg, NULL, 0, 0, 0, 0, UCHAR );

  ptrTheSeg = &theSeg;
  ptrNextSeg = &nextSeg;

  BAL_InitImage( &theFuse, NULL, 0, 0, 0, 0, UCHAR );
  BAL_InitImage( &nextFuse, NULL, 0, 0, 0, 0, UCHAR );

  ptrTheFuse = &theFuse;
  ptrNextFuse = &nextFuse;


  for ( i=first; i<=last; i++ ) {

    if ( (int)i >= segmentationName->n_data || (int)i >= cellProperties->n_data )
      continue;

    if ( _verbose_ >= 2 )
      fprintf( stderr, "... lineage of image #%3lu\n", i );

    /* read first segmentation image
     * or get it from previous reading
     */
    if ( i == first || ptrNextSeg == (bal_image*)NULL ) {
      if ( BAL_ReadImage( &theSeg, segmentationName->data[i], 0 ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to read segmentation image '%s'\n", proc, segmentationName->data[i] );
        chunk->ret = -1;
        return( (void*)NULL );
      }
      ptrTheSeg = &theSeg;
    }
    else {
      ptrTmp = ptrTheSeg;
      ptrTheSeg = ptrNextSeg;
      ptrNextSeg = ptrTmp;
    }

    /* at this point:
     * - ptrTheSeg is allocated
     */

    /* lineage
     */

    if ( _isPropertyInList( propertyList, _LINEAGE_ ) == 1 ) {
      if ( cellProperties->data[i].next_acquisition_time > cellProperties->data[i].acquisition_time ) {


        /*** transformation computation (if required)
         ***/

        ptrResampleTrsf = (bal_transformation *)NULL;

        if ( fusionName->n_data > 0 ) {

          if ( _verbose_ >= 2 )
            fprintf( stderr, "    linear registration of image #%3lu\n", i );

           /* read first fusion image
            * or get it from previous reading
            */
          if ( i == first || ptrNextFuse == (bal_image*)NULL ) {
            if ( BAL_ReadImage( &theFuse, fusionName->data[i], normalisation ) != 1 ) {
              BAL_FreeImage( ptrTheSeg );
              if ( _verbose_ )
                fprintf( stderr, "%s: unable to read fusion image '%s'\n", proc, fusionName->data[i] );
              chunk->ret = -1;
              return( (void*)NULL );
            }
            ptrTheFuse = &theFuse;
          }
          else {
            ptrTmp = ptrTheFuse;
            ptrTheFuse = ptrNextFuse;
            ptrNextFuse = ptrTmp;
          }

          /* at this point:
           * - ptrTheSeg is allocated
           * - ptrTheFuse is allocated
           */

          /* read next fusion image
           */
          if ( i == first ) ptrTmp = &nextFuse;
          else ptrTmp = ptrNextFuse;

          if ( BAL_ReadImage( ptrTmp, fusionName->data[i+1], 0 ) != 1 ) {
            BAL_FreeImage( ptrTheFuse );
            BAL_FreeImage( ptrTheSeg );
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to read fusion image '%s'\n", proc, fusionName->data[i+1] );
            chunk->ret = -1;
            return( (void*)NULL );
          }

          /* linear registration of both image
           * reference = fusionName->data[i]
           * floating = fusionName->data[i+1]
           */

          ptrAffineTrsf = BAL_PyramidalBlockMatching( ptrTheFuse, ptrNextFuse,
                                                             (bal_transformation *)NULL,
                                                             (bal_transformation *)NULL,
                                                             affineRegistration );
          if ( ptrAffineTrsf == (bal_transformation*)NULL ) {
            BAL_FreeImage( ptrTheFuse );
            BAL_FreeImage( ptrNextFuse );
            BAL_FreeImage( ptrTheSeg );
            if ( _verbose_ )
              fprintf( stderr, "%s: error during linear registration (ref=%s)\n", proc, fusionName->data[i] );
            chunk->ret = -1;
            return( (void*)NULL );
          }

          if ( 0 ) {
            char trsfname[256];
            sprintf( trsfname, "affine%03lu.trsf", i );
            (void)BAL_WriteTransformation( ptrAffineTrsf, trsfname );
          }

          if ( nonlinearRegistration->max_iterations.lowest > 0 || nonlinearRegistration->max_iterations.highest > 0 ) {
            /* non-linear registration
             */
            if ( _verbose_ >= 2 )
              fprintf( stderr, "    non-linear registration of image #%3lu\n", i );

            ptrNonlinearTrsf = BAL_PyramidalBlockMatching( ptrTheFuse, ptrNextFuse,
                                                           ptrAffineTrsf,
                                                           (bal_transformation *)NULL,
                                                           nonlinearRegistration );
            if ( ptrNonlinearTrsf == (bal_transformation*)NULL ) {
              BAL_FreeImage( ptrTheFuse );
              BAL_FreeImage( ptrNextFuse );
              BAL_FreeImage( ptrTheSeg );
              BAL_FreeTransformation( ptrAffineTrsf );
              vtfree( ptrAffineTrsf );
              if ( _verbose_ )
                fprintf( stderr, "%s: error during non-linear registration (ref=%s)\n", proc, fusionName->data[i] );
              chunk->ret = -1;
              return( (void*)NULL );
            }
            /* composition
             */
            BAL_InitTransformation( &theComposedTrsf );
            if ( BAL_AllocTransformationComposition( &theComposedTrsf,
                                                     ptrAffineTrsf,
                                                     ptrNonlinearTrsf,
                                                     ptrTheFuse ) != 1 ) {
              BAL_FreeImage( ptrTheFuse );
              BAL_FreeImage( ptrNextFuse );
              BAL_FreeImage( ptrTheSeg );
              BAL_FreeTransformation( ptrAffineTrsf );
              vtfree( ptrAffineTrsf );
              BAL_FreeTransformation( ptrNonlinearTrsf );
              vtfree( ptrNonlinearTrsf );
              if ( _verbose_ )
                fprintf( stderr, "%s: error during transformation allocation (ref=%s)\n", proc, fusionName->data[i] );
              chunk->ret = -1;
              return( (void*)NULL );
            }
            if ( BAL_TransformationComposition( &theComposedTrsf,
                                                ptrAffineTrsf,
                                                ptrNonlinearTrsf ) != 1 ) {
              BAL_FreeImage( ptrTheFuse );
              BAL_FreeImage( ptrNextFuse );
              BAL_FreeImage( ptrTheSeg );
              BAL_FreeTransformation( ptrAffineTrsf );
              vtfree( ptrAffineTrsf );
              BAL_FreeTransformation( ptrNonlinearTrsf );
              vtfree( ptrNonlinearTrsf );
              BAL_FreeTransformation( &theComposedTrsf );
              if ( _verbose_ )
                fprintf( stderr, "%s: error during transformation composition (ref=%s)\n", proc, fusionName->data[i] );
              chunk->ret = -1;
              return( (void*)NULL );
            }
            BAL_FreeTransformation( ptrAffineTrsf );
            vtfree( ptrAffineTrsf );
            BAL_FreeTransformation( ptrNonlinearTrsf );
            vtfree( ptrNonlinearTrsf );
            ptrResampleTrsf = &theComposedTrsf;
            ptrResampleTrsfIsAllocated = 0;
          }
          else {
              ptrResampleTrsf = ptrAffineTrsf;
              ptrResampleTrsfIsAllocated = 1;
          }

          if ( 0 ) {
            char img_name[256];
            bal_image fooImg;
            if ( nonlinearRegistration->max_iterations.lowest > 0 || nonlinearRegistration->max_iterations.highest > 0 )
              sprintf( img_name, "resampled-non-linear-%03lu.mha", i+1 );
            else
              sprintf( img_name, "resampled-linear-%03lu.mha", i+1 );
            (void)BAL_AllocImageFromImage( &fooImg, img_name, ptrTheFuse, ptrNextFuse->type );
            (void)BAL_ResampleImage( ptrNextFuse, &fooImg, ptrResampleTrsf, LINEAR );
            (void)BAL_WriteImage( &fooImg, img_name );
          }

          BAL_FreeImage( ptrTheFuse );
        } /* if ( fusionName->n_data > 0 ) */

        /* at this point:
         * - ptrTheSeg is allocated
         * if fusion image are given
         * - ptrNextFuse is allocated
         * - ptrResampleTrsf is allocated
         * - the pointer ptrAffineTrsf may be allocated
         */

        /*** read next segmentation image
         ***/

        if ( i == first ) ptrTmp = &nextSeg;
        else ptrTmp = ptrNextSeg;

        if ( BAL_ReadImage( ptrTmp, segmentationName->data[i+1], 0 ) != 1 ) {
          if ( fusionName->n_data > 0 ) {
            BAL_FreeImage( ptrNextFuse );
            BAL_FreeTransformation( ptrResampleTrsf );
            if ( ptrResampleTrsfIsAllocated ) vtfree( ptrResampleTrsf );
          }
          BAL_FreeImage( ptrTheSeg );
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to read image '%s'\n", proc, segmentationName->data[i+1] );
          chunk->ret = -1;
          return( (void*)NULL );
        }

        /*** resample next segmentation image (if required)
         ***/

        ptrLineageSeg = (bal_image*)NULL;

        if ( fusionName->n_data > 0 ) {

            if ( _verbose_ >= 2 )
                fprintf( stderr, "    resampling of image #%3lu\n", i+1 );

            if ( BAL_AllocImageFromImage( &resampleSeg, "resampledSegmentation.mha", ptrTheSeg, ptrNextSeg->type ) != 1 ) {
            BAL_FreeImage( ptrNextFuse );
            BAL_FreeTransformation( ptrResampleTrsf );
            if ( ptrResampleTrsfIsAllocated ) vtfree( ptrResampleTrsf );
            BAL_FreeImage( ptrNextSeg );
            BAL_FreeImage( ptrTheSeg );
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to allocate resampled segmentation image\n", proc );
            chunk->ret = -1;
            return( (void*)NULL );
          }
          /* if cell images are to be smoothed
           * BAL_ResampleCellImage() (from blockmatching/bal-transformation-cell.c) should be used
           */
          if ( BAL_ResampleImage( ptrNextSeg, &resampleSeg, ptrResampleTrsf, NEAREST ) != 1 ) {
            BAL_FreeImage( &resampleSeg );
            BAL_FreeImage( ptrNextFuse );
            BAL_FreeTransformation( ptrResampleTrsf );
            if ( ptrResampleTrsfIsAllocated ) vtfree( ptrResampleTrsf );
            BAL_FreeImage( ptrNextSeg );
            BAL_FreeImage( ptrTheSeg );
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to resample segmentation image\n", proc );
            chunk->ret = -1;
            return( (void*)NULL );
          }
          ptrLineageSeg = &resampleSeg;
        }
        else {
          ptrLineageSeg = ptrNextSeg;
        }



        /* at this point:
         * - ptrTheSeg is allocated
         * - ptrNextSeg is allocated
         * if fusion image are given
         * - ptrNextFuse is allocated (kept for not re-reading the image)
         * - resampleSeg is allocated
         * - ptrResampleTrsf is allocated
         * - the pointer ptrAffineTrsf may be allocated
         */

        /*** compute forward intersection
         ***/

        if ( _verbose_ >= 2 && fusionName->n_data > 0 )
            fprintf( stderr, "    lineage of image #%3lu\n", i );

        if ( BAL_CellImageForwardIntersection( ptrTheSeg, ptrLineageSeg,
                                   &(cellProperties->data[i]) ) != 1 ) {
            if ( fusionName->n_data > 0 ) {
              BAL_FreeImage( &resampleSeg );
              BAL_FreeImage( ptrNextFuse );
              BAL_FreeTransformation( ptrResampleTrsf );
              if ( ptrResampleTrsfIsAllocated ) vtfree( ptrResampleTrsf );
            }
            BAL_FreeImage( ptrTheSeg );
            BAL_FreeImage( ptrNextSeg );
            if ( _verbose_ )
              fprintf( stderr, "%s: some error occurs during processing images ['%s', '%s']\n",
                       proc, segmentationName->data[i], segmentationName->data[i+1] );
            chunk->ret = -1;
            return( (void*)NULL );
        }

        /* release memory
         */
        if ( fusionName->n_data > 0 ) {
          BAL_FreeImage( &resampleSeg );
          BAL_FreeTransformation( ptrResampleTrsf );
          if ( ptrResampleTrsfIsAllocated ) vtfree( ptrResampleTrsf );
        }

        /* end of lineage calculation (two conditions)
         */
      }
    }

    /* release current image
     */
    BAL_FreeImage( ptrTheSeg );

    /* at this point:
     * - ptrNextSeg is allocated
     * if fusion image are given
     * - ptrNextFuse is allocated (kept for not re-reading the image)
     */
  }

  BAL_FreeImage( ptrNextSeg );
  if ( fusionName->n_data > 0 ) {
    BAL_FreeImage( ptrNextFuse );
  }

  chunk->ret = 1;
  return( (void*)NULL );
}





int API_INTERMEDIARY_sequenceCellProperties( char *theformat_segmentation,
                                             char *theformat_fusion,
                                             int firstindex, int lastindex,
                                             char *output_xml_name,
                                             char *output_diagnosis_name,
                                             char *param_str_1, char *param_str_2 )
{
  char *proc = "API_INTERMEDIARY_sequenceCellProperties";
  lineCmdParamCellProperties par;
  stringList theSegmentationList;
  stringList theFusionList;
  typeCellSequence cellProperties;
  typeChunks chunks;
  int maxchunks;
  _propertiesParam propertiesParam;
  _lineageParam lineageParam;
  int i;
  FILE *f;



  /* parameter initialization
   */
  API_InitParam_cellProperties( &par );

  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_cellProperties( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_cellProperties( param_str_2, &par );

  if ( firstindex < 0 || lastindex < 0 ) {
    API_FreeParam_cellProperties( &par );
    if ( _verbose_)
      fprintf( stderr, "%s: index range [%d,%d] has negative values\n", proc, firstindex, lastindex );
    return( -1 );
  }

  if ( par.output.n_data == 0 ) BAL_DefaultPropertyList( &(par.output) );



  /* for update purposes (to get the same index for names and cell image properties)
   * the numbering of names should begin with 0
   */
  initStringList( &theSegmentationList );
  if ( buildStringListFromFormat( theformat_segmentation, 0, lastindex, &theSegmentationList ) != 1 ) {
      API_FreeParam_cellProperties( &par );
      if ( _verbose_)
        fprintf( stderr, "%s: unable to build segmentation image list\n", proc);
      return( -1 );
  }

  if ( _debug_ >= 1 )
      printStringList( stderr, &theSegmentationList, (char*)NULL );

  initStringList( &theFusionList );
  if ( theformat_fusion != (char*)NULL && theformat_fusion[0] != '\0' ) {
    if ( buildStringListFromFormat( theformat_fusion, 0, lastindex, &theFusionList ) != 1 ) {
      freeStringList( &theSegmentationList );
      API_FreeParam_cellProperties( &par );
      if ( _verbose_)
        fprintf( stderr, "%s: unable to build fusion image list\n", proc);
      return( -1 );
    }
    BAL_SetVerboseInBalBlockMatching( 0 );
    BAL_SetVerboseInBalImage( 0 );
  }


  /************************************************************
   *
   *  here is the stuff
   *
   ************************************************************/

  BAL_InitCellSequence( &cellProperties );
  if ( BAL_AllocCellSequence( &cellProperties, lastindex+1 ) != 1 ) {
    freeStringList( &theFusionList );
    freeStringList( &theSegmentationList );
    API_FreeParam_cellProperties( &par );
    if ( _verbose_)
      fprintf( stderr, "%s: allocation of property list failed\n", proc);
    return( -1 );
  }

  /* this has to be rethought,
   * in case there is a gap larger than 1 between two time points
   */
  for ( i=0; i<=lastindex; i++ ) {
      if ( i > 0 )
        cellProperties.data[i].previous_acquisition_time = i-1;
      else
        cellProperties.data[i].previous_acquisition_time = 0;
      cellProperties.data[i].acquisition_time = i;
      if ( i < lastindex )
        cellProperties.data[i].next_acquisition_time = i+1;
      else
        cellProperties.data[i].next_acquisition_time = 0;
  }

  /* there are lastindex-firstindex+1 images to be processed,
   * do not forget to set the minimal numbers of elements in one chunk to one
   * (processing one single image per processor is ok!)
   *
   * - image names are indexed from 0 to lastindex
   *   but only the range [ firstindex, lastindex ] has to be processed
   * - cell images properties are numbered from 0 to lastindex
   *   the 'firstIndex' first properties are then empty lists
   */



  /****** computation of still image properties
   * => we can parallelized on images
   ****** computation of lineages
   * => we can parallelized if there is no registration
   * (registration already occurs in registration itself, but also
   *  in resampling, image filtering, ...)
   */

  maxchunks = getMaxChunks( );
  if ( par.chunks_for_properties > 0 ) {
    setMaxChunks( par.chunks_for_properties );
  }

  initChunks( &chunks );
  setMinElementsInChunks( 1 );

  if ( buildChunks( &chunks, firstindex, lastindex, proc ) != 1 ) {
      BAL_FreeCellSequence( &cellProperties );
      freeStringList( &theFusionList );
      freeStringList( &theSegmentationList );
      API_FreeParam_cellProperties( &par );
      if ( _verbose_)
        fprintf( stderr, "%s: unable to compute chunks\n", proc);
      return( -1 );
  }

  propertiesParam.segmentationName = &theSegmentationList;
  propertiesParam.propertyList = &(par.output);
  propertiesParam.surfaceEstimationType = par.surfaceEstimationType;
  propertiesParam.cellProperties = &cellProperties;

  for ( i=0; i<chunks.n_allocated_chunks; i++ )
    chunks.data[i].parameters = (void*)(&propertiesParam);

  if ( processChunks( &_cellSequenceProperties, &chunks, proc ) != 1 ) {
    freeChunks( &chunks );
    freeStringList( &theFusionList );
    BAL_FreeCellSequence( &cellProperties );
    freeStringList( &theSegmentationList );
    API_FreeParam_cellProperties( &par );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute cell properties\n", proc );
    return( -1 );
  }

  freeChunks( &chunks );

  if ( par.chunks_for_properties > 0 ) {
    setMaxChunks( maxchunks );
  }

  /****** computation of lineage
   */

  if ( _isPropertyInList( &(par.output), _LINEAGE_ ) == 1
       || _isPropertyInList( &(par.output), _FORWARD_NEIGHBORS_ ) == 1
       || _isPropertyInList( &(par.output), _BACKWARD_NEIGHBORS_ ) == 1 ) {

    maxchunks = getMaxChunks( );

    initChunks( &chunks );
    setMinElementsInChunks( 1 );

    if ( theFusionList.n_data > 1 ) {
      /* il y a les images de fusion, donc il y aura du recalage
       * a faire.
       * on ne construit qu'un seul chunk avec toutes les images
       */
      if ( buildChunks( &chunks, 0, 0, proc ) != 1 ) {
        BAL_FreeCellSequence( &cellProperties );
        freeStringList( &theFusionList );
        freeStringList( &theSegmentationList );
        API_FreeParam_cellProperties( &par );
        if ( _verbose_)
          fprintf( stderr, "%s: unable to compute chunks\n", proc);
        return( -1 );
      }
      chunks.data[0].first = firstindex;
      chunks.data[0].last = lastindex;
    }
    else {
     /* pas de fusion, pas de recalage,
      * on peut paralleliser
      * on modifie le nombre maximum de chunks si demande
      */
     if ( par.chunks_for_properties > 0 ) {
        setMaxChunks( par.chunks_for_properties );
      }
      if ( buildChunks( &chunks, firstindex, lastindex, proc ) != 1 ) {
        BAL_FreeCellSequence( &cellProperties );
        freeStringList( &theFusionList );
        freeStringList( &theSegmentationList );
        API_FreeParam_cellProperties( &par );
        if ( _verbose_)
          fprintf( stderr, "%s: unable to compute chunks\n", proc);
        return( -1 );
      }
    }

    lineageParam.fusionName = &theFusionList;
    lineageParam.segmentationName = &theSegmentationList;
    lineageParam.propertyList = &(par.output);
    lineageParam.cellProperties = &cellProperties;
    lineageParam.normalisation = par.normalisation;
    lineageParam.affineRegistration = &(par.affineRegistration);
    lineageParam.nonlinearRegistration = &(par.nonlinearRegistration);

    for ( i=0; i<chunks.n_allocated_chunks; i++ )
      chunks.data[i].parameters = (void*)(&lineageParam);

    if ( processChunks( &_cellSequenceLineage, &chunks, proc ) != 1 ) {
      freeChunks( &chunks );
      freeStringList( &theFusionList );
      BAL_FreeCellSequence( &cellProperties );
      freeStringList( &theSegmentationList );
      API_FreeParam_cellProperties( &par );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute cell lineage\n", proc );
      return( -1 );
    }

    for ( i=0; i<chunks.n_allocated_chunks; i++ ) {
      if ( chunks.data[i].ret == -1 ) {
        if ( _isPropertyInList( &(par.output), _LINEAGE_ ) == 1 ) {
          fprintf( stderr, "%s: unable to compute cell lineage for chunk #%d\n", proc, i );
          fprintf( stderr, "\t remove 'lineage' from properties to be printed out\n" );
          _removePropertyFromList( &(par.output), _LINEAGE_ );
        }
      }
    }

    freeChunks( &chunks );

    if ( theFusionList.n_data <= 0 && par.chunks_for_properties > 0 ) {
      setMaxChunks( maxchunks );
    }


    /* some post-processing
     */

    if ( BAL_CellSequenceBackwardIntersection( &cellProperties ) != 1 ) {
      BAL_FreeCellSequence( &cellProperties );
      API_FreeParam_cellProperties( &par );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when extracting backward intersections\n", proc );
      return( -1 );
    }
  }


  /* release some memory
   */
  freeStringList( &theFusionList );
  freeStringList( &theSegmentationList );


  /* ...
   */

  if ( 0 ) {
      f = fopen( "sequence_whole.txt", "w" );
      BAL_FprintfCellSequence( f, &cellProperties );
      fclose( f );
  }


  /* diagnosis
   */

  if ( output_diagnosis_name != (char*)NULL && output_diagnosis_name[0] != '\0' ) {
    if ( _writeCellSequenceDiagnosis( output_diagnosis_name, &cellProperties, &(par.output) ) != 1 ) {
        BAL_FreeCellSequence( &cellProperties );
        API_FreeParam_cellProperties( &par );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when computing diagnosis\n", proc );
        return( -1 );
    }
  }

  /* output writing
   */

  if ( par.allow_stdin_stdout != 0 &&
       (output_xml_name == (char*)NULL || output_xml_name[0] == '\0'
        || (output_xml_name[0] == '>' && output_xml_name[1] == '\0')) ) {
    f = stdout;
  }
  else {
    f = fopen( output_xml_name, "w" );
    if ( f == (FILE*)NULL ) {
      BAL_FreeCellSequence( &cellProperties );
      API_FreeParam_cellProperties( &par );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when opening '%s'\n", proc, output_xml_name );
      return( -1 );
    }
  }

  XML_FprintfCellSequenceProperty( f, &cellProperties, &(par.output) );

  if ( par.allow_stdin_stdout != 0 &&
       (output_xml_name == (char*)NULL || output_xml_name[0] == '\0'
        || (output_xml_name[0] == '>' && output_xml_name[1] == '\0')) ) {
    ;
  }
  else {
    fclose( f );
  }

  BAL_FreeCellSequence( &cellProperties );
  API_FreeParam_cellProperties( &par );

  return( 1 );
}





int API_INTERMEDIARY_sequenceCellPropertiesUpdate( char *theformat_segmentation,
                                                   char *theformat_fusion,
                                                   char *input_xml_name,
                                                   char *output_xml_name,
                                                   char *output_diagnosis_name,
                                                   char *param_str_1, char *param_str_2 )
{
  char *proc = "API_INTERMEDIARY_sequenceCellPropertiesUpdate";
  lineCmdParamCellProperties par;
  stringList theSegmentationList;
  stringList theFusionList;
  typeCellSequence cellProperties;
  typePropertyList readProperties;
  typePropertyList writeProperties;
  int lastindex;
  FILE *f;

  typePropertyList lineageDaughterProperty;
  int i;
  typeChunks chunks;
  _propertiesParam propertiesParam;
  _lineageParam lineageParam;


  /* parameter initialization
   */
  API_InitParam_cellProperties( &par );

  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_cellProperties( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_cellProperties( param_str_2, &par );



  /* read input properties
   */

  BAL_InitCellSequence( &cellProperties );
  BAL_InitPropertyList( &readProperties );

  f = fopen( input_xml_name, "r" );
  if ( f == (FILE*)NULL ) {
    BAL_FreeCellSequence( &cellProperties );
    freeStringList( &theSegmentationList );
    API_FreeParam_cellProperties( &par );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when opening '%s'\n", proc, input_xml_name );
    return( -1 );
  }

  if ( XML_FscanfCellSequenceProperty( f, &cellProperties, &readProperties ) != 1 ) {
      fclose( f );
      BAL_FreeCellSequence( &cellProperties );
      freeStringList( &theSegmentationList );
      API_FreeParam_cellProperties( &par );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when reading '%s'\n", proc, input_xml_name );
      return( -1 );
  }

  fclose( f );

  lastindex = cellProperties.n_data-1;

  if ( 0 ) {
      f = fopen( "sequence_read.txt", "w" );
      BAL_FprintfPropertyList( f, &readProperties, (char*)NULL );
      BAL_FprintfCellSequence( f, &cellProperties );
      fclose( f );
  }


  /* this has to be rethought,
   * in case there is a gap larger than 1 between two time points
   */
  for ( i=0; i<=lastindex; i++ ) {
      if ( i > 0 )
        cellProperties.data[i].previous_acquisition_time = i-1;
      else
        cellProperties.data[i].previous_acquisition_time = 0;
      cellProperties.data[i].acquisition_time = i;
      if ( i < lastindex )
        cellProperties.data[i].next_acquisition_time = i+1;
      else
        cellProperties.data[i].next_acquisition_time = 0;
  }



  /* for update purposes (to get the same index for names and cell image properties)
   * the numbering of names should begin with 0
   */
  initStringList( &theSegmentationList );
  if ( buildStringListFromFormat( theformat_segmentation, 0, lastindex, &theSegmentationList ) != 1 ) {
      API_FreeParam_cellProperties( &par );
      if ( _verbose_)
        fprintf( stderr, "%s: unable to build input image list\n", proc);
      return( -1 );
  }

  if ( _debug_ >= 1 )
      printStringList( stderr, &theSegmentationList, (char*)NULL );

  initStringList( &theFusionList );
  if ( theformat_fusion != (char*)NULL && theformat_fusion[0] != '\0' ) {
    if ( buildStringListFromFormat( theformat_fusion, 0, lastindex, &theFusionList ) != 1 ) {
      freeStringList( &theSegmentationList );
      API_FreeParam_cellProperties( &par );
      if ( _verbose_)
        fprintf( stderr, "%s: unable to build fusion image list\n", proc);
      return( -1 );
    }
    BAL_SetVerboseInBalBlockMatching( 0 );
    BAL_SetVerboseInBalImage( 0 );
  }



  /************************************************************
   *
   *  here is the stuff
   *
   ************************************************************/




  /* re-init selected images of the sequence
   */
  if ( BAL_InitCellSequencePartial( &cellProperties, &readProperties, &(par.updateList) ) != 1 ) {
    BAL_FreeCellSequence( &cellProperties );
    freeStringList( &theSegmentationList );
    API_FreeParam_cellProperties( &par );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when partially initializing sequence cell properties\n", proc );
    return( -1 );
  }

  if ( 0 ) {
      f = fopen( "sequence_reinit.txt", "w" );
      BAL_FprintfCellSequence( f, &cellProperties );
      fclose( f );
  }


  /* recompute cell properties
   */
  BAL_InitPropertyList( &lineageDaughterProperty );
  if ( _isPropertyInList( &readProperties, _LINEAGE_ ) ) {
    lineageDaughterProperty.data[0] = _LINEAGE_;
    lineageDaughterProperty.n_data = 1;
  }

  initChunks( &chunks );
  setMinElementsInChunks( 1 );
  if ( buildChunks( &chunks, 0, 0, proc ) != 1 ) {
      BAL_FreeCellSequence( &cellProperties );
      freeStringList( &theSegmentationList );
      API_FreeParam_cellProperties( &par );
      if ( _verbose_)
        fprintf( stderr, "%s: unable to compute chunks\n", proc);
      return( -1 );
  }

  propertiesParam.segmentationName = &theSegmentationList;
  propertiesParam.propertyList = &(par.output);
  propertiesParam.surfaceEstimationType = par.surfaceEstimationType;
  propertiesParam.cellProperties = &cellProperties;

  lineageParam.fusionName = &theFusionList;
  lineageParam.segmentationName = &theSegmentationList;
  lineageParam.propertyList = &(par.output);
  lineageParam.cellProperties = &cellProperties;
  lineageParam.normalisation = par.normalisation;
  lineageParam.affineRegistration = &(par.affineRegistration);
  lineageParam.nonlinearRegistration = &(par.nonlinearRegistration);




  for ( i=0; i<par.updateList.n_data; i++ ) {

      /***** properties computation of #i
       */
      propertiesParam.propertyList = &readProperties;
      chunks.data[0].first = par.updateList.data[i];
      chunks.data[0].last = par.updateList.data[i];
      chunks.data[0].parameters = (void*)(&propertiesParam);

      if ( _verbose_ >= 2 )
          fprintf( stderr, "%s: computation of image #%d properties\n",
                   proc, par.updateList.data[i] );

      if ( processChunks( &_cellSequenceProperties, &chunks, proc ) != 1 ) {
        freeChunks( &chunks );
        freeStringList( &theFusionList );
        BAL_FreeCellSequence( &cellProperties );
        freeStringList( &theSegmentationList );
        API_FreeParam_cellProperties( &par );
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to compute cell properties\n", proc );
        return( -1 );
      }


      if ( _isPropertyInList( &readProperties, _LINEAGE_ ) ) {

        /***** lineage computation of #i
         */
        lineageParam.propertyList = &(par.output);
        chunks.data[0].parameters = (void*)(&lineageParam);

        if ( _verbose_ >= 2 )
            fprintf( stderr, "%s: computation of image #%d lineage\n",
                     proc, par.updateList.data[i] );

        if ( processChunks( &_cellSequenceLineage, &chunks, proc ) != 1 ) {
          freeChunks( &chunks );
          freeStringList( &theFusionList );
          BAL_FreeCellSequence( &cellProperties );
          freeStringList( &theSegmentationList );
          API_FreeParam_cellProperties( &par );
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to compute cell properties\n", proc );
          return( -1 );
        }

        /****** lineage computation of #i-1
         * forward lineage (list of daughter cells) is to be recomputed
         * for image par.updateList.data[i]-1
         */
        if ( par.updateList.data[i] > 0
             && _isUpdateInList( &(par.updateList), par.updateList.data[i]-1 ) == 0
             && cellProperties.data[ par.updateList.data[i]-1 ].n_data > 0 ) {

          lineageParam.propertyList = &lineageDaughterProperty;
          chunks.data[0].first = par.updateList.data[i] - 1;
          chunks.data[0].last = par.updateList.data[i] - 1;
          chunks.data[0].parameters = (void*)(&lineageParam);

          if ( _verbose_ >= 2 )
              fprintf( stderr, "%s: computation of image #%d lineage\n",
                       proc, par.updateList.data[i]-1 );

          if ( processChunks( &_cellSequenceLineage, &chunks, proc ) != 1 ) {
            freeChunks( &chunks );
            freeStringList( &theFusionList );
            BAL_FreeCellSequence( &cellProperties );
            freeStringList( &theSegmentationList );
            API_FreeParam_cellProperties( &par );
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to compute cell properties\n", proc );
            return( -1 );
          }
        }

      }
  }

  freeChunks( &chunks );
  freeStringList( &theFusionList );
  freeStringList( &theSegmentationList );



  /* some post-processing
   */

  for ( i=0; i<par.updateList.n_data; i++ ) {
    /* update of backward lineage for image par.updateList.data[i]
     * from image par.updateList.data[i]-1
     */
    if ( par.updateList.data[i] > 0 && cellProperties.data[ par.updateList.data[i]+1 ].n_data > 0 ) {
      if ( BAL_CellImageBackwardIntersection( &(cellProperties.data[ par.updateList.data[i]-1 ]),
                                         &(cellProperties.data[ par.updateList.data[i] ]) ) != 1 ) {
        BAL_FreeCellSequence( &cellProperties );
        API_FreeParam_cellProperties( &par );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when extracting backward lineage of image #%d\n", proc, par.updateList.data[i] );
        return( -1 );
      }

    }
    /* update of backward lineage for image par.updateList.data[i]+1
     * from image par.updateList.data[i]
     */
    if ( par.updateList.data[i] < cellProperties.n_data-1 && _isUpdateInList( &(par.updateList), par.updateList.data[i]+1 ) == 0 ) {
      if ( BAL_CellImageBackwardIntersection( &(cellProperties.data[ par.updateList.data[i] ]),
                                         &(cellProperties.data[ par.updateList.data[i]+1 ]) ) != 1 ) {
        BAL_FreeCellSequence( &cellProperties );
        API_FreeParam_cellProperties( &par );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when extracting backward lineage of image #%d\n", proc, par.updateList.data[i]+1 );
        return( -1 );
      }
    }

  }


  /* ...
   */

  if ( 0 ) {
      f = fopen( "sequence_update.txt", "w" );
      BAL_FprintfCellSequence( f, &cellProperties );
      fclose( f );
  }

  /* diagnosis
   */

  if ( output_diagnosis_name != (char*)NULL && output_diagnosis_name[0] != '\0' ) {
    if ( _writeCellSequenceDiagnosis( output_diagnosis_name, &cellProperties, &readProperties ) != 1 ) {
        BAL_FreeCellSequence( &cellProperties );
        API_FreeParam_cellProperties( &par );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when computing diagnosis\n", proc );
        return( -1 );
    }
  }

  /* output writing
   */

  f = (FILE*)NULL;
  BAL_InitPropertyList( &writeProperties );
  BAL_GetWritePropertyList( &writeProperties, &readProperties, &(par.output) );

  if ( par.allow_stdin_stdout != 0 &&
       (output_xml_name == (char*)NULL || output_xml_name[0] == '\0'
        || (output_xml_name[0] == '>' && output_xml_name[1] == '\0')) ) {
    f = stdout;
  }
  else {
    f = fopen( output_xml_name, "w" );
    if ( f == (FILE*)NULL ) {
      BAL_FreeCellSequence( &cellProperties );
      API_FreeParam_cellProperties( &par );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when opening '%s'\n", proc, output_xml_name );
      return( -1 );
    }
  }

  if ( f != (FILE*)NULL ) {
    if ( _verbose_ >= 2 )
      BAL_FprintfPropertyList( stderr, &writeProperties, "* properties to be written" );
    XML_FprintfCellSequenceProperty( f, &cellProperties, &writeProperties );
    if ( par.allow_stdin_stdout != 0 &&
         (output_xml_name == (char*)NULL || output_xml_name[0] == '\0'
          || (output_xml_name[0] == '>' && output_xml_name[1] == '\0')) ) {
      ;
    }
    else {
      fclose( f );
    }
  }

  BAL_FreeCellSequence( &cellProperties );
  API_FreeParam_cellProperties( &par );

  return( 1 );
}





int API_INTERMEDIARY_imageCellProperties( char *theim_name,
                                          char *output_xml_name,
                                          char *param_str_1, char *param_str_2 )
{
  char *proc = "API_INTERMEDIARY_imageCellProperties";
  lineCmdParamCellProperties par;
  bal_image image;
  typeCellImage cellImage;
  FILE *f;

  /* parameter initialization
   */
  API_InitParam_cellProperties( &par );

  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_cellProperties( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_cellProperties( param_str_2, &par );

  /************************************************************
   *
   *  here is the stuff
   *
   ************************************************************/

  if ( BAL_ReadImage( &image, theim_name, 0 ) != 1 ) {
    API_FreeParam_cellProperties( &par );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to read input image '%s'\n", proc, theim_name );
    return( -1 );
  }

  /* output creation
   */
  BAL_InitCellImage( &cellImage );
  cellImage.acquisition_time = par.acquisition_time;

  /* API call
   */

  if ( API_cellImageProperties( &image, &cellImage, param_str_1, param_str_2 ) != 1 ) {
      BAL_FreeCellImage( &cellImage );
      BAL_FreeImage( &image );
      API_FreeParam_cellProperties( &par );
      if ( _verbose_ )
        fprintf( stderr, "%s: some error occurs during processing\n", proc );
      return( -1 );
  }

  BAL_FreeImage( &image );

  /* output writing
   */

  if ( par.allow_stdin_stdout != 0 &&
       (output_xml_name == (char*)NULL || output_xml_name[0] == '\0'
        || (output_xml_name[0] == '>' && output_xml_name[1] == '\0')) ) {
    f = stdout;
  }
  else {
    f = fopen( output_xml_name, "w" );
    if ( f == (FILE*)NULL ) {
      BAL_FreeCellImage( &cellImage );
      API_FreeParam_cellProperties( &par );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when opening '%s'\n", proc, output_xml_name );
      return( -1 );
    }
  }

  XML_FprintfCellImageProperty( f, &cellImage, &(par.output) );

  if ( par.allow_stdin_stdout != 0 &&
       (output_xml_name == (char*)NULL || output_xml_name[0] == '\0'
        || (output_xml_name[0] == '>' && output_xml_name[1] == '\0')) ) {
    ;
  }
  else {
    fclose( f );
  }

  /* memory freeing
   */
  BAL_FreeCellImage( &cellImage );
  API_FreeParam_cellProperties( &par );


  return( 1 );
}





int API_INTERMEDIARY_properties( char *input_xml_name,
                                 char *output_xml_name,
                                 char *output_diagnosis_name,
                                 char *param_str_1, char *param_str_2 )
{
  char *proc = "API_INTERMEDIARY_properties";
  lineCmdParamCellProperties par;
  typeCellSequence cellProperties;
  typePropertyList readProperties;
  typePropertyList writeProperties;
  FILE *f = (FILE*)NULL;

  /* parameter initialization
   */
  API_InitParam_cellProperties( &par );

  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_cellProperties( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_cellProperties( param_str_2, &par );

  /* reading the input xml file
   */
  BAL_InitCellSequence( &cellProperties );
  BAL_InitPropertyList( &readProperties );

  f = fopen( input_xml_name, "r" );
  if ( f == (FILE*)NULL ) {
    BAL_FreeCellSequence( &cellProperties );
    API_FreeParam_cellProperties( &par );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when opening '%s' for reading\n", proc, input_xml_name );
    return( -1 );
  }

  if ( XML_FscanfCellSequenceProperty( f, &cellProperties, &readProperties ) != 1 ) {
      fclose( f );
      BAL_FreeCellSequence( &cellProperties );
      API_FreeParam_cellProperties( &par );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when reading '%s'\n", proc, input_xml_name );
      return( -1 );
  }

  fclose( f );


  /* ...
   */
  if ( 0 ) {
      f = fopen( "sequence_properties.txt", "w" );
      BAL_FprintfCellSequence( f, &cellProperties );
      fclose( f );
  }

  /* diagnosis
   */

  if ( output_diagnosis_name != (char*)NULL && output_diagnosis_name[0] != '\0' ) {
    if ( _writeCellSequenceDiagnosis( output_diagnosis_name, &cellProperties, &readProperties ) != 1 ) {
        BAL_FreeCellSequence( &cellProperties );
        API_FreeParam_cellProperties( &par );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when computing diagnosis\n", proc );
        return( -1 );
    }
  }

  /* output writing
   */

  f = (FILE*)NULL;
  BAL_InitPropertyList( &writeProperties );
  BAL_GetWritePropertyList( &writeProperties, &readProperties, &(par.output) );

  if ( par.allow_stdin_stdout != 0 &&
       (output_xml_name == (char*)NULL || output_xml_name[0] == '\0'
        || (output_xml_name[0] == '>' && output_xml_name[1] == '\0')) ) {
    f = stdout;
  }
  else {
    if ( output_xml_name != (char*)NULL && output_xml_name[0] != '\0' ) {
      f = fopen( output_xml_name, "w" );
      if ( f == (FILE*)NULL ) {
        BAL_FreeCellSequence( &cellProperties );
        API_FreeParam_cellProperties( &par );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when opening '%s' for writing\n", proc, output_xml_name );
        return( -1 );
      }
    }
  }

  if ( f != (FILE*)NULL ) {
    if ( _verbose_ >= 2 )
        BAL_FprintfPropertyList( stderr, &writeProperties, "* properties to be written" );
    XML_FprintfCellSequenceProperty( f, &cellProperties, &writeProperties );
    if ( par.allow_stdin_stdout != 0 &&
         (output_xml_name == (char*)NULL || output_xml_name[0] == '\0'
          || (output_xml_name[0] == '>' && output_xml_name[1] == '\0')) ) {
      ;
    }
    else {
      fclose( f );
    }
  }

  BAL_FreeCellSequence( &cellProperties );
  API_FreeParam_cellProperties( &par );

  return( 1 );
}




/************************************************************
 *
 *
 *
 ************************************************************/





int API_cellImageProperties( bal_image *image, typeCellImage *CellImage, char *param_str_1, char *param_str_2 )
{
  char *proc = "API_cellImageProperties";
  lineCmdParamCellProperties par;



  /* parameter initialization
   */
  API_InitParam_cellProperties( &par );

  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_cellProperties( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_cellProperties( param_str_2, &par );


  /************************************************************
   *
   *  here is the stuff
   *
   ************************************************************/

  if ( BAL_CellImageProperties( image, CellImage,
                                &(par.output), par.surfaceEstimationType ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when processing\n", proc );
    return( -1 );
  }

  return( 1 );
}






/************************************************************
 *
 * static functions
 *
 ************************************************************/



static char **_Str2Array( int *argc, char *str )
{
  char *proc = "_Str2Array";
  int n = 0;
  char *s = str;
  char **array, **a;

  if ( s == (char*)NULL || strlen( s ) == 0 ) {
    if ( _verbose_ >= 2 )
      fprintf( stderr, "%s: empty input string\n", proc );
    *argc = 0;
    return( (char**)NULL );
  }

  /* go to the first valid character
   */
  while ( *s == ' ' || *s == '\n' || *s == '\t' )
    s++;

  if ( *s == '\0' ) {
    if ( _verbose_ >= 2 )
      fprintf( stderr, "%s: weird, input string contains only separation characters\n", proc );
    *argc = 0;
    return( (char**)NULL );
  }

  /* count the number of strings
   */
  for ( n = 0; *s != '\0'; ) {
    n ++;
    while ( *s != ' ' && *s != '\n' && *s != '\t' && *s != '\0' )
      s ++;
    while ( *s == ' ' || *s == '\n' || *s == '\t' )
      s ++;
  }

  if ( _verbose_ >= 5 )
    fprintf( stderr, "%s: found %d strings\n", proc, n );

  /* the value of the strings will be duplicated
   * so that the input string can be freed
   */
  array = (char**)vtmalloc( n * sizeof(char*) + (strlen(str)+1) * sizeof(char), "array", proc );
  if ( array == (char**)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    *argc = 0;
    return( (char**)NULL );
  }

  a = array;
  a += n;
  s = (char*)a;
  (void)strncpy( s, str, strlen( str ) );
  s[ strlen( str ) ] = '\0';

  while ( *s == ' ' || *s == '\n' || *s == '\t' ) {
    *s = '\0';
    s++;
  }

  for ( n = 0; *s != '\0'; ) {
    array[n] = s;
    n ++;
    while ( *s != ' ' && *s != '\n' && *s != '\t' && *s != '\0' )
      s ++;
    while ( *s == ' ' || *s == '\n' || *s == '\t' ) {
      *s = '\0';
      s ++;
    }
  }

  *argc = n;
  return( array );
}





/************************************************************
 *
 * help / documentation
 *
 ************************************************************/



static char *usage = "[image-in]\n\
 [-fusion-format %s]\n\
 [-segmentation-format|-seg-format|-format %s]\n\
 [-f[irst] %d -l[ast] %d]\n\
 [-i|-input|-input-xml %s] [[-o|-output|-output-xml] %s]\n\
 [-diagnosis %s]\n\
 [-a|-acquisition-time %d]\n\
 [-update %d %d %d ... %d]\n\
 [-feature|-property lineage|forward-neighbors|backward-neighbors|\n\
    ...|h_min|volume|surface|sigma|label_in_time|...\n\
    ...|barycenter|principal-value|principal-vector|fate|...\n\
    ...|all-cells|name|history|contact]\n\
 [-s|-surface|-surface-estimation 6n|6-neighbors|windreich|lindblad]\n\
 [-max-chunks-properties %d]\n\
 [-affine-pyramid-lowest-level|-affine-py-ll %d]\n\
 [-affine-pyramid-highest-level|-affine-py-hl %d]\n\
 [-no-non-linear-registration|-no-non-linear]\n\
 [-non-linear-registration|-non-linear]\n\
 [-non-linear-max-iteration[s]|-non-linear-max-iter|-non-linear-iteration[s] %d]\n\
 [-non-linear-pyramid-lowest-level|-non-linear-py-ll %d]\n\
 [-non-linear-pyramid-highest-level|-non-linear-py-hl %d]\n\
 [-parallel|-no-parallel] [-max-chunks %d]\n\
 [-parallelism-type|-parallel-type default|none|openmp|omp|pthread|thread]\n\
 [-omp-scheduling|-omps default|static|dynamic-one|dynamic|guided]\n\
 [-verbose|-v] [-no-verbose|-noverbose|-nv]\n\
 [-debug|-D] [-no-debug|-nodebug]\n\
 [-allow-pipe|-pipe] [-no-allow-pipe|-no-pipe|-nopipe]\n\
 [-print-parameters|-param]\n\
 [-print-time|-time] [-no-time|-notime]\n\
 [-trace-memory|-memory] [-no-memory|-nomemory]\n\
 [-help|-h]";



static char *detail = "\
 # file names \n\
 image-in: single input image\n\
 -fusion-format %s : \n\
    format 'a la printf' of fusion images (used to compute co-registrations)\n\
    must contain one '%d'\n\
 -segmentation-format|-seg-format|-format %s : \n\
    format 'a la printf' of segmentation images to be processed\n\
    must contain one '%d'\n\
 -f[irst] %d            : first value of the index in the format\n\
 -l[ast] %d             : last value of the index in the format\n\
 -input|-input-xml %s   : input xml file\n\
 [-output|-output-xml] file-out\n\
 -diagnosis %s: summary of some diagnosis\n\
# image information \n\
 -a|-acquisition-time %d: image label or acquisition time\n\
    (used for single image processing), used to create cell label\n\
 -update %d %d %d ... %d: indexes of images to be updated\n\
     properties of the cells belonging to these images will be updated\n\
     as well as lineage of previous images\n\
# properties\n\
 -feature|-property lineage|forward-neighbors|backward-neighbors|\n\
     ...|h_min|volume|surface|sigma|label_in_time|...\n\
     ...|barycenter|principal-value|principal-vector|fate|...\n\
     ...|all-cells|name|history|contact:\n\
     features to be written in the output file (most of them are not implemented yet)\n\
 -s|-surface|-surface-estimation %s : contact surface estimation method\n\
    6n|6-neighbors: outer 6-neighbors of the cell, Astec historical method\n\
    windreich:\n\
       Voxel-based surface area estimation: from theory to practice\n\
       G. Windreich, N. Kiryati, and G. Lohmann\n\
       Pattern Recognition, Volume 36, Issue 11, November 2003, Pages 2531-2541\n\
    lindblad: (default)\n\
       Surface area estimation of digitized 3D objects using weighted local configurations,\n\
       Joakim Lindblad, Image and Vision Computing, 23(2):111-122, 2005\n\
# properties-related parallelism parameter\n\
  -max-chunks-properties %d:\n\
    single-image based properties are computed with a parallelism on data,\n\
    i.e. several images are processed in parallel. Since this may cause disk\n\
    access troubles, this parameter allows to control the number of files\n\
    opened simulteanously\n\
# linear registration parameters\n\
 -affine-pyramid-lowest-level|-affine-py-ll %d: \n\
    pyramid lowest level (0 = original dimension, default is 3)\n\
 -affine-pyramid-highest-level|-affine-py-hl %\n\
    pyramid highest level (default is 5)\n\
# non-linear registration parameters\n\
 -no-non-linear-registration|-no-non-linear :\n\
    do not perform non-linear registration (decrease computational cost)\
 -non-linear-registration|-non-linear: \n\
    perform non-linear registration (default)\n\
 -non-linear-max-iteration[s]|-non-linear-max-iter|-non-linear-iteration[s] %d: \n\
    maximum number of iterations ate each level (default is 10)\n\
 -non-linear-pyramid-lowest-level|-non-linear-py-ll %d: \n\
    pyramid lowest level (0 = original dimension, default is 3)\n\
 -non-linear-pyramid-highest-level|-non-linear-py-hl %d: \n\
    pyramid highest level (default is 5)\n\
# parallelism parameters\n\
 -parallel|-no-parallel:\n\
 -max-chunks %d:\n\
 -parallelism-type|-parallel-type default|none|openmp|omp|pthread|thread:\n\
 -omp-scheduling|-omps default|static|dynamic-one|dynamic|guided:\n\
# general parameters \n\
  -verbose|-v: increase verboseness\n\
    parameters being read several time, use '-nv -v -v ...'\n\
    to set the verboseness level\n\
  -no-verbose|-noverbose|-nv: no verboseness at all\n\
  -debug|-D: increase debug level\n\
  -no-debug|-nodebug: no debug indication\n\
  -allow-pipe|-pipe: allow the use of stdin/stdout (with '-')\n\
  -no-allow-pipe|-no-pipe|-nopipe: do not allow the use of stdin/stdout\n\
  -print-parameters|-param:\n\
  -print-time|-time:\n\
  -no-time|-notime:\n\
  -trace-memory|-memory:\n\
  -no-memory|-nomemory:\n\
  -h: print option list\n\
  -help: print option list + details\n\
";





char *API_Help_cellProperties( int h )
{
    if ( h == 0 )
        return( usage );
    return( detail );
}





void API_ErrorParse_cellProperties( char *program, char *str, int flag )
{
    if ( flag >= 0 ) {
        if ( program != (char*)NULL )
           (void)fprintf(stderr,"Usage: %s %s\n", program, usage);
        else
            (void)fprintf(stderr,"Command line options: %s\n", usage);
    }
    if ( flag == 1 ) {
      (void)fprintf( stderr, "--------------------------------------------------\n" );
      (void)fprintf(stderr,"%s",detail);
      (void)fprintf( stderr, "--------------------------------------------------\n" );
    }
    if ( str != (char*)NULL )
      (void)fprintf(stderr,"Error: %s\n",str);
    exit( 1 );
}





/************************************************************
 *
 * parameters management
 *
 ************************************************************/



void API_InitParam_cellProperties( lineCmdParamCellProperties *p )
{
    (void)strncpy( p->input_segmentation_name, "\0", 1 );

    (void)strncpy( p->input_fusion_format, "\0", 1 );
    (void)strncpy( p->input_segmentation_format, "\0", 1 );
    p->firstindex = -1;
    p->lastindex = -1;

    (void)strncpy( p->input_xml_name, "\0", 1 );
    (void)strncpy( p->output_xml_name, "\0", 1 );

    (void)strncpy( p->output_diagnosis_name, "\0", 1 );

    p->acquisition_time = 0;
    BAL_InitUpdateList( &(p->updateList) );

    BAL_InitPropertyList( &(p->output) );
    p->surfaceEstimationType = _LINDBLAD_;

    p->chunks_for_properties = -1;

    p->normalisation = 1;

    BAL_InitBlockMatchingPyramidalParameters( &(p->affineRegistration) );
    p->affineRegistration.transformation_type = AFFINE_3D;
    BAL_InitBlockMatchingPyramidalParametersForLinearTransformation( &(p->affineRegistration) );
    p->affineRegistration.pyramid_lowest_level = 3;
    p->affineRegistration.pyramid_highest_level = 5;
    p->affineRegistration.estimator.lowest.type = TYPE_WLTS;
    p->affineRegistration.estimator.lowest.retained_fraction = 0.55;

    BAL_InitBlockMatchingPyramidalParameters( &(p->nonlinearRegistration) );
    p->nonlinearRegistration.transformation_type = VECTORFIELD_3D;
    BAL_InitBlockMatchingPyramidalParametersForVectorfieldTransformation( &(p->nonlinearRegistration) );
    p->nonlinearRegistration.pyramid_lowest_level = 3;
    p->nonlinearRegistration.pyramid_highest_level = 5;
    p->nonlinearRegistration.estimator.lowest.type = TYPE_WLTS;
    p->nonlinearRegistration.estimator.lowest.retained_fraction = 0.55;
    p->nonlinearRegistration.estimator.lowest.sigma.x = 2.0;
    p->nonlinearRegistration.estimator.lowest.sigma.y = 2.0;
    p->nonlinearRegistration.estimator.lowest.sigma.z = 2.0;
    p->nonlinearRegistration.estimator.highest.sigma = p->nonlinearRegistration.estimator.lowest.sigma;
    p->nonlinearRegistration.elastic_regularization_sigma.lowest.x = 2.0;
    p->nonlinearRegistration.elastic_regularization_sigma.lowest.y = 2.0;
    p->nonlinearRegistration.elastic_regularization_sigma.lowest.z = 2.0;
    p->nonlinearRegistration.elastic_regularization_sigma.highest = p->nonlinearRegistration.elastic_regularization_sigma.lowest;
    p->nonlinearRegistration.pyramid_gaussian_filtering = 1;

    p->allow_stdin_stdout = 0;
    p->print_lineCmdParam = 0;
    p->print_time = 0;
    p->trace_allocations = 0;
}



void API_FreeParam_cellProperties( lineCmdParamCellProperties *p )
{
  BAL_FreeUpdateList( &(p->updateList) );
  API_InitParam_cellProperties( p );
}



void API_PrintParam_cellProperties( FILE *theFile, char *program,
                                  lineCmdParamCellProperties *p, char *str )
{
  FILE *f = theFile;
  if ( theFile == (FILE*)NULL ) f = stderr;

  fprintf( f, "==================================================\n" );
  fprintf( f, "= in line command parameters" );
  if ( program != (char*)NULL )
    fprintf( f, " for '%s'", program );
  if ( str != (char*)NULL )
    fprintf( f, "= %s\n", str );
  fprintf( f, "\n"  );
  fprintf( f, "==================================================\n" );


  fprintf( f, "# image names\n" );

  fprintf( f, "- input_segmentation_name is " );
  if ( p->input_segmentation_name != (char*)NULL && p->input_segmentation_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_segmentation_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- input_fusion_format is " );
  if ( p->input_fusion_format != (char*)NULL && p->input_fusion_format[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_fusion_format );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- input_segmentation_format is " );
  if ( p->input_segmentation_format != (char*)NULL && p->input_segmentation_format[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_segmentation_format );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- firstindex  = %d\n", p->firstindex );
  fprintf( f, "- lastindex  = %d\n", p->lastindex );

  fprintf( f, "- input_xml_name is " );
  if ( p->input_xml_name != (char*)NULL && p->input_xml_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_xml_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- output_xml_name is " );
  if ( p->output_xml_name != (char*)NULL && p->output_xml_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_xml_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- output_diagnosis_name is " );
  if ( p->output_diagnosis_name != (char*)NULL && p->output_diagnosis_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_diagnosis_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- acquisition_time  = %d\n", p->acquisition_time );
  fprintf( f, "- updateList = " );
  BAL_FprintfUpdateList( f, &(p->updateList) );

  BAL_FprintfPropertyList( f, &(p->output), "- output." );
  fprintf( f, "- surfaceEstimationType = " );
  BAL_FprintfEnumSurfaceEstimation( f, p->surfaceEstimationType );
  fprintf( f, "\n" );

  fprintf( f, "# dedicated parallelism parameters\n" );
  fprintf( f, "- p->chunks_for_properties = %d\n", p->chunks_for_properties );

  fprintf( f, "# registration parameters\n" );
  fprintf( f, "- p->normalisation = %d\n", p->normalisation );
  fprintf( f, "# parameters for affine registration\n" );
  BAL_PrintBlockMatchingPyramidalParameters( f, &(p->affineRegistration) );
  if ( 0 ) {
  fprintf( f, "# parameters for non-linear registration\n" );
  BAL_PrintBlockMatchingPyramidalParameters( f, &(p->nonlinearRegistration) );
  }

  fprintf( f, "# general parameters\n" );
  fprintf( f, "- allows stdin/stdout  = %d\n", p->allow_stdin_stdout );
  fprintf( f, "- print parameters     = %d\n", p->print_lineCmdParam );
  fprintf( f, "- print time           = %d\n", p->print_time );
  fprintf( f, "- p->trace_allocations = %d\n", p->trace_allocations );

  fprintf( f, "==================================================\n" );
}





/************************************************************
 *
 * parameters parsing
 *
 ************************************************************/



static void _API_ParseParam_cellProperties( char *str, lineCmdParamCellProperties *p )
{
  char *proc = "_API_ParseParam_cellProperties";
  char **argv;
  int i, argc;

  if ( str == (char*)NULL || strlen(str) == 0 )
      return;

  argv = _Str2Array( &argc, str );
  if ( argv == (char**)NULL || argc == 0 ) {
      if ( _debug_ ) {
          fprintf( stderr, "%s: weird, no arguments were found\n", proc );
      }
      return;
  }

  if ( _debug_ > 4 ) {
      fprintf( stderr, "%s: translation from\n", proc );
      fprintf( stderr, "   '%s'\n", str );
      fprintf( stderr, "into\n" );
      for ( i=0; i<argc; i++ )
          fprintf( stderr, "   argv[%2d] = '%s'\n", i, argv[i] );
  }

  API_ParseParam_cellProperties( 0, argc, argv, p );

  vtfree( argv );
}





static int _n_call_parse_ = 0;

void API_ParseParam_cellProperties( int firstargc, int argc, char *argv[],
                                  lineCmdParamCellProperties *p )
{
  int i;
  int inputisread = 0;
  int outputisread = 0;
  char text[STRINGLENGTH];
  int updateValue, status;
  int maxchunks;

  _n_call_parse_ ++;

  /* option line parsing
   */
  for ( i=firstargc; i<argc; i++ ) {

      /* strings beginning with '-'
       */
      if ( argv[i][0] == '-' ) {
          if ( argv[i][1] == '\0' ) {
            if ( inputisread == 0 ) {
              if ( p->allow_stdin_stdout ) {
                (void)strcpy( p->input_segmentation_name,  "<" );  /* standard input */
                inputisread = 1;
              }
              else {
                API_ErrorParse_cellProperties( (char*)NULL, "parsing '-', stdin is not allowed\n", 0 );
              }
            }
            else if ( outputisread == 0 ) {
              if ( p->allow_stdin_stdout ) {
                (void)strcpy( p->output_xml_name,  ">" );  /* standard output */
                outputisread = 1;
              }
              else {
                API_ErrorParse_cellProperties( (char*)NULL, "parsing '-', stdout is not allowed\n", 0 );
              }
            }
            else {
              API_ErrorParse_cellProperties( (char*)NULL, "too many file names, parsing '-' ...\n", 0 );
            }
          }

          /* file names
           */

          else if ( strcmp ( argv[i], "-fusion-format" ) == 0 ) {
            i++;
            if ( i >= argc) API_ErrorParse_cellProperties( (char*)NULL, "parsing -fusion-format...\n", 0 );
            (void)strcpy( p->input_fusion_format, argv[i] );
          }

          else if ( strcmp ( argv[i], "-segmentation-format" ) == 0
                    || strcmp ( argv[i], "-seg-format" ) == 0
                    || strcmp ( argv[i], "-format" ) == 0 ) {
            i++;
            if ( i >= argc) API_ErrorParse_cellProperties( (char*)NULL, "parsing -segmentation-format...\n", 0 );
            (void)strcpy( p->input_segmentation_format, argv[i] );
            inputisread = 1;
          }


          else if ( (strcmp ( argv[i], "-f" ) == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "-first" ) == 0 && argv[i][6] == '\0') ) {
            i++;
            if ( i >= argc) API_ErrorParse_cellProperties( (char*)NULL, "parsing -first ...\n", 0 );
            status = sscanf( argv[i], "%d", &(p->firstindex) );
            if ( status <= 0 ) API_ErrorParse_cellProperties( (char*)NULL, "parsing -first ...", 0 );
          }
          else if ( (strcmp ( argv[i], "-l" ) == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "-last" ) == 0 && argv[i][5] == '\0') ) {
            i++;
            if ( i >= argc) API_ErrorParse_cellProperties( (char*)NULL, "parsing -last ...\n", 0 );
            status = sscanf( argv[i], "%d", &(p->lastindex) );
            if ( status <= 0 ) API_ErrorParse_cellProperties( (char*)NULL,"parsing -last ...", 0 );
          }

          else if ( (strcmp ( argv[i], "-i" ) == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "-input" ) == 0 && argv[i][6] == '\0')
                    || strcmp ( argv[i], "-input-xml" ) == 0 ) {
            i++;
            if ( i >= argc) API_ErrorParse_cellProperties( (char*)NULL, "parsing -input-xml ...\n", 0 );
            (void)strcpy( p->input_xml_name, argv[i] );
            outputisread = 1;
          }

          else if ( (strcmp ( argv[i], "-o" ) == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "-output" ) == 0 && argv[i][7] == '\0')
                    || strcmp ( argv[i], "-output-xml" ) == 0 ) {
            i++;
            if ( i >= argc) API_ErrorParse_cellProperties( (char*)NULL, "parsing -output-xml ...\n", 0 );
            (void)strcpy( p->output_xml_name, argv[i] );
            outputisread = 1;
          }

          else if ( strcmp ( argv[i], "-diagnosis" ) == 0 ) {
            i++;
            if ( i >= argc) API_ErrorParse_cellProperties( (char*)NULL, "parsing -diagnosis ...\n", 0 );
            (void)strcpy( p->output_diagnosis_name, argv[i] );
          }

          /* informations
           */

          else if ( (strcmp ( argv[i], "-a" ) == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "-acquisition-time" ) == 0) ) {
              i++;
              if ( i >= argc) API_ErrorParse_cellProperties( (char*)NULL, "parsing -acquisition-time ...\n", 0 );
              status = sscanf( argv[i], "%d", &(p->acquisition_time) );
              if ( status <= 0 ) API_ErrorParse_cellProperties( (char*)NULL,"parsing -acquisition-time ...", 0 );
          }

          else if ( (strcmp ( argv[i], "-u" ) == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "-update" ) == 0 && argv[i][7] == '\0') ) {
            i++;
            if ( i >= argc) API_ErrorParse_cellProperties( (char*)NULL, "parsing -update ...\n", 0 );
            do {
              status = sscanf( argv[i], "%d", &updateValue );
              if ( status == 1 ) {
                if ( BAL_AddUpdateToList( &(p->updateList), updateValue ) != 1 ) {
                  API_ErrorParse_cellProperties( (char*)NULL, "adding index to be updated ...\n", 0 );
                }
                i ++;
              }
            } while ( status == 1 && i < argc );
            if ( status <= 0 ) i--;
          }

          /* properties
           */

          else if ( strcmp ( argv[i], "-feature" ) == 0
                    || strcmp ( argv[i], "-property" ) == 0 ) {
              i++;
              if ( i >= argc) API_ErrorParse_cellProperties( (char*)NULL, "parsing -feature ...\n", 0 );
              while( 1 ) {
                if ( strcmp ( argv[i], "lineage" ) == 0  ) {
                  if ( p->output.n_data == p->output.n_allocated_data )
                      API_ErrorParse_cellProperties( (char*)NULL, "-feature, too many features...\n", 0 );
                  p->output.data[ p->output.n_data++ ] = _LINEAGE_;
                }
                else if ( strcmp ( argv[i], "forward-neighbors" ) == 0  ) {
                  if ( p->output.n_data == p->output.n_allocated_data )
                      API_ErrorParse_cellProperties( (char*)NULL, "-feature, too many features...\n", 0 );
                  p->output.data[ p->output.n_data++ ] = _FORWARD_NEIGHBORS_;
                }
                else if ( strcmp ( argv[i], "backward-neighbors" ) == 0  ) {
                  if ( p->output.n_data == p->output.n_allocated_data )
                      API_ErrorParse_cellProperties( (char*)NULL, "-feature, too many features...\n", 0 );
                  p->output.data[ p->output.n_data++ ] = _BACKWARD_NEIGHBORS_;
                }
                else if ( strcmp ( argv[i], "h_min" ) == 0  ) {
                  if ( p->output.n_data == p->output.n_allocated_data )
                    API_ErrorParse_cellProperties( (char*)NULL, "-feature, too many features...\n", 0 );
                  p->output.data[ p->output.n_data++ ] = _H_MIN_;
                }
                else if ( strcmp ( argv[i], "volume" ) == 0  ) {
                  if ( p->output.n_data == p->output.n_allocated_data )
                    API_ErrorParse_cellProperties( (char*)NULL, "-feature, too many features...\n", 0 );
                  p->output.data[ p->output.n_data++ ] = _VOLUME_;
                }
                else if ( strcmp ( argv[i], "surface" ) == 0  ) {
                  if ( p->output.n_data == p->output.n_allocated_data )
                    API_ErrorParse_cellProperties( (char*)NULL, "-feature, too many features...\n", 0 );
                  p->output.data[ p->output.n_data++ ] = _SURFACE_;
                }
                else if ( strcmp ( argv[i], "sigma" ) == 0  ) {
                  if ( p->output.n_data == p->output.n_allocated_data )
                    API_ErrorParse_cellProperties( (char*)NULL, "-feature, too many features...\n", 0 );
                  p->output.data[ p->output.n_data++ ] = _SIGMA_;
                }
                else if ( strcmp ( argv[i], "label_in_time" ) == 0  ) {
                  if ( p->output.n_data == p->output.n_allocated_data )
                    API_ErrorParse_cellProperties( (char*)NULL, "-feature, too many features...\n", 0 );
                  p->output.data[ p->output.n_data++ ] = _LABEL_IN_TIME_;
                }
                else if ( strcmp ( argv[i], "barycenter" ) == 0  ) {
                  if ( p->output.n_data == p->output.n_allocated_data )
                    API_ErrorParse_cellProperties( (char*)NULL, "-feature, too many features...\n", 0 );
                  p->output.data[ p->output.n_data++ ] = _BARYCENTER_;
                }
                else if ( strcmp ( argv[i], "principal-value" ) == 0  ) {
                  if ( p->output.n_data == p->output.n_allocated_data )
                    API_ErrorParse_cellProperties( (char*)NULL, "-feature, too many features...\n", 0 );
                  p->output.data[ p->output.n_data++ ] = _PRINCIPAL_VALUE_;
                }
                else if ( strcmp ( argv[i], "principal-vector" ) == 0  ) {
                  if ( p->output.n_data == p->output.n_allocated_data )
                    API_ErrorParse_cellProperties( (char*)NULL, "-feature, too many features...\n", 0 );
                  p->output.data[ p->output.n_data++ ] = _PRINCIPAL_VECTOR_;
                }
                else if ( strcmp ( argv[i], "fate" ) == 0  ) {
                  if ( p->output.n_data == p->output.n_allocated_data )
                    API_ErrorParse_cellProperties( (char*)NULL, "-feature, too many features...\n", 0 );
                  p->output.data[ p->output.n_data++ ] = _FATE_;
                }
                else if ( strcmp ( argv[i], "all-cells" ) == 0  ) {
                  if ( p->output.n_data == p->output.n_allocated_data )
                    API_ErrorParse_cellProperties( (char*)NULL, "-feature, too many features...\n", 0 );
                  p->output.data[ p->output.n_data++ ] = _ALL_CELLS_;
                }
                else if ( strcmp ( argv[i], "name" ) == 0  ) {
                  if ( p->output.n_data == p->output.n_allocated_data )
                    API_ErrorParse_cellProperties( (char*)NULL, "-feature, too many features...\n", 0 );
                  p->output.data[ p->output.n_data++ ] = _NAME_;
                }
                else if ( strcmp ( argv[i], "history" ) == 0  ) {
                  if ( p->output.n_data == p->output.n_allocated_data )
                    API_ErrorParse_cellProperties( (char*)NULL, "-feature, too many features...\n", 0 );
                  p->output.data[ p->output.n_data++ ] = _HISTORY_;
                }
                else if ( strcmp ( argv[i], "contact" ) == 0
                          || strcmp ( argv[i], "contact-surface" ) == 0 ) {
                  if ( p->output.n_data == p->output.n_allocated_data )
                    API_ErrorParse_cellProperties( (char*)NULL, "-feature, too many features...\n", 0 );
                  p->output.data[ p->output.n_data++ ] = _CONTACT_SURFACE_;
                }
                else {
                  i--;
                  break;
                }
                if ( i == argc - 1 ) break;
                i++;
              }
          }

          else if ( (strcmp ( argv[i], "-s" ) == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "-surface" ) == 0)
                    || (strcmp ( argv[i], "-surface-estimation" ) == 0) ) {
              i++;
              if ( i >= argc) API_ErrorParse_cellProperties( (char*)NULL, "parsing -surface-estimation ...\n", 0 );
              if ( (strcmp ( argv[i], "6n" ) == 0 && argv[i][2] == '\0')
                   || strcmp ( argv[i], "6-neighbors" ) == 0  ) {
                  p->surfaceEstimationType = _OUTER_6NEIGHBORS_;
              }
              else if ( strcmp ( argv[i], "lindblad" ) == 0 ) {
                  p->surfaceEstimationType = _LINDBLAD_;
              }
              else if ( strcmp ( argv[i], "windreich" ) == 0 ) {
                  p->surfaceEstimationType = _WINDREICH_;
              }
              else {
                  fprintf( stderr, "unknown surface estimation type: '%s'\n", argv[i] );
                  API_ErrorParse_cellProperties( (char*)NULL,"parsing -surface-estimation ...", 0 );
              }
          }
          /* dedicated parallelism parameters
           */
          else if ( strcmp ( argv[i], "-max-chunks-properties" ) == 0 ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_cellProperties( (char*)NULL, "parsing -max-chunks-properties ...\n", 0 );
             status = sscanf( argv[i], "%d", &(p->chunks_for_properties) );
             if ( status <= 0 ) API_ErrorParse_cellProperties( (char*)NULL, "parsing -max-chunks-properties ...\n", 0 );
          }

          /* linear registration parameters
           */
          else if ( strcmp ( argv[i], "-affine-pyramid-lowest-level" ) == 0
                    || strcmp ( argv[i], "-affine-py-ll" ) == 0 ) {
              i++;
              if ( i >= argc) API_ErrorParse_cellProperties( (char*)NULL, "parsing -affine-pyramid-lowest-level ...\n", 0 );
              status = sscanf( argv[i], "%d", &(p->affineRegistration.pyramid_lowest_level) );
              if ( status <= 0 ) API_ErrorParse_cellProperties( (char*)NULL,"parsing -affine-pyramid-lowest-level ...", 0 );
          }
          else if ( strcmp ( argv[i], "-affine-pyramid-highest-level" ) == 0
                    || strcmp ( argv[i], "-affine-py-hl" ) == 0 ) {
              i++;
              if ( i >= argc) API_ErrorParse_cellProperties( (char*)NULL, "parsing -affine-pyramid-highest-level ...\n", 0 );
              status = sscanf( argv[i], "%d", &(p->affineRegistration.pyramid_highest_level) );
              if ( status <= 0 ) API_ErrorParse_cellProperties( (char*)NULL,"parsing -affine-pyramid-highest-level ...", 0 );
          }

          /* non-linear registration parameters
           */
          else if ( strcmp ( argv[i], "-no-non-linear-registration" ) == 0
                    || strcmp ( argv[i], "-no-non-linear" ) == 0 ) {
              p->nonlinearRegistration.max_iterations.lowest = p->nonlinearRegistration.max_iterations.highest = 0;
          }
          else if ( strcmp ( argv[i], "-non-linear-registration" ) == 0
                    || strcmp ( argv[i], "-non-linear" ) == 0 ) {
              p->nonlinearRegistration.max_iterations.lowest = p->nonlinearRegistration.max_iterations.highest = 10;
          }
          else if ( strcmp ( argv[i], "-non-linear-max-iterations" ) == 0
                    || strcmp ( argv[i], "-non-linear-max-iteration" ) == 0
                    || strcmp ( argv[i], "-non-linear-max-iter" ) == 0
                    || strcmp ( argv[i], "-non-linear-iterations" ) == 0
                    || strcmp ( argv[i], "-non-linear-iteration" ) == 0 ) {
              i++;
              if ( i >= argc) API_ErrorParse_cellProperties( (char*)NULL, "parsing -non-linear-max-iterations ...\n", 0 );
              status = sscanf( argv[i], "%d", &(p->nonlinearRegistration.max_iterations.highest) );
              if ( status <= 0 ) API_ErrorParse_cellProperties( (char*)NULL,"parsing -non-linear-max-iterations ...", 0 );
              p->nonlinearRegistration.max_iterations.lowest = p->nonlinearRegistration.max_iterations.highest;
          }

          else if ( strcmp ( argv[i], "-non-linear-pyramid-lowest-level" ) == 0
                    || strcmp ( argv[i], "-non-linear-py-ll" ) == 0 ) {
              i++;
              if ( i >= argc) API_ErrorParse_cellProperties( (char*)NULL, "parsing -non-linear-pyramid-lowest-level ...\n", 0 );
              status = sscanf( argv[i], "%d", &(p->nonlinearRegistration.pyramid_lowest_level) );
              if ( status <= 0 ) API_ErrorParse_cellProperties( (char*)NULL,"parsing -non-linear-pyramid-lowest-level ...", 0 );
          }
          else if ( strcmp ( argv[i], "-non-linear-pyramid-highest-level" ) == 0
                    || strcmp ( argv[i], "-non-linear-py-hl" ) == 0 ) {
              i++;
              if ( i >= argc) API_ErrorParse_cellProperties( (char*)NULL, "parsing -non-linear-pyramid-highest-level ...\n", 0 );
              status = sscanf( argv[i], "%d", &(p->nonlinearRegistration.pyramid_highest_level) );
              if ( status <= 0 ) API_ErrorParse_cellProperties( (char*)NULL,"parsing -non-linear-pyramid-highest-level ...", 0 );
          }

          /* parallelism parameters
           */
          else if ( strcmp ( argv[i], "-parallel" ) == 0 ) {
             setParallelism( _DEFAULT_PARALLELISM_ );
          }

          else if ( strcmp ( argv[i], "-no-parallel" ) == 0 ) {
             setParallelism( _NO_PARALLELISM_ );
          }

          else if ( strcmp ( argv[i], "-parallelism-type" ) == 0 ||
                      strcmp ( argv[i], "-parallel-type" ) == 0 ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_cellProperties( (char*)NULL, "parsing -parallelism-type ...\n", 0 );
             if ( strcmp ( argv[i], "default" ) == 0 ) {
               setParallelism( _DEFAULT_PARALLELISM_ );
             }
             else if ( strcmp ( argv[i], "none" ) == 0 ) {
               setParallelism( _NO_PARALLELISM_ );
             }
             else if ( strcmp ( argv[i], "openmp" ) == 0 || strcmp ( argv[i], "omp" ) == 0 ) {
               setParallelism( _OMP_PARALLELISM_ );
             }
             else if ( strcmp ( argv[i], "pthread" ) == 0 || strcmp ( argv[i], "thread" ) == 0 ) {
               setParallelism( _PTHREAD_PARALLELISM_ );
             }
             else {
               fprintf( stderr, "unknown parallelism type: '%s'\n", argv[i] );
               API_ErrorParse_cellProperties( (char*)NULL, "parsing -parallelism-type ...\n", 0 );
             }
          }

          else if ( strcmp ( argv[i], "-max-chunks" ) == 0 && argv[i][11] == '\0' ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_cellProperties( (char*)NULL, "parsing -max-chunks ...\n", 0 );
             status = sscanf( argv[i], "%d", &maxchunks );
             if ( status <= 0 ) API_ErrorParse_cellProperties( (char*)NULL, "parsing -max-chunks ...\n", 0 );
             if ( maxchunks >= 1 ) setMaxChunks( maxchunks );
          }

          else if ( strcmp ( argv[i], "-omp-scheduling" ) == 0 ||
                   ( strcmp ( argv[i], "-omps" ) == 0 && argv[i][5] == '\0') ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_cellProperties( (char*)NULL, "parsing -omp-scheduling, no argument\n", 0 );
             if ( strcmp ( argv[i], "default" ) == 0 ) {
               setOmpScheduling( _DEFAULT_OMP_SCHEDULING_ );
             }
             else if ( strcmp ( argv[i], "static" ) == 0 ) {
               setOmpScheduling( _STATIC_OMP_SCHEDULING_ );
             }
             else if ( strcmp ( argv[i], "dynamic-one" ) == 0 ) {
               setOmpScheduling( _DYNAMIC_ONE_OMP_SCHEDULING_ );
             }
             else if ( strcmp ( argv[i], "dynamic" ) == 0 ) {
               setOmpScheduling( _DYNAMIC_OMP_SCHEDULING_ );
             }
             else if ( strcmp ( argv[i], "guided" ) == 0 ) {
               setOmpScheduling( _GUIDED_OMP_SCHEDULING_ );
             }
             else {
               fprintf( stderr, "unknown omp scheduling type: '%s'\n", argv[i] );
               API_ErrorParse_cellProperties( (char*)NULL, "parsing -omp-scheduling ...\n", 0 );
             }
          }

          /* general parameters
           */
          else if ( (strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0')
                    || (strcmp ( argv[i], "--help" ) == 0 && argv[i][6] == '\0') ) {
             API_ErrorParse_cellProperties( (char*)NULL, (char*)NULL, 1);
          }
          else if ( (strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0') ) {
             API_ErrorParse_cellProperties( (char*)NULL, (char*)NULL, 0);
          }
          else if ( strcmp ( argv[i], "-verbose" ) == 0
                    || (strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0') ) {
            if ( _n_call_parse_ == 1 ) {
              if ( _verbose_ <= 0 ) _verbose_ = 1;
              else                  _verbose_ ++;
              BAL_IncrementVerboseInBalCellProperties( );
            }
          }
          else if ( strcmp ( argv[i], "-no-verbose" ) == 0
                    || strcmp ( argv[i], "-noverbose" ) == 0
                    || (strcmp ( argv[i], "-nv" ) == 0 && argv[i][3] == '\0') ) {
              _verbose_ = 0;
              BAL_SetVerboseInBalCellProperties( 0 );
          }
          else if ( (strcmp ( argv[i], "-debug" ) == 0 && argv[i][6] == '\0')
                    || (strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0') ) {
            if ( _n_call_parse_ == 1 ) {
              if ( _debug_ <= 0 ) _debug_ = 1;
              else                _debug_ ++;
            }
            incrementDebugInChunks( );
          }
          else if ( (strcmp ( argv[i], "-no-debug" ) == 0 && argv[i][9] == '\0')
                    || (strcmp ( argv[i], "-nodebug" ) == 0 && argv[i][8] == '\0') ) {
              _debug_ = 0;
              setDebugInChunks( 0 );
          }

          else if ( strcmp ( argv[i], "-allow-pipe" ) == 0
                    || (strcmp ( argv[i], "-pipe" ) == 0 && argv[i][5] == '\0') ) {
            p->allow_stdin_stdout = 1;
          }

          else if ( strcmp ( argv[i], "-no-allow-pipe" ) == 0
                    || (strcmp ( argv[i], "-no-pipe" ) == 0 && argv[i][8] == '\0')
                    || (strcmp ( argv[i], "-nopipe" ) == 0 && argv[i][7] == '\0') ) {
            p->allow_stdin_stdout = 0;
          }

          else if ( strcmp ( argv[i], "-print-parameters" ) == 0
                    || (strcmp ( argv[i], "-param" ) == 0 && argv[i][6] == '\0') ) {
             p->print_lineCmdParam = 1;
          }

          else if ( strcmp ( argv[i], "-print-time" ) == 0
                     || (strcmp ( argv[i], "-time" ) == 0 && argv[i][5] == '\0') ) {
             p->print_time = 1;
          }
          else if ( (strcmp ( argv[i], "-notime" ) == 0 && argv[i][7] == '\0')
                      || (strcmp ( argv[i], "-no-time" ) == 0 && argv[i][8] == '\0') ) {
             p->print_time = 0;
          }

          else if ( strcmp ( argv[i], "-trace-memory" ) == 0
                     || (strcmp ( argv[i], "-memory" ) == 0 && argv[i][7] == '\0') ) {
             if ( _n_call_parse_ == 1 ) {
               incrementTraceInVtMalloc( );
               if ( p->trace_allocations  <= 0 ) p->trace_allocations  = 1;
               else                              p->trace_allocations  ++;
             }
             if ( 0 ) setParallelism( _NO_PARALLELISM_ );
          }
          else if ( (strcmp ( argv[i], "-nomemory" ) == 0 && argv[i][9] == '\0')
                      || (strcmp ( argv[i], "-no-memory" ) == 0 && argv[i][10] == '\0') ) {
             setTraceInVtMalloc( 0 );
          }

          /* unknown option
           */
          else {
              sprintf(text,"unknown option %s\n",argv[i]);
              API_ErrorParse_cellProperties( (char*)NULL, text, 0);
          }
      }

      /* strings beginning with a character different from '-'
       */
      else {
          if ( strlen( argv[i] ) >= STRINGLENGTH ) {
              fprintf( stderr, "... parsing '%s'\n", argv[i] );
              API_ErrorParse_cellProperties( (char*)NULL, "too long file name ...\n", 0 );
          }
          else if ( inputisread == 0 ) {
              (void)strcpy( p->input_segmentation_name, argv[i] );
              inputisread = 1;
          }
          else if ( outputisread == 0 ) {
              (void)strcpy( p->output_xml_name, argv[i] );
              outputisread = 1;
          }
          else {
              fprintf( stderr, "... parsing '%s'\n", argv[i] );
              API_ErrorParse_cellProperties( (char*)NULL, "too many file names ...\n", 0 );
          }
      }
  }

  /* if not enough file names
   */
  if ( inputisread == 0 ) {
    if ( p->allow_stdin_stdout ) {
      (void)strcpy( p->input_segmentation_name,  "<" );  /* standard input */
      inputisread = 1;
    }
  }
  if ( outputisread == 0 ) {
    if ( p->allow_stdin_stdout ) {
      (void)strcpy( p->output_xml_name,  ">" );  /* standard output */
      outputisread = 1;
    }
  }

}
