/*************************************************************************
 * tree.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2016, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mar 22 mar 2016 12:02:49 CET
 *
 * ADDITIONS, CHANGES
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include <vtmalloc.h>

#include <api-tree.h>

#include <vt_common.h>





static int _verbose_ = 1;





/* static function definitions
 */

static int _writeOutputImage( typeVoxelTree *iTree,
                              char *name, bufferType type, vt_image *image,
                              enumVoxelTreeAttribute attribute );

static char *_Array2Str( int argc, char *argv[] );
static char *_BaseName( char *p );
static double _GetTime();
static double _GetClock();






int main( int argc, char *argv[] )
{
  lineCmdParamTree par;
  vt_image image;
  vt_image imanchor;
  vt_image *ptranchor = (vt_image *)NULL;
  typeVoxelTree internalTree;
  typeTree tree;
  char *lineoptions;
  typeTreeOutputParam outputParam;


  double time_init = _GetTime();
  double time_exit;
  double clock_init = _GetClock();
  double clock_exit;


  /* parameter initialization
   */
  API_InitParam_tree( &par );



  /* parameter parsing
   */
  if ( argc <= 1 )
      API_ErrorParse_tree( _BaseName( argv[0] ), (char*)NULL, 0 );
  API_ParseParam_tree( 1, argc, argv, &par );
  


  /* input image reading
   */
  VT_Image( &image );
  if ( VT_ReadInrimage( &image, par.input_name ) != 1 ) {
    API_ErrorParse_tree( _BaseName( argv[0] ), "unable to read input image ...\n", 0 );
  }

  if ( par.input_inv == 1 )  VT_InverseImage( &image );
  if ( par.input_swap == 1 ) VT_SwapImage( &image );

  VT_Image( &imanchor );
  if ( par.anchor_name != (char*)NULL && par.anchor_name[0] != '\0' ) {
    if ( VT_ReadInrimage( &imanchor , par.anchor_name ) != 1 ) {
      VT_FreeImage( &image );
      API_ErrorParse_tree( _BaseName( argv[0] ), "unable to read anchor image ...\n", 0 );
    }
    ptranchor = &imanchor;
  }


  /* compute internal tree
   */
  initVoxelTree( &internalTree );

  /* API call
   */

  lineoptions = _Array2Str( argc, argv );
  if ( lineoptions == (char*)NULL ) {
      if ( ptranchor != (vt_image*)NULL ) VT_FreeImage( &imanchor );
      VT_FreeImage( &image );
      API_ErrorParse_tree( _BaseName( argv[0] ), "unable to translate command line options ...\n", 0 );
  }
  if ( API_tree( &image, ptranchor, &internalTree, lineoptions, (char*)NULL ) != 1 ) {
      vtfree( lineoptions );
      freeVoxelTree( &internalTree );
      if ( ptranchor != (vt_image*)NULL ) VT_FreeImage( &imanchor );
      VT_FreeImage( &image );
      API_ErrorParse_tree( _BaseName( argv[0] ), "some error occurs during processing ...\n", -1 );
  }
  vtfree( lineoptions );
  if ( ptranchor != (vt_image*)NULL ) VT_FreeImage( &imanchor );



  /* output images
   */

  if ( _writeOutputImage( &internalTree, par.output_name, UCHAR, &image, _NO_ATTRIBUTE_ ) != 1 ) {
    freeVoxelTree( &internalTree );
    VT_FreeImage( &image );
    API_ErrorParse_tree( _BaseName( argv[0] ), "unable to write output image ...\n", -1 );
  }

  if ( _writeOutputImage( &internalTree, par.output_branches_name, UCHAR, &image, _BRANCH_LABEL_ ) != 1 ) {
    freeVoxelTree( &internalTree );
    VT_FreeImage( &image );
    API_ErrorParse_tree( _BaseName( argv[0] ), "unable to write branch output image ...\n", -1 );
  }

  if ( _writeOutputImage( &internalTree, par.output_components_name, UCHAR, &image, _COMPONENT_LABEL_ ) != 1 ) {
    freeVoxelTree( &internalTree );
    VT_FreeImage( &image );
    API_ErrorParse_tree( _BaseName( argv[0] ), "unable to write component output image ...\n", -1 );
  }

  if ( _writeOutputImage( &internalTree, par.output_neighbors_name, UCHAR, &image, _POINT_NEIGHBORS_ ) != 1 ) {
    freeVoxelTree( &internalTree );
    VT_FreeImage( &image );
    API_ErrorParse_tree( _BaseName( argv[0] ), "unable to write neighbor output image ...\n", -1 );
  }

  if ( _writeOutputImage( &internalTree, par.output_order_name, UCHAR, &image, _BRANCH_ORDER_ ) != 1 ) {
    freeVoxelTree( &internalTree );
    VT_FreeImage( &image );
    API_ErrorParse_tree( _BaseName( argv[0] ), "unable to write order output image ...\n", -1 );
  }



  /* memory freeing
   */
  VT_FreeImage( &image );



  /* output text files
   */
  if ( writeComponentLengthFromVoxelTree( par.text_components_length,
                                          &internalTree, par.write_text_details ) != 1 ) {
    freeVoxelTree( &internalTree );
    API_ErrorParse_tree( _BaseName( argv[0] ), "unable to write component length text file ...\n", -1 );
  }

  if ( writeBranchLengthFromVoxelTree( par.text_branches_length,
                                          &internalTree, par.write_text_details ) != 1 ) {
    freeVoxelTree( &internalTree );
    API_ErrorParse_tree( _BaseName( argv[0] ), "unable to write branch length text file ...\n", -1 );
  }



  /* tree structure files
   */
  initTree( &tree );
  if ( voxelTreeToTree( &internalTree, &tree ) != 1 ) {
    freeTree( &tree );
    freeVoxelTree( &internalTree );
    API_ErrorParse_tree( _BaseName( argv[0] ), "unable to build tree ...\n", -1 );
  }

  freeVoxelTree( &internalTree );

  initTreeOutputParam( &outputParam );

  outputParam.color = _COLOR_ANCHOR_;
  if ( writeVTKLegacyFile( par.vtktree_name, &tree, &outputParam, (char *)NULL ) != 1 ) {
    freeTree( &tree );
    API_ErrorParse_tree( _BaseName( argv[0] ), "unable to write vtk file ...\n", -1 );
  }

  outputParam.color = _COLOR_BRANCH_LABEL_;
  if ( writeVTKLegacyFile( par.vtktree_branches_name, &tree, &outputParam, (char *)NULL ) != 1 ) {
    freeTree( &tree );
    API_ErrorParse_tree( _BaseName( argv[0] ), "unable to write vtk branch file ...\n", -1 );
  }

  outputParam.color = _COLOR_COMPONENT_LABEL_;
  if ( writeVTKLegacyFile( par.vtktree_components_name, &tree, &outputParam, (char *)NULL ) != 1 ) {
    freeTree( &tree );
    API_ErrorParse_tree( _BaseName( argv[0] ), "unable to write vtk component file ...\n", -1 );
  }

  outputParam.color = _COLOR_BRANCH_ORDER_;
  if ( writeVTKLegacyFile( par.vtktree_order_name, &tree, &outputParam, (char *)NULL ) != 1 ) {
    freeTree( &tree );
    API_ErrorParse_tree( _BaseName( argv[0] ), "unable to write vtk order file ...\n", -1 );
  }

  freeTree( &tree );







/*

  if ( par.vtktree_name[0] != '\0' || par.length_name[0] != '\0' ) {
    initTree( &tree );
    lineoptions = _Array2Str( argc, argv );
    if ( lineoptions == (char*)NULL ) {
        if ( resim != &image ) VT_FreeImage( &imres );
        if ( ptranchor != (vt_image*)NULL ) VT_FreeImage( &imanchor );
        VT_FreeImage( &image );
        API_ErrorParse_tree( _BaseName( argv[0] ), "unable to translate command line options ...\n", 0 );
    }

    if ( API_symbolicTree( resim, ptranchor, &tree, lineoptions, (char*)NULL ) != 1 ) {
        vtfree( lineoptions );
        if ( resim != &image ) VT_FreeImage( &imres );
        if ( ptranchor != (vt_image*)NULL ) VT_FreeImage( &imanchor );
        VT_FreeImage( &image );
        API_ErrorParse_tree( _BaseName( argv[0] ), "some error occurs during processing ...\n", -1 );
    }

    if ( par.vtktree_name[0] != '\0' ) {
        if ( writeVTKLegacyFile( par.vtktree_name, &tree, &(par.treeOutputParam), (char*)NULL ) != 1 ) {
            freeTree( &tree );
            vtfree( lineoptions );
            if ( resim != &image ) VT_FreeImage( &imres );
            if ( ptranchor != (vt_image*)NULL ) VT_FreeImage( &imanchor );
            VT_FreeImage( &image );
            API_ErrorParse_tree( _BaseName( argv[0] ), "some error occurs during processing ...\n", -1 );
        }
    }

    if ( par.length_name[0] != '\0' ) {
        if ( writeLengthFile( par.length_name, &tree, &(par.treeOutputParam) ) != 1 ) {
            freeTree( &tree );
            vtfree( lineoptions );
            if ( resim != &image ) VT_FreeImage( &imres );
            if ( ptranchor != (vt_image*)NULL ) VT_FreeImage( &imanchor );
            VT_FreeImage( &image );
            API_ErrorParse_tree( _BaseName( argv[0] ), "some error occurs during processing ...\n", -1 );
        }
    }
    freeTree( &tree );
    vtfree( lineoptions );
  }

  */
  




  time_exit = _GetTime();
  clock_exit = _GetClock();

  if ( par.print_time ) { 
    fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName( argv[0] ), time_exit - time_init );
    fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
    fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );
  }


  return( 0 );
}





/************************************************************
 *
 *
 *
 ************************************************************/

int _writeOutputImage( typeVoxelTree *internalTree,
                       char *name, bufferType type, vt_image *image,
                       enumVoxelTreeAttribute attribute )
{
  char *proc = "_writeOutputImage";
  vt_image imres;
  int theDim[3];

  theDim[0] = image->dim.x;
  theDim[1] = image->dim.y;
  theDim[2] = image->dim.z;

  if ( name == (char*)NULL || name[0] == '\0' || name[0] == '>' )
    return( 1 );

  VT_Image( &imres );
  VT_InitFromImage( &imres, image, name, type );
  if ( VT_AllocImage( &imres ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when allocating output image\n", proc );
    return( -1 );
  }

  if ( voxelTreeToImage( internalTree, imres.buf, type, theDim, attribute ) != 1 ) {
    VT_FreeImage( &imres );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when translating tree into image\n", proc );
    return( -1 );
  }

  if ( VT_WriteInrimage( &imres ) == -1 ) {
    VT_FreeImage( &imres );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when writing output image\n", proc );
    return( -1 );
  }

  VT_FreeImage( &imres );

  return( 1 );
}


/************************************************************
 *
 * static functions
 *
 ************************************************************/



static char *_Array2Str( int argc, char *argv[] )
{
  char *proc = "_Array2Str";
  int i, l;
  char *s, *t;

  if ( argc <= 1 || argv == (char**)NULL ) {
    if ( _verbose_ >= 2 )
      fprintf( stderr, "%s: no options in argv[]\n", proc );
    return( (char*)NULL );
  }

  /* there are argc-1 strings
   * compute the sum of string lengths from 1 to argc-1
   * + number of interval between successive strings (argc-2)
   * + 1 to add a trailing '\0'
   */
  for ( l=argc-1, i=1; i<argc; i++ ) {
    l += strlen( argv[i] );
  }

  s = (char*)vtmalloc( l * sizeof( char ), "s", proc );
  if ( s == (char*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( (char*)NULL );
  }

  for ( t=s, i=1; i<argc; i++ ) {
    (void)strncpy( t, argv[i], strlen( argv[i] ) );
    t += strlen( argv[i] );
    if ( i < argc-1 ) {
      *t = ' ';
      t++;
    }
    else {
      *t = '\0';
    }
  }

  return( s );
}



static char *_BaseName( char *p )
{
  int l;
  if ( p == (char*)NULL ) return( (char*)NULL );
  l = strlen( p ) - 1;
  while ( l >= 0 && p[l] != '/' ) l--;
  if ( l < 0 ) l = 0;
  if ( p[l] == '/' ) l++;
  return( &(p[l]) );
}



static double _GetTime()
{
  struct timeval tv;
  gettimeofday(&tv, (void *)0);
  return ( (double) tv.tv_sec + tv.tv_usec*1e-6 );
}



static double _GetClock() 
{
  return ( (double) clock() / (double)CLOCKS_PER_SEC );
}
