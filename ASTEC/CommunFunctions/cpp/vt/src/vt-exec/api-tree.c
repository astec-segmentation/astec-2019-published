/*************************************************************************
 * api-tree.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2016, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mar 22 mar 2016 12:03:24 CET
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <chunks.h>
#include <topological-component-image.h>
#include <topological-voxel-tree-main.h>
#include <topological-voxel-tree-order.h>
#include <threshold.h>
#include <vtmalloc.h>

#include <vt_common.h>

#include <api-tree.h>






static int _verbose_ = 1;
static int _debug_ = 0;


static void _API_ParseParam_tree( char *str, lineCmdParamTree *p );



/************************************************************
 *
 * main APIs
 *
 ************************************************************/



int API_tree( vt_image *image, vt_image *anchor,
              typeVoxelTree *iTree,
              char *param_str_1, char *param_str_2 )
{
  char *proc = "API_tree";
  lineCmdParamTree par;
  int theDim[3];
  float theSize[3];
  void *theBuf = (void*)NULL;
  bufferType theType = TYPE_UNKNOWN;
  void *ancBuf = (void*)NULL;
  bufferType ancType =TYPE_UNKNOWN;



  /* parameter initialization
   */
  API_InitParam_tree( &par );

  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_tree( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_tree( param_str_2, &par );

  if ( par.print_lineCmdParam )
      API_PrintParam_tree( stderr, proc, &par, (char*)NULL );

  /************************************************************
   *
   *  here is the stuff
   *
   ************************************************************/

  theDim[0] = image->dim.x;
  theDim[1] = image->dim.y;
  theDim[2] = image->dim.z;

  theSize[0] = image->siz.x;
  theSize[1] = image->siz.y;
  theSize[2] = image->siz.z;

  theBuf = image->buf;
  theType = image->type;

  if ( anchor != (vt_image*)NULL ) {
    ancBuf = anchor->buf;
    ancType = anchor->type;
  }

  if ( buildVoxelTree( theBuf, theType, ancBuf, ancType, theDim, theSize,
                       iTree,
                       par.realBranchLength,  par.voxelBranchLength, par.maxEndEdges,
                       par.maxorder ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when building voxel tree\n", proc );
    return( -1 );
  }

  return( 1 );
}





int API_symbolicTree( vt_image *image, vt_image *anchor,
                      typeTree *tree, char *param_str_1, char *param_str_2 )
{
  char *proc = "API_symbolicTree";
  lineCmdParamTree par;
  int theDim[3];
  typeComponentImage treeImage;
  void *ancBuf = (void*)NULL;
  bufferType ancType =TYPE_UNKNOWN;

  /* parameter initialization
   */
  API_InitParam_tree( &par );

  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_tree( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_tree( param_str_2, &par );

  if ( par.print_lineCmdParam )
      API_PrintParam_tree( stderr, proc, &par, (char*)NULL );

  /************************************************************
   *
   *  here is the stuff
   *
   ************************************************************/

  theDim[0] = image->dim.x;
  theDim[1] = image->dim.y;
  theDim[2] = image->dim.z;

  if ( anchor != (vt_image*)NULL ) {
    ancBuf = anchor->buf;
    ancType = anchor->type;
  }


  initComponentImage( &treeImage );
  if ( allocComponentImage( &treeImage, USHORT, theDim ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error for tree image\n", proc );
    return( -1 );
  }

  treeImage.voxelSize[0] = image->siz.x;
  treeImage.voxelSize[1] = image->siz.y;
  treeImage.voxelSize[2] = image->siz.z;
  if ( imageToComponentImage( image->buf, image->type, &treeImage ) != 1 ) {
    freeComponentImage( &treeImage );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when labeling components in tree\n", proc );
    return( -1 );
  }

  if ( treeImageToTree( &treeImage, ancBuf, ancType, tree ) != 1 ) {
    freeComponentImage( &treeImage );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when building tree\n", proc );
    return( -1 );
  }

  freeComponentImage( &treeImage );
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



static char *usage = "[image-in] [[-output-image|-oi] image-out]\n\
 [-anchor-image|-anchor %s]\n\
 [-pruning-real-length|-pruning-length|-prl %f]\n\
 [-pruning-voxel-length|-pvl %d]\n\
 [-pruning-end-branches|-pruning-end-edges|-maximal-end-branches|-pe %d]\n\
 [-maximal-order %d]\n\
 [-output-image-branches|-oib %s]\n\
 [-output-image-components|-oic %s]\n\
 [-output-image-neighbors|-oin %s]\n\
 [-output-image-order|-oio %s]\n\
 [-output-vtk|-ov|-vtk-output|-vtk %s]\n\
 [-output-vtk-branches|-ovb %s]\n\
 [-output-vtk-components|-ovn %s]\n\
 [-output-vtk-order|-ovo %s]\n\
 [-write-length-details|-wld]\n\
 [-no-write-length-details|-nwld]\n\
 [-length-components|-length-output|-length|-lc %s]\n\
 [-length-branches|-lb %s]\n\
 [-parallel|-no-parallel] [-max-chunks %d]\n\
 [-parallelism-type|-parallel-type default|none|openmp|omp|pthread|thread]\n\
 [-omp-scheduling|-omps default|static|dynamic-one|dynamic|guided]\n\
 [-inv] [-swap] [output-image-type | -type s8|u8|s16|u16...]\n\
 [-verbose|-v] [-nv|-noverbose] [-debug|-D] [-nodebug]\n\
 [-print-parameters|-param]\n\
 [-print-time|-time] [-notime]\n\
 [-help|-h]";

/*
 *  [-tree-type all-branches|inner-branches|end-branches|all|inner|end]\n\
         [-vtk-tree-files single|multiple|s|m]\n\
         [-vtk-color marked-points|marked]\n\
         [-vtk-label all-edges|all|inner-edges|inner|end-edges|end]\n\
                 */

static char *detail = "\
 if 'image-in' is equal to '-', stdin will be used\n\
 if 'image-out' is not specified or equal to '-', stdout will be used\n\
 if both are not specified, stdin and stdout will be used\n\
# ...\n\
 build output from a skeleton image\n\
# input files\n\
 -anchor-image|-anchor %s: points that can not be deleted\n\
# processing parameters\n\
 -pruning-real-length|-pruning-length|-prl %f\n\
   remove iteratively all end branches whose length (in real units)\n\
   is shorter that the given value\n\
 -pruning-voxel-length|-pvl %d\n\
   remove iteratively all end components whose length (in voxel)\n\
   is shorter than the given value. If '-prune-length' is also\n\
   given, both conditions stand\n\
 -pruning-end-branches|-pruning-end-edges|-maximal-end-branches|-pe %d:\n\
  maximal number of end branches to be kept.\n\
  Removes end components (the smallest first) until this number is reached\n\
 -maximal-order %d:\n\
  maximal order for components to be kept\n\
  1 corresponds to anchor tree\n\
  2 corresponds to longest branches connected to the anchor tree\n\
  etc.\n\
# output image files\n\
 -output-image-branches|-oib %s:\n\
   one label per branche; branches of order n are recursively defined\n\
   as the longest path in a subtree, issued from order n-1 branches\n\
   with order 1 being the anchor tree\n\
 -output-image-components|-oic %s:\n\
   one label par component. A component is either a curve/edge\n\
   or a junction. Junctions may be formed by several points.\n\
 -output-image-neighbors|-oin %s:\n\
   counts the number of neighbors around each point.\n\
   For trees made of curves:\n\
   o end points have 1 neighbor\n\
   o curve points have 2 neighbors\n\
   o junction points have more than 2 neighbors\n\
 -output-image-order|-oio %s:\n\
   each point is assigned the order of the branch it belongs\n\
   1 is for the anchor tree points\n\
# output vtk files\n\
 -output-vtk|-ov|-vtk-output|-vtk %s:\n\
   vtk file (eg for paraview vizualisation)\n\
 -output-vtk-branches|-ovb %s:\n\
   vtk file, points are colored by branch label\n\
 -output-vtk-components|-ovn %s:\n\
   vtk file, points are colored by component label\n\
 -output-vtk-order|-ovo %s:\n\
   vtk file, points are colored by branch order\n\
# text output files\n\
 -write-length-details|-wld:\n\
   write component/branch label in addition to length\n\
 -no-write-length-details|-nwld:\n\
 -length-components|-length-output|-length|-lc %s]\n\
   write component lengths\n\
 -length-branches|-lb %s:\n\
   write branch lengths\n\
# parallelism parameters\n\
 -parallel|-no-parallel:\n\
 -max-chunks %d:\n\
 -parallelism-type|-parallel-type default|none|openmp|omp|pthread|thread:\n\
 -omp-scheduling|-omps default|static|dynamic-one|dynamic|guided:\n\
# general image related parameters\n\
  -inv: inverse 'image-in'\n\
  -swap: swap 'image-in' (if encoded on 2 or 4 bytes)\n\
   output-image-type: -o 1    : unsigned char\n\
                      -o 2    : unsigned short int\n\
                      -o 2 -s : short int\n\
                      -o 4 -s : int\n\
                      -r      : float\n\
  -type s8|u8|s16|u16|... \n\
   default is type of input image\n\
# general parameters \n\
  -verbose|-v: increase verboseness\n\
    parameters being read several time, use '-nv -v -v ...'\n\
    to set the verboseness level\n\
  -noverbose|-nv: no verboseness at all\n\
  -debug|-D: increase debug level\n\
  -nodebug: no debug indication\n\
  -print-parameters|-param:\n\
  -print-time|-time:\n\
  -no-time|-notime:\n\
  -h: print option list\n\
  -help: print option list + details\n\
";


 /*
  *
    -image-type|-output neighbors|components|uniform:\n\
      type of output/labeling for the image output.\n\
      neighbors: count the number of neighbors\n\
      components: label each basic component (ie edges and junction)\n\
      uniform:\n\
      branches: label each branch (may be formed by several edges)\n\

    -vtk-output|-vtk %s\n\
     vtk file (eg for paraview vizualisation)\n\
    -length-output|-length %s\n\
      output file containing the branch lengths (one per line)\n\
    -length-with-label|-lwl\n\
      each line of the above line is of the form 'label: length'\n\
    -tree-type all-branches|inner-branches|end-branches|all|inner|end\n\
      selection of the tree branches to be outputed\n\
    -vtk-tree-files single|multiple\n\
      a file per branch or a file for all branches\n\
    -vtk-color marked-points|marked\n\
      marked points (issued from the anchor image) are colored at 1\n\
      and the other at 0\n\
    -vtk-label all-edges|all|inner-edges|inner|end-edges|end\n\
      gives a different label/color to edges\n\
    */



char *API_Help_tree( int h )
{
    if ( h == 0 )
        return( usage );
    return( detail );
}





void API_ErrorParse_tree( char *program, char *str, int flag )
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



void API_InitParam_tree( lineCmdParamTree *p )
{
    (void)strncpy( p->input_name, "\0", 1 );
    (void)strncpy( p->anchor_name, "\0", 1 );

    p->realBranchLength = -1.0;
    p->voxelBranchLength = -1;
    p->maxEndEdges = -1;
    p->maxorder = -1;


    (void)strncpy( p->output_name, "\0", 1 );
    (void)strncpy( p->output_branches_name, "\0", 1 );
    (void)strncpy( p->output_components_name, "\0", 1 );
    (void)strncpy( p->output_neighbors_name, "\0", 1 );
    (void)strncpy( p->output_order_name, "\0", 1 );

    (void)strncpy( p->vtktree_name, "\0", 1 );
    (void)strncpy( p->vtktree_branches_name, "\0", 1 );
    (void)strncpy( p->vtktree_components_name, "\0", 1 );
    (void)strncpy( p->vtktree_order_name, "\0", 1 );

    p->write_text_details = 1;
    (void)strncpy( p->text_components_length, "\0", 1 );
    (void)strncpy( p->text_branches_length, "\0", 1 );

    /*
    p->imageOutputType = _TREE_IMAGE_UNIFORM_;

    p->realBranchLength = -1.0;
    p->voxelBranchLength = -1;
    p->maxEndEdges = -1;

    initTreeOutputParam( &(p->treeOutputParam) );
    */

    p->input_inv = 0;
    p->input_swap = 0;
    p->output_type = TYPE_UNKNOWN;

    p->print_lineCmdParam = 0;
    p->print_time = 0;
}





void API_PrintParam_tree( FILE *theFile, char *program,
                                  lineCmdParamTree *p, char *str )
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


  fprintf( f, "# input names\n" );

  fprintf( f, "- input image is " );
  if ( p->input_name != (char*)NULL && p->input_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- anchor image is " );
  if ( p->anchor_name != (char*)NULL && p->anchor_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->anchor_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "# tree pruning parameters\n" );
  fprintf( f, "- minimal branch length (real units) = %f\n", p->realBranchLength );
  fprintf( f, "- minimal branch length (voxels) = %d\n", p->voxelBranchLength );
  fprintf( f, "- maximal number of end branches = %d\n", p->maxEndEdges );
  fprintf( f, "- maximal order = %d\n", p->maxorder );


  fprintf( f, "# output image names\n" );

  fprintf( f, "- output image name is " );
  if ( p->output_name != (char*)NULL && p->output_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- output branch image name is " );
  if ( p->output_branches_name != (char*)NULL && p->output_branches_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_branches_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- output component image name is " );
  if ( p->output_components_name != (char*)NULL && p->output_components_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_components_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- output neighbor image name is " );
  if ( p->output_neighbors_name != (char*)NULL && p->output_neighbors_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_neighbors_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- output order image name is " );
  if ( p->output_order_name != (char*)NULL && p->output_order_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_order_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "# output vtk tree names\n" );

  fprintf( f, "- vtk tree file is " );
  if ( p->vtktree_name != (char*)NULL && p->vtktree_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->vtktree_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- branch-colored vtk tree file is " );
  if ( p->vtktree_branches_name != (char*)NULL && p->vtktree_branches_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->vtktree_branches_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- component-colored vtk tree file is " );
  if ( p->vtktree_components_name != (char*)NULL && p->vtktree_components_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->vtktree_components_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "-  order-colored tree file is " );
  if ( p->vtktree_order_name != (char*)NULL && p->vtktree_order_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->vtktree_order_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "# output text names\n" );

  fprintf( f, "- write_text_details = %d\n ", p->write_text_details );

  fprintf( f, "- component length file is " );
  if ( p->text_components_length != (char*)NULL && p->text_components_length[0] != '\0' )
    fprintf( f, "'%s'\n", p->text_components_length );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- branch length file is " );
  if ( p->text_branches_length != (char*)NULL && p->text_branches_length[0] != '\0' )
    fprintf( f, "'%s'\n", p->text_branches_length );
  else
    fprintf( f, "'NULL'\n" );


  fprintf( f, "# ...\n" );

  /*
  fprintf( f, "- image output type is " );
  switch ( p->imageOutputType ) {
  default :     fprintf( f, "unknown\n" ); break;
  case _TREE_IMAGE_NEIGHBORS_ :  fprintf( f, "_TREE_IMAGE_NEIGHBORS_\n" ); break;
  case _TREE_IMAGE_COMPONENTS_ : fprintf( f, "_TREE_IMAGE_COMPONENTS_\n" ); break;
  case _TREE_IMAGE_UNIFORM_ :    fprintf( f, "_TREE_IMAGE_UNIFORM_\n" ); break;
  }

  fprintf( f, "# tree output parameters\n" );
  fprintfTreeOutputParam( f, &(p->treeOutputParam), "treeOutputParam" );

  */

  fprintf( f, "# ...\n" );

  fprintf( f, "# general image related parameters\n" );

  fprintf( f, "- input image inverse = %d\n", p->input_inv );
  fprintf( f, "- input image swap = %d\n", p->input_swap );
  fprintf( f, "- output image type = " );
  switch ( p->output_type ) {
  default :     fprintf( f, "TYPE_UNKNOWN\n" ); break;
  case SCHAR :  fprintf( f, "SCHAR\n" ); break;
  case UCHAR :  fprintf( f, "UCHAR\n" ); break;
  case SSHORT : fprintf( f, "SSHORT\n" ); break;
  case USHORT : fprintf( f, "USHORT\n" ); break;
  case UINT :   fprintf( f, "UINT\n" ); break;
  case SINT :   fprintf( f, "INT\n" ); break;
  case ULINT :  fprintf( f, "ULINT\n" ); break;
  case FLOAT :  fprintf( f, "FLOAT\n" ); break;
  case DOUBLE : fprintf( f, "DOUBLE\n" ); break;
  }

  fprintf( f, "==================================================\n" );
}





/************************************************************
 *
 * parameters parsing
 *
 ************************************************************/



static void _API_ParseParam_tree( char *str, lineCmdParamTree *p )
{
  char *proc = "_API_ParseParam_tree";
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

  API_ParseParam_tree( 0, argc, argv, p );

  vtfree( argv );
}





static int _n_call_parse_ = 0;

void API_ParseParam_tree( int firstargc, int argc, char *argv[],
                                  lineCmdParamTree *p )
{
  int i;
  int inputisread = 0;
  int outputisread = 0;
  char text[STRINGLENGTH];
  int status;
  int maxchunks;
  int o=0, s=0, r=0;

  _n_call_parse_ ++;

  /* option line parsing
   */
  for ( i=firstargc; i<argc; i++ ) {

      /* strings beginning with '-'
       */
      if ( argv[i][0] == '-' ) {
          if ( argv[i][1] == '\0' ) {
            if ( inputisread == 0 ) {
              (void)strcpy( p->input_name,  "<" );  /* standart input */
              inputisread = 1;
            }
            else if ( outputisread == 0 ) {
              (void)strcpy( p->output_name,  ">" );  /* standart output */
              outputisread = 1;
            }
            else {
              API_ErrorParse_tree( (char*)NULL, "too many file names, parsing '-' ...\n", 0 );
            }
          }

          /* anchor image
           */
          else if ( strcmp ( argv[i], "-anchor-image" ) == 0 ||
                    (strcmp ( argv[i], "-anchor" ) == 0 && argv[i][7] == '\0') ) {
            i++;
            if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -anchor-image ...\n", 0 );
            if ( strlen( argv[i] ) >= STRINGLENGTH ) {
                fprintf( stderr, "... parsing '%s'\n", argv[i] );
                API_ErrorParse_tree( (char*)NULL, "too long file name ...\n", 0 );
            }
            (void)strcpy( p->anchor_name, argv[i] );
          }

          /* pruning parameters
           */
          else if ( strcmp ( argv[i], "-pruning-length" ) == 0 ||
                    strcmp ( argv[i], "-pruning-real-length" ) == 0 ||
                    (strcmp ( argv[i], "-prl" ) == 0 && argv[i][4] == '\0') ) {
              i ++;
              if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -pruning-real-length ...\n", 0 );
              status = sscanf( argv[i], "%f", &(p->realBranchLength) );
              if ( status <= 0 ) API_ErrorParse_tree( (char*)NULL, "parsing -pruning-real-length ...\n", 0 );
          }
          else if ( strcmp ( argv[i], "-pruning-voxel-length" ) == 0 ||
                    (strcmp ( argv[i], "-pvl" ) == 0 && argv[i][4] == '\0') ) {
              i ++;
              if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -pruning-voxel-length ...\n", 0 );
              status = sscanf( argv[i], "%d", &(p->voxelBranchLength) );
              if ( status <= 0 ) API_ErrorParse_tree( (char*)NULL, "parsing -pruning-voxel-length ...\n", 0 );
          }
          else if ( strcmp ( argv[i], "-pruning-end-branches" ) == 0 ||
                    strcmp ( argv[i], "-pruning-end-edges" ) == 0 ||
                    strcmp ( argv[i], "-maximal-end-branches" ) == 0 ||
                    (strcmp ( argv[i], "-pe" ) == 0 && argv[i][3] == '\0') ) {
              i ++;
              if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -pruning-end-branches ...\n", 0 );
              status = sscanf( argv[i], "%d", &(p->maxEndEdges) );
              if ( status <= 0 ) API_ErrorParse_tree( (char*)NULL, "parsing -pruning-end-branches ...\n", 0 );
          }

          else if ( strcmp ( argv[i], "-maximal-order" ) == 0 ) {
              i ++;
              if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -maximal-order ...\n", 0 );
              status = sscanf( argv[i], "%d", &(p->maxorder) );
              if ( status <= 0 ) API_ErrorParse_tree( (char*)NULL, "parsing -maximal-order ...\n", 0 );
          }

          /* output images
           */
          else if ( (strcmp ( argv[i], "-output-image" ) == 0 && argv[i][13] == '\0') ||
                    (strcmp ( argv[i], "-oi" ) == 0 && argv[i][3] == '\0') ) {
            i++;
            if ( i >= argc)
              API_ErrorParse_tree( (char*)NULL, "parsing -output-image ...\n", 0 );
            if ( outputisread )
              API_ErrorParse_tree( (char*)NULL, "parsing -output-image ...\n", 0 );
            if ( strlen( argv[i] ) >= STRINGLENGTH ) {
                fprintf( stderr, "... parsing '%s'\n", argv[i] );
                API_ErrorParse_tree( (char*)NULL, "too long file name ...\n", 0 );
            }
            (void)strcpy( p->output_name, argv[i] );
          }

          else if ( (strcmp ( argv[i], "-output-image-branches" ) == 0 && argv[i][22] == '\0') ||
                    (strcmp ( argv[i], "-oib" ) == 0 && argv[i][4] == '\0') ) {
            i++;
            if ( i >= argc)
              API_ErrorParse_tree( (char*)NULL, "parsing -output-image-branches ...\n", 0 );
            if ( strlen( argv[i] ) >= STRINGLENGTH ) {
                fprintf( stderr, "... parsing '%s'\n", argv[i] );
                API_ErrorParse_tree( (char*)NULL, "too long file name ...\n", 0 );
            }
            (void)strcpy( p->output_branches_name, argv[i] );
          }

          else if ( (strcmp ( argv[i], "-output-image-components" ) == 0 && argv[i][24] == '\0') ||
                    (strcmp ( argv[i], "-oic" ) == 0 && argv[i][4] == '\0') ) {
            i++;
            if ( i >= argc)
              API_ErrorParse_tree( (char*)NULL, "parsing -output-image-components ...\n", 0 );
            if ( strlen( argv[i] ) >= STRINGLENGTH ) {
                fprintf( stderr, "... parsing '%s'\n", argv[i] );
                API_ErrorParse_tree( (char*)NULL, "too long file name ...\n", 0 );
            }
            (void)strcpy( p->output_components_name, argv[i] );
          }

          else if ( (strcmp ( argv[i], "-output-image-neighbors" ) == 0 && argv[i][23] == '\0') ||
                    (strcmp ( argv[i], "-oin" ) == 0 && argv[i][4] == '\0') ) {
            i++;
            if ( i >= argc)
              API_ErrorParse_tree( (char*)NULL, "parsing -output-image-neighbors ...\n", 0 );
            if ( strlen( argv[i] ) >= STRINGLENGTH ) {
                fprintf( stderr, "... parsing '%s'\n", argv[i] );
                API_ErrorParse_tree( (char*)NULL, "too long file name ...\n", 0 );
            }
            (void)strcpy( p->output_neighbors_name, argv[i] );
          }

          else if ( (strcmp ( argv[i], "-output-image-order" ) == 0 && argv[i][19] == '\0') ||
                    (strcmp ( argv[i], "-oio" ) == 0 && argv[i][3] == '\0') ) {
            i++;
            if ( i >= argc)
              API_ErrorParse_tree( (char*)NULL, "parsing -output-image-order ...\n", 0 );
            if ( strlen( argv[i] ) >= STRINGLENGTH ) {
                fprintf( stderr, "... parsing '%s'\n", argv[i] );
                API_ErrorParse_tree( (char*)NULL, "too long file name ...\n", 0 );
            }
            (void)strcpy( p->output_order_name, argv[i] );
          }




          /* output tree structures
           */
          else if ( (strcmp ( argv[i], "-output-vtk" ) == 0 && argv[i][11] == '\0') ||
                    (strcmp ( argv[i], "-ov" ) == 0 && argv[i][3] == '\0') ||
                    strcmp ( argv[i], "-vtk-output" ) == 0 ||
                    (strcmp ( argv[i], "-vtk" ) == 0 && argv[i][4] == '\0') ) {
            i++;
            if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -output-vtk ...\n", 0 );
            (void)strcpy( p->vtktree_name, argv[i] );
          }

          else if ( (strcmp ( argv[i], "-output-vtk-branches" ) == 0 && argv[i][20] == '\0') ||
                    (strcmp ( argv[i], "-ovb" ) == 0 && argv[i][4] == '\0') ) {
            i++;
            if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -output-vtk-branches ...\n", 0 );
            (void)strcpy( p->vtktree_branches_name, argv[i] );
          }

          else if ( (strcmp ( argv[i], "-output-vtk-components" ) == 0 && argv[i][22] == '\0') ||
                    (strcmp ( argv[i], "-ovc" ) == 0 && argv[i][4] == '\0') ) {
            i++;
            if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -output-vtk-components ...\n", 0 );
            (void)strcpy( p->vtktree_components_name, argv[i] );
          }

          else if ( (strcmp ( argv[i], "-output-vtk-order" ) == 0 && argv[i][17] == '\0') ||
                    (strcmp ( argv[i], "-ovo" ) == 0 && argv[i][4] == '\0') ) {
            i++;
            if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -output-vtk-order ...\n", 0 );
            (void)strcpy( p->vtktree_order_name, argv[i] );
          }


          /* output text files
           */

          else if ( strcmp ( argv[i], "-write-length-details" ) == 0 ||
                    (strcmp ( argv[i], "-wld" ) == 0 && argv[i][4] == '\0') ) {
            p->write_text_details = 1;
          }
          else if ( strcmp ( argv[i], "-no-write-length-details" ) == 0 ||
                    (strcmp ( argv[i], "-nwld" ) == 0 && argv[i][4] == '\0') ) {
            p->write_text_details = 0;
          }

          else if ( strcmp ( argv[i], "-length-components" ) == 0 ||
                    strcmp ( argv[i], "-length-output" ) == 0 ||
                    (strcmp ( argv[i], "-length" ) == 0 && argv[i][7] == '\0') ||
                    (strcmp ( argv[i], "-lc" ) == 0 && argv[i][3] == '\0') ) {
            i++;
            if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -length-components ...\n", 0 );
            (void)strcpy( p->text_components_length, argv[i] );
          }

          else if ( strcmp ( argv[i], "-length-branches" ) == 0 ||
                    (strcmp ( argv[i], "-lb" ) == 0 && argv[i][3] == '\0') ) {
            i++;
            if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -length-branches ...\n", 0 );
            (void)strcpy( p->text_branches_length, argv[i] );
          }


          /* ...
           */


          /*
          else if ( strcmp ( argv[i], "-tree-type" ) == 0 ) {
              i++;
              if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -tree-type ...\n", 0 );
              if ( strcmp ( argv[i], "all-branches" ) == 0
                   || strcmp ( argv[i], "all-edges" ) == 0
                   || (strcmp ( argv[i], "all" ) == 0 && argv[i][3] == '\0') ) {
                p->treeOutputParam.edgeToBeWritten = _UNDEFINED_EDGE_;
              }
              else if ( strcmp ( argv[i], "inner-branches" ) == 0
                        || strcmp ( argv[i], "inner-edges" ) == 0
                        || (strcmp ( argv[i], "inner" ) == 0 && argv[i][5] == '\0') ) {
                p->treeOutputParam.edgeToBeWritten = _INNER_EDGE_;
              }
              else if ( strcmp ( argv[i], "end-branches" ) == 0
                        || strcmp ( argv[i], "end-edges" ) == 0
                        || (strcmp ( argv[i], "end" ) == 0 && argv[i][3] == '\0') ) {
                p->treeOutputParam.edgeToBeWritten = _END_EDGE_;
              }
              else {
                fprintf( stderr, "unknown tree output type: '%s'\n", argv[i] );
                API_ErrorParse_tree( (char*)NULL, "parsing -tree-type ...\n", 0 );
              }
          }

          else if ( strcmp ( argv[i], "-vtk-tree-files" ) == 0 ) {
              i++;
              if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -vtk-tree-files ...\n", 0 );
              if ( (strcmp ( argv[i], "single" ) == 0 && argv[i][6] == '\0')
                   || (strcmp ( argv[i], "s" ) == 0 && argv[i][1] == '\0') ) {
                p->treeOutputParam.filetype = _SINGLE_FILE_;
              }
              else if ( (strcmp ( argv[i], "multiple" ) == 0 && argv[i][8] == '\0')
                        || (strcmp ( argv[i], "m" ) == 0 && argv[i][1] == '\0') ) {
                  p->treeOutputParam.filetype = _MULTIPLE_FILES_;
              }
              else {
                fprintf( stderr, "unknown file type: '%s'\n", argv[i] );
                API_ErrorParse_tree( (char*)NULL, "parsing -vtk-tree-files ...\n", 0 );
              }
          }

          else if ( strcmp ( argv[i], "-vtk-color" ) == 0 ) {
              i++;
              if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -vtk-color ...\n", 0 );
              if ( (strcmp ( argv[i], "marked-points" ) == 0 && argv[i][13] == '\0')
                   || (strcmp ( argv[i], "marked" ) == 0 && argv[i][6] == '\0') ) {
                p->treeOutputParam.color = _COLOR_MARKED_POINTS_;
              }
              else {
                fprintf( stderr, "unknown color type: '%s'\n", argv[i] );
                API_ErrorParse_tree( (char*)NULL, "parsing -vtk-color ...\n", 0 );
              }
          }

          else if ( strcmp ( argv[i], "-vtk-label" ) == 0 ) {
              i++;
              if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -vtk-label ...\n", 0 );
              if ( strcmp ( argv[i], "all-branches" ) == 0
                   || (strcmp ( argv[i], "all-edges" ) == 0 && argv[i][9] == '\0')
                   || (strcmp ( argv[i], "all" ) == 0 && argv[i][3] == '\0') ) {
                p->treeOutputParam.color = _LABEL_ALL_EDGES_;
              }
              else if ( strcmp ( argv[i], "inner-branches" ) == 0
                        || (strcmp ( argv[i], "inner-edges" ) == 0 && argv[i][9] == '\0')
                        || (strcmp ( argv[i], "inner" ) == 0 && argv[i][3] == '\0') ) {
                if ( p->treeOutputParam.color == _LABEL_END_EDGES_ )
                  p->treeOutputParam.color = _LABEL_ALL_EDGES_;
                else
                  p->treeOutputParam.color = _LABEL_INNER_EDGES_;
              }
              else if ( strcmp ( argv[i], "end-branches" ) == 0
                        || (strcmp ( argv[i], "end-edges" ) == 0 && argv[i][9] == '\0')
                        || (strcmp ( argv[i], "end" ) == 0 && argv[i][3] == '\0') ) {
                if ( p->treeOutputParam.color == _LABEL_INNER_EDGES_ )
                  p->treeOutputParam.color = _LABEL_ALL_EDGES_;
                else
                  p->treeOutputParam.color = _LABEL_END_EDGES_;
              }
              else{
                fprintf( stderr, "unknown label type: '%s'\n", argv[i] );
                API_ErrorParse_tree( (char*)NULL, "parsing -vtk-label ...\n", 0 );
              }
          }
          */


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
             if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -parallelism-type ...\n", 0 );
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
               API_ErrorParse_tree( (char*)NULL, "parsing -parallelism-type ...\n", 0 );
             }
          }

          else if ( strcmp ( argv[i], "-max-chunks" ) == 0 ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -max-chunks ...\n", 0 );
             status = sscanf( argv[i], "%d", &maxchunks );
             if ( status <= 0 ) API_ErrorParse_tree( (char*)NULL, "parsing -max-chunks ...\n", 0 );
             if ( maxchunks >= 1 ) setMaxChunks( maxchunks );
          }

          else if ( strcmp ( argv[i], "-omp-scheduling" ) == 0 ||
                   ( strcmp ( argv[i], "-omps" ) == 0 && argv[i][5] == '\0') ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -omp-scheduling, no argument\n", 0 );
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
               API_ErrorParse_tree( (char*)NULL, "parsing -omp-scheduling ...\n", 0 );
             }
          }

          /* general image related parameters
           */
          else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
             p->input_inv = 1;
          }
          else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
             p->input_swap = 1;
          }

          else if ( strcmp ( argv[i], "-r" ) == 0 && argv[i][2] == '\0' ) {
             r = 1;
          }
          else if ( strcmp ( argv[i], "-s" ) == 0 && argv[i][2] == '\0' ) {
             s = 1;
          }
          else if ( strcmp ( argv[i], "-o" ) == 0 && argv[i][2] == '\0' ) {
             i += 1;
             if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -o...\n", 0 );
             status = sscanf( argv[i],"%d",&o );
             if ( status <= 0 ) API_ErrorParse_tree( (char*)NULL, "parsing -o...\n", 0 );
          }
          else if ( strcmp ( argv[i], "-type" ) == 0 && argv[i][5] == '\0' ) {
            i += 1;
            if ( i >= argc)    API_ErrorParse_tree( (char*)NULL, "parsing -type...\n", 0 );
            if ( strcmp ( argv[i], "s8" ) == 0 && argv[i][2] == '\0' ) {
               p->output_type = SCHAR;
            }
            else if ( strcmp ( argv[i], "u8" ) == 0 && argv[i][2] == '\0' ) {
               p->output_type = UCHAR;
            }
            else if ( strcmp ( argv[i], "s16" ) == 0 && argv[i][3] == '\0' ) {
              p->output_type = SSHORT;
            }
            else if ( strcmp ( argv[i], "u16" ) == 0 && argv[i][3] == '\0' ) {
              p->output_type = USHORT;
            }
            else if ( strcmp ( argv[i], "s32" ) == 0 && argv[i][3] == '\0' ) {
              p->output_type = SINT;
            }
            else if ( strcmp ( argv[i], "u32" ) == 0 && argv[i][3] == '\0' ) {
              p->output_type = UINT;
            }
            else if ( strcmp ( argv[i], "s64" ) == 0 && argv[i][3] == '\0' ) {
              p->output_type = SLINT;
            }
            else if ( strcmp ( argv[i], "u64" ) == 0 && argv[i][3] == '\0' ) {
              p->output_type = ULINT;
            }
            else if ( strcmp ( argv[i], "r32" ) == 0 && argv[i][3] == '\0' ) {
              p->output_type = FLOAT;
            }
            else if ( strcmp ( argv[i], "r64" ) == 0 && argv[i][3] == '\0' ) {
              p->output_type = DOUBLE;
            }
            else {
              API_ErrorParse_tree( (char*)NULL, "parsing -type...\n", 0 );
            }
          }

          /* general parameters
           */
          else if ( (strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0')
                    || (strcmp ( argv[i], "--help" ) == 0 && argv[i][6] == '\0') ) {
             API_ErrorParse_tree( (char*)NULL, (char*)NULL, 1);
          }
          else if ( (strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0') ) {
             API_ErrorParse_tree( (char*)NULL, (char*)NULL, 0);
          }
          else if ( strcmp ( argv[i], "-verbose" ) == 0
                    || (strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0') ) {
            if ( _n_call_parse_ == 1 ) {
              if ( _verbose_ <= 0 ) _verbose_ = 1;
              else                  _verbose_ ++;
              if ( _VT_VERBOSE_ <= 0 ) _VT_VERBOSE_ = 1;
              else                     _VT_VERBOSE_ ++;
              incrementVerboseInTopologicalComponentImage( );
              incrementVerboseInTopologicalGenericTree( );
              incrementVerboseInTopologicalVoxelTree( );
              incrementVerboseInTopologicalVoxelTreeMain( );
              incrementVerboseInTopologicalVoxelTreeOrder( );
              incrementVerboseInTopologicalVoxelTreePruning( );
            }
          }
          else if ( strcmp ( argv[i], "-no-verbose" ) == 0
                    || strcmp ( argv[i], "-noverbose" ) == 0
                    || (strcmp ( argv[i], "-nv" ) == 0 && argv[i][3] == '\0') ) {
              _verbose_ = 0;
              _VT_VERBOSE_ = 0;
              setVerboseInTopologicalComponentImage( 0 );
              setVerboseInTopologicalGenericTree( 0 );
              setVerboseInTopologicalVoxelTree( 0 );
              setVerboseInTopologicalVoxelTreeMain( 0 );
              setVerboseInTopologicalVoxelTreeOrder( 0 );
              setVerboseInTopologicalVoxelTreePruning( 0 );
          }
          else if ( (strcmp ( argv[i], "-debug" ) == 0 && argv[i][6] == '\0')
                    || (strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0') ) {
            if ( _n_call_parse_ == 1 ) {
              if ( _debug_ <= 0 ) _debug_ = 1;
              else                _debug_ ++;
              if ( _VT_DEBUG_ <= 0 ) _VT_DEBUG_ = 1;
              else                   _VT_DEBUG_ ++;
              incrementDebugInTopologicalComponentImage( );
              incrementDebugInTopologicalGenericTree( );
              incrementDebugInTopologicalVoxelTree( );
              incrementDebugInTopologicalVoxelTreeMain( );
              incrementDebugInTopologicalVoxelTreeOrder( );
              incrementDebugInTopologicalVoxelTreePruning( );
            }
          }
          else if ( (strcmp ( argv[i], "-no-debug" ) == 0 && argv[i][9] == '\0')
                    || (strcmp ( argv[i], "-nodebug" ) == 0 && argv[i][8] == '\0') ) {
              _debug_ = 0;
              _VT_DEBUG_ = 0;
              setDebugInTopologicalComponentImage( 0 );
              setDebugInTopologicalGenericTree( 0 );
              setDebugInTopologicalVoxelTree( 0 );
              setDebugInTopologicalVoxelTreeMain( 0 );
              setDebugInTopologicalVoxelTreeOrder( 0 );
              setDebugInTopologicalVoxelTreePruning( 0 );
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

          /* unknown option
           */
          else {
              sprintf(text,"unknown option %s\n",argv[i]);
              API_ErrorParse_tree( (char*)NULL, text, 0);
          }
      }

      /* strings beginning with a character different from '-'
       */
      else {
          if ( strlen( argv[i] ) >= STRINGLENGTH ) {
              fprintf( stderr, "... parsing '%s'\n", argv[i] );
              API_ErrorParse_tree( (char*)NULL, "too long file name ...\n", 0 );
          }
          else if ( inputisread == 0 ) {
              (void)strcpy( p->input_name, argv[i] );
              inputisread = 1;
          }
          else if ( outputisread == 0 ) {
              (void)strcpy( p->output_name, argv[i] );
              outputisread = 1;
          }
          else {
              fprintf( stderr, "... parsing '%s'\n", argv[i] );
              API_ErrorParse_tree( (char*)NULL, "too many file names, parsing '-'...\n", 0 );
          }
      }
  }

  /* output image type
   */
  if ( (o != 0) || (s != 0) || (r != 0) ) {
    if ( (o == 1) && (s == 1) && (r == 0) ) p->output_type = SCHAR;
    else if ( (o == 1) && (s == 0) && (r == 0) ) p->output_type = UCHAR;
    else if ( (o == 2) && (s == 0) && (r == 0) ) p->output_type = USHORT;
    else if ( (o == 2) && (s == 1) && (r == 0) ) p->output_type = SSHORT;
    else if ( (o == 4) && (s == 1) && (r == 0) ) p->output_type = SINT;
    else if ( (o == 0) && (s == 0) && (r == 1) ) p->output_type = FLOAT;
    else {
        API_ErrorParse_tree( (char*)NULL, "unable to determine output image type ...\n", 0 );
    }
  }

}
