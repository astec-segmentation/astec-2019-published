/*************************************************************************
 * api-adhocFuse.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2015, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Ven  2 fév 2018 10:18:50 CET
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <chunks.h>
#include <vtmalloc.h>

#include <vt_common.h>

#include <vt_adhocFuse.h>

#include <api-adhocFuse.h>






static int _verbose_ = 1;
static int _debug_ = 0;


static void _API_ParseParam_adhocFuse( char *str, lineCmdParamAdhocFuse *p );



/************************************************************
 *
 * main API
 *
 ************************************************************/



int API_adhocFuse( vt_image *intImage, vt_image *recImage,
                   vt_image *minMskImage, vt_image *maxMskImage,
                   vt_image *segImage,
                   vt_image *resMinImage, vt_image *resMaxImage,
                   vt_image *resIntImage, vt_image *resFusImage,
                   char *param_str_1, char *param_str_2 )
{
  char *proc = "API_adhocFuse";
  lineCmdParamAdhocFuse par;

  vt_image tmpIntImage;
  vt_image *ptrIntImage = (vt_image*)NULL;

  int globalMin = -1;
  int globalMax = -1;

  int theDim[3];



  /* parameter initialization
   */
  API_InitParam_adhocFuse( &par );

  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_adhocFuse( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_adhocFuse( param_str_2, &par );

  if ( par.print_lineCmdParam )
      API_PrintParam_adhocFuse( stderr, proc, &par, (char*)NULL );



  /************************************************************
   *
   *  here is the stuff
   *
   ************************************************************/



  /* allocation of resampled intensity image (if required)
   */
  VT_Image( &tmpIntImage );
  if ( resIntImage == (vt_image*)NULL || resIntImage->buf == (void*)NULL ) {
    VT_InitFromImage( &tmpIntImage, intImage, (char*)NULL, recImage->type );
    if ( VT_AllocImage( &tmpIntImage ) != 1 ) {
      if ( _verbose_ ) {
        fprintf( stderr, "%s: auxiliary intensity image allocation failed\n", proc );
      }
      return( -1 );
    }
    ptrIntImage = &tmpIntImage;
  }
  else {
    ptrIntImage = resIntImage;
  }



  if ( AdhocFuse_RescaleIntensityImage( intImage, minMskImage, maxMskImage, segImage,
                                        ptrIntImage, resMinImage, resMaxImage,
                                        &globalMin, &globalMax,
                                        par.percentileMin, par.percentileMax,
                                        par.methodMin, par.methodMax,
                                        par.sigma ) != 1 ) {
      if ( ptrIntImage == &tmpIntImage ) VT_FreeImage( &tmpIntImage );
      if ( _verbose_ ) {
        fprintf( stderr, "%s: error when renormalizing image\n", proc );
      }
      return( -1 );
  }


  /* fusion
   */

  theDim[0] = intImage->dim.x * intImage->dim.v;
  theDim[1] = intImage->dim.y;
  theDim[2] = intImage->dim.z;

  if ( resFusImage != (vt_image*)NULL && resFusImage->buf != (void*)NULL ) {
    if ( maxImages( ptrIntImage->buf, ptrIntImage->type, recImage->buf, recImage->type,
                    resFusImage->buf, resFusImage->type, theDim ) != 1 ) {
      if ( ptrIntImage == &tmpIntImage ) VT_FreeImage( &tmpIntImage );
      if ( _verbose_ ) {
        fprintf( stderr, "%s: error when fusing images\n", proc );
      }
      return( -1 );
    }
  }

  /* ... */

  if ( ptrIntImage == &tmpIntImage ) VT_FreeImage( &tmpIntImage );

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



static char *usage = "-intensity-image|-ii %s\n\
 -reconstructed-image|-ri %s\n\
 [-minimum-mask-image|-min-mi %s]\n\
 [-maximum-mask-image|-max-mi %s]\n\
 [-segmentation-image|-si %s]\n\
 [-result-minimum-image|-rmini %s]\n\
 [-result-maximum-image|-rmaxi %s]\n\
 [-result-intensity-image|-rii %s]\n\
 [-result-fused-image|-rfi %s]\n\
 [-min-percentile|-min-p %f]\n\
 [-max-percentile|-max-p %f]\n\
 [-min-method|-method-min global|cell|cellborder|cellinterior|voxel]\n\
 [-max-method|-method-max global|cell|cellborder|cellinterior|voxel]\n\
 [-sigma %f]\n\
 [-parallel|-no-parallel] [-max-chunks %d]\n\
 [-parallelism-type|-parallel-type default|none|openmp|omp|pthread|thread]\n\
 [-omp-scheduling|-omps default|static|dynamic-one|dynamic|guided]\n\
 [-inv] [-swap] [output-image-type | -type s8|u8|s16|u16...]\n\
 [-verbose|-v] [-no-verbose|-noverbose|-nv]\n\
 [-debug|-D] [-no-debug|-nodebug]\n\
 [-allow-pipe|-pipe] [-no-allow-pipe|-no-pipe|-nopipe]\n\
 [-print-parameters|-param]\n\
 [-print-time|-time] [-no-time|-notime]\n\
 [-trace-memory|-memory] [-no-memory|-nomemory]\n\
 [-help|-h]";



static char *detail = "\
  -intensity-image|-ii %s\n\
  -reconstructed-image|-ri %s\n\
  -minimum-mask-image|-min-mi %s\n\
  -maximum-mask-image|-max-mi %s\n\
  -segmentation-image|-si %s\n\
  -result-minimum-image|-rmini %s\n\
  -result-maximum-image|-rmaxi %s\n\
  -result-intensity-image|-rii %s\n\
  -result-fused-image|-rfi %s\n\
# ...\n\
  -min-percentile|-min-p %f\n\
  -max-percentile|-max-p %f\n\
  -min-method|-method-min global|cell|cellborder|cellinterior|voxel\n\
  -max-method|-method-max global|cell|cellborder|cellinterior|voxel\n\
  -sigma %f\n\
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





char *API_Help_adhocFuse( int h )
{
    if ( h == 0 )
        return( usage );
    return( detail );
}





void API_ErrorParse_adhocFuse( char *program, char *str, int flag )
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



void API_InitParam_adhocFuse( lineCmdParamAdhocFuse *p )
{
    (void)strncpy( p->input_intensity_name, "\0", 1 );
    (void)strncpy( p->input_reconstructed_name, "\0", 1 );
    (void)strncpy( p->input_minimum_mask_name, "\0", 1 );
    (void)strncpy( p->input_maximum_mask_name, "\0", 1 );
    (void)strncpy( p->input_segmentation_name, "\0", 1 );

    (void)strncpy( p->output_maximum_name, "\0", 1 );
    (void)strncpy( p->output_minimum_name, "\0", 1 );
    (void)strncpy( p->output_intensity_name, "\0", 1 );
    (void)strncpy( p->output_fusion_name, "\0", 1 );

    p->input_inv = 0;
    p->input_swap = 0;
    p->output_type = TYPE_UNKNOWN;

    p->percentileMin = 0.01;
    p->percentileMax = 0.99;

    p->methodMin = _GLOBAL_;
    p->methodMax = _GLOBAL_;

    p->sigma = 0.0;

    p->allow_stdin_stdout = 1;
    p->print_lineCmdParam = 0;
    p->print_time = 0;
    p->trace_allocations = 0;
}





void API_PrintParam_adhocFuse( FILE *theFile, char *program,
                                  lineCmdParamAdhocFuse *p, char *str )
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

  fprintf( f, "- p->input_intensity_name=" );
  if ( p->input_intensity_name != (char*)NULL && p->input_intensity_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_intensity_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- p->input_reconstructed_name=" );
  if ( p->input_reconstructed_name != (char*)NULL && p->input_reconstructed_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_reconstructed_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- p->input_minimum_mask_name=" );
  if ( p->input_minimum_mask_name != (char*)NULL && p->input_minimum_mask_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_minimum_mask_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- p->input_maximum_mask_name=" );
  if ( p->input_maximum_mask_name != (char*)NULL && p->input_maximum_mask_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_maximum_mask_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- p->input_segmentation_name=" );
  if ( p->input_segmentation_name != (char*)NULL && p->input_segmentation_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_segmentation_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- p->output_minimum_name=" );
  if ( p->output_minimum_name != (char*)NULL && p->output_minimum_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_minimum_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- p->output_maximum_name=" );
  if ( p->output_maximum_name != (char*)NULL && p->output_maximum_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_maximum_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- p->output_intensity_name=" );
  if ( p->output_intensity_name != (char*)NULL && p->output_intensity_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_intensity_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- p->output_fusion_name=" );
  if ( p->output_fusion_name != (char*)NULL && p->output_fusion_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_fusion_name );
  else
    fprintf( f, "'NULL'\n" );


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

  fprintf( f, "# ...\n" );
  fprintf( f, "- p->percentileMin = %f\n", p->percentileMin );
  fprintf( f, "- p->percentileMax = %f\n", p->percentileMax );
  fprintf( f, "- p->methodMin = " );
  switch( p->methodMin ) {
  default : fprintf( f, "unknown\n" ); break;
  case _GLOBAL_ : fprintf( f, "_GLOBAL_\n" ); break;
  case _CELL_ : fprintf( f, "_CELL_\n" ); break;
  case _CELL_BORDER_ : fprintf( f, "_CELL_BORDER_\n" ); break;
  case _CELL_INTERIOR_ : fprintf( f, "_CELL_INTERIOR_\n" ); break;
  case _VOXEL_ : fprintf( f, "_VOXEL_\n" ); break;
  }
  fprintf( f, "- p->methodMax = " );
  switch( p->methodMax ) {
  default : fprintf( f, "unknown\n" ); break;
  case _GLOBAL_ : fprintf( f, "_GLOBAL_\n" ); break;
  case _CELL_ : fprintf( f, "_CELL_\n" ); break;
  case _CELL_BORDER_ : fprintf( f, "_CELL_BORDER_\n" ); break;
  case _CELL_INTERIOR_ : fprintf( f, "_CELL_INTERIOR_\n" ); break;
  case _VOXEL_ : fprintf( f, "_VOXEL_\n" ); break;
  }
  fprintf( f, "- sigma = %f\n", p->sigma );

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



static void _API_ParseParam_adhocFuse( char *str, lineCmdParamAdhocFuse *p )
{
  char *proc = "_API_ParseParam_adhocFuse";
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

  API_ParseParam_adhocFuse( 0, argc, argv, p );

  vtfree( argv );
}





static int _n_call_parse_ = 0;

void API_ParseParam_adhocFuse( int firstargc, int argc, char *argv[],
                                  lineCmdParamAdhocFuse *p )
{
  int i;
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
            API_ErrorParse_adhocFuse( (char*)NULL, "parsing '-' ...\n", 0 );
          }

          /* file names
           */
          else if ( strcmp ( argv[i], "-intensity-image" ) == 0
                    || (strcmp ( argv[i], "-ii" ) == 0 && argv[i][3] == '\0') ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_adhocFuse( (char*)NULL, "parsing -intensity-image ...\n", 0 );
            (void)strcpy( p->input_intensity_name, argv[i] );
          }

          else if ( strcmp ( argv[i], "-reconstructed-image" ) == 0
                    || (strcmp ( argv[i], "-ri" ) == 0 && argv[i][3] == '\0') ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_adhocFuse( (char*)NULL, "parsing -reconstructed-image ...\n", 0 );
            (void)strcpy( p->input_reconstructed_name, argv[i] );
          }

          else if ( strcmp ( argv[i], "-minimum-mask-image" ) == 0
                    || (strcmp ( argv[i], "-min-mi" ) == 0 && argv[i][7] == '\0') ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_adhocFuse( (char*)NULL, "parsing -minimum-mask-image ...\n", 0 );
            (void)strcpy( p->input_minimum_mask_name, argv[i] );
          }

          else if ( strcmp ( argv[i], "-maximum-mask-image" ) == 0
                    || (strcmp ( argv[i], "-max-mi" ) == 0 && argv[i][7] == '\0') ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_adhocFuse( (char*)NULL, "parsing -maximum-mask-image ...\n", 0 );
            (void)strcpy( p->input_maximum_mask_name, argv[i] );
          }

          else if ( strcmp ( argv[i], "-segmentation-image" ) == 0
                    || (strcmp ( argv[i], "-si" ) == 0 && argv[i][3] == '\0') ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_adhocFuse( (char*)NULL, "parsing -segmentation-image ...\n", 0 );
            (void)strcpy( p->input_segmentation_name, argv[i] );
          }

          else if ( strcmp ( argv[i], "-result-minimum-image" ) == 0
                    || (strcmp ( argv[i], "-rmini" ) == 0 && argv[i][6] == '\0') ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_adhocFuse( (char*)NULL, "parsing -result-minimum-image ...\n", 0 );
            (void)strcpy( p->output_minimum_name, argv[i] );
          }

          else if ( strcmp ( argv[i], "-result-maximum-image" ) == 0
                    || (strcmp ( argv[i], "-rmaxi" ) == 0 && argv[i][6] == '\0') ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_adhocFuse( (char*)NULL, "parsing -result-maximum-image ...\n", 0 );
            (void)strcpy( p->output_maximum_name, argv[i] );
          }

          else if ( strcmp ( argv[i], "-result-intensity-image" ) == 0
                    || (strcmp ( argv[i], "-rii" ) == 0 && argv[i][4] == '\0') ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_adhocFuse( (char*)NULL, "parsing -result-intensity-image ...\n", 0 );
            (void)strcpy( p->output_intensity_name, argv[i] );
          }

          else if ( strcmp ( argv[i], "-result-fused-image" ) == 0
                    || (strcmp ( argv[i], "-rfi" ) == 0 && argv[i][4] == '\0') ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_adhocFuse( (char*)NULL, "parsing -result-fused-image ...\n", 0 );
            (void)strcpy( p->output_fusion_name, argv[i] );
          }



          /* ...
           */
          else if ( strcmp ( argv[i], "-min-percentile" ) == 0
                    || (strcmp ( argv[i], "-min-p" ) == 0 && argv[i][6] == '\0') ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_adhocFuse( (char*)NULL, "parsing -min-percentile ...\n", 0 );
             status = sscanf( argv[i], "%f", &(p->percentileMin) );
             if ( status <= 0 ) API_ErrorParse_adhocFuse( (char*)NULL, "parsing -min-percentile ...\n", 0 );
          }
          else if ( strcmp ( argv[i], "-max-percentile" ) == 0
                    || (strcmp ( argv[i], "-max-p" ) == 0 && argv[i][6] == '\0') ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_adhocFuse( (char*)NULL, "parsing -max-percentile ...\n", 0 );
             status = sscanf( argv[i], "%f", &(p->percentileMax) );
             if ( status <= 0 ) API_ErrorParse_adhocFuse( (char*)NULL, "parsing -max-percentile ...\n", 0 );
          }
          else if ( (strcmp ( argv[i], "-min-method" ) == 0 && argv[i][11] == '\0') ||
                    (strcmp ( argv[i], "-method-min" ) == 0 && argv[i][11] == '\0') ) {
              i ++;
              if ( i >= argc)    API_ErrorParse_adhocFuse( (char*)NULL, "parsing -min-method ...\n", 0 );
              if ( strcmp ( argv[i], "global" ) == 0 ) {
                  p->methodMin = _GLOBAL_;
              }
              else if ( strcmp ( argv[i], "cell" ) == 0 && argv[i][4] == '\0' ) {
                  p->methodMin = _CELL_;
              }
              else if ( strcmp ( argv[i], "cellborder" ) == 0 && argv[i][10] == '\0' ) {
                  p->methodMin = _CELL_BORDER_;
              }
              else if ( strcmp ( argv[i], "cellinterior" ) == 0 && argv[i][12] == '\0' ) {
                  p->methodMin = _CELL_INTERIOR_;
              }
              else if ( strcmp ( argv[i], "voxel" ) == 0 && argv[i][5] == '\0' ) {
                  p->methodMin = _VOXEL_;
              }
              else {
                  fprintf( stderr, "unknown method: '%s'\n", argv[i] );
                  API_ErrorParse_adhocFuse( (char*)NULL, "parsing -min-method ...\n", 0 );
              }
          }
          else if ( (strcmp ( argv[i], "-max-method" ) == 0 && argv[i][11] == '\0') ||
                    (strcmp ( argv[i], "-method-max" ) == 0 && argv[i][11] == '\0') ) {
              i ++;
              if ( i >= argc)    API_ErrorParse_adhocFuse( (char*)NULL, "parsing -max-method ...\n", 0 );
              if ( strcmp ( argv[i], "global" ) == 0 ) {
                  p->methodMax = _GLOBAL_;
              }
              else if ( strcmp ( argv[i], "cell" ) == 0 && argv[i][4] == '\0' ) {
                  p->methodMax = _CELL_;
              }
              else if ( strcmp ( argv[i], "cellborder" ) == 0 && argv[i][10] == '\0' ) {
                  p->methodMax = _CELL_BORDER_;
              }
              else if ( strcmp ( argv[i], "cellinterior" ) == 0 && argv[i][12] == '\0' ) {
                  p->methodMax = _CELL_INTERIOR_;
              }
              else if ( strcmp ( argv[i], "voxel" ) == 0 && argv[i][5] == '\0' ) {
                  p->methodMax = _VOXEL_;
              }
              else {
                  fprintf( stderr, "unknown method: '%s'\n", argv[i] );
                  API_ErrorParse_adhocFuse( (char*)NULL, "parsing -max-method ...\n", 0 );
              }
          }

          else if ( strcmp ( argv[i], "-sigma" ) == 0 ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_adhocFuse( (char*)NULL, "parsing -sigma ...\n", 0 );
             status = sscanf( argv[i], "%f", &(p->sigma) );
             if ( status <= 0 ) API_ErrorParse_adhocFuse( (char*)NULL, "parsing -sigma ...\n", 0 );
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
             if ( i >= argc)    API_ErrorParse_adhocFuse( (char*)NULL, "parsing -parallelism-type ...\n", 0 );
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
               API_ErrorParse_adhocFuse( (char*)NULL, "parsing -parallelism-type ...\n", 0 );
             }
          }

          else if ( strcmp ( argv[i], "-max-chunks" ) == 0 ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_adhocFuse( (char*)NULL, "parsing -max-chunks ...\n", 0 );
             status = sscanf( argv[i], "%d", &maxchunks );
             if ( status <= 0 ) API_ErrorParse_adhocFuse( (char*)NULL, "parsing -max-chunks ...\n", 0 );
             if ( maxchunks >= 1 ) setMaxChunks( maxchunks );
          }

          else if ( strcmp ( argv[i], "-omp-scheduling" ) == 0 ||
                   ( strcmp ( argv[i], "-omps" ) == 0 && argv[i][5] == '\0') ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_adhocFuse( (char*)NULL, "parsing -omp-scheduling, no argument\n", 0 );
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
               API_ErrorParse_adhocFuse( (char*)NULL, "parsing -omp-scheduling ...\n", 0 );
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
             if ( i >= argc)    API_ErrorParse_adhocFuse( (char*)NULL, "parsing -o...\n", 0 );
             status = sscanf( argv[i],"%d",&o );
             if ( status <= 0 ) API_ErrorParse_adhocFuse( (char*)NULL, "parsing -o...\n", 0 );
          }
          else if ( strcmp ( argv[i], "-type" ) == 0 && argv[i][5] == '\0' ) {
            i += 1;
            if ( i >= argc)    API_ErrorParse_adhocFuse( (char*)NULL, "parsing -type...\n", 0 );
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
              API_ErrorParse_adhocFuse( (char*)NULL, "parsing -type...\n", 0 );
            }
          }

          /* general parameters
           */
          else if ( (strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0')
                    || (strcmp ( argv[i], "--help" ) == 0 && argv[i][6] == '\0') ) {
             API_ErrorParse_adhocFuse( (char*)NULL, (char*)NULL, 1);
          }
          else if ( (strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0') ) {
             API_ErrorParse_adhocFuse( (char*)NULL, (char*)NULL, 0);
          }
          else if ( strcmp ( argv[i], "-verbose" ) == 0
                    || (strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0') ) {
            if ( _n_call_parse_ == 1 ) {
              if ( _verbose_ <= 0 ) _verbose_ = 1;
              else                  _verbose_ ++;
              if ( _VT_VERBOSE_ <= 0 ) _VT_VERBOSE_ = 1;
              else                     _VT_VERBOSE_ ++;
            }
          }
          else if ( strcmp ( argv[i], "-no-verbose" ) == 0
                    || strcmp ( argv[i], "-noverbose" ) == 0
                    || (strcmp ( argv[i], "-nv" ) == 0 && argv[i][3] == '\0') ) {
              _verbose_ = 0;
              _VT_VERBOSE_ = 0;
          }
          else if ( (strcmp ( argv[i], "-debug" ) == 0 && argv[i][6] == '\0')
                    || (strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0') ) {
            if ( _n_call_parse_ == 1 ) {
              if ( _debug_ <= 0 ) _debug_ = 1;
              else                _debug_ ++;
              if ( _VT_DEBUG_ <= 0 ) _VT_DEBUG_ = 1;
              else                   _VT_DEBUG_ ++;
            }
          }
          else if ( (strcmp ( argv[i], "-no-debug" ) == 0 && argv[i][9] == '\0')
                    || (strcmp ( argv[i], "-nodebug" ) == 0 && argv[i][8] == '\0') ) {
              _debug_ = 0;
              _VT_DEBUG_ = 0;
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
              API_ErrorParse_adhocFuse( (char*)NULL, text, 0);
          }
      }

      /* strings beginning with a character different from '-'
       */
      else {
        fprintf( stderr, "... parsing '%s'\n", argv[i] );
        API_ErrorParse_adhocFuse( (char*)NULL, "unknown option ...\n", 0 );
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
        API_ErrorParse_adhocFuse( (char*)NULL, "unable to determine output image type ...\n", 0 );
    }
  }

}
