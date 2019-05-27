/*************************************************************************
 * api-adhocFuse.h -
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

#ifndef _api_adhocfuse_h_
#define _api_adhocfuse_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <pixel-operation.h>

#include <vt_typedefs.h>
#include <vt_image.h>

#include <vt_adhocFuse.h>







typedef struct lineCmdParamAdhocFuse {

  /* image names and output type
   */
  char input_intensity_name[STRINGLENGTH];
  char input_reconstructed_name[STRINGLENGTH];
  char input_minimum_mask_name[STRINGLENGTH];
  char input_maximum_mask_name[STRINGLENGTH];
  char input_segmentation_name[STRINGLENGTH];

  char output_minimum_name[STRINGLENGTH];
  char output_maximum_name[STRINGLENGTH];
  char output_intensity_name[STRINGLENGTH];
  char output_fusion_name[STRINGLENGTH];

  int input_inv;
  int input_swap;
  ImageType output_type;

  /* specific arguments
   */
  float percentileMin;
  float percentileMax;
  enumExtremaMethod methodMin;
  enumExtremaMethod methodMax;

  float sigma;

  /* ... */

  /* general parameters
   */
  int allow_stdin_stdout;
  int print_lineCmdParam;
  int print_time;
  int trace_allocations;

} lineCmdParamAdhocFuse;



extern int API_adhocFuse( vt_image *intImage, vt_image *recImage,
                          vt_image *minMskImage, vt_image *maxMskImage,
                          vt_image *segImage,
                          vt_image *resMinImage, vt_image *resMaxImage,
                          vt_image *resIntmage, vt_image *recFusImage,
                          char *param_str_1,
                          char *param_str_2 );



extern char *API_Help_adhocFuse( int h );

extern void API_ErrorParse_adhocFuse( char *program, char *str, int flag );

extern void API_InitParam_adhocFuse( lineCmdParamAdhocFuse *par );

extern void API_PrintParam_adhocFuse( FILE *theFile, char *program,
                                         lineCmdParamAdhocFuse *par,
                                         char *str );

extern void API_ParseParam_adhocFuse( int firstargc, int argc, char *argv[],
                                 lineCmdParamAdhocFuse *p );



#ifdef __cplusplus
}
#endif

#endif
