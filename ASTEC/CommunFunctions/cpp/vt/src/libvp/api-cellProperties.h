/*************************************************************************
 * api-cellProperties.h -
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

#ifndef _api_cellproperties_h_
#define _api_cellproperties_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <vt_typedefs.h>

#include <bal-image.h>
#include <bal-blockmatching-param.h>

#include <bal-cellProperties.h>





typedef struct lineCmdParamCellProperties {

  /* image names and output type
   */
  char input_segmentation_name[STRINGLENGTH];

  char input_fusion_format[STRINGLENGTH];
  char input_segmentation_format[STRINGLENGTH];
  int firstindex;
  int lastindex;

  char input_xml_name[STRINGLENGTH];
  char output_xml_name[STRINGLENGTH];

  char output_diagnosis_name[STRINGLENGTH];

  /* specific arguments
   */

  int acquisition_time;
  typeUpdateList updateList;

  typePropertyList output;
  enumSurfaceEstimation surfaceEstimationType;

  /* dedicated parallelism parameters
   */
  int chunks_for_properties;

  /* parameters for image registration
   */
  int normalisation;
  bal_blockmatching_pyramidal_param affineRegistration;
  bal_blockmatching_pyramidal_param nonlinearRegistration;


  /* ... */

  /* general parameters
   */
  int allow_stdin_stdout;
  int print_lineCmdParam;
  int print_time;
  int trace_allocations;

} lineCmdParamCellProperties;


extern int API_INTERMEDIARY_sequenceCellProperties( char *theformat_segmentation,
                                                    char *theformat_fusion,
                                                    int firstindex, int lastindex,
                                                    char *output_xml_name,
                                                    char *output_diagnosis_name,
                                                    char *param_str_1, char *param_str_2 );

extern int API_INTERMEDIARY_sequenceCellPropertiesUpdate( char *theformat_segmentation,
                                                          char *theformat_fusion,
                                                          char *intput_xml_name,
                                                          char *output_xml_name,
                                                          char *output_diagnosis_name,
                                                          char *param_str_1, char *param_str_2 );

extern int API_INTERMEDIARY_imageCellProperties( char *theim_name,
                                                 char *output_xml_name,
                                                 char *param_str_1,
                                                 char *param_str_2 );

extern int API_INTERMEDIARY_properties( char *input_xml_name,
                                        char *output_xml_name,
                                        char *output_diagnosis_name,
                                        char *param_str_1, char *param_str_2 );


extern int API_cellImageProperties( bal_image *image,
                             typeCellImage *cellImage,
                             char *param_str_1,
                             char *param_str_2 );



extern char *API_Help_cellProperties( int h );

extern void API_ErrorParse_cellProperties( char *program, char *str, int flag );

extern void API_InitParam_cellProperties( lineCmdParamCellProperties *par );

extern void API_FreeParam_cellProperties( lineCmdParamCellProperties *par );

extern void API_PrintParam_cellProperties( FILE *theFile, char *program,
                                         lineCmdParamCellProperties *par,
                                         char *str );

extern void API_ParseParam_cellProperties( int firstargc, int argc, char *argv[],
                                 lineCmdParamCellProperties *p );



#ifdef __cplusplus
}
#endif

#endif
