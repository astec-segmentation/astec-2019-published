/*************************************************************************
 * api-membrane.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2016, all rights reserved
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * Ven 2 dec 2016 13:36:43 CEST
 *
 * ADDITIONS, CHANGES
 *
 */

#ifndef _api_Membrane_h_
#define _api_Membrane_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <vt_typedefs.h>
#include <vt_image.h>


#include <vt_common.h>
#include <mt_membrane3D.h>
#include <connexe.h>
#include <time.h>

#include <chunks.h>

typedef enum {
  MULTI_SCALE,
  SINGLE_SCALE
} enumComputation;



typedef struct lineCmdParamMembrane {

  /* image names and output type
   */
  /*char input_name[STRINGLENGTH];
  char output_name[STRINGLENGTH];
  int input_inv;
  int input_swap;
  ImageType output_type;
*/
  /* specific arguments
   */

  /* ... */

  /* general parameters
   */
  /*
  int allow_stdin_stdout;
  */
  int print_lineCmdParam;
  int print_time;
  int trace_allocations;


  vt_names names;
  int type;

  int writeImages;
  enumComputation typeComputation;
  enumMode mode;

  enumStructureColor structureColor;

  double scale1;
  double scale2;
  int nbscales;

  double zfact;

  int dimension;
  int hsp;

  int flagMask;

} lineCmdParamMembrane;



extern int API_membrane( vt_image *image,
                             vt_3Dimres *imsRes,
                             vt_image *mask,
                             char *param_str_1,
                             char *param_str_2 );

extern int API_extrema(vt_3Dimres *imsRes,
                       vt_image *mask,
                       int flagMask,
                       char *param_str_1,
                       char *param_str_2 );


extern char *API_Help_membrane( int h );

extern void API_ErrorParse_membrane( char *program, char *str, int flag );

extern void API_InitParam_membrane( lineCmdParamMembrane *par );

/*extern void API_PrintParam_membrane( FILE *theFile, char *program,
                                         lineCmdParamMembrane *par,
                                         char *str );
*/
extern void API_ParseParam_membrane( int firstargc, int argc, char *argv[],
                                 lineCmdParamMembrane *p );



#ifdef __cplusplus
}
#endif

#endif
