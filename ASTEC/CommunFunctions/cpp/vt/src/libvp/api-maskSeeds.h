/*************************************************************************
 * api-maskSeeds.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mer 25 jul 2018 17:43:47 CEST
 *
 * ADDITIONS, CHANGES
 *
 */

#ifndef _api_maskseeds_h_
#define _api_maskseeds_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <morphotools.h>

#include <vt_typedefs.h>
#include <vt_neighborhood.h>
#include <vt_image.h>



typedef struct lineCmdParamMaskSeeds {

  /* image names and output type
   */
  char input_seed_name[STRINGLENGTH];
  char input_cell_name[STRINGLENGTH];
  char output_name[STRINGLENGTH];



  /* general parameters
   */
  int print_lineCmdParam;
  int print_time;

} lineCmdParamMaskSeeds;



extern int API_maskSeeds( vt_image *imseed, vt_image *imcell, vt_image *imres,
                           char *param_str_1, char *param_str_2 );

extern char *API_Help_maskSeeds( int h );

extern void API_ErrorParse_maskSeeds( char *program, char *str, int flag );

extern void API_InitParam_maskSeeds( lineCmdParamMaskSeeds *par );

extern void API_PrintParam_maskSeeds( FILE *theFile, char *program,
                                         lineCmdParamMaskSeeds *par,
                                         char *str );

extern void API_ParseParam_maskSeeds( int firstargc, int argc, char *argv[],
                                 lineCmdParamMaskSeeds *p );



#ifdef __cplusplus
}
#endif

#endif
