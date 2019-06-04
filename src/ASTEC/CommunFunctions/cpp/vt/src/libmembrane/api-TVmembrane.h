/*************************************************************************
 * api-TVmembrane.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2016, all rights reserved
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * Ven 2 dec 2016 17:22:42 CEST
 *
 * ADDITIONS, CHANGES
 *
  */

#ifndef _api_TVmembrane_h_
#define _api_TVmembrane_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <vt_typedefs.h>
#include <vt_image.h>


#include <vt_common.h>
#include <mt_membrane3D.h>
#include <parcelling.h>

#include <time.h>
#include <sys/time.h>

#include <chunks.h>

typedef enum enumSampleMode {
  _UNKNOWN_ /* unknown mode */,
  _RANDOM_SAMPLING_,
  _REGULAR_SAMPLING_
} enumSampleMode;


typedef struct lineCmdParamTVmembrane {

  /* image names and output type
   */
  char inputBinary[STRINGLENGTH];
  char inputMask[STRINGLENGTH];

  char outputGeneric[STRINGLENGTH];

  char outputSampling[STRINGLENGTH];

  char outputAngles[STRINGLENGTH];
  char outputEigenvalues[STRINGLENGTH];
  char outputTensor[STRINGLENGTH];
  char outputBinary[STRINGLENGTH];

  char outputExtrema[STRINGLENGTH];

  char outputSuffix[STRINGLENGTH];


  /* sampling options
   */
  double sample;
  int power;
  enumSampleMode sampleMode;
  int itermax;

  /* tensor voting options
   */
  double scale;
  int flagReal;
  double zfact;
  hessianMode initHessian;
  int niter;   /*  pour sparse voting */

  enumTVmode TVmode;
  int nangles; /*  nanglesiter */
  int nsticks; /*  nombre de sticks pour le calcul des plate fields */

  /* shape extraction
   */
  shapeExtraction shape;


  /* general options
   */
  int writeImages;
  int dimension;

  /* general parameters
   */
  int print_lineCmdParam;
  int print_time;
  int trace_allocations;

} lineCmdParamTVmembrane;


extern int API_Sampling(vt_image *imageBin,
                 char *param_str_1,
                 char *param_str_2 );


extern int API_TVmembrane( vt_image **imagesIn,
                           vt_3Dtensor *theTensor,
                           char *param_str_1,
                           char *param_str_2 );

extern int API_TVmembrane2D( vt_image **imagesIn,
                           mt_2Dtensor *theTensor2D,
                           char *param_str_1,
                           char *param_str_2 );

extern int API_ShapeExtraction( vt_3Dtensor *theTensor3D, char *param_str_1, char *param_str_2 );
extern int API_ShapeExtraction2D( mt_2Dtensor *theTensor2D, char *param_str_1, char *param_str_2 );


extern char *API_Help_TVmembrane( int h );

extern void API_ErrorParse_TVmembrane( char *program, char *str, int flag );

extern void API_InitParam_TVmembrane( lineCmdParamTVmembrane *par );

extern void API_PrintParam_TVmembrane( FILE *theFile, char *program,
                                         lineCmdParamTVmembrane *par,
                                         char *str );

extern void API_ParseParam_TVmembrane( int firstargc, int argc, char *argv[],
                                 lineCmdParamTVmembrane *p );



#ifdef __cplusplus
}
#endif

#endif
