/*************************************************************************
 * multipleTrsfs.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mar 17 dec 2013 22:44:32 CET
 *
 * ADDITIONS, CHANGES
 *
 */


#include <sys/time.h> /* gettimeofday() */
#include <time.h> /* clock() */
#include <string.h>

#include <string-tools.h>

#include <bal-transformation-compose.h>
#include <bal-transformation-inversion.h>
#include <bal-transformation-averaging.h>
#include <bal-transformation-propagation.h>


static int _verbose_ = 1;
static int _debug_ = 0;



typedef enum {
  _AVERAGE_,
  _PROPAGATION_
} typeMultipleTrsfsMethod;

typedef enum {
  _MEMORY_,
  _STREAMING_
} typeMultipleTrsfsComputation;


typedef struct local_parameter {

  char *thetransformation_format;
  int firstindex;
  int lastindex;

  int referenceindex;

  int nfloatingbefore;
  int nfloatingafter;

  char *restransformation_format;
  enumTypeTransfo transformation_type;

  typeMultipleTrsfsMethod method;

  int max_iterations;

  int use_inverse;

  enumTransformationAveraging operation;
  bal_estimator estimator;
  double sigma;

  typeMultipleTrsfsComputation computation;
  int print_time;

} local_parameter;





/*------- Definition des fonctions statiques ----------*/
static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p );
static void _InitParam( local_parameter *par );
static double _GetTime();
static double _GetClock();




static char *program = NULL;

static char *usage = "[format-in] -f[irst] %d -l[ast] %d]\n\
 [-ref|-reference %d]\n\
 [-nfb|-nfloatingbefore %d] [-nfa|-nfloatingafter %d]\n\
 [[-res] format-out]\n\
 [-transformation-type|-transformation|-trsf-type %s]\n\
 [-method average|propagation]\n\
 [-max-iteration|-max-iter|-iteration|-i %d]\n\
 [-use-inverse|-inverse]\n\
 [-smoothing-sigma|-sigma %lf]\n\
 [-mean|-robust-mean]\n\
 [-estimator-type|-estimator|-es-type wlts|lts|wls|ls]\n\
 [-lts-cut|-lts-fraction %lf] [-lts-deviation %f] [-lts-iterations %d]\n\
 [-fluid-sigma|-lts-sigma %lf %lf %lf]\n\
 [-lts-cut|-lts-fraction %lf]\n\
 [-streaming | -memory]\n\
 [-v] [-help]";

static char *detail = "\
\n\
given the transformations between couple of images (typically successive\n\
images in a series), compute the transformations for every images with\n\
respect to the same reference image\n\
\n\
[format-in]    # format 'a la printf' of transformations to be processed\n\
               # must contain two '%d'\n\
               # the first one if for the floating image\n\
               # the second one if for the reference image\n\
[-first %d]    # first value of the index in the format\n\
[-last %d]     # last value of the index in the format\n\
[-reference %d] # index of the reference image for result transformations\n\
[-nfloatingbefore %d] # relative left half-interval for floating images\n\
[-nfloatingafter %d]  # relative right half-interval for floating images\n\
   # if 'before' or 'after' are non-null, the interval of transformation\n\
   # for a given reference 'ref' is [ref-before, ref-1] U [ref+1, ref+after]\n\
[-res] format-out   # format 'a la printf' for output transformations\n\
[-transformation-type|-transformation|-trsf-type %s] # transformation type\n\
  translation2D, translation3D, translation-scaling2D, translation-scaling3D,\n\
  rigid2D, rigid3D, rigid, similitude2D, similitude3D, similitude,\n\
  affine2D, affine3D, affine, vectorfield2D, vectorfield3D, vectorfield, vector\n\
[-method average|propagation] # \n\
  average: compute transformation estimation by averaging composition\n\
           of transformation\n\
  propagation: compose transformations from the reference one\n\
[-max-iteration|-max-iter|-iteration|-i %d] # number of iterations\n\
[-use-inverse|-inverse] # use also inverse transformations\n\
[-smoothing-sigma|-sigma %lf] # std deviation for transformation smoothing\n\
[-mean]        # equivalent to '-estimator-type ls'\n\
[-robust-mean] # equivalent to '-estimator-type lts'\n\
[-estimator-type|-estimator|-es-type %s] # transformation estimator\n\
  wlts: weighted least trimmed squares\n\
  lts: least trimmed squares\n\
  wls: weighted least squares\n\
  ls: least squares\n\
  weighted methods are not implemented yet\n\
[-lts-cut|-lts-fraction %lf] # for trimmed estimations (robust mean),\n\
  fraction of pairs that are kept\n\
[-lts-deviation %lf] # for trimmed estimations, defines the threshold to discard\n\
  pairings, ie 'average + this_value * standard_deviation'\n\
[-lts-iterations %d] # for trimmed estimations, the maximal number of iterations\n\
[-fluid-sigma|-lts-sigma %lf %lf %lf] # sigma for fluid regularization,\n\
  ie field interpolation and regularization for pairings (only for vector field)\n\
[-streaming] # computation is done by reading one transformation after the other\n\
[-memory]    # computation is done by loading all transformation in memory\n\
 -v : mode verbose\n\
\n";







int main( int argc, char *argv[] )
{
  local_parameter p;

  stringArray transformationFileArray;
  bal_transformationArray transformationArray;
  int ntransformations = 0;
  bal_transformation invTrsf;

  int i, j;
  char imname[STRINGLENGTH];

  stringList resTransformationFileList;
  bal_transformationList resTransformationList; 

  double time_init = _GetTime();
  double time_exit;
  double clock_init = _GetClock();
  double clock_exit;


  /***************************************************
   *
   * parsing parameters
   *
   ***************************************************/
  program = argv[0];

  /* no arguments
   */
  if ( argc == 1 ) _ErrorParse( NULL, 0 );


  /* parsing parameters 
   */
  _InitParam( &p );
  _Parse( argc, argv, &p );
  

  
  ntransformations = p.lastindex - p.firstindex + 1;
  if ( p.referenceindex < p.firstindex || p.lastindex < p.referenceindex )
    p.referenceindex = p.firstindex;



  /* transformation file names
     string[j][i] is printf( format, i, j )
     since transformations are denoted floating_on_reference and allows to resample
     the floating image onto the reference
     string[j][i] refers to the transformation T_{i<-j}
   */
  initStringArray( &transformationFileArray );

  if ( p.nfloatingbefore <= 0 && p.nfloatingafter <= 0 ) {
      /* build all names
       */
      if ( _debug_ >= 3 ) {
          fprintf( stderr, "%s: build all names\n", argv[0] );
      }
      if ( buildStringArrayFromFormat( p.thetransformation_format,
                                       p.firstindex, p.lastindex,
                                       p.firstindex, p.lastindex,
                                       &transformationFileArray ) != 1 ) {
        _ErrorParse( "unable to build input transformation list from format\n", 0);
      }
  }
  else {
      if ( _debug_ >= 3 ) {
          fprintf( stderr, "%s: build selection of names\n", argv[0] );
      }
      if ( allocStringArray( &transformationFileArray,
                             p.lastindex-p.firstindex+1, p.lastindex-p.firstindex+1 ) != 1 ) {
          _ErrorParse( "unable to allocate input transformation list\n", 0);
      }

      for ( j=p.firstindex; j<=p.lastindex; j++ ) {

          if ( p.nfloatingbefore >= 1 ) {
              for ( i=1; i<=p.nfloatingbefore; i++ ) {
                  if ( j-i < p.firstindex || j-i > p.lastindex ) continue;
                  sprintf( imname, p.thetransformation_format, j-i, j );
                  transformationFileArray.array[j-p.firstindex][j-i-p.firstindex] = (char*)malloc( (strlen(imname)+1) * sizeof(char) );
                  if ( transformationFileArray.array[j-p.firstindex][j-i-p.firstindex] == (char*)NULL ) {
                      freeStringArray( &transformationFileArray );
                      if ( _verbose_ )
                          fprintf( stderr, "... unable to allocate string[%d][%d]\n", j, j-i );
                      _ErrorParse( "unable to allocate input transformation list\n", 0);
                  }
                  (void)strncpy( transformationFileArray.array[j-p.firstindex][j-i-p.firstindex], imname, strlen(imname)+1 );
              }
          }

          sprintf( imname, p.thetransformation_format, j, j );
          transformationFileArray.array[j-p.firstindex][j-p.firstindex] = (char*)malloc( (strlen(imname)+1) * sizeof(char) );
          if ( transformationFileArray.array[j-p.firstindex][j-p.firstindex] == (char*)NULL ) {
              freeStringArray( &transformationFileArray );
              if ( _verbose_ )
                  fprintf( stderr, "... unable to allocate string[%d][%d]\n", j, j );
              _ErrorParse( "unable to allocate input transformation list\n", 0);
          }
          (void)strncpy( transformationFileArray.array[j-p.firstindex][j-p.firstindex], imname, strlen(imname)+1 );

          if ( p.nfloatingafter >= 1 ) {
              for ( i=1; i<=p.nfloatingafter; i++ ) {
                  if ( j+i < p.firstindex || j+i > p.lastindex ) continue;
                  sprintf( imname, p.thetransformation_format, j+i, j );
                  transformationFileArray.array[j-p.firstindex][j+i-p.firstindex] = (char*)malloc( (strlen(imname)+1) * sizeof(char) );
                  if ( transformationFileArray.array[j-p.firstindex][j+i-p.firstindex] == (char*)NULL ) {
                      freeStringArray( &transformationFileArray );
                      if ( _verbose_ )
                          fprintf( stderr, "... unable to allocate string[%d][%d]\n", j, j+i );
                      _ErrorParse( "unable to allocate input transformation list\n", 0);
                  }
                  (void)strncpy( transformationFileArray.array[j-p.firstindex][j+i-p.firstindex], imname, strlen(imname)+1 );
              }
          }
      }

  }

  
  if ( _debug_ >= 4 ) {
    printStringArray( stderr, &transformationFileArray, "Input transformations" );
  }




  /* reading transformations
   */
  BAL_InitTransformationArray( &transformationArray );

  if ( BAL_AllocTransformationArray( &transformationArray, 
                                     ntransformations, ntransformations ) != 1 ) {
    freeStringArray( &transformationFileArray );
    _ErrorParse( "error when allocating input transformations\n", 0);
  }

  
  if ( BAL_ReadTransformationArray( &transformationArray, &transformationFileArray ) != 1 ) {
    BAL_FreeTransformationArray( &transformationArray );
    freeStringArray( &transformationFileArray );
    _ErrorParse( "error when reading input transformations\n", 0);
  }

  if ( _debug_ >= 3 )
    BAL_TestTransformationArray( &transformationArray, &transformationFileArray );


  freeStringArray( &transformationFileArray );


  
  

  /* result transformation file names
   */

  initStringList( &resTransformationFileList );

  if ( p.restransformation_format != (char*)NULL ) {
    if ( buildStringListFromFormat( p.restransformation_format, p.firstindex, p.lastindex, 
                                    &resTransformationFileList ) != 1 ) {
      BAL_FreeTransformationArray( &transformationArray );
      _ErrorParse( "unable to build input transformation list from format\n", 0);
    }
  }
  else {
    BAL_FreeTransformationArray( &transformationArray );
    _ErrorParse( "no format for the result transformations\n", 0);
  }

  
  if ( _debug_ >= 4 )
    printStringList( stderr, &resTransformationFileList, "Result transformations" );






  /* initializing result transformations
     should we check whether read transformations are compatible with
     given type ?
  */
  BAL_InitTransformationList( &resTransformationList );

  switch ( p.transformation_type ) {
  default :
    freeStringList( &resTransformationFileList );
    BAL_FreeTransformationArray( &transformationArray );
    _ErrorParse( "such transformation type not handled yet\n", 0);
    return( 1 );
    
  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :
    
    if ( BAL_FullAllocTransformationList( &resTransformationList, 
                                          ntransformations,
                                          p.transformation_type, (bal_image *)NULL ) != 1 ) {
      freeStringList( &resTransformationFileList );
      BAL_FreeTransformationArray( &transformationArray );
      _ErrorParse( "unable to allocate result transformations (matrices)\n", 0);
    }
    break;
    
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
    
    /* should provide a template
     */
     freeStringList( &resTransformationFileList );
     BAL_FreeTransformationArray( &transformationArray );
     _ErrorParse( "such transformation type not handled yet (vector field)\n", 0);
  }



  switch ( p.method ) {
  default :
    BAL_FreeTransformationArray( &transformationArray );
    BAL_FreeTransformationList( &resTransformationList );
    freeStringList( &resTransformationFileList );
    _ErrorParse( "unknowm computation method\n", 0);
    break;
  case _AVERAGE_ :
    if ( BAL_EstimateTransformationByAveraging( &transformationArray,
                                             &resTransformationList,
                                             &(p.estimator),
                                             p.sigma,
                                             p.max_iterations,
                                             p.use_inverse ) != 1 ) {
      BAL_FreeTransformationArray( &transformationArray );
      BAL_FreeTransformationList( &resTransformationList );
      freeStringList( &resTransformationFileList );
      _ErrorParse( "unable to compute result transformations\n", 0);
    }
    break;
  case _PROPAGATION_ :
    if ( BAL_EstimateTransformationByPropagation( &transformationArray,
                                                  &resTransformationList,
                                                  p.referenceindex-p.firstindex ) != 1 ) {
      BAL_FreeTransformationArray( &transformationArray );
      BAL_FreeTransformationList( &resTransformationList );
      freeStringList( &resTransformationFileList );
      _ErrorParse( "not implemented yet\n", 0);
      _ErrorParse( "unable to compute result transformations\n", 0);
    }
    break;
  }






  BAL_FreeTransformationArray( &transformationArray );


  /* we've got T_{i<-r} transformation, with an unknown reference r
   * for the _AVERAGE_ case
   * we will then compute T_{i<-r} o T_{ref<-r}^{-1}
   */
  if ( p.firstindex <= p.referenceindex && p.referenceindex <= p.lastindex ) {
    switch ( p.method ) {
    default :
      break;
    case _AVERAGE_ :
      BAL_InitTransformation( &invTrsf );
      if ( BAL_AllocTransformation( &invTrsf,
                                    resTransformationList.data[p.referenceindex - p.firstindex].type,
                                    &(resTransformationList.data[p.referenceindex - p.firstindex].vx) ) != 1) {
        BAL_FreeTransformationList( &resTransformationList );
        freeStringList( &resTransformationFileList );
        _ErrorParse( "unable to allocate inverse reference result transformations\n", 0);
      }
      if ( BAL_InverseTransformation( &(resTransformationList.data[p.referenceindex - p.firstindex]),
           &invTrsf ) != 1 ) {
        BAL_FreeTransformation( &invTrsf );
        BAL_FreeTransformationList( &resTransformationList );
        freeStringList( &resTransformationFileList );
        _ErrorParse( "unable to invert reference result transformations\n", 0);
      }
      for ( i=p.firstindex; i<=p.lastindex; i++ ) {
        if ( i == p.referenceindex ) {
          BAL_SetTransformationToIdentity( &(resTransformationList.data[p.referenceindex - p.firstindex]) );
        }
        else {
          if ( BAL_TransformationComposition( &(resTransformationList.data[i - p.firstindex]),
                                              &(resTransformationList.data[i - p.firstindex]),
                                              &invTrsf ) != 1 ) {
            BAL_FreeTransformation( &invTrsf );
            BAL_FreeTransformationList( &resTransformationList );
            freeStringList( &resTransformationFileList );
            _ErrorParse( "unable to compose transformations\n", 0);
          }
        }
      }
      BAL_FreeTransformation( &invTrsf );
      break;
    }
  }

  if ( BAL_WriteTransformationList( &resTransformationList, &resTransformationFileList ) != 1 ) {
    BAL_FreeTransformationList( &resTransformationList );
    freeStringList( &resTransformationFileList );
    _ErrorParse( "unable to write result transformations\n", 0);
  }



  BAL_FreeTransformationList( &resTransformationList );
  freeStringList( &resTransformationFileList );









  

  
  time_exit = _GetTime();
  clock_exit = _GetClock();

  if (  p.print_time ) {
    fprintf( stderr, "%s: elapsed time = %f\n", program, time_exit - time_init );
    fprintf( stderr, "%s: elapsed time = %f\n", program, clock_exit - clock_init );
  }

  return( 0 );
}








static void _Parse( int argc, char *argv[], local_parameter *p )
{
  int i, status;
  int inputisread = 0;
  int outputisread = 0;

  program = argv[0];

  for ( i=1; i<argc; i++ ) {
  
    if ( argv[i][0] == '-' ) {

      /* file names
       */

      if ( strcmp ( argv[i], "-format" ) == 0 ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -format...\n", 0 );
        if ( p->thetransformation_format != (char*)NULL || inputisread > 0 )
          _ErrorParse( "parsing -format: input has already been parsed ...\n", 0 );
        p->thetransformation_format = argv[i];
        inputisread = 1;
      }

      else if ( (strcmp ( argv[i], "-f" ) == 0 && argv[i][2] == '\0') 
                || (strcmp ( argv[i], "-first" ) == 0 && argv[i][6] == '\0') ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -first ...\n", 0 );
        status = sscanf( argv[i], "%d", &(p->firstindex) );
        if ( status <= 0 ) _ErrorParse( "parsing -first ...", 0 );
      }
      else if ( (strcmp ( argv[i], "-l" ) == 0 && argv[i][2] == '\0') 
                || (strcmp ( argv[i], "-last" ) == 0 && argv[i][5] == '\0') ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -last ...\n", 0 );
        status = sscanf( argv[i], "%d", &(p->lastindex) );
        if ( status <= 0 ) _ErrorParse( "parsing -last ...", 0 );
      }

      else if ( (strcmp ( argv[i], "-ref" ) == 0 && argv[i][4] == '\0')
                || (strcmp ( argv[i], "-reference" ) == 0) ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -reference ...\n", 0 );
        status = sscanf( argv[i], "%d", &(p->referenceindex) );
        if ( status <= 0 ) _ErrorParse( "parsing -reference ...", 0 );
      }

      else if ( (strcmp ( argv[i], "-nfb" ) == 0 && argv[i][4] == '\0')
                || (strcmp ( argv[i], "-nfloatingbefore" ) == 0) ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -nfloatingbefore ...\n", 0 );
        status = sscanf( argv[i], "%d", &(p->nfloatingbefore) );
        if ( status <= 0 ) _ErrorParse( "parsing -nfloatingbefore ...", 0 );
      }
      else if ( (strcmp ( argv[i], "-nfa" ) == 0 && argv[i][4] == '\0')
                || (strcmp ( argv[i], "-nfloatingafter" ) == 0) ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -nfloatingafter ...\n", 0 );
        status = sscanf( argv[i], "%d", &(p->nfloatingafter) );
        if ( status <= 0 ) _ErrorParse( "parsing -nfloatingafter ...", 0 );
      }


      else if ( strcmp ( argv[i], "-res" ) == 0  && argv[i][4] == '\0' ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -res...\n", 0 );
        if ( p->restransformation_format != (char*)NULL || outputisread > 0 )
          _ErrorParse( "parsing -res: output has already been parsed ...\n", 0 );
        p->restransformation_format = argv[i];
        outputisread = 1;
      }

      /* transformation type
       */
      else if ( strcmp ( argv[i], "-transformation-type" ) == 0 
                || strcmp ( argv[i], "-transformation" ) == 0
                || strcmp ( argv[i], "-trsf-type" ) == 0 ) {
        i ++;
        if ( i >= argc)    _ErrorParse( "-transformation-type", 0 );
        if ( strcmp ( argv[i], "translation2D" ) == 0 ) {
          p->transformation_type = TRANSLATION_2D;
        }
        else if ( strcmp ( argv[i], "translation3D" ) == 0 ) {
          p->transformation_type = TRANSLATION_3D;
        }
        else if ( strcmp ( argv[i], "translation" ) == 0 && argv[i][11] == '\0') {
          p->transformation_type = TRANSLATION_3D;
        }
        else if ( strcmp ( argv[i], "translation-scaling2D" ) == 0 ) {
          p->transformation_type = TRANSLATION_SCALING_2D;
        }
        else if ( strcmp ( argv[i], "translation-scaling3D" ) == 0 ) {
          p->transformation_type = TRANSLATION_SCALING_3D;
        }
        else if ( strcmp ( argv[i], "rigid2D" ) == 0 ) {
          p->transformation_type = RIGID_2D;
        }
        else if ( strcmp ( argv[i], "rigid3D" ) == 0 ) {
          p->transformation_type = RIGID_3D;
        }
        else if ( (strcmp ( argv[i], "rigid" ) == 0 && argv[i][5] == '\0') ) {
          p->transformation_type = RIGID_3D;
        }
        else if ( strcmp ( argv[i], "similitude2D" ) == 0 ) {
          p->transformation_type = SIMILITUDE_2D;
        }
        else if ( strcmp ( argv[i], "similitude3D" ) == 0 ) {
          p->transformation_type = SIMILITUDE_3D;
        }
        else if ( strcmp ( argv[i], "similitude" ) == 0 ) {
          p->transformation_type = SIMILITUDE_3D;
        }
        else if ( strcmp ( argv[i], "affine2D" ) == 0 ) {
          p->transformation_type = AFFINE_2D;
        }
        else if ( strcmp ( argv[i], "affine3D" ) == 0 ) {
          p->transformation_type = AFFINE_3D;
        }
        else if ( strcmp ( argv[i], "affine" ) == 0 ) {
          p->transformation_type = AFFINE_3D;
        }
        /*
          else if ( strcmp ( argv[i], "spline" ) == 0 ) {
          p->transformation_type = SPLINE;
          }
        */
        else if ( strcmp ( argv[i], "vectorfield" ) == 0
                  || (strcmp ( argv[i], "vector" ) == 0 && argv[i][6] == '\0') ) {
          p->transformation_type = VECTORFIELD_3D;
        }
        else if ( strcmp ( argv[i], "vectorfield3D" ) == 0
                  || strcmp ( argv[i], "vector3D" ) == 0 ) {
          p->transformation_type = VECTORFIELD_3D;
        }
        else if ( strcmp ( argv[i], "vectorfield2D" ) == 0
                  || strcmp ( argv[i], "vector2D" ) == 0 ) {
          p->transformation_type = VECTORFIELD_2D;
        }
        else {
          fprintf( stderr, "unknown transformation type: '%s'\n", argv[i] );
          _ErrorParse( "-transformation-type", 0 );
        }
      }


      /* method
       */
      else if ( strcmp ( argv[i], "-method" ) == 0 ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -method ...\n", 0 );
        if ( strcmp ( argv[i], "average" ) == 0 ) {
          p->method = _AVERAGE_;
        }
        else if ( strcmp ( argv[i], "propagation" ) == 0 ) {
          p->method = _PROPAGATION_;
        }
        else {
          _ErrorParse( "parsing -method: unknown method...\n", 0 );
        }
      }

      /* computation tuning
       */
      else if ( (strcmp ( argv[i], "-i" ) == 0 && argv[i][2] == '\0') 
                || (strcmp ( argv[i], "-iteration" ) == 0 && argv[i][10] == '\0')
                || (strcmp ( argv[i], "-max-iteration" ) == 0)
                || (strcmp ( argv[i], "-max-iter" ) == 0 && argv[i][9] == '\0') ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -max-iteration ...\n", 0 );
        status = sscanf( argv[i], "%d", &(p->max_iterations) );
        if ( status <= 0 ) _ErrorParse( "parsing -max-iteration ...", 0 );
      }

      /*
       */
       else if ( strcmp ( argv[i], "-inverse" ) == 0 
                 || strcmp ( argv[i], "-use-inverse" ) == 0 ) {
        p->use_inverse = 1;
      }

      /* post-estimation smoothing
       */
      else if ( (strcmp ( argv[i], "-sigma" ) == 0 && argv[i][6] == '\0') 
                || strcmp ( argv[i], "-smoothing-sigma" ) == 0 ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -smoothing-sigma ...\n", 0 );
        status = sscanf( argv[i], "%lf", &(p->sigma) );
        if ( status <= 0 ) _ErrorParse( "parsing -smoothing-sigma ...", 0 );
      }

      /* averaging operation
       */

      else if ( strcmp ( argv[i], "-mean" ) == 0 ) {
        p->estimator.type = TYPE_LS;
      }
      else if ( strcmp ( argv[i], "-robust-mean" ) == 0 
                || strcmp ( argv[i], "-rmean" ) == 0 ) {
        p->estimator.type = TYPE_LTS;
      }

      else if ( strcmp ( argv[i], "-estimator-type") == 0 
                || strcmp ( argv[i], "-estimator") == 0
                || strcmp ( argv[i], "-es-type") == 0 ) {
        i ++;
        if ( i >= argc)    _ErrorParse( "-estimator-type", 0 );
        if ( (strcmp ( argv[i], "ltsw" ) == 0 && argv[i][4] == '\0')
             || (strcmp ( argv[i], "wlts" ) == 0 && argv[i][4] == '\0') ) {
          p->estimator.type = TYPE_WLTS;
        }
        else if ( strcmp ( argv[i], "lts" ) == 0 && argv[i][3] == '\0' ) {
          p->estimator.type = TYPE_LTS;
        }
        else if ( (strcmp ( argv[i], "lsw" ) == 0 && argv[i][3] == '\0')
                  || (strcmp ( argv[i], "wls" ) == 0 && argv[i][3] == '\0') ) {
          p->estimator.type = TYPE_WLS;
        }
        else if ( strcmp ( argv[i], "ls" ) == 0 && argv[i][2] == '\0' ) {
          p->estimator.type = TYPE_LS;
        }
        else {
          fprintf( stderr, "unknown estimator type: '%s'\n", argv[i] );
          _ErrorParse( "-estimator-type", 0 );
        }
      }

      else if ( strcmp ( argv[i], "-lts-fraction" ) == 0 
                || strcmp ( argv[i], "-lts-cut" ) == 0) {
        i ++;
        if ( i >= argc)    _ErrorParse( "-lts-fraction", 0 );
        status = sscanf( argv[i], "%lf", &(p->estimator.retained_fraction) );
        if ( status <= 0 ) _ErrorParse( "-lts-fraction", 0 );
      }

      else if ( strcmp ( argv[i], "-lts-deviation" ) == 0 ) {
        i ++;
        if ( i >= argc)    _ErrorParse( "-lts-deviation", 0 );
        status = sscanf( argv[i], "%lf", &(p->estimator.standard_deviation_threshold) );
        if ( status <= 0 ) _ErrorParse( "-lts-deviation", 0 );
      }

      else if ( strcmp ( argv[i], "-lts-iterations" ) == 0 ) {
        i ++;
        if ( i >= argc)    _ErrorParse( "-lts-iterations", 0 );
        status = sscanf( argv[i], "%d", &(p->estimator.max_iterations) );
        if ( status <= 0 ) _ErrorParse( "-lts-iterations", 0 );
      }
    
      else if ( (strcmp (argv[i], "-fluid-sigma" ) == 0 && argv[i][12] == '\0')
                || (strcmp (argv[i], "-lts-sigma" ) == 0 && argv[i][10] == '\0') ) {
        i ++;
        if ( i >= argc)    _ErrorParse( "parsing -lts-sigma %lf", 0 );
        status = sscanf( argv[i], "%lf", &(p->estimator.sigma.x) );
        if ( status <= 0 ) _ErrorParse( "parsing -lts-sigma %lf", 0 );
        i ++;
        if ( i >= argc) {
          p->estimator.sigma.y = p->estimator.sigma.x;
          p->estimator.sigma.z = p->estimator.sigma.x;
        }
        else {
          status = sscanf( argv[i], "%lf", &(p->estimator.sigma.y) );
          if ( status <= 0 ) {
            i--;
            p->estimator.sigma.y = p->estimator.sigma.x;
            p->estimator.sigma.z = p->estimator.sigma.x;
          }
          else {
            i ++;
            if ( i >= argc) p->estimator.sigma.z = 0;
            else {
              status = sscanf( argv[i], "%lf", &(p->estimator.sigma.z) );
              if ( status <= 0 ) {
                i--;
                p->estimator.sigma.z = 0;
              }
            }
          }
        }
      }


      /* kind of computation
       */
      else if ( (strcmp ( argv[i], "-memory" ) == 0 && argv[i][7] == '\0') ) {
        p->computation = _MEMORY_;
      }
      else if ( (strcmp ( argv[i], "-streaming" ) == 0 && argv[i][10] == '\0') ) {
        p->computation = _STREAMING_;
      }

      /* general arguments*/
       else if ( (strcmp ( argv[i], "-time" ) == 0 && argv[i][5] == '\0') ) {
        p->print_time = 1;
      }
      else if ( (strcmp ( argv[i], "-notime" ) == 0 && argv[i][7] == '\0')  
                || (strcmp ( argv[i], "-no-time" ) == 0 && argv[i][8] == '\0') ) {
        p->print_time = 0;
      }

      else if ( (strcmp ( argv[i], "-verbose") == 0 && argv[i][8] == '\0')
                || (strcmp ( argv[i], "-v") == 0 && argv[i][2] == '\0') ) {
        if ( _verbose_ <= 0 ) _verbose_ = 1;
        else _verbose_ ++;
        BAL_IncrementVerboseInBalTransformationAveraging( );
        BAL_IncrementVerboseInBalTransformationPropagation( );
        BAL_IncrementVerboseInBalTransformation(  );
        BAL_IncrementVerboseInBalTransformationInversion( );
      }
      
      else if ( strcmp ( argv[i], "--help" ) == 0 
                || ( strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0' )
                || ( strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0' )
                || ( strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0' ) ) {
        _ErrorParse( NULL, 1 );
      }
      
      else if ( (strcmp ( argv[i], "-debug" ) == 0 && argv[i][6] == '\0')
                || (strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0') ) {
        if ( _debug_ >= 1 )
            _debug_ ++;
        else
            _debug_ = 1;
        BAL_IncrementDebugInBalTransformationPropagation( );
        BAL_IncrementDebugInBalTransformationAveraging(  );
      }
      else if ( (strcmp ( argv[i], "-no-debug" ) == 0 && argv[i][9] == '\0')
                || (strcmp ( argv[i], "-nodebug" ) == 0 && argv[i][8] == '\0') ) {
          _debug_ = 0;
          BAL_IncrementDebugInBalTransformationPropagation( 0 );
          BAL_IncrementDebugInBalTransformationAveraging( 0 );
      }

      /* unknown option
       */
      else {
        fprintf(stderr,"unknown option: '%s'\n",argv[i]);
      }
    }
    
    /*--- saisie des noms d'images ---*/
    else if ( argv[i][0] != 0 ) {
      if ( inputisread == 0 ) {
        p->thetransformation_format = argv[i];
        inputisread = 1;
      }
      else if ( outputisread == 0 ) {
        p->restransformation_format = argv[i];
        outputisread = 1;
      }
      else 
        fprintf(stderr,"too many file names: '%s'\n",argv[i]);
    }

  }
  
}





static void _ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  if ( str != (char*)NULL )
    (void)fprintf(stderr,"Erreur : %s",str);
  exit( 1 );
}





static void _InitParam( local_parameter *p )
{
  p->thetransformation_format = (char*)NULL;
  p->firstindex = 0;
  p->lastindex = 0;

  p->referenceindex = -1;

  p->nfloatingbefore = -1;
  p->nfloatingafter = -1;

  p->restransformation_format = (char*)NULL;

  p->transformation_type = UNDEF_TRANSFORMATION;


  p->method = _AVERAGE_;

  /* specific parameters for averaging operation
   */
  p->max_iterations = 0;

  p->use_inverse = 0;

  p->operation = _MEAN_;

  BAL_InitEstimator( &(p->estimator) );
  p->estimator.type = TYPE_LS;
  p->estimator.retained_fraction = 0.7500000;
  p->sigma = -1.0;

  p->computation = _MEMORY_;

  p->print_time = 0;
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

