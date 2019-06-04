/*************************************************************************
 * meanTrsfs.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Lun 18 nov 2013 17:09:53 CET
 *
 * ADDITIONS, CHANGES
 *
 */


#include <sys/time.h> /* gettimeofday() */
#include <time.h> /* clock() */
#include <string.h>

#include <string-tools.h>

#include <bal-transformation-averaging.h>

static int _verbose_ = 0;
static int _debug_ = 0;





typedef enum {
  _MEMORY_,
  _STREAMING_
} typeComputation;


typedef struct local_parameter {

  char *transformationlist_name;
  char *restransformation_name;

  char *nameformat;
  int firstindex;
  int lastindex;

  enumTypeTransfo transformation_type;

  enumTransformationAveraging operation;

  bal_estimator estimator;

  typeComputation computation;

  int print_time;

} local_parameter;





/*------- Definition des fonctions statiques ----------*/
static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p );
static void _InitParam( local_parameter *par );
static double _GetTime();
static double _GetClock();




static char *program = NULL;

static char *usage = "[[-transformation-list|-list] %s]\n\
 [-format %s -f[irst] %d -l[ast] %d]\n\
 [[-res] image-out]\n\
 [-transformation-type|-transformation|-trsf-type %s]\n\
 [-mean|-robust-mean]\n\
 [-estimator-type|-estimator|-es-type wlts|lts|wls|ls]\n\
 [-lts-cut|-lts-fraction %lf] [-lts-deviation %f] [-lts-iterations %d]\n\
 [-fluid-sigma|-lts-sigma %lf %lf %lf]\n\
 [-lts-cut|-lts-fraction %lf]\n\
 [-streaming | -memory]\n\
 [-v] [-help]";

static char *detail = "\
[-transformation-list|-list %s] # text file = list of transformations to be processed\n\
[-format %s]   # format 'a la printf' of transformations to be processed, must contain a '%d'\n\
[-first %d]    # first value of the index in the format\n\
[-last %d]     # last value of the index in the format\n\
[-res %s]      # output transformation\n\
[-transformation-type|-transformation|-trsf-type %s] # transformation type\n\
  translation2D, translation3D, translation-scaling2D, translation-scaling3D,\n\
  rigid2D, rigid3D, rigid, similitude2D, similitude3D, similitude,\n\
  affine2D, affine3D, affine, vectorfield2D, vectorfield3D, vectorfield, vector\n\
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

  stringList transformationFileList;
  bal_transformation resTrsf;
  bal_transformationList transformationList; 

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
  


  /* transformation file names
   */

  initStringList( &transformationFileList );

  if ( p.nameformat != (char*)NULL ) {
    if ( buildStringListFromFormat( p.nameformat, p.firstindex, p.lastindex, &transformationFileList ) != 1 ) {
      _ErrorParse( "unable to build input transformation list from format\n", 0);
    }
  }
  else if ( buildStringListFromFile( p.transformationlist_name, &transformationFileList ) != 1 ) {
    _ErrorParse( "unable to read input transformation list from file\n", 0);
  }

  if ( _debug_ ) printStringList( stderr, &transformationFileList, "Input transformations" );

  
  
  /* reading transformations
   */
  BAL_InitTransformationList( &transformationList );
  
  if ( BAL_AllocTransformationList( &transformationList, transformationFileList.n_data ) != 1 ) {
    freeStringList( &transformationFileList );
    _ErrorParse( "error when allocating input transformations\n", 0);
  }
  
  if (  BAL_ReadTransformationList( &transformationList, &transformationFileList ) != 1 ) {
    BAL_FreeTransformationList( &transformationList );
    freeStringList( &transformationFileList );
    _ErrorParse( "error when reading input transformations\n", 0);
  }

  if ( transformationList.n_trsfs < transformationList.n_allocated_trsfs ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: %d transformations were read, %d (%d) were expected\n",
	       program, transformationList.n_trsfs, transformationList.n_allocated_trsfs,
         transformationFileList.n_data );
    }
  }


  freeStringList( &transformationFileList );

  

  /* initializing result transformation
     should we check whether read transformations are compatible with
     given type ?
   */
  BAL_InitTransformation( &resTrsf );

  switch ( p.transformation_type ) {
  default :
    BAL_FreeTransformationList( &transformationList );
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
    
    if ( BAL_AllocTransformation( &resTrsf, p.transformation_type, (bal_image *)NULL ) != 1 ) {
      BAL_FreeTransformationList( &transformationList );
      _ErrorParse( "unable to allocate result transformation (matrix)\n", 0);
    }
    break;
    
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
    
    /* should provide a template
     */
     BAL_FreeTransformationList( &transformationList );
      _ErrorParse( "such transformation type not handled yet (vector field)\n", 0);
  }



  /* ...
   */
  if ( BAL_ComputeAverageTransformation( &transformationList,  &resTrsf, &(p.estimator) ) != 1 ) {
    BAL_FreeTransformation( &resTrsf );
    BAL_FreeTransformationList( &transformationList );
    _ErrorParse( "unable to compute average transformation\n", 0);
  }

  BAL_FreeTransformationList( &transformationList );


  if ( BAL_WriteTransformation( &resTrsf, p.restransformation_name ) != 1 ) {
    BAL_FreeTransformation( &resTrsf );
    _ErrorParse( "unable to write result transformation\n", 0);
  }

  BAL_FreeTransformation( &resTrsf );

  /* end
   */
  
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

      if ( strcmp ( argv[i], "-transformation-list" ) == 0  
		|| (strcmp ( argv[i], "-list" ) == 0 && argv[i][5] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -transformation-list...\n", 0 );
	if ( p->transformationlist_name != (char*)NULL || inputisread > 0 ) 
	  _ErrorParse( "parsing -transformation-list: input has already been parsed ...\n", 0 );
	p->transformationlist_name = argv[i];
	inputisread = 1;
      }

      else if ( strcmp ( argv[i], "-format" ) == 0 ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -format...\n", 0 );
	if ( p->transformationlist_name != (char*)NULL || inputisread > 0 ) 
	  _ErrorParse( "parsing -format: input has already been parsed ...\n", 0 );
	p->nameformat = argv[i];
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

      else if ( strcmp ( argv[i], "-res" ) == 0  && argv[i][4] == '\0' ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -res...\n", 0 );
	if ( p->restransformation_name != (char*)NULL || outputisread > 0 ) 
	  _ErrorParse( "parsing -res: output has already been parsed ...\n", 0 );
	p->restransformation_name = argv[i];
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


      /* operation
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
	BAL_IncrementVerboseInBalTransformationAveraging(  );
	BAL_IncrementVerboseInBalTransformation(  );
      }
      
      else if ( strcmp ( argv[i], "--help" ) == 0 
		|| ( strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0' )
		|| ( strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0' )
		|| ( strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0' ) ) {
	_ErrorParse( NULL, 1 );
      }
      
       else if ( (strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0') ) {
	_debug_ = 1;
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
	p->transformationlist_name = argv[i];
	inputisread = 1;
      }
      else if ( outputisread == 0 ) {
	p->restransformation_name = argv[i];
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
  p->transformationlist_name = (char*)NULL;
  p->restransformation_name = (char*)NULL;
    
  p->nameformat = (char*)NULL;
  p->firstindex = 0;
  p->lastindex = 0;

  p->transformation_type = UNDEF_TRANSFORMATION;

  p->operation = _MEAN_;

  BAL_InitEstimator( &(p->estimator) );
  p->estimator.type = TYPE_LS;
  p->estimator.retained_fraction = 0.7500000;

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

