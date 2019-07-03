/*************************************************************************
 * changeMultipleTrsfs.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mar 21 jan 2014 22:43:26 CET
 *
 * ADDITIONS, CHANGES
 *
 */


#include <sys/time.h> /* gettimeofday() */
#include <time.h> /* clock() */
#include <string.h>

#include <string-tools.h>

#include <bal-transformation-tools.h>
#include <bal-transformation-list-tools.h>
#include <bal-transformation-averaging.h>

static int _verbose_ = 1;
static int _debug_ = 0;








typedef struct local_parameter {

  char *thetransformation_format;
  int firstindex;
  int lastindex;
  int refindex;
  char *thetransformation_list;


  char *theimage_name;
  char *theimage_format;
  char *theimage_list;

  int threshold;

  bal_doublePoint template_voxel;

  char *restransformation_format;
  char *restransformation_list;

  /* build the resulting template from
   * a given image
   */
  int reftemplate;
  char *res_template_name;

  int isotropic;
  bal_doublePoint res_template_voxel;
  int margin[3];
  int extension[3];

  enumTypeTransfo transformation_type;

  int print_time;

} local_parameter;





/*------- Definition des fonctions statiques ----------*/
static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p );
static void _InitParam( local_parameter *par );
static double _GetTime();
static double _GetClock();




static char *program = NULL;

static char *usage = "[-trsf-format|-format] format-in\n\
 -f[irst] %d -l[ast] %d\n\
 [-index-reference|-r[ef|eference] %d]\n\
 [-trsf-list %s]\n\
 [-template|-image %s]\n\
 [-template-format|-image-format %s]\n\
 [-template-list|-image-list %s]\n\
 [-threshold|-t %d]\n\
 [[-res-trsf-format|-res-format|-res] format-out] \n\
 [-res-trsf-list|-res-list %s] \n\
 [-res-index-template|-index-template %d]\n\
 [-result-template|-res-template|-res-t %s] \n\
 [-result-isotropic-voxel|-result-isotropic|-res-iso [%lf]]\n\
 [-result-template-voxel|-result-voxel-size|-res-voxel|-res-pixel|-rvs %lf %lf [%lf]]\n\
 [-margin %d]\n\
 [-x|-y|-z]\n\
 [-transformation-type|-transformation|-trsf-type %s]\n\
 [-v] [-help]";

static char *detail = "\
\n\
given a list of transformation towards a reference,\n\
compute a new template that contains all transformed images as well\n\
as the new transformations\n\
\n\
[-trsf-format|-format] format-in  # format 'a la printf' of transformation files\n\
             # to be processed. It must contain one '%d'\n\
             # depicts transformations of the form T_{i<-ref}\n\
             # (ie allows to resample image I_i in geometry of I_ref)\n\
[-first %d]  # first value of the index in the format\n\
[-last %d]   # last value of the index in the format\n\
[-index-reference|-reference %d] # index of the reference transformation\n\
             # the corresponding image will only be translated\n\
             # if none, only translations will change in transformations\n\
             # (ie reference is not changed)\n\
[-trsf-list %s] # list of transformations to be processed\n\
[-template|-t|-image %s] # template image corresponding to the input images\n\
             # one template for all\n\
[-template-format|-image-format %s] # input templates/images\n\
             # one template per transformation\n\
[-template-list|-image-list %s] # list of templates/images\n\
             # one template per transformation\n\
[-threshold|-t %d] # threshold on input templates/images to compute the\n\
             # useful bounding box (else it is the entire image)\n\
             # points whose values are greater or equal to the\n\
             # threshold are kept in the template\n\
[-template-voxel|-voxel-size|-voxel|-pixel|-vs %lf %lf [%lf]] #\n\
             # change the voxel sizes of the input template\n\
             # (when a single input template is given)\n\
[-res-trsf-format|-res-format|-res] format-out\n\
             # format 'a la printf' for output transformations\n\
             # will allow to resample input image into the resulting\n\
             # template (thus still of the form T_{i<-ref})\n\
             # reference is changed if '-index-reference' is used\n\
[-res-trsf-list|-res-list %s] # list of result transformations\n\
[-res-index-template|-index-template %d] # if given, the output template\n\
             # is a transformed image\n\
             # used with '-template' : the transformed single template image \n\
             # used with '-template-format' : the transformed template image \n\
             #   of given index\n\
[[-result-template|-res-template|-res-t  %s] # \n\
             # output template image corresponding to the output transformations\n\
[-result-isotropic-voxel|-res-iso [%lf]] # make voxels isotropic for the output template\n\
             # (if no value is given, uses the smallest voxel size from the template(s))\n\
[-result-template-voxel|-result-voxel-size|-res-voxel|-res-pixel|-rvs %lf %lf [%lf]] #\n\
             # result template image voxel sizes\n\
[-margin %d]   # add a margin (in voxels)\n\
[-x|-y|-z]     # dimension to be 'extended', if none is specified, all are extended\n\
[-transformation-type|-transformation|-trsf-type %s] # transformation type\n\
  translation2D, translation3D, translation-scaling2D, translation-scaling3D,\n\
  rigid2D, rigid3D, rigid, similitude2D, similitude3D, similitude,\n\
  affine2D, affine3D, affine, vectorfield2D, vectorfield3D, vectorfield, vector\n\
 -v : mode verbose\n\
\n";







int main( int argc, char *argv[] )
{
  local_parameter p;

  stringList theTransformationFileList;
  bal_transformationList theTransformationList; 

  stringList theImageFileList;
  bufferType typeImage = TYPE_UNKNOWN;

  typeBoxList boxList;

  int i, ntransformations = 0;

  stringList resTransformationFileList;
  bal_transformationList resTransformationList; 

  bal_image theTemplate;
  bal_image resTemplate;
  bufferType typeTemplate = UCHAR;
  char nameTemplate[STRINGLENGTH];
  float vsize;

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
  

  

  if ( p.firstindex == p.lastindex ) {
    p.reftemplate = p.firstindex;
  }
  nameTemplate[0] = '\0';




  /***************************************************
   *
   * inputs
   *
   ***************************************************/


  /* transformation file names
     string[j][i] is printf( format, i, j )
     since transformations are denoted floating_on_reference and allows to resample
     the floating image onto the reference
     string[j][i] refers to the transformation T_{i<-j}
   */
  initStringList( &theTransformationFileList );

  if ( p.thetransformation_format != (char*)NULL ) {
    if ( buildStringListFromFormat( p.thetransformation_format,
                                    p.firstindex, p.lastindex,
                                    &theTransformationFileList ) != 1 ) {
      _ErrorParse( "unable to build input transformation list from format\n", 0);
    }
  }

  if ( p.thetransformation_list != (char*)NULL ) {
    if ( buildStringListFromFile( p.thetransformation_list,
                                  &theTransformationFileList ) != 1 ) {
      _ErrorParse( "unable to build input transformation list from file\n", 0);
    }

  }

  ntransformations = theTransformationFileList.n_data;


  if ( _debug_ >= 3 ) 
    printStringList( stderr, &theTransformationFileList, "Input transformations" );





  /* reading transformations
   */
  BAL_InitTransformationList( &theTransformationList );

  if ( BAL_AllocTransformationList( &theTransformationList, 
                                     ntransformations ) != 1 ) {
    freeStringList( &theTransformationFileList );
    _ErrorParse( "error when allocating input transformations\n", 0);
  }

  
  if ( BAL_ReadTransformationList( &theTransformationList, &theTransformationFileList ) != 1 ) {
    BAL_FreeTransformationList( &theTransformationList );
    freeStringList( &theTransformationFileList );
    _ErrorParse( "error when reading input transformations\n", 0);
  }

  freeStringList( &theTransformationFileList );

  

  /* box list
   */
  _initBoxList( &boxList );

  if ( _allocBoxList( &boxList, ntransformations ) != 1 ) {
    BAL_FreeTransformationList( &theTransformationList );
    _ErrorParse( "error when allocating box list\n", 0);
  }

  /* a same template for all transformations
   */
  if ( p.theimage_name != (char*)NULL ) {
    if ( _fillBoxesWithImageName( &boxList, &typeImage, p.theimage_name,
                                  &theTransformationList, (int)-100000,
                                  p.refindex - p.firstindex,
                                  &(p.template_voxel) ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to fill boxes with image '%s\n", argv[0], p.theimage_name );
      _freeBoxList( &boxList );
      BAL_FreeTransformationList( &theTransformationList );
      _ErrorParse( "unable to fill boxes with image\n", 0);
    }
    /* the resulting template will be
     * the transformed input image template
     */
    if ( p.firstindex <= p.reftemplate && p.reftemplate <= p.lastindex ) {
      typeTemplate = typeImage;
      strncpy( nameTemplate, p.theimage_name, STRINGLENGTH );
    }
  }

  /* a template (an image) per transformation
   */
  else if ( p.theimage_format != (char*)NULL || p.theimage_list != (char*)NULL ) {
    initStringList( &theImageFileList );
    if ( p.theimage_format != (char*)NULL ) {
      if ( buildStringListFromFormat( p.theimage_format,
                                      p.firstindex, p.lastindex,
                                      &theImageFileList ) != 1 ) {
        _freeBoxList( &boxList );
        BAL_FreeTransformationList( &theTransformationList );
        _ErrorParse( "unable to build input image list from format\n", 0);
      }
    }
    if ( p.theimage_list != (char*)NULL ) {
      if ( buildStringListFromFile( p.theimage_list,
                                    &theImageFileList ) != 1 ) {
        freeStringList( &theImageFileList );
        _freeBoxList( &boxList );
        BAL_FreeTransformationList( &theTransformationList );
        _ErrorParse( "unable to build input image list from file\n", 0);
      }
    }

    if ( theImageFileList.n_data != theTransformationList.n_trsfs ) {
      freeStringList( &theImageFileList );
      _freeBoxList( &boxList );
      BAL_FreeTransformationList( &theTransformationList );
      _ErrorParse( "image and transformation lists have different numbers of elements\n", 0);
    }

    if ( _fillBoxesWithImageNames( &boxList, &typeImage, &theImageFileList,
                                   &theTransformationList, p.threshold,
                                   p.refindex - p.firstindex,
                                   p.reftemplate - p.firstindex ) != 1 ) {
      if ( _verbose_ )
          fprintf( stderr, "%s: unable to fill boxes with image names\n", argv[0] );
       freeStringList( &theImageFileList );
       _freeBoxList( &boxList );
       BAL_FreeTransformationList( &theTransformationList );
       _ErrorParse( "unable to fill boxes with image format\n", 0);
    }

    freeStringList( &theImageFileList );
  }
  else {
    _freeBoxList( &boxList );
    BAL_FreeTransformationList( &theTransformationList );
    _ErrorParse( "no image template(s)\n", 0);
  }



  if ( _debug_ ) {
      _fprintfBoxList( stderr, &boxList, "after transformation" );
  }
 




  /***************************************************
   *
   * outputs
   *
   ***************************************************/
  

  /* initializing result transformations
     should we check whether read transformations are compatible with
     given type ?
  */
  BAL_InitTransformationList( &resTransformationList );

  switch ( p.transformation_type ) {
  default :
    _freeBoxList( &boxList );
    BAL_FreeTransformationList( &theTransformationList );
    _ErrorParse( "such transformation type not handled yet (or no type was given)\n", 0);
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
      _freeBoxList( &boxList );
      BAL_FreeTransformationList( &theTransformationList );
      _ErrorParse( "unable to allocate result transformations (matrices)\n", 0);
    }
    break;
    
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
    
    /* should provide a template
     */
    _freeBoxList( &boxList );
    BAL_FreeTransformationList( &theTransformationList );
    _ErrorParse( "such transformation type not handled yet (vector field)\n", 0);
  }

  /* result template
   */

  BAL_InitImage( &resTemplate, NULL, 0, 0, 0, 0, typeTemplate );

  if ( p.res_template_voxel.x > 0 && p.res_template_voxel.y > 0 && p.res_template_voxel.z > 0 ) {
      resTemplate.vx = p.res_template_voxel.x;
      resTemplate.vy = p.res_template_voxel.y;
      resTemplate.vz = p.res_template_voxel.z;
  }
  else if ( p.isotropic ) {
    vsize = boxList.data[0].voxel.x;
    for ( i=0; i<boxList.n_data; i++ ) {
      if ( vsize > boxList.data[i].voxel.x ) vsize = boxList.data[i].voxel.x;
      if ( vsize > boxList.data[i].voxel.y ) vsize = boxList.data[i].voxel.y;
      if ( vsize > boxList.data[i].voxel.z ) vsize = boxList.data[i].voxel.z;
    }
    resTemplate.vx = vsize;
    resTemplate.vy = vsize;
    resTemplate.vz = vsize;
  }
  else {
    resTemplate.vx = boxList.data[0].voxel.x;
    resTemplate.vy = boxList.data[0].voxel.y;
    resTemplate.vz = boxList.data[0].voxel.z;
  }

  BAL_SetImageVoxelSizes( &resTemplate, resTemplate.vx, resTemplate.vy, resTemplate.vz );



  /***************************************************
   *
   * 
   *
   ***************************************************/

  if ( p.extension[0] == 0 && p.extension[1] == 0 && p.extension[2] == 0 )
      p.extension[0] = p.extension[1] = p.extension[2] = 1;

  if ( BAL_ChangeTransformationList( &theTransformationList, &boxList,
                                     &resTransformationList, &resTemplate,
                                     p.refindex - p.firstindex,
                                     p.margin, p.extension ) != 1 ) {
    BAL_FreeImage( &resTemplate );
    BAL_FreeTransformationList( &resTransformationList );
    _freeBoxList( &boxList );
    BAL_FreeTransformationList( &theTransformationList );
    _ErrorParse( "unable to recompute transformations\n", 0);
  }


  _freeBoxList( &boxList );
  BAL_FreeTransformationList( &theTransformationList );




  
   /* writing result transformations
   */

  initStringList( &resTransformationFileList );

  if ( p.restransformation_format != (char*)NULL || p.thetransformation_list != (char*)NULL ) {

    if ( p.restransformation_format != (char*)NULL ) {
      if ( buildStringListFromFormat( p.restransformation_format, p.firstindex, p.lastindex,
                                      &resTransformationFileList ) != 1 ) {
        BAL_FreeImage( &resTemplate );
        BAL_FreeTransformationList( &resTransformationList );
        _ErrorParse( "unable to build output transformation list from format\n", 0);
      }
    }

    if ( p.restransformation_list != (char*)NULL ) {
      if ( buildStringListFromFile( p.restransformation_list,
                                      &resTransformationFileList ) != 1 ) {
        freeStringList( &theImageFileList );
        BAL_FreeImage( &resTemplate );
        BAL_FreeTransformationList( &resTransformationList );
        _ErrorParse( "unable to build output transformation list from file\n", 0);
      }
    }


    if ( _debug_ >= 3 )
      printStringList( stderr, &resTransformationFileList, "Result transformations" );

    if ( resTransformationFileList.n_data != resTransformationList.n_trsfs ) {
      freeStringList( &resTransformationFileList );
      BAL_FreeImage( &resTemplate );
      BAL_FreeTransformationList( &resTransformationList );
      _ErrorParse( "image and transformation lists have different numbers of elements\n", 0);
    }


    if ( BAL_WriteTransformationList( &resTransformationList, &resTransformationFileList ) != 1 ) {
      freeStringList( &resTransformationFileList );
      BAL_FreeImage( &resTemplate );
      BAL_FreeTransformationList( &resTransformationList );
      _ErrorParse( "unable to write result transformations\n", 0);
    }

    freeStringList( &resTransformationFileList );

  }



  /* writing result template image
   */
  if ( p.res_template_name != (char*)NULL ) {

    if ( BAL_AllocImage( &resTemplate ) != 1 ) {
      BAL_FreeTransformationList( &resTransformationList );
      _ErrorParse( "unable to allocate result template image\n", 0);
    }

    /* the resulting template image is a transformed image
     */
    if ( nameTemplate[0] != '\0' ) {

      if ( BAL_ReadImage( &theTemplate, nameTemplate, 0 ) != 1 ) {
        BAL_FreeImage( &resTemplate );
        BAL_FreeTransformationList( &resTransformationList );
        if ( _verbose_ )
          fprintf( stderr, " ... unable to read image '%s'\n", nameTemplate );
        _ErrorParse( "unable to read input template image\n", 0);
      }

      if ( BAL_ResampleImage( &theTemplate, &resTemplate,
                              &(resTransformationList.data[p.reftemplate - p.firstindex]),
                              LINEAR ) != 1 ) {
        BAL_FreeImage( &resTemplate );
        BAL_FreeTransformationList( &resTransformationList );
        _ErrorParse( "unable to transform input template image\n", 0);
      }
    }

    if ( BAL_WriteImage( &resTemplate, p.res_template_name ) != 1 ) {
      _ErrorParse( "unable to write result template image\n", 0);
    }
    BAL_FreeImage( &resTemplate );
  }

  BAL_FreeTransformationList( &resTransformationList );





  

  
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

      if ( strcmp ( argv[i], "-trsf-format" ) == 0
           || (strcmp ( argv[i], "-format" ) == 0 && argv[i][7] == '\0') ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -trsf-format...\n", 0 );
        if ( p->thetransformation_format != (char*)NULL || inputisread > 0 )
          _ErrorParse( "parsing -trsf-format: input has already been parsed ...\n", 0 );
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
      else if ( (strcmp ( argv[i], "-r" ) == 0 && argv[i][2] == '\0') 
                || (strcmp ( argv[i], "-ref" ) == 0 && argv[i][4] == '\0')
                || (strcmp ( argv[i], "-reference" ) == 0)
                || (strcmp ( argv[i], "-index-reference" ) == 0) ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -index-reference ...\n", 0 );
        status = sscanf( argv[i], "%d", &(p->refindex) );
        if ( status <= 0 ) _ErrorParse( "parsing -index-reference ...", 0 );
      }

      else if ( strcmp ( argv[i], "-trsf-list" ) == 0 ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -trsf-list...\n", 0 );
        p->thetransformation_list = argv[i];
      }

      else if ( (strcmp ( argv[i], "-template") == 0 && argv[i][9] == '\0')
                || (strcmp ( argv[i], "-image") == 0 && argv[i][6] == '\0') ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -template", 0 );
        p->theimage_name = argv[i];
      }

      else if ( strcmp ( argv[i], "-template-format") == 0
                || strcmp ( argv[i], "-image-format") == 0 ) {
          i++;
          if ( i >= argc) _ErrorParse( "parsing -template-format", 0 );
          p->theimage_format = argv[i];
      }

      else if ( strcmp ( argv[i], "-template-list") == 0
                || strcmp ( argv[i], "-image-list") == 0 ) {
          i++;
          if ( i >= argc) _ErrorParse( "parsing -template-list", 0 );
          p->theimage_list = argv[i];
      }

      else if ( (strcmp ( argv[i], "-threshold") == 0 && argv[i][10] == '\0')
                || (strcmp ( argv[i], "-t") == 0 && argv[i][2] == '\0') ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -threshold", 0 );
        status = sscanf( argv[i], "%d", &(p->threshold) );
        if ( status <= 0 ) _ErrorParse( "parsing -threshold", 0 );
      }

      else if ( strcmp ( argv[i], "-template-voxel") == 0
                || strcmp ( argv[i], "-voxel-size") == 0
                || (strcmp (argv[i], "-voxel" ) == 0 && argv[i][6] == '\0')
                || (strcmp (argv[i], "-pixel" ) == 0 && argv[i][6] == '\0')
                || (strcmp (argv[i], "-vs" ) == 0 && argv[i][3] == '\0') ) {
        i ++;
        if ( i >= argc)    _ErrorParse( "parsing -template-voxel %lf", 0 );
        status = sscanf( argv[i], "%lf", &(p->template_voxel.x) );
        if ( status <= 0 ) _ErrorParse( "parsing -template-voxel %lf", 0 );
        i ++;
        if ( i >= argc)    _ErrorParse( "parsing -template-voxel %lf %lf", 0 );
        status = sscanf( argv[i], "%lf", &(p->template_voxel.y) );
        if ( status <= 0 ) _ErrorParse( "parsing -template-voxel %lf %lf", 0 );
        i ++;
        if ( i >= argc) p->template_voxel.z = 1;
        else {
          status = sscanf( argv[i], "%lf", &(p->template_voxel.z) );
          if ( status <= 0 ) {
            i--;
            p->template_voxel.z = 1;
          }
        }
      }

      /* outputs
       */

      else if ( strcmp ( argv[i], "-res-trsf-format" ) == 0
                || strcmp ( argv[i], "-res-format" ) == 0
                || (strcmp ( argv[i], "-res" ) == 0  && argv[i][4] == '\0') ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -res-trsf-format...\n", 0 );
        if ( p->restransformation_format != (char*)NULL || outputisread > 0 )
          _ErrorParse( "parsing -res-trsf-format: output has already been parsed ...\n", 0 );
        p->restransformation_format = argv[i];
        outputisread = 1;
      }

      else if ( strcmp ( argv[i], "-res-trsf-list" ) == 0
                || strcmp ( argv[i], "-res-list" ) == 0 ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -res-trsf-list...\n", 0 );
        p->restransformation_list = argv[i];
      }

      else if ( (strcmp ( argv[i], "-index-template" ) ==0)
                || (strcmp ( argv[i], "-res-index-template" ) == 0)
                || (strcmp ( argv[i], "-result-index-template" ) == 0) ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -result-index-template ...\n", 0 );
        status = sscanf( argv[i], "%d", &(p->reftemplate) );
        if ( status <= 0 ) _ErrorParse( "parsing -result-index-template ...", 0 );
      }

      else if ( strcmp ( argv[i], "-result-template") == 0
                || strcmp ( argv[i], "-res-template") == 0
                || (strcmp ( argv[i], "-res-t") == 0 && argv[i][6] == '\0') ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -res-template", 0 );
        p->res_template_name = argv[i];
      }

      /* computation tuning
       */


      else if ( strcmp ( argv[i], "-result-template-voxel") == 0
                || strcmp ( argv[i], "-result-voxel-size") == 0
                || (strcmp (argv[i], "-res-voxel" ) == 0 && argv[i][10] == '\0')
                || (strcmp (argv[i], "-res-pixel" ) == 0 && argv[i][10] == '\0')
                || (strcmp (argv[i], "-rvs" ) == 0 && argv[i][4] == '\0') ) {
        i ++;
        if ( i >= argc)    _ErrorParse( "parsing -result-template-voxel %lf", 0 );
        status = sscanf( argv[i], "%lf", &(p->res_template_voxel.x) );
        if ( status <= 0 ) _ErrorParse( "parsing -result-template-voxel %lf", 0 );
        i ++;
        if ( i >= argc)    _ErrorParse( "parsing -result-template-voxel %lf %lf", 0 );
        status = sscanf( argv[i], "%lf", &(p->res_template_voxel.y) );
        if ( status <= 0 ) _ErrorParse( "parsing -result-template-voxel %lf %lf", 0 );
        i ++;
        if ( i >= argc) p->res_template_voxel.z = 1;
        else {
          status = sscanf( argv[i], "%lf", &(p->res_template_voxel.z) );
          if ( status <= 0 ) {
            i--;
            p->res_template_voxel.z = 1;
          }
        }
      }

      else if ( (strcmp ( argv[i], "-res-iso" ) == 0 && argv[i][8] == '\0')
                || (strcmp ( argv[i], "-result-isotropic" ) == 0 && argv[i][17] == '\0')
                || strcmp ( argv[i], "-result-isotropic-voxel" ) == 0 ) {
        if ( i+1 < argc ) {
          status = sscanf( argv[i+1], "%lf", &(p->res_template_voxel.x) );
          if ( status <= 0 ) {
              p->res_template_voxel.x = -1.0;
          }
          else {
              i++;
              p->res_template_voxel.z = p->res_template_voxel.y = p->res_template_voxel.x;
          }
        }
        p->isotropic = 1;
      }

      else if ( strcmp ( argv[i], "-margin" ) == 0 ) {
        i++;
        if ( i >= argc) _ErrorParse( "parsing -margin ...\n", 0 );
        status = sscanf( argv[i], "%d", &(p->margin[0]) );
        if ( status <= 0 ) _ErrorParse( "parsing -margin ...", 0 );
        p->margin[1] = p->margin[2] = p->margin[0];
      }

      else if ( ( strcmp ( argv[i], "-x" ) == 0 && argv[i][2] == '\0' ) ) {
          p->extension[0] = 1;
      }
      else if ( ( strcmp ( argv[i], "-y" ) == 0 && argv[i][2] == '\0' ) ) {
          p->extension[1] = 1;
      }
      else if ( ( strcmp ( argv[i], "-z" ) == 0 && argv[i][2] == '\0' ) ) {
          p->extension[2] = 1;
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
        BAL_IncrementVerboseInBalTransformationListTools( );
        BAL_IncrementVerboseInBalTransformationAveraging( );
        BAL_IncrementVerboseInBalTransformation( );
      }
      else if ( (strcmp ( argv[i], "-no-verbose") == 0 && argv[i][11] == '\0')
                || (strcmp ( argv[i], "-noverbose") == 0 && argv[i][10] == '\0')
                || (strcmp ( argv[i], "-nv") == 0 && argv[i][3] == '\0') ) {
        _verbose_ = 0;
        BAL_SetVerboseInBalTransformationListTools( 0 );
        BAL_SetVerboseInBalTransformationAveraging( 0 );
        BAL_SetVerboseInBalTransformation( 0 );
      }
      
      else if ( strcmp ( argv[i], "--help" ) == 0 
                || ( strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0' ) ) {
        _ErrorParse( NULL, 1 );
      }

      else if ( ( strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0' )
                || ( strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0' ) ) {
        _ErrorParse( NULL, 0 );
      }

       else if ( (strcmp ( argv[i], "-debug" ) == 0 && argv[i][6] == '\0')
                 || (strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0') ) {
        if ( _debug_ <= 0 ) _debug_ = 1;
        else _debug_ ++;
        BAL_IncrementDebugInBalTransformationListTools( );
      }
      else if ( (strcmp ( argv[i], "-no-debug" ) == 0 && argv[i][9] == '\0')
                || (strcmp ( argv[i], "-nodebug" ) == 0 && argv[i][8] == '\0') ) {
       _debug_ = 0;
       BAL_SetDebugInBalTransformationListTools( 0 );
     }

      /* unknown option
       */
      else {
        fprintf(stderr,"unknown option: '%s'\n",argv[i]);
      }
    }
    
    /*--- saisie des noms d'images ---*/
    else if ( argv[i][0] != 0 ) {
      if ( inputisread == 0 && p->thetransformation_list == (char*)NULL ) {
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
  (void)fprintf(stderr,"Usage : %s %s\n", program, usage);
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
  p->refindex = -1;
  p->thetransformation_list = (char*)NULL;

  p->reftemplate = -1;

  p->theimage_name = (char*)NULL;
  p->theimage_format = (char*)NULL;
  p->theimage_list = (char*)NULL;

  p->threshold = -100000;

  p->template_voxel.x = -1.0;
  p->template_voxel.y = -1.0;
  p->template_voxel.z = -1.0;

  p->restransformation_format = (char*)NULL;
  p->restransformation_list = (char*)NULL;

  p->transformation_type = UNDEF_TRANSFORMATION;

  p->res_template_name = (char*)NULL;

  p->isotropic = 0;

  p->res_template_voxel.x = -1.0;
  p->res_template_voxel.y = -1.0;
  p->res_template_voxel.z = -1.0;

  p->margin[0] = 0;
  p->margin[1] = 0;
  p->margin[2] = 0;

  p->extension[0] = 0;
  p->extension[1] = 0;
  p->extension[2] = 0;

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

