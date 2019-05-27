/*************************************************************************
 * TVmembrane.c -
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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include <pixel-operation.h>
#include <vtmalloc.h>

#include <api-TVmembrane.h>

/*
 * #include <vt_common.h>
 */




static int _verbose_ = 1;
static int _debug_ = 0;



/* static function definitions
 */

static vt_image *_readImage( char *name );
static vt_image *_readHessianImage( char *binaryname, char *suffix );

static char *_Array2Str( int argc, char *argv[] );
static char *_BaseName( char *p );
static double _GetTime();
static double _GetClock();






int main( int argc, char *argv[] )
{
  lineCmdParamTVmembrane par;
  char *lineoptions;


  int flag_3D = 1;
  vt_image *imageBin = (vt_image*)NULL;
  vt_image *imageTht = (vt_image*)NULL;
  vt_image *imagePhi = (vt_image*)NULL;
  vt_image *imageMask = (vt_image*)NULL;
  vt_image* imagesIn[3] = { (vt_image*)NULL, (vt_image*)NULL, (vt_image*)NULL };
  vt_3Dtensor theTensor;
  mt_2Dtensor theTensor2D;

  char name[DOUBLESTRINGLENGTH];
  char *ptrSuffix;
  char *inrSuffix = "inr";

  double time_init = _GetTime();
  double time_exit;
  double clock_init = _GetClock();
  double clock_exit;


  /* parameter initialization
   */
  API_InitParam_TVmembrane( &par );



  /* parameter parsing
   */
  if ( argc <= 1 )
      API_ErrorParse_TVmembrane( _BaseName( argv[0] ), (char*)NULL, 0 );
  API_ParseParam_TVmembrane( 1, argc, argv, &par );
  
  ptrSuffix = ( par.outputSuffix != (char*)NULL && par.outputSuffix[0] != '\0' ) ? par.outputSuffix : inrSuffix;


  /* input images reading
   */

  if ( par.initHessian == NONE ) {

    imageBin = _VT_Inrimage( par.inputBinary );
    if ( imageBin == (vt_image*)NULL ) {
      fprintf(stderr, "%s unreadable\n", par.inputBinary );
      API_ErrorParse_TVmembrane( _BaseName( argv[0] ), "unable to read input binary image  ...\n", 0 );
    }
    if ( _verbose_ )
      fprintf( stdout, "Input image is '%s'\n", par.inputBinary );

    if ( par.dimension == 2 || imageBin->dim.z == 1 )
      flag_3D = 0;
  }

  else {

    imageBin = _readImage( par.inputBinary );
    if ( imageBin == (vt_image*)NULL ) {
        fprintf( stderr, "binary image '%s' unreadable\n", par.inputBinary );
        API_ErrorParse_TVmembrane( _BaseName( argv[0] ), "unable to read input binary image ...\n", 0 );
    }
    if ( _verbose_ )
      fprintf( stdout, "Binary input image is '%s'\n", par.inputBinary );

    if ( par.dimension == 2 || imageBin->dim.z == 1 )
      flag_3D = 0;

    imageTht = _readHessianImage( par.inputBinary, ".theta" );
    if ( imageTht == (vt_image*)NULL ) {
        VT_FreeImage( imageBin );
        VT_Free( (void**)&imageBin );
        fprintf( stderr, "theta image '%s' unreadable\n", name);
        API_ErrorParse_TVmembrane( _BaseName( argv[0] ), "unable to read input theta image ...\n", 0 );
    }

    if ( flag_3D == 1 ) {
      imagePhi = _readHessianImage( par.inputBinary, ".phi" );
      if ( imageTht == (vt_image*)NULL ) {
          VT_FreeImage( imageTht );
          VT_Free( (void**)&imageTht );
          VT_FreeImage( imageBin );
          VT_Free( (void**)&imageBin );
          fprintf( stderr, "theta image '%s' unreadable\n", name);
          API_ErrorParse_TVmembrane( _BaseName( argv[0] ), "unable to read input phi image ...\n", 0 );
      }
    }
  }


  /* is there a mask image ?
   */
  if ( par.inputMask != (char*)NULL && par.inputMask[0] != '\0' ) {

      int theDim[3] = { (int)imageBin->dim.x,
                      (int)imageBin->dim.y,
                      (int)imageBin->dim.z };

      imageMask = _readImage( par.inputMask );
      if ( imageMask == (vt_image*)NULL ) {
        if ( imagePhi != (vt_image*)NULL ) {
            VT_FreeImage( imagePhi );
            VT_Free( (void**)&imagePhi );
        }
        if ( imageTht != (vt_image*)NULL ) {
            VT_FreeImage( imageTht );
            VT_Free( (void**)&imageTht );
        }
        VT_FreeImage( imageBin );
        VT_Free( (void**)&imageBin );
        fprintf(stderr, "mask image '%s' unreadable\n", par.inputMask );
        API_ErrorParse_TVmembrane( _BaseName( argv[0] ), "unable to read input mask image ...\n", 0 );
      }

      fprintf(stdout, "Mask image is '%s'\n", par.inputMask);

      /* mask application
       */
      if ( maskImage( imageBin->buf, imageBin->type,
                      imageMask->buf, imageMask->type,
                      imageBin->buf, imageBin->type, theDim ) != 1 ) {
          VT_FreeImage( imageMask );
          VT_Free( (void**)&imageMask );
          if ( imagePhi != (vt_image*)NULL ) {
              VT_FreeImage( imagePhi );
              VT_Free( (void**)&imagePhi );
          }
          if ( imageTht != (vt_image*)NULL ) {
              VT_FreeImage( imageTht );
              VT_Free( (void**)&imageTht );
          }
          VT_FreeImage( imageBin );
          VT_Free( (void**)&imageBin );
          API_ErrorParse_TVmembrane( _BaseName( argv[0] ), "error when masking image ...\n", 0 );
      }

      VT_FreeImage( imageMask );
      VT_Free( (void**)&imageMask );
  }





  /* API calls
   */

  lineoptions = _Array2Str( argc, argv );
  if ( lineoptions == (char*)NULL ) {
      if ( imagePhi != (vt_image*)NULL ) {
          VT_FreeImage( imagePhi );
          VT_Free( (void**)&imagePhi );
      }
      if ( imageTht != (vt_image*)NULL ) {
          VT_FreeImage( imageTht );
          VT_Free( (void**)&imageTht );
      }
      VT_FreeImage( imageBin );
      VT_Free( (void**)&imageBin );
      API_ErrorParse_TVmembrane( _BaseName( argv[0] ), "unable to translate command line options ...\n", 0 );
  }



  /* 1. Echantillonnage de l'image
   */

  if ( API_Sampling( imageBin, lineoptions, (char*)NULL ) != 1 ) {
      vtfree( lineoptions );
      if ( imagePhi != (vt_image*)NULL ) {
          VT_FreeImage( imagePhi );
          VT_Free( (void**)&imagePhi );
      }
      if ( imageTht != (vt_image*)NULL ) {
          VT_FreeImage( imageTht );
          VT_Free( (void**)&imageTht );
      }
      VT_FreeImage( imageBin );
      VT_Free( (void**)&imageBin );
      API_ErrorParse_TVmembrane( _BaseName( argv[0] ), "some error occurs during sampling ...\n", -1 );
  }

  name[0] = '\0';
  if ( par.outputSampling != (char*)NULL && par.outputSampling[0] != '\0' ) {
          sprintf( name, "%s", par.outputSampling );
  }
  else if ( par.writeImages == 1 ) {
    if ( par.outputGeneric != (char*)NULL && par.outputGeneric[0] != '\0' ) {
      sprintf( name, "%s.sample.%s", par.outputGeneric, ptrSuffix );
    }
    else {
        sprintf( name, "sample.%s", ptrSuffix );
    }
  }
  if ( name[0] != '\0' ) {
    if ( VT_WriteInrimageWithName( imageBin, name ) != 1 ) {
        vtfree( lineoptions );
        if ( imagePhi != (vt_image*)NULL ) {
            VT_FreeImage( imagePhi );
            VT_Free( (void**)&imagePhi );
        }
        if ( imageTht != (vt_image*)NULL ) {
            VT_FreeImage( imageTht );
            VT_Free( (void**)&imageTht );
        }
        VT_FreeImage( imageBin );
        VT_Free( (void**)&imageBin );
        API_ErrorParse_TVmembrane( _BaseName( argv[0] ), "error when writing sampled image ...\n", -1 );
    }
  }



  /* 2. Tensor voting
   */

  /* Ecriture des adresses dans l'input imagesIn de la fonction de TV
   */
  imagesIn[0] = imageBin;
  imagesIn[1] = imageTht;
  if (flag_3D == 1) imagesIn[2] = imagePhi;

  /*
   * ALLOCATION IMAGE TENSEUR
   */

  if ( flag_3D == 1 ) {
     if ( VT_Alloc3DtensorFromImage( &theTensor, imageBin, par.outputGeneric ) != 1 ) {
         vtfree( lineoptions );
         if ( imagePhi != (vt_image*)NULL ) {
             VT_FreeImage( imagePhi );
             VT_Free( (void**)&imagePhi );
         }
         if ( imageTht != (vt_image*)NULL ) {
             VT_FreeImage( imageTht );
             VT_Free( (void**)&imageTht );
         }
         VT_FreeImage( imageBin );
         VT_Free( (void**)&imageBin );
         API_ErrorParse_TVmembrane( _BaseName( argv[0] ), "unable to allocate tensor image ...\n", -1 );
     }
  }
  else {
      if ( VT_Alloc2DTensorFromImage( &theTensor2D, imageBin, par.outputGeneric ) != 1 ) {
          vtfree( lineoptions );
          if ( imageTht != (vt_image*)NULL ) {
              VT_FreeImage( imageTht );
              VT_Free( (void**)&imageTht );
          }
          VT_FreeImage( imageBin );
          VT_Free( (void**)&imageBin );
          API_ErrorParse_TVmembrane( _BaseName( argv[0] ), "unable to allocate tensor image ...\n", -1 );
      }
  }

  /* API call
   * tensor voting
   */

  /* Calcul de l'image de tenseurs theTensor2D (cas 2D) ou theTensor (cas 3D)  */
  if ( flag_3D == 1 ) {
      if ( API_TVmembrane( imagesIn, &theTensor, lineoptions, (char*)NULL ) != 1 ) {
          VT_Free3Dtensor(&theTensor);
          vtfree( lineoptions );
          if ( imagePhi != (vt_image*)NULL ) {
              VT_FreeImage( imagePhi );
              VT_Free( (void**)&imagePhi );
          }
          if ( imageTht != (vt_image*)NULL ) {
              VT_FreeImage( imageTht );
              VT_Free( (void**)&imageTht );
          }
          VT_FreeImage( imageBin );
          VT_Free( (void**)&imageBin );
          API_ErrorParse_TVmembrane( _BaseName( argv[0] ), "some error occurs during processing ...\n", -1 );
    }
  }
  else {
      if ( API_TVmembrane2D( imagesIn, &theTensor2D, lineoptions, (char*)NULL ) != 1) {
          VT_Free2DTensor(&theTensor2D);
          vtfree( lineoptions );
          if ( imageTht != (vt_image*)NULL ) {
              VT_FreeImage( imageTht );
              VT_Free( (void**)&imageTht );
          }
          VT_FreeImage( imageBin );
          VT_Free( (void**)&imageBin );
          API_ErrorParse_TVmembrane( _BaseName( argv[0] ), "some error occurs during processing ...\n", -1 );
    }
  }



  /* we can release some images
   */
  if ( imagePhi != (vt_image*)NULL ) {
      VT_FreeImage( imagePhi );
      VT_Free( (void**)&imagePhi );
  }
  if ( imageTht != (vt_image*)NULL ) {
      VT_FreeImage( imageTht );
      VT_Free( (void**)&imageTht );
  }
  VT_FreeImage( imageBin );
  VT_Free( (void**)&imageBin );
  imagesIn[0] = (vt_image*)NULL;
  imagesIn[1] = (vt_image*)NULL;
  imagesIn[2] = (vt_image*)NULL;



  /* output image writing
   */

  if ( _VT_VERBOSE_ )
    fprintf(stdout, "Ecriture des images...\n");
  if ( flag_3D == 1 ) {
    if ( par.outputGeneric != (char*)NULL && par.outputGeneric[0] != '\0' ) {
      if ( VT_Write3DtensorWithNames( &theTensor, par.outputGeneric, par.outputGeneric,
                                      par.outputGeneric, par.outputGeneric, par.outputSuffix ) != 1 ) {
          VT_Free3Dtensor(&theTensor);
          vtfree( lineoptions );
          API_ErrorParse_TVmembrane( _BaseName( argv[0] ), "some error occurs during writing images ...\n", -1 );
      }
    }
    else {
        if ( VT_Write3DtensorWithNames( &theTensor, par.outputTensor, par.outputEigenvalues,
                                        par.outputAngles, par.outputBinary, par.outputSuffix ) != 1 ) {
            VT_Free3Dtensor(&theTensor);
            vtfree( lineoptions );
            API_ErrorParse_TVmembrane( _BaseName( argv[0] ), "some error occurs during writing images ...\n", -1 );
        }
    }
  }
  else {
      if ( par.outputGeneric != (char*)NULL && par.outputGeneric[0] != '\0' ) {
        if ( VT_Write2DtensorWithNames( &theTensor2D, par.outputGeneric, par.outputGeneric,
                                        par.outputGeneric, par.outputGeneric, par.outputSuffix ) != 1 ) {
            VT_Free2DTensor(&theTensor2D);
            vtfree( lineoptions );
            API_ErrorParse_TVmembrane( _BaseName( argv[0] ), "some error occurs during writing images ...\n", -1 );
        }
      }
      else {
          if ( VT_Write2DtensorWithNames( &theTensor2D, par.outputTensor, par.outputEigenvalues,
                                          par.outputAngles, par.outputBinary, par.outputSuffix ) != 1 ) {
              VT_Free2DTensor(&theTensor2D);
              vtfree( lineoptions );
              API_ErrorParse_TVmembrane( _BaseName( argv[0] ), "some error occurs during writing images ...\n", -1 );
          }
      }
  }




  /* 3. shape maxima extraction (API qui fait allocation, extraction, ecriture et desallocation)
   *
   */

  if ( par.shape != NO_EXTRACTION ) {
       if ( flag_3D == 1 ) {
           if ( API_ShapeExtraction( &theTensor, lineoptions, (char*)NULL ) != 1 ) {
               VT_Free3Dtensor(&theTensor);
               vtfree( lineoptions );
               API_ErrorParse_TVmembrane( _BaseName( argv[0] ), "some error occurs during shape extraction ...\n", -1 );
           }
       }
       else {
           if ( API_ShapeExtraction2D( &theTensor2D, lineoptions, (char*)NULL ) != 1 ) {
               VT_Free3Dtensor(&theTensor);
               vtfree( lineoptions );
               API_ErrorParse_TVmembrane( _BaseName( argv[0] ), "some error occurs during shape extraction ...\n", -1 );
           }
       }
  }

  /* memory freeing
   */
  vtfree( lineoptions );

  if (flag_3D == 1)
    VT_Free3Dtensor(&theTensor);
  else
    VT_Free2DTensor(&theTensor2D);

  if ( par.trace_allocations ) {
    fprintfVtMallocTrace( stderr );
    clearVtMalloc();
  }

  /* temps de calcul
   */

  time_exit = _GetTime();
  clock_exit = _GetClock();

  if ( par.print_time ) {

    fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName( argv[0] ), time_exit - time_init );
    fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
    fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );

    if ( par.outputGeneric != (char*)NULL && par.outputGeneric[0] != '\0' ) {
      sprintf( name, "%s.timesampling.txt", par.outputGeneric );
      FILE *fichier = fopen (name , "a" );
      if ( fichier == NULL ) {
        perror (name);
      }
      else {
          int i;
          fprintf(fichier, "\n");
          for(i=0;i<argc;i++) fprintf( fichier, "%s ", argv[i]);
          fprintf(fichier, "\n");
          fprintf(fichier, "%s: elapsed (real) time = %f\n", _BaseName( argv[0] ), time_exit - time_init );
          fprintf(fichier, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
          fprintf(fichier, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );
          fclose(fichier);
      }
    }
  }

  return( 0 );
}





/************************************************************
 *
 * static functions
 *
 ************************************************************/


static char *suffixes[] = { ".inr", ".inr.gz",
                            ".hdr",
                            ".mha", ".mha.gz",
                            ".tif",
                            "" };

static vt_image *_readImage( char *name )
{
  char *proc = "_readImage";
  char tmp[STRINGLENGTH];
  vt_image *im;
  int i;

  if ( _debug_ )
    fprintf( stderr, "%s: try to read '%s' ... ", proc, name );

  im = _VT_Inrimage( name );

  if ( _debug_ ) {
    if ( im == (vt_image *)NULL ) fprintf( stderr, "failure\n" );
    else fprintf( stderr, "success\n" );
  }

  for ( i=0; strlen(suffixes[i]) > 0 && im == (vt_image *)NULL; i++ ) {

    sprintf( tmp, "%s%s", name, suffixes[i] );

    if ( _debug_ )
      fprintf( stderr, "%s: try to read '%s' ... ", proc, tmp );

    im = _VT_Inrimage( tmp );

    if ( _debug_ ) {
      if ( im == (vt_image *)NULL ) fprintf( stderr, "failure\n" );
      else fprintf( stderr, "success\n" );
    }
  }
  return( im );
}




static vt_image *_readHessianImage( char *binaryname, char *suffix )
{
  char *proc = "_readHessianImage";
  char tmp[STRINGLENGTH];
  int i, j, lengthbinary, l;


  /* remove image suffix if any
   */
  lengthbinary = strlen( binaryname );

  for ( i=0, tmp[0] = '\0'; strlen(suffixes[i]) > 0 && tmp[0] =='\0'; i++ ) {
    l = strlen(suffixes[i]);
    if ( strncmp( &(binaryname[lengthbinary-l]), suffixes[i], l ) == 0 ) {
      if ( _debug_ )
        fprintf( stderr, "%s: recognize suffix '%s' ... \n", proc, suffixes[i] );
      (void)strncpy( tmp, binaryname, lengthbinary-l );
      tmp[lengthbinary-l]='\0';
    }
  }

  if ( tmp[0] =='\0' ) {
    (void)strcpy( tmp, binaryname );
  }

  /* remove stuff from the last '.', if any
   */
  l = strlen( tmp );

  for ( j=-1, i=l-1; i>=0 && j == -1; i-- ) {
    if ( tmp[i] == '.' ) j = i;
  }
  if (j>=0)
    for ( i=j; i<l; i++ ) tmp[i] = '\0';

  /* prefix is recognized
   */
  if ( _debug_ )
      fprintf(stdout, "%s: read file with prefix '%s'...\n", proc, tmp );

  /* read image
   */
  if ( suffix != (char*)NULL && suffix[0] != '\0' )
    strcat( tmp, suffix );

  return( _readImage( tmp ) );
}





static char *_Array2Str( int argc, char *argv[] )
{
  char *proc = "_Array2Str";
  int i, l;
  char *s, *t;

  if ( argc <= 1 || argv == (char**)NULL ) {
    if ( _verbose_ >= 2 )
      fprintf( stderr, "%s: no options in argv[]\n", proc );
    return( (char*)NULL );
  }

  /* there are argc-1 strings
   * compute the sum of string lengths from 1 to argc-1
   * + number of interval between successive strings (argc-2)
   * + 1 to add a trailing '\0'
   */
  for ( l=argc-1, i=1; i<argc; i++ ) {
    l += strlen( argv[i] );
  }

  s = (char*)vtmalloc( l * sizeof( char ), "s", proc );
  if ( s == (char*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( (char*)NULL );
  }

  for ( t=s, i=1; i<argc; i++ ) {
    (void)strncpy( t, argv[i], strlen( argv[i] ) );
    t += strlen( argv[i] );
    if ( i < argc-1 ) {
      *t = ' ';
      t++;
    }
    else {
      *t = '\0';
    }
  }

  return( s );
}



static char *_BaseName( char *p )
{
  int l;
  if ( p == (char*)NULL ) return( (char*)NULL );
  l = strlen( p ) - 1;
  while ( l >= 0 && p[l] != '/' ) l--;
  if ( l < 0 ) l = 0;
  if ( p[l] == '/' ) l++;
  return( &(p[l]) );
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
