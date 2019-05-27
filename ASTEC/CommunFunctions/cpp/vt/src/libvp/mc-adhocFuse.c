/*************************************************************************
 * mc-adhocFuse.c - template for executable creation
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
#include <time.h>
#include <sys/time.h>

#include <vtmalloc.h>

#include <api-adhocFuse.h>

#include <vt_common.h>





static int _verbose_ = 1;





/* static function definitions
 */

static char *_Array2Str( int argc, char *argv[] );
static char *_BaseName( char *p );
static double _GetTime();
static double _GetClock();






int main( int argc, char *argv[] )
{
  lineCmdParamAdhocFuse par;
  vt_image *inputIntImage = (vt_image*)NULL;
  vt_image *inputRecImage = (vt_image*)NULL;
  vt_image *inputMinMskImage = (vt_image*)NULL;
  vt_image *ptrInputMinMskImage = (vt_image*)NULL;
  vt_image *inputMaxMskImage = (vt_image*)NULL;
  vt_image *ptrInputMaxMskImage = (vt_image*)NULL;
  vt_image *inputSegImage = (vt_image*)NULL;

  ImageType outputIntImageType = UCHAR;

  vt_image outputMinImage;
  vt_image outputMaxImage;
  vt_image outputIntImage;
  vt_image outputFusImage;
  vt_image *ptrOutputMinImage = (vt_image*)NULL;
  vt_image *ptrOutputMaxImage = (vt_image*)NULL;
  vt_image *ptrOutputIntImage = (vt_image*)NULL;
  vt_image *ptrOutputFusImage = (vt_image*)NULL;

  char *lineoptions;

  double time_init = _GetTime();
  double time_exit;
  double clock_init = _GetClock();
  double clock_exit;


  /* parameter initialization
   */
  API_InitParam_adhocFuse( &par );
  VT_Image( &outputIntImage );
  VT_Image( &outputFusImage );



  /* parameter parsing
   */
  if ( argc <= 1 )
      API_ErrorParse_adhocFuse( _BaseName( argv[0] ), (char*)NULL, 0 );
  API_ParseParam_adhocFuse( 1, argc, argv, &par );
  
  if ( par.output_type != TYPE_UNKNOWN )
      outputIntImageType = par.output_type;




  /* input images reading
   */
  if ( par.input_intensity_name == (char*)NULL || par.input_intensity_name[0] == '\0' ) {
    API_ErrorParse_adhocFuse( _BaseName( argv[0] ), "no input intensity image ...\n", 0 );
  }

  inputIntImage = _VT_Inrimage( par.input_intensity_name );
  if ( inputIntImage == (vt_image*)NULL ) {
      API_ErrorParse_adhocFuse( _BaseName( argv[0] ), "unable to read input intensity image ...\n", 0 );
  }

  if ( par.input_reconstructed_name != (char*)NULL && par.input_reconstructed_name[0] != '\0' ) {
    inputRecImage = _VT_Inrimage( par.input_reconstructed_name );
    if ( inputRecImage == (vt_image*)NULL ) {
      VT_FreeImage( inputIntImage );
      VT_Free( (void**)&inputIntImage );
      API_ErrorParse_adhocFuse( _BaseName( argv[0] ), "unable to read input reconstruction image ...\n", 0 );
    }
    outputIntImageType = inputRecImage->type;
  }



  if ( par.input_minimum_mask_name != (char*)NULL && par.input_minimum_mask_name[0] != '\0' ) {
    inputMinMskImage = _VT_Inrimage( par.input_minimum_mask_name );
    if ( inputMinMskImage == (vt_image*)NULL ) {
      if ( inputRecImage != (vt_image*)NULL ) {
        VT_FreeImage( inputRecImage );
        VT_Free( (void**)&inputRecImage );
      }
      VT_FreeImage( inputIntImage );
      VT_Free( (void**)&inputIntImage );
      API_ErrorParse_adhocFuse( _BaseName( argv[0] ), "unable to read input minimum mask image ...\n", 0 );
    }
    ptrInputMinMskImage = inputMinMskImage;
  }

  if ( par.input_maximum_mask_name != (char*)NULL && par.input_maximum_mask_name[0] != '\0' ) {
    if ( strlen( par.input_minimum_mask_name ) == strlen( par.input_maximum_mask_name )
         && strncmp( par.input_maximum_mask_name, par.input_minimum_mask_name, strlen( par.input_minimum_mask_name ) ) == 0 ) {
       ptrInputMaxMskImage = inputMinMskImage;
    }
    else {
      inputMaxMskImage = _VT_Inrimage( par.input_maximum_mask_name );
      if ( inputMaxMskImage == (vt_image*)NULL ) {
        if ( inputMinMskImage != (vt_image*)NULL ) {
          VT_FreeImage( inputMinMskImage );
          VT_Free( (void**)&inputMinMskImage );
        }
        if ( inputRecImage != (vt_image*)NULL ) {
          VT_FreeImage( inputRecImage );
          VT_Free( (void**)&inputRecImage );
        }
        VT_FreeImage( inputIntImage );
        VT_Free( (void**)&inputIntImage );
        API_ErrorParse_adhocFuse( _BaseName( argv[0] ), "unable to read input maximum mask image ...\n", 0 );
      }
      ptrInputMaxMskImage = inputMaxMskImage;
    }
  }

  if ( par.input_segmentation_name != (char*)NULL && par.input_segmentation_name[0] != '\0' ) {
    inputSegImage = _VT_Inrimage( par.input_segmentation_name );
    if ( inputSegImage == (vt_image*)NULL ) {
      if ( inputMaxMskImage != (vt_image*)NULL ) {
        VT_FreeImage( inputMaxMskImage );
        VT_Free( (void**)&inputMaxMskImage );
      }
      if ( inputMinMskImage != (vt_image*)NULL ) {
        VT_FreeImage( inputMinMskImage );
        VT_Free( (void**)&inputMinMskImage );
      }
      if ( inputRecImage != (vt_image*)NULL ) {
        VT_FreeImage( inputRecImage );
        VT_Free( (void**)&inputRecImage );
      }
      VT_FreeImage( inputIntImage );
      VT_Free( (void**)&inputIntImage );
      API_ErrorParse_adhocFuse( _BaseName( argv[0] ), "unable to read input segmentation image ...\n", 0 );
    }
  }


  /* output image creation
   */

  VT_Image( &outputIntImage );
  ptrOutputIntImage = (vt_image*)NULL;

  if ( par.output_minimum_name != (char*)NULL && par.output_minimum_name[0] != '\0' ) {
    VT_InitFromImage( &outputMinImage, inputIntImage, par.output_minimum_name, inputIntImage->type );

    if ( VT_AllocImage( &outputMinImage ) != 1 ) {
      if ( inputSegImage != (vt_image*)NULL ) {
          VT_FreeImage( inputSegImage );
          VT_Free( (void**)&inputSegImage );
      }
      if ( inputMaxMskImage != (vt_image*)NULL ) {
        VT_FreeImage( inputMaxMskImage );
        VT_Free( (void**)&inputMaxMskImage );
      }
      if ( inputMinMskImage != (vt_image*)NULL ) {
        VT_FreeImage( inputMinMskImage );
        VT_Free( (void**)&inputMinMskImage );
      }
      if ( inputRecImage != (vt_image*)NULL ) {
        VT_FreeImage( inputRecImage );
        VT_Free( (void**)&inputRecImage );
      }
      VT_FreeImage( inputIntImage );
      VT_Free( (void**)&inputIntImage );
      API_ErrorParse_adhocFuse( _BaseName( argv[0] ), "unable to allocate output minimum image ...\n", 0 );
    }
    ptrOutputMinImage = &outputMinImage;
  }

  if ( par.output_maximum_name != (char*)NULL && par.output_maximum_name[0] != '\0' ) {
    VT_InitFromImage( &outputMaxImage, inputIntImage, par.output_maximum_name, inputIntImage->type );

    if ( VT_AllocImage( &outputMaxImage ) != 1 ) {
      if ( ptrOutputMinImage != (vt_image*)NULL ) {
        VT_FreeImage( ptrOutputMinImage );
        ptrOutputMinImage = (vt_image*)NULL;
      }
      if ( inputSegImage != (vt_image*)NULL ) {
          VT_FreeImage( inputSegImage );
          VT_Free( (void**)&inputSegImage );
      }
      if ( inputMaxMskImage != (vt_image*)NULL ) {
        VT_FreeImage( inputMaxMskImage );
        VT_Free( (void**)&inputMaxMskImage );
      }
      if ( inputMinMskImage != (vt_image*)NULL ) {
        VT_FreeImage( inputMinMskImage );
        VT_Free( (void**)&inputMinMskImage );
      }
      if ( inputRecImage != (vt_image*)NULL ) {
        VT_FreeImage( inputRecImage );
        VT_Free( (void**)&inputRecImage );
      }
      VT_FreeImage( inputIntImage );
      VT_Free( (void**)&inputIntImage );
      API_ErrorParse_adhocFuse( _BaseName( argv[0] ), "unable to allocate output maximum image ...\n", 0 );
    }
    ptrOutputMaxImage = &outputMaxImage;
  }

  if ( par.output_intensity_name != (char*)NULL && par.output_intensity_name[0] != '\0' ) {
    VT_InitFromImage( &outputIntImage, inputIntImage, par.output_intensity_name, outputIntImageType );
    if ( par.output_type != TYPE_UNKNOWN ) outputIntImage.type = par.output_type;

    if ( VT_AllocImage( &outputIntImage ) != 1 ) {
      if ( ptrOutputMaxImage != (vt_image*)NULL ) {
        VT_FreeImage( ptrOutputMaxImage );
        ptrOutputMaxImage = (vt_image*)NULL;
      }
      if ( ptrOutputMinImage != (vt_image*)NULL ) {
        VT_FreeImage( ptrOutputMinImage );
        ptrOutputMinImage = (vt_image*)NULL;
      }
      if ( inputSegImage != (vt_image*)NULL ) {
          VT_FreeImage( inputSegImage );
          VT_Free( (void**)&inputSegImage );
      }
      if ( inputMaxMskImage != (vt_image*)NULL ) {
        VT_FreeImage( inputMaxMskImage );
        VT_Free( (void**)&inputMaxMskImage );
      }
      if ( inputMinMskImage != (vt_image*)NULL ) {
        VT_FreeImage( inputMinMskImage );
        VT_Free( (void**)&inputMinMskImage );
      }
      if ( inputRecImage != (vt_image*)NULL ) {
        VT_FreeImage( inputRecImage );
        VT_Free( (void**)&inputRecImage );
      }
      VT_FreeImage( inputIntImage );
      VT_Free( (void**)&inputIntImage );
      API_ErrorParse_adhocFuse( _BaseName( argv[0] ), "unable to allocate output renormalized image ...\n", 0 );
    }
    ptrOutputIntImage = &outputIntImage;
  }

  VT_Image( &outputFusImage );
  ptrOutputFusImage = (vt_image*)NULL;

  if ( inputRecImage != (vt_image*)NULL
       && par.output_fusion_name != (char*)NULL
       && par.output_fusion_name[0] != '\0' ) {
    VT_InitFromImage( &outputFusImage, inputRecImage, par.output_fusion_name, inputRecImage->type );

    if ( VT_AllocImage( &outputFusImage ) != 1 ) {
      if ( ptrOutputIntImage != (vt_image*)NULL ) {
        VT_FreeImage( ptrOutputIntImage );
        ptrOutputIntImage = (vt_image*)NULL;
      }
      if ( ptrOutputMaxImage != (vt_image*)NULL ) {
        VT_FreeImage( ptrOutputMaxImage );
        ptrOutputMaxImage = (vt_image*)NULL;
      }
      if ( ptrOutputMinImage != (vt_image*)NULL ) {
        VT_FreeImage( ptrOutputMinImage );
        ptrOutputMinImage = (vt_image*)NULL;
      }
      if ( inputSegImage != (vt_image*)NULL ) {
          VT_FreeImage( inputSegImage );
          VT_Free( (void**)&inputSegImage );
      }
      if ( inputMaxMskImage != (vt_image*)NULL ) {
        VT_FreeImage( inputMaxMskImage );
        VT_Free( (void**)&inputMaxMskImage );
      }
      if ( inputMinMskImage != (vt_image*)NULL ) {
        VT_FreeImage( inputMinMskImage );
        VT_Free( (void**)&inputMinMskImage );
      }
      if ( inputRecImage != (vt_image*)NULL ) {
        VT_FreeImage( inputRecImage );
        VT_Free( (void**)&inputRecImage );
      }
      VT_FreeImage( inputIntImage );
      VT_Free( (void**)&inputIntImage );
      API_ErrorParse_adhocFuse( _BaseName( argv[0] ), "unable to allocate output merged image ...\n", 0 );
    }
    ptrOutputFusImage = &outputFusImage;
  }





  /* API call
   */

  lineoptions = _Array2Str( argc, argv );
  if ( lineoptions == (char*)NULL ) {
    if ( ptrOutputFusImage != (vt_image*)NULL ) {
      VT_FreeImage( ptrOutputFusImage );
      ptrOutputFusImage = (vt_image*)NULL;
    }
    if ( ptrOutputIntImage != (vt_image*)NULL ) {
      VT_FreeImage( ptrOutputIntImage );
      ptrOutputIntImage = (vt_image*)NULL;
    }
    if ( ptrOutputMaxImage != (vt_image*)NULL ) {
      VT_FreeImage( ptrOutputMaxImage );
      ptrOutputMaxImage = (vt_image*)NULL;
    }
    if ( ptrOutputMinImage != (vt_image*)NULL ) {
      VT_FreeImage( ptrOutputMinImage );
      ptrOutputMinImage = (vt_image*)NULL;
    }
    if ( inputSegImage != (vt_image*)NULL ) {
        VT_FreeImage( inputSegImage );
        VT_Free( (void**)&inputSegImage );
    }
    if ( inputMaxMskImage != (vt_image*)NULL ) {
      VT_FreeImage( inputMaxMskImage );
      VT_Free( (void**)&inputMaxMskImage );
    }
    if ( inputMinMskImage != (vt_image*)NULL ) {
      VT_FreeImage( inputMinMskImage );
      VT_Free( (void**)&inputMinMskImage );
    }
    if ( inputRecImage != (vt_image*)NULL ) {
      VT_FreeImage( inputRecImage );
      VT_Free( (void**)&inputRecImage );
    }
   VT_FreeImage( inputIntImage );
    VT_Free( (void**)&inputIntImage );
    API_ErrorParse_adhocFuse( _BaseName( argv[0] ), "unable to translate command line options ...\n", 0 );
  }



  if ( API_adhocFuse( inputIntImage, inputRecImage,
                      ptrInputMinMskImage, ptrInputMaxMskImage,
                      inputSegImage,
                      ptrOutputMinImage, ptrOutputMaxImage,
                      ptrOutputIntImage, ptrOutputFusImage, lineoptions, (char*)NULL ) != 1 ) {
      vtfree( lineoptions );
      if ( ptrOutputFusImage != (vt_image*)NULL ) {
        VT_FreeImage( ptrOutputFusImage );
        ptrOutputFusImage = (vt_image*)NULL;
      }
      if ( ptrOutputIntImage != (vt_image*)NULL ) {
        VT_FreeImage( ptrOutputIntImage );
        ptrOutputIntImage = (vt_image*)NULL;
      }
      if ( ptrOutputMaxImage != (vt_image*)NULL ) {
        VT_FreeImage( ptrOutputMaxImage );
        ptrOutputMaxImage = (vt_image*)NULL;
      }
      if ( ptrOutputMinImage != (vt_image*)NULL ) {
        VT_FreeImage( ptrOutputMinImage );
        ptrOutputMinImage = (vt_image*)NULL;
      }
      if ( inputSegImage != (vt_image*)NULL ) {
          VT_FreeImage( inputSegImage );
          VT_Free( (void**)&inputSegImage );
      }
      if ( inputMaxMskImage != (vt_image*)NULL ) {
        VT_FreeImage( inputMaxMskImage );
        VT_Free( (void**)&inputMaxMskImage );
      }
      if ( inputMinMskImage != (vt_image*)NULL ) {
        VT_FreeImage( inputMinMskImage );
        VT_Free( (void**)&inputMinMskImage );
      }
      if ( inputRecImage != (vt_image*)NULL ) {
        VT_FreeImage( inputRecImage );
        VT_Free( (void**)&inputRecImage );
      }
      VT_FreeImage( inputIntImage );
      VT_Free( (void**)&inputIntImage );
      API_ErrorParse_adhocFuse( _BaseName( argv[0] ), "some error occurs during processing ...\n", -1 );
  }
  vtfree( lineoptions );



  /* freeing input images
   */
  if ( inputSegImage != (vt_image*)NULL ) {
      VT_FreeImage( inputSegImage );
      VT_Free( (void**)&inputSegImage );
  }
  if ( inputMaxMskImage != (vt_image*)NULL ) {
    VT_FreeImage( inputMaxMskImage );
    VT_Free( (void**)&inputMaxMskImage );
  }
  if ( inputMinMskImage != (vt_image*)NULL ) {
    VT_FreeImage( inputMinMskImage );
    VT_Free( (void**)&inputMinMskImage );
  }
  if ( inputRecImage != (vt_image*)NULL ) {
    VT_FreeImage( inputRecImage );
    VT_Free( (void**)&inputRecImage );
  }
  VT_FreeImage( inputIntImage );
  VT_Free( (void**)&inputIntImage );



  /* output image writing
   */

  if ( ptrOutputMinImage != (vt_image*)NULL ) {
      if ( VT_WriteInrimage( ptrOutputMinImage ) == -1 ) {
          if ( ptrOutputFusImage != (vt_image*)NULL ) {
            VT_FreeImage( ptrOutputFusImage );
            ptrOutputFusImage = (vt_image*)NULL;
          }
          if ( ptrOutputIntImage != (vt_image*)NULL ) {
            VT_FreeImage( ptrOutputIntImage );
            ptrOutputIntImage = (vt_image*)NULL;
          }
          if ( ptrOutputMaxImage != (vt_image*)NULL ) {
            VT_FreeImage( ptrOutputMaxImage );
            ptrOutputMaxImage = (vt_image*)NULL;
          }
          if ( ptrOutputMinImage != (vt_image*)NULL ) {
            VT_FreeImage( ptrOutputMinImage );
            ptrOutputMinImage = (vt_image*)NULL;
          }
          API_ErrorParse_adhocFuse( _BaseName( argv[0] ), "unable to write output minimum image ...\n", -1 );
      }
  }

  if ( ptrOutputMinImage != (vt_image*)NULL ) {
    VT_FreeImage( ptrOutputMinImage );
    ptrOutputMinImage = (vt_image*)NULL;
  }

  if ( ptrOutputMaxImage != (vt_image*)NULL ) {
      if ( VT_WriteInrimage( ptrOutputMaxImage ) == -1 ) {
          if ( ptrOutputFusImage != (vt_image*)NULL ) {
            VT_FreeImage( ptrOutputFusImage );
            ptrOutputFusImage = (vt_image*)NULL;
          }
          if ( ptrOutputIntImage != (vt_image*)NULL ) {
            VT_FreeImage( ptrOutputIntImage );
            ptrOutputIntImage = (vt_image*)NULL;
          }
          if ( ptrOutputMaxImage != (vt_image*)NULL ) {
            VT_FreeImage( ptrOutputMaxImage );
            ptrOutputMaxImage = (vt_image*)NULL;
          }
          API_ErrorParse_adhocFuse( _BaseName( argv[0] ), "unable to write output maximum image ...\n", -1 );
      }
  }

  if ( ptrOutputMaxImage != (vt_image*)NULL ) {
    VT_FreeImage( ptrOutputMaxImage );
    ptrOutputMaxImage = (vt_image*)NULL;
  }

  if ( ptrOutputIntImage != (vt_image*)NULL ) {
    if ( VT_WriteInrimage( ptrOutputIntImage ) == -1 ) {
        if ( ptrOutputFusImage != (vt_image*)NULL ) {
          VT_FreeImage( ptrOutputFusImage );
          ptrOutputFusImage = (vt_image*)NULL;
        }
        if ( ptrOutputIntImage != (vt_image*)NULL ) {
          VT_FreeImage( ptrOutputIntImage );
          ptrOutputIntImage = (vt_image*)NULL;
        }
        API_ErrorParse_adhocFuse( _BaseName( argv[0] ), "unable to write output intensity image ...\n", -1 );
    }
  }

  if ( ptrOutputIntImage != (vt_image*)NULL ) {
    VT_FreeImage( ptrOutputIntImage );
    ptrOutputIntImage = (vt_image*)NULL;
  }

  if ( ptrOutputFusImage != (vt_image*)NULL ) {
    if ( VT_WriteInrimage( ptrOutputFusImage ) == -1 ) {
        if ( ptrOutputFusImage != (vt_image*)NULL ) {
          VT_FreeImage( ptrOutputFusImage );
          ptrOutputFusImage = (vt_image*)NULL;
        }
        API_ErrorParse_adhocFuse( _BaseName( argv[0] ), "unable to write output fused image ...\n", -1 );
    }
  }

  if ( ptrOutputFusImage != (vt_image*)NULL ) {
    VT_FreeImage( ptrOutputFusImage );
    ptrOutputFusImage = (vt_image*)NULL;
  }


  /* ... */



  if ( par.trace_allocations ) {
    fprintfVtMallocTrace( stderr );
    clearVtMalloc();
  }

  time_exit = _GetTime();
  clock_exit = _GetClock();

  if ( par.print_time ) { 
    fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName( argv[0] ), time_exit - time_init );
    fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
    fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );
  }


  return( 0 );
}





/************************************************************
 *
 * static functions
 *
 ************************************************************/



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
