/*************************************************************************
 * api-TVmembrane.c -
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

#include <chunks.h>
#include <vtmalloc.h>

#include <vt_common.h>

#include <api-TVmembrane.h>






static int _verbose_ = 1;
static int _debug_ = 0;


static void _API_ParseParam_TVmembrane( char *str, lineCmdParamTVmembrane *p );



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



/*
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
*/





static int _count_points( void *theBuf, bufferType theBufType, int *theDim )
{
  char *proc = "_count_points";
  int i, v = theDim[0]*theDim[1]*theDim[2];
  int n = 0;

  switch ( theBufType ) {

  default :
    if ( _VT_VERBOSE_ ) {
      fprintf( stderr, "%s: such label type not handled yet\n", proc );
    }
    return( -1 );

  case UCHAR :
    {
      u8 *buf = (u8*)theBuf;
      for ( i=0; i<v; i++ )
        if ( buf[i] ) n++;
    }
    break;

  case USHORT :
    {
      u16 *buf = (u16*)theBuf;
      for ( i=0; i<v; i++ )
        if ( buf[i] ) n++;
    }
    break;

  case FLOAT :
    {
      r32 *buf = (r32*)theBuf;
      for ( i=0; i<v; i++ )
        if ( buf[i] != 0 ) n++;
    }
    break;

  }

  return( n );
}







/************************************************************
 *
 * main APIs
 *
 ************************************************************/


int API_Sampling(vt_image *imageBin,
                 char *param_str_1,
                 char *param_str_2 )

{

    char *proc = "API_Sampling";
    lineCmdParamTVmembrane par;

    double time_sampling;
    double clock_sampling;
    int flag_3D=1;

    int nparcels;
    int **theSeeds = NULL;
    int *seeds = NULL;
    int theDim[3];
    int i,j,k,n;
    char Ftime[DOUBLESTRINGLENGTH];

    /* parameter initialization
     */
    API_InitParam_TVmembrane( &par );

    /* parameter parsing
     */
    if ( param_str_1 != (char*)NULL )
        _API_ParseParam_TVmembrane( param_str_1, &par );
    if ( param_str_2 != (char*)NULL )
        _API_ParseParam_TVmembrane( param_str_2, &par );

    if ( par.print_lineCmdParam )
        API_PrintParam_TVmembrane( stderr, proc, &par, (char*)NULL );

    /************************************************************
     *
     *  here is the stuff
     *
     ************************************************************/

    if ( par.dimension == 2 || imageBin->dim.z == 1 )
      flag_3D = 0;


    /* Echantillonnage de l'image */
    if (par.power == 1)
    {
      par.sample = exp((par.sample>=0 ? - par.sample : par.sample)*log(10));
    }

    if ( par.sample < 1.0 )
    {
      time_sampling = _GetTime();
      clock_sampling = _GetClock();
      if ( _VT_VERBOSE_ )
        fprintf(stdout, "Sampling step : sampling coefficient = %f\n", par.sample);
      if (par.sample > 0 && par.sample < 1)
      {
        if (flag_3D == 1)
        {
          if (par.sampleMode == _RANDOM_SAMPLING_)
          {
            if (MT_SampleBin(imageBin, par.sample) != 1 )
            {
              API_ErrorParse_TVmembrane( (char*)NULL, "unexpected error while sampling binary image\n", 0);
            }
          }
          if (par.sampleMode == _REGULAR_SAMPLING_)
          {

            theDim[0] = imageBin->dim.x;
            theDim[1] = imageBin->dim.y;
            theDim[2] = imageBin->dim.z;

            nparcels =(int) (((double)_count_points(imageBin->buf, imageBin->type, theDim)) * par.sample);

            /* fprintf(stdout, "nparcels=%d\n", nparcels); */

            seeds = (int*)vtmalloc( 3*nparcels * sizeof( int ), "seeds", proc );
            theSeeds = (int**)vtmalloc( nparcels * sizeof( int* ), "theSeeds", proc );
            if ( seeds == NULL || theSeeds == NULL ) {
              if ( seeds != NULL ) vtfree( seeds );
              if ( theSeeds != NULL ) vtfree( theSeeds );
              API_ErrorParse_TVmembrane( (char*)NULL, "error when allocating seeds arrays\n", 0);
            }

            for ( n=0;n<nparcels;n++ ) {
              theSeeds[n] = &(seeds[3*n]);
            }

            if ( parcelling( imageBin->buf, imageBin->type,
                             theSeeds, nparcels,
                             NULL, TYPE_UNKNOWN,
                             NULL, TYPE_UNKNOWN,
                             theDim, 0, NULL ) != 1 ) {
              if ( seeds != NULL ) vtfree( seeds );
              if ( theSeeds != NULL ) vtfree( theSeeds );
              API_ErrorParse_TVmembrane( (char*)NULL, "error when processing\n", 0);
            }
            switch(imageBin->type)
            {
            case UCHAR :
              {
                u8 ***array = (unsigned char ***)imageBin->array;
                for (i=0;i<theDim[0];i++)
                for (j=0;j<theDim[1];j++)
                for (k=0;k<theDim[2];k++)
                {
                  array[k][j][i]= (unsigned char) 0 ;
                }
                for (n=0;n<nparcels;n++)
                {
                  array[theSeeds[n][2]][theSeeds[n][1]][theSeeds[n][0]]=(unsigned char) 255;
                }
                break;
              }
            case FLOAT :
              {
                r32 ***array = (float ***)imageBin->array;
                for (n=0;n<nparcels;n++)
                {
                  array[theSeeds[n][2]][theSeeds[n][1]][theSeeds[n][0]]*= -1;
                }
                for (i=0;i<theDim[0];i++)
                for (j=0;j<theDim[1];j++)
                for (k=0;k<theDim[2];k++)
                {
                  array[k][j][i]= (array[k][j][i]>0) ? 0.0 : -array[k][j][i];
                }
                break;
              }
            default:
              if ( seeds != NULL ) vtfree( seeds );
              if ( theSeeds != NULL ) vtfree( theSeeds );
              API_ErrorParse_TVmembrane( (char*)NULL, "error: image type not implemented yet\n", 0);
            }
            if ( seeds != NULL ) vtfree( seeds );
            if ( theSeeds != NULL ) vtfree( theSeeds );
          }
        }
        else
        {
          if (MT_SampleBin2D(imageBin, par.sample) != 1 )
          {
            API_ErrorParse_TVmembrane( (char*)NULL, "unexpected error while sampling binary image\n", 0);
          }
        }
      }
      else
      {
        API_ErrorParse_TVmembrane( (char*)NULL, "incorrect sampling parameter : 0 < expected value <= 1\n", 0);
      }

      double time_exit = _GetTime();
      double clock_exit = _GetClock();
      if ( par.print_time ) {

        fprintf( stderr, "%s: elapsed (real) time = %f\n", proc, time_exit - time_sampling );
        fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_sampling );
        fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_sampling)/(time_exit - time_sampling) );

        if ( par.outputGeneric != (char*)NULL && par.outputGeneric[0] != '\0' ) {
          sprintf( Ftime, "%s.timesampling.txt", par.outputGeneric );
          FILE *fichier = fopen (Ftime, "a" );
          if ( fichier == NULL ) {
            perror (Ftime);
          }
          else {
              fprintf(fichier, "\n");
              fprintf( fichier, "%s ", param_str_1);
              fprintf( fichier, "%s ", param_str_2);
              fprintf(fichier, "\n");
              fprintf(fichier, "%s: elapsed (real) time = %f\n", proc, time_exit - time_sampling );
              fprintf(fichier, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_sampling );
              fprintf(fichier, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_sampling)/(time_exit - time_sampling) );
              fclose(fichier);
          }
        }
      }


    }
    return(1);
}





int API_TVmembrane( vt_image **imagesIn, vt_3Dtensor *theTensor, char *param_str_1, char *param_str_2 )
{
  char *proc = "API_TVmembrane";
  lineCmdParamTVmembrane par;
  double _siz=1.0;
  double zfact;
  /* parameter initialization
   */
  API_InitParam_TVmembrane( &par );

  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_TVmembrane( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_TVmembrane( param_str_2, &par );

  if ( par.print_lineCmdParam )
      API_PrintParam_TVmembrane( stderr, proc, &par, (char*)NULL );

  /************************************************************
   *
   *  here is the stuff
   *
   ************************************************************/
  zfact = par.zfact;
  if (par.flagReal == 1)
  {
      _siz=imagesIn[0]->siz.x;
      zfact = imagesIn[0]->siz.x/imagesIn[0]->siz.z;
      if (_VT_VERBOSE_)
      {
          fprintf(stdout, "Voting scale %f given in real coordinates\n", par.scale );
      }
  }
  else
  {
      if (_VT_VERBOSE_)
      {
          fprintf(stdout, "Voting scale %f given in voxel coordinates\n", par.scale );
      }

  }
  if(_VT_VERBOSE_)
      fprintf(stdout, "Entree dans la fonction de tensor voting\n");

  if ( MT_Compute3DTensorVoting( theTensor, imagesIn, par.scale/_siz, zfact, par.niter,
                                 par.nangles,par.nsticks, par.TVmode, par.initHessian,
                                 par.outputGeneric, par.writeImages) != 1 )
  {
      API_ErrorParse_TVmembrane( (char*)NULL, "problem while computing the tensor voting\n", 0 );
      return(-1);
  }


  if(_VT_VERBOSE_)
    fprintf(stdout, "Retour dans la fonction principale\n");


  return( 1 );
}





int API_TVmembrane2D( vt_image **imagesIn, mt_2Dtensor *theTensor2D, char *param_str_1, char *param_str_2 )
{
  char *proc = "API_TVmembrane";
  lineCmdParamTVmembrane par;
  double _siz=1.0;
  /* parameter initialization
   */
  API_InitParam_TVmembrane( &par );

  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_TVmembrane( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_TVmembrane( param_str_2, &par );

  if ( par.print_lineCmdParam )
      API_PrintParam_TVmembrane( stderr, proc, &par, (char*)NULL );

  /************************************************************
   *
   *  here is the stuff
   *
   ************************************************************/
  if (par.flagReal == 1)
  {
      _siz=imagesIn[0]->siz.x;
      if (_VT_VERBOSE_)
      {
          fprintf(stdout, "Voting scale %f given in real coordinates\n", par.scale );
      }
  }
  else
      if (_VT_VERBOSE_)
        fprintf(stdout, "Voting scale %f given in voxel coordinates\n", par.scale );


  if(_VT_VERBOSE_)
      fprintf(stdout, "Entree dans la fonction de tensor voting\n");


  if ( MT_Compute2DTensorVoting( theTensor2D, imagesIn, par.scale/_siz, par.niter,
                                par.nangles, par.TVmode, par.initHessian,
                                par.outputGeneric, par.writeImages) != 1 )
  {
      API_ErrorParse_TVmembrane( (char*)NULL, "problem while computing the tensor voting 2D\n", 0 );
      return(-1);
  }


  if(_VT_VERBOSE_)
    fprintf(stdout, "Retour dans la fonction principale\n");


  return( 1 );
}





int API_ShapeExtraction( vt_3Dtensor *theTensor3D, char *param_str_1, char *param_str_2 )
{

    char *proc = "API_ShapeExtraction";
    lineCmdParamTVmembrane par;
    char name[DOUBLESTRINGLENGTH];
    char *ptrSuffix;
    char *inrSuffix = "inr";

    vt_image theExtrema;
    double zfact;

    /* parameter initialization
     */
    API_InitParam_TVmembrane( &par );

    /* parameter parsing
     */
    if ( param_str_1 != (char*)NULL )
        _API_ParseParam_TVmembrane( param_str_1, &par );
    if ( param_str_2 != (char*)NULL )
        _API_ParseParam_TVmembrane( param_str_2, &par );

    if ( par.print_lineCmdParam )
        API_PrintParam_TVmembrane( stderr, proc, &par, (char*)NULL );

    if ( par.shape == NO_EXTRACTION ) {
        return(1);
    }

    /************************************************************
     *
     *  here is the stuff
     *
     ************************************************************/
    zfact = par.zfact;
    if ( par.flagReal == 1 ) {
        zfact = theTensor3D->imvp1.siz.z/theTensor3D->imvp1.siz.x;
        if ( _verbose_ ) {
            fprintf(stdout, "Voting scale %f given in real coordinates\n", par.scale );
        }
    }

    VT_InitFromImage( &theExtrema, &(theTensor3D->imxx), (char*)NULL, FLOAT );
    if ( VT_AllocImage( &theExtrema ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate extrema image\n", proc );
      return( -1 );
    }

    ptrSuffix = ( par.outputSuffix != (char*)NULL && par.outputSuffix[0] != '\0' ) ? par.outputSuffix : inrSuffix;



    /* 1. Surfaces extraction
     */

    if ( par.shape == MEMBRANE || par.shape == MV || par.shape == MB ||
         par.shape == MVB ) {

        MT_ComputeTensorSurfaceExtrema( theTensor3D, &theExtrema, zfact );

        /* build the name
         */
        if ( par.outputGeneric != (char*)NULL && par.outputGeneric[0] != '\0' ) {
          sprintf( name, "%s.plane.%s", par.outputGeneric, ptrSuffix );
        }
        else if ( par.outputExtrema != (char*)NULL && par.outputExtrema[0] != '\0' ) {
          if ( par.shape == MEMBRANE )
            sprintf( name, "%s", par.outputExtrema );
          else
            sprintf( name, "%s.plane.%s", par.outputExtrema, ptrSuffix );
        }
        else {
          sprintf( name, "plane.%s", ptrSuffix );
        }

        if ( _verbose_ )
          fprintf(stdout, "%s: writing plane extraction into '%s' ...\n", proc, name );

        if ( VT_WriteInrimageWithName( &theExtrema, name ) != 1 ) {
          VT_FreeImage( &theExtrema );
          if ( _verbose_ )
            fprintf( stderr, "%s: error when writing '%s'\n", proc, name );
          return( -1 );
        }
    }



    /* 2. Lines extraction
     */
    if ( par.shape == VESSEL || par.shape == MV || par.shape == VB ||
         par.shape == MVB ) {

          MT_ComputeTensorLineExtrema( theTensor3D, &theExtrema );

          /* build the name
           */
          if ( par.outputGeneric != (char*)NULL && par.outputGeneric[0] != '\0' ) {
            sprintf( name, "%s.line.inr", par.outputGeneric );
          }
          else if ( par.outputExtrema != (char*)NULL && par.outputExtrema[0] != '\0' ) {
            if ( par.shape == VESSEL )
              sprintf( name, "%s", par.outputExtrema );
            else
              sprintf( name, "%s.line.%s", par.outputExtrema, ptrSuffix );
          }
          else {
            sprintf( name, "line.%s", ptrSuffix );
          }

          if ( _verbose_ )
            fprintf(stdout, "%s: writing line extraction into '%s' ...\n", proc, name );

          if ( VT_WriteInrimageWithName( &theExtrema, name ) != 1 ) {
            VT_FreeImage( &theExtrema );
            if ( _verbose_ )
              fprintf( stderr, "%s: error when writing '%s'\n", proc, name );
            return( -1 );
          }
    }


    /* 3. Balls extraction
     */
    if ( par.shape == BALL || par.shape == MB || par.shape == VB ||
         par.shape == MVB ) {

        MT_ComputeTensorBallExtrema( theTensor3D, &theExtrema );

        /* build the name
         */
        if ( par.outputGeneric != (char*)NULL && par.outputGeneric[0] != '\0' ) {
          sprintf( name, "%s.ball.%s", par.outputGeneric, ptrSuffix );
        }
        else if ( par.outputExtrema != (char*)NULL && par.outputExtrema[0] != '\0' ) {
          if ( par.shape == BALL )
            sprintf( name, "%s", par.outputExtrema );
          else
            sprintf( name, "%s.ball.%s", par.outputExtrema, ptrSuffix );
        }
        else {
          sprintf( name, "ball.%s", ptrSuffix );
        }

        if ( _verbose_ )
          fprintf(stdout, "%s: writing ball extraction into '%s' ...\n", proc, name );

        if ( VT_WriteInrimageWithName( &theExtrema, name ) != 1 ) {
          VT_FreeImage( &theExtrema );
          if ( _verbose_ )
            fprintf( stderr, "%s: error when writing '%s'\n", proc, name );
          return( -1 );
        }
    }

    VT_FreeImage(&theExtrema);

    return(1);
}




int API_ShapeExtraction2D( mt_2Dtensor *theTensor2D, char *param_str_1, char *param_str_2 )
{

    char *proc = "API_ShapeExtraction2D";
    lineCmdParamTVmembrane par;
    char name[DOUBLESTRINGLENGTH];
    char *ptrSuffix;
    char *inrSuffix = "inr";

    vt_image theExtrema;


    /* parameter initialization
     */
    API_InitParam_TVmembrane( &par );

    /* parameter parsing
     */
    if ( param_str_1 != (char*)NULL )
        _API_ParseParam_TVmembrane( param_str_1, &par );
    if ( param_str_2 != (char*)NULL )
        _API_ParseParam_TVmembrane( param_str_2, &par );

    if ( par.print_lineCmdParam )
        API_PrintParam_TVmembrane( stderr, proc, &par, (char*)NULL );

    if ( par.shape == NO_EXTRACTION ) {
        return(1);
    }

    /************************************************************
     *
     *  here is the stuff
     *
     ************************************************************/

    VT_InitFromImage( &theExtrema, &(theTensor2D->imxx), (char*)NULL, FLOAT );
    if ( VT_AllocImage( &theExtrema ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate extrema image\n", proc );
      return( -1 );
    }

    ptrSuffix = ( par.outputSuffix != (char*)NULL && par.outputSuffix[0] != '\0' ) ? par.outputSuffix : inrSuffix;



    /* 1. Surfaces=Lines extraction
     */

    if ( par.shape == MEMBRANE || par.shape == MV || par.shape == MB ||
         par.shape == MVB || par.shape == VESSEL || par.shape == VB ) {

        MT_Compute2DTensorLineicExtrema( theTensor2D, &theExtrema );

        /* build the name
         */
        if ( par.outputGeneric != (char*)NULL && par.outputGeneric[0] != '\0' ) {
          sprintf( name, "%s.line.%s", par.outputGeneric, ptrSuffix );
        }
        else if ( par.outputExtrema != (char*)NULL && par.outputExtrema[0] != '\0' ) {
          if ( par.shape == MEMBRANE || par.shape == VESSEL )
            sprintf( name, "%s", par.outputExtrema );
          else
            sprintf( name, "%s.line.%s", par.outputExtrema, ptrSuffix );
        }
        else {
          sprintf( name, "line.%s", ptrSuffix );
        }

        if ( _verbose_ )
          fprintf(stdout, "%s: writing line extraction into '%s' ...\n", proc, name );

        if ( VT_WriteInrimageWithName( &theExtrema, name ) != 1 ) {
          VT_FreeImage( &theExtrema );
          if ( _verbose_ )
            fprintf( stderr, "%s: error when writing '%s'\n", proc, name );
          return( -1 );
        }
    }



    /* 2. Balls extraction
     */
    if ( par.shape == BALL || par.shape == MB || par.shape == VB ||
         par.shape == MVB ) {

        MT_Compute2DTensorBallExtrema( theTensor2D, &theExtrema );

        /* build the name
         */
        if ( par.outputGeneric != (char*)NULL && par.outputGeneric[0] != '\0' ) {
          sprintf( name, "%s.ball.%s", par.outputGeneric, ptrSuffix );
        }
        else if ( par.outputExtrema != (char*)NULL && par.outputExtrema[0] != '\0' ) {
          if ( par.shape == BALL )
            sprintf( name, "%s", par.outputExtrema );
          else
            sprintf( name, "%s.ball.%s", par.outputExtrema, ptrSuffix );
        }
        else {
          sprintf( name, "ball.%s", ptrSuffix );
        }

        if ( _verbose_ )
          fprintf(stdout, "%s: writing ball extraction into '%s' ...\n", proc, name );

        if ( VT_WriteInrimageWithName( &theExtrema, name ) != 1 ) {
          VT_FreeImage( &theExtrema );
          if ( _verbose_ )
            fprintf( stderr, "%s: error when writing '%s'\n", proc, name );
          return( -1 );
        }
    }

    VT_FreeImage(&theExtrema);

    return(1);
}











/************************************************************
 *
 * static functions
 *
 ************************************************************/



static char **_Str2Array( int *argc, char *str )
{
  char *proc = "_Str2Array";
  int n = 0;
  char *s = str;
  char **array, **a;

  if ( s == (char*)NULL || strlen( s ) == 0 ) {
    if ( _verbose_ >= 2 )
      fprintf( stderr, "%s: empty input string\n", proc );
    *argc = 0;
    return( (char**)NULL );
  }

  /* go to the first valid character
   */
  while ( *s == ' ' || *s == '\n' || *s == '\t' )
    s++;

  if ( *s == '\0' ) {
    if ( _verbose_ >= 2 )
      fprintf( stderr, "%s: weird, input string contains only separation characters\n", proc );
    *argc = 0;
    return( (char**)NULL );
  }

  /* count the number of strings
   */
  for ( n = 0; *s != '\0'; ) {
    n ++;
    while ( *s != ' ' && *s != '\n' && *s != '\t' && *s != '\0' )
      s ++;
    while ( *s == ' ' || *s == '\n' || *s == '\t' )
      s ++;
  }

  if ( _verbose_ >= 5 )
    fprintf( stderr, "%s: found %d strings\n", proc, n );

  /* the value of the strings will be duplicated
   * so that the input string can be freed
   */
  array = (char**)vtmalloc( n * sizeof(char*) + (strlen(str)+1) * sizeof(char), "array", proc );
  if ( array == (char**)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    *argc = 0;
    return( (char**)NULL );
  }

  a = array;
  a += n;
  s = (char*)a;
  (void)strncpy( s, str, strlen( str ) );
  s[ strlen( str ) ] = '\0';

  while ( *s == ' ' || *s == '\n' || *s == '\t' ) {
    *s = '\0';
    s++;
  }

  for ( n = 0; *s != '\0'; ) {
    array[n] = s;
    n ++;
    while ( *s != ' ' && *s != '\n' && *s != '\t' && *s != '\0' )
      s ++;
    while ( *s == ' ' || *s == '\n' || *s == '\t' ) {
      *s = '\0';
      s ++;
    }
  }

  *argc = n;
  return( array );
}





/************************************************************
 *
 * help / documentation
 *
 ************************************************************/



static char *usage = "image-in \n\
\t [-output-generic-prefix|-output-generic] prefix-out\n\
\t [-output-sampling sampled-image-out]\n\
\t [-output-tensor tensor-prefix-out]\n\
\t [-output-eigenvalues eigenvalues-prefix-out]\n\
\t [-output-angles angles-prefix-out]\n\
\t [-output-binary binary-prefix-out]\n\
\t [-output-extrema extrema-prefix-out]\n\
\t [-output-suffix suffix-out]\n\
\t [-mask %s]\n\
\t [-sample %lf] [-power]\n\
\t [-random] [-regular]\n\
\t [-sampling-iterations|-iterations|-si|-i %d]\n\
\t [-scale %lf] [-real]\n\
\t [-zfact %lf]\n\
\t [-hessian [plane|line]]\n\
\t [-niter %d]\n\
\t [-cftv|-tvclassic] [-nangles %d] [-nsticks %d]\n\
\t [-shape m|v|b|mv|mb|vb|mvb]\n\
\t [-write-images|-wi]\n\
\t [-2D]\n\
\t [-images-parallel-voting %d]\n\
\t [-parallel|-no-parallel] [-max-chunks %d]\n\
\t [-parallelism-type|-parallel-type default|none|openmp|omp|pthread|thread]\n\
\t [-omp-scheduling|-omps default|static|dynamic-one|dynamic|guided]\n\
\t [-verbose|-v] [-no-verbose|-noverbose|-nv]\n\
\t [-debug|-D] [-no-debug|-nodebug]\n\
\t [-print-parameters|-param]\n\
\t [-print-time|-time] [-no-time|-notime]\n\
\t [-trace-memory|-memory] [-no-memory|-nomemory]\n\
\t [-help|--help|-h|--h]\n";

/*[-cftv|-tvclassic] [-nangles %d] [-nsticks %d]*/

static char *detail = "\
si 'image-in' est '-', on prendra stdin (incompatible avec le mode '-hessian')\n\
  [-output-generic-prefix|-output-generic] generic-prefix-out\n\
    generic suffix for all the output below. Hence all the outputs are written:\n\
  -output-sampling sampled-image-out\n\
    specifies the sampled image file name\n\
    else it is generic-prefix-out.sample.suffix-out (requires -write-images)\n\
  -output-tensor tensor-prefix-out\n\
    tensor image file names are:\n\
    tensor-prefix-out.[imxx|imyy|imxy].suffix-out for both 2D and 3D cases\n\
    and tensor-prefix-out.[imzz|imxz|imyz].suffix-out (3D case only)\n\
  -output-eigenvalues eigenvalues-prefix-out\n\
    eigenvalues image file names are:\n\
    eigenvalues-prefix-out.[imvp1|imvp2].suffix-out for both 2D and 3D cases\n\
    and eigenvalues-prefix-out.imvp3.suffix-out (3D case only)\n\
  -output-angles angles-prefix-out\n\
    angles image file names are:\n\
    angles-prefix-out.[imtheta1|imtheta2].suffix-out for both 2D and 3D cases\n\
    and angles-prefix-out.[imtheta3|imphi1|imphi2|imphi3].suffix-out (3D case only)\n\
  -output-binary binary-prefix-out\n\
    binary image file name is:\n\
    binary-prefix-out.[iszero].suffix-out for both 2D and 3D cases\n\
  -output-extrema extrema-prefix-out\n\
    if the extrema to be extracted are from one type only\n\
    (ie -shape m|v|b), it specifies the output extrema image file name\n\
    else it specifies its prefix, and names are build with\n\
    extrema-prefix-out.[plane|line|ball].suffix-out\n\
  -output-suffix suffix-out\n\
    allows to specify the suffix (hence the format of ouput images) when\n\
    only prefixes are given. Default is 'inr';\n\
#\n\
  -mask : mask image applied to the input image to restrict the domain of tensor\n\
    voting tokens (must be same dimensions with input image)\n\
# sampling options (done before tensor voting)\n\
  -sample %lf : echantillonne les votants de l'image initiale selon le coefficient (defaut : 1.0)\n\
  -power : l'argument de -sample est alors converti en puissance pour \n\
    l'echantillonnage : pourcentage=10^{- sample-coefficient}\n\
    eg: '-sample 2 -power' is equivalent to '-sample 0.01' and means that 1%% of points are kept\n\
  -random: random sampling (default)\n\
  -regular: regular sampling (a la Voronoi)\n\
  -sampling-iterations|-iterations|-si|-i %d:\n\
    maximal iterations number for regular sampling\n\
# tensor voting options\n\
  -scale %lf : echelle du tensor voting\n\
  -real : parametre d'echelle en coordonnees reelles si l'option est activee\n\
    (defaut : echelle en voxels)\n\
  -zfact %lf : facteur de resolution selon l'axe des z (defaut = 1, ie resolution isotrope).\n\
    Option a n'utiliser qu'en cas d'echelle en voxels (donc incompatible avec l'option -real).\n\
  -hessian : recuperation des informations directionnelles donnees par\n\
    la matrice hessienne calculee lors de la 1ere binarisation\n\
    (donner le prefixe ou le binaire en argument, necessite imtheta (+ imphi en 3D)\n\
    avec le meme prefixe que le binaire, formats supportes .inr et .hdr)\n\
  -hessian plane|line : plane si les directions sont des normales\n\
    aux structures (surfaces, defaut) ou line pour des tangentes aux structures\n\
    (courbes)\n\
  -niter %d : nombre d'iterations de SPARSE voting (en l'absence d'initialisation hessienne)\n\
  -cftv : champs de tenseurs de type Closed Form Tensor Voting\n\
  -tvclassic : champs de tenseurs de type Medioni\n\
  -nangles %d : nombre d'angles pour la discretisation des champs de \n\
    tenseurs\n\
  -nsticks %d : nombre d'angles pour calculer le champ de tenseurs plate\n\
# shape extraction (done after tensor voting, if specified)\n\
  -shape m|v|b|mv|mb|vb|mvb : extrait les maxima de {membranes, vaisseaux, balles}\n\
# general options\n\
  -write-images|-wi\n\
  -2D: slice per slice computation\n\
# parallelism parameters (specific)\n\
  -images-parallel-voting %d:\n\
    (3D) tensor voting is parallelized by tensor component (xx, yy, zz, xy, xz, yz)\n\
    meaning there is only 6 threads. This guarantees the same behavior than a\n\
    sequential version. This is the default behavior.\n\
    Additional parallelism can be done for the computation\n\
    of each component, requiring the allocation of auxiliary images.\n\
    '-images-parallel-voting n' means that 'n' additional images will be allocated\n\
    per component (a total of 6*n images), and that 'n+1' threads will be used for\n\
    each component.\n\
# parallelism parameters (general)\n\
 -parallel|-no-parallel:\n\
 -max-chunks %d:\n\
 -parallelism-type|-parallel-type default|none|openmp|omp|pthread|thread:\n\
 -omp-scheduling|-omps default|static|dynamic-one|dynamic|guided:\n\
# general parameters \n\
  -verbose|-v: increase verboseness\n\
    parameters being read several time, use '-nv -v -v ...'\n\
    to set the verboseness level\n\
  -no-verbose|-noverbose|-nv: no verboseness at all\n\
  -debug|-D: increase debug level\n\
  -no-debug|-nodebug: no debug indication\n\
  -print-parameters|-param:\n\
  -print-time|-time:\n\
  -no-time|-notime:\n\
  -trace-memory|-memory:\n\
  -no-memory|-nomemory:\n\
  -h: print option list\n\
  -help: print option list + details\n";






char *API_Help_TVmembrane( int h )
{
    if ( h == 0 )
        return( usage );
    return( detail );
}





void API_ErrorParse_TVmembrane( char *program, char *str, int flag )
{
    if ( flag >= 0 ) {
        if ( program != (char*)NULL )
           (void)fprintf(stderr,"Usage: %s %s\n", program, usage);
        else
            (void)fprintf(stderr,"Command line options: %s\n", usage);
    }
    if ( flag == 1 ) {
      (void)fprintf( stderr, "--------------------------------------------------\n" );
      (void)fprintf(stderr,"%s",detail);
      (void)fprintf( stderr, "--------------------------------------------------\n" );
    }
    if ( str != (char*)NULL )
      (void)fprintf(stderr,"Error: %s\n",str);
    exit( 1 );
}





/************************************************************
 *
 * parameters management
 *
 ************************************************************/



void API_InitParam_TVmembrane( lineCmdParamTVmembrane *par )
{
  /* names
   */

  par->inputBinary[0] = '\0';
  par->inputMask[0] = '\0';

  par->outputGeneric[0] = '\0';

  par->outputSampling[0] = '\0';

  par->outputAngles[0] = '\0';
  par->outputEigenvalues[0] = '\0';
  par->outputTensor[0] = '\0';
  par->outputBinary[0] = '\0';

  par->outputExtrema[0] = '\0';

  par->outputSuffix[0] = '\0';

  par->sample = 1.0;
  par->power = 0;
  par->sampleMode = _RANDOM_SAMPLING_;
  par->itermax = 0;

  par->scale = 1;
  par->flagReal = 0;
  par->zfact = 1;
  par->initHessian = NONE;
  par->niter = 1;

  par->TVmode = TVCLASSIC;
  par->nangles = 4;
  par->niter = 1;

  par->shape = NO_EXTRACTION;

  par->writeImages = 0;
  par->dimension = 3;

  par->print_lineCmdParam = 0;
  par->print_time = 0;
  par->trace_allocations = 0;
}





void API_PrintParam_TVmembrane( FILE *theFile, char *program,
                                  lineCmdParamTVmembrane *par, char *str )
{
  FILE *f = theFile;
  if ( theFile == (FILE*)NULL ) f = stderr;

  fprintf( f, "==================================================\n" );
  fprintf( f, "= in line command parameters" );
  if ( program != (char*)NULL )
    fprintf( f, " for '%s'", program );
  if ( str != (char*)NULL )
    fprintf( f, "= %s\n", str );
  fprintf( f, "\n"  );
  fprintf( f, "==================================================\n" );


  fprintf( f, "# image names\n" );

  fprintf( f, "- inputBinary =        " );
  if ( par->inputBinary != (char*)NULL && par->inputBinary[0] != '\0' )
    fprintf( f, "'%s'\n", par->inputBinary );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- inputMask =          " );
  if ( par->inputMask != (char*)NULL && par->inputMask[0] != '\0' )
    fprintf( f, "'%s'\n", par->inputMask );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- outputGeneric =      " );
  if ( par->outputGeneric != (char*)NULL && par->outputGeneric[0] != '\0' )
    fprintf( f, "'%s'\n", par->outputGeneric );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- outputSampling =     " );
  if ( par->outputSampling != (char*)NULL && par->outputSampling[0] != '\0' )
    fprintf( f, "'%s'\n", par->outputSampling );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- outputAngles =       " );
  if ( par->outputAngles != (char*)NULL && par->outputAngles[0] != '\0' )
    fprintf( f, "'%s'\n", par->outputAngles );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- outputEigenvalues =  " );
  if ( par->outputEigenvalues != (char*)NULL && par->outputEigenvalues[0] != '\0' )
    fprintf( f, "'%s'\n", par->outputEigenvalues );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- outputTensor =      " );
  if ( par->outputTensor != (char*)NULL && par->outputTensor[0] != '\0' )
    fprintf( f, "'%s'\n", par->outputTensor );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- outputBinary =      " );
  if ( par->outputBinary != (char*)NULL && par->outputBinary[0] != '\0' )
    fprintf( f, "'%s'\n", par->outputBinary );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- outputExtrema =      " );
  if ( par->outputExtrema != (char*)NULL && par->outputExtrema[0] != '\0' )
    fprintf( f, "'%s'\n", par->outputExtrema );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- outputSuffix =      " );
  if ( par->outputSuffix != (char*)NULL && par->outputSuffix[0] != '\0' )
    fprintf( f, "'%s'\n", par->outputSuffix );
  else
    fprintf( f, "'NULL'\n" );


  fprintf( f, "==================================================\n" );

  fprintf( f, "# Sampling  parameters\n" );

  fprintf( f, "- sample = %f\n", par->sample );
  fprintf( f, "- power = %d\n", par->power );
  fprintf( f, "- sampleMode = " );
  switch( par->sampleMode ) {
  default : fprintf( f, "unknown mode\n" ); break;
  case _UNKNOWN_ :          fprintf( f, "_UNKNOWN_\n" ); break;
  case _RANDOM_SAMPLING_ :  fprintf( f, "_RANDOM_SAMPLING_\n" ); break;
  case _REGULAR_SAMPLING_ : fprintf( f, "_REGULAR_SAMPLING_\n" ); break;
  }
  fprintf( f, "- itermax = %d\n", par->itermax );

  fprintf( f, "# Tensor voting parameters\n" );

  fprintf( f, "- scale = %f\n", par->scale );
  fprintf( f, "- flagReal = %d\n", par->flagReal );
  fprintf( f, "- zfact = %f\n", par->zfact );
  fprintf( f, "- initHessian = " );
  switch( par->initHessian ) {
  default : fprintf( f, "unknown mode\n" ); break;
  case NONE :  fprintf( f, "NONE\n" ); break;
  case PLANE : fprintf( f, "PLANE\n" ); break;
  case LINE :  fprintf( f, "LINE\n" ); break;
  }
  fprintf( f, "- niter = %d\n", par->niter );

  fprintf( f, "- TVmode = " );
  switch( par->TVmode ) {
  default : fprintf( f, "unknown mode\n" ); break;
  case TVCLASSIC : fprintf( f, "TVCLASSIC\n" ); break;
  case CFTV :      fprintf( f, "CFTV\n" ); break;
  }
  fprintf( f, "- nangles = %d\n", par->nangles );
  fprintf( f, "- nsticks = %d\n", par->nsticks );

  fprintf( f, "# Shape extraction parameters\n" );
  fprintf( f, "- shape = " );
  switch( par->shape ) {
  default : fprintf( f, "unknown mode\n" ); break;
  case NO_EXTRACTION : fprintf( f, "NO_EXTRACTION\n" ); break;
  case MEMBRANE :      fprintf( f, "MEMBRANE\n" ); break;
  case VESSEL :        fprintf( f, "VESSEL\n" ); break;
  case BALL : fprintf( f, "BALL\n" ); break;
  case MV :   fprintf( f, "MV\n" ); break;
  case MB :   fprintf( f, "MV\n" ); break;
  case VB :   fprintf( f, "VB\n" ); break;
  case MVB :  fprintf( f, "MVB\n" ); break;
  }

  fprintf( f, "# General parameters\n" );
  fprintf( f, "- writeImages = %d\n", par->writeImages );
  fprintf( f, "- dimension = %d\n", par->dimension );

  fprintf( f, "==================================================\n" );

  fprintf( f, "# general image related parameters\n" );

  fprintf( f, "# general parameters\n" );
  fprintf( f, "- print parameters     = %d\n", par->print_lineCmdParam );
  fprintf( f, "- print time           = %d\n", par->print_time );
  fprintf( f, "- par->trace_allocations = %d\n", par->trace_allocations );

  fprintf( f, "==================================================\n" );
}





/************************************************************
 *
 * parameters parsing
 *
 ************************************************************/



static void _API_ParseParam_TVmembrane( char *str, lineCmdParamTVmembrane *par )
{
  char *proc = "_API_ParseParam_TVmembrane";
  char **argv;
  int i, argc;

  if ( str == (char*)NULL || strlen(str) == 0 )
      return;

  argv = _Str2Array( &argc, str );
  if ( argv == (char**)NULL || argc == 0 ) {
      if ( _debug_ ) {
          fprintf( stderr, "%s: weird, no arguments were found\n", proc );
      }
      return;
  }

  if ( _debug_ > 4 ) {
      fprintf( stderr, "%s: translation from\n", proc );
      fprintf( stderr, "   '%s'\n", str );
      fprintf( stderr, "into\n" );
      for ( i=0; i<argc; i++ )
          fprintf( stderr, "   argv[%2d] = '%s'\n", i, argv[i] );
  }

  API_ParseParam_TVmembrane( 0, argc, argv, par );

  vtfree( argv );
}





static int _n_call_parse_ = 0;

void API_ParseParam_TVmembrane( int firstargc, int argc, char *argv[],
                                  lineCmdParamTVmembrane *par )
{
  int i=firstargc;
  int inputisread = 0;
  int outputisread = 0;
  char text[STRINGLENGTH];
  int status;
  int maxchunks;
  int additionalImages;
  /*int o=0, s=0, r=0;*/

  _n_call_parse_ ++;

  /* option line parsing
   */

  while ( i < argc ) {
    if ( argv[i][0] == '-' ) {

      if ( argv[i][1] == '\0' ) {
          /*--- standard input ---*/
          if ( inputisread == 0 ) {
            strcpy( par->inputBinary, "<" );
            inputisread = 1;
          }
          else {
            API_ErrorParse_TVmembrane( (char*)NULL, "too many file names, parsing '-' ...\n", 0 );
          }
      }

      /*--- names
       */
      else if ( strcmp ( argv[i], "-output-generic-prefix" ) == 0
                || strcmp ( argv[i], "-output-generic" ) == 0 ) {
        i += 1;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -output-generic-prefix...\n", 0 );
        if ( outputisread != 0 )
          API_ErrorParse_TVmembrane( (char*)NULL, "parsing -output-generic-prefix: output already read ...\n", 0 );
         strncpy( par->outputGeneric, argv[i], STRINGLENGTH );
         outputisread = 1;
      }
      else if ( strcmp ( argv[i], "-output-sampling" ) == 0 ) {
        i += 1;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -output-sampling ...\n", 0 );
        strncpy( par->outputSampling, argv[i], STRINGLENGTH );
      }
      else if ( strcmp ( argv[i], "-output-tensor" ) == 0 ) {
        i += 1;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -output-tensor ...\n", 0 );
        strncpy( par->outputTensor, argv[i], STRINGLENGTH );
      }
      else if ( strcmp ( argv[i], "-output-eigenvalues" ) == 0 ) {
        i += 1;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -output-eigenvalues ...\n", 0 );
        strncpy( par->outputEigenvalues, argv[i], STRINGLENGTH );
      }
      else if ( strcmp ( argv[i], "-output-angles" ) == 0 ) {
        i += 1;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -output-angles ...\n", 0 );
        strncpy( par->outputAngles, argv[i], STRINGLENGTH );
      }
      else if ( strcmp ( argv[i], "-output-binary" ) == 0 ) {
        i += 1;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -output-binary ...\n", 0 );
        strncpy( par->outputBinary, argv[i], STRINGLENGTH );
      }
      else if ( strcmp ( argv[i], "-output-extrema" ) == 0 ) {
        i += 1;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -output-extrema ...\n", 0 );
        strncpy( par->outputExtrema, argv[i], STRINGLENGTH );
      }
      else if ( strcmp ( argv[i], "-output-suffix" ) == 0 ) {
        i += 1;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -output-suffix ...\n", 0 );
        strncpy( par->outputSuffix, argv[i], STRINGLENGTH );
      }


      else if ( strcmp ( argv[i], "-mask" ) == 0 ) {
        i += 1;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -mask...\n", 0 );
        strncpy( par->inputMask, argv[i], STRINGLENGTH );
      }

      /*--- sampling options
       */
      else if ( strcmp ( argv[i], "-sample" ) == 0 ) {
        i += 1;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -sample...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->sample) );
        if ( status <= 0 ) API_ErrorParse_TVmembrane( (char*)NULL, "parsing -sample...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-power" ) == 0 && argv[i][6] == '\0'  ) {
        par->power = 1;
      }

      else if ( strcmp ( argv[i], "-random" ) == 0 ) {
        par->sampleMode = _RANDOM_SAMPLING_;
      }

      else if ( strcmp ( argv[i], "-regular" ) == 0 ) {
        par->sampleMode = _REGULAR_SAMPLING_;
      }

      else if ( strcmp ( argv[i], "-sampling-iterations" ) == 0
                || strcmp ( argv[i], "-iterations" ) == 0
                || (strcmp ( argv[i], "-si" ) == 0 && argv[i][3] == '\0')
                || (strcmp ( argv[i], "-i" ) == 0 && argv[i][2] == '\0') ) {
        i += 1;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -sampling-iterations...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->itermax) );
        if ( status <= 0 ) API_ErrorParse_TVmembrane( (char*)NULL, "parsing -sampling-iterations...\n", 0 );
        parcelling_setNumberOfIterations( par->itermax );
      }

      /*--- tensor voting options
       */
      else if ( strcmp ( argv[i], "-scale" ) == 0 ) {
        i += 1;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -scale...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->scale) );
        if ( status <= 0 ) API_ErrorParse_TVmembrane( (char*)NULL, "parsing -scale...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-real" ) == 0 ) {
        par->flagReal = 1;
      }
      else if ( strcmp ( argv[i], "-zfact" ) == 0 ) {
        i += 1;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -zfact...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->zfact) );
        if ( status <= 0 ) API_ErrorParse_TVmembrane( (char*)NULL, "parsing -zfact...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-hessian" ) == 0 ) {
        par->initHessian = NONE;
        i += 1;
        if ( i < argc) {
          if ( strcmp(argv[i],"plane") == 0 )
                        par->initHessian = PLANE;
          else if ( strcmp(argv[i],"line") == 0 )
                        par->initHessian = LINE;
        }
        if (par->initHessian == NONE) {
              par->initHessian = PLANE;
              i -= 1;
        }
      }

      else if ( strcmp ( argv[i], "-niter" ) == 0 ) {
        i += 1;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -niter...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->niter) );
        if ( status <= 0 ) API_ErrorParse_TVmembrane( (char*)NULL, "parsing -niter...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-cftv" ) == 0 && argv[i][5] == '\0' ) {
        par->TVmode = CFTV;
      }
      else if ( strcmp ( argv[i], "-tvclassic" ) == 0 && argv[i][10] == '\0' ) {
        par->TVmode = TVCLASSIC;
      }

      else if ( strcmp ( argv[i], "-nangles" ) == 0 ) {
        i += 1;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -nangles...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->nangles) );
        if ( status <= 0 ) API_ErrorParse_TVmembrane( (char*)NULL, "parsing -nangles...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-nsticks" ) == 0 ) {
        i += 1;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -nsticks...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->nsticks) );
        if ( status <= 0 ) API_ErrorParse_TVmembrane( (char*)NULL, "parsing -nsticks...\n", 0 );
      }

      /* shape extraction
       */
      else if ( strcmp ( argv[i], "-shape" ) == 0 && argv[i][6] == '\0' ) {
        i += 1;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -shape...\n", 0 );
        if ( strcmp(argv[i],"m") == 0 && argv[i][1] == '\0' )
                  par->shape = MEMBRANE;
        else if ( strcmp(argv[i],"v") == 0 && argv[i][1] == '\0' )
                  par->shape = VESSEL;
        else if ( strcmp(argv[i],"b") == 0 && argv[i][1] == '\0' )
                  par->shape = BALL;
        else if ( strcmp(argv[i],"mv") == 0 && argv[i][2] == '\0' )
                  par->shape = MV;
        else if ( strcmp(argv[i],"mb") == 0 && argv[i][2] == '\0' )
                  par->shape = MB;
        else if ( strcmp(argv[i],"vb") == 0 && argv[i][2] == '\0' )
                  par->shape = VB;
        else if ( strcmp(argv[i],"mvb") == 0 && argv[i][3] == '\0' )
                  par->shape = MVB;
        else
          API_ErrorParse_TVmembrane( (char*)NULL, "parsing -shape...\n", 0 );

      }

      /*--- general parameters
       */
      else if ( strcmp ( argv[i], "-write-images" ) == 0
                || (strcmp ( argv[i], "-wi" ) == 0 && argv[i][3] == '\0') ) {
        par->writeImages = 1;
      }

      else if ( strcmp ( argv[i], "-2D" ) == 0 && argv[i][3] == '\0' ) {
        par->dimension = 2;
      }

      /* parallelism (specific)
       */
      else if ( strcmp ( argv[i], "-images-parallel-voting" ) == 0 ) {
          i ++;
          if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -images-parallel-voting", 0 );
          status = sscanf( argv[i], "%d", &additionalImages );
          if ( status <= 0 ) API_ErrorParse_TVmembrane( (char*)NULL, "parsing -images-parallel-voting", 0 );
          if ( additionalImages >= 0 )
              MT_SetAdditionalImagesForParallelVoting( additionalImages );
      }

      /* parallelism (general)
       */
      else if ( strcmp ( argv[i], "-parallel" ) == 0 && argv[i][9] == '\0' ) {
        setParallelism( _DEFAULT_PARALLELISM_ );
      }

      else if ( strcmp ( argv[i], "-no-parallel" ) == 0 ) {
         setParallelism( _NO_PARALLELISM_ );
      }

      else if ( strcmp ( argv[i], "-parallelism-type" ) == 0 ||
                 strcmp ( argv[i], "-parallel-type" ) == 0 ) {
        i ++;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -parallelism-type", 0 );
        if ( strcmp ( argv[i], "default" ) == 0 ) {
          setParallelism( _DEFAULT_PARALLELISM_ );
        }
        else if ( strcmp ( argv[i], "none" ) == 0 ) {
          setParallelism( _NO_PARALLELISM_ );
        }
        else if ( strcmp ( argv[i], "openmp" ) == 0 || strcmp ( argv[i], "omp" ) == 0 ) {
          setParallelism( _OMP_PARALLELISM_ );
        }
        else if ( strcmp ( argv[i], "pthread" ) == 0 || strcmp ( argv[i], "thread" ) == 0 ) {
          setParallelism( _PTHREAD_PARALLELISM_ );
        }
      }


      else if ( strcmp ( argv[i], "-max-chunks" ) == 0 ) {
        i ++;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -max-chunks", 0 );
        status = sscanf( argv[i], "%d", &maxchunks );
        if ( status <= 0 ) API_ErrorParse_TVmembrane( (char*)NULL, "parsing -max-chunks", 0 );
        if ( maxchunks >= 1 ) setMaxChunks( maxchunks );
      }

      else if ( strcmp ( argv[i], "-omp-scheduling" ) == 0 ||
                ( strcmp ( argv[i], "-omps" ) == 0 && argv[i][3] == '\0') ) {
        i ++;
        if ( i >= argc)    API_ErrorParse_TVmembrane( (char*)NULL, "parsing -parallel-scheduling", 0 );
        if ( strcmp ( argv[i], "default" ) == 0 ) {
          setOmpScheduling( _DEFAULT_OMP_SCHEDULING_ );
        }
        else if ( strcmp ( argv[i], "static" ) == 0 ) {
          setOmpScheduling( _STATIC_OMP_SCHEDULING_ );
        }
        else if ( strcmp ( argv[i], "dynamic-one" ) == 0 ) {
          setOmpScheduling( _DYNAMIC_ONE_OMP_SCHEDULING_ );
        }
        else if ( strcmp ( argv[i], "dynamic" ) == 0 ) {
          setOmpScheduling( _DYNAMIC_OMP_SCHEDULING_ );
        }
        else if ( strcmp ( argv[i], "guided" ) == 0 ) {
          setOmpScheduling( _GUIDED_OMP_SCHEDULING_ );
        }
        else {
          fprintf( stderr, "unknown omp scheduling type: '%s'\n", argv[i] );
          API_ErrorParse_TVmembrane( (char*)NULL, "parsing -omp-scheduling ...\n", 0 );
        }
      }


      /* --- general parameters
       */

      else if ( strcmp ( argv[i], "-verbose" ) == 0
                || (strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0') ) {
        if ( _n_call_parse_ == 1 ) {
          if ( _verbose_ <= 0 ) _verbose_ = 1;
          else                  _verbose_ ++;
          if ( _VT_VERBOSE_ <= 0 ) _VT_VERBOSE_ = 1;
          else                     _VT_VERBOSE_ ++;
          MT_IncrementVerboseInMtMembrane2D();
          MT_IncrementVerboseInMtMembrane3D();
        }
      }
      else if ( strcmp ( argv[i], "-no-verbose" ) == 0
                || strcmp ( argv[i], "-noverbose" ) == 0
                || (strcmp ( argv[i], "-nv" ) == 0 && argv[i][3] == '\0') ) {
          _verbose_ = 0;
          _VT_VERBOSE_ = 0;
          MT_SetVerboseInMtMembrane2D( 0 );
          MT_SetVerboseInMtMembrane3D( 0 );
      }

      else if ( (strcmp ( argv[i], "-debug" ) == 0 && argv[i][6] == '\0')
                || (strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0') ) {
        if ( _n_call_parse_ == 1 ) {
          if ( _debug_ <= 0 ) _debug_ = 1;
          else                _debug_ ++;
          if ( _VT_DEBUG_ <= 0 ) _VT_DEBUG_ = 1;
          else                   _VT_DEBUG_ ++;
        }
      }
      else if ( (strcmp ( argv[i], "-no-debug" ) == 0 && argv[i][9] == '\0')
                || (strcmp ( argv[i], "-nodebug" ) == 0 && argv[i][8] == '\0') ) {
          _debug_ = 0;
          _VT_DEBUG_ = 0;
      }

      else if ( strcmp ( argv[i], "-print-parameters" ) == 0
                || (strcmp ( argv[i], "-param" ) == 0 && argv[i][6] == '\0') ) {
         par->print_lineCmdParam = 1;
      }

      else if ( strcmp ( argv[i], "-print-time" ) == 0
                 || (strcmp ( argv[i], "-time" ) == 0 && argv[i][5] == '\0') ) {
         par->print_time = 1;
      }
      else if ( (strcmp ( argv[i], "-notime" ) == 0 && argv[i][7] == '\0')
                  || (strcmp ( argv[i], "-no-time" ) == 0 && argv[i][8] == '\0') ) {
         par->print_time = 0;
      }

      else if ( strcmp ( argv[i], "-trace-memory" ) == 0
                 || (strcmp ( argv[i], "-memory" ) == 0 && argv[i][7] == '\0') ) {
         if ( _n_call_parse_ == 1 ) {
           incrementTraceInVtMalloc( );
           if ( par->trace_allocations  <= 0 ) par->trace_allocations  = 1;
           else                              par->trace_allocations  ++;
         }
      }
      else if ( (strcmp ( argv[i], "-nomemory" ) == 0 && argv[i][9] == '\0')
                  || (strcmp ( argv[i], "-no-memory" ) == 0 && argv[i][10] == '\0') ) {
         setTraceInVtMalloc( 0 );
      }

      else if ( (strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0')
                || (strcmp ( argv[i], "--help" ) == 0 && argv[i][6] == '\0') ) {
         API_ErrorParse_TVmembrane( (char*)NULL, (char*)NULL, 1);
      }
      else if ( (strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0')
                || (strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0') ) {
         API_ErrorParse_TVmembrane( (char*)NULL, (char*)NULL, 0);
      }

      /*--- option inconnue ---*/
      else {
        sprintf(text,"unknown option %s\n",argv[i]);
        API_ErrorParse_TVmembrane( (char*)NULL, text, 0 );

      }
    }
    /*--- saisie des noms d'images ---*/
    else if ( argv[i][0] != 0 ) {
        if ( strlen( argv[i] ) >= STRINGLENGTH ) {
            fprintf( stderr, "... parsing '%s'\n", argv[i] );
            API_ErrorParse_TVmembrane( (char*)NULL, "too long file name ...\n", 0 );
        }
        else if ( inputisread == 0 ) {
            (void)strcpy( par->inputBinary, argv[i] );
            inputisread = 1;
        }
        else if ( outputisread == 0 ) {
            (void)strcpy( par->outputGeneric, argv[i] );
            outputisread = 1;
        }
        else {
            fprintf( stderr, "... parsing '%s'\n", argv[i] );
            API_ErrorParse_TVmembrane( (char*)NULL, "too many file names ...\n", 0 );
        }
    }
    i += 1;
  }

}
