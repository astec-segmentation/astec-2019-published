/*************************************************************************
 * api-membrane.c -
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <chunks.h>
#include <vtmalloc.h>

#include <vt_common.h>

#include <api-membrane.h>






static int _verbose_ = 1;
static int _debug_ = 0;


static void _API_ParseParam_membrane( char *str, lineCmdParamMembrane *p );



/************************************************************
 *
 * main API
 *
 ************************************************************/



int API_membrane( vt_image *image, vt_3Dimres *imsRes, vt_image *mask, char *param_str_1, char *param_str_2 )
{
  char *proc = "API_membrane";
  lineCmdParamMembrane par;


  int flag_3D=1;


  /* parameter initialization
   */
  API_InitParam_membrane( &par );

  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_membrane( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_membrane( param_str_2, &par );

  /*
   * if ( par.print_lineCmdParam )
      API_PrintParam_membrane( stderr, proc, &par, (char*)NULL );
  */

  /************************************************************
   *
   *  here is the stuff
   *
   ************************************************************/



  if ( par.dimension == 2 || image->dim.z == 1 )
    flag_3D = 0;

  switch ( par.typeComputation ) {

  case SINGLE_SCALE :

      par.scale2 = par.scale1;
      par.nbscales = 1;
      /* Falls through. */

  case MULTI_SCALE :

    if ( flag_3D ) {
      /*  3D case */
      if ( MT_Compute3DMultiScale( image, imsRes, mask, par.flagMask, par.scale1,
                           par.scale2, par.nbscales, par.zfact,
                           par.structureColor, par.mode, par.hsp ) != 1 ) {
          API_ErrorParse_membrane( proc, "some error occurs during processing ...\n", -1 );
      }
    }
    else {
      /*  2D case */
      if ( MT_Compute2DMultiScale( image, imsRes, mask, par.flagMask, par.scale1,
                                   par.scale2, par.nbscales,
                                   par.structureColor ) != 1 ) {
          API_ErrorParse_membrane( proc, "some error occurs during processing ...\n", -1 );
      }

    }

    break;

  default :
    {
      API_ErrorParse_membrane( proc, "incorrect typeComputation selected\n", -1 );
    }
  }

  return( 1 );
}

int API_extrema(vt_3Dimres *imsRes, vt_image *mask, int flagMask, char *param_str_1, char *param_str_2 )
{
    /*
    char *proc = "API_extrema";
     */
    lineCmdParamMembrane par;

    int flag_3D;

    /* parameter initialization
     */
    API_InitParam_membrane( &par );

    /* parameter parsing
     */
    if ( param_str_1 != (char*)NULL )
        _API_ParseParam_membrane( param_str_1, &par );
    if ( param_str_2 != (char*)NULL )
        _API_ParseParam_membrane( param_str_2, &par );

    /*
     * if ( par.print_lineCmdParam )
        API_PrintParam_membrane( stderr, proc, &par, (char*)NULL );
    */

    /************************************************************
     *
     *  here is the stuff
     *
     ************************************************************/


    if ( par.dimension == 2 || imsRes->imRep.dim.z == 1 )
      flag_3D = 0;

    if (flag_3D) {
        MT_Compute3DExtrema( imsRes, mask, flagMask, par.zfact, &(imsRes->imTheta) );
    }
    else {
        MT_Compute2DExtrema( imsRes, mask, flagMask, &(imsRes->imTheta) );
    }
    sprintf( imsRes->imTheta.name, "%s.ext.inr", par.names.out );

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

/*

static char *usage = "[image-in] [image-out]\n\
 [-parallel|-no-parallel] [-max-chunks %d]\n\
 [-parallelism-type|-parallel-type default|none|openmp|omp|pthread|thread]\n\
 [-omp-scheduling|-omps default|static|dynamic-one|dynamic|guided]\n\
 [-inv] [-swap] [output-image-type | -type s8|u8|s16|u16...]\n\
 [-verbose|-v] [-no-verbose|-noverbose|-nv]\n\
 [-debug|-D] [-no-debug|-nodebug]\n\
 [-allow-pipe|-pipe] [-no-allow-pipe|-no-pipe|-nopipe]\n\
 [-print-parameters|-param]\n\
 [-print-time|-time] [-no-time|-notime]\n\
 [-trace-memory|-memory] [-no-memory|-nomemory]\n\
 [-help|-h]";



static char *detail = "\
 if 'image-in' is equal to '-', stdin will be used\n\
 if 'image-out' is not specified or equal to '-', stdout will be used\n\
 if both are not specified, stdin and stdout will be used\n\
# ...\n\
# parallelism parameters\n\
 -parallel|-no-parallel:\n\
 -max-chunks %d:\n\
 -parallelism-type|-parallel-type default|none|openmp|omp|pthread|thread:\n\
 -omp-scheduling|-omps default|static|dynamic-one|dynamic|guided:\n\
# general image related parameters\n\
  -inv: inverse 'image-in'\n\
  -swap: swap 'image-in' (if encoded on 2 or 4 bytes)\n\
   output-image-type: -o 1    : unsigned char\n\
                      -o 2    : unsigned short int\n\
                      -o 2 -s : short int\n\
                      -o 4 -s : int\n\
                      -r      : float\n\
  -type s8|u8|s16|u16|... \n\
   default is type of input image\n\
# general parameters \n\
  -verbose|-v: increase verboseness\n\
    parameters being read several time, use '-nv -v -v ...'\n\
    to set the verboseness level\n\
  -no-verbose|-noverbose|-nv: no verboseness at all\n\
  -debug|-D: increase debug level\n\
  -no-debug|-nodebug: no debug indication\n\
  -allow-pipe|-pipe: allow the use of stdin/stdout (with '-')\n\
  -no-allow-pipe|-no-pipe|-nopipe: do not allow the use of stdin/stdout\n\
  -print-parameters|-param:\n\
  -print-time|-time:\n\
  -no-time|-notime:\n\
  -trace-memory|-memory:\n\
  -no-memory|-nomemory:\n\
  -h: print option list\n\
  -help: print option list + details\n\
";
*/


static char *usage = "[image-in] [prefix-images-out] [-mask image-mask]\n\
\t [-init %lf] [-last %lf] [-real] [-nb %d | -single] [-wi] [-2D]\n\
\t [-wi] [-black|-noir|-white|-blanc] [-modelbased|-acme]\n\
\t [-hsp %d]\n\
\t [-parallel|-no-parallel] [-max-chunks %d]\n\
\t [-parallelism-type|-parallel-type default|none|openmp|omp|pthread|thread]\n\
\t [-omp-scheduling|-omps default|static|dynamic-one|dynamic|guided]\n\
\t [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t 'prefix-images-out' : indiquer le prefixe des images de sortie (ie sans extension)\n\
\t -mask image-mask : image binaire (u8 ou u16) telle que le calcul de la fonction de reponse n'est realise que dans ce masque (voxels non-nuls)\n\
\t -init %lf : echelle (sigma) fine\n\
\t -last %lf : echelle (sigma) grossiere\n\
\t -real : afin d'exprimer les echelles en coordonnees reelles, ie en considerant la resolution de l'image\n\
\t (note : sans cette option, on travaille en echelle voxellique, ce qui implique de travailler sur des images isotropes)\n\
\t -nb   %d  : nombre d'echelles (interpolation logarithmique)\n\
\t -single : une seule echelle de calcul (specifiee dans -init)\n\
\t -wi : ecrit toutes les images intermediaires (extension '.inr')\n\
\t       .theta, .phi : angles definissant le vecteur 3D\n\
\t       .rep : reponse 3D multi-echelle\n\
\t       .scale : echelle de la reponse maximale\n\
\t       .ext : extrema (plans centraux) : SEULE IMAGE SAUVEE PAR DEFAUT\n\
\t [-modelbased|-acme] : selectionne le mode de calcul de la reponse (defaut : modelbased ~ [Krissian2000]-like)\n\
\t -hsp  %d  : half size plane (pour le cas modelbased : integration sur \n\
\t des petits plans plutot que seulement sur 2 points) (option inactive par defaut)\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";
/*\t -zfact %lf: rapport d'echelle entre l'axe z et les axes x et y\n\*/




char *API_Help_membrane( int h )
{
    if ( h == 0 )
        return( usage );
    return( detail );
}





void API_ErrorParse_membrane( char *program, char *str, int flag )
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



void API_InitParam_membrane( lineCmdParamMembrane *par )
{
    /*
    (void)strncpy( p->input_name, "\0", 1 );
    (void)strncpy( p->output_name, "\0", 1 );
    p->input_inv = 0;
    p->input_swap = 0;
    p->output_type = TYPE_UNKNOWN;

    p->allow_stdin_stdout = 1;
    */
    par->print_lineCmdParam = 0;
    par->print_time = 0;
    par->trace_allocations = 0;


    VT_Names( &(par->names) );
    par->type = TYPE_UNKNOWN;

    par->writeImages = 0;
    par->typeComputation = MULTI_SCALE;
    par->mode = MODELBASED;

    par->structureColor = _WHITE_;

    par->scale1 = 1.0;
    par->scale2 = 1.0;
    par->nbscales = 1;
    par->zfact = 1.0;

    par->dimension = 3;
    par->hsp = 0;

    par->flagMask = 0;

}



/*

void API_PrintParam_membrane( FILE *theFile, char *program,
                                  lineCmdParamMembrane *p, char *str )
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

  fprintf( f, "- input image is " );
  if ( p->input_name != (char*)NULL && p->input_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- output image is " );
  if ( p->output_name != (char*)NULL && p->output_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_name );
  else
    fprintf( f, "'NULL'\n" );


  fprintf( f, "# ...\n" );

  fprintf( f, "# general image related parameters\n" );

  fprintf( f, "- input image inverse = %d\n", p->input_inv );
  fprintf( f, "- input image swap = %d\n", p->input_swap );
  fprintf( f, "- output image type = " );
  switch ( p->output_type ) {
  default :     fprintf( f, "TYPE_UNKNOWN\n" ); break;
  case SCHAR :  fprintf( f, "SCHAR\n" ); break;
  case UCHAR :  fprintf( f, "UCHAR\n" ); break;
  case SSHORT : fprintf( f, "SSHORT\n" ); break;
  case USHORT : fprintf( f, "USHORT\n" ); break;
  case UINT :   fprintf( f, "UINT\n" ); break;
  case SINT :   fprintf( f, "INT\n" ); break;
  case ULINT :  fprintf( f, "ULINT\n" ); break;
  case FLOAT :  fprintf( f, "FLOAT\n" ); break;
  case DOUBLE : fprintf( f, "DOUBLE\n" ); break;
  }

  fprintf( f, "# general parameters\n" );
  fprintf( f, "- allows stdin/stdout  = %d\n", p->allow_stdin_stdout );
  fprintf( f, "- print parameters     = %d\n", p->print_lineCmdParam );
  fprintf( f, "- print time           = %d\n", p->print_time );
  fprintf( f, "- p->trace_allocations = %d\n", p->trace_allocations );

  fprintf( f, "==================================================\n" );
}

*/



/************************************************************
 *
 * parameters parsing
 *
 ************************************************************/



static void _API_ParseParam_membrane( char *str, lineCmdParamMembrane *p )
{
  char *proc = "_API_ParseParam_membrane";
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

  API_ParseParam_membrane( 0, argc, argv, p );

  vtfree( argv );
}




static char program[STRINGLENGTH];

static int _n_call_parse_ = 0;

void API_ParseParam_membrane( int firstargc, int argc, char *argv[],
                                  lineCmdParamMembrane *par )
{
    int i, nb, status;
    int maxchunks;
    char text[STRINGLENGTH];

    _n_call_parse_ ++;

    if ( VT_CopyName( program, argv[0] ) != 1 )
      VT_Error("Error while copying program name", (char*)NULL);
    if ( argc == 1 )    API_ErrorParse_membrane( (char*)NULL, (char*)NULL, 0);
        /*MT_ErrorParse("\n", 0 );*/

    /*--- lecture des parametres ---*/
    i = firstargc; nb = 0;
    while ( i < argc ) {
      if ( argv[i][0] == '-' ) {
        if ( argv[i][1] == '\0' ) {
          if ( nb == 0 ) {
            /*--- standart input ---*/
            strcpy( par->names.in, "<" );
            nb += 1;
          }
        }
        /*--- arguments generaux ---*/
        else if ( strcmp ( argv[i], "-help" ) == 0 ) {
          /*MT_ErrorParse("\n", 1);*/
          API_ErrorParse_membrane( (char*)NULL, (char*)NULL, 1);
        }
        else if ( strcmp ( argv[i], "-v" ) == 0 ) {
          _VT_VERBOSE_ = 1;
          VT_IncrementVerboseInVtTube3D(  );
        }
        else if ( strcmp ( argv[i], "-D" ) == 0 ) {
          _VT_DEBUG_ = 1;
          VT_IncrementDebugInVtTube3D(  );
        }

        /*--- traitement eventuel de l'image d'entree ---*/
        else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
          par->names.inv = 1;
        }
        else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
          par->names.swap = 1;
        }

        else if ( strcmp ( argv[i], "-mask" ) == 0 ) {
            i += 1;
            if ( i >= argc)    API_ErrorParse_membrane( (char*)NULL, "parsing -mask...\n", 0);
            strcpy( par->names.ext, argv[i] );
            par->flagMask=1;
        }




        else if ( strcmp ( argv[i], "-single" ) == 0 ) {
          par->typeComputation = SINGLE_SCALE;
        }

        else if ( strcmp ( argv[i], "-acme" ) == 0 ) {
          par->mode = ACME;
        }

        else if ( strcmp ( argv[i], "-modelbased" ) == 0 ) {
          par->mode = MODELBASED;
        }


        else if ( strcmp ( argv[i], "-wi" ) == 0 ) {
          par->writeImages = 1;
        }



        else if ( strcmp ( argv[i], "-black" ) == 0
                  || strcmp ( argv[i], "-noir" ) == 0 ) {
          par->structureColor = _BLACK_;
        }
        else if ( strcmp ( argv[i], "-white" ) == 0
                  || strcmp ( argv[i], "-blanc" ) == 0 ) {
          par->structureColor = _WHITE_;
        }



        else if ( strcmp ( argv[i], "-2D" ) == 0 ) {
          par->dimension = 2;
        }



        /* Parametres de calcul */
        else if ( strcmp ( argv[i], "-init" ) == 0 ) {
          i += 1;
          if ( i >= argc)    API_ErrorParse_membrane( (char*)NULL, "parsing -init...\n", 0);
          status = sscanf( argv[i],"%lf",&(par->scale1) );
          if ( status <= 0 ) API_ErrorParse_membrane( (char*)NULL, "parsing -init...\n", 0);
        }
        else if ( strcmp ( argv[i], "-last" ) == 0 ) {
          i += 1;
          if ( i >= argc)    API_ErrorParse_membrane( (char*)NULL, "parsing -last...\n", 0);
          status = sscanf( argv[i],"%lf",&(par->scale2) );
          if ( status <= 0 ) API_ErrorParse_membrane( (char*)NULL, "parsing -last...\n", 0);
        }
        else if ( strcmp ( argv[i], "-nb" ) == 0 ) {
          i += 1;
          if ( i >= argc)    API_ErrorParse_membrane( (char*)NULL, "parsing -nb...\n", 0);
          status = sscanf( argv[i],"%d",&(par->nbscales) );
          if ( status <= 0 ) API_ErrorParse_membrane( (char*)NULL, "parsing -nb...\n", 0);
        }
        else if ( strcmp ( argv[i], "-zfact" ) == 0 ) {
          i += 1;
          if ( i >= argc)    API_ErrorParse_membrane( (char*)NULL, "parsing -zfact...\n", 0);
          status = sscanf( argv[i],"%lf",&(par->zfact) );
          if ( status <= 0 ) API_ErrorParse_membrane( (char*)NULL, "parsing -zfact...\n", 0);
        }
        else if ( strcmp ( argv[i], "-hsp" ) == 0 ) {
          i += 1;
          if ( i >= argc)    API_ErrorParse_membrane( (char*)NULL, "parsing -hsp...\n", 0);
          status = sscanf( argv[i],"%d",&(par->hsp) );
          if ( status <= 0 ) API_ErrorParse_membrane( (char*)NULL, "parsing -hsp...\n", 0);
        }
        else if ( strcmp ( argv[i], "-real" ) == 0 ) {
          par->zfact = -1;
        }

        /* parallelism
         */
      else if ( strcmp ( argv[i], "-parallel" ) == 0 ) {
        setParallelism( _DEFAULT_PARALLELISM_ );
      }

      else if ( strcmp ( argv[i], "-no-parallel" ) == 0 ) {
         setParallelism( _NO_PARALLELISM_ );
      }

      else if ( strcmp ( argv[i], "-parallelism-type" ) == 0 ||
                 strcmp ( argv[i], "-parallel-type" ) == 0 ) {
        i ++;
        if ( i >= argc)    API_ErrorParse_membrane( (char*)NULL, "parsing -parallelism-type...\n", 0);
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
        if ( i >= argc)    API_ErrorParse_membrane( (char*)NULL, "parsing -max-chunks...\n", 0);
        status = sscanf( argv[i], "%d", &maxchunks );
        if ( status <= 0 ) API_ErrorParse_membrane( (char*)NULL, "parsing -max-chunks...\n", 0);
        if ( maxchunks >= 1 ) setMaxChunks( maxchunks );
      }

      else if ( strcmp ( argv[i], "-omp-scheduling" ) == 0 ||
                ( strcmp ( argv[i], "-omps" ) == 0 && argv[i][3] == '\0') ) {
        i ++;
        if ( i >= argc)    API_ErrorParse_membrane( (char*)NULL, "parsing -omp-scheduling...\n", 0);
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
          API_ErrorParse_membrane( (char*)NULL, "parsing -omp-scheduling...\n", 0);
        }
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
           if ( 0 ) setParallelism( _NO_PARALLELISM_ );
        }
        else if ( (strcmp ( argv[i], "-nomemory" ) == 0 && argv[i][9] == '\0')
                    || (strcmp ( argv[i], "-no-memory" ) == 0 && argv[i][10] == '\0') ) {
           setTraceInVtMalloc( 0 );
        }


        /*--- option inconnue ---*/
        else {
          sprintf(text,"unknown option %s\n",argv[i]);
          API_ErrorParse_membrane( (char*)NULL, text, 0);
        }
      }
      /*--- saisie des noms d'images ---*/
      else if ( argv[i][0] != 0 ) {
        if ( nb == 0 ) {
          strncpy( par->names.in, argv[i], STRINGLENGTH );
          nb += 1;
        }
        else if ( nb == 1 ) {
          strncpy( par->names.out, argv[i], STRINGLENGTH );
          nb += 1;
        }
        else
          API_ErrorParse_membrane( (char*)NULL, "too much file names when parsing\n", 0);
      }
      i += 1;
    }

    /*--- s'il n'y a pas assez de noms ... ---*/
    if (nb == 0) {
      strcpy( par->names.in,  "<" );  /* standart input */
      strcpy( par->names.out, ">" );  /* standart output */
    }
    if (nb == 1)
      strcpy( par->names.out, ">" );  /* standart output */

}
