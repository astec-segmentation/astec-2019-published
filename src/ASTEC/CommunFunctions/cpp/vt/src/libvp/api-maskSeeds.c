/*************************************************************************
 * api-maskSeeds.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mer 25 jul 2018 18:06:10 CEST
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <api-maskSeeds.h>

#include <vt_maskSeeds.h>




static int _verbose_ = 1;
static int _debug_ = 0;







static void _API_ParseParam_maskSeeds( char *str, lineCmdParamMaskSeeds *p );



/************************************************************
 *
 * main API
 *
 ************************************************************/



int API_maskSeeds( vt_image *imseed, vt_image *imcell, vt_image *imres,
                    char *param_str_1, char *param_str_2 )
{
  char *proc = "API_maskSeeds";
  lineCmdParamMaskSeeds par;





  /* parameter initialization
   */
  API_InitParam_maskSeeds( &par );

  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_maskSeeds( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_maskSeeds( param_str_2, &par );

  if ( par.print_lineCmdParam )
      API_PrintParam_maskSeeds( stderr, proc, &par, (char*)NULL );

  /************************************************************
   *
   *  here is the stuff
   *
   ************************************************************/
  if ( VT_SelectSeedsinCells( imseed, imcell, imres ) != 1 ) {
      if ( _verbose_ ) {
        fprintf( stderr, "%s: error when processing\n", proc );
      }
      return( -1 );
  }



  return( 1 );
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
  array = (char**)malloc( n * sizeof(char*) + (strlen(str)+1) * sizeof(char) );
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



static char *usage = "-seed-image|-si %s -cell-image|-ci %s \n\
 [-output-image|-oi|-o] %s\n\
 [-verbose|-v] [-nv|-noverbose] [-debug|-D] [-nodebug]\n\
 [-print-parameters|-param]\n\
 [-print-time|-time] [-notime]\n\
 [-help|-h]";



static char *detail = "\
 select the 'seeds' (ie labeled connected components) from the seed image\n\
 that are totally included in single cell (intersect only one label)\n\
 from the cell image. Other seeds are removed.\n\
 -seed-image|-si %s\n\
 -cell-image|-ci %s\n\
 [-output-image|-oi|-o] %s\n\
# general parameters \n\
  -verbose|-v: increase verboseness\n\
    parameters being read several time, use '-nv -v -v ...'\n\
    to set the verboseness level\n\
  -noverbose|-nv: no verboseness at all\n\
  -debug|-D: increase debug level\n\
  -nodebug: no debug indication\n\
  -print-parameters|-param:\n\
  -print-time|-time:\n\
  -no-time|-notime:\n\
  -h: print option list\n\
  -help: print option list + details\n\
";





char *API_Help_maskSeeds( int h )
{
  if ( h == 0 )
    return( usage );
  return( detail );
}





void API_ErrorParse_maskSeeds( char *program, char *str, int flag )
{
  if ( program != (char*)NULL )
     (void)fprintf(stderr,"Usage: %s %s\n", program, usage);
  else
      (void)fprintf(stderr,"Command line options: %s\n", usage);
  if ( flag == 1 ) {
    (void)fprintf( stderr, "--------------------------------------------------\n" );
    (void)fprintf(stderr,"%s",detail);
    (void)fprintf( stderr, "--------------------------------------------------\n" );
  }
  if ( str != (char*)NULL )
    (void)fprintf(stderr,"Error: %s",str);
  exit( 1 );
}





/************************************************************
 *
 * parameters management
 *
 ************************************************************/



void API_InitParam_maskSeeds( lineCmdParamMaskSeeds *p )
{
    (void)strncpy( p->input_seed_name, "\0", 1 );
    (void)strncpy( p->input_cell_name, "\0", 1 );
    (void)strncpy( p->output_name, "\0", 1 );

    p->print_lineCmdParam = 0;
    p->print_time = 0;
}





void API_PrintParam_maskSeeds( FILE *theFile, char *program,
                                  lineCmdParamMaskSeeds *p, char *str )
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

  fprintf( f, "- input seed image is " );
  if ( p->input_seed_name != (char*)NULL && p->input_seed_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_seed_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- input cell image is " );
  if ( p->input_cell_name != (char*)NULL && p->input_cell_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_cell_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- output image is " );
  if ( p->output_name != (char*)NULL && p->output_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "==================================================\n" );
}





/************************************************************
 *
 * parameters parsing
 *
 ************************************************************/



static void _API_ParseParam_maskSeeds( char *str, lineCmdParamMaskSeeds *p )
{
  char *proc = "_API_ParseParam_maskSeeds";
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

  if ( _debug_ ) {
      fprintf( stderr, "%s: translation from\n", proc );
      fprintf( stderr, "   '%s'\n", str );
      fprintf( stderr, "into\n" );
      for ( i=0; i<argc; i++ )
          fprintf( stderr, "   argv[%2d] = '%s'\n", i, argv[i] );
  }

  API_ParseParam_maskSeeds( 0, argc, argv, p );

  free( argv );
}





static int _n_call_parse_ = 0;

void API_ParseParam_maskSeeds( int firstargc, int argc, char *argv[],
                                  lineCmdParamMaskSeeds *p )
{
  int i;
  int outputisread = 0;
  char text[STRINGLENGTH];

  _n_call_parse_ ++;

  /* option line parsing
   */
  for ( i=firstargc; i<argc; i++ ) {

      /* strings beginning with '-'
       */
      if ( argv[i][0] == '-' ) {
          if ( argv[i][1] == '\0' ) {
              API_ErrorParse_maskSeeds( (char*)NULL, "parsing '-' ...\n", 0 );
          }

          /* file names
           */
          else if ( strcmp ( argv[i], "-seed-image" ) == 0
                    || (strcmp ( argv[i], "-si" ) == 0 && argv[i][3] == '\0') ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_maskSeeds( (char*)NULL, "parsing -seed-image ...\n", 0 );
            (void)strcpy( p->input_seed_name, argv[i] );
          }

          else if ( strcmp ( argv[i], "-cell-image" ) == 0
                    || (strcmp ( argv[i], "-ci" ) == 0 && argv[i][3] == '\0') ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_maskSeeds( (char*)NULL, "parsing -cell-image ...\n", 0 );
            (void)strcpy( p->input_cell_name, argv[i] );
          }

          else if ( strcmp ( argv[i], "-output-image" ) == 0
                    || (strcmp ( argv[i], "-oi" ) == 0 && argv[i][3] == '\0')
                    || (strcmp ( argv[i], "-o" ) == 0 && argv[i][2] == '\0') ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_maskSeeds( (char*)NULL, "parsing -output-image ...\n", 0 );
            (void)strcpy( p->output_name, argv[i] );
            outputisread = 1;
          }

          /* general parameters
           */
          else if ( (strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0')
                    || (strcmp ( argv[i], "--help" ) == 0 && argv[i][6] == '\0') ) {
             API_ErrorParse_maskSeeds( (char*)NULL, (char*)NULL, 1);
          }
          else if ( (strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0') ) {
             API_ErrorParse_maskSeeds( (char*)NULL, (char*)NULL, 0);
          }
          else if ( strcmp ( argv[i], "-verbose" ) == 0
                    || (strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0') ) {
            if ( _n_call_parse_ == 1 ) {
              if ( _verbose_ <= 0 ) _verbose_ = 1;
              else                  _verbose_ ++;
              if ( _VT_VERBOSE_ <= 0 ) _VT_VERBOSE_ = 1;
              else                     _VT_VERBOSE_ ++;
            }
          }
          else if ( strcmp ( argv[i], "-noverbose" ) == 0
                    || (strcmp ( argv[i], "-nv" ) == 0 && argv[i][3] == '\0') ) {
              _verbose_ = 0;
              _VT_VERBOSE_ = 0;
          }
          else if ( (strcmp ( argv[i], "-debug" ) == 0 && argv[i][6] == '\0')
                    || (strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0') ) {
            if ( _n_call_parse_ == 1 ) {
              _debug_ = 1;
              _VT_DEBUG_ = 1;
            }
          }
          else if ( (strcmp ( argv[i], "-no-debug" ) == 0 && argv[i][9] == '\0')
                    || (strcmp ( argv[i], "-nodebug" ) == 0 && argv[i][8] == '\0') ) {
              _debug_ = 0;
              _VT_DEBUG_ = 0;
          }

          else if ( strcmp ( argv[i], "-print-parameters" ) == 0
                    || (strcmp ( argv[i], "-param" ) == 0 && argv[i][6] == '\0') ) {
             p->print_lineCmdParam = 1;
          }

          else if ( strcmp ( argv[i], "-print-time" ) == 0
                     || (strcmp ( argv[i], "-time" ) == 0 && argv[i][5] == '\0') ) {
             p->print_time = 1;
          }
          else if ( (strcmp ( argv[i], "-notime" ) == 0 && argv[i][7] == '\0')
                      || (strcmp ( argv[i], "-no-time" ) == 0 && argv[i][8] == '\0') ) {
             p->print_time = 0;
          }

          /* unknown option
           */
          else {
              sprintf(text,"unknown option %s\n",argv[i]);
              API_ErrorParse_maskSeeds( (char*)NULL, text, 0);
          }
      }

      /* strings beginning with a character different from '-'
       */
      else {
          if ( strlen( argv[i] ) >= STRINGLENGTH ) {
              fprintf( stderr, "... parsing '%s'\n", argv[i] );
              API_ErrorParse_maskSeeds( (char*)NULL, "too long file name ...\n", 0 );
          }
          else if ( outputisread == 0 ) {
              (void)strcpy( p->output_name, argv[i] );
              outputisread = 1;
          }
          else {
              fprintf( stderr, "... parsing '%s'\n", argv[i] );
              API_ErrorParse_maskSeeds( (char*)NULL, "too many file names, parsing '-'...\n", 0 );
          }
      }
  }

}
