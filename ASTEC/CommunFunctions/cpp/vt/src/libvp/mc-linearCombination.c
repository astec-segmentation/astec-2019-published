/*************************************************************************
 * mc-linearCombination.c -
 *
 * Copyright (c) INRIA 2019, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Mer  6 mar 2019 12:28:41 CET
 *
 * ADDITIONS, CHANGES
 *
 */

#include <sys/time.h> /* gettimeofday() */
#include <time.h> /* clock() */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <string-tools.h>
#include <typedefs.h>

#include <vt_names.h>
#include <vt_error.h>

#include <api-linearCombination.h>

static int _time_ = 1;

typedef struct local_par {
  char *resName;
  bufferType type;
} local_par;

/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par, stringList *imageNames, stringList *weightNames );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
static char *_BaseName( char *p );
static double VT_GetTime();
static double VT_GetClock();



static char *usage = "\n\
 -res %s\n\
 -weight[s] %s ... %s -image[s] %s ... %s \n\
 [-time|-notime]\n\
 [-v|-nv] [-D] [-help] [encoding-type]";

static char *detail = "\
###\n\
[-v]        # be more verbose\n\
[-D]        # some debug information (if any)\n\
[encoding-type] # for the ouput image\n\
  -o 1    : unsigned char\n\
  -o 2    : unsigned short int\n\
  -o 2 -s : short int\n\
  -o 4 -s : int\n\
  -r      : float\n\
  -type s8|u8|s16|u16|...\n\
  default is same type than 'image-in's\n\
\n";

static char program[STRINGLENGTH];




int main( int argc, char *argv[] )
{
  stringList imageNames;
  stringList weightNames;
  local_par par;

  double time_init, time_exit;
  double clock_init, clock_exit;

  time_init = VT_GetTime();
  clock_init = VT_GetClock();


  initStringList( &imageNames );
  initStringList( &weightNames );

  VT_InitParam( &par );
  VT_Parse( argc, argv, &par, &imageNames, &weightNames );

  if ( 1 ) {
    printStringList( stderr, &imageNames, "image name list" );
    printStringList( stderr, &weightNames, "weight name list" );
  }

  if ( API_linearCombination( &weightNames, &imageNames, par.resName, par.type ) != 1 ) {
    VT_ErrorParse( "error when computing\n", 0 );
  }


  time_exit = VT_GetTime();
  clock_exit = VT_GetClock();

  if ( _time_ ) {
    fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName( argv[0] ), time_exit - time_init );
    fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
    fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );
  }


  return( 0 );
}





static void VT_Parse( int argc,
                      char *argv[],
                      local_par *par, stringList *imageNames, stringList *weightNames )
{
  int i, status;
  int o=0, s=0, r=0;
  char text[STRINGLENGTH];

  if ( VT_CopyName( program, argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);
  if ( argc == 1 ) VT_ErrorParse("\n", 0 );

  /*--- lecture des parametres ---*/
  i = 1;
  while ( i < argc ) {
    if ( argv[i][0] == '-' ) {

      /*--- arguments generaux ---*/
      if ( strcmp ( argv[i], "-help" ) == 0 ) {
        VT_ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0' ) {
        _VT_VERBOSE_ = 1;
      }
      else if ( strcmp ( argv[i], "-nv" ) == 0 && argv[i][3] == '\0' ) {
        _VT_VERBOSE_ = 0;
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0' ) {
        _VT_DEBUG_ = 1;
      }
      else if ( strcmp ( argv[i], "-time" ) == 0 && argv[i][5] == '\0' ) {
        _time_ = 1;
      }
      else if ( (strcmp ( argv[i], "-notime" ) == 0 && argv[i][7] == '\0')
                || (strcmp ( argv[i], "-no-time" ) == 0 && argv[i][8] == '\0') ) {
        _time_ = 0;
      }

      else if ( strcmp ( argv[i], "-weight" ) == 0
                || strcmp ( argv[i], "-weights" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -weights...\n", 0 );
        while ( i < argc && argv[i][0] != '-' ) {
          if ( addStringToList( argv[i], weightNames) != 1 ) {
            fprintf( stderr, "error when adding '%s' to weights\n", argv[i] );
            VT_ErrorParse( "parsing -weights...\n", 0 );
          }
          i++;
        }
        i--;
      }

      else if ( strcmp ( argv[i], "-image" ) == 0
                || strcmp ( argv[i], "-images" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -images...\n", 0 );
        while ( i < argc && argv[i][0] != '-' ) {
          if ( addStringToList( argv[i], imageNames) != 1 ) {
            fprintf( stderr, "error when adding '%s' to images\n", argv[i] );
            VT_ErrorParse( "parsing -images...\n", 0 );
          }
          i++;
        }
        i--;
      }

      else if ( strcmp ( argv[i], "-res" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -res...\n", 0 );
        par->resName = argv[i];
      }


      /*--- lecture du type de l'image de sortie ---*/

      else if ( strcmp ( argv[i], "-r" ) == 0 ) {
        r = 1;
      }
      else if ( strcmp ( argv[i], "-s" ) == 0 ) {
        s = 1;
      }
      else if ( strcmp ( argv[i], "-o" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -o...\n", 0 );
        status = sscanf( argv[i],"%d",&o );
        if ( status <= 0 ) VT_ErrorParse( "parsing -o...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-type" ) == 0 && argv[i][5] == '\0' ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -type...\n", 0 );
        if ( strcmp ( argv[i], "s8" ) == 0 && argv[i][2] == '\0' ) {
           par->type = SCHAR;
        }
        else if ( strcmp ( argv[i], "u8" ) == 0 && argv[i][2] == '\0' ) {
           par->type = UCHAR;
        }
        else if ( strcmp ( argv[i], "s16" ) == 0 && argv[i][3] == '\0' ) {
          par->type = SSHORT;
        }
        else if ( strcmp ( argv[i], "u16" ) == 0 && argv[i][3] == '\0' ) {
          par->type = USHORT;
        }
        else if ( strcmp ( argv[i], "s32" ) == 0 && argv[i][3] == '\0' ) {
          par->type = SINT;
        }
        else if ( strcmp ( argv[i], "u32" ) == 0 && argv[i][3] == '\0' ) {
          par->type = UINT;
        }
        else if ( strcmp ( argv[i], "s64" ) == 0 && argv[i][3] == '\0' ) {
          par->type = SLINT;
        }
        else if ( strcmp ( argv[i], "u64" ) == 0 && argv[i][3] == '\0' ) {
          par->type = ULINT;
        }
        else if ( strcmp ( argv[i], "r32" ) == 0 && argv[i][3] == '\0' ) {
          par->type = FLOAT;
        }
        else if ( strcmp ( argv[i], "r64" ) == 0 && argv[i][3] == '\0' ) {
          par->type = DOUBLE;
        }
        else {
          VT_ErrorParse( "parsing -type: unknown type...\n", 0 );
        }
      }

      /*--- option inconnue ---*/
      else {
        sprintf(text,"unknown option %s\n",argv[i]);
        VT_ErrorParse(text, 0);
      }
    }
    /*--- saisie des noms d'images ---*/
    else {
      sprintf(text,"unknown option %s\n",argv[i]);
      VT_ErrorParse(text, 0);
    }
    i += 1;
  }

  /*--- type de l'image resultat ---*/
  if ( (o == 1) && (s == 1) && (r == 0) ) par->type = SCHAR;
  if ( (o == 1) && (s == 0) && (r == 0) ) par->type = UCHAR;
  if ( (o == 2) && (s == 1) && (r == 0) ) par->type = SSHORT;
  if ( (o == 2) && (s == 0) && (r == 0) ) par->type = USHORT;
  if ( (o == 4) && (s == 1) && (r == 0) ) par->type = SINT;
  if ( (o == 4) && (s == 0) && (r == 0) ) par->type = UINT;
  if ( (o == 0 || o == 4) && (s == 0) && (r == 1) ) par->type = FLOAT;
  if ( (o == 8) && (s == 0) && (r == 1) ) par->type = DOUBLE;

  /* if ( par->type == TYPE_UNKNOWN ) VT_Warning("no specified type", program); */
}



static void VT_ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"Erreur : %s",str);
  exit( 1 );
}



static void VT_InitParam( local_par *par )
{
  par->resName = (char*)NULL;
  par->type = TYPE_UNKNOWN;
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


static double VT_GetTime()
{
  struct timeval tv;
  gettimeofday(&tv, (void *)0);
  return ( (double) tv.tv_sec + tv.tv_usec*1e-6 );
}

static double VT_GetClock()
{
  return ( (double) clock() / (double)CLOCKS_PER_SEC );
}
