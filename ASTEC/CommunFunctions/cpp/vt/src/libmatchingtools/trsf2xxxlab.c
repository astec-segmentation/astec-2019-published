/*************************************************************************
 * trsf2xxxlab.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Sun Feb 10 16:59:14 CET 2013
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */

#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>

#include <fcntl.h>
#include <unistd.h>

#include <string.h>
#include <sys/time.h> /* gettimeofday() */
#include <time.h> /* clock() */
#include <math.h>

#include <histogram.h>

#include <bal-stddef.h>
#include <bal-image.h>
#include <bal-transformation.h>
#include <bal-transformation-tools.h>



static int _time_ = 0;
static int _clock_ = 0;

static int _verbose_ = 1;
static int _debug_ = 0;






typedef enum {
  _ANGLE_,
  _TX_,
  _TY_,
  _TZ_,
  _A11_,
  _A12_,
  _A13_,
  _A21_,
  _A22_,
  _A23_,
  _A31_,
  _A32_,
  _A33_
} enumTrsfParameter;


typedef struct local_par {
  char *theformat;
  char *resname;

  int first;
  int last;

  enumHistogramFile xxxlab;

  enumTrsfParameter parameter;

} local_par;




/*------- Definition des fonctions statiques ----------*/
static void _Parse( int argc, char *argv[], local_par *par );
static void _ErrorParse( char *str, int l );
static void _InitParam( local_par *par );
static char *_BaseName( char *p );
static double _GetTime();
static double _GetClock();




static char *usage = "[format-in] [file-out]\n\
 [-first %d] [-last %d]\n\
 [-parameter angle | a11|a12|a13|a21|a22|a23|a31|a32|a33 | tx|ty|tz]\n\
 [-matlab|-scilab]\n\
 [-v] [-D] [-help]";

static char *detail = "\
\t -v : mode verbose\n\
\t -D : mode debug\n\
\n";

static char *program = NULL;






int main( int argc, char *argv[] )
{
  local_par par;
  double time_init, time_exit;
  double clock_init, clock_exit;

  int i;
  char filename[1024];
  bal_transformation theTrsf; 
  double *parameter = NULL;
  int *index = NULL;
  int n, nvalues = 0;
  double da;
  double min, max;

  time_init = _GetTime();
  clock_init = _GetClock();


  /*--- initialisation des parametres ---*/
  _InitParam( &par );
  
  /*--- lecture des parametres ---*/
  _Parse( argc, argv, &par );
  




  if ( par.theformat == (char*)NULL ) {
    fprintf( stderr, "%s: no input format\n", program );
    return( 0 );
  }
  if ( par.last - par.first + 1 <= 0 ) {
    fprintf( stderr, "%s: last < first \n", program );
    return( 0 );
  }
    
  parameter = (double*)malloc( (par.last - par.first + 1) * sizeof(double) );
  if ( parameter == (double*)NULL ) {
    fprintf( stderr, "%s: unable to allocate array (1)\n", program );
    return( 0 );
  }
  index = (int*)malloc( (par.last - par.first + 1) * sizeof(int) );
  if ( index == (int*)NULL ) {
    free( parameter );
    fprintf( stderr, "%s: unable to allocate array (2)\n", program );
    return( 0 );
  }
  
  for ( nvalues=0, i=par.first; i<=par.last; i++ ) {
    sprintf( filename, par.theformat, i );
    if ( _debug_ )
      fprintf( stderr, "#%4d: processing '%s'\n", i, filename );
    if ( BAL_ReadTransformation( &theTrsf, filename ) != 1 ) {
      if ( _verbose_ )
         fprintf( stderr, "%s: unable to read transformation '%s'\n", program, filename );
      continue;
    }

    switch ( par.parameter ) {

    default :
      
      free( index );
      free( parameter );
      fprintf( stderr, "%s: such parameter not handled yet\n", program );
      return( 0 );

    case _ANGLE_ :
      if ( BAL_RotationAngle( &theTrsf, &da ) != 1 ) {
        free( index );
        free( parameter );
        fprintf( stderr, "%s: unable to compute angle from transformation '%s'\n", program, filename );
        return( 0 );
      }
      parameter[nvalues] = da;
      break;
      
    case _TX_ :
      parameter[nvalues] =  theTrsf.mat.m[3];
      break;

    case _TY_ :
      parameter[nvalues] =  theTrsf.mat.m[7];
      break;

    case _TZ_ :
      parameter[nvalues] =  theTrsf.mat.m[11];
      break;

    case _A11_ :
      parameter[nvalues] =  theTrsf.mat.m[0];
      break;

    case _A12_ :
      parameter[nvalues] =  theTrsf.mat.m[1];
      break;

    case _A13_ :
      parameter[nvalues] =  theTrsf.mat.m[2];
      break;

    case _A21_ :
      parameter[nvalues] =  theTrsf.mat.m[4];
      break;

    case _A22_ :
      parameter[nvalues] =  theTrsf.mat.m[5];
      break;

    case _A23_ :
      parameter[nvalues] =  theTrsf.mat.m[6];
      break;

    case _A31_ :
      parameter[nvalues] =  theTrsf.mat.m[8];
      break;

    case _A32_ :
      parameter[nvalues] =  theTrsf.mat.m[9];
      break;

    case _A33_ :
      parameter[nvalues] =  theTrsf.mat.m[10];
      break;

    }

    fprintf( stderr, "#%4d: %f\n", i, parameter[nvalues] );
    index[nvalues] = i;
    nvalues ++;
  }

  if ( nvalues <= 0 ) {
    free( index );
    free( parameter );
    fprintf( stderr, "%s: no read values\n", program );
    return( 0 );
  }

  if ( par.resname == (char*)NULL ) {
    free( index );
    free( parameter );
    fprintf( stderr, "%s: no output name\n", program );
    return( 0 );
  }


  min = max = parameter[0];
  for ( n=1; n<nvalues; n++ ) {
    if ( min > parameter[n] ) min = parameter[n];
    if ( max < parameter[n] ) max = parameter[n];
  }


  {
    int fd;
    FILE *f;

    sprintf( filename, "%s.raw", par.resname );

    fd = open(filename, O_CREAT | O_TRUNC | O_WRONLY, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
    if ( fd == -1 ) {
      free( index );
      free( parameter );
      fprintf( stderr, "%s: unable to open '%s' for writing\n", program, filename );
      return( 0 );
    }

    if ( write( fd, index, (nvalues) * sizeof(int) ) == -1 ) {
      fprintf( stderr, "%s: unable to write in '%s'\n", program, filename );
      free( index );
      free( parameter );
      close( fd );
      return( 0 );
    }

    if ( write( fd, parameter, (nvalues) * sizeof(double) ) == -1 ) {
      fprintf( stderr, "%s: unable to write in '%s'\n", program, filename );
      free( index );
      free( parameter );
      close( fd );
      return( 0 );
    }

    close( fd );
    free( index );
    free( parameter );


    switch( par.xxxlab ) {

    default :
    case _SCILAB_ :
      
      sprintf( filename, "%s.sce", par.resname );
      f = fopen( filename, "w" );
      
      if ( f == (FILE*)NULL ) {
        fprintf( stderr, "%s: unable to open '%s' for writing\n", program, filename );
        return( 0 );
      }

      fprintf( f, "\n" );
      fprintf( f, "f=mopen('%s.raw','r');\n", _BaseName( par.resname ) );
      fprintf( f, "I%s=mget( %d, 'i', f);\n", _BaseName( par.resname ), nvalues );
      fprintf( f, "A%s=mget( %d, 'd', f);\n", _BaseName( par.resname ), nvalues );
      fprintf( f, "mclose(f);\n" );
      fprintf( f, "\n" );
      
      fprintf( f, "\n" );
      fprintf( f, "figure;\n" );
      fprintf( f, "plot( I%s, A%s );\n", _BaseName( par.resname ), _BaseName( par.resname ) );
      fprintf( f, "\n" );
      fprintf( f, "// ymin= %lf\n", min );
      fprintf( f, "// ymax= %lf\n", max );
      fprintf( f, "\n" );

      fprintf( f, "\n" );
      fprintf( f, "// xs2jpg(gcf(),'FIG%s.jpg');\n", _BaseName( par.resname )  );
      fprintf( f, "xs2png(gcf(),'FIG%s.png');\n", _BaseName( par.resname )  );
      fprintf( f, "\n" );

      fclose( f );
      
      break;

    case _MATLAB_ :
      
      sprintf( filename, "%s.m", par.resname );
      f = fopen( filename, "w" );

      if ( f == (FILE*)NULL ) {
        fprintf( stderr, "%s: unable to open '%s' for writing\n", program, filename );
        return( 0 );
      }

      fprintf( f, "\n" );
      fprintf( f, "f=fopen('%s.raw','r');\n", _BaseName( par.resname ) );
      fprintf( f, "I%s=fread( f, [%d], int%lu);\n", _BaseName( par.resname ),
               nvalues, 8*sizeof(int) );
      fprintf( f, "A%s=fread( f, [%d], float%lu);\n", _BaseName( par.resname ),
               nvalues, 8*sizeof(double) );
      fprintf( f, "fclose(f);\n" );
      fprintf( f, "\n" );
      
      fprintf( f, "\n" );
      fprintf( f, "figure;\n" );
      fprintf( f, "plot( I%s, A%s );\n", _BaseName( par.resname ), _BaseName( par.resname ) );
      fprintf( f, "\n" );
      fclose( f );

      break;

    }

  }

  
  time_exit = _GetTime();
  clock_exit = _GetClock();

  if ( _time_ ) 
    fprintf( stderr, "%s: elapsed time = %f\n", program, time_exit - time_init );

  if ( _clock_ ) 
    fprintf( stderr, "%s: elapsed time = %f\n", program, clock_exit - clock_init );

  return( 0 );
}








static void _Parse( int argc, 
                      char *argv[],
                      local_par *par )
{
  int i, status;
  
  program = argv[0];

  if ( argc == 1 ) _ErrorParse("\n", 0 );
  
  /*--- lecture des parametres ---*/
  i = 1;
  while ( i < argc ) {
    if ( argv[i][0] == '-' ) {
      if ( (strcmp ( argv[i], "-f" ) == 0 && argv[i][2] == '\0')
           || strcmp ( argv[i], "-first" ) == 0 ) {
        i += 1;
        if ( i >= argc)    _ErrorParse( "parsing -first...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->first) );
        if ( status <= 0 ) _ErrorParse( "parsing -first...\n", 0 );
      }

      else if ( (strcmp ( argv[i], "-l" ) == 0 && argv[i][2] == '\0')
                || strcmp ( argv[i], "-last" ) == 0 ) {
        i += 1;
        if ( i >= argc)    _ErrorParse( "parsing -last...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->last) );
        if ( status <= 0 ) _ErrorParse( "parsing -last...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-matlab" ) == 0 ) {
        par->xxxlab = _MATLAB_;
      }

      else if ( strcmp ( argv[i], "-scilab" ) == 0 ) {
        par->xxxlab = _SCILAB_;
      }

      else if ( strcmp ( argv[i], "-parameter" ) == 0 
                || (strcmp ( argv[i], "-param" ) == 0 && argv[i][6] == '\0')
                || (strcmp ( argv[i], "-p" ) == 0 && argv[i][2] == '\0') ) {
        i += 1;
        if ( i >= argc)    _ErrorParse( "parsing -parameter...\n", 0 );
        if ( strcmp ( argv[i], "angle" ) == 0 ) {
          par->parameter = _ANGLE_;
        }
        else if ( strcmp ( argv[i], "tx" ) == 0 ) {
          par->parameter = _TX_;
        }
        else if ( strcmp ( argv[i], "ty" ) == 0 ) {
          par->parameter = _TY_;
        }
        else if ( strcmp ( argv[i], "tz" ) == 0 ) {
          par->parameter = _TZ_;
        }
        else if ( strcmp ( argv[i], "a11" ) == 0 ) {
          par->parameter = _A11_;
        }
        else if ( strcmp ( argv[i], "a12" ) == 0 ) {
          par->parameter = _A12_;
        }
        else if ( strcmp ( argv[i], "a13" ) == 0 ) {
          par->parameter = _A13_;
        }
        else if ( strcmp ( argv[i], "a21" ) == 0 ) {
          par->parameter = _A21_;
        }
        else if ( strcmp ( argv[i], "a22" ) == 0 ) {
          par->parameter = _A22_;
        }
        else if ( strcmp ( argv[i], "a23" ) == 0 ) {
          par->parameter = _A23_;
        }
        else if ( strcmp ( argv[i], "a31" ) == 0 ) {
          par->parameter = _A31_;
        }
        else if ( strcmp ( argv[i], "a32" ) == 0 ) {
          par->parameter = _A32_;
        }
        else if ( strcmp ( argv[i], "a33" ) == 0 ) {
          par->parameter = _A33_;
        }
        else {
          _ErrorParse( "parsing -parameter: unknown parameter...\n", 0 );
        }
      }


      /*--- arguments generaux ---*/
      else if ( strcmp ( argv[i], "-help" ) == 0 ) {
        _ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 ) {
        _verbose_ = 1;
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 ) {
        _debug_ = 1;
      }
      

      /*--- option inconnue ---*/
      else {
        fprintf(stderr,"unknown option: '%s'\n",argv[i]);
      }
    }
    /*--- saisie des noms d'images ---*/
    else {
      if ( par->theformat == NULL ) {
        par->theformat = argv[i];
      }
      else if ( par->resname == NULL ) {
        par->resname = argv[i];
      }
      else {
        fprintf(stderr,"too many file names: '%s'\n",argv[i]);
      }
    }
    i += 1;
  }
  

}






static void _ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"Erreur : %s",str);
  exit( 1 );
}








static void _InitParam( local_par *par )
{
  par->theformat = (char*)NULL;
  par->resname = (char*)NULL;
  par->first = 0;
  par->last = 0;

  par->xxxlab = _SCILAB_;
  par->parameter = _ANGLE_;
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
