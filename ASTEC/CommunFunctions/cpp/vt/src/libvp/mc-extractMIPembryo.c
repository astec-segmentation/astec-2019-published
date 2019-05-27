/*************************************************************************
 * mc-extractMIPembryo.c -
 * OBSOLETE, use execTemplate.c instead
 *
 * Copyright (c) INRIA 1999, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mar 19 jui 2018 07:54:13 CEST
 *
 * ADDITIONS, CHANGES
 *
 */

#include <sys/time.h> /* gettimeofday() */
#include <time.h> /* clock() */

#include <connexe.h>

#include <vt_common.h>
#include <vt_mip.h>


static int _time_ = 1;
static int _clock_ = 1;

typedef struct local_par {
  vt_names names;
  bufferType type;
} local_par;




/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
static double VT_GetTime();
static double VT_GetClock();




static char *usage = "[image-in] [image-out]\n\
 [-time|-notime] [-clock|-noclock]\n\
 [-inv] [-swap] [-v|-nv] [-D] [-help] [encoding-type]";

static char *detail = "\
if 'image-in' is equal to '-', we consider stdin\n\
if 'image-out' is not specified, we consider stdout\n\
if both are not specified, we consider stdin and stdout\n\
###\n\
[-inv]      # inverse 'image-in'\n\
[-swap]     # swap bytes of 'image-in' (if encoded on 2 bytes)\n\
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
  local_par par;
  vt_image *image;
  vt_image immip;
  vt_image imbin;
  vt_3m m;
  int theDim[3];
  double time_init, time_exit;
  double clock_init, clock_exit;

  time_init = VT_GetTime();
  clock_init = VT_GetClock();


  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);
  
  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.names.inv == 1 )  VT_InverseImage( image );
  if ( par.names.swap == 1 ) VT_SwapImage( image );
  
  /*--- initialisation de l'image mip ---*/
  VT_Image( &immip );
  VT_InitImage( &immip, par.names.out, image->dim.x, image->dim.y, (int)1, image->type );
  immip.siz.x = image->siz.x;
  immip.siz.y = image->siz.y;
  if ( VT_AllocImage( &immip ) != 1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse( "unable to allocate mip image\n", 0 );
  }

  if ( VT_MIP( image, &immip, (vt_image*)NULL, (vt_image*)NULL ) != 1 ) {
      VT_FreeImage( &immip );
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse( "error when computing mip image\n", 0 );
  }

  VT_FreeImage( image );
  VT_Free( (void**)&image );

  if ( 0 ) {
      if ( VT_WriteInrimage( &immip ) == -1 ) {
        VT_FreeImage( &immip );
        VT_ErrorParse("unable to write output image\n", 0);
      }
  }

  /* moyenne
   */
  if ( VT_3m( &immip, (vt_image*)NULL, &m ) == -1 ) {
      VT_FreeImage( &immip );
      VT_ErrorParse( "error when computing mean\n", 0 );
  }

  /* composantes connexes
   */
  VT_InitFromImage( &imbin, &immip, par.names.out, UCHAR );
  if ( VT_AllocImage( &imbin ) != 1 ) {
    VT_FreeImage( &immip );
    VT_ErrorParse( "unable to allocate output image\n", 0 );
  }

  theDim[0] = imbin.dim.x;
  theDim[1] = imbin.dim.y;
  theDim[2] = imbin.dim.z;

  if ( CountConnectedComponentsWithAllParams( immip.buf, immip.type,
                                              imbin.buf, imbin.type,
                                              theDim, m.moy, 4, -1, 1, 1 ) != 1 ) {
      VT_FreeImage( &imbin );
      VT_FreeImage( &immip );
      VT_ErrorParse( "error when extracting connected component\n", 0 );
  }

  VT_FreeImage( &immip );

  if ( VT_WriteInrimage( &imbin ) == -1 ) {
    VT_FreeImage( &immip );
    VT_ErrorParse("unable to write output image\n", 0);
  }
  
  VT_FreeImage( &immip );

  time_exit = VT_GetTime();
  clock_exit = VT_GetClock();

  if ( _time_ ) 
    fprintf( stderr, "%s: elapsed time = %f\n", program, time_exit - time_init );

  if ( _clock_ ) 
    fprintf( stderr, "%s: elapsed time = %f\n", program, clock_exit - clock_init );

  return( 0 );
}








static void VT_Parse( int argc, 
		      char *argv[], 
		      local_par *par )
{
  int i, nb, status;
  int o=0, s=0, r=0;
  char text[STRINGLENGTH];
  
  if ( VT_CopyName( program, argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);
  if ( argc == 1 ) VT_ErrorParse("\n", 0 );
  
  /*--- lecture des parametres ---*/
  i = 1; nb = 0;
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
	VT_ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0' ) {
	_VT_VERBOSE_ = 1;
        incrementVerboseInConnexe( );
      }
      else if ( strcmp ( argv[i], "-nv" ) == 0 && argv[i][3] == '\0' ) {
	_VT_VERBOSE_ = 0;
        setVerboseInConnexe( 0 );
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0' ) {
	_VT_DEBUG_ = 1;
        incrementDebugInConnex( );
      }
      else if ( strcmp ( argv[i], "-time" ) == 0 && argv[i][5] == '\0' ) {
	_time_ = 1;
      }
      else if ( strcmp ( argv[i], "-notime" ) == 0 && argv[i][7] == '\0' ) {
	_time_ = 0;
      }
      else if ( strcmp ( argv[i], "-clock" ) == 0 && argv[i][6] == '\0' ) {
	_clock_ = 1;
      }
      else if ( strcmp ( argv[i], "-noclock" ) == 0 && argv[i][8] == '\0' ) {
	_clock_ = 0;
      }


      /*--- traitement eventuel de l'image d'entree ---*/
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
	par->names.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
	par->names.swap = 1;
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
	VT_ErrorParse("too much file names when parsing\n", 0 );
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
  VT_Names( &(par->names) );
  par->type = TYPE_UNKNOWN;
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
