/*************************************************************************
 * tube3D.c -
 *
 * $Id: vessel.c,v 1.2 2006/05/16 09:33:34 greg Exp $
 *
 * Copyright (c) INRIA 2000
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * 
 *
 * ADDITIONS, CHANGES
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>

#include <chunks.h>

#include <vt_common.h>
#include <vt_tube3D.h>
#include <vt_tube2Dmatlab.h>



typedef enum {
  MULTI_SCALE,
  MATLAB_VECTORS
} enumComputation;



typedef struct local_par {
  vt_names names;
  int type;

  int writeImages;
  enumComputation typeComputation;

  enumStructureColor structureColor;

  double scale1;
  double scale2;
  int nbscales;

  int dimension;
  methodType mode;

  int print_time;

} local_par;





/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );

static double _GetTime();
static double _GetClock();
static char *_BaseName( char *p );


static char *usage = "[image-in] [image-out]\n\
 [-init %lf] [-last %lf] [-nb %d] [-wi] [-2D]\n\
 [-krissian|-frangi] [-matlab] [-wi] [-black|-noir|-white|-blanc]\n\
 [-parallel|-no-parallel] [-max-chunks %d]\n\
 [-parallelism-type|-parallel-type default|none|openmp|omp|pthread|thread]\n\
 [-omp-scheduling|-omps default|static|dynamic-one|dynamic|guided]\n\
 [-time] [-notime]\n\
 [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -init %lf : echelle (sigma) fine\n\
\t -last %lf : echelle (sigma) grossiere\n\
\t -nb   %d  : nombre d'echelles (interpolation logarithmique)\n\
\t -wi : ecrit toutes les images intermediaires\n\
\t       theta, phi : angles definissant le vecteur 3D\n\
\t       rep : reponse 3D multi-echelle\n\
\t       scale : echelle de la reponse maximale\n\
\t -krissian : Krissian multiscale response based on Hessian eigenvectors and\n\
\t             Gradients (default mode)\n\
\t -frangi : Frangi multiscale response based on Hessian eigenvalues\n\
\t -matlab : vecteurs 3D en matlab\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image;
  vt_3Dimres imsRes;
  int flag_3D = 1;
  
  double time_init = _GetTime();
  double time_exit;
  double clock_init = _GetClock();
  double clock_exit;


  /* pas de parallelisme
     par defaut
  */
  setMaxChunks( 1 );

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
  

  if ( par.dimension == 2 || image->dim.z == 1 )
    flag_3D = 0;
  
  if ( VT_Alloc3DImres( &imsRes, image, par.names.out, flag_3D ) != 1 ) {
    VT_ErrorParse("unable to allocate response images\n", 0 );
  }







  switch ( par.typeComputation ) {
  default :
  case MULTI_SCALE :

    if ( flag_3D ) {

      if ( VT_Compute3DMultiScale( image, &imsRes, par.scale1,
				   par.scale2, par.nbscales, 
				   par.structureColor, par.mode ) != 1 ) {
	VT_ErrorParse("unable to compute response\n", 0 );
      }
      if ( par.writeImages ) VT_Write3DImres( &imsRes, flag_3D );
      
      VT_Compute3DExtrema( &imsRes, &(imsRes.imTheta) );
      sprintf( imsRes.imTheta.name, "%s.ext.inr", par.names.out );
      VT_WriteInrimage( &(imsRes.imTheta) );

    }
    else {

      if ( VT_Compute2DMultiScale( image, &imsRes, par.scale1,
				   par.scale2, par.nbscales, 
				   par.structureColor, par.mode ) != 1 ) {
	VT_ErrorParse("unable to compute response\n", 0 );
      }
      if ( par.writeImages ) VT_Write3DImres( &imsRes, flag_3D );

      VT_Compute2DExtrema( &imsRes, &(imsRes.imTheta) );
      sprintf( imsRes.imTheta.name, "%s.ext.inr", par.names.out );
      VT_WriteInrimage( &(imsRes.imTheta) );
    }
    

    break;



  case MATLAB_VECTORS :
    {
      char name[DOUBLESTRINGLENGTH];
      FILE *f, *fopen();
      int i, j, fd;


      
      i = strlen( par.names.out )-1;
      for ( ; i >= 0 && par.names.out[i] != '/' ; i-- ) 
	;

      sprintf( name, "%s.raw", par.names.out );
      fd = creat( name, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );

      sprintf( name, "%s.m", par.names.out );
      f = fopen( name, "w" );
      
      fprintf( f, "\n" );
      fprintf( f, "\n" );
      fprintf( f, "\n" );
      fprintf( f, "%%\n" );
      fprintf( f, "%%" );
      for ( j=0; j<argc; j++ ) fprintf( f, " %s", argv[j] );
      fprintf( f, "\n" );
      fprintf( f, "%%\n" );
      fprintf( f, "\n" );
      fprintf( f, "\n" );
      fprintf( f, "\n" );

      fprintf( f, "echo off\n" );
      fprintf( f, "fid = fopen('%s.raw', 'r' );\n", &(par.names.out[i+1]) );

      fprintf( f, "figure;\n" );
      fprintf( f, "hold on;\n" );
      VT_3DDrawImage( image, fd, f );

      if ( VT_Compute3DMultiScale( image, &imsRes, par.scale1,
				 par.scale2, par.nbscales, 
				   par.structureColor, par.mode ) != 1 ) {
	VT_ErrorParse("unable to compute response\n", 0 );
      }

      VT_3DDrawWeightedVectors( &(imsRes.imRep), &(imsRes.imTheta), 
				&(imsRes.imPhi), fd, f, 2.0 );


      fprintf( f, "hold off;\n" );
      fprintf( f, "fclose( fid );\n" );
      close( fd );
      fclose( f );
    }
  }



  
  /*--- liberations memoires ---*/
  VT_Free3DImres( &imsRes );
  VT_FreeImage( image );
  VT_Free( (void**)&image );



  time_exit = _GetTime();
  clock_exit = _GetClock();

  if ( par.print_time ) { 
    fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName( program ), time_exit - time_init );
    fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
    fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );
  }


  return( 1 );
}








static void VT_Parse( int argc, 
		      char *argv[], 
		      local_par *par )
{
  int i, nb, status;
  int maxchunks;
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



      
      else if ( strcmp ( argv[i], "-matlab" ) == 0 ) {
	par->typeComputation = MATLAB_VECTORS;
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


      else if ( strcmp ( argv[i], "-krissian" ) == 0 ) {
	par->mode = KRISSIAN;
      }

      else if ( strcmp ( argv[i], "-frangi" ) == 0 ) {
	par->mode = FRANGI;
      }




      else if ( strcmp ( argv[i], "-init" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -init...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->scale1) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -init...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-last" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -last...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->scale2) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -last...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-nb" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -nb...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->nbscales) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -nb...\n", 0 );
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
	if ( i >= argc)    VT_ErrorParse( "-parallelism-type", 0 );
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
	else {
	  fprintf( stderr, "unknown parallelism type: '%s'\n", argv[i] );
	  VT_ErrorParse( "-parallelism-type", 0 );
	}
      }
      

      else if ( strcmp ( argv[i], "-max-chunks" ) == 0 ) {
	i ++;
	if ( i >= argc)    VT_ErrorParse( "-max-chunks", 0 );
	status = sscanf( argv[i], "%d", &maxchunks );
	if ( status <= 0 ) VT_ErrorParse( "-max-chunks", 0 );
	if ( maxchunks >= 1 ) setMaxChunks( maxchunks );
      }
      
      else if ( strcmp ( argv[i], "-omp-scheduling" ) == 0 || 
	      ( strcmp ( argv[i], "-omps" ) == 0 && argv[i][3] == '\0') ) {
	i ++;
	if ( i >= argc)    VT_ErrorParse( "-parallel-scheduling", 0 );
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
	  VT_ErrorParse( "-omp-scheduling", 0 );
	}
      }



      else if ( (strcmp ( argv[i], "-time" ) == 0 && argv[i][5] == '\0') ) {
	par->print_time = 1;
      }
      else if ( (strcmp ( argv[i], "-notime" ) == 0 && argv[i][7] == '\0') 
		|| (strcmp ( argv[i], "-no-time" ) == 0 && argv[i][8] == '\0') ) {
	par->print_time = 0;
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

  par->writeImages = 0;
  par->typeComputation = 0;

  par->structureColor = _WHITE_;

  par->scale1 = 1.0;
  par->scale2 = 1.0;
  par->nbscales = 1;

  par->dimension = 3;
  par->mode = KRISSIAN;

  par->print_time = 0;
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
