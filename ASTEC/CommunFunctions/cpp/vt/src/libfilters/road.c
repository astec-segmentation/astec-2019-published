/*************************************************************************
 * road.c -
 *
 * $Id: road.c,v 1.9 2006/05/16 09:33:34 greg Exp $
 *
 * Copyright (c) INRIA 1999
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

#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>

#include <vt_common.h>
#include <vt_tube2D.h>
#include <vt_tube2Dmatlab.h>


typedef enum {
  ONE_SCALE,
  MATLAB_VECTORS
} enumComputation;



typedef struct local_par {
  vt_names names;
  int type;

  int writeImages;
  enumComputation typeComputation;

  enumStructureColor structureColor;
  

  double theCoefficient;

  double theTensorLargeCoefficient;
  double theTensorSmallCoefficient;
  double theTensorMultiCoefficient;


} local_par;





/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[image-in] [image-out]\n\
\t [-matlab] [-wi] [-black|-noir|-white|-blanc]\n\
\t [-sigma %lf] [-tsigma %lf] [-mult-tsigma|-mts %lf]\n\
\t [-large-tsigma|-lts %lf] [-small-tsigma|-sts %lf]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n\
\t options-de-type : -o 1    : unsigned char\n\
\t                   -o 2    : unsigned short int\n\
\t                   -o 2 -s : short int\n\
\t                   -o 4 -s : int\n\
\t                   -r      : float\n\
\t si aucune de ces options n'est presente, on prend le type de 'image-in'\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image;

  float theCoeffs[3];
  vt_2Dimages theIms;
  vt_2Dtensor theTensor;
  vt_image theExtrema;
  


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
  
  

  switch( par.structureColor ) {
  case _BLACK_ :
    fprintf( stderr, "recherche de structures noires\n" );
    break;
  case _WHITE_ :
    fprintf( stderr, "recherche de structures blanches\n" );
    break;
  }

  theCoeffs[0] = theCoeffs[1] = theCoeffs[2] = par.theCoefficient;




  switch ( par.typeComputation ) {
  default :
  case ONE_SCALE :

    VT_Alloc2Dimages( &theIms, image, par.names.out );

    VT_Filter2Dimages( image, &theIms, theCoeffs );
    VT_Compute2DEigenVectors( &(theIms.hessien), &(theIms.ime) );
    VT_FilterOn2DEigenValues( &theIms, par.structureColor );
    VT_Compute2DResponse( &theIms, par.structureColor, 1.0, par.theCoefficient );
    VT_Compute2DMaskedExtrema( &(theIms.imr), &(theIms.hessien.imtheta1), 
			 &(theIms.ime) );

    {
      char name[DOUBLESTRINGLENGTH];
      sprintf( name, "%s.tensor", par.names.out );
      VT_Alloc2Dtensor( &theTensor, image, name );
      sprintf( name, "%s.tensor.extrema.inr", par.names.out );
      VT_InitFromImage( &theExtrema, image, name, UCHAR );
      VT_AllocImage( &theExtrema );
    }

    /*
    VT_2DTensorGaussianVoting( &theTensor, &(theIms.imr), 
			       &(theIms.hessien.imtheta2), 
			       &(theIms.ime), par.theTensorCoefficient );
    */
    VT_2DTensorVoting( &theTensor, &(theIms.imr), 
		       &(theIms.hessien.imtheta2), 
		       &(theIms.ime), 
		       par.theTensorLargeCoefficient, 
		       par.theTensorSmallCoefficient,
		       par.theTensorMultiCoefficient );

    VT_Compute2DEigenVectors( &theTensor, (vt_image *)NULL );

    VT_Compute2DMaskedExtrema( &(theTensor.imvp1), &(theTensor.imtheta2), 
			 &theExtrema );

    VT_Write2Dtensor( &theTensor );
    VT_Write2Dimages( &theIms );
    VT_WriteInrimage( &theExtrema );

    VT_Free2Dimages( &theIms );

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
      
      
      VT_Alloc2Dimages( &theIms, image, par.names.out );
      VT_Filter2Dimages( image, &theIms, theCoeffs );
      VT_Compute2DEigenVectors( &(theIms.hessien), &(theIms.ime) );
      VT_Compute2DResponse( &theIms, par.structureColor, 1.0, par.theCoefficient );


      VT_2DDrawImage( image, 0, fd, f );
      VT_2DDrawWeightedVectors( &(theIms.hessien.imvp1), &(theIms.hessien.imtheta2), 
				fd, f ); 
      fprintf( f, "title('originale + (grde |val prop.|)*(petit |vec prop.|)');\n" );
      fprintf( f, "hold off;\n" );





      fprintf( f, "figure;\n" );
      fprintf( f, "hold on;\n" );
      VT_2DDrawImage( image, fd, 0, f );
      VT_2DDrawWeightedVectors( &(theIms.imr), &(theIms.hessien.imtheta2), 
					fd, f ); 
      fprintf( f, "title('originale + (reponse)*(petit |vec prop.|)');\n" );
      fprintf( f, "hold off;\n" );









      fprintf( f, "figure;\n" );
      fprintf( f, "hold on;\n" );

      sprintf( name, "%s.tensor", par.names.out );
      VT_Alloc2Dtensor( &theTensor, image, name );

      /*
      VT_2DTensorGaussianVoting( &theTensor, &(theIms.imr), 
	                         &(theIms.hessien.imtheta2), 
				 &(theIms.ime), par.theTensorCoefficient );
      */
      VT_2DTensorVoting( &theTensor, &(theIms.imr), 
			 &(theIms.hessien.imtheta2), 
			 &(theIms.ime), 
			 par.theTensorLargeCoefficient, 
			 par.theTensorSmallCoefficient,
			 par.theTensorMultiCoefficient );

      VT_Compute2DEigenVectors( &theTensor, (vt_image *)NULL );
      VT_2DDrawImage( &(theTensor.imvp1), fd, 0, f );
      VT_2DDrawWeightedVectors( &(theTensor.imvp1), &(theTensor.imtheta1), fd, f );

      fprintf( f, "title('grde |val prop.| tenseur| + (id.)*(grd |vec prop.|)');\n" );
      fprintf( f, "hold off;\n" );




      fprintf( f, "figure;\n" );
      fprintf( f, "hold on;\n" );
      VT_2DDrawImage( &(theTensor.imvp1), fd, 0, f );
      VT_2DDrawWeightedVectors( &(theTensor.imvp1), &(theTensor.imtheta2), fd, f );

      fprintf( f, "title('grde |val prop.| tenseur| + (id.)*(petit |vec prop.|)');\n" );
      fprintf( f, "hold off;\n" );




      if ( par.writeImages ) VT_Write2Dimages( &theIms );
      VT_Free2Dimages( &theIms );

      fprintf( f, "fclose( fid );\n" );
      close( fd );
      fclose( f );
      if ( 0 ) {
	fprintf( stderr, "Pour voir l'image originale dans la meme geometrie\n" );
	fprintf( stderr, "Tran -db %s | Tran -sy | xviewer & \n", image->name );
      }
    }
  }



  
  /*--- liberations memoires ---*/
  VT_FreeImage( image );
  VT_Free( (void**)&image );
  return( 1 );
}








static void VT_Parse( int argc, 
		      char *argv[], 
		      local_par *par )
{
  int i, nb, status;
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
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 ) {
	_VT_DEBUG_ = 1;
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



      
      else if ( strcmp ( argv[i], "-black" ) == 0 || strcmp ( argv[i], "-noir" ) == 0 ) {
	par->structureColor = _BLACK_;
      }
      else if ( strcmp ( argv[i], "-white" ) == 0 || strcmp ( argv[i], "-blanc" ) == 0 ) {
	par->structureColor = _WHITE_;
      }





      else if ( strcmp ( argv[i], "-sigma" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -sigma...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->theCoefficient) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -sigma...\n", 0 );
      }



      else if ( strcmp ( argv[i], "-tsigma" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -tsigma...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->theTensorLargeCoefficient) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -tsigma...\n", 0 );
	par->theTensorSmallCoefficient = par->theTensorLargeCoefficient;
      }
      
      else if ( strcmp ( argv[i], "-mult-tsigma" ) == 0 || 
		strcmp ( argv[i], "-mts" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -mult-tsigma...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->theTensorMultiCoefficient) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -mult-tsigma...\n", 0 );
      }
      
      else if ( strcmp ( argv[i], "-large-tsigma" ) == 0 || 
		strcmp ( argv[i], "-lts" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -large-tsigma...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->theTensorLargeCoefficient) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -large-tsigma...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-small-tsigma" ) == 0  || 
		strcmp ( argv[i], "-sts" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -small-tsigma...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->theTensorSmallCoefficient) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -small-tsigma...\n", 0 );
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

  par->theCoefficient = 1.0;
  par->theTensorLargeCoefficient = 1.0;
  par->theTensorSmallCoefficient = 1.0;
  par->theTensorMultiCoefficient = 3.0;
}
