/*************************************************************************
 * minimum.c -
 *
 * $Id: test-contours.c,v 1.2 2002/09/05 17:15:06 greg Exp $
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
#include <vt_isocontours.h>

#include <vt_tube2Dmatlab.h>
#include <vt_contoursMatlab.h>


typedef struct local_par {
  int slice;
  double threshold;
  vt_names names;
  int type;
} local_par;




/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[image-in] [image-out]\n\
\t [-slice %d] [-th %lf] [-matlab %s]\n\
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
\t si aucune de ces options n'est presente, on prend le type de 'image-in'\n\
\n\
 $Revision: 1.2 $ $Date: 2002/09/05 17:15:06 $ $Author: greg $\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image, imres;

  typeStructure structure;
  typeSlice *slice;
  int n, ncontours = 0;
  int z;

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



  
  /*--- initialisation de l'image resultat ---*/
  VT_Image( &imres );
  VT_InitFromImage( &imres, image, par.names.out, UCHAR );

  if ( par.type != TYPE_UNKNOWN ) imres.type = par.type;
  if ( VT_AllocImage( &imres ) != 1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }

  VT_ComputeNumberOfContoursPerCell( image, &imres, par.threshold );

  /*--- ecriture de l'image resultat ---*/
  if ( VT_WriteInrimage( &imres ) == -1 ) {
    VT_FreeImage( image );
    VT_FreeImage( &imres );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to write output image\n", 0);
  }

  VT_FreeImage( &imres );



  ncontours = 0;
  initStructure( &structure );

  for ( z=0; z<(int)image->dim.z; z++ ) {

    slice = VT_Compute2DIsoContours( image, z, par.threshold );
    if ( slice == NULL ) {
      printf( "bof\n" );
    }
    
    if ( 0 ) printf( "found %d contours in slice %d\n", slice->n, z  );

    ncontours += slice->n;

    if ( addSliceToStructure( slice, &structure ) != 1 ) {
      VT_ErrorParse("error when adding slice to structure\n", 0);
    }

  }
  
  printf( "found %d contours\n", ncontours );
  if ( 0 ) printStructure( &structure, 2 );

  /*
  for (j=0; j<listOfContours.n; j++ ) {
    printf( "%3d #=%d z=%f\n",j, listOfContours.theCts[j].n, 
  }
  */
  
  if ( par.names.ext[0] != '\0' ) {
    char name[DOUBLESTRINGLENGTH];
    FILE *f, *fopen();
    int startname, fd;
    int j;

    startname = strlen( par.names.ext )-1;
    for ( ; startname >= 0 && par.names.ext[startname] != '/' ; startname-- ) 
      ;

    sprintf( name, "%s.raw", par.names.ext );
    fd = creat( name, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
    
    sprintf( name, "%s.m", par.names.ext );
    f = fopen( name, "w" );

    fprintf( f, "\n\n\n" );
    fprintf( f, "%%\n" );
    fprintf( f, "%%" );
    for ( j=0; j<argc; j++ ) fprintf( f, " %s", argv[j] );
    fprintf( f, "\n" );
    fprintf( f, "%%\n" );
    fprintf( f, "\n\n\n" );

    fprintf( f, "echo off\n" );
    fprintf( f, "fid = fopen('%s.raw', 'r' );\n", &(par.names.ext[startname+1]) );
    
    fprintf( f, "figure;\n" );
    fprintf( f, "hold on;\n" );
    
    VT_2DDrawImage( image, par.slice, fd, f );

    for (j=0; j<structure.n; j++ ) {
      if ( structure.theSlices[j]->z > par.slice - 0.1 &&
	   structure.theSlices[j]->z < par.slice + 0.1 )
	for ( n=0; n<structure.theSlices[j]->n; n++ ) {
	  MAT_DrawContour2D( structure.theSlices[j]->theContours[n], fd, f, n );
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


      else if ( strcmp ( argv[i], "-th" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -th...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->threshold) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -th...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-slice" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -slice...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->slice) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -slice...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-matlab" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -matlab...\n", 0 );
	strncpy( par->names.ext, argv[i], STRINGLENGTH ); 
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
  if ( (o == 1) && (s == 1) && (r == 0) )  par->type = SCHAR;
  if ( (o == 1) && (s == 0) && (r == 0) ) par->type = UCHAR;
  if ( (o == 2) && (s == 0) && (r == 0) ) par->type = USHORT;
  if ( (o == 2) && (s == 1) && (r == 0) )  par->type = SSHORT;
  if ( (o == 4) && (s == 1) && (r == 0) )  par->type = SINT;
  if ( (o == 0) && (s == 0) && (r == 1) )  par->type = FLOAT;
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
  par->threshold = 100;
  par->slice = 0;
  par->type = TYPE_UNKNOWN;
}
