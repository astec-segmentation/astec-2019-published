/*************************************************************************
 * minimum.c -
 *
 * $Id: minimum.c,v 1.5 2000/08/16 16:31:56 greg Exp $
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

#include <vt_common.h>
#include <recbuffer.h>

typedef struct local_par {
  vt_names names;
  int type;
} local_par;




/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[image-in] [image-out]\n\
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
 $Revision: 1.5 $ $Date: 2000/08/16 16:31:56 $ $Author: greg $\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image, imres;
  
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
  VT_InitFromImage( &imres, image, par.names.out, FLOAT );

  if ( VT_AllocImage( &imres ) != 1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }
  
  {
    int dim[3];
    int cont[3] = { 20, 20, 20 };
    float coeffs[3] = { 1.0, 1.0, 1.0 };
    dim[0] = image->dim.x;
    dim[1] = image->dim.y;
    dim[2] = image->dim.z;
    GradientHessianGradient( image->buf, image->type,
			     imres.buf, imres.type,
			     dim, cont, coeffs, GAUSSIAN_DERICHE );
  }
  /*--- ecriture de l'image resultat ---*/
  if ( VT_WriteInrimage( &imres ) == -1 ) {
    VT_FreeImage( image );
    VT_FreeImage( &imres );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to write output image\n", 0);
  }
  
  /*--- liberations memoires ---*/
  VT_FreeImage( image );
  VT_FreeImage( &imres );
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
  if ( (o == 4) && (s == 1) && (r == 0) )  par->type = INT;
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
  par->type = TYPE_UNKNOWN;
}
