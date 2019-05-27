/*************************************************************************
 * addPlane.c -
 *
 * $Id: addPlane.c,v 1.0 2014/07/29 16:26:00 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2014/07/29
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>
#include <vt_symmetryPlane.h>


typedef struct local_par {

  vt_names names;
  double a;
  double b;
  double c;
  double d;
  int vert;
  unsigned char l;
} local_par;

/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );


static char *usage = "[image-in] [image-out] \n\
\t [-n %lf %lf %lf %lf] [-a %lf -b %lf -c %lf -d %lf] [-x %d]\n\
\t [-label | -l %d] [-v] [-D] [-help]";


static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -n %lf %lf %lf %lf : parametres du plan\n\
\t (plan d'equation n[0]*x+n[1]*y+n[2]*z+n[3]=0)\n\
\t -a %lf -b %lf -c %lf -d %lf : parametres du plan\n\
\t (plan d'equation a*x+b*y+c*z+d=0)\n\
\t -x %d : position du plan si plan vertical (coordonnee en X)\n\
\t -label %d : valeur des pixels du plan (default : 255)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];



int main( int argc, char *argv[] )
{
  /* char name[STRINGLENGTH]; */
  local_par par;
  vt_image *imageIn;
  /* vt_image imres; */
  double n[4];
  
  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );

  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );

  /*--- Lecture image d'entree ---*/
  imageIn = _VT_Inrimage( par.names.in );
  if ( imageIn == (vt_image*)NULL )
  {
	VT_ErrorParse("unable to read input image\n", 0);
  }
  
  /*--- Ajout du plan ---*/
  n[0]=par.a;
  n[1]=par.b;
  n[2]=par.c;
  n[3]=par.d;
  if (par.vert==1) {
	if( VT_AddPlaneToImage(imageIn, n, &(par.l)) != 1) {
	  VT_FreeImage(imageIn);
	  VT_ErrorParse("unable to add plane to image\n", 0);
    }
  }
  else
  {
    if( VT_AddVPlaneToImage(imageIn, par.a, &(par.l)) != 1) {
	  VT_FreeImage(imageIn);
	  VT_ErrorParse("unable to add plane to image\n", 0);
	}
  }



  /*--- Ecriture du resultat ---*/
  if ( VT_WriteInrimageWithName( imageIn, par.names.out ) == -1 ) {
	VT_FreeImage( imageIn );
  if(_VT_VERBOSE_)
	VT_ErrorParse("unable to write image\n", 0);
  }
  VT_FreeImage( imageIn );

  return( 0 );
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
      
      /*--- traitement eventuel de l'image d'entree ---*/



      /*--- Parametres d'input/output ---*/

      /* Parametres de calcul */

      else if ( strcmp ( argv[i], "-n" ) == 0 ) {
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -n...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->a) );
		if ( status <= 0 ) VT_ErrorParse( "parsing -n...\n", 0 );
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -n...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->b) );
		if ( status <= 0 ) VT_ErrorParse( "parsing -n...\n", 0 );
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -n...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->c) );
		if ( status <= 0 ) VT_ErrorParse( "parsing -n...\n", 0 );
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -n...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->d) );
		if ( status <= 0 ) VT_ErrorParse( "parsing -n...\n", 0 );
		par->vert=1;
      }
	  
      else if ( strcmp ( argv[i], "-a" ) == 0 || strcmp ( argv[i], "-x" ) == 0 ) {
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -a...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->a) );
		if ( status <= 0 ) VT_ErrorParse( "parsing -a...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-b" ) == 0 ) {
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -b...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->b) );
		if ( status <= 0 ) VT_ErrorParse( "parsing -b...\n", 0 );
		par->vert=1;
      }

      else if ( strcmp ( argv[i], "-c" ) == 0 ) {
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -c...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->c) );
		if ( status <= 0 ) VT_ErrorParse( "parsing -c...\n", 0 );
		par->vert=1;
      }

      else if ( strcmp ( argv[i], "-d" ) == 0 ) {
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -d...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->d) );
		if ( status <= 0 ) VT_ErrorParse( "parsing -d...\n", 0 );
		par->vert=1;
      }


      else if ( (strcmp ( argv[i], "-l" ) == 0) || 
				(strcmp ( argv[i], "-label" ) == 0) ) {
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -label...\n", 0 );
		int tmp;
		status = sscanf( argv[i],"%d",&tmp );
		if ( status <= 0 ) VT_ErrorParse( "parsing -label...\n", 0 );
		par->l = (unsigned char) tmp;
      }
/*
      else if ( strcmp ( argv[i], "-bin" ) == 0 ) {
        par->bin = 1;
      }
*/


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
  {
    strcpy( par->names.out, ">" );  /* standart output */
  }

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
  par->a=0.0;
  par->b=0.0;
  par->c=0.0;
  par->d=0.0;
  par->vert=0;
  par->l=(unsigned char) 255;
}


