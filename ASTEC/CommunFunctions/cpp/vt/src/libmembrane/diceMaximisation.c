/*************************************************************************
 * diceMaximisation.c -
 *
 * $Id: diceMaximisation.c,v 1.0 2014/08/07 10:37:00 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2013/07/30
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vtmalloc.h>

#include <vt_common.h>
#include <mt_dice.h>
#include <time.h>


typedef struct local_par {

  vt_names names;
  int plane;
  int dimension;
  int vertical;
  double a;
  double b;
  double c;
  double d;
  int bin;
  double delta;
  int vox;
} local_par;



/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );





static char *usage = "[image-in] [filename-out] [-symmetry|s|n %f [%f %f %f]]\n\
\t [-delta %d] [-plane|p %s] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-ext' est absent, on prendra stdin\n\
\t si 'filename-out' est absent, on prendra stdout\n\
\t si les trois sont absents, on prendra stdin et stdout\n\
\t -symmetry %f : calcule le dice entre les volumes a gauche et a droite du \n\
\t                plan d'equation x = %f\n\
\t -n %f %f %f %f : calcule le dice entre les volumes de part et \n\
\t                d'autre du plan d'equation %f*x+%f*y+%f*z+%f = 0\n\
\t -voxel : equation du plan en coordonnees voxelliques\n\
\t -delta %d : fenetre de calcul du maximum de dice\n\
\t -plane %s : sauve dans %s l'image du plan\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *imageIn;
  int flag_3D = 1;
  size_t i;
  unsigned char *theU8;
  unsigned short int *theUSHORT;
  double n[4];
  int cst=0;
  double dice;
  int delta;
  int d;
  clock_t start, stop;
  double elapsed;
  double *dices;
  double dicemax=0;
  double dmax=0;

  start = clock();


  /*--- initialisation des parametres ---*/
  MT_InitParam( &par );

  /*--- lecture des parametres ---*/
  MT_Parse( argc, argv, &par );

  /*--- lecture de l'image d'entree ---*/
  imageIn = _VT_Inrimage( par.names.in );
  if ( imageIn == (vt_image*)NULL )
    MT_ErrorParse("unable to read input image\n", 0);

  
  if ( par.dimension == 2 || imageIn->dim.z == 1 )
    flag_3D = 0;

  
  /* if(_VT_VERBOSE_) */
  /*   fprintf(stdout, "Calculs...\n"); */

  /*--- calculs ---*/

  if(1 || par.bin == 1)
  {
	switch (imageIn->type) {
	case UCHAR:
	case SCHAR:
	  theU8=(unsigned char *)imageIn->buf;
	  for (i=0 ; i<(imageIn->dim.x*imageIn->dim.y*imageIn->dim.z) ; i++)
	    theU8[i] = (theU8[i]!=(unsigned char)0) ? (unsigned char) 1 : (unsigned char) 0;
	  break;
	case USHORT :
	case SSHORT :
	  theUSHORT=(unsigned short int *)imageIn->buf;
	  for (i=0 ; i<(imageIn->dim.x*imageIn->dim.y*imageIn->dim.z) ; i++)
	    theUSHORT[i] = (theUSHORT[i]!=(unsigned short int)0) ? (unsigned short int) 1 : (unsigned short int) 0;
	  break;
	default:
      VT_FreeImage( imageIn  );
      VT_Free( (void**)&imageIn );
      MT_ErrorParse("input and extra images must have the same size\n", 0);
    }
  }
  
  delta=(int)fabs(par.delta);
  dices = vtmalloc( (2*delta+1)*sizeof(double), "dices", argv[0] );

  if ( flag_3D ) {
    /*  3D case */
    if (par.vertical == 1) {
	  cst=(int)par.a;
      if(par.vox) {
          cst /= imageIn->siz.x;
      }
	  i=0;
	  for (d=-delta; d<=delta; d++) {
		if ( MT_ComputeSymmetryDice3DBis( imageIn, cst+d, &dice) != 1) {
		  VT_FreeImage( imageIn  );
		  VT_Free( (void**)&imageIn );
		  vtfree(dices);
		  dices=NULL;
		  MT_ErrorParse("unable to compute response\n", 0 );
	    }
	    if (dicemax<dice) {
		  dicemax=dice;
		  dmax=(int)cst+d;
	    }
		dices[i++] = dice;
	  }
	}
	else {
	  n[0]=par.a;
	  n[1]=par.b;
	  n[2]=par.c;
	  n[3]=par.d;
      if(par.vox == 0) {
          n[0]*=imageIn->siz.x;
          n[1]*=imageIn->siz.y;
          n[2]*=imageIn->siz.z;
          double tmp = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
          for (i=0; i<4; i++) n[i]/=tmp;
      }
	  double m[4];
	  m[0]=n[0]; m[1]=n[1]; m[2]=n[2];
	  i=0;
	  for (d=-delta; d<=delta; d++) {
		m[3]=n[3]+d;
		if( MT_ComputeSymmetryDiceOblic3DBis( imageIn, m, &dice) != 1) {
		  VT_FreeImage( imageIn  );
		  VT_Free( (void**)&imageIn );
		  vtfree(dices);
		  dices=NULL;
		  MT_ErrorParse("unable to compute response\n", 0 );
		}
	    if (dicemax<dice) {
		  dicemax=dice;
		  dmax=m[3];
	    }
		dices[i++] = dice;
	  }
	}
    if(par.vox == 0) {
        n[0]/=imageIn->siz.x;
        n[1]/=imageIn->siz.y;
        n[2]/=imageIn->siz.z;
        double tmp = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        for (i=0; i<4; i++) n[i]/=tmp;
        dmax/=tmp;
    }
	if (par.names.out[0] != 0) {
	  FILE *fichier = fopen (par.names.out, "w" );
	  if (fichier == NULL)
	  {
		perror (par.names.out);
	  }
	  if (par.vertical==1)
		fprintf(fichier, "diceold: %f \ndold: %d\n", dices[delta], (int)cst);
	  else
		fprintf(fichier, "diceold: %f \ndold: %f\n", dices[delta], n[3]);
	  fprintf(fichier, "dicemax: %f \ndmax: %f\n", dicemax, dmax);
	  if (par.vertical == 0)
		fprintf(fichier, "new n: %f %f %f %f\n", n[0], n[1], n[2], dmax);
	  fclose(fichier);
	}
	if (_VT_VERBOSE_ || par.names.out[0] == 0) {
	  if (par.vertical==1)
		fprintf(stdout, "diceold: %f \ndold: %d\n", dices[delta], (int)cst);
	  else
		fprintf(stdout, "diceold: %f \ndold: %f\n", dices[delta], n[3]);
	  fprintf(stdout, "dicemax: %f\ndmax: %f\n", dicemax, dmax);
	  if (par.vertical == 0)
		fprintf(stdout, "new n: %f %f %f %f\n", n[0], n[1], n[2], dmax);
	}
	/* fprintf(stdout, "dices = ["); */
	/* for( i=0; i<=2*delta; i++) */
	/* 	fprintf(stdout, " %f\t", dices[i]); */
	/* fprintf(stdout, "]\n"); */
	

    vtfree(dices);
    dices=NULL;

  }
  else {
    /*  2D case */
	  VT_FreeImage( imageIn  );
	  VT_Free( (void**)&imageIn );
      MT_ErrorParse("2D case not handled yet\n", 0 );
  }

  /*--- liberations memoires ---*/
  VT_FreeImage( imageIn  );
  VT_Free( (void**)&imageIn );

  stop = clock();
  elapsed = (double)(stop-start)/CLOCKS_PER_SEC;

  if(_VT_VERBOSE_)
    fprintf(stdout, "Elapsed time : \t%.1fs\n", elapsed);

  return( 0 );
}








static void MT_Parse( int argc,
                      char *argv[],
                      local_par *par )
{
  int i, nb;
  char text[STRINGLENGTH];
  int status;

  if ( VT_CopyName( program, argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);
  if ( argc == 1 ) MT_ErrorParse("\n", 0 );

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
        MT_ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 ) {
        _VT_VERBOSE_ = 1;
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 ) {
        _VT_DEBUG_ = 1;
      }

      /*--- traitement eventuel de l'image d'entree ---*/



      /*--- Parametres d'input/output ---*/

      else if ( strcmp ( argv[i], "-plane" ) == 0 || strcmp ( argv[i], "-p" ) == 0 ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -plane...\n", 0 );
		strncpy( par->names.ext, argv[i], STRINGLENGTH );      
		par->plane = 1;
	  }




      else if ( strcmp ( argv[i], "-2D" ) == 0 ) {
        par->dimension = 2;
      }



      /* Parametres de calcul */

      else if ( strcmp ( argv[i], "-bin" ) == 0 ) {
        par->bin = 1;
      }

      else if ( strcmp ( argv[i], "-voxel" ) == 0 ) {
        par->vox = 1;
      }


      else if ( strcmp ( argv[i], "-symmetry" ) == 0 ||
		strcmp(argv[i], "-s" ) == 0 || strcmp ( argv[i], "-n" ) == 0 ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -symmetry...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->a) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -symmetry...\n", 0 );
		par->vertical = 1;
		i += 1;
		if (i < argc) {
		  status = sscanf( argv[i],"%lf",&(par->b) );
		  if ( status <= 0 ) i -= 1;
		  else {
			par->vertical=0;
			i += 1;
			if ( i >= argc)    MT_ErrorParse( "parsing -symmetry...\n", 0 );
			status = sscanf( argv[i],"%lf",&(par->c) );
			if ( status <= 0 ) MT_ErrorParse( "parsing -symmetry...\n", 0 );
			i += 1;
			if ( i >= argc)    MT_ErrorParse( "parsing -symmetry...\n", 0 );
			status = sscanf( argv[i],"%lf",&(par->d) );
			if ( status <= 0 ) MT_ErrorParse( "parsing -symmetry...\n", 0 );
		  }
		}
		else i -= 1;
      }

      else if ( strcmp ( argv[i], "-delta" ) == 0) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -delta...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->delta) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -delta...\n", 0 );
	  }
      
      /*--- option inconnue ---*/
      else {
        sprintf(text,"unknown option %s\n",argv[i]);
        MT_ErrorParse(text, 0);
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
        MT_ErrorParse("too much file names when parsing\n", 0 );
    }
    i += 1;
  }

  /*--- s'il n'y a pas assez de noms ... ---*/
  if (nb == 0) {
    strcpy( par->names.in,  "<" );  /* standart input */
    strcpy( par->names.out, "" );  /* standart output */
  }
  if (nb == 1)
  {
    strcpy( par->names.out, "" );  /* standart output */
  }
  if (nb > 2)
    MT_ErrorParse("too much file names when parsing\n", 0 );
}







static void MT_ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"Erreur : %s",str);
  exit( 1 );
}








static void MT_InitParam( local_par *par )
{
  VT_Names( &(par->names) );

  par->vertical = -1;
  par->dimension = 3;
  par->bin = 0;
  par->plane = 0;
  par->a=0;
  par->b=0;
  par->c=0;
  par->d=0;
  par->delta=15;
  par->vox=0;
}
