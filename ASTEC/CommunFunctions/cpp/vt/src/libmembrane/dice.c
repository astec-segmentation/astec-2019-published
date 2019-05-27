/*************************************************************************
 * dice.c -
 *
 * $Id: dice.c,v 1.0 2013/07/30 11:05:00 gael Exp $
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

#include <vt_common.h>
#include <mt_dice.h>
#include <time.h>


typedef struct local_par {

  vt_names names;
  int seed;
  int dimension;
  int writeImages;
  int sym;
  double n[4];
  int obl;
  int bin;
} local_par;



/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );





static char *usage = "[image-in] [image-ext] [filename-out] [-symmetry|s|n %f [%f %f %f]]\n\
\t [-graines|-seed|-intersection] [-bin] [-wi] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-ext' est absent, on prendra stdin\n\
\t si 'filename-out' est absent, on prendra stdout\n\
\t si les trois sont absents, on prendra stdin et stdout\n\
\t -[graines|seed] : ecrit le tableau de confusion a la place des Dice\n\
\t -intersection : ecrit le tableau des intersections label/label\n\
\t -symmetry %d : calcule le dice entre les volumes a gauche et a droite du \n\
\t                plan d'equation x = %d\n\
\t -n %f %f %f %f : calcule le dice entre les volumes de part et d'autre du \n\
\t                plan d'equation n[0]*x+n[1]*y+n[2]*z+n[3] = 0\n\
\t -bin : calcul du dice sur les images binaires (sans les labels)\n\
\t -wi : ecrit toutes les images intermediaires\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  char filename[256];
  char *temp;
  local_par par;
  vt_image *imageIn;
  vt_image *imageExt;
  int flag_3D = 1;
  size_t i;
  unsigned char *theU8;
  unsigned short int *theUSHORT;
  
  clock_t start, stop;
  double elapsed;

  start = clock();


  /*--- initialisation des parametres ---*/
  MT_InitParam( &par );

  /*--- lecture des parametres ---*/
  MT_Parse( argc, argv, &par );

  /*--- lecture de l'image d'entree ---*/
  imageIn = _VT_Inrimage( par.names.in );
  if ( imageIn == (vt_image*)NULL )
    MT_ErrorParse("unable to read input image\n", 0);

  /*--- lecture de l'image extra ---*/
  if(par.sym<=0 && par.obl == 0) {
	imageExt = _VT_Inrimage( par.names.ext );
	if ( imageExt == (vt_image*)NULL )
	{
      VT_FreeImage( imageIn  );
      VT_Free( (void**)&imageIn );
      MT_ErrorParse("unable to read input image\n", 0);
	}
  }
  else {
	sprintf(par.names.out, "%s", par.names.ext);
  }
  
  temp = strcpy(filename, par.names.out);
  if(temp == NULL)
  {
	if(par.sym<0 && par.obl == 0) {
      VT_FreeImage( imageExt );
      VT_Free( (void**)&imageExt );
	}
    VT_FreeImage( imageIn  );
    VT_Free( (void**)&imageIn );
    MT_ErrorParse("problem while using strcpy\n", 0);
  }

  if (par.sym < 0 && par.obl == 0 && (imageIn->dim.x != imageExt->dim.x || imageIn->dim.y != imageExt->dim.y ||
      imageIn->dim.z != imageExt->dim.z))
  {
    VT_FreeImage( imageExt );
    VT_Free( (void**)&imageExt );
    VT_FreeImage( imageIn  );
    VT_Free( (void**)&imageIn );
    MT_ErrorParse("input and extra images must have the same size\n", 0);
  }
  
  if ( par.dimension == 2 || imageIn->dim.z == 1 )
    flag_3D = 0;

  
  if(_VT_VERBOSE_)
    fprintf(stdout, "Calculs...\n");

  /*--- calculs ---*/

  if(par.bin == 1)
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
	  if(par.sym<0 && par.obl == 0) {
		VT_FreeImage( imageExt );
		VT_Free( (void**)&imageExt );
	  }
      VT_FreeImage( imageIn  );
      VT_Free( (void**)&imageIn );
      MT_ErrorParse("input and extra images must have the same size\n", 0);
	}

	if (par.sym<0 && par.obl == 0)
	  switch (imageExt->type) {
	  case UCHAR:
	  case SCHAR:
		theU8=(unsigned char *)imageExt->buf;
		for (i=0 ; i<(imageExt->dim.x*imageExt->dim.y*imageExt->dim.z) ; i++)
	      theU8[i] = (theU8[i]!=(unsigned char)0) ? (unsigned char) 1 : (unsigned char) 0;
		break;
	  case USHORT :
	  case SSHORT :
		theUSHORT=(unsigned short int *)imageExt->buf;
		for (i=0 ; i<(imageExt->dim.x*imageExt->dim.y*imageExt->dim.z) ; i++)
	      theUSHORT[i] = (theUSHORT[i]!=(unsigned short int)0) ? (unsigned short int) 1 : (unsigned short int) 0;
		break;
	  default:
		VT_FreeImage( imageExt );
		VT_Free( (void**)&imageExt );
		VT_FreeImage( imageIn  );
		VT_Free( (void**)&imageIn );
		MT_ErrorParse("input and extra images must have the same size\n", 0);
	  }
  }

  if ( flag_3D ) {
    /*  3D case */
    if (par.sym<0 && par.obl == 0) {
	  if (par.seed == 0)  {
        if ( MT_ComputeDice3D( *imageIn, *imageExt, filename) != 1 ) {
		  VT_FreeImage( imageExt );
		  VT_FreeImage( imageIn  );
		  VT_Free( (void**)&imageExt );
		  VT_Free( (void**)&imageIn );
		  MT_ErrorParse("unable to compute response\n", 0 );
		}
	  }
	  else {
        if (par.seed == 1)  {
          if ( MT_ComputeConfusionTab3D( *imageIn, *imageExt, filename) != 1 ) {
            VT_FreeImage( imageExt );
            VT_FreeImage( imageIn  );
            VT_Free( (void**)&imageExt );
            VT_Free( (void**)&imageIn );
            MT_ErrorParse("unable to compute response\n", 0 );
          }
        }
        else {
          if ( MT_ComputeIntersection3D( *imageIn, *imageExt, filename) != 1 ) {
            VT_FreeImage( imageExt );
            VT_FreeImage( imageIn  );
            VT_Free( (void**)&imageExt );
            VT_Free( (void**)&imageIn );
            MT_ErrorParse("unable to compute response\n", 0 );
          }
        }
      }
	}
	else {
	  if (par.obl == 0) {
		if(  MT_ComputeSymmetryDice3D( imageIn, par.sym, filename) != 1) {
		VT_FreeImage( imageIn  );
		VT_Free( (void**)&imageIn );
        MT_ErrorParse("unable to compute response\n", 0 );
	    }
	  }
	  else {
		if (MT_ComputeSymmetryDiceOblic3D( imageIn, par.n, filename) != 1) {
		  VT_FreeImage( imageIn  );
		  VT_Free( (void**)&imageIn );
          MT_ErrorParse("unable to compute response\n", 0 );
	    }
	  }
	}
  }
  else {
    /*  2D case */
    if ( MT_ComputeDice2D( *imageIn, *imageExt, filename) != 1 ) {
	  if (par.obl == 0 && par.sym < 0) {
		VT_FreeImage( imageExt );
		VT_Free( (void**)&imageExt );
	  }
	  VT_FreeImage( imageIn  );
	  VT_Free( (void**)&imageIn );
      MT_ErrorParse("unable to compute response\n", 0 );
    }

  }

  /*--- liberations memoires ---*/
  if (par.obl == 0 && par.sym < 0 ) {
	VT_FreeImage( imageExt );
    VT_Free( (void**)&imageExt );
  }
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

      else if ( strcmp ( argv[i], "-wi" ) == 0 ) {
        par->writeImages = 1;
      }




      else if ( strcmp ( argv[i], "-2D" ) == 0 ) {
        par->dimension = 2;
      }



      /* Parametres de calcul */

      else if ( strcmp ( argv[i], "-bin" ) == 0 ) {
        par->bin = 1;
      }



      else if ( strcmp ( argv[i], "-graines" ) == 0 ||
          strcmp ( argv[i], "-seed" ) == 0 ||strcmp ( argv[i], "-seeds" ) == 0
          || strcmp ( argv[i], "-graines" ) == 0) {
        par->seed = 1;
      }

      else if ( strcmp ( argv[i], "-intersection" ) == 0 ) {
        par->seed = 2;
      }

      else if ( strcmp ( argv[i], "-symmetry" ) == 0 ||
		(strcmp ( argv[i], "-s" ) == 0 ) ||(strcmp ( argv[i], "-n" ) == 0 )) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -symmetry...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->n[0]) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -symmetry...\n", 0 );
		i += 1;
		if ( i >= argc) {
		  i -= 1;
		  par->sym = (int) par->n[0];
		}
		else {
		  status = sscanf( argv[i],"%lf",&(par->n[1]) );
		  if ( status <= 0 ) {
			i -= 1;
			par->sym = (int) par->n[0];
		  }
		  else {
			i += 1;
			if ( i >= argc)    MT_ErrorParse( "parsing -n...\n", 0 );
			status = sscanf( argv[i],"%lf",&(par->n[2]) );
			if ( status <= 0 ) MT_ErrorParse( "parsing -n...\n", 0 );
			i += 1;
			if ( i >= argc)    MT_ErrorParse( "parsing -n...\n", 0 );
			status = sscanf( argv[i],"%lf",&(par->n[3]) );
			if ( status <= 0 ) MT_ErrorParse( "parsing -n...\n", 0 );
			par->obl=1;
		  }
		}
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
        strncpy( par->names.ext, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else if ( nb == 2 ) {
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
    strcpy( par->names.ext, "<" );  /* standart extra */
    strcpy( par->names.out, ">" );  /* standart output */
  }
  if (nb == 1)
  {
    strcpy( par->names.ext, "<" );  /* standart extra */
    strcpy( par->names.out, ">" );  /* standart output */
  }
  if (nb == 2)
    strcpy( par->names.out, ">" );  /* standart output */
  if (nb > 3 )
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

  par->writeImages = 0;
  par->seed = 0; /*  0: dice ; 1: seed/label confusion table ; 2: intersection table */
  par->dimension = 3;
  par->sym = -1;
  par->obl=0;
  par->bin = 0;
}
