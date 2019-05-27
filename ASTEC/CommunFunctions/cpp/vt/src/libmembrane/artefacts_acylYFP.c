/*************************************************************************
 * artefacts_acylYFP.c -
 *
 * $Id: artefacts_acylYFP.c,v 1.0 2014/04/03 15:53:00 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2014/04/03
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>
#include <time.h>


typedef enum {
  _DIFF_,
  _UNKNOWN_
} artefactMode;




typedef struct local_par {

  vt_names names;
  char Fname[STRINGLENGTH];
  /*   int dimension; */
  int writeImages;
  artefactMode mode;

} local_par;



/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );





static char *usage = "[image-in] [image-out] [-f %s] [-m diff]\n\
\t  [-wi] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -f %s : fichier de sortie log\n\
\t -m [diff] : mode de detection d'artefact\n\
\t -wi : ecrit toutes les images intermediaires\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  char name[DOUBLESTRINGLENGTH];
  local_par par;
  vt_image *imageIn;
  vt_image imout;
  vt_image *imageOut;
  int dim[3];
  bufferType t;
  int i,j,k;
  
  unsigned char ***arrayInU8;
  unsigned short int ***arrayInUSHORT;
  unsigned int ***arrayInUINT;
  float ***arrayInFLOAT;
  unsigned char ***arrayOutU8;
  unsigned short int ***arrayOutUSHORT;
  unsigned int ***arrayOutUINT;
  float ***arrayOutFLOAT;

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


  /*--- Image de sortie ---*/
  sprintf( name, "%s", par.names.out );
  t =  imageIn->type;
  dim[0]=imageIn->dim.x;
  dim[1]=imageIn->dim.y;
  dim[2]=imageIn->dim.z;

  VT_Image( &imout );
  VT_InitVImage( &imout, name, 1, dim[0], dim[1], dim[2], t );
		     
  imout.siz.x = imageIn->siz.x;
  imout.siz.y = imageIn->siz.y;
  imout.siz.z = imageIn->siz.z;

  if ( VT_AllocImage( &imout ) != 1 ) {
	VT_FreeImage( imageIn );
	MT_ErrorParse("unable to allocate output image\n", 0);
  }


  imageOut = &imout;


  /*--- calculs ---*/

  if(_VT_VERBOSE_)
    fprintf(stdout, "Calculs...\n");

  switch (t) {
  case UCHAR:
  case SCHAR:
	arrayInU8 = (unsigned char ***)imageIn->array;
	arrayOutU8 = (unsigned char ***)imageOut->array;
	switch ( par.mode ) {
	case _DIFF_:
	  for (k=0 ; k<dim[2]-1 ; k++) {
		for (j=0;j<dim[1];j++)
		for (i=0;i<dim[0];i++) {
		  int d = (int) arrayInU8[k+1][j][i] - (int)arrayInU8[k][j][i];
		  arrayOutU8[k][j][i]=(unsigned char) ((d>=0) ? d : -d);
		}
	  }
	  break;
	default:
	  VT_FreeImage( imageIn);
	  VT_FreeImage( imageOut);
	  MT_ErrorParse("Artefact mode not handled yet\n", 0);	  
	}
    break;
  case USHORT:
  case SSHORT:
  	arrayInUSHORT = (unsigned short int ***)imageIn->array;
	arrayOutUSHORT = (unsigned short int ***)imageOut->array;
	switch ( par.mode ) {
	case _DIFF_:
	  for (k=0 ; k<dim[2]-1 ; k++) {
		for (j=0;j<dim[1];j++)
		for (i=0;i<dim[0];i++) {
		  int d = (int)arrayInUSHORT[k+1][j][i]-(int)arrayInUSHORT[k][j][i];
		  arrayOutUSHORT[k][j][i]=(unsigned short int) ((d>=0) ? d : -d);
		}
	  }
	  break;
	default:
	  VT_FreeImage( imageIn);
	  VT_FreeImage( imageOut);
	  MT_ErrorParse("Artefact mode not handled yet\n", 0);	  
	}
	break;
  case UINT:
  case SINT:
	arrayInUINT = (unsigned int ***)imageIn->array;
	arrayOutUINT = (unsigned int ***)imageOut->array;
	switch ( par.mode ) {
	case _DIFF_:
	  for (k=0 ; k<dim[2]-1 ; k++) {
		for (j=0;j<dim[1];j++)
		for (i=0;i<dim[0];i++) {
		  int d = (int)arrayInUINT[k+1][j][i]-(int)arrayInUINT[k][j][i];
		  arrayOutUINT[k][j][i]=(unsigned int) ((d>=0) ? d : -d);
		}
	  }
	  break;
	default:
	  VT_FreeImage( imageIn);
	  VT_FreeImage( imageOut);
	  MT_ErrorParse("Artefact mode not handled yet\n", 0);	  
	}
    break;
  case FLOAT:
	arrayInFLOAT = (float ***)imageIn->array;
	arrayOutFLOAT = (float ***)imageOut->array;
	switch ( par.mode ) {
	case _DIFF_:
	  for (k=0 ; k<dim[2]-1 ; k++) {
		for (j=0;j<dim[1];j++)
		for (i=0;i<dim[0];i++) {
		  float d = arrayInFLOAT[k+1][j][i]-arrayInFLOAT[k][j][i];
		  arrayOutFLOAT[k][j][i]=((d>=0) ? d : -d);
		}
	  }
	  break;
	default:
	  VT_FreeImage( imageIn);
	  VT_FreeImage( imageOut);
	  MT_ErrorParse("Artefact mode not handled yet\n", 0);	  
	}  
    break;
  default:
	VT_FreeImage( imageIn);
	VT_FreeImage( imageOut);
	MT_ErrorParse("Image type not handled yet\n", 0);
  }


  if (par.writeImages == 1)
    if ( VT_WriteInrimage( imageOut ) == -1 ) {
      VT_FreeImage( imageIn );
	  VT_FreeImage( imageOut );
	  MT_ErrorParse("unable to write output image\n", 0);
    }
  
  /*--- liberations memoires ---*/
  VT_FreeImage( imageIn  );
  VT_FreeImage( imageOut );

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



	  /*
      else if ( strcmp ( argv[i], "-2D" ) == 0 ) {
        par->dimension = 2;
      }
      */ 



      /* Parametres de calcul */

      else if ( strcmp ( argv[i], "-f" ) == 0 ) {
		  		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -f...\n", 0 );
        strncpy( par->Fname, argv[i], STRINGLENGTH );
      }

      else if ( strcmp ( argv[i], "-m" ) == 0 ) {
		  		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -m...\n", 0 );
        if ( strcmp ( argv[i], "diff" ) == 0 ) par->mode = _DIFF_;
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
    strcpy( par->names.out, ">" );  /* standart output */
  }
  if (nb == 1)
  {
    strcpy( par->names.out, ">" );  /* standart output */
  }

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
  par->Fname[0] = '\0';
  par->mode = _UNKNOWN_;
}
