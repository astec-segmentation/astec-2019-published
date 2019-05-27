/*************************************************************************
 * mergeLabels.c -
 *
 * $Id: mergeLabels.c,v 1.0 2014/05/28 13:59:51 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2014/05/28
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>

typedef struct local_par {
  vt_names names;
  char labels[STRINGLENGTH];
  int bckgrd;
} local_par;



/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );

static int _verbose_ = 0;





static char *usage = "[image-in] [image-ext] [image-out]\n\
\t [-labels | -l %s] [-bckgrd | -b %d] [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-ext' est absent, on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les trois sont absents, on prendra stdin et stdout\n\
\t -bckgrd | -b %d : valeur du background dans image-in\n\
\t -labels | -l %d : labels de image-ext a inserer dans image-in\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];

static int 	indexLabel(unsigned int valeur,unsigned int *labelsExt, int n);






int main( int argc, char *argv[] )
{

  local_par par;
  vt_image *imageIn;
  vt_image *imageExt;
  vt_image imres;
  char name[DOUBLESTRINGLENGTH];
  int z,y,x;
  int i;
  
  
  int valeur;

  unsigned char ***inU8=NULL, ***extU8=NULL;
  unsigned short int ***inU16=NULL, ***extU16=NULL;
  unsigned short int ***res;
  FILE* fi;
  char mot[256];
  int n_lab;
  
  int status; 
  int ind;
  
  unsigned int *labelsExt;




  /*--- initialisation des parametres ---*/
  MT_InitParam( &par );

  /*--- lecture des parametres ---*/
  MT_Parse( argc, argv, &par );

  /*--- lecture de l'image d'entree ---*/
  imageIn = _VT_Inrimage( par.names.in );
  if ( imageIn == (vt_image*)NULL )
    MT_ErrorParse("unable to read input image-in\n", 0);

  imageExt = _VT_Inrimage( par.names.ext );
  if ( imageExt == (vt_image*)NULL ) {
    VT_FreeImage( imageIn );
    MT_ErrorParse("unable to read input image-ext\n", 0);
  }
  
  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.names.inv == 1 )  VT_InverseImage( imageIn );
  if ( par.names.inv == 1 )  VT_InverseImage( imageExt );
  if ( par.names.swap == 1 ) VT_SwapImage( imageIn );
  if ( par.names.swap == 1 ) VT_SwapImage( imageExt );



  if (imageIn->dim.x != imageExt->dim.x || imageIn->dim.y != imageExt->dim.y ||
	  imageIn->dim.z != imageExt->dim.z )
  {
    VT_FreeImage( imageIn );
    VT_FreeImage( imageExt );
    MT_ErrorParse("Input image sizes do not match\n", 0);
  }
	
  /* Alloc image sortie */

  sprintf( name, "%s", par.names.out );
  VT_InitImage( &imres, name, imageIn->dim.x, imageIn->dim.y,
              imageIn->dim.z, (int)USHORT );
  if ( VT_AllocImage( &(imres) ) != 1 ) {
    VT_FreeImage( imageIn );
    VT_FreeImage( imageExt );
    MT_ErrorParse("problem while allocating imres\n", 0 );
  }

  /*--- Lecture labels ---*/

  fi=fopen(par.labels,"r");
  labelsExt=malloc(256*sizeof(unsigned int));
  i=0;
  while (fgets(mot, 256, fi) != NULL) {
	if ( i && (i % 256) ) {
	  labelsExt=realloc(labelsExt,(i+256)*sizeof(unsigned int));
	  if (labelsExt==NULL)
	  { 
		VT_FreeImage( imageIn );
		VT_FreeImage( imageExt );
		VT_FreeImage( &imres );
		MT_ErrorParse("problem while allocating table for labels\n", 0 );
	  }
	}
	status = sscanf( mot,"%u",&(labelsExt[i++]) );
	if ( status <= 0 ) {
	  fprintf(stderr, "Erreur de lecture (iteration %d), mot = %s\n", i, mot );
	  MT_ErrorParse( "label file corrupted...\n", 0 );
	}
  }
  n_lab=i;
  
  if(_verbose_) {
	fprintf(stdout, "labelsExt = { ");
	for (i=0; i<n_lab ; i++)
	  fprintf(stdout, "%u ", labelsExt[i]);
	fprintf(stdout, "}\n");
  }

  

  /*--- calculs ---*/





  switch (imageIn->type) {
    case SCHAR:
    case UCHAR:
      inU8=(unsigned char ***)imageIn->array;
      break;
    case SSHORT:
    case USHORT:
      inU16=(unsigned short int ***)imageIn->array;
      break;
    case TYPE_UNKNOWN:
    default:
	  VT_FreeImage( imageIn );
	  VT_FreeImage( imageExt );
	  VT_FreeImage( &imres );
	  free(labelsExt);
	  labelsExt=NULL;
      VT_Error("image type unknown or not supported for this program",program);
      return( 0 );
  }
  
  switch (imageExt->type) {
    case SCHAR:
    case UCHAR:
      extU8=(unsigned char ***)imageExt->array;
      break;
    case SSHORT:
    case USHORT:
      extU16=(unsigned short int ***)imageExt->array;
      break;
    case TYPE_UNKNOWN:
    default:
	  VT_FreeImage( imageIn );
	  VT_FreeImage( imageExt );
	  VT_FreeImage( &imres );
	  free(labelsExt);
	  labelsExt=NULL;
      VT_Error("image type unknown or not supported for this program",program);
      return( 0 );
  }
  
  res=(unsigned short int ***)imres.array;


  /* Calcul du label maximum dans imageIn
   */
  unsigned short int labInMax=0;
  for(z=0;z<(int)imageIn->dim.z;z++)
  for(y=0;y<(int)imageIn->dim.y;y++)
  for(x=0;x<(int)imageIn->dim.x;x++)
  {
    switch (imageIn->type) {
    case SCHAR:
    case UCHAR:
	  if ( labInMax < (unsigned short int) inU8[z][y][x] )
		labInMax = (unsigned short int) inU8[z][y][x]; 
      break;
	case SSHORT:
	case USHORT:
	  if ( labInMax < (unsigned short int) inU16[z][y][x] )
		labInMax = (unsigned short int) inU16[z][y][x]; 
	  break;
	default:
	  VT_FreeImage( imageIn );
	  VT_FreeImage( imageExt );
	  VT_FreeImage( &imres );
	  free(labelsExt);
	  labelsExt=NULL;
      VT_Error("image type unknown or not supported for this program",program);
	  return(0);
	}
  }

  for(z=0;z<(int)imageIn->dim.z;z++)
  for(y=0;y<(int)imageIn->dim.y;y++)
  for(x=0;x<(int)imageIn->dim.x;x++)
  {
    switch (imageIn->type) {
    case SCHAR:
    case UCHAR:
	  if (inU8[z][y][x] != (unsigned char) par.bckgrd )
	  {
		res[z][y][x]=(unsigned short int) inU8[z][y][x];
		continue;
	  }
      break;
	case SSHORT:
	case USHORT:
	  if (inU16[z][y][x] != (unsigned short int) par.bckgrd )
	  {
		res[z][y][x]=(unsigned short int) inU16[z][y][x];
		continue;
	  }	
	  break;
	default:
	  VT_FreeImage( imageIn );
	  VT_FreeImage( imageExt );
	  VT_FreeImage( &imres );
	  free(labelsExt);
	  labelsExt=NULL;
      VT_Error("image type unknown or not supported for this program",program);
	  return(0);
	}
	switch (imageExt->type) {
	case SCHAR:
	case UCHAR:
	  valeur=(unsigned int) extU8[z][y][x];
	  break;
	case SSHORT:
	case USHORT:
	  valeur=(unsigned int) extU16[z][y][x];	
	  break;
	default:
	  VT_FreeImage( imageIn );
	  VT_FreeImage( imageExt );
	  VT_FreeImage( &imres );
	  free(labelsExt);
	  labelsExt=NULL;
      VT_Error("image type unknown or not supported for this program",program);
	  return(0);
	}
	ind=indexLabel(valeur,labelsExt, n_lab);
	if (ind<0)
	{
	  res[z][y][x]=(unsigned short int) par.bckgrd;
	  continue;
	}
	res[z][y][x]=(unsigned short int) labInMax+ind+1;
  }


  VT_WriteInrimage( &(imres) );

  if (_verbose_)
    fprintf(stdout, "Labels ajoutes >= %u\n", labInMax+1);

  /*--- liberations memoires ---*/
  VT_FreeImage( imageIn );
  VT_FreeImage( imageExt );
  VT_Free( (void**)&imageIn );
  VT_Free( (void**)&imageExt );
  VT_FreeImage( &imres );
  free(labelsExt);
  labelsExt=NULL;
  
  
  
  return( 0 );
}








static void MT_Parse( int argc,
                      char *argv[],
                      local_par *par )
{
  int i, nb, status;
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
        _verbose_++;
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


      /* Parametres de calcul */

      else if ( strcmp ( argv[i], "-labels" ) == 0 || 
				strcmp ( argv[i], "-l" ) == 0 ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -labels...\n", 0 );
        strncpy( par->labels, argv[i], STRINGLENGTH );
      }

      else if ( strcmp ( argv[i], "-bckgrd" ) == 0 || 
				strcmp ( argv[i], "-b" ) == 0 ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -bckgrd...\n", 0 );
		status = sscanf( argv[i],"%d",&(par->bckgrd) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -bckgrd...\n", 0 );
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
  par->bckgrd = 1;
  strcpy(par->labels, "" );
}


static int 	indexLabel(unsigned int valeur,unsigned int *labelsExt, int n)
{
  int res=n-1;
  while (res>=0 && labelsExt[res] != valeur)
    res--;
  
  return(res);
}
