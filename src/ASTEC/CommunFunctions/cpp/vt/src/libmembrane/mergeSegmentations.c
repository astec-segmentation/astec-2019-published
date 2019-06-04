/*************************************************************************
 * mergeSegmentations.c -
 *
 * $Id: mergeSegmentations.c,v 1.0 2014/04/02 11:52:00 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2014/04/02
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>
#include <mt_mergeSegmentations.h>
#include <time.h>


typedef struct local_par {

  vt_names names;
  int dilatation;
  /*   int dimension; */
  int writeImages;

} local_par;



/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );





static char *usage = "[image-in] [image-ext] [image-out]\n\
\t [-dilatation | -d %d] [-wi] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-ext' est absent, on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les trois sont absents, on prendra stdin et stdout\n\
\t -[dilatation | d] %d : parametre pour l'element structurant de la dilatation\n\
\t -wi : ecrit toutes les images intermediaires\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  /* char filename[256]; */
  /* char *temp; */
  char name[STRINGLENGTH], tmp[STRINGLENGTH];
  local_par par;
  vt_image *imageIn;
  vt_image *imageExt;
  vt_image imout;
  vt_image *imageOut;
  vt_image *imageMask;
  vt_image mask;
  int dim[3];
  /* int flag_3D = 1; */
  bufferType type = USHORT;
  int i,j;
  bufferType t;

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
  imageExt = _VT_Inrimage( par.names.ext );
  if ( imageExt == (vt_image*)NULL )
  {
    VT_FreeImage( imageIn  );
    VT_Free( (void**)&imageIn );
    MT_ErrorParse("unable to read ext image\n", 0);
  }

  /*--- Image de sortie ---*/
  sprintf( name, "%s", par.names.out );
  type =  USHORT;
  dim[0]=imageIn->dim.x;
  dim[1]=imageIn->dim.y;
  dim[2]=imageIn->dim.z;

  VT_Image( &imout );
  VT_InitVImage( &imout, name, 1, dim[0], dim[1], dim[2], type );
  /* VT_InitVImage( &imres, name, 1,  */
  /*      imtemplate->dim.x, imtemplate->dim.y, imtemplate->dim.z,  */
  /*      imtemplate->type ); */
		     
  imout.siz.x = imageIn->siz.x;
  imout.siz.y = imageIn->siz.y;
  imout.siz.z = imageIn->siz.z;

  if ( VT_AllocImage( &imout ) != 1 ) {
	VT_FreeImage( imageIn );
    VT_FreeImage( imageExt );
	MT_ErrorParse("unable to allocate output image\n", 0);
  }


  imageOut = &imout;


  /*--- calculs ---*/

  if(_VT_VERBOSE_)
    fprintf(stdout, "Calculs...\n");

  /*--- Image Mask ---*/
  if(_VT_VERBOSE_)
    fprintf(stdout, "Img Mask...\n");
  /*  allocation  */
  sprintf(name, "%s", imageIn->name);
  for (i=0 ; name[i]!='\0'; i++) {}
  while(i>0 && name[i--]!='.') {}
  if(i==0) 
  {
	fprintf(stderr, "Error when renaming mask image\n");
	return(0);
  }
  for(j=0;j<=i;j++) 
	tmp[j]=name[j];
  tmp[j++]='_';tmp[j++]='m';tmp[j++]='a';tmp[j++]='s';tmp[j++]='k';
  while (name[++i] != '\0')
	tmp[j++]=name[i];
  tmp[j]='\0';
   
  dim[0]=imageIn->dim.x;
  dim[1]=imageIn->dim.y;
  dim[2]=imageIn->dim.z;  

  t=UCHAR;

  VT_Image( &mask );
  VT_InitVImage( &mask, tmp, 1, dim[0], dim[1], dim[2], t );

  mask.siz.x = imageIn->siz.x;
  mask.siz.y = imageIn->siz.y;
  mask.siz.z = imageIn->siz.z;

  if ( VT_AllocImage( &mask ) != 1 ) {
    fprintf(stderr, "Error while allocating mask image\n");
	return(0);
  }

  imageMask=&mask;

  if( MT_BackgroundMask( *imageIn, imageMask, (int) 1) != 1) {
	VT_FreeImage( imageIn );
	VT_FreeImage( imageExt );
	VT_FreeImage( imageOut );
	VT_FreeImage( imageMask );
  }

  if(_VT_VERBOSE_)
    fprintf(stdout, "Writing Mask...\n");

  if (par.writeImages == 1)
    if ( VT_WriteInrimage( imageMask ) == -1 ) {
      VT_FreeImage( imageMask );
      VT_FreeImage( imageIn );
	  VT_FreeImage( imageExt );
	  VT_FreeImage( imageOut );
	  MT_ErrorParse("unable to write mask image\n", 0);
    }
  
  /*--- Listage des regions d'imageExt dans le fond de imageIn ---*/

  /* TODO */

  /*--- liberations memoires ---*/
  VT_FreeImage( imageExt );
  VT_FreeImage( imageIn  );
  VT_FreeImage( imageMask );
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

      else if ( strcmp ( argv[i], "-dilatation" ) == 0 ||
          strcmp ( argv[i], "-d" ) == 0 ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -dilatation...\n", 0 );
		status = sscanf( argv[i],"%d",&(par->dilatation) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -dilatation...\n", 0 );
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

  par->writeImages = 0;
  par->dilatation = 0;
  /* par->dimension = 3; */

}
