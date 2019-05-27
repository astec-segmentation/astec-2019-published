/*************************************************************************
 * extractLabels.c -
 *
 * $Id: extractLabels.c,v 1.0 2014/05/28 13:59:51 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2014/10/22
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>
#define LABELSLENGTH 1000

typedef struct local_par {
  vt_names names;
  int labels[LABELSLENGTH];
  int nlabels;
  int inv;
  int set;
} local_par;



/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );

static int _verbose_ = 0;





static char *usage = "[image-in] [image-mask] [image-out] [-set]\n\
\t [-labels | -l %d [%d [...]]] [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-mask' est absent, on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les trois sont absents, on prendra stdin et stdout\n\
\t -set : attribue aux labels la valeur correspondant d'image-mask (255 si plusieurs valeurs)\n\
\t -labels | -l %d [...] : labels a extraire dans image-in (si pas image-mask)\n\
\t -inv : inverse 'image-mask' ou extrait les labels non cites dans -labels\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];


int main( int argc, char *argv[] )
{

  local_par par;
  vt_image *imageIn;
  vt_image *imageExt=NULL;
  vt_image imres;
  int *LAB;
  int max;
  size_t i;
  int j;
  int nLAB;
  

  unsigned char *inU8=NULL, *extU8=NULL;
  unsigned short int *inU16=NULL;
  unsigned char *resU8=NULL;
  unsigned short int *resU16=NULL;





  /*--- initialisation des parametres ---*/
  MT_InitParam( &par );

  /*--- lecture des parametres ---*/
  MT_Parse( argc, argv, &par );
  if (par.nlabels>0 && par.set > 0)
      MT_ErrorParse("conflicting options -set and -labels\n", 0);
  /*--- lecture de l'image d'entree ---*/
  imageIn = _VT_Inrimage( par.names.in );
  if ( imageIn == (vt_image*)NULL )
    MT_ErrorParse("unable to read input image-in\n", 0);

  if(par.nlabels == 0) {
    imageExt = _VT_Inrimage( par.names.ext );
    if ( imageExt == (vt_image*)NULL ) {
      VT_FreeImage( imageIn );
      MT_ErrorParse("unable to read input image-ext\n", 0);
    }
  }

  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.names.swap == 1 ) VT_SwapImage( imageIn );
  if ( par.names.swap == 1 && par.nlabels==0) VT_SwapImage( imageExt );


  if (par.nlabels == 0 && (imageIn->dim.x != imageExt->dim.x || imageIn->dim.y != imageExt->dim.y ||
      imageIn->dim.z != imageExt->dim.z ))
  {
    VT_FreeImage( imageIn );
    VT_FreeImage( imageExt );
    MT_ErrorParse("Input image sizes do not match\n", 0);
  }
	
  /* Alloc image sortie */

  if (par.set==0) VT_InitImage( &imres, par.names.out, imageIn->dim.x, imageIn->dim.y,
              imageIn->dim.z, imageIn->type );
  if (par.set==1) VT_InitImage( &imres, par.names.out, imageIn->dim.x, imageIn->dim.y,
              imageIn->dim.z, imageExt->type );
  if ( VT_AllocImage( &(imres) ) != 1 ) {
    VT_FreeImage( imageIn );
    if (par.nlabels==0) VT_FreeImage( imageExt );
    MT_ErrorParse("problem while allocating imres\n", 0 );
  }

  imres.siz.x=imageIn->siz.x;
  imres.siz.y=imageIn->siz.y;
  imres.siz.z=imageIn->siz.z;


  /* Switch for buffers */

  switch (imageIn->type) {
  case UCHAR:
  case SCHAR:
      inU8=(unsigned char*) imageIn->buf;
      break;
  case USHORT:
  case SSHORT:
      inU16=(unsigned short int*) imageIn->buf;
      break;
  default:
      VT_FreeImage( imageIn );
      if (par.nlabels==0) VT_FreeImage( imageExt );
      VT_FreeImage( &imres );
      MT_ErrorParse("Such image type not handled yet\n", 0);
  }

  if (par.nlabels==0){
    switch (imageExt->type) {
    case UCHAR:
    case SCHAR:
      extU8=(unsigned char*) imageExt->buf;
      break;
    default:
      VT_FreeImage( imageIn );
      VT_FreeImage( imageExt );
      VT_FreeImage( &imres );
      MT_ErrorParse("Such image type not handled yets (mask image)\n", 0);
    }
  }

  switch (imres.type) {
  case UCHAR:
  case SCHAR:
      resU8=(unsigned char*) imres.buf;
      break;
  case USHORT:
  case SSHORT:
      resU16=(unsigned short int*) imres.buf;
      break;
  default:
      VT_FreeImage( imageIn );
      if (par.nlabels==0) VT_FreeImage( imageExt );
      VT_FreeImage( &imres );
      MT_ErrorParse("Such image type not handled yet\n", 0);
  }

  /* Compute the max element of imageIn */

  max=0;
  switch(imageIn->type){
  case UCHAR:
  case SCHAR:
      max=255;
      break;
  case SSHORT:
  case USHORT:
      for (i=0 ; i<imageIn->dim.x*imageIn->dim.y*imageIn->dim.z ; i++)
        if (inU16[i]>max) max=inU16[i];
      break;
  default:
      VT_FreeImage( imageIn );
      if(par.nlabels == 0 ) VT_FreeImage( imageExt );
      VT_FreeImage( &imres );
      MT_ErrorParse("Such image type not handled yets (mask image)\n", 0);
  }
  nLAB=max+1;
  LAB=malloc(nLAB*sizeof(int));

  for (j=0;j<nLAB;j++)
      if(par.inv==0) LAB[j]=0;
      else LAB[j]=1;

  /*--- Calcul des labels a conserver ---*/

  if(par.nlabels==0) {
    for (i=0; i<imageIn->dim.x*imageIn->dim.y*imageIn->dim.z ; i++)
    {
        if(extU8[i]==0) continue;
        switch (imageIn->type){
        case UCHAR:
        case SCHAR:
            if (par.set==0) {
              if (par.inv==0) LAB[(int)inU8[i]]=1;
              else LAB[(int)inU8[i]]=0;
            }
            else {
              if((int)inU8[i]==0) continue;
              if (LAB[(int)inU8[i]]==(int)extU8[i] || LAB[(int)inU8[i]]==0)
                LAB[(int)inU8[i]]=(int)extU8[i];
              else
                LAB[(int)inU8[i]]=(int)255;
            }
            break;
        case USHORT:
        case SSHORT:
            if (par.set==0) {
              if (par.inv==0) LAB[(int)inU16[i]]=1;
              else LAB[(int)inU16[i]]=0;
            }
            else {
                if((int)inU16[i]==0) continue;
                if (LAB[(int)inU16[i]]==(int)extU8[i] || LAB[(int)inU16[i]]==0)
                  LAB[(int)inU16[i]]=(int)extU8[i];
                else
                  LAB[(int)inU16[i]]=(int)255;
            }
            break;
        default:
            free(LAB);
            VT_FreeImage( imageIn );
            VT_FreeImage( imageExt );
            VT_FreeImage( &imres );
            MT_ErrorParse("Such image type not handled yets (mask image)\n", 0);
        }
    }
    VT_FreeImage( imageExt );
  }
  else {
      for(j=0; j<par.nlabels; j++)
          if( par.inv == 0) LAB[par.labels[j]]=1;
          else LAB[par.labels[j]]=0;
  }

  /*--- Calcul de l'image resultat ---*/

  switch (imageIn->type) {
    case SCHAR:
    case UCHAR:
      for (i=0 ; i<imageIn->dim.x*imageIn->dim.y*imageIn->dim.z; i++)
          if(par.set==0)
          {
            resU8[i] = (LAB[(int)inU8[i]]==0)? 0 : inU8[i];
          }
          else
          {
            switch (imres.type) {
            case UCHAR:
            case SCHAR:
                resU8[i] = (unsigned char) LAB[(int)inU8[i]];
                break;
            case USHORT:
            case SSHORT:
                resU16[i] = (unsigned short int) LAB[(int)inU8[i]];
                break;
            default:
                free(LAB);
                VT_FreeImage( imageIn );
                VT_FreeImage( &imres );
                MT_ErrorParse("Such image type not handled yets (imres)\n", 0);
            }
          }
      break;
    case SSHORT:
    case USHORT:
      for (i=0 ; i<imageIn->dim.x*imageIn->dim.y*imageIn->dim.z; i++)
          if(par.set==0)
            resU16[i] = (LAB[(int)inU16[i]]==0)? 0 : inU16[i];
          else
          {
              switch (imres.type) {
              case UCHAR:
              case SCHAR:
                  resU8[i] = (unsigned char) LAB[(int)inU16[i]];
                  break;
              case USHORT:
              case SSHORT:
                  resU16[i] = (unsigned short int) LAB[(int)inU16[i]];
                  break;
              default:
                  free(LAB);
                  VT_FreeImage( imageIn );
                  VT_FreeImage( &imres );
                  MT_ErrorParse("Such image type not handled yets (imres)\n", 0);
              }
          }
      break;
    default:
	  VT_FreeImage( imageIn );
	  VT_FreeImage( &imres );
      free(LAB);
      VT_Error("image type unknown or not supported for this program",program);
      return( 0 );
  }

  /* Ecriture de imres */

  VT_WriteInrimage( &(imres) );

  /*--- liberations memoires ---*/
  VT_FreeImage( imageIn );
  VT_FreeImage( &imres );
  free(LAB);
  
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

      /*--- parametres de calcul ---*/
      else if ( strcmp ( argv[i], "-set" ) == 0 ) {
        par->set = 1;
      }

      /*--- traitement eventuel de l'image d'entree ---*/
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
        par->names.inv = 1;
        par->inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
        par->names.swap = 1;
      }


      /* Parametres de calcul */


      else if ( strcmp ( argv[i], "-l" ) == 0 || strcmp ( argv[i], "-labels" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -l...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->labels[par->nlabels]) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -l...\n", 0 );
        par->nlabels += 1;
        i += 1;
        while ( i < argc && par->nlabels < LABELSLENGTH)  {
          status = sscanf( argv[i],"%d",&(par->labels[par->nlabels]) );
          if ( status <= 0 ) break;
          par->nlabels += 1;
          i += 1;
        }
        i--;
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
  if (nb == 2) {
    strcpy( par->names.out, ">" );  /* standart output */
    if (par->nlabels>0)
      strncpy( par->names.out, par->names.ext, STRINGLENGTH);
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

  par->nlabels=0;
  par->inv=0;
  par->set=0;
}

