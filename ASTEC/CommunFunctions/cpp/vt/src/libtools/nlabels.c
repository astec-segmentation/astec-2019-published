/*************************************************************************
 * nlabels.c - retourne le nombre de labels dans une image CHAR, SHORT ou INT
 *
 *
 * $Id: nlabels.c,v 1.0 2014/14/01 15:00:00 gael Exp $
 *
 * DESCRIPTION:
 *
 *
 *
 *
 *
 * AUTHOR:
 * Gael Michelin
 *
 * CREATION DATE:
 * Jan 2014
 *
 * Copyright Gregoire Malandain, INRIA
 *
 *
 * ADDITIONS, CHANGES:
 *
 *
 */

#include <vt_common.h>


typedef struct local_par {
  vt_names names;
} local_par;

/*------- Definition des fonctions statiques ----------*/

static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );



static char *usage = "[image-in]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
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
  int size;
  int i,j,k,n=0;
  int MAX;
  u8 V8[256];
  u16 *V16;
  u32 *V32;
  u16 *tmp16;
  u32 *tmp32;
  u16 val16;
  u32 val32;
  u8 *buf8;
  u16 *buf16;
  u32 *buf32;

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

  size = image->dim.x * image->dim.y * image->dim.z;

  switch ( image->type ) {
  case UCHAR :
  case SCHAR :
    if (_VT_VERBOSE_)
      fprintf(stdout, "TYPE : CHAR\n");
    buf8 = (u8*)image->buf;
    for (i=0 ; i<256 ; i++) V8[i]=(unsigned char)0;
    for (i=0 ; i<size ; i++)
      V8[buf8[i]] = (unsigned char) 255;
    for (i=0 ; i<256 ; i++) if (V8[i] != (unsigned char) 0) n++;
    break;
  case USHORT :
  case SSHORT :
    if (_VT_VERBOSE_)
      fprintf(stdout, "TYPE : SHORT\n");
    MAX=256;
    tmp16 = malloc(MAX*sizeof(u16));
    if (tmp16 != NULL)
    {
      V16=tmp16;
    }
    else
    {
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse( "error while reallocating label vector", 0 );
    }
    for (i=0 ; i<MAX ; i++) V16[i]=-1;
    buf16 = (u16*)image->buf;
    for (i=0 ; i<size ; i++)
    {
      if (n > 0)
      {
        if (buf16[i] == val16)
          continue;
        j=0;
        k=j-1;
        while (j<n && k!=j)
        {
          k=j;
          if (buf16[i] == V16[j])
          {
            val16=V16[j];
            continue;
          }
          j++;
        }
        if (j < n)
          continue;
      }
      /* Ici, on a affaire a un nouveau label...
       */
      if (n==MAX)
      {
        tmp16=realloc(V16, (MAX+256)*sizeof(u16));
        if (tmp16 != NULL)
        {
          MAX+=256;
          V16=tmp16;
        }
        else
        {
          free(V16);
          VT_FreeImage( image );
          VT_Free( (void**)&image );
          VT_ErrorParse( "error while reallocating label vector", 0 );
        }
      }
      V16[n++] = buf16[i];
      val16=buf16[i];
    }
    free(V16);
    break;
  case UINT :
  case SINT :
    if (_VT_VERBOSE_)
      fprintf(stdout, "TYPE : INT\n");
    MAX=256;
    tmp32 = malloc(MAX*sizeof(u32));
    if (tmp32 != NULL)
    {
      V32=tmp32;
    }
    else
    {
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse( "error while reallocating label vector", 0 );
    }
    for (i=0 ; i<MAX ; i++) V32[i]=-1;
    buf32 = (u32*)image->buf;
    for (i=0 ; i<size ; i++)
    {
      if (n > 0)
      {
        if (buf32[i] == val32)
          continue;
        j=0;
        k=j-1;
        while (j<n && k!=j)
        {
          k=j;
          if (buf32[i] == V32[j])
          {
            val32=V32[j];
            continue;
          }
          j++;
        }
        if (j < n)
          continue;
      }
      /* Ici, on a affaire a un nouveau label...
       */
      if (n==MAX)
      {
        tmp32=realloc(V32, (MAX+256)*sizeof(u32));
        if (tmp32 != NULL)
        {
          MAX+=256;
          V32=tmp32;
        }
        else
        {
          free(V32);
          VT_FreeImage( image );
          VT_Free( (void**)&image );
          VT_ErrorParse( "error while reallocating label vector", 0 );
        }
      }
      V32[n++] = buf32[i];
      val32=buf32[i];
    }
    free(V32);
    break;
  default :
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse( "image type does not correspond to label type", 0 );
  }

  VT_FreeImage( image );
  VT_Free( (void**)&image );
  /* if (_VT_VERBOSE_) */
  fprintf(stdout, "%d\n", n);
  return( n );
}











static void VT_Parse( int argc, char *argv[], local_par *par )
{
  int i, nb;
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
}
