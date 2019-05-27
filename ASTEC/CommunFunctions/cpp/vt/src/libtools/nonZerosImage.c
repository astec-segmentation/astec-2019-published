/*************************************************************************
 * nonZerosImage.c - retourne 1 si l'image n'est pas constituee que de zeros, 0 sinon
 *
 *
 * $Id: nonZerosImage.c,v 1.0 2017/14/03 14:33:00 gael Exp $
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
 * Mar 2017
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
\t -D : mode debug\n";

static char program[STRINGLENGTH];

int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image;
  int size;
  int i;
  u8 *buf8;
  u16 *buf16;
  u32 *buf32;
  u64 *buf64;
  r32 *bufFlt;
  r64 *bufDbl;

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
    if (_VT_DEBUG_)
      fprintf(stdout, "TYPE : CHAR\n");
    buf8 = (u8*)image->buf;
    for (i=0 ; i<size ; i++)
      if (buf8[i] != (unsigned char) 0)
      {
          VT_FreeImage( image );
          VT_Free( (void**)&image );
          if (_VT_VERBOSE_)
              fprintf(stdout, "1\n");
          return 1;
      }
    break;
  case USHORT :
  case SSHORT :
    if (_VT_DEBUG_)
      fprintf(stdout, "TYPE : SHORT\n");
    buf16 = (u16*)image->buf;
    for (i=0 ; i<size ; i++)
      if (buf16[i] != (unsigned short int) 0)
      {
          VT_FreeImage( image );
          VT_Free( (void**)&image );
          if (_VT_VERBOSE_)
              fprintf(stdout, "1\n");
          return 1;
      }
    break;
  case UINT :
  case SINT :
    if (_VT_DEBUG_)
      fprintf(stdout, "TYPE : INT\n");
    buf32 = (u32*)image->buf;
    for (i=0 ; i<size ; i++)
      if (buf32[i] != (unsigned int) 0)
      {
          VT_FreeImage( image );
          VT_Free( (void**)&image );
          if (_VT_VERBOSE_)
              fprintf(stdout, "1\n");
          return 1;
      }
    break;
  case ULINT :
  case SLINT :
    if (_VT_DEBUG_)
      fprintf(stdout, "TYPE : LONG\n");
    buf64 = (u64*)image->buf;
    for (i=0 ; i<size ; i++)
      if (buf64[i] != (unsigned long int) 0)
      {
          VT_FreeImage( image );
          VT_Free( (void**)&image );
          if (_VT_VERBOSE_)
              fprintf(stdout, "1\n");
          return 1;
      }
    break;
  case FLOAT:
      if (_VT_DEBUG_)
        fprintf(stdout, "TYPE : FLOAT\n");
      bufFlt = (r32*)image->buf;
      for (i=0 ; i<size ; i++)
        if (bufFlt[i] != (float) 0)
        {
            VT_FreeImage( image );
            VT_Free( (void**)&image );
            if (_VT_VERBOSE_)
                fprintf(stdout, "1\n");
            return 1;
        }
      break;
  case DOUBLE:
      if (_VT_DEBUG_)
        fprintf(stdout, "TYPE : DOUBLE\n");
      bufDbl = (r64*)image->buf;
      for (i=0 ; i<size ; i++)
        if (bufDbl[i] != (double) 0)
        {
            VT_FreeImage( image );
            VT_Free( (void**)&image );
            if (_VT_VERBOSE_)
                fprintf(stdout, "1\n");
            return 1;
        }
      break;
  default :
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse( "image type does not correspond to label type", 0 );
  }

  VT_FreeImage( image );
  VT_Free( (void**)&image );
  if (_VT_VERBOSE_)
      fprintf(stdout, "1\n");
  return( 0 );
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
