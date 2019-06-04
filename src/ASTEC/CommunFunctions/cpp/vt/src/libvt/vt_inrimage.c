/*************************************************************************
 * vt_inrimage.c -
 *
 * $Id: vt_inrimage.c,v 1.10 2006/04/14 08:39:32 greg Exp $
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


#include <string.h>

#include <vtmalloc.h>

#include <vt_inrimage.h>

typedef enum {
  INRTYPE_UNKNOWN=0,
  INRTYPE_FLOAT=1,
  INRTYPE_UNSIGNEDFIXED=2,
  INRTYPE_SIGNEDFIXED=3
} InrimageType;

/*------- Definition des fonctions statiques ----------*/
#include <ImageIO.h> 








/************************************************************
 *
 * reading procedure
 *
 ************************************************************/



static int _imageToVtImage( _image *inr, vt_image *image )
{
  char *proc = "_imageToVtImage";
  ImageType type = TYPE_UNKNOWN;
  int i, j;

  /* conversion
   */
  switch ( inr->wordKind ) {
  case WK_FIXED :
    switch( inr->sign ) {
    case SGN_UNSIGNED :
      if ( inr->wdim  == sizeof( unsigned char ) ) {
        type = UCHAR;
      }
      else if ( inr->wdim  == sizeof( unsigned short int ) ) {
        type = USHORT;
      }
      else if ( inr->wdim  == sizeof( unsigned int ) ) {
        type = UINT;
      }
      else {
        if ( _VT_VERBOSE_ )
          fprintf( stderr, "%s: unknown WK_FIXED UNSIGNED word dim (%d) for image '%s'\n ", proc, inr->wdim, image->name );
        return 0;
      }
      break;
    case SGN_SIGNED :
      if ( inr->wdim  == sizeof( char ) ) {
        type = SCHAR;
      }
      else if ( inr->wdim  == sizeof( short int ) ) {
        type = SSHORT;
      }
      else if ( inr->wdim  == sizeof( int ) ) {
        type = SINT;
      }
      else {
        if ( _VT_VERBOSE_ )
          fprintf( stderr, "%s: unknown WK_FIXED SIGNED word dim (%d) for image '%s'\n ", proc, inr->wdim, image->name );
        return 0;
      }
      break;
    default :
      if ( _VT_VERBOSE_ )
        fprintf( stderr, "%s: unknown wordSign for image '%s'\n ", proc, image->name );
      return 0;    
    }
    break;
  case WK_FLOAT :
    if ( inr->wdim  == sizeof( float ) ) {
      type = FLOAT;
    }
    else if ( inr->wdim  == sizeof( double ) ) {
      type = DOUBLE;
    }
    else {
      if ( _VT_VERBOSE_ )
        fprintf( stderr, "%s: unknown WK_FLOAT word dim for image '%s'\n ", proc, image->name );
      return 0;
    }
    break;
  default :
    if ( _VT_VERBOSE_ )
      fprintf( stderr, "%s: unknown wordKind for image '%s'\n ", proc, image->name );
    return 0;
  }
  
  /* copie des attributs
   */
  image->type = type;

  image->dim.x = inr->xdim; 
  image->dim.y = inr->ydim; 
  image->dim.z = inr->zdim; 
  image->dim.v = inr->vdim; 

  image->siz.x = inr->vx;
  image->siz.y = inr->vy;
  image->siz.z = inr->vz;

  if ( inr->t_is_set ) {
    image->off.x = inr->tx;
    image->off.y = inr->ty;
    image->off.z = inr->tz;
    image->off_is_set = 1;
  }

  if ( inr->q_is_set ) {
    image->quat.x = inr->qb;
    image->quat.y = inr->qc;
    image->quat.z = inr->qd;
    image->quat_is_set = 1;
  }

  if ( inr->qform_code ) {
    for ( i=0; i<4; i++ )
    for ( j=0; j<4; j++ )
      image->qform_rfv[i][j] = inr->qform_toreal[i][j];
    image->geometry = _VT_QFORM_GEOMETRY_;
    image->qform_code = inr->qform_code;
  }
  else {
    for ( i=0; i<4; i++ )
    for ( j=0; j<4; j++ )
      image->qform_rfv[i][j] = 0.0;
    image->qform_rfv[0][0] = inr->vx;
    image->qform_rfv[1][1] = inr->vy;
    image->qform_rfv[2][2] = inr->vz;
    image->qform_rfv[3][3] = 1.0;
    image->geometry = _VT_DEFAULT_GEOMETRY_;
  }

  if ( inr->sform_code ) {
    for ( i=0; i<4; i++ )
    for ( j=0; j<4; j++ )
      image->sform_rfv[i][j] = inr->sform_toreal[i][j];
    image->sform_code = inr->sform_code;
  }


  return 1;
}




/* Lecture d'une image.

   Cette fonction essaye de lire le fichier de nom name 
   (si c'est "<", c'est le standard input). Si l'en-tete
   est correctement lu, une structure image est allouee
   et on lit l'image que l'on range dans le buffer de 
   la structure. En cas d'erreur, la structure est desallouee.
   Ne reconnait qu'INRIMAGE-4.

RETURN
  Retourne (vt_image*)NULL en cas d'erreur.
*/


int VT_ReadInrimage( vt_image *image, char *name )
{
  char *proc = "VT_ReadInrimage";
  _image *inr;

  VT_FreeImage( image );
  /* lecture
   */
  inr = _readImage( name );
  if( !inr ) {
    if ( _VT_VERBOSE_ )
      fprintf( stderr, "%s: unable to read '%s'\n ", proc, name );
    return 0;
  }

  strcpy( image->name, name );

  if ( _imageToVtImage( inr, image ) != 1 ) {
    _freeImage( inr );
    if ( _VT_VERBOSE_ )
      fprintf( stderr, "%s: unable to decode header of '%s'\n ", proc, name );
    return 0;
  }

  image->buf = inr->data;
  if ( VT_AllocArrayImage( image ) != 1 ) {
    if ( _VT_VERBOSE_ )
      fprintf( stderr, "%s: unable to build image array\n", proc );
    _freeImage(inr);
    return 0;
  }

  /* deallocates the _image structure
   */
  _freeImageStructure( inr );

  return 1;
}





/* Lecture d'une image avec creation de la structure associee.

   Cette fonction essaye de lire le fichier de nom name
   (si c'est "<", c'est le standard input). Si l'en-tete
   est correctement lu, une structure image est allouee
   et on lit l'image que l'on range dans le buffer de
   la structure. En cas d'erreur, la structure est desallouee.
   Ne reconnait qu'INRIMAGE-4.

RETURN
  Retourne (vt_image*)NULL en cas d'erreur.
*/
vt_image* _VT_Inrimage( char *name /* file name containing the inrimage image */ )
{
  char *proc = "_VT_Inrimage";
  vt_image *image;


  image = (vt_image*)VT_Malloc( (unsigned)sizeof(vt_image) );
  if (image == (vt_image*)NULL) {
    VT_Error("allocation of image structure failed", proc );
    return( (vt_image*)NULL );
  }
  VT_Image( image );


  if ( VT_ReadInrimage( image, name ) != 1 ) {
    VT_Error("Unable to read file", proc );
    VT_Free( (void**)&image );
    return( (vt_image*)NULL );
  }

  return( image );
}









/************************************************************
 *
 * writing procedure
 *
 ************************************************************/



static int _vtImageToImage( vt_image *image,  _image *inr )
{
  char *proc = "_vtImageToImage";
  int i, j;

  if ( image == (vt_image*)NULL ) {
    if ( _VT_VERBOSE_ )
      fprintf( stderr, "%s: input image is a NULL pointer\n", proc );
    return( 0 );
  }

  switch ( image->type ) {
  case SCHAR :
    inr->wordKind = WK_FIXED;
    inr->sign     = SGN_SIGNED;
    inr->wdim     = sizeof( char );
    break;
  case UCHAR :
    inr->wordKind = WK_FIXED;
    inr->sign     = SGN_UNSIGNED;
    inr->wdim     = sizeof( char );
    break;
  case SSHORT :
    inr->wordKind = WK_FIXED;
    inr->sign     = SGN_SIGNED;
    inr->wdim     = sizeof( short int );
    break;
  case USHORT :
    inr->wordKind = WK_FIXED;
    inr->sign     = SGN_UNSIGNED;
    inr->wdim     = sizeof( short int );
    break;
  case FLOAT :
    inr->wordKind = WK_FLOAT;
    inr->sign     = SGN_UNKNOWN;
    inr->wdim     = sizeof( float );
    break;
  case DOUBLE :
    inr->wordKind = WK_FLOAT;
    inr->sign     = SGN_UNKNOWN;
    inr->wdim     = sizeof( double );
    break;
  case SINT :
    inr->wordKind = WK_FIXED;
    inr->sign     = SGN_SIGNED;
    inr->wdim     = sizeof( int );
    break;
  case UINT :
    inr->wordKind = WK_FIXED;
    inr->sign     = SGN_UNSIGNED;
    inr->wdim     = sizeof( int );
    break;
  default :
    if ( _VT_VERBOSE_ )
      fprintf( stderr, "%s: unknown type for image '%s'\n ", proc, image->name );
    return 0;
  }

  inr->xdim = image->dim.x;
  inr->ydim = image->dim.y;
  inr->zdim = image->dim.z;
  inr->vdim = image->dim.v;

  inr->vx = image->siz.x;
  inr->vy = image->siz.y;
  inr->vz = image->siz.z;

  if ( image->off_is_set ) {
    inr->tx = image->off.x;
    inr->ty = image->off.y;
    inr->tz = image->off.z;
    inr->t_is_set = 1;
  }

  if ( image->quat_is_set ) {
    inr->qb = image->quat.x;
    inr->qc = image->quat.y;
    inr->qd = image->quat.z;
    inr->q_is_set = 1;
  }

  if ( image->geometry == _VT_QFORM_GEOMETRY_ || image->qform_code ) {
    for ( i=0; i<4; i++ )
    for ( j=0; j<4; j++ )
      inr->qform_toreal[i][j] = image->qform_rfv[i][j];
    inr->qform_code = image->qform_code;
  }

  if ( image->sform_code ) {
    for ( i=0; i<4; i++ )
    for ( j=0; j<4; j++ )
      inr->sform_toreal[i][j] = image->sform_rfv[i][j];
    inr->sform_code = image->sform_code;
  }

  inr->data = image->buf;

  return( 1 );
}








int VT_WriteInrimageWithName( vt_image *image, char *name )
{
  char *proc = "VT_WriteInrimageWithName";
  _image *inr = NULL;

  inr = _initImage();
  if ( _vtImageToImage( image, inr ) != 1 ) {
    if ( _VT_VERBOSE_ || _VT_DEBUG_ )
      fprintf( stderr, "%s: unable to translate vt_image\n", proc );
    vtfree( inr );
    return 0;
  }

  if ( _writeImage( inr, name ) == -1 ) {
    if ( _VT_VERBOSE_ )
      fprintf( stderr, "%s: unable to write image '%s'", proc, image->name );
    vtfree( inr );
    return 0;
  }

  vtfree( inr );
  return 1;
}


/* Fonction d'ecriture d'une image inrimage.

RETURN
   Retourne 0 en cas de probleme.
*/


int VT_WriteInrimage( vt_image *image )
{
  return( VT_WriteInrimageWithName( image, image->name ) );
}





int VT_WriteInrimageHeaderWithName( vt_image *image, char *name )
{
  char *proc = "VT_WriteInrimageHeaderWithName";
  _image *inr = NULL;

  inr = _initImage();
  if ( _vtImageToImage( image, inr ) != 1 ) {
    if ( _VT_VERBOSE_ || _VT_DEBUG_ )
      fprintf( stderr, "%s: unable to translate vt_image\n", proc );
    vtfree( inr );
    return 0;
  }

  if ( _writeImageHeader( inr, name ) != 0 ) {
    if ( _VT_VERBOSE_ || _VT_DEBUG_ )
      fprintf( stderr, "%s: unable to write image header '%s'\n", proc, image->name );
    vtfree( inr );
    return 0;
  }

  vtfree( inr );

  return 1;
}


/* Fonction d'ecriture d'une image inrimage.

RETURN
   Retourne 0 en cas de probleme.
*/


int VT_WriteInrimageHeader( vt_image *image )
{
  return( VT_WriteInrimageHeaderWithName( image, image->name ) );
}













