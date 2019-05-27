/****************************************************
 * extImage.c - 
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Fri Feb  8 10:32:01 CET 2013
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */



#include <extract.h>
#include <vt_common.h>






typedef struct local_par {
  vt_names names;

  int origin_definition;
  vt_4vpt origin;
  vt_4vpt dim;
  vt_4vpt slice;

  int extractRed;
  int extractGreen;
  int extractBlue;
  
} local_par;



/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );



static char *usage = "[image-in] [image-out]\n\
 [-origin %d %d [%d] | [-ix|-xorigin %d] [-iy|-yorigin %d] [-iz|-zorigin %d] [-vorigin %d]]\n\
 [-0 | -1]\n\
 [-dim %d %d [%d] | [-x|-xdim %d] [-y|-ydim %d] [-z|-zdim %d] [-vdim %d]]\n\
 [-template %s]\n\
 [-xy %d | -yz %d | -xz %d] [-xyz %d | -red | -green | -blue]\n\
 [-v] [-help] [-D]";



static char *detail = "\
\n\
   Extract a sub-image from 'image-in'.\n\
\n\
 [-origin %d %d [%d]]   # origin of the extracted image in the input image\n\
                          default is (0,0,0)\n\
 [-ix|-xorigin %d]      # X coordinates of the origin\n\
 [-iy|-yorigin %d]      # Y coordinates of the origin\n\
 [-iz|-zorigin %d]      # Z coordinates of the origin\n\
 [-vorigin %d]          # V coordinates of the origin\n\
 [-0]                   # first coordinate is 0 (default)\n\
 [-1]                   # first coordinate is 1\n\
                        # WARNING: all coordinates have to be given\n\
                        # (else the result will be erroneous)\n\
 [-dim %d %d [%d]]      # extracted image dimensions\n\
 [-x|-xdim %d]          # X dimension of the extracted image\n\
 [-y|-ydim %d]          # Y dimension of the extracted image\n\
 [-z|-zdim %d]          # Z dimension of the extracted image\n\
 [-vdim %d]             # V dimension of the extracted image\n\
 [-template %s]         # template image for the (x,y,z) dimensions\n\
                          of the extracted image\n\
 [-xy %d]   # extract the XY slice #%d\n\
            # this is equivalent to '-origin 0 0 %d -dim dimx dimy 1'\n\
 [-xz %d]   # extract the XZ slice #%d\n\
            # the extracted image has sizes (dimx, dimz, 1)\n\
            # this is different of '-origin 0 %d 0 -dim dimx 1 dimz'\n\
 [-yz %d]   # extract the YZ slice #%d\n\
            # the extracted image has sizes (dimy, dimz, 1)\n\
            # this is different of '-origin %d 0 0 -dim 1 dimy dimz'\n\
 [-xyz %d]  # extract the vectorial component #%d\n\
 [-red]     # extract the red component\n\
            # this is equivalent to '-vorigin 0 -vdim 1'\n\
 [-green]   # extract the green component\n\
            # this is equivalent to '-vorigin 1 -vdim 1'\n\
 [-blue]    # extract the blue component\n\
            # this is equivalent to '-vorigin 2 -vdim 1'\n\
\n";



static char program[STRINGLENGTH];



int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image;
  vt_image imres;
  vt_image *ptrImres = (vt_image *)NULL;
  vt_image imvres;
  vt_image *ptrImvres = (vt_image *)NULL;
  vt_image *imtemplate;
  
  int theDim[3];
  int resDim[3] = {-1, -1, -1};
  int resDimV = 0;
  int theLeftCorner[3] = {0, 0, 0};
  int theRightCorner[3] = {0, 0, 0};
  int theLeftV = 0;
  int theRightV;

  

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
    

  
  /*
   * if the indices do not begin with 0
   */
  if ( par.origin_definition > 0 ) {
    par.origin.v -= par.origin_definition;
    par.origin.x -= par.origin_definition;
    par.origin.y -= par.origin_definition;
    par.origin.z -= par.origin_definition;
    
    par.slice.v -= par.origin_definition;
    par.slice.x -= par.origin_definition;
    par.slice.y -= par.origin_definition;
    par.slice.z -= par.origin_definition;
  }

  if ( par.extractRed ) par.slice.v = 0;
  else if ( par.extractGreen ) par.slice.v = 1;
  else if ( par.extractBlue ) par.slice.v = 2;




  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);
  
  theDim[0] = image->dim.x * image->dim.v;
  theDim[1] = image->dim.y;
  theDim[2] = image->dim.z;



  /* initializing result image
     - with reference image, if any
     - with parameters, if any
  */

  VT_Image( &imres );

  if ( par.names.ext[0] != '\0' ) {
    imtemplate = _VT_Inrimage( par.names.ext );
    if ( imtemplate == (vt_image*)NULL ) {
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse("unable to read template image\n", 0);
    }
    
    theLeftCorner[0] = par.origin.x;
    theLeftCorner[1] = par.origin.y;
    theLeftCorner[2] = par.origin.z;

    theRightCorner[0] = par.origin.x + imtemplate->dim.x;
    theRightCorner[1] = par.origin.y + imtemplate->dim.y;
    theRightCorner[2] = par.origin.z + imtemplate->dim.z;

    VT_FreeImage( imtemplate );
    VT_Free( (void**)&imtemplate );

    if ( theLeftCorner[0] < 0 ) theLeftCorner[0] = 0;
    if ( theLeftCorner[1] < 0 ) theLeftCorner[1] = 0;
    if ( theLeftCorner[2] < 0 ) theLeftCorner[2] = 0;
    
    if ( theRightCorner[0] > (int)image->dim.x ) theRightCorner[0] = image->dim.x;
    if ( theRightCorner[1] > (int)image->dim.y ) theRightCorner[1] = image->dim.y;
    if ( theRightCorner[2] > (int)image->dim.z ) theRightCorner[2] = image->dim.z;

    ptrImres = &imres;
    VT_InitVImage( &imres, par.names.out, image->dim.v, 
                   theRightCorner[0] - theLeftCorner[0],
                   theRightCorner[1] - theLeftCorner[1],
                   theRightCorner[2] - theLeftCorner[2],
                   image->type );
    
    imres.siz.x = image->siz.x;
    imres.siz.y = image->siz.y;
    imres.siz.z = image->siz.z;

    resDim[0] = imres.dim.x * image->dim.v;
    resDim[1] = imres.dim.y;
    resDim[2] = imres.dim.z;
  }
      


  /* initialisation with dimension parameters
   */
  else if ( par.dim.x > 0 && par.dim.y > 0 ) {

    theLeftCorner[0] = par.origin.x;
    theLeftCorner[1] = par.origin.y;

    theRightCorner[0] = par.origin.x + par.dim.x;
    theRightCorner[1] = par.origin.y + par.dim.y;

    if ( theLeftCorner[0] < 0 ) theLeftCorner[0] = 0;
    if ( theLeftCorner[1] < 0 ) theLeftCorner[1] = 0;
    
    if ( theRightCorner[0] > (int)image->dim.x ) theRightCorner[0] = image->dim.x;
    if ( theRightCorner[1] > (int)image->dim.y ) theRightCorner[1] = image->dim.y;

    if ( par.dim.z > 0 ) {
      
      theLeftCorner[2] = par.origin.z;
      theRightCorner[2] = par.origin.z + par.dim.z;

      if ( theLeftCorner[2] < 0 ) theLeftCorner[2] = 0;
      if ( theRightCorner[2] > (int)image->dim.z ) theRightCorner[2] = image->dim.z;
      
    }
    else {

      theLeftCorner[2] = 0;
      theRightCorner[2] = 1;

    }
    
    ptrImres = &imres;
    VT_InitVImage( &imres, par.names.out, image->dim.v,
                   theRightCorner[0] - theLeftCorner[0],
                   theRightCorner[1] - theLeftCorner[1],
                   theRightCorner[2] - theLeftCorner[2],
                   image->type );
    
    imres.siz.x = image->siz.x;
    imres.siz.y = image->siz.y;
    imres.siz.z = image->siz.z;

    resDim[0] = imres.dim.x * image->dim.v;
    resDim[1] = imres.dim.y;
    resDim[2] = imres.dim.z;
  }



  /* initialisation with slice information
   */
  else if ( par.slice.z >= 0 && par.slice.z < (int)image->dim.z ) {

    ptrImres = &imres;
    VT_InitVImage( &imres, par.names.out, image->dim.v,
                   image->dim.x, image->dim.y, 1,
                   image->type );
    imres.siz.x = image->siz.x;
    imres.siz.y = image->siz.y;
    imres.siz.z = image->siz.z;

    resDim[0] = image->dim.x * image->dim.v;
    resDim[1] = image->dim.y;
    resDim[2] = 1;

    theLeftCorner[2] = par.slice.z;

  }

  else if ( par.slice.y >= 0 && par.slice.y < (int)image->dim.y ) {

    ptrImres = &imres;
    VT_InitVImage( &imres, par.names.out, image->dim.v,
                   image->dim.x, image->dim.z, 1,
                   image->type );
    imres.siz.x = image->siz.x;
    imres.siz.y = image->siz.z;
    imres.siz.z = image->siz.y;

    resDim[0] = image->dim.x * image->dim.v;
    resDim[1] = 1;
    resDim[2] = image->dim.z;

    theLeftCorner[1] = par.slice.y;

  }

  else if ( par.slice.x >= 0 && par.slice.x < (int)image->dim.x ) {

    ptrImres = &imres;
    VT_InitVImage( &imres, par.names.out, image->dim.v,
                   image->dim.y, image->dim.z, 1,
                   image->type );
    imres.siz.x = image->siz.y;
    imres.siz.y = image->siz.z;
    imres.siz.z = image->siz.x;

    resDim[0] = 1 * image->dim.v;
    resDim[1] = image->dim.y;
    resDim[2] = image->dim.z;

    theLeftCorner[0] = par.slice.x;

  }

  /* no dimension parameters
   */

  else {
    /* no X, Y or Z cropping
     */
    ptrImres = image;

    resDim[0] = image->dim.x * image->dim.v;
    resDim[1] = image->dim.y;
    resDim[2] = image->dim.z;

    theLeftCorner[0] = 0;
    theLeftCorner[1] = 0;
    theLeftCorner[2] = 0;
  }



  /***************************************************
   *
   * X, Y or Z cropping
   *
   ***************************************************/
  
  if ( ptrImres != image ) {
    if ( VT_AllocImage( ptrImres ) != 1 ) {
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse("unable to allocate output image (1)\n", 0);
  }

    if (  _VT_VERBOSE_ ||  _VT_DEBUG_ ) {
      fprintf( stderr, "... will extract a sub-image of dimension [%d %d %d] from origin (%d %d %d)\n",
               resDim[0], resDim[1], resDim[2], theLeftCorner[0],  theLeftCorner[1],  theLeftCorner[2] );
    }

    theLeftCorner[0] *= image->dim.v;
    ExtractFromBuffer( image->buf, theDim, imres.buf, resDim, theLeftCorner, image->type );

    VT_FreeImage( image );
    VT_Free( (void**)&image );
  }

  if ( ptrImres->dim.v == 1 ) {
    if ( VT_WriteInrimage( ptrImres ) == -1 ) {
      VT_FreeImage( ptrImres );
      if ( ptrImres == image ) VT_Free( (void**)&image );
      VT_ErrorParse("unable to write output image (1)\n", 0);
    }
    return( 0 );
  }

  /***************************************************
   *
   * V cropping
   *
   ***************************************************/

  theDim[0] = resDim[0] = ptrImres->dim.x;
  theDim[1] = resDim[1] = ptrImres->dim.y;
  theDim[2] = resDim[2] = ptrImres->dim.z;

  theLeftCorner[0] = 0;
  theLeftCorner[1] = 0;
  theLeftCorner[2] = 0;

  if ( par.dim.v > 0 ) {

    theLeftV = par.origin.v;
    theRightV = par.origin.v + par.dim.v;
    if ( theLeftV < 0 ) theLeftV = 0;
    if ( theRightV > (int)ptrImres->dim.v ) theRightV = ptrImres->dim.v;
    resDimV = theRightV-theLeftV;

    ptrImvres = &imvres;
    VT_InitVImage( &imvres, par.names.out, resDimV,
                   ptrImres->dim.x, ptrImres->dim.y, ptrImres->dim.z, ptrImres->type );
    imvres.siz.x = ptrImres->siz.x;
    imvres.siz.y = ptrImres->siz.y;
    imvres.siz.z = ptrImres->siz.z;
  }
  else if ( par.slice.v >= 0 || par.slice.v < (int)ptrImres->dim.v ) {
    theLeftV = par.slice.v;
    resDimV = 1;

    ptrImvres = &imvres;
    VT_InitVImage( &imvres, par.names.out, resDimV,
                   ptrImres->dim.x, ptrImres->dim.y, ptrImres->dim.z, ptrImres->type );
    imvres.siz.x = ptrImres->siz.x;
    imvres.siz.y = ptrImres->siz.y;
    imvres.siz.z = ptrImres->siz.z;
  }
  else {
    ptrImvres = ptrImres;
  }

  if ( ptrImvres != ptrImres ) {
    if ( VT_AllocImage( ptrImvres) != 1 ) {
      VT_FreeImage( ptrImres );
      if ( ptrImres == image ) VT_Free( (void**)&image );
      VT_ErrorParse("unable to allocate output image (2)\n", 0);
    }

    if (  _VT_VERBOSE_ ||  _VT_DEBUG_ ) {
      fprintf( stderr, "... will extract a sub-image of dimension [%d] from origin (%d)\n",
               resDimV, theLeftV );
    }

    VectorialExtractFromBuffer( ptrImres->buf, theDim, ptrImres->dim.v, ptrImvres->buf, resDim, resDimV, theLeftCorner, theLeftV, image->type );

    VT_FreeImage( ptrImres );
    if ( ptrImres == image ) VT_Free( (void**)&image );
  }
    
  /*--- ecriture de l'image resultat ---*/
  if ( VT_WriteInrimage( ptrImvres ) == -1 ) {
      VT_FreeImage( ptrImvres );
      if ( ptrImvres == image ) VT_Free( (void**)&image );
      VT_ErrorParse("unable to write output image(2)\n", 0);
  }

  /*--- liberations memoires ---*/
  VT_FreeImage( ptrImvres );
  if ( ptrImvres == image ) VT_Free( (void**)&image );

  return( 0 );
}



static void VT_Parse( int argc, char *argv[], local_par *par )
{
  int i, nb, status;
  char text[STRINGLENGTH];
  float t;

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

      /* some general options
       */

      else if ( strcmp ( argv[i], "--help" ) == 0 
                || ( strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0' )
                || ( strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0' )
                || ( strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0' ) ) {
        VT_ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-verbose" ) == 0
                || (strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0') ) {
        if ( _VT_VERBOSE_ <= 0 )
          _VT_VERBOSE_ = 1;
        else
          _VT_VERBOSE_ ++;
      }
      else if ( strcmp ( argv[i], "-debug" ) == 0
                || (strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0') ) {
        if ( _VT_DEBUG_  <= 0 )
          _VT_DEBUG_ = 1;
        else
          _VT_DEBUG_ ++;
      }

      /* origin of the subblock
       */
      
      else if ( strcmp ( argv[i], "-1" ) == 0 && argv[i][2] == '\0' ) {
        par->origin_definition = 1;
      }
      else if ( strcmp ( argv[i], "-0" ) == 0 && argv[i][2] == '\0' ) {
        par->origin_definition = 0;
      }
      
      else if ( (strcmp (argv[i], "-origin" ) == 0 && argv[i][7] == '\0')
                || (strcmp (argv[i], "-o" ) == 0 && argv[i][2] == '\0') ) {
        i ++;
        if ( i >= argc)    VT_ErrorParse( "parsing -origin %d\n", 0 );
        status = sscanf( argv[i], "%f", &t );
        if ( status <= 0 ) VT_ErrorParse( "parsing -origin %d\n", 0 );
        par->origin.x = (t > 0.0) ? (int)(t+0.5) : (int)(t-0.5);
        i ++;
        if ( i >= argc)    VT_ErrorParse( "parsing -origin %d %d\n", 0 );
        status = sscanf( argv[i], "%f", &t );
        if ( status <= 0 ) VT_ErrorParse( "parsing -origin %d %d\n", 0 );
        par->origin.y = (t > 0.0) ? (int)(t+0.5) : (int)(t-0.5);
        i ++;
        if ( i >= argc) par->origin.z = (par->origin_definition == 1) ? 1 : 0;
        else {
          status = sscanf( argv[i], "%f", &t );
          if ( status <= 0 ) {
            i--;
            par->origin.z = (par->origin_definition == 1) ? 1 : 0;
          }
          else {
            par->origin.z = (t > 0.0) ? (int)(t+0.5) : (int)(t-0.5);
          }
        }
      }

      else if ( (strcmp (argv[i], "-xorigin" ) == 0 && argv[i][8] == '\0')
                || (strcmp ( argv[i], "-ix" ) == 0 && argv[i][3] == '\0') ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -xorigin...\n", 0 );
        status = sscanf( argv[i], "%f", &t );
        if ( status <= 0 ) VT_ErrorParse( "parsing -xorigin", 0 );
        par->origin.x = (t > 0.0) ? (int)(t+0.5) : (int)(t-0.5);
      }
      else if ( (strcmp (argv[i], "-yorigin" ) == 0 && argv[i][8] == '\0')
                || (strcmp ( argv[i], "-iy" ) == 0 && argv[i][3] == '\0') ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -yorigin...\n", 0 );
        status = sscanf( argv[i], "%f", &t );
        if ( status <= 0 ) VT_ErrorParse( "parsing -yorigin", 0 );
        par->origin.y = (t > 0.0) ? (int)(t+0.5) : (int)(t-0.5);
      }
      else if ( (strcmp (argv[i], "-zorigin" ) == 0 && argv[i][8] == '\0')
                || (strcmp ( argv[i], "-iz" ) == 0 && argv[i][3] == '\0') ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -zorigin...\n", 0 );
        status = sscanf( argv[i], "%f", &t );
        if ( status <= 0 ) VT_ErrorParse( "parsing -zorigin...\n", 0 );
        par->origin.z = (t > 0.0) ? (int)(t+0.5) : (int)(t-0.5);
      }
      else if ( (strcmp (argv[i], "-vorigin" ) == 0 && argv[i][8] == '\0') ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -vorigin...\n", 0 );
        status = sscanf( argv[i], "%f", &t );
        if ( status <= 0 ) VT_ErrorParse( "parsing -vorigin...\n", 0 );
        par->origin.v = (t > 0.0) ? (int)(t+0.5) : (int)(t-0.5);
      }

      /* dimensions of the subblock
       */

      else if ( strcmp (argv[i], "-dim" ) == 0 && argv[i][4] == '\0' ) {
        i ++;
        if ( i >= argc)    VT_ErrorParse( "parsing -dim %d\n", 0 );
        status = sscanf( argv[i], "%d", &(par->dim.x) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -dim %d\n", 0 );
        i ++;
        if ( i >= argc)    VT_ErrorParse( "parsing -dim %d %d\n", 0 );
        status = sscanf( argv[i], "%d", &(par->dim.y) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -dim %d %d\n", 0 );
        i ++;
        if ( i >= argc) par->dim.z = 1;
        else {
          status = sscanf( argv[i], "%d", &(par->dim.z) );
          if ( status <= 0 ) {
            i--;
            par->dim.z = 1;
          }
        }
      }
      
      else if ( (strcmp ( argv[i], "-xdim") == 0 && argv[i][5] == '\0')
                || (strcmp ( argv[i], "-x") == 0 && argv[i][2] == '\0') ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -xdim ...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->dim.x) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -xdim ...\n", 0 );
      }
      else if ( (strcmp ( argv[i], "-ydim") == 0 && argv[i][5] == '\0')
                || (strcmp ( argv[i], "-y") == 0 && argv[i][2] == '\0') ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -ydim ...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->dim.y) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -ydim ...\n", 0 );
      }
      else if ( (strcmp ( argv[i], "-zdim") == 0 && argv[i][5] == '\0')
                || (strcmp ( argv[i], "-z") == 0 && argv[i][2] == '\0') ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -zdim ...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->dim.z) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -zdim ...\n", 0 );
      }
      else if ( (strcmp ( argv[i], "-vdim") == 0 && argv[i][5] == '\0') ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -vdim ...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->dim.v) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -vdim ...\n", 0 );
      }

      /* slices
       */

      else if ( strncmp ( argv[i], "-yz", 3 ) == 0 && argv[i][3] == '\0' ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -yz...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->slice.x) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -yz...\n", 0 );
      }
      else if ( strncmp ( argv[i], "-xz", 3 ) == 0 && argv[i][3] == '\0' ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -xz...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->slice.y) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -xz...\n", 0 );
      }
      else if ( strncmp ( argv[i], "-xy", 3 ) == 0 && argv[i][3] == '\0' ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -xy...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->slice.z) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -xy...\n", 0 );
      }
      else if ( strncmp ( argv[i], "-xyz", 4 ) == 0 && argv[i][4] == '\0' ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -xyz...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->slice.v) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -xyz...\n", 0 );
      }

      else if ( strncmp ( argv[i], "-red", 4 ) == 0 && argv[i][4] == '\0' ) {
        par->extractRed  = 1;
      }
      else if ( strncmp ( argv[i], "-green", 6 ) == 0 && argv[i][6] == '\0' ) {
        par->extractGreen  = 1;
      }
      else if ( strncmp ( argv[i], "-blue", 5 ) == 0 && argv[i][5] == '\0' ) {
        par->extractBlue  = 1;
      }

      /* template
       */

      else if ( strcmp ( argv[i], "-template") == 0
           || (strcmp ( argv[i], "-t") == 0 && argv[i][2] == '\0')
           || (strcmp ( argv[i], "-dims") == 0 && argv[i][5] == '\0') ) {
        i++;
        if ( i >= argc) VT_ErrorParse( "parsing -template\n", 0 );
        strncpy( par->names.ext, argv[i], STRINGLENGTH );
      }

      /*--- option inconnue ---*/
      else {
        sprintf(text,"unknown option %s\n",argv[i]);
        VT_ErrorParse(text, 0);
      }
    }
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

  par->origin_definition = 0;
  par->origin.v = par->origin.x = par->origin.y = par->origin.z = 0;
  par->dim.v = par->dim.x = par->dim.y = par->dim.z = -1;
  par->slice.v = par->slice.x = par->slice.y = par->slice.z = -1;

  par->extractRed = 0;
  par->extractGreen = 0;
  par->extractBlue = 0;
}
