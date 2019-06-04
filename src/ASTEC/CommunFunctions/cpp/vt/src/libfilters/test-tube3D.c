/*************************************************************************
 * minimum.c -
 *
 * $Id: test-tube3D.c,v 1.5 2001/06/08 15:10:50 greg Exp $
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

#include <vt_common.h>

#include <tube3D.h>
#include <tube3Dutil.h>
#include <connexe.h>

typedef struct local_par {
  vt_names names;
  int type;
} local_par;




/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[image-in] [image-out]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n\
\t options-de-type : -o 1    : unsigned char\n\
\t                   -o 2    : unsigned short int\n\
\t                   -o 2 -s : short int\n\
\t                   -o 4 -s : int\n\
\t                   -r      : float\n\
\t si aucune de ces options n'est presente, on prend le type de 'image-in'\n\
\n\
 $Revision: 1.5 $ $Date: 2001/06/08 15:10:50 $ $Author: greg $\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image;
  vt_image  imrep, imthe, imphi, imext;
  vt_image  imbin, imthi, imbrb, imcln;
  int theDim[3];
  double scale5[5] = { 1.0, 1.5, 2.0, 2.5, 3.0 };
  double scale3[5] = { 1.0, 2.0, 3.0 };
  int nb_scales = 3;

  typeFiberParameter *fiber_params = NULL;
  int n, nb_fibers;

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
  
  /*--- initialisation de l'image resultat ---*/
  VT_InitFromImage( &imrep, image, par.names.out, FLOAT );
  VT_InitFromImage( &imthe, image, par.names.out, FLOAT );
  VT_InitFromImage( &imphi, image, par.names.out, FLOAT );
  VT_InitFromImage( &imext, image, par.names.out, FLOAT );
  VT_InitFromImage( &imbin, image, par.names.out, UCHAR );
  VT_InitFromImage( &imthi, image, par.names.out, UCHAR );
  VT_InitFromImage( &imbrb, image, par.names.out, UCHAR );
  VT_InitFromImage( &imcln, image, par.names.out, UCHAR );
  sprintf( imrep.name, "%s.rep.inr", par.names.out );
  sprintf( imthe.name, "%s.the.inr", par.names.out );
  sprintf( imphi.name, "%s.phi.inr", par.names.out );
  sprintf( imext.name, "%s.ext.inr", par.names.out );
  sprintf( imbin.name, "%s.bin.inr", par.names.out );
  sprintf( imthi.name, "%s.thi.inr", par.names.out );
  sprintf( imbrb.name, "%s.brb.inr", par.names.out );
  sprintf( imcln.name, "%s.cln.inr", par.names.out );

  if ( VT_AllocImage( &imrep ) != 1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to allocate output image (1) \n", 0);
  }
  if ( VT_AllocImage( &imthe ) != 1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to allocate output image (2) \n", 0);
  }
  if ( VT_AllocImage( &imphi ) != 1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to allocate output image (3) \n", 0);
  }

  _VerboseInTube3D();
  _VerboseInTube3D();
  _VerboseInTube3Dutil();
  _VerboseInTube3Dutil();

  theDim[0] = image->dim.x;
  theDim[1] = image->dim.y;
  theDim[2] = image->dim.z;

  fprintf( stderr, " filtering ...\n" );

  switch ( nb_scales ) {
  case 5 :
    if ( compute_3D_line_response_at_multiple_scales( image->buf, image->type,
						      imrep.buf, imthe.buf, imphi.buf, theDim, scale5, 5 ) != 1 ) {
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse("error in computation\n", 0);
    }
    break;
  default :      
  case 3 :
    if ( compute_3D_line_response_at_multiple_scales( image->buf, image->type,
						      imrep.buf, imthe.buf, imphi.buf, theDim, scale3, 3 ) != 1 ) {
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse("error in computation\n", 0);
    }
  }
  
  VT_FreeImage( image );
  VT_Free( (void**)&image );


  if ( VT_AllocImage( &imext ) != 1 ) {
    VT_ErrorParse("unable to allocate output image (4) \n", 0);
  }
  
  fprintf( stderr, " extracting extrema ...\n" );
  if ( compute_3D_line_response_extrema( imext.buf, imrep.buf, imthe.buf, imphi.buf, theDim) != 1 ) {
    VT_ErrorParse("error in extrema computation\n", 0);
  
  }

  (void)VT_WriteInrimage( &imrep );
  (void)VT_WriteInrimage( &imthe );
  (void)VT_WriteInrimage( &imphi );
  VT_FreeImage( &imrep );
  VT_FreeImage( &imthe );
  VT_FreeImage( &imphi );


  if ( VT_AllocImage( &imbin ) != 1 ) {
    VT_ErrorParse("unable to allocate output image (5) \n", 0);
  }



  fprintf( stderr, " hysteresis thresholding ...\n" );
  
  switch ( theDim[0] ) {
  default :
  case  64 :
    if ( HysteresisThresholdingWithAllParams( imext.buf, FLOAT,
					      imbin.buf, UCHAR,
					      theDim, 
					      6.0, 12.0, /* low and high threshold */
					      26, /* connectivity */
					      50, /* minimal size of component */
					      1,  /* at least one point above the high threshold
						     in each connected component */
					      0,  /* all components */
					      1 /* binary output */ ) == -1 ) {
      VT_ErrorParse("error in thresholding\n", 0);
    }
    break;
      case  32 :
    if ( HysteresisThresholdingWithAllParams( imext.buf, FLOAT,
					      imbin.buf, UCHAR,
					      theDim, 
					      6.0, 12.0, /* low and high threshold */
					      26, /* connectivity */
					      10, /* minimal size of component */
					      1,  /* at least one point above the high threshold
						     in each connected component */
					      0,  /* all components */
					      1 /* binary output */ ) == -1 ) {
      VT_ErrorParse("error in thresholding\n", 0);
    }
    break;

  }


  (void)VT_WriteInrimage( &imext );
  VT_FreeImage( &imext );




  /*  _VT_BDD_GREG_CURVES 
      libvt/vt_bdd_amincir.c
      vt_amliste.c:vt_pt_amincir *_VT_ThinPtList( vt_image *image, int dim, int *nb )
   */
  if ( VT_AllocImage( &imthi ) != 1 ) {
    VT_ErrorParse("unable to allocate output image (6) \n", 0);
  }

  fprintf( stderr, " thinning ...\n" );
  if ( thin_3D_thick_lines( imbin.buf, imthi.buf, theDim ) != 1 ) {
    VT_ErrorParse("error in thinning\n", 0);
  }

  (void)VT_WriteInrimage( &imbin );
  VT_FreeImage( &imbin );




  if ( VT_AllocImage( &imbrb ) != 1 ) {
    VT_ErrorParse("unable to allocate output image (6) \n", 0);
  }
  if ( remove_small_simple_curves( imthi.buf, imbrb.buf, theDim, 10 ) != 1  ) {
    VT_ErrorParse("error in removing small branches\n", 0);
  }
  




  (void)VT_WriteInrimage( &imthi );
  VT_FreeImage( &imthi );


  if ( VT_AllocImage( &imcln ) != 1 ) {
    VT_ErrorParse("unable to allocate output image (7) \n", 0);
  }
  if ( remove_non_simple_components( imbrb.buf, imcln.buf, theDim, 10 ) != 1  ) {
    VT_ErrorParse("error in extracting valid components\n", 0);
  }





  fiber_params = compute_fibers_parameters( imcln.buf, theDim, 10, &nb_fibers );
  if ( fiber_params == NULL ) {
    VT_ErrorParse("error in computing parameters\n", 0);
  }
  
  printf( "\n" );
  for ( n=1; n <= nb_fibers; n++ ) {
    printf( "(1) FIBRE #%3d: size=%3d (%2d %2d %2d)-(%2d %2d %2d)\n",
	    n, fiber_params[n].size,
	    fiber_params[n].extremity_1[0], fiber_params[n].extremity_1[1], 
	    fiber_params[n].extremity_1[2], fiber_params[n].extremity_2[0], 
	    fiber_params[n].extremity_2[1], fiber_params[n].extremity_2[2] );
  }
  printf( "\n" );

  free( fiber_params );

  
  


  fiber_params = compute_fibers_parameters( imbrb.buf, theDim, 10, &nb_fibers );
  if ( fiber_params == NULL ) {
    VT_ErrorParse("error in computing parameters (2)\n", 0);
  }
  
  printf( "\n" );
  for ( n=1; n <= nb_fibers; n++ ) {
    printf( "(2) FIBRE #%3d: size=%3d (%2d %2d %2d)-(%2d %2d %2d)\n",
	    n, fiber_params[n].size,
	    fiber_params[n].extremity_1[0], fiber_params[n].extremity_1[1], 
	    fiber_params[n].extremity_1[2], fiber_params[n].extremity_2[0], 
	    fiber_params[n].extremity_2[1], fiber_params[n].extremity_2[2] );
  }
  printf( "\n" );

  free( fiber_params );

  
  






  (void)VT_WriteInrimage( &imbrb );
  VT_FreeImage( &imbrb );

  (void)VT_WriteInrimage( &imcln );
  VT_FreeImage( &imcln );







  return( 1 );
}








static void VT_Parse( int argc, 
		      char *argv[], 
		      local_par *par )
{
  int i, nb, status;
  int o=0, s=0, r=0;
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
      /*--- lecture du type de l'image de sortie ---*/
      else if ( strcmp ( argv[i], "-r" ) == 0 ) {
	r = 1;
      }
      else if ( strcmp ( argv[i], "-s" ) == 0 ) {
	s = 1;
      }
      else if ( strcmp ( argv[i], "-o" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -o...\n", 0 );
	status = sscanf( argv[i],"%d",&o );
	if ( status <= 0 ) VT_ErrorParse( "parsing -o...\n", 0 );
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
  
  /*--- type de l'image resultat ---*/
  if ( (o == 1) && (s == 1) && (r == 0) )  par->type = SCHAR;
  if ( (o == 1) && (s == 0) && (r == 0) ) par->type = UCHAR;
  if ( (o == 2) && (s == 0) && (r == 0) ) par->type = USHORT;
  if ( (o == 2) && (s == 1) && (r == 0) )  par->type = SSHORT;
  if ( (o == 4) && (s == 1) && (r == 0) )  par->type = SINT;
  if ( (o == 0) && (s == 0) && (r == 1) )  par->type = FLOAT;
  /* if ( par->type == TYPE_UNKNOWN ) VT_Warning("no specified type", program); */
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
  par->type = TYPE_UNKNOWN;
}
