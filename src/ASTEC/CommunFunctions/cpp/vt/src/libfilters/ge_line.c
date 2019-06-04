#include <vt_common.h>

#include <vt_geline.h>
#include <vt_geline2D.h>
#include <vt_geline3D.h>
#include <vt_gelrec3D.h>

#include <connexe.h>

#define _EXTREMA_   1
#define _REPONSE_   2
#define _ECHELLE_   3
#define _EXTREMA2_  4
#define _SEUILLAGE_ 5
#define _RECONS_    6

#define _BIN_   0
#define _GREY_  1
#define _SCALE_ 2
#define _CSTE_  3

typedef struct local_par {
    /*--- dimension ---*/
    int dim;
    /*--- seuillage ---*/
    float mul_max_sh;
    float mul_max_sb;
    int size_cc;
    /*--- type de sortie ---*/
    int type_output;
    int type_recons;
    float rayon_recons;
    /*--- parametres de calcul ---*/
    vt_line par_line;
    vt_names names;
} local_par;



/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );



static char *usage = "[image-in] [image-out]\n\
\t [-first %f] [-last %f] [-nb %d] [-black | -white] [-csb %f] [-csh %f] [-tcc %d]\n\
\t [-reponse | -echelle | -extrema | -extrema2 | -seuillage | -recons [-grey | -scale | -rad %f]]\n\
\t [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -first %f : premier coefficient de filtrage (d'echelle)\n\
\t -last %f  : dernier coefficient de filtrage (d'echelle)\n\
\t -nb %d    : nombre d'echelles\n\
\t -reponse  : output = image des reponses maximales\n\
\t -echelle  : output = image des echelles correpondant aux reponses maximales\n\
\t -extrema  : output = image des reponses maximales egalement extrema\n\
\t -extrema2 : output = extrema calcules sur l'image des reponses maximales\n\
\t -recons   : output = reconstruction\n\
\t -csb %f : coefficient multiplicateur du maximum pour le seuil bas (1/6)\n\
\t -csh %f : coefficient multiplicateur du maximum pour le seuil haut (1/3)\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];





int main( int argc, char *argv[] )
{
        local_par par;
	vt_image *image;
	vt_images images;
	vt_resline res;
	r32 *fbuf;
	u8 *ubuf;
	u8 *sbuf;
	register int i, v;
	vt_3m m;

	/*--- lecture des parametres ---*/
	VT_Parse( argc, argv, &par );

	/*--- lecture de l'image d'entree ---*/
	image = _VT_Inrimage( par.names.in );
	if ( image == (vt_image*)NULL ) 
		VT_ErrorParse("unable to read input image\n", 0);
	if ( image->dim.z == 1 ) par.dim = VT_2D;

	/*--- operations eventuelles sur l'image d'entree ---*/
	if ( par.names.inv == 1 )  VT_InverseImage( image );
	if ( par.names.swap == 1 ) VT_SwapImage( image );

	/*--- initialisation de l'image resultat ---*/
	VT_InitRes( &res );
	if ( VT_AllocRes( &res, image, par.dim ) != 1 ) {
	    VT_FreeImage( image );
	    VT_Free( (void**)&image );
	    VT_ErrorParse("unable to allocate output image\n", 0);
        }

	/*--- calcul ---*/
	VT_InitImages( &images );
	if ( VT_AllocImages( &images, image, par.dim ) != 1 ) {
	  VT_FreeImage( image );
	  VT_Free( (void**)&image );
	  VT_FreeRes( &res, par.dim );
	  VT_ErrorParse("unable to allocate images\n", 0);
	}

	if ( par.dim == VT_2D ) {
	    if ( VT_ComputeLine2D( &res, &images, image, &(par.par_line) ) != 1 ) {
		VT_FreeImage( image );
		VT_Free( (void**)&image );
		VT_FreeRes( &res, par.dim );
		VT_FreeImages( &images, par.dim );
		VT_ErrorParse("unable to compute result\n", 0);
	    }
	} else {
	    if ( VT_ComputeLine3D( &res, &images, image, &(par.par_line) ) != 1 ) {
		VT_FreeImage( image );
		VT_Free( (void**)&image );
		VT_FreeRes( &res, par.dim );
		VT_FreeImages( &images, par.dim );
		VT_ErrorParse("unable to compute result\n", 0);
	    }
	}

	/*--- calcul de nouveaux extrema ? ---*/
	switch ( par.type_output ) {	 
	case _SEUILLAGE_ :
	case _RECONS_ :
	case _EXTREMA2_ :
	    if ( par.dim == VT_2D ) {
	      if ( VT_ExtractMaxRad2D( &res, &images, &(par.par_line) ) != 1 ) {
		VT_FreeImage( image );
		VT_Free( (void**)&image );
		VT_FreeRes( &res, par.dim );
		VT_ErrorParse("unable to compute extrema\n", 0);
	      }
	    } else {
		if ( VT_ExtractMaxRad3D( &res, &images, &(par.par_line) ) != 1 ) {
		    VT_FreeImage( image );
		    VT_Free( (void**)&image );
		    VT_FreeRes( &res, par.dim );
		    VT_ErrorParse("unable to compute extrema\n", 0);
		}
	    }
	}

	/*--- liberation des images intermediaires ---*/
	VT_FreeImages( &images, par.dim );

	/*--- ecriture de l'image resultat ---*/
	switch ( par.type_output ) {	    
	case _ECHELLE_ :
	    strcpy( res.imscale.name, par.names.out );
	    if ( VT_WriteInrimage( &(res.imscale) ) == -1 ) {
                VT_FreeImage( image );
                VT_Free( (void**)&image );
		VT_FreeRes( &res, par.dim );
                VT_ErrorParse("unable to write output image\n", 0);
	    }
	    break;
	case _REPONSE_ :
	    strcpy( res.imres.name, par.names.out );
	    if ( VT_WriteInrimage( &(res.imres) ) == -1 ) {
                VT_FreeImage( image );
                VT_Free( (void**)&image );
		VT_FreeRes( &res, par.dim );
                VT_ErrorParse("unable to write output image\n", 0);
	    }
	    break;
	case _SEUILLAGE_ :
	case _RECONS_ :
	    /*--- on masque l'image des reponses avec celle des extrema ---*/
	    fbuf = (r32*)(res.imres.buf);
	    ubuf = (u8*)(res.imext.buf);
	    v = image->dim.x * image->dim.y * image->dim.z;
	    for ( i = 0; i < v; i ++, ubuf ++, fbuf ++) {
		if ( (*ubuf) == (u8)0 ) *fbuf = 0.0;
	    }
	    /*--- calcul des min, moy et max des valeurs des extrema ---*/
	    if ( VT_3m( &(res.imres), &m ) != 1 ) {
		VT_FreeImage( image );
                VT_Free( (void**)&image );
		VT_FreeRes( &res, par.dim );
                VT_ErrorParse("unable to compute maximum value of extrema\n", 0);
	    }
	    /*--- seuillage : on passe dans EpidaureLib ---*/
	    {
	      int dim[3];
	      float lowThres = (float)( m.max * par.mul_max_sb );
	      float higThres = (float)( m.max * par.mul_max_sh );
	      char message[256];

	      dim[0] = image->dim.x;
	      dim[1] = image->dim.y;
	      dim[2] = image->dim.z;

	      Connexe_SetConnectivity( 26 );
	      Connexe_SetMinimumSizeOfComponents( par.size_cc );

	      sprintf( message," seuillage avec seuil bas = %f, seuil haut = %f, tcc = %d", 
		       (float)( lowThres ), (float)( higThres ), par.size_cc );
	      VT_Message( message, "" );
	      if ( HysteresisThresholding( res.imres.buf, FLOAT,
					   res.imext.buf, UCHAR, dim,
					   lowThres, higThres ) <= 0 ) {
		VT_FreeImage( image );
		VT_Free( (void**)&image );
		VT_FreeRes( &res, par.dim );
		VT_ErrorParse("unable to threshold image\n", 0);
	      }
	    }
	    /*--- seuillage fini ---*/
	    if ( par.type_output == _RECONS_ ) {
	      /*--- on masque l'image des echelles avec le resultat du seuillage ---*/
	      ubuf = (u8*)(res.imext.buf);
	      sbuf = (u8*)(res.imscale.buf);
	      v = image->dim.x * image->dim.y * image->dim.z;
	      for ( i = 0; i < v; i ++, ubuf ++, sbuf ++) {
		if ( (*ubuf) == (u8)0 ) *sbuf = (u8)0;
	      }
	      /*--- reconstruction ---*/
	      if ( par.dim == VT_2D ) {
		switch( par.type_recons ) {
		case _GREY_ :
		  if ( VT_GreyReconstruct2D( &res ) != 1 ) {
		    VT_FreeImage( image );
		    VT_Free( (void**)&image );
		    VT_FreeRes( &res, par.dim );
		    VT_ErrorParse("unable to (grey) reconstruct image\n", 0);
		  }
		  break;
		case _BIN_ :
		default :
		  if ( VT_Reconstruct2D( &(res.imext), &(res.imscale), &(par.par_line) ) != 1 ) {
		    VT_FreeImage( image );
			VT_Free( (void**)&image );
			VT_FreeRes( &res, par.dim );
			VT_ErrorParse("unable to reconstruct image\n", 0);
		  }
		}
	      } else {
		switch( par.type_recons ) {
		case _SCALE_ :
		  if ( VT_ScaleReconstruct3D( &res, &(par.par_line) ) != 1 ) {
		    VT_FreeImage( image );
		    VT_Free( (void**)&image );
		    VT_FreeRes( &res, par.dim );
		    VT_ErrorParse("unable to (scale) reconstruct image\n", 0);
		  }
		  break;
		case _CSTE_ :
		  if ( VT_CsteReconstruct3D( &res, (double)(par.rayon_recons) ) != 1 ) {
		    VT_FreeImage( image );
		    VT_Free( (void**)&image );
		    VT_FreeRes( &res, par.dim );
		    VT_ErrorParse("unable to (cste) reconstruct image\n", 0);
		  }
		  break;
		case _GREY_ :
		  if ( VT_GreyReconstruct3D( &res ) != 1 ) {
		    VT_FreeImage( image );
		    VT_Free( (void**)&image );
		    VT_FreeRes( &res, par.dim );
		    VT_ErrorParse("unable to (grey) reconstruct image\n", 0);
		  }
		  break;
		case _BIN_ :
		default :
		  if ( VT_Reconstruct3D( &(res.imext), &(res.imscale), &(par.par_line) ) != 1 ) {
		    VT_FreeImage( image );
		    VT_Free( (void**)&image );
		    VT_FreeRes( &res, par.dim );
		    VT_ErrorParse("unable to reconstruct image\n", 0);
		  }
		}
	      }
	    }
	    /*--- ecriture resultat ---*/
	    strcpy( res.imext.name, par.names.out );
	    if ( VT_WriteInrimage( &(res.imext) ) == -1 ) {
                VT_FreeImage( image );
                VT_Free( (void**)&image );
		VT_FreeRes( &res, par.dim );
                VT_ErrorParse("unable to write output image\n", 0);
	    }
	    break;
	case _EXTREMA2_ :
	case _EXTREMA_ :
	default :
	    fbuf = (r32*)(res.imres.buf);
	    ubuf = (u8*)(res.imext.buf);
	    v = image->dim.x * image->dim.y * image->dim.z;
	    for ( i = 0; i < v; i ++, ubuf ++, fbuf ++) {
		if ( (*ubuf) == (u8)0 ) *fbuf = 0.0;
	    }
	    strcpy( res.imres.name, par.names.out );
	    if ( VT_WriteInrimage( &(res.imres) ) == -1 ) {
                VT_FreeImage( image );
                VT_Free( (void**)&image );
		VT_FreeRes( &res, par.dim );
                VT_ErrorParse("unable to write output image\n", 0);
	    }
	}
	/*--- liberations memoires ---*/
	VT_FreeImage( image );
	VT_Free( (void**)&image );
        VT_FreeRes( &res, par.dim );
	return( 1 );
}





static void VT_Parse( int argc, char *argv[], local_par *par )
{
    int i, nb, status;
    char text[STRINGLENGTH];
    
    if ( VT_CopyName( program, argv[0] ) != 1 )
	VT_Error("Error while copying program name", (char*)NULL);
    if ( argc == 1 ) VT_ErrorParse("\n", 0 );
    
    /*--- initialisation des parametres ---*/
    VT_InitParam( par );

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
	    /*--- parametres ---*/
	    else if ( strcmp ( argv[i], "-first" ) == 0 ) {
		i += 1;
                if ( i >= argc)    VT_ErrorParse( "parsing -first...\n", 0 );
                status = sscanf( argv[i],"%f",&(par->par_line.first_coeff) );
                if ( status <= 0 ) VT_ErrorParse( "parsing -first...\n", 0 );
	    }
	    else if ( strcmp ( argv[i], "-last" ) == 0 ) {
		i += 1;
                if ( i >= argc)    VT_ErrorParse( "parsing -last...\n", 0 );
                status = sscanf( argv[i],"%f",&(par->par_line.last_coeff) );
                if ( status <= 0 ) VT_ErrorParse( "parsing -last...\n", 0 );
	    }
	    else if ( strcmp ( argv[i], "-nb" ) == 0 ) {
		i += 1;
                if ( i >= argc)    VT_ErrorParse( "parsing -nb...\n", 0 );
                status = sscanf( argv[i],"%d",&(par->par_line.nb_coeff) );
                if ( status <= 0 ) VT_ErrorParse( "parsing -nb...\n", 0 );
	    }
	    /*--- coefficient du seuillage ---*/
	    else if ( strcmp ( argv[i], "-csb" ) == 0 ) {
                i += 1;
                if ( i >= argc)    VT_ErrorParse( "parsing -csb...\n", 0 );
                status = sscanf( argv[i],"%f",&(par->mul_max_sb) );
                if ( status <= 0 ) VT_ErrorParse( "parsing -csb...\n", 0 );
            }
	    else if ( strcmp ( argv[i], "-csh" ) == 0 ) {
                i += 1;
                if ( i >= argc)    VT_ErrorParse( "parsing -csh...\n", 0 );
                status = sscanf( argv[i],"%f",&(par->mul_max_sh) );
                if ( status <= 0 ) VT_ErrorParse( "parsing -csh...\n", 0 );
            } 
	    else if ( strcmp ( argv[i], "-tcc" ) == 0 ) {
                i += 1;
                if ( i >= argc)    VT_ErrorParse( "parsing -tcc...\n", 0 );
                status = sscanf( argv[i],"%d",&(par->size_cc) );
                if ( status <= 0 ) VT_ErrorParse( "parsing -tcc...\n", 0 );
            } 
	    /*--- type des structures ---*/
	    else if ( strcmp ( argv[i], "-white" ) == 0 ) {
		par->par_line.type_structures =  VT_WHITE;
	    }
	    else if ( strcmp ( argv[i], "-black" ) == 0 ) {
		par->par_line.type_structures =  VT_BLACK;
	    }
	    /*--- type de sortie ---*/
	    else if ( strcmp ( argv[i], "-extrema" ) == 0 ) {
		par->type_output = _EXTREMA_;
	    }
	    else if ( strcmp ( argv[i], "-extrema2" ) == 0 ) {
		par->type_output = _EXTREMA2_;
	    }
	    else if ( strcmp ( argv[i], "-reponse" ) == 0 ) {
		par->type_output = _REPONSE_;
	    }
	    else if ( strcmp ( argv[i], "-echelle" ) == 0 ) {
		par->type_output = _ECHELLE_;
	    }
	    else if ( strcmp ( argv[i], "-seuillage" ) == 0 ) {
		par->type_output = _SEUILLAGE_;
	    }
	    else if ( strcmp ( argv[i], "-recons" ) == 0 ) {
		par->type_output = _RECONS_;
	    }
	    else if ( strcmp ( argv[i], "-scale" ) == 0 ) {
		par->type_recons = _SCALE_;
	    }
	    else if ( strcmp ( argv[i], "-rad" ) == 0 ) {
                i += 1;
                if ( i >= argc)    VT_ErrorParse( "parsing -rad...\n", 0 );
                status = sscanf( argv[i],"%f",&(par->rayon_recons) );
                if ( status <= 0 ) VT_ErrorParse( "parsing -rad...\n", 0 );
		par->type_recons = _CSTE_;
            } 
	    else if ( strcmp ( argv[i], "-grey" ) == 0 ) {
		par->type_recons = _GREY_;
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
    par->dim = VT_3D;
    par->mul_max_sh = 0.33333333333333333333; /* 1/3 */
    par->mul_max_sb = 0.16666666666666666666; /* 1/6 */
    par->size_cc = 1;
    par->type_output = _EXTREMA_;
    par->type_recons = _BIN_;
    par->rayon_recons = 1.0;
    VT_Line( &(par->par_line) );
    VT_Names( &(par->names) );
}
