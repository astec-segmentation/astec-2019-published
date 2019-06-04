/*************************************************************************
 * mc-removeLinec -
 *
 * $$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Jeu  5 dec 2013 18:53:48 CET
 *
 * ADDITIONS, CHANGES
 *
 */

#include <sys/time.h> /* gettimeofday() */
#include <time.h> /* clock() */

#include <chunks.h>

#include <vt_common.h>

#include <vt_removeLine.h>





static int _time_ = 1;
static int _clock_ = 1;

typedef struct local_par {
  vt_names names;

  char *input_correction_name;
  char *output_correction_name;
  
  typeRemoveLineParameter p;

  int type;
} local_par;




/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
static char *_BaseName( char *p );
static double VT_GetTime();
static double VT_GetClock();




static char *usage = "[image-in] [image-out]\n\
 [-output-corrections|-oc %s]\n\
 [-input-corrections|-ic %s]\n\
 [-method|-m global|g | regional|r | local|l]\n\
 [-contrastFraction|-c %f]\n\
 [-yRejection|-y %f]\n\
 [-xzKept|-xz %f]\n\
 [-automated]\n\
 [-inv] [-swap] [-v] [-D]\n\
 [-parallel|-no-parallel]\n\
 [-help] [options-de-type]";

static char *detail = "\
-output-corrections|-oc %s: save the computed corrections\n\
-input-corrections|-ic %s: does not compute the corrections\n\
  but applies the read corrections\n\
-method %s:\n\
  local|l: local analysis\n\
       # XZ slices selected (see below) are corrected w.r.t.\n\
       # neighboring slices\n\
  regional|r: regional analysis\n\
       # increase or decrease of intensity are detected/selected\n\
       # (through the coefficient describes below)\n\
       # so that plateaus of higher intensity can be identified\n\
       # if a negative value is given as coefficient,\n\
       # an 'optimal' value of coefficient is searched\n\
  global|g: global analysis\n\
-contrastFraction %lf # ratio in ]0,1[ = percentage / 100 \n\
       # a XZ slice is selected for removal if the percentage of points\n\
       # with a value larger than its neighbors in the upper (or the lower)\n\
       # adjacent XZ slice is larger than the given percentage.\n\
       # The lower the ratio, the larger the number of slices to be corrected.\n\
       # This value has then to be larger than 0.5 (e.g. 0.6)\n\
 si 'image-in' est '-', on prendra stdin\n\
 si 'image-out' est absent, on prendra stdout\n\
 si les deux sont absents, on prendra stdin et stdout\n\
 -inv : inverse 'image-in'\n\
 -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
 -v : mode verbose\n\
 -D : mode debug\n\
 options-de-type : -o 1    : unsigned char\n\
                   -o 2    : unsigned short int\n\
                   -o 2 -s : short int\n\
                   -o 4 -s : int\n\
                   -r      : float\n\
 si aucune de ces options n'est presente, on prend le type de 'image-in'\n\
\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image, imres;
  typeCorrectionList theList;
  double time_init, time_exit;
  double clock_init, clock_exit;

  time_init = VT_GetTime();
  clock_init = VT_GetClock();


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
  VT_Image( &imres );
  VT_InitFromImage( &imres, image, par.names.out, image->type );
   if ( VT_AllocImage( &imres ) != 1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }
  
  VT_InitCorrectionList( &theList );


  /*---  ---*/

  if ( par.input_correction_name != (char*)NULL ) {
      if ( VT_ReadCorrectionList( par.input_correction_name, &theList) != 1 ) {
          VT_FreeImage( &imres );
          VT_FreeImage( image );
          VT_Free( (void**)&image );
          VT_FreeCorrectionList( &theList );
          VT_ErrorParse( "unable to read input corrections\n", 0 );
      }
      if ( VT_CorrectLines( image, &imres, &theList ) != 1 ) {
          VT_FreeImage( &imres );
          VT_FreeImage( image );
          VT_Free( (void**)&image );
          VT_FreeCorrectionList( &theList );
          VT_ErrorParse( "unable to apply input corrections\n", 0 );
      }
  }
  else {
      if ( VT_RemoveLines( image, &imres, &theList, &(par.p) ) != 0 ) {
          VT_FreeImage( &imres );
          VT_FreeImage( image );
          VT_Free( (void**)&image );
          VT_FreeCorrectionList( &theList );
          VT_ErrorParse( "unable to remove line(s)\n", 0 );
      }
  }

 VT_FreeImage( image );
 VT_Free( (void**)&image );

 /*--- ecriture de l'image resultat ---*/
 if ( VT_WriteInrimage( &imres ) == -1 ) {
   VT_FreeImage( &imres );
   VT_FreeCorrectionList( &theList );
   VT_ErrorParse("unable to write output image\n", 0);
 }

 VT_FreeImage( &imres );

 if ( par.output_correction_name != (char*)NULL ) {
     if ( VT_WriteCorrectionList( par.output_correction_name, &theList ) != 1 ) {
         VT_FreeCorrectionList( &theList );
         VT_ErrorParse("unable to write output correction list\n", 0);
     }
 }
  
 VT_FreeCorrectionList( &theList );

  time_exit = VT_GetTime();
  clock_exit = VT_GetClock();

  if ( _time_ || _clock_ ) {
    fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName( argv[0] ), time_exit - time_init );
    fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
    fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );
  }

  return( 0 );
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

      else if ( strcmp ( argv[i], "-output-corrections" ) == 0
                || (strcmp ( argv[i], "-oc" ) == 0 && argv[i][3] == '\0') ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -output-corrections ...\n", 0 );
        par->output_correction_name = argv[i];
      }

      else if ( strcmp ( argv[i], "-input-corrections" ) == 0
                || (strcmp ( argv[i], "-ic" ) == 0 && argv[i][3] == '\0') ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -input-corrections ...\n", 0 );
        par->input_correction_name = argv[i];
      }

      else if ( strcmp ( argv[i], "-method" ) == 0
                || (strcmp ( argv[i], "-m" ) == 0 && argv[i][2] == '\0') ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -method ...\n", 0 );
        if ( strcmp ( argv[i], "global" ) == 0
             || (strcmp ( argv[i], "g" ) == 0 && argv[i][1] == '\0') )
            par->p.method = _GLOBAL_;
        else if ( strcmp ( argv[i], "regional" ) == 0
             || (strcmp ( argv[i], "r" ) == 0 && argv[i][1] == '\0') )
            par->p.method = _REGIONAL_;
        else if ( strcmp ( argv[i], "local" ) == 0
                  || (strcmp ( argv[i], "l" ) == 0 && argv[i][1] == '\0') )
                 par->p.method = _LOCAL_;
        else
            VT_ErrorParse( "parsing -method: unknown method ...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-contrastFraction" ) == 0
                || (strcmp ( argv[i], "-c" ) == 0 && argv[i][2] == '\0') ) {
	i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -contrastFraction ...\n", 0 );
        status = sscanf( argv[i],"%f", &(par->p.contrastSignificantFraction) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -contrastFraction ...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-yRejection" ) == 0
                || (strcmp ( argv[i], "-y" ) == 0 && argv[i][2] == '\0') ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -yRejection ...\n", 0 );
        status = sscanf( argv[i],"%f", &(par->p.yRejectedFraction) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -yRejection ...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-xzKept" ) == 0
                || (strcmp ( argv[i], "-xz" ) == 0 && argv[i][3] == '\0') ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -xzKept ...\n", 0 );
        status = sscanf( argv[i],"%f", &(par->p.xzKeptFraction) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -xzKept ...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-automated" ) == 0 ) {
        par->p.automatedChoices = 1;
      }

      /*--- arguments generaux ---*/
      else if ( strcmp ( argv[i], "-help" ) == 0 ) {
	VT_ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 ) {
	_VT_VERBOSE_ = 1;
        VT_IncrementVerboseInVtRemoveLine();
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 ) {
	_VT_DEBUG_ = 1;
      }

      else if ( strcmp ( argv[i], "-parallel" ) == 0 ) {
         setParallelism( _DEFAULT_PARALLELISM_ );
      }

      else if ( strcmp ( argv[i], "-no-parallel" ) == 0 ) {
         setParallelism( _NO_PARALLELISM_ );
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

  par->input_correction_name = (char*)NULL;
  par->output_correction_name = (char*)NULL;

  initRemoveLineParameter( &(par->p) );

  par->type = TYPE_UNKNOWN;
}



static char *_BaseName( char *p )
{
  int l;
  if ( p == (char*)NULL ) return( (char*)NULL );
  l = strlen( p ) - 1;
  while ( l >= 0 && p[l] != '/' ) l--;
  if ( l < 0 ) l = 0;
  if ( p[l] == '/' ) l++;
  return( &(p[l]) );
}



static double VT_GetTime() 
{
  struct timeval tv;
  gettimeofday(&tv, (void *)0);
  return ( (double) tv.tv_sec + tv.tv_usec*1e-6 );
}

static double VT_GetClock() 
{
  return ( (double) clock() / (double)CLOCKS_PER_SEC );
}
