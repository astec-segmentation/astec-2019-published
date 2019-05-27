/*************************************************************************
 * progressiveThreshold.c -
 *
 * $Id: progressiveThreshold.c,v 1.0 2017/10/27 10:22:00 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2017/10/27
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>
#include <mt_anisotropicHist.h>

#include <chunks.h>


typedef struct local_par {

  vt_names names;
  int wi;
  double x_beg;
  double y_beg;
  double z_beg;
  double x_end;
  double y_end;
  double z_end;

} local_par;



/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );





static char *usage = "[image-in] [image-out] [-th-begin[-x|y|z] %lf] [-th-end[-x|y|z] %lf] \n\
\t [-wi] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -th-begin %lf: valeur du threshold pour z=0\n\
\t -th-end %lf: valeur du threshold pour z=end\n\
\t -th-[begin|end]-[x|y|z] %lf: definit un threshold axial en [x|y|z] pour z=0 (begin) ou z=end (end)\n\
\t (Note: le threshold en un z quelconque est calcule en tant que ratio entre le threshold a z=0 et celui a z=end)\n\
\t -wi : ecrit les images reponses projetees sur X, Y, Z + eventuelles binarisations\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";


static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_Rep_Angles imageIn;
  VT_InitRep_Angles(&imageIn);

  vt_image *imRep=NULL, *imTht=NULL, *imPhi=NULL;
  vt_image imRepAnisotropes[3];
  vt_image imBin;
  char name[DOUBLESTRINGLENGTH], prefix[STRINGLENGTH];
  int i;

  /*--- initialisation des parametres ---*/
  MT_InitParam( &par );

  /*--- lecture des parametres ---*/
  MT_Parse( argc, argv, &par );

  /*--- lecture des images d'entree ---*/

  if(1)
  {
    sprintf( prefix, "%s", par.names.in );
    i=0;
    while (prefix[i]!='\0')
    {
      i++;
    }
    if (i>4 && prefix[i-4]=='.' &&
        ((prefix[i-3]=='h' && prefix[i-2]=='d' && prefix[i-1]=='r') ||
         (prefix[i-3]=='i' && prefix[i-2]=='n' && prefix[i-1]=='r') ||
         (prefix[i-3]=='m' && prefix[i-2]=='h' && prefix[i-1]=='a') ))
    {
      prefix[i-4]='\0';
      sprintf( name, "%s", par.names.in );
    }
    else
    {
      sprintf( name, "%s.hdr", prefix );
    }

    imRep = _VT_Inrimage( name );
    if ( imRep == (vt_image*)NULL )
    {
      /* format .hdr */
      i=0;
      while (name[i]!='\0')
      {
        i++;
      }
      name[i++] = '.';
      name[i++] = 'h';
      name[i++] = 'd';
      name[i++] = 'r';
      name[i]  = '\0';

      imRep = _VT_Inrimage( name );
      if ( imRep == (vt_image*)NULL )
      {
        /* format .inr */
        name[i-3] = 'i';
        name[i-2] = 'n';

        imRep = _VT_Inrimage( name );
        if ( imRep == (vt_image*)NULL )
        {
          fprintf(stderr, "%s unreadable\n", par.names.in );
          MT_ErrorParse("unable to read input imageA\n", 0);
        }
      }
    }
    if (_VT_VERBOSE_)
      fprintf(stdout, "Image d'entree : %s\n", name);
    imageIn.imrep = imRep;

    sprintf( name, "%s.theta.hdr", prefix );
    imTht = _VT_Inrimage( name );
    if ( imTht == (vt_image*)NULL )
    {
      i=0;
      while (name[i]!='\0')
      {
        i++;
      }
      name[i-3] = 'i';
      name[i-2] = 'n';
      imTht = _VT_Inrimage( name );
      if ( imTht == (vt_image*)NULL )
      {
        int j=0;
        i=0;
        while (prefix[i]!='\0')
        {
          if (prefix[i]=='.') j = i;
          i++;
        }
        if (j>0)
        {
          prefix[j]='\0';
          sprintf( name, "%s.theta.inr", prefix );

          imTht = _VT_Inrimage( name );
          if ( imTht == (vt_image*)NULL )
          {
            i=0;
            while (name[i]!='\0')
            {
              i++;
            }
            name[i-3] = 'h';
            name[i-2] = 'd';
            imTht = _VT_Inrimage( name );
            if ( imTht == (vt_image*)NULL )
            {
                VT_FreeImage( imRep );
                fprintf(stderr, "%s unreadable\n", name);
                MT_ErrorParse("unable to read input imageC\n", 0);
            }
          }
        }
      }
    }
    if (_VT_VERBOSE_)
      fprintf(stdout, "Image d'entree theta : %s\n", name);
    imageIn.imtheta = imTht;

    sprintf( name, "%s.phi.inr", prefix );
    imPhi = _VT_Inrimage( name );
    if ( imPhi == (vt_image*)NULL )
    {
      i=0;
      while (name[i]!='\0')
      {
        i++;
      }
      name[i-3] = 'h';
      name[i-2] = 'd';
      imPhi = _VT_Inrimage( name );
      if ( imPhi == (vt_image*)NULL )
      {
        VT_FreeImage( imTht );
        VT_FreeImage( imRep );
        fprintf(stderr, "%s unreadable\n", name);
        MT_ErrorParse("unable to read input imageD\n", 0);
      }
    }
    if (_VT_VERBOSE_)
      fprintf(stdout, "Image d'entree phi : %s\n", name);
    imageIn.imphi = imPhi;

  }



  /*--- calculs ---*/

  if(MT_ComputeRepAnisotropic(imageIn, imRepAnisotropes) == 0) {
      VT_FreeRep_Angles ( &imageIn );
      VT_FreeImage(imRepAnisotropes);
      VT_FreeImage(imRepAnisotropes+1);
      VT_FreeImage(imRepAnisotropes+2);
      MT_ErrorParse("unable to compute the axial contributions\n", 0);
  }

  if(par.wi==1)
  {
    if(_VT_VERBOSE_)
      fprintf(stdout, "Ecriture des projections sur X, Y, Z...\n");

    VT_WriteInrimage( imRepAnisotropes);
    VT_WriteInrimage( imRepAnisotropes+1);
    VT_WriteInrimage( imRepAnisotropes+2);
  }


  double Temp[1];
  if (MT_BinariseProgressive(imageIn, imRepAnisotropes,
      &imBin, par.names.out, par.x_beg, par.y_beg, par.z_beg, par.x_end, par.y_end, par.z_end,
      _BIN_DOT2_, Temp) < 0)
  {
    VT_FreeRep_Angles ( &imageIn );
    VT_FreeImage(imRepAnisotropes);
    VT_FreeImage(imRepAnisotropes+1);
    VT_FreeImage(imRepAnisotropes+2);
    MT_ErrorParse("unable to binarise the image\n", 0);
  }

  if(_VT_VERBOSE_)
    fprintf(stdout, "Ecriture de %s ...\n", par.names.out);
  VT_WriteInrimage( &imBin);



  /*--- liberations memoires ---*/
  VT_FreeRep_Angles ( &imageIn );
  VT_FreeImage(imRepAnisotropes);
  VT_FreeImage(imRepAnisotropes+1);
  VT_FreeImage(imRepAnisotropes+2);
  VT_FreeImage(&imBin);

  return( 0 );
}








static void MT_Parse( int argc,
                      char *argv[],
                      local_par *par )
{
  int i, nb, status;
  int maxchunks;
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
        par->wi = 1;
      }


      /* Parametres de calcul */

      else if ( strcmp ( argv[i], "-th-begin" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -th-begin...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->x_beg) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -th-begin...\n", 0 );
        par->y_beg=par->x_beg;
        par->z_beg=par->x_beg;
      }

      else if ( strcmp ( argv[i], "-th-end" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -th-end...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->x_end) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -th-end...\n", 0 );
        par->y_end=par->x_end;
        par->z_end=par->x_end;
      }

      else if ( strcmp ( argv[i], "-th-begin-x" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -th-begin-x...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->x_beg) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -th-begin-x...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-th-begin-y" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -th-begin-y...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->y_beg) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -th-begin-y...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-th-begin-z" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -th-begin-z...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->z_beg) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -th-begin-z...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-th-end-x" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -th-end-x...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->x_end) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -th-end-x...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-th-end-y" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -th-end-y...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->y_end) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -th-end-y...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-th-end-z" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -th-end-z...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->z_end) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -th-end-z...\n", 0 );
      }


      /* parallelism
       */
    else if ( strcmp ( argv[i], "-parallel" ) == 0 ) {
      setParallelism( _DEFAULT_PARALLELISM_ );
    }

    else if ( strcmp ( argv[i], "-no-parallel" ) == 0 ) {
       setParallelism( _NO_PARALLELISM_ );
    }

    else if ( strcmp ( argv[i], "-parallelism-type" ) == 0 ||
               strcmp ( argv[i], "-parallel-type" ) == 0 ) {
      i ++;
      if ( i >= argc)    MT_ErrorParse( "-parallelism-type", 0 );
      if ( strcmp ( argv[i], "default" ) == 0 ) {
        setParallelism( _DEFAULT_PARALLELISM_ );
      }
      else if ( strcmp ( argv[i], "none" ) == 0 ) {
        setParallelism( _NO_PARALLELISM_ );
      }
      else if ( strcmp ( argv[i], "openmp" ) == 0 || strcmp ( argv[i], "omp" ) == 0 ) {
        setParallelism( _OMP_PARALLELISM_ );
      }
      else if ( strcmp ( argv[i], "pthread" ) == 0 || strcmp ( argv[i], "thread" ) == 0 ) {
        setParallelism( _PTHREAD_PARALLELISM_ );
      }
    }


    else if ( strcmp ( argv[i], "-max-chunks" ) == 0 ) {
      i ++;
      if ( i >= argc)    MT_ErrorParse( "-max-chunks", 0 );
      status = sscanf( argv[i], "%d", &maxchunks );
      if ( status <= 0 ) MT_ErrorParse( "-max-chunks", 0 );
      if ( maxchunks >= 1 ) setMaxChunks( maxchunks );
    }

    else if ( strcmp ( argv[i], "-omp-scheduling" ) == 0 ||
              ( strcmp ( argv[i], "-omps" ) == 0 && argv[i][3] == '\0') ) {
      i ++;
      if ( i >= argc)    MT_ErrorParse( "-parallel-scheduling", 0 );
      if ( strcmp ( argv[i], "default" ) == 0 ) {
        setOmpScheduling( _DEFAULT_OMP_SCHEDULING_ );
      }
      else if ( strcmp ( argv[i], "static" ) == 0 ) {
        setOmpScheduling( _STATIC_OMP_SCHEDULING_ );
      }
      else if ( strcmp ( argv[i], "dynamic-one" ) == 0 ) {
        setOmpScheduling( _DYNAMIC_ONE_OMP_SCHEDULING_ );
      }
      else if ( strcmp ( argv[i], "dynamic" ) == 0 ) {
        setOmpScheduling( _DYNAMIC_OMP_SCHEDULING_ );
      }
      else if ( strcmp ( argv[i], "guided" ) == 0 ) {
        setOmpScheduling( _GUIDED_OMP_SCHEDULING_ );
      }
      else {
        fprintf( stderr, "unknown omp scheduling type: '%s'\n", argv[i] );
        MT_ErrorParse( "-omp-scheduling", 0 );
      }
    }


      /*--- option inconnue ---*/
      else {
        sprintf(text,"unknown option %s\n",argv[i]);
        MT_ErrorParse(text, 0);
      }
    }
    /*--- saisie des noms d'image in / fichier out ---*/
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
        MT_ErrorParse("too much file names when parsing\n", 0 );
    }
    i += 1;
  }

  /*--- s'il n'y a pas assez de noms ... ---*/
  if (nb == 0) {
    strcpy( par->names.in,  "<" );  /* standart input */
    strcpy( par->names.out, ">" );  /* standart output */
  }
  if (nb == 1)
  {
    strcpy( par->names.out, ">" );  /* standart output */
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

  par->wi = 0;

  par->x_beg = 0;
  par->y_beg = 0;
  par->z_beg = 0;
  par->x_end = 0;
  par->y_end = 0;
  par->z_end = 0;
}

