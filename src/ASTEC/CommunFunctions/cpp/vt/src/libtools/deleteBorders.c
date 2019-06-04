/*************************************************************************
 * membrane.c -
 *
 * $Id: membrane.c,v 1.0 2013/08/12 10:34:51 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2013/08/12
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>

#define TABLENGTH 2000

typedef struct local_par {
  vt_names names;
  vt_names nomgraines;
  int areGraines;
} local_par;



/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );

static int _verbose_ = 0;





static char *usage = "[image-in] [image-out] [-seeds|-graines %s %s]\n\
\t [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t [-graines|-seeds] %s %s : passe a zero les graines eliminees (param : nom graines + nouveau nom)\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];







int main( int argc, char *argv[] )
{

  local_par par;
  vt_image *image;
  vt_image *graines;
  vt_image imres;
  vt_image Graines;
  /*  int flag_3D = 1; */
  char name[DOUBLESTRINGLENGTH];
  sprintf( name, "%s", par.names.out );


  /*--- initialisation des parametres ---*/
  MT_InitParam( &par );

  /*--- lecture des parametres ---*/
  MT_Parse( argc, argv, &par );

  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL )
    MT_ErrorParse("unable to read input image\n", 0);

  if (par.areGraines==1)
  {
    graines = _VT_Inrimage( par.nomgraines.in );
    if ( graines == (vt_image*)NULL )
    {
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      MT_ErrorParse("unable to read input image\n", 0);
    }
  }

  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.names.inv == 1 )  VT_InverseImage( image );
  if ( par.names.swap == 1 ) VT_SwapImage( image );


  /* if ( image->dim.z == 1 )
   *  flag_3D = 0;
   */



  /* Alloc image sortie */

  int t = (int) image->type;
  sprintf( name, "%s", par.names.out );
  VT_InitImage( &imres, name, image->dim.x, image->dim.y,
              image->dim.z, t );
  if ( VT_AllocImage( &(imres) ) != 1 ) {
    MT_ErrorParse("problem while allocating imres\n", 0 );
    VT_FreeImage( image );
    VT_FreeImage( graines );
    return( -1 );
  }

  imres.siz.x=image->siz.x;
  imres.siz.y=image->siz.y;
  imres.siz.z=image->siz.z;

  if(par.areGraines==1)
  {
    sprintf( name, "%s", par.nomgraines.out );
    VT_InitImage( &Graines, name, graines->dim.x, graines->dim.y,
              graines->dim.z, t );
    if ( VT_AllocImage( &(Graines) ) != 1 ) {
      MT_ErrorParse("problem while allocating graines_res\n", 0 );
      VT_FreeImage( image );
      VT_FreeImage( graines );
      VT_FreeImage( &imres );
      return( -1 );
    }

    Graines.siz.x=graines->siz.x;
    Graines.siz.y=graines->siz.y;
    Graines.siz.z=graines->siz.z;
  }

  /*--- calculs ---*/

  unsigned short int ***imgarray=NULL;
  unsigned short int ***resarray=NULL;
  unsigned short int ***garray=NULL;
  unsigned short int ***Garray=NULL;

  unsigned char ***imgarrayu8=NULL;
  unsigned char ***resarrayu8=NULL;
  unsigned char ***garrayu8=NULL;
  unsigned char ***Garrayu8=NULL;

  switch (image->type) {
    case SCHAR:
    case UCHAR:
      imgarrayu8= (unsigned char ***)image->array;
      resarrayu8= (unsigned char ***)imres.array;
      garrayu8= (unsigned char ***)graines->array;
      Garrayu8= (unsigned char ***)Graines.array;
      break;
    case SSHORT:
    case USHORT:
      imgarray= (unsigned short int ***)image->array;
      resarray= (unsigned short int ***)imres.array;
      garray= (unsigned short int ***)graines->array;
      Garray= (unsigned short int ***)Graines.array;
      break;
    case TYPE_UNKNOWN:
    default:
      VT_Error("image type unknown or not supported for this program",program);
      return( 0 );
  }

  int z,y,x;
  unsigned char toDeleteu8[TABLENGTH];
  unsigned short int toDelete[TABLENGTH];
  int cpt=0, c;
  unsigned char valu8;
  unsigned short int val;
  for(y=0; y<(int)image->dim.y; y++)
  {
    for(x=0; x<(int)image->dim.x; x++)
    {
        switch ( image->type ) {
          case SCHAR :
          case UCHAR :
            valu8 = (imgarrayu8[0][y][x]); /* z = zmin */
            for(c=0;c<cpt;c++)
              if (toDeleteu8[c] == valu8 )
              {
                c = -1;
                break;
              }
            if (c == cpt)
              toDeleteu8[cpt++]=valu8;

            valu8 = (imgarrayu8[image->dim.z-1][y][x]); /* z = zmax */
            for(c=0;c<cpt;c++)
              if (toDeleteu8[c] == valu8)
              {
                c = -1;
                break;
              }
            if (c == cpt)
              toDeleteu8[cpt++]=valu8;
            break;
          case SSHORT :
          case USHORT :
            val = (imgarray[0][y][x]);      /* z = zmin */
            for(c=0;c<cpt;c++)
              if (toDelete[c] == val )
              {
                c = -1;
                break;
              }
            if (c == cpt)
              toDelete[cpt++]=val;

            val = (imgarray[image->dim.z-1][y][x]); /* z = zmax */
            for(c=0;c<cpt;c++)
              if (toDelete[c] == val)
              {
                c = -1;
                break;
              }
            if (c == cpt)
              toDelete[cpt++]=val;
            break;
          default:
            break;
        }
    }
  }

  for(z=1;z<(int)image->dim.z-1;z++)
    for(x=0;x<(int)image->dim.x;x++)
    {
        switch ( image->type ) {
          case SCHAR :
          case UCHAR :
            valu8 = imgarrayu8[z][0][x];      /* y = ymin */
            for(c=0;c<cpt;c++)
              if (toDeleteu8[c] == valu8 )
              {
                c = -1;
                break;
              }
            if (c == cpt)
              toDeleteu8[cpt++]=valu8;

            valu8 = imgarrayu8[z][image->dim.y-1][x]; /* y = ymax */
            for(c=0;c<cpt;c++)
              if (toDeleteu8[c] == valu8)
              {
                c = -1;
                break;
              }
            if (c == cpt)
              toDeleteu8[cpt++]=valu8;
            break;
          case SSHORT :
          case USHORT :

            val = imgarray[z][0][x];      /* y = ymin */
            for(c=0;c<cpt;c++)
              if (toDelete[c] == val )
              {
                c = -1;
                break;
              }
            if (c == cpt)
              toDelete[cpt++]=val;

            val = imgarray[z][image->dim.y-1][x]; /* y = ymax */
            for(c=0;c<cpt;c++)
              if (toDelete[c] == val)
              {
                c = -1;
                break;
              }
            if (c == cpt)
              toDelete[cpt++]=val;
            break;
          default:
            break;
        }
    }
  for(z=1;z<(int)image->dim.z-1;z++)
    for(y=1;y<(int)image->dim.y-1;y++)
    {
        switch ( image->type ) {
          case SCHAR :
          case UCHAR :
            valu8 = imgarrayu8[z][y][0];      /* x = xmin */
             for(c=0;c<cpt;c++)
               if (toDeleteu8[c] == valu8 )
               {
                 c = -1;
                 break;
               }
             if (c == cpt)
               toDeleteu8[cpt++]=valu8;

             valu8 = imgarrayu8[z][y][image->dim.x-1]; /* x = xmax */
             for(c=0;c<cpt;c++)
               if (toDeleteu8[c] == valu8)
               {
                 c = -1;
                 break;
               }
             if (c == cpt)
               toDeleteu8[cpt++]=valu8;
             break;
          case SSHORT :
          case USHORT :
            val = imgarray[z][y][0];      /* x = xmin */
             for(c=0;c<cpt;c++)
               if (toDelete[c] == val )
               {
                 c = -1;
                 break;
               }
             if (c == cpt)
               toDelete[cpt++]=val;

             val = imgarray[z][y][image->dim.x-1]; /* x = xmax */
             for(c=0;c<cpt;c++)
               if (toDelete[c] == val)
               {
                 c = -1;
                 break;
               }
             if (c == cpt)
               toDelete[cpt++]=val;
             break;
          default:
            break;
        }
    }

  if (_verbose_)
  {
    fprintf(stdout, "cpt = %d\n", cpt);

    fprintf(stdout, "toDelete = [");
    for (c = 0;c< cpt; c++)
      switch (image->type) {
        case SCHAR:
        case UCHAR:
          fprintf(stdout, "%d, ", (int)toDeleteu8[c]);
          break;
        case SSHORT:
        case USHORT:
          fprintf(stdout, "%d, ", (int)toDelete[c]);
          break;
        default:
          break;
      }
    fprintf(stdout, "]\n");
  }

  for (z = 0; z<(int)image->dim.z; z++)
    for (y = 0; y<(int)image->dim.y; y++)
      for (x = 0; x<(int)image->dim.x; x++)
      {
        switch ( image->type ) {
          case SCHAR :
          case UCHAR :
            valu8 = imgarrayu8[z][y][x];
            for(c = 0; c < cpt; c++)
              if(valu8 == toDeleteu8[c])
                break;
            if (c == cpt)
            {
              resarrayu8[z][y][x]=valu8;
              if (par.areGraines == 1)
                Garrayu8[z][y][x]=garrayu8[z][y][x];
            }
            else
            {
              resarrayu8[z][y][x]='\0';
              if (par.areGraines == 1)
                Garrayu8[z][y][x]='\0';
            }
            break;
          case SSHORT :
          case USHORT :
            val = imgarray[z][y][x];
            for(c = 0; c < cpt; c++)
              if(val == toDelete[c])
                break;
            if (c == cpt)
            {
              resarray[z][y][x]=val;
              if (par.areGraines == 1)
                Garray[z][y][x]=garray[z][y][x];
            }
            else
            {
              resarray[z][y][x]=0;
              if (par.areGraines == 1)
                Garray[z][y][x]=0;
            }
            break;
          default:
            break;
        }
      }

  VT_WriteInrimage( &(imres) );
  if(par.areGraines==1)
    VT_WriteInrimage( &(Graines) );


  /*--- liberations memoires ---*/
  VT_FreeImage( image );
  VT_Free( (void**)&image );
  VT_FreeImage( &imres );
  if(par.areGraines==1)
  {
    VT_FreeImage( graines );
    VT_Free( (void**)&graines );
    VT_FreeImage( &Graines );
  }

  return( cpt ); /* returns the number of deleted labels */
}








static void MT_Parse( int argc,
                      char *argv[],
                      local_par *par )
{
  int i, nb;
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

      /*--- traitement eventuel de l'image d'entree ---*/
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
        par->names.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
        par->names.swap = 1;
      }


      /* Parametres de calcul */

      else if ( strcmp ( argv[i], "-seeds" ) == 0  ||
          strcmp ( argv[i], "-graines" ) == 0 ||
          strcmp ( argv[i], "-seed" ) == 0 ||
          strcmp ( argv[i], "-graine" ) == 0 ) {
        par->areGraines = 1;
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -graines...\n", 0 );
        strncpy( par->nomgraines.in, argv[i], STRINGLENGTH );
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -graines...\n", 0 );
        strncpy( par->nomgraines.out, argv[i], STRINGLENGTH );
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
  VT_Names( &(par->nomgraines) );
  par->areGraines = 0;
}
