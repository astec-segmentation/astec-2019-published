/*************************************************************************
 * fuseLabels.c -
 *
 * $Id: fuseLabels.c,v 1.0 2014/09/15 17:25:51 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2014/09/15
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>

#define TABLENGTH 2000

typedef struct local_par {
  vt_names names;
  float gauches[TABLENGTH];
  float droites[TABLENGTH];
  int N;
  float set[TABLENGTH];
  int S;
  float vector[TABLENGTH][2];
  int V;
  int flagLut;
  int flagExclusive;
} local_par;



/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );


static int VT_AllocLutWithFileName( int ** lut, int *lut_size, char* file_name );
static int VT_FillLutWithFileName( int **lut, int lut_size, char* file_name, int flagExclusive);

static int _verbose_ = 0;



static int findInVec(float val, float *vec, int n);

static char *usage = "[image-in] [image-out]\n\
\t [-p |-pair %d %d [%d %d [...]]] [-s|set %d [%d...]] [-sv|setvector %d %d [...]] \n\
\t [-lut %s [-e[xclusive]]]\n\
\t[-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -p %d %d [...] : attribue la 1ere valeur a tous les pixels ayant la 2eme valeur\n\
\t -s %d [...] : met a zero l'ensemble des pixels dont la valeur est precisee\n\
\t -sv %d [...] : met a zero l'ensemble des pixels dont la valeur est entre les 2 bornes donnees\n\
\t -lut %s : transmet en parametre un fichier de look-up-table dans lequel chaque ligne respecte le format : <ancien_label> <nouveau_label> (/!\\ contraintes : uniquement des entiers positifs ; uniquement des images entree/sortie encodees sur le meme nombre de bits /!\\)\n\
\t -e : met a zero l'ensemble des pixels ayant une valeur non precisee par la lut\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{

  local_par par;
  vt_image *image;
  /* vt_image imres;
   * int flag_3D = 1;
   */
  int z,y,x;
  int i;
  unsigned short int d;
  float val;

  int *lut=(int*)NULL;
  int lut_size = 0;

  /*--- initialisation des parametres ---*/
  MT_InitParam( &par );

  /*--- lecture des parametres ---*/
  MT_Parse( argc, argv, &par );

  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL )
    MT_ErrorParse("unable to read input image\n", 0);

  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.names.inv == 1 )  VT_InverseImage( image );
  if ( par.names.swap == 1 ) VT_SwapImage( image );


  /* if ( image->dim.z == 1 )
   * flag_3D = 0;
   */

  /* Chargement de la lut */
  if (par.flagLut)
  {
      /*--- allocation de la lut ---*/
      if ( VT_AllocLutWithFileName(&lut, &lut_size, par.names.ext) != 1) {
          MT_ErrorParse("unable to allocate lut for lut file\n", 0);
          return(0);
      }
      if (lut_size == 0)
      {
          /* Aucune modification */
          VT_WriteInrimageWithName( image, par.names.out );
          VT_FreeImage( image );
          VT_Free( (void**)&image );
          return(1);
      }

      /*--- remplissage de la lut ---*/
      /*fprintf(stdout, "Filling LUT of size %d...\n", lut_size);*/
      if (VT_FillLutWithFileName(&lut, lut_size, par.names.ext, par.flagExclusive ) != 1)
      {
          VT_FreeImage(image);
          VT_Free( (void**)&image );
          if (lut != NULL)
          {
            free(lut);
            lut = NULL;
          }
          MT_ErrorParse("unable to fill lut for lut file\n", 0);
          return(0);
      }
      /*fprintf(stdout, "Done : lut_size = %d\n",lut_size);*/

  }
  /*
  if (par.flagLut)
  {
      int index;
      for (index=0 ; index<lut_size ; index ++)
      {
          fprintf(stdout, "lut[%d] = %d\n",index, lut[index]);
      }

  }
  */
  /*--- calculs ---*/

  unsigned char ***arrayu8=NULL;
  /* char ***array8; */
  unsigned short int ***arrayu16=NULL;
  /* short int ***array16; */
  float ***arrayflt=NULL;


  /*      array8=(char ***)image->array;
        if (_VT_VERBOSE_)   fprintf(stdout, "Type = char \n");
        break;
  */
  /*      array16=(short int ***)image->array;
        if (_VT_VERBOSE_) fprintf(stdout, "Type = short int \n");
        break;
  */
  switch (image->type) {
    case SCHAR:
    case UCHAR:
      arrayu8=(unsigned char ***)image->array;
      break;
  case SSHORT:
  case USHORT:
    arrayu16=(unsigned short int ***)image->array;
    break;
  case FLOAT:
    arrayflt=(float ***)image->array;
    break;
  case TYPE_UNKNOWN:
    default:
      VT_Error("image type unknown or not supported for this program",program);
      return( 0 );
  }

  /*        i=findInVec((int) array8[z][y][x], par.droites, par.N);
          if (i<0) continue;
          array8[z][y][x]=(char)par.gauches[i];
          break;
  */
  /*        i=findInVec((int) array16[z][y][x], par.droites, par.N);
          if (i<0) continue;
          array16[z][y][x]=(short int)par.gauches[i];
          cpt++;
        break;
  */


  int Set;
  int ind;
  int new_val;

  for(z=0;z<(int)image->dim.z;z++)
  for(y=0;y<(int)image->dim.y;y++)
  for(x=0;x<(int)image->dim.x;x++) {
    if (par.flagLut) {
        switch (image->type) {
        case SCHAR:
        case UCHAR:
            ind=(int) arrayu8[z][y][x];
            break;
        case SSHORT:
        case USHORT:
            ind=(int) arrayu16[z][y][x];
            break;
        case FLOAT:
            /*ind=(int) arrayflt[z][y][x];
            break;*/
        case TYPE_UNKNOWN:
        default:
            fprintf(stdout, "toto\n");
            VT_Error("image type unknown or not supported for this program",program);
            VT_FreeImage( image );
            VT_Free( (void**)&image );
           return( 0 );
        }
        if ( ind>=lut_size)
        {
            if (par.flagExclusive)
              new_val=0;
            else
              continue;
        } else{
            if (ind<0) {
                fprintf(stdout, "Negative original value error: old_val[%d][%d][%d]=%d\n",z,y,x,ind);
                VT_Error("old label found in old image does not satisfy positivity constraint, exiting...\n",program);
                VT_FreeImage( image );
                VT_Free( (void**)&image );
                if (lut != NULL)
                {
                    free(lut);
                    lut=NULL;
                }
               return( 0 );
            }
            else {
                new_val=lut[ind];
            }
        }
        if (new_val <0)
            continue;

        switch (image->type) {
        case SCHAR:
        case UCHAR:
            arrayu8[z][y][x] = (unsigned char) new_val;
            break;
        case SSHORT:
        case USHORT:
            arrayu16[z][y][x] = (unsigned short int) new_val;
            break;
        case FLOAT:
            /*arrayflt[z][y][x] = (float) new_val;
            break;*/
        case TYPE_UNKNOWN:
        default:
            fprintf(stdout, "toto2\n");
            VT_Error("image type unknown or not supported for this program",program);
            VT_FreeImage( image );
            VT_Free( (void**)&image );
            if (lut != NULL)
            {
                free(lut);
                lut=NULL;
            }
           return( 0 );
        }

    }
    else {
        Set=0;
        switch (image->type) {
        case SCHAR:
        case UCHAR:
          val = (float) arrayu8[z][y][x];
          for (ind=0; ind<par.V; ind++) {
            if (val>=par.vector[ind][0] && val<=par.vector[ind][1] )
            {
              arrayu8[z][y][x]=(unsigned char)0;
              Set=1;
              break;
            }
          }
          if(Set==1) continue;
          for (ind=0; ind<par.S; ind++) {
              i=findInVec(val, par.set, par.S);
              if (i<0) continue;
              arrayu8[z][y][x]=(unsigned char)0;
              Set=1;
          }
          if(Set==1) continue;
          i=findInVec(val, par.droites, par.N);
          if (i<0) continue;
          arrayu8[z][y][x]=(unsigned char)par.gauches[i];
          break;
        case SSHORT:
        case USHORT:
            d=(unsigned short int) arrayu16[z][y][x];
            val = (float) d;
            for (ind=0; ind<par.V; ind++) {
              if (val>=par.vector[ind][0] && val<=par.vector[ind][1] )
              {
                arrayu16[z][y][x]=(unsigned short int)0;
                Set=1;
                break;
              }
            }
            if(Set==1) continue;
            for (ind=0; ind<par.S; ind++) {
                i=findInVec(val, par.set, par.S);
                if (i<0) continue;
                arrayu16[z][y][x]=(unsigned short int)0;
                Set=1;
            }
            if(Set==1) continue;
            i=findInVec(val, par.droites, par.N);
            if (i<0) {continue;}
            arrayu16[z][y][x]=(unsigned short int)par.gauches[i];
          break;
        case FLOAT:
            val = (float) arrayflt[z][y][x];
            for (ind=0; ind<par.V; ind++) {
              if (val>=par.vector[ind][0] && val<=par.vector[ind][1] )
              {
                arrayflt[z][y][x]=0.0;
                Set=1;
                break;
              }
            }
            if(Set==1) continue;
            for (ind=0; ind<par.S; ind++) {
                i=findInVec(val, par.set, par.S);
                if (i<0) continue;
                arrayflt[z][y][x]=0.0;
                Set=1;
            }
            if(Set==1) continue;
            i=findInVec(val, par.droites, par.N);
            if (i<0) {continue;}
            arrayflt[z][y][x]=par.gauches[i];
          break;
        case TYPE_UNKNOWN:
        default:
          VT_Error("image type unknown or not supported for this program",program);
          VT_FreeImage( image );
          VT_Free( (void**)&image );

          return( 0 );
        }
    }
  }

  fprintf(stdout, "Write out at %s\n", par.names.out);

  VT_WriteInrimageWithName( image, par.names.out );

  /*--- liberations memoires ---*/
  if (lut != NULL)
  {
      free(lut);
      lut=NULL;
  }
  VT_FreeImage( image );
  VT_Free( (void**)&image );

  return( 0 );
}








static void MT_Parse( int argc,
                      char *argv[],
                      local_par *par )
{
  int i, nb;
  char text[STRINGLENGTH];
  int status;

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

      else if ( strcmp ( argv[i], "-p" ) == 0 || strcmp ( argv[i], "-pair" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -p...\n", 0 );
        status = sscanf( argv[i],"%f",&(par->gauches[par->N]) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -p...\n", 0 );
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -p...\n", 0 );
        status = sscanf( argv[i],"%f",&(par->droites[par->N]) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -p...\n", 0 );
        par->N += 1;
        i += 1;
        while ( i < argc)  {
          status = sscanf( argv[i],"%f",&(par->gauches[par->N]) );
          if ( status <= 0 ) {i--;  break;}
          i += 1;
          if ( i >= argc)    MT_ErrorParse( "parsing -p...\n", 0 );
          status = sscanf( argv[i],"%f",&(par->droites[par->N]) );
          if ( status <= 0 ) MT_ErrorParse( "parsing -p...\n", 0 );
          par->N += 1;
          i += 1;
        }
      }

      else if ( strcmp ( argv[i], "-sv" ) == 0 || strcmp ( argv[i], "-v" ) == 0 || strcmp ( argv[i], "-setvector" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -sv...\n", 0 );
        status = sscanf( argv[i],"%f",&(par->vector[par->V][0]) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -sv...\n", 0 );
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -sv...\n", 0 );
        status = sscanf( argv[i],"%f",&(par->vector[par->V][1]) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -sv...\n", 0 );
        par->V += 1;
        i += 1;
        while ( i < argc)  {
          status = sscanf( argv[i],"%f",&(par->vector[par->V][0]) );
          if ( status <= 0 ) {i--;  break;}
          i += 1;
          if ( i >= argc)    MT_ErrorParse( "parsing -sv...\n", 0 );
          status = sscanf( argv[i],"%f",&(par->vector[par->V][1]) );
          if ( status <= 0 ) MT_ErrorParse( "parsing -sv...\n", 0 );
          par->V += 1;
          i += 1;
        }
      }

      else if ( strcmp ( argv[i], "-s" ) == 0 || strcmp ( argv[i], "-set" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -s...\n", 0 );
        status = sscanf( argv[i],"%f",&(par->set[par->S]) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -s...\n", 0 );
        par->S += 1;
        i += 1;
        while ( i < argc)  {
          status = sscanf( argv[i],"%f",&(par->set[par->S]) );
          if ( status <= 0 ) {i--;  break;}
          par->S += 1;
          i += 1;
        }
      }

      else if ( strcmp ( argv[i], "-lut" ) == 0 ) {
        par->flagLut = 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -lut...\n", 0 );
        i += 1;
        strncpy( par->names.ext, argv[i], STRINGLENGTH );
      }

      else if ( strcmp ( argv[i], "-e" ) == 0 || strcmp(argv[i], "-exclusive") == 0) {
        par->flagExclusive = 1;
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
      else {
        fprintf(stderr, "N=%d\nS=%d\nV=%d\n", par->N, par->S, par->V);
        MT_ErrorParse("too much file names when parsing\n", 0 );
      }
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
  par->N=0;
  par->S=0;
  par->V=0;
  par->flagLut=0;
  par->flagExclusive=0;
}


static int findInVec(float val, float *vec, int n)
{
  if (val<0) return(-2);
  int i;
  for (i=0; i<n ; i++)
    if (vec[i]==val) return(i);
  return (-1);
}



static int VT_AllocLutWithFileName( int ** lut_pointer, int *lut_size, char* file_name )
{
    char *proc="VT_AllocLutWithFileName";

    int *lut = NULL;

    FILE *f;
    char line[512];
    int l1, l2;
    int o;
    int n = -1;


    if ((f = fopen (file_name, "r")) == NULL) {
      fprintf( stderr, "%s: unable to open '%s' for reading\n", proc, file_name );
      return( -1 );
    }


    while ( fgets(line, 512, f) != NULL ) {
        o = sscanf( line, "%d %d\n", &l1, &l2);
        /*fprintf(stdout, "line %s : %d %d\n", line, l1, l2);*/
        if ( o == 2 ) {
            if (l1<0) {
              if (_verbose_)
                fprintf( stderr, "%s: strictly negative old label value (%d) detected in '%s' ; Exiting...\n", proc, l1, file_name );
              fclose(f);
              return( -1 );
            }
            n = (l1>n) ? l1 : n;
        }
    }

    fclose( f );
    if (n<0)
    {
        if (_verbose_)
            fprintf( stderr, "%s: no positive old values detected in '%s'\n", proc, file_name );
        return( 1 );
    }

    n++;


    lut = (int*)malloc(n * sizeof(int));
    if (lut == NULL) {
      fprintf( stderr, "%s: unable to allocate the lut ; Exiting...\n", proc);
      return( -1 );
    }

    *lut_pointer=lut;

    *lut_size = n;
    return(1);
}


static int VT_FillLutWithFileName( int **lut_pointer, int lut_size, char* file_name, int flagExclusive )
{
    char *proc="VT_AllocLutWithFileName";

    /*int *h = NULL;*/

    FILE *f;
    char line[512];
    int l1, l2;
    int o;
    int i;
    int *lut=(int*) *lut_pointer;
    /*int max=0;*/

    /*fprintf(stdout, "init lut\n");*/
    int defaultvalue=(flagExclusive) ? 0 : -1;

    for (i = 0; i<lut_size ; i++)
        lut[i]=defaultvalue;

    /*fprintf(stdout, "read lut file\n");*/
    if ((f = fopen (file_name, "r")) == NULL) {
      fprintf( stderr, "%s: unable to open '%s' for reading\n", proc, file_name );
      return( -1 );
    }


    /*fprintf(stdout, "set lut\n");*/
    while ( fgets(line, 512, f) != NULL ) {
        o = sscanf( line, "%d %d\n", &l1, &l2);
        if ( o == 2 ) {
            if (l1<0 || l1 >=lut_size){
              if (_verbose_)
                fprintf( stderr, "%s: strictly negative or too high old label value (%d) detected in '%s' ; Exiting...\n", proc, l1, file_name );
              fclose(f);
              return( -1 );
            }
            if (l2 < 0)
            {
                if (_verbose_)
                  fprintf( stderr, "%s: strictly negative new label value (%d) detected in '%s' ; Exiting...\n", proc, l2, file_name );
                fclose(f);
                return( -1 );
            }
            lut[l1] = l2;
            /*max = (l2>max) ? l2 : max;*/
        }
    }

    fclose( f );
    return(1);
}
