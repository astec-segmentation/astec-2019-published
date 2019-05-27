/*************************************************************************
 * directionHistogramMaxima.c -
 *
 * $Id: directionHistogramMaxima.c,v 1.0 2017/08/08 15:10:00 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2017/08/08
 *
 * ADDITIONS, CHANGES
 *
 */

#include <time.h>

#include <vtmalloc.h>

#include <vt_common.h>
#include <vt_symmetryPlane.h>


#define SEUIL 1e-4
#define NMAX 256

typedef struct local_par {

  vt_names names;
  int nmax;
  int max[NMAX];
  double elevation_fraction;
} local_par;



/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );


static char *usage = "[image-in] [file-out] [-max %d [%d [...]]] [-frac %lf]\n\
\t [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -max %d [%d [...]] : direction extraite la %d-eme plus grande (0=plus grand)\n\
\t -frac %lf : directions extraites : celles dont le maximum est superieur ou egal a l'elevation max \n\
\t          multipliee par la fraction donnee (de valeur >=0 (0.5 = defaut) et inferieure a 1)\n\
\t Note: les deux options ne sont pas prevues pour etre utilisees en simultane\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
    /*   char *temp; */
    local_par par;
    vt_image *imageSphere=NULL;
    /*int dim[3];*/
    FILE *fileout;

    /*  liste des maxima de l'histogramme directionnel */
    double *nx, *ny, *nz;
    double *values;
    int taille, i;

    /*--- initialisation des parametres ---*/
    MT_InitParam( &par );

    /*--- lecture des parametres ---*/
    MT_Parse( argc, argv, &par );

    imageSphere = _VT_Inrimage( par.names.in );
    if ( imageSphere == (vt_image*)NULL )
    {
      fprintf(stderr, "%s unreadable\n", par.names.in);
      MT_ErrorParse("unable to read sphere image \n", 0);
    }
    /*dim[0]=imageSphere->dim.x;
    dim[1]=imageSphere->dim.y;
    dim[2]=imageSphere->dim.z;
    */
    /*if (sphereToNormal(imageSphere->array, imageSphere->type, dim, n) != 1)
    {
      VT_FreeImage( imageSphere );
      fprintf(stderr, "%s unreadable\n", par.names.in);
      MT_ErrorParse("unable to compute normal from sphere image \n", 0);
    }*/

    if (VT_ExtractMaxima(imageSphere, &nx, &ny, &nz, &values, &taille) != 1)
    {

      VT_FreeImage( imageSphere );
      MT_ErrorParse("maxima extraction failed \n", 0);
    }
    if(_VT_VERBOSE_) {
      fprintf(stdout, "Diagramme des directions : %d maxima locaux\n",taille);
      for (i=0; i<taille;  i++)
            fprintf(stdout, "normale %d\t = { %f\t %f\t %f }\t\tvalue = %f\n",i,nx[i],ny[i],nz[i],values[i]);
    }

    VT_FreeImage( imageSphere );

    if (par.nmax > 0)
    {
      for (i=0; i<par.nmax; i++)
        if(par.max[i] >= taille)
        {
          vtfree(nx);
          vtfree(ny);
          vtfree(nz);
          vtfree(values);
          MT_ErrorParse("not enough maxima found regarding to the parsed parameters\n", 0);
        }
    }

    if(_VT_VERBOSE_) {
      fprintf(stdout, "Diagramme des directions : %d maxima locaux\n",taille);
      for (i=0; i<taille;  i++)
            fprintf(stdout, "normale %d\t = { %f\t %f\t %f }\t\tvalue = %f\n",i,nx[i],ny[i],nz[i],values[i]);
    }


    if (strcmp(par.names.out, ">") == 0)
        fileout=stdout;
    else
        fileout = fopen( par.names.out, "w");

    if (par.nmax > 0)
    {
        for (i=0;i<par.nmax;i++)
            fprintf(fileout, "%f %f %f\n",nx[par.max[i]],ny[par.max[i]],nz[par.max[i]]);
    }
    else {
        for (i=0; i<taille;  i++) {
            if (values[i]<values[0]*par.elevation_fraction)
                break;
            fprintf(fileout, "%f %f %f\n",nx[i],ny[i],nz[i]);
        }
    }

    if (strcmp(par.names.out, ">") != 0)
        fclose( fileout );

    /*--- liberations memoires ---*/
    vtfree(nx);
    vtfree(ny);
    vtfree(nz);
    vtfree(values);

    return(0);

}















static void MT_Parse( int argc,
                      char *argv[],
                      local_par *par )
{
  int i, nb, status;
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

      /* Parametres de calcul */

      else if ( strcmp ( argv[i], "-frac" ) == 0 ) {
                i += 1;
                if ( i >= argc)    MT_ErrorParse( "parsing -frac...\n", 0 );
                status = sscanf( argv[i],"%lf",&(par->elevation_fraction) );
                if ( status <= 0 ) MT_ErrorParse( "parsing -frac...\n", 0 );
      }


      else if ( strcmp ( argv[i], "-max" ) == 0 ) {
                i += 1;
                if ( i >= argc)    MT_ErrorParse( "parsing -max...\n", 0 );
                status = sscanf( argv[i],"%d",&(par->max[par->nmax]) );
                if ( status <= 0 ) MT_ErrorParse( "parsing -max...\n", 0 );
                par->nmax = 1;
                while(1) {
                    i += 1;
                    if (i >= argc) {
                        i -= 1;
                        break;
                    }
                    int tmp;
                    status = sscanf( argv[i],"%d",&tmp );
                    if (status <= 0) {
                        i -= 1;
                        break;
                    }
                    if (par->nmax==NMAX) MT_ErrorParse("parsing -max : maximum number of maximums reached...\n",0);
                    par->max[par->nmax++]=tmp;
                }
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
      strcpy( par->names.out,  ">" );  /* standart output */
  }
  if (nb == 1) {
    strcpy( par->names.out,  ">" );  /* standart output */
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
  par->elevation_fraction=0.5;
  par->nmax=0;
}
