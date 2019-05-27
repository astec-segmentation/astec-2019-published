/*************************************************************************
 * testTVsample.c -
 *
 * $Id: testTVsample.c,v 1.0 2013/12/18 17:54:00 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2013/12/18
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>
#include <transfo.h>
#include <convert.h>

#include <mt_membrane3D.h>

/* #define NMAXSCALES 100; */

typedef struct local_par {
  int logfile;

  int niter;
  int nangles;
  int nsticks;

  double scale;
  /*   double scales[100]; */
  double zfact;
  vt_names names;
  int initHessian;

  double sample;
  int power;
  MT_SAMPLINGMODE samplingmode;


  double pt[3];
  double normal[3];

  double value;

  vt_4vpt dim;
  vt_fpt voxel;
  ImageType type;

} local_par;





/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
static int _verbose_ = 1;





static int drawPlane( void *inputBuf, void *resultBuf, int *theDim, bufferType type, int *pt, double *normal, double value );
static int SampleBinRegular(vt_image *theIm, double sample, int *thePt, double *ray);


static char *usage = "[outprefix] [-niter] [-nangles %d] [-nsticks %d] [-scale %lf] [-zfact %lf] [-hessian]\n\
\t [-sample %lf [-power]] [-random] [-regular]\n\
\t [-dim %d %d [%d] | [-x %d] [-y %d] [-z %d]] [-v %d] [-template %s]\n\
\t [-voxel | -pixel | -vs %f %f [%f] | [-vx %f] [-vy %f] [-vz %f] ]\n\
\t [-o|-b|-bytes %d [-r|-f] [-s]]\n\
\t [-v] [-D] [-log] [-help]";

static char *detail = "\
\t si 'outprefix' est absent, on prendra stdout\n\
\t -niter %d : nombre d'iterations\n\t si cette option n'est pas presente, on fait une seule iteration\n\
\t -nangles %d : nombre d'angles pour discretiser la demi-sphere\n\
\t -nsticks %d : nombre d'angles pour discretiser le plan normal au champ de tenseur plate\n\
\t -scale %lf : echelle du vote\n\
\t -zfact %lf : facteur de resolution selon l'axe des z\n\
\t -hessian : initialisation par des tensors sticks en fct de la normale au plan (pas d'iteration)\n\
\t -sample %lf : echantillonne les votants de l'image initiale selon le \
coefficient \n\
\t -power : l'argument de -sample est alors converti en puissance pour l'echantillonnage : %age=10^-argument\n\
\t -regular : echantillonne de maniere reguliere \n\
\t -random : echantillonne de maniere aleatoire \n\
\t [-dim %d %d [%d]]      # output image dimensions\n\
\t [-x %d]                # X dimension of the ouput image\n\
\t [-y %d]                # Y dimension of the ouput image\n\
\t [-z %d]                # Z dimension of the ouput image\n\
\t [-template %s]         # template image for the dimensions\n\
\t                          of the output image\n\
\t [-voxel %f %f [%f]]    # output image voxel sizes\n\
\t -v : mode verbose\n\
\t -D : mode debug\n\
\t -log : ecrit un logfile .log\n\
\t -help : help\n\
\t [-o %d]            # number of bytes for encoding\n\
\t [-r]                   # real\n\
\t [-f]                   # fixed (ie integer)\n\
\t [-s]                   # signed\n\
\t  e.g. '-o 1'    = unsigned char\n\
\t       '-o 2 -s' = signed short int\n\
\t       '-o 4 -r' = float\n\
 $Revision: 1.1 $ $Date: 2013/12/18 08:55:56 $ $Author: gael $\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  char Flog[DOUBLESTRINGLENGTH];
  FILE *flog=NULL;

  int writeImages=0;

  vt_image *image, imres, *theIm, imtht, *theTht, imphi, *thePhi;
  vt_image *imtemplate;
  int thePt[3];
  int i;
  double theNormal[3];
  char name[DOUBLESTRINGLENGTH];
  vt_image* imagesIn[3];
  double ray=0;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );

  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );

  /*   fprintf(stdout, "imIn=%s...\n", par.names.in ); */
  /*   fprintf(stdout, "imOut=%s...\n", par.names.out ); */

  /* LOGFILE */
  if ( par.logfile != 0 )
  {
    sprintf(Flog, "%s.log", par.names.in);
    flog = fopen (Flog, "a" );
  }
  if (flog == NULL)
  {
    perror (Flog);
  }


  /* creation de l'image
   */

  image = _VT_Inrimage( par.names.in );

  if ( image == (vt_image*)NULL ) {
    /* initializing result image
       - with reference image, if any
       - with parameters, if any
    */

    VT_Image( &imres );

    if ( par.names.ext[0] != '\0' ) {
      imtemplate = _VT_Inrimage( par.names.ext );
      if ( imtemplate == (vt_image*)NULL ) {
        VT_ErrorParse("unable to read template image \n", 0);
      }
      sprintf( name, "%s.hdr", par.names.in );

      VT_InitVImage( &imres, name, imtemplate->dim.v,
                     imtemplate->dim.x, imtemplate->dim.y, imtemplate->dim.z,
                     imtemplate->type );
      imres.siz.x = imtemplate->siz.x;
      imres.siz.y = imtemplate->siz.y;
      imres.siz.z = imtemplate->siz.z;

      VT_FreeImage( imtemplate );
      VT_Free( (void**)&imtemplate );

      if ( par.dim.v > 0 ) imres.dim.v = par.dim.v;
      if ( par.dim.x > 0 ) imres.dim.x = par.dim.x;
      if ( par.dim.y > 0 ) imres.dim.y = par.dim.y;
      if ( par.dim.z > 0 ) imres.dim.z = par.dim.z;

      if ( par.voxel.x > 0.0 ) imres.siz.x = par.voxel.x;
      if ( par.voxel.y > 0.0 ) imres.siz.y = par.voxel.y;
      if ( par.voxel.z > 0.0 ) imres.siz.z = par.voxel.z;

      if ( par.type != TYPE_UNKNOWN ) imres.type = par.type;
    }

    else {

      if ( par.dim.x < 0 && par.dim.y < 0 ) {
        VT_ErrorParse( "negative X and Y dimensions\n", 0);
      }

      if ( par.dim.v < 0 ) par.dim.v = 1;
      if ( par.dim.x < 0 ) par.dim.x = 1;
      if ( par.dim.y < 0 ) par.dim.y = 1;
      if ( par.dim.z < 0 ) par.dim.z = 1;

      if ( par.voxel.x < 0.0 ) par.voxel.x = 1.0;
      if ( par.voxel.y < 0.0 ) par.voxel.y = 1.0;
      if ( par.voxel.z < 0.0 ) par.voxel.z = 1.0;

      if ( par.type == TYPE_UNKNOWN ) par.type = UCHAR;
      sprintf( name, "%s.hdr", par.names.in );
      VT_InitVImage( &imres, name, par.dim.v,
                     par.dim.x, par.dim.y, par.dim.z,
                     par.type );

      imres.siz.x = par.voxel.x;
      imres.siz.y = par.voxel.y;
      imres.siz.z = par.voxel.z;
    }


    if ( VT_AllocImage( &imres ) != 1 ) {
      VT_ErrorParse("unable to allocate output image\n", 0);
    }

    /* VT_AllocImage met l'image a zero
       VT_FillImage( &imres, 0.0 );
    */

    theIm = &imres;

  }
  else {
    theIm = image;
    if ( par.names.out[0] != '\0' ) {
      sprintf( name, "%s.hdr", par.names.in );
      if ( 0 )
        fprintf( stderr, "%s: set '%s' as output image name\n", program, name );
      if ( VT_CopyName( image->name, name ) == 0 ) {
        VT_ErrorParse("unable to copy image name\n", 0);
      }
    }
  }

  /* ... */
  int theDim[3];
  theDim[0] = theIm->dim.x;
  theDim[1] = theIm->dim.y;
  theDim[2] = theIm->dim.z;



  thePt[0] = (par.pt[0] > 0) ? (int)(par.pt[0]+0.5) : (int)(par.pt[0]-0.5);
  thePt[1] = (par.pt[1] > 0) ? (int)(par.pt[1]+0.5) : (int)(par.pt[1]-0.5);
  thePt[2] = (par.pt[2] > 0) ? (int)(par.pt[2]+0.5) : (int)(par.pt[2]-0.5);

  theNormal[0] = par.normal[0];
  theNormal[1] = par.normal[1];
  theNormal[2] = par.normal[2];

  if ( drawPlane( theIm->buf, theIm->buf, theDim, theIm->type,
                 thePt, theNormal, par.value ) != 1 ) {
    VT_FreeImage( theIm );
    if ( image != (vt_image*)NULL )
      VT_Free( (void**)&image );
    VT_ErrorParse( "error when drawing rectangle\n", 0 );
  }

  if(par.initHessian == 1)
  {
    double theta, phi;

    UnitVectorToSphericalAngles( theNormal, &theta, &phi );
    sprintf( name, "%s.theta.inr", par.names.in );

    VT_InitFromImage( &imtht, theIm, name, (int)FLOAT);

    sprintf( name, "%s.phi.inr", par.names.in );

    VT_InitVImage( &imphi, name, theIm->dim.v,
                   theIm->dim.x, theIm->dim.y, theIm->dim.z,
                   FLOAT);


    if ( VT_AllocImage( &imtht ) != 1 ) {
      VT_FreeImage( theIm );
      if ( image != (vt_image*)NULL )
        VT_Free( (void**)&image );
      VT_ErrorParse("unable to allocate output image\n", 0);
    }


    theTht = &imtht;

    if ( VT_AllocImage( &imphi ) != 1 ) {
      VT_FreeImage( theIm );
      VT_FreeImage( theTht );
      if ( image != (vt_image*)NULL )
        VT_Free( (void**)&image );
      VT_ErrorParse("unable to allocate output image\n", 0);
    }

    thePhi = &imphi;

    theTht->siz.x = par.voxel.x;
    theTht->siz.y = par.voxel.y;
    theTht->siz.z = par.voxel.z;
    thePhi->siz.x = par.voxel.x;
    thePhi->siz.y = par.voxel.y;
    thePhi->siz.z = par.voxel.z;

    float *bufT, *bufP;
    bufT=(float*)theTht->buf;
    bufP=(float*)thePhi->buf;
    for (i=0; i<theDim[0]*theDim[1]*theDim[2]; i++)
    {
      bufT[i]=(float)theta;
      bufP[i]=(float)phi;
    }

    imagesIn[1]=theTht;
    imagesIn[2]=thePhi;
  }

  /*--- ecriture de l'image resultat ---*/
  if ( VT_WriteInrimage( theIm ) == -1 ) {
    VT_FreeImage( theIm );
    if ( image != (vt_image*)NULL )
      VT_Free( (void**)&image );
    if(par.initHessian==1)
    {
      VT_FreeImage( theTht );
      VT_FreeImage( thePhi );
    }
    VT_ErrorParse("unable to write output image\n", 0);
  }

  /* Echantillonnage de l'image */
  if (par.power == 1)
  {
    par.sample = exp((par.sample>=0 ? - par.sample : par.sample)*log(10));
  }

  if (par.logfile != 0 && flog != NULL)
  {
    fprintf(flog, "sampling : %f\n", par.sample);
  }


  if ( par.sample != 1.0 )
  {
    if ( _VT_VERBOSE_ )
      fprintf(stdout, "Sampling step : sampling coefficient = %f\n", par.sample);
    if (par.sample > 0 && par.sample < 1)
    {

      /*--- Image sampling ---*/

      switch (par.samplingmode) {
      case _RANDOM_ :
        if (par.logfile != 0 && flog != NULL)
        {
          fprintf(flog, "mode     : %s\n", "_RANDOM_");
        }
        if (MT_SampleBin(theIm, par.sample) != 1 )
        {
          if (par.initHessian == 0)
          {
            VT_FreeImage( theTht );
            VT_FreeImage( thePhi );
          }
          VT_FreeImage( theIm );
          if ( image != (vt_image*)NULL )
            VT_Free( (void**)&image );
          VT_ErrorParse("unexpected error while sampling binary image\n", 0);
        }
        break;
      case _REGULAR_ :
        if (par.logfile != 0 && flog != NULL)
        {
          fprintf(flog, "mode     : %s\n", "_REGULAR_");
        }
        if(theNormal[0]==theNormal[1] && theNormal[0]==0 && theNormal[2]==1)
        {
          if (par.logfile != 0 && flog != NULL)
          {
            fprintf(flog, "center    : %d  %d  %d\n", thePt[0],thePt[1],thePt[2]);
          }
          if (SampleBinRegular(theIm, par.sample, thePt, &ray) != 1 )
          {
            if (par.initHessian == 0)
            {
              VT_FreeImage( theTht );
              VT_FreeImage( thePhi );
            }
            VT_FreeImage( theIm );
            if ( image != (vt_image*)NULL )
              VT_Free( (void**)&image );
            VT_ErrorParse("unexpected error while sampling binary image\n", 0);
          }
          ray /= sqrt(3);
          if (par.logfile != 0 && flog != NULL)
          {
            fprintf(flog, "ray      : %f\n", ray);
          }
        }
        else
          fprintf(stderr, "_REGULAR_ sampling mode still unimplemented\n");
        break;
      default :
        VT_FreeImage( theIm );
        if ( image != (vt_image*)NULL )
          VT_Free( (void**)&image );
        if(par.initHessian==1)
        {
          VT_FreeImage( theTht );
          VT_FreeImage( thePhi );
        }
        VT_ErrorParse("unable to write output image\n", 0);
      }
    }
    if (1 )
    {
      char tmp[DOUBLESTRINGLENGTH];
      sprintf(tmp, "%s.sample.inr", par.names.in);
      /* fprintf(stdout, "Ecriture de %s ...\n", tmp); */
      VT_WriteInrimageWithName(theIm, tmp);
    }
    else
    {
      if (par.initHessian == 0)
      {
        VT_FreeImage( theTht );
        VT_FreeImage( thePhi );
      }
      VT_FreeImage( theIm );
      if ( image != (vt_image*)NULL )
        VT_Free( (void**)&image );
      VT_ErrorParse("incorrect sampling parameter : 0 < expected value <= 1\n", 0);
    }
  }
  imagesIn[0]=theIm;

  if (par.logfile != 0 && flog != NULL && par.samplingmode == _REGULAR_)
  {
    u8 ***theArray;
    theArray = (u8***)(theIm->array);
    fprintf(flog, "(x0 y0 z0) = %d  %d  %d\n", thePt[0], thePt[1], thePt[2]);
    fprintf(flog, "theArray[z0][y0][x0] = %d\n", (int)theArray[thePt[2]][thePt[1]][thePt[0]]);
    fprintf(flog, "(x0 y0+ray z0) = %d  %d  %d\n", thePt[0], thePt[1]+(int)ray, thePt[2]);
    fprintf(flog, "theArray[z0][y0+ray][x0] = %d\n", (int)theArray[thePt[2]][thePt[1]+(int)ray][thePt[0]]);
  }


  /*--- debut du test ---*/

  vt_3Dtensor theTensor;
  enumTVmode mode = TVCLASSIC;


  /*--- Tensor Voting step ---*/
  if (par.logfile != 0 && flog != NULL)
  {
    fprintf(flog, "TVscale  : %f\n", par.scale);
    fprintf(flog, "hessian  : %s\n", (par.initHessian == 0) ? "0" : "1");
  }

  if (MT_Compute3DTensorVoting(&theTensor, imagesIn, par.scale, par.zfact, par.niter,
    par.nangles,par.nsticks, mode, par.initHessian,
    par.names.in, writeImages) != 1)
  {
    if (par.initHessian == 1)
    {
      VT_FreeImage( theTht );
      VT_FreeImage( thePhi );
    }
    VT_FreeImage( theIm );
    if ( image != (vt_image*)NULL )
      VT_Free( (void**)&image );
    VT_ErrorParse("problem while computing the tensor voting\n", 0 );
    return(-1);
  }


  if (par.logfile != 0 && flog != NULL && par.samplingmode == _REGULAR_)
  {
    vt_image *theTVimg;
    r32 ***theSurf;

    theTVimg = &(theTensor.imvp3);
    theSurf = (r32***)(theTVimg->array);
    fprintf(flog, "theVP3[z0][y0][x0] = %f\n", theSurf[thePt[2]][thePt[1]][thePt[0]]);
    fprintf(flog, "theVP3[z0][y0+ray][x0] = %f\n", theSurf[thePt[2]][thePt[1]+(int)ray][thePt[0]]);
    /*     theTVimg = &(theTensor.imvp2); */
    /*     theSurf = (r32***)(theTVimg->array); */
    /*     fprintf(flog, "theVP2[z0][y0][x0] = %f\n", theSurf[thePt[2]][thePt[1]][thePt[0]]); */
    /*     theTVimg = &(theTensor.imvp1); */
    /*     theSurf = (r32***)(theTVimg->array); */
    /*     fprintf(flog, "theVP1[z0][y0][x0] = %f\n", theSurf[thePt[2]][thePt[1]][thePt[0]]); */
  }

  fprintf(stdout, "Retour dans la fonction principale\n");

  fprintf(stdout, "Ecriture des images...\n");
  VT_Write3Dtensor( &theTensor );


  /*--- liberations memoires ---*/
  VT_Free3Dtensor(&theTensor);
  VT_FreeImage( theIm );
  if ( image != (vt_image*)NULL )
    VT_Free( (void**)&image );
  if (par.initHessian == 1)
  {
    VT_FreeImage( theTht );
    VT_FreeImage( thePhi );
  }

  fprintf(stdout, "Fin des operations\n");

  return( 0 );
}








static void VT_Parse( int argc,
                      char *argv[],
                      local_par *par )
{
  int i, nb, status;
  int  o=0, s=0, r=0;
  char text[STRINGLENGTH];

  if ( VT_CopyName( program, argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);
  if ( argc == 1 ) VT_ErrorParse("\n", 0 );

  /*--- lecture des parametres ---*/
  i = 1; nb = 0;
  while ( i < argc ) {
    if ( argv[i][0] == '-' ) {
      if ( argv[i][1] == '\0' ) {
        VT_ErrorParse("parsing - \n", 1);
      }
      /*--- arguments generaux ---*/
      else if ( strcmp ( argv[i], "-help" ) == 0 ) {
        VT_ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 ) {
        _VT_VERBOSE_ = 1;
        VT_IncrementVerboseInVtTube3D(  );
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 ) {
        _VT_DEBUG_ = 1;
        VT_IncrementDebugInVtTube3D(  );
      }

      else if ( strcmp ( argv[i], "-log" ) == 0 ) {
        par->logfile= 1;
      }

      /*--- arguments specifiques ---*/



      else if ( strcmp ( argv[i], "-niter" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -niter...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->niter) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -niter...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-nangles" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -nangles...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->nangles) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -nangles...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-nsticks" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -nsticks...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->nsticks) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -nsticks...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-scale" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -scale...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->scale) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -scale...\n", 0 );
      }
/*
      else if ( (strcmp ( argv[i], "-scales" ) == 0) ) {
        i += 1;
        int j=0;
        if ( i >= argc)    VT_ErrorParse( "parsing -scales ...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->scales[0]) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -scales ...\n", 0 );
        do {
          j ++;
          if (j>=100)
            VT_ErrorParse( "parsing -scales ... too much scales\n", 0 );
          i++;

          if ( i < argc) {
            status = sscanf( argv[i], "%lf", &(par->membraneSurFond[j]) );
            if ( status <= 0 ) {
              i--;
              par->membraneSurFond[j] = -1;
            }
          }
          else {
            status=0;
          }
        } while(status>0);
        par->nMSF = j;
      }
*/



      else if ( strcmp ( argv[i], "-zfact" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -zfact...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->zfact) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -zfact...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-hessian" ) == 0 ) {
        par->initHessian = 1;
      }

      else if ( strcmp ( argv[i], "-sample" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -sample...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->sample) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -sample...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-power" ) == 0 ) {
        par->power = 1;
      }

      else if ( strcmp ( argv[i], "-regular" ) == 0 ) {
        par->samplingmode = _REGULAR_;
      }

      else if ( strcmp ( argv[i], "-random" ) == 0 ) {
        par->samplingmode = _RANDOM_;
      }


      /* image value
       */
      else if ( strcmp ( argv[i], "-value" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -value...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->value) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -value...\n", 0 );
      }

      /*--- dimension de l'image ---*/
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
        if ( i < argc ) {
          status = sscanf( argv[i], "%d", &(par->dim.z) );
          if ( status <= 0 ) {
            i--;
          }
        }
      }

      else if ( strcmp ( argv[i], "-x") == 0 && argv[i][2] == '\0' ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -x...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->dim.x) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -x...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-y" ) == 0  && argv[i][2] == '\0' ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -y...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->dim.y) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -y...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-z" ) == 0  && argv[i][2] == '\0' ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -z...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->dim.z) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -z...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-v" ) == 0   && argv[i][2] == '\0' ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -v...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->dim.v) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -v...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-voxel" ) == 0
                || strcmp ( argv[i], "-pixel" ) == 0
                || (strcmp ( argv[i], "-vs" ) == 0  && argv[i][3] == '\0') ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -voxel..\n", 0 );
        status = sscanf( argv[i],"%f",&(par->voxel.x) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -voxel...\n", 0 );
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -voxel..\n", 0 );
        status = sscanf( argv[i],"%f",&(par->voxel.y) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -voxel...\n", 0 );
        i += 1;
        if ( i < argc ) {
          status = sscanf( argv[i],"%f",&(par->voxel.z) );
          if ( status <= 0 ) {
            i--;
          }
        }
      }

      else if ( strcmp ( argv[i], "-vx" ) == 0   && argv[i][3] == '\0' ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -vx...\n", 0 );
        status = sscanf( argv[i], "%f", &(par->voxel.x) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -vx...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-vy" ) == 0   && argv[i][3] == '\0' ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -vy...\n", 0 );
        status = sscanf( argv[i], "%f", &(par->voxel.y) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -vy...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-vz" ) == 0   && argv[i][3] == '\0' ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -vz...\n", 0 );
        status = sscanf( argv[i], "%f", &(par->voxel.z) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -vz...\n", 0 );
      }



      /* template
       */

      else if ( strcmp ( argv[i], "-template") == 0
           || (strcmp ( argv[i], "-t") == 0 && argv[i][2] == '\0')
           || (strcmp ( argv[i], "-dims") == 0 && argv[i][5] == '\0') ) {
        i++;
        /* fprintf(stdout, "COUCOU\n\t\t%s\t\t%s\n\n", argv[i-1], argv[i]); */
        if ( i >= argc) VT_ErrorParse( "parsing -template\n", 0 );
        strncpy( par->names.ext, argv[i], STRINGLENGTH );
      }


      /* image encoding type
       */

      else if ( strcmp ( argv[i], "-r" ) == 0 ) {
        r = 1;
      }
      else if ( strcmp ( argv[i], "-s" ) == 0 ) {
        s = 1;
      }
      else if ( strcmp ( argv[i], "-bytes" ) == 0
                || (strcmp ( argv[i], "-b" ) == 0 && argv[i][2] == '\0')
                || (strcmp ( argv[i], "-o" ) == 0 && argv[i][2] == '\0') ) {
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
    else if ( argv[i][0] != 0 ) {
      if ( nb == 0 ) {
        strncpy( par->names.in, argv[i], STRINGLENGTH );
        nb += 1;
      } else if ( nb == 1 ) {
        strncpy( par->names.out, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else
        VT_ErrorParse("too much file names when parsing\n", 0 );
    }
    i += 1;
  }
  if (nb == 0) {
    strcpy( par->names.out, ">" );  /* standart output */
  }

  /*--- type de l'image resultat ---*/
  if ( r == 1 ) {
    switch( o ) {
    default :
      VT_ErrorParse( "such byte size not handled for floating types\n", 0 );
      break;
    case 0 :
    case 4 :
      par->type = FLOAT; break;
    case 8 :
      par->type = FLOAT; break;
    }
  }
  else {
    switch( o ) {
    default :
      VT_ErrorParse( "such byte size not handled for integer types\n", 0 );
      break;
    case 0 :
    case 1 :
      par->type = ( s == 1 ) ? SCHAR : UCHAR;
      break;
    case 2 :
      par->type = ( s == 1 ) ? SSHORT : USHORT;
      break;
    case 4 :
      par->type = ( s == 1 ) ? SINT : UINT;
      break;
    case 8 :
      if ( s == 1 )
        VT_ErrorParse( "signed long int not handled yet\n", 0 );
      else
        par->type = ULINT;
      break;
    }
  }

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

  par->logfile = 0;

  par->niter = 2;
  par->nangles = 3;
  par->nsticks = 5;
  par->scale = 2.0;
  /*   par->scales[0] = par->scale; */
  /*  par->scales[1] = 0; */
  par->zfact = 1;
  par->initHessian = 0;
  par->sample=1;
  par->power=0;
  par->samplingmode=_RANDOM_;

  par->value = 255;

  par->dim.x = 256;
  par->dim.y = 256;
  par->dim.z = 64;
  par->dim.v = 1;
  par->voxel.x = par->voxel.y = par->voxel.z = 1.0;
  par->type = UCHAR;

  par->pt[0] = 127.0;
  par->pt[1] = 127.0;
  par->pt[2] = 31.0;

  par->normal[0] = 0;
  par->normal[1] = 0;
  par->normal[2] = 1;

}


static int drawPlane( void *inputBuf,
                   void *resultBuf,
                   int *theDim,
                   bufferType type,
                   int *pt,
                   double *normal,
                   double value )
{
  char *proc = "drawPlane";
  size_t size;
  int x, y, z, i;

  size = (size_t)theDim[0] * (size_t)theDim[1] * (size_t)theDim[2];
  if ( ConvertBuffer( inputBuf, type, resultBuf, type, size ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to copy buffer\n", proc );
    return( 0 );
  }

  switch( type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    return( 0 );

  case UCHAR :
    {
      u8 *theBuf = (u8*)resultBuf;
      u8 v;
      if ( value <= 0.0 ) v = 0;
      else if ( value >= 255.0 ) v = 255;
      else v = (int)(value + 0.5);
      if (normal[0]==normal[1] && normal[0]==0 && normal[2]==1)
      {
        z=pt[2];
        i=z*theDim[0]*theDim[1];
        for(y=0;y<theDim[1];y++)
        for(x=0;x<theDim[0];x++, i++)
          theBuf[i] = v;
      }
      else
      {
        double sum=sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
        normal[0]/=sum;
        normal[1]/=sum;
        normal[2]/=sum;
        for ( i=0, z=0; z<theDim[2]; z++ )
        for ( y=0; y<theDim[1]; y++ )
        for ( x=0; x<theDim[0]; x++, i++ ) {
          if ( (int) (normal[0]*(x-pt[0]) + normal[1]*(y-pt[1]) + normal[2]*(z-pt[2]) + 0.5 ) == 0   ) {
            theBuf[i] = v;
          }
        }
      }
    }
    break;

  }
  return( 1 );
}

static int SampleBinRegular(vt_image *theIm, double sample, int *thePt, double *Ray)
{
  char *proc = "SampleBinRegular";
  int theDim[3];
  double xp, yp;
  double xm, ym;
  int i, j;
  int n1=0;

  bufferType type;

  double x0;
  double y0=(double)thePt[1];
  int z0=thePt[2];

  double ray;
  ray = sqrt(2/(sqrt(3)*sample));
  *Ray = ray;
  /* fprintf(stdout, "thePt = [%d  %d  %d]\n", thePt[0],thePt[1],thePt[2]); */

  theDim[0] = theIm->dim.x;
  theDim[1] = theIm->dim.y;
  theDim[2] = theIm->dim.z;
  type = theIm->type;

  switch (type){
  case UCHAR :
  {
    u8 ***theArray;
    theArray = (u8***)(theIm->array);
    /* for(k=0;k<theDim[2];k++) */
    for(j=0;j<theDim[1];j++)
    for(i=0;i<theDim[0];i++)
    {
      theArray[z0][j][i]=(u8)0;
    }
    yp=ym=y0;
    ym-=sqrt(3)*ray/2;
    x0=(double)thePt[0];
    /* fprintf(stdout, "(x0 y0 z0) = [%d  %d  %d]\n", (int)((double)x0+0.5),(int)((double)y0+0.5), z0); */
    while(yp<(double)theDim[1]-0.5)
    {
      xp=xm=x0;
      while(xp<(double)theDim[0]-0.5)
      {
        if (theArray[z0][(int)(yp+0.5)][(int)(xp+0.5)] == (u8)0)
        {
          n1++;
          theArray[z0][(int)(yp+0.5)][(int)(xp+0.5)] = (u8)255;
          /* if (n1 == 1) */
	  /* fprintf(stdout, "(x y z) = [%d  %d  %d]\n", (int)(xp+0.5), (int)(yp+0.5), z0); */
        }
        xp+=ray;
      }
      while(xm>=-0.5)
      {
        if (theArray[z0][(int)(yp+0.5)][(int)(xm+0.5)] == (u8)0)
        {
          n1++;
          theArray[z0][(int)(yp+0.5)][(int)(xm+0.5)] = (u8)255;
        }
        xm-=ray;
      }
      x0=x0-ray/2;
      if(x0<-0.5)
        x0+=ray;
      yp+=sqrt(3)*ray/2;
    }
    x0=(double)thePt[0]-ray/2;
    if(x0<-0.5)
      x0+=ray;
    while(ym>=-0.5)
    {
      xp=xm=x0;
      while(xp<(double)theDim[0]-0.5)
      {
        if (theArray[z0][(int)(ym+0.5)][(int)(xp+0.5)] == (u8)0)
        {
          n1++;
          theArray[z0][(int)(ym+0.5)][(int)(xp+0.5)] = (u8)255;
        }
        xp+=ray;
      }
      while(xm>=-0.5)
      {
        if (theArray[z0][(int)(ym+0.5)][(int)(xm+0.5)] == (u8)0)
        {
          n1++;
          theArray[z0][(int)(ym+0.5)][(int)(xm+0.5)] = (u8)255;
        }
        xm-=ray;
      }
      x0=x0-ray/2;
      if(x0<-0.5)
        x0+=ray;
      ym-=sqrt(3)*ray/2;
    }
    /* fprintf(stdout, "theArray[31][127][127] = %d\n", (int)theArray[31][127][127]); */
    /* fprintf(stdout, "theArray[31][128][127] = %d\n", (int)theArray[31][128][127]); */
    break;
  }
  default:
    fprintf( stderr, "%s: such image type not handled yet\n", proc );
    return( 0 );
  }
  return( 1 );
}

