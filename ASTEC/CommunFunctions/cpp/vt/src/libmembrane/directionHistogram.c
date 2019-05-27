/*************************************************************************
 * directionHistogram.c -
 *
 * $Id: directionHistogram.c,v 1.0 2014/03/18 14:12:00 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2014/03/18
 *
 * ADDITIONS, CHANGES
 *
 */

#include <time.h>


#include <transfo.h>
#include <vtmalloc.h>

#include <vt_common.h>
#include <vt_morpho.h>
#include <morphotools.h>
#include <chamferdistance.h>

#define PI 3.141592654

typedef struct local_par {

  vt_names names;
  double rayon;
  double sigma;
  /*  int writeImages; */
  int bin;
  int angle;
  double coef;
} local_par;



/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );

/*
static int drawSphereCrown(void *BufOut,
		int *theDim,
		bufferType type,
		double *center,
		double radius,
		double value );
*/
		
static int drawSphereCrown2( void *BufOut,
		int *theDim,
		bufferType type,
		double *center,
		double radius,
		double value );

static double gaussianFun(double x, double sigma2x2);

static int addToSphere(void ***array, 	/* array de la sphere-histogramme */
		       int *dim,		/* size of the array */
		       bufferType type,/* type de l'array */
		       double *c, 		/* position du centre de la sphere */
		       int angle,		/* mode de lissage (ps ou angle) */
		       double sigma, 	/* ecart-type du lissage de l'histo */
		       double coef,	/* ponderation du "vote" */
		       double tht, 	/* angle 1 */
		       double phi		/* angle 2 */
						);

static int Neighborhood2Int ( Neighborhood N )
{
  int connectivity = 26;
  switch ( N ) {
  case N04 :
    connectivity = 4; break;
  case N06 :
    connectivity = 6; break;
  case N08 :
    connectivity = 8; break;
  case N10 :
    connectivity = 10; break;
  case N18 :
    connectivity = 18; break;
  case N26 :
    connectivity = 26; break;
  }
  return( connectivity );
}


static char *usage = "[image-in] [image-out]\n\
\t [-sigma %lf] [-coef %lf] [-rayon %lf] [-bin] [-angle | -noangle] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les trois sont absents, on prendra stdin et stdout\n\
\t -sigma %lf : ecart-type pour l'accumulation des votes (defaut : PI/32)\n\
\t -coef %lf : coefficient multiplicatif s'appliquant au parametre sigma, \n\
\t             utile pour proceder a des ajustements par rapport au parametrage par defaut (defaut : 1.0)\n\
\t -rayon %lf : rayon (precision) de la sphere recevant les votes\n\
\t -bin : aucune ponderation des contributions\n\
\t -angle : calcule les contributions en fct de l'ecart d'angle\n\
\t -noangle : calcule les contributions en fct de l'ecart de prod. scalaire\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  char name[DOUBLESTRINGLENGTH];
  local_par par;
  vt_image *imageIn;
  vt_image *imageTht;
  vt_image *imagePhi=NULL;
  vt_image imres;
  vt_image *theIm;
  int dim[3];
  double pt[3];
  bufferType t;
  /*  vt_image* imagesIn[3]; */
  int marker = 0;
  int flag_3D = 1;
  /*  int type = (int)UCHAR; */
  char prefix[STRINGLENGTH];
  int i, j, k;
  float ***theTht, ***thePhi=NULL;
  float ***theIn32=NULL;
  unsigned char ***theInU8=NULL;
  double coef, tht, phi;
  
  clock_t start, stop;
  double elapsed;

  start = clock();


  /*--- initialisation des parametres ---*/
  MT_InitParam( &par );

  /*--- lecture des parametres ---*/
  MT_Parse( argc, argv, &par );

  /*--- lecture des images d'entree ---*/
  sprintf( prefix, "%s", par.names.in );
  i=0;
  while (prefix[i]!='\0')
  {
    i++;
  }
  if (i>4 && prefix[i-4]=='.' &&
      ((prefix[i-3]=='h' && prefix[i-2]=='d' && prefix[i-1]=='r') ||
       (prefix[i-3]=='i' && prefix[i-2]=='n' && prefix[i-1]=='r')))
  {
    prefix[i-4]='\0';
    sprintf( name, "%s", par.names.in );
    marker = 1;
  }
  else
  {
    sprintf( name, "%s.inr", prefix );
  }
  imageIn = _VT_Inrimage( name );
  if ( imageIn == (vt_image*)NULL )
  {
	if ( marker == 1 )
	  MT_ErrorParse("unable to read input image A\n", 0);	  
    sprintf( name, "%s.hdr", prefix );	
	imageIn = _VT_Inrimage( name );
	if ( imageIn == (vt_image*)NULL )
	  MT_ErrorParse("unable to read input image A\n", 0);	  
  }
  if(_VT_VERBOSE_)
	fprintf(stdout, "Image A : %s\n", name);
  switch ( imageIn->type ) {
  case UCHAR:
    theInU8 = (unsigned char ***)imageIn->array;
    break;
  case FLOAT:
    theIn32 = (float ***)imageIn->array;
    break;
  default:
    MT_ErrorParse("image Phi type not handled yet\n", 0);		   
  }
	  
  if ( imageIn->dim.z == 1 )
    flag_3D = 0;

  sprintf( name, "%s.theta.inr", prefix );
  imageTht = _VT_Inrimage( name );
  if ( imageTht == (vt_image*)NULL )
  {
     fprintf(stderr, "%s unreadable\n", name);
     sprintf( name, "%s.theta.hdr", prefix );
     imageTht = _VT_Inrimage( name );
     if ( imageTht == (vt_image*)NULL )
     {
		fprintf(stderr, "%s unreadable\n", name);
	    j=0;
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
		  imageTht = _VT_Inrimage( name );
		  if ( imageTht == (vt_image*)NULL )
		  {
			fprintf(stderr, "%s unreadable\n", name);
			sprintf( name, "%s.theta.hdr", prefix );
			imageTht = _VT_Inrimage( name );
			if ( imageTht == (vt_image*)NULL )
			{
			  VT_FreeImage( imageIn );
			  fprintf(stderr, "%s unreadable\n", name);
			  MT_ErrorParse("unable to read input image B\n", 0);		   
			}
		  }
	    }
	}
  }	
  if(_VT_VERBOSE_)
	fprintf(stdout, "Image d'entree theta : %s\n", name);
  if (imageTht->type != FLOAT) {
	VT_FreeImage( imageIn );
	VT_FreeImage( imageTht );
	MT_ErrorParse("image Theta type not handled yet\n", 0);		   
    
  }
  theTht = (float ***)imageTht->array;

  if (flag_3D == 1)
  {
	sprintf( name, "%s.phi.inr", prefix );
	imagePhi = _VT_Inrimage( name );
	if ( imagePhi == (vt_image*)NULL )
	{
	  fprintf(stderr, "%s unreadable\n", name);
	  sprintf( name, "%s.phi.hdr", prefix );
	  imagePhi = _VT_Inrimage( name );
	  if ( imagePhi == (vt_image*)NULL )
	  {
		VT_FreeImage( imageIn );
		VT_FreeImage( imageTht );
		fprintf(stderr, "%s unreadable\n", name);
		MT_ErrorParse("unable to read input image C\n", 0);		   
	  }
	}
	if(_VT_VERBOSE_)
	  fprintf(stdout, "Image d'entree phi : %s\n", name);
	if (imagePhi->type != FLOAT) {
	  VT_FreeImage( imageIn );
	  VT_FreeImage( imageTht );
	  VT_FreeImage( imagePhi );
	  MT_ErrorParse("image Phi type not handled yet\n", 0);		   
	}
	thePhi = (float ***)imagePhi->array;
  }


  /*--- Spherical histogram initialization ---*/
  sprintf( name, "%s", par.names.out );
  t =  FLOAT;
  dim[0]=(int) par.rayon*2+1;
  dim[1]=dim[0];
  if (flag_3D==1)
	dim[2]=dim[0];
  else
	dim[2]=1;
  pt[0] = par.rayon;
  pt[1] = par.rayon;
  if (flag_3D==1)
	pt[2] = par.rayon;
  else
	pt[2] = 0;
	
  VT_Image( &imres );
  VT_InitVImage( &imres, name, 1, dim[0], dim[1], dim[2], t );
  /*VT_InitVImage( &imres, name, 1,  */
  /*     imtemplate->dim.x, imtemplate->dim.y, imtemplate->dim.z,  */
  /*     imtemplate->type ); */
		     
  imres.siz.x = 1.0;
  imres.siz.y = 1.0;
  imres.siz.z = 1.0;

  if ( VT_AllocImage( &imres ) != 1 ) {
	VT_FreeImage( imageIn );
    VT_FreeImage( imageTht );
	if ( flag_3D != 0 )
		VT_FreeImage( imagePhi );
	MT_ErrorParse("unable to allocate output image\n", 0);
  }


  theIm = &imres;

  if(_VT_VERBOSE_)
    fprintf(stdout, "Draw sphere...\n");
  /*  if ( drawSphereCrown( theIm->buf, dim, t, pt, par.rayon, (double)0 ) != 1 ) { */
  if ( drawSphereCrown2( theIm->buf, dim, t, pt, par.rayon, (double)0 ) != 1 ) {
    VT_FreeImage( &imres );
    VT_FreeImage( imageIn );
    VT_FreeImage( imageTht );
	if ( flag_3D != 0 )
		VT_FreeImage( imagePhi );
    MT_ErrorParse( "error when drawing sphere\n", 0 );
  }

  if(_VT_VERBOSE_)
    fprintf(stdout, "Calculs...\n");

  /*--- calculs ---*/

  if ( flag_3D == 1 ) {
    /* 3D case */
    for (k=0;k<(int)imageIn->dim.z;k++)
    for (j=0;j<(int)imageIn->dim.y;j++)
    for (i=0;i<(int)imageIn->dim.x;i++)
	{
	  switch (imageIn->type) {
	  case UCHAR:
	    coef=(double) (theInU8[k][j][i]); 
	    break;
	  case FLOAT:
	    coef=(double) (theIn32[k][j][i]); 
	    break;
	  default:
		VT_FreeImage( &imres );
		VT_FreeImage( imageIn );
		VT_FreeImage( imageTht );
		VT_FreeImage( imagePhi );
	    MT_ErrorParse( "input image type not handled yet\n", 0 );
	  }
	  tht=theTht[k][j][i];
	  phi=thePhi[k][j][i];
	  if (par.bin == 1) 
		coef = (coef > 0.0) ? 1.0 : 0.0;
    if ( coef > 0.0 && addToSphere(theIm->array, dim, t, pt, par.angle, par.sigma*par.coef, coef, tht, phi) == 0 )
	  {
		VT_FreeImage( &imres );
		VT_FreeImage( imageIn );
		VT_FreeImage( imageTht );
		VT_FreeImage( imagePhi );
	    MT_ErrorParse( "Error while computing histogram\n", 0 );
	  }
	  
	}
    
  }
  else {
    /* 2D case : TODO */
	VT_FreeImage( imageIn );
	VT_FreeImage( imageTht  );
	VT_FreeImage( &imres  );
    MT_ErrorParse( "2D case not handled yet\n", 0 );
  }

  /*--- liberations memoires ---*/
  VT_FreeImage( imageIn );
  VT_FreeImage( imageTht  );
  if (flag_3D ==1)
	VT_FreeImage( imagePhi  );
  
  
  /*--- ecriture de l'image resultat ---*/
  if ( VT_WriteInrimage( &imres ) == -1 ) {
    VT_FreeImage( &imres );
    MT_ErrorParse("unable to write output image\n", 0);
  }
  VT_FreeImage( &imres );


  stop = clock();
  elapsed = (double)(stop-start)/CLOCKS_PER_SEC;

  if(_VT_VERBOSE_)
    fprintf(stdout, "Elapsed time : \t%.1fs\n", elapsed);

  return( 0 );
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
/*
      else if ( strcmp ( argv[i], "-wi" ) == 0 ) {
        par->writeImages = 1;
      }




      else if ( strcmp ( argv[i], "-2D" ) == 0 ) {
        par->dimension = 2;
      }
*/


      /* Parametres de calcul */

      else if ( strcmp ( argv[i], "-radius" ) == 0 ||
		(strcmp ( argv[i], "-rayon" ) == 0 ) ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -rayon...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->rayon) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -rayon...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-sigma" ) == 0 ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -sigma...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->sigma) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -sigma...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-bin" ) == 0 ) {
        par->bin = 1;
      }

      else if ( strcmp ( argv[i], "-angle" ) == 0 ) {
        par->angle = 1;
      }

      else if ( strcmp ( argv[i], "-noangle" ) == 0 ) {
        par->angle = 0;
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

  /*  par->writeImages = 0; */
  par->rayon = 10.0;
  par->angle = 1;
  par->sigma = PI/32;
  par->bin = 0;
  par->coef=1.0;
}

int drawSphereCrown2( void *BufOut,
		int *theDim,
		bufferType type,
		double *center,
		double radius,
		double value )
{
  char *proc = "drawSphereCrown2";
  int x, y, z, i;
  double r2;
  double radius2 = radius*radius;
  Neighborhood local_connexite;
  typeStructuringElement SE;
  float *bufTmp;
  float *theBuf;
  
  switch( type ) {
  default :
    if ( _VT_VERBOSE_ )
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    return( 0 );
  case FLOAT :
    {
      theBuf = (float*)BufOut;
      float v;
      v = (float)value;
      for ( i=0, z=0; z<theDim[2]; z++ )
      for ( y=0; y<theDim[1]; y++ )
      for ( x=0; x<theDim[0]; x++, i++ ) {
	r2 = (x-center[0])*(x-center[0])
	  + (y-center[1])*(y-center[1])
	  + (z-center[2])*(z-center[2]);
	if ( (r2 <= radius2) )
	  theBuf[i] = v;
	else
	  theBuf[i] = -1;
      }
    }

    bufTmp = (float*)vtmalloc( theDim[0]*theDim[1]*theDim[2]*sizeof(float), "bufTmp", proc );
    if (bufTmp==(float*) NULL) 
	  return(0);
    break;
  }

  initStructuringElement( &SE );

  SE.nbIterations = 1;
  local_connexite = N06;

  SE.connectivity = Neighborhood2Int( local_connexite );
  SE.radius = 0;

  

  if ( morphologicalErosion( theBuf, bufTmp, type, theDim, &SE ) != 1 ) {
    freeStructuringElement( &SE );
    return(0);
  }

  for (i=0; i<theDim[0]*theDim[1]*theDim[2]; i++)
	theBuf[i] = (theBuf[i]==bufTmp[i]) ? -1.0 : 0.0;

  freeStructuringElement( &SE );
  vtfree(bufTmp);
  return(1);
}

/*
int drawSphereCrown( void *BufOut,
		int *theDim,
		bufferType type,
		double *center,
		double radius,
		double value )
{
  char *proc = "drawSphereCrown";
  int x, y, z, i;
  double r2;
  double radius2 = radius*radius;
  double radius2m = (radius-1)*(radius-1);


  switch( type ) {
  default :
    if ( _VT_VERBOSE_ )
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    return( 0 );
    
  case UCHAR :
    {
      u8 *theBuf8 = (u8*)BufOut;
      u8 v;
      if ( value <= 0.0 ) v = 0;
      else if ( value >= 255.0 ) v = (unsigned char)255;
      else v = (int)(value + 0.5);

      for ( i=0, z=0; z<theDim[2]; z++ )
      for ( y=0; y<theDim[1]; y++ )
      for ( x=0; x<theDim[0]; x++, i++ ) {
	r2 = (x-center[0])*(x-center[0])
	  + (y-center[1])*(y-center[1])
	  + (z-center[2])*(z-center[2]);
	if ( (r2 <= radius2 ) && (r2 > radius2m))
	  theBuf8[i] = v;
      }
    }
    break;
  case FLOAT :
    {
      float *theBuf = (float*)BufOut;
      float v;
      v = (float)value;

      for ( i=0, z=0; z<theDim[2]; z++ )
      for ( y=0; y<theDim[1]; y++ )
      for ( x=0; x<theDim[0]; x++, i++ ) {
	r2 = (x-center[0])*(x-center[0])
	  + (y-center[1])*(y-center[1])
	  + (z-center[2])*(z-center[2]);
	if ( (r2 <= radius2) && (r2 > radius2m))
	  theBuf[i] = v;
	else
	  theBuf[i] = -1;
      }
    }
    break;
  }
  return( 1 );
}
*/

int addToSphere(void ***array, 			/* array de la sphere-histogramme */
		int *dim,		/* size of the array */
		bufferType type,/* array type */
		double *c, 		/* position du centre de la sphere */
		int angle,		/* mode de lissage (ps ou angle) */
		double sigma, 	/* ecart-type du lissage de l'histo */
		double coef,	/* ponderation du "vote" */
		double tht, 	/* angle 1 */
		double phi		/* angle 2 */
						)
{
  /*fprintf(stdout, "coef = %f\ttht = %f\tphi = %f\n", coef, tht, phi); */
  char *proc = "addToSphere";
  double n[3]; /* vecteur correspondant a la normale a la structure */
  double v[3]; /* vecteur normal a la sphere-histogramme */
  int i, j, k;
  double s;
  double sum=0;
  double ps;
  double sigma2x2=2*sigma*sigma;
  double val;
  double alpha;
  
  unsigned char ***theArrayu8;
  float ***theArray32;
  
  switch (type) {
  case UCHAR:
	theArrayu8=(unsigned char ***)array;
	/*fprintf(stdout, "UCHAR\n"); */
	break;
  case FLOAT:
	theArray32=(float ***)array;
	/*fprintf(stdout, "FLOAT\n"); */
	break;
  default:
    if ( _VT_VERBOSE_ )
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    return( 0 );
  }
  
  for (k=0; k<dim[2]; k++)
  for (j=0; j<dim[1]; j++)
  for (i=0; i<dim[0]; i++) {
	switch (type) {
	case UCHAR:
	  val=(double)theArrayu8[k][j][i];
	  break;
	case FLOAT:
	  val=(double)theArray32[k][j][i];
	  break;
	default:
	  if ( _VT_VERBOSE_ )
        fprintf( stderr, "%s: such image type not handled yet\n", proc );
      return( 0 );
	}
	/*fprintf( stdout, "%1.0f\t",val); */
	if (val < 0)
	  continue;
	/*fprintf( stdout, "\n"); */
	v[0]=(double) i - c[0];
	v[1]=(double) j - c[1];
	v[2]=(double) k - c[2];
	s=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	v[0]/=s;
	v[1]/=s;
	v[2]/=s;
	SphericalAnglesToUnitVector(tht, phi, n);
	ps=n[0]*v[0]+n[1]*v[1]+n[2]*v[2];
	/*fprintf(stdout, "ps = %f\t", ps); */
	if ( angle == 0 ) 
	  sum += gaussianFun(1-fabs(ps), sigma2x2);	
	else 
	{
	  alpha=acos(ps);
	  alpha=(alpha>PI/2) ? alpha-PI : alpha;
	  sum += gaussianFun(alpha, sigma2x2);	
	}
	/*fprintf(stdout, "sum = %f\n", sum); */
  }
  
  for (k=0; k<dim[2]; k++)
  for (j=0; j<dim[1]; j++)
  for (i=0; i<dim[0]; i++) {
	switch (type) {
	case UCHAR:
	  val=(double)theArrayu8[k][j][i];
	  break;
	case FLOAT:
	  val=(double)theArray32[k][j][i];
	  break;
	default:
	  if ( _VT_VERBOSE_ )
        fprintf( stderr, "%s: such image type not handled yet\n", proc );
      return( 0 );
	}
	if (val < 0)
	  continue;
	v[0]=(double) i - c[0];
	v[1]=(double) j - c[1];
	v[2]=(double) k - c[2];
	s=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	v[0]/=s;
	v[1]/=s;
	v[2]/=s;
	SphericalAnglesToUnitVector(tht, phi, n);
	ps=n[0]*v[0]+n[1]*v[1]+n[2]*v[2];
	
	switch (type) {
	case UCHAR:
	  if (angle == 0)
		theArrayu8[k][j][i]+=(unsigned char)(coef*gaussianFun(1-fabs(ps),sigma2x2)/sum);
	  else 
	  {
		alpha = acos(ps);
		alpha = (alpha>PI/2) ? alpha-PI : alpha;
		theArrayu8[k][j][i]+=(unsigned char)(coef*gaussianFun(alpha,sigma2x2)/sum);
	  }
	  break;
	case FLOAT:
	  if (angle == 0)
		theArray32[k][j][i]+=(float)(coef*gaussianFun(1-fabs(ps),sigma2x2)/sum);
	  else 
	  {
		alpha = acos(ps);
		alpha = (alpha>PI/2) ? alpha-PI : alpha;
		theArray32[k][j][i]+=(float)(coef*gaussianFun(alpha,sigma2x2)/sum);
	  }
	  break;
	default:
	  if ( _VT_VERBOSE_ )
        fprintf( stderr, "%s: such image type not handled yet\n", proc );
      return( 0 );
	}
  }
  
  return( 1 );
}

double gaussianFun(double x, double sigma2x2)
{
  double r=exp(-x*x/sigma2x2);
  return(r);
}
