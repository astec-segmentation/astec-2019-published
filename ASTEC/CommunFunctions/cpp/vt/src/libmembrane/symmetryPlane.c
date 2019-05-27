/*************************************************************************
 * symmetryPlane.c -
 *
 * $Id: symmetryPlane.c,v 1.0 2014/03/24 16:44:00 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2014/03/24
 *
 * ADDITIONS, CHANGES
 *
 */

#include <time.h>

#include <vtmalloc.h>

#include <vt_common.h>
#include <vt_symmetryPlane.h>




#define SEUIL 1e-4

typedef struct local_par {

  vt_names names;
  vt_names name_plane;
  double sigma;
  int distribution;
  int bin;
  double angles[2];
  double n[4];
  int markerAngle;
  int markerN;
  int markerSphere;
  int max;
  int plane;
  int iter;
  double sigma_p;
  double sigma_d;
  double deltaL;
  double Lmin;
  int planeEq;
  int init;
  int real;
  int realIn;
} local_par;



/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );


static char *usage = "[image-in] [-trsf %s] [-angles %lf %lf]\n\
\t [-n %lf %lf %lf [%lf [-vox-in]]] -sphere %s [-max %d]] [-sigma %lf] [-bin] \n\
\t [-plane | -plan %s] [-weq %s [-voxel]] [-delta %lf] [-iter %d] [-d %f] [-dmin %lf]\n\
\t [-p %f] [-distribution %s] [-real-d|-voxel-d|-vox-d] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -trsf %s : calcule la transformation permettant d'aligner le plan avec les\n\
\t            axes de l'image\n\
\t -angles %lf %lf : angles de la direction consideree (theta + phi)\n\
\t -n %lf %lf %lf : composantes de la direction consideree\n\
\t -n %lf %lf %lf %lf : fixe le plan de depart pour l'algorithme least square\n\
\t -sigma %lf : ecart-type pour la ponderation selon l'ecart angulaire (defaut : PI/64)\n\
\t -sphere %s : direction deduite du fichier histogramme directionnel %s\n\
\t -max %d : direction extraite la %d-eme plus grande (0=plus grand: default)\n\
\t -bin : aucune ponderation de la distribution\n\
\t -plane %s : calcule la position du plan de symetrie\n\
\t -weq %s : ecrit dans %s l'equation du plan [-voxel : equation en coor. reelles]\n\
\t -iter %d : nombre d'iterations maxi\n\
\t -d %lf : option de ponderation des points par une gaussienne d'e-t %lf \n\
\t          par rapport a la distance au plan de l'iteration precedante \n\
\t -delta %lf : pas de diminution de l'e-t des distances (default: 2.0)\n\
\t -dmin %d : e-t minimal pour la distance au plan pour le least square (default: 5.0)\n\
\t -real-d : distances in real coordinates (default)\n\
\t -vox[el]-d : distances in voxel coordinates\n\
\t -p %lf : option de ponderation des points par une gaussienne d'e-t %lf \n\
\t          par rapport a l'ecart angulaire au plan de l'iteration precedante (default: PI/64)\n\
\t -distribution %s : sauve l'image des distributions de l'angle principal\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  char name[DOUBLESTRINGLENGTH];
  /*   char *temp; */
  local_par par;
  vt_image *imageIn;
  vt_image *imageTht=NULL;
  vt_image *imagePhi=NULL;
  vt_image imres;
  vt_image implane;
  vt_image *theIm;
  int dim[3];
  /*   double pt[3]; */
  bufferType t;
  /*   vt_image* imagesIn[3]; */
  int marker = 0;
  int flag_3D = 1;
  /*   int type = (int)UCHAR; */
  char prefix[STRINGLENGTH];
  int i, j, k;
  float ***theTht=NULL, ***thePhi=NULL;
  float ***theIn32=NULL;
  float ***theOut;
  unsigned char ***thePlane=NULL;
  unsigned char ***theInU8=NULL;
  double n[4], m[3], nold[4];
  double sigma2x2;
  double coef, tht=0, phi=0, alpha;
  
  /*  GM: correction of P declaration */
  double P[3], c[2], tmp;
  double cst;
  double delta;
  double deltaL;
  int ind, ind_max;
  double *H=NULL;
  double E, Eold, L;
  int iter;
  vt_weighted_vector *list;
  int nb;

  vt_fpt siz;

  /*  liste des maxima de l'histogramme directionnel */
  double *nx, *ny, *nz;
  double *values;
  int taille;
  double Lmin;
  
  clock_t start, stop;
  double elapsed;

  start = clock();


  /*--- initialisation des parametres ---*/
  MT_InitParam( &par );

  /*--- lecture des parametres ---*/
  MT_Parse( argc, argv, &par );

  if(par.markerAngle+par.markerN+par.markerSphere != 1)
  {
	MT_ErrorParse("Wrong direction parameters...\n", 0);
  }

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
    MT_ErrorParse("image In type not handled yet\n", 0);		   
  }
	  
  if ( imageIn->dim.z == 1 )
    flag_3D = 0;

  if (par.sigma_p > 0) { 
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
  }

  /*--- Calcul du vecteur normal en input ---*/
  if (_VT_VERBOSE_)
	fprintf(stdout, "Computing direction from input data...\n");
  if (par.markerN==1)
  {
	double sum=sqrt(par.n[0]*par.n[0]+par.n[1]*par.n[1]+par.n[2]*par.n[2]);
	n[0]=par.n[0]/sum;
	n[1]=par.n[1]/sum;
	n[2]=par.n[2]/sum;
	if (par.init == 0) n[3]=par.n[3]/sum;
  }
  if (par.markerAngle==1)
  {
	SphericalAnglesToUnitVector(par.angles[0], par.angles[1], n);
  }
  if(par.markerSphere==1)
  {
	sprintf( name, "%s", par.names.ext );
	vt_image *imageSphere = _VT_Inrimage( name );
	if ( imageSphere == (vt_image*)NULL )
	{
	  VT_FreeImage( imageIn );
	  if (par.sigma_p > 0) { 
		VT_FreeImage( imageTht );
		if (flag_3D == 1)
		  VT_FreeImage( imagePhi );
	  }
	  fprintf(stderr, "%s unreadable\n", name);
	  MT_ErrorParse("unable to read sphere image \n", 0);		   
	}
	dim[0]=imageSphere->dim.x;
	dim[1]=imageSphere->dim.y;
	dim[2]=imageSphere->dim.z;
	
	if (sphereToNormal(imageSphere->array, imageSphere->type, dim, n) != 1)
	{
	  VT_FreeImage( imageSphere );
	  VT_FreeImage( imageIn );
	  if (par.sigma_p > 0) { 
		VT_FreeImage( imageTht );
		if (flag_3D == 1)
		  VT_FreeImage( imagePhi );
	  }
	  fprintf(stderr, "%s unreadable\n", name);
	  MT_ErrorParse("unable to compute normal from sphere image \n", 0);
	}
		  
	if (VT_ExtractMaxima(imageSphere, &nx, &ny, &nz, &values, &taille) != 1)
	{
	  
	  VT_FreeImage( imageSphere );
	  VT_FreeImage( imageIn );
	  if (par.sigma_p > 0) { 
		VT_FreeImage( imageTht );
		if (flag_3D == 1)
		  VT_FreeImage( imagePhi );
	  }
	  MT_ErrorParse("maxima extraction failed \n", 0);		   
	}
	if(_VT_VERBOSE_) {
	  fprintf(stdout, "Diagramme des directions : %d maxima locaux\n",taille);	
	  for (i=0; i<taille;  i++) 
		fprintf(stdout, "normale %d\t = { %f\t %f\t %f }\t\tvalue = %f\n",i,nx[i],ny[i],nz[i],values[i]);
	}
	
	VT_FreeImage( imageSphere );
	
	
	if (par.max != 0) 
	{
	  if(par.max >= taille)
	  {
		VT_FreeImage( imageIn );
		if (par.sigma_p > 0) { 
		  VT_FreeImage( imageTht );
		  if (flag_3D == 1)
			VT_FreeImage( imagePhi );
		}
		MT_ErrorParse("not enough maxima found regarding to parsed parameters\n", 0);		   
	  }
	  n[0]=nx[par.max];
	  n[1]=ny[par.max];
	  n[2]=nz[par.max];
	}
    vtfree(nx);
    vtfree(ny);
    vtfree(nz);
    vtfree(values);
  }
  
  
  /*  Pour la matrice de transformation */
  double dimx=(double)imageIn->dim.x;
  double dimy=(double)imageIn->dim.y;
  double dimz=(double)imageIn->dim.z;

  /*--- Output initialization ---*/
  siz.x=imageIn->siz.x;
  siz.y=imageIn->siz.y;
  siz.z=imageIn->siz.z;

  sprintf( name, "%s", par.names.out );
  t =  FLOAT;
  VT_Image( &imres );
  /* VT_InitVImage( &imres, name, 1, dim[0], dim[1], dim[2], t ); */
  VT_InitVImage( &imres, name, 1, 
		     imageIn->dim.x, imageIn->dim.y, imageIn->dim.z, t );
		     
  imres.siz.x = siz.x;
  imres.siz.y = siz.y;
  imres.siz.z = siz.z;

  if ( VT_AllocImage( &imres ) != 1 ) {
	VT_FreeImage( imageIn );
	if (par.sigma_p > 0) { 
	  VT_FreeImage( imageTht );
	  if ( flag_3D != 0 )
		VT_FreeImage( imagePhi );
	}
	MT_ErrorParse("unable to allocate output image\n", 0);
  }


  theIm = &imres;
  theOut = (float***) theIm->array;

  if ( par.plane == 1)
  {
	sprintf( name, "%s", par.name_plane.ext );
	t =  UCHAR;
	VT_Image( &implane );
	VT_InitVImage( &implane, name, 1, 
		     imageIn->dim.x, imageIn->dim.y, imageIn->dim.z, t );
		     
    implane.siz.x = imageIn->siz.x;
    implane.siz.y = imageIn->siz.y;
    implane.siz.z = imageIn->siz.z;

	if ( VT_AllocImage( &implane ) != 1 ) {
	  VT_FreeImage( &imres );
	  VT_FreeImage( imageIn );
	  if (par.sigma_p > 0) { 
		VT_FreeImage( imageTht );
		if ( flag_3D != 0 )
		  VT_FreeImage( imagePhi );
	  }
	  MT_ErrorParse("unable to allocate output image\n", 0);
	}


	theIm = &implane;
	thePlane = (unsigned char***) theIm->array;
  }

  /*--- calculs ---*/

  /* sigma2x2=2*par.sigma*par.sigma; */
  sigma2x2=2*par.sigma_p*par.sigma_p;
  if (1 || par.plane == 1)
  {
	P[0]=0.0; P[1]=0.0; P[2]=0.0;
	c[0]=-(n[0]*P[0]+n[1]*P[1]+n[2]*P[2]); c[1]=c[0];
	
    P[0]=0; P[1]=0; P[2]=(double)imageIn->dim.z*siz.z;
	tmp=-(n[0]*P[0]+n[1]*P[1]+n[2]*P[2]); 
	if(tmp<c[0]) c[0]=tmp;
	if(tmp>c[1]) c[1]=tmp;

    P[0]=0; P[1]=(double)imageIn->dim.y*siz.y; P[2]=0;
	tmp=-(n[0]*P[0]+n[1]*P[1]+n[2]*P[2]); 
	if(tmp<c[0]) c[0]=tmp;
	if(tmp>c[1]) c[1]=tmp;

    P[0]=(double)imageIn->dim.x*siz.x; P[1]=0; P[2]=0;
	tmp=-(n[0]*P[0]+n[1]*P[1]+n[2]*P[2]); 
	if(tmp<c[0]) c[0]=tmp;
	if(tmp>c[1]) c[1]=tmp;

    P[0]=0; P[1]=(double)imageIn->dim.y*siz.y; P[2]=(double)imageIn->dim.z*siz.z;
	tmp=-(n[0]*P[0]+n[1]*P[1]+n[2]*P[2]); 
	if(tmp<c[0]) c[0]=tmp;
	if(tmp>c[1]) c[1]=tmp;

    P[0]=(double)imageIn->dim.x*siz.x; P[1]=0; P[2]=(double)imageIn->dim.z*siz.z;
	tmp=-(n[0]*P[0]+n[1]*P[1]+n[2]*P[2]); 
	if(tmp<c[0]) c[0]=tmp;
	if(tmp>c[1]) c[1]=tmp;

    P[0]=(double)imageIn->dim.x*siz.x; P[1]=(double)imageIn->dim.y*siz.y; P[2]=0;
	tmp=-(n[0]*P[0]+n[1]*P[1]+n[2]*P[2]); 
	if(tmp<c[0]) c[0]=tmp;
	if(tmp>c[1]) c[1]=tmp;

    P[0]=(double)imageIn->dim.x*siz.x; P[1]=(double)imageIn->dim.y*siz.y; P[2]=(double)imageIn->dim.z*siz.z;
	tmp=-(n[0]*P[0]+n[1]*P[1]+n[2]*P[2]); 
	if(tmp<c[0]) c[0]=tmp;
	if(tmp>c[1]) c[1]=tmp;
	
    delta=(siz.x<siz.y) ? ((siz.x<siz.z) ? 0.5*siz.x : 0.5*siz.z) : ((siz.y<siz.z) ? 0.5*siz.y : 0.5*siz.z);
	
    ind_max=(int)(1.5+(c[1]-c[0])/delta);
	if (par.init == 1) {
      H = vtmalloc( sizeof(double)*ind_max, "H", argv[0] );
	  if(H==(double*)NULL) 
	  {
		VT_FreeImage( &imres );
		if (par.plane == 1) 
		  VT_FreeImage( &implane );
		VT_FreeImage( imageIn );
		if (par.sigma_p > 0) { 
		  VT_FreeImage( imageTht );
		  if ( flag_3D != 0 )
			VT_FreeImage( imagePhi );
	    }
		MT_ErrorParse("unable to allocate Histogram vector\n", 0);
	  }
	  for (ind=0;ind<ind_max;ind++) H[ind]=0;
	}
  }
  
  if ( flag_3D == 1 ) {
    /*  3D case */
    for (k=0;k<(int)imageIn->dim.z;k++)
    for (j=0;j<(int)imageIn->dim.y;j++)
    for (i=0;i<(int)imageIn->dim.x;i++)
	{
	  if(par.plane==1)
	  {
		thePlane[k][j][i]=(unsigned char)0;
	  }
      if(par.distribution == 1 || par.init == 1)   {
		switch (imageIn->type) {
		case UCHAR:
	      coef=(double) (theInU8[k][j][i]); 
	      break;
		case FLOAT:
	      coef=(double) (theIn32[k][j][i]); 
	      break;
		default:
		  VT_FreeImage( &imres );
		  if (par.plane == 1) 
		    VT_FreeImage( &implane );
		  VT_FreeImage( imageIn );
		  if (par.sigma_p > 0) { 
			VT_FreeImage( imageTht );
			VT_FreeImage( imagePhi );
		  }
	      MT_ErrorParse( "input image type not handled yet\n", 0 );
		}
		if (par.sigma_p > 0)
		{
		  tht=theTht[k][j][i];
		  phi=thePhi[k][j][i];
		}
		if (par.bin == 1) 
		  coef = (coef > 0.0) ? 1.0 : 0.0;
		if ( coef > 0.0 )/* && addToSphere(theIm->array, dim, t, par.rayon, pt, par.angle, par.sigma, coef, tht, phi) == 0 )  */
		{
		  if (par.sigma_p > 0) {
			SphericalAnglesToUnitVector(tht, phi, m);
			alpha=acos(n[0]*m[0]+n[1]*m[1]+n[2]*m[2]);
			alpha=(alpha>PI/2) ? alpha-PI : alpha;
			theOut[k][j][i]=coef*gaussianFun(alpha,sigma2x2);
		  }
		  else
			theOut[k][j][i]=coef;
		  /* VT_FreeImage( &imres ); */
		  /* VT_FreeImage( imageIn ); */
		  /* VT_FreeImage( imageTht ); */
		  /* VT_FreeImage( imagePhi ); */
		  /* MT_ErrorParse( "Error while computing distribution\n", 0 ); */
        }
        if(par.init == 1)   {
          cst=-(i*siz.x*n[0]+j*siz.y*n[1]+k*siz.z*n[2]);
		  ind=(int)(0.5+(cst-c[0])/delta);
		  H[ind]+=theOut[k][j][i];
		}
	  }
	  
	}
	if(par.init == 1 )
	{
	  /* fprintf(stdout, "\n\n"); */
      for (ind=1; ind<ind_max; ind++) {
	/* fprintf(stdout, "%f ", H[ind]); */
		H[ind]+=H[ind-1];
	  }
      /* fprintf(stdout, "\n\n"); */
      ind=0;
      while(ind<ind_max && H[ind]<0.5*H[ind_max-1]) ind++;
      if(ind==ind_max)
      {
		free(H);
		H=NULL;
		VT_FreeImage( &imres );
		VT_FreeImage( &implane );
		VT_FreeImage( imageIn );
		if (par.sigma_p > 0) { 
		  VT_FreeImage( imageTht  );
		  if (flag_3D ==1)
			VT_FreeImage( imagePhi  );
	    }
		MT_ErrorParse("Unexpected error while searching plane position\n", 0);
	  }
	  vtfree(H);
	  H=NULL;

	  n[3]=c[0]+delta*ind;
	  /* n[3]=-215; */
	  
	  
	  /* if (par.plane==1) */
	  /* for (k=0;k<imageIn->dim.z;k++) */
	  /* for (j=0;j<imageIn->dim.y;j++) */
	  /* for (i=0;i<imageIn->dim.x;i++)  */
	  /* { */
	  /* if(fabs(n[3]+i*n[0]+j*n[1]+k*n[2])<1) */
	  /* thePlane[k][j][i]=(unsigned char)1; */
	  /* } */
	  
	}
	if (par.plane==1) {
	  int tmp=1;
	  if(VT_AddPlaneToImage( &implane, n, &tmp ) != 1)
	  {
		VT_FreeImage( &implane );
		VT_FreeImage( imageIn );
		if (par.sigma_p > 0) { 
		  VT_FreeImage( imageTht  );
		  VT_FreeImage( imagePhi  );
		}
		VT_FreeImage( &imres  );
		VT_FreeWeightedVector( list);
		MT_ErrorParse("unable to write plane image\n", 0);
	  }
	}
	if(_VT_VERBOSE_ )
      fprintf(stdout, "Plane #%d:\t a x + b y + c z + d = 0 (real coordinates)\n\t a = %f\n\t b = %f\n\t c = %f\n\t d = %f\n\n",
					  0, n[0], n[1], n[2], n[3]);

  }
  else {
    /*  2D case : TODO */
	VT_FreeImage( imageIn );
	if (par.sigma_p > 0) { 
	  VT_FreeImage( imageTht  );
	}
	VT_FreeImage( &imres  );
	if (par.plane == 1) 
	  VT_FreeImage( &implane );
	MT_ErrorParse( "2D case not handled yet\n", 0 );
  }

  if( par.distribution==1 && par.sigma_d >= par.Lmin)
  {
    /*  GM: correction of dim declaration */
      int dim[3];
      dim[0]=imageIn->dim.x;
      dim[1]=imageIn->dim.y;
      dim[2]=imageIn->dim.z;
      VT_DistributionCoef(theOut, dim, siz, n, par.sigma_d);
  }

  if(1)
  {
	if(flag_3D==1)
	{
      double _siz = (siz.x<siz.y) ? ((siz.x<siz.z) ? siz.x : siz.z) : ((siz.y<siz.z) ? siz.y : siz.z);
      if (par.realIn == 0) _siz=1;
      L=par.sigma_d/_siz; /*  constant value ? */
      deltaL = par.deltaL/_siz;
      Lmin=par.Lmin/_siz;
	  int cpt=2;


      do
	  {
		if(L<Lmin) continue;
		
		if(_VT_VERBOSE_)
		  fprintf(stdout, "\nPlane with band of semi-width %f...\n\n", L);
	  
		if (par.sigma_p>0) {
		  if (VT_AllocateAndBuildWeightedVector( imageIn, imageTht, imagePhi, 	
                                3*L*_siz, n, &list,  &nb) != 1)
		  {
			VT_FreeImage( &implane );
			VT_FreeImage( imageIn );
			if (par.sigma_p > 0) { 
			  VT_FreeImage( imageTht  );
			  VT_FreeImage( imagePhi  );
			}
			VT_FreeImage( &imres  );

			MT_ErrorParse("unable to write plane image\n", 0);
		  }
	    }
	    else {
		  if (VT_AllocateAndBuildWeightedVectorNoAngle( imageIn, 
                                3*L*_siz, n, &list,  &nb) != 1)
		  {
			VT_FreeImage( &implane );
			VT_FreeImage( imageIn );
			if (par.sigma_p > 0) { 
			  VT_FreeImage( imageTht  );
			  VT_FreeImage( imagePhi  );
			}
			VT_FreeImage( &imres  );
			MT_ErrorParse("unable to write plane image\n", 0);
		  }
		}

        if (_VT_VERBOSE_)
        {
            fprintf(stdout,"      <vt_weighted_vector> list size = %d\n", nb);
        }

		iter = 1;
		if(_VT_VERBOSE_ && 0)
		  fprintf(stdout, "Plane %d:\t a x + b y + c z + d = 0\n\t a = %f\n\t b = %f\n\t c = %f\n\t d = %f\n\n",
				iter-1, n[0], n[1], n[2], n[3]);
	  
        E=VT_MeasureSquareError(list, nb, n, L*_siz, par.sigma);
        if (_VT_VERBOSE_)
        {
            fprintf(stdout,"      <vt_weighted_vector> list square error = %f\n", E);
        }

		do
		{
		  Eold=E;
          int _i;
          for (_i=0;_i<4;_i++)nold[_i]=n[_i];
	  /*  plan courant P_i : n[0]*x+n[1]*y+n[2]*z+cst=0 */

				
		  if(_VT_VERBOSE_ && 0)
			fprintf(stdout, "Total least square iteration # %d\n",iter);
		
		  iter++;

		
          if (VT_totalLeastSquarePlane(list,nb, L*_siz, par.sigma_p, n, &E) != 1)
		  {
			VT_FreeImage( imageIn );
			if (par.sigma_p > 0) { 
			  VT_FreeImage( imageTht  );
			  if (flag_3D ==1)
				VT_FreeImage( imagePhi  );
			}
			VT_FreeImage( &imres  );
			
			if (par.plane == 1 ) 
			  VT_FreeImage( &implane );
			VT_FreeWeightedVector( list);
			MT_ErrorParse( "unable to compute the total least square solution\n", 0 );
		  }
		  if(_VT_VERBOSE_ && 0)
            fprintf(stdout, "Plane %d:\t a x + b y + c z + d = 0 (real coordinates)\n\t a = %f\n\t b = %f\n\t c = %f\n\t d = %f\n\t Error = %f\n\t Ratio = %f\n\n",
			  iter-1, n[0], n[1], n[2], n[3], E, E/Eold);
			
        } while(iter<=par.iter && (fabs((nold[3]-n[3]))>SEUIL ||
                                   sqrt((nold[0]-n[0])*(nold[0]-n[0])+(nold[1]-n[1])*(nold[1]-n[1])+(nold[2]-n[2])*(nold[2]-n[2]))>SEUIL));
		/* while(iter<=par.iter && fabs(1-(E/Eold))>SEUIL); */
		

		VT_FreeWeightedVector( list);

		if (par.plane != 0) {
		  if(VT_AddPlaneToImage( &implane, n, &cpt ) != 1)
		  {
			VT_FreeImage( &implane );
			VT_FreeImage( imageIn );
			if (par.sigma_p > 0) { 
			  VT_FreeImage( imageTht  );
			  VT_FreeImage( imagePhi  );
			}
			VT_FreeImage( &imres  );
			VT_FreeWeightedVector( list);
            MT_ErrorParse("unable to write plane image\n", 0);
		  }
		}
		
	    L-=deltaL;

   		if (_VT_VERBOSE_ || L<Lmin)
         fprintf(stdout,"Plane #%d:\t a x + b y + c z + d = 0 (real coordinates)\n\t a = %f\n\t b = %f\n\t c = %f\n\t d = %f\n\t convergence after %d iterations \n\n",
			cpt-1, n[0], n[1], n[2], n[3], iter);

	    cpt++;

      } while (L>=Lmin);

	} 
	else
	{
	  /*  2D case : TODO */
	  VT_FreeImage( imageIn );
	  VT_FreeImage( imageTht  );
	  VT_FreeImage( &imres  );
	  if (par.plane == 1 ) 
		VT_FreeImage( &implane );
      MT_ErrorParse( "2D case not handled yet\n", 0 );
	}
  }




  /*--- liberations memoires ---*/
  VT_FreeImage( imageIn );
  if (par.sigma_p > 0) { 
    VT_FreeImage( imageTht  );
    if (flag_3D ==1)
	  VT_FreeImage( imagePhi  );
  }
  
  /*--- ecriture de l'image resultat ---*/
  if ( par.distribution == 1) {
	if ( VT_WriteInrimage( &imres ) == -1 ) {
      VT_FreeImage( &imres );
	  if (par.plane == 1 ) 
		VT_FreeImage( &implane );
      if(_VT_VERBOSE_)
		MT_ErrorParse("unable to write output image\n", 0);
	}
  }
  VT_FreeImage( &imres );

  if (par.plane == 1 ) {
/*	for (k=0;k<imageIn->dim.z;k++)
	for (j=0;j<imageIn->dim.y;j++)
	for (i=0;i<imageIn->dim.x;i++) 
	{
	if(fabs(n[3]+i*n[0]+j*n[1]+k*n[2])<=1)
	  thePlane[k][j][i]=(unsigned char)255;
	else
	  thePlane[k][j][i]=(unsigned char)0;
  }
*/
	if ( VT_WriteInrimage( &implane ) == -1 ) {
      VT_FreeImage( &implane );
	  if(_VT_VERBOSE_)
		MT_ErrorParse("unable to write plane image\n", 0);
	}
	VT_FreeImage( &implane );


  }
  
  /*--- Calcul de la transformation pour l'alignement du plan avec l'axe x ---*/
  
  /*  Matrice de rotation : */
  double R[9];	/*  Matrice de rotation de nouvelles vers anciennes coordonnees */
  double uy=-n[2], uz=n[1], C, s; /*  ux=0, vecteur u = vecteur de rotation */
  s=sqrt(uy*uy+uz*uz);
  if(s<1e-6)
  {
	R[0]=R[4]=R[8]=1.0;
	R[1]=R[2]=R[3]=R[5]=R[6]=R[7]=0.0;
  }
  else {
	uy/=s;
	uz/=s;
	C=n[0];
	s=sqrt(1-C*C);
	/*
	R=[                                                       ...
	ux^2+(1-ux^2)*c   ux*uy*(1-c)-uz*s    ux*uz*(1-c)+uy*s  ; ...
	ux*uy*(1-c)+uz*s  uy^2+(1-uy^2)*c     uy*uz*(1-c)-ux*s  ; ...
	ux*uz*(1-c)-uy*s  uy*uz*(1-c)+ux*s    uz^2+(1-uz^2)*c   ; ...
	];
	*/
	R[0]=C;
	R[1]=-uz*s;
	R[2]=uy*s;
	R[3]=uz*s;
	R[4]=uy*uy+(1-uy*uy)*C;
	R[5]=uy*uz*(1-C);
	R[6]=-uy*s;
	R[7]=uy*uz*(1-C);
	R[8]=uz*uz+(1-uz*uz)*C;  
  }
  
  /*  Liste des huit coins de l'image d'origine */
  double O[8][3];
  double rmin[3], rmax[3], TEMP;
  double trans[3];
  
  O[0][0]=0;		O[0][1]=0;		O[0][2]=0;
  O[1][0]=dimx-1;	O[1][1]=0;		O[1][2]=0;
  O[2][0]=0;		O[2][1]=dimy-1;	O[2][2]=0;
  O[3][0]=0;		O[3][1]=0;		O[3][2]=dimz-1;
  O[4][0]=dimx-1;	O[4][1]=dimy-1;	O[4][2]=0;
  O[5][0]=dimx-1;	O[5][1]=0;		O[5][2]=dimz-1;
  O[6][0]=0;		O[6][1]=dimy-1;	O[6][2]=dimz-1;
  O[7][0]=dimx-1;	O[7][1]=dimy-1;	O[7][2]=dimz-1;

  for (i=0; i<8; i++) {
      O[i][0]*=siz.x;
      O[i][1]*=siz.y;
      O[i][2]*=siz.z;
  }

  for (j=0; j<3; j++) {
    /* rmin[j]=RR[j]*O[0][0]+RR[1+j*3]*O[0][1]+RR[2+j*3]*O[0][2]; */
    /* rmax[j]=rmin[j]; */
	rmin[j]=0;
	rmax[j]=0;
  }
  
  for (i=1;i<8;i++)
  {
	for (j=0; j<3; j++) {
	  TEMP=R[j]*O[i][0]+R[3+j]*O[i][1]+R[6+j]*O[i][2];
	  rmin[j]=(rmin[j]<=TEMP) ? rmin[j] : TEMP;
	  rmax[j]=(rmax[j]>=TEMP) ? rmax[j] : TEMP;
	}
  }
  
  
  /*  Calcul des nouvelles dimensions */
  dimx=(int)((rmax[0]-rmin[0])/siz.x)+1;
  dimy=(int)((rmax[1]-rmin[1])/siz.y)+1;
  dimz=(int)((rmax[2]-rmin[2])/siz.z)+1;
  
  /*  Calcul de la translation inverse */
  for (i = 0; i<3; i++)
	trans[i]=R[i*3]*rmin[0]+R[1+i*3]*rmin[1]+R[2+i*3]*rmin[2];
  
  /*  Calcul du nouveau plan */
  int indmax=0;
  double nmax=n[0];
  double newx;
  for (i=1; i<3; i++)
	if (fabs(n[i])>fabs(nmax)) {
	  nmax=n[i]; 
	  indmax=i;
	}
  
  if (indmax==0)
	newx=R[0]*(-n[3]/nmax-trans[0])-R[3]*trans[1]-R[6]*trans[2];
  if (indmax==1)
	newx=-R[0]*trans[0]+R[3]*(-n[3]/nmax-trans[1])-R[6]*trans[2];
  if (indmax==2)
	newx=-R[0]*trans[0]-R[3]*trans[1]+R[6]*(-n[3]/nmax-trans[2]);
  
  
  
  /*  Ecriture des informations */
  if (_VT_VERBOSE_)
  {
    fprintf(stdout, "Plane equation");
    if (par.real == 1)
    {
        fprintf(stdout, " (real coordinates):\n\n");
        fprintf(stdout, "n: \t%f \t%f \t%f \t%f\n",
          n[0], n[1], n[2], n[3]);
    }
    /* else */
    {
        fprintf(stdout, " (voxel coordinates):\n\n");
        double tmp[3];
        tmp[0]=n[0]*siz.x;
        tmp[1]=n[1]*siz.y;
        tmp[2]=n[2]*siz.z;
        double somme=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
        fprintf(stdout, "n: \t%f \t%f \t%f \t%f\n",
          tmp[0]/somme, tmp[1]/somme, tmp[2]/somme, n[3]/somme);
    }
  }
  if (par.planeEq == 1) {
	FILE* fichier = NULL;
	fichier = fopen(par.name_plane.in, "w");
	if( fichier == NULL)
	{
	  MT_ErrorParse("Erreur pour l'ecriture de l'equation dans un fichier\n", 0);
	}
	fprintf(fichier, "#\n# \tPlane equation:\n#\n");
    if (par.real == 1)
        fprintf(fichier, "n: \t%f \t%f \t%f \t%f\n",
		  n[0], n[1], n[2], n[3]);
    else {
        double tmp[3];
        tmp[0]=n[0]*siz.x;
        tmp[1]=n[1]*siz.y;
        tmp[2]=n[2]*siz.z;
        double somme=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
        fprintf(fichier, "n: \t%f \t%f \t%f \t%f\n",
          tmp[0]/somme, tmp[1]/somme, tmp[2]/somme, n[3]/somme);
    }
	fclose(fichier);
  }
  
  
  
  if (_VT_VERBOSE_) {
	fprintf(stdout, "\nTransformation matrix:\n\n");
	fprintf(stdout, "%f \t%f \t%f \t%f\n%f \t%f \t%f \t%f\n%f \t%f \t%f \t%f\n0 \t0 \t0 \t1\n",
		  R[0], R[1], R[2], trans[0], 
		  R[3], R[4], R[5], trans[1], 
		  R[6], R[7], R[8], trans[2]);
	fprintf(stdout, "\nNew dimensions:\n\n");
	fprintf(stdout, "\n%d \t%d \t%d\n", (int)dimx+1, (int)dimy+1, (int)dimz+1);
    fprintf(stdout, "\nNew plane position:\n\n");
    fprintf(stdout, "\nX = %f\n\n", newx);
    fprintf(stdout, "\nNew voxel size:\n\n");
    fprintf(stdout, "\nvs = %f %f %f\n\n", siz.x, siz.y, siz.z);
  }
  
  if (par.name_plane.out[0] != '\0') {
	FILE* fichier = NULL;
	fichier = fopen(par.name_plane.out, "w");
	if( fichier == NULL)
	{
	  MT_ErrorParse("Erreur pour l'ecriture de la transformation dans un fichier\n", 0);
	}
	fprintf(fichier, "#\n# \tTransformation matrix:\n#\n");
	fprintf(fichier, "%f \t%f \t%f \t%f\n%f \t%f \t%f \t%f\n%f \t%f \t%f \t%f\n0 \t0 \t0 \t1\n",
		  R[0], R[1], R[2], trans[0], 
		  R[3], R[4], R[5], trans[1], 
		  R[6], R[7], R[8], trans[2]);
	fprintf(fichier, "#\n# \tNew dimensions: \t");
	fprintf(fichier, "\t%d \t%d \t%d\n", (int)dimx+1, (int)dimy+1, (int)dimz+1);
    fprintf(fichier, "#\n# \tNew plane position in X: \t");
    fprintf(fichier, "\t%f\n#\n", newx);
    fprintf(fichier, "#\n# \tNew voxel size: \t");
    fprintf(fichier, "\t%f %f %f\n#\n", siz.x, siz.y, siz.z);

	fclose(fichier);
  }
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

      /* Parametres de calcul */

      else if ( strcmp ( argv[i], "-trsf" ) == 0 || 
				strcmp ( argv[i], "-transformation" ) == 0 ) {
        i +=1;
        if (i >=argc)	   MT_ErrorParse( "parsing -trsf...\n", 0);
        strncpy( par->name_plane.out, argv[i], STRINGLENGTH );
      }

      else if ( strcmp ( argv[i], "-distribution" ) == 0 ) {
		  par->distribution = 1;
        i +=1;
        if (i >=argc)	   MT_ErrorParse( "parsing -distribution...\n", 0);
        par->distribution=1;
        strncpy( par->names.out, argv[i], STRINGLENGTH );
      }

      else if ( strcmp ( argv[i], "-plane" ) == 0 || 
				strcmp ( argv[i], "-plan" ) == 0 ) {
        par->plane = 1;
        i +=1;
        if (i >=argc)	   MT_ErrorParse( "parsing -plane...\n", 0);
        strncpy( par->name_plane.ext, argv[i], STRINGLENGTH );
      }

      else if ( strcmp ( argv[i], "-weq" ) == 0 ) {
        par->planeEq = 1;
        i +=1;
        if (i >=argc)	   MT_ErrorParse( "parsing -weq...\n", 0);
        strncpy( par->name_plane.in, argv[i], STRINGLENGTH );
      }

      else if ( strcmp ( argv[i], "-voxel" ) == 0 ) {
        par->real = 0;
      }

      else if ( strcmp ( argv[i], "-voxel-d" ) == 0 || strcmp ( argv[i], "-vox-d" ) == 0 ) {
        par->realIn = 0;
      }

      else if ( strcmp ( argv[i], "-real-d" ) == 0 ) {
        par->realIn = 1;
      }

      else if ( strcmp ( argv[i], "-delta" ) == 0 ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -delta...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->deltaL) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -delta...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-p" ) == 0 ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -p...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->sigma_p) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -p...\n", 0 );
      }
	  
      else if ( strcmp ( argv[i], "-d" ) == 0 ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -d...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->sigma_d) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -d...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-dmin" ) == 0 ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -dmin...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->Lmin) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -dmin...\n", 0 );
      }


      else if ( strcmp ( argv[i], "-sigma" ) == 0 ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -sigma...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->sigma) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -sigma...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-iter" ) == 0 ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -iter...\n", 0 );
		status = sscanf( argv[i],"%d",&(par->iter) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -iter...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-bin" ) == 0 ) {
        par->bin = 1;
      }

      else if ( strcmp ( argv[i], "-angles" ) == 0 ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -angles...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->angles[0]) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -angles...\n", 0 );
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -angles : 2 angles expected...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->angles[1]) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -angles : 2 angles expected...\n", 0 );
		par->markerAngle=1;
      }

      else if ( strcmp ( argv[i], "-n" ) == 0 ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -n...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->n[0]) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -n...\n", 0 );
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -n (2nd arg)...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->n[1]) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -n (2nd arg)...\n", 0 );
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -n (3rd arg)...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->n[2]) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -n (3rd arg)...\n", 0 );
		par->markerN=1;
		i += 1;
		if ( i >= argc)    i -= 1;
		status = sscanf( argv[i],"%lf",&(par->n[3]) );
		if ( status <= 0 ) i -= 1;
		else par->init = 0;

      }

      else if ( strcmp ( argv[i], "-sphere" ) == 0 ||
		(strcmp ( argv[i], "-sphere" ) == 0 ) ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -sphere...\n", 0 );
		strncpy( par->names.ext, argv[i], STRINGLENGTH );
		par->markerSphere=1;
      }

      else if ( strcmp ( argv[i], "-max" ) == 0 ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -max...\n", 0 );
		status = sscanf( argv[i],"%d",&(par->max) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -max...\n", 0 );
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
  VT_Names( &(par->name_plane) );

  par->distribution = 0;
  par->angles[0]=0;
  par->angles[1]=0;
  par->markerAngle=0;
  par->n[0]=0;
  par->n[1]=0;
  par->n[2]=0;
  par->markerN=0;
  par->markerSphere=0;
  par->sigma = PI/64;
  par->bin = 0;
  par->plane = 0;
  par->iter=100;
  /*par->sigma_p=0.0;*/
  par->sigma_p=PI/64;
  par->sigma_d=10*0.3;
  par->Lmin=5*0.3;
  par->deltaL = 2.0*0.3;
  par->max=0;
  par->planeEq=0;
  par->init = 1;
  par->realIn = 1;
  par->real = 1;
}



