/*
 * vt_directionDistribution.c -
 *
 * $Id: vt_directionDistribution.c,v 1.0 2014/07/23 11:03:34 gael Exp $
 *
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2014/07/23
 *
 * ADDITIONS, CHANGES
 *
 *
 */


#include <convert.h>
#include <regionalmax.h>
#include <vtmalloc.h>


#include <vt_common.h>

#include <vt_directionDistribution.h>



static int _verbose_ = 0;


static double EPS=1e-16;


int VT_AllocateAndBuildWeightedVector(	vt_image* imageIn,	/* nuage de points */
					vt_image* imageTht, 	/* angle spherique theta */
					vt_image* imagePhi, 	/* angle spherique phi */
					double L, 		/* demi-largeur du sous-espace */
					double *n, 		/* parametres du plan initial */
					vt_weighted_vector **image,	/* liste a allouer */
					int *nb	/* taille de la liste */
								)
{
  char *proc = "VT_AllocateAndBuildWeightedVector";
  
  int N;
  int i,j,k;
  int dimx, dimy, dimz;
  
  unsigned char ***theU8, *bufU8;
  float ***theFLT, *bufFLT, ***theTht, ***thePhi;
  double coef; 
  
  vt_weighted_vector *list;
  vt_weighted_vector *theVector;
  
  
  dimx=imageIn->dim.x;
  dimy=imageIn->dim.y;
  dimz=imageIn->dim.z;
  
  switch(imageIn->type) {
  case UCHAR : 
    theU8 = (unsigned char ***)imageIn->array;
    bufU8 = (unsigned char *)imageIn->buf;
	break;
  case FLOAT : 
    theFLT = (float ***)imageIn->array;
    bufFLT = (float *)imageIn->buf;
	break;
  default : 
	fprintf(stderr, "%s : such image type not handled yet\n",proc);
	return(0);
  }

  switch(imageTht->type) {
  case FLOAT : 
    theTht = (float ***)imageTht->array;
	break;
  default : 
	fprintf(stderr, "%s : such image type not handled yet\n",proc);
	return(0);
  }

  switch(imagePhi->type) {
  case FLOAT : 
    thePhi = (float ***)imagePhi->array;
	break;
  default : 
	fprintf(stderr, "%s : such image type not handled yet\n",proc);
	return(0);
  }

  N=0;
  for (i = 0; i < dimx*dimy*dimz ; i++) 
  {
    switch ( imageIn->type) {
    case UCHAR :
      if( bufU8[i] > (unsigned char) 0) N++;
      break;
    case FLOAT:
	  if( bufFLT[i]> 0.0 ) N++;
	  break;
	default:
	  fprintf(stderr, "%s : such image type not handled yet\n",proc);
	  return(0);
    }
  }
  
  list = vtmalloc( N*sizeof(vt_weighted_vector), "list", proc );
  if (list == (vt_weighted_vector*) NULL)
  {
	fprintf(stderr, "%s : such image type not handled yet\n",proc);
	return(0);
  }
  /* fprintf(stdout, "ALLOC: list = %d", list); */
  N=0;
  for (k = 0; k<dimz ; k++)
  for (j = 0; j<dimy ; j++)
  for (i = 0; i<dimx ; i++)
  {
	if (fabs(n[0]*i+n[1]*j+n[2]*k+n[3])>L) continue;
	switch (imageIn->type) {
	  case UCHAR : 
		coef=(double)theU8[k][j][i];
		break;
	  case FLOAT :
		coef=(double)theFLT[k][j][i];
		break;
	  default :
		return(0);
	}
	if (coef < EPS) continue;	/* restriction au nuage de points */
	theVector=list+N;
	theVector->x = i;
	theVector->y = j;
	theVector->z = k;
	theVector->coef = coef;
	theVector->weight = 0;
	SphericalAnglesToUnitVector(theTht[k][j][i], thePhi[k][j][i], theVector->n);
	N++;
  }
  list= vtrealloc( list, N*sizeof(vt_weighted_vector), "list", proc );
  *nb=N;
  
  /* for (N=0;N<3; N++) { */
  /*   theVector=list+N; */
  /*   fprintf(stdout, "liste[%d] : {%d, %d, %d}\n", N, theVector->x, theVector->y, theVector->z); */
  /*   fprintf(stdout, "            {%f, %f, %f}\n", theVector->n[0],theVector->n[1], theVector->n[2]); */
  /* } */
  *image=list;
  return(1);
}

void VT_FreeWeightedVector(vt_weighted_vector *v)
{
  vtfree(v);
  v=(vt_weighted_vector *)NULL;
  return;
}





int VT_totalLeastSquarePlane(	vt_weighted_vector* image, 	/* nuage de points */
				int nb, 			/* nb de points du nuage */
				double L, 			/* demi-largeur du sous-espace */
				double sigma,		/* ecart-type pour la ponderation du nuage en fct des ecarts d'angles */
				double *n, 			/* parametres du plan initial */
				double *E)			/* residu du least square */
{
  char *proc = "VT_totalLeastSquarePlane";
  /* fprintf(stdout, "ABCD: list: %d\n", image); */

  /* fprintf(stdout, "nb = %d\n", nb); */
  
  /* vt_image temp; */
  /*   float ***theTEMP; */
  vt_weighted_vector *vec;

  int i;
  /*   unsigned char ***theU8; */
  /*   float ***theFLT; */
  /*   float ***theTht; */
  /*   float ***thePhi; */
  
  double coef;
  double sum_coef;

  double A[3];
  
  double M[9];
  double valprop[3], vecprop[9];
  
  double norm;

  double *m;
  double sigma2x2=2*sigma*sigma;
  double L2x2=2*L*L;
  
  double alpha, a, b, c;
  
  norm=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  if (norm < EPS)
  {
	fprintf(stderr, "%s: unexpected norm value of n[0..2]\n",proc);
	return(0);
  }
  n[0]/=norm;
  n[1]/=norm;
  n[2]/=norm;
  n[3]/=norm;
  
  
/*
  switch (imageIn->type) {
  case UCHAR :
	theU8=(unsigned char ***)imageIn->array;
	break;
  case FLOAT : 
	theFLT=(float ***)imageIn->array;
	break;  
  default:
	fprintf(stderr, "%s: such imageIn type not handled yet\n", proc);
    return(0);
  }

  switch (imageTht->type) {
  case FLOAT : 
	theTht=(float ***)imageTht->array;
	break;  
  default:
	fprintf(stderr, "%s: such imageTht type not handled yet\n", proc);
    return(0);
  }

  switch (imagePhi->type) {
  case FLOAT : 
	thePhi=(float ***)imagePhi->array;
	break;  
  default:
	fprintf(stderr, "%s: such imagePhi type not handled yet\n", proc);
    return(0);
  }
*/
  
/*   char *name="TEMP.inr"; */
/*   VT_Image( &temp ); */
/*   VT_InitVImage( &temp, name, 1,  */
/* 		     imageIn->dim.x, imageIn->dim.y, imageIn->dim.z, FLOAT ); */
		     
/*   temp.siz.x = 1.0; */
/*   temp.siz.y = 1.0; */
/*   temp.siz.z = 1.0; */

/*   if ( VT_AllocImage( &temp ) != 1 ) { */
/* 	fprintf(stderr,"%s: unable to allocate temp image\n", proc); */
/* 	return(0); */
/*   } */

/*   theTEMP=(float***)temp.array; */


  /*--- Calcul du point A barycentre du nuage de points ---*/
  A[0]=0.0; A[1]=0.0; A[2]=0.0;
  sum_coef=0.0;


  for (i=0; i<nb ; i++)
  {  
	
	vec=image+i;
	
	if (fabs(n[0]*vec->x+n[1]*vec->y+n[2]*vec->z+n[3])>2.45*L) /*  restriction au sous espace */
	{
	  vec->weight = -1;
	  continue; 
	}
	coef=vec->coef;

	m = vec->n;

	alpha=acos(n[0]*m[0]+n[1]*m[1]+n[2]*m[2]);
	alpha=(alpha>PI/2) ? alpha-PI : alpha;
	coef*=gaussianFun(alpha,sigma2x2);
	coef*=gaussianFun(n[0]*vec->x+n[1]*vec->y+n[2]*vec->z+n[3], L2x2);

	
	sum_coef+=coef;
	A[0]+=coef*vec->x;
	A[1]+=coef*vec->y;
	A[2]+=coef*vec->z;
	
	vec->weight=coef;
	
  }
  A[0]/=sum_coef;
  A[1]/=sum_coef;
  A[2]/=sum_coef;
  
  if(_verbose_)
	fprintf(stdout, "A = [ %f,\t %f,\t %f ]\n",A[0], A[1], A[2]);

  /*--- Calcul de la matrice M de la forme quadratique a minimiser ---*/
  M[0]=0.0; M[1]=0.0; M[2]=0.0; M[4]=0.0; M[5]=0.0; M[8]=0.0;
  

  for(i=0;i<nb;i++)
  {

	vec=image+i;
	m=vec->n;

	coef=vec->weight;
	
	if(coef < 0) continue;
	
	a=vec->x-A[0];
	b=vec->y-A[1];
	c=vec->z-A[2];
	
	M[0]+=coef*a*a;
	M[1]+=coef*a*b;
	M[2]+=coef*a*c;
	M[4]+=coef*b*b;
	M[5]+=coef*b*c;
	M[8]+=coef*c*c;
  }

  M[3]=M[1];
  M[6]=M[2];
  M[7]=M[5];


  /*--- Calcul de l'eigensystem de la matrice M ---*/
  if ( _ComputeEigensOfSymetricSquareMatrix( M, valprop, vecprop, 3 )
       != 1 ) {
    if (_verbose_)
      fprintf(stderr, "%s: error in computing eigens\n", proc);
    return( -1 );
  }
  if ( _SortEigensInAbsIncreasingOrder( valprop, vecprop, 3 ) != 1 ) {

    if (_verbose_)
      fprintf(stderr, "%s: error in sorting eigens\n", proc);
    return( -1 );
  }

  /*--- Ecriture des resultats ---*/
  *E=valprop[0];
  n[0]=vecprop[0]; 
  n[1]=vecprop[3]; 
  n[2]=vecprop[6]; 
  n[3]=-(n[0]*A[0]+n[1]*A[1]+n[2]*A[2]);
  if(_verbose_)
	fprintf(stdout, "N = [ %f,\t %f,\t %f ]\t\tc = %f\t\tE = %f\n", 
			n[0],n[1],n[2],n[3],*E);
  
  /* VT_WriteInrimage( &temp ); */
  /* VT_FreeImage(&temp); */
  return(1);
}

double VT_MeasureSquareError(vt_weighted_vector *v, int nb, double *n, double L, double sigma_p)
{
  double E=0;
  double L2x2=2*L*L;
  double sigma2x2=2*sigma_p*sigma_p;
  int i;
  double *m;
  double alpha;
  double coef;
  vt_weighted_vector *vec;
  
  for (i=0 ; i<nb ; i++)
  {
	vec = v+i;

	if (fabs(n[0]*vec->x+n[1]*vec->y+n[2]*vec->z+n[3])>2.45*L) continue;

	m = vec->n;

	alpha=acos(n[0]*m[0]+n[1]*m[1]+n[2]*m[2]);
	alpha=(alpha>PI/2) ? alpha-PI : alpha;
	coef=vec->coef;
	coef*=gaussianFun(alpha,sigma2x2);
	coef*=gaussianFun(n[0]*vec->x+n[1]*vec->y+n[2]*vec->z+n[3], L2x2);

	E+=coef;
  }
  
  return (E);
}

double gaussianFun(double x, double sigma2x2)
{
  if (sigma2x2<EPS) return(x);
  return(exp(-x*x/sigma2x2));
}

int sphereToNormal(void*** array, bufferType t, int* dim, double* n )
{
  char *proc = "sphereToNormal";

  int i,j,k;
  unsigned char ***theU8;
  float ***theFLOAT;
  int mid;
  
  double max=0;
  double val;
  
  switch (t) {
  case UCHAR:
    theU8 = (unsigned char ***)array;
    break;
  case FLOAT:
    theFLOAT = (float ***)array;  
    break;
  default:
    if ( _VT_VERBOSE_ )
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    return( 0 );
  }
  
  mid=(int)(0.5*(dim[0]-1));
  
  for (k=0;k<dim[2];k++)
  for (j=0;j<dim[1];j++)
  for (i=0;i<dim[0];i++)
  {
	switch (t) {
	case UCHAR: 
	  val=(double)theU8[k][j][i];
	  break;
	case FLOAT:
	  val=(double)theFLOAT[k][j][i];
	  break;
	default:
	  if ( _VT_VERBOSE_ )
		fprintf( stderr, "%s: such image type not handled yet\n", proc );
	return( 0 );
	}
	if (val<=max)
	  continue;
	n[0]=i-mid;
	n[1]=j-mid;
	n[2]=k-mid;
	max=val;
  }
  val=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  n[0]=n[0]/val;
  n[1]=n[1]/val;
  n[2]=n[2]/val;
  if (_VT_VERBOSE_)
	fprintf(stdout, "Normal vector: \nnx\t%f\nny\t%f\nnz\t%f\n", n[0],n[1],n[2]);
  return(1);
}



int VT_AddPlaneToImage(vt_image* img, double *n, void * value)
{
  char *proc="VT_AddPlaneToImage";
  int i,j,k;
  unsigned char ***theU8, *vU8;
  unsigned short int ***theUSHORT, *vUSHORT;
  float ***theFLT, *vFLT;
  
  switch (img->type) {
  case UCHAR :
	theU8=(unsigned char ***)img->array;
    vU8 = (unsigned char *)value;
    /* fprintf(stdout, "img->type = UCHAR\t\tvalue=%d\n", (int)(*vU8)); */
    break;
  case USHORT :
    theUSHORT=(unsigned short int ***)img->array;
    vUSHORT = (unsigned short int *)value;
    /* fprintf(stdout, "img->type = USHORT\t\tvalue=%d\n", (int)(*vUSHORT)); */
    break;
  case FLOAT :
    theFLT=(float ***)img->array;
    vFLT = (float *)value;
    /* fprintf(stdout, "img->type = FLOAT\t\tvalue=%f\n", (*vFLT)); */
    break;
  default:
    fprintf(stderr, "%s : such image type not handled yet\n", proc);
    return(0);
  }
  
  for (k=0;k<img->dim.z;k++)
  for (j=0;j<img->dim.y;j++)
  for (i=0;i<img->dim.x;i++) 
  {
	if(fabs(n[3]+i*n[0]+j*n[1]+k*n[2])<=1)
	  switch (img->type) {
	  case UCHAR :
		theU8[k][j][i] = *vU8;
		break;
	  case USHORT :
		theUSHORT[k][j][i] = *vUSHORT;
		break;
	  case FLOAT :
		theFLT[k][j][i] = *vFLT;
		break;
	  default :
		fprintf(stderr, "%s : such image type not handled yet\n", proc);
	    return(0);
	  }
  }
  return(1);
}


int VT_ExtractMaxima(vt_image* histo, double** normalx, double** normaly, double** normalz, double **values, int *N)
{
  char *proc="VT_ExtractMaxima";
  char name[256];
  double *nx, *ny, *nz, coef;
  double *v;
  int nb, i,j,k, ind, IND, dim[3];
  vt_image imres;
  double val;
  unsigned char ***theU8, ***theV;
  float ***theFLT;
  unsigned short int ***theUSHORT;
  
     
  dim[0] = histo->dim.x;
  dim[1] = histo->dim.y;
  dim[2] = histo->dim.z;

  /*--- initialisation de l'image resultat ---*/
  VT_Image( &imres );
  sprintf(name, "histo_norma.inr");
  VT_InitFromImage( &imres, histo, name, UCHAR );
  
  if ( VT_AllocImage( &imres ) != 1 ) {
    fprintf(stderr, "%s : unable to allocate the normalized histogram image\n", proc);
  }

  /*--- normalisation + maxima regionaux ---*/
  switch (histo->type) {
  case FLOAT :
	theFLT=(float ***)histo->array;
    if ( VT_NormaImage( histo, &imres ) != 1 ) {
      VT_FreeImage( &imres );
      fprintf(stderr, "%s : unable to normalize the histogram\n", proc);
    }
    regionalmax( imres.buf, imres.buf, UCHAR, dim, 1, 1 );
	break;
  case USHORT :
	theUSHORT=(unsigned short int ***) histo->array;
    if ( VT_NormaImage( histo, &imres ) != 1 ) {
      VT_FreeImage( &imres );
      fprintf(stderr, "%s : unable to normalize the histogram\n", proc);
    }
    regionalmax( imres.buf, imres.buf, UCHAR, dim, 1, 1 );
	break;
  case UCHAR :
	theU8=(unsigned char***) histo->array;
    regionalmax( histo->buf, imres.buf, UCHAR, dim, 1, 1 );
	break;
  default :
	VT_FreeImage(&imres);
	fprintf(stderr, "%s : such image type not handled yet\n", proc);
	return(0);
  }
  theV = (unsigned char ***) imres.array;

  /* fprintf(stdout, "Normalisation et regionalmax effectues (r = %d)...\n", r); */
  
  /*--- Recuperation des maxima ---*/
  nb=0;
  for (k = 0; k<=(int)(0.5*dim[2]) ; k++)
  for ( j=0; j<dim[1]; j++)
  for (i=0; i<dim[0]; i++)
  {
    /*  restriction aux directions de la demi sphere */
	if(k == (int)(0.5*dim[2]) ) {
	  if (j < (int)(0.5*dim[1])) continue;
	  if (j == (int)(0.5*dim[1]) && i < (int)(0.5*dim[0])) continue;
	}

	if (theV[k][j][i] == (unsigned char) 0) continue ;

	switch (histo->type) {
	case UCHAR: 
	  val=(double)theU8[k][j][i];
	  break;
	case USHORT :
	  val=(double)theUSHORT[k][j][i];
	  break;
	case FLOAT :
	  val=(double)theFLT[k][j][i];
	  break;
	default:
	  VT_FreeImage(&imres);
	  if (nb>0) {
		free(v);
		free(nx);
		free(ny);
		free(nz);
	  }
	  fprintf(stderr, "%s: such image type not handled yet\n", proc);
	  return(0);
	}
	  
	ind=0; 
	while(ind<nb && v[ind]>val) ind++;
	nb++;
	if (nb==1)
	{
          v = vtmalloc( nb*sizeof(double), "v", proc );
          nx = vtmalloc( nb*sizeof(double), "nx", proc );
          ny = vtmalloc( nb*sizeof(double), "ny", proc );
          nz = vtmalloc( nb*sizeof(double), "nz", proc );
	}
	else
	{
          v= vtrealloc( v, nb*sizeof(double), "v", proc );
          nx= vtrealloc( nx, nb*sizeof(double), "nx", proc );
          ny= vtrealloc( ny, nb*sizeof(double), "ny", proc );
          nz= vtrealloc( nz, nb*sizeof(double), "nz", proc );
	}
	if (v == NULL || nx == NULL || ny == NULL || nz == NULL)
	{
	  VT_FreeImage(&imres);
	  if (v != NULL) vtfree(v);
	  if (nx != NULL) vtfree(nx);
	  if (ny != NULL) vtfree(ny);
	  if (nz != NULL) vtfree(nz);
          fprintf(stderr, "%s : unable to re-allocate a vector\n", proc);
	  return(0);
	}

	IND=nb-1;
	while(IND>ind) 
	{
	  v[IND]=v[IND-1];
	  nx[IND]=nx[IND-1];
	  ny[IND]=ny[IND-1];
	  nz[IND]=nz[IND-1];
	  IND--;
	}
	v[ind]=val;
	coef=sqrt((double)((i-(int)(dim[0]/2))*(i-(int)(dim[0]/2)) + 
				(j-(int)(dim[1]/2))*(j-(int)(dim[1]/2)) + 
				(k-(int)(dim[2]/2))*(k-(int)(dim[2]/2))));
	nx[ind]=((double)(i-(int)(dim[0]/2)))/coef;
	ny[ind]=((double)(j-(int)(dim[1]/2)))/coef;
	nz[ind]=((double)(k-(int)(dim[2]/2)))/coef;
  }

  *N=nb;
  *values=v;
  *normalx=nx;
  *normaly=ny;
  *normalz=nz;

  return(1);
}
