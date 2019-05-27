/*
 * mt_membane2D.c -
 *
 * $Id: mt_membrane2D.c,v 1.0 2013/06/20 11:28:34 gael Exp $
 *
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 * 
 * CREATION DATE: 
 * 2013/06/20
 *
 * ADDITIONS, CHANGES
 *
 *
 */
 
#include <eigens.h>
#include <vtmalloc.h>

#include <mt_membrane2D.h>


static int _borderLength_[3] = { 10, 10, 10 };
static int _verbose_ = 1;


#define FLTZERO 1e-8




/* rand_(a) retourne une valeur entiere aleatoire comprise entre 0 et a-1
 */
int rand_(int a){
    return(rand()%a);
}





void MT_SetVerboseInMtMembrane2D( int v )
{
  _verbose_ = v;
}

void MT_IncrementVerboseInMtMembrane2D(  )
{
  _verbose_ ++;
}

void MT_DecrementVerboseInMtMembrane2D(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}




/************************************************************
 *
 * tensors
 *
 *
 ************************************************************/


int  VT_Alloc2DTensorFromImage( mt_2Dtensor *par, vt_image *im,
          char *genericname )
{
  char name[256];

  sprintf( name, "%s.imxx.inr", genericname );
  VT_InitFromImage( &(par->imxx), im, name,  (int)FLOAT );
  sprintf( name, "%s.imxy.inr", genericname );
  VT_InitFromImage( &(par->imxy), im, name,  (int)FLOAT );
  sprintf( name, "%s.imyy.inr", genericname );
  VT_InitFromImage( &(par->imyy), im, name,  (int)FLOAT );


  sprintf( name, "%s.imvp1.inr", genericname );
  VT_InitFromImage( &(par->imvp1),  im, name,  (int)FLOAT );
  sprintf( name, "%s.imvp2.inr", genericname );
  VT_InitFromImage( &(par->imvp2),  im, name,  (int)FLOAT );

  sprintf( name, "%s.imtheta1.inr", genericname );
  VT_InitFromImage( &(par->imtheta1),  im, name,  (int)FLOAT );
  sprintf( name, "%s.imtheta2.inr", genericname );
  VT_InitFromImage( &(par->imtheta2),  im, name,  (int)FLOAT );

  sprintf( name, "%s.iszero.inr", genericname );
  VT_InitFromImage( &(par->iszero),  im, name,  (int)UCHAR );

  if ( VT_AllocImage( &(par->imxx) ) != 1 ) {
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imxy) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imyy) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imvp1) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imvp2) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imvp1) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imtheta1) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imtheta2) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    VT_FreeImage( &(par->imtheta1) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->iszero) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    VT_FreeImage( &(par->imtheta1) );
    VT_FreeImage( &(par->imtheta2) );
    return( -1 );
  }


  return( 1 );


}

int  VT_Alloc2DTensor( mt_2Dtensor *par,
    char *genericname /* image generic name */,
        int dimx /* X dimension */,
        int dimy /* Y dimension */,
        int dimz /* Z dimension */,
        int type /* image type  */ )
{
  char name[256];

  sprintf( name, "%s.imxx.inr", genericname );
  VT_InitVImage( &(par->imxx), name, (int)1, dimx, dimy, dimz, type );
  sprintf( name, "%s.imxy.inr", genericname );
  VT_InitVImage( &(par->imxy), name, (int)1, dimx, dimy, dimz, type );
  sprintf( name, "%s.imyy.inr", genericname );
  VT_InitVImage( &(par->imyy), name, (int)1, dimx, dimy, dimz, type );



  sprintf( name, "%s.imvp1.inr", genericname );
  VT_InitVImage( &(par->imvp1), name, (int)1, dimx, dimy, dimz, type );
  sprintf( name, "%s.imvp2.inr", genericname );
  VT_InitVImage( &(par->imvp2), name, (int)1, dimx, dimy, dimz, type );

  sprintf( name, "%s.imtheta1.inr", genericname );
  VT_InitVImage( &(par->imtheta1), name, (int)1, dimx, dimy, dimz, type );
  sprintf( name, "%s.imtheta2.inr", genericname );
  VT_InitVImage( &(par->imtheta2), name, (int)1, dimx, dimy, dimz, type );

  sprintf( name, "%s.iszero.inr", genericname );
  VT_InitVImage( &(par->iszero), name, (int)1, dimx, dimy, dimz, (int)UCHAR );

  if ( VT_AllocImage( &(par->imxx) ) != 1 ) {
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imxy) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imyy) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imvp1) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imvp2) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imvp1) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imtheta1) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imtheta2) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    VT_FreeImage( &(par->imtheta1) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->iszero) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    VT_FreeImage( &(par->imtheta1) );
    VT_FreeImage( &(par->imtheta2) );
    return( -1 );
  }

  return( 1 );

}

void VT_Free2DTensor ( mt_2Dtensor *par )
{
  VT_FreeImage( &(par->imxx) );
  VT_FreeImage( &(par->imxy) );
  VT_FreeImage( &(par->imyy) );
  VT_FreeImage( &(par->imvp1) );
  VT_FreeImage( &(par->imvp2) );
  VT_FreeImage( &(par->imtheta1) );
  VT_FreeImage( &(par->imtheta2) );
  VT_FreeImage( &(par->iszero) );
}

void VT_Write2DTensor( mt_2Dtensor *par )
{
  VT_WriteInrimage( &(par->imxx) );
  VT_WriteInrimage( &(par->imxy) );
  VT_WriteInrimage( &(par->imyy) );
  VT_WriteInrimage( &(par->imvp1) );
  VT_WriteInrimage( &(par->imvp2) );
  VT_WriteInrimage( &(par->imtheta1) );
  VT_WriteInrimage( &(par->imtheta2) );
  VT_WriteInrimage( &(par->iszero) );

}





int VT_Write2DtensorWithNames( mt_2Dtensor *par,
                               char *tensorname,
                               char *eigenvaluename,
                               char *anglename,
                               char *binaryname,
                               char *suffix )
{
  char *proc = "VT_Write2DtensorWithNames";
  char name[STRINGLENGTH];
  char *ptrSuffix;
  char *inrSuffix = "inr";

  ptrSuffix = ( suffix != (char*)NULL && suffix[0] != '\0' ) ? suffix : inrSuffix;

  if ( tensorname != (char*)NULL && tensorname[0] != '\0' ) {
    sprintf( name, "%s.imxx.%s", tensorname, ptrSuffix );
    if ( VT_WriteInrimageWithName( &(par->imxx), name ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to write imxx image\n", proc );
      return( -1 );
    }
    sprintf( name, "%s.imyy.%s", tensorname, ptrSuffix );
    if ( VT_WriteInrimageWithName( &(par->imyy), name ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to write imyy image\n", proc );
      return( -1 );
    }
    sprintf( name, "%s.imxy.%s", tensorname, ptrSuffix );
    if ( VT_WriteInrimageWithName( &(par->imxy), name ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to write imxy image\n", proc );
      return( -1 );
    }
  }

  if ( eigenvaluename != (char*)NULL && eigenvaluename[0] != '\0' ) {
    sprintf( name, "%s.imvp1.%s", eigenvaluename, ptrSuffix );
    if ( VT_WriteInrimageWithName( &(par->imvp1), name ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to write imvp1 image\n", proc );
      return( -1 );
    }
    sprintf( name, "%s.imvp2.%s", eigenvaluename, ptrSuffix );
    if ( VT_WriteInrimageWithName( &(par->imvp2), name ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to write imvp2 image\n", proc );
      return( -1 );
    }
  }

  if ( anglename != (char*)NULL && anglename[0] != '\0' ) {
    sprintf( name, "%s.imtheta1.%s", anglename, ptrSuffix );
    if ( VT_WriteInrimageWithName( &(par->imtheta1), name ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to write imtheta1 image\n", proc );
      return( -1 );
    }
    sprintf( name, "%s.imtheta2.%s", anglename, ptrSuffix );
    if ( VT_WriteInrimageWithName( &(par->imtheta2), name ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to write imtheta2 image\n", proc );
      return( -1 );
    }
  }

  if ( binaryname != (char*)NULL && binaryname[0] != '\0' ) {
    sprintf( name, "%s.iszero.%s", binaryname, ptrSuffix );
    if ( VT_WriteInrimageWithName( &(par->iszero), name ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to write iszero image\n", proc );
      return( -1 );
    }
  }

  return( 1 );
}





int VT_Write2DTensorWithName( mt_2Dtensor *par, char *genericname)
{
    char *proc = "VT_Write2DTensorWithName";
    if ( VT_Write2DtensorWithNames( par, genericname, genericname, genericname, genericname, "inr" ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to write images\n", proc );
        return( -1 );
    }
    return( 1 );
}











/************************************************************
 *
 *
 *
 *
 ************************************************************/


static void _init2DTensorFromAngle(mt_2Dtensor *theTensor)
{
  int x,y,z;
  int dimx = theTensor->imxx.dim.x;
  int dimy = theTensor->imxx.dim.y;
  int dimz = theTensor->imxx.dim.z;

  float ***theXX = (float ***)theTensor->imxx.array;
  float ***theYY = (float ***)theTensor->imyy.array;
  float ***theXY = (float ***)theTensor->imxy.array;
  float ***theTheta2 = (float ***)theTensor->imtheta2.array;
  float ***theVP2 = (float ***)theTensor->imvp2.array;

  unsigned char ***zeros = (unsigned char ***)theTensor->iszero.array;

  double v2[2], theta;

  for (z=0; z<dimz; z++)
  for (y=0; y<dimy; y++)
  for (x=0; x<dimx; x++)
  {
    if (zeros[z][y][x] == 1)
      continue;

    theta = (double)theTheta2[z][y][x];
    v2[0]=cos(theta); v2[1]=sin(theta);

    theXX[z][y][x] = (float)v2[0]*v2[0]*theVP2[z][y][x];
    theYY[z][y][x] = (float)v2[1]*v2[1]*theVP2[z][y][x];
    theXY[z][y][x] = (float)v2[0]*v2[1]*theVP2[z][y][x];
  }
  return;
}







/************************************************************
 *
 * angles
 *
 *
 ************************************************************/


static int _compute2DAngles(double **angles, int Nangles)
{
  char *proc = "_compute2DAngles";
  int i;
  double pi = 3.14159265358979323846;


  double theta;

  if( Nangles <0)
  {
    if (_verbose_)
      fprintf(stderr, "%s: error : Nangles should be positive integer\n", proc);
    return(-1);
  }


  (*angles) = vtmalloc(sizeof(double)*Nangles, "(*angles)", proc );
  if((*angles)==NULL)
  {
    if (_verbose_)
      fprintf(stderr, "%s: error in allocating angles\n", proc);
    return(-1);
  }

  for (i=0; i<Nangles; i++)
  {
    theta = i*pi/Nangles;
    (*angles)[i] = theta;
  }

  return( Nangles );
}





static int _nearest2DAngle(double theta, double *angles, int Nangles)
{
  int n, N=0;
  double dot, DOT=0;

  double v[2], v0[2];
  v0[0]=cos(theta); v0[1]=sin(theta);

  for (n=0;n<Nangles;n++)
  {
    v[0]=cos(angles[n]); v[1]=sin(angles[n]);

    dot = fabs(v[0]*v0[0]+v[1]*v0[1]);
    if (dot<DOT)
      continue;
    DOT = dot;
    N = n;
    if (1-DOT < FLTZERO)
      break;
  }

  return (N);
}






/************************************************************
 *
 * static procedures: filtering
 *
 *
 ************************************************************/



static int _filter2Dimages( vt_image *im, vt_2Dimages *par, float *theCoeffs )
{
  char *proc = "_filter2Dimages";
  filterType theFilter = GAUSSIAN_DERICHE;
  derivativeOrder theDerivatives[3] = { NODERIVATIVE, NODERIVATIVE, NODERIVATIVE };
  int theDim[3];

  theDim[0] = im->dim.x;
  theDim[1] = im->dim.y;
  theDim[2] = im->dim.z;

  {
    size_t x, y, z;
    s8 ***theE = (s8 ***)par->ime.array;
    float ***theR = (float ***)par->imr.array;

    for ( z=0; z<par->ime.dim.z; z++ )
    for ( y=0; y<par->ime.dim.y; y++ )
    for ( x=0; x<par->ime.dim.x; x++ ) {
      theE[z][y][x] = 0;
      theR[z][y][x] = 0.0;
    }
  }


  theDerivatives[0] = NODERIVATIVE;
  theDerivatives[1] = DERIVATIVE_0;

  if ( RecursiveFilterOnBuffer( im->buf, im->type,
                                par->imx.buf, par->imx.type,
                                theDim, _borderLength_,
                                theDerivatives, theCoeffs, theFilter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computing derivatives (.,0)\n", proc );
    return( -1 );
  }
  theDerivatives[1] = NODERIVATIVE;
  theDerivatives[0] = DERIVATIVE_2;
  if ( RecursiveFilterOnBuffer( par->imx.buf, par->imx.type,
                                par->hessien.imxx.buf, par->hessien.imxx.type,
                                theDim, _borderLength_,
                                theDerivatives, theCoeffs, theFilter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computing derivatives (2,.)\n", proc );
    return( -1 );
  }
  theDerivatives[0] = DERIVATIVE_1;
  if ( RecursiveFilterOnBuffer( par->imx.buf, par->imx.type,
                                par->imx.buf, par->imx.type,
                                theDim, _borderLength_,
                                theDerivatives, theCoeffs, theFilter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computing derivatives (1,.)\n", proc );
    return( -1 );
  }




  theDerivatives[0] = DERIVATIVE_0;
  theDerivatives[1] = NODERIVATIVE;

  if ( RecursiveFilterOnBuffer( im->buf, im->type,
                                par->imy.buf, par->imy.type,
                                theDim, _borderLength_,
                                theDerivatives, theCoeffs, theFilter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computing derivatives (0,.)\n", proc );
    return( -1 );
  }
  theDerivatives[0] = NODERIVATIVE;
  theDerivatives[1] = DERIVATIVE_2;
  if ( RecursiveFilterOnBuffer( par->imy.buf, par->imy.type,
                                par->hessien.imyy.buf, par->hessien.imyy.type,
                                theDim, _borderLength_,
                                theDerivatives, theCoeffs, theFilter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computing derivatives (.,2)\n", proc );
    return( -1 );
  }
  theDerivatives[1] = DERIVATIVE_1;
  if ( RecursiveFilterOnBuffer( par->imy.buf, par->imy.type,
                                par->imy.buf, par->imy.type,
                                theDim, _borderLength_,
                                theDerivatives, theCoeffs, theFilter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computing derivatives (.,1)\n", proc );
    return( -1 );
  }

  theDerivatives[0] = DERIVATIVE_1;
  theDerivatives[1] = DERIVATIVE_1;

  if ( RecursiveFilterOnBuffer( im->buf, im->type,
                                par->hessien.imxy.buf, par->hessien.imxy.type,
                                theDim, _borderLength_,
                                theDerivatives, theCoeffs, theFilter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computing derivatives (1,1)\n", proc );
    return( -1 );
  }

  return( 1 );
}





static void _compute2DResponse( vt_2Dimages *par, enumStructureColor color,
                           double tau, double sigma )
{
  double coeff = tau;
  size_t x, y, z;

  s8 ***theE = (s8 ***)par->ime.array;
  float ***theR = (float ***)par->imr.array;
  float ***theTHETA1 = (float ***)par->hessien.imtheta1.array;

  double dx, dy, gradx, grady;
  double ps1, ps2;

  if ( coeff <= 0.0 ) coeff = sqrt(3.0);




  for ( z=0; z<par->imx.dim.z; z++ )
  for ( y=0; y<par->imx.dim.y; y++ )
  for ( x=0; x<par->imx.dim.x; x++ ) {
    if ( theE[z][y][x] < 0 ) {
      theR[z][y][x] = 0.0;
      continue;
    }
    /* on va en M + c * sigma * V
       blanc : ps1 < 0
       noir  : ps1 > 0
     */
    dx = x + coeff * sigma * cos( theTHETA1[z][y][x] );
    dy = y + coeff * sigma * sin( theTHETA1[z][y][x] );
    gradx = _GetInterpolated2DValue( &(par->imx), dx, dy, z );
    grady = _GetInterpolated2DValue( &(par->imy), dx, dy, z );
    ps1 = cos( theTHETA1[z][y][x] ) * gradx + sin( theTHETA1[z][y][x] ) * grady;


    /* on va en M - c * sigma * V
       blanc : ps2 > 0
       noir  : ps2 < 0
     */
    dx = x - coeff * sigma * cos( theTHETA1[z][y][x] );
    dy = y - coeff * sigma * sin( theTHETA1[z][y][x] );
    gradx = _GetInterpolated2DValue( &(par->imx), dx, dy, z );
    grady = _GetInterpolated2DValue( &(par->imy), dx, dy, z );
    ps2 = cos( theTHETA1[z][y][x] ) * gradx + sin( theTHETA1[z][y][x] ) * grady;



    /* on cherche une reponse positive
     */
    switch( color ) {
    case _WHITE_ :
      /*
        theR[z][y][x] = ps2 - ps1;
      */
      if ( ps2 > (-ps1) ) {
        theR[z][y][x] = -ps1;
      } else {
        theR[z][y][x] = ps2;
      }
      break;
    case _BLACK_ :
      /*
        theR[z][y][x] = ps1 - ps2;
      */
       if ( ps1 > (-ps2) ) {
        theR[z][y][x] = -ps2;
      } else {
        theR[z][y][x] = ps1;
      }

      break;
    }

    if ( theR[z][y][x] < 0.0 ) {
      theR[z][y][x] = 0.0;
      theE[z][y][x] = -128;
    }

    if ( theR[z][y][x] < 0 ) printf( "%f ", theR[z][y][x] );

  }
}










/************************************************************
 *
 * static procedures: mt_2Dtensor
 *
 *
 ************************************************************/


static int _set2DTensor(mt_2Dtensor *tensorImg, int x, int y, int slice,
          enumVote mode)
{
  char *proc = "_set2DTensor";
  float ***theXX = (float***)(tensorImg->imxx.array);
  float ***theYY = (float***)(tensorImg->imyy.array);
  float ***theXY = (float***)(tensorImg->imxy.array);

  float ***theVP1 = (float***)(tensorImg->imvp1.array);
  float ***theVP2 = (float***)(tensorImg->imvp2.array);

  float ***theTht1 = (float***)(tensorImg->imtheta1.array);
  float ***theTht2 = (float***)(tensorImg->imtheta2.array);

  unsigned char ***theZeros = (unsigned char ***) (tensorImg->iszero.array);

  double hessien[4];
  double valprop[2], vecprop[4];
  double v[2], theta;

  float oldVPsum, newVPsum;

  hessien[0] = (double)theXX[slice][y][x];
  hessien[1] = hessien[2] = (double)theXY[slice][y][x];
  hessien[3] = (double)theYY[slice][y][x];

  if ( _ComputeEigensOfSymetricSquareMatrix( hessien, valprop, vecprop, 2 )
       != 1 ) {
    if (_verbose_)
      fprintf(stderr, "%s: error in computing eigens\n", proc);
    return( -1 );
  }
  if ( _SortEigensInAbsIncreasingOrder( valprop, vecprop, 2 ) != 1 ) {

    if (_verbose_)
      fprintf(stderr, "%s: error in sorting eigens\n", proc);
    return( -1 );
  }

  /* Les valeurs et vecteurs propres sont tries par fabs croissante */
  if (mode == SPARSE) {
    oldVPsum=theVP1[slice][y][x]+theVP2[slice][y][x];
    newVPsum=(float)(valprop[0]+valprop[1]);
  }

  theVP1[slice][y][x]=(float)valprop[0];
  v[0]=vecprop[0]; v[1]=vecprop[2];
  if (fabs(v[0])<FLTZERO)
    theta = asin(v[1]);
  else
    theta = atan(v[1]/v[0]);

  theTht1[slice][y][x]=(float)theta;

  theVP2[slice][y][x]=(float)valprop[1];
  v[0]=vecprop[1]; v[1]=vecprop[3];
  if (fabs(v[0])<FLTZERO)
    theta = asin(v[1]);
  else
    theta = atan(v[1]/v[0]);
  theTht2[slice][y][x]=(float)theta;

  if (valprop[1]>FLTZERO)
  {
    theZeros[slice][y][x] = (unsigned char)0;
    if (mode == SPARSE)
    {
        /* On "normalise" : sum(vp avant vote)=sum(vp apres vote) */
        theVP1[slice][y][x]=theVP1[slice][y][x]*oldVPsum/newVPsum;
        theVP2[slice][y][x]=theVP2[slice][y][x]*oldVPsum/newVPsum;

        theXX[slice][y][x]=theXX[slice][y][x]*oldVPsum/newVPsum;
        theYY[slice][y][x]=theYY[slice][y][x]*oldVPsum/newVPsum;
        theXY[slice][y][x]=theXY[slice][y][x]*oldVPsum/newVPsum;
    }

  }
  else
  {
      theZeros[slice][y][x] = (unsigned char)1;
  }

  return(1);
}










/************************************************************
 *
 *
 *
 *
 ************************************************************/



int MT_Compute2DMultiScale( vt_image *theIm,
                            vt_3Dimres *imsRes,
                            vt_image *mask,
                            int flagMask,
                            double scale1,
                            double scale2,
                            int nbscales,
                            enumStructureColor color )
{
  char *proc = "MT_Compute2DMultiScale";
  int s, nbs = nbscales;
  double scale;

  vt_2Dimages theIms;

  float theCoeff[3];

  size_t x, y, z;
  float ***maxRep = (float***)imsRes->imRep.array;
  float ***maxTht = (float***)imsRes->imTheta.array;
  float ***maxScl = (float***)imsRes->imScale.array;

  float ***theRep;
  float ***theTht;

  float tau = sqrt(3);


  unsigned char ***theMaskU8 = (unsigned char ***)NULL;
  unsigned short int ***theMaskU16 = (unsigned short int ***)NULL;



  if ( VT_Alloc2Dimages( &theIms, theIm, proc ) != 1 ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to allocate 2D auxiliary images\n", proc );
      return( -1 );
    }
  }

  theRep = (float***)theIms.imr.array;
  theTht = (float***)theIms.hessien.imtheta1.array;



  if ( scale2 <= scale1 ) nbs = 1;

  if (flagMask == 1) {
      switch (mask->type) {
      default:
      case UCHAR:
      case SCHAR:
          theMaskU8 = (unsigned char ***) mask->array;
          break;
      case USHORT:
      case SSHORT:
          theMaskU16 = (unsigned short int ***) mask->array;
          break;
      }
  }

  for ( s = 0; s < nbs; s ++ ) {

    if ( s == 0 ) {
      scale = scale1;
    } else {
      scale = exp( ((double)(nbs-1-s) * log(scale1) + (double)(s) * log(scale2)) /
                   (double)(nbs-1) );
    }

    if ( _verbose_ ) {
      fprintf( stderr, "... processing scale #%d, sigma = %f\n", s, scale );
    }



    theCoeff[0] = theCoeff[1] = theCoeff[2] = scale;


    if ( _filter2Dimages( theIm, &theIms, theCoeff ) != 1 ) {
      VT_Free2Dimages( &theIms );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to 2D filter\n", proc );
            return( -1 );

    }

    VT_Compute2DEigenVectors( &(theIms.hessien), &(theIms.ime) );
    VT_FilterOn2DEigenValues( &theIms, color );
    _compute2DResponse( &theIms, color, tau, scale );

    /* on recopie le resultat dans les images resultats
     */

    if ( s == 0 ) {
      for ( z = 0; z < theIm->dim.z; z ++ )
      for ( y = 0; y < theIm->dim.y; y ++ )
      for ( x = 0; x < theIm->dim.x; x ++ ) {
          if (flagMask == 1)
              switch (mask->type) {
              default:
              case UCHAR:
              case SCHAR:
                  if (theMaskU8[z][y][x] == (unsigned char) 0) {
                      theRep[z][y][x] = 0.0;
                  }
                  break;
              case USHORT:
              case SSHORT:
                  if (theMaskU16[z][y][x] == (unsigned short int) 0) {
                      theRep[z][y][x] = 0.0;
                  }
                  break;
              }
        if ( theRep[z][y][x] > 0.0 ) {
          maxRep[z][y][x] = scale * theRep[z][y][x];
          maxTht[z][y][x] = theTht[z][y][x];
          maxScl[z][y][x] = scale * theIm->siz.x;
        } else {
          maxRep[z][y][x] = maxTht[z][y][x] = 0.0;
          maxScl[z][y][x] = 0.0;
        }
      }
    } else {
      for ( z = 0; z < theIm->dim.z; z ++ )
      for ( y = 0; y < theIm->dim.y; y ++ )
      for ( x = 0; x < theIm->dim.x; x ++ ) {
        if ( scale * theRep[z][y][x] > maxRep[z][y][x] ) {
          maxRep[z][y][x] = scale * theRep[z][y][x];
          maxTht[z][y][x] = theTht[z][y][x];
          maxScl[z][y][x] = scale * theIm->siz.x;
        }
      }
    }

  } /* for ( s = 0; s < nbs; s ++ ) */

  VT_Free2Dimages( &theIms );
  return( 1 );
}





void MT_Compute2DExtrema( vt_3Dimres *imsRes, vt_image *mask, int flagMask,
                          vt_image *imExt )
{
  float ***theExt = (float***)(imExt->array);
  float ***theRep = (float***)(imsRes->imRep.array);
  float ***theTht = (float***)(imsRes->imTheta.array);

  unsigned char ***theMaskU8 = (unsigned char ***)NULL;
  unsigned short int ***theMaskU16 = (unsigned short int ***)NULL;


  size_t x, y, z;
  double rep, r1, r2, theta, dx, dy;


  if (flagMask == 1)
  {
      switch (mask->type) {
      default:
      case UCHAR:
      case SCHAR:
          theMaskU8 = (unsigned char ***) mask->array;
          break;
      case USHORT:
      case SSHORT:
          theMaskU16 = (unsigned short int***) mask->array;
          break;
      }
  }

  for ( z=0; z<imExt->dim.z; z++ )
  for ( y=0; y<imExt->dim.y; y++ )
  for ( x=0; x<imExt->dim.x; x++ ) {

      if (flagMask == 1)
      {
          switch (mask->type) {
          default:
          case UCHAR:
          case SCHAR:
              if (theMaskU8[z][y][x] == (unsigned char) 0 ||
                      (x+1 < imExt->dim.x && theMaskU8[z][y][x+1] == (unsigned char) 0) ||
                      (x >= 1 && theMaskU8[z][y][x-1] == (unsigned char) 0) ||
                      (y+1 < imExt->dim.y && theMaskU8[z][y+1][x] == (unsigned char) 0) ||
                      (y >= 1 && theMaskU8[z][y-1][x] == (unsigned char) 0) ||
                      (x >= 1 && y >= 1 && theMaskU8[z][y-1][x-1] == (unsigned char) 0) || (x >= 1 && y+1 < imExt->dim.y && theMaskU8[z][y+1][x-1] == (unsigned char) 0) ||
                      (x+1 < imExt->dim.x && y >= 1 && theMaskU8[z][y-1][x+1] == (unsigned char) 0) || (x+1 < imExt->dim.x && y+1 < imExt->dim.y && theMaskU8[z][y+1][x+1] == (unsigned char) 0)  )
                  continue;
              break;
          case USHORT:
          case SSHORT:
              if (theMaskU16[z][y][x] == (unsigned short int) 0 ||
                      (x+1 < imExt->dim.x && theMaskU16[z][y][x+1] == (unsigned short int) 0) ||
                      (x >= 1 && theMaskU16[z][y][x-1] == (unsigned short int) 0) ||
                      (y+1 < imExt->dim.y && theMaskU16[z][y+1][x] == (unsigned short int) 0) ||
                      (y >= 1 && theMaskU16[z][y-1][x] == (unsigned short int) 0) ||
                      (x >= 1 && y >= 1 && theMaskU16[z][y-1][x-1] == (unsigned short int) 0) || (x >= 1 && y+1 < imExt->dim.y && theMaskU16[z][y+1][x-1] == (unsigned short int) 0) ||
                      (x+1 < imExt->dim.x && y >= 1 && theMaskU16[z][y-1][x+1] == (unsigned short int) 0) || (x+1 < imExt->dim.x && y+1 < imExt->dim.y && theMaskU16[z][y+1][x+1] == (unsigned short int) 0)  )
                  continue;
              break;
          }
      }



    theta = (double)theTht[z][y][x];
    rep = theRep[z][y][x];
    theExt[z][y][x] = 0.0;

    dx = x + cos( theta );
    dy = y + sin( theta );
    r1 = _GetInterpolated2DValue( &(imsRes->imRep), dx, dy, z );

    if ( rep < r1 ) continue;

    dx = x - cos( theta );
    dy = y - sin( theta );
    r2 = _GetInterpolated2DValue( &(imsRes->imRep), dx, dy, z );

    if ( rep < r2 ) continue;

    theExt[z][y][x] = rep;
  }

}










/************************************************************
 *
 * /////////////// TENSOR VOTING ///////////////
 *
 ************************************************************/



int MT_SampleBin2D(vt_image *imageBin, double sample)
{
  char *proc = "MT_SampleBin2D";

  int nind;
  int *ind;
  size_t i, z;
  int n0;
  int tmp;
  int j;

  unsigned char* bufU8=NULL;
  float* bufFLT=NULL;

  if (sample < 0.0 || sample > 1.0)
  {
    fprintf(stderr, "%s : bad sample argument, expected to be between 0 and 1 (sample=%f)", proc, sample);
    return(-1);
  }

  for(z=0;z<imageBin->dim.z;z++)
  {
    nind=0;
    switch (imageBin->type)
    {
    case UCHAR :
      bufU8=(unsigned char*)imageBin->buf + z * imageBin->dim.x * imageBin->dim.y;
      for (i=0; i<(imageBin->dim.x*imageBin->dim.y); i++)
        if(bufU8[i] != '\0')
          nind ++;
      break;
    case FLOAT :
      bufFLT=(float*)imageBin->buf + z * imageBin->dim.x * imageBin->dim.y;
      for (i=0; i<(imageBin->dim.x*imageBin->dim.y); i++)
        if(bufFLT[i] > 0.0)
          nind ++;
      break;
    default:
      fprintf(stderr, "%s : image extension not implemented yet\n", proc);
      return(-1);
    }

    ind = vtmalloc( nind*sizeof(int), "ind", proc );
    j=0;

    switch (imageBin->type)
    {
    case UCHAR :
      for (i=0; i<(imageBin->dim.x*imageBin->dim.y); i++)
        if(bufU8[i] != '\0')
          ind[j++] = i;
      break;
    case FLOAT :
      for (i=0; i<(imageBin->dim.x*imageBin->dim.y); i++)
        if(bufFLT[i] > 0.0)
          ind[j++] = i;
      break;
    default:
      vtfree(ind);
      ind=(int*)NULL;
      fprintf(stderr, "%s : image extension not implemented yet\n", proc);
      return(-1);
    }

    n0 = (int) ( ((double) nind)*(1.0-sample) );

    /* fprintf(stdout, "nind=%d\tn0=%d\n", nind, n0);
     */

    srand(time(NULL));

    for (i=0; i<(size_t)n0; i++)
    {
      j=rand_(nind-i);
      tmp=ind[j];
       ind[j]=ind[nind-i-1];
      ind[nind-i-1]=tmp;
      switch (imageBin->type)
      {
      case UCHAR :
        bufU8[tmp]='\0';
        break;
      case FLOAT :
        bufFLT[tmp]=0.0;
        break;
      default:
        vtfree(ind);
        ind=(int*)NULL;
        fprintf(stderr, "%s : image extension not implemented yet\n", proc);
        return(-1);
      }
    }

    if(0 && _verbose_)
    {
      fprintf(stdout, "indices supprimes : [   ");
      for (i=0; i<(size_t)n0; i++)
        fprintf(stdout, "%d   ", ind[nind-n0+i]);
      fprintf(stdout, "]\n");
    }

    vtfree(ind);
    ind=(int*)NULL;
  }
  return (1);
}





static void _reinit2DTensorImage(mt_2Dtensor *theTensor)
{
  /* Reinitialise seulement les champs XX, YY et XY */
  int x,y,z;
  int dimx = theTensor->imxx.dim.x;
  int dimy = theTensor->imxx.dim.y;
  int dimz = theTensor->imxx.dim.z;

  float ***theXX = (float ***)theTensor->imxx.array;
  float ***theYY = (float ***)theTensor->imyy.array;
  float ***theXY = (float ***)theTensor->imxy.array;

  for (z=0; z<dimz; z++)
  for (y=0; y<dimy; y++)
  for (x=0; x<dimx; x++)
  {
    theXX[z][y][x] = 0;
    theYY[z][y][x] = 0;
    theXY[z][y][x] = 0;
  }
  return;

}





static int _add2DField(mt_2Dtensor *tensorImg, mt_2Dtensor field, double coef,
          int x, int y, int slice, enumVote typeVote)
{
  if (coef < FLTZERO)
    return(1);


  float ***fieldXX = (float***)(field.imxx.array);
  float ***fieldYY = (float***)(field.imyy.array);
  float ***fieldXY = (float***)(field.imxy.array);

  unsigned char ***fieldZeros = (unsigned char***)(field.iszero.array);


  float ***theXX = (float***)(tensorImg->imxx.array);
  float ***theYY = (float***)(tensorImg->imyy.array);
  float ***theXY = (float***)(tensorImg->imxy.array);

  unsigned char ***theZeros = (unsigned char***)(tensorImg->iszero.array);


  int xf0, yf0;
  int x0, y0, x1, y1;
  int i, j;
  int ifield, jfield;

  int sf[2]={(int)field.imxx.dim.x,(int)field.imxx.dim.y};

  int hsf[2];

  int xmax = tensorImg->imxx.dim.x;
  int ymax = tensorImg->imxx.dim.y;

  hsf[0]=sf[0]/2;
  hsf[1]=sf[1]/2;


  if (hsf[0]<=x)
  {
    xf0 = 0;
    x0 = x-hsf[0];
  }
  else
  {
    x0 = 0;
    xf0 = hsf[0]-x;
  }
  if (hsf[1]<=y)
  {
    yf0 = 0;
    y0 = y-hsf[1];
  }
  else
  {
    y0 = 0;
    yf0 = hsf[1]-y;
  }

  /* xf0 = (hsf[0]<=x) ? 0 : hsf[0]-x; */
  /* yf0 = (hsf[1]<=y) ? 0 : hsf[1]-y; */

  /* xf1 = (x+hsf[0]<xmax) ? sf[0] : (xmax-1-x)+hsf[0]+1; */
  /* yf1 = (y+hsf[1]<ymax) ? sf[1] : (ymax-1-y)+hsf[1]+1; */

  x1 = (x+hsf[0]<xmax) ? x+hsf[0]+1 : xmax;
  y1 = (y+hsf[1]<ymax) ? y+hsf[1]+1 : ymax;

  jfield = yf0;
  for (j=y0; j<y1; j++)
  {
    ifield = xf0;
    for (i=x0; i<x1;i++)
    {

      if (fieldZeros[0][jfield][ifield]==0)
        if (typeVote == DENSE || theZeros[slice][j][i]==0)
        {
            theZeros[slice][j][i]=0;
            theXX[slice][j][i]+=(float)coef*fieldXX[0][jfield][ifield];
            theYY[slice][j][i]+=(float)coef*fieldYY[0][jfield][ifield];
            theXY[slice][j][i]+=(float)coef*fieldXY[0][jfield][ifield];
        }
      ifield++;
    }
    jfield++;
  }

  return(1);

}





int MT_Compute2DTensorVoting(mt_2Dtensor *theTensor,
                             vt_image **imBin,
                             double scale, /* echelle en voxels */
                             int Niter,
                             int nangles,
                             enumTVmode TVmode,
                             hessianMode initHessian,
                             char *parName,
                             int writeImages)
{


  char *proc = "MT_Compute2DTensorVoting";
  int hsf;

  int dimFields[2];


  double *angles;
  int Nangles;

  double theCoeff[2];
  theCoeff[0] = theCoeff[1] = scale;
  hsf = ceil(sqrt(-pow(scale,2)*log(0.01)));
  dimFields[0] = dimFields[1] = 2*hsf+1;


  int x,y,slice;
  int dimx = imBin[0]->dim.x;
  int dimy = imBin[0]->dim.y;
  int nslice = imBin[0]->dim.z;



  float ***theXX;
  float ***theYY;
  float ***theXY;

  float ***theVP1;
  float ***theVP2;

  float ***theTheta2;

  unsigned char ***theZEROS;

  float ***imTht=NULL;

  float ***bin=NULL;
  unsigned char ***bin0;

  switch (imBin[0]->type)
  {
  case UCHAR :
    bin0 = (unsigned char***)imBin[0]->array;
    bin = vtmalloc( nslice*sizeof(float**), "bin", proc );
    for (slice=0;slice<nslice;slice++)
    {
      bin[slice] = vtmalloc( dimy*sizeof(float*), "bin[slice]", proc );
      for (y=0;y<dimy;y++)
      {
        bin[slice][y] = vtmalloc( dimx*sizeof(float), "bin[slice][y]", proc );
        for (x=0;x<dimx;x++)
          bin[slice][y][x] = bin0[slice][y][x] > 0 ? (float) 1 : (float) 0;
      }
    }
    break;
  case FLOAT :
    bin = (float***)imBin[0]->array;
    break;
  default:
    break;
  }




  if (initHessian != NONE)
  {
    imTht = (float***)imBin[1]->array;
  }

  int n,m,i;

  mt_2Dtensor *sfields;
  mt_2Dtensor bfield;

  double lambda1, lambda2;

  enumVote typeVote;
  char name[256];

  double theta;
  float val;

  /* Calcul des angles repartis de facon homogene sur le cercle unite */
  if (_verbose_)
    fprintf(stdout, "%s: computing angles...\n", proc);

  Nangles = _compute2DAngles(&angles, nangles);
  if (Nangles <= 0) {
    if (_verbose_)
      fprintf(stderr, "%s: error in computing angles\n", proc);
    return( -1 );
  }



  /* Calcul des champs de tenseurs Stick/Plate/Ball */
  if (_verbose_)
    fprintf(stdout, "%s: computing stick fields...\n", proc);

  if (MT_Compute2DStickFields( &sfields, dimFields,
          angles, theCoeff, Nangles, TVmode) != 1)
  {
    vtfree(angles);
    angles = NULL;
    switch (imBin[0]->type)
    {
    case UCHAR :
      for(slice=0;slice<nslice;slice++)
      {
        for(y=0;y<dimy;y++)
        {
          vtfree(bin[slice][y]);
          bin[slice][y]=NULL;
        }
        vtfree(bin[slice]);
        bin[slice]=NULL;
      }
      vtfree(bin);
      bin=NULL;
      break;
    default:
      break;
    }
    return (-1);
  }

  if (initHessian == NONE)
  {
    if (_verbose_)
      fprintf(stdout, "%s: computing ball fields...\n", proc);

    if (MT_Compute2DBallFieldFromStick( &bfield, sfields, dimFields,
            Nangles ) != 1)
    {
      for (n=0;n<Nangles;n++){
        VT_Free2DTensor(sfields+n);
      }
      switch (imBin[0]->type)
      {
      case UCHAR :
        for(slice=0;slice<nslice;slice++)
        {
          for(y=0;y<dimy;y++)
          {
            vtfree(bin[slice][y]);
            bin[slice][y]=NULL;
          }
          vtfree(bin[slice]);
          bin[slice]=NULL;
        }
        vtfree(bin);
        bin=NULL;
        break;
      default:
        break;
      }
      vtfree(sfields);
      sfields=NULL;
      vtfree(angles);
      angles = NULL;
      return (-1);
    }
  }

  /*  Tensor image allocation */
  /*
  if (_verbose_)
    fprintf(stdout, "%s: allocating tensor image...\n", proc);

  sprintf( name, "%s", parName );
  if (VT_Alloc2DTensorFromImage( theTensor, (imBin[0]), name ) !=1)
  {
    if (_verbose_)
      fprintf(stderr, "%s: allocation of tensor image failed\n", proc);

    for (n=0;n<Nangles;n++){
      VT_Free2DTensor(sfields+n);
    }
    VT_Free2DTensor(&bfield);
    vtfree(sfields);
    sfields=NULL;
    vtfree(angles);
    angles = NULL;
    switch (imBin[0]->type)
    {
    case UCHAR :
      for(slice=0;slice<nslice;slice++)
      {
        for(y=0;y<dimy;y++)
        {
          vtfree(bin[slice][y]);
          bin[slice][y]=NULL;
        }
        vtfree(bin[slice]);
        bin[slice]=NULL;
      }
      vtfree(bin);
      bin=NULL;
      break;
    default:
      break;
    }

    return (-1);

  }
  */

  theXX = (float ***)theTensor->imxx.array;
  theYY = (float ***)theTensor->imyy.array;
  theXY = (float ***)theTensor->imxy.array;

  theVP1 = (float ***)theTensor->imvp1.array;
  theVP2 = (float ***)theTensor->imvp2.array;

  theTheta2 = (float ***)theTensor->imtheta2.array;

  theZEROS = (unsigned char ***)theTensor->iszero.array;


  /* Initialisation de l'image tenseur a partir d'imBin */
  if (_verbose_)
    fprintf(stdout,"%s: initialisation du tenseur a partir d'imBin...\n", proc);

  n=0;
  m=0;

  for (slice=0;slice<nslice;slice++)
  {
    if (_verbose_ && (nslice<=30 || slice%(nslice/30)==0))
      fprintf(stdout,"\tFrame #%d/%d (%d%%)... \n",
          slice+1, slice, (int)(100*(float)slice/(float)nslice));
    for (y=0;y<dimy;y++)
    for (x=0;x<dimx;x++)
    {
      val = bin[slice][y][x];
      if ( val < FLTZERO)
      {
        /*  le voxel(x,y,z) ne vote pas ni ne recoit de votes  */
        theXX[slice][y][x] = 0;
        theYY[slice][y][x] = 0;
        theXY[slice][y][x] = 0;
        theZEROS[slice][y][x] = (unsigned char)1;
        n++;
        continue;
      }
      if (initHessian == NONE)
      {
      /* Tenseur Ball */
      theXX[slice][y][x] = val;
      theYY[slice][y][x] = val;
      theXY[slice][y][x] = 0;

      theVP1[slice][y][x] = val;
      theVP2[slice][y][x] = val;
      }
      else
      {
        /* Tenseur Stick */

        theTheta2[slice][y][x] = imTht[slice][y][x];

        theVP1[slice][y][x] = 0;
        theVP2[slice][y][x] = val;
      }
      theZEROS[slice][y][x] = (unsigned char)0;
      m++;
    }
  }

  if (initHessian != NONE)
    _init2DTensorFromAngle(theTensor);


  if (_verbose_)
  {
    fprintf(stdout, "%s : nombre de zeros = %d \tnombre de non-zeros = %d \n",
        proc, n,m);
  }

  /* Niter etapes de votes SPARSE + 1 etape de votes DENSE */

  typeVote = SPARSE;

  if (initHessian != NONE)
    Niter=0;

  for (i=0;i<=Niter;i++)
  {

    if (writeImages && parName != (char*)NULL && parName[0] != '\0' )
      {
      if (_verbose_)
      {
        sprintf( name, "%s.eparse.%d", parName,i );
        fprintf(stdout,
        "Ecriture du champ de tenseurs eparse #%d dans %s.inr...\n", i, name);
      }
      VT_Write2DTensorWithName( theTensor, name );
    }
    if (_verbose_)
      fprintf(stdout, "%s: etape de vote : iteration #%d/%d...\n",
        proc, i+1, Niter+1);

    /* Vote DENSE si derniere iteration */
    if (i == Niter)
    {
        typeVote = DENSE;
    }


    if (_verbose_) {
      if ( typeVote == SPARSE) fprintf(stdout, "\tSPARSE voting...\n");
      else fprintf(stdout, "\tDENSE voting...\n");
        }

    _reinit2DTensorImage(theTensor);

    /* Cumul des votes (seuls imxx,yy,zz,xy,xz,yz et iszero sont affectes) */
    for (slice=0;slice<nslice;slice++)
    {
      if (_verbose_)
        if (nslice<=30 || slice%(nslice/30)==0)
          fprintf(stdout, "\t\tTreating frame #%d/%d (%d%%)...\n",
              slice+1, nslice, (int)(100*(float)slice/(float)nslice));

      for (y=0;y<dimy;y++)
      for (x=0;x<dimx;x++)
      {
        if ( theZEROS[slice][y][x] != 0)
        {
          /*  le voxel(x,y,z) ne vote pas  */
          continue;
        }

        lambda1 = theVP1[slice][y][x];
        lambda2 = theVP2[slice][y][x];

        if (lambda2-lambda1 > FLTZERO)  /* STICKS */
        {
          theta = theTheta2[slice][y][x];
          n = _nearest2DAngle(theta, angles, Nangles);
          _add2DField(theTensor, sfields[n], lambda2-lambda1,x,y,slice, typeVote);
        }

        if (lambda1 > FLTZERO)          /* BALLS */
          _add2DField(theTensor, bfield, lambda1,x,y,slice, typeVote);
      }
    }

    /* Conversion de la hessienne en vp, angles, etc... */
    if (_verbose_)
      fprintf(stdout, "\tConversion des tenseurs en valeurs/vecteurs propres, \
          etc...\n");

    for (slice=0;slice<nslice;slice++)
    for (y=0;y<dimy;y++)
    for (x=0;x<dimx;x++)
    {
      if ( theZEROS[slice][y][x] == 1)
      {
        /*  le voxel(x,y,z) contient un tenseur nul */
        continue;
      }

      if (_set2DTensor(theTensor, x,y,slice, typeVote) !=1)
      {
        if (_verbose_)
          fprintf(stderr, "%s: probleme lors du calcul des parametres des \
              tenseurs", proc);
        for (n=0;n<Nangles;n++){
          VT_Free2DTensor(sfields+n);
        }
        vtfree(sfields);
        sfields=NULL;
        VT_Free2DTensor(&bfield);
        vtfree(angles);
        angles = NULL;
        VT_Free2DTensor(theTensor);
        switch (imBin[0]->type)
        {
        case UCHAR :
          for(slice=0;slice<nslice;slice++)
          {
            for(y=0;y<dimy;y++)
            {
              vtfree(bin[slice][y]);
              bin[slice][y]=NULL;
            }
            vtfree(bin[slice]);
            bin[slice]=NULL;
          }
          vtfree(bin);
          bin=NULL;
          break;
        default:
          break;
        }
        return(-1);
      }

    }

  } /* Niter */


  if (_verbose_)
    fprintf(stdout, "%s: Liberation de la memoire allouee...\n", proc);

  for (n=0;n<Nangles;n++){
    VT_Free2DTensor(sfields+n);
  }
  vtfree(sfields);
  sfields=NULL;
  if (initHessian == NONE)
    VT_Free2DTensor(&bfield);
  vtfree(angles);
  angles = NULL;

  switch (imBin[0]->type)
  {
  case UCHAR :
    for(slice=0;slice<nslice;slice++)
    {
      for(y=0;y<dimy;y++)
      {
        vtfree(bin[slice][y]);
        bin[slice][y]=NULL;
      }
      vtfree(bin[slice]);
      bin[slice]=NULL;
    }
    vtfree(bin);
    bin=NULL;
    break;
  default:
    break;
  }


  return(1);
}









/************************************************************
 *
 *
 *
 *
 ************************************************************/



void MT_Compute2DTensorLineicExtrema( mt_2Dtensor *imTensor,
    vt_image *imExt )
{
  size_t x, y, z;

  float ***theExt = (float***)(imExt->array);

  float ***theVP1 = (float***)(imTensor->imvp1.array);
  float ***theVP2 = (float***)(imTensor->imvp2.array);
  float ***theTht2 = (float***)(imTensor->imtheta2.array);
  unsigned char ***theZeros = (unsigned char***)(imTensor->iszero.array);



  double v2[2];

  double rep, rn[2];

  double theta2;

  double vx=0.0, vy=0.0;

  
  double rx, ry;
  double dx, dy;
  int j, ix, iy;
  double dxdy;
  double v4,v5;

  /* just to avoid compilation warnings
   */
  rn[0] = rn[1] = 0.0;

  for ( z=0; z<imExt->dim.z; z++ )
  for ( y=0; y<imExt->dim.y; y++ )
  for ( x=0; x<imExt->dim.x; x++ ) {

    theta2 = (double)theTht2[z][y][x];
    
    theExt[z][y][x] = 0.0;
    if ( theZeros[z][y][x] != 0 ) continue;

    rep = (double)(theVP2[z][y][x]-theVP1[z][y][x]);



    /* v2 est le vecteur normal a la structure lineique
       v1 est tangent a la structure
    */
    /* v1[0]=cos(theta1); v1[1]=sin(theta1); */
    v2[0]=cos(theta2); v2[1]=sin(theta2);


    for (j=0;j<2;j++) {
      switch(j) {
      default:
      case 0 :
  vx = v2[0];    vy = v2[1]; break;
      case 1 :
  vx = -v2[0];  vy = -v2[1]; break;
      }

      rx = x + vx;
      ix = (int)rx;
      if ( rx <= 0.0 || ix >= (int)imExt->dim.x-1 ) {
        j = 3;
        continue;
      }
      ry = y + vy;
      iy = (int)ry;
      if ( ry <= 0.0 || iy >= (int)imExt->dim.y-1 ) {
        j = 3;
        continue;
      }

      dx = rx - ix;
      dy = ry - iy;

      dxdy = dx*dy;

      v4 = 1-dx;
      v5 = 1-dy;

      rn[j] = 0;
      rn[j]+= dxdy       *(theVP2[z][iy+1][ix+1]-theVP1[z][iy+1][ix+1]);
      rn[j]+= dx*v5      *(theVP2[z][iy+1][ix  ]-theVP1[z][iy+1][ix  ]);
      rn[j]+= v4*dy      *(theVP2[z][iy  ][ix+1]-theVP1[z][iy  ][ix+1]);
      rn[j]+= v4*v5      *(theVP2[z][iy  ][ix  ]-theVP1[z][iy  ][ix  ]);

      if ( rn[j] > rep ) {
        j = 3;
        continue;
      }

    } /* fin de la boucle suivant la normale a la membrane */

    if ( rn[0] == rep && rn[1] == rep ) {
      continue;
    }

    if ( j < 3 ) theExt[z][y][x] = rep;
  } /* x,y,z */


}





void MT_Compute2DTensorBallExtrema( mt_2Dtensor *imTensor,
    vt_image *imExt )
{
  size_t x, y, z;

  float ***theExt = (float***)(imExt->array);
  float ***theVP1 = (float***)(imTensor->imvp1.array);
  
  float rep;
  
  int i;
  
  int dx[]={-1,  0,  1,
			-1,      1,
			-1,  0,  1};
  
  int dy[]={-1, -1, -1,
			 0,      0,
			 1,  1,  1};
  
  for ( z=0; z<imExt->dim.z  ; z++ ) 
  for ( y=1; y<imExt->dim.y-1; y++ ) 
  for ( x=1; x<imExt->dim.x-1; x++ ) {
	rep=theVP1[z][y][x];
	
	if (rep<=0.0) continue;
	
	
	for (i=0; i<8; i++) {
	  if ( rep > theVP1[z][y+dy[i]][x+dx[i]] )
		continue;
	  i=9;
	}
	if ( i == 26 )
	  theExt[z][y][x] = rep;
  }
}





static int _compute2DStickField( mt_2Dtensor *sfield, int *dimFields,
          double angle, double *theCoeffs, enumTVmode mode)
{
  int x,y;

  int mid[2] = {dimFields[0]/2, dimFields[1]/2};
  double r[2];
  double n[2];
  double a[2];
  double rr[3]; /* rr = r*transpose(r) matrice symetrique 2x2 = 3 elts */
  double v[2];
  double sinteta;
  double fact=0;
  double alpha;
  double s, kappa, l;
  double SUM = 0;

  double scale = theCoeffs[0];
  double pi = 3.14159265358979323846;
  double c = -16*log(0.1)*(scale-1)/pow(pi,2);

  float ***theXX = (float***)(sfield->imxx.array);


  float ***theYY = (float***)(sfield->imyy.array);
  float ***theXY = (float***)(sfield->imxy.array);

  float ***theVP1 = (float***)(sfield->imvp1.array);
  float ***theVP2 = (float***)(sfield->imvp2.array);


  float ***theTht1 = (float***)(sfield->imtheta1.array);
  float ***theTht2 = (float***)(sfield->imtheta2.array);

  unsigned char ***theZeros = (unsigned char***)(sfield->iszero.array);

  n[0]=cos(angle); n[1]=sin(angle);

  for (y=0;y<dimFields[1];y++)
  for (x=0;x<dimFields[0];x++) {

    r[0]=(double) x-mid[0];
    r[1]=(double) y-mid[1];
    l = sqrt(pow(r[0],2) + pow(r[1],2));

    if (l == 0)
    {
      /*  Ici le point qui recoit le vote est le point votant */
      theXX[0][y][x] = n[0]*n[0];
      theYY[0][y][x] = n[1]*n[1];
      theXY[0][y][x] = n[0]*n[1];
      theVP1[0][y][x]= 0;
      theVP2[0][y][x]= 1;
      theTht1[0][y][x]= angle+pi/(double)2;
      theTht2[0][y][x]= angle;
      theZeros[0][y][x]= 0;
      SUM += 1;
      continue;
    }

    r[0]=r[0]/l;
    r[1]=r[1]/l;


    if (mode == TVCLASSIC && fabs(r[0]*n[0]+r[1]*n[1])>sqrt(2)/2)
    {
      /* angle inferieur a pi/4 entre la normale en le point votant et
         le vecteur reliant les deux points -> on ne vote pas */
      theXX[0][y][x] = 0;
      theYY[0][y][x] = 0;
      theXY[0][y][x] = 0;
      theVP1[0][y][x]= 0;
      theVP2[0][y][x]= 0;
      theTht1[0][y][x]= angle+pi/(double)2;
      theTht2[0][y][x]= angle;
      theZeros[0][y][x]= 1;
      continue;
    }

    /* Calcul du vecteur v normal a l'arc de cercle reliant le point
       votant et le point recevant le vote */
    rr[0]=r[0]*r[0]; /* rxx */
    rr[1]=r[1]*r[1]; /* ryy */
    rr[2]=r[0]*r[1]; /* rxy */

    a[0] = rr[0]*n[0]+rr[2]*n[1];
    a[1] = rr[2]*n[0]+rr[1]*n[1];

    v[0] = n[0]-2*a[0];
    v[1] = n[1]-2*a[1];

    /* Calcul du facteur de decroissance en le pixel [x,y] */
    sinteta = r[0]*n[0]+r[1]*n[1];

    switch (mode) {
    case TVCLASSIC :
      s = (sinteta == 0) ? l : l*asin(sinteta)/sinteta;
      kappa = 2*sinteta/l;
      fact = exp(-(pow(s,2)+c*pow(kappa,2))/pow(scale,2));
      break;
    case CFTV :
      fact = exp(-pow(l/scale,2))*(1-pow(sinteta,2));
      break;
    }

    if (fact < 0.01) fact = 0;

    if (fabs(v[0])<FLTZERO)
      alpha = asin(v[1]);
    else
      alpha = atan(v[1]/v[0]);

    theXX[0][y][x] = v[0]*v[0]*fact;
    theYY[0][y][x] = v[1]*v[1]*fact;
    theXY[0][y][x] = v[0]*v[1]*fact;

    theVP1[0][y][x]= 0;
    theVP2[0][y][x]= fact;

    theTht1[0][y][x]= alpha+pi/(double)2;
    theTht2[0][y][x]= alpha;

    theZeros[0][y][x]= (fact>0) ? 0 : 1;

    SUM += fact;
  } /* for(x,y,z) */
  for (y=0;y<dimFields[1];y++)
  for (x=0;x<dimFields[0];x++)
  {
    theVP2[0][y][x] = theVP2[0][y][x]/(float)SUM;
    theXX[0][y][x] = theXX[0][y][x]/(float)SUM;
    theYY[0][y][x] = theYY[0][y][x]/(float)SUM;
    theXY[0][y][x] = theXY[0][y][x]/(float)SUM;
  }

  return( 1 );
}





int MT_Compute2DStickFields( mt_2Dtensor **sfields, int *dimFields,
          double *angles, double *theCoeffs, int Nangles, enumTVmode mode)
{
  int n,m;
  char *proc = "MT_Compute2DStickFields";
  double angle;
  char name[256];
  int t = (int) FLOAT;
  mt_2Dtensor *field;

  (*sfields) =  vtmalloc( Nangles*sizeof(mt_2Dtensor), "(*sfields)", proc );

  if((*sfields)==NULL)
  {
    if (_verbose_)
      fprintf(stderr, "%s: allocation de (*sfield) echouee", proc);
    return(-1);
  }

  for (n=0; n<Nangles; n++) {
    sprintf( name, "%s.%d", "stickfield", n );
    field = (*sfields+n);
    if ( VT_Alloc2DTensor(field, name, dimFields[0], dimFields[1], 1, t ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate 2D tensor images\n", proc );
      for (m=0;m<n;m++){
        field = (*sfields+m);
        VT_Free2DTensor( field );
      }
      vtfree((*sfields));
      sfields[0]=NULL;
      return( -1 );

    }
  }

  for (n=0; n<Nangles; n++) {
    angle=angles[n];
    field = (*sfields+n);
    if (_compute2DStickField(field, dimFields, angle, theCoeffs, mode) != 1)
    {
      if (_verbose_)
        fprintf(stderr, "%s: error in computing stick field #%d\n", proc, n);
      for (m=0;m<Nangles;m++)
      {
        field = (*sfields+m);
        VT_Free2DTensor( field );
      }
      vtfree((*sfields));
      (*sfields)=NULL;
      return( -1 );
    }
  }

  return( 1 );

}





int MT_Compute2DBallFieldFromStick( mt_2Dtensor *bfield,
          mt_2Dtensor *sfields, int *dimFields,
          int Nangles)
{
  char *proc = "MT_Compute2DBallFieldFromStick";
  int n, x, y;
  double hessien[4];
  double valprop[2];
  double vecprop[4];
  double v[2];
  double SUM = 0;
  double theta;

  float ***sXX;
  float ***sYY;
  float ***sXY;

  float ***theXX;
  float ***theYY;
  float ***theXY;

  float ***theVP1;
  float ***theVP2;

  float ***theTht1;
  float ***theTht2;

  unsigned char ***theZeros;

  char name[256];
  int t = (int) FLOAT;

  sprintf( name, "%s", "ballfield" );
  if ( VT_Alloc2DTensor( bfield, name, dimFields[0], dimFields[1], 1, t ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate 2D tensor images\n", proc );
    return( -1 );

  }

  theXX = (float***)(bfield->imxx.array);
  theYY = (float***)(bfield->imyy.array);
  theXY = (float***)(bfield->imxy.array);

  theVP1 = (float***)(bfield->imvp1.array);
  theVP2 = (float***)(bfield->imvp2.array);

  theTht1 = (float***)(bfield->imtheta1.array);
  theTht2 = (float***)(bfield->imtheta2.array);

  theZeros = (unsigned char***)(bfield->iszero.array);

  for (y=0;y<dimFields[1];y++)
  for (x=0;x<dimFields[0];x++) {
    theXX[0][y][x] = 0;
    theYY[0][y][x] = 0;
    theXY[0][y][x] = 0;
    theZeros[0][y][x]=1;
  }

  for (n=0;n<Nangles;n++)
  {
    mt_2Dtensor* sfield = sfields+n;
    sXX = (float ***)sfield->imxx.array;
    sYY = (float ***)sfield->imyy.array;
    sXY = (float ***)sfield->imxy.array;

    for (y=0;y<dimFields[1];y++)
    for (x=0;x<dimFields[0];x++) {

      theXX[0][y][x] = theXX[0][y][x] + sXX[0][y][x];
      theYY[0][y][x] = theYY[0][y][x] + sYY[0][y][x];
      theXY[0][y][x] = theXY[0][y][x] + sXY[0][y][x];
    }
  }

  for (y=0;y<dimFields[1];y++)
  for (x=0;x<dimFields[0];x++) {
    hessien[0] = theXX[0][y][x];
    hessien[1] = hessien[2] = theXY[0][y][x];
    hessien[3] = theYY[0][y][x];

    if ( _ComputeEigensOfSymetricSquareMatrix( hessien, valprop, vecprop, 2 )
        != 1 ) {
      if (_verbose_)
        fprintf(stderr, "%s: error in computing eigens\n", proc);
      VT_Free2DTensor(bfield);
      return( -1 );
    }
    if ( _SortEigensInAbsIncreasingOrder( valprop, vecprop, 2 ) != 1 ) {

      if (_verbose_)
        fprintf(stderr, "%s: error in sorting eigens\n", proc);
      VT_Free2DTensor(bfield);
      return( -1 );
    }

    /* Les valeurs et vecteurs propres sont tries par fabs croissante */
    theVP1[0][y][x]=valprop[0];
    v[0]=vecprop[0]; v[1]=vecprop[2];
    if (fabs(v[0])<FLTZERO)
      theta = asin(v[1]);
    else
      theta = atan(v[1]/v[0]);
    theTht1[0][y][x]=theta;

    theVP2[0][y][x]=valprop[1];
    v[0]=vecprop[1]; v[1]=vecprop[3];
    if (fabs(v[0])<FLTZERO)
      theta = asin(v[1]);
    else
      theta = atan(v[1]/v[0]);
    theTht2[0][y][x]=theta;

    if (valprop[1]>0)
    {
      theZeros[0][y][x] = 0;
      SUM += valprop[1];
    }
  }

  for (y=0;y<dimFields[1];y++)
  for (x=0;x<dimFields[0];x++)
  {
    theVP1[0][y][x] = theVP1[0][y][x]/(float)SUM;
    theVP2[0][y][x] = theVP2[0][y][x]/(float)SUM;
    theXX[0][y][x] = theXX[0][y][x]/(float)SUM;
    theYY[0][y][x] = theYY[0][y][x]/(float)SUM;
    theXY[0][y][x] = theXY[0][y][x]/(float)SUM;

  }
  return( 1 );

}


