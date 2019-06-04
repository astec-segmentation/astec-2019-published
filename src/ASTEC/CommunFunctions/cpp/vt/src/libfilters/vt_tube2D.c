/*************************************************************************
 * vt_tube2D.c -
 *
 * $Id: vt_tube2D.c,v 1.8 2006/05/16 09:33:34 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 *
 * CREATION DATE:
 * Mon Jul 10 17:03:57 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <vt_tube2D.h>

#include <linearFiltering.h>

static int _borderLength_[3] = { 10, 10, 10 };
static int _verbose_ = 1;













typedef struct {
  int x;
  int y;
  double c;
} type2DConvolutionPoint;



/*
  largeSigma : dans la direction du vecteur V = ( cos(T), sin(T) )
  smallSigma : dans la direction orthogonale du vecteur

  autour d'un point M, on garde les points tels que
  ((mM . V)/largeSigma)^2 + (|mM - (mM.V) V|/smallSigma)^2 <= multSigma^2

*/

int VT_2DTensorVoting( vt_2Dtensor *par,
               vt_image *imWeight,
               vt_image *imTheta,
               vt_image *imMask,
               double largeSigma,
               double smallSigma,
               double multSigma )
{
  char *proc = "VT_2DTensorVoting";
  int x, y, z;
  int i, j;
  int n, nb;
  int r;
  double sum;

  s8 ***theMask = (s8 ***)NULL;
  float ***theWeight = (float ***)imWeight->array;
  float ***theTHETA = (float ***)imTheta->array;

  float ***theXX = (float ***)par->imxx.array;
  float ***theXY = (float ***)par->imxy.array;
  float ***theYY = (float ***)par->imyy.array;

  type2DConvolutionPoint *convolutionArray = (type2DConvolutionPoint *)NULL;

  double c2 = multSigma * multSigma;
  double ls2 = largeSigma * largeSigma;
  double ss2 = smallSigma * smallSigma;

  double txx, txy, tyy, vx, vy;
  double largeC, smallC;
  double c, e;


  if ( largeSigma - largeSigma/1000 < smallSigma
       && smallSigma < largeSigma + largeSigma/1000 )
    return( VT_2DTensorGaussianVoting( par, imWeight, imTheta, imMask, smallSigma ) );

  if ( largeSigma > smallSigma ) r = (int)(multSigma*largeSigma) + 1;
  else                   r = (int)(multSigma*smallSigma) + 1;

  if ( r <= 1 ) {
    fprintf( stderr, "%s: too small window (r=%d, sigmas=[%f,%f], coeff=%f)\n",
         proc, r, largeSigma, smallSigma, multSigma );
    return( -1 );
  }



  convolutionArray = (type2DConvolutionPoint *)VT_Malloc( (2*r+1)*(2*r+1) *
                              sizeof( type2DConvolutionPoint ) );


  if ( imMask != (vt_image*)NULL ) {
    if ( imMask->type == SCHAR ) theMask = (s8 ***)imMask->array;
  }


  for ( z=0; z<(int)par->imxx.dim.z; z++ )
  for ( y=0; y<(int)par->imxx.dim.y; y++ )
  for ( x=0; x<(int)par->imxx.dim.x; x++ ) {
    theXX[z][y][x] = theXY[z][y][x] = theYY[z][y][x] = 0.0;
  }


  for ( z=0; z<(int)par->imxx.dim.z; z++ )
  for ( y=0; y<(int)par->imxx.dim.y; y++ ) {

    if ( _verbose_ )
      fprintf( stderr, "\r process slice %4d/%lu   row %4d/%lu",
           z+1, par->imxx.dim.z, y+1, par->imxx.dim.y );

    for ( x=0; x<(int)par->imxx.dim.x; x++ ) {

      if ( theMask != (s8 ***)NULL ) {
    if ( theMask[z][y][x] < 0 ) continue;
      }
      if ( theWeight[z][y][x] <= 0.0 ) continue;

      vx = cos( theTHETA[z][y][x] );
      vy = sin( theTHETA[z][y][x] );

      nb = 0;
      sum = 0;
      for ( j = -r; j <= r; j ++ )
      for ( i = -r; i <= r; i ++ ) {
    /* M = (x,y) m = (x+i,y+j) Mm = (i,j)
       largeC = mM.V
       smallC = |mM - (mM.V) V|^2
       mM - (mM.V) V = ( x-i ) - largeC * ( vx )
                       ( y-j )            ( vy )
       mM - (mM.V) V = ( x-i - largeC * vx )
                       ( y-j - largeC * vy )
    */
    largeC = i*vx + j*vy;
    smallC = (i - largeC*vx)*(i - largeC*vx) +
             (j - largeC*vy)*(j - largeC*vy);
    e = largeC*largeC/ls2 + smallC/ss2;
    if ( e > c2 ) continue;

    /* on a un point du masque
       on le compte dans la somme;
    */
    c = exp( -e / 2.0 );
    sum += c;

    if ( x+i < 0 || x+i >= (int)par->imxx.dim.x ||
         y+j < 0 || y+j >= (int)par->imxx.dim.y ) continue;

    convolutionArray[ nb ].x = x+i;
    convolutionArray[ nb ].y = y+j;
    convolutionArray[ nb ].c = c;
    nb ++;
      }

      /* on a un masque de convolution dont les coefficients
     ont pour somme 'sum'
     dont les indices varient de 0 a 'nb'
      */
      txx = theWeight[z][y][x] * vx * vx;
      txy = theWeight[z][y][x] * vx * vy;
      tyy = theWeight[z][y][x] * vy * vy;

      /* on applique la convolution
       */
      for ( n=0; n<nb; n++ ) {
    i = convolutionArray[ n ].x;
    j = convolutionArray[ n ].y;
    c = convolutionArray[ n ].c / sum;
    theXX[z][j][i] += c * txx;
    theXY[z][j][i] += c * txy;
    theYY[z][j][i] += c * tyy;
      }

    }
  }

  if ( _verbose_ )
    fprintf( stderr, " ... done\n" );


  VT_Free( (void**)&convolutionArray );
  return( 1 );
}














int VT_2DTensorGaussianVoting( vt_2Dtensor *par,
                   vt_image *imWeight,
                   vt_image *imTheta,
                   vt_image *imMask,
                   double sigma )
{
  char *proc = "VT_2DTensorGaussianVoting";
  int x, y, z;


  s8 ***theMask = (s8 ***)NULL;
  float ***theWeight = (float ***)imWeight->array;
  float ***theTHETA = (float ***)imTheta->array;

  float ***theXX = (float ***)par->imxx.array;
  float ***theXY = (float ***)par->imxy.array;
  float ***theYY = (float ***)par->imyy.array;

  typeFilteringCoefficients filter[3];

  int theDim[3];


  initFilteringCoefficients( &(filter[0]) );
  initFilteringCoefficients( &(filter[1]) );
  initFilteringCoefficients( &(filter[2]) );
  filter[0].type = GABOR_YOUNG_2002;
  filter[1].type = GABOR_YOUNG_2002;
  filter[2].type = GABOR_YOUNG_2002;
  filter[0].coefficient = sigma;
  filter[1].coefficient = sigma;
  filter[2].coefficient = sigma;
  filter[0].derivative = DERIVATIVE_0;
  filter[1].derivative = DERIVATIVE_0;
  filter[2].derivative = NODERIVATIVE;



  if ( imMask != (vt_image*)NULL ) {
    if ( imMask->type == SCHAR ) theMask = (s8 ***)imMask->array;
  }

  for ( z=0; z<(int)par->imxx.dim.z; z++ )
  for ( y=0; y<(int)par->imxx.dim.y; y++ )
  for ( x=0; x<(int)par->imxx.dim.x; x++ ) {
    if ( theMask != (s8 ***)NULL ) {
      if ( theMask[z][y][x] < 0 ) {
    theXX[z][y][x] = theXY[z][y][x] = theYY[z][y][x] = 0.0;
    continue;
      }
    }
    if ( theWeight[z][y][x] <= 0.0 ) {
      theXX[z][y][x] = theXY[z][y][x] = theYY[z][y][x] = 0.0;
    } else {
      theXX[z][y][x] = theWeight[z][y][x]*cos( theTHETA[z][y][x] )*cos( theTHETA[z][y][x] );
      theXY[z][y][x] = theWeight[z][y][x]*cos( theTHETA[z][y][x] )*sin( theTHETA[z][y][x] );
      theYY[z][y][x] = theWeight[z][y][x]*sin( theTHETA[z][y][x] )*sin( theTHETA[z][y][x] );
    }
  }


  if ( sigma <= 0.0 ) return ( 1 );


  theDim[0] = par->imxx.dim.x;
  theDim[1] = par->imxx.dim.y;
  theDim[2] = par->imxx.dim.z;

  if ( separableLinearFiltering( par->imxx.buf, par->imxx.type,
				 par->imxx.buf, par->imxx.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computating derivatives (XX)\n", proc );
    return( -1 );
  }
  
  if ( separableLinearFiltering( par->imxy.buf, par->imxy.type,
				 par->imxy.buf, par->imxy.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computating derivatives (XY)\n", proc );
    return( -1 );
  }

  if ( separableLinearFiltering( par->imyy.buf, par->imyy.type,
				 par->imyy.buf, par->imyy.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computating derivatives (YY)\n", proc );
    return( -1 );
  }

  return( 1 );
}












int VT_Compute2DMultiScale( vt_image *theIm,
                vt_3Dimres *imsRes,
                double scale1,
                double scale2,
                int nbscales,
                enumStructureColor color,
                methodType mode )
{
  char *proc = "VT_Compute2DMultiScale";
  int s, nbs = nbscales;
  double scale;

  vt_2Dimages theIms;

  float theCoeff[3];

  int x, y, z;
  float ***maxRep = (float***)imsRes->imRep.array;
  float ***maxTht = (float***)imsRes->imTheta.array;
  float ***maxScl = (float***)imsRes->imScale.array;

  float ***theRep;
  float ***theTht;






  if ( VT_Alloc2Dimages( &theIms, theIm, proc ) != 1 ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to allocate 2D auxiliary images\n", proc );
      return( -1 );
    }
  }

  theRep = (float***)theIms.imr.array;
  theTht = (float***)theIms.hessien.imtheta1.array;



  if ( scale2 <= scale1 ) nbs = 1;

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

    if (mode == KRISSIAN) {
      if ( VT_Filter2Dimages( theIm, &theIms, theCoeff ) != 1 ) {
        VT_Free2Dimages( &theIms );
        if ( _verbose_ ) {
          fprintf( stderr, "%s: unable to 2D filter\n", proc );
          return( -1 );
        }
      }
    }
    else {
      if ( VT_Filter2Dimages2ndOrder( theIm, &theIms, theCoeff ) != 1 ) {
        VT_Free2Dimages( &theIms );
        if ( _verbose_ ) {
          fprintf( stderr, "%s: unable to 2D filter\n", proc );
          return( -1 );
        }
      }
    }

    VT_Compute2DEigenVectors( &(theIms.hessien), &(theIms.ime) );
    VT_FilterOn2DEigenValues( &theIms, color );
    switch (mode) {
    default:
    case KRISSIAN:
      VT_Compute2DResponse( &theIms, color, 1.0, scale );
      break;
    case FRANGI:
      VT_Compute2DResponseFrangi( &theIms, scale );
      break;
    }
    /* on recopie le resultat dans les images resultats
     */

    if ( s == 0 ) {
      for ( z = 0; z < (int)theIm->dim.z; z ++ )
      for ( y = 0; y < (int)theIm->dim.y; y ++ )
      for ( x = 0; x < (int)theIm->dim.x; x ++ ) {
    if ( theRep[z][y][x] > 0.0 ) {
      maxRep[z][y][x] = (mode==KRISSIAN ? scale : 1.0) * theRep[z][y][x];
      maxTht[z][y][x] = theTht[z][y][x];
      maxScl[z][y][x] = scale * theIm->siz.x;
    } else {
      maxRep[z][y][x] = maxTht[z][y][x] = 0.0;
      maxScl[z][y][x] = 0.0;
    }
      }
    } else {
      for ( z = 0; z < (int)theIm->dim.z; z ++ )
      for ( y = 0; y < (int)theIm->dim.y; y ++ )
      for ( x = 0; x < (int)theIm->dim.x; x ++ ) {
    if ( (mode==KRISSIAN ? scale : 1.0) * theRep[z][y][x] > maxRep[z][y][x] ) {
      maxRep[z][y][x] = (mode==KRISSIAN ? scale : 1.0) * theRep[z][y][x];
      maxTht[z][y][x] = theTht[z][y][x];
      maxScl[z][y][x] = scale * theIm->siz.x;
    }
      }
    }

  } /* for ( s = 0; s < nbs; s ++ ) */

  VT_Free2Dimages( &theIms );
  return( 1 );
}







void VT_Compute2DExtrema( vt_3Dimres *imsRes,
              vt_image *imExt )
{
  float ***theExt = (float***)(imExt->array);
  float ***theRep = (float***)(imsRes->imRep.array);
  float ***theTht = (float***)(imsRes->imTheta.array);

  int x, y, z;
  double rep, r1, r2, theta, dx, dy;



  for ( z=0; z<(int)imExt->dim.z; z++ )
  for ( y=0; y<(int)imExt->dim.y; y++ )
  for ( x=0; x<(int)imExt->dim.x; x++ ) {

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






void VT_Compute2DMaskedExtrema( vt_image *imRes,
              vt_image *imTheta,
              vt_image *imMask )
{
  char *proc = "VT_Compute2MaskedDExtrema";
  float ***theR = (float ***)imRes->array;
  float ***theTHETA = (float ***)imTheta->array;

  int x, y, z;
  double dx, dy;
  double r1, r2;

  switch ( imMask->type ) {
  default :
    fprintf( stderr, "%s: such type not handled in switch\n", proc );
    return;
  case SCHAR :
    {
      s8 ***theE = (s8 ***)imMask->array;
      for ( z=0; z<(int)imRes->dim.z; z++ )
      for ( y=0; y<(int)imRes->dim.y; y++ )
      for ( x=0; x<(int)imRes->dim.x; x++ ) {
    if ( theE[z][y][x] < 0 ) continue;
    dx = x + cos( theTHETA[z][y][x] );
    dy = y + sin( theTHETA[z][y][x] );
    r1 = _GetInterpolated2DValue( imRes, dx, dy, z );
    if ( theR[z][y][x] < r1 ) continue;

    dx = x - cos( theTHETA[z][y][x] );
    dy = y - sin( theTHETA[z][y][x] );
    r2 = _GetInterpolated2DValue( imRes, dx, dy, z );
    if ( theR[z][y][x] < r2 ) continue;

    theE[z][y][x] = 127;
      }
    }
    break;
  case UCHAR :
    {
      u8 ***theE = (u8 ***)imMask->array;
      for ( z=0; z<(int)imRes->dim.z; z++ )
      for ( y=0; y<(int)imRes->dim.y; y++ )
      for ( x=0; x<(int)imRes->dim.x; x++ ) {
    theE[z][y][x] = 0;
    if (  theR[z][y][x] <= 0.0 ) continue;

    dx = x + cos( theTHETA[z][y][x] );
    dy = y + sin( theTHETA[z][y][x] );
    r1 = _GetInterpolated2DValue( imRes, dx, dy, z );
    if ( theR[z][y][x] <= r1 ) continue;


    dx = x - cos( theTHETA[z][y][x] );
    dy = y - sin( theTHETA[z][y][x] );
    r2 = _GetInterpolated2DValue( imRes, dx, dy, z );
    if ( theR[z][y][x] <= r2 ) continue;

    theE[z][y][x] = 255;
      }
    }
    break;
  }

}













void VT_Compute2DResponse( vt_2Dimages *par, enumStructureColor color,
               double theta, double sigma )
{
  double coeff = theta;
  int x, y, z;

  s8 ***theE = (s8 ***)par->ime.array;
  float ***theR = (float ***)par->imr.array;
  float ***theTHETA1 = (float ***)par->hessien.imtheta1.array;

  double dx, dy, gradx, grady;
  double ps1, ps2;


  if ( coeff <= 0.0 ) coeff = 1.0;




  for ( z=0; z<(int)par->imx.dim.z; z++ )
  for ( y=0; y<(int)par->imx.dim.y; y++ )
  for ( x=0; x<(int)par->imx.dim.x; x++ ) {
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



void VT_Compute2DResponseFrangi( vt_2Dimages *par, double sigma )
{
  int x, y, z;

  s8 ***theE = (s8 ***)par->ime.array;
  float ***theR = (float ***)par->imr.array;
  float ***theVP1 = (float ***)par->hessien.imvp1.array;
  float ***theVP2 = (float ***)par->hessien.imvp2.array;

  vt_image imtmp;
  float ***theS;
  int alloc;

  double vp1, vp2;
  double beta=0.5;
  double beta2=2*beta*beta;
  double c=1.0;
  double c2=2*c*c;
  double R2,S2;
  double Hmax=0;

  VT_Image( &imtmp );
  VT_InitFromImage( &imtmp, &(par->imr), "tmp.inr", FLOAT);
  if ( VT_AllocImage( &imtmp ) != 1 ) {
    alloc=0;
    fprintf(stderr, "unable to allocate temporary image, taking default value c=1.0\n");
  }
  else {
    alloc=1;
    theS=(float ***)imtmp.array;
  }



  for ( z=0; z<(int)par->imx.dim.z; z++ )
  for ( y=0; y<(int)par->imx.dim.y; y++ )
  for ( x=0; x<(int)par->imx.dim.x; x++ ) {
    if ( theE[z][y][x] < 0 ) {
      theR[z][y][x] = 0.0;
      continue;
    }

      /* Detecteur de Frangi
       *
       */
      vp1=theVP1[z][y][x]; vp1=vp1*vp1; /* vp1 = val.propre 1 au carre */
      vp2=theVP2[z][y][x]; vp2=vp2*vp2; /* vp2 = val.propre 2 au carre */
      R2=vp2/vp1;
      S2=sigma*sigma*sigma*sigma*(vp1+vp2);
      if (alloc==0)
        theR[z][y][x] = exp(-R2/(beta2)) * (1-exp(-S2/c2));
      else {
        theR[z][y][x] = exp(-R2/(beta2));
        theS[z][y][x] = S2;
        if (Hmax < S2) Hmax = S2;
      }
  }

  if (alloc==1) {
    /* c=sqrt(Hmax)/2; */
    c2=Hmax/2;
    /* fprintf(stderr, "Parameter c=%f\n", sqrt(Hmax)/2); */
    for ( z=0; z<(int)par->imx.dim.z; z++ )
    for ( y=0; y<(int)par->imx.dim.y; y++ )
    for ( x=0; x<(int)par->imx.dim.x; x++ ) {
      if ( theE[z][y][x] < 0 ) {
        continue;
      }
      S2=theS[z][y][x];
      theR[z][y][x] *= 1-exp(-S2/c2);
    }

    VT_FreeImage( &imtmp );
  }

}













void VT_FilterOn2DEigenValues( vt_2Dimages *par, enumStructureColor color )
{
  int x, y, z;

  s8 ***theE = (s8 ***)par->ime.array;
  float ***theVP1 = (float ***)par->hessien.imvp1.array;
  float ***theVP2 = (float ***)par->hessien.imvp2.array;
  float ***theTHETA1 = (float ***)par->hessien.imtheta1.array;
  float ***theTHETA2 = (float ***)par->hessien.imtheta2.array;


  /* case blanc : la plus grande valeur propre doit etre negative
   */
  switch( color ) {

  case _WHITE_ :

    for ( z=0; z<(int)par->imx.dim.z; z++ )
    for ( y=0; y<(int)par->imx.dim.y; y++ )
    for ( x=0; x<(int)par->imx.dim.x; x++ ) {
      /*
      if ( ( fabs(theVP1[0][y][x]) > fabs(theVP2[0][y][x]) && theVP1[0][y][x] > 0.0 )
       || ( fabs(theVP1[0][y][x]) < fabs(theVP2[0][y][x]) && theVP2[0][y][x] > 0.0 ) ) {
      */
      if ( theVP1[z][y][x] > 0.0 ) {
    theE[z][y][x] = -128;
    theVP1[z][y][x] = theVP2[z][y][x] = 0.0;
    theTHETA1[z][y][x] = theTHETA2[z][y][x] = 0.0;
      }
    }
    break;

  case _BLACK_ :

    for ( z=0; z<(int)par->imx.dim.z; z++ )
    for ( y=0; y<(int)par->imx.dim.y; y++ )
    for ( x=0; x<(int)par->imx.dim.x; x++ ) {
      /*
      if ( ( fabs(theVP1[0][y][x]) > fabs(theVP2[0][y][x]) && theVP1[0][y][x] < 0.0 )
       || ( fabs(theVP1[0][y][x]) < fabs(theVP2[0][y][x]) && theVP2[0][y][x] < 0.0 ) ) {
      */
      if ( theVP1[z][y][x] < 0.0 ) {
    theE[z][y][x] = -128;
    theVP1[z][y][x] = theVP2[z][y][x] = 0.0;
    theTHETA1[z][y][x] = theTHETA2[z][y][x] = 0.0;
      }
    }
  }
}













/*
 * VP1 contient la plus valeur propre (en valeur absolue)
 * fabs( theVP1[z][y][x] ) >= fabs( theVP2[z][y][x] )
 *
 */


void VT_Compute2DEigenVectors( vt_2Dtensor *par,
                   vt_image *imMask )
{
  int x, y, z;
  double sum, prod, delta;
  double vp1, vp2;

  float ***theXX = (float ***)par->imxx.array;
  float ***theXY = (float ***)par->imxy.array;
  float ***theYY = (float ***)par->imyy.array;

  float ***theVP1 = (float ***)par->imvp1.array;
  float ***theVP2 = (float ***)par->imvp2.array;
  float ***theTHETA1 = (float ***)par->imtheta1.array;
  float ***theTHETA2 = (float ***)par->imtheta2.array;

  s8 ***theE = (s8 ***)NULL;




  if ( imMask != (vt_image *)NULL )
    theE = (s8 ***)imMask->array;

  for ( z=0; z<(int)par->imxx.dim.z; z++ )
  for ( y=0; y<(int)par->imxx.dim.y; y++ )
  for ( x=0; x<(int)par->imxx.dim.x; x++ ) {

    if ( theE != (s8 ***)NULL ) theE[z][y][x] = 0;

    sum  = theXX[z][y][x] + theYY[z][y][x];
    prod = theXX[z][y][x] * theYY[z][y][x] - theXY[z][y][x] * theXY[z][y][x];
    delta = sum*sum - 4.0 * prod;


    if ( delta >= 0.0 ) {

      vp1 = ( sum - sqrt( delta ) ) / 2.0;
      vp2 = ( sum + sqrt( delta ) ) / 2.0;

      /* float  atan2f   (float  y, float  x);
     The atan2() and atan2f() functions compute the principal value of the arc
     tangent of y/x, in the interval [-pi,pi] radians.  The sign of atan2() and
     atan2f() is determined by the sign of y.  The value of atan2(y,x) is com-
     puted as follows where f is the number of fraction bits associated with the
     data type.
     _____________________________________________________
     Value of Input Arguments    Angle Returned
     _____________________________________________________
     x = 0 or y/x > 2**(f+1)     pi/2 * (sign y)
     x > 0 and y/x <= 2**(f+1)   atan(y/x)
     x < 0 and y/x <= 2**(f+1)   pi * (sign y) + atan(y/x)
     _____________________________________________________

      */

      if ( fabs( vp1 ) > fabs( vp2 ) ) {
    theVP1[z][y][x] = vp1;
    theVP2[z][y][x] = vp2;
    theTHETA1[z][y][x] = atan2( (double)(vp1 - theXX[z][y][x]),
                    (double)theXY[z][y][x] );
    theTHETA2[z][y][x] = atan2( (double)(vp2 - theXX[z][y][x]),
                    (double)theXY[z][y][x] );
      } else {
    theVP1[z][y][x] = vp2;
    theVP2[z][y][x] = vp1;
    theTHETA1[z][y][x] = atan2( (double)(vp2 - theXX[z][y][x]),
                    (double)theXY[z][y][x] );
    theTHETA2[z][y][x] = atan2( (double)(vp1 - theXX[z][y][x]),
                     (double)theXY[z][y][x] );
      }



    } else {

      if ( theE != (s8 ***)NULL ) theE[z][y][x] = -128;
      theVP1[z][y][x] = theVP2[z][y][x] = 0.0;
      theTHETA1[z][y][x] = theTHETA2[z][y][x] = 0.0;

    }
  }
}















int VT_Filter2Dimages( vt_image *im, vt_2Dimages *par, float *theCoeffs )
{
  char *proc = "VT_Filter2Dimages";
  typeFilteringCoefficients filter[3];
  int theDim[3];

  theDim[0] = im->dim.x;
  theDim[1] = im->dim.y;
  theDim[2] = im->dim.z;


  /* initialisation des coefficients
   */

  initFilteringCoefficients( &(filter[0]) );
  initFilteringCoefficients( &(filter[1]) );
  initFilteringCoefficients( &(filter[2]) );
  filter[0].type = GABOR_YOUNG_2002;
  filter[1].type = GABOR_YOUNG_2002;
  filter[2].type = GABOR_YOUNG_2002;
  filter[0].coefficient = theCoeffs[0];
  filter[1].coefficient = theCoeffs[1];
  filter[2].coefficient = theCoeffs[2];
  filter[0].derivative = NODERIVATIVE;
  filter[1].derivative = NODERIVATIVE;
  filter[2].derivative = NODERIVATIVE;

  
  {
    int x, y, z;
    s8 ***theE = (s8 ***)par->ime.array;
    float ***theR = (float ***)par->imr.array;

    for ( z=0; z<(int)par->ime.dim.z; z++ )
    for ( y=0; y<(int)par->ime.dim.y; y++ )
    for ( x=0; x<(int)par->ime.dim.x; x++ ) {
      theE[z][y][x] = 0;
      theR[z][y][x] = 0.0;
    }
  }



  filter[0].derivative = NODERIVATIVE;
  filter[1].derivative = DERIVATIVE_0;

  if ( separableLinearFiltering( im->buf, im->type,
				 par->imx.buf, par->imx.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computating derivatives (.,0)\n", proc );
    return( -1 );
  }

  filter[0].derivative = DERIVATIVE_2;
  filter[1].derivative = NODERIVATIVE;

  if ( separableLinearFiltering( par->imx.buf, par->imx.type,
				 par->hessien.imxx.buf, par->hessien.imxx.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computating derivatives (2,.)\n", proc );
    return( -1 );
  }

  filter[0].derivative = DERIVATIVE_1;
  filter[1].derivative = NODERIVATIVE;

  if ( separableLinearFiltering( par->imx.buf, par->imx.type,
				 par->imx.buf, par->imx.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computating derivatives (1,.)\n", proc );
    return( -1 );
  }



  filter[0].derivative = DERIVATIVE_0;
  filter[1].derivative = NODERIVATIVE;

  if ( separableLinearFiltering( im->buf, im->type,
				 par->imy.buf, par->imy.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computating derivatives (0,.)\n", proc );
    return( -1 );
  }

  filter[0].derivative = NODERIVATIVE;
  filter[1].derivative = DERIVATIVE_2;

  if ( separableLinearFiltering( par->imy.buf, par->imy.type,
				 par->hessien.imyy.buf, par->hessien.imyy.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computating derivatives (.,2)\n", proc );
    return( -1 );
  }

  filter[0].derivative = NODERIVATIVE;
  filter[1].derivative = DERIVATIVE_1;

  if ( separableLinearFiltering( par->imy.buf, par->imy.type,
				 par->imy.buf, par->imy.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computating derivatives (.,1)\n", proc );
    return( -1 );
  }



  filter[0].derivative = DERIVATIVE_1;
  filter[1].derivative = DERIVATIVE_1;

  if ( separableLinearFiltering( im->buf, im->type,
				 par->hessien.imxy.buf, par->hessien.imxy.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computating derivatives (1,1)\n", proc );
    return( -1 );
  }

  return( 1 );
}



int VT_Filter2Dimages2ndOrder( vt_image *im, vt_2Dimages *par, float *theCoeffs )
{
  char *proc = "VT_Filter2Dimages2ndOrder";
  typeFilteringCoefficients filter[3];
  int theDim[3];

  theDim[0] = im->dim.x;
  theDim[1] = im->dim.y;
  theDim[2] = im->dim.z;



  /* initialisation des coefficients
   */

  initFilteringCoefficients( &(filter[0]) );
  initFilteringCoefficients( &(filter[1]) );
  initFilteringCoefficients( &(filter[2]) );
  filter[0].type = GABOR_YOUNG_2002;
  filter[1].type = GABOR_YOUNG_2002;
  filter[2].type = GABOR_YOUNG_2002;
  filter[0].coefficient = theCoeffs[0];
  filter[1].coefficient = theCoeffs[1];
  filter[2].coefficient = theCoeffs[2];
  filter[0].derivative = NODERIVATIVE;
  filter[1].derivative = NODERIVATIVE;
  filter[2].derivative = NODERIVATIVE;


  {
    int x, y, z;
    s8 ***theE = (s8 ***)par->ime.array;
    float ***theR = (float ***)par->imr.array;

    for ( z=0; z<(int)par->ime.dim.z; z++ )
    for ( y=0; y<(int)par->ime.dim.y; y++ )
    for ( x=0; x<(int)par->ime.dim.x; x++ ) {
      theE[z][y][x] = 0;
      theR[z][y][x] = 0.0;
    }
  }


  filter[0].derivative = DERIVATIVE_2;
  filter[1].derivative = DERIVATIVE_0;

  if ( separableLinearFiltering( im->buf, im->type,
				 par->hessien.imxx.buf, par->hessien.imxx.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computating derivatives (2,0)\n", proc );
    return( -1 );
  }
  filter[0].derivative = DERIVATIVE_0;
  filter[1].derivative = DERIVATIVE_2;

  if ( separableLinearFiltering( im->buf, im->type,
				 par->hessien.imyy.buf, par->hessien.imyy.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computating derivatives (0,2)\n", proc );
    return( -1 );
  }
  filter[0].derivative = DERIVATIVE_1;
  filter[1].derivative = DERIVATIVE_1;

  if ( separableLinearFiltering( im->buf, im->type,
				 par->hessien.imxy.buf, par->hessien.imxy.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error in computating derivatives (1,1)\n", proc );
    return( -1 );
  }

  return( 1 );
}






















int VT_Alloc2Dtensor ( vt_2Dtensor *par, vt_image *im, char *genericname )
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

  return( 1 );
}



int VT_Alloc2Dimages ( vt_2Dimages *par, vt_image *im, char *genericname )
{
  char name[256];

  sprintf( name, "%s.imx.inr", genericname );
  VT_InitFromImage( &(par->imx),  im, name, (int)FLOAT );
  sprintf( name, "%s.imy.inr", genericname );
  VT_InitFromImage( &(par->imy),  im, name,  (int)FLOAT );


  sprintf( name, "%s.imr.inr", genericname );
  VT_InitFromImage( &(par->imr),  im, name,  (int)FLOAT );
  sprintf( name, "%s.ime.inr", genericname );
  VT_InitFromImage( &(par->ime),  im, name,  (int)SCHAR );

  if ( VT_AllocImage( &(par->imx) ) != 1 ) {
    return( -1 );
    }
  if ( VT_AllocImage( &(par->imy) ) != 1 ) {
    VT_FreeImage( &(par->imx) );
    return( -1 );
  }
  if ( VT_Alloc2Dtensor( &(par->hessien), im, genericname ) != 1 ) {
    VT_FreeImage( &(par->imx) );
    VT_FreeImage( &(par->imy) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imr) ) != 1 ) {
    VT_FreeImage( &(par->imx) );
    VT_FreeImage( &(par->imy) );
    VT_Free2Dtensor( &(par->hessien) );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->ime) ) != 1 ) {
    VT_FreeImage( &(par->imx) );
    VT_FreeImage( &(par->imy) );
    VT_Free2Dtensor( &(par->hessien) );
    VT_FreeImage( &(par->imr) );
    return( -1 );
  }


  return( 1 );
}





void VT_Free2Dtensor( vt_2Dtensor *par )
{
  VT_FreeImage( &(par->imxx) );
  VT_FreeImage( &(par->imxy) );
  VT_FreeImage( &(par->imyy) );

  VT_FreeImage( &(par->imvp1) );
  VT_FreeImage( &(par->imvp2) );
  VT_FreeImage( &(par->imtheta1) );
  VT_FreeImage( &(par->imtheta2) );
}



void VT_Free2Dimages( vt_2Dimages *par )
{
  VT_FreeImage( &(par->imx) );
  VT_FreeImage( &(par->imy) );

  VT_Free2Dtensor( &(par->hessien) );

  VT_FreeImage( &(par->imr) );
  VT_FreeImage( &(par->ime) );

}


void VT_Write2Dtensor( vt_2Dtensor *par )
{
  VT_WriteInrimage( &(par->imxx) );
  VT_WriteInrimage( &(par->imxy) );
  VT_WriteInrimage( &(par->imyy) );

  VT_WriteInrimage( &(par->imvp1) );
  VT_WriteInrimage( &(par->imvp2) );
  VT_WriteInrimage( &(par->imtheta1) );
  VT_WriteInrimage( &(par->imtheta2) );
}


void VT_Write2Dimages( vt_2Dimages *par )
{
  VT_WriteInrimage( &(par->imx) );
  VT_WriteInrimage( &(par->imy) );

  VT_Write2Dtensor( &(par->hessien) );

  VT_WriteInrimage( &(par->imr) );
  VT_WriteInrimage( &(par->ime) );

}





int VT_Alloc3DImres( vt_3Dimres *par, vt_image *im, char *genericname,
             int alloc_phi )
{
  char *proc = "VT_Alloc3DImres";
  char name[256];

  sprintf( name, "%s.theta.inr", genericname );
  VT_InitFromImage( &(par->imTheta), im, name, (int)FLOAT );
  sprintf( name, "%s.phi.inr", genericname );
  VT_InitFromImage( &(par->imPhi), im, name, (int)FLOAT );
  sprintf( name, "%s.rep.inr", genericname );
  VT_InitFromImage( &(par->imRep), im, name, (int)FLOAT );
  sprintf( name, "%s.scale.inr", genericname );
  VT_InitFromImage( &(par->imScale), im, name, (int)FLOAT );

  if ( VT_AllocImage( &(par->imTheta) ) != 1 ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to allocate theta image\n", proc );
    }
    return( -1 );
  }
  if ( alloc_phi ) {
    if ( VT_AllocImage( &(par->imPhi) ) != 1 ) {
      if ( _verbose_ ) {
    fprintf( stderr, "%s: unable to allocate phi image\n", proc );
      }
      VT_FreeImage( &(par->imTheta) );
      return( -1 );
    }
  }
  if ( VT_AllocImage( &(par->imRep) ) != 1 ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to allocate response image\n", proc );
    }
    VT_FreeImage( &(par->imPhi) );
    VT_FreeImage( &(par->imTheta) );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imScale) ) != 1 ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to allocate scale image\n", proc );
    }
    VT_FreeImage( &(par->imPhi) );
    VT_FreeImage( &(par->imTheta) );
    VT_FreeImage( &(par->imRep  ) );
    return( -1 );
  }

  return( 1 );
}




void VT_Write3DImres( vt_3Dimres *par, int write_phi )
{
  VT_WriteInrimage( &(par->imTheta) );
  if ( write_phi )
    VT_WriteInrimage( &(par->imPhi) );
  VT_WriteInrimage( &(par->imRep) );
  VT_WriteInrimage( &(par->imScale) );
}


void VT_Write3DImresAngles( vt_3Dimres *par, int write_phi )
{
  VT_WriteInrimage( &(par->imTheta) );
  if ( write_phi )
    VT_WriteInrimage( &(par->imPhi) );
}



void VT_Free3DImres( vt_3Dimres *par )
{
  VT_FreeImage( &(par->imTheta) );
  VT_FreeImage( &(par->imPhi) );
  VT_FreeImage( &(par->imRep) );
  VT_FreeImage( &(par->imScale) );
}
