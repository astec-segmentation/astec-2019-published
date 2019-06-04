/*************************************************************************
 * vt_tube3D.c -
 *
 * $Id: vt_tube3D.c,v 1.12 2006/05/16 09:33:34 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Mon Sep 11 18:20:24 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <vt_tube3D.h>
#include <transfo.h>
#include <eigens.h>
#include <convert.h>

#include <linearFiltering.h>

static int _borderLength_[3] = { 20, 20, 20 };
static int _verbose_ = 1;
static int _debug_ = 1;



int VT_GetVerboseInVtTube3D( )
{
  return( _verbose_ );
}

void VT_SetVerboseInVtTube3D( int v )
{
  _verbose_ = v;
}

void VT_IncrementVerboseInVtTube3D(  )
{
  _verbose_ ++;
}

void VT_DecrementVerboseInVtTube3D(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}

int VT_GetDebugInVtTube3D( )
{
  return( _debug_ );
}

void VT_SetDebugInVtTube3D( int v )
{
  _debug_ = v;
}

void VT_IncrementDebugInVtTube3D(  )
{
  _debug_ ++;
}

void VT_DecrementDebugInVtTube3D(  )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}





int VT_Compute3DMultiScale( vt_image *theIm,
			    vt_3Dimres *imsRes,
			    double scale1,
			    double scale2,
			    int nbscales,
			    enumStructureColor color,
			    methodType mode ) 
{
  char *proc = "VT_Compute3DMultiScale";
  int s, nbs = nbscales;
  double scale;

  vt_3Dimages ims3D;
  vt_2Dimauxs ims2D;

  double theta = 1.0;

  int    nbCirclePoints;
  typeCirclePoint *thePts = (typeCirclePoint*)NULL;

  float theCoeff[3];

  int slice;

  typeResponseInSlice aux2D;

  int x, y, z;
  float ***maxRep = (float***)imsRes->imRep.array;
  float ***maxTht = (float***)imsRes->imTheta.array;
  float ***maxPhi = (float***)imsRes->imPhi.array;
  float ***maxScl = (float***)imsRes->imScale.array;

  float ***theRep;
  float ***theTht;
  float ***thePhi;





  if ( VT_Alloc3Dimages( &ims3D, theIm, proc ) != 1 ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to allocate 3D auxiliary images\n", proc );
      return( -1 );
    }
  }


  theRep = (float***)ims3D.imzz.array;
  theTht = (float***)ims3D.tmp0.array;
  thePhi = (float***)ims3D.tmp1.array;


  if ( VT_Alloc2Dimauxs( &ims2D, theIm, proc ) != 1 ) {
    VT_Free3DImages( &ims3D );
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to allocate 2D auxiliary images\n", proc );
      return( -1 );
    }
  }






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
      /* Calcul des images 3D suivantes
        tmp0 : . . 0
        tmp1 : . . 1
        imx  : 1 0 0
        imy  : 0 1 0
        imz  : 0 0 1
        imzz : 0 0 2
      */
      if ( VT_Filter3Dimages( theIm, &ims3D, theCoeff ) != 1 ) {
		VT_Free2Dimauxs( &ims2D );
		VT_Free3DImages( &ims3D );
		if ( _verbose_ ) {
		  fprintf( stderr, "%s: unable to 3D filter\n", proc );
		  return( -1 );
		}
      }
	}
	else {
      /* Calcul des images 3D suivantes
        tmp0 : . . 0
        tmp1 : . . 1
        imzz : 0 0 2
	  */
	  if ( VT_Filter3Dimages2ndOrder( theIm, &ims3D, theCoeff ) != 1 ) {
		VT_Free2Dimauxs( &ims2D );
		VT_Free3DImages( &ims3D );
		if ( _verbose_ ) {
		  fprintf( stderr, "%s: unable to 3D filter\n", proc );
		  return( -1 );
		}
      }
	}

    
    thePts = VT_BuildCircle( theta*scale, &nbCirclePoints );
    if ( thePts == (typeCirclePoint*)NULL ) {
      VT_Free2Dimauxs( &ims2D );
      VT_Free3DImages( &ims3D );
      if ( _verbose_ ) {
	fprintf( stderr, "%s: unable to allocate points list\n", proc );
	return( -1 );
      }
    }





    for ( slice = 0; slice < (int)theIm->dim.z; slice ++ ) {


      /* Calcul des images 2D suivantes
	 imxx : 2 0 0
	 imyy : 0 2 0
	 imxy : 1 1 0
	 imxz : 1 0 1
	 imyz : 0 1 1
      */
      if ( VT_Filter2Dimauxs( &ims3D, &ims2D, theCoeff, slice ) != 1 ) {
	free( thePts );
	VT_Free2Dimauxs( &ims2D );
	VT_Free3DImages( &ims3D );
	if ( _verbose_ ) {
	  fprintf( stderr, "%s: unable to 2D filter\n", proc );
	  return( -1 );
	}
      }


      aux2D.theXX = ((float***)(ims2D.imxx.array))[0];
      aux2D.theYY = ((float***)(ims2D.imyy.array))[0];
      aux2D.theXY = ((float***)(ims2D.imxy.array))[0];

      aux2D.theXZ = ((float***)(ims2D.imxz.array))[0];
      aux2D.theYZ = ((float***)(ims2D.imyz.array))[0];

      aux2D.theZZ = ((float***)(ims3D.imzz.array))[slice];

      aux2D.theX  = (float***)(ims3D.imx.array);
      aux2D.theY  = (float***)(ims3D.imy.array);
      aux2D.theZ  = (float***)(ims3D.imz.array);
      
      aux2D.theRep   = theRep[slice];
      aux2D.theTheta = theTht[slice];
      aux2D.thePhi   = thePhi[slice];

      VT_Compute3DresponseInSlice( &aux2D, thePts, nbCirclePoints,
				   theta*scale, theIm->dim.x, theIm->dim.y, 
				   theIm->dim.z, slice, color, mode );
      

    } /* for ( slice = 0; slice < theIm->dim.z; slice ++ ) */

    
    free( thePts );
    thePts = (typeCirclePoint*)NULL;

    /* on recopie le resultat dans les images resultats
     */
    
    if ( s == 0 ) {
      for ( z = 0; z < (int)theIm->dim.z; z ++ )
      for ( y = 0; y < (int)theIm->dim.y; y ++ )
      for ( x = 0; x < (int)theIm->dim.x; x ++ ) {
	if ( theRep[z][y][x] > 0.0 ) {
	  maxRep[z][y][x] = scale * theRep[z][y][x];
	  maxTht[z][y][x] = theTht[z][y][x];
	  maxPhi[z][y][x] = thePhi[z][y][x];
	  maxScl[z][y][x] = scale * theIm->siz.x;
	} else {
	  maxRep[z][y][x] = maxTht[z][y][x] = maxPhi[z][y][x] = 0.0;
	  maxScl[z][y][x] = 0.0;
	}
      }
    } else {
      for ( z = 0; z < (int)theIm->dim.z; z ++ )
      for ( y = 0; y < (int)theIm->dim.y; y ++ )
      for ( x = 0; x < (int)theIm->dim.x; x ++ ) {
	if (mode==KRISSIAN) {
	  if ( scale * theRep[z][y][x] > maxRep[z][y][x] ) {
		maxRep[z][y][x] = scale * theRep[z][y][x];
		maxTht[z][y][x] = theTht[z][y][x];
		maxPhi[z][y][x] = thePhi[z][y][x];
		maxScl[z][y][x] = scale * theIm->siz.x;
	  }
	}
	else { /*FRANGI*/
	  if ( theRep[z][y][x] > maxRep[z][y][x] ) {
		maxRep[z][y][x] = theRep[z][y][x];
		maxTht[z][y][x] = theTht[z][y][x];
		maxPhi[z][y][x] = thePhi[z][y][x];
		maxScl[z][y][x] = scale * theIm->siz.x;
	  }
	}
      }
    }

  } /* for ( s = 0; s < nbs; s ++ ) */


  VT_Free3DImages( &ims3D );
  VT_Free2Dimauxs( &ims2D );
  return( 1 );
}
			   



























/* calcul des extrema 3D pour des tubes
 */

void VT_Compute3DExtrema( vt_3Dimres *imsRes,
			  vt_image *imExt )
{
  int x, y, z;
  
  float ***theExt = (float***)(imExt->array);
  float ***theRep = (float***)(imsRes->imRep.array);
  float ***theTht = (float***)(imsRes->imTheta.array);
  float ***thePhi = (float***)(imsRes->imPhi.array);

  double v1[3], v2[3], v3[3];
  
  double rep, r;

  double theta, phi;

  double vx=0.0, vy=0.0, vz=0.0;
  
  double c = sqrt(2.0)/2.0;

  double rx, ry, rz;
  double dx, dy, dz;
  int i, ix, iy, iz;
  double dxdy, dxdz, dydz, dxdydz;
  double v4, v5, v6;





  for ( z=0; z<(int)imExt->dim.z; z++ )
  for ( y=0; y<(int)imExt->dim.y; y++ )
  for ( x=0; x<(int)imExt->dim.x; x++ ) {

    theta = (double)theTht[z][y][x];
    phi   = (double)thePhi[z][y][x];
    
    theExt[z][y][x] = 0.0;
    if ( theRep[z][y][x] <= 0.0 ) continue;

    rep = theRep[z][y][x];



    /* v1 est le vecteur directeur du vaisseau
       v2 et v3 sont deux vecteurs orthogonaux
    */
    SphericalAnglesToUnitsVectors( theta, phi, v1, v2, v3 );
    
    for ( i=0; i<8; i++ ) {

      switch( i ) {
      default :
      case 0 :
	vx = v2[0];    vy = v2[1];    vz = v2[2];    break;
      case 1 :
	vx = -v2[0];   vy = -v2[1];   vz = -v2[2];   break;
      case 2 :
	vx = v3[0];    vy = v3[1];    vz = v3[2];    break;
      case 3 :
	vx = -v3[0];   vy = -v3[1];   vz = -v3[2];   break;
      case 4 :
	vx =  c * v2[0] + c * v3[0];
	vy =  c * v2[1] + c * v3[1];
	vz =  c * v2[2] + c * v3[2];
	break;
      case 5 :
	vx = -c * v2[0] + c * v3[0];
	vy = -c * v2[1] + c * v3[1];
	vz = -c * v2[2] + c * v3[2];
	break;
       case 6 :
	vx =  c * v2[0] - c * v3[0];
	vy =  c * v2[1] - c * v3[1];
	vz =  c * v2[2] - c * v3[2];
	break;
      case 7 :
	vx = -c * v2[0] - c * v3[0];
	vy = -c * v2[1] - c * v3[1];
	vz = -c * v2[2] - c * v3[2];
	break;
      }

      rx = x + vx;
      ix = (int)rx;
      if ( rx <= 0.0 || ix >= (int)imExt->dim.x-1 ) {
	i = 10;
	continue;
      }
      ry = y + vy;
      iy = (int)ry;
      if ( ry <= 0.0 || iy >= (int)imExt->dim.y-1 ) {
	i = 10;
	continue;
      }
      rz = z + vz;
      iz = (int)rz;
      if ( rz <= 0.0 || iz >= (int)imExt->dim.z-1 ) {
	i = 10;
	continue;
      }


      /* si on veut que tous les points autour
	 aient une reponse positive
	 C'est tres (trop?) restrictif
      */
      if ( 0 ) {
	if ( theRep[iz+1][iy+1][ix+1] <= 0.0 ||
	     theRep[iz+1][iy+1][ix  ] <= 0.0 ||
	     theRep[iz+1][iy  ][ix+1] <= 0.0 ||
	     theRep[iz+1][iy  ][ix  ] <= 0.0 ||
	     theRep[iz  ][iy+1][ix+1] <= 0.0 ||
	     theRep[iz  ][iy+1][ix  ] <= 0.0 ||
	     theRep[iz  ][iy  ][ix+1] <= 0.0 ||
	     theRep[iz  ][iy  ][ix  ] <= 0.0 ) {
	  i = 10;
	  continue;
	}
      }
      
      

      dx = rx - ix;
      dy = ry - iy;
      dz = rz - iz;

      dxdy = dx*dy;
      dxdz = dx*dz;
      dydz = dy*dz;
      dxdydz = dxdy*dz;

      v6 = dxdz-dxdydz;
      v5 = dxdy-dxdydz;
      v4 = dx-dxdy-v6;

      r = 0;
      r += dxdydz        * theRep[iz+1][iy+1][ix+1];
      r += (dydz-dxdydz) * theRep[iz+1][iy+1][ix  ];
      r += v6            * theRep[iz+1][iy  ][ix+1];
      r += (dz-dydz-v6)  * theRep[iz+1][iy  ][ix  ];
      r += v5            * theRep[iz  ][iy+1][ix+1];
      r += (dy-dydz-v5)  * theRep[iz  ][iy+1][ix  ];
      r += v4            * theRep[iz  ][iy  ][ix+1];
      r += (1-dy-dz+dydz-v4) * theRep[iz  ][iy][ix];
      

      
      if ( r >= rep ) {
	i = 10;
	continue;
      }
      
    } /* fin de la boucle sur les 4 directions
       */
    if ( i < 10 ) theExt[z][y][x] = rep;
  }
}




















void VT_Compute3DresponseInSlice( typeResponseInSlice *aux,
				  typeCirclePoint *thePts,
				  int nbPts,
				  double radius,
				  int dimx,
				  int dimy,
				  int dimz,
				  int slice,
				  enumStructureColor color,
				  methodType mode )
{
  float **theXX = aux->theXX;
  float **theYY = aux->theYY;
  float **theZZ = aux->theZZ;
  float **theXY = aux->theXY;
  float **theXZ = aux->theXZ;
  float **theYZ = aux->theYZ;
  
  float ***theX = aux->theX;
  float ***theY = aux->theY;
  float ***theZ = aux->theZ;

  float **theRep   = aux->theRep;
  float **theTheta = aux->theTheta;
  float **thePhi   = aux->thePhi;

  int x, y;

  double hessien[9];
  double valprop[3];
  double vecprop[9];
  double v[3];


  int n;
  double vx, vy, vz;
  double gx, gy, gz;
  double rx, ry, rz;
  double dx, dy, dz;
  int ix, iy, iz;
  double dxdy, dxdz, dydz, dxdydz;
  double v4, v5, v6;

  int nbPosPts;
  int nbNegPts;
  int nbValPts;
  /*
  int minValPts = (int)(0.5 + 2.0 * (double)nbPts / 3.0);
  */
  double sumValPts;
  double sumPosPts;
  double sumNegPts;

  double alpha, beta, c, Ra, Rb, S;


  double theta, phi;

  vt_image imtmp;
  float ***theS = (float***)NULL;
  double Hmax=0;
  int alloc = 0;

  if (mode==FRANGI) {
	alpha = 0.5;
	beta  = 0.5;
	c     = 1;
    VT_Image( &imtmp );
    VT_InitImage( &imtmp, "tmp.inr", dimx, dimy, 1, FLOAT);
    if ( VT_AllocImage( &imtmp ) != 1 ) {
      alloc=0;
      fprintf(stderr, "unable to allocate temporary image, taking default value c=1.0\n");
    }
    else {
      alloc=1;
      theS=(float ***)imtmp.array;
    }
  }

  for ( y = 0; y < dimy; y ++ )
  for ( x = 0; x < dimx; x ++ ) {



    /* calcul des valeurs et vecteurs propres
     */
    hessien[0] = theXX[y][x];
    hessien[1] = hessien[3] = theXY[y][x];
    hessien[2] = hessien[6] = theXZ[y][x];
    hessien[4] = theYY[y][x];
    hessien[5] = hessien[7] = theYZ[y][x];
    hessien[8] = theZZ[y][x];

	
    
    theRep[y][x] = theTheta[y][x] = thePhi[y][x] = 0;



    if ( _ComputeEigensOfSymetricSquareMatrix( hessien, valprop, vecprop, 3 ) != 1 ) {
      theRep[y][x] = theTheta[y][x] = thePhi[y][x] = 0.0;
      continue;
    }
    if ( _SortEigensInAbsIncreasingOrder( valprop, vecprop, 3 ) != 1 ) {
      theRep[y][x] = theTheta[y][x] = thePhi[y][x] = 0.0;
      continue;
    }




    switch( color ) {
    case _WHITE_ :

      /* le point appartient a un tube blanc sur fond
	 noir si les deux plus grandes sont negatives
	 et beaucoup plus grandes que la plus petite
	 
	 On a fabs(valprop[0]) <= fabs(valprop[1]) <= fabs(valprop[2])
	 
	 On veut :
	 1.        valprop[1] < 0         &&   valprop[2] < 0 
	 2.        fabs(valprop[0]) * 2   <=   fabs(valprop[1]) 	(si KRISSIAN)
	  */
      
      if ( valprop[1] >= 0 || valprop[2] >= 0 ) {
	theRep[y][x] = theTheta[y][x] = thePhi[y][x] = 0.0;
	continue;
      }

      break;

    case _BLACK_ :

      /* le point appartient a un tube noir sur fond
	 blanc si les deux plus grandes sont positives
	 et beaucoup plus grandes que la plus petite
	 
	 On a fabs(valprop[0]) <= fabs(valprop[1]) <= fabs(valprop[2])
	 
	 On veut :
	 1.        valprop[1] > 0         &&   valprop[2] > 0
	 2.        fabs(valprop[0]) * 2   <=   fabs(valprop[1])		(si KRISSIAN)
	 
      */
      if ( valprop[1] <= 0 || valprop[2] <= 0 ) {
	theRep[y][x] = theTheta[y][x] = thePhi[y][x] = 0.0;
	continue;
      }
      
      break;

    }


	if (mode==KRISSIAN) {	/*KRISSIAN mode*/

      if ( fabs(valprop[0]) * 2 > fabs(valprop[1]) )  {
		theRep[y][x] = theTheta[y][x] = thePhi[y][x] = 0.0;
		continue;
      }



      /* OK, ici on regarde le point 
         le vecteur de la direction du vaisseau est celui
         associe a la plus petite valeur propre
         ie vecprop[0,3,6]
         les deux autres sont vecprop[1,4,7] et vecprop[2,5,8]
      */

      nbPosPts = nbNegPts = nbValPts = 0;
      sumValPts = sumPosPts = sumNegPts = 0;

	  for ( n=0; n<nbPts; n++ ) {

		thePts[n].f = 0;
		thePts[n].v = 0;

		vx = thePts[n].c2 * vecprop[1] + thePts[n].c3 * vecprop[2];
		rx = x + radius * vx;
		ix = (int)rx;
		if ( rx <= 0.0 || ix >= dimx-1 ) {
		  continue;
		}

		vy = thePts[n].c2 * vecprop[4] + thePts[n].c3 * vecprop[5];
		ry = y + radius * vy;
		iy = (int)ry;
		if ( ry <= 0.0 || iy >= dimy-1 ) {
		  continue;
		}

		vz = thePts[n].c2 * vecprop[7] + thePts[n].c3 * vecprop[8];
		rz = slice + radius * vz;
		iz = (int)rz;
		if ( rz <= 0.0 || iz >= dimz-1 ) {
		  continue;
		}



		dx = rx - ix;
		dy = ry - iy;
		dz = rz - iz;

		dxdy = dx*dy;
		dxdz = dx*dz;
		dydz = dy*dz;
		dxdydz = dxdy*dz;

		v6 = dxdz-dxdydz;
		v5 = dxdy-dxdydz;
		v4 = dx-dxdy-v6;
      
		gx = 0;
		gx += dxdydz        * theX[iz+1][iy+1][ix+1];
		gx += (dydz-dxdydz) * theX[iz+1][iy+1][ix  ];
		gx += v6            * theX[iz+1][iy  ][ix+1];
		gx += (dz-dydz-v6)  * theX[iz+1][iy  ][ix  ];
		gx += v5            * theX[iz  ][iy+1][ix+1];
		gx += (dy-dydz-v5)  * theX[iz  ][iy+1][ix  ];
		gx += v4            * theX[iz  ][iy  ][ix+1];
		gx += (1-dy-dz+dydz-v4) * theX[iz  ][iy][ix];

		gy = 0;
		gy += dxdydz        * theY[iz+1][iy+1][ix+1];
		gy += (dydz-dxdydz) * theY[iz+1][iy+1][ix  ];
		gy += v6            * theY[iz+1][iy  ][ix+1];
		gy += (dz-dydz-v6)  * theY[iz+1][iy  ][ix  ];
		gy += v5            * theY[iz  ][iy+1][ix+1];
		gy += (dy-dydz-v5)  * theY[iz  ][iy+1][ix  ];
		gy += v4            * theY[iz  ][iy  ][ix+1];
		gy += (1-dy-dz+dydz-v4) * theY[iz  ][iy][ix];

		gz = 0;
		gz += dxdydz        * theZ[iz+1][iy+1][ix+1];
		gz += (dydz-dxdydz) * theZ[iz+1][iy+1][ix  ];
		gz += v6            * theZ[iz+1][iy  ][ix+1];
		gz += (dz-dydz-v6)  * theZ[iz+1][iy  ][ix  ];
		gz += v5            * theZ[iz  ][iy+1][ix+1];
		gz += (dy-dydz-v5)  * theZ[iz  ][iy+1][ix  ];
		gz += v4            * theZ[iz  ][iy  ][ix+1];
		gz += (1-dy-dz+dydz-v4) * theZ[iz  ][iy][ix];

		thePts[n].f = 1;
		/* cas _WHITE_
         */
		thePts[n].v = - (gx * vx + gy * vy + gz * vz);
		nbValPts ++;
		sumValPts += thePts[n].v;

        if ( thePts[n].v > 0 ) {
		  nbPosPts ++;
		  sumPosPts += thePts[n].v;
        }
		if ( thePts[n].v < 0 ) {
		  nbNegPts ++;
		  sumNegPts -= thePts[n].v;
		}
      
	  }

      /* on inverse pour le cas noir 
       */
      if ( color == _BLACK_ ) {
		int nb;
		double sum;
      
		for ( n=0; n<nbPts; n++ ) {
		  thePts[n].v = (- thePts[n].v);
		}
		sumValPts = (-sumValPts);

		sum = sumPosPts;
		sumPosPts = sumNegPts;
		sumNegPts = sum;

		nb = nbPosPts;
		nbPosPts = nbNegPts;
		nbNegPts = nb;

      }


      /* on a la liste des points
		on veut :
		1. suffisamment de points valides
           ie plus de 2/3 du nombre total de points
		2. suffisamment de points 'positifs'
           ie plus des 2/3 des points 'valides'
		3. une somme positive des points valides
		*/
      /*
      if ( nbValPts <= 0 || nbValPts < minValPts ) continue;
      if ( 3*nbPosPts < 2*nbValPts ) continue;
      */

      if ( nbPosPts <= nbNegPts ) continue;
      if ( sumValPts <= 0.0 ) continue;

      v[0] = vecprop[0];
      v[1] = vecprop[3];
	  v[2] = vecprop[6];
    
      /* comme reponse on prend la moyenne
       */

      theRep[y][x] = sumValPts / (double)nbValPts;
      UnitVectorToSphericalAngles( v, &theta, &phi );
      theTheta[y][x] = theta;
      thePhi[y][x]   = phi;

	}
	else { /*FRANGI mode*/

	  /* Frangi : 
	   * 	theRep = (1-exp(-A)) * exp(-B) * (1-exp(-C))
	   * where
	   * 	A=Ra*Ra/(2*alpha*alpha)
	   * 	B=Rb*Rb/(2*beta*beta)
	   * 	C=S*S/(2*c*c)
	   * with
	   * 	alpha, beta are thresholds fixed to 0.5
	   * 	c is threshold fixed to 1 (or 0.5*max{norm(Hess)} on the image)
	   * and
	   * 	Ra distinguishes between plane-like and line-like structures
	   * 	Rb distinguishes between blob-like and line-like structures
	   * 	S is a measure of the second order structureness
	   */
	  
	  Ra=fabs(valprop[1]/valprop[2]);
	  Rb=fabs(valprop[0])/sqrt(fabs(valprop[1]*valprop[2]));
	  S=sqrt(valprop[0]*valprop[0]+valprop[1]*valprop[1]+valprop[2]*valprop[2]);
	  
	  /* Multiscale consideration :
	   * 2nd order derivatives supposed to be multiplied by scale*scale
	   * Here scale==radius and this multiplication only impacts value of S
	   */
	  S=S*radius*radius; 
      if (alloc==0) {
	  
        theRep[y][x]=(1-exp(-Ra*Ra/(2*alpha*alpha))) * exp(-Rb*Rb/(2*beta*beta)) *
			 (1-exp(-S*S/(2*c*c)));

      } else {
        theRep[y][x]=(1-exp(-Ra*Ra/(2*alpha*alpha))) * exp(-Rb*Rb/(2*beta*beta));
        theS[0][y][x]=S;
        if(Hmax<S) Hmax=S;
      }

      v[0] = vecprop[0];
      v[1] = vecprop[3];
	  v[2] = vecprop[6];

      UnitVectorToSphericalAngles( v, &theta, &phi );
      theTheta[y][x] = theta;
      thePhi[y][x]   = phi;

	}
  }

  if (alloc==1) {
    c=Hmax/2;
    for ( y = 0; y < dimy; y ++ )
    for ( x = 0; x < dimx; x ++ ) {
      if (theRep[y][x] == 0.0)
        continue;
      S=theS[0][y][x];
      theRep[y][x] *= (1-exp(-S*S/(2*c*c)));
    }
    VT_FreeImage(&imtmp);
  }
}


















typeCirclePoint *VT_BuildCircle( double radius, int *n )
{
  char *proc = "VT_BuildCircle";
  int    nbCirclePoints = (int)(2.0 * 3.1415926536 * radius) + 1;
  typeCirclePoint *thePts = (typeCirclePoint*)NULL;
  int p;
  
  /* pour calculer la reponse,
     on se place sur un cercle de rayon radius=theta*sigma
     que l'on dicretise en N=(int)(2*pi*radius) + 1
     pour obtenir environ 1 point par unite curviligne
     soit v2 et v3 les "vecteurs" du cercle
     les points se calculent alors
     par cos(i*2*pi/N) v2 + sin(i*2*pi/N) v3
     pour i = 0 ... N-1
  */

  if ( nbCirclePoints < 4 ) nbCirclePoints = 4;

  if ( nbCirclePoints % 4 > 0 ) nbCirclePoints += 4 - nbCirclePoints % 4;

  *n = 0;
  thePts = (typeCirclePoint*)VT_Malloc( nbCirclePoints*sizeof(typeCirclePoint) );
  if( thePts == (typeCirclePoint*)NULL ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: allocation failed\n", proc );
    }
    return( (typeCirclePoint*)NULL );
  }
  *n = nbCirclePoints;


  thePts[0].c2 = 1.0;  
  thePts[0].c3 = 0.0;


  for ( p = 1; p < nbCirclePoints; p ++ ) {
    thePts[p].c2 = cos( p * 2.0 * 3.1415926536 / (double)nbCirclePoints );
    thePts[p].c3 = sin( p * 2.0 * 3.1415926536 / (double)nbCirclePoints );
  }

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "%s: build %d points with radius = %f\n", proc, nbCirclePoints, radius );
  }

  return( thePts );
}



























int VT_Filter3Dimages( vt_image *im, vt_3Dimages *par, float *theCoeffs )
{
  char *proc = "VT_Filter3Dimages";
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



  /* calcul
   */

  filter[0].derivative = NODERIVATIVE;
  filter[1].derivative = NODERIVATIVE;
  filter[2].derivative = DERIVATIVE_0;

  if ( _verbose_ >= 2 ) fprintf( stderr, " 3D filtering (.,.,0)\n" );

  if ( _verbose_ >= 2 ) {
    VT_PrintImage( stderr, im, "input image" );
    VT_PrintImage( stderr, &(par->tmp0), "par->tmp0" );
  }

  if ( separableLinearFiltering( im->buf, im->type,
				 par->tmp0.buf, par->tmp0.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (.,.,0)\n", proc );
    return( -1 );
  }



  filter[0].derivative = DERIVATIVE_0;
  filter[1].derivative = DERIVATIVE_1;
  filter[2].derivative = NODERIVATIVE;

  if( _verbose_ >= 2 ) fprintf( stderr, " 3D filtering (0,1,.)\n" );


  if ( separableLinearFiltering( par->tmp0.buf, par->tmp0.type,
				 par->imy.buf, par->imy.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (0,1,0)\n", proc );
    return( -1 );
  }



  filter[0].derivative = DERIVATIVE_1;
  filter[1].derivative = DERIVATIVE_0;
  filter[2].derivative = NODERIVATIVE;

  if( _verbose_ >= 2 ) fprintf( stderr, " 3D filtering (1,0,.)\n" );

  if ( separableLinearFiltering( par->tmp0.buf, par->tmp0.type,
				 par->imx.buf, par->imx.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (1,0,0)\n", proc );
    return( -1 );
  }



  filter[0].derivative = NODERIVATIVE;
  filter[1].derivative = NODERIVATIVE;
  filter[2].derivative = DERIVATIVE_1;

  if ( separableLinearFiltering( im->buf, im->type,
				 par->tmp1.buf, par->tmp1.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (.,.,1)\n", proc );
    return( -1 );
  }



  filter[0].derivative = DERIVATIVE_0;
  filter[1].derivative = DERIVATIVE_0;
  filter[2].derivative = NODERIVATIVE;

  if ( separableLinearFiltering( par->tmp1.buf, par->tmp1.type,
				 par->imz.buf, par->imz.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (0,0,1)\n", proc );
    return( -1 );
  }



  filter[0].derivative = DERIVATIVE_0;
  filter[1].derivative = DERIVATIVE_0;
  filter[2].derivative = DERIVATIVE_2;

  if ( separableLinearFiltering( im->buf, im->type,
				 par->imzz.buf, par->imzz.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (0,0,2)\n", proc );
    return( -1 );
  }



  return( 1 );
}




int VT_Filter3Dimages2ndOrder( vt_image *im, vt_3Dimages *par, float *theCoeffs )
{
  char *proc = "VT_Filter3Dimages2ndOrder";
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



  /* calcul
   */

  filter[0].derivative = NODERIVATIVE;
  filter[1].derivative = NODERIVATIVE;
  filter[2].derivative = DERIVATIVE_0;

  if ( _verbose_ >= 2 ) fprintf( stderr, " 3D filtering (.,.,0)\n" );

  if ( _verbose_ >= 2 ) {
    VT_PrintImage( stderr, im, "input image" );
    VT_PrintImage( stderr, &(par->tmp0), "par->tmp0" );
  }


  if ( separableLinearFiltering( im->buf, im->type,
				 par->tmp0.buf, par->tmp0.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (.,.,0)\n", proc );
    return( -1 );
  }



  filter[0].derivative = NODERIVATIVE;
  filter[1].derivative = NODERIVATIVE;
  filter[2].derivative = DERIVATIVE_1;

  if ( separableLinearFiltering( im->buf, im->type,
				 par->tmp1.buf, par->tmp1.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (.,.,1)\n", proc );
    return( -1 );
  }



  filter[0].derivative = DERIVATIVE_0;
  filter[1].derivative = DERIVATIVE_0;
  filter[2].derivative = DERIVATIVE_2;

  if ( separableLinearFiltering( im->buf, im->type,
				 par->imzz.buf, par->imzz.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (0,0,2)\n", proc );
    return( -1 );
  }



  return( 1 );
}



















int VT_Alloc3Dimages( vt_3Dimages *par, vt_image *im, char *genericname ) 
{
  char name[256];

  sprintf( name, "%s.imx.inr", genericname );
  VT_InitFromImage( &(par->imx), im, name, (int)FLOAT );  
  sprintf( name, "%s.imy.inr", genericname );
  VT_InitFromImage( &(par->imy), im, name, (int)FLOAT );  
  sprintf( name, "%s.imz.inr", genericname );
  VT_InitFromImage( &(par->imz), im, name, (int)FLOAT );


  sprintf( name, "%s.tmp0.inr", genericname );
  VT_InitFromImage( &(par->tmp0), im, name, (int)FLOAT );
  sprintf( name, "%s.tmp1.inr", genericname );
  VT_InitFromImage( &(par->tmp1), im, name, (int)FLOAT );

  
  sprintf( name, "%s.imzz.inr", genericname );
  VT_InitFromImage( &(par->imzz), im, name, (int)FLOAT );

  
  if ( VT_AllocImage( &(par->imx) ) != 1 ) {
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imy) ) != 1 ) {
    VT_FreeImage( &(par->imx) );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imz) ) != 1 ) {
    VT_FreeImage( &(par->imy) );
    VT_FreeImage( &(par->imx) );
    return( -1 );
  }


  if ( VT_AllocImage( &(par->tmp0) ) != 1 ) {
    VT_FreeImage( &(par->imz) );
    VT_FreeImage( &(par->imy) );
    VT_FreeImage( &(par->imx) );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->tmp1) ) != 1 ) {
    VT_FreeImage( &(par->tmp0) );
    VT_FreeImage( &(par->imz) );
    VT_FreeImage( &(par->imy) );
    VT_FreeImage( &(par->imx) );
    return( -1 );
  }


  if ( VT_AllocImage( &(par->imzz) ) != 1 ) {
    VT_FreeImage( &(par->tmp1) );
    VT_FreeImage( &(par->tmp0) );
    VT_FreeImage( &(par->imz) );
    VT_FreeImage( &(par->imy) );
    VT_FreeImage( &(par->imx) );
    return( -1 );
  }
  
  return( 1 );
}








void VT_Write3DImages( vt_3Dimages *par )
{
  VT_WriteInrimage( &(par->imx) );
  VT_WriteInrimage( &(par->imy) );
  VT_WriteInrimage( &(par->imz) );

  VT_WriteInrimage( &(par->tmp0) );
  VT_WriteInrimage( &(par->tmp1) );

  VT_WriteInrimage( &(par->imzz) );
}





void VT_Free3DImages( vt_3Dimages *par )
{
  VT_FreeImage( &(par->imx) );
  VT_FreeImage( &(par->imy) );
  VT_FreeImage( &(par->imz) );

  VT_FreeImage( &(par->tmp0) );
  VT_FreeImage( &(par->tmp1) );

  VT_FreeImage( &(par->imzz) );
}















int VT_Filter2Dimauxs( vt_3Dimages *ims, vt_2Dimauxs *par,
		       float *theCoeffs, int slice )
{
  char *proc = "VT_Filter2Dimauxs";
  typeFilteringCoefficients filter[3];
  int theDim[3];
  
  float ***theTmp0 = (float ***)(ims->tmp0.array);
  float ***theTmp1 = (float ***)(ims->tmp1.array);

  if ( ims->tmp0.type != FLOAT || ims->tmp1.type != FLOAT ) {
    return( -1 );
  }

  theDim[0] = ims->imzz.dim.x;
  theDim[1] = ims->imzz.dim.y;
  theDim[2] = 1;



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



  /* calcul
   */

  filter[0].derivative = DERIVATIVE_2;
  filter[1].derivative = DERIVATIVE_0;
  filter[2].derivative = NODERIVATIVE;

  if ( separableLinearFiltering( (void*)(&(theTmp0[slice][0][0])), ims->tmp0.type,
				 par->imxx.buf, par->imxx.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (2,0,0) in slice %d\n", 
	       proc, slice );
    return( -1 );
  }
  
  
  filter[0].derivative = DERIVATIVE_0;
  filter[1].derivative = DERIVATIVE_2;
  filter[2].derivative = NODERIVATIVE;

  if ( separableLinearFiltering( (void*)(&(theTmp0[slice][0][0])), ims->tmp0.type,
				 par->imyy.buf, par->imyy.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (0,2,0) in slice %d\n", 
	       proc, slice );
    return( -1 );
  }
  
  
  filter[0].derivative = DERIVATIVE_1;
  filter[1].derivative = DERIVATIVE_1;
  filter[2].derivative = NODERIVATIVE;
  
  if ( separableLinearFiltering( (void*)(&(theTmp0[slice][0][0])), ims->tmp0.type,
				 par->imxy.buf, par->imxy.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (1,1,0) in slice %d\n", 
	       proc, slice );
    return( -1 );
  }
  
  
  filter[0].derivative = DERIVATIVE_1;
  filter[1].derivative = DERIVATIVE_0;
  filter[2].derivative = NODERIVATIVE;

  if ( separableLinearFiltering( (void*)(&(theTmp1[slice][0][0])), ims->tmp1.type,
				 par->imxz.buf, par->imxz.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (1,0,1) in slice %d\n", 
	       proc, slice );
    return( -1 );
  }
  
  
  filter[0].derivative = DERIVATIVE_0;
  filter[1].derivative = DERIVATIVE_1;
  filter[2].derivative = NODERIVATIVE;
  
  if ( separableLinearFiltering( (void*)(&(theTmp1[slice][0][0])), ims->tmp1.type,
				 par->imyz.buf, par->imyz.type,
				 theDim, _borderLength_, filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (0,1,1) in slice %d\n", 
	       proc, slice );
    return( -1 );
  }
  
  
  
  return( 1 );
}









int VT_Alloc2Dimauxs( vt_2Dimauxs *par, vt_image *im, char *genericname ) 
{
  char name[256];

  sprintf( name, "%s.imxx.inr", genericname );
  VT_InitFromImage( &(par->imxx), im, name, (int)FLOAT );  
  par->imxx.dim.z = 1;
  sprintf( name, "%s.imyy.inr", genericname );
  VT_InitFromImage( &(par->imyy), im, name, (int)FLOAT );  
  par->imyy.dim.z = 1;
  sprintf( name, "%s.imxy.inr", genericname );
  VT_InitFromImage( &(par->imxy), im, name, (int)FLOAT );  
  par->imxy.dim.z = 1;



  sprintf( name, "%s.imxz.inr", genericname );
  VT_InitFromImage( &(par->imxz), im, name, (int)FLOAT );  
  par->imxz.dim.z = 1;
  sprintf( name, "%s.imyz.inr", genericname );
  VT_InitFromImage( &(par->imyz), im, name, (int)FLOAT );  
  par->imyz.dim.z = 1;




  if ( VT_AllocImage( &(par->imxx) ) != 1 ) {
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imyy) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imxy) ) != 1 ) {
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imxx) );
    return( -1 );
  }



  if ( VT_AllocImage( &(par->imxz) ) != 1 ) {
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imxx) );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imyz) ) != 1 ) {
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imxx) );
    return( -1 );
  }




  return( 1 );
}












void VT_Write2Dimauxs( vt_2Dimauxs *par )
{
  VT_WriteInrimage( &(par->imxx) );
  VT_WriteInrimage( &(par->imyy) );
  VT_WriteInrimage( &(par->imxy) );

  VT_WriteInrimage( &(par->imxz) );
  VT_WriteInrimage( &(par->imyz) );
}



void VT_Free2Dimauxs( vt_2Dimauxs *par )
{
  VT_FreeImage( &(par->imxx) );
  VT_FreeImage( &(par->imyy) );
  VT_FreeImage( &(par->imxy) );

  VT_FreeImage( &(par->imxz) );
  VT_FreeImage( &(par->imyz) );
}























/* calcule le repere (dans le monde reel)
   permettant de reechantillonner une coupe orthogonale
   a partir de la direction du vaisseau (dans le volume image)
   
   La direction du vaisseau doit etre la troisieme direction
   (qui sera associee au Z de la coupe reechantillonnee)
 */


void VT_ComputeNextRealFrame( double thetaInImage,
			      double phiInImage,
			      double *ptSize,
			      double *newRealFrame,
			      double *oldRealFrame )
{
  double t, p;
  double v1[3], v2[3], v3[3];
  double n;

  /* on cherche la direction du vaisseau 
     dans les coordonnees image du volume
  */
  SphericalAnglesToUnitVector( thetaInImage, phiInImage, v3 );
  /* on transforme cette direction en direction 'reelle'
     en multipliant par les dimensions du point
  */
  v3[0] *= ptSize[0];
  v3[1] *= ptSize[1];
  v3[2] *= ptSize[2];

  
  /* on renormalise 
  */
  n = sqrt( v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2]);
  v3[0] /= n;
  v3[1] /= n;
  v3[2] /= n;
  
  /* si il y a une ancienne base :
     - on verifie l'orientation de v3 
     - on calcule les nouveaux vecteurs de la base par produit vectoriel
       v2(new) = v3(new) * v1(old)
       v1(new) = v2(new) * v3(new)
  */
  if ( oldRealFrame != (double*)NULL ) {
    if ( v3[0]*oldRealFrame[2]+v3[1]*oldRealFrame[5]+v3[2]*oldRealFrame[8] < 0 ) {
      v3[0] = -v3[0];
      v3[1] = -v3[1];
      v3[2] = -v3[2];
    }
    newRealFrame[2] = v3[0];
    newRealFrame[5] = v3[1];
    newRealFrame[8] = v3[2];
    
    newRealFrame[1] = v3[1] * oldRealFrame[6] - v3[2] * oldRealFrame[3];
    newRealFrame[4] = v3[2] * oldRealFrame[0] - v3[0] * oldRealFrame[6];
    newRealFrame[7] = v3[0] * oldRealFrame[3] - v3[1] * oldRealFrame[0];
    n = sqrt( newRealFrame[1]*newRealFrame[1] +
	      newRealFrame[4]*newRealFrame[4] + newRealFrame[7]*newRealFrame[7] );
    newRealFrame[1] /= n;
    newRealFrame[4] /= n;
    newRealFrame[7] /= n;

    newRealFrame[0] = newRealFrame[4] * v3[2] - newRealFrame[7] * v3[1];
    newRealFrame[3] = newRealFrame[7] * v3[0] - newRealFrame[1] * v3[2];
    newRealFrame[6] = newRealFrame[1] * v3[1] - newRealFrame[4] * v3[0];
    n = sqrt( newRealFrame[0]*newRealFrame[0] +
	      newRealFrame[3]*newRealFrame[3] + newRealFrame[6]*newRealFrame[6] );
    newRealFrame[0] /= n;
    newRealFrame[3] /= n;
    newRealFrame[6] /= n;

    return;
  }

  /* s'il n'y a pas d'ancienne base
     on recalcule theta et phi et on on deduit la base
  */
  UnitVectorToSphericalAngles( v3, &t, &p );
  SphericalAnglesToUnitsVectors( t, p, v3, v1, v2 );
  
  newRealFrame[0] = v1[0];
  newRealFrame[3] = v1[1];
  newRealFrame[6] = v1[2];

  newRealFrame[1] = v2[0];
  newRealFrame[4] = v2[1];
  newRealFrame[7] = v2[2];

  newRealFrame[2] = v3[0];
  newRealFrame[5] = v3[1];
  newRealFrame[8] = v3[2];

}









void VT_ComputeRealFrame( double *frame,
			  double *previousFrame )
{
  double t, p;
  double v1[3], v2[3], v3[3];
  double n;

  /* la direction du vaisseau est le dernier vecteur
     du repere, il est cense etre normalise
  */
  v3[0] = frame[2];
  v3[1] = frame[5];
  v3[2] = frame[8];
  if ( 0 ) {
    n = sqrt( v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2]);
    v3[0] /= n;
    v3[1] /= n;
    v3[2] /= n;
  }


  /* si il y a une ancienne base :
     - on verifie l'orientation de v3 
     - on calcule les nouveaux vecteurs de la base par produit vectoriel
       v2(new) = v3(new) * v1(old)
       v1(new) = v2(new) * v3(new)
  */
  if ( previousFrame != (double*)NULL ) {
    if ( v3[0]*previousFrame[2]+v3[1]*previousFrame[5]+v3[2]*previousFrame[8] < 0 ) {
      v3[0] = -v3[0];
      v3[1] = -v3[1];
      v3[2] = -v3[2];
    }
    frame[2] = v3[0];
    frame[5] = v3[1];
    frame[8] = v3[2];
    
    frame[1] = v3[1] * previousFrame[6] - v3[2] * previousFrame[3];
    frame[4] = v3[2] * previousFrame[0] - v3[0] * previousFrame[6];
    frame[7] = v3[0] * previousFrame[3] - v3[1] * previousFrame[0];
    n = sqrt( frame[1]*frame[1] +
	      frame[4]*frame[4] + frame[7]*frame[7] );
    frame[1] /= n;
    frame[4] /= n;
    frame[7] /= n;

    frame[0] = frame[4] * v3[2] - frame[7] * v3[1];
    frame[3] = frame[7] * v3[0] - frame[1] * v3[2];
    frame[6] = frame[1] * v3[1] - frame[4] * v3[0];
    n = sqrt( frame[0]*frame[0] +
	      frame[3]*frame[3] + frame[6]*frame[6] );
    frame[0] /= n;
    frame[3] /= n;
    frame[6] /= n;

    return;
  }

  /* s'il n'y a pas d'ancienne base
     on recalcule theta et phi et on on deduit la base
  */
  UnitVectorToSphericalAngles( v3, &t, &p );
  SphericalAnglesToUnitsVectors( t, p, v3, v1, v2 );
  
  frame[0] = v1[0];
  frame[3] = v1[1];
  frame[6] = v1[2];

  frame[1] = v2[0];
  frame[4] = v2[1];
  frame[7] = v2[2];

  frame[2] = v3[0];
  frame[5] = v3[1];
  frame[8] = v3[2];
}









void VT_ComputeRealFrameWithPivot( double *frame,
				   double *previousFrame, 
				   double tolerance )
{
  double z1[3];
  double z2[3];
  double nz;
  double r[3];
  double nr;
  double n;
  double sinus, cosinus;
  double ROT[9];

  if ( previousFrame == (double*)NULL ) {
    VT_ComputeRealFrame( frame, NULL );
    return;
  }


  z1[0] = previousFrame[2];
  z1[1] = previousFrame[5];
  z1[2] = previousFrame[8];
  nz = sqrt( z1[0] * z1[0] + z1[1] * z1[1] + z1[2] * z1[2] );
  z1[0] /= nz;
  z1[1] /= nz;
  z1[2] /= nz;

  z2[0] = frame[2];
  z2[1] = frame[5];
  z2[2] = frame[8];
  nz = sqrt( z2[0] * z2[0] + z2[1] * z2[1] + z2[2] * z2[2] );
  z2[0] /= nz;
  z2[1] /= nz;
  z2[2] /= nz;

  /* on verifie le sens des vecteurs
   */
  if ( z1[0]*z2[0] + z1[1]*z2[1] + z1[2]*z2[2] < 0.0 ) {
    z2[0] = -z2[0];
    z2[1] = -z2[1];
    z2[2] = -z2[2];
  }

  /* on a z_frame
   */
  frame[2] = z2[0];
  frame[5] = z2[1];
  frame[8] = z2[2];


  /* l'axe de rotation est r = z1 vectoriel z2
   */
  r[0] = z1[1] * z2[2] - z1[2] * z2[1];
  r[1] = z1[2] * z2[0] - z1[0] * z2[2];
  r[2] = z1[0] * z2[1] - z1[1] * z2[0];
  nr = sqrt( r[0] * r[0] + r[1] * r[1] + r[2] * r[2] );
  sinus = nr;

  
  /* plans paralleles ou presque,
     les deux bases doivent etre quasiment identiques
     par construction (id cas generique)
     y_frame = z_frame  vectoriel  x_previousframe
     x_frame = y_frame  vectoriel  z_frame
  */
  if ( nr < tolerance ) {

    frame[1] = z2[1] * previousFrame[6] - z2[2] * previousFrame[3];
    frame[4] = z2[2] * previousFrame[0] - z2[0] * previousFrame[6];
    frame[7] = z2[0] * previousFrame[3] - z2[1] * previousFrame[0];
    n = sqrt( frame[1]*frame[1] + frame[4]*frame[4] + frame[7]*frame[7] );
    frame[1] /= n;
    frame[4] /= n;
    frame[7] /= n;

    frame[0] = frame[4] * frame[8] - frame[7] * frame[5];
    frame[3] = frame[7] * frame[2] - frame[1] * frame[8];
    frame[6] = frame[1] * frame[5] - frame[4] * frame[2];
    n = sqrt( frame[0]*frame[0] + frame[3]*frame[3] + frame[6]*frame[6] );
    frame[0] /= n;
    frame[3] /= n;
    frame[6] /= n;

    return;
  }

  r[0] /= nr;
  r[1] /= nr;
  r[2] /= nr;
  cosinus = z1[0] * z2[0] + z1[1] * z2[1] + z1[2] * z2[2];
  
  /* on construit explicitement la matrice de rotation
   */
  ROT[0] = 1.0 - (1.0-cosinus) * (r[1]*r[1] + r[2]*r[2]);
  ROT[4] = 1.0 - (1.0-cosinus) * (r[0]*r[0] + r[2]*r[2]);
  ROT[8] = 1.0 - (1.0-cosinus) * (r[0]*r[0] + r[1]*r[1]);

  ROT[1] = - sinus * r[2] + (1.0-cosinus) * r[0] * r[1];
  ROT[2] =   sinus * r[1] + (1.0-cosinus) * r[0] * r[2];
  ROT[3] =   sinus * r[2] + (1.0-cosinus) * r[0] * r[1];
  ROT[5] = - sinus * r[0] + (1.0-cosinus) * r[1] * r[2];
  ROT[6] = - sinus * r[1] + (1.0-cosinus) * r[0] * r[2];
  ROT[7] =   sinus * r[0] + (1.0-cosinus) * r[1] * r[2];

  /* on l'applique
   */
  frame[1] = ROT[0] * previousFrame[1] + ROT[1] * previousFrame[4] + ROT[2] * previousFrame[7];
  frame[4] = ROT[3] * previousFrame[1] + ROT[4] * previousFrame[4] + ROT[5] * previousFrame[7];
  frame[7] = ROT[6] * previousFrame[1] + ROT[7] * previousFrame[4] + ROT[8] * previousFrame[7];
  n = sqrt( frame[1]*frame[1] + frame[4]*frame[4] + frame[7]*frame[7] );
  frame[1] /= n;
  frame[4] /= n;
  frame[7] /= n;

  frame[0] = ROT[0] * previousFrame[0] + ROT[1] * previousFrame[3] + ROT[2] * previousFrame[6];
  frame[3] = ROT[3] * previousFrame[0] + ROT[4] * previousFrame[3] + ROT[5] * previousFrame[6];
  frame[6] = ROT[6] * previousFrame[0] + ROT[7] * previousFrame[3] + ROT[8] * previousFrame[6];
  n = sqrt( frame[0]*frame[0] + frame[3]*frame[3] + frame[6]*frame[6] );
  frame[0] /= n;
  frame[3] /= n;
  frame[6] /= n;

  /*
  if ( 0 ) {
    double t[3];
    double p;
    t[0] = ROT[0] * z1[0] + ROT[1] * z1[1] + ROT[2] * z1[2];
    t[1] = ROT[3] * z1[0] + ROT[4] * z1[1] + ROT[5] * z1[2];
    t[2] = ROT[6] * z1[0] + ROT[7] * z1[1] + ROT[8] * z1[2];
    p = z2[0] * t[0] + z2[1] * t[1] + z2[2] * t[2];
    fprintf( stdout, "1 - %f = %g\n", p, 1.0 - p );
    fprintf( stdout, "\n" );
  }
  */
}









static void VT_ComputeNextVoxelFrame( double thetaInImage,
				      double phiInImage,
				      double *ptInImage,
				      double *ptSize,
				      double *sliceSize,
				      double *sliceCtr,
				      double *newRealFrame,
				      double *oldRealFrame, 
				      double *mat )
{
  double frame[9];


  /* calcul de la matrice de rotation dans le
     monde 'reel'
  */
  VT_ComputeNextRealFrame( thetaInImage, phiInImage, ptSize,
			   frame, oldRealFrame );

  if ( newRealFrame != (double*)NULL ) {
    (void)memcpy( newRealFrame, frame, 9*sizeof( double ) );
  }
  
  
  mat[ 0] = sliceSize[0] * frame[0];
  mat[ 1] = sliceSize[1] * frame[1];
  mat[ 2] = sliceSize[2] * frame[2];

  mat[ 4] = sliceSize[0] * frame[3];
  mat[ 5] = sliceSize[1] * frame[4];
  mat[ 6] = sliceSize[2] * frame[5];

  mat[ 8] = sliceSize[0] * frame[6];
  mat[ 9] = sliceSize[1] * frame[7];
  mat[10] = sliceSize[2] * frame[8];

  mat[12] = mat[13] = mat[14] = 0.0;
  mat[15] = 1;
  
  /* le centre de la coupe (z=0) doit donner le point de l'image
     on incorpore deja le dernier changement d'echelle
   */
  mat[ 3] = ptInImage[0] - (mat[0]*sliceCtr[0] + mat[1]*sliceCtr[1])/ptSize[0];
  mat[ 7] = ptInImage[1] - (mat[4]*sliceCtr[0] + mat[5]*sliceCtr[1])/ptSize[1];
  mat[11] = ptInImage[2] - (mat[8]*sliceCtr[0] + mat[9]*sliceCtr[1])/ptSize[2];
  
  
  /* dernier changement d'echelle
   */
  mat[ 0] /= ptSize[0];
  mat[ 1] /= ptSize[0];
  mat[ 2] /= ptSize[0];

  mat[ 4] /= ptSize[1];
  mat[ 5] /= ptSize[1];
  mat[ 6] /= ptSize[1];

  mat[ 8] /= ptSize[2];
  mat[ 9] /= ptSize[2];
  mat[10] /= ptSize[2];


}





int VT_ReechTriLinSlice( vt_image *theSlice,
			 int slice,
			 vt_image *theIm,
			 double *ptInImage,
			 double thetaInImage,
			 double phiInImage,
			 double *newRealFrame,
			 double *oldRealFrame,
			 double *newMat )
{
  char *proc ="VT_ReechTriLinSlice";
  double mat[16];
  double ptSize[3], sliceSize[3];
  double sliceCtr[2];
  int theDim[3], sliceDim[3];
  int offset;

  if ( theIm->type != theSlice->type ) {
    VT_Error( "images have different types", proc );
    return( -1 );
  }

  ptSize[0] = theIm->siz.x;
  ptSize[1] = theIm->siz.y;
  ptSize[2] = theIm->siz.z;

  sliceSize[0] = theSlice->siz.x;
  sliceSize[1] = theSlice->siz.y;
  sliceSize[2] = theSlice->siz.z;

  sliceCtr[0] = (double)(theSlice->dim.x - 1)/2.0;
  sliceCtr[1] = (double)(theSlice->dim.y - 1)/2.0;


  VT_ComputeNextVoxelFrame( thetaInImage, phiInImage, ptInImage,
			    ptSize,
			    sliceSize, sliceCtr, 
			    newRealFrame, oldRealFrame, mat );



  


  if ( newMat != (double*)NULL ) {
    (void)memcpy( newMat, mat, 16*sizeof( double ) );
  }


  /* la matrice permet maintenant de passer de theSlice
     a theIm (pt(im) = mat * pt(slice))
  */

  theDim[0] = theIm->dim.x;
  theDim[1] = theIm->dim.y;
  theDim[2] = theIm->dim.z;

  sliceDim[0] = theSlice->dim.x;
  sliceDim[1] = theSlice->dim.y;
  sliceDim[2] = 1;

  offset = slice * theSlice->dim.x * theSlice->dim.y;

  switch( theIm->type ) {
  case UCHAR :
    Reech3DTriLin4x4_u8( theIm->buf, theDim, 
			 (void*)&((u8*)(theSlice->buf))[offset], 
			 sliceDim, mat );
    break;
  case SCHAR :
    Reech3DTriLin4x4_s8( theIm->buf, theDim, 
			 (void*)&((s8*)(theSlice->buf))[offset], 
			 sliceDim, mat );
    break;
  case USHORT :
    Reech3DTriLin4x4_u16( theIm->buf, theDim, 
			  (void*)&((u16*)(theSlice->buf))[offset], 
			  sliceDim, mat );
    break;
  case SSHORT :
    Reech3DTriLin4x4_s16( theIm->buf, theDim, 
			  (void*)&((s16*)(theSlice->buf))[offset], 
			  sliceDim, mat );
    break;
  default :
    VT_Error( "unable to deal with such image type", proc );
    return( -1 );
  }
  return( 1 );
}





int VT_ReechTriLinSlice2( vt_image *theSlice,
			 int slice,
			 vt_image *theIm,
			 double *mat )
{
  char *proc ="VT_ReechTriLinSlice2";
  int theDim[3], sliceDim[3];
  int offset;

  if ( theIm->type != theSlice->type ) {
    VT_Error( "images have different types", proc );
    return( -1 );
  }

  /* la matrice permet maintenant de passer de theSlice
     a theIm (pt(im) = mat * pt(slice))
  */

  theDim[0] = theIm->dim.x;
  theDim[1] = theIm->dim.y;
  theDim[2] = theIm->dim.z;

  sliceDim[0] = theSlice->dim.x;
  sliceDim[1] = theSlice->dim.y;
  sliceDim[2] = 1;

  offset = slice * theSlice->dim.x * theSlice->dim.y;

  switch( theIm->type ) {
  case UCHAR :
    Reech3DTriLin4x4_u8( theIm->buf, theDim, 
			 (void*)&((u8*)(theSlice->buf))[offset], 
			 sliceDim, mat );
    break;
  case SCHAR :
    Reech3DTriLin4x4_s8( theIm->buf, theDim, 
			 (void*)&((s8*)(theSlice->buf))[offset], 
			 sliceDim, mat );
    break;
  case USHORT :
    Reech3DTriLin4x4_u16( theIm->buf, theDim, 
			  (void*)&((u16*)(theSlice->buf))[offset], 
			  sliceDim, mat );
    break;
  case SSHORT :
    Reech3DTriLin4x4_s16( theIm->buf, theDim, 
			  (void*)&((s16*)(theSlice->buf))[offset], 
			  sliceDim, mat );
    break;
  default :
    VT_Error( "unable to deal with such image type", proc );
    return( -1 );
  }
  return( 1 );
}





























int VT_AllocCSplinesInSlice( typeCSplinesInSlice *par,
			     vt_image *theIm )
{
  char *proc = "VT_AllocCSplinesInSlice";
  char *generic  = "cspline";
  char name[256];
  
  sprintf( name, "%s.x.inr", generic );
  VT_InitImage( &(par->imx), name, theIm->dim.x, theIm->dim.y, 1, FLOAT );
  sprintf( name, "%s.y.inr", generic );
  VT_InitImage( &(par->imy), name, theIm->dim.x, theIm->dim.y, 1, FLOAT );
  sprintf( name, "%s.z.inr", generic );
  VT_InitImage( &(par->imz), name, theIm->dim.x, theIm->dim.y, 1, FLOAT );

  sprintf( name, "%s.xx.inr", generic );
  VT_InitImage( &(par->imxx), name, theIm->dim.x, theIm->dim.y, 1, FLOAT );
  sprintf( name, "%s.xy.inr", generic );
  VT_InitImage( &(par->imxy), name, theIm->dim.x, theIm->dim.y, 1, FLOAT );
  sprintf( name, "%s.xz.inr", generic );
  VT_InitImage( &(par->imxz), name, theIm->dim.x, theIm->dim.y, 1, FLOAT );

  sprintf( name, "%s.yy.inr", generic );
  VT_InitImage( &(par->imyy), name, theIm->dim.x, theIm->dim.y, 1, FLOAT );
  sprintf( name, "%s.yz.inr", generic );
  VT_InitImage( &(par->imyz), name, theIm->dim.x, theIm->dim.y, 1, FLOAT );

  sprintf( name, "%s.zz.inr", generic );
  VT_InitImage( &(par->imzz), name, theIm->dim.x, theIm->dim.y, 1, FLOAT );

  sprintf( name, "%s.d1.inr", generic );
  VT_InitImage( &(par->d1), name, theIm->dim.x, theIm->dim.y, 1, FLOAT );

  sprintf( name, "%s.d2.inr", generic );
  VT_InitImage( &(par->d2), name, theIm->dim.x, theIm->dim.y, 1, FLOAT );



  if ( VT_AllocImage( &(par->imx ) ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate X slice\n", proc );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imy ) ) != 1 ) {
    VT_FreeImage( &(par->imx ) );
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate Y slice\n", proc );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imz ) ) != 1 ) {
    VT_FreeImage( &(par->imx ) );
    VT_FreeImage( &(par->imy ) );
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate Z slice\n", proc );
    return( -1 );
  }



  if ( VT_AllocImage( &(par->imxx ) ) != 1 ) {
    VT_FreeImage( &(par->imx ) );
    VT_FreeImage( &(par->imy ) );
    VT_FreeImage( &(par->imz ) );
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate XX slice\n", proc );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imxy ) ) != 1 ) {
    VT_FreeImage( &(par->imx ) );
    VT_FreeImage( &(par->imy ) );
    VT_FreeImage( &(par->imz ) );
    VT_FreeImage( &(par->imxx ) );
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate XY slice\n", proc );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imxz ) ) != 1 ) {
    VT_FreeImage( &(par->imx ) );
    VT_FreeImage( &(par->imy ) );
    VT_FreeImage( &(par->imz ) );
    VT_FreeImage( &(par->imxx ) );
    VT_FreeImage( &(par->imxy ) );
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate XZ slice\n", proc );
    return( -1 );
  }



  if ( VT_AllocImage( &(par->imyy ) ) != 1 ) {
    VT_FreeImage( &(par->imx ) );
    VT_FreeImage( &(par->imy ) );
    VT_FreeImage( &(par->imz ) );
    VT_FreeImage( &(par->imxx ) );
    VT_FreeImage( &(par->imxy ) );
    VT_FreeImage( &(par->imxz ) );
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate YY slice\n", proc );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imyz ) ) != 1 ) {
    VT_FreeImage( &(par->imx ) );
    VT_FreeImage( &(par->imy ) );
    VT_FreeImage( &(par->imz ) );
    VT_FreeImage( &(par->imxx ) );
    VT_FreeImage( &(par->imxy ) );
    VT_FreeImage( &(par->imxz ) );
    VT_FreeImage( &(par->imyy ) );
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate YZ slice\n", proc );
    return( -1 );
  }



  if ( VT_AllocImage( &(par->imzz ) ) != 1 ) {
    VT_FreeImage( &(par->imx ) );
    VT_FreeImage( &(par->imy ) );
    VT_FreeImage( &(par->imz ) );
    VT_FreeImage( &(par->imxx ) );
    VT_FreeImage( &(par->imxy ) );
    VT_FreeImage( &(par->imxz ) );
    VT_FreeImage( &(par->imyy ) );
    VT_FreeImage( &(par->imyz ) );
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate ZZ slice\n", proc );
    return( -1 );
  }


  if ( VT_AllocImage( &(par->d1 ) ) != 1 ) {
    VT_FreeImage( &(par->imx ) );
    VT_FreeImage( &(par->imy ) );
    VT_FreeImage( &(par->imz ) );
    VT_FreeImage( &(par->imxx ) );
    VT_FreeImage( &(par->imxy ) );
    VT_FreeImage( &(par->imxz ) );
    VT_FreeImage( &(par->imyy ) );
    VT_FreeImage( &(par->imyz ) );
    VT_FreeImage( &(par->imzz ) );
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate D1 slice\n", proc );
    return( -1 );
  }


  if ( VT_AllocImage( &(par->d2 ) ) != 1 ) {
    VT_FreeImage( &(par->imx ) );
    VT_FreeImage( &(par->imy ) );
    VT_FreeImage( &(par->imz ) );
    VT_FreeImage( &(par->imxx ) );
    VT_FreeImage( &(par->imxy ) );
    VT_FreeImage( &(par->imxz ) );
    VT_FreeImage( &(par->imyy ) );
    VT_FreeImage( &(par->imyz ) );
    VT_FreeImage( &(par->imzz ) );
    VT_FreeImage( &(par->d1 ) );
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate D2 slice\n", proc );
    return( -1 );
  }


  return( 1 );
}



void VT_FreeCSplinesInSlice( typeCSplinesInSlice *par )
{
  VT_FreeImage( &(par->imx ) );
  VT_FreeImage( &(par->imy ) );
  VT_FreeImage( &(par->imz ) );

  VT_FreeImage( &(par->imxx ) );
  VT_FreeImage( &(par->imxy ) );
  VT_FreeImage( &(par->imxz ) );
  VT_FreeImage( &(par->imyy ) );
  VT_FreeImage( &(par->imyz ) );
  VT_FreeImage( &(par->imzz ) );

  VT_FreeImage( &(par->d1 ) );
  VT_FreeImage( &(par->d2 ) );
}










int VT_ReechCSplineSlice( vt_image *theSlice,
			  vt_image *theGrad,
			  vt_image *theLapl,
			  int slice,
			  typeCSplineCoefficients *theCoeff,
			  typeCSplinesInSlice *aux,
			  vt_image *theIm,
			  double *ptInImage,
			  double thetaInImage,
			  double phiInImage,
			  double *newRealFrame,
			  double *oldRealFrame,
			  double *newMat )
{
  char *proc ="VT_ReechCSplineSlice";
  double mat[16];
  double ptSize[3], sliceSize[3];
  double sliceCtr[2];
  int sliceDim[3];
  int offset;
  int deriv[3] = {0,0,0};

  int i;
  r32 *theX = (r32*)aux->imx.buf;
  r32 *theY = (r32*)aux->imy.buf;
  r32 *theZ = (r32*)aux->imz.buf;
  r32 *theXX = (r32*)aux->imxx.buf;
  r32 *theXY = (r32*)aux->imxy.buf;
  r32 *theXZ = (r32*)aux->imxz.buf;
  r32 *theYY = (r32*)aux->imyy.buf;
  r32 *theYZ = (r32*)aux->imyz.buf;
  r32 *theZZ = (r32*)aux->imzz.buf;
  r32 *theD1 = (r32*)aux->d1.buf;
  r32 *theD2 = (r32*)aux->d2.buf;
 
  double dx, dy, v[3], dxx, dxy, dyy;

  ptSize[0] = theIm->siz.x;
  ptSize[1] = theIm->siz.y;
  ptSize[2] = theIm->siz.z;

  sliceSize[0] = theSlice->siz.x;
  sliceSize[1] = theSlice->siz.y;
  sliceSize[2] = theSlice->siz.z;
  
  sliceCtr[0] = (double)(theSlice->dim.x - 1)/2.0;
  sliceCtr[1] = (double)(theSlice->dim.y - 1)/2.0;
  

  VT_ComputeNextVoxelFrame( thetaInImage, phiInImage, ptInImage,
			    ptSize,
			    sliceSize, sliceCtr, 
			    newRealFrame, oldRealFrame, mat );

  

  


  if ( newMat != (double*)NULL ) {
    (void)memcpy( newMat, mat, 16*sizeof( double ) );
  }
  

  /* la matrice permet maintenant de passer de theSlice
     a theIm (pt(im) = mat * pt(slice))
  */

  /*theDim[0] = theIm->dim.x;
  theDim[1] = theIm->dim.y;
  theDim[2] = theIm->dim.z;
    */
  sliceDim[0] = theSlice->dim.x;
  sliceDim[1] = theSlice->dim.y;
  sliceDim[2] = 1;



  offset = slice * theSlice->dim.x * theSlice->dim.y;


  if ( Reech3DCSpline4x4WithCoefficients( theCoeff, (r32*)(aux->d1.buf), 
					  sliceDim, mat, 0, deriv ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to compute (0,0,0) slice\n", proc );
    return( -1 );
  }
  





  switch( theSlice->type ) {
  default :
    if ( _verbose_ ) 
      fprintf( stderr, "%s: such output type not handled in switch (1)\n", proc );
    return( -1 );
  case UCHAR :
    if ( ConvertBuffer( aux->d1.buf, FLOAT, (void*)&((u8*)(theSlice->buf))[offset],
			theSlice->type, theSlice->dim.x*theSlice->dim.y ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to convert buffer\n", proc );
      return( -1 );
    }
    break;
  case USHORT :
    if ( ConvertBuffer( aux->d1.buf, FLOAT, (void*)&((u16*)(theSlice->buf))[offset],
		   theSlice->type, theSlice->dim.x*theSlice->dim.y ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to convert buffer\n", proc );
      return( -1 );
    }
    break;
  case SSHORT :
    if ( ConvertBuffer( aux->d1.buf, FLOAT, (void*)&((s16*)(theSlice->buf))[offset],
		   theSlice->type, theSlice->dim.x*theSlice->dim.y ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to convert buffer\n", proc );
      return( -1 );
    }
    break;
  case FLOAT :
    if ( ConvertBuffer( aux->d1.buf, FLOAT, (void*)&((r32*)(theSlice->buf))[offset],
		   theSlice->type, theSlice->dim.x*theSlice->dim.y ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to convert buffer\n", proc );
      return( -1 );
    }
    break;
  }






  if ( theGrad != (vt_image*)NULL && theLapl == (vt_image*)NULL ) {

    deriv[0] = 1;   deriv[1] = 0;   deriv[2] = 0;
    if ( Reech3DCSpline4x4WithCoefficients( theCoeff, (r32*)(aux->imx.buf), 
					    sliceDim, mat, 0, deriv ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to compute (1,0,0) slice\n", proc );
      return( -1 );
    }

    deriv[0] = 0;   deriv[1] = 1;   deriv[2] = 0;
    if ( Reech3DCSpline4x4WithCoefficients( theCoeff, (r32*)(aux->imy.buf), 
					    sliceDim, mat, 0, deriv ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to compute (0,1,0) slice\n", proc );
      return( -1 );
    }

    deriv[0] = 0;   deriv[1] = 0;   deriv[2] = 1;
    if ( Reech3DCSpline4x4WithCoefficients( theCoeff, (r32*)(aux->imz.buf), 
					    sliceDim, mat, 0, deriv ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to compute (0,0,1) slice\n", proc );
      return( -1 );
    }
    

    /* calcul du gradient 2D dans la coupe
    */
    for ( i=0; i<(int)(theSlice->dim.x*theSlice->dim.y); i++ ) {
      dx = theX[i]*newRealFrame[ 0] + 
	theY[i]*newRealFrame[ 3] + theZ[i]*newRealFrame[ 6];
      dy = theX[i]*newRealFrame[ 1] + 
	theY[i]*newRealFrame[ 4] + theZ[i]*newRealFrame[ 7];
      theD1[i] = sqrt( dx*dx + dy*dy );
    }


    switch( theGrad->type ) {
    default :
      if ( _verbose_ ) 
	fprintf( stderr, "%s: such output type not handled in switch (1)\n", proc );
      return( -1 );
    case UCHAR :
      if ( ConvertBuffer( aux->d1.buf, FLOAT, (void*)&((u8*)(theGrad->buf))[offset],
			  theGrad->type, theGrad->dim.x*theGrad->dim.y ) != 1 ) {
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: unable to convert buffer\n", proc );
	return( -1 );
      }
      break;
    case USHORT :
      if ( ConvertBuffer( aux->d1.buf, FLOAT, (void*)&((u16*)(theGrad->buf))[offset],
		     theGrad->type, theGrad->dim.x*theGrad->dim.y ) != 1 ) {
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: unable to convert buffer\n", proc );
	return( -1 );
      }
      break;
    case FLOAT :
      if ( ConvertBuffer( aux->d1.buf, FLOAT, (void*)&((r32*)(theGrad->buf))[offset],
		     theGrad->type, theGrad->dim.x*theGrad->dim.y ) != 1 ) {
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: unable to convert buffer\n", proc );
	return( -1 );
      }
      break;
    }

  } /* if ( theGrad != (vt_image*)NULL ) */




  if ( theGrad != (vt_image*)NULL && theLapl != (vt_image*)NULL ) {

    deriv[0] = 1;   deriv[1] = 0;   deriv[2] = 0;
    if ( Reech3DCSpline4x4WithCoefficients( theCoeff, (r32*)(aux->imx.buf), 
					    sliceDim, mat, 0, deriv ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to compute (1,0,0) slice\n", proc );
      return( -1 );
    }

    deriv[0] = 0;   deriv[1] = 1;   deriv[2] = 0;
    if ( Reech3DCSpline4x4WithCoefficients( theCoeff, (r32*)(aux->imy.buf), 
					    sliceDim, mat, 0, deriv ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to compute (0,1,0) slice\n", proc );
      return( -1 );
    }

    deriv[0] = 0;   deriv[1] = 0;   deriv[2] = 1;
    if ( Reech3DCSpline4x4WithCoefficients( theCoeff, (r32*)(aux->imz.buf), 
					    sliceDim, mat, 0, deriv ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to compute (0,0,1) slice\n", proc );
      return( -1 );
    }

    deriv[0] = 2;   deriv[1] = 0;   deriv[2] = 0;
    if ( Reech3DCSpline4x4WithCoefficients( theCoeff, (r32*)(aux->imxx.buf), 
					    sliceDim, mat, 0, deriv ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to compute (2,0,0) slice\n", proc );
      return( -1 );
    }

    deriv[0] = 1;   deriv[1] = 1;   deriv[2] = 0;
    if ( Reech3DCSpline4x4WithCoefficients( theCoeff, (r32*)(aux->imxy.buf), 
					    sliceDim, mat, 0, deriv ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to compute (1,1,0) slice\n", proc );
      return( -1 );
    }

    deriv[0] = 1;   deriv[1] = 0;   deriv[2] = 1;
    if ( Reech3DCSpline4x4WithCoefficients( theCoeff, (r32*)(aux->imxz.buf), 
					    sliceDim, mat, 0, deriv ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to compute (1,0,1) slice\n", proc );
      return( -1 );
    }

    deriv[0] = 0;   deriv[1] = 2;   deriv[2] = 0;
    if ( Reech3DCSpline4x4WithCoefficients( theCoeff, (r32*)(aux->imyy.buf), 
					    sliceDim, mat, 0, deriv ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to compute (0,2,0) slice\n", proc );
      return( -1 );
    }

    deriv[0] = 0;   deriv[1] = 1;   deriv[2] = 1;
    if ( Reech3DCSpline4x4WithCoefficients( theCoeff, (r32*)(aux->imyz.buf), 
					    sliceDim, mat, 0, deriv ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to compute (0,1,1) slice\n", proc );
      return( -1 );
    }

    deriv[0] = 0;   deriv[1] = 0;   deriv[2] = 2;
    if ( Reech3DCSpline4x4WithCoefficients( theCoeff, (r32*)(aux->imzz.buf), 
					    sliceDim, mat, 0, deriv ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to compute (0,0,2) slice\n", proc );
      return( -1 );
    }
    
    /* calcul du gradient x Hessien . gradient
    */
    for ( i=0; i<(int)(theSlice->dim.x*theSlice->dim.y); i++ ) {

      dx  = theX[i]*newRealFrame[ 0] + 
	    theY[i]*newRealFrame[ 3] + theZ[i]*newRealFrame[ 6];
      dy  = theX[i]*newRealFrame[ 1] + 
	    theY[i]*newRealFrame[ 4] + theZ[i]*newRealFrame[ 7];

      theD1[i] =  dx*dx + dy*dy ;

      v[0] = theXX[i]*newRealFrame[ 0] +
	     theXY[i]*newRealFrame[ 3] + theXZ[i]*newRealFrame[ 6]; 
      v[1] = theXY[i]*newRealFrame[ 0] +
	     theYY[i]*newRealFrame[ 3] + theYZ[i]*newRealFrame[ 6]; 
      v[2] = theXZ[i]*newRealFrame[ 0] +
	     theYZ[i]*newRealFrame[ 3] + theZZ[i]*newRealFrame[ 6]; 

      dxx = v[0]*newRealFrame[ 0] + v[1]*newRealFrame[ 3] + 
	    v[2]*newRealFrame[ 6];
      dxy = v[0]*newRealFrame[ 1] + v[1]*newRealFrame[ 4] + 
	    v[2]*newRealFrame[ 7];

      v[0] = theXX[i]*newRealFrame[ 1] +
	     theXY[i]*newRealFrame[ 4] + theXZ[i]*newRealFrame[ 7]; 
      v[1] = theXY[i]*newRealFrame[ 1] +
	     theYY[i]*newRealFrame[ 4] + theYZ[i]*newRealFrame[ 7]; 
      v[2] = theXZ[i]*newRealFrame[ 1] +
	     theYZ[i]*newRealFrame[ 4] + theZZ[i]*newRealFrame[ 7]; 

      dyy = v[0]*newRealFrame[ 1] + v[1]*newRealFrame[ 4] + 
	    v[2]*newRealFrame[ 7];
      
      theD2[i] = dx * (dxx * dx + dxy * dy) +
	         dy * (dxy * dx + dyy * dy);

      if ( theD1[i] > 1e-10 ) theD2[i] /= theD1[i];
      theD1[i] = sqrt( theD1[i] );
    }

    switch( theGrad->type ) {
    default :
      if ( _verbose_ ) 
	fprintf( stderr, "%s: such output type not handled in switch (1)\n", proc );
      return( -1 );
    case UCHAR :
      if ( ConvertBuffer( aux->d1.buf, FLOAT, (void*)&((u8*)(theGrad->buf))[offset],
		     theGrad->type, theGrad->dim.x*theGrad->dim.y ) != 1 ) {
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: unable to convert buffer\n", proc );
	return( -1 );
      }
      break;
    case USHORT :
      if ( ConvertBuffer( aux->d1.buf, FLOAT, (void*)&((u16*)(theGrad->buf))[offset],
		     theGrad->type, theGrad->dim.x*theGrad->dim.y ) != 1 ) {
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: unable to convert buffer\n", proc );
	return( -1 );
      }
      break;
    case FLOAT :
      if ( ConvertBuffer( aux->d1.buf, FLOAT, (void*)&((r32*)(theGrad->buf))[offset],
		     theGrad->type, theGrad->dim.x*theGrad->dim.y ) != 1 ) {
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: unable to convert buffer\n", proc );
	return( -1 );
      }
      break;
    }

    switch( theLapl->type ) {
    default :
      if ( _verbose_ ) 
	fprintf( stderr, "%s: such output type not handled in switch (1)\n", proc );
      return( -1 );
    case UCHAR :
      if ( ConvertBuffer( aux->d2.buf, FLOAT, (void*)&((u8*)(theLapl->buf))[offset],
		     theLapl->type, theLapl->dim.x*theLapl->dim.y ) != 1 ) {
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: unable to convert buffer\n", proc );
	return( -1 );
      }
      break;
    case USHORT :
      if ( ConvertBuffer( aux->d2.buf, FLOAT, (void*)&((u16*)(theLapl->buf))[offset],
		     theLapl->type, theLapl->dim.x*theLapl->dim.y ) != 1 ) {
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: unable to convert buffer\n", proc );
	return( -1 );
      }
      break;
    case FLOAT :
      if ( ConvertBuffer( aux->d2.buf, FLOAT, (void*)&((r32*)(theLapl->buf))[offset],
		     theLapl->type, theLapl->dim.x*theLapl->dim.y ) != 1 ) {
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: unable to convert buffer\n", proc );
	return( -1 );
      }
      break;
    }

  }


  return( 1 );
}





int VT_ReechCSplineSlice2( vt_image *theSlice,
			  int slice,
			  typeCSplineCoefficients *theCoeff,
			  typeCSplinesInSlice *aux,
			  double *mat )
{
  char *proc ="VT_ReechCSplineSlice2";
  int sliceDim[3];
  int offset;
  int deriv[3] = {0,0,0};


  sliceDim[0] = theSlice->dim.x;
  sliceDim[1] = theSlice->dim.y;
  sliceDim[2] = 1;

  offset = slice * theSlice->dim.x * theSlice->dim.y;


  if ( Reech3DCSpline4x4WithCoefficients( theCoeff, (r32*)(aux->d1.buf), 
					  sliceDim, mat, 0, deriv ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to compute (0,0,0) slice\n", proc );
    return( -1 );
  }


  switch( theSlice->type ) {
  default :
    if ( _verbose_ ) 
      fprintf( stderr, "%s: such output type not handled in switch (1)\n", proc );
    return( -1 );
  case UCHAR :
    if ( ConvertBuffer( aux->d1.buf, FLOAT, (void*)&((u8*)(theSlice->buf))[offset],
			theSlice->type, theSlice->dim.x*theSlice->dim.y ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to convert buffer\n", proc );
      return( -1 );
    }
    break;
  case USHORT :
    if ( ConvertBuffer( aux->d1.buf, FLOAT, (void*)&((u16*)(theSlice->buf))[offset],
		   theSlice->type, theSlice->dim.x*theSlice->dim.y ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to convert buffer\n", proc );
      return( -1 );
    }
    break;
  case SSHORT :
    if ( ConvertBuffer( aux->d1.buf, FLOAT, (void*)&((s16*)(theSlice->buf))[offset],
		   theSlice->type, theSlice->dim.x*theSlice->dim.y ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to convert buffer\n", proc );
      return( -1 );
    }
    break;
  case FLOAT :
    if ( ConvertBuffer( aux->d1.buf, FLOAT, (void*)&((r32*)(theSlice->buf))[offset],
		   theSlice->type, theSlice->dim.x*theSlice->dim.y ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to convert buffer\n", proc );
      return( -1 );
    }
    break;
  }

  return( 1 );
}

