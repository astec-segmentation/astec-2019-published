/*************************************************************************
 * vt_tube3Delf.c -
 *
 * $Id: vt_tube3Delf.c,v 1.5 2000/10/05 16:01:36 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Tue Sep 12 14:26:12 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <vt_tube3Delf.h>

static int _verbose_ = 1;


int VT_3DResponseAtOneScale( vt_image *theIm,
			     vt_3Dimages *ims, 
			     vt_2Dimauxs *aux, 
			     float theSigma,
                             float theta __attribute__ ((unused)) )
{
  char *proc = "VT_3DResponseAtOneScale";
  float theCoeffs[3] = { theSigma, theSigma, theSigma }; 
  int slice;

  /* int    nbCirclePoints; */
  /* typeCirclePoint *thePts = (typeCirclePoint*)NULL; */


  /* pour calculer la reponse,
     on se place sur un cercle de rayon radius=theta*sigma
     que l'on dicretise en N=(int)(2*pi*radius) + 1
     pour obtenir environ 1 point par unite curviligne
     soit v2 et v3 les "vecteurs" du cercle
     les points se calculent alors
     par cos(i*2*pi/N) v2 + sin(i*2*pi/N) v3
     pour i = 0 ... N-1
  */
  
  /* thePts = VT_BuildCircle( theta*theSigma, &nbCirclePoints ); */

  




  /* le filtrage 3D pour les images :
     le gradient,
     les derivees partielles premiere et seconde en Z
   */
  if ( VT_Filter3Dimages( theIm, ims, theCoeffs ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to process 3D filters\n", proc );
    return( -1 );
  }
  






  for ( slice = 0; slice < (int)theIm->dim.z; slice ++ ) {
    

    /* le filtrage 2D pour les coupes
     */
    if ( VT_Filter2Dimauxs( ims, aux, theCoeffs, slice ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to process 2D filters in slice #%d\n", 
		 proc, slice );
      return( -1 );
    }
    




  }




  return( 1 );
}

