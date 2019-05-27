/*************************************************************************
 * tube3D.h -
 *
 * $Id: tube3D.h,v 1.5 2001/06/08 15:10:52 greg Exp $
 *
 * Copyright (c) INRIA 2000
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Mon Sep 11 18:21:59 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _tube3D_h_
#define _tube3D_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
  /* #include <malloc.h> */
#include <math.h>
#include <string.h>


#include <typedefs.h>
#include <recbuffer.h>
#include <eigens.h>
#include <issimple3D.h>



/* calcul d'une reponse multi-echelles a partir
   des reponses pre-calculees a chaque echelle

    RETURN:
   -1 en cas d'erreur
    1 si tout va bien

*/
extern int combine_3D_line_response_at_several_scales( float *resResponse,
						float *resTheta,
						float *resPhi,
						float **bufResponse,
						float **bufTheta,
						float **bufPhi,
						int *theDim,
						int nbScales );






/* calcul de la reponse a plusieurs echelles du
   detecteur de structures tubulaires 3D

    RETURN:
   -1 en cas d'erreur
    1 si tout va bien
  
*/
extern int compute_3D_line_response_at_several_scales( void *theBuf,
						bufferType typeBuf,
						float **bufResponse,
						float **bufTheta,
						float **bufPhi,
						int *theDim,
						double *scale,
						int nbScales );





/* calcul incremental de la reponse multi-echelle du
   detecteur de structures tubulaires 3D.

   DESCRIPTION:
   Pour la premiere echelle (ie scale[0])
   le resultat est directement calcule dans les images
   resultat. 
   Pour les autres echelles (ie scale[i], i de 1 a nbScales-1)
   le resultat est calcule dans des images intermediaires
   et le resultat global est mis a jour.

   Ainsi, le buffer *bufResponse contiendra la reponse multi-echelle
   maximale, et les buffers *bufTheta et *bufPhi la direction
   presumee de la structure tubulaire calculee a l'echelle 
   correspondante.

   RETURN:
   -1 en cas d'erreur
      toutefois, il se peut que le calcul de 
      la premiere echelle ait eu lieu, et 
      que le resultat soit range dans les images 
      resultat
    1 si tout va bien

*/
extern int compute_3D_line_response_at_multiple_scales( void *theBuf,
						 bufferType typeBuf,
						 float *bufResponse,
						 float *bufTheta,
						 float *bufPhi,
						 int *theDim,
						 double *scale,
						 int nbScales );









/* calcul de la reponse a une seule echelle du
   detecteur de structures tubulaires 3D

   DESCRIPTION:
   Le principe est de calculer en tout point le hessien
   (la matrice des derivees secondes) de l'image.
   L'etude des valeurs et vecteurs propres de cette
   matrice nous renseigne sur la structure locale 
   de l'image. En particulier, les structures tubulaires
   blanches sur fond noir sont caracterisees par
   2 valeurs propres negatives plus grandes en valeur absolue
   que la troisieme.
   On verifie donc que les 2 plus grandes valeurs propres
   en valeur absolue sont negatives, et que la plus petite
   valeur propre en valeur absolue est "beaucoup" (un facteur 2)
   plus petite que la valeur mediane.

   Le vecteur propre associe a cette troisieme valeur absolue 
   (la plus petite en valeur absolue donc) donne la direction
   presumee de la structure tubulaire observee.

   La reponse est la moyenne, sur un cercle autour de la 
   direction presumee de la structure tubulaire, d'une reponse
   a un filtre detecteur de contour, c'est-a-dire au produit 
   scalaire du gradient et de la direction radiale.

   On garde evidemment la valeur de la reponse (multipliee 
   par l'echelle pour repondre aux besoins multi-echelles)
   ainsi que la direction presumee de la structure tubulaire
   grace a deux angles (theta, phi). Le vecteur directeur etant
   alors :
   ( x )   ( cos theta sin phi )
   ( y ) = ( sin theta sin phi )
   ( z )   ( cos phi )
   

   PARAMETRES:
   void *theBuf : buffer de l'image 3D d'entree
   bufferType typeBuf : type de l'image d'entree
      defini dans typedefs.h
   float *bufResponse : buffer de la reponse 
      du detecteur. Doit deja etre alloue.
   float *bufTheta, float *bufPhi : buffers 3D
      donnant la direction de la structure 
      tubulaire. Doivent deja etre alloues.
   int *theDim : tableau de 3 entiers donnant les
      dimensions selon X, Y et Z des buffers 3D.
      Chacun des buffers doit donc etre de taille
      theDim[0]*theDim[1]*theDim[2]
   double scale : echelle de detection.

   RETURN:
   -1 en cas d'erreur
    1 si tout va bien
*/

extern int compute_3D_line_response_at_one_scale( void *theBuf,
					   bufferType typeBuf,
					   float *bufResponse,
					   float *bufTheta,
					   float *bufPhi,
					   int *theDim,
					   double scale );






/* calcul des extrema de la reponse d'un detecteur 3D de 
   structures tubulaires.
   
   DESCRIPTION:
   Pour chaque point de l'image des reponses, on calcule 
   les 2 vecteurs (supposes) orthogonaux grace au deux angles 
   (theta, phi) v1 et v2
   (  sin theta )         ( cos theta cos phi )
   ( -cos theta )   et    ( sin theta cos phi )
   ( 0 )                  ( -sin phi )
   
   Un point est extremal si sa reponse est plus grande (strictement)
   a la reponse interpolee lineairement aux 8 points suivants :
   M +/- vi (i=1,2) , M +/- a1 v1 +/- a2 v2 (ai = sqrt(2)/2)

   Les points non extremaux sont mis a 0.
   

   RETURN:
   -1 en cas d'erreur
    1 si tout va bien
*/

extern int compute_3D_line_response_extrema( float *bufExtrema, 
				      float *bufResponse,
				      float *bufTheta,
				      float *bufPhi,
				      int *theDim );




/* aminicissement adapte aux structures filaires.

   DESCRIPTION:
   Cette procedure efface iterativement les points simples.

   Pour chaque iteration, on examine chacune des 6
   directions pour laquelle on recherche des epaisseurs 
   (selon la 18-connexite) :
   - direction #1 : Z decroissant 
     on veut [x][y][z+1] appartenant au fond
     et [x+i][y+j][z-1] appartenant a l'objet 
     (pour au moins un couple i,j dans {-1,0,1})
   - direction #2 : Z croissant
   - direction #3 : Y decroissant 
   - direction #4 : Y croissant 
   - direction #5 : X decroissant 
   - direction #6 : X croissant 

   Pour chaque direction :
   - on commence par marquer les points effacables
     * points simples dans l'image
     * ET points simples dans l'image sans les points effacables
     Cela permet de ne pas effacer trop de points
   - certains de ces points simples sont marques comme points
     de fin (deviennent ineffacables) ce sont ceux qui n'ont qu'un
     seul voisin dans l'image (en comptant les points effacables)
   - lorsque tous les points sont parcourus, on efface tous les
     point effacables

   On arrete lorsqu'aucun point n'a ete efface pendant une iteration
   (c-a-d l'examen des 6 directions).


   RETURN:
   -1 en cas d'erreur
    1 si tout va bien
*/

extern int thin_3D_thick_lines( unsigned char *input_buffer,
				unsigned char *output_buffer,
				int *theDim );



extern void _VerboseInTube3D();
extern void _NoVerboseInTube3D();


#ifdef __cplusplus
}
#endif

#endif /* _tube3D_h_ */

