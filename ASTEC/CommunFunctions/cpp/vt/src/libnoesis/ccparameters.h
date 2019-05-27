/*************************************************************************
 * ccparameters.h - extraction de parametres sur des parties numerotees
 *
 * $Id: ccparameters.h,v 1.2 2001/04/03 10:27:21 greg Exp $
 *
 * DESCRIPTION: 
 *
 *
 *
 *
 *
 * AUTHOR:
 * Gregoire Malandain
 *
 * CREATION DATE:
 * Sun Feb 11 10:37:45 MET 2001
 *
 * Copyright Gregoire Malandain, INRIA
 *
 *
 * ADDITIONS, CHANGES:
 *
 *
 */

#ifndef _ccparameters_h_
#define _ccparameters_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <typedefs.h>

typedef struct {
  int volume;
  double maxDiameter;
  double medDiameter;
  double minDiameter;
} typeParameter;


extern void _VerboseInCcParameters();
extern void _NoVerboseInCcParameters();



/* calcul de parametres sur des parties etiquettees
   
   on calcule les parametres sur les parties numerotees
   de 1 a nbLabels (inclus). On renvoie un tableau de 
   structures de type typeParameter de taille (nbLabels+1).
   on adresse donc la structure correspondant a chaque
   partie par son etiquette. Exemple, la structure correspondant
   a la partie d'indice i sera l'element d'indice i.

   La dimension des voxels est prise en compte si
   le parametre correspondant (voxelSizes) n'est pas NULL.
   voxelSizes[0] : taille du point selon X
   voxelSizes[1] : taille du point selon Y
   voxelSizes[2] : taille du point selon Z


   La desallocation du tableau (par la procedure standard free)
   est a la charge de l'appelant.

   RETURN:
   renvoie NULL en cas d'erreur.
*/
   
   

extern typeParameter *ComputeParameterFromLabels( unsigned short int *theBuf,
					   const int *theDim,
					   const double *voxelSizes,
					   const int nbLabels );



extern void _SetFeretDiameterCalculusToComputation();
extern void _SetFeretDiameterCalculusToEstimation();


extern void _VerboseInCcParameters();
extern void _NoVerboseInCcParameters();

#ifdef __cplusplus
}
#endif

#endif /* _ccparameters_h_ */
