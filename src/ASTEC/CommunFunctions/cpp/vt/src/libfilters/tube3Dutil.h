/*************************************************************************
 * tube3Dutil.h -
 *
 * $Id: tube3Dutil.h,v 1.2 2001/06/08 15:10:52 greg Exp $
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

#ifndef _tube3Dutil_h_
#define _tube3Dutil_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
  /* #include <malloc.h> */
#include <string.h>


#include <typedefs.h>
#include <connexe.h>
#include <issimple3D.h>




/* Elimination des barbules.

   DESCRIPTION:
   On presuppose que l'image d'entree ne contient que des courbes amincies.
   De ce fait, la caracterisation topologique de tout point
   peut se faire en comptant le nombre de voisins d'un point :
   - 0 : point isole
   - 1 : extremite de courbe
   - 2 : point de courbe 
   - >2 : jonction entre courbes

   Dans un premier temps, on reconnait les composantes connexes des
   jonctions (4eme type de point) et des courbes pures (3 premiers
   types de points). 
   
   On elimine iterativement les courbes pures terminales, qui sont 
   definies comme des courbes contenant au plus 'maximal_length' 
   points et connectees a une seule jonction. C'est une definition 
   abusive (au sens topologique) car il faudrait aussi verifier que 
   cette courbe contient une extremite de courbe : notre definition
   permet aussi d'effacer des boucles.
   
   Apres avoir elimine une courbe terminale,
   - on met a jour la jonction attenante (des points peuvent etre 
     effaces)
   - on met a jour les courbes (2 courbes peuvent fusionnees)

   On arrete lorsqu'il n'y a plus de courbes terminales de longueur
   inferieure ou egale a 'maximal_length'.

   PARAMETRES:
   unsigned char *input_buffer : buffer de l'image 3D d'entree
   unsigned char *output_buffer : buffer de l'image 3D de sortie
      (doit deja etre alloue)
   int *theDim : tableau de 3 entiers donnant les
      dimensions selon X, Y et Z des buffers 3D.
      Chacun des buffers doit donc etre de taille
      theDim[0]*theDim[1]*theDim[2]
   int maximal_length : longueur maximale des branches a effacer
      (on ne garde que des branches terminales de longueur
      strictement superieure)

   RETURN:
   -1 en cas d'erreur
    1 si tout va bien
 */
extern int remove_small_simple_curves( unsigned char *input_buffer,
				unsigned char *output_buffer,
				int *theDim,
				int maximal_length );







/* Elimination des composantes qui ne sont pas des courbes pures.
   
   DESCRIPTION:
   On presuppose que l'image d'entree ne contient que des courbes amincies.
   De ce fait, la caracterisation topologique de tout point
   peut se faire en comptant le nombre de voisins d'un point :
   - 0 : point isole
   - 1 : extremite de courbe
   - 2 : point de courbe 
   - >2 : jonction entre courbes

   On elimine les composantes connexes qui contiennent 
   au moins un point de jonction ou qui n'ont pas 2 points
   extremites ou
   qui sont de taille
   trop petite.

   RETURN:
   -1 en cas d'erreur
    1 si tout va bien
 */
extern int remove_non_simple_components( unsigned char *input_buffer,
				  unsigned char *output_buffer,
				  int *theDim,
				  int minimal_size );





typedef struct {
  int size;
  int extremity_1[3];
  int extremity_2[3];
} typeFiberParameter;

/* calcul de parametres sur des courbes 

   On calcule les parametres suivants :
   - taille de la composante
   - coordonnees des 2 extremites de la composantes
   pour chacune des courbes. Les parametres sont ranges
   dans un tableau de structures de taille (*nb_fibers)+1
   ou (*nb_fibers) est le nombre de courbes pures ouvertes
   trouvees : les parametres des (*nb_fibers) courbes
   se trouvent dans les elements d'indice 1 a (*nb_fibers)
   (l'element d'indice 0 n'est pas utilise)
   
   RETURN :
   - pointeur sur un tableau de structures de type typeFiberParameter
     et de taille (*nb_fibers)+1 en cas de succes
   - NULL sinon
 */

extern typeFiberParameter * compute_fibers_parameters( unsigned char *input_buffer,
						int *theDim,
						int minimal_length,
						int *nb_fibers );

extern void _VerboseInTube3Dutil();
extern void _NoVerboseInTube3Dutil();



#ifdef __cplusplus
}
#endif

#endif /* _tube3Dutil_h_ */

