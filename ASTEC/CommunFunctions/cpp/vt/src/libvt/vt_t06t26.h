#ifndef _vt_t06t26_h_
#define _vt_t06t26_h_

#ifdef __cplusplus
extern "C" {
#endif

/* Evalue t06 et t26 dans un 26-voisinage.

DESCRIPTION:
   - t26 est le nombre de 26-composantes connexes
   de l'objet dans un 26-voisinage.
   Les points consideres sont ceux dont la valeur est 
   differente de 0. Le point central n'est pas compte
   (il est considere comme appartenant au fond).

   - t06 est le nombre de 6-composantes connexes 
   du fond dans un 18-voisinage.
   Les points consideres sont ceux dont la valeur est 
   egale a 0. Le point central n'est pas compte
   (il est considere comme appartenant a l'objet).
   On ne compte que les 6-composantes qui sont 
   6-adjacentes au point central.

   - Le voisinage est numerote comme suit
        0  1  2  -   9 10 11  -  18 19 20
        3  4  5  -  12 13 14  -  21 22 23
        6  7  8  -  15 16 17  -  24 25 26

*/
extern void Compute_T06_and_T26( int neighbors[27], int *t06, int *t26 );

#ifdef __cplusplus
}
#endif

#endif /* _vt_t06t26_h_ */
