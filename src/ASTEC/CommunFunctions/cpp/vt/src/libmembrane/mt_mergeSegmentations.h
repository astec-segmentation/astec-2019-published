/*************************************************************************
 * mt_mergeSegmentations.h -
 *
 * $Id: mt_mergeSegmentations.h,v 1.6 2014/04/02 11:30:00 gael Exp $
 *
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2014/04/02
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _mt_mergeSeg_h_
#define _mt_mergeSeg_h_


#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include <vt_image.h>


extern int MT_BackgroundMask( vt_image imageIn, 
							  vt_image *imageMask,
							  int value);

extern int MT_ExtOnBackground(vt_image imageIn,
							  vt_image imageExt, 
							  vt_image imageMask,
			      int dilatation,			/* argument de dilatation pr chacun des labels de imageExt */
			      int *labelsExt,			/* liste des labels de imageExt qui intersectent la region non nulle d'imageMask */
			      int *volumes,				/* volumes des labels de labelsExt */
			      int *volumesIntersections,/* volume des intersections entre les labelsExt et le background */
			      int *nbVoisins,			/* nombre de labels de imageIn en contact avec les labelsExt */
			      int *n					/* nombre de labelsExt comptes... (size des tableaux) */
							  );



#ifdef __cplusplus
}
#endif

#endif /* _mt_mergeSeg_h_ */
