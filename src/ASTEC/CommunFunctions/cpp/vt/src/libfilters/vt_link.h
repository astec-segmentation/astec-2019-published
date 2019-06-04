/*************************************************************************
 * vt_link.h -
 *
 * $Id: vt_link.h,v 1.1 2000/10/20 11:15:53 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Wed Oct 11 12:00:17 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _vt_link_h_
#define _vt_link_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>


extern int *VT_ExtractLinkedCurve( vt_image *theIm,
			    int *firstPt,
			    int *nbPts, int sb, int sh );


extern int VT_ExtractFirstPoint( vt_image *theIm,
				 int *firstPt, int sb, int sh );



#ifdef __cplusplus
}
#endif

#endif /* _vt_link_h_ */

