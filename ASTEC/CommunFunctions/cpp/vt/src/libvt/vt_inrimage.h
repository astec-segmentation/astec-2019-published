/*************************************************************************
 * vt_inrimage.h -
 *
 * $Id: vt_inrimage.h,v 1.2 2000/09/07 07:42:33 greg Exp $
 *
 * Copyright (c) INRIA 1999
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * 
 *
 * ADDITIONS, CHANGES
 *
 */


#ifndef _vt_inrimage_h_
#define _vt_inrimage_h_


#ifdef __cplusplus
extern "C" {
#endif


#include <vt_image.h>



/* reading procedure
 */
extern int VT_ReadInrimage( vt_image *image, char *name );
extern vt_image *_VT_Inrimage( char *name );



/* writing procedure
 */
extern int VT_WriteInrimageWithName( vt_image *image, char *name );
extern int VT_WriteInrimage( vt_image *image );

extern int VT_WriteInrimageHeaderWithName( vt_image *image, char *name );
extern int VT_WriteInrimageHeader( vt_image *image );

#ifdef __cplusplus
}
#endif

#endif /* _vt_inrimage_h_ */
