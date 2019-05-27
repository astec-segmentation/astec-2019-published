
#ifndef _vt_contoursMatlab_h_
#define _vt_contoursMatlab_h_



#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <vt_contours.h>

extern void MAT_DrawContour2D( typeContour2D *c,
			       int rawDataFileDesc,
			       FILE *commandFile,
			       int d );

extern void MAT_DrawContour2Din3D( typeContour2D *c,
				   double z,
				   double *mat,
				   int rawDataFileDesc,
				   FILE *commandFile,
				   int d, char *options );

extern void MAT_DrawContour3D( typeContour3D *c,
			       double *mat,
			       int rawDataFileDesc,
			       FILE *commandFile,
			       int d );

#ifdef __cplusplus
}
#endif

#endif /* _vt_contoursMatlab_h_ */
