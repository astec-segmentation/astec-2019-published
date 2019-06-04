#ifndef _vt_isocontours_h_
#define _vt_isocontours_h_



#ifdef __cplusplus
extern "C" {
#endif


#include <vt_common.h>
#include <vt_contours.h>

typeListOfContours3D *VT_Compute2DIsoContoursInXYSlice( vt_image *theIm, 
					     int z, /* slice to be processed */
					     double threshold );
typeListOfContours3D *VT_Compute2DIsoContoursInXZSlice( vt_image *theIm, 
					     int z, /* slice to be processed */
					     double threshold );
typeListOfContours3D *VT_Compute2DIsoContoursInYZSlice( vt_image *theIm, 
					     int z, /* slice to be processed */
					     double threshold );

typeSlice *VT_Compute2DIsoContours( vt_image *theIm, 
				    int z, /* slice to be processed */
				    double threshold );

int VT_ExtractContour2DAndAddToSlice( vt_image *theIm, 
				      vt_image *theMarks,
				      int startx, int starty, int startz,
				      double threshold, 
				      int EDGE,
				      typeSlice *slice );

typeContour2D *VT_ExtractContour2D( vt_image *theIm, 
				    vt_image *theMarks,
				    int x, int y, int z,
				    double threshold, 
				    int DIR );


extern int _PointInContour( typeContour2D *c,
		     double x, double y );






extern int VT_ComputeNumberOfContoursPerCell( vt_image *theIm,
					      vt_image *resIm,
					      double threshold );

#ifdef __cplusplus
}
#endif

#endif /* _vt_isocontours_h_ */
