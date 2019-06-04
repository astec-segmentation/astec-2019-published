
#ifndef _vt_contours_h_
#define _vt_contours_h_



#ifdef __cplusplus
extern "C" {
#endif

#include <string.h>


typedef enum {
  _OPEN_,
  _CLOSED_
} enumTopologicalType;


typedef struct {
  double x;
  double y;
} typePoint2D;

typedef struct {
  int n;
  int nalloc;
  typePoint2D *thePts;
  enumTopologicalType topology;
} typeContour2D;

typedef struct {
  double x;
  double y;
  double z;
} typePoint3D;

typedef struct {
  int n;
  int nalloc;
  typePoint3D *thePts;
  enumTopologicalType topology;
} typeContour3D;

#define _POINTS_ 100

void initContour2D( typeContour2D *c );
int  allocContour2D( typeContour2D *c, int n );
void freeContour2D( typeContour2D *c );
void printContour2D( typeContour2D *c, int p );
int  addPointToContour2D( double x, double y, typeContour2D *c );
typeContour2D *copyContour2D( typeContour2D *c );

double surfaceContour2D( typeContour2D *c, double x, double y );

void initContour3D( typeContour3D *c );
int  allocContour3D( typeContour3D *c, int n );
void freeContour3D( typeContour3D *c );
void printContour3D( typeContour3D *c, int p );
int  addPointToContour3D( double x, double y, double z, typeContour3D *c );
typeContour3D *copyContour3D( typeContour3D *c );



typedef struct {
  double z;
  double *mat;
  int n;
  int nalloc;
  typeContour2D **theContours;
} typeSlice;

#define _CONTOURS2D_ 5

void initSlice( typeSlice *s );
int  addContour2DToSlice( typeContour2D *c, typeSlice *s );
void freeSlice( typeSlice *s );
void printSlice( typeSlice *s, int p );



typedef struct {
  int n;
  int nalloc;
  typeSlice **theSlices;
} typeStructure;

#define _SLICES_ 10

void initStructure( typeStructure *s );
int  addSliceToStructure( typeSlice *s,
			  typeStructure *structure );
int  addContour2DToStructure( typeContour2D *c,
			     typeStructure *structure, double z );
void freeStructure( typeStructure *s );
void printStructure( typeStructure *s, int p );





typedef struct {
  int n;
  int nalloc;
  typeContour3D **theContours;
} typeListOfContours3D;

#define _CONTOURS3D_ 10

void initListOfContours3D( typeListOfContours3D *l );
int  addContour3DToListOfContours3D( typeContour3D *c, 
				    typeListOfContours3D *l );
typeListOfContours3D* copyListOfContours3D( typeListOfContours3D *l );
void freeListOfContours3D( typeListOfContours3D *l );

int matlab_writeListOfContours3D( FILE *mfile, int rfile, typeListOfContours3D *l,
				  double *size,
				  char *desc, char *options );



int readStructure( typeStructure *s,
		   char *filename );
int writeStructure( typeStructure *s,
		    char *filename );

#ifdef __cplusplus
}
#endif

#endif /* _vt_contours_h_ */
