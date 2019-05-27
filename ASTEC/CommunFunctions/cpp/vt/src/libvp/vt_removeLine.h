#ifndef _vt_removeLine_h_
#define _vt_removeLine_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <vt_image.h>



extern void VT_SetVerboseInVtRemoveLine( int v );
extern void VT_IncrementVerboseInVtRemoveLine(  );
extern void VT_DecrementVerboseInVtRemoveLine(  );

typedef enum enumRemoveLineMethod {
  _LOCAL_,
  _REGIONAL_,
  _GLOBAL_
} enumRemoveLineMethod;



typedef struct typeRemoveLineParameter {
  enumRemoveLineMethod method;
  float contrastSignificantFraction;
  float yRejectedFraction;
  float xzKeptFraction;
  int automatedChoices;
} typeRemoveLineParameter;

extern void initRemoveLineParameter( typeRemoveLineParameter *p );
extern void fprintfRemoveLineParameter( FILE *f, typeRemoveLineParameter *p );



typedef struct typeCorrection {
    int y;
    float a;
    float b;
} typeCorrection;

typedef struct typeCorrectionList {
  typeCorrection *data;
  int n_data;
  int n_allocated_data;
} typeCorrectionList;


extern void VT_InitCorrectionList( typeCorrectionList *l );
extern void VT_FreeCorrectionList( typeCorrectionList *l );
extern int VT_ReadCorrectionList( char *name, typeCorrectionList *l );
extern int VT_WriteCorrectionList( char *name, typeCorrectionList *l );


extern int VT_CorrectLines( vt_image *theIm, vt_image *resIm,
                            typeCorrectionList *l );


extern int VT_RemoveLines( vt_image *theIm, vt_image *resIm,
                           typeCorrectionList *l, typeRemoveLineParameter *p );


#ifdef __cplusplus
}
#endif

#endif
