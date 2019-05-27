#ifndef _vt_geline_h_
#define _vt_geline_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <vt_common.h>

#include <vt_recfilters.h>

typedef enum  {
  VT_BLACK = 1,
  VT_WHITE = 2
} BlackOrWhite;

typedef struct vt_line {
    /*--- multi-echelle ---*/
    float first_coeff;
    float last_coeff;
    int nb_coeff;
    /*--- type des vaisseaux ---*/
    BlackOrWhite type_structures;
    /*--- filtrage ---*/
    vt_recfilters par_filt;
} vt_line;

typedef struct vt_resline {
  /*--- declarations pour 2D et 3D ---*/
  vt_image imres;
  vt_image imext;
  vt_image imscale;
  vt_image imdirx;
  vt_image imdiry;
  vt_image imradp;
  vt_image imradn;
  /*--- declarations pour 3D seulement ----*/
  vt_image imdirz;
  vt_image imdir2x;
  vt_image imdir2y;
  vt_image imdir2z;
  vt_image imrad2p;
  vt_image imrad2n;
} vt_resline;

typedef struct vt_images {
  /*--- declarations pour 2D et 3D ---*/
  vt_image imx;
  vt_image imy;
  vt_image imxx;
  vt_image imxy;
  vt_image imyy;
  vt_image imr;
  vt_image ime;
  /*--- declarations pour 3D seulement ----*/
  vt_image imz;
  vt_image imxz;
  vt_image imyz;
  vt_image imzz;
} vt_images;



extern void VT_Line( vt_line *par );
extern int  VT_ExtractMaxima3D( vt_resline *res );
extern int  VT_ExtractMaxRad3D( vt_resline *res, vt_images *ims, vt_line *par );
extern int  VT_ExtractMaxima2D( vt_resline *res );
extern int  VT_ExtractMaxRad2D( vt_resline *res, vt_images *ims, vt_line *par );
extern int  VT_Reconstruct2D( vt_image *ext, vt_image *scl, vt_line *par );
extern int  VT_GreyReconstruct2D( vt_resline *res );
extern void VT_InitRes( vt_resline *par );
extern int  VT_AllocRes( vt_resline *par, vt_image *im, int dim );
extern void VT_FreeRes( vt_resline *par, int dim );
extern void VT_InitImages( vt_images *par );
extern int  VT_AllocImages( vt_images *par, vt_image *im, int dim );
extern void VT_FreeImages( vt_images *par, int dim );



#ifdef __cplusplus
}
#endif

#endif /* _vt_geline_h_ */
