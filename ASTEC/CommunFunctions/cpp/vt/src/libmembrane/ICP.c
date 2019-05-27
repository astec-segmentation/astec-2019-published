/*************************************************************************
 * ICP.c -
 *
 * $Id: planeRegistration.c,v 1.0 2015/07/22 10:14:51 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2015/07/22
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vtmalloc.h>

#include <vt_common.h>
#include <bal-transformation.h>

#include <vt_planeRegistration.h>
#include <vt_tree.h>



#define IMAX 500
#define DELTA 0.01

typedef struct local_par {
  vt_names names; /*  image In, image Ext and initial transformation */
  vt_names namesOut; /*  trsf, residuals and residuals-img */
  vt_names nameTemplate; /*  template */
  int flag_affine;
  int flag_bary;
  int flagInitTrsf;
  bal_estimator estimator;
  int imax;
  int flagTemplate;
  int flagImOut;
  int flagImIn;
  int flagImExt;
  int flagLabels;
  double percentile;
} local_par;


/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );

static int _verbose_ = 0;


static char *usage = "[-ref %s | -ref-points %s] [-flo %s | -flo-points %s] [-init-trsf %s]\n\
\t [-trsf %s] [-residuals|-residuals-out %s] [-residuals-img %s [-template %s]] [-r %f]\n\
\t [-rigid | -affine | -translation | -[sim|similitude]] [-labels] [-bary | -barycentre | -barycentres]\n\
\t [-estimator-type|-estimator|-es-type wlts|lts|wls|ls]\n\
\t [-lts-fraction %lf] [-lts-deviation %f] [-lts-iterations %d]\n\
\t [-fluid-sigma|-lts-sigma[-ll|-hl] %lf %lf %lf] [-imax %d]\n\
\t [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t -flo %s : image-in (binary or labellized image)\n\
\t -ref %s : image-ext (binary or labellized image)\n\
\t -flo-points %s : file containing the floating point list (to be registered)\n\
\t -ref-points %s : file containing the reference point list (on which floating ones are registered)\n\
\t -init-trsf %s : fichier d'initialisation de la transformation T_ref<-flo a calculer\n\
\t -trsf %s : fichier dans lequel est enregistree la transformation T_ref<-flo calculee\n\
\t -residuals %s : fichier dans lequel est enregistree la liste (non triee) des residus\n\
\t -residuals-img %s : les voxels de l'image sortie valent les residus des voxels correspondants de l'image d'entree\n\
\t -template %s : image template pour l'option -residuals-img\n\
\t -r %f : extrait la valeur de residu du percentil correspondant (nombre entre 0 et 1)\n\
\t -labels : point-sets to be registered are paired with voxels with the same value\n\
\t -bary : barycentre-to-barycentre transformation initialisation\n\
\t ### transformation type ###\n\
\t -rigid : computes rigid transformation (set as default)\n\
\t -affine : computes affine transformation\n\
\t -translation : computes translation transformation\n\
\t -similitude : computes similitude transformation\n\
\t ### transformation estimation ###\n\
\t [-estimator-type|-estimator|-es-type %s] # transformation estimator\n\
\t wlts: weighted least trimmed squares\n\
\t lts: least trimmed squares\n\
\t wls: weighted least squares\n\
\t ls: least squares\n\
\t [-lts-fraction %lf] # for trimmed estimations, fraction of pairs that are kept\n\
\t [-lts-deviation %lf] # for trimmed estimations, defines the threshold to discard\n\
\t pairings, ie 'average + this_value * standard_deviation'\n\
\t [-lts-iterations %d] # for trimmed estimations, the maximal number of iterations\n\
\t [-fluid-sigma|-lts-sigma] %lf %lf %lf] # sigma for fluid regularization,\n\
\t ie field interpolation and regularization for pairings (only for vector field)\n\
\t [-imax] %d] # maximal number of iterations (default : 500)\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];

int main( int argc, char *argv[] )
{
  fprintf(stdout, "Warning : convention for output transformation should be reviewed. Expected output trsf is T_flo<-ref instead of T_ref<-flo.\n");
  local_par par;
  vt_image *imageIn;
  vt_image *imageExt;
  vt_image *imageTemplate;
  vt_image imResiduals;
  /* int flag_3D = 1; */
  /* char name[256]; */

  vt_pointList p_in;
  vt_pointList p_ext;
  vt_pointList *p_ext_list = NULL;

  int nClasses = 255;
  int i;
  vt_pointList cornersFlt;
  vt_pointList cornersTrsf;
  vt_pointList cornersOld;

  vt_point *pointTrsf, *pointOld;

  int sizx=0,sizy=0,sizz=0;
  int dimx=0,dimy=0,dimz=0;

  double *residuals=NULL;
  int *ind=NULL;

  FILE* fichier;

  /*--- initialisation des parametres ---*/
  MT_InitParam( &par );

  /*--- lecture des parametres ---*/
  MT_Parse( argc, argv, &par );

  if (par.flagImIn) /*  -in */
  {
      /*--- lecture de l'image d'entree in ---*/
      imageIn = _VT_Inrimage( par.names.in );
      if ( imageIn == (vt_image*)NULL )
        MT_ErrorParse("unable to read input image in\n", 0);
      sizx=imageIn->siz.x;
      sizy=imageIn->siz.y;
      sizz=imageIn->siz.z;
      dimx=imageIn->dim.x;
      dimy=imageIn->dim.y;
      dimz=imageIn->dim.z;

      /*--- operations eventuelles sur l'image d'entree ---*/
      if ( par.names.inv == 1 )  VT_InverseImage( imageIn );
      if ( par.names.swap == 1 ) VT_SwapImage( imageIn );


      if (VT_AllocPointListWithImage(&p_in, imageIn) != 1) {
          VT_FreeImage( imageIn );
          VT_Free( (void**)&imageIn);
          MT_ErrorParse("unable to allocate pointset image for image-in\n", 0);
          return(0);
      }

      if(VT_ExtractPointListWithImage(*imageIn, &p_in) != 1)
      {
          VT_FreeImage( imageIn );
          VT_Free( (void**)&imageIn);
          VT_FreePointList(&p_in);
          MT_ErrorParse("unable to compute pointset image for image-in\n", 0);
          return(0);
      }

      /*--- liberations memoires partielles ---*/
      VT_FreeImage(imageIn);
  }
  else  /*  -in-points */
  {
      /*--- allocation de la liste de points in ---*/
      if ( VT_AllocPointListWithFileName(&p_in, par.names.in) != 1) {
          MT_ErrorParse("unable to allocate pointset image for file-in\n", 0);
          return(0);
      }

      /*--- extraction de la liste de points in ---*/
      if(VT_ExtractPointListWithFileName(par.names.in, &p_in) != 1)
      {
          VT_FreePointList(&p_in);
          MT_ErrorParse("unable to compute pointset image for file-in\n", 0);
          return(0);
      }

  }


  if (par.flagImExt) /*  -ext */
  {
      /*--- lecture de l'image d'entree ext ---*/
      imageExt = _VT_Inrimage( par.names.ext );
      if ( imageExt == (vt_image*)NULL )
        MT_ErrorParse("unable to read input image ext\n", 0);

      /*--- operations eventuelles sur l'image d'entree ---*/
      if ( par.names.inv == 1 )  VT_InverseImage( imageExt );
      if ( par.names.swap == 1 ) VT_SwapImage( imageExt );

      if (par.flagLabels) {
      /*
      if (VT_AllocPointListWithImage(&p_ext, imageExt) != 1) {
          VT_FreeImage( imageExt );
          VT_Free( (void**)&imageExt);
          VT_FreePointList(&p_in);
          MT_ErrorParse("unable to allocate pointset image for image-ext\n", 0);
          return(0);
      }

      if(VT_ExtractPointListWithImage(*imageExt, &p_ext) != 1)
      {
          VT_FreeImage( imageExt );
          VT_Free( (void**)&imageExt);
          VT_FreePointList(&p_in);
          VT_FreePointList(&p_ext);
          MT_ErrorParse("unable to compute pointset image for image-ext\n", 0);
          return(0);
      }
      */
          if (VT_AllocPointListsWithImage(&p_ext_list, &nClasses, imageExt) != 1) {
              VT_FreeImage( imageExt );
              VT_Free( (void**)&imageExt);
              VT_FreePointList(&p_in);
              MT_ErrorParse("unable to allocate point lists image for image-ext\n", 0);
              return(0);
          }

          /* fprintf(stdout, "Point lists allocated ; nClasses = %d\n", nClasses); */

          if(VT_ExtractPointListsWithImage(*imageExt, &p_ext_list, nClasses) != 1)
          {
              VT_FreeImage( imageExt );
              VT_Free( (void**)&imageExt);
              VT_FreePointList(&p_in);
              for (i = 0 ; i < nClasses ; i++) {
                  VT_FreePointList(&(p_ext_list[i]));
              }
              vtfree(p_ext_list);
              p_ext_list = NULL;
              MT_ErrorParse("unable to compute point lists with image-ext\n", 0);
              return(0);
          }
          /* fprintf(stdout, "Point lists extracted\n"); */

      }
      else {
          if (VT_AllocPointListWithImage(&p_ext, imageExt) != 1) {
              VT_FreeImage( imageExt );
              VT_Free( (void**)&imageExt);
              VT_FreePointList(&p_in);
              MT_ErrorParse("unable to allocate pointset image for image-ext\n", 0);
              return(0);
          }

          if(VT_ExtractPointListWithImage(*imageExt, &p_ext) != 1)
          {
              VT_FreeImage( imageExt );
              VT_Free( (void**)&imageExt);
              VT_FreePointList(&p_in);
              VT_FreePointList(&p_ext);
              MT_ErrorParse("unable to compute pointset image for image-ext\n", 0);
              return(0);
          }
      }


      /*--- liberations memoires partielles ---*/
      VT_FreeImage(imageExt);
  }
  else /*  -ext-points */
  {
      if (par.flagLabels) {
          /*--- allocation de la liste de points ext ---*/
          if ( VT_AllocPointListsWithFileName(&p_ext_list, &nClasses, par.names.ext) != 1) {
              VT_FreePointList(&p_in);
              MT_ErrorParse("unable to allocate pointset image for file-ext\n", 0);
              return(0);
          }

          /*--- extraction de la liste de points ext ---*/
          if(VT_ExtractPointListsWithFileName(par.names.ext, &p_ext_list, nClasses) != 1)
          {
              VT_FreePointList(&p_in);
              for (i = 0 ; i < nClasses ; i++) {
                  VT_FreePointList(&(p_ext_list[i]));
              }
              vtfree(p_ext_list);
              p_ext_list = NULL;
              MT_ErrorParse("unable to compute pointset image for file-ext\n", 0);
              return(0);
          }
      }
      else {
          /*--- allocation de la liste de points ext ---*/
          if ( VT_AllocPointListWithFileName(&p_ext, par.names.ext) != 1) {
              VT_FreePointList(&p_in);
              MT_ErrorParse("unable to allocate pointset image for file-ext\n", 0);
              return(0);
          }

          /*--- extraction de la liste de points ext ---*/
          if(VT_ExtractPointListWithFileName(par.names.ext, &p_ext) != 1)
          {
              VT_FreePointList(&p_in);
              VT_FreePointList(&p_ext);
              MT_ErrorParse("unable to compute pointset image for file-ext\n", 0);
              return(0);
          }
      }

  }

  /*
  int _i;
  for (_i = 1000 ; _i < 1030 ; _i++)
  {
    pointOld = &(p_in.list[_i]);
    fprintf(stdout, "PointIn[%d].{x,y,z}    = { %f, %f, %f }\n", _i, pointOld->x, pointOld->y, pointOld->z);
    fprintf(stdout, "           .value}     = %d\n", pointOld->value);
    fprintf(stdout, "           .attribute} = %d\n", pointOld->attribute);
    fprintf(stdout, "           .weight}    = %f\n\n", pointOld->weight);
  }
  return(0);
  */

  /* if ( dimz == 1 ) */
  /*   flag_3D = 0; */

  /*--- calculs ---*/

  tree_node *node_ext;
  tree_node **nodes_ext;
  bal_transformation theTrsf;
  double T_ext_in[16];
  double residual;
  i = 0;

  /*  CRITERE D'ARRET : mouvement des 8 coins de l'image petit a l'echelle d'un voxel (ou de la bbox si liste de points lue dans un fichier) */

  VT_AllocPointList(&cornersFlt, 8);
  VT_AllocPointList(&cornersTrsf, 8);
  VT_AllocPointList(&cornersOld, 8);

  double xmin, xmax, ymin, ymax, zmin, zmax;
  if (par.flagImIn) {
    xmin=0;
    ymin=0;
    zmin=0;
    xmax=dimx*sizx;
    ymax=dimy*sizy;
    zmax=dimz*sizz;
  }
  else {
    int _i;
    vt_point *P;
    if(p_in.n<1) {
        fprintf(stderr, "point list in cardinal < 1, exiting\n");
        return(0);
    }
    P = &(p_in.list[0]);
    xmin = P->x;
    ymin = P->y;
    zmin = P->z;
    xmax = P->x;
    ymax = P->y;
    zmax = P->z;

    for (_i = 1 ; _i<p_in.n ; _i++)
    {
        P = &(p_in.list[_i]);
        if (xmin>P->x) xmin = P->x;
        if (ymin>P->y) ymin = P->y;
        if (zmin>P->z) zmin = P->z;
        if (xmax<P->x) xmax = P->x;
        if (ymax<P->y) ymax = P->y;
        if (zmax<P->z) zmax = P->z;
    }
  }

  pointTrsf = &(cornersFlt.list[0]);
  pointTrsf->x = xmin;  pointTrsf->y = ymin;    pointTrsf->z = zmin;    pointTrsf->attribute = 0;    pointTrsf->weight = 0; pointTrsf->value = 255;
  pointTrsf = &(cornersFlt.list[1]);
  pointTrsf->x = xmax;  pointTrsf->y = ymin;    pointTrsf->z = zmin;    pointTrsf->attribute = 0;    pointTrsf->weight = 0; pointTrsf->value = 255;
  pointTrsf = &(cornersFlt.list[2]);
  pointTrsf->x = xmin;  pointTrsf->y = ymax;    pointTrsf->z = zmin;    pointTrsf->attribute = 0;    pointTrsf->weight = 0; pointTrsf->value = 255;
  pointTrsf = &(cornersFlt.list[3]);
  pointTrsf->x = xmax;  pointTrsf->y = ymax;    pointTrsf->z = zmin;    pointTrsf->attribute = 0;    pointTrsf->weight = 0; pointTrsf->value = 255;
  pointTrsf = &(cornersFlt.list[4]);
  pointTrsf->x = xmin;  pointTrsf->y = ymin;    pointTrsf->z = zmax;    pointTrsf->attribute = 0;    pointTrsf->weight = 0; pointTrsf->value = 255;
  pointTrsf = &(cornersFlt.list[5]);
  pointTrsf->x = xmax;  pointTrsf->y = ymin;    pointTrsf->z = zmax;    pointTrsf->attribute = 0;    pointTrsf->weight = 0; pointTrsf->value = 255;
  pointTrsf = &(cornersFlt.list[6]);
  pointTrsf->x = xmin;  pointTrsf->y = ymax;    pointTrsf->z = zmax;    pointTrsf->attribute = 0;    pointTrsf->weight = 0; pointTrsf->value = 255;
  pointTrsf = &(cornersFlt.list[7]);
  pointTrsf->x = xmax;  pointTrsf->y = ymax;    pointTrsf->z = zmax;    pointTrsf->attribute = 0;    pointTrsf->weight = 0; pointTrsf->value = 255;

  for (i=0 ; i<8 ; i++) {
      pointTrsf = &(cornersFlt.list[i]);
      pointOld = &(cornersOld.list[i]);
      pointOld->x = pointTrsf->x; pointOld->y = pointTrsf->y; pointOld->z = pointTrsf->z;
      pointOld->attribute = pointTrsf->attribute; pointOld->weight = pointTrsf->weight; pointOld->value = pointTrsf->value;
  }

  /*  Construct KD-Tree(s) for points ext */

  if (par.flagLabels)
  {
      vt_kdtrees(&nodes_ext, p_ext_list, nClasses);
  }
  else
  {
      vt_kdtree(&node_ext, p_ext, 0);
  }

  /*for(i=0;i<nClasses;i++) {
      fprintf(stdout, "%d : \n", i);
      printTree(nodes_ext[i]);
  }*/

  if(par.flagInitTrsf)
  {
      BAL_InitTransformation( &theTrsf );
      /*  Lecture trsf initiale */
      if ( BAL_ReadTransformation( &theTrsf, par.names.out ) != 1 ) {
        VT_FreePointList(&p_in);
        if(par.flagLabels) {
           for (i = 0 ; i< nClasses ; i++)
             clearTree(&(nodes_ext[i]));
           vtfree(nodes_ext);
           nodes_ext=NULL;
        }
        else {
          clearTree(&node_ext);
        }

        fprintf( stderr, "%s: unable to read 'real' transformation '%s'\n", program, par.names.out );
        return (-1);
      }
      double *m = theTrsf.mat.m;
      for(i=0; i<16; i++)
        T_ext_in[i] = m[i];
  }
  else
  {
      for (i=0; i<16;i++)
        T_ext_in[i]=0.0;
      for(i=0;i<4;i++)
        T_ext_in[i*5]=1.0;
  }

  /* Option -bary : initial translation to get barycenters superimpose */
  if ( par.flag_bary )
  {
      double B_in_trsf[3], B_ext[3];
      double w;

      VT_PointListBarycentreInitTrsf(p_in, T_ext_in, B_in_trsf);

      if(par.flagLabels) {
          double B_tmp[3];
          double w_tmp;
          int j;
          VT_TreeBarycentre(nodes_ext[0], B_ext, &w);
          for (i=1 ; i<nClasses ; i++){
            VT_TreeBarycentre(nodes_ext[i], B_tmp, &w_tmp);
            for (j=0 ; j<3 ; j++) B_ext[j]=(B_ext[j]*w + B_tmp[j]*w_tmp)/(w + w_tmp);
            w+=w_tmp;
          }
      }
      else {
          VT_TreeBarycentre(nodes_ext[i], B_ext, &w);
      }
      if (_verbose_)
      fprintf(stdout, "B_in_trsf = (  %f   %f   %f  ), weight_in  = %f\nB_ext     = (  %f   %f   %f  ), weight_ext = %f\n",
              B_in_trsf[0], B_in_trsf[1], B_in_trsf[2], (float)p_in.n, B_ext[0], B_ext[1], B_ext[2], w);
      for (i=0 ; i<3 ; i++) T_ext_in[4*i+3] += B_ext[i] - B_in_trsf[i];
  }


  int cond = 1;

  i = 0;

  while(i<par.imax && cond) {
    i++;

    if (10<par.imax && (i % (par.imax/10) )== 0 )
    {
        fprintf(stdout, "ICP : %d %% \n", (int)(100*i/par.imax));
    }
    if (par.flagLabels) {


        if (VT_TreesComputeCP(p_in, p_ext_list, nodes_ext, nClasses, T_ext_in, T_ext_in, &residual, par.flag_affine, par.estimator) != 1) {
            VT_FreePointList(&p_in);
            for (i = 0 ; i< nClasses ; i++) {
              p_ext = p_ext_list[i];
              if (p_ext.list) {
                  VT_FreePointList(&(p_ext_list[i]));
              }
              clearTree(&(nodes_ext[i]));
            }
            vtfree(nodes_ext);
            nodes_ext=NULL;
            vtfree(p_ext_list);
            p_ext_list = NULL;

            MT_ErrorParse("unable to compute iteration #%d of ICP... Aborting\n", 0);
            return(0);
        }

    }
    else {
      if (VT_TreeComputeCP(p_in, p_ext, node_ext, T_ext_in, T_ext_in, &residual, par.flag_affine, par.estimator) != 1) {
        VT_FreePointList(&p_in);
        VT_FreePointList(&p_ext);
        clearTree(&node_ext);

        MT_ErrorParse("unable to compute iteration #%d of ICP... Aborting\n", 0);
        return(0);
      }
    }
    /*  Conditions d'arret : les coins de l'image flottante se deplacent de moins d'1/10 de voxel */

    VT_TransformPointList(T_ext_in, cornersFlt, &cornersTrsf);

    int j;

    int convergence = 1;

    for (j=0 ; j<8 ; j++) {

        pointTrsf = &(cornersTrsf.list[j]);
        pointOld = &(cornersOld.list[j]);

        double deltax = fabs(pointTrsf->x-pointOld->x);
        double deltay = fabs(pointTrsf->y-pointOld->y);
        double deltaz = fabs(pointTrsf->z-pointOld->z);

        /* if (0) */
        /*     fprintf(stdout, "deltax = %f\tdeltay = %f\tdeltaz = %f\n", deltax, deltay, deltaz); */

        if( deltax>sizx*DELTA ||
            deltay>sizy*DELTA ||
            deltaz>sizz*DELTA )
            convergence = 0;

        pointOld->x = pointTrsf->x; pointOld->y = pointTrsf->y; pointOld->z = pointTrsf->z;
        pointOld->attribute = pointTrsf->attribute; pointOld->weight = pointTrsf->weight;
    }

    cond = (convergence == 1) ? 0 : 1;

  }

  fprintf(stdout, "Final number of iterations : %d\n", i);

  VT_FreePointList(&cornersFlt);
  VT_FreePointList(&cornersTrsf);
  VT_FreePointList(&cornersOld);

  /* write outputs */

  /*if(par.flagImOut && applyTrsf(
                     par.names.ext,
                     par.namesOut.out,
                     p.real_transformation_name,
                     p.voxel_transformation_name,
                     p.result_real_transformation_name,
                     p.result_voxel_transformation_name,
                     p.template_image_name,
                     p.dim,
                     p.voxel,
                     p.resize,
                     p.interpolation,
                     p.type,
                     _debug_,
                     _verbose_
                     ))
  {
    fprintf( stderr, "%s: Failure.\n",program);
    return(0);
  }*/


  fichier = NULL;
  if(par.namesOut.in[0] != '\0') /*  trsf */
  {
      fichier = fopen(par.namesOut.in, "w");
      if( fichier == NULL)
      {
        MT_ErrorParse("Erreur pour l'ecriture du resultat dans un fichier\n", 0);
      }
  }
  else
  {
      fichier = stdout;
  }
  fprintf(fichier, "#\n# \tTransformation matrix (ext<-in):\n#\n");
  fprintf(fichier, "%f \t%f \t%f \t%f\n%f \t%f \t%f \t%f\n%f \t%f \t%f \t%f\n%f \t%f \t%f \t%f\n",
          T_ext_in[0],T_ext_in[1],T_ext_in[2],T_ext_in[3],
          T_ext_in[4],T_ext_in[5],T_ext_in[6],T_ext_in[7],
          T_ext_in[8],T_ext_in[9],T_ext_in[10],T_ext_in[11],
          T_ext_in[12],T_ext_in[13],T_ext_in[14],T_ext_in[15]);

  fprintf(fichier, "#\n");
  fprintf(fichier, "\n");
  fclose(fichier);


  /* Quantification : residuals vector computation */

  if (par.flagLabels) {
      if (VT_LabelsResidualsVector(p_in, nodes_ext, T_ext_in, &residuals) != 1)
      {
        if(p_in.list)
          VT_FreePointList(&p_in);
        for (i = 0 ; i < nClasses ; i++) {
            p_ext = p_ext_list[i];
            if (p_ext.list)
              VT_FreePointList(&(p_ext_list[i]));
            clearTree(&(nodes_ext[i]));
        }
        vtfree(nodes_ext);
        vtfree(p_ext_list);
        p_ext_list = NULL;
        MT_ErrorParse("Error while computing residuals\n", 0);
      }
  }
  else {
    if (VT_ResidualsVector(p_in, node_ext, T_ext_in, &residuals) != 1)
    {
      if(p_in.list)
        VT_FreePointList(&p_in);
      if (p_ext.list) {
        VT_FreePointList(&p_ext);
      }
      clearTree(&node_ext);
      MT_ErrorParse("Error while computing residuals\n", 0);
    }
  }

  /* Image of residuals */
  if (par.flagImOut)
  {
      if (par.flagTemplate)
      {
        /*--- lecture de l'image template ---*/
        imageTemplate = _VT_Inrimage( par.nameTemplate.in );
        if ( imageTemplate == (vt_image*)NULL )
        {
          vtfree(residuals); residuals=NULL;
          if(p_in.list)
            VT_FreePointList(&p_in);
          if(par.flagLabels) {
            for (i = 0 ; i < nClasses ; i++) {
                p_ext = p_ext_list[i];
                if (p_ext.list)
                  VT_FreePointList(&(p_ext_list[i]));
                clearTree(&(nodes_ext[i]));
            }
            vtfree(nodes_ext);
            vtfree(p_ext_list);
            p_ext_list = NULL;
          }
          else {
            if (p_ext.list) {
              VT_FreePointList(&p_ext);
            }
            clearTree(&node_ext);
          }
          MT_ErrorParse("unable to read input image template\n", 0);
        }


        /*--- init et allocation image de sortie ---*/

        VT_InitImage( &imResiduals, par.namesOut.out, imageTemplate->dim.x, imageTemplate->dim.y,
                      imageTemplate->dim.z, FLOAT );

        if ( VT_AllocImage( &(imResiduals) ) != 1 ) {
          VT_FreeImage(imageTemplate);
          if(p_in.list)
            VT_FreePointList(&p_in);
          if(par.flagLabels) {
            for (i = 0 ; i < nClasses ; i++) {
                p_ext = p_ext_list[i];
                if (p_ext.list)
                  VT_FreePointList(&(p_ext_list[i]));
                clearTree(&(nodes_ext[i]));
            }
            vtfree(nodes_ext);
            vtfree(p_ext_list);
            p_ext_list = NULL;
          }
          else {
            if (p_ext.list) {
              VT_FreePointList(&p_ext);
            }
            clearTree(&node_ext);
          }
          MT_ErrorParse("problem while allocating imResiduals\n", 0 );
        }

        imResiduals.siz.x=imageTemplate->siz.x;
        imResiduals.siz.y=imageTemplate->siz.y;
        imResiduals.siz.z=imageTemplate->siz.z;

        /*--- liberations memoires partielles ---*/
        VT_FreeImage(imageTemplate);
     }
     else
     {
        dimx=(int)xmax+1;
        dimy=(int)ymax+1;
        dimz=(int)zmax+1;

        /*--- init et allocation image de sortie ---*/

        VT_InitImage( &imResiduals, par.namesOut.out, dimx, dimy,
                      dimz, FLOAT );

        if ( VT_AllocImage( &(imResiduals) ) != 1 ) {
          if(p_in.list)
            VT_FreePointList(&p_in);
          if(par.flagLabels) {
            for (i = 0 ; i < nClasses ; i++) {
                p_ext = p_ext_list[i];
                if (p_ext.list)
                  VT_FreePointList(&(p_ext_list[i]));
                clearTree(&(nodes_ext[i]));
            }
            vtfree(nodes_ext);
            vtfree(p_ext_list);
            p_ext_list = NULL;
          }
          else {
            if (p_ext.list) {
              VT_FreePointList(&p_ext);
            }
            clearTree(&node_ext);
          }
          MT_ErrorParse("problem while allocating imResiduals\n", 0 );
        }

        imResiduals.siz.x=1;
        imResiduals.siz.y=1;
        imResiduals.siz.z=1;
     }

     VT_SetResidualImage(&imResiduals, p_in, residuals);

     VT_WriteInrimage(&imResiduals);
     VT_FreeImage(&imResiduals);
  }

  /* residuals vector sorting */

  if (VT_TriRapideVector(residuals, p_in.n, &ind) != 1)
  {
      vtfree(residuals); residuals = NULL;
      if(p_in.list)
        VT_FreePointList(&p_in);
      if(par.flagLabels) {
        for (i = 0 ; i < nClasses ; i++) {
            p_ext = p_ext_list[i];
            if (p_ext.list)
              VT_FreePointList(&(p_ext_list[i]));
            clearTree(&(nodes_ext[i]));
        }
        vtfree(nodes_ext);
        vtfree(p_ext_list);
        p_ext_list = NULL;
      }
      else {
        if (p_ext.list) {
          VT_FreePointList(&p_ext);
        }
        clearTree(&node_ext);
      }
      MT_ErrorParse("Error while sorting residuals\n", 0);
  }

  /* sorted residuals writting in output txt file */

  fichier = NULL;
  if(par.namesOut.ext[0] != '\0') /*  residuals */
  {
      fichier = fopen(par.namesOut.ext, "w");
      if( fichier == NULL)
      {
        MT_ErrorParse("Erreur pour l'ecriture des residus dans un fichier\n", 0);
      }
      for (i = 0 ; i<p_in.n ; i++)
        fprintf(fichier, "%f\n", residuals[ind[i]]);
  }
  else
  {
      fichier = stdout;
  }

  if (par.percentile>=0.0 ) {
    if (par.percentile<=1.0) {
      i = (int) (par.percentile * (p_in.n -1));
      if (i>=p_in.n) {
          fprintf(stderr, "Error : value %d is over the vector bound %d\n", i, p_in.n-1);
      }
      else {
        fprintf(fichier, "# Residuals %f percentile value : %f\n", par.percentile*100, residuals[ind[i]]);
      }
    }
    else
    {
        fprintf(stderr, "Unexpected option '-r' value %f : should be >= 0 and <= 1\n", par.percentile);
    }
  }

  fclose(fichier);


  /* free memory*/
  vtfree(ind); ind = NULL;
  vtfree(residuals); residuals = NULL;
  if(p_in.list)
    VT_FreePointList(&p_in);
  if(par.flagLabels) {
    for (i = 0 ; i < nClasses ; i++) {
        p_ext = p_ext_list[i];
        if (p_ext.list)
          VT_FreePointList(&(p_ext_list[i]));
        clearTree(&(nodes_ext[i]));
    }
    vtfree(nodes_ext);
    vtfree(p_ext_list);
    p_ext_list = NULL;
  }
  else {
    if (p_ext.list) {
      VT_FreePointList(&p_ext);
    }
    clearTree(&node_ext);
  }

  return( 0 );
}














static void MT_Parse( int argc,
                      char *argv[],
                      local_par *par )
{
  int i, nb;
  char text[STRINGLENGTH];
  int status;

  if ( VT_CopyName( program, argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);
  if ( argc == 1 ) MT_ErrorParse("\n", 0 );

  /*--- lecture des parametres ---*/
  i = 1;
  nb = 0;
  while ( i < argc ) {
    if ( argv[i][0] == '-' ) {

      /*--- arguments generaux ---*/
      if ( strcmp ( argv[i], "-help" ) == 0 ) {
        MT_ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 ) {
        _VT_VERBOSE_ = 1;
        _verbose_++;
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 ) {
        _VT_DEBUG_ = 1;
      }

      /*--- traitement eventuel de l'image d'entree ---*/
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
          par->names.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
          par->names.swap = 1;
      }


      /* Parametres de calcul */

      else if ( strcmp ( argv[i], "-imax" ) == 0 ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-imax", 0 );
        status = sscanf( argv[i], "%d", &(par->imax) );
        if ( status <= 0 ) MT_ErrorParse( "-imax", 0 );
      }

      else if ( strcmp ( argv[i], "-init-trsf" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -init-trsf...\n",0);
        sprintf(par->names.out, "%s", argv[i]);
        par->flagInitTrsf = 1;
      }

      else if ( strcmp ( argv[i], "-trsf" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -trsf...\n",0);
        sprintf(par->namesOut.in, "%s", argv[i]);
      }

      else if ( strcmp ( argv[i], "-residuals-out" ) == 0 || strcmp ( argv[i], "-residuals" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -residuals...\n",0);
        sprintf(par->namesOut.ext, "%s", argv[i]);
      }

      else if ( strcmp ( argv[i], "-r" ) == 0
            || strcmp ( argv[i], "-percentile" ) == 0) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-r", 0 );
        status = sscanf( argv[i], "%lf", &(par->percentile) );
        if ( status <= 0 ) MT_ErrorParse( "-r", 0 );
      }

      else if ( strcmp ( argv[i], "-residuals-img" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -residuals-img...\n",0);
        sprintf(par->namesOut.out, "%s", argv[i]);
        par->flagImOut=1;
      }

      else if ( strcmp ( argv[i], "-template" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -template...\n",0);
        sprintf(par->nameTemplate.in, "%s", argv[i]);
        par->flagTemplate=1;
      }

      else if ( strcmp ( argv[i], "-rigid" ) == 0 ){
        par->flag_affine=0;
      }

      else if ( strcmp ( argv[i], "-affine" ) == 0 ){
        par->flag_affine=1;
      }

      else if ( strcmp ( argv[i], "-translation" ) == 0 ){
        par->flag_affine=2;
      }

      else if ( strcmp ( argv[i], "-similitude" ) == 0 || strcmp ( argv[i], "-sim" ) == 0 ){
        par->flag_affine=3;
      }

      else if ( strcmp ( argv[i], "-bary" ) == 0 || strcmp ( argv[i], "-barycentre" ) == 0 || strcmp ( argv[i], "-barycentres" ) == 0 ){
        par->flag_bary=1;
      }

      /*ESTIMATOR PARAMETERS*/

      else if ( strcmp ( argv[i], "-estimator-type") == 0
            || strcmp ( argv[i], "-estimator") == 0
            || strcmp ( argv[i], "-es-type") == 0 ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-estimator-type", 0 );
        if ( (strcmp ( argv[i], "ltsw" ) == 0 && argv[i][4] == '\0')
         || (strcmp ( argv[i], "wlts" ) == 0 && argv[i][4] == '\0') ) {
      par->estimator.type = TYPE_WLTS;
        }
        else if ( strcmp ( argv[i], "lts" ) == 0 && argv[i][3] == '\0' ) {
      par->estimator.type = TYPE_LTS;
        }
        else if ( (strcmp ( argv[i], "lsw" ) == 0 && argv[i][3] == '\0')
          || (strcmp ( argv[i], "wls" ) == 0 && argv[i][3] == '\0') ) {
      par->estimator.type = TYPE_WLS;
        }
        else if ( strcmp ( argv[i], "ls" ) == 0 && argv[i][2] == '\0' ) {
      par->estimator.type = TYPE_LS;
        }
        else {
      fprintf( stderr, "unknown estimator type: '%s'\n", argv[i] );
      MT_ErrorParse( "-estimator-type", 0 );
        }
      }

      else if ( strcmp ( argv[i], "-lts-fraction" ) == 0
            || strcmp ( argv[i], "-lts-cut" ) == 0) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-lts-fraction", 0 );
        status = sscanf( argv[i], "%lf", &(par->estimator.retained_fraction) );
        if ( status <= 0 ) MT_ErrorParse( "-lts-fraction", 0 );
      }

      else if ( strcmp ( argv[i], "-lts-deviation" ) == 0 ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-lts-deviation", 0 );
        status = sscanf( argv[i], "%lf", &(par->estimator.standard_deviation_threshold) );
        if ( status <= 0 ) MT_ErrorParse( "-lts-deviation", 0 );
      }

      else if ( strcmp ( argv[i], "-lts-iterations" ) == 0 ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-lts-iterations", 0 );
        status = sscanf( argv[i], "%d", &(par->estimator.max_iterations) );
        if ( status <= 0 ) MT_ErrorParse( "-lts-iterations", 0 );
      }

      else if ( (strcmp (argv[i], "-fluid-sigma" ) == 0 && argv[i][12] == '\0')
            || (strcmp (argv[i], "-lts-sigma" ) == 0 && argv[i][10] == '\0') ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "parsing -lts-sigma %lf", 0 );
        status = sscanf( argv[i], "%lf", &(par->estimator.sigma.x) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -lts-sigma %lf", 0 );
        i ++;
        if ( i >= argc) {
      par->estimator.sigma.y = par->estimator.sigma.x;
      par->estimator.sigma.z = par->estimator.sigma.x;
        }
        else {
      status = sscanf( argv[i], "%lf", &(par->estimator.sigma.y) );
      if ( status <= 0 ) {
        i--;
        par->estimator.sigma.y = par->estimator.sigma.x;
        par->estimator.sigma.z = par->estimator.sigma.x;
      }
      else {
        i ++;
        if ( i >= argc) par->estimator.sigma.z = 0;
        else {
          status = sscanf( argv[i], "%lf", &(par->estimator.sigma.z) );
          if ( status <= 0 ) {
            i--;
            par->estimator.sigma.z = 0;
          }
        }
      }
        }
      }

      else if ( strcmp ( argv[i], "-labels" ) == 0 ) {
        par->flagLabels = 1;
      }


      else if ( strcmp ( argv[i], "-in" ) == 0 || strcmp ( argv[i], "-flo" ) == 0 ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-flo", 0 );
        strncpy( par->names.in, argv[i], STRINGLENGTH );
        nb++;
      }

      else if ( strcmp ( argv[i], "-ext" ) == 0 || strcmp ( argv[i], "-ref" ) == 0 ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-ref", 0 );
        strncpy( par->names.ext, argv[i], STRINGLENGTH );
        nb++;
      }


      else if ( strcmp ( argv[i], "-in-points" ) == 0 || strcmp ( argv[i], "-flo-points" ) == 0 ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-flo-points", 0 );
        strncpy( par->names.in, argv[i], STRINGLENGTH );
        par->flagImIn=0;
        nb++;
      }

      /*else if ( strcmp ( argv[i], "-ext" ) == 0 || strcmp ( argv[i], "-ref" ) == 0 ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-ref", 0 );
        strncpy( par->names.ext, argv[i], STRINGLENGTH );
        nb++;
      }*/

      else if ( strcmp ( argv[i], "-ext-points" ) == 0 || strcmp ( argv[i], "-ref-points" ) == 0 ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-ext-points", 0 );
        strncpy( par->names.ext, argv[i], STRINGLENGTH );
        par->flagImExt=0;
        nb++;
      }

      /*--- option inconnue ---*/
      else {
        sprintf(text,"unknown option %s\n",argv[i]);
        MT_ErrorParse(text, 0);
      }
    }
    /*--- saisie des noms d'images ---*/
    else if ( argv[i][0] != 0 ) {
      if ( nb == 0 ) {
        strncpy( par->names.in, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else if ( nb == 1 ) {
        strncpy( par->names.ext, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else {
        MT_ErrorParse("too much file names when parsing\n", 0 );
      }
    }
    i += 1;
  }

  /*--- s'il n'y a pas assez de noms ... ---*/
  if (nb == 0) {
      MT_ErrorParse("not enough file names when parsing\n", 0 );
  }
  if (nb == 1)
    strcpy( par->names.ext, ">" );  /* standart ext */


}




static void MT_ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"Erreur : %s",str);
  exit( 1 );
}








static void MT_InitParam( local_par *par )
{
  VT_Names( &(par->names) );
  VT_Names( &(par->namesOut) );
  VT_Names( &(par->nameTemplate) );
  par->flag_affine=0;
  BAL_InitEstimator( &(par->estimator) );
  par->estimator.type = TYPE_LS;
  par->imax = IMAX;
  par->flagTemplate=0;
  par->flagImOut=0;
  par->flagImIn=1;
  par->flagImExt=1;
  par->flagInitTrsf=0;
  par->flagLabels=0;
  par->percentile=-1.0;
  par->flag_bary=0;
}

