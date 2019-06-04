/*************************************************************************
 * planeRegistration.c -
 *
 * $Id: planeRegistration.c,v 1.0 2014/09/17 13:37:51 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2014/09/17
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vtmalloc.h>

#include <vt_common.h>
#include <vt_planeRegistration.h>
#include <mt_dice.h>

/*#include <bal-estimator.h> */

#define TABLENGTH 2000
#define PI 3.141592654

typedef struct local_par {
  vt_names name_gradients;
  vt_names name_labels;
  vt_names name_planes;
  vt_names name_dices;
  int flag_name_planes_in;
  int flag_name_planes_ext;
  vt_names name_res;
  vt_names name_residuals;
  double p_in[4];
  double p_ext[4];
  int flag_planes_in;
  int flag_planes_ext;
  int flag_affine;
  bal_estimator estimator;
  double delta;
  int force;
  int npairs;
  int pairs[2][TABLENGTH];
  int dices_in;
  int dices_out;
  int background_in;
  int background_ext;
} local_par;


/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );

static int _verbose_ = 0;


static char *usage = "[-label-[ref|in] %s] [-label-[flo|ext] %s] [-background-[ref|in] %d] [-background-[flo|ext] %d]\n\
\t [-gradient-[ref|in] %s] [-gradient-[flo|ext] %s] [-delta %lf]\n\
\t [-pairs|pair|p] %d %d [...] | %s]\n\
\t [-p-[ref|in] %f %f %f %f] [-p-[flo|ext] %f %f %f %f] \n\
\t [-rigid | -affine] \n\
\t [-estimator-type|-estimator|-es-type wlts|lts|wls|ls]\n\
\t [-lts-fraction %lf] [-lts-deviation %f] [-lts-iterations %d]\n\
\t [-fluid-sigma|-lts-sigma[-ll|-hl] %lf %lf %lf]\n\
\t [-force-trsf direct|reversal|%d]\n\
\t [-pairs-out %s] [-pairs-in %s] [-trsf %s] [-dices %s] [-residuals %s]\n\
\t [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t -label-[ref|flo] %s : nom des fichiers image ou texte de labels ref|flo\n\
\t N.B.: le cas echeant, le fichier texte attendu doit respecter le format suivant : \n\
\t 1ere ligne (optionnel) : en-tête au format \"Voxelsize %f %f %f\" stipulant le voxelsize d'origine de l'image de laquelle proviennent les coordonnees barycentriques des labels\n\
\t liste des barycentres  : une ligne par barycentre ; chaque ligne respecte le format \"%d %f %f %f\" ou \"%d %f %f %f %f\" correspondant au label, coordonnee en x, en y, en z et (le cas echeant) au volume (ou poids) du label \n\
\t lignes de commentaires : commencent avec le caractere \"#\"\n\
\t -background-[ref|flo] %d : label de fond a ignorer dans l'image [ref|flo] correspondante (defaut = 0)\n\
\t -gradient-[ref|flo] %s : nom des images de gradients ref|flo (barycentre des labels si non specifie)\n\
\t -delta %lf : ignore les labels d'une distance > au parametre\n\
\t -pairs %d %d [...] | %s : associe les paires de labels specifiees (ref-flo [ref-flo [...]])\n\
\t -p-ref %f %f %f %f : equation du plan de l'image-ref (coordonnees reelles)\n\
\t -p-flo %f %f %f %f : equation du plan de l'image-flo (coordonnees reelles)\n\
\t -trsf %s : fichier dans lequel est enregistree la transformation T_flo<-ref calculee\n\
\t -pairs-in %s : fichier dans lequel sont enregistres les appariements utilises pour le calcul\n\
\t -pairs-out %s : fichier dans lequel sont enregistres les appariements trouves apres calcul\n\
\t -force-trsf : applique l'algorithme avec les transformations uniquement [directes|retournees],\n\
\t ou seulement sur l'iteration donnee\n\
\t -residuals %s : enregistre les valeurs de residus pour les paires de labels mis en correspondance\n\
\t -dices %s : calcule l'indice de similarite de Dice entre chaque paire de labels en correspondance ; \n\
\t             OPTION VALABLE UNIQUEMENT AVEC DES PARAMETRES -label-[ref|flo] DE TYPE IMAGE \n\
\t ### transformation type ###\n\
\t -rigid : computes rigid transformation (set as default)\n\
\t -affine : computes affine transformation (set as default)\n\
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
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{

  local_par par;
  vt_image *images_gradients[2];
  vt_image *images_labels[2];

  double G_in[3];
  double G_ext[3];
  double tmp;
  double R_in[9], R_ext[9], R_sym[9];
  double t_in[3], t_ext[3], t_sym[3];
  double T_theta_ext[16], T_ext_theta[16], T_sym[16], T_ext_in[16];
  double T_in_ext_min[16];
  double r, rmin=-1, *residus=NULL;
  double sd, sdmin=0, *std_deviations=NULL;
  int *pairs[2];
  int npairs;



  int i, imin=0;

  double theta, *thetas=NULL, delta;
  int n;

  vt_image_barycentres b_in;
  vt_image_barycentres b_ext;
  vt_image_barycentres b_in_trsf;
  vt_image_barycentres b_ext_trsf;

  /*--- initialisation des parametres ---*/
  MT_InitParam( &par );

  /*--- lecture des parametres ---*/
  MT_Parse( argc, argv, &par );

  if (par.npairs==0 && par.delta>0){
    if ((par.flag_planes_in == 0 && par.flag_name_planes_in == 0) || (par.flag_planes_ext == 0 && par.flag_name_planes_ext == 0 ))
      MT_ErrorParse("missing input plane specifications\n", 0);
    if ((par.flag_planes_in == 1 && par.flag_name_planes_in == 1) || (par.flag_planes_ext == 1 && par.flag_name_planes_ext == 1 ))
      MT_ErrorParse("conflicting input plane specifications\n", 0);
  }


  /* SI ENTREE = IMAGES */
  /*--- lecture des images labels d'entree ---*/
  images_labels[0] = _VT_Inrimage( par.name_labels.in );
  if ( images_labels[0] == (vt_image*)NULL ) {
    /* SI ENTREE = FICHIERS TEXTES avec lignes au format %d %lf %lf %lf */
      if (VT_SetImageBarycentresWithList(&b_in, par.name_labels.in, par.background_in) != 1) {
          MT_ErrorParse("unable to read input image or text labels ref\n", 0);
          return(0);
      }
  }
  else {
      /*--- operations eventuelles sur l'image d'entree ---*/
      if ( par.name_labels.inv == 1 )  {VT_InverseImage( images_labels[0]); }
      if ( par.name_labels.swap == 1 ) {VT_SwapImage( images_labels[0]); }


      /*--- operations sur les images de labels ---*/
      if (VT_AllocImageBarycentresWithImage(&b_in, *images_labels[0], par.background_in) != 1) {
          VT_FreeImage( images_labels[0] );
          VT_Free( (void**)&images_labels[0]);
          MT_ErrorParse("unable to allocate barycentre image for label-ref\n", 0);
          return(0);
      }

      if(VT_ComputeLabelBarycentres(*images_labels[0], &b_in, par.background_in) != 1)
      {
          VT_FreeImage( images_labels[0] );
          VT_Free( (void**)&images_labels[0]);
          VT_FreeImageBarycentres(&b_in);
          MT_ErrorParse("unable to compute barycentre image for label-ref\n", 0);
          return(0);
      }

      /*--- liberations memoires partielles ---*/
      VT_FreeImage( images_labels[0] );  VT_Free( (void**)&images_labels[0]);
  }


  if(par.delta>0)
      VT_ImageBarycentresCentralLayersExtraction(&b_in, par.p_in, par.delta);

  /*VT_PrintImageBarycentres(b_in);*/



  if(_VT_VERBOSE_) {
      vt_barycentre theBary;
      fprintf(stdout, "label-ref : %s\n", par.name_labels.in);
      fprintf(stdout, "\t- Number of labels : %d\n", b_in.n);
      fprintf(stdout, "\t- Voxel size : %f %f %f\n", b_in.vx, b_in.vy, b_in.vz);
      if( _verbose_>2 && b_in.n>0)
      {
          theBary=b_in.barycentres[0];
          fprintf(stdout, "\t- Gravity centers  : Pt=[ %f  %f  %f ]\t Label=%d\t Weight=%d\n", theBary.x, theBary.y, theBary.z, theBary.label, (int)theBary.weight);
          for (i=1; i<b_in.n; i++)
          {
              theBary=b_in.barycentres[i];
              fprintf(stdout, "\t                     Pt=[ %f  %f  %f ]\t Label=%d\t Weight=%d\n", theBary.x, theBary.y, theBary.z, theBary.label, (int)theBary.weight);
          }
      }
  }

  /* SI ENTREE = IMAGES */
  /*--- lecture des images labels d'entree ---*/
  images_labels[1] = _VT_Inrimage( par.name_labels.ext );
  if ( images_labels[1] == (vt_image*)NULL ) {
    /* SI ENTREE = FICHIERS TEXTES avec lignes au format %d %lf %lf %lf */
      if (VT_SetImageBarycentresWithList(&b_ext, par.name_labels.ext, par.background_ext) != 1) {
          if ( images_labels[0] != (vt_image*)NULL ) {
              VT_FreeImage( images_labels[0] );
              VT_Free( (void**)&images_labels[0]);
          }
          MT_ErrorParse("unable to read input image or text labels flo\n", 0);
          return(0);
      }
  }
  else{
      /*--- operations eventuelles sur l'image d'entree ---*/
      if ( par.name_labels.inv == 1 )  {VT_InverseImage( images_labels[1]); }
      if ( par.name_labels.swap == 1 ) {VT_SwapImage( images_labels[1]); }


      if (VT_AllocImageBarycentresWithImage(&b_ext, *images_labels[1], par.background_ext) != 1) {
          VT_FreeImage( images_labels[1] );
          VT_Free( (void**)&images_labels[1]);
          VT_FreeImageBarycentres(&b_in);
          return(0);
      }

      if(VT_ComputeLabelBarycentres(*images_labels[1], &b_ext, par.background_ext) != 1)
      {
          VT_FreeImage( images_labels[1] );
          VT_Free( (void**)&images_labels[1]);
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_ext);
          MT_ErrorParse("unable to compute barycentre image for label-flo\n", 0);
          return(0);
      }

      /*--- liberations memoires partielles ---*/
      VT_FreeImage( images_labels[1] );  VT_Free( (void**)&images_labels[1]);
  }

  if(par.delta>0)
      VT_ImageBarycentresCentralLayersExtraction(&b_ext, par.p_ext, par.delta);

  if(_VT_VERBOSE_) {
      vt_barycentre theBary;
      fprintf(stdout, "label-flo : %s\n", par.name_labels.ext);
      fprintf(stdout, "\t- Number of labels : %d\n", b_ext.n);
      fprintf(stdout, "\t- Voxel size : %f %f %f\n", b_ext.vx, b_ext.vy, b_ext.vz);
      if(_verbose_>2 && b_ext.n>0)
      {
          theBary=b_ext.barycentres[0];
          fprintf(stdout, "\t- Gravity centers  : Pt=[ %f  %f  %f ]\t Label=%d\t Weight=%d\n", theBary.x, theBary.y, theBary.z, theBary.label, (int)theBary.weight);
          for (i=1; i<b_ext.n; i++)
          {
              theBary=b_ext.barycentres[i];
              fprintf(stdout, "\t                     Pt=[ %f  %f  %f ]\t Label=%d\t Weight=%d\n", theBary.x, theBary.y, theBary.z, theBary.label, (int)theBary.weight);
          }
      }
  }

  /* Allocation des images des barycentres transformes */

  if (VT_AllocImageBarycentres(&b_in_trsf, b_in.n) != 1) {
    VT_FreeImageBarycentres(&b_in);
    VT_FreeImageBarycentres(&b_ext);
    MT_ErrorParse("unable to allocate image of barycentres b_in_trsf...\n", 0);
  }

  if(par.npairs==0) {

      if (par.name_gradients.in[0]!='\0') {
        /*--- lecture de l'image gradient in ---*/
        images_gradients[0] = _VT_Inrimage( par.name_gradients.in );
        if ( images_gradients[0] == (vt_image*)NULL )
          MT_ErrorParse("unable to read input image gradient ref\n", 0);

        /*--- operations eventuelles sur l'image d'entree ---*/
        if ( par.name_gradients.inv == 1 )  VT_InverseImage( images_gradients[0] );
        if ( par.name_gradients.swap == 1 ) VT_SwapImage( images_gradients[0] );

        /*--- suite des calculs ---*/

        if (VT_Barycentre(images_gradients[0], G_in) != 1) {
          VT_FreeImage( images_gradients[0] );
          VT_Free( (void**)&images_gradients[0]);
          MT_ErrorParse("unable to compute gradient-ref barycentre\n", 0);
        }



        /*--- liberations memoires partielles ---*/
        VT_FreeImage( images_gradients[0] );
        VT_Free( (void**)&images_gradients[0]);
      }
      else {
        VT_BaryImgLabels(b_in, G_in);
      }

      if(par.name_gradients.ext[0]!='\0') {
        /*--- lecture de l'image gradient ext ---*/
        images_gradients[1] = _VT_Inrimage( par.name_gradients.ext );
        if ( images_gradients[1] == (vt_image*)NULL ) {
          MT_ErrorParse("unable to read input image gradient flo\n", 0);
        }

        /*--- operations eventuelles sur l'image d'entree ---*/
        if ( par.name_gradients.inv == 1 )  VT_InverseImage( images_gradients[1] );
        if ( par.name_gradients.swap == 1 ) VT_SwapImage( images_gradients[1] );

        if(VT_Barycentre(images_gradients[1], G_ext) != 1) {
          VT_FreeImage( images_gradients[1] );
          VT_Free( (void**)&images_gradients[1]);
          MT_ErrorParse("unable to compute gradient-flo barycentre\n", 0);
        }

        /*--- liberations memoires partielles ---*/
        VT_FreeImage( images_gradients[1] );
        VT_Free( (void**)&images_gradients[1]);
      }
      else {
        VT_BaryImgLabels(b_ext, G_ext);
      }

      if (_VT_VERBOSE_) {
        fprintf(stdout, "Barycentre_ref  = [ %f   %f   %f ]\n", G_in[0],G_in[1],G_in[2]);
        fprintf(stdout, "Barycentre_flo = [ %f   %f   %f ]\n", G_ext[0],G_ext[1],G_ext[2]);
      }

      /* Barycentre projections */

      tmp=par.p_in[0]*G_in[0]+par.p_in[1]*G_in[1]+par.p_in[2]*G_in[2]+par.p_in[3];
      G_in[0] -= par.p_in[0]*tmp;
      G_in[1] -= par.p_in[1]*tmp;
      G_in[2] -= par.p_in[2]*tmp;

      tmp=par.p_ext[0]*G_ext[0]+par.p_ext[1]*G_ext[1]+par.p_ext[2]*G_ext[2]+par.p_ext[3];
      G_ext[0] -= par.p_ext[0]*tmp;
      G_ext[1] -= par.p_ext[1]*tmp;
      G_ext[2] -= par.p_ext[2]*tmp;


      if (_VT_VERBOSE_) {
        fprintf(stdout, "Barycentre_ref projected  = [ %f   %f   %f ]\n", G_in[0],G_in[1],G_in[2]);
        fprintf(stdout, "Barycentre_flo projected = [ %f   %f   %f ]\n", G_ext[0],G_ext[1],G_ext[2]);
      }

      /*VT_PrintImageBarycentres(b_in); */
      /*VT_PrintImageBarycentres(b_ext); */

      /* Transformation Matrix */
      /* Matrice de retournement */
      VT_RotationMatrix180(par.p_ext, R_sym);
      VT_Translation(G_ext, G_ext, R_sym, t_sym);

      T_sym[0]=R_sym[0]; T_sym[1]=R_sym[1]; T_sym[2]=R_sym[2];  T_sym[3]=t_sym[0];
      T_sym[4]=R_sym[3]; T_sym[5]=R_sym[4]; T_sym[6]=R_sym[5];  T_sym[7]=t_sym[1];
      T_sym[8]=R_sym[6]; T_sym[9]=R_sym[7]; T_sym[10]=R_sym[8]; T_sym[11]=t_sym[2];
      T_sym[12]=0;       T_sym[13]=0;       T_sym[14]=0;        T_sym[15]=1;

      if (1 && _verbose_>3) {
        fprintf(stdout, "#Transformation matrix (sym<->flo):\n");
        fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_sym[0],T_sym[1],T_sym[2],T_sym[3]);
        fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_sym[4],T_sym[5],T_sym[6],T_sym[7]);
        fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_sym[8],T_sym[9],T_sym[10],T_sym[11]);
        fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_sym[12],T_sym[13],T_sym[14],T_sym[15]);
      }

      /* Matrice de transformation */
      VT_RotationMatrix(par.p_in, par.p_ext, R_ext);
      VT_Translation(G_in, G_ext, R_ext, t_ext);

      T_ext_theta[0]=R_ext[0]; T_ext_theta[1]=R_ext[1]; T_ext_theta[2]=R_ext[2];  T_ext_theta[3]=t_ext[0];
      T_ext_theta[4]=R_ext[3]; T_ext_theta[5]=R_ext[4]; T_ext_theta[6]=R_ext[5];  T_ext_theta[7]=t_ext[1];
      T_ext_theta[8]=R_ext[6]; T_ext_theta[9]=R_ext[7]; T_ext_theta[10]=R_ext[8]; T_ext_theta[11]=t_ext[2];
      T_ext_theta[12]=0;       T_ext_theta[13]=0;       T_ext_theta[14]=0;        T_ext_theta[15]=1;

      if (1 && _verbose_>3) {
        fprintf(stdout, "#Transformation matrix (flo<-theta):\n");
        fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_ext_theta[0],T_ext_theta[1],T_ext_theta[2],T_ext_theta[3]);
        fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_ext_theta[4],T_ext_theta[5],T_ext_theta[6],T_ext_theta[7]);
        fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_ext_theta[8],T_ext_theta[9],T_ext_theta[10],T_ext_theta[11]);
        fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_ext_theta[12],T_ext_theta[13],T_ext_theta[14],T_ext_theta[15]);
      }

      VT_RotationMatrix(par.p_ext, par.p_in, R_in);
      VT_Translation(G_ext, G_in, R_in, t_in);

      T_theta_ext[0]=R_in[0]; T_theta_ext[1]=R_in[1]; T_theta_ext[2]=R_in[2];  T_theta_ext[3]=t_in[0];
      T_theta_ext[4]=R_in[3]; T_theta_ext[5]=R_in[4]; T_theta_ext[6]=R_in[5];  T_theta_ext[7]=t_in[1];
      T_theta_ext[8]=R_in[6]; T_theta_ext[9]=R_in[7]; T_theta_ext[10]=R_in[8]; T_theta_ext[11]=t_in[2];
      T_theta_ext[12]=0;      T_theta_ext[13]=0;      T_theta_ext[14]=0;       T_theta_ext[15]=1;

      if (1 && _verbose_>3) {
        fprintf(stdout, "#Transformation matrix (theta<-flo):\n");
        fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_theta_ext[0],T_theta_ext[1],T_theta_ext[2],T_theta_ext[3]);
        fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_theta_ext[4],T_theta_ext[5],T_theta_ext[6],T_theta_ext[7]);
        fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_theta_ext[8],T_theta_ext[9],T_theta_ext[10],T_theta_ext[11]);
        fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_theta_ext[12],T_theta_ext[13],T_theta_ext[14],T_theta_ext[15]);
       }


      /* Allocation des images des barycentres transformes */

      if (VT_AllocImageBarycentres(&b_ext_trsf, b_ext.n) != 1) {
        VT_FreeImageBarycentres(&b_in);
        VT_FreeImageBarycentres(&b_ext);
        VT_FreeImageBarycentres(&b_in_trsf);
        MT_ErrorParse("unable to allocate image of barycentres b_in_trsf...\n", 0);
      }
      /* Liste d'appariements : boucle sur les angles */

      delta=PI/64;
      n=(int)2*PI/delta;

      thetas = vtmalloc( 2*n*sizeof(double), "thetas", argv[0] );
      if (thetas==(double*)NULL) {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext);
          VT_FreeImageBarycentres(&b_ext_trsf);
          MT_ErrorParse("unable to allocate thetas vector", 0);
      }
      residus = vtmalloc( 2*n*sizeof(double), "residus", argv[0] );
      if (residus==(double*)NULL) {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext);
          VT_FreeImageBarycentres(&b_ext_trsf);
          vtfree(thetas); thetas=NULL;
          MT_ErrorParse("unable to allocate residus vector", 0);
      }
      std_deviations = vtmalloc( 2*n*sizeof(double), "std_deviations", argv[0] );
      if (std_deviations==(double*)NULL) {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext);
          VT_FreeImageBarycentres(&b_ext_trsf);
          vtfree(thetas); thetas=NULL;
          vtfree(residus); residus=NULL;
          MT_ErrorParse("unable to allocate std_deviations vector", 0);
      }

      for (i=0 ; i<2*n ; i++) {

        if (par.force==-2 && i%2 == 1) continue;
        if (par.force==-3 && i%2 == 0) continue;
        if (par.force>=0 && i!=par.force) continue;

        double R_theta[9];          /* Matrice de rotation de in d'angle theta autour de par.p_in[0..2] */
        double t_theta[3];          /* Translation de in vers image tournee d'angle theta */

        double T_in_theta[16];         /*POINT_theta = T_theta * POINT_in */

        double T_in_theta_ext[16];  /*POINT_theta = T_in_theta_ext * POINT_ext */

        double T_in_theta_ext_sym[16];/* Combination with the possible embryo rotation symmetry transformation */

        double T_tmp[16];

        theta = ((double)((int)(0.5*i)))*delta;
        thetas[i]=theta;

        VT_RotationMatrixWithAngle(par.p_in, theta, R_theta);
        VT_Translation(G_in, G_in, R_theta, t_theta);
        T_in_theta[0]=R_theta[0]; T_in_theta[1]=R_theta[1]; T_in_theta[2]=R_theta[2];  T_in_theta[3]=t_theta[0];
        T_in_theta[4]=R_theta[3]; T_in_theta[5]=R_theta[4]; T_in_theta[6]=R_theta[5];  T_in_theta[7]=t_theta[1];
        T_in_theta[8]=R_theta[6]; T_in_theta[9]=R_theta[7]; T_in_theta[10]=R_theta[8]; T_in_theta[11]=t_theta[2];
        T_in_theta[12]=0;      T_in_theta[13]=0;      T_in_theta[14]=0;       T_in_theta[15]=1;

        if (1 && _verbose_>4) {
          fprintf(stdout, "#Transformation matrix (ref<-theta = %f * PI):\n", theta/PI);
          fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_in_theta[0],T_in_theta[1],T_in_theta[2],T_in_theta[3]);
          fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_in_theta[4],T_in_theta[5],T_in_theta[6],T_in_theta[7]);
          fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_in_theta[8],T_in_theta[9],T_in_theta[10],T_in_theta[11]);
          fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_in_theta[12],T_in_theta[13],T_in_theta[14],T_in_theta[15]);
         }


        T_in_theta_ext[0] = T_in_theta[0] * T_theta_ext[0] + T_in_theta[1] * T_theta_ext[4] + T_in_theta[2] * T_theta_ext[8] + T_in_theta[3] * T_theta_ext[12];
        T_in_theta_ext[1] = T_in_theta[0] * T_theta_ext[1] + T_in_theta[1] * T_theta_ext[5] + T_in_theta[2] * T_theta_ext[9] + T_in_theta[3] * T_theta_ext[13];
        T_in_theta_ext[2] = T_in_theta[0] * T_theta_ext[2] + T_in_theta[1] * T_theta_ext[6] + T_in_theta[2] * T_theta_ext[10]+ T_in_theta[3] * T_theta_ext[14];
        T_in_theta_ext[3] = T_in_theta[0] * T_theta_ext[3] + T_in_theta[1] * T_theta_ext[7] + T_in_theta[2] * T_theta_ext[11]+ T_in_theta[3] * T_theta_ext[15];

        T_in_theta_ext[4] = T_in_theta[4] * T_theta_ext[0] + T_in_theta[5] * T_theta_ext[4] + T_in_theta[6] * T_theta_ext[8] + T_in_theta[7] * T_theta_ext[12];
        T_in_theta_ext[5] = T_in_theta[4] * T_theta_ext[1] + T_in_theta[5] * T_theta_ext[5] + T_in_theta[6] * T_theta_ext[9] + T_in_theta[7] * T_theta_ext[13];
        T_in_theta_ext[6] = T_in_theta[4] * T_theta_ext[2] + T_in_theta[5] * T_theta_ext[6] + T_in_theta[6] * T_theta_ext[10]+ T_in_theta[7] * T_theta_ext[14];
        T_in_theta_ext[7] = T_in_theta[4] * T_theta_ext[3] + T_in_theta[5] * T_theta_ext[7] + T_in_theta[6] * T_theta_ext[11]+ T_in_theta[7] * T_theta_ext[15];

        T_in_theta_ext[8] = T_in_theta[8] * T_theta_ext[0] + T_in_theta[9] * T_theta_ext[4] + T_in_theta[10]* T_theta_ext[8] + T_in_theta[11]* T_theta_ext[12];
        T_in_theta_ext[9] = T_in_theta[8] * T_theta_ext[1] + T_in_theta[9] * T_theta_ext[5] + T_in_theta[10]* T_theta_ext[9] + T_in_theta[11]* T_theta_ext[13];
        T_in_theta_ext[10]= T_in_theta[8] * T_theta_ext[2] + T_in_theta[9] * T_theta_ext[6] + T_in_theta[10]* T_theta_ext[10]+ T_in_theta[11]* T_theta_ext[14];
        T_in_theta_ext[11]= T_in_theta[8] * T_theta_ext[3] + T_in_theta[9] * T_theta_ext[7] + T_in_theta[10]* T_theta_ext[11]+ T_in_theta[11]* T_theta_ext[15];

        T_in_theta_ext[12]= T_in_theta[12]* T_theta_ext[0] + T_in_theta[13]* T_theta_ext[4] + T_in_theta[14]* T_theta_ext[8] + T_in_theta[15]* T_theta_ext[12];
        T_in_theta_ext[13]= T_in_theta[12]* T_theta_ext[1] + T_in_theta[13]* T_theta_ext[5] + T_in_theta[14]* T_theta_ext[9] + T_in_theta[15]* T_theta_ext[13];
        T_in_theta_ext[14]= T_in_theta[12]* T_theta_ext[2] + T_in_theta[13]* T_theta_ext[6] + T_in_theta[14]* T_theta_ext[10]+ T_in_theta[15]* T_theta_ext[14];
        T_in_theta_ext[15]= T_in_theta[12]* T_theta_ext[3] + T_in_theta[13]* T_theta_ext[7] + T_in_theta[14]* T_theta_ext[11]+ T_in_theta[15]* T_theta_ext[15];

        if (1 && _verbose_>4) {
          fprintf(stdout, "#Transformation matrix (ref<-flo, theta = %f * PI):\n", theta/PI);
          fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_in_theta_ext[0],T_in_theta_ext[1],T_in_theta_ext[2],T_in_theta_ext[3]);
          fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_in_theta_ext[4],T_in_theta_ext[5],T_in_theta_ext[6],T_in_theta_ext[7]);
          fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_in_theta_ext[8],T_in_theta_ext[9],T_in_theta_ext[10],T_in_theta_ext[11]);
          fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_in_theta_ext[12],T_in_theta_ext[13],T_in_theta_ext[14],T_in_theta_ext[15]);
        }

        if (i%2==1)
        {
            T_in_theta_ext_sym[0] = T_in_theta_ext[0] * T_sym[0] + T_in_theta_ext[1] * T_sym[4] + T_in_theta_ext[2] * T_sym[8] + T_in_theta_ext[3] * T_sym[12];
            T_in_theta_ext_sym[1] = T_in_theta_ext[0] * T_sym[1] + T_in_theta_ext[1] * T_sym[5] + T_in_theta_ext[2] * T_sym[9] + T_in_theta_ext[3] * T_sym[13];
            T_in_theta_ext_sym[2] = T_in_theta_ext[0] * T_sym[2] + T_in_theta_ext[1] * T_sym[6] + T_in_theta_ext[2] * T_sym[10]+ T_in_theta_ext[3] * T_sym[14];
            T_in_theta_ext_sym[3] = T_in_theta_ext[0] * T_sym[3] + T_in_theta_ext[1] * T_sym[7] + T_in_theta_ext[2] * T_sym[11]+ T_in_theta_ext[3] * T_sym[15];

            T_in_theta_ext_sym[4] = T_in_theta_ext[4] * T_sym[0] + T_in_theta_ext[5] * T_sym[4] + T_in_theta_ext[6] * T_sym[8] + T_in_theta_ext[7] * T_sym[12];
            T_in_theta_ext_sym[5] = T_in_theta_ext[4] * T_sym[1] + T_in_theta_ext[5] * T_sym[5] + T_in_theta_ext[6] * T_sym[9] + T_in_theta_ext[7] * T_sym[13];
            T_in_theta_ext_sym[6] = T_in_theta_ext[4] * T_sym[2] + T_in_theta_ext[5] * T_sym[6] + T_in_theta_ext[6] * T_sym[10]+ T_in_theta_ext[7] * T_sym[14];
            T_in_theta_ext_sym[7] = T_in_theta_ext[4] * T_sym[3] + T_in_theta_ext[5] * T_sym[7] + T_in_theta_ext[6] * T_sym[11]+ T_in_theta_ext[7] * T_sym[15];

            T_in_theta_ext_sym[8] = T_in_theta_ext[8] * T_sym[0] + T_in_theta_ext[9] * T_sym[4] + T_in_theta_ext[10]* T_sym[8] + T_in_theta_ext[11]* T_sym[12];
            T_in_theta_ext_sym[9] = T_in_theta_ext[8] * T_sym[1] + T_in_theta_ext[9] * T_sym[5] + T_in_theta_ext[10]* T_sym[9] + T_in_theta_ext[11]* T_sym[13];
            T_in_theta_ext_sym[10]= T_in_theta_ext[8] * T_sym[2] + T_in_theta_ext[9] * T_sym[6] + T_in_theta_ext[10]* T_sym[10]+ T_in_theta_ext[11]* T_sym[14];
            T_in_theta_ext_sym[11]= T_in_theta_ext[8] * T_sym[3] + T_in_theta_ext[9] * T_sym[7] + T_in_theta_ext[10]* T_sym[11]+ T_in_theta_ext[11]* T_sym[15];

            T_in_theta_ext_sym[12]= T_in_theta_ext[12]* T_sym[0] + T_in_theta_ext[13]* T_sym[4] + T_in_theta_ext[14]* T_sym[8] + T_in_theta_ext[15]* T_sym[12];
            T_in_theta_ext_sym[13]= T_in_theta_ext[12]* T_sym[1] + T_in_theta_ext[13]* T_sym[5] + T_in_theta_ext[14]* T_sym[9] + T_in_theta_ext[15]* T_sym[13];
            T_in_theta_ext_sym[14]= T_in_theta_ext[12]* T_sym[2] + T_in_theta_ext[13]* T_sym[6] + T_in_theta_ext[14]* T_sym[10]+ T_in_theta_ext[15]* T_sym[14];
            T_in_theta_ext_sym[15]= T_in_theta_ext[12]* T_sym[3] + T_in_theta_ext[13]* T_sym[7] + T_in_theta_ext[14]* T_sym[11]+ T_in_theta_ext[15]* T_sym[15];
        }
        else {
            int j;
            for (j=0; j<16; j++)
                T_in_theta_ext_sym[j]=T_in_theta_ext[j];
        }


        if (VT_TransformImageBarycentres(T_in_theta_ext_sym, b_ext, &b_ext_trsf) != 1) {
            VT_FreeImageBarycentres(&b_in);
            VT_FreeImageBarycentres(&b_in_trsf);
            VT_FreeImageBarycentres(&b_ext);
            VT_FreeImageBarycentres(&b_ext_trsf);
            vtfree(residus); residus=NULL;
            vtfree(std_deviations); std_deviations=NULL;
            vtfree(thetas);  thetas=NULL;
          MT_ErrorParse("unable to process the image of barycentres transformation\n",0);
        }
        if(_verbose_>5) {
            fprintf(stdout, "\nBarycentres ref :\n");
            VT_PrintImageBarycentres(b_in);
            fprintf(stdout, "\nBarycentres flo :\n");
            VT_PrintImageBarycentres(b_ext);
            fprintf(stdout, "\nBarycentres flo trsf :\n");
            VT_PrintImageBarycentres(b_ext_trsf);
        }

        if (VT_ComputePairsOfPoints(b_in, b_ext_trsf, pairs, &npairs, &r) != 1)
        {
            VT_FreeImageBarycentres(&b_in);
            VT_FreeImageBarycentres(&b_in_trsf);
            VT_FreeImageBarycentres(&b_ext);
            VT_FreeImageBarycentres(&b_ext_trsf);
            vtfree(thetas); thetas=NULL;
            vtfree(residus); residus=NULL;
            vtfree(std_deviations); std_deviations=NULL;
            MT_ErrorParse("unable to compute the pairs of points to be associated\n",0);
        }

        if (_verbose_ > 4) {
          fprintf(stdout, "npairs = %d ; r = %lf\n", npairs, r);
          if (_verbose_ > 5) {
            int i_tmp;
            for (i_tmp=0 ; i_tmp<npairs ; i_tmp++)
              fprintf(stdout, "( %d, %d )\n", pairs[0][i_tmp], pairs[1][i_tmp]);
          }
        }


        /* CALCUL DE T_tmp = TRANSFO DE IN_THETA VERS EXT */
        if( VT_ComputePointTrsf(pairs, npairs, b_in, b_ext, &r, &sd, T_tmp, par.flag_affine, par.estimator) != 1) {
            VT_FreeImageBarycentres(&b_in);
            VT_FreeImageBarycentres(&b_in_trsf);
            VT_FreeImageBarycentres(&b_ext);
            VT_FreeImageBarycentres(&b_ext_trsf);
            vtfree(residus); residus=NULL;
            vtfree(std_deviations); std_deviations=NULL;
            vtfree(thetas);  thetas=NULL;
            vtfree(pairs[0]); pairs[0]=NULL;
            vtfree(pairs[1]); pairs[1]=NULL;
            MT_ErrorParse("unable to compute the points to points transformation...\n", 0);
        }


        if(_verbose_>2)
        {

            fprintf(stdout, "i=%d\ttheta=%f\n", i, theta);
            fprintf(stdout, "residu = %f\n", r);

            fprintf(stdout, "Nombre d'appariements: %d\n", npairs);

            if (_verbose_>3) {
              if (i%2==0) {
		/*            fprintf(stdout, "#Transformation matrix (coordonnees de ext vers theta = %f * PI):\n", theta/PI); */
                fprintf(stdout, "#Transformation matrix (ref<-flo, theta = %f * PI ):\n", theta/PI);
              }
              else {
                fprintf(stdout, "#Transformation matrix (ref<-flo_S, theta = %f * PI ):\n", theta/PI);
              }
              fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_in_theta_ext_sym[0],T_in_theta_ext_sym[1],T_in_theta_ext_sym[2],T_in_theta_ext_sym[3]);
              fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_in_theta_ext_sym[4],T_in_theta_ext_sym[5],T_in_theta_ext_sym[6],T_in_theta_ext_sym[7]);
              fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_in_theta_ext_sym[8],T_in_theta_ext_sym[9],T_in_theta_ext_sym[10],T_in_theta_ext_sym[11]);
              fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_in_theta_ext_sym[12],T_in_theta_ext_sym[13],T_in_theta_ext_sym[14],T_in_theta_ext_sym[15]);
            }

            if (_verbose_>3) {
              if (i%2==0) {
		/*            fprintf(stdout, "#Transformation matrix (coordonnees de ext vers theta = %f * PI):\n", theta/PI); */
                fprintf(stdout, "#Transformation matrix (flo<-ref, theta = %f * PI ):\n", theta/PI);
              }
              else {
                fprintf(stdout, "#Transformation matrix (flo_S<-ref, theta = %f * PI ):\n", theta/PI);
              }
              fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_tmp[0],T_tmp[1],T_tmp[2],T_tmp[3]);
              fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_tmp[4],T_tmp[5],T_tmp[6],T_tmp[7]);
              fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_tmp[8],T_tmp[9],T_tmp[10],T_tmp[11]);
              fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_tmp[12],T_tmp[13],T_tmp[14],T_tmp[15]);
            }

            /*if (_verbose_>3) { */
            /*    if (i%2==0) fprintf(stdout, "Appariements (theta = %f * PI):", theta/PI); */
            /*    else fprintf(stdout, "Appariements (retourne, theta = %f * PI):", theta/PI); */
            /*    for (j=0; j<npairs  ; j++) { */
            /*        if(j%2==0) fprintf(stdout, "\n"); */
            /*        theB_in=b_in.barycentres[pairs[0][j]]; */
            /*        theB_ext=b_ext.barycentres[pairs[1][j]]; */
            /*        fprintf(stdout, "( %d  %d ) \t", theB_in.label, theB_ext.label); */
            /*    } */
            /*    fprintf(stdout, "\n"); */
            /*} */
            fprintf(stdout, "\n");

        }
        residus[i]=r;
        std_deviations[i]=sd;

        if((((r>=0 && r<rmin) || rmin<0 ) && 1)||0) {
            int j;
            for (j=0; j<16; j++) {
                T_ext_in[j]=T_tmp[j];
                T_in_ext_min[j]=T_in_theta_ext_sym[j];
            }
            rmin=r;
            sdmin=sd;
            imin=i;
            /*fprintf(stderr, "test tata\n");*/

            if(par.name_res.in[0] != '\0')
            {
                int j;
                /*fprintf(stdout, "npairs = %d\n", npairs); */
                FILE* fichier = NULL;
                fprintf(stderr, "par.name_res.in : %s\n", par.name_res.in);
                fichier = fopen(par.name_res.in, "w");
                if( fichier == NULL)
                {
                    VT_FreeImageBarycentres(&b_in);
                    VT_FreeImageBarycentres(&b_in_trsf);
                    VT_FreeImageBarycentres(&b_ext);
                    VT_FreeImageBarycentres(&b_ext_trsf);
                    vtfree(residus); residus=NULL;
                    vtfree(thetas);  thetas=NULL;
                    vtfree(pairs[0]); pairs[0]=NULL;
                    vtfree(pairs[1]); pairs[1]=NULL;
                    MT_ErrorParse("Erreur pour l'ecriture de la liste d'appariements dans un fichier\n", 0);
                }
                fprintf(stderr, "file opened\n");

                fprintf(fichier, "#Liste d'appariements (ref-flo) trouves:\n#\n");
                vt_barycentre theB_in;
                vt_barycentre theB_ext;
                for (j=0; j<npairs; j++) {
                    theB_in=b_in.barycentres[pairs[0][j]];
                    theB_ext=b_ext.barycentres[pairs[1][j]];
                    fprintf(fichier, "%d \t%d\n", theB_in.label, theB_ext.label);
                }
                fclose(fichier);
            }

        }

        vtfree(pairs[0]); pairs[0]=NULL;
        vtfree(pairs[1]); pairs[1]=NULL;


      }
      if (_verbose_>1) {
        fprintf(stdout, "#Transformation matrix min (ref<-flo) before point Transformation computation:\n");
        fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_in_ext_min[0],T_in_ext_min[1],T_in_ext_min[2],T_in_ext_min[3]);
        fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_in_ext_min[4],T_in_ext_min[5],T_in_ext_min[6],T_in_ext_min[7]);
        fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_in_ext_min[8],T_in_ext_min[9],T_in_ext_min[10],T_in_ext_min[11]);
        fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_in_ext_min[12],T_in_ext_min[13],T_in_ext_min[14],T_in_ext_min[15]);
      }

      /*if(_verbose_>2) {
          int j;
          vt_barycentre theB_in, theB_ext;
          fprintf(stdout, "%d appariements utilises pour le calcul de la transformation :\n", npairs);
          for (j=0; j<npairs  ; j++) {
            theB_in=b_in.barycentres[pairs[0][j]];
            theB_ext=b_ext.barycentres[pairs[1][j]];
            fprintf(stdout, "# %d / %d : ( %d  %d )\n", j+1, npairs, theB_in.label, theB_ext.label);
          }
      }*/


  }
  else {
    /* CALCUL DE T_ext_in = TRANSFO DE IN_THETA VERS EXT */
      pairs[0] = vtmalloc(par.npairs*sizeof(int), "pairs[0]", argv[0] );
      if(pairs[0]==NULL) {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext);
          MT_ErrorParse("unable to allocate a vector...\n", 0);
      }
      pairs[1] = vtmalloc(par.npairs*sizeof(int), "pairs[1]", argv[0] );
      if(pairs[1]==NULL) {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext);
          vtfree(pairs[0]); pairs[0]=NULL;
          MT_ErrorParse("unable to allocate a vector...\n", 0);
      }
      for (i=0;i<par.npairs;i++) {
          pairs[0][i]=VT_FindLabel(b_in, par.pairs[0][i] );
          pairs[1][i]=VT_FindLabel(b_ext, par.pairs[1][i]);
          if (pairs[0][i]<0 || pairs[1][i]<0 )
          {
              VT_FreeImageBarycentres(&b_in);
              VT_FreeImageBarycentres(&b_in_trsf);
              VT_FreeImageBarycentres(&b_ext);
              vtfree(pairs[0]); pairs[0]=NULL;
              vtfree(pairs[1]); pairs[1]=NULL;
              MT_ErrorParse("label not find, exiting the process...\n", 0);
          }
      }
      /* GM: _VT_VERBOSE_ is a boolean
       * */
      if (_VT_VERBOSE_ ) fprintf(stdout, "npairs=%d\nPairs : \n", par.npairs);
      if (_VT_VERBOSE_ )
        for (i=0;i<par.npairs;i++) {
          fprintf(stdout, "( %d \t%d )\n", par.pairs[0][i], par.pairs[1][i]);
        }
      if( VT_ComputePointTrsf(pairs, par.npairs, b_in, b_ext, &r, &sd, T_ext_in, par.flag_affine, par.estimator) != 1) {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext);
          MT_ErrorParse("unable to compute the points to points transformation...\n", 0);
          vtfree(pairs[0]); pairs[0]=NULL;
          vtfree(pairs[1]); pairs[1]=NULL;
      }

      /*fprintf(stderr, "test toto\n");*/
      if(par.name_res.in[0] != '\0')
      {
	/*fprintf(stdout, "npairs = %d\n", npairs); */
          FILE* fichier = NULL;
          fprintf(stderr, "par.name_res.in : %s\n", par.name_res.in);
          fichier = fopen(par.name_res.in, "w");
          if( fichier == NULL)
          {
              VT_FreeImageBarycentres(&b_in);
              VT_FreeImageBarycentres(&b_in_trsf);
              VT_FreeImageBarycentres(&b_ext);
              MT_ErrorParse("unable to compute the points to points transformation...\n", 0);
              vtfree(pairs[0]); pairs[0]=NULL;
              vtfree(pairs[1]); pairs[1]=NULL;
              MT_ErrorParse("Erreur pour l'ecriture de la liste d'appariements dans un fichier\n", 0);
          }
          fprintf(stderr, "file opened\n");
          fprintf(fichier, "#Liste d'appariements (ref-flo) entres :\n#\n");
          vt_barycentre theB_in;
          vt_barycentre theB_ext;
          for (i=0; i<par.npairs; i++) {
              theB_in=b_in.barycentres[pairs[0][i]];
              theB_ext=b_ext.barycentres[pairs[1][i]];
              fprintf(fichier, "%d \t%d\n", theB_in.label, theB_ext.label);
          }
          fclose(fichier);
      }

      vtfree(pairs[0]); pairs[0]=NULL;
      vtfree(pairs[1]); pairs[1]=NULL;
  }


  if(_VT_VERBOSE_) {
      fprintf(stdout, "Resultat de l'appariement :\n");
      if(_verbose_ > 1) {
        if (par.npairs==0){
            fprintf(stdout, "Residus:\n[  ");
            for (i=0; i<2*n; i++) fprintf(stdout, "%f  ", residus[i]);
            fprintf(stdout, "]\n");
            fprintf(stdout, "Ecart-type des Residus:\n[  ");
            for (i=0; i<2*n; i++) fprintf(stdout, "%f  ", std_deviations[i]);
            fprintf(stdout, "]\n");
        }
      }
      if(par.npairs==0) {
          fprintf(stdout, "imin = %d\n", imin);
          if (imin%2==0) fprintf(stdout, "#Transformation matrix (flo<-ref, obtained with theta = %f * PI ):\n", thetas[imin]/PI);
          /*else fprintf(stdout, "#Transformation matrix (ext_R<-in, theta = %f * PI ):\n", thetas[imin]/PI); */
      }
      else fprintf(stdout, "#Transformation matrix (flo<-ref):\n");
      fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_ext_in[0],T_ext_in[1],T_ext_in[2],T_ext_in[3]);
      fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_ext_in[4],T_ext_in[5],T_ext_in[6],T_ext_in[7]);
      fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_ext_in[8],T_ext_in[9],T_ext_in[10],T_ext_in[11]);
      fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_ext_in[12],T_ext_in[13],T_ext_in[14],T_ext_in[15]);
      if (par.npairs==0) fprintf(stdout, "residu = %f\nstd deviation = %f\n", rmin, sdmin);
      else fprintf(stdout, "residu = %f\nstd deviation = %f\n", r, sd);
  }
  if(par.name_residuals.out[0] != '\0')
  {
      FILE* fichier = NULL;
      fichier = fopen(par.name_residuals.out, "w");
      if( fichier == NULL)
      {
        MT_ErrorParse("Erreur pour l'ecriture du resultat dans un fichier\n", 0);
      }
      if (par.npairs==0){
          fprintf(fichier, "%%Residus:\nresiduals=[  ");
          for (i=0; i<2*n; i++) fprintf(fichier, "%f  ", residus[i]);
          fprintf(fichier, "];\n");
          fprintf(fichier, "%%Ecart-type des Residus:\nstd_dev=[  ");
          for (i=0; i<2*n; i++) fprintf(fichier, "%f  ", std_deviations[i]);
          fprintf(fichier, "];\n");
          fprintf(fichier, "%%Iter min:\nimin=%d ;\n", imin+1); /* +1 pour matlab. */
      }
      else {
          fprintf(fichier, "%%Residu:\nresidual=%f ;\n", rmin);
          fprintf(fichier, "%%Ecart-type du Residu:\nstd_dev=%f ;\n", sdmin);
      }

      fclose(fichier);
  }
  if(par.name_res.out[0] != '\0')
  {
      FILE* fichier = NULL;
      fichier = fopen(par.name_res.out, "w");
      if( fichier == NULL)
      {
        MT_ErrorParse("Erreur pour l'ecriture du resultat dans un fichier\n", 0);
      }
      fprintf(fichier, "#\n# \tTransformation matrix (flo<-ref):\n#\n");
      fprintf(fichier, "%f \t%f \t%f \t%f\n%f \t%f \t%f \t%f\n%f \t%f \t%f \t%f\n%f \t%f \t%f \t%f\n",
              T_ext_in[0],T_ext_in[1],T_ext_in[2],T_ext_in[3],
              T_ext_in[4],T_ext_in[5],T_ext_in[6],T_ext_in[7],
              T_ext_in[8],T_ext_in[9],T_ext_in[10],T_ext_in[11],
              T_ext_in[12],T_ext_in[13],T_ext_in[14],T_ext_in[15]);
      if (par.npairs==0) fprintf(fichier, "#Residu = %f\n#Std deviation = %f\n", rmin, sdmin);
      else fprintf(fichier, "#Residu = %f\n#Std deviation = %f\n", r, sd);
      fprintf(fichier, "#\n");
      fprintf(fichier, "\n");
      fclose(fichier);
  }

  if (par.name_res.ext[0] != '\0') {

      FILE* fichier = NULL;
      fichier = fopen(par.name_res.ext, "w");
      if( fichier == NULL)
      {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_ext);
          VT_FreeImageBarycentres(&b_in_trsf);
          if(par.npairs==0) {
              VT_FreeImageBarycentres(&b_ext_trsf);
              vtfree(residus); residus=NULL;
              vtfree(std_deviations); std_deviations=NULL;
              vtfree(thetas);  thetas=NULL;
          }

          MT_ErrorParse("Erreur pour l'ecriture de la liste d'appariements dans un fichier\n", 0);
      }

      if (VT_TransformImageBarycentres(T_ext_in, b_in, &b_in_trsf) != 1) {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext);
          if(par.npairs==0) {
              VT_FreeImageBarycentres(&b_ext_trsf);
              vtfree(residus); residus=NULL;
              vtfree(std_deviations); std_deviations=NULL;
              vtfree(thetas);  thetas=NULL;
          }
        MT_ErrorParse("unable to process the image of barycentres transformation\n",0);
      }

      if (VT_ComputePairsOfPoints(b_in_trsf, b_ext, pairs, &npairs, &r) != 1)
      {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext);
          if(par.npairs==0) {
              VT_FreeImageBarycentres(&b_ext_trsf);
              vtfree(thetas); thetas=NULL;
              vtfree(residus); residus=NULL;
              vtfree(std_deviations); std_deviations=NULL;
          }
          MT_ErrorParse("unable to compute the pairs of points to be associated\n",0);
      }

      fprintf(fichier, "#Liste d'appariements (ref-flo) recalcules :\n#\n");
      vt_barycentre theB_in;
      vt_barycentre theB_ext;



      for (i=0; i<npairs; i++) {
          theB_in=b_in.barycentres[pairs[0][i]];
          theB_ext=b_ext.barycentres[pairs[1][i]];
          fprintf(fichier, "%d \t%d\n", theB_in.label, theB_ext.label);
      }


      /*fprintf(fichier, "#Residu correspondant: %f\n", r); */
      fclose(fichier);
      vtfree(pairs[0]); pairs[0]=NULL;
      vtfree(pairs[1]); pairs[1]=NULL;

  }




  if (par.npairs == 0 && par.dices_in == 1)
  {
      FILE* fichier = NULL;
      fichier = fopen(par.name_dices.in, "w");
      double *dices;

      if( fichier == NULL)
      {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_ext);
          VT_FreeImageBarycentres(&b_in_trsf);
          if(par.npairs==0) {
              VT_FreeImageBarycentres(&b_ext_trsf);
              vtfree(residus); residus=NULL;
              vtfree(std_deviations); std_deviations=NULL;
              vtfree(thetas);  thetas=NULL;
          }

          MT_ErrorParse("Erreur pour l'ecriture de la liste des Dices dans un fichier\n", 0);
      }

      if (VT_TransformImageBarycentres(T_in_ext_min, b_ext, &b_ext_trsf) != 1) {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext);
          VT_FreeImageBarycentres(&b_ext_trsf);
          vtfree(residus); residus=NULL;
          vtfree(std_deviations); std_deviations=NULL;
          vtfree(thetas);  thetas=NULL;
          MT_ErrorParse("unable to process the image of barycentres transformation\n",0);
      }

      if (VT_ComputePairsOfPoints(b_in, b_ext_trsf, pairs, &npairs, &r) != 1)
      {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext);
          VT_FreeImageBarycentres(&b_ext_trsf);
          vtfree(thetas); thetas=NULL;
          vtfree(residus); residus=NULL;
          vtfree(std_deviations); std_deviations=NULL;
          MT_ErrorParse("unable to compute the pairs of points to be associated\n",0);
      }

      /*--- lecture des images labels d'entree ---*/
      images_labels[0] = _VT_Inrimage( par.name_labels.in );
      if ( images_labels[0] == (vt_image*)NULL ) {
        vtfree(pairs[0]); pairs[0]=NULL;
        vtfree(pairs[1]); pairs[1]=NULL;
        VT_FreeImageBarycentres(&b_in);
        VT_FreeImageBarycentres(&b_ext);
        VT_FreeImageBarycentres(&b_in_trsf);
        VT_FreeImageBarycentres(&b_ext_trsf);
        vtfree(residus); residus=NULL;
        vtfree(std_deviations); std_deviations=NULL;
        vtfree(thetas);  thetas=NULL;
        MT_ErrorParse("unable to read input image labels ref\n", 0);
      }
      images_labels[1] = _VT_Inrimage( par.name_labels.ext );
      if ( images_labels[1] == (vt_image*)NULL ) {
        VT_FreeImage(images_labels[0]);
        vtfree(pairs[0]); pairs[0]=NULL;
        vtfree(pairs[1]); pairs[1]=NULL;
        VT_FreeImageBarycentres(&b_in);
        VT_FreeImageBarycentres(&b_ext);
        VT_FreeImageBarycentres(&b_in_trsf);
        VT_FreeImageBarycentres(&b_ext_trsf);
        vtfree(residus); residus=NULL;
        vtfree(std_deviations); std_deviations=NULL;
        vtfree(thetas);  thetas=NULL;
        MT_ErrorParse("unable to read input image labels flo\n", 0);
      }

      /*--- operations eventuelles sur l'image d'entree ---*/
      if ( par.name_labels.inv == 1 )  {VT_InverseImage( images_labels[0]); VT_InverseImage( images_labels[1]); }
      if ( par.name_labels.swap == 1 ) {VT_SwapImage( images_labels[0]); VT_SwapImage( images_labels[1]); }



      dices = vtmalloc(npairs*sizeof(double), "dices", argv[0] );
      if(dices==NULL)
      {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_ext);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext_trsf);
          vtfree(residus); residus=NULL;
          vtfree(std_deviations); std_deviations=NULL;
          vtfree(thetas);  thetas=NULL;
          VT_FreeImage(images_labels[0]);
          VT_FreeImage(images_labels[1]);
          vtfree(pairs[0]); pairs[0]=NULL;
          vtfree(pairs[1]); pairs[1]=NULL;
          MT_ErrorParse("unable to allocate Dices vector\n", 0);
      }

      int *labels[2];
      labels[0] = vtmalloc(npairs*sizeof(int), "labels[0]", argv[0] );
      if(dices==NULL)
      {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_ext);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext_trsf);
          vtfree(residus); residus=NULL;
          vtfree(std_deviations); std_deviations=NULL;
          vtfree(thetas);  thetas=NULL;
          VT_FreeImage(images_labels[0]);
          VT_FreeImage(images_labels[1]);
          vtfree(pairs[0]); pairs[0]=NULL;
          vtfree(pairs[1]); pairs[1]=NULL;
          vtfree(dices);
          MT_ErrorParse("unable to allocate Dices vector\n", 0);
      }
      labels[1] = vtmalloc(npairs*sizeof(int), "labels[1]", argv[0] );
      if(dices==NULL)
      {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_ext);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext_trsf);
          vtfree(residus); residus=NULL;
          vtfree(std_deviations); std_deviations=NULL;
          vtfree(thetas);  thetas=NULL;
          VT_FreeImage(images_labels[0]);
          VT_FreeImage(images_labels[1]);
          vtfree(pairs[0]); pairs[0]=NULL;
          vtfree(pairs[1]); pairs[1]=NULL;
          vtfree(dices);
          vtfree(labels[0]);
          MT_ErrorParse("unable to allocate Dices vector\n", 0);
      }

      vt_barycentre theB_in;
      vt_barycentre theB_ext;

      for (i=0; i<npairs; i++) {
          theB_in=b_in.barycentres[pairs[0][i]];
          theB_ext=b_ext.barycentres[pairs[1][i]];
          labels[0][i]=theB_in.label;
          labels[1][i]=theB_ext.label;
      }


      VT_ComputePairDices(images_labels[1], images_labels[0], T_in_ext_min, labels[1], labels[0], npairs, dices);


      fprintf(fichier, "#Labels Ref \tLabels Flo \tDices\n");
      double mean_dices=0;
      double std_dev_dices=0;
      for (i=0; i<npairs; i++) {
          fprintf(fichier, "%d \t%d \t%f\n", labels[0][i], labels[1][i], dices[i]);
          mean_dices+=dices[i];
          std_dev_dices+=dices[i]*dices[i];
      }
      mean_dices/=npairs;
      std_dev_dices=sqrt((std_dev_dices)/npairs-mean_dices*mean_dices);
      fprintf(fichier, "#Mean Dice: %f\n", mean_dices);
      fprintf(fichier, "#Standard deviation: %f\n", std_dev_dices);


      fclose(fichier);
      VT_FreeImage(images_labels[0]);
      VT_FreeImage(images_labels[1]);
      vtfree(pairs[0]); pairs[0]=NULL;
      vtfree(pairs[1]); pairs[1]=NULL;
      vtfree(dices);
      vtfree(labels[0]);
      vtfree(labels[1]);
  }




  if (par.dices_out==1)
  {
      FILE* fichier = NULL;
      fichier = fopen(par.name_dices.out, "w");
      double *dices;

      if( fichier == NULL)
      {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_ext);
          VT_FreeImageBarycentres(&b_in_trsf);
          if(par.npairs==0) {
              VT_FreeImageBarycentres(&b_ext_trsf);
              vtfree(residus); residus=NULL;
              vtfree(std_deviations); std_deviations=NULL;
              vtfree(thetas);  thetas=NULL;
          }

          MT_ErrorParse("Erreur pour l'ecriture de la liste des Dices dans un fichier\n", 0);
      }

      if (VT_TransformImageBarycentres(T_ext_in, b_in, &b_in_trsf) != 1) {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext);
          if(par.npairs==0) {
              VT_FreeImageBarycentres(&b_ext_trsf);
              vtfree(residus); residus=NULL;
              vtfree(std_deviations); std_deviations=NULL;
              vtfree(thetas);  thetas=NULL;
          }
        MT_ErrorParse("unable to process the image of barycentres transformation\n",0);
      }

      if (VT_ComputePairsOfPoints(b_in_trsf, b_ext, pairs, &npairs, &r) != 1)
      {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext);
          if(par.npairs==0) {
              VT_FreeImageBarycentres(&b_ext_trsf);
              vtfree(thetas); thetas=NULL;
              vtfree(residus); residus=NULL;
              vtfree(std_deviations); std_deviations=NULL;
          }
          MT_ErrorParse("unable to compute the pairs of points to be associated\n",0);
      }

      /*--- lecture des images labels d'entree ---*/
      images_labels[0] = _VT_Inrimage( par.name_labels.in );
      if ( images_labels[0] == (vt_image*)NULL ) {
        vtfree(pairs[0]); pairs[0]=NULL;
        vtfree(pairs[1]); pairs[1]=NULL;
        VT_FreeImageBarycentres(&b_in);
        VT_FreeImageBarycentres(&b_ext);
        VT_FreeImageBarycentres(&b_in_trsf);
        if(par.npairs==0) {
            VT_FreeImageBarycentres(&b_ext_trsf);
            vtfree(residus); residus=NULL;
            vtfree(std_deviations); std_deviations=NULL;
            vtfree(thetas);  thetas=NULL;
        }
        MT_ErrorParse("unable to read input image labels ref\n", 0);
      }
      images_labels[1] = _VT_Inrimage( par.name_labels.ext );
      if ( images_labels[1] == (vt_image*)NULL ) {
        VT_FreeImage(images_labels[0]);
        vtfree(pairs[0]); pairs[0]=NULL;
        vtfree(pairs[1]); pairs[1]=NULL;
        VT_FreeImageBarycentres(&b_in);
        VT_FreeImageBarycentres(&b_ext);
        VT_FreeImageBarycentres(&b_in_trsf);
        if(par.npairs==0) {
            VT_FreeImageBarycentres(&b_ext_trsf);
            vtfree(residus); residus=NULL;
            vtfree(std_deviations); std_deviations=NULL;
            vtfree(thetas);  thetas=NULL;
        }
        MT_ErrorParse("unable to read input image labels flo\n", 0);
      }

      /*--- operations eventuelles sur l'image d'entree ---*/
      if ( par.name_labels.inv == 1 )  {VT_InverseImage( images_labels[0]); VT_InverseImage( images_labels[1]); }
      if ( par.name_labels.swap == 1 ) {VT_SwapImage( images_labels[0]); VT_SwapImage( images_labels[1]); }



      dices = vtmalloc(npairs*sizeof(double), "dices", argv[0] );
      if(dices==NULL)
      {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_ext);
          VT_FreeImageBarycentres(&b_in_trsf);
          if(par.npairs==0) {
              VT_FreeImageBarycentres(&b_ext_trsf);
              vtfree(residus); residus=NULL;
              vtfree(std_deviations); std_deviations=NULL;
              vtfree(thetas);  thetas=NULL;
          }
          VT_FreeImage(images_labels[0]);
          VT_FreeImage(images_labels[1]);
          vtfree(pairs[0]); pairs[0]=NULL;
          vtfree(pairs[1]); pairs[1]=NULL;
          MT_ErrorParse("unable to allocate Dices vector\n", 0);
      }

      int *labels[2];
      labels[0] = vtmalloc(npairs*sizeof(int), "labels[0]", argv[0] );
      if(dices==NULL)
      {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_ext);
          VT_FreeImageBarycentres(&b_in_trsf);
          if(par.npairs==0) {
              VT_FreeImageBarycentres(&b_ext_trsf);
              vtfree(residus); residus=NULL;
              vtfree(std_deviations); std_deviations=NULL;
              vtfree(thetas);  thetas=NULL;
          }
          VT_FreeImage(images_labels[0]);
          VT_FreeImage(images_labels[1]);
          vtfree(pairs[0]); pairs[0]=NULL;
          vtfree(pairs[1]); pairs[1]=NULL;
          vtfree(dices);
          MT_ErrorParse("unable to allocate Dices vector\n", 0);
      }
      labels[1] = vtmalloc(npairs*sizeof(int), "labels[1]", argv[0] );
      if(dices==NULL)
      {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_ext);
          VT_FreeImageBarycentres(&b_in_trsf);
          if(par.npairs==0) {
              VT_FreeImageBarycentres(&b_ext_trsf);
              vtfree(residus); residus=NULL;
              vtfree(std_deviations); std_deviations=NULL;
              vtfree(thetas);  thetas=NULL;
          }
          VT_FreeImage(images_labels[0]);
          VT_FreeImage(images_labels[1]);
          vtfree(pairs[0]); pairs[0]=NULL;
          vtfree(pairs[1]); pairs[1]=NULL;
          vtfree(dices);
          vtfree(labels[0]);
          MT_ErrorParse("unable to allocate Dices vector\n", 0);
      }

      vt_barycentre theB_in;
      vt_barycentre theB_ext;

      for (i=0; i<npairs; i++) {
          theB_in=b_in.barycentres[pairs[0][i]];
          theB_ext=b_ext.barycentres[pairs[1][i]];
          labels[0][i]=theB_in.label;
          labels[1][i]=theB_ext.label;
      }


      VT_ComputePairDices(images_labels[0], images_labels[1], T_ext_in, labels[0], labels[1], npairs, dices);


      fprintf(fichier, "#Labels Ref \tLabels Flo \tDices\n");
      double mean_dices=0;
      double std_dev_dices=0;
      for (i=0; i<npairs; i++) {
          fprintf(fichier, "%d \t%d \t%f\n", labels[0][i], labels[1][i], dices[i]);
          mean_dices+=dices[i];
          std_dev_dices+=dices[i]*dices[i];
      }
      mean_dices/=npairs;
      std_dev_dices=sqrt((std_dev_dices)/npairs-mean_dices*mean_dices);
      fprintf(fichier, "#Mean Dice: %f\n", mean_dices);
      fprintf(fichier, "#Standard deviation: %f\n", std_dev_dices);

      fclose(fichier);
      VT_FreeImage(images_labels[0]);
      VT_FreeImage(images_labels[1]);
      vtfree(pairs[0]); pairs[0]=NULL;
      vtfree(pairs[1]); pairs[1]=NULL;
      vtfree(dices);
      vtfree(labels[0]);
      vtfree(labels[1]);
  }

  /*--- liberations memoires ---*/
  VT_FreeImageBarycentres(&b_in);
  VT_FreeImageBarycentres(&b_ext);
  VT_FreeImageBarycentres(&b_in_trsf);
  if(par.npairs==0) {
      VT_FreeImageBarycentres(&b_ext_trsf);
      vtfree(residus); residus=NULL;
      vtfree(std_deviations); std_deviations=NULL;
      vtfree(thetas);  thetas=NULL;
  }

  return( 0 );
}








static void MT_Parse( int argc,
                      char *argv[],
                      local_par *par )
{
  int i;
  char text[STRINGLENGTH];
  int status;

  if ( VT_CopyName( program, argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);
  if ( argc == 1 ) MT_ErrorParse("\n", 0 );

  /*--- lecture des parametres ---*/
  i = 1;
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
          par->name_gradients.inv = 1;
          par->name_labels.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
          par->name_gradients.swap = 1;
          par->name_labels.swap = 1;
      }


      /* Parametres de calcul */

      else if ( strcmp ( argv[i], "-dices" ) == 0 || strcmp ( argv[i], "-dices-out" ) == 0 ){
        par->dices_out=1;
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -dices-out...\n",0);
        sprintf(par->name_dices.out, "%s", argv[i]);
      }

      else if ( strcmp ( argv[i], "-dices-in" ) == 0 ){
        par->dices_in=1;
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -dices-in...\n",0);
        sprintf(par->name_dices.in, "%s", argv[i]);
      }

      else if ( strcmp ( argv[i], "-gradient-in" ) == 0 || strcmp ( argv[i], "-gradient-ref" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -gradient-ref...\n",0);
        sprintf(par->name_gradients.in, "%s", argv[i]);
      }

      else if ( strcmp ( argv[i], "-gradient-ext" ) == 0 || strcmp ( argv[i], "-gradient-flo" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -gradient-flo...\n",0);
        sprintf(par->name_gradients.ext, "%s", argv[i]);
      }
      else if ( strcmp ( argv[i], "-label-in" ) == 0 || strcmp ( argv[i], "-label-ref" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -label-ref...\n",0);
        sprintf(par->name_labels.in, "%s", argv[i]);
      }
      else if ( strcmp ( argv[i], "-label-ext" ) == 0 || strcmp ( argv[i], "-label-flo" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -label-flo...\n",0);
        sprintf(par->name_labels.ext, "%s", argv[i]);
      }


      else if ( strcmp ( argv[i], "-delta" ) == 0 ) {
          i += 1;
          if ( i >= argc)    MT_ErrorParse( "parsing -delta...\n", 0 );
          status = sscanf( argv[i],"%lf",&(par->delta) );
          if ( status <= 0 ) MT_ErrorParse( "parsing -delta...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-trsf" ) == 0 || strcmp ( argv[i], "-file-out" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -trsf...\n",0);
        sprintf(par->name_res.out, "%s", argv[i]);
      }

      else if ( strcmp ( argv[i], "-pairs-in" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -pairs-in...\n",0);
        sprintf(par->name_res.in, "%s", argv[i]);
      }

      else if ( strcmp ( argv[i], "-pairs-out" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -pairs-out...\n",0);
        sprintf(par->name_res.ext, "%s", argv[i]);
      }

      else if ( strcmp ( argv[i], "-pairs" ) == 0 ||
                strcmp ( argv[i], "-pair" ) == 0  ||
                strcmp ( argv[i], "-p" ) == 0){
          i += 1;
          if ( i >= argc)    MT_ErrorParse( "parsing -p...\n", 0 );

          FILE *f;
          if ((f = fopen (argv[i], "r")) != NULL) {
	    /* File */
              char line[512];
              int f1, f2;
              int n=0;
              while ( fgets(line, 512, f) != NULL ) {
                  if ( ( sscanf( line, "%d %d\n", &f1, &f2) == 2 ) ) {
                      n+=1;
                  }
              }

              if ( n==0 ) {
                  MT_ErrorParse( "parsing -p... incompatible file.\n", 0);
              }
              if (n >= TABLENGTH)
              {
                  MT_ErrorParse( "parsing -p... too many pairs.\n", 0);
              }
              fclose( f );

              if ((f = fopen (argv[i], "r")) == NULL) {
                  MT_ErrorParse( "parsing -p... file disappeared.\n", 0);
              }


              if ( fgets(line, 512, f) == NULL ) {
                  MT_ErrorParse( "parsing -p... an error occured while reading file.\n", 0);
              }

              while ( ( sscanf( line, "%d %d\n", &f1, &f2 ) < 2 )  ) {
                if ( fgets(line, 512, f) == NULL ) {
                    MT_ErrorParse( "parsing -p... an error occured while reading file.\n", 0);
                }
              }
              par->pairs[0][0] = f1;
              par->pairs[1][0] = f2;
              par->npairs = 1;

              /* get the n-1 other rows
               */
              while (par->npairs<n)  {
                if ( fgets(line, 512, f) == NULL ) {
                    MT_ErrorParse( "parsing -p... an error occured while reading file.\n", 0);
                }
                if ( sscanf( line, "%d %d\n", &f1, &f2 ) != 2 ) {
                    MT_ErrorParse( "parsing -p... an error occured while reading file.\n", 0);
                }
                par->pairs[0][par->npairs] = f1;
                par->pairs[1][par->npairs] = f2;
                par->npairs++;
              }

              fclose( f );
          }
          else {
	    /*Not File */
              status = sscanf( argv[i],"%d",&(par->pairs[0][par->npairs]) );
              if ( status <= 0 ) MT_ErrorParse( "parsing -p...\n", 0 );
              i += 1;
              if ( i >= argc)    MT_ErrorParse( "parsing -p...\n", 0 );
              status = sscanf( argv[i],"%d",&(par->pairs[1][par->npairs]) );
              if ( status <= 0 ) MT_ErrorParse( "parsing -p...\n", 0 );
              par->npairs += 1;
              i += 1;
              while ( i < argc)  {
                if (par->npairs>=TABLENGTH)
                    MT_ErrorParse( "parsing -p : too many pairs...\n", 0 );
                status = sscanf( argv[i],"%d",&(par->pairs[0][par->npairs]) );
                if ( status <= 0 ) {i--;  break;}
                i += 1;
                if ( i >= argc)    MT_ErrorParse( "parsing -p...\n", 0 );
                status = sscanf( argv[i],"%d",&(par->pairs[1][par->npairs]) );
                if ( status <= 0 ) MT_ErrorParse( "parsing -p...\n", 0 );
                par->npairs += 1;
                i += 1;
              }
          }
      }

      else if ( strcmp ( argv[i], "-residuals" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -residuals...\n",0);
        sprintf(par->name_residuals.out, "%s", argv[i]);
      }



      else if ( strcmp ( argv[i], "-p-in-file" ) == 0 || strcmp ( argv[i], "-p-ref-file" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -p-ref-file...\n",0);
        sprintf(par->name_planes.in, "%s", argv[i]);
        par->flag_name_planes_in=1;
      }
      else if ( strcmp ( argv[i], "-p-ext-file" ) == 0 || strcmp ( argv[i], "-p-flo-file" ) == 0 ){
        i += 1;
        if ( i>= argc ) MT_ErrorParse( "parsing -p-flo-file...\n",0);
        sprintf(par->name_planes.ext, "%s", argv[i]);
        par->flag_name_planes_ext=1;
      }

      else if ( strcmp ( argv[i], "-p-in" ) == 0 || strcmp ( argv[i], "-p-ref" ) == 0 ) {
          i += 1;
          if ( i >= argc)    MT_ErrorParse( "parsing -p-ref...\n", 0 );
          status = sscanf( argv[i],"%lf",&(par->p_in[0]) );
          if ( status <= 0 ) MT_ErrorParse( "parsing -p-ref...\n", 0 );
          i += 1;
          if ( i >= argc)    MT_ErrorParse( "parsing -p-ref...\n", 0 );
          status = sscanf( argv[i],"%lf",&(par->p_in[1]) );
          if ( status <= 0 ) MT_ErrorParse( "parsing -p-ref...\n", 0 );
          i += 1;
          if ( i >= argc)    MT_ErrorParse( "parsing -p-ref...\n", 0 );
          status = sscanf( argv[i],"%lf",&(par->p_in[2]) );
          if ( status <= 0 ) MT_ErrorParse( "parsing -p-ref...\n", 0 );
          i += 1;
          if ( i >= argc)    MT_ErrorParse( "parsing -p-ref...\n", 0 );
          status = sscanf( argv[i],"%lf",&(par->p_in[3]) );
          if ( status <= 0 ) MT_ErrorParse( "parsing -p-ref...\n", 0 );
          par->flag_planes_in=1;
      }


      else if ( strcmp ( argv[i], "-p-ext" ) == 0 || strcmp ( argv[i], "-p-flo" ) == 0 ) {
          i += 1;
          if ( i >= argc)    MT_ErrorParse( "parsing -p-flo...\n", 0 );
          status = sscanf( argv[i],"%lf",&(par->p_ext[0]) );
          if ( status <= 0 ) MT_ErrorParse( "parsing -p-flo...\n", 0 );
          i += 1;
          if ( i >= argc)    MT_ErrorParse( "parsing -p-flo...\n", 0 );
          status = sscanf( argv[i],"%lf",&(par->p_ext[1]) );
          if ( status <= 0 ) MT_ErrorParse( "parsing -p-flo...\n", 0 );
          i += 1;
          if ( i >= argc)    MT_ErrorParse( "parsing -p-flo...\n", 0 );
          status = sscanf( argv[i],"%lfd",&(par->p_ext[2]) );
          if ( status <= 0 ) MT_ErrorParse( "parsing -p-flo...\n", 0 );
          i += 1;
          if ( i >= argc)    MT_ErrorParse( "parsing -p-flo...\n", 0 );
          status = sscanf( argv[i],"%lf",&(par->p_ext[3]) );
          if ( status <= 0 ) MT_ErrorParse( "parsing -p-flo...\n", 0 );
          par->flag_planes_ext=1;
      }

      else if ( strcmp ( argv[i], "-rigid" ) == 0 ){
        par->flag_affine=0;
      }

      else if ( strcmp ( argv[i], "-affine" ) == 0 ){
        par->flag_affine=1;
      }

      else if ( strcmp ( argv[i], "-force-trsf" ) == 0 ){
          i += 1;
          if ( i >= argc)    MT_ErrorParse( "parsing -force-trsf...\n", 0 );
          status = sscanf( argv[i],"%d",&(par->force) );
          if ( status <= 0 ) {
              if ( strcmp ( argv[i], "direct" ) == 0 )
                  par->force = -2;
              else if ( strcmp ( argv[i], "reversal" ) == 0 || strcmp ( argv[i], "mirror" ) == 0 || strcmp ( argv[i], "miroir" ) == 0 )
                  par->force = -3;
              else MT_ErrorParse( "parsing -force-trsf...\n", 0 );
          }
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

      else if ( strcmp ( argv[i], "-background-in" ) == 0 || strcmp ( argv[i], "-background-ref" ) == 0 ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( argv[i-1], 0 );
        status = sscanf( argv[i], "%d", &(par->background_in) );
        if ( status <= 0 ) MT_ErrorParse( argv[i-1], 0 );
      }

      else if ( strcmp ( argv[i], "-background-ext" ) == 0 || strcmp ( argv[i], "-background-flo" ) == 0) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( argv[i-1], 0 );
        status = sscanf( argv[i], "%d", &(par->background_ext) );
        if ( status <= 0 ) MT_ErrorParse( argv[i-1], 0 );
      }



      /*--- option inconnue ---*/
      else {
        sprintf(text,"unknown option %s\n",argv[i]);
        MT_ErrorParse(text, 0);
      }
    }
    /*--- saisie des noms d'images ---*/
    /*else if ( argv[i][0] != 0 ) {
      if ( nb == 0 ) {
        strncpy( par->names.in, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else if ( nb == 1 ) {
        strncpy( par->names.out, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else {
        fprintf(stderr, "N=%d\nS=%d\nV=%d\n", par->N, par->S, par->V);
        MT_ErrorParse("too much file names when parsing\n", 0 );
      }
    }*/
    i += 1;
  }

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
  VT_Names( &(par->name_gradients) );
  VT_Names( &(par->name_labels) );
  VT_Names( &(par->name_planes) );
  VT_Names( &(par->name_res) );
  VT_Names( &(par->name_dices) );
  VT_Names( &(par->name_residuals) );
  par->flag_name_planes_in=0;
  par->flag_name_planes_ext=0;
  par->flag_planes_in=0;
  par->flag_planes_ext=0;
  par->flag_affine=0;
  BAL_InitEstimator( &(par->estimator) );
  par->estimator.type = TYPE_LS;
  par->delta=0;
  par->force=-1;
  par->npairs=0;
  par->dices_out=0;
  par->dices_in=0;
  par->background_in=0;
  par->background_ext=0;
}

