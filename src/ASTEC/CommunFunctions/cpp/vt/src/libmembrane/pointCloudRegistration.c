/*************************************************************************
 * pointCloudRegistration.c -
 *
 * $Id: pointCloudRegistration.c,v 1.0 2017/07/27 13:37:51 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2017/07/27
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
/*#define PI 3.141592654*/

typedef struct local_par {
  vt_names name_labels;
  vt_names name_dices;
  vt_names name_res;
  vt_names name_residuals;
  int flag_affine;
  bal_estimator estimator;
  double delta;
  int npairs;
  int pairs[2][TABLENGTH];
  int dices_in;
  int dices_out;
  int background_in;
  int background_ext;
  int flag_skip_unfound;
} local_par;


/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );

static int _verbose_ = 0;


static char *usage = "[-label-[ref|in] %s] [-label-[flo|ext] %s] [-background-[ref|in] %d] [-background-[flo|ext] %d]\n\
\t [-pairs|pair|p] %d %d [...] | %s] [-skip[-not-found]]\n\
\t [-rigid | -affine] \n\
\t [-estimator-type|-estimator|-es-type wlts|lts|wls|ls]\n\
\t [-lts-fraction %lf] [-lts-deviation %f] [-lts-iterations %d]\n\
\t [-fluid-sigma|-lts-sigma[-ll|-hl] %lf %lf %lf]\n\
\t [-pairs-out %s] [-pairs-in %s] [-trsf %s] [-dices %s] [-residuals %s]\n\
\t [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t -label-[ref|flo] %s : nom des fichiers image ou texte de labels ref|flo\n\
\t N.B.: le cas echeant, le fichier texte attendu doit respecter le format suivant : \n\
\t 1ere ligne (optionnel) : en-tête au format \"Voxelsize %f %f %f\" stipulant le voxelsize d'origine de l'image de laquelle proviennent les coordonnees barycentriques des labels\n\
\t liste des barycentres  : une ligne par barycentre ; chaque ligne respecte le format \"%d %f %f %f\" ou \"%d %f %f %f %f\" correspondant au label, coordonnee en x, en y, en z et (le cas echeant) au volume (ou poids) du label \n\
\t lignes de commentaires : commencent avec le caractere \"#\"\n\
\t -background-[ref|flo] %d : label de fond a ignorer dans l'image [ref|flo] correspondante (defaut = 0)\n\
\t -pairs %d %d [...] | %s : associe les paires de labels specifiees (ref-flo [ref-flo [...]])\n\
\t -skip-not-found : poursuit le calcul du recalage malgre d'eventuelles associations de labels inexistantes dans les labels d'entree\n\
\t -trsf %s : fichier dans lequel est enregistree la transformation T_flo<-ref calculee\n\
\t -pairs-in %s : fichier dans lequel sont enregistres les appariements utilises pour le calcul\n\
\t -pairs-out %s : fichier dans lequel sont enregistres les appariements trouves apres calcul\n\
\t ou seulement sur l'iteration donnee\n\
\t -residuals %s : enregistre les valeurs de residus pour les paires de labels mis en correspondance\n\
\t -dices %s : calcule l'indice de similarite de Dice entre chaque paire de labels en correspondance ; \n\
\t             OPTION VALABLE UNIQUEMENT AVEC DES PARAMETRES -label-[ref|flo] DE TYPE IMAGE \n\
\t ### transformation type ###\n\
\t -rigid : computes rigid transformation (set as default)\n\
\t -affine : computes affine transformation \n\
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
  vt_image *images_labels[2];

  double T_ext_in[16];
  double r;
  double sd;
  int *pairs[2];
  int npairs;



  int i, j;


  vt_image_barycentres b_in;
  vt_image_barycentres b_ext;
  vt_image_barycentres b_in_trsf;

  /*--- initialisation des parametres ---*/
  MT_InitParam( &par );

  /*--- lecture des parametres ---*/
  MT_Parse( argc, argv, &par );


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

  /*return(0);*/

  /* Allocation des images des barycentres transformes */

  if (VT_AllocImageBarycentres(&b_in_trsf, b_in.n) != 1) {
    VT_FreeImageBarycentres(&b_in);
    VT_FreeImageBarycentres(&b_ext);
    MT_ErrorParse("unable to allocate image of barycentres b_in_trsf...\n", 0);
  }

  if(par.npairs==0) {
    MT_ErrorParse("Error: there is no given labels association. Exiting.", 0);
  }
  else {
    /* CALCUL DE T_ext_in = TRANSFO DE IN_THETA VERS EXT */
      pairs[0] = vtmalloc( par.npairs*sizeof(int), "", argv[0] );
      if(pairs[0]==NULL) {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext);
          MT_ErrorParse("unable to allocate a vector...\n", 0);
      }
      pairs[1] = vtmalloc( par.npairs*sizeof(int), "", argv[0] );
      if(pairs[1]==NULL) {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext);
          vtfree(pairs[0]); pairs[0]=NULL;
          MT_ErrorParse("unable to allocate a vector...\n", 0);
      }
      j=0;
      npairs=par.npairs;
      for (i=0;i<par.npairs;i++) {
          pairs[0][j]=VT_FindLabel(b_in, par.pairs[0][i] );
          pairs[1][j]=VT_FindLabel(b_ext, par.pairs[1][i]);
          if (pairs[0][j]<0 || pairs[1][j]<0 )
          {
              if (par.flag_skip_unfound) {
                npairs--;
                if (_verbose_)
                    fprintf(stdout, "Warning label association (%d - %d) not found.\n", par.pairs[0][i], par.pairs[1][i]);
              }
              else{
                VT_FreeImageBarycentres(&b_in);
                VT_FreeImageBarycentres(&b_in_trsf);
                VT_FreeImageBarycentres(&b_ext);
                vtfree(pairs[0]); pairs[0]=NULL;
                vtfree(pairs[1]); pairs[1]=NULL;
                MT_ErrorParse("label not find, exiting the process... The user may try again using the option '-skip'.\n", 0);
              }
          }
          else {
              j++;
          }
      }
      /* GM: _VT_VERBOSE_ is a boolean
       * */
      if (_VT_VERBOSE_ ) fprintf(stdout, "npairs=%d\nPairs : \n", npairs);
      /* gm : a supprimer...
      if (_VT_VERBOSE_ )
        for (i=0;i<par.npairs;i++) {
          fprintf(stdout, "( %d \t%d )\n", par.pairs[0][i], par.pairs[1][i]);
        }
      */
      if( VT_ComputePointTrsf(pairs, npairs, b_in, b_ext, &r, &sd, T_ext_in, par.flag_affine, par.estimator) != 1) {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext);
          MT_ErrorParse("unable to compute the points to points transformation...\n", 0);
          vtfree(pairs[0]); pairs[0]=NULL;
          vtfree(pairs[1]); pairs[1]=NULL;
      }

      if(par.name_res.in[0] != '\0')
      {
        /*fprintf(stdout, "npairs = %d\n", npairs); */
          FILE* fichier = NULL;
          /*fprintf(stderr, "par.name_res.in : %s\n", par.name_res.in);*/
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
          /*fprintf(stderr, "file opened\n");*/
          fprintf(fichier, "#Liste d'appariements (ref-flo) entres :\n#\n");
          vt_barycentre theB_in;
          vt_barycentre theB_ext;
          for (i=0; i<npairs; i++) {
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
      fprintf(stdout, "#Transformation matrix (flo<-ref):\n");
      fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_ext_in[0],T_ext_in[1],T_ext_in[2],T_ext_in[3]);
      fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_ext_in[4],T_ext_in[5],T_ext_in[6],T_ext_in[7]);
      fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_ext_in[8],T_ext_in[9],T_ext_in[10],T_ext_in[11]);
      fprintf(stdout, "\t%f\t %f\t %f\t %f\n", T_ext_in[12],T_ext_in[13],T_ext_in[14],T_ext_in[15]);
      fprintf(stdout, "residu = %f\nstd deviation = %f\n", r, sd);
  }
  if(par.name_residuals.out[0] != '\0')
  {
      FILE* fichier = NULL;
      fichier = fopen(par.name_residuals.out, "w");
      if( fichier == NULL)
      {
        MT_ErrorParse("Erreur pour l'ecriture du resultat dans un fichier\n", 0);
      }
      fprintf(fichier, "%%Residu:\nresidual=%f ;\n", r);
      fprintf(fichier, "%%Ecart-type du Residu:\nstd_dev=%f ;\n", sd);

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
      fprintf(fichier, "#Residu = %f\n#Std deviation = %f\n", r, sd);
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

          MT_ErrorParse("Erreur pour l'ecriture de la liste d'appariements dans un fichier\n", 0);
      }

      if (VT_TransformImageBarycentres(T_ext_in, b_in, &b_in_trsf) != 1) {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext);
          MT_ErrorParse("unable to process the image of barycentres transformation\n",0);
      }

      if (VT_ComputePairsOfPoints(b_in_trsf, b_ext, pairs, &npairs, &r) != 1)
      {
          VT_FreeImageBarycentres(&b_in);
          VT_FreeImageBarycentres(&b_in_trsf);
          VT_FreeImageBarycentres(&b_ext);
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




  /*--- liberations memoires ---*/
  VT_FreeImageBarycentres(&b_in);
  VT_FreeImageBarycentres(&b_ext);
  VT_FreeImageBarycentres(&b_in_trsf);

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
          par->name_labels.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
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


      else if  (strcmp(argv[i], "-skip") == 0 || strcmp(argv[i], "-skip-not-found") == 0 ) {
        par->flag_skip_unfound = 1;
      }


      else if ( strcmp ( argv[i], "-rigid" ) == 0 ){
        par->flag_affine=0;
      }

      else if ( strcmp ( argv[i], "-affine" ) == 0 ){
        par->flag_affine=1;
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
  VT_Names( &(par->name_labels) );
  VT_Names( &(par->name_res) );
  VT_Names( &(par->name_dices) );
  VT_Names( &(par->name_residuals) );
  par->flag_affine=0;
  BAL_InitEstimator( &(par->estimator) );
  par->estimator.type = TYPE_LS;
  par->delta=0;
  par->npairs=0;
  par->dices_out=0;
  par->dices_in=0;
  par->background_in=0;
  par->background_ext=0;
  par->flag_skip_unfound=0;
}

