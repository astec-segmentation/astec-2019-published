/*************************************************************************
 * labelBarycenters.c -
 *
 * $Id: labelBarycenters.c,v 1.0 2017/08/01 14:47:51 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2017/08/01
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>

typedef enum _UNIT {
 REAL,
 VOXEL
} _UNIT;

typedef struct local_par {
    vt_names names;
    int flag_voxelsize;
    int flag_volume;
    _UNIT unit;
} local_par;

typedef struct {
    double x;
    double y;
    double z;
    int weight;
    int label;
} vt_barycentre;

typedef struct vt_image_barycentres {
    vt_barycentre *barycentres;
    int n;
    double vx;
    double vy;
    double vz;
} vt_image_barycentres;

/*------- Definition des fonctions statiques ----------*/

static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
static int labelBarycenters(vt_image labels, vt_image_barycentres* b, int background);
static int ncell(vt_image img, int background);
static int allocImageBarycentresWithImage(vt_image_barycentres* b, vt_image labels, int background);
static int initImageBarycentres(vt_image_barycentres* b, vt_image image, int background);
static void freeImageBarycentres(vt_image_barycentres* b);
static int findLabel(vt_image_barycentres b, int val);
static int partition(vt_image_barycentres b, int *tableau, int deb, int fin);
static void tri_rapide_bis(vt_image_barycentres b, int *tableau, int debut,int fin);
static void tri_rapide(vt_image_barycentres b, int *tableau);


static int _verbose_ = 0;

static char *usage = "[image-in] [file-out] \n\
\t [-print-voxelsize|-print-vs|-pvs] [-print-volume|-print-vol|-pv]\n\
\t [-unit voxel|real]\n\
\t [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'file-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t N.B. : ce programme calcul les barycentres des differents labels d'une image encodee sur 8 ou 16 bits seulement\n\
\t -print-voxelsize : genere une en-tete de format \"Voxelsize %f %f %f\" renseignant les resolutions d'image vx, vy, vz\n\
\t -print-volume : ajoute en fin de chaque ligne du fichier de sortie l'information du volume (nb de pixels) de chaque label\n\
\t -unit : calcule les barycentres en coordonnees reelles (real, par defaut) ou voxelliques (voxel)\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];




int main( int argc, char *argv[] )
{
    local_par par;
    vt_image *image;
    vt_image_barycentres b;
    vt_barycentre tmpBary, *theBary;
    FILE* fileout;
    int background = 0;
    int i;
    int *index;

    /*--- initialisation des parametres ---*/
    VT_InitParam( &par );

    /*--- lecture des parametres ---*/
    VT_Parse( argc, argv, &par );

    /*--- lecture de l'image d'entree ---*/
    image = _VT_Inrimage( par.names.in );
    if ( image == (vt_image*)NULL )
      VT_ErrorParse("unable to read input image\n", 0);

    /*--- operations eventuelles sur l'image d'entree ---*/
    if ( par.names.inv == 1 )  VT_InverseImage( image );
    if ( par.names.swap == 1 ) VT_SwapImage( image );


    /*--- operations sur les images de labels ---*/
    if (allocImageBarycentresWithImage(&b, *image, background) != 1) {
        VT_FreeImage( image );
        VT_Free( (void**)&image);
        VT_ErrorParse("unable to allocate barycentre image\n", 0);
        return(0);
    }

    if(labelBarycenters(*image, &b, background) != 1)
    {
        VT_FreeImage( image );
        VT_Free( (void**)&image);
        freeImageBarycentres(&b);
        VT_ErrorParse("unable to compute barycentre image\n", 0);
        return(0);
    }

    /* barycenters are computed in voxel units
     */
    switch ( par.unit ) {
    default :
    case REAL :
        for (i=0;i<b.n ; i++) {
            theBary = &(b.barycentres[i]);
            tmpBary.x = image->qform_rfv[0][0] * theBary->x
                    + image->qform_rfv[0][1] * theBary->y
                    + image->qform_rfv[0][2] * theBary->z
                    + image->qform_rfv[0][3];
            tmpBary.y = image->qform_rfv[1][0] * theBary->x
                    + image->qform_rfv[1][1] * theBary->y
                    + image->qform_rfv[1][2] * theBary->z
                    + image->qform_rfv[1][3];
            tmpBary.z = image->qform_rfv[2][0] * theBary->x
                    + image->qform_rfv[2][1] * theBary->y
                    + image->qform_rfv[2][2] * theBary->z
                    + image->qform_rfv[2][3];
            theBary->x = tmpBary.x;
            theBary->y = tmpBary.y;
            theBary->z = tmpBary.z;
        }
        break;
    case VOXEL :
        break;
    }


    /*--- liberations memoires partielles ---*/
    VT_FreeImage( image );
    VT_Free( (void**)&image);

    /*--- tri ---*/
    index=malloc(sizeof(int)*b.n);
    tri_rapide(b, index);


    /*--- ecriture des barycentres ---*/
    if (strcmp(par.names.out, ">") == 0)
        fileout=stdout;
    else
        fileout = fopen( par.names.out, "w");

    if (par.flag_voxelsize)
        fprintf(fileout, "Voxelsize %f %f %f\n", b.vx, b.vy, b.vz);
    if (_verbose_){
        if (par.flag_volume)
            fprintf(fileout, "# LABEL X Y Z VOLUME\n");
        else
            fprintf(fileout, "# LABEL X Y Z\n");
    }
    for (i = 0 ; i < b.n ; i++)
    {
        theBary = &(b.barycentres[index[i]]);
        if (par.flag_volume)
            fprintf(fileout, "%d %f %f %f %d\n", theBary->label, theBary->x, theBary->y, theBary->z, theBary->weight);
        else
            fprintf(fileout, "%d %f %f %f\n", theBary->label, theBary->x, theBary->y, theBary->z);
    }

    if (strcmp(par.names.out, ">") != 0)
        fclose( fileout );

    /*--- liberations memoires ---*/
    freeImageBarycentres(&b);
    free(index);

    return(0);
}







static void VT_Parse( int argc,
                      char *argv[],
                      local_par *par )
{
    int i, nb;
    char text[STRINGLENGTH];

    if ( VT_CopyName( program, argv[0] ) != 1 )
      VT_Error("Error while copying program name", (char*)NULL);
    if ( argc == 1 ) VT_ErrorParse("\n", 0 );

    /*--- lecture des parametres ---*/
    i = 1; nb = 0;
    while ( i < argc ) {
      if ( argv[i][0] == '-' ) {
        if ( argv[i][1] == '\0' ) {
          if ( nb == 0 ) {
            /*--- standart input ---*/
            strcpy( par->names.in, "<" );
            nb += 1;
          }
        }
        /*--- arguments generaux ---*/
        else if ( strcmp ( argv[i], "-help" ) == 0 ) {
          VT_ErrorParse("\n", 1);
        }
        else if ( strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0') {
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

        else if ( strcmp ( argv[i], "-print-voxelsize" ) == 0
                  || strcmp ( argv[i], "-print-vs" ) == 0
                  || strcmp ( argv[i], "-vs" ) == 0
                  || strcmp ( argv[i], "-voxelsize" ) == 0
                  || (strcmp ( argv[i], "-pvs" ) == 0 && argv[i][4] == '\0') ) {
          par->flag_voxelsize = 1;
        }
        else if ( strcmp ( argv[i], "-print-volume" ) == 0
                  || strcmp ( argv[i], "-print-vol" ) == 0
                  || strcmp ( argv[i], "-volume" ) == 0
                  || strcmp ( argv[i], "-vol" ) == 0
                  || (strcmp ( argv[i], "-pv" ) == 0 && argv[i][3] == '\0') ) {
          par->flag_volume = 1;
        }

        else if ( strcmp ( argv[i], "-unit" ) == 0 ) {
          i++;
          if ( i >= argc)
            VT_ErrorParse( "no argument for -unit", 0);
          if ( strcmp ( argv[i], "voxel" ) == 0 ) {
            par->unit = VOXEL;
          }
          else if ( strcmp ( argv[i], "real" ) == 0 ) {
              par->unit = REAL;
          }
          else {
            VT_ErrorParse( "unknown unit", 0);
          }
        }

        /*--- option inconnue ---*/
      else {
        sprintf(text,"unknown option %s\n",argv[i]);
        VT_ErrorParse(text, 0);
      }
    }

    /*--- saisie des noms d'images ---*/
    else if ( argv[i][0] != 0 ) {
      if ( nb == 0 ) {
        strncpy( par->names.in, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else if ( nb == 1 ) {
        strncpy( par->names.out, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else
        VT_ErrorParse("too much file names when parsing\n", 0 );
    }
    i += 1;
  }

  /*--- s'il n'y a pas assez de noms ... ---*/
  if (nb == 0) {
    strcpy( par->names.in,  "<" );  /* standart input */
    strcpy( par->names.out, ">" );  /* standart output */
  }
  if (nb == 1)
    strcpy( par->names.out, ">" );  /* standart output */

}



static void VT_ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"Erreur : %s",str);
  exit( 1 );
}








static void VT_InitParam( local_par *par )
{
  VT_Names( &(par->names) );
  par->flag_volume = 0;
  par->flag_voxelsize = 0;
  par->unit = REAL;
}


int ncell(vt_image img, int background)
{
    char *proc="ncell";

    unsigned char *bufu8=NULL;
    unsigned short int *bufu16=NULL;

    size_t i;
    int j, val, valold;
    int n=0;
    int *labels;
    int delta=256;

    int add;

    switch (img.type) {
    default:
        if(1 || _verbose_)
            fprintf(stderr, "%s: such image type not handled yet\n", proc);
        return(-1);
    case UCHAR:
    case SCHAR:
        bufu8=(unsigned char *)img.buf;
        break;
    case USHORT:
    case SSHORT:
        bufu16=(unsigned short int *)img.buf;
        break;
    }



    labels=malloc(delta*sizeof(int));
    if(labels==(int*)NULL)
    {
        if(1 || _verbose_)
            fprintf(stderr, "%s: problem while allocating label vector\n", proc);
        return(-1);
    }


    valold=background;
    for (i=0 ; i<img.dim.x*img.dim.y*img.dim.z ; i++) {
        switch (img.type){
        default:
            if(1 || _verbose_)
                fprintf(stderr, "%s: such image type not handled yet\n", proc);
            free(labels);
            return(-1);
        case UCHAR:
        case SCHAR:
            val=(int)bufu8[i];
            break;
        case USHORT:
        case SSHORT:
            val=(int)bufu16[i];
            break;
        }

        if (val==valold || val==background) continue;
        valold=val;

        add=1;
        for (j=0 ; j<n ; j++) {
            if (val!=labels[j]) continue;
            add=0;
            break;
        }
        if (add==0) continue;

        if (n%delta == 0)
        {
            labels=realloc(labels, (delta+n)*sizeof(int));
            if(labels==(int*)NULL)
            {
                if(1 || _verbose_)
                    fprintf(stderr, "%s: problem while reallocating label vector\n", proc);
                return(-1);
            }
        }
        labels[n++]=val;
    }

    free(labels);
    return (n);
}



int allocImageBarycentresWithImage(vt_image_barycentres* b, vt_image labels, int background)
{
    char *proc="allocImageBarycentresWithImage";
    int n=ncell(labels, background);
    if (n<=0) {
        fprintf(stderr, "%s: negative number of cells\n", proc);
        return(0);
    }
    b->barycentres=malloc(n*sizeof(vt_barycentre));
    if (b->barycentres == (vt_barycentre*)NULL) {
        fprintf(stderr, "%s: malloc failed\n", proc);
        return(0);
    }
    b->n = n;
    return(initImageBarycentres(b, labels, background));
}




int initImageBarycentres(vt_image_barycentres* b, vt_image image, int background)
{
    char *proc="initImageBarycentres";
    unsigned char *bufu8=NULL;
    unsigned short int *bufu16=NULL;

    size_t i;
    int j, val, valold;
    int n=0;
    int *labels;
    int delta=256;

    int add;

    b->vx=image.siz.x;
    b->vy=image.siz.y;
    b->vz=image.siz.z;

    switch (image.type) {
    default:
        if(1 || _verbose_)
            fprintf(stderr, "%s: such image type not handled yet\n", proc);
        return(-1);
    case UCHAR:
    case SCHAR:
        bufu8=(unsigned char *)image.buf;
        break;
    case USHORT:
    case SSHORT:
        bufu16=(unsigned short int *)image.buf;
        break;
    }



    labels=malloc(delta*sizeof(int));
    if(labels==(int*)NULL)
    {
        if(1 || _verbose_)
            fprintf(stderr, "%s: problem while allocating label vector\n", proc);
        return(-1);
    }

    valold=background;
    for (i=0 ; i<image.dim.x*image.dim.y*image.dim.z ; i++) {
        switch (image.type){
        default:
            if(1 || _verbose_)
                fprintf(stderr, "%s: such image type not handled yet\n", proc);
            free(labels);
            return(-1);
        case UCHAR:
        case SCHAR:
            val=(int)bufu8[i];
            break;
        case USHORT:
        case SSHORT:
            val=(int)bufu16[i];
            break;
        }

        if (val==valold || val==background) continue;
        valold=val;

        add=1;
        for (j=0 ; j<n ; j++) {
            if (val!=labels[j]) continue;
            add=0;
            break;
        }
        if (add==0) continue;

        if (n%delta == 0)
        {
            labels=realloc(labels, (delta+n)*sizeof(int));
            if(labels==(int*)NULL)
            {
                if(1 || _verbose_)
                    fprintf(stderr, "%s: problem while reallocating label vector\n", proc);
                return(-1);
            }
        }
        labels[n++]=val;
    }
    if(n != b->n)
    {
        if (1 || _verbose_)
            fprintf(stderr, "%s: incompatible vt_barycentre_image field n = %d and vt_image number of cells = %d\n", proc, b->n, n);
        free(labels);
        return(0);
    }

    /*  To improve */
    for (i=0; i<(size_t)b->n; i++) {
        vt_barycentre *theBary = &(b->barycentres[i]);
        theBary->label=labels[i];
        theBary->weight=0;
        theBary->x=0;
        theBary->y=0;
        theBary->z=0;
    }

    free(labels);
    return (1);
}




int labelBarycenters(vt_image labels, vt_image_barycentres* b, int background)
{
    char *proc="labelBarycenters";

    unsigned char ***arrayu8=NULL;
    unsigned short int ***arrayu16=NULL;

    size_t i, j, k;
    int val, valold;
    int ind;

    vt_barycentre *theBary;

    switch (labels.type) {
    default:
        if(_verbose_)
            fprintf(stderr, "%s: such image type not handled yet\n", proc);
        return(-1);
    case UCHAR:
    case SCHAR:
        arrayu8=(unsigned char ***)labels.array;
        break;
    case USHORT:
    case SSHORT:
        arrayu16=(unsigned short int ***)labels.array;
        break;
    }

    valold=background;
    ind=0;
    for (k=0 ; k<labels.dim.z ; k++)
    for (j=0 ; j<labels.dim.y ; j++)
    for (i=0 ; i<labels.dim.x ; i++) {
        switch (labels.type){
        default:
            if(_verbose_)
                fprintf(stderr, "%s: such image type not handled yet\n", proc);
            return(-1);
        case UCHAR:
        case SCHAR:
            val=(int)arrayu8[k][j][i];
            break;
        case USHORT:
        case SSHORT:
            val=(int)arrayu16[k][j][i];
            break;
        }
        if(val==background) continue;
        if(val!=valold) {
            ind=findLabel(*b, val);
            if(ind<0)
            {
                if(_verbose_)
                    fprintf(stderr, "%s: did not find the index of label # %d\n", proc, val);
                return(-1);
            }
        }
        valold=val;
        theBary=&(b->barycentres[ind]);
        if(theBary->label != val) {
            if(_verbose_)
                fprintf(stderr, "%s: unexpected returned index # %d for value %d\n", proc, ind, val);
            return(-1);
        }
        theBary->weight += 1;
        theBary->x += i;
        theBary->y += j;
        theBary->z += k;
    }

    for (ind=0 ; ind<b->n ; ind++) {
        theBary=&(b->barycentres[ind]);
        if (theBary->weight==0) continue;
        theBary->x /= theBary->weight;
        theBary->y /= theBary->weight;
        theBary->z /= theBary->weight;
        /*theBary->weight *= b->vx * b->vy * b->vz;*/
    }

    return(1);
}

int findLabel(vt_image_barycentres b, int val)
{
    int i;
    vt_barycentre theBary;
    for (i=0 ; i<b.n ; i++)
    {
        theBary=b.barycentres[i];
        if (val==theBary.label)
            return(i);
    }
    return (-1);
}


void freeImageBarycentres(vt_image_barycentres* b)
{
    free(b->barycentres);
    b=(vt_image_barycentres*) NULL;
}





int partition(vt_image_barycentres b, int *tableau, int deb, int fin)
{
    int compt=deb;
    int pivot=b.barycentres[tableau[deb]].label;
    int i;
    int v;

    for(i=deb+1;i<=fin;i++)
    {
        if(b.barycentres[tableau[i]].label<pivot)
        {
            compt++;
            v=tableau[i];
            tableau[i]=tableau[compt];
            tableau[compt]=v;
        }
    }
    v=tableau[deb];
    tableau[deb]=tableau[compt];
    tableau[compt]=v;
    return(compt);
}

void tri_rapide_bis(vt_image_barycentres b, int *tableau, int debut,int fin)
{
    if(debut<fin)
    {
        int pivot=partition(b, tableau,debut,fin);
        tri_rapide_bis(b, tableau,debut,pivot-1);
        tri_rapide_bis(b, tableau,pivot+1,fin);
    }
}

void tri_rapide(vt_image_barycentres b, int *tableau)
{
    int i;
    for (i = 0 ; i<b.n ; i++)
        tableau[i]=i;
     tri_rapide_bis(b, tableau, 0, b.n-1);
}
