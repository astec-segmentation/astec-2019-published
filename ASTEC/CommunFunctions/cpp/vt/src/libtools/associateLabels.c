/*************************************************************************
 * associateLabels.c -
 *
 * $Id: associateLabels.c,v 1.0 2014/09/25 13:14:51 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2014/09/25
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>

typedef struct local_par {
    vt_names name_labels;
    vt_names name_outs;
    int flag_force_labels;
    int flag_zeros;
    int flag_ref;
} local_par;


typedef struct {
  int l,c;
  int *m;
} _MATRIX;


/*------- Definition des fonctions statiques ----------*/
static int _read_mat(char *name, _MATRIX *m);
static int _read_mat_3(char *name, _MATRIX *m);
static int _write_mat(char *name, _MATRIX *m);
static void _init_mat( _MATRIX *mat );
static int _alloc_mat ( _MATRIX *mat, int l, int c );
static void _free_mat ( _MATRIX *mat );
static void _compute_table ( _MATRIX *mat, int flag_ref );
static int _ind( _MATRIX mat, int c, int val);

static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );

static int _verbose_ = 0;

static char *usage = "[image-[in|ref]] [image-[ext-flo]] [labels-file]\n\
\t [image-[in|ref]-out] [image-[ext|flo]-out] [labels-file-table]\n\
\t [-z] [-f] [-ref] [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t -z : les zeros en entree sont conserves en sortie\n\
\t -f : force la valeur de sortie des labels (tableau d'entree n*3)\n\
\t -ref : utilise les labels de l'image in comme labels de reference\n\
\t labels-file-table : fichier de correspondances label-ref / label-flo (tel qu'en sortie de planeRegistration ou pointCloudRegistration -pairs-out %s)\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];



int main( int argc, char *argv[] )
{

  local_par par;
  char name[DOUBLESTRINGLENGTH];
  vt_image *image_in;
  vt_image *image_ext;

  _MATRIX mat;
  int i, t, val, ind=0, oldval;
  unsigned char *bufU8_in=NULL, *inU8=NULL;
  unsigned short int *bufU16_in=NULL, *inU16=NULL;

  unsigned char *bufU8_ext=NULL, *extU8=NULL;
  unsigned short int *bufU16_ext=NULL, *extU16=NULL;

  vt_image out_in;
  vt_image out_ext;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );

  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );

  /*--- lecture des images labels d'entree ---*/
  image_in = _VT_Inrimage( par.name_labels.in );
  if ( image_in  == (vt_image*)NULL )
    VT_ErrorParse("unable to read image label-in\n", 0);

  image_ext = _VT_Inrimage( par.name_labels.ext );
  if ( image_ext  == (vt_image*)NULL )
    VT_ErrorParse("unable to read image label-ext\n", 0);

  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.name_labels.inv == 1 )  {
      VT_InverseImage( image_in);
      VT_InverseImage( image_ext);
  }

  if ( par.name_labels.swap == 1 ) {
      VT_SwapImage( image_in);
      VT_SwapImage( image_ext);
  }

  switch (image_in->type) {
  default:
      VT_FreeImage(image_in);
      VT_FreeImage(image_ext);
      VT_ErrorParse("such image-in type not handled yet\n", 0);
      return( 1 );
  case UCHAR:
  case SCHAR:
      inU8 = (unsigned char *)image_in->buf;
      break;
  case USHORT:
  case SSHORT:
      inU16 = (unsigned short int *)image_in->buf;
      break;
  }

  switch (image_ext->type) {
  default:
      VT_FreeImage(image_in);
      VT_FreeImage(image_ext);
      VT_ErrorParse("such image-ext type not handled yet\n", 0);
      return( 1 );
  case UCHAR:
  case SCHAR:
      extU8 = (unsigned char *)image_ext->buf;
      break;
  case USHORT:
  case SSHORT:
      extU16 = (unsigned short int *)image_ext->buf;
      break;
  }

  /*--- lecture du fichier des labels ---*/


  _init_mat(&mat);
  if(par.flag_force_labels==0){
    if (_read_mat(par.name_labels.out, &mat)!=1)
    {
      VT_FreeImage(image_in);
      VT_FreeImage(image_ext);
      VT_ErrorParse("unable to read input label file\n", 0);
    }

    _compute_table(&mat, par.flag_ref);

  }
  else {
    if (_read_mat_3(par.name_labels.out, &mat)!=1)
    {
      VT_FreeImage(image_in);
      VT_FreeImage(image_ext);
      VT_ErrorParse("unable to read input label file\n", 0);
    }
  }

  if (par.name_outs.out[0] != '\0') {
    if(_write_mat(par.name_outs.out, &mat) != 1)  {
      VT_FreeImage(image_in);
      VT_FreeImage(image_ext);
      _free_mat(&mat);
      VT_ErrorParse("unable to write label table file\n", 0);
    }
  }

  /*--- allocation des images de sortie ---*/
  if (mat.l<255) t =  UCHAR;
  else if (mat.l<32767) t = USHORT;
  else {
      VT_FreeImage(image_in);
      VT_FreeImage(image_ext);
      _free_mat(&mat);
      VT_ErrorParse("too many labels\n", 0);
  }

  sprintf( name, "%s", par.name_outs.in );
  VT_Image( &out_in );
  VT_InitVImage( &out_in, name, 1,
             image_in->dim.x, image_in->dim.y, image_in->dim.z, t );

  out_in.siz.x = image_in->siz.x;
  out_in.siz.y = image_in->siz.y;
  out_in.siz.z = image_in->siz.z;

  if ( VT_AllocImage( &out_in ) != 1 ) {
      VT_FreeImage(image_in);
      VT_FreeImage(image_ext);
      _free_mat(&mat);
      VT_ErrorParse("unable to allocate output image-in\n", 0);
  }

  sprintf( name, "%s", par.name_outs.ext );
  VT_Image( &out_ext );
  VT_InitVImage( &out_ext, name, 1,
             image_ext->dim.x, image_ext->dim.y, image_ext->dim.z, t );

  out_ext.siz.x = image_ext->siz.x;
  out_ext.siz.y = image_ext->siz.y;
  out_ext.siz.z = image_ext->siz.z;

  if ( VT_AllocImage( &out_ext) != 1 ) {
      VT_FreeImage(image_in);
      VT_FreeImage(image_ext);
      VT_FreeImage(&out_in);
      _free_mat(&mat);
      VT_ErrorParse("unable to allocate output image-ext\n", 0);
  }

  switch (t) {
  default:
      VT_FreeImage(image_in);
      VT_FreeImage(image_ext);
      VT_FreeImage(&out_in);
      VT_FreeImage(&out_ext);
      _free_mat(&mat);
      VT_ErrorParse("error with out image types\n", 0);
      return( 1 );
  case UCHAR:
      bufU8_in = (unsigned char *)out_in.buf;
      bufU8_ext = (unsigned char *)out_ext.buf;
      break;
  case USHORT:
      bufU16_in = (unsigned short int *)out_in.buf;
      bufU16_ext = (unsigned short int *)out_ext.buf;
      break;
  }

  /*--- Ecriture dans les images ---*/
  oldval=0;
  for (i=0 ; i<(int)(image_in->dim.x*image_in->dim.y*image_in->dim.z) ; i++ ) {
      switch (image_in->type)  {
      case UCHAR:
      case SCHAR:
          val=(int)inU8[i];
          break;
      case USHORT:
      case SSHORT:
          val=(int)inU16[i];
          break;
      default:
          VT_FreeImage(image_in);
          VT_FreeImage(image_ext);
          VT_FreeImage(&out_in);
          VT_FreeImage(&out_ext);
          _free_mat(&mat);
          VT_ErrorParse("error with out image types\n", 0);
      }
      if(i==0 || val != oldval) {
          if (val!=0 || par.flag_zeros != 1)
            ind=_ind(mat, 0, val);
          oldval=val;
      }
      if(val==0 && par.flag_zeros == 1) {
          switch (t) {
          case UCHAR:
              bufU8_in[i]=(unsigned char) 0;
              break;
          case USHORT:
              bufU16_in[i]=(unsigned short int) 0;
              break;
          default:
              VT_FreeImage(image_in);
              VT_FreeImage(image_ext);
              VT_FreeImage(&out_in);
              VT_FreeImage(&out_ext);
              _free_mat(&mat);
              VT_ErrorParse("error with out image types\n", 0);
          }
          continue;
      }
      switch (t) {
      case UCHAR:
          if(ind<0)
            bufU8_in[i]=(unsigned char) 255;
          else
            bufU8_in[i]=(unsigned char) mat.m[2+ind*3];
          break;
      case USHORT:
          if(ind<0)
            bufU16_in[i]=(unsigned short int) 32767;
          else
            bufU16_in[i]=(unsigned short int) mat.m[2+ind*3];
          break;
      default:
          VT_FreeImage(image_in);
          VT_FreeImage(image_ext);
          VT_FreeImage(&out_in);
          VT_FreeImage(&out_ext);
          _free_mat(&mat);
          VT_ErrorParse("error with out image types\n", 0);
      }
  }
  VT_FreeImage(image_in);


  for (i=0 ; i<(int)(image_ext->dim.x*image_ext->dim.y*image_ext->dim.z) ; i++ ) {
      switch (image_ext->type)  {
      case UCHAR:
      case SCHAR:
          val=(int)extU8[i];
          break;
      case USHORT:
      case SSHORT:
          val=(int)extU16[i];
          break;
      default:
          VT_FreeImage(image_ext);
          VT_FreeImage(&out_in);
          VT_FreeImage(&out_ext);
          _free_mat(&mat);
          VT_ErrorParse("error with out image types\n", 0);
      }
      if(i==0 || val != oldval) {
          if (val!=0 || par.flag_zeros != 1)
            ind=_ind(mat, 1, val);
          oldval=val;
      }
      if(val==0 && par.flag_zeros == 1) {
          switch (t) {
          case UCHAR:
              bufU8_ext[i]=(unsigned char) 0;
              break;
          case USHORT:
              bufU16_ext[i]=(unsigned short int) 0;
              break;
          default:
              VT_FreeImage(image_in);
              VT_FreeImage(image_ext);
              VT_FreeImage(&out_in);
              VT_FreeImage(&out_ext);
              _free_mat(&mat);
              VT_ErrorParse("error with out image types\n", 0);
          }
          continue;
      }
      switch (t) {
      case UCHAR:
          if(ind<0)
            bufU8_ext[i]=(unsigned char) 255;
          else
            bufU8_ext[i]=(unsigned char) mat.m[2+ind*3];
          break;
      case USHORT:
          if(ind<0)
            bufU16_ext[i]=(unsigned short int) 32767;
          else
            bufU16_ext[i]=(unsigned short int) mat.m[2+ind*3];
          break;
      default:
          VT_FreeImage(image_ext);
          VT_FreeImage(&out_in);
          VT_FreeImage(&out_ext);
          _free_mat(&mat);
          VT_ErrorParse("error with out image types\n", 0);
      }
  }
  VT_FreeImage(image_ext);
  _free_mat(&mat);


  if ( VT_WriteInrimage( &out_in ) == -1 ) {
    VT_FreeImage( &out_in );
    VT_FreeImage( &out_ext );
    VT_ErrorParse("unable to write output image-in\n", 0);
  }
  VT_FreeImage( &out_in );

  if ( VT_WriteInrimage( &out_ext ) == -1 ) {

    VT_FreeImage( &out_ext );
    VT_ErrorParse("unable to write output image-ext\n", 0);
  }
  VT_FreeImage( &out_ext );

  return( 0 );
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
          strcpy( par->name_labels.in, "<" );
          nb += 1;
        }
      }
      /*--- arguments generaux ---*/
      else if ( strcmp ( argv[i], "-help" ) == 0 ) {
        VT_ErrorParse("\n", 1);
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

      else if ( strcmp ( argv[i], "-z" ) == 0 ) {
        par->flag_zeros = 1;
      }

      else if ( strcmp ( argv[i], "-f" ) == 0 ) {
        par->flag_force_labels = 1;
      }

      else if ( strcmp ( argv[i], "-ref" ) == 0 ) {
        par->flag_ref = 1;
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
        strncpy( par->name_labels.in, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else if ( nb == 1 ) {
        strncpy( par->name_labels.ext, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else if ( nb == 2 ) {
        strncpy( par->name_labels.out, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else if ( nb == 3 ) {
        strncpy( par->name_outs.in, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else if ( nb == 4 ) {
        strncpy( par->name_outs.ext, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else if ( nb == 5 ) {
        strncpy( par->name_outs.out, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else {
        VT_ErrorParse("too much file names when parsing\n", 0 );
      }
    }
    i += 1;
  }

  /*--- s'il n'y a pas assez de noms ... ---*/
  if (nb < 4) {
      VT_ErrorParse("not enough file names when parsing\n", 0 );
  }


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
  VT_Names( &(par->name_labels) );
  VT_Names( &(par->name_outs) );
  par->flag_zeros = 0;
  par->flag_force_labels = 0;
  par->flag_ref = 0;
}




void _init_mat( _MATRIX *mat )
{
  mat->m = (int*)NULL;
  mat->l = 0;
  mat->c = 0;
}


int _alloc_mat ( _MATRIX *mat, int l, int c )
{

  mat->l = 0;
  mat->c = 0;

  mat->m =  (int *) calloc( (l*c), sizeof(int) );
  if ( mat->m == (int*)NULL ) return( -1 );

  mat->l = l;
  mat->c = c;

  return ( 1 );
}

void _free_mat ( _MATRIX *mat )
{
  if ( mat->m != (int*)NULL )
    free ( mat->m );

  mat->m = NULL;
  mat->l = 0;
  mat->c = 0;
}



int _read_mat( char *name, _MATRIX *m )
{
  int r;
  FILE *f;
  char line[512];
  int f1, f2;

  int n;

  if ((f = fopen (name, "r")) == NULL) {
    if ( _verbose_ )
      fprintf( stderr, "_read_mat: unable to open '%s' for reading\n", name );
    return( -1 );
  }

  /*
  if ( fgets(line, 512, f) == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "_read_mat: an error occurs when reading '%s'\n", name );
    return( -1 );
  }
  */
  /* counting the number of rows
   */


  n=0;

  while ( fgets(line, 512, f) != NULL ) {
      if ( ( sscanf( line, "%d %d\n", &f1, &f2) == 2 ) ) {
          n+=1;
      }
  }

  if ( n==0 ) {
    if ( _verbose_ )
      fprintf( stderr, "_read_mat: 0 row detected in '%s'\n", name );
    return( -1 );
  }
  fclose( f );




  if (_alloc_mat ( m, n, 3) != 1) {
      fprintf( stderr, "_read_mat: error while matrix allocation\n" );
      return(-1);
  }

  /* get the 2 first values, ie the first row
   */



  if ((f = fopen (name, "r")) == NULL) {
    if ( _verbose_ )
      fprintf( stderr, "_read_mat: unable to open '%s' for reading\n", name );
    return( -1 );
  }


  if ( fgets(line, 512, f) == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "_read_mat: an error occurs when reading '%s'\n", name );
    return( -1 );
  }

  while ( ( sscanf( line, "%d %d\n", &f1, &f2 ) < 2 )  ) {
    if ( fgets(line, 512, f) == NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "_read_mat: an error occurs when reading '%s'\n", name );
      return( -1 );
    }
  }


  m->m[0] = f1; m->m[1] = f2;

  r = 1;

  /* get the n-1 other rows
   */

  for ( r=1; r<n; r ++ ) {
    if ( fgets(line, 512, f) == NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "_read_mat: an error occurs when reading '%s'\n", name );
      return( -1 );
    }
    if ( sscanf( line, "%d %d\n", &f1, &f2 ) != 2 ) {
      if ( _verbose_ )
        fprintf( stderr, "_read_mat: line for row #%d does not contain 2 values in '%s'\n", r+1, name );
      return( -1 );
    }
    m->m[r*3+0] = f1;
    m->m[r*3+1] = f2;
  }

  fclose( f );
  return( 1 );
}




int _read_mat_3( char *name, _MATRIX *m )
{
  int r;
  FILE *f;
  char line[512];
  int f1, f2, f3;

  int n;

  if ((f = fopen (name, "r")) == NULL) {
    if ( _verbose_ )
      fprintf( stderr, "_read_mat_3: unable to open '%s' for reading\n", name );
    return( -1 );
  }

  /*
  if ( fgets(line, 512, f) == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "_read_mat: an error occurs when reading '%s'\n", name );
    return( -1 );
  }
  */
  /* counting the number of rows
   */


  n=0;

  while ( fgets(line, 512, f) != NULL ) {
      if ( ( sscanf( line, "%d %d %d\n", &f1, &f2, &f3) == 3 ) ) {
          n+=1;
      }
  }

  if ( n==0 ) {
    if ( _verbose_ )
      fprintf( stderr, "_read_mat_3: 0 row detected in '%s'\n", name );
    return( -1 );
  }
  fclose( f );




  if (_alloc_mat ( m, n, 3) != 1) {
      fprintf( stderr, "_read_mat_3: error while matrix allocation\n" );
      return(-1);
  }

  /* get the 2 first values, ie the first row
   */



  if ((f = fopen (name, "r")) == NULL) {
    if ( _verbose_ )
      fprintf( stderr, "_read_mat_3: unable to open '%s' for reading\n", name );
    return( -1 );
  }


  if ( fgets(line, 512, f) == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "_read_mat_3: an error occurs when reading '%s'\n", name );
    return( -1 );
  }

  while ( ( sscanf( line, "%d %d %d\n", &f1, &f2, &f3 ) < 3 )  ) {
    if ( fgets(line, 512, f) == NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "_read_mat_3: an error occurs when reading '%s'\n", name );
      return( -1 );
    }
  }


  m->m[0] = f1; m->m[1] = f2; m->m[2] = f3;

  r = 1;

  /* get the n-1 other rows
   */

  for ( r=1; r<n; r ++ ) {
    if ( fgets(line, 512, f) == NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "_read_mat_3: an error occurs when reading '%s'\n", name );
      return( -1 );
    }
    if ( sscanf( line, "%d %d %d\n", &f1, &f2, &f3 ) != 3 ) {
      if ( _verbose_ )
        fprintf( stderr, "_read_mat_3: line for row #%d does not contain 3 values in '%s'\n", r+1, name );
      return( -1 );
    }
    m->m[r*3+0] = f1;
    m->m[r*3+1] = f2;
    m->m[r*3+2] = f3;
  }

  fclose( f );
  return( 1 );
}



int _write_mat( char *name, _MATRIX *m )
{
  int i;
  FILE *f;
  /* char *format = "%d %d\n";
  */
  char *format = "%d \t%d \t---> %d\n";

  if ( name == NULL || name[0] == '\0'
       || (strncmp( name, "stderr", 6 ) == 0 && name[6] == '\0') ) {
    for (i = 0; i < m->l; i++) {
      fprintf( stderr, format,
           m->m[i*2], m->m[1+i*2] );
    }
  }
  else if ( strncmp( name, "stdout", 6 ) == 0 && name[6] == '\0' ) {
    for (i = 0; i < m->l; i++) {
      fprintf( stdout, format,
           m->m[i*2], m->m[1+i*2] );
    }
  }
  else {
    f = fopen( name, "w" );
    if ( f == NULL ) {
      if ( _verbose_ )
    fprintf( stderr, "_write_mat: unable to open '%s' for writing\n", name );
      return( -1 );
    }
    fprintf( f, "In \tExt \t---> New label\n" );
    for (i = 0; i < m->l; i++) {
      fprintf( f, format,
           m->m[i*3], m->m[1+i*3], m->m[2+i*3] );
    }
    fclose( f );
  }
  return( 1 );
}


void _compute_table ( _MATRIX *mat, int flag_ref )
{
    int i;

    for (i=0 ; i<mat->l ; i++) {
        if (flag_ref==0)
          mat->m[(i+1)*mat->c-1]=i+1;
        else
          mat->m[(i+1)*mat->c-1]=mat->m[i*mat->c];
    }
}

int _ind( _MATRIX mat, int c, int val)
{
    int i;
    for (i=0; i<mat.l; i++)
        if(mat.m[c+i*mat.c]==val)
            return(i);
    return(-1);
}
