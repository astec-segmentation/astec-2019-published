/*************************************************************************
 * labelBorders.c -
 *
 * $Id: labelBorders.c,v 1.0 2013/08/12 10:34:51 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2013/08/12
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>

#define TABLENGTH 2000

typedef struct local_par {
  vt_names names;
  int gauches[TABLENGTH];
  int droites[TABLENGTH];
  int N;
  int c;
} local_par;



/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );

static int _verbose_ = 0;



static int findInVec(int val, int *vec, int n);

static char *usage = "[image-in] [image-out]\n\
\t [-p |-pair %d %d [%d %d [...]]] [-cross | -c] [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -p %d %d [...] : restreint les bords aux interfaces entre les paires de labels stipulees\n\
\t -cross : restreint les bords aux interfaces entre 3 ou plus labels\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{

  local_par par;
  vt_image *image;
  vt_image imres;
  char name[DOUBLESTRINGLENGTH];
  int z,y,x;
  int i, j, g, d, a[3], b, l, cpt;
  int X[]={1, 0, 1, 0, 1, 0, 1};
  int Y[]={0, 1, 1, 0, 0, 1, 1};
  int Z[]={0, 0, 0, 1, 1, 1, 1};
  int ind_X[]={1, 3, 5};
  int ind_Y[]={0, 3, 4};
  int ind_Z[]={0, 1, 2};

  /*--- initialisation des parametres ---*/
  MT_InitParam( &par );

  /*--- lecture des parametres ---*/
  MT_Parse( argc, argv, &par );

  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL )
    MT_ErrorParse("unable to read input image\n", 0);

  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.names.inv == 1 )  VT_InverseImage( image );
  if ( par.names.swap == 1 ) VT_SwapImage( image );





  /* Alloc image sortie */

  sprintf( name, "%s", par.names.out );
  VT_InitImage( &imres, name, image->dim.x, image->dim.y,
              image->dim.z, (int)UCHAR );
  if ( VT_AllocImage( &(imres) ) != 1 ) {
    MT_ErrorParse("problem while allocating imres\n", 0 );
    VT_FreeImage( image );
    return( -1 );
  }
  imres.siz.x=image->siz.x;
  imres.siz.y=image->siz.y;
  imres.siz.z=image->siz.z;

  /*--- calculs ---*/

  unsigned char ***arrayu8=NULL;
  unsigned short int ***arrayu16=NULL;
  unsigned char ***arrayres;

  switch (image->type) {
    case SCHAR:
    case UCHAR:
      arrayu8=(unsigned char ***)image->array;
      break;
    case SSHORT:
    case USHORT:
      arrayu16=(unsigned short int ***)image->array;
      break;
    case TYPE_UNKNOWN:
    default:
      VT_Error("image type unknown or not supported for this program",program);
      return( 0 );
  }
  arrayres=(unsigned char ***)imres.array;

  for(z=0;z<(int)image->dim.z;z++)
  for(y=0;y<(int)image->dim.y;y++)
  for(x=0;x<(int)image->dim.x;x++)
    arrayres[z][y][x]=(unsigned char)0;

  switch (image->type) {
    case SCHAR:
    case UCHAR:
      for(z=0;z<(int)image->dim.z-1;z++)
      for(y=0;y<(int)image->dim.y-1;y++)
      for(x=0;x<(int)image->dim.x-1;x++) {
		if (par.N==0 && par.c==0) {
          if(arrayu8[z][y][x]!=arrayu8[z][y][x+1]) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y][x+1]=(unsigned char)255;
          }
          if(arrayu8[z][y][x]!=arrayu8[z][y+1][x]) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y+1][x]=(unsigned char)255;
		  }
          if(arrayu8[z][y][x]!=arrayu8[z+1][y][x]) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z+1][y][x]=(unsigned char)255;
          }
		  continue;
		}
		if(par.N>0 && par.c==0) {
		  g=arrayu8[z][y][x];
		  i=findInVec((int)g, par.gauches, par.N);
		  if (i<0) {
			i=findInVec((int)g, par.droites, par.N);
			if (i<0) continue;
			d=par.gauches[i];
		  }
		  else {
			d=par.droites[i];
		  }
          if(arrayu8[z][y][x+1]==(unsigned char)d) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y][x+1]=(unsigned char)255;
          }
          if(arrayu8[z][y+1][x]==(unsigned char)d) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y+1][x]=(unsigned char)255;
		  }
          if(arrayu8[z+1][y][x]==(unsigned char)d) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z+1][y][x]=(unsigned char)255;
          }
		  continue;
		}
		if(par.N==0 && par.c==1)
		{
		  a[0]=(int)arrayu8[z][y][x];
		  cpt=1;
		  for (i=0; i<7; i++) {
			if (cpt>=3) break;
			b=(int)arrayu8[z+Z[i]][y+Y[i]][x+X[i]];
			l=1; 
			for (j=0; j<cpt; j++) 
			  if (b==a[j]) l=0;
			if (l==1) {
			  a[cpt]=b;
			  cpt+=1;
			}
		  }
		  if(cpt>=3) {
			arrayres[z][y][x]=(unsigned char)255;
			for (i=0; i<7; i++) {
			  arrayres[z+Z[i]][y+Y[i]][x+X[i]]=(unsigned char)255;
			}
		  }
		  continue;
		}
		MT_ErrorParse("Uncompatible parameters while parsing.\n", 0);
      }
      z=image->dim.z-1;
      for(y=0;y<(int)image->dim.y-1;y++)
      for(x=0;x<(int)image->dim.x-1;x++)  {
		if (par.N == 0 && par.c==0) {
          if(arrayu8[z][y][x]!=arrayu8[z][y][x+1]) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y][x+1]=(unsigned char)255;
          }
          if(arrayu8[z][y][x]!=arrayu8[z][y+1][x]) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y+1][x]=(unsigned char)255;
		  }
		  continue;
        }
        if(par.N>0 && par.c==0) {
		  g=arrayu8[z][y][x];
		  i=findInVec((int)g, par.gauches, par.N);
		  if (i<0) {
			i=findInVec((int)g, par.droites, par.N);
			if (i<0) continue;
			d=par.gauches[i];
		  }
		  else {
			d=par.droites[i];
		  }
          if(arrayu8[z][y][x+1]==(unsigned char)d) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y][x+1]=(unsigned char)255;
          }
          if(arrayu8[z][y+1][x]==(unsigned char)d) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y+1][x]=(unsigned char)255;
		  }
		  continue;
		}
		if(par.N==0 && par.c==1)
		{
		  a[0]=(int)arrayu8[z][y][x];
		  cpt=1;
		  for (i=0; i<3; i++) {
			if (cpt>=3) break;
			b=(int)arrayu8[z+Z[ind_Z[i]]][y+Y[ind_Z[i]]][x+X[ind_Z[i]]];
			l=1; 
			for (j=0; j<cpt; j++) 
			  if (b==a[j]) l=0;
			if (l==1) {
			  a[cpt]=b;
			  cpt+=1;
			}
		  }
		  if(cpt>=3) {
			arrayres[z][y][x]=(unsigned char)255;
			for (i=0; i<3; i++) {
			  arrayres[z+Z[ind_Z[i]]][y+Y[ind_Z[i]]][x+X[ind_Z[i]]]=(unsigned char)255;
			}
		  }
		  continue;
		}
		MT_ErrorParse("Uncompatible parameters while parsing.\n", 0);
      }
      y=image->dim.y-1;
      for(z=0;z<(int)image->dim.z-1;z++)
      for(x=0;x<(int)image->dim.x-1;x++)  {
		if (par.N == 0 && par.c==0) {
          if(arrayu8[z][y][x]!=arrayu8[z][y][x+1]) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y][x+1]=(unsigned char)255;
          }
          if(arrayu8[z][y][x]!=arrayu8[z+1][y][x]) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z+1][y][x]=(unsigned char)255;
          }
		  continue;
		}
		if(par.N>0 && par.c==0) {
		  g=arrayu8[z][y][x];
		  i=findInVec((int)g, par.gauches, par.N);
		  if (i<0) {
			i=findInVec((int)g, par.droites, par.N);
			if (i<0) continue;
			d=par.gauches[i];
		  }
		  else {
			d=par.droites[i];
		  }
          if(arrayu8[z][y][x+1]==(unsigned char)d) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y][x+1]=(unsigned char)255;
          }
          if(arrayu8[z+1][y][x]==(unsigned char)d) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z+1][y][x]=(unsigned char)255;
          }
		  continue;
		}
		if(par.N==0 && par.c==1)
		{
		  a[0]=(int)arrayu8[z][y][x];
		  cpt=1;
		  for (i=0; i<3; i++) {
			if (cpt>=3) break;
			b=(int)arrayu8[z+Z[ind_Y[i]]][y+Y[ind_Y[i]]][x+X[ind_Y[i]]];
			l=1; 
			for (j=0; j<cpt; j++) 
			  if (b==a[j]) l=0;
			if (l==1) {
			  a[cpt]=b;
			  cpt+=1;
			}
		  }
		  if(cpt>=3) {
			arrayres[z][y][x]=(unsigned char)255;
			for (i=0; i<3; i++) {
			  arrayres[z+Z[ind_Y[i]]][y+Y[ind_Y[i]]][x+X[ind_Y[i]]]=(unsigned char)255;
			}
		  }
		  continue;
		}
		MT_ErrorParse("Uncompatible parameters while parsing.\n", 0);
      }
      x=image->dim.x-1;
      for(z=0;z<(int)image->dim.z-1;z++)
      for(y=0;y<(int)image->dim.y-1;y++)  {
		if (par.N == 0 && par.c==0) {
		  if(arrayu8[z][y][x]!=arrayu8[z][y+1][x]) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y+1][x]=(unsigned char)255;
          }
          if(arrayu8[z][y][x]!=arrayu8[z+1][y][x]) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z+1][y][x]=(unsigned char)255;
          }
		  continue;
		}
		if(par.N>0 && par.c==0) {
		  g=arrayu8[z][y][x];
		  i=findInVec((int)g, par.gauches, par.N);
		  if (i<0) {
			i=findInVec((int)g, par.droites, par.N);
			if (i<0) continue;
			d=par.gauches[i];
		  }
		  else {
			d=par.droites[i];
		  }
          if(arrayu8[z][y+1][x]==(unsigned char)d) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y+1][x]=(unsigned char)255;
		  }
          if(arrayu8[z+1][y][x]==(unsigned char)d) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z+1][y][x]=(unsigned char)255;
          }
		  continue;
		}
		if(par.N==0 && par.c==1)
		{
		  a[0]=(int)arrayu8[z][y][x];
		  cpt=1;
		  for (i=0; i<3; i++) {
			if (cpt>=3) break;
			b=(int)arrayu8[z+Z[ind_X[i]]][y+Y[ind_X[i]]][x+X[ind_X[i]]];
			l=1; 
			for (j=0; j<cpt; j++) 
			  if (b==a[j]) l=0;
			if (l==1) {
			  a[cpt]=b;
			  cpt+=1;
			}
		  }
		  if(cpt>=3) {
			arrayres[z][y][x]=(unsigned char)255;
			for (i=0; i<3; i++) {
			  arrayres[z+Z[ind_X[i]]][y+Y[ind_X[i]]][x+X[ind_X[i]]]=(unsigned char)255;
			}
		  }
		  continue;
		}
		MT_ErrorParse("Uncompatible parameters while parsing.\n", 0);
      }
      break;
    case SSHORT:
    case USHORT:
      for(z=0;z<(int)image->dim.z-1;z++)
      for(y=0;y<(int)image->dim.y-1;y++)
      for(x=0;x<(int)image->dim.x-1;x++) {
		if (par.N == 0 && par.c==0) {
          if(arrayu16[z][y][x]!=arrayu16[z][y][x+1]) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y][x+1]=(unsigned char)255;
          }
          if(arrayu16[z][y][x]!=arrayu16[z][y+1][x]) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y+1][x]=(unsigned char)255;
          }
          if(arrayu16[z][y][x]!=arrayu16[z+1][y][x]) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z+1][y][x]=(unsigned char)255;
		  }
		  continue;
		}
		if(par.N>0 && par.c==0) {
		  g=arrayu16[z][y][x];
		  i=findInVec((int)g, par.gauches, par.N);
		  if (i<0) {
			i=findInVec((int)g, par.droites, par.N);
			if (i<0) continue;
			d=par.gauches[i];
		  }
		  else {
			d=par.droites[i];
		  }
          if(arrayu16[z][y][x+1]==(unsigned short int)d) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y][x+1]=(unsigned char)255;
          }
          if(arrayu16[z][y+1][x]==(unsigned short int)d) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y+1][x]=(unsigned char)255;
		  }
          if(arrayu16[z+1][y][x]==(unsigned short int)d) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z+1][y][x]=(unsigned char)255;
          }
		  continue;
		}
		if(par.N==0 && par.c==1)
		{
		  a[0]=(int)arrayu16[z][y][x];
		  cpt=1;
		  for (i=0; i<7; i++) {
			if (cpt>=3) break;
			b=(int)arrayu16[z+Z[i]][y+Y[i]][x+X[i]];
			l=1; 
			for (j=0; j<cpt; j++) 
			  if (b==a[j]) l=0;
			if (l==1) {
			  a[cpt]=b;
			  cpt+=1;
			}
		  }
		  if(cpt>=3) {
			arrayres[z][y][x]=(unsigned char)255;
			for (i=0; i<7; i++) {
			  arrayres[z+Z[i]][y+Y[i]][x+X[i]]=(unsigned char)255;
			}
		  }
		  continue;
		}
		MT_ErrorParse("Uncompatible parameters while parsing.\n", 0);
      }
      z=image->dim.z-1;
      for(y=0;y<(int)image->dim.y-1;y++)
      for(x=0;x<(int)image->dim.x-1;x++)  {
		if (par.N == 0 && par.c==0) {
          if(arrayu16[z][y][x]!=arrayu16[z][y][x+1]) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y][x+1]=(unsigned char)255;
          }
          if(arrayu16[z][y][x]!=arrayu16[z][y+1][x]) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y+1][x]=(unsigned char)255;
          }
		  continue;
		}
		if(par.N>0 && par.c==0) {
		  g=arrayu16[z][y][x];
		  i=findInVec((int)g, par.gauches, par.N);
		  if (i<0) {
			i=findInVec((int)g, par.droites, par.N);
			if (i<0) continue;
			d=par.gauches[i];
		  }
		  else {
			d=par.droites[i];
		  }
          if(arrayu16[z][y][x+1]==(unsigned short int)d) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y][x+1]=(unsigned char)255;
          }
          if(arrayu16[z][y+1][x]==(unsigned short int)d) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y+1][x]=(unsigned char)255;
		  }
		  continue;
		}
		if(par.N==0 && par.c==1)
		{
		  a[0]=(int)arrayu16[z][y][x];
		  cpt=1;
		  for (i=0; i<3; i++) {
			if (cpt>=3) break;
			b=(int)arrayu16[z+Z[ind_Z[i]]][y+Y[ind_Z[i]]][x+X[ind_Z[i]]];
			l=1; 
			for (j=0; j<cpt; j++) 
			  if (b==a[j]) l=0;
			if (l==1) {
			  a[cpt]=b;
			  cpt+=1;
			}
		  }
		  if(cpt>=3) {
			arrayres[z][y][x]=(unsigned char)255;
			for (i=0; i<3; i++) {
			  arrayres[z+Z[ind_Z[i]]][y+Y[ind_Z[i]]][x+X[ind_Z[i]]]=(unsigned char)255;
			}
		  }
		  continue;
		}
		MT_ErrorParse("Uncompatible parameters while parsing.\n", 0);
      }
      y=image->dim.y-1;
      for(z=0;z<(int)image->dim.z-1;z++)
      for(x=0;x<(int)image->dim.x-1;x++)  {
        if( par.N == 0 && par.c==0 ) {
          if(arrayu16[z][y][x]!=arrayu16[z][y][x+1]) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y][x+1]=(unsigned char)255;
          }
          if(arrayu16[z][y][x]!=arrayu16[z+1][y][x]) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z+1][y][x]=(unsigned char)255;
          }
		  continue;
		}
		if(par.N>0 && par.c==0) {
		  g=arrayu16[z][y][x];
		  i=findInVec((int)g, par.gauches, par.N);
		  if (i<0) {
			i=findInVec((int)g, par.droites, par.N);
			if (i<0) continue;
			d=par.gauches[i];
		  }
		  else {
			d=par.droites[i];
		  }
          if(arrayu16[z][y][x+1]==(unsigned short int)d) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y][x+1]=(unsigned char)255;
          }
          if(arrayu16[z+1][y][x]==(unsigned short int)d) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z+1][y][x]=(unsigned char)255;
          }
		  continue;
		}
		if(par.N==0 && par.c==1)
		{
		  a[0]=(int)arrayu16[z][y][x];
		  cpt=1;
		  for (i=0; i<3; i++) {
			if (cpt>=3) break;
			b=(int)arrayu16[z+Z[ind_Y[i]]][y+Y[ind_Y[i]]][x+X[ind_Y[i]]];
			l=1; 
			for (j=0; j<cpt; j++) 
			  if (b==a[j]) l=0;
			if (l==1) {
			  a[cpt]=b;
			  cpt+=1;
			}
		  }
		  if(cpt>=3) {
			arrayres[z][y][x]=(unsigned char)255;
			for (i=0; i<3; i++) {
			  arrayres[z+Z[ind_Y[i]]][y+Y[ind_Y[i]]][x+X[ind_Y[i]]]=(unsigned char)255;
			}
		  }
		  continue;
		}
		MT_ErrorParse("Uncompatible parameters while parsing.\n", 0);
      }
      x=image->dim.x-1;
      for(z=0;z<(int)image->dim.z-1;z++)
      for(y=0;y<(int)image->dim.y-1;y++)  {
		if(par.N == 0 && par.c==0) {
          if(arrayu16[z][y][x]!=arrayu16[z][y+1][x]) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y+1][x]=(unsigned char)255;
          }
          if(arrayu16[z][y][x]!=arrayu16[z+1][y][x]) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z+1][y][x]=(unsigned char)255;
          }
		  continue;
		}
		if(par.N>0 && par.c==0) {
		  g=arrayu16[z][y][x];
		  i=findInVec((int)g, par.gauches, par.N);
		  if (i<0) {
			i=findInVec((int)g, par.droites, par.N);
			if (i<0) continue;
			d=par.gauches[i];
		  }
		  else {
			d=par.droites[i];
		  }
          if(arrayu16[z][y+1][x]==(unsigned short int)d) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z][y+1][x]=(unsigned char)255;
		  }
          if(arrayu16[z+1][y][x]==(unsigned short int)d) {
			arrayres[z][y][x]=(unsigned char)255;
			arrayres[z+1][y][x]=(unsigned char)255;
          }
		  continue;
		}
		if(par.N==0 && par.c==1)
		{
		  a[0]=(int)arrayu16[z][y][x];
		  cpt=1;
		  for (i=0; i<3; i++) {
			if (cpt>=3) break;
			b=(int)arrayu16[z+Z[ind_X[i]]][y+Y[ind_X[i]]][x+X[ind_X[i]]];
			l=1; 
			for (j=0; j<cpt; j++) 
			  if (b==a[j]) l=0;
			if (l==1) {
			  a[cpt]=b;
			  cpt+=1;
			}
		  }
		  if(cpt>=3) {
			arrayres[z][y][x]=(unsigned char)255;
			for (i=0; i<3; i++) {
			  arrayres[z+Z[ind_X[i]]][y+Y[ind_X[i]]][x+X[ind_X[i]]]=(unsigned char)255;
			}
		  }
		  continue;
		}
		MT_ErrorParse("Uncompatible parameters while parsing.\n", 0);
      }
      break;
    case TYPE_UNKNOWN:
    default:
      VT_Error("image type unknown or not supported for this program",program);
      return( 0 );
  }

  VT_WriteInrimage( &(imres) );

  /*--- liberations memoires ---*/
  VT_FreeImage( image );
  VT_Free( (void**)&image );
  VT_FreeImage( &imres );

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

      else if ( strcmp ( argv[i], "-p" ) == 0 || strcmp ( argv[i], "-pair" ) == 0 ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -p...\n", 0 );
		status = sscanf( argv[i],"%d",&(par->gauches[par->N]) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -p...\n", 0 );
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -p...\n", 0 );
		status = sscanf( argv[i],"%d",&(par->droites[par->N]) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -p...\n", 0 );
		par->N += 1;
		i += 1;
		while ( i < argc)  {
		  status = sscanf( argv[i],"%d",&(par->gauches[par->N]) );
		  if ( status <= 0 ) break;
		  i += 1;
		  if ( i >= argc)    MT_ErrorParse( "parsing -p...\n", 0 );
		  status = sscanf( argv[i],"%d",&(par->droites[par->N]) );
		  if ( status <= 0 ) MT_ErrorParse( "parsing -p...\n", 0 );
		  par->N += 1;
		  i += 1;
	    }
      }

      else if ( strcmp ( argv[i], "-cross" ) == 0 || strcmp ( argv[i], "-c" ) == 0 ) {
        par->c = 1;
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
        strncpy( par->names.out, argv[i], STRINGLENGTH );
        nb += 1;
      }
      else
        MT_ErrorParse("too much file names when parsing\n", 0 );
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
  par->N=0;
  par->c=0;
}


static int findInVec(int val, int *vec, int n)
{
  int i;
  for (i=0; i<n ; i++)
    if (vec[i]==val) return(i);
  return (-1);
}
