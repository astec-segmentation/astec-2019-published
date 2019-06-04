/*****************************************************************************
 * extraction_arbre.cpp
 *
 * par Manon Linder
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "vt_common.h"
#include "D_graphe.h"
#include "G_classif.h"

typedef struct local_par 
{
  vt_names names;
  
  int writeImages;
  int res;
  int tronc;
  
} local_par;

using namespace std;

/*------- Definition des fonctions statiques ----------*/

static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( string str, int l );
static void VT_InitParam( local_par *par );

static string usage = "[image-skeleton-in] [image-Distance-in] [graphe-out]\n\
\t [-graph] [-tronc] \n\
\t [-help]";

static string detaillocal = "\
\t [image-skeleton-in] : Image du squelette caracterise \n\
\t [image-Distance-in] : Carte des distances \n\
\t [graphe-out] : Nom du fichier .txt contenant le graphe \n\
\t [-graph] : Construction du graphe \n\
\t [-tronc] : Recherche du tronc";

static char program[STRINGLENGTH];

/*------------------------------------------------------------------------------------------------------------*/

int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *imageSquelette;
  vt_image *imageDistance;
  vt_image imtmp;
	unsigned int rmax = 0;
	unsigned int a, b, c;
	unsigned char ***sk_array;
	unsigned char ***sk_array2;
	unsigned char ***dm_array; 
	char name[256];
	unsigned char *bufSquelette;
	unsigned char *bufTmp;
	int madimz;
	int madimy;
	int madimx;
	vector < int > mesdim;	
	
  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
 
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree 1 : Labels ---*/
  imageSquelette = _VT_Inrimage( par.names.in );
  if ( imageSquelette == (vt_image*)NULL ) 
  {
		cout << "unable to read input image Labels\n"<< endl;
		return(0);
	}

  /*--- lecture de l'image d'entree 2 : Distance ---*/
  imageDistance = _VT_Inrimage( par.names.ext );
  if ( imageDistance == (vt_image*)NULL ) 
  {
		cout << "unable to read input image Distance\n"<< endl;
		return(0);
	}

	/*--- verification meme taille pour les deux images d'entree ---*/
	if ( imageSquelette->dim.x != imageDistance->dim.x ||
       imageSquelette->dim.y != imageDistance->dim.y ||
       imageSquelette->dim.z != imageDistance->dim.z ) 
	{
		cout << "input images must have same dimensions\n"<< endl;
		return(0);
	}

/*------------------------------------------------------------------------------------------------------------*/
  /* 
	* ETAPE 1 : EXTRAIRE LES INFORMATIONS UTILES 
	*
	*/

	/*--- verification image LABELS ---*/	
	if (imageSquelette->type == UCHAR)
	{ 
		sk_array=(unsigned char***)imageSquelette->array;
	}
	else
	{
		cout << "input image Labels type not handled yet\n"<< endl;
		return(0);
	}
	
	/*--- verification image DISTANCE ---*/	
	if (imageDistance->type == UCHAR)
	{
		dm_array=(unsigned char***)imageDistance->array;
	}
	else
	{
		cout << "input image distance type not handled yet\n"<< endl;
		return(0);
	}

	/*--- Copie de l'image de squelette dans imtmp ---*/
	sprintf( name, "lipoutou.inr" );
	VT_Image( &imtmp );	//pre-initialisation
  VT_InitFromImage( &imtmp, imageSquelette, name, imageSquelette->type);
  if(VT_AllocImage( &imtmp ) != 1)
  {
  	VT_FreeImage(imageSquelette);
  	VT_FreeImage(imageDistance);
  	cout << "allocation of lipoutou.inr failed" << endl;
  	return( 0 );
  }	 
  switch (imageSquelette->type) 
  {
  	case UCHAR :
  		bufSquelette = (unsigned char *)imageSquelette->buf;
  		bufTmp = (unsigned char *)imtmp.buf;
			sk_array2=(unsigned char***)imtmp.array;
  		break;
 		default :
  		VT_FreeImage(imageSquelette);
  		VT_FreeImage(imageDistance);  	
  		VT_FreeImage(&imtmp);
  		cout << "type of lipoutou.inr not handled yet" << endl;
  		return( 0 );
  }
	for (a = 0 ; a < imtmp.dim.x*imtmp.dim.y*imtmp.dim.z ; a++) 
	{
		bufTmp[a] = bufSquelette[a];
	}
  
/*------------------------------------------------------------------------------------------------------------*/ 
	/* 
	 * ETAPE 2 : EXTRACTION ARBRE SIMPLE PAR LE TRONC
	 *
	 */
	   
	/*--- Point depart du tronc = extremite dont rayon le plus grand ---*/
  unsigned int z = 0;
  unsigned int y = 0;
  unsigned int x = 0;

	madimz = imageSquelette->dim.z;
  madimy = imageSquelette->dim.y;
  madimx = imageSquelette->dim.x;
  mesdim.push_back(madimz);
  mesdim.push_back(madimy);
  mesdim.push_back(madimx);
  
  for (a=0; a < imageSquelette->dim.z; a++)
	for (b=0; b < imageSquelette->dim.y; b++)
	for (c=0; c < imageSquelette->dim.x; c++)
	{
		if (sk_array[a][b][c]==(unsigned char)200 && dm_array[a][b][c]>=rmax) //Recherche extremites
		{
			rmax = dm_array[a][b][c];
			z=a;
			y=b;
			x=c;
		}
	}

  D_coord_ori Mescoord_ori(z,y,x,sk_array[z][y][x]);
	cout << "**** First Pixel ****" << endl;
	Mescoord_ori.AfficheOPixel();
	
	/*--- Creation d'un nouveau graphe ---*/
	D_graphe mon_graphe;
	D_node the_special_node;

	mon_graphe.SkeletonProcess(0, the_special_node, 1, Mescoord_ori, sk_array, dm_array, mesdim);
	
	if (par.res == 1)
	{
		cout << " Trouvez la branche qui correspond au tronc du reseau puis fermez la fenetre" << endl;
	
		mon_graphe.visualization3DLabel();
	
		cout << " Entrez le numero de la branche : ";
	
		int trunk = 0 ;
	
		cin >> trunk;
	
		// Vérification point origine
		if (sk_array2[z][y][x]!=(unsigned char)200)
		{
			cout << " Mauvaise Branche " << endl;
			VT_FreeImage( imageSquelette );
  		VT_FreeImage( imageDistance );
  		VT_FreeImage( &imtmp );
			return (0);
		}
		if (trunk != 0)
		{
			vector<D_node> l_node;
			l_node = mon_graphe.getListeNode();
			z = l_node[trunk].getValuePixel(l_node[trunk].getSizeNode()-1).getCoordZ();
			y = l_node[trunk].getValuePixel(l_node[trunk].getSizeNode()-1).getCoordY();
			x = l_node[trunk].getValuePixel(l_node[trunk].getSizeNode()-1).getCoordX();

			if (sk_array2[z][y][x]!=(unsigned char)200)
			{
				cout << " Mauvaise Branche " << endl;
				VT_FreeImage( imageSquelette );
  			VT_FreeImage( imageDistance );
  			VT_FreeImage( &imtmp );
				return (0);
			}
			else
			{
				D_coord_ori mon_depart(z,y,x,sk_array2[z][y][x]);
				cout << "**** The Pixel ****" << endl;
				mon_depart.AfficheOPixel();
				D_graphe mon_graphe;
				D_node the_special_node;
				mon_graphe.SkeletonProcess(0, the_special_node, 1, mon_depart, sk_array2, dm_array, mesdim);
				mon_graphe.saveGraphe(par.names.out);
			}
		}
		else
		{
			cout << "**** The Pixel ****" << endl;
			Mescoord_ori.AfficheOPixel();
			D_graphe mon_graphe;
			D_node the_special_node;
			mon_graphe.SkeletonProcess(0, the_special_node, 1, Mescoord_ori, sk_array2, dm_array, mesdim);
			mon_graphe.saveGraphe(par.names.out);
		}
	}
	
	if (par.tronc == 1)
	{
		D_coord_ori mon_depart;
		mon_depart = mon_graphe.SearchTrueDepartureMixAll(sk_array2);
	
		cout << "**** A Pixel ****" << endl;
		mon_depart.AfficheOPixel();
		mon_graphe.visualizationTag();
	
	
		//if(mon_depart.getOCoordZ() != z && mon_depart.getOCoordY() != y && mon_depart.getOCoordX()!= x)
		//{
		//	D_graphe mon_graphe;
		//	
		//	D_node the_special_node;
		//
		//	mon_graphe.SkeletonProcess(0, the_special_node, 1, mon_depart, sk_array2, dm_array);
		//	mon_graphe.saveGraphe();
		//}
		//else
		//{
		//	mon_graphe.saveGraphe();
		//}
	}

	
			
/*------------------------------------------------------------------------------------------------------------*/ 
	

  /*--- liberations memoires ---*/
  VT_FreeImage( imageSquelette );
  VT_FreeImage( imageDistance );
  VT_FreeImage( &imtmp );
  return( 1 );
}

/*------------------------------------------------------------------------------------------------------------*/

static void VT_Parse( int argc, char *argv[], local_par *par )
{
	int i, nb;
	char text[STRINGLENGTH];
  
	if ( VT_CopyName( program, argv[0] ) != 1 )
	{
   cout << "Error while copying program name " << (char*)NULL << endl;
	}
	if ( argc == 1 ) VT_ErrorParse("\n", 0 );
  
  /*--- lecture des parametres ---*/
	i = 1; nb = 0;
	while ( i < argc ) 
	{
    if ( argv[i][0] == '-' ) 
    {
      if ( argv[i][1] == '\0' ) 
      {
				if ( nb == 0 ) 
				{
					/*--- standart input ---*/
					strcpy( par->names.in, "<" );
					nb += 1;
				}
      }
      /*--- arguments generaux ---*/
      else if ( strcmp ( argv[i], "-help" ) == 0 ) 
      {
				VT_ErrorParse("\n", 1);
      }
      
      /*--- options ---*/
			else if ( strcmp ( argv[i], "-graph" ) == 0 ) 
      {
				par->res = 1;
			}
			else if ( strcmp ( argv[i], "-tronc" ) == 0 ) 
      {
				par->tronc = 1;
			}
      
      /*--- option inconnue ---*/
      else 
      {
				sprintf(text,"unknown option %s\n",argv[i]);
				VT_ErrorParse(text, 0);
      }
    }
    
		/*--- saisie des noms d'images ---*/
    else if ( argv[i][0] != 0 ) 
    {
      if ( nb == 0 ) 
      { 
				strncpy( par->names.in, argv[i], STRINGLENGTH );  
				nb += 1;
      }
      else if ( nb == 1 ) 
      {
				strncpy( par->names.ext, argv[i], STRINGLENGTH );  
				nb += 1;
      }
      else if ( nb == 2 ) 
      {
				strncpy( par->names.out, argv[i], STRINGLENGTH );  
				nb += 1;
      }
      else 
				VT_ErrorParse("too much file names when parsing\n", 0 );
    }
    i += 1;
	}
  
	/*--- s'il n'y a pas assez de noms ... ---*/
	if (nb == 0) 
	{
    strcpy( par->names.in,  "<" );  /* standart input */
    strcpy( par->names.out, ">" );  /* standart output */
  }
  if (nb == 1)
    strcpy( par->names.out, ">" );  /* standart output */
}


static void VT_ErrorParse( string str, int flag )
{
	cerr << "Usage : "<< program << " " << usage << endl;
  if ( flag == 1 )
  {
  	cerr << detaillocal << endl;
  }
  cerr << "Erreur : "<< str << endl;
  exit(0);
}


static void VT_InitParam( local_par *par )
{
  VT_Names( &(par->names) );

  par->writeImages = 0;
  par->res = 0;
  par->tronc = 0;
}

