/*****************************************************************************
 * extraction_arbre.c
 *
 * par Manon Linder
 *
 *
 */

#include "D_graphe.h"
#include "G_classif.h"


using namespace std;

typedef struct local_par {
  vt_names names;
  int vox_size;
	int matrice_connexite;
	int di_ord;
	int lo_ord;
	int aro_ord;
	int abocassot_ord;
	int abo_ord;
	int all_res;
	int viz_class;
	int viz_tag;

} local_par;

/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( string str, int l );
static void VT_InitParam( local_par *par );


static string usage = "[graphe-in] [result-out] [-vox %d]\n\
\t [-cm] [-do] [-lo] [-aro] [-abocassot] [-abo] [-all]\n\
\t [-viz] [-vitag]\n\
\t [-help]";

static string detaillocal = "\
\t [graphe-in] : Graphe en .txt \n\
\t [result-out] : Fichier des resultats en .txt \n\
\t [-vox %d] : taille du voxel en micro-m \n\
\t -cm : matrice de connexite\n\
\t -do : diametre par ordre\n\
\t -lo : longueur par ordre\n\
\t -aro : ratio aire par ordre\n\
\t -abocassot : ratio aire par ordre suivant la methode de Cassot\n\
\t -abo : angle bifurcation par ordre\n\
\t -all : enregistre tous les resultats\n\
\t -viz : visualisation\n\
\t -vitag : visualisation tag";

static char program[STRINGLENGTH];

////////////////////////////////////////////////////////////////////////////////////////////////
int main( int argc, char *argv[] )
{
	local_par par;
	
	/*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
	if (par.vox_size == 0)
	{
		cout << "Il manque la taille des voxels"<< endl;
		return(0);
	}
	
	/*--- Creation d'un nouveau graphe ---*/
	G_classif mon_graphe;
	
	mon_graphe.openGraphe(par.names.in);

	mon_graphe.getConnexiteType(par.names.out);

	mon_graphe.Strahler_Diameter(par.names.out, par.vox_size);

	if(par.matrice_connexite==1 || par.all_res == 1)
	{
		mon_graphe.Connectivity_matrix(par.names.out);
	}
	
	if(par.di_ord==1 || par.all_res == 1)
	{
		mon_graphe.Diameter_order("Save", par.names.out);
	}
	
	if(par.lo_ord==1 || par.all_res == 1)
	{
		mon_graphe.Lenght_order(par.names.out);
	}
	
	if(par.aro_ord==1 || par.all_res == 1)
	{
		mon_graphe.Area_ratio_order(par.names.out);
	}
	
	if(par.abocassot_ord==1)
	{
		mon_graphe.Abo_Cassot(par.names.out, par.vox_size);
	}
	
	if(par.abo_ord==1 || par.all_res == 1)
	{
		mon_graphe.Angle_Bifurcation(par.names.out, par.vox_size);
	}
	
	if(par.viz_class==1)
	{
		mon_graphe.visualization();
	}
	
	if(par.viz_tag==1)
	{
		mon_graphe.visualizationTag();
	}	
	
	return(1);
}

////////////////////////////////////////////////////////////////////////////////////////////////

static void VT_Parse( int argc, char *argv[], local_par *par )
{
	int i, nb, status;
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
			else if ( strcmp ( argv[i], "-cm" ) == 0 ) 
      {
				par->matrice_connexite = 1;
			}
			else if ( strcmp ( argv[i], "-do" ) == 0 ) 
      {
				par->di_ord = 1;
			}
			else if ( strcmp ( argv[i], "-lo" ) == 0 ) 
      {
				par->lo_ord = 1;
			}
			else if ( strcmp ( argv[i], "-aro" ) == 0 ) 
      {
				par->aro_ord = 1;
			}
			else if ( strcmp ( argv[i], "-abocassot" ) == 0 ) 
      {
				par->abocassot_ord = 1;
			}
			else if ( strcmp ( argv[i], "-abo" ) == 0 ) 
      {
				par->abo_ord = 1;
			}
			else if ( strcmp ( argv[i], "-all" ) == 0 ) 
      {
				par->all_res = 1;
			}
			else if ( strcmp ( argv[i], "-viz" ) == 0 ) 
      {
				par->viz_class = 1;
			}
			else if ( strcmp ( argv[i], "-vitag" ) == 0 ) 
      {
				par->viz_tag = 1;
			}
			
			else if ( strcmp ( argv[i], "-cm" ) == 0 ) 
      {
				par->matrice_connexite = 1;
			}
			
		  else if ( strcmp ( argv[i], "-vox" ) == 0 ) 
		  {
				i += 1;
				if ( i >= argc)
				{
					VT_ErrorParse( "parsing -nb...\n", 0 );
				}
				status = sscanf( argv[i],"%d",&(par->vox_size) );
				if ( status <= 0 ) 
				{
					VT_ErrorParse( "parsing -nb...\n", 0 );
				}
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
  par->vox_size = 0;
	par->matrice_connexite = 0;
	par->di_ord = 0;
	par->lo_ord = 0;
	par->aro_ord = 0;
	par->abocassot_ord = 0;
	par->abo_ord = 0;
	par->all_res = 0;
	par->viz_class = 0;
	par->viz_tag = 0;
}
