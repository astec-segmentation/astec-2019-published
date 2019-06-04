/*****************************************************************************
* G_classif.cpp
*
* par Manon Linder
*
*
*/

#include "G_classif.h"


using namespace std;

/*------- CLASS G_CLASSIF ----------*/

G_classif::G_classif() : D_graphe(), connexite(0), Maxorder(0), ordre(0), diametre_moy_node(0), VAR_diametre_moy_node(0), diametre_moy_order(0), SD_diametre_moy_order(0), connectivity_matrix(0), lenght_moy_node(0), lenght_moy_order(0), area_node(0), 	total_ratio(0), major_ratio(0), minor_ratio(0), ass_ratio(0), angle_node(0), angle_m1(0), angle_m2(0), angle_12(0), angle_beta(0), angle_omega(0)
{

}

// Statistiques
void G_classif::getConnexiteType(char* name)
{
	getListeJct();
	int compte = 0;
	std::vector<int> con_0(0);
	std::vector<int> con_1(0);
	std::vector<int> con_2(0);
	std::vector<int> con_3(0);
	std::vector<int> con_4(0);
	std::vector<int> con_5(0);
	
	for(unsigned int j = 0; j < l_node_jct.size(); j++)
	{		
		for(int i = 0; i < getNumberEdge(); i++)
		{
			if (getCoupleNode(i).getValueNode1() == l_node_jct[j])
			{
				compte++;
			}
		}
		switch (compte)
		{
			case 1:
				con_0.push_back(l_node_jct[j]);
				break;
			case 2:
				con_1.push_back(l_node_jct[j]);
				break;
			case 3:
				con_2.push_back(l_node_jct[j]);
				break;
			case 4:
				con_3.push_back(l_node_jct[j]);
				break;
			case 5:
				con_4.push_back(l_node_jct[j]);
				break;
			case 6:
				con_5.push_back(l_node_jct[j]);
				break;
		}
		compte = 0;
	}
	connexite.push_back(con_0);
	connexite.push_back(con_1);
	connexite.push_back(con_2);
	connexite.push_back(con_3);
	connexite.push_back(con_4);
	connexite.push_back(con_5);
	
	/*cout << "-- CONNEXITE --" << endl;
	cout << "Faux noeud fin segment : " << connexite[0].size() << endl;
	cout << "Faux noeud entre segments : " << connexite[1].size() << endl;
	cout << "Bifurcation : " << connexite[2].size() << endl;
	cout << "Trifurcation : " << connexite[3].size() << endl;
	cout << "Quadrifurcation : " << connexite[4].size() << endl;
	cout << "Pentafurcation : " << connexite[5].size() << endl;*/

	//Sauvegarde
	ofstream monFlux(name);
	
	if (monFlux)
	{
		// CONNEXITE //
		monFlux << "RESULTATS : " << name << endl;
		monFlux << "-- TYPE de jonction --" << endl;
		monFlux << "Faux noeud fin segment : " << connexite[0].size() << endl;
		monFlux << "Faux noeud entre segments : " << connexite[1].size() << endl;
		monFlux << "Bifurcation : " << connexite[2].size() << endl;
		monFlux << "Trifurcation : " << connexite[3].size() << endl;
		monFlux << "Quadrifurcation : " << connexite[4].size() << endl;
		monFlux << "Pentafurcation : " << connexite[5].size() << endl;
	}
	else
	{
		cout << "ERREUR : Impossible d'ouvrir le fichier." << endl;
	}
}

void G_classif::Morphometry_node(double voxel)
{
	double somme_diametre = 0 ;
	double somme_diametre_carre = 0 ;
	double somme_longueur = 0 ;
	double longueur_interm = voxel;
	double aire = 0;
	const double mon_pi = 3.14159265358979323846;
	double size_voxel;
	size_voxel = voxel ; // en micro metre
	
	for(int i = 0; i < getNumberNode(); i++)
	{
		for(int j = 0 ; j < l_node[i].getSizeNode(); j++)
		{
			somme_diametre += 2*l_node[i].getValuePixel(j).getRadius() * size_voxel;
			somme_diametre_carre += pow(2*l_node[i].getValuePixel(j).getRadius()*size_voxel, 2); // puissance
		}
		
		for(int j = 0 ; j < l_node[i].getSizeNode()-1; j++)
		{
			longueur_interm = sqrt(pow(l_node[i].getValuePixel(j).getCoordZ()-l_node[i].getValuePixel(j+1).getCoordZ() ,2) + pow(l_node[i].getValuePixel(j).getCoordY()-l_node[i].getValuePixel(j+1).getCoordY() ,2) + pow(l_node[i].getValuePixel(j).getCoordX()-l_node[i].getValuePixel(j+1).getCoordX() ,2));
			
			somme_longueur += longueur_interm * size_voxel;
		}
		//aire de la section, cercle
		aire = mon_pi*pow((somme_diametre/l_node[i].getSizeNode()),2)*pow(size_voxel,2)/4;
		
		diametre_moy_node.push_back(somme_diametre/l_node[i].getSizeNode());
		
		VAR_diametre_moy_node.push_back((somme_diametre_carre/l_node[i].getSizeNode())-pow((somme_diametre/l_node[i].getSizeNode()),2));
		
		lenght_moy_node.push_back(somme_longueur);
		area_node.push_back(aire);
		
		somme_diametre = 0;
		somme_diametre_carre = 0 ;
		somme_longueur = voxel ;
		aire = 0;
	}
}

/*******************************************************************************************/
void G_classif::Strahler(char* name)
{
	vector<int> branche;
	vector<int> filles;
	vector<int> filles2;
	int order_max = 0;
	getListeJct();
	
	for (int i = getHauteurMaxA(); i >= 1; i--) // Pour chaque generation
	{	
		for(unsigned int j = 0; j < l_node_jct.size(); j++)
		{
			if (getHauteurJunction(j) == i) // Pour chaque jonction de la generation
			{
				for(int k = 0; k < getNumberEdge(); k++)
				{
					if (getCoupleNode(k).getValueNode1() == l_node_jct[j]) // Recherche des branches associees a jonction
					{
						branche.push_back(getCoupleNode(k).getValueNode2());
					}
				}
			
				for (unsigned int g=1; g < branche.size()-1; g++)
				{
					filles.push_back(l_node[branche[g]].getOrderNode()); // Enregistrement ordre fille dans vector
				}

				for (unsigned int h=2; h < branche.size(); h++)
				{
					filles2.push_back(l_node[branche[h]].getOrderNode()); // Enregistrement ordre fille dans vector2
				}
				if (std::equal(filles.begin(), filles.end(), filles2.begin()))
				{
					l_node[branche[0]].setOrderNode(l_node[branche[1]].getOrderNode()+1);
				}
				else
				{
					for (unsigned int l=1 ; l < branche.size();l++)
					{
						if (l_node[branche[l]].getOrderNode()>=order_max)
						{
							order_max = l_node[branche[l]].getOrderNode();
							l_node[branche[0]].setOrderNode(order_max);
						}
					}
				}
				branche.clear();
				filles.clear();
				filles2.clear();
				order_max = 0;
			}
		}
	}
	
	// Ordre max dans classif
	getListeDoubleSeg();
	Maxorder = l_node[0].getOrderNode();
	
	//Sauvegarde
	ofstream monFlux(name,ios::app);
	
	if (monFlux)
	{
		ordre.resize(Maxorder+1,0);
		for(unsigned int m = 0; m < l_node_dbseg.size(); m++)
		{
			ordre[l_node[getValueDoubleSeg(m)].getOrderNode()] += 1;
		}
		monFlux << "-- ORDRE STRAHLER : Segment par ordre --" << endl;
		for (unsigned int p = 0; p < ordre.size(); p++)
		{
			monFlux << "Ordre "<< p << " : " << ordre[p] << endl;
		}
	}
	else
	{
		cout << "ERREUR : Impossible d'ouvrir le fichier." << endl;
	}
}

void G_classif::Diameter_order(std::string save, char* name)
{
	getListeDoubleSeg();
	
	vector<double> diametre_moy_nodeCarre;
	vector<double> diametreCarre_moy_node;
	diametre_moy_nodeCarre.resize(Maxorder+1,0);
	diametreCarre_moy_node.resize(Maxorder+1,0);
	
	diametre_moy_order.clear();
	SD_diametre_moy_order.clear();
	diametre_moy_order.resize(Maxorder+1,0);
	SD_diametre_moy_order.resize(Maxorder+1,0);

	std::vector<int> compte(Maxorder+1);
	
	if (save == "NotSave")
	{
		for(unsigned int m = 0; m < l_node_dbseg.size(); m++)
		{
			diametre_moy_order[l_node[getValueDoubleSeg(m)].getOrderNode()] += diametre_moy_node[getValueDoubleSeg(m)];
		
			diametreCarre_moy_node[l_node[getValueDoubleSeg(m)].getOrderNode()] += pow(diametre_moy_node[getValueDoubleSeg(m)],2);
		
			compte[l_node[getValueDoubleSeg(m)].getOrderNode()] += 1;
		}
		for (int i = 0; i < Maxorder+1; i++)
		{
			diametre_moy_order[i] = diametre_moy_order[i]/compte[i];
			diametre_moy_nodeCarre[i] = pow(diametre_moy_order[i],2);
			diametreCarre_moy_node[i] = diametreCarre_moy_node[i]/compte[i];

			SD_diametre_moy_order[i] = sqrt(diametreCarre_moy_node[i]-diametre_moy_nodeCarre[i]);
		}
	}
	
	//Sauvegarde
	if (save == "Save")
	{
		for(unsigned int m = 0; m < l_node_dbseg.size(); m++)
		{
			diametre_moy_order[l_node[getValueDoubleSeg(m)].getOrderNode()] += diametre_moy_node[getValueDoubleSeg(m)];
		
			diametreCarre_moy_node[l_node[getValueDoubleSeg(m)].getOrderNode()] += pow(diametre_moy_node[getValueDoubleSeg(m)],2);
		
			compte[l_node[getValueDoubleSeg(m)].getOrderNode()] += 1;
		}
		for (int i = 0; i < Maxorder+1; i++)
		{
			diametre_moy_order[i] = diametre_moy_order[i]/compte[i];
			diametre_moy_nodeCarre[i] = pow(diametre_moy_order[i],2);
			diametreCarre_moy_node[i] = diametreCarre_moy_node[i]/compte[i];

			SD_diametre_moy_order[i] = sqrt(diametreCarre_moy_node[i]-diametre_moy_nodeCarre[i]);
		}

		ofstream monFlux(name,ios::app);
	
		if (monFlux)
		{
			monFlux << "-- DIAMETRE MOYEN en micromètre & SD par ordre --" << endl;
			for (int p = 0; p < Maxorder+1; p++)
			{
				monFlux << "Ordre "<< p << " : " << diametre_moy_order[p] << " +/- " << SD_diametre_moy_order[p] << endl;
			}
		}
		else
		{
			cout << "ERREUR : Impossible d'ouvrir le fichier." << endl;
		}
	}
}


void G_classif::Strahler_Diameter(char* name, double voxel)
{
	Strahler(name);
	Morphometry_node(voxel);
	vector<int> branche;
	vector<int> filles;
	vector<int> filles2;
	int order_max = 0;
	vector<double> comparaison_ordre(Maxorder+1);
	vector<double> last_ordre;
	int nochange = 1;
	int bcl = 0;
	
	while (nochange == 1)
	{
		nochange = 0;
		Diameter_order("NotSave", name);
		
		last_ordre.clear();
		last_ordre.resize(Maxorder+1,0);
		
		for(int t = 0; t < Maxorder+1; t++)
		{
			last_ordre[t] = ordre[t];
		}
		
		for (int b = 0; b < Maxorder; b++)
		{
			comparaison_ordre[b] = (diametre_moy_order[b]+SD_diametre_moy_order[b]+diametre_moy_order[b+1]-SD_diametre_moy_order[b+1])/2;
		}
	
		for (int i = getHauteurMaxA(); i >= 1; i--) // Pour chaque generation
		{	
			for(unsigned int j = 0; j < l_node_jct.size(); j++)
			{
				if (getHauteurJunction(j) == i) // Pour chaque jonction de la generation
				{
					for(int k = 0; k < getNumberEdge(); k++)
					{
						if (getCoupleNode(k).getValueNode1() == l_node_jct[j]) // Recherche des branches associees a jonction
						{
							branche.push_back(getCoupleNode(k).getValueNode2());
						}
					}
			
					for (unsigned int g=1; g < branche.size()-1; g++)
					{
						filles.push_back(l_node[branche[g]].getOrderNode()); // Enregistrement ordre fille dans vector
					}

					for (unsigned int h=2; h < branche.size(); h++)
					{
						filles2.push_back(l_node[branche[h]].getOrderNode()); // Enregistrement ordre fille dans vector2
					}
					if (std::equal(filles.begin(), filles.end(), filles2.begin()))
					{
						if (diametre_moy_node[branche[0]] > comparaison_ordre[l_node[branche[1]].getOrderNode()])
						{
							l_node[branche[0]].setOrderNode(l_node[branche[1]].getOrderNode()+1);
						}
						else
						{
							l_node[branche[0]].setOrderNode(l_node[branche[1]].getOrderNode());
						}
					}
					else
					{
						for (unsigned int l=1 ; l < branche.size();l++)
						{
							if (l_node[branche[l]].getOrderNode()>=order_max)
							{
								order_max = l_node[branche[l]].getOrderNode();
								l_node[branche[0]].setOrderNode(order_max);
							}
						}
					}
					branche.clear();
					filles.clear();
					filles2.clear();
					order_max = 0;
				}
			}
		}
		
		// Ordre max dans classif
		getListeDoubleSeg();
		Maxorder = l_node[0].getOrderNode();
		
		// Creation nouveaux ordres
		ordre.clear();
		ordre.resize(Maxorder+1,0);
		for(unsigned int m = 0; m < l_node_dbseg.size(); m++)
		{
			ordre[l_node[getValueDoubleSeg(m)].getOrderNode()] += 1;
		}
		
		// Comparaison avec l'ancien
		for(int t = 0; t < Maxorder+1; t++)
		{
			if (last_ordre[t] != ordre[t])
			{
				nochange = 1;
			}
		}
		bcl+=1;
	}
	
	//Sauvegarde
	ofstream monFlux(name,ios::app);
	
	if (monFlux)
	{
		vector<int> ordred(Maxorder+1,0);			
		for(unsigned int m = 0; m < l_node_dbseg.size(); m++)
		{
			ordred[l_node[getValueDoubleSeg(m)].getOrderNode()] += 1;
		}
		monFlux << "-- ORDRE STRAHLER DIAMETRE en micrometre : Segment par ordre -- Iteration : " << bcl << endl;
		for (unsigned int p = 0; p < ordred.size(); p++)
		{
			monFlux << "Ordre "<< p << " : " << ordred[p] << endl;
		}
	}
	else
	{
		cout << "ERREUR : Impossible d'ouvrir le fichier." << endl;
	}
	
}

/*******************************************************************************************/

void G_classif::Connectivity_matrix(char* name)
{
	vector<int> branche;
	vector<int> compte(Maxorder+1,0);
	connectivity_matrix.resize(Maxorder+1);
	getListeJct();
	
	for (int i = 0; i < Maxorder+1; i++)
	{
		connectivity_matrix[i].resize(Maxorder+1,0);
	}
	
	for(unsigned int j = 0; j < l_node_jct.size(); j++)
	{
		for(int k = 0; k < getNumberEdge(); k++)
		{
			if (getCoupleNode(k).getValueNode1() == l_node_jct[j])
			{
				branche.push_back(getCoupleNode(k).getValueNode2());
			}
		}
		for (unsigned int p=1; p < branche.size(); p++)
		{
			connectivity_matrix[l_node[branche[0]].getOrderNode()][l_node[branche[p]].getOrderNode()]+=1;
			compte[l_node[branche[0]].getOrderNode()]+=1;
		}
		branche.clear();
	}
	
	for (int m = 0; m < Maxorder+1; m++)
	{
		for (int n = 0; n < Maxorder+1; n++)
		{
			connectivity_matrix[m][n] = connectivity_matrix[m][n]/compte[m];
		}
	}
	
	//Sauvegarde
	ofstream monFlux(name,ios::app);
	
	if (monFlux)
	{
		monFlux << "-- MATRICE DE CONNEXITE par ordre --" << endl;
		for (int p = 0; p < Maxorder+1; p++)
		{
			monFlux << "Parents ordre " << p << " : " << endl;
			for (int o = 0; o < Maxorder+1; o++)
			{
				monFlux << "Filles ordre "<< o << " : " << connectivity_matrix[p][o] << endl;
			}
		}
	}
	else
	{
		cout << "ERREUR : Impossible d'ouvrir le fichier." << endl;
	}
	
}



/*******************************************************************************************/

void G_classif::Lenght_order(char* name)
{
	std::vector<int> compte(Maxorder+1);
	lenght_moy_order.resize(Maxorder+1,0);
	getListeDoubleSeg();
	
	for(unsigned int m = 1; m < l_node_dbseg.size(); m++)
	{
		lenght_moy_order[l_node[getValueDoubleSeg(m)].getOrderNode()] += lenght_moy_node[getValueDoubleSeg(m)];
				
		compte[l_node[getValueDoubleSeg(m)].getOrderNode()] += 1;
	}
	for (int i = 0; i < Maxorder+1; i++)
	{
		lenght_moy_order[i] = lenght_moy_order[i]/compte[i];
	}
	
	
	//Sauvegarde
	ofstream monFlux(name,ios::app);
	
	if (monFlux)
	{
		monFlux << "-- LONGUEUR MOYENNE en micromètre par ordre --" << endl;
		for (int p = 0; p < Maxorder+1; p++)
		{
			monFlux << "Ordre "<< p << " : " << lenght_moy_order[p] << endl;
		}
	}
	else
	{
		cout << "ERREUR : Impossible d'ouvrir le fichier." << endl;
	}

}

/*******************************************************************************************/
void G_classif::Abo_Cassot(char* name, double voxel)
{
	Barycentre();
	
	angle_m1.clear();
	angle_m2.clear();
	angle_12.clear();
	angle_beta.clear();
	angle_omega.clear();
	
	int fille2;
	int size = 0 ;
        int numero = -1;
	double rayonjct;
	double longueur_interm = voxel;
	double a_cos;
	vector<int> branche;
	
	double CoordZ;
	double CoordY;
	double CoordX;
	
	vector<double> vect_m;
	vector<double> vect_f1;
	vector<double> vect_f2;
	
	double norme_m;
        double norme_f1 = 0.0;
        double norme_f2 = 0.0;
	
	vector<double> produit_vect;
	double norme_prodvect;
	double pmix;				
	double comortho;					
	double spb;			
	vector<double> aplanbif;
	double renorm;
					
	vector<double> prod_scal;
	vector<double> angle_interm;
	double size_voxel;
	size_voxel = voxel ;
	const double mon_pi = 3.14159265358979323846;
	
	std::vector<int> compte(Maxorder+1);
	angle_m1.resize(Maxorder+1,0);
	angle_m2.resize(Maxorder+1,0);
	angle_12.resize(Maxorder+1,0);
	angle_beta.resize(Maxorder+1,0);
	angle_omega.resize(Maxorder+1,0);
	
	// Pour chaque bifurcation				
	for (unsigned int i = 0; i < connexite[2].size(); i++)
	{
		if (connexite[2][i] != 1)
		{
			size = 0;
			branche.push_back(connexite[2][i]); // id jct
			for(int k = 0; k < getNumberEdge(); k++)
			{
				if (getCoupleNode(k).getValueNode1() == connexite[2][i])
				{
					branche.push_back(getCoupleNode(k).getValueNode2()); // id mere, fille1, fille2
				}
			}
		
			if (diametre_moy_node[branche[3]] > diametre_moy_node[branche[2]]) 
			{
				fille2 = branche[2];
				branche[2]=branche[3];
				branche[3]=fille2;
			}
		
			rayonjct = l_node[branche[0]].getValuePixel(0).getRadius()*size_voxel;
			
			// Numero du dernier point a l'interieur du cercle rayon jonction
			// Cas branche mère
			for(int m = l_node[branche[1]].getSizeNode()-1; m >= 0; m--)
			{
				longueur_interm = sqrt(pow(l_node[branche[1]].getValuePixel(m).getCoordZ()-l_node[branche[0]].getValuePixel(0).getCoordZ() ,2) + pow(l_node[branche[1]].getValuePixel(m).getCoordY()-l_node[branche[0]].getValuePixel(0).getCoordY() ,2) + pow(l_node[branche[1]].getValuePixel(m).getCoordX()-l_node[branche[0]].getValuePixel(0).getCoordX() ,2))*size_voxel;
				if (longueur_interm <= rayonjct)
				{
					size +=1;				
				}
			}
		
			//cout << "size " << size << endl;
		
			if (size == l_node[branche[1]].getSizeNode())
			{
				numero = 0 ;
			}
			if (size < l_node[branche[1]].getSizeNode())
			{
				if (size==0 && l_node[branche[1]].getSizeNode() > 1)
				{
					numero = l_node[branche[1]].getSizeNode()-2 ;
				}
				else 
				{
					numero = l_node[branche[1]].getSizeNode()-1-size ;
				}
			}

			CoordZ = l_node[branche[1]].getValuePixel(numero).getCoordZ();
			CoordY = l_node[branche[1]].getValuePixel(numero).getCoordY();
			CoordX = l_node[branche[1]].getValuePixel(numero).getCoordX();
		
			// Calcul Vecteur + norme
			vect_m.push_back(CoordZ-l_node[branche[0]].getValuePixel(0).getCoordZ());
			vect_m.push_back(CoordY-l_node[branche[0]].getValuePixel(0).getCoordY());
			vect_m.push_back(CoordX-l_node[branche[0]].getValuePixel(0).getCoordX());
			norme_m = sqrt(pow(vect_m[0],2)+pow(vect_m[1],2)+pow(vect_m[2],2));
						
			// Cas branches filles
			for (int f = 2; f<4 ; f++)
			{
				size = 0;
			
				//cout << "**************************" << endl;
				//cout << branche[f] << endl;
				//cout << "rayon "<< rayonjct << endl;
				for(int m = 0; m < l_node[branche[f]].getSizeNode(); m++)
				{
					longueur_interm = sqrt(pow(l_node[branche[f]].getValuePixel(m).getCoordZ()-l_node[branche[0]].getValuePixel(0).getCoordZ() ,2) + pow(l_node[branche[f]].getValuePixel(m).getCoordY()-l_node[branche[0]].getValuePixel(0).getCoordY() ,2) + pow(l_node[branche[f]].getValuePixel(m).getCoordX()-l_node[branche[0]].getValuePixel(0).getCoordX() ,2))*size_voxel;
					if (longueur_interm <= rayonjct)
					{
						size +=1;				
					}
				}
				//cout << "size " << size << endl;
				if (size == l_node[branche[f]].getSizeNode())
				{
					numero = l_node[branche[f]].getSizeNode()-1 ;
				}
				if (size < l_node[branche[f]].getSizeNode())
				{
					if (size==0 && l_node[branche[f]].getSizeNode() >1)
					{
						numero = 1 ;
					}
					else 
					{
						numero = size ;
					}
				}

				CoordZ = l_node[branche[f]].getValuePixel(numero).getCoordZ();
				CoordY = l_node[branche[f]].getValuePixel(numero).getCoordY();
				CoordX = l_node[branche[f]].getValuePixel(numero).getCoordX();
				//cout << CoordZ << " " << CoordY << " " << CoordX << endl;
				// Calcul Vecteur + norme
				if (f==2)
				{
					vect_f1.push_back(CoordZ-l_node[branche[0]].getValuePixel(0).getCoordZ());
					vect_f1.push_back(CoordY-l_node[branche[0]].getValuePixel(0).getCoordY());
					vect_f1.push_back(CoordX-l_node[branche[0]].getValuePixel(0).getCoordX());
					norme_f1 = sqrt(pow(vect_f1[0],2)+pow(vect_f1[1],2)+pow(vect_f1[2],2));
				}
				if (f==3)
				{
					vect_f2.push_back(CoordZ-l_node[branche[0]].getValuePixel(0).getCoordZ());
					vect_f2.push_back(CoordY-l_node[branche[0]].getValuePixel(0).getCoordY());
					vect_f2.push_back(CoordX-l_node[branche[0]].getValuePixel(0).getCoordX());
					norme_f2 = sqrt(pow(vect_f2[0],2)+pow(vect_f2[1],2)+pow(vect_f2[2],2));
				}
			}
		
			// Normalisation des Vecteurs
			for(unsigned int l=0; l < vect_m.size(); l++)
			{
				vect_m[l] = vect_m[l]/norme_m;
				vect_f1[l] = vect_f1[l]/norme_f1;
				vect_f2[l] = vect_f2[l]/norme_f2;
			}
		
			// Produit scalaire
			// CAS M1
			prod_scal.push_back(vect_m[0]*vect_f1[0] + vect_m[1]*vect_f1[1] + vect_m[2]*vect_f1[2]);
		
			// CAS M2
			prod_scal.push_back(vect_m[0]*vect_f2[0] + vect_m[1]*vect_f2[1] + vect_m[2]*vect_f2[2]);
		
			// CAS 12
			prod_scal.push_back(vect_f1[0]*vect_f2[0] + vect_f1[1]*vect_f2[1] + vect_f1[2]*vect_f2[2]);

			// Composante du produit vectoriel
			produit_vect.push_back(vect_f1[1]*vect_f2[2] - vect_f1[2]*vect_f2[1]);
			produit_vect.push_back(vect_f1[2]*vect_f2[0] - vect_f1[0]*vect_f2[2]);
			produit_vect.push_back(vect_f1[0]*vect_f2[1] - vect_f1[1]*vect_f2[0]);
			norme_prodvect = sqrt(pow(produit_vect[0],2)+pow(produit_vect[1],2)+pow(produit_vect[2],2));
		
			// Produit mixte
			pmix = vect_m[0]*produit_vect[0] + vect_m[1]*produit_vect[1] + vect_m[2]*produit_vect[2];
				
			if (norme_prodvect > 0)
			{
				// Normalisation prod_vect
				produit_vect[0] = produit_vect[0]/norme_prodvect;
				produit_vect[1] = produit_vect[1]/norme_prodvect;
				produit_vect[2] = produit_vect[2]/norme_prodvect;
			}	
			
			// Normalisation
			comortho = pmix/norme_prodvect;
			spb = abs(comortho);
			
			if (spb < 1)
			{
				aplanbif.push_back((180/mon_pi)*atan(spb/sqrt(1-pow(spb,2))));			
			}
			else
			{
				aplanbif.push_back(90);
			}
			
			// Projection sur plan des branches filles de la branche mere
			for(unsigned int l=0; l < vect_m.size(); l++)
			{
				vect_m[l] = vect_m[l] - comortho*produit_vect[l];
			}
			renorm = sqrt(pow(vect_m[0],2)+pow(vect_m[1],2)+pow(vect_m[2],2));
			if (renorm>0)
			{
				vect_m[0] = vect_m[0]/renorm;
				vect_m[1] = vect_m[1]/renorm;
				vect_m[2] = vect_m[2]/renorm;
			}
			
			// Calcul angles
			// CAS M1
			a_cos = vect_m[0]*vect_f1[0] + vect_m[1]*vect_f1[1] + vect_m[2]*vect_f1[2];
			if (a_cos < -1)
			{
				a_cos = -1;
			}
			if (a_cos > 1)
			{
				a_cos = 1;
			}
			angle_interm.push_back(acos(a_cos));
			//cout << (180/mon_pi)*acos(a_cos) << endl;
			// CAS M2
			a_cos = vect_m[0]*vect_f2[0] + vect_m[1]*vect_f2[1] + vect_m[2]*vect_f2[2];
			if (a_cos < -1)
			{
				a_cos = -1;
			}
			if (a_cos > 1)
			{
				a_cos = 1;
			}
			angle_interm.push_back(acos(a_cos));
			//	cout << (180/mon_pi)*acos(a_cos) << endl;
			// CAS 12
			a_cos = vect_f1[0]*vect_f2[0] + vect_f1[1]*vect_f2[1] + vect_f1[2]*vect_f2[2];
			if (a_cos < -1)
			{
				a_cos = -1;
			}
			if (a_cos > 1)
			{
				a_cos = 1;
			}
			angle_interm.push_back(acos(a_cos));
			//cout << (180/mon_pi)*acos(a_cos) << endl;
			// CAS Assymetrie Beta
			angle_interm.push_back(angle_interm[0]-angle_interm[1]);
		
			// CAS Planarite Omega
			angle_interm.push_back(angle_interm[0]+angle_interm[1]+angle_interm[2]-mon_pi);

			angle_node.push_back(angle_interm);
			
			angle_m1[l_node[branche[1]].getOrderNode()] += angle_interm[0];
			angle_m2[l_node[branche[1]].getOrderNode()] += angle_interm[1];
			angle_12[l_node[branche[1]].getOrderNode()] += angle_interm[2];
			angle_beta[l_node[branche[1]].getOrderNode()] += angle_interm[3];
			angle_omega[l_node[branche[1]].getOrderNode()] += angle_interm[4];
		
			compte[l_node[branche[1]].getOrderNode()] += 1;
		}
		// Nettoyage
		branche.clear();
		vect_m.clear();
		vect_f1.clear();
		vect_f2.clear();
		prod_scal.clear();
		angle_interm.clear();
		produit_vect.clear();
		aplanbif.clear();
	}
	
	for (int j = 0; j < Maxorder+1; j++)
	{
		angle_m1[j] = angle_m1[j]/compte[j];
		angle_m2[j] = angle_m2[j]/compte[j];
		angle_12[j] = angle_12[j]/compte[j];
		angle_beta[j] = angle_beta[j]/compte[j];
		angle_omega[j] = angle_omega[j]/compte[j];	
	}

	compte.clear();	
	
	//Sauvegarde
	ofstream monFlux(name,ios::app);

	if (monFlux)
	{
		monFlux << "-- ANGLES CASSOT en degré par ordre --" << endl;
		monFlux << "Angles M1" << endl;			
		for (int p = 0; p < Maxorder+1; p++)
		{
			monFlux << "Ordre "<< p << " : " << angle_m1[p]*180/mon_pi << endl;
		}
		monFlux << "Angles M2" << endl;		
		for (int p = 0; p < Maxorder+1; p++)
		{
			monFlux << "Ordre "<< p << " : " << angle_m2[p]*180/mon_pi << endl;
		}
		monFlux << "Angles 12" << endl;		
		for (int p = 0; p < Maxorder+1; p++)
		{
			monFlux << "Ordre "<< p << " : " << angle_12[p]*180/mon_pi << endl;
		}
		monFlux << "Angles Beta Assymetrie" << endl;		
		for (int p = 0; p < Maxorder+1; p++)
		{
			monFlux << "Ordre "<< p << " : " << angle_beta[p]*180/mon_pi << endl;
		}
		monFlux << "Angles Omega Planarite" << endl;		
		for (int p = 0; p < Maxorder+1; p++)
		{
			monFlux << "Ordre "<< p << " : " << angle_omega[p]*180/mon_pi << endl;
		}
	}
	else
	{
		cout << "ERREUR : Impossible d'ouvrir le fichier." << endl;
	}

}





/*******************************************************************************************/
void G_classif::Angle_Bifurcation(char* name, double voxel)
{
	Barycentre();
	
	angle_m1.clear();
	angle_m2.clear();
	angle_12.clear();
	angle_beta.clear();
	angle_omega.clear();
	
	int fille2 = 0;
	int size = 0 ;
	int numero = 0;
	double rayonjct = 0;
	double longueur_interm = voxel;
	double a_cos = 0;
	vector<int> branche;
	
	double CoordZ=0;
	double CoordY=0;
	double CoordX=0;
	
	vector<double> vect_m;
	vector<double> vect_f1;
	vector<double> vect_f2;
	
	double norme_m=0;
	double norme_f1=0;
	double norme_f2=0;
				
	vector<double> prod_scal;
	vector<double> angle_interm;
	double size_voxel;
	size_voxel = voxel ;
	const double mon_pi = 3.14159265358979323846;
	
	std::vector<int> compte(Maxorder+1);
	angle_m1.resize(Maxorder+1,0);
	angle_m2.resize(Maxorder+1,0);
	angle_12.resize(Maxorder+1,0);
	angle_beta.resize(Maxorder+1,0);
	angle_omega.resize(Maxorder+1,0);
	
	// Pour chaque bifurcation				
	for (unsigned int i = 0; i < connexite[2].size(); i++)
	{
		if (connexite[2][i] != 1)
		{
			size = 0;
			branche.push_back(connexite[2][i]); // id jct
			for(int k = 0; k < getNumberEdge(); k++)
			{
				if (getCoupleNode(k).getValueNode1() == connexite[2][i])
				{
					branche.push_back(getCoupleNode(k).getValueNode2()); // id mere, fille1, fille2
				}
			}
		
			if (diametre_moy_node[branche[3]] > diametre_moy_node[branche[2]]) 
			{
				fille2 = branche[2];
				branche[2]=branche[3];
				branche[3]=fille2;
			}
		
			rayonjct = l_node[branche[0]].getValuePixel(0).getRadius()*size_voxel;
			
			// Numero du dernier point a l'interieur du cercle rayon jonction
			// Cas branche mère
			for(int m = l_node[branche[1]].getSizeNode()-1; m >= 0; m--)
			{
				longueur_interm = sqrt(pow(l_node[branche[1]].getValuePixel(m).getCoordZ()-l_node[branche[0]].getValuePixel(0).getCoordZ() ,2) + pow(l_node[branche[1]].getValuePixel(m).getCoordY()-l_node[branche[0]].getValuePixel(0).getCoordY() ,2) + pow(l_node[branche[1]].getValuePixel(m).getCoordX()-l_node[branche[0]].getValuePixel(0).getCoordX() ,2))*size_voxel;
				if (longueur_interm <= rayonjct)
				{
					size +=1;				
				}
			}
		
			//cout << "size " << size << endl;
		
			if (size == l_node[branche[1]].getSizeNode())
			{
				numero = 0 ;
			}
			if (size < l_node[branche[1]].getSizeNode())
			{
				if (size==0 && l_node[branche[1]].getSizeNode() > 1)
				{
					numero = l_node[branche[1]].getSizeNode()-2 ;
				}
				else 
				{
					numero = l_node[branche[1]].getSizeNode()-1-size ;
				}
			}

			CoordZ = l_node[branche[1]].getValuePixel(numero).getCoordZ();
			CoordY = l_node[branche[1]].getValuePixel(numero).getCoordY();
			CoordX = l_node[branche[1]].getValuePixel(numero).getCoordX();
		
			// Calcul Vecteur + norme
			vect_m.push_back(CoordZ-l_node[branche[0]].getValuePixel(0).getCoordZ());
			vect_m.push_back(CoordY-l_node[branche[0]].getValuePixel(0).getCoordY());
			vect_m.push_back(CoordX-l_node[branche[0]].getValuePixel(0).getCoordX());
			norme_m = sqrt(pow(vect_m[0],2)+pow(vect_m[1],2)+pow(vect_m[2],2));
						
			// Cas branches filles
			for (int f = 2; f<4 ; f++)
			{
				size = 0;
			
				//cout << "**************************" << endl;
				//cout << branche[f] << endl;
				//cout << "rayon "<< rayonjct << endl;
				for(int m = 0; m < l_node[branche[f]].getSizeNode(); m++)
				{
					longueur_interm = sqrt(pow(l_node[branche[f]].getValuePixel(m).getCoordZ()-l_node[branche[0]].getValuePixel(0).getCoordZ() ,2) + pow(l_node[branche[f]].getValuePixel(m).getCoordY()-l_node[branche[0]].getValuePixel(0).getCoordY() ,2) + pow(l_node[branche[f]].getValuePixel(m).getCoordX()-l_node[branche[0]].getValuePixel(0).getCoordX() ,2))*size_voxel;
					if (longueur_interm <= rayonjct)
					{
						size +=1;				
					}
				}
				//cout << "size " << size << endl;
				if (size == l_node[branche[f]].getSizeNode())
				{
					numero = l_node[branche[f]].getSizeNode()-1 ;
				}
				if (size < l_node[branche[f]].getSizeNode())
				{
					if (size==0 && l_node[branche[f]].getSizeNode() >1)
					{
						numero = 1 ;
					}
					else 
					{
						numero = size ;
					}
				}

				CoordZ = l_node[branche[f]].getValuePixel(numero).getCoordZ();
				CoordY = l_node[branche[f]].getValuePixel(numero).getCoordY();
				CoordX = l_node[branche[f]].getValuePixel(numero).getCoordX();
				//cout << CoordZ << " " << CoordY << " " << CoordX << endl;
				// Calcul Vecteur + norme
				if (f==2)
				{
					vect_f1.push_back(CoordZ-l_node[branche[0]].getValuePixel(0).getCoordZ());
					vect_f1.push_back(CoordY-l_node[branche[0]].getValuePixel(0).getCoordY());
					vect_f1.push_back(CoordX-l_node[branche[0]].getValuePixel(0).getCoordX());
					norme_f1 = sqrt(pow(vect_f1[0],2)+pow(vect_f1[1],2)+pow(vect_f1[2],2));
				}
				if (f==3)
				{
					vect_f2.push_back(CoordZ-l_node[branche[0]].getValuePixel(0).getCoordZ());
					vect_f2.push_back(CoordY-l_node[branche[0]].getValuePixel(0).getCoordY());
					vect_f2.push_back(CoordX-l_node[branche[0]].getValuePixel(0).getCoordX());
					norme_f2 = sqrt(pow(vect_f2[0],2)+pow(vect_f2[1],2)+pow(vect_f2[2],2));
				}
			}
			
			// Produit scalaire
			// CAS M1
			prod_scal.push_back(vect_m[0]*vect_f1[0] + vect_m[1]*vect_f1[1] + vect_m[2]*vect_f1[2]);
		
			// CAS M2
			prod_scal.push_back(vect_m[0]*vect_f2[0] + vect_m[1]*vect_f2[1] + vect_m[2]*vect_f2[2]);
		
			// CAS 12
			prod_scal.push_back(vect_f1[0]*vect_f2[0] + vect_f1[1]*vect_f2[1] + vect_f1[2]*vect_f2[2]);
		
		
			// Angles
			// CAS M1
			a_cos = prod_scal[0]/(norme_m*norme_f1);
			if (a_cos < -1)
			{
				a_cos = -1;
			}
			if (a_cos > 1)
			{
				a_cos = 1;
			}
			angle_interm.push_back(acos(a_cos));

			// CAS M2
			a_cos = prod_scal[1]/(norme_m*norme_f2);
			if (a_cos < -1)
			{
				a_cos = -1;
			}
			if (a_cos > 1)
			{
				a_cos = 1;
			}
			angle_interm.push_back(acos(a_cos));
		
			// CAS 12
			a_cos = prod_scal[2]/(norme_f1*norme_f2);
			if (a_cos < -1)
			{
				a_cos = -1;
			}
			if (a_cos > 1)
			{
				a_cos = 1;
			}
			angle_interm.push_back(acos(a_cos));
				
			// CAS Assymetrie Beta
			angle_interm.push_back(angle_interm[0]-angle_interm[1]);
		
			// CAS Planarite Omega
			angle_interm.push_back(angle_interm[0]+angle_interm[1]+angle_interm[2]-mon_pi);

			angle_node.push_back(angle_interm);
		
			angle_m1[l_node[branche[1]].getOrderNode()] += angle_interm[0];
			angle_m2[l_node[branche[1]].getOrderNode()] += angle_interm[1];
			angle_12[l_node[branche[1]].getOrderNode()] += angle_interm[2];
			angle_beta[l_node[branche[1]].getOrderNode()] += angle_interm[3];
			angle_omega[l_node[branche[1]].getOrderNode()] += angle_interm[4];
		
			compte[l_node[branche[1]].getOrderNode()] += 1;


		}
		// Nettoyage
		branche.clear();
		vect_m.clear();
		vect_f1.clear();
		vect_f2.clear();
		prod_scal.clear();
		angle_interm.clear();
	}
	
	for (int j = 0; j < Maxorder+1; j++)
	{
		angle_m1[j] = angle_m1[j]/compte[j];
		angle_m2[j] = angle_m2[j]/compte[j];
		angle_12[j] = angle_12[j]/compte[j];
		angle_beta[j] = angle_beta[j]/compte[j];
		angle_omega[j] = angle_omega[j]/compte[j];	
	}

	compte.clear();

	//Sauvegarde
	ofstream monFlux(name,ios::app);
	
	if (monFlux)
	{
		monFlux << "-- ANGLES en degré par ordre --" << endl;
		monFlux << "Angles M1" << endl;			
		for (int p = 0; p < Maxorder+1; p++)
		{
			monFlux << "Ordre "<< p << " : " << angle_m1[p]*180/mon_pi << endl;
		}
		monFlux << "Angles M2" << endl;		
		for (int p = 0; p < Maxorder+1; p++)
		{
			monFlux << "Ordre "<< p << " : " << angle_m2[p]*180/mon_pi << endl;
		}
		monFlux << "Angles 12" << endl;		
		for (int p = 0; p < Maxorder+1; p++)
		{
			monFlux << "Ordre "<< p << " : " << angle_12[p]*180/mon_pi << endl;
		}
		monFlux << "Angles Beta Assymetrie" << endl;		
		for (int p = 0; p < Maxorder+1; p++)
		{
			monFlux << "Ordre "<< p << " : " << angle_beta[p]*180/mon_pi << endl;
		}
		monFlux << "Angles Omega Planarite" << endl;		
		for (int p = 0; p < Maxorder+1; p++)
		{
			monFlux << "Ordre "<< p << " : " << angle_omega[p]*180/mon_pi << endl;
		}
	}
	else
	{
		cout << "ERREUR : Impossible d'ouvrir le fichier." << endl;
	}

}

/*******************************************************************************************/

void G_classif::Area_ratio_order(char* name)
{
	vector<int> branche;
	double fille2;
	int fille2bis;
	vector<double> mes_aires;
	std::vector<int> compte(Maxorder+1);
	total_ratio.resize(Maxorder+1,0);
	major_ratio.resize(Maxorder+1,0);
	minor_ratio.resize(Maxorder+1,0);
	ass_ratio.resize(Maxorder+1,0);
				
	for (unsigned int i = 0; i < connexite[2].size(); i++)
	{
		if (connexite[2][i] != 1)
		{
			for(int k = 0; k < getNumberEdge(); k++)
			{
				if (getCoupleNode(k).getValueNode1() == connexite[2][i])
				{
					branche.push_back(getCoupleNode(k).getValueNode2());
					mes_aires.push_back(area_node[getCoupleNode(k).getValueNode2()]);
				}
			}
			if (mes_aires[2]>mes_aires[1])
			{
				fille2=mes_aires[1];
				mes_aires[1]=mes_aires[2];
				mes_aires[2]=fille2;
				fille2bis=branche[1];
				branche[1]=branche[2];
				branche[2]=fille2bis;
			}	
			total_ratio[l_node[branche[0]].getOrderNode()] += (mes_aires[1] + mes_aires[2])/mes_aires[0];
			major_ratio[l_node[branche[0]].getOrderNode()] += mes_aires[1]/mes_aires[0];
			minor_ratio[l_node[branche[0]].getOrderNode()] += mes_aires[2]/mes_aires[0];
			ass_ratio[l_node[branche[0]].getOrderNode()] += mes_aires[2]/mes_aires[1];
			compte[l_node[branche[0]].getOrderNode()] += 1;
			branche.clear();
			mes_aires.clear();
		}
	}
	for (int j = 0; j < Maxorder+1; j++)
	{
		total_ratio[j] = total_ratio[j]/compte[j];
		major_ratio[j] = major_ratio[j]/compte[j];
		minor_ratio[j] = minor_ratio[j]/compte[j];
		ass_ratio[j] = ass_ratio[j]/compte[j];
	}
	
	//Sauvegarde
	ofstream monFlux(name,ios::app);
	
	if (monFlux)
	{
		monFlux << "-- RATIOS DES AIRES par ordre --" << endl;
		monFlux << "Total Ratio" << endl;			
		for (int p = 0; p < Maxorder+1; p++)
		{
			monFlux << "Ordre "<< p << " : " << total_ratio[p] << endl;
		}
		monFlux << "Major Ratio" << endl;		
		for (int p = 0; p < Maxorder+1; p++)
		{
			monFlux << "Ordre "<< p << " : " << major_ratio[p] << endl;
		}
		monFlux << "Minor Ratio" << endl;		
		for (int p = 0; p < Maxorder+1; p++)
		{
			monFlux << "Ordre "<< p << " : " << minor_ratio[p] << endl;
		}
		monFlux << "Assymetry Ratio" << endl;		
		for (int p = 0; p < Maxorder+1; p++)
		{
			monFlux << "Ordre "<< p << " : " << ass_ratio[p] << endl;
		}
	}
	else
	{
		cout << "ERREUR : Impossible d'ouvrir le fichier." << endl;
	}
	
}

/*******************************************************************************************/

void G_classif::visualization()
{
	// renderer renders into the render window.
  vtkSmartPointer<vtkRenderer> ren1 = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
  renWin->AddRenderer(ren1);
  vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  iren->SetRenderWindow(renWin);
	ren1->SetBackground(1, 1, 1);
	
	for(int i = 0; i < getNumberNode(); i++) // Boucle pour chaque point de chaque noeud
	{
		if (l_node[i].getTypeNode() != "junction") //Si segment couleur selon ordre
		{
			for(int j = 0 ; j < l_node[i].getSizeNode(); j++)
			{
				vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
  			sphereSource->SetThetaResolution(50);
  			sphereSource->SetPhiResolution(50);
  			sphereSource->SetCenter(l_node[i].getValuePixel(j).getCoordZ(), l_node[i].getValuePixel(j).getCoordY(), l_node[i].getValuePixel(j).getCoordX());
  			sphereSource->SetRadius(l_node[i].getValuePixel(j).getRadius());
  			sphereSource->Update();
 
  			vtkSmartPointer<vtkPolyDataMapper> sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  			sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
  			vtkSmartPointer<vtkActor> sphereActor = vtkSmartPointer<vtkActor>::New();
 				
 				switch(l_node[i].getOrderNode())
 				{
					case 0:
						sphereActor->GetProperty()->SetColor(0,0,0.5);
						break;
					case 1:
						sphereActor->GetProperty()->SetColor(0,0,1);
						break;
					case 2:
						sphereActor->GetProperty()->SetColor(0,0.5,0);
						break;
					case 3:
						sphereActor->GetProperty()->SetColor(1,0.8,0);
						break;
					case 4:
						sphereActor->GetProperty()->SetColor(1,0,0);
						break;
					case 5:
						sphereActor->GetProperty()->SetColor(0.5,0,0);
						break;
				}
				sphereActor->SetMapper(sphereMapper);
				ren1->AddActor(sphereActor);
			}
		}
		
		if (l_node[i].getTypeNode() == "junction") //Si junction couleur selon ordre segment precedent
		{
			for(int j = 0 ; j < l_node[i].getSizeNode(); j++)
			{
				vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
  			sphereSource->SetThetaResolution(50);
  			sphereSource->SetPhiResolution(50);
  			sphereSource->SetCenter(l_node[i].getValuePixel(j).getCoordZ(), l_node[i].getValuePixel(j).getCoordY(), l_node[i].getValuePixel(j).getCoordX());
  			sphereSource->SetRadius(l_node[i].getValuePixel(j).getRadius());
  			sphereSource->Update();
 
  			vtkSmartPointer<vtkPolyDataMapper> sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  			sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
  			vtkSmartPointer<vtkActor> sphereActor = vtkSmartPointer<vtkActor>::New();
 				
 				switch(l_node[i-1].getOrderNode())
 				{
					case 0:
						sphereActor->GetProperty()->SetColor(0,0,0.5);
						break;
					case 1:
						sphereActor->GetProperty()->SetColor(0,0,1);
						break;
					case 2:
						sphereActor->GetProperty()->SetColor(0,0.5,0);
						break;
					case 3:
						sphereActor->GetProperty()->SetColor(1,0.8,0);
						break;
					case 4:
						sphereActor->GetProperty()->SetColor(1,0,0);
						break;
					case 5:
						sphereActor->GetProperty()->SetColor(0.5,0,0);
						break;
				}
				sphereActor->SetMapper(sphereMapper);
				ren1->AddActor(sphereActor);
			}
		}
	}
	// This starts the event loop and invokes an initial render.
  iren->Initialize();
  iren->Start();
}

void G_classif::visualizationJct()
{
	// renderer renders into the render window.
  vtkSmartPointer<vtkRenderer> ren1 = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
  renWin->AddRenderer(ren1);
  vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  iren->SetRenderWindow(renWin);
	ren1->SetBackground(1, 1, 1);
	
	for(int i = 0; i < getNumberNode(); i++) // Boucle pour chaque point de chaque noeud
	{
		if (l_node[i].getTypeNode() == "junction") //Si junction couleur selon ordre segment precedent
		{
			for(int j = 0 ; j < l_node[i].getSizeNode(); j++)
			{
				vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
  			sphereSource->SetThetaResolution(50);
  			sphereSource->SetPhiResolution(50);
  			sphereSource->SetCenter(l_node[i].getValuePixel(j).getCoordX(), l_node[i].getValuePixel(j).getCoordY(), l_node[i].getValuePixel(j).getCoordZ());
  			sphereSource->SetRadius(1);
  			sphereSource->Update();
 
  			vtkSmartPointer<vtkPolyDataMapper> sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  			sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
  			vtkSmartPointer<vtkActor> sphereActor = vtkSmartPointer<vtkActor>::New();
 				
 				switch(l_node[i-1].getOrderNode())
 				{
					case 0:
						sphereActor->GetProperty()->SetColor(0,0,0.5);
						break;
					case 1:
						sphereActor->GetProperty()->SetColor(0,0,1);
						break;
					case 2:
						sphereActor->GetProperty()->SetColor(0,0.5,0);
						break;
					case 3:
						sphereActor->GetProperty()->SetColor(1,0.8,0);
						break;
					case 4:
						sphereActor->GetProperty()->SetColor(1,0,0);
						break;
					case 5:
						sphereActor->GetProperty()->SetColor(0.5,0,0);
						break;
				}
				sphereActor->SetMapper(sphereMapper);
				ren1->AddActor(sphereActor);
			}
		}
	}
	// This starts the event loop and invokes an initial render.
  iren->Initialize();
  iren->Start();
}

void G_classif::visualizationRadius()
{
	// renderer renders into the render window.
  vtkSmartPointer<vtkRenderer> ren1 = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
  renWin->AddRenderer(ren1);
  vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  iren->SetRenderWindow(renWin);
	ren1->SetBackground(1, 1, 1);
	
	for(int i = 0; i < getNumberNode(); i++) // Boucle pour chaque point de chaque noeud
	{
		if (l_node[i].getTypeNode() != "junction") //Si segment couleur selon ordre
		{
			for(int j = 0 ; j < l_node[i].getSizeNode(); j++)
			{
				vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
  			sphereSource->SetThetaResolution(50);
  			sphereSource->SetPhiResolution(50);
  			sphereSource->SetCenter(l_node[i].getValuePixel(j).getCoordX(), l_node[i].getValuePixel(j).getCoordY(), l_node[i].getValuePixel(j).getCoordZ());
  			sphereSource->SetRadius(1);
  			sphereSource->Update();
 
  			vtkSmartPointer<vtkPolyDataMapper> sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  			sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
  			vtkSmartPointer<vtkActor> sphereActor = vtkSmartPointer<vtkActor>::New();
 				
 				switch(l_node[i].getOrderNode())
 				{
					case 0:
						sphereActor->GetProperty()->SetColor(0,0,0.5);
						break;
					case 1:
						sphereActor->GetProperty()->SetColor(0,0,1);
						break;
					case 2:
						sphereActor->GetProperty()->SetColor(0,0.5,0);
						break;
					case 3:
						sphereActor->GetProperty()->SetColor(1,0.8,0);
						break;
					case 4:
						sphereActor->GetProperty()->SetColor(1,0,0);
						break;
					case 5:
						sphereActor->GetProperty()->SetColor(0.5,0,0);
						break;
				}
				sphereActor->SetMapper(sphereMapper);
				ren1->AddActor(sphereActor);
			}
		}
	}
	
	// This starts the event loop and invokes an initial render.
  iren->Initialize();
  iren->Start();
}
