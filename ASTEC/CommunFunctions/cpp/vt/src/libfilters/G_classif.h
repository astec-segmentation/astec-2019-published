/*****************************************************************************
* G_classif.h
*
* par Manon Linder
*
*
*/


#ifndef _G_classif_h_
#define _G_classif_h_


#include "D_graphe.h"
#include <cmath>

// Visu 3D sphere
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include <vtkXMLPolyDataWriter.h>

class G_classif : public D_graphe
{
	public : // Methodes
	
	G_classif();
	
	// Statistiques simples
	void getConnexiteType(char* name);
	void Morphometry_node(double voxel);
		
	// Fonctions de classification de l'arbre
	void Strahler(char* name);
	void Diameter_order(std::string save, char* name);
	void Strahler_Diameter(char* name, double voxel);
	
	// Matrice de connectixite
	void Connectivity_matrix(char* name);
	
	// Statistiques par ordre
	void Lenght_order(char* name);
	
	// Fonctions calcul angles
	void Abo_Cassot(char* name, double voxel);
	void Angle_Bifurcation(char* name, double voxel);
	
	// Fonctions calcul aires
	void Area_ratio_order(char* name);
	
	// Fonctions de visualisation de l'arbre avec VTK
	void visualization();	
	void visualizationJct();
	void visualizationRadius();
	
	private :
	std::vector< std::vector <int> > connexite;
	int Maxorder;
	std::vector<int> ordre;
	std::vector<double> diametre_moy_node;
	std::vector<double> VAR_diametre_moy_node;
	std::vector<double> diametre_moy_order;
	std::vector<double> SD_diametre_moy_order;
	std::vector< std::vector <double> > connectivity_matrix;
	std::vector<double> lenght_moy_node;
	std::vector<double> lenght_moy_order;
	std::vector<double> area_node;
	std::vector<double> total_ratio;
	std::vector<double> major_ratio;
	std::vector<double> minor_ratio;
	std::vector<double> ass_ratio;
	std::vector< std::vector <double> > angle_node;
	std::vector<double> angle_m1;
	std::vector<double> angle_m2;
	std::vector<double> angle_12;
	std::vector<double> angle_beta;
	std::vector<double> angle_omega;
};

#endif /* _G_classif_h_ */



