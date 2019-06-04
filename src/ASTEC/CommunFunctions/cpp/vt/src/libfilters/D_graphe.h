/*****************************************************************************
* D_graphe.h
*
* par Manon Linder
*
*
*/


#ifndef _D_graphe_h_
#define _D_graphe_h_


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <cstring>
#include <vector> 

#include "vt_common.h"

// Vtk en commun
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

// Visu 2DGraph
#include <vtkCircularLayoutStrategy.h>
#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkGraphLayoutView.h>
#include <vtkIntArray.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkPoints.h>
#include <vtkIntArray.h>
#include <vtkLookupTable.h>
#include <vtkViewTheme.h>
#include <vtkCamera.h>

// Visu 3D
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPoints.h>
#include <vtkLine.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkVectorText.h>
#include <vtkActor.h>
#include <vtkFollower.h>


// Visu Tag
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include <vtkXMLPolyDataWriter.h>

/*------- CLASS D_COORD_ORI ----------*/

class D_coord_ori
{
	public :
	
	D_coord_ori();
	
	D_coord_ori(double ocoordz, double ocoordy, double ocoordx, int olabel);
	double getOCoordX() const;
	double getOCoordY() const;
	double getOCoordZ() const;
	int getOLabel() const; 
	void AfficheOPixel() const;
	
	protected :
	double o_coordz;
	double o_coordy;
	double o_coordx;
	int o_label;
};


/*------- CLASS D_PIXEL ----------*/

class D_pixel
{
	public :
	
	D_pixel();
	
	D_pixel(double coordz, double coordy, double coordx, double radius);
	
	void setValuePixel(double coordz, double coordy, double coordx, double radius);
	double getCoordX() const;
	double getCoordY() const;
	double getCoordZ() const; 
	double getRadius() const;
	void AffichePixel() const;
	
	protected :
	double p_coordz;
	double p_coordy;
	double p_coordx;
	double p_radius;
};


/*------- CLASS D_NODE ----------*/

class D_node
{
	public :
	
	D_node();
	
	void setTypeNode(std::string type);
	void resetNode(D_pixel const& new_pixel);
	std::string getTypeNode();
	void AfficheTypeNode(); 
	void InsertNextPixel(D_pixel const& new_pixel);
	void InsertNextPixelOri(D_coord_ori const& new_pixel);
	int getSizeNode() const;
	int getSizeONode() const;
	D_pixel getValuePixel(int no_pixel);
	D_coord_ori getValueOPixel(int no_pixel);
	void AfficheTheNode() const;
	void AfficheTheONode() const;
	void setHauteurNoeud(int mahauteur);
	int getHauteurNoeud();
	int getOrderNode();
	void setOrderNode(int nomordre);
	
	protected :
	
	std::vector<D_pixel> l_pixel;
	std::vector<D_coord_ori> l_opixel;
	std::string TypeNode;
	int no_ordre;
	int hauteur_noeud;
};


/*------- CLASS D_EDGE ----------*/

class D_edge
{
	public :
	
	D_edge();
	
	void InsertCoupleNode(int id_node1, int id_node2, int wght);
	int getValueNode1();
	int getValueNode2();
	int getWeight();
	void AfficheCoupleNode() const;
	
	protected :
	
	std::vector<int> e_node;
};


/*------- CLASS D_GRAPHE ----------*/
class D_graphe
{ 
	public : // Methodes
	
	D_graphe();
	
	// Fonctions relatives au graphe
	void InsertNextNode(D_node const& new_node);
	void InsertNextEdge(D_edge const& new_edge);
	int getIDNodeLastSegment();
	D_edge getCoupleNode(int no_edge);
	int getNumberEdge();
	int getNumberNode();
	int getNumberSegment();
	int getNumberSegTerminal();
	int getNumberJunction();
	int getNumberTotPixel();
	int getNumberTotPixelJunction();
	int getNumberTotPixelSegment();
	std::vector<int> getListeJct();
	std::vector<int> getListeTerm();
	std::vector<int> getListeSeg();
	std::vector<int> getListeDoubleSeg();
	int getValueDoubleSeg(int nb_i);
	int getValueJct(int nb_i);
	void setHauteurMaxA(int maHauteurA);
	int getHauteurMaxA();
	int getHauteurJunction(int no_jct);
	void AfficheNode();
	void AfficheEdge();
	void Liste_Adjacence();
	std::vector<D_node> getListeNode();
	
	// Fonctions du parcours du squelette
	
	void SkeletonProcess(int mahauteurA, D_node mynode, int no_jct, D_coord_ori mescoord, unsigned char ***sk_array, unsigned char ***dm_array, std::vector < int > mesdims);

	D_node SearchNeighbour(D_node mynode, D_pixel mespix, unsigned char ***sk_array, unsigned char ***dm_array);
	
	// Recherche du vrai Node 0 et de son point d'origine
	void Diametre_moyen();
	void Barycentre();
	
	void SelectPropag();
	void SelectAngle();
	void SelectDiam();
	
	D_coord_ori SearchTrueDeparturePropag(unsigned char ***sk_array);
	D_coord_ori SearchTrueDepartureAngle(unsigned char ***sk_array);
	D_coord_ori SearchTrueDepartureDiam(unsigned char ***sk_array);
	D_coord_ori SearchTrueDepartureDiamBary(unsigned char ***sk_array);
	D_coord_ori SearchTrueDepartureMixAll(unsigned char ***sk_array);	
	
	// Visualisation
	void visualization2DGrapheLabel();
	void visualization3DLabel();
	void visualizationTag();

	// Fonctions enregistrement du graphe
	void saveGraphe(char* name);
	void openGraphe(char* name);
	
	protected :// Attributs
	std::vector<D_node> l_node;
	std::vector<D_edge> l_edgeN;
	int id_lastsgnode;
	int NumberSegment;
	int NumberSegTerminal;
	int NumberJunction;
	int NumberTotPixel;
	int NumberTotPixelJunction;
	int NumberTotPixelSegment;
	std::vector<int> l_node_jct;
	std::vector<int> l_node_term;
	std::vector<int> l_node_seg;
	std::vector<int> l_node_dbseg;
	int HauteurA;
	std::vector<double> diametre_moyen;
	std::vector<int> tag;
	std::vector< std::vector <int> > l_Ad;
};

#endif /* _D_graphe_h_ */

