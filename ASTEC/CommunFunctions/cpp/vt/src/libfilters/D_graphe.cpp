/*****************************************************************************
* G_graphe.cpp
*
* par Manon Linder
*
*
*/

#include "D_graphe.h"


using namespace std;

/*------- CLASS D_COORD_ORI ----------*/

D_coord_ori::D_coord_ori() : o_coordz(0), o_coordy(0), o_coordx(0), o_label(0)
{
}
	
D_coord_ori::D_coord_ori(double ocoordz, double ocoordy, double ocoordx, int olabel) : o_coordz(ocoordz), o_coordy(ocoordy), o_coordx(ocoordx), o_label(olabel)
{
}

double D_coord_ori::getOCoordZ() const
{
	return o_coordz;
}

double D_coord_ori::getOCoordY() const
{
	return o_coordy;
}

double D_coord_ori::getOCoordX() const
{
	return o_coordx;
}

int D_coord_ori::getOLabel() const
{
	return o_label;
}

void D_coord_ori::AfficheOPixel() const
{
	cout << "OCoordZ :" << o_coordz << endl << "OCoordY :" << o_coordy << endl << "OCoordX :" << o_coordx << endl << "OLabel : " << o_label << endl;
}

/*------- CLASS D_PIXEL ----------*/

D_pixel::D_pixel() : p_coordz(0), p_coordy(0), p_coordx(0), p_radius(0)
{
}

D_pixel::D_pixel(double coordz, double coordy, double coordx, double radius) : p_coordz(coordz), p_coordy(coordy), p_coordx(coordx), p_radius(radius)
{
}

void D_pixel::setValuePixel(double coordz, double coordy, double coordx, double radius)
{
	p_coordz = coordz;
	p_coordy = coordy;
	p_coordx = coordx;
	p_radius = radius;
}

double D_pixel::getCoordZ() const
{
	return p_coordz;
}

double D_pixel::getCoordY() const
{
	return p_coordy;
}

double D_pixel::getCoordX() const
{
	return p_coordx;
}

double D_pixel::getRadius() const
{
	return p_radius;
}

void D_pixel::AffichePixel() const
{
	cout << "CoordZ :" << p_coordz << endl << "CoordY :" << p_coordy << endl << "CoordX :" << p_coordx << endl << "Radius :" << p_radius << endl;
}

/*------- CLASS D_NODE ----------*/

D_node::D_node() : l_pixel(0), l_opixel(0), TypeNode("tronc"), no_ordre(0), hauteur_noeud(0)
{
}

void D_node::setTypeNode(string type)
{
	TypeNode = type;
}

void D_node::resetNode(D_pixel const& new_pixel)
{
	l_pixel.clear();
	l_pixel.push_back(new_pixel);
}

string D_node::getTypeNode()
{
	return TypeNode;
}

void D_node::AfficheTypeNode()
{
	cout << TypeNode << endl;
}

void D_node::InsertNextPixel(D_pixel const& new_pixel)
{
	l_pixel.push_back(new_pixel);
}

void D_node::InsertNextPixelOri(D_coord_ori const& new_pixel)
{
	l_opixel.push_back(new_pixel);
}

int D_node::getSizeNode() const
{
	return l_pixel.size();
}

int D_node::getSizeONode() const
{
	return l_opixel.size();
}
	
D_pixel D_node::getValuePixel(int no_pixel)
{
	return l_pixel[no_pixel];
}	

D_coord_ori D_node::getValueOPixel(int no_pixel)
{
	return l_opixel[no_pixel];
}	

void D_node::AfficheTheNode() const
{
	for( unsigned int i = 0; i < l_pixel.size(); i++ )
	{
		cout << "Pixel Numero [" << i << "] =\tz : " << l_pixel[i].getCoordZ() << "\ty : " << l_pixel[i].getCoordY()  << "\tx : " <<l_pixel[i].getCoordX() << "\tradius : " << l_pixel[i].getRadius() << endl;
	}
}

void D_node::AfficheTheONode() const
{
	for( unsigned int i = 0; i < l_opixel.size(); i++ )
	{
		cout << "Pixel Numero [" << i << "] =\tz :" << l_opixel[i].getOCoordZ() << "\ty :" << l_opixel[i].getOCoordY()  << "\tx :" <<l_opixel[i].getOCoordX() << endl;
	}
}

void D_node::setHauteurNoeud(int mahauteur)
{
	hauteur_noeud = mahauteur;
}

int D_node::getHauteurNoeud()
{
	return hauteur_noeud;
}

int D_node::getOrderNode()
{
	return no_ordre;
}

void D_node::setOrderNode(int nomordre)
{
	no_ordre=nomordre;
}

/*------- CLASS D_EDGE ----------*/

D_edge::D_edge() : e_node(0)
{
}

void D_edge::InsertCoupleNode(int id_node1, int id_node2, int wght)
{
	e_node.push_back(id_node1);
	e_node.push_back(id_node2);
	e_node.push_back(wght);
}

int D_edge::getValueNode1()
{
	return e_node[0];
}

int D_edge::getValueNode2()
{
	return e_node[1];
}

int D_edge::getWeight()
{
	return e_node[2];
}

void D_edge::AfficheCoupleNode() const
{
	cout << "Couple Node =" << endl;
	for( unsigned int i = 0; i < 2; i++ )
	{
		cout << "node[" << i << "] = " << e_node[i] << endl;
	}
	cout << "Weight " << e_node[2] << endl;
}

/*------- CLASS D_GRAPHE ----------*/

D_graphe::D_graphe() : l_node(0), l_edgeN(0), id_lastsgnode(0), NumberSegment(0), NumberSegTerminal(0), NumberJunction(0), NumberTotPixel(0), NumberTotPixelJunction(0), NumberTotPixelSegment(0), l_node_jct(0), l_node_term(0), l_node_seg(0), l_node_dbseg(0), HauteurA(0), diametre_moyen(0), tag(0), l_Ad(0)
{
}	 

// Fonctions relatives au graphe
void D_graphe::InsertNextNode(D_node const& new_node)
{
	l_node.push_back(new_node);
}

void D_graphe::InsertNextEdge(D_edge const& new_edge)
{
	l_edgeN.push_back(new_edge);
}

int D_graphe::getIDNodeLastSegment()
{
	for(int i = 0; i < getNumberNode(); i++)
	{
		if (l_node[i].getTypeNode() == "segment")
		{
			id_lastsgnode = i;
		}
	}
	return id_lastsgnode;
}

D_edge D_graphe::getCoupleNode(int no_edge)
{
	return l_edgeN[no_edge];
}

int D_graphe::getNumberEdge()
{
	return l_edgeN.size();
}

int D_graphe::getNumberNode()
{
	return l_node.size();
}

int D_graphe::getNumberSegment()
{
	NumberSegment = 0;
	for(int i = 0; i < getNumberNode(); i++)
	{
		if (l_node[i].getTypeNode() == "segment")
		{
			NumberSegment++;
		}
	}
	return NumberSegment;
}

int D_graphe::getNumberSegTerminal()
{
	NumberSegTerminal = 0;
	for(int i = 0; i < getNumberNode(); i++)
	{
		if (l_node[i].getTypeNode() == "terminal")
		{
			NumberSegTerminal++;
		}
	}
	return NumberSegTerminal;
}

int D_graphe::getNumberJunction()
{
	NumberJunction = 0;
	for(int i = 0; i < getNumberNode(); i++)
	{
		if (l_node[i].getTypeNode() == "junction")
		{
			NumberJunction++;
		}
	}
	return NumberJunction;
}

int D_graphe::getNumberTotPixel()
{
	NumberTotPixel = 0;
	for(int i = 0; i < getNumberNode(); i++)
	{
		NumberTotPixel+=l_node[i].getSizeNode();
	}
	return NumberTotPixel;
}

int D_graphe::getNumberTotPixelJunction()
{
	NumberTotPixelJunction = 0;
	for(int i = 0; i < getNumberNode(); i++)
	{
		if (l_node[i].getTypeNode() == "junction")
		{
			NumberTotPixelJunction+=l_node[i].getSizeNode();
		}
	}
	return NumberTotPixelJunction;
}

int D_graphe::getNumberTotPixelSegment()
{
	NumberTotPixelSegment = 0;
	for(int i = 0; i < getNumberNode(); i++)
	{
		if (l_node[i].getTypeNode() != "junction")
		{
			NumberTotPixelSegment+=l_node[i].getSizeNode();
		}
	}
	return NumberTotPixelSegment;
}


std::vector<int> D_graphe::getListeJct()
{ 
	l_node_jct.clear();
	for(int i = 0; i < getNumberNode(); i++)
	{
		if (l_node[i].getTypeNode() == "junction")
		{
			l_node_jct.push_back(i);
		}
	}
	return l_node_jct;
}

std::vector<int> D_graphe::getListeTerm()
{ 
	l_node_term.clear();
	for(int i = 0; i < getNumberNode(); i++)
	{
		if (l_node[i].getTypeNode() == "terminal")
		{
			l_node_term.push_back(i);
		}
	}
	return l_node_term;
}

std::vector<int> D_graphe::getListeSeg()
{ 
	l_node_seg.clear();
	for(int i = 0; i < getNumberNode(); i++)
	{
		if (l_node[i].getTypeNode() == "segment")
		{
			l_node_seg.push_back(i);
		}
	}
	return l_node_seg;
}

std::vector<int> D_graphe::getListeDoubleSeg()
{ 
	l_node_dbseg.clear();
	for(int i = 0; i < getNumberNode(); i++)
	{
		if (l_node[i].getTypeNode() != "junction")
		{
			l_node_dbseg.push_back(i);
		}
	}
	return l_node_dbseg;
}

int D_graphe::getValueDoubleSeg(int nb_i)
{
	return l_node_dbseg[nb_i];
}

int D_graphe::getValueJct(int nb_i)
{
	return l_node_jct[nb_i];
}

void D_graphe::setHauteurMaxA(int maHauteurA)
{
	if (maHauteurA > HauteurA)
	{
		HauteurA=maHauteurA;
	}
}

int D_graphe::getHauteurMaxA()
{
	return HauteurA;
}

int D_graphe::getHauteurJunction(int no_jct)
{
	return l_node[l_node_jct[no_jct]].getHauteurNoeud();
}

void D_graphe::AfficheNode()
{
	for(int i = 0; i < getNumberNode(); i++)
	{
		cout <<i<<endl;
		l_node[i].AfficheTypeNode();
		if (l_node[i].getTypeNode() == "junction")
		{
			cout << l_node[i].getHauteurNoeud() << endl;
		}
		cout << l_node[i].getSizeNode() <<endl;
		l_node[i].AfficheTheNode();
	}
}

void D_graphe::AfficheEdge()
{
	for(int i = 0; i < getNumberEdge(); i++)
	{
		cout <<i<<endl;
		l_edgeN[i].AfficheCoupleNode();
	}
}

void D_graphe::Liste_Adjacence()
{
	getListeJct();
	vector<int> l_interm;
	
	for(unsigned int i = 0; i < l_node_jct.size(); i++)
	{
		l_interm.push_back(l_node_jct[i]); // Numero jonction
		for(int k = 0; k < getNumberEdge(); k++)
		{
			if (getCoupleNode(k).getValueNode1() == l_node_jct[i])
			{
				l_interm.push_back(getCoupleNode(k).getValueNode2()); // Branches m et filles
			}
		}
		l_Ad.push_back(l_interm);
		l_interm.clear();
	}
}

std::vector<D_node> D_graphe::getListeNode()
{
	return l_node;
}

/*****************************************************************************************************************/
// Fonctions du parcours du squelette

D_node D_graphe::SearchNeighbour(D_node mynode, D_pixel mespix, unsigned char ***sk_array, unsigned char ***dm_array)
{
	int x, y, z;
	int zscan;
	int yscan;
	int xscan;

	/*--- Recuperation des pixels d'origine ---*/
	zscan = mespix.getCoordZ();
	yscan = mespix.getCoordY();
	xscan = mespix.getCoordX();
	
	for (z=zscan-1; z<=zscan+1; z++)
	for (y=yscan-1; y<=yscan+1; y++)
	for (x=xscan-1; x<=xscan+1; x++)
	{
		sk_array[zscan][yscan][xscan] = (unsigned char)0;	// valeur pixel a 0 // petit poucet
		if (sk_array[z][y][x] == (unsigned char)255)
		{
			D_pixel the_pix(z, y, x, dm_array[z][y][x]);
			sk_array[z][y][x] = (unsigned char)0;
			mynode.InsertNextPixel(the_pix);
		}
	}
	return(mynode);
}

void D_graphe::SkeletonProcess(int maHauteurA, D_node the_node, int no_jct, D_coord_ori mescoord, unsigned char ***sk_array, unsigned char ***dm_array, vector < int > mesdims)
{
	int zo,yo,xo;
	int zscan;
	int yscan;
	int xscan;
	int x, y, z;
	int l_sg;
	
	//mescoord.AfficheOPixel();
	/*--- Recuperation des pixels d'origine ---*/
	zo = mescoord.getOCoordZ();
	yo = mescoord.getOCoordY();
	xo = mescoord.getOCoordX();
	
	if (getNumberNode()>2)
	{	
		if (the_node.getTypeNode() == "junction")
		{
			// Creation nouveau noeud
			D_node the_node;
		}
	}
	
	/*--- PARCOURS SEGMENT ---*/
	// Enregistrement donnees pixel
	D_pixel the_pix(zo, yo, xo, (unsigned int)dm_array[zo][yo][xo]);
	D_pixel the_pix_next_event;
	
	// Insertion dans noeud actuel
	the_node.InsertNextPixel(the_pix);
	
	/*--- Parcours segment ---*/
	xscan = xo;
	yscan = yo;
	zscan = zo;
	sk_array[zscan][yscan][xscan] = (unsigned char)100;	// petit poucet
	
	while(sk_array[zscan][yscan][xscan] == (unsigned char)100)
	{
		sk_array[zscan][yscan][xscan] = (unsigned char)0;	// petit poucet	
		for (z=zscan-1; z<=zscan+1; z++)
		for (y=yscan-1; y<=yscan+1; y++)
		for (x=xscan-1; x<=xscan+1; x++)
		{
			if (sk_array[z][y][x] == (unsigned char)100) // voisin=suite segment
			{ 
				D_pixel the_pix(z, y, x, dm_array[z][y][x]);
				the_node.InsertNextPixel(the_pix);
				//the_node.AfficheTheNode();
				// passage au pixel suivant
				xscan=x; 
				yscan=y;
				zscan=z;
				// Sortie boucle for
				x=xscan+4; 
				y=yscan+4;
				z=zscan+4;
			}
			if( z < mesdims[0] && y < mesdims[1] && x < mesdims[2] )
			{
				if (sk_array[z][y][x] == (unsigned char)200 || sk_array[z][y][x] == (unsigned char)255)
				{
					// passage au pixel suivant
					the_pix_next_event.setValuePixel(z, y, x, dm_array[z][y][x]);
					// Sortie boucle for
					x=xscan+4; 
					y=yscan+4;
					z=zscan+4;
				}
			}
		}
	}
	//the_node.AfficheTheNode();
	
	xscan=the_pix_next_event.getCoordX(); 
	yscan=the_pix_next_event.getCoordY();
	zscan=the_pix_next_event.getCoordZ();
	
	//cout << (unsigned int)sk_array[zscan][yscan][xscan] << endl;
	
	/*--- CAS segment TERMINAL---*/			
	if (sk_array[zscan][yscan][xscan] == (unsigned char)200) // voisin = point terminal
	{ 	
		
		D_pixel the_pix(zscan, yscan, xscan, dm_array[zscan][yscan][xscan]);
		the_node.InsertNextPixel(the_pix);
		the_node.setTypeNode("terminal");
		//the_node.AfficheTypeNode();
		//the_node.AfficheTheNode();
				
		// Insertion du noeud dans le graphe
		InsertNextNode(the_node);
		//cout << "Number node " << getNumberNode() << endl;
		if (getNumberNode() > 2)
		{			
			// Creation arc entre ce noeud et le precedent
			D_edge the_edge;
			the_edge.InsertCoupleNode(no_jct,getNumberNode()-1, 1);
			//the_edge.AfficheCoupleNode();
			InsertNextEdge(the_edge);
		}
	}
	
	/*--- CAS segment PUIS jonction ---*/					
	if (sk_array[zscan][yscan][xscan] == (unsigned char)255) // voisin = point jct
	{
		// SCAN DES POINTS DU "NOEUD JONCTION"		
		// Creation noeud intermediaire
		D_node the_node_interm;		
		
		// Enregistrement donnees pixel
		D_pixel the_pix(zscan, yscan, xscan, dm_array[zscan][yscan][xscan]);

		// Insertion dans noeud intermediaire
		the_node_interm.InsertNextPixel(the_pix);
		//the_node_interm.AfficheTheNode();
		// valeur pixel a 0 // petit poucet
		sk_array[zscan][yscan][xscan] = (unsigned char)0;	

		// Recherche point jonction autour du 1er point de jonction
		the_node_interm = SearchNeighbour(the_node_interm, the_pix, sk_array, dm_array);

		int j = 0;
		while(j<3)
		{
			// Recuperation du nombre de point dans d_points
			int Nbre_pt_junction = the_node_interm.getSizeNode();
			
			// Recuperation des coord dans d_pixels
			for(int i = 0; i < Nbre_pt_junction; i++)
			{
				the_node_interm = SearchNeighbour(the_node_interm, the_node_interm.getValuePixel(i), sk_array, dm_array);
			}
			j++;
		}
		//cout << "***************************"<< endl;
		//cout << "Jonction" << endl;
		//the_node_interm.AfficheTheNode();
		//cout << endl;
		
		// SCAN DES POINTS APRES LE "NOEUD JONCTION"
		/*--- Stockage points origines ---*/
		// Recuperation du nombre de point jonction dans d_points
		int Nbre_pt_junction2 = the_node_interm.getSizeNode();
		
		// Creation pseudo-noeud stockage points originies
		D_node futur_seg;	
		D_node futur_term;		
		// Recuperation des coord des points d'origine dans d_pixels
		for(int i = 0; i < Nbre_pt_junction2; i++)
		{
			zscan = the_node_interm.getValuePixel(i).getCoordZ();
			yscan = the_node_interm.getValuePixel(i).getCoordY();
			xscan = the_node_interm.getValuePixel(i).getCoordX();
			sk_array[zscan][yscan][xscan] = (unsigned char)0;
			
			for (z=zscan-1; z<=zscan+1; z++)
			for (y=yscan-1; y<=yscan+1; y++)
			for (x=xscan-1; x<=xscan+1; x++)
			{
				if (sk_array[z][y][x] == (unsigned char)100)
				{
					D_coord_ori mescoordfuturseg(z,y,x,sk_array[z][y][x]);
					//mescoordfuturseg.AfficheOPixel();
					futur_seg.InsertNextPixelOri(mescoordfuturseg);
				}
				if (sk_array[z][y][x] == (unsigned char)200) // cas pixel terminal
				{
					D_coord_ori mescoordfuturterm(z,y,x,sk_array[z][y][x]);
					//mescoordfuturterm.AfficheOPixel();
					futur_term.InsertNextPixelOri(mescoordfuturterm);
				}
			}
		}
		//cout << "Segment" << endl;
		//futur_seg.AfficheTheONode();
		//cout << "Terminal" << endl;
		//futur_term.AfficheTheONode();
		
		// CAS 1 : plusieurs points SEGMENTS -> Pas de remaniement de l'arbre
		if (futur_seg.getSizeONode() > 1 && futur_term.getSizeONode() >= 0)
		{
			// FERMETURE DU SEGMENT
			the_node.setTypeNode("segment");
			//the_node.AfficheTypeNode();
			//the_node.AfficheTheNode();

			// Insertion du noeud dans le graphe
			InsertNextNode(the_node);
			//cout << "Number node " << getNumberNode() << endl;
			if (getNumberNode() > 2)
			{			
				// Creation arc entre ce noeud et le precedent
				D_edge the_edge;
				the_edge.InsertCoupleNode(no_jct,getNumberNode()-1,1);
				//the_edge.AfficheCoupleNode();
				InsertNextEdge(the_edge);
			}
			// FERMETURE DE LA JONCTION
			maHauteurA+=1;
			setHauteurMaxA(maHauteurA);
			the_node_interm.setHauteurNoeud(maHauteurA);
			the_node_interm.setTypeNode("junction");
			//the_node_interm.AfficheTypeNode();
			//the_node_interm.AfficheTheNode();
			// Insertion du noeud dans le graphe
			InsertNextNode(the_node_interm);
			//cout << "Number node " << getNumberNode() << endl;
	
			if (getNumberNode() == 2)
			{
				D_edge the_edge;
				the_edge.InsertCoupleNode(1,0,1);
				//the_edge.AfficheCoupleNode();
				InsertNextEdge(the_edge);
			}
		
			if (getNumberNode() > 2)
			{			
				// Creation arc entre ce noeud et le precedent
				D_edge the_edge;
				no_jct = getNumberNode()-1;
				l_sg = getIDNodeLastSegment();
				the_edge.InsertCoupleNode(no_jct,l_sg,1);
				//the_edge.AfficheCoupleNode();
				InsertNextEdge(the_edge);
			}
			// RECURSIVITE
			for(int i = 0; i < Nbre_pt_junction2; i++)
			{
				zscan = the_node_interm.getValuePixel(i).getCoordZ();
				yscan = the_node_interm.getValuePixel(i).getCoordY();
				xscan = the_node_interm.getValuePixel(i).getCoordX();
				sk_array[zscan][yscan][xscan] = (unsigned char)0;
			
				for (z=zscan-1; z<=zscan+1; z++)
				for (y=yscan-1; y<=yscan+1; y++)
				for (x=xscan-1; x<=xscan+1; x++)
				{
					if (sk_array[z][y][x] == (unsigned char)100)
					{
						D_node the_node;
						D_coord_ori mescoord(z,y,x,sk_array[z][y][x]);
						SkeletonProcess(maHauteurA, the_node, no_jct, mescoord, sk_array, dm_array, mesdims);
					}
				}
			}
		}
		
		// CAS 2 : un point SEGMENTS et points TERMINAUX -> Remaniement de l'arbre : Fusion segment + jonction + segment
		if (futur_seg.getSizeONode() == 1 && futur_term.getSizeONode() >= 0)
		{
			
			// INSERTION POINT JONCTION DANS SEGMENT
			for(int i = 0; i < Nbre_pt_junction2; i++)
			{
				zscan = the_node_interm.getValuePixel(i).getCoordZ();
				yscan = the_node_interm.getValuePixel(i).getCoordY();
				xscan = the_node_interm.getValuePixel(i).getCoordX();
				// Enregistrement donnees pixel
				D_pixel the_pix(zscan, yscan, xscan, dm_array[zscan][yscan][xscan]);

				// Insertion dans noeud segment
				the_node.InsertNextPixel(the_pix);
			}
			the_node.setTypeNode("segment");
			// RECURSIVITE
			setHauteurMaxA(maHauteurA);
			D_coord_ori mescoordf(futur_seg.getValueOPixel(0).getOCoordZ(),futur_seg.getValueOPixel(0).getOCoordY(),futur_seg.getValueOPixel(0).getOCoordX(),futur_seg.getValueOPixel(0).getOLabel());
			SkeletonProcess(maHauteurA, the_node, no_jct, mescoordf, sk_array, dm_array, mesdims);			
		}
		
		// CAS 3 : que des points TERMINAUX -> Fusion segment + jonction
		if (futur_seg.getSizeONode() == 0 && futur_term.getSizeONode() > 0)
		{
			
			// INSERTION POINT JONCTION DANS SEGMENT
			for(int i = 0; i < Nbre_pt_junction2; i++)
			{
				zscan = the_node_interm.getValuePixel(i).getCoordZ();
				yscan = the_node_interm.getValuePixel(i).getCoordY();
				xscan = the_node_interm.getValuePixel(i).getCoordX();
				// Enregistrement donnees pixel
				D_pixel the_pix(zscan, yscan, xscan, dm_array[zscan][yscan][xscan]);

				// Insertion dans noeud segment
				the_node.InsertNextPixel(the_pix);
			}

			// FERMETURE DU SEGMENT
			the_node.setTypeNode("terminal");
			//the_node.AfficheTypeNode();
			//the_node.AfficheTheNode();

			// Insertion du noeud dans le graphe
			InsertNextNode(the_node);
			//cout << "Number node " << getNumberNode() << endl;
			if (getNumberNode() > 2)
			{			
				// Creation arc entre ce noeud et le precedent
				D_edge the_edge;
				the_edge.InsertCoupleNode(no_jct,getNumberNode()-1,1);
				//the_edge.AfficheCoupleNode();
				InsertNextEdge(the_edge);
			}
		}
	}
}

/*******************************************************************************************/
// Recherche du vrai Node 0 et de son point d'origine
void D_graphe::Diametre_moyen()
{
	double somme_diametre = 0 ;
	double size_voxel = 16 ; // en micro metre
	
	for(int i = 0; i < getNumberNode(); i++)
	{
		for(int j = 0 ; j < l_node[i].getSizeNode(); j++)
		{
			somme_diametre += 2*l_node[i].getValuePixel(j).getRadius() * size_voxel;
		}
		diametre_moyen.push_back(somme_diametre/l_node[i].getSizeNode());
		somme_diametre = 0;
	}
}

void D_graphe::Barycentre()
{
	getListeJct();
	
	for(unsigned int i = 0; i < l_node_jct.size(); i++)
	{
		double z = 0;
		double y = 0;
		double x = 0;
		double rad = 0;
		int compte = 0;
		if (l_node[l_node_jct[i]].getSizeNode() > 1)
		{
			for(int m = 0; m < l_node[l_node_jct[i]].getSizeNode(); m++)
			{
				z += l_node[l_node_jct[i]].getValuePixel(m).getCoordZ();
				y += l_node[l_node_jct[i]].getValuePixel(m).getCoordY();
				x += l_node[l_node_jct[i]].getValuePixel(m).getCoordX();
				rad += l_node[l_node_jct[i]].getValuePixel(m).getRadius();
					
				compte+=1;
			}
			z = z/compte;
			y = y/compte;
			x = x/compte;
			rad = rad/compte;
			D_pixel the_pix(z, y, x, rad);
			l_node[l_node_jct[i]].resetNode(the_pix);
		}
	}
}

void D_graphe::SelectPropag()
{
	vector<int> branche;
	vector<int> jonction;
        int untag_br = -1;
        int next_jct = -1;
	unsigned int compte;
	double moyen_br;
	
	for(int i = 0; i < getNumberNode(); i++)
	{
		if (l_node[i].getTypeNode() == "junction")
		{
			// Recuperation branches de la jonction
			for(int k = 0; k < getNumberEdge(); k++)
			{
				if (getCoupleNode(k).getValueNode1() == i)
				{
					branche.push_back(getCoupleNode(k).getValueNode2());
				}
			}
			// Compte le nombre de branche tag = 1
			for (unsigned int j = 0 ; j < branche.size() ; j++)
			{
				if (tag[branche[j]] == 1)
				{
					compte +=1;
					moyen_br += diametre_moyen[branche[j]];
				}
				if (tag[branche[j]] == 0)
				{
					untag_br = branche[j];
					for(int k = 0; k < getNumberEdge(); k++)
					{
						if (getCoupleNode(k).getValueNode2() == untag_br)
						{
							jonction.push_back(getCoupleNode(k).getValueNode1());
						}
					}
					for (unsigned int n = 0 ; n < jonction.size() ; n++)
					{
						if (jonction[n] != i)
						{
							next_jct = jonction[n];
						}
					}
				}
			}
			// Si 1 branche pas tag parmi toutes
			if (compte == branche.size()-1)
			{
				moyen_br = moyen_br/compte;
				if (/*diametre_moyen[untag_br] >= moyen_br || */diametre_moyen[next_jct] >= diametre_moyen[i])
				{
					tag[untag_br] = 1;
					tag[i] = 1;
					SelectPropag();
				}
			}
			branche.clear();
			jonction.clear();
			compte = 0 ;
			moyen_br = 0;
		}
	}
}

void D_graphe::SelectAngle()
{
	vector<int> branche;
	vector<int> jonction;

	double CoordZ;
	double CoordY;
	double CoordX;
	
	vector<double> vect_m;
	vector<double> vect_f1;
	vector<double> vect_f2;
	
	double norme_m;
        double norme_f1 = 0.0;
        double norme_f2 = 0.0;
				
	vector<double> prod_scal;
	vector<double> angle_interm;
	//const double mon_pi = 3.14159265358979323846;
	
	tag.clear();
	tag.resize(getNumberNode());
	
	for (unsigned int i = 0 ; i < l_node_jct.size() ; i++)
	{
		branche.push_back(i);
		for(int k = 0; k < getNumberEdge(); k++)
		{
			if (getCoupleNode(k).getValueNode1() == l_node_jct[i])
			{
					branche.push_back(getCoupleNode(k).getValueNode2());
			}
		}

		/// Calcul d'angle
		if (branche.size() == 4)
		{
			// Cas branche mere
			CoordZ = l_node[branche[1]].getValuePixel(l_node[branche[1]].getSizeNode()-1).getCoordZ();
			CoordY = l_node[branche[1]].getValuePixel(l_node[branche[1]].getSizeNode()-1).getCoordY();
			CoordX = l_node[branche[1]].getValuePixel(l_node[branche[1]].getSizeNode()-1).getCoordX();
		
			// Calcul Vecteur + norme
			vect_m.push_back(CoordZ-l_node[branche[0]].getValuePixel(0).getCoordZ());
			vect_m.push_back(CoordY-l_node[branche[0]].getValuePixel(0).getCoordY());
			vect_m.push_back(CoordX-l_node[branche[0]].getValuePixel(0).getCoordX());
			norme_m = sqrt(pow(vect_m[0],2)+pow(vect_m[1],2)+pow(vect_m[2],2));
						
			// Cas branches filles
			for (int f = 2; f<4 ; f++)
			{
				CoordZ = l_node[branche[f]].getValuePixel(l_node[branche[f]].getSizeNode()-1).getCoordZ();
				CoordY = l_node[branche[f]].getValuePixel(l_node[branche[f]].getSizeNode()-1).getCoordY();
				CoordX = l_node[branche[f]].getValuePixel(l_node[branche[f]].getSizeNode()-1).getCoordX();

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
			angle_interm.push_back(acos(prod_scal[0]/(norme_m*norme_f1)));
		
			// CAS M2
			angle_interm.push_back(acos(prod_scal[1]/(norme_m*norme_f2)));
		
			// CAS 12
			angle_interm.push_back(acos(prod_scal[2]/(norme_f1*norme_f2)));
				
			// Branche mere reellement mere
			if(angle_interm[2]<angle_interm[1] && angle_interm[2]<angle_interm[0])
			{
				//cout << "m " << endl;
				tag[branche[1]]=1;
				tag[branche[2]]=0;
				tag[branche[3]]=0;			
			}
			// Branche 1 = branche mere
			if(angle_interm[1]<angle_interm[2] && angle_interm[1]<angle_interm[0])
			{
				//cout << "f2 " << endl;
				tag[branche[1]]=0;
				tag[branche[2]]=0;
				tag[branche[3]]=1;	
			}
			// Branche 2 = branche mere
			if(angle_interm[0]<angle_interm[1] && angle_interm[0]<angle_interm[2])
			{
				//cout << "f1 " << endl;
				tag[branche[1]]=0;
				tag[branche[2]]=1;
				tag[branche[3]]=0;	
			}
		}
		branche.clear();
		vect_m.clear();
		vect_f1.clear();
		vect_f2.clear();
		prod_scal.clear();
		angle_interm.clear();
	}

	for(int i = 0; i < getNumberNode(); i++)
	{		
		if(tag[i]==1 && l_node[i].getTypeNode() != "terminal")
		{
			tag[i]=0;
		}
	}
}


void D_graphe::SelectDiam()
{
	tag.clear();
	tag.resize(getNumberNode());
	vector<int> branche;
	vector<int> jonction;
	vector<int> ssjonction;
	vector <vector<double> > diametre3;
	vector < double > intermed;
        double calculjct = 0.0;
	
	Diametre_moyen();
	getListeJct(); 
	
	// Moyenne des jonctions
	for(unsigned int j = 0; j < l_node_jct.size(); j++)
	{
		calculjct += diametre_moyen[l_node_jct[j]];
		
	}
	calculjct = calculjct/l_node_jct.size();

	// Selection des jonctions dont moyenne propre sup a moyenne jct et qui possede seg term
	for(unsigned int j = 0; j < l_node_jct.size(); j++)
	{
		if (diametre_moyen[l_node_jct[j]]>(calculjct*3.3/2))
		{
			ssjonction.push_back(l_node_jct[j]);
		}
	}

	for(unsigned int s=0; s < ssjonction.size();s++)
	{
		for(int k = 0; k < getNumberEdge(); k++)
		{
			if (getCoupleNode(k).getValueNode1() == ssjonction[s])
			{
				if(l_node[getCoupleNode(k).getValueNode2()].getTypeNode() == "terminal")
				{
					tag[getCoupleNode(k).getValueNode2()]=1;
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////
D_coord_ori D_graphe::SearchTrueDeparturePropag(unsigned char ***sk_array)
{
	vector<int> tag_jonction;
	std::vector<int> list_propag;
	//// Pre Marquage
	tag.clear();
	tag.resize(getNumberNode());

	tag[0] = 1;
	
	for(int i = 0; i < getNumberNode(); i++)
	{
		if (l_node[i].getTypeNode() == "terminal")
		{
			tag[i] = 1;
		}
	}
	SelectPropag();
	
	for(int i = 0; i < getNumberNode(); i++)
	{
		if (tag[i] == 0 && l_node[i].getTypeNode() == "junction")
		{
			tag_jonction.push_back(i);
		}
	}
	for(unsigned int s=0; s < tag_jonction.size();s++)
	{
		for(int k = 0; k < getNumberEdge(); k++)
		{
			if (getCoupleNode(k).getValueNode1() == tag_jonction[s])
			{
				if(l_node[getCoupleNode(k).getValueNode2()].getTypeNode() == "terminal")
				{
					list_propag.push_back(getCoupleNode(k).getValueNode2());
				}
			}
		}
	}
	
	int mon_pix;
	int z;
	int y;
	int x;
	int ma_br =0	;
	if(ma_br == 0)
	{
		mon_pix = 0;
	}
	else 
	{
		mon_pix = l_node[ma_br].getSizeNode()-1;
	}
	z = l_node[ma_br].getValuePixel(mon_pix).getCoordZ();
	y = l_node[ma_br].getValuePixel(mon_pix).getCoordY();
	x = l_node[ma_br].getValuePixel(mon_pix).getCoordX();

	
	D_coord_ori Mescoord_ori(z,y,x,sk_array[z][y][x]);
	return(Mescoord_ori);
}


D_coord_ori D_graphe::SearchTrueDepartureAngle(unsigned char ***sk_array)
{
	SelectAngle();
	std::vector<int> list_angle;
	
	for(int i = 0; i < getNumberNode(); i++)
	{		
		if(tag[i]==1)
		{
			list_angle.push_back(i);
		}
	}
	
	int mon_pix;
	int z;
	int y;
	int x;
	int ma_br =0	;
	if(ma_br == 0)
	{
		mon_pix = 0;
	}
	else 
	{
		mon_pix = l_node[ma_br].getSizeNode()-1;
	}
	z = l_node[ma_br].getValuePixel(mon_pix).getCoordZ();
	y = l_node[ma_br].getValuePixel(mon_pix).getCoordY();
	x = l_node[ma_br].getValuePixel(mon_pix).getCoordX();

	
	D_coord_ori Mescoord_ori(z,y,x,sk_array[z][y][x]);
	return(Mescoord_ori);
}

D_coord_ori D_graphe::SearchTrueDepartureDiam(unsigned char ***sk_array)
{
	std::vector<int> list_diam;
	std::vector<int> list_jdiam;

									
	SelectDiam();
	for(int i = 0; i < getNumberNode(); i++)
	{		
		if(tag[i]==1)
		{
			list_diam.push_back(i);
		}
	}
	
	for(unsigned int s=0; s < list_diam.size();s++)
	{
		for(int k = 0; k < getNumberEdge(); k++)
		{
			if (getCoupleNode(k).getValueNode2() == list_diam[s])
			{
				list_jdiam.push_back(getCoupleNode(k).getValueNode1());
			}
		}
	}
	
	int mon_pix;
	int z;
	int y;
	int x;
	int ma_br =0	;
	if(ma_br == 0)
	{
		mon_pix = 0;
	}
	else 
	{
		mon_pix = l_node[ma_br].getSizeNode()-1;
	}
	z = l_node[ma_br].getValuePixel(mon_pix).getCoordZ();
	y = l_node[ma_br].getValuePixel(mon_pix).getCoordY();
	x = l_node[ma_br].getValuePixel(mon_pix).getCoordX();

	
	D_coord_ori Mescoord_ori(z,y,x,sk_array[z][y][x]);
	return(Mescoord_ori);
}


D_coord_ori D_graphe::SearchTrueDepartureDiamBary(unsigned char ***sk_array)
{
	std::vector<int> list_diam;	
	std::vector<int> list_jdiam;
	std::vector<double> list_longueur;	
	
	SelectDiam();
	
	for(int i = 0; i < getNumberNode(); i++)
	{		
		if(tag[i]==1)
		{
			list_diam.push_back(i);
		}
	}
	
	for(unsigned int s=0; s < list_diam.size();s++)
	{
		for(int k = 0; k < getNumberEdge(); k++)
		{
			if (getCoupleNode(k).getValueNode2() == list_diam[s])
			{
				list_jdiam.push_back(getCoupleNode(k).getValueNode1());
			}
		}
	}
		
	// Barycentre des jonctions
	double zb = 0;
	double yb = 0;
	double xb = 0;
	double radb = 0;
	int compte = 0;
	for(int i = 0; i < getNumberNode(); i++)
	{
		if (l_node[i].getTypeNode() == "junction")
		{
			zb += l_node[i].getValuePixel(0).getCoordZ();
			yb += l_node[i].getValuePixel(0).getCoordY();
			xb += l_node[i].getValuePixel(0).getCoordX();
			radb += l_node[i].getValuePixel(0).getRadius();
			
			compte+=1;
		}
	}
	zb = zb/compte;
	yb = yb/compte;
	xb = xb/compte;
	radb = radb/compte;
	D_pixel the_pix(zb, yb, xb, radb);
	
	double longeur_interm = 0;
	
	// Calcul distance par rapport au barycentre
	for(unsigned int i = 0; i < list_diam.size(); i++)
	{
		longeur_interm = sqrt(pow(the_pix.getCoordZ()-l_node[list_jdiam[i]].getValuePixel(0).getCoordZ() ,2) + pow(the_pix.getCoordY()-l_node[list_jdiam[i]].getValuePixel(0).getCoordY() ,2) + pow(the_pix.getCoordX()-l_node[list_jdiam[i]].getValuePixel(0).getCoordX() ,2));
		/*cout << longeur_interm << endl;*/
		list_longueur.push_back(longeur_interm);
	}
	int mon_max =0 ;
	int mon_it = 0;
	
	for (unsigned int i = 0; i < list_longueur.size(); i++)
	{
		if (list_longueur[i]>mon_max)
		{
			mon_max = list_longueur[i];
			mon_it = i;
		}
	}
	int mon_itit = 0 ;
	for(unsigned int s=0; s < list_diam.size();s++)
	{
		for(int k = 0; k < getNumberEdge(); k++)
		{
			if (getCoupleNode(k).getValueNode1() == list_jdiam[mon_it])
			{	
				
				if(l_node[getCoupleNode(k).getValueNode2()].getTypeNode() == "terminal")
				{
					mon_itit = getCoupleNode(k).getValueNode2();
				}
			}
		}
	}

	tag.clear();
	tag.resize(getNumberNode());

	tag[mon_itit] = 1;
	
	int mon_pix;
	int z;
	int y;
	int x;

	if(mon_itit == 0)
	{
		mon_pix = 0;
	}
	else 
	{
		mon_pix = l_node[mon_itit].getSizeNode()-1;
	}
	z = l_node[mon_itit].getValuePixel(mon_pix).getCoordZ();
	y = l_node[mon_itit].getValuePixel(mon_pix).getCoordY();
	x = l_node[mon_itit].getValuePixel(mon_pix).getCoordX();

	
	D_coord_ori Mescoord_ori(z,y,x,sk_array[z][y][x]);
	return(Mescoord_ori);
}

D_coord_ori D_graphe::SearchTrueDepartureMixAll(unsigned char ***sk_array)
{
	Diametre_moyen();
	getListeJct();
	std::vector<int> list_angle;
	std::vector<int> list_propag;
	std::vector<int> list_diam;
		std::vector<int> list_jdiam;
		std::vector<int> list_merge;
		std::vector<int> list_merge2;
				std::vector<int> list_final;
								std::vector<double> list_longueur;
	vector<int> tag_jonction;
	vector<int> term_jonction;
		int mycount = 0;
		
		
	SelectAngle();
	for(int i = 0; i < getNumberNode(); i++)
	{		
		if(tag[i]==1)
		{
			list_angle.push_back(i);
		}
	}
	
	//// Pre Marquage
	tag.clear();
	tag.resize(getNumberNode());

	tag[0] = 1;
	
	for(int i = 0; i < getNumberNode(); i++)
	{
		if (l_node[i].getTypeNode() == "terminal")
		{
			tag[i] = 1;
		}
	}
	SelectPropag();
	
	for(int i = 0; i < getNumberNode(); i++)
	{
		if (tag[i] == 0 && l_node[i].getTypeNode() == "junction")
		{
			tag_jonction.push_back(i);
		}
	}
	for(unsigned int s=0; s < tag_jonction.size();s++)
	{
		for(int k = 0; k < getNumberEdge(); k++)
		{
			if (getCoupleNode(k).getValueNode1() == tag_jonction[s])
			{
				if(l_node[getCoupleNode(k).getValueNode2()].getTypeNode() == "terminal")
				{
					list_propag.push_back(getCoupleNode(k).getValueNode2());
				}
			}
		}
	}
	
	SelectDiam();
	for(int i = 0; i < getNumberNode(); i++)
	{		
		if(tag[i]==1)
		{
			list_diam.push_back(i);
		}
	}
	
	for(unsigned int s=0; s < list_diam.size();s++)
	{
		for(int k = 0; k < getNumberEdge(); k++)
		{
			if (getCoupleNode(k).getValueNode2() == list_diam[s])
			{
				list_jdiam.push_back(getCoupleNode(k).getValueNode1());
			}
		}
	}
	
	// Barycentre des jonctions
	double zb = 0;
	double yb = 0;
	double xb = 0;
	double radb = 0;
	int compte = 0;
	for(int i = 0; i < getNumberNode(); i++)
	{
		if (l_node[i].getTypeNode() == "junction")
		{
			zb += l_node[i].getValuePixel(0).getCoordZ();
			yb += l_node[i].getValuePixel(0).getCoordY();
			xb += l_node[i].getValuePixel(0).getCoordX();
			radb += l_node[i].getValuePixel(0).getRadius();
			
			compte+=1;
		}
	}
	zb = zb/compte;
	yb = yb/compte;
	xb = xb/compte;
	radb = radb/compte;
	D_pixel the_pix(zb, yb, xb, radb);
	
	double longeur_interm = 0;
	
	// Calcul distance par rapport au barycentre
	for(unsigned int i = 0; i < list_diam.size(); i++)
	{
		longeur_interm = sqrt(pow(the_pix.getCoordZ()-l_node[list_jdiam[i]].getValuePixel(0).getCoordZ() ,2) + pow(the_pix.getCoordY()-l_node[list_jdiam[i]].getValuePixel(0).getCoordY() ,2) + pow(the_pix.getCoordX()-l_node[list_jdiam[i]].getValuePixel(0).getCoordX() ,2));
		/*cout << longeur_interm << endl;*/
		list_longueur.push_back(longeur_interm);
	}
	int mon_max =0 ;
	int mon_it = 0;
	
	for (unsigned int i = 0; i < list_longueur.size(); i++)
	{
		if (list_longueur[i]>mon_max)
		{
			mon_max = list_longueur[i];
			mon_it = i; 
		}
	}
	int mon_itit = 0 ;
	for(unsigned int s=0; s < list_diam.size();s++)
	{
		for(int k = 0; k < getNumberEdge(); k++)
		{
			if (getCoupleNode(k).getValueNode1() == list_jdiam[mon_it])
			{	
				if(l_node[getCoupleNode(k).getValueNode2()].getTypeNode() == "terminal")
				{
					mon_itit = getCoupleNode(k).getValueNode2();
				}
			}
		}
	}
	
	// Points terminaux les plus grands
	getListeTerm();
	std::vector<int> list_points;
	int pointerm_max;
	
	pointerm_max = l_node[0].getValuePixel(0).getRadius();
	
	
	for(unsigned int j = 0; j < l_node_term.size(); j++)
	{
		if (l_node[l_node_term[j]].getValuePixel(l_node[l_node_term[j]].getSizeNode()-1).getRadius() > pointerm_max)
		{
			pointerm_max = l_node[l_node_term[j]].getValuePixel(l_node[l_node_term[j]].getSizeNode()-1).getRadius();
		}
	}
	if (l_node[0].getValuePixel(0).getRadius() == pointerm_max)
	{
		list_points.push_back(0);
	}
	
	for(unsigned int j = 0; j < l_node_term.size(); j++)
	{
		if (l_node[l_node_term[j]].getValuePixel(l_node[l_node_term[j]].getSizeNode()-1).getRadius() == pointerm_max)
		{
			list_points.push_back(l_node_term[j]);
		}
	}
	
	//// Fusion liste
	list_merge.resize(list_angle.size() + list_propag.size());
	/*cout << list_angle.size() << endl ;
	cout << list_propag.size() << endl ;
	cout << list_diam.size() << endl ;
	cout << list_merge.size() << endl;*/
	merge(list_angle.begin(),list_angle.end(),list_propag.begin(),list_propag.end(),list_merge.begin());
	
	list_merge2.resize(list_merge.size() + list_points.size());
	merge(list_merge.begin(),list_merge.end(),list_points.begin(),list_points.end(),list_merge2.begin());
	
	list_final.resize(list_merge2.size() + list_diam.size() );
	//list_final.resize(list_angle.size() + list_diam.size() );
	//cout << list_final.size() << endl;
	merge(list_merge2.begin(),list_merge2.end(),list_diam.begin(),list_diam.end(),list_final.begin());
	
	tag.clear();
	tag.resize(getNumberNode());	
	
	for(int i = 0; i < getNumberNode(); i++)
	{ 
		mycount = count(list_final.begin(), list_final.end(), i);
		//cout << " Node : " << i << " compte : " << mycount << endl;
		if (i == mon_itit)
		{
			mycount = mycount+1;
		}
		if (mycount>=3 /*&& i == mon_itit*/)
		{
			tag[i] = 1 ;
		}
		//else if(mycount>=2)
		//{
		//	tag[i] = 1 ;
		//}
	}
	
	int mon_pix;
	int z;
	int y;
	int x;
	int ma_br =0	;
	if(ma_br == 0)
	{
		mon_pix = 0;
	}
	else 
	{
		mon_pix = l_node[ma_br].getSizeNode()-1;
	}
	z = l_node[ma_br].getValuePixel(mon_pix).getCoordZ();
	y = l_node[ma_br].getValuePixel(mon_pix).getCoordY();
	x = l_node[ma_br].getValuePixel(mon_pix).getCoordX();

	
	D_coord_ori Mescoord_ori(z,y,x,sk_array[z][y][x]);
	return(Mescoord_ori);
}


/*******************************************************************************************/
void  D_graphe::visualization2DGrapheLabel()
{
	Liste_Adjacence();
	vector <int> ma_jct;
	vector <vector <int> > iterator;
	
	vtkSmartPointer<vtkMutableUndirectedGraph> g = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
  
  // Tableau pour les labels
  vtkSmartPointer<vtkIntArray> vertexIDs = vtkSmartPointer<vtkIntArray>::New();
  vertexIDs->SetNumberOfComponents(1);
  vertexIDs->SetName("VertexIDs");

  // Tableau pour les coordonnees de points
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  // Tableau pour les couleurs
  vtkSmartPointer<vtkIntArray> edgeColors = vtkSmartPointer<vtkIntArray>::New();
  edgeColors->SetNumberOfComponents(1);
  edgeColors->SetName("Color");
 
  vtkSmartPointer<vtkLookupTable> lookupTable = vtkSmartPointer<vtkLookupTable>::New();
  lookupTable->SetNumberOfTableValues(2);
  lookupTable->SetTableValue(0, 1.0, 0.0, 0.0); // red
  lookupTable->SetTableValue(1, 0.0, 0.0, 1.0); // blue
  lookupTable->Build();
	
	
	// Premier noeud
	vtkIdType v0 = g->AddVertex();
	vertexIDs->InsertNextValue(0);
	points->InsertNextPoint(l_node[0].getValuePixel(0).getCoordX(), l_node[0].getValuePixel(0).getCoordY(), l_node[0].getValuePixel(0).getCoordZ());
	
	// Premiere jonction
	vtkIdType v1 = g->AddVertex();
	vertexIDs->InsertNextValue(1);
	points->InsertNextPoint(l_node[1].getValuePixel(0).getCoordX(), l_node[1].getValuePixel(0).getCoordY(), l_node[1].getValuePixel(0).getCoordZ());
	
	// Premier arc
	g->AddEdge(v0, v1);
	edgeColors->InsertNextValue(1);
	
 
 	// Numero de jonction
	ma_jct.push_back(1);
	// Numero ID dans graphe
	ma_jct.push_back(v1);
	// Sauvegarde
	iterator.push_back(ma_jct);
	ma_jct.clear();
        int laj = -1;
	

	for(unsigned int i = 0; i < l_Ad.size(); i++)
	{
		//cout << "*********"<< endl;
		//cout << "JONCTION " << l_Ad[i][0] << endl;
		
		for (unsigned int f = 2; f < l_Ad[i].size(); f++)
		{
			if (l_node[l_Ad[i][f]].getTypeNode() == "terminal")
			{
				//cout << "Segment Terminal " << l_Ad[i][f] << endl;
				vtkIdType v2 = g->AddVertex();
				//cout << " v2 " << v2 << endl;
				vertexIDs->InsertNextValue(l_Ad[i][f]);				
				points->InsertNextPoint(l_node[l_Ad[i][f]].getValuePixel(l_node[l_Ad[i][f]].getSizeNode()-1).getCoordX(), l_node[l_Ad[i][f]].getValuePixel(l_node[l_Ad[i][f]].getSizeNode()-1).getCoordY(), l_node[l_Ad[i][f]].getValuePixel(l_node[l_Ad[i][f]].getSizeNode()-1).getCoordZ());

				for(unsigned int p = 0; p < iterator.size(); p++)
				{
					if(iterator[p][0] == l_Ad[i][0])
					{
						laj = iterator[p][1];
						//cout << iterator[p][0] << " Last jonction pdt : " << laj << endl;
					}
				}
				
				g->AddEdge(laj, v2);
				edgeColors->InsertNextValue(1);
			}
			if (l_node[l_Ad[i][f]].getTypeNode() == "segment")
			{
				//cout << "Segment " << l_Ad[i][f] << endl;
				for(unsigned int m = 0; m < l_Ad.size(); m++)
				{
					if(l_Ad[i][f] == l_Ad[m][1])
					{ 
						vtkIdType v2 = g->AddVertex();
						//cout << "Last v2 : " << v2 << endl;
						vertexIDs->InsertNextValue(l_Ad[m][0]);
						points->InsertNextPoint(l_node[l_Ad[m][0]].getValuePixel(0).getCoordX(), l_node[l_Ad[m][0]].getValuePixel(0).getCoordY(), l_node[l_Ad[m][0]].getValuePixel(0).getCoordZ());
						for(unsigned int p = 0; p < iterator.size(); p++)
						{
							if(iterator[p][0] == l_Ad[i][0])
							{
								laj = iterator[p][1];
								//cout << iterator[p][0] << " Last jonction pdt : " << laj << endl;
							}
						}
						g->AddEdge(laj, v2);
						//cout << l_Ad[m][0] << " " << v2<< endl;
						ma_jct.push_back(l_Ad[m][0]);
						ma_jct.push_back(v2);
						iterator.push_back(ma_jct);
						ma_jct.clear();
						edgeColors->InsertNextValue(0);
					}
				}
			}
		}
	}

  // Ajouter les informations dans le graphe
  g->GetEdgeData()->AddArray(edgeColors);
  g->GetVertexData()->AddArray(vertexIDs);
  g->SetPoints(points);
  
  // Visualisation du graphe
  vtkSmartPointer<vtkGraphLayoutView> graphLayoutView = vtkSmartPointer<vtkGraphLayoutView>::New();
	graphLayoutView->SetLayoutStrategyToSimple2D();
	//graphLayoutView->SetLayoutStrategyToForceDirected();
	//graphLayoutView->SetLayoutStrategyToClustering2D();
	//graphLayoutView->SetLayoutStrategyToFast2D();
	//graphLayoutView->SetLayoutStrategyToCircular();
  graphLayoutView->AddRepresentationFromInput(g);
	graphLayoutView->SetEdgeColorArrayName("Color");
  graphLayoutView->ColorEdgesOn();
  graphLayoutView->SetVertexLabelVisibility(true);
  graphLayoutView->SetVertexLabelArrayName("VertexIDs"); //default is "labels"

  graphLayoutView->ResetCamera();
  graphLayoutView->Render();
  graphLayoutView->GetInteractor()->Start();
   
}

void  D_graphe::visualization3DLabel()
{
	Liste_Adjacence();
	
	vector <int> node;
	
	// Renderer et Windows
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  
 	renderer->SetBackground(0, 0, 0);

	// Couleur segments
	unsigned char red[3] = {255, 0, 0};
  unsigned char blue[3] = {0, 0, 255};
  vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName("Colors");	
	
	// Structure de sauvegarde
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	
	// Premier noeud
	double origin[3] = {l_node[0].getValuePixel(0).getCoordX(), l_node[0].getValuePixel(0).getCoordY(), l_node[0].getValuePixel(0).getCoordZ()};
	
	// Premiere jonction
	double nj[3] = {l_node[1].getValuePixel(0).getCoordX(), l_node[1].getValuePixel(0).getCoordY(), l_node[1].getValuePixel(0).getCoordZ()};
	
	// Sauvegarder les points dans vtkPoints
	points->InsertNextPoint(origin);
	node.push_back(0);
	points->InsertNextPoint(nj);
	node.push_back(1);
	
	// Creer la premiere ligne
	vtkSmartPointer<vtkLine> aline = vtkSmartPointer<vtkLine>::New();
	aline->GetPointIds()->SetId(0,0);
	aline->GetPointIds()->SetId(1,1);	
	colors->InsertNextTupleValue(red);
	lines->InsertNextCell(aline);

	// Labelissation du noeud 0
	stringstream sstr;
	sstr << node[0];
	string str1 = sstr.str();
	char * cstr = new char [str1.length()+1];
  strcpy (cstr, str1.c_str());
  char * p = strtok(cstr," ");

  vtkSmartPointer<vtkVectorText> textnode = vtkSmartPointer<vtkVectorText>::New();
	textnode->SetText(p);
	textnode->Update();
	
	vtkSmartPointer<vtkPolyDataMapper> mappertext = vtkSmartPointer<vtkPolyDataMapper>::New();
	mappertext->SetInputConnection(textnode->GetOutputPort());
	
	vtkSmartPointer<vtkFollower> actortext = vtkSmartPointer<vtkFollower>::New();
	actortext->SetMapper(mappertext);
	actortext->SetScale(5, 5, 5);
	actortext->AddPosition(origin);
	
  renderer->AddActor(actortext);
	renderer->ResetCameraClippingRange();
	actortext->SetCamera(renderer->GetActiveCamera());	

        int laj = -1;
	
	for(unsigned int i = 0; i < l_Ad.size(); i++)
	{
		//cout << "*********"<< endl;
		//cout << "JONCTION " << l_Ad[i][0] << endl;
		
		for (unsigned int f = 2; f < l_Ad[i].size(); f++)
		{
			if (l_node[l_Ad[i][f]].getTypeNode() == "terminal")
			{
				//cout << "Segment Terminal " << l_Ad[i][f] << endl;			
				double nj[3] = {l_node[l_Ad[i][f]].getValuePixel(l_node[l_Ad[i][f]].getSizeNode()-1).getCoordX(), l_node[l_Ad[i][f]].getValuePixel(l_node[l_Ad[i][f]].getSizeNode()-1).getCoordY(), l_node[l_Ad[i][f]].getValuePixel(l_node[l_Ad[i][f]].getSizeNode()-1).getCoordZ()};
				points->InsertNextPoint(nj);
				node.push_back(l_Ad[i][f]);

				for(unsigned int p = 0; p < node.size(); p++)
				{
					if(node[p] == l_Ad[i][0])
					{
						laj = p;
						//cout << "id jct " << node[p] << " : " << laj << endl;
					}
				}
				vtkSmartPointer<vtkLine> aline = vtkSmartPointer<vtkLine>::New();
				aline->GetPointIds()->SetId(0,laj);
				aline->GetPointIds()->SetId(1,node.size()-1);
				colors->InsertNextTupleValue(red);
				//cout << "arc " << laj << " -> " << node.size()-1 << endl;
				//cout << "soit " << node[laj] << " -> " << node[node.size()-1] << endl;
				lines->InsertNextCell(aline);
				
				// Labelissation de la branche terminale
				stringstream sstr;
				sstr << node[node.size()-1];
				string str1 = sstr.str();
				char * cstr = new char [str1.length()+1];
 				strcpy (cstr, str1.c_str());
 			 	char * p = strtok(cstr," ");
	
  			vtkSmartPointer<vtkVectorText> textnode = vtkSmartPointer<vtkVectorText>::New();
				textnode->SetText(p);
				textnode->Update();
	
				vtkSmartPointer<vtkPolyDataMapper> mappertext = vtkSmartPointer<vtkPolyDataMapper>::New();
				mappertext->SetInputConnection(textnode->GetOutputPort());
	
				vtkSmartPointer<vtkFollower> actortext = vtkSmartPointer<vtkFollower>::New();
				actortext->SetMapper(mappertext);
				actortext->SetScale(5, 5, 5);
				actortext->AddPosition(nj);
	
			  renderer->AddActor(actortext);
			  renderer->ResetCameraClippingRange();
				actortext->SetCamera(renderer->GetActiveCamera());	
			}
			if (l_node[l_Ad[i][f]].getTypeNode() == "segment")
			{
				//cout << "Segment " << l_Ad[i][f] << endl;
				for(unsigned int m = 0; m < l_Ad.size(); m++)
				{
					if(l_Ad[i][f] == l_Ad[m][1])
					{ 
						double nj[3] = {l_node[l_Ad[m][0]].getValuePixel(0).getCoordX(), l_node[l_Ad[m][0]].getValuePixel(0).getCoordY(), l_node[l_Ad[m][0]].getValuePixel(0).getCoordZ()};
						points->InsertNextPoint(nj);
						node.push_back(l_Ad[m][0]);
						
						for(unsigned int p = 0; p < node.size(); p++)
						{
							if(node[p] == l_Ad[i][0])
							{
								laj = p;
								//cout << "id jct " << node[p] << " : " << laj << endl;
							}
						}
						vtkSmartPointer<vtkLine> aline = vtkSmartPointer<vtkLine>::New();
						aline->GetPointIds()->SetId(0,laj);
						aline->GetPointIds()->SetId(1,node.size()-1);
						colors->InsertNextTupleValue(blue);
						//cout << "arc " << laj << " -> " << node.size()-1 << endl;
						//cout << "soit " << node[laj] << " -> " << node[node.size()-1] << endl;
						lines->InsertNextCell(aline);
					}
				}
			}
		}
	}
 	
	// Creation d'un polydata pour tout sauvegarder dedans
	vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();
	linesPolyData->SetPoints(points);
	linesPolyData->SetLines(lines);
	linesPolyData->GetCellData()->SetScalars(colors);

	// Points a chaque point
	vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
	vertexGlyphFilter->AddInputData(linesPolyData);
	vertexGlyphFilter->Update();

	// Mapper
	vtkSmartPointer<vtkPolyDataMapper> mapperligne = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapperligne->SetInputData(linesPolyData);

	vtkSmartPointer<vtkPolyDataMapper> mapperglyph = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapperglyph->SetInputConnection(vertexGlyphFilter->GetOutputPort());

  
  // Actor
	vtkSmartPointer<vtkActor> actorligne = vtkSmartPointer<vtkActor>::New();
	actorligne->SetMapper(mapperligne);
	actorligne->GetProperty()->SetLineWidth(4);

	vtkSmartPointer<vtkActor> actorglyph = vtkSmartPointer<vtkActor>::New();
	actorglyph->SetMapper(mapperglyph);
	actorglyph->GetProperty()->SetPointSize(4);
  actorglyph->GetProperty()->SetColor(1,0,0);
  

	renderer->AddActor(actorligne);
	renderer->AddActor(actorglyph);
	
	renderer->ResetCamera();	
	
	renderWindow->Render();
	renderWindowInteractor->Start();

}

void  D_graphe::visualizationTag()
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
		if (tag[i] == 0) //Si segment couleur selon ordre
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
 				
 				sphereActor->GetProperty()->SetColor(0,0,0.5);
				sphereActor->SetMapper(sphereMapper);
				ren1->AddActor(sphereActor);
			}
		}
		if (tag[i] == 1) //Si segment couleur selon ordre
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
 				
 				sphereActor->GetProperty()->SetColor(1,0.8,0);
				sphereActor->SetMapper(sphereMapper);
				ren1->AddActor(sphereActor);
			}
		}
	}
	// This starts the event loop and invokes an initial render.
  iren->Initialize();
  iren->Start();
}
/*******************************************************************************************/
// Fonctions enregistrement du graphe

void  D_graphe::saveGraphe(char* name)
{
	ofstream monFlux(name);
	
	if (monFlux)
	{
		monFlux << "GRAPHE" << endl;
		monFlux << "Nombre_Arc " <<	getNumberEdge() << endl;
		monFlux << "Nombre_Noeud " <<	getNumberNode() << endl;
		monFlux << "Nombre_Jonction " <<	getNumberJunction() << endl;
		monFlux << "Nombre_Segment " <<	getNumberSegment() << endl;
		monFlux << "Nombre_Segment_Terminal " <<	getNumberSegTerminal() << endl;
		monFlux << "Nombre_Total_de_Pixel_de_Jonction_du_graphe " << getNumberTotPixelJunction() << endl;
		monFlux	<< "Nombre_Total_de_Pixel_de_Segment_du_graphe " << getNumberTotPixelSegment() << endl;
		monFlux << "Hauteur_Maximale_Arbre " << getHauteurMaxA() << endl;
		monFlux << endl;
		
		monFlux << "NOEUDS" << endl;
		for(int i = 0; i < getNumberNode(); i++)
		{
			monFlux << "Node " << i << endl;
			monFlux << "Type " << l_node[i].getTypeNode() << endl;
			monFlux << "Size " << l_node[i].getSizeNode() << endl;
			if (l_node[i].getTypeNode() == "junction")
			{
				monFlux << "Hauteur " << l_node[i].getHauteurNoeud() << endl;
			}
			for(int m = 0; m < l_node[i].getSizeNode(); m++)
			{
				monFlux <<"Pix " << l_node[i].getValuePixel(m).getCoordZ() << " " << l_node[i].getValuePixel(m).getCoordY()  << " " <<l_node[i].getValuePixel(m).getCoordX() << " " << l_node[i].getValuePixel(m).getRadius() << endl;
			}
			monFlux << "End " << endl;
			monFlux <<endl;
		}
		monFlux << "ARCS" << endl;
		for(int i = 0; i < getNumberEdge(); i++)
		{
			monFlux << i << " " << l_edgeN[i].getValueNode1() << " " << l_edgeN[i].getValueNode2() << " " << l_edgeN[i].getWeight() << endl;
		}
		monFlux << endl;
	}
	else
	{
		cout << "ERREUR : Impossible d'ouvrir le fichier." << endl;
	}
}

void D_graphe::openGraphe(char* name)
{
	ifstream monFlux(name, ios::in);
	
	if (monFlux)
	{
		string qlq ;
		int nqlq; 
		int Nombre_Arc = 0 ;
		int Nombre_Noeud = 0 ;
		int Nombre_Jonction = 0 ;
		int Nombre_Segment = 0 ;	
		int Nombre_Segment_Terminal = 0 ;	
		int Nbr_Tot_Pix_Jct = 0 ;	
		int Nbr_Tot_Pix_Seg = 0 ;		
		vector<int> size(0);
		int j = 0;
		int i = 0;
		double x, y;
		int c1, c2, wght;
		int radius;
		string ligne;
		istringstream ss;
    
    while (getline(monFlux, ligne)) 
		{ 
    	if (ligne == "GRAPHE")
   		{
   			while (getline(monFlux, ligne) && j < 9) 
   			{
   				ss.str(ligne);
   				ss >> qlq >> nqlq;
   				ss.clear();
   				if (qlq == "Nombre_Arc")
   				{
   					Nombre_Arc = nqlq;
   					l_edgeN.resize(Nombre_Arc);
   					j+=1;
   				}
   				if (qlq == "Nombre_Noeud")
   				{
   					Nombre_Noeud = nqlq;
   					l_node.resize(Nombre_Noeud);
   					j+=1;
   				}
   				if (qlq == "Nombre_Jonction")
   				{
   					Nombre_Jonction = nqlq;
   					j+=1;
   				}
   				if (qlq == "Nombre_Segment")
   				{
   					Nombre_Segment = nqlq;
   					j+=1;
   				}
   				if (qlq == "Nombre_Segment_Terminal")
   				{
   					Nombre_Segment_Terminal = nqlq;
   					j+=1;
   				}
   				if (qlq == "Nombre_Total_de_Pixel_de_Jonction_du_graphe")
   				{
   					Nbr_Tot_Pix_Jct = nqlq;
   					j+=1;
   				}
   				if (qlq == "Nombre_Total_de_Pixel_de_Segment_du_graphe")
   				{
   					Nbr_Tot_Pix_Seg = nqlq;
   					j+=1;
   				}
   				if (qlq == "Hauteur_Maximale_Arbre")
   				{
   					setHauteurMaxA(nqlq);
   					j+=1;
   				}
				}
			}				
			
			if (ligne == "NOEUDS")
   		{	
   			while (getline(monFlux, ligne) && ligne != "ARCS") 
				{ 
					ss.str(ligne);
					ss >> qlq >> nqlq;
					ss.clear();
					if (qlq == "Node")
					{
						i = nqlq;
					}
					if (qlq == "Type")
					{
						ss >> qlq;
						ss.clear();
						l_node[i].setTypeNode(qlq);
					}
					if (qlq == "Size")
					{
						ss >> nqlq;
						ss.clear();
						size.push_back(nqlq);
					}
					if (qlq == "Hauteur")
					{
						ss >> nqlq;
						ss.clear();
						l_node[i].setHauteurNoeud(nqlq);
					}	
					if (qlq == "Pix")
					{
						ss >> y >> x >> radius;
						ss.clear();
						D_pixel the_pix(nqlq, y, x, radius);
						l_node[i].InsertNextPixel(the_pix);
					}
				}
   		}
			if (ligne == "ARCS")
   		{	
   			while (getline(monFlux, ligne)) 
   			{
   				ss.str(ligne);
					ss >> nqlq >> c1 >> c2 >> wght;
					ss.clear();
					l_edgeN[nqlq].InsertCoupleNode(c1,c2,wght);	
   			}
   		}	
   	}
   	
   	monFlux.close();
   	
   	// Verification de la reconstruction du graphe
   	if (Nombre_Jonction != getNumberJunction() || Nombre_Segment != getNumberSegment() || Nombre_Segment_Terminal != getNumberSegTerminal() || Nbr_Tot_Pix_Jct != getNumberTotPixelJunction() || Nbr_Tot_Pix_Seg != getNumberTotPixelSegment())
		{
			cout << "ERREUR : Reconstruction du graphe incorrecte" << endl;
		}	
		for (unsigned int i = 0; i < size.size() ; i++)
		{
			if (size[i] != l_node[i].getSizeNode())
			{
				cout << "ERREUR : Reconstruction du graphe incorrecte" << endl;
			}
		}	
	}
	else
	{
		cout << "ERREUR : Impossible d'ouvrir le fichier." << endl;
	}
}




























