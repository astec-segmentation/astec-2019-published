
#include <vt_geline2D.h>

static int  VT_ComputeOneLine2D( vt_images *ims, vt_image *theim, vt_recfilters *par, vt_line *lpar );
static int  VT_FilterImages( vt_images *ims, vt_image *im, vt_recfilters *par );





int VT_ComputeLine2D( vt_resline *res, vt_images *ims, vt_image *theim, vt_line *par )
{
    char *proc="VT_ComputeLine2D";
    register int x, y, z;
    int dimx, dimy, dimz, i;
    double coeff, step;
    vt_recfilters local_par;
    r32 ***bufr, ***bufx, ***bufy;
    u8 ***bufe, ***bufs;
    r32 ***rbuf, ***xbuf, ***ybuf;
    u8 ***ebuf;
    float r;

    /*--- tests sur resultats ---*/
    if ( VT_Test2Image( &(res->imres),   theim, proc ) == -1 ) return( -1 );
    if ( VT_Test2Image( &(res->imext),   theim, proc ) == -1 ) return( -1 );
    if ( VT_Test2Image( &(res->imscale), theim, proc ) == -1 ) return( -1 );
    if ( VT_Test2Image( &(res->imdirx),  theim, proc ) == -1 ) return( -1 );
    if ( VT_Test2Image( &(res->imdiry),  theim, proc ) == -1 ) return( -1 );
    if ( (res->imres.type != FLOAT) || (res->imext.type != UCHAR) || (res->imscale.type != UCHAR) ||
	 (res->imdirx.type != FLOAT) || (res->imdiry.type != FLOAT) ) {
	VT_Error( "bad type for result images", proc );
	return( -1 );
    }
    /*--- tests sur images (a faire) ---*/
    /*--- tests sur parametres ---*/
    if ( (par->first_coeff <= 0.0) || (par->nb_coeff <= 0) || ((par->nb_coeff > 1) && (par->last_coeff <= par->first_coeff)) ) {
	VT_Error( "bad coefficients for multiscale analysis", proc );
	return( -1 );
    }
    if ( par->nb_coeff > 255 ) {
	VT_Error( "too many scales", proc );
	return( -1 );
    }

    /*--- initialisation ---*/
    dimx = theim->dim.x;
    dimy = theim->dim.y;
    dimz = theim->dim.z;
    
    /*--- buffers ---*/
    bufr = (r32 ***)(res->imres.array);
    bufe = (u8 ***)(res->imext.array);
    bufs = (u8 ***)(res->imscale.array);
    bufx = (r32 ***)(res->imdirx.array);
    bufy = (r32 ***)(res->imdiry.array);
    rbuf = (r32 ***)(ims->imr.array);
    ebuf = (u8 ***)(ims->ime.array);
    xbuf = (r32 ***)(ims->imxx.array);
    ybuf = (r32 ***)(ims->imyy.array);

    /*--- calcul pour les echelles ---*/
    step = 0.0;
    if ( par->nb_coeff > 1 ) 
	step = (double)(par->last_coeff - par->first_coeff) / (double)(par->nb_coeff - 1);
    else 
      par->last_coeff = par->first_coeff;
    coeff = par->last_coeff;
    local_par = par->par_filt;
    for (i = 0; i < par->nb_coeff; i++, coeff -= step ) {
	local_par.value_coefficient.x = (float)coeff;
	local_par.value_coefficient.y = (float)coeff;
	local_par.value_coefficient.z = (float)coeff;
	/*--- calcul pour une echelle ---*/
	if ( _VT_VERBOSE_ == 1 ) {
	    char message[256];
	    sprintf( message,"filtrage avec coefficient = %f", (float)coeff );
	    VT_Message( message, proc );
	}
	if ( VT_ComputeOneLine2D( ims, theim, &local_par, par ) != 1 ) {
	    VT_Error( "unable to compute one scale", proc );
	    return( -1 );
	}

	/*--- resultats :
          solution 1 : collecter les extrema a chaque echelle dans l'image
	               extrema resultat. Ce n'est pas satisfaisant, car les
		       extrema se delocalisent a chaque echelle et cela 
		       finit par faire des lignes epaisses.
	  solution 2 : pour un point, prendre la reponse maximale selon
	               les echelles et regarder si cette reponse etait
		       un extrema a cette echelle.                         ---*/
	
	if ( i == 0 ) {
	    for ( z = 0; z < dimz; z ++ )
	    for ( y = 0; y < dimy; y ++ )
	    for ( x = 0; x < dimx; x ++ ) {
		bufr[z][y][x] = rbuf[z][y][x];
		bufe[z][y][x] = ebuf[z][y][x];
		bufs[z][y][x] = (u8)(par->nb_coeff);
		bufx[z][y][x] = xbuf[z][y][x];
		bufy[z][y][x] = ybuf[z][y][x];
	    }
	} else {
	    for ( z = 0; z < dimz; z ++ )
	    for ( y = 0; y < dimy; y ++ )
	    for ( x = 0; x < dimx; x ++ ) {
		r = rbuf[z][y][x];
		if ( r > bufr[z][y][x] ) {
		    bufr[z][y][x] = r;
		    bufe[z][y][x] = ebuf[z][y][x];
		    bufs[z][y][x] = (u8)(par->nb_coeff - i);
		    bufx[z][y][x] = xbuf[z][y][x];
		    bufy[z][y][x] = ybuf[z][y][x];
		}
	    }
	}
    }

    /*--- fin ---*/
    return( 1 );
}





static int VT_ComputeOneLine2D( vt_images *ims, vt_image *theim, vt_recfilters *par, vt_line *lpar )
{
    char *proc="VT_ComputeOneLine2D";
    register int x, y, z;
    int dimx, dimy, dimz, dimx1, dimy1;
    r32 ***bufx, ***bufy, ***bufxx, ***bufxy, ***bufyy;
    r32 ***bufr;
    u8 ***bufe;
    double det, sigma, aux=0.0;
    double l, v[2];
    int xi, yi;
    double xr, yr, dx, dy, coeff[2][2];
    double gx, gy, r1, r2;

    /*--- test ---*/
    if ( VT_Test2Image( &(ims->imr), theim, proc ) == -1 ) return( -1 );
    if ( VT_Test2Image( &(ims->ime), theim, proc ) == -1 ) return( -1 );
    if ( (ims->imr.type != FLOAT) || (ims->ime.type != UCHAR) ) {
	VT_Error( "bad type for result images", proc );
	return( -1 );
    }

    /*--- initialisation ---*/
    dimx = theim->dim.x;   dimx1 = dimx - 1;
    dimy = theim->dim.y;   dimy1 = dimy - 1;
    dimz = theim->dim.z;    
    sigma = par->value_coefficient.x;

    /*--- filtrage ---*/
    if ( VT_FilterImages( ims, theim, par ) != 1 ) {
	VT_Error( "unable to compute derivatives", proc );
	return( -1 );
    }

    /*--- buffers ---*/
    bufx  = (r32 ***)(ims->imx.array);
    bufy  = (r32 ***)(ims->imy.array);
    bufxx = (r32 ***)(ims->imxx.array);
    bufxy = (r32 ***)(ims->imxy.array);
    bufyy = (r32 ***)(ims->imyy.array);
    bufr  = (r32 ***)(ims->imr.array);
    bufe  = (u8 ***)(ims->ime.array);

    /*--- calcul de la reponse : on remplit bufr i.e. ims->imr ---*/
    for ( z = 0; z < dimz; z ++ )
    for ( y = 0; y < dimy; y ++ )
    for ( x = 0; x < dimx; x ++ ) {

	/*--- initialisation ---*/
	bufr[z][y][x] = 0.0;

	/*--- calcul de l'orientation :
	      1. calcul du determinant du polynome caracteristique du hessien
	      2. calcul de la plus grande valeur propre (en valeur absolue)
	         c'est la valeur propre stable qui donne la direction
		 orthogonale du vaisseau.
		 si cette valeur propre est positive, c'est un vaisseau
		 noir sur fond blanc, sinon c'est l'inverse.
	      3. calcul et normalisation du vecteur propre associe             ---*/
	aux  = bufxx[z][y][x] + bufyy[z][y][x];
	det  = aux * aux - 4.0 * ( bufxx[z][y][x] * bufyy[z][y][x] - bufxy[z][y][x] * bufxy[z][y][x] );
	if ( det < 0.0 ) continue;
	det = VT_Sqrt( det );
	/*--- calcul de la valeur propre et du vecteur propre associe ---*/
	if ( aux < 0.0 ) l = ( aux - det ) / 2.0;
	else             l = ( aux + det ) / 2.0;
	/*--- est-ce une structure a conserver ? ---*/
	if ( (l > 0.0) && (lpar->type_structures == VT_WHITE) ) continue;
	if ( (l < 0.0) && (lpar->type_structures == VT_BLACK) ) continue;
	/*--- calcul du vecteur propre associe ---*/
 	v[0] = bufxy[z][y][x];    v[1] = l - bufxx[z][y][x];
 	det = VT_Sqrt( v[0]*v[0] + v[1]*v[1] );
	if ( det < 0.0000001 ) continue;
 	v[0] /= det;              v[1] /= det;
	/*--- on range un vecteur propre dans (bufxx, bufyy) ---*/
	bufxx[z][y][x] = v[0];
	bufyy[z][y][x] = v[1];

	/*--- calcul de la reponse du filtre gradient selon v en deux points :
	      M_{1,2} = M +/- sigma * v
	      on interpole les derivees selon X et Y en ces points et on
	      fait le produit scalaire avec v.                              ---*/
	xr = (double)(x) + sigma * v[0];   yr = (double)(y) + sigma * v[1];
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) ) continue;
	xi = (int)xr;                       yi = (int)yr;
	dx = xr - (double)xi;               dy = yr - (double)yi;
	/*--- interpolation de la 1ere derivee ---*/
	coeff[0][0] = (1.0 - dx) * (1.0 - dy); /* xi,   yi   */
	coeff[0][1] = dx         * (1.0 - dy); /* xi+1, yi   */
	coeff[1][0] = (1.0 - dx) * dy;         /* xi,   yi+1 */
	coeff[1][1] = dx         * dy;         /* xi+1, yi+1 */
	gx  = coeff[0][0] * bufx[z][yi][xi]   + coeff[0][1] * bufx[z][yi][xi+1];
	gx += coeff[1][0] * bufx[z][yi+1][xi] + coeff[1][1] * bufx[z][yi+1][xi+1];
	gy  = coeff[0][0] * bufy[z][yi][xi]   + coeff[0][1] * bufy[z][yi][xi+1];
	gy += coeff[1][0] * bufy[z][yi+1][xi] + coeff[1][1] * bufy[z][yi+1][xi+1];
	r1 = gx * v[0] + gy * v[1];
	/*--- 2nde derivee ---*/
	xr = (double)(x) - sigma * v[0];   yr = (double)(y) - sigma * v[1];
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) ) continue;
	xi = (int)xr;                       yi = (int)yr;
	dx = xr - (double)xi;               dy = yr - (double)yi;
	/*--- interpolation de la 2nde derivee : on recalcule les
	      coefficients pour des raison numeriques             ---*/
	coeff[0][0] = (1.0 - dx) * (1.0 - dy); /* xi,   yi   */
	coeff[0][1] = dx         * (1.0 - dy); /* xi+1, yi   */
	coeff[1][0] = (1.0 - dx) * dy;         /* xi,   yi+1 */
	coeff[1][1] = dx         * dy;         /* xi+1, yi+1 */
	gx  = coeff[0][0] * bufx[z][yi][xi]   + coeff[0][1] * bufx[z][yi][xi+1];
	gx += coeff[1][0] * bufx[z][yi+1][xi] + coeff[1][1] * bufx[z][yi+1][xi+1];
	gy  = coeff[0][0] * bufy[z][yi][xi]   + coeff[0][1] * bufy[z][yi][xi+1];
	gy += coeff[1][0] * bufy[z][yi+1][xi] + coeff[1][1] * bufy[z][yi+1][xi+1];
	r2 = gx * v[0] + gy * v[1];
	/*--- on verifie que les deux reponses ont bien des signes oppposes ---*/
	if ( r1*r2 >= 0.0 ) continue;
	/*--- on change les signes pour comparer ---*/
	if ( r1 < 0.0 ) r1 = (- r1);
	if ( r2 < 0.0 ) r2 = (- r2);
	if ( r2 > r1 ) bufr[z][y][x] = r1;
	else           bufr[z][y][x] = r2;

    }

    /*--- calcul de l'extremalite : on remplit bufr i.e. ims->ime.
          on change egalement les buffers bufx et bufy
	  qui representaient les derivees selon x et y de
	  l'image : maintenant ce sont les coordonnees du vecteur
	  propre associee a la valeur propre.                      ---*/
    for ( z = 0; z < dimz; z ++ )
    for ( y = 0; y < dimy; y ++ )
    for ( x = 0; x < dimx; x ++ ) {

	/*--- initialisation ---*/
	bufe[z][y][x] = (u8)0;
	
	if (bufr[z][y][x] == 0.0) continue;

	/*--- calcul de la reponse du filtre non lineaire en deux points :
	      r1 et r2. Le point courant est extrema si
	      bufr[z][y][x] >= r1 et bufr[z][y][x] >= r2                   ---*/
	xr = (double)x + bufxx[z][y][x];   yr = (double)y + bufyy[z][y][x];
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) ) continue;
	xi = (int)xr;               yi = (int)yr;
	dx = xr - (double)xi;       dy = yr - (double)yi;
	/*--- interpolation de la 1ere reponse ---*/
	coeff[0][0] = (1.0 - dx) * (1.0 - dy); /* xi,   yi   */
	coeff[0][1] = dx         * (1.0 - dy); /* xi+1, yi   */
	coeff[1][0] = (1.0 - dx) * dy;         /* xi,   yi+1 */
	coeff[1][1] = dx         * dy;         /* xi+1, yi+1 */
	r1  = coeff[0][0] * bufr[z][yi][xi]   + coeff[0][1] * bufr[z][yi][xi+1];
	r1 += coeff[1][0] * bufr[z][yi+1][xi] + coeff[1][1] * bufr[z][yi+1][xi+1];
	if ( r1 > bufr[z][y][x] ) continue;
	/*--- 2nde reponse ---*/
	xr = (double)x - bufxx[z][y][x];   yr = (double)y - bufyy[z][y][x];
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) ) continue;
	xi = (int)xr;               yi = (int)yr;
	dx = xr - (double)xi;       dy = yr - (double)yi;
	/*--- interpolation de la 2nde reponse ---*/
	coeff[0][0] = (1.0 - dx) * (1.0 - dy); /* xi,   yi   */
	coeff[0][1] = dx         * (1.0 - dy); /* xi+1, yi   */
	coeff[1][0] = (1.0 - dx) * dy;         /* xi,   yi+1 */
	coeff[1][1] = dx         * dy;         /* xi+1, yi+1 */
	r2  = coeff[0][0] * bufr[z][yi][xi]   + coeff[0][1] * bufr[z][yi][xi+1];
	r2 += coeff[1][0] * bufr[z][yi+1][xi] + coeff[1][1] * bufr[z][yi+1][xi+1];
	if ( r2 > bufr[z][y][x] ) continue;

	bufe[z][y][x] = (u8)255;
	
    }

    return( 1 );

}





/*============================================================*/

static int VT_FilterImages( vt_images *ims, vt_image *im, vt_recfilters *par )
{
    vt_recfilters local_par;
    char *proc="VT_FilterImages";

    local_par = *par;
    local_par.derivative.z = VT_NODERIVATIVE;

    /*--- derivee selon X ---*/
    local_par.derivative.x = VT_DERIVATIVE_1;
    local_par.derivative.y = VT_DERIVATIVE_0;
    if ( VT_RecFilterOnImage( im, &(ims->imx), &local_par ) != 1 ) {
	VT_Error( "unable to compute dI/dx", proc );
	return( -1 );
    }

    /*--- derivee selon Y ---*/
    local_par.derivative.x = VT_DERIVATIVE_0;
    local_par.derivative.y = VT_DERIVATIVE_1;
    if ( VT_RecFilterOnImage( im, &(ims->imy), &local_par ) != 1 ) {
	VT_Error( "unable to compute dI/dy", proc );
	return( -1 );
    }

    /*--- derivee selon XX ---*/
    local_par.derivative.x = VT_DERIVATIVE_2;
    local_par.derivative.y = VT_DERIVATIVE_0;
    if ( VT_RecFilterOnImage( im, &(ims->imxx), &local_par ) != 1 ) {
	VT_Error( "unable to compute d^2I/dx^2", proc );
	return( -1 );
    }

    /*--- derivee selon XY ---*/
    local_par.derivative.x = VT_DERIVATIVE_1;
    local_par.derivative.y = VT_DERIVATIVE_1;
    if ( VT_RecFilterOnImage( im, &(ims->imxy), &local_par ) != 1 ) {
	VT_Error( "unable to compute d^2I/dxdy", proc );
	return( -1 );
    }

    /*--- derivee selon YY ---*/
    local_par.derivative.x = VT_DERIVATIVE_0;
    local_par.derivative.y = VT_DERIVATIVE_2;
    if ( VT_RecFilterOnImage( im, &(ims->imyy), &local_par ) != 1 ) {
	VT_Error( "unable to compute d^2I/dy^2", proc );
	return( -1 );
    }
    
    return( 1 );
}
