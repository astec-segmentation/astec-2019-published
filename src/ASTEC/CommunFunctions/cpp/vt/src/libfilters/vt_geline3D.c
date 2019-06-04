
#include <vt_geline3D.h>

#define _LOCAL_DEBUG_ 0

static int  VT_ComputeOneLine3D( vt_images *ims, vt_image *theim, vt_recfilters *par, vt_line *lpar );
static int  VT_FilterImages( vt_images *ims, vt_image *im, vt_recfilters *par );





int VT_ComputeLine3D( vt_resline *res, vt_images *ims, vt_image *theim, vt_line *par )
{
    char *proc="VT_ComputeLine3D";
    register int x, y, z;
    int dimx, dimy, dimz, i;
    double coeff, step;
    vt_recfilters local_par;
    r32 ***bufr, ***bufx, ***bufy, ***bufz;
    r32 ***bufx2, ***bufy2, ***bufz2;
    u8 ***bufe, ***bufs;
    r32 ***rbuf, ***xbuf, ***ybuf, ***zbuf;
    r32 ***x2buf, ***y2buf, ***z2buf;
    u8 ***ebuf;
    float r;

    /*--- tests sur resultats ---*/
    if ( VT_Test2Image( &(res->imres),   theim, proc ) == -1 ) return( -1 );
    if ( VT_Test2Image( &(res->imext),   theim, proc ) == -1 ) return( -1 );
    if ( VT_Test2Image( &(res->imscale), theim, proc ) == -1 ) return( -1 );
    if ( VT_Test2Image( &(res->imdirx),  theim, proc ) == -1 ) return( -1 );
    if ( VT_Test2Image( &(res->imdiry),  theim, proc ) == -1 ) return( -1 );
    if ( VT_Test2Image( &(res->imdirz),  theim, proc ) == -1 ) return( -1 );
    if ( VT_Test2Image( &(res->imdir2x),  theim, proc ) == -1 ) return( -1 );
    if ( VT_Test2Image( &(res->imdir2y),  theim, proc ) == -1 ) return( -1 );
    if ( VT_Test2Image( &(res->imdir2z),  theim, proc ) == -1 ) return( -1 );
    if ( (res->imres.type != FLOAT) || (res->imext.type != UCHAR) || (res->imscale.type != UCHAR) ) {
	VT_Error( "bad type for result images", proc );
	return( -1 );
    }
    if ( (res->imdirx.type != FLOAT) || (res->imdiry.type != FLOAT) || (res->imdirz.type != FLOAT) ) {
	VT_Error( "bad type for result images (2)", proc );
	return( -1 );
    }
    if ( (res->imdir2x.type != FLOAT) || (res->imdir2y.type != FLOAT) || (res->imdir2z.type != FLOAT) ) {
	VT_Error( "bad type for result images (2)", proc );
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
    bufz = (r32 ***)(res->imdirz.array);
    bufx2 = (r32 ***)(res->imdir2x.array);
    bufy2 = (r32 ***)(res->imdir2y.array);
    bufz2 = (r32 ***)(res->imdir2z.array);
    rbuf = (r32 ***)(ims->imr.array);
    ebuf = (u8 ***)(ims->ime.array);
    xbuf = (r32 ***)(ims->imxx.array);
    ybuf = (r32 ***)(ims->imyy.array);
    zbuf = (r32 ***)(ims->imzz.array);
    x2buf = (r32 ***)(ims->imyz.array);
    y2buf = (r32 ***)(ims->imxy.array);
    z2buf = (r32 ***)(ims->imxz.array);

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
	if ( VT_ComputeOneLine3D( ims, theim, &local_par, par ) != 1 ) {
	    VT_Error( "unable to compute one scale", proc );
	    return( -1 );
	}

	/*
	if ( _LOCAL_DEBUG_ == 1 ) {
	    (void)VT_WriteInrimage( &(ims->imx) );
	    (void)VT_WriteInrimage( &(ims->imy) );
	    (void)VT_WriteInrimage( &(ims->imxx) );
	    (void)VT_WriteInrimage( &(ims->imxy) );
	    (void)VT_WriteInrimage( &(ims->imyy) );
	    (void)VT_WriteInrimage( &(ims->imr) );
	    (void)VT_WriteInrimage( &(ims->ime) );
	    (void)VT_WriteInrimage( &(ims->imz) );
	    (void)VT_WriteInrimage( &(ims->imxz) );
	    (void)VT_WriteInrimage( &(ims->imyz) );
	    (void)VT_WriteInrimage( &(ims->imzz) );
	}
	*/

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
		bufz[z][y][x] = zbuf[z][y][x];
		bufx2[z][y][x] = x2buf[z][y][x];
		bufy2[z][y][x] = y2buf[z][y][x];
		bufz2[z][y][x] = z2buf[z][y][x];
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
		    bufz[z][y][x] = zbuf[z][y][x];
		    bufx2[z][y][x] = x2buf[z][y][x];
		    bufy2[z][y][x] = y2buf[z][y][x];
		    bufz2[z][y][x] = z2buf[z][y][x];
		}
	    }
	}
    }

    /*--- fin ---*/
    return( 1 );
}





static int VT_ComputeOneLine3D( vt_images *ims, vt_image *theim, vt_recfilters *par, vt_line *lpar )
{
    char *proc="VT_ComputeOneLine3D";
    register int x, y, z;
    int dimx, dimy, dimz, dimx1, dimy1, dimz1;
    r32 ***bufx, ***bufy, ***bufz, ***bufxx, ***bufxy, ***bufxz, ***bufyy, ***bufyz, ***bufzz;
    r32 ***bufr;
    u8 ***bufe;
    double sigma;
    double hessien[9], val_propre[3], vec_propre[9];
    int xi, yi, zi;
    double xr, yr, zr, dx, dy, dz, coeff[2][2][2];
    double gx, gy, gz, r[2][2], r1, r2;
    u16 ***ubuf = (u16 ***)theim->array;

    /*--- test ---*/
    if ( VT_Test2Image( &(ims->imr), theim, proc ) == -1 ) return( -1 );
    if ( VT_Test2Image( &(ims->ime), theim, proc ) == -1 ) return( -1 );
    if ( (ims->imr.type != FLOAT) || (ims->ime.type != UCHAR) ) {
	VT_Error( "bad type for result images", proc );
	return( -1 );
    }
    /*
    if ( theim->type != USHORT ) {
      VT_Error( "bad type for input image (should be unsigned short)", proc );
      return( -1 );
    }
    */

    /*--- initialisation ---*/
    dimx = theim->dim.x;   dimx1 = dimx - 1;
    dimy = theim->dim.y;   dimy1 = dimy - 1;
    dimz = theim->dim.z;   dimz1 = dimz - 1;
    sigma = par->value_coefficient.x;

    /*--- filtrage ---*/
    if ( VT_FilterImages( ims, theim, par ) != 1 ) {
	VT_Error( "unable to compute derivatives", proc );
	return( -1 );
    }
    
    /*--- buffers ---*/
    bufx  = (r32 ***)(ims->imx.array);
    bufy  = (r32 ***)(ims->imy.array);
    bufz  = (r32 ***)(ims->imz.array);
    bufxx = (r32 ***)(ims->imxx.array);
    bufxy = (r32 ***)(ims->imxy.array);
    bufxz = (r32 ***)(ims->imxz.array);
    bufyy = (r32 ***)(ims->imyy.array);
    bufyz = (r32 ***)(ims->imyz.array);
    bufzz = (r32 ***)(ims->imzz.array);
    bufr  = (r32 ***)(ims->imr.array);
    bufe  = (u8 ***)(ims->ime.array);

    /*--- calcul de la reponse : on remplit bufr i.e. ims->imr ---*/
    for ( z = 0; z < dimz; z ++ )
    for ( y = 0; y < dimy; y ++ )
    for ( x = 0; x < dimx; x ++ ) {

	/*--- initialisation ---*/
	bufr[z][y][x] = (r32)0.0;
	r[0][0] = r[0][1] = r[1][0] = r[1][1] = 0.0;

	if ( ubuf[z][y][x] == 0 ) continue;

	/*--- calcul de l'orientation :
	      1. calcul du determinant du polynome caracteristique du hessien
	      2. calcul de la plus grande valeur propre (en valeur absolue)
	         c'est la valeur propre stable qui donne la direction
		 orthogonale du vaisseau.
		 si cette valeur propre est positive, c'est un vaisseau
		 noir sur fond blanc, sinon c'est l'inverse.
	      3. calcul et normalisation du vecteur propre associe             ---*/
	hessien[0] = bufxx[z][y][x];
	hessien[1] = hessien[3] = bufxy[z][y][x];
	hessien[2] = hessien[6] = bufxz[z][y][x];
	hessien[4] = bufyy[z][y][x];
	hessien[5] = hessien[7] = bufyz[z][y][x];
	hessien[8] = bufzz[z][y][x];
	/*--- calcul des valeurs et vecteurs propres :
	      les valeurs propres sont ordonnees par ordre croissant 
	      les vecteurs propres (colonnes) sont normalises       ---*/
	
	if ( 0 ) {
	  fprintf( stderr, "fatal error (a corriger)\n" );
          exit ( 2 );
	} else {
	  if ( E_DMVVpropresMatSym( hessien, val_propre, vec_propre, (int)3 ) == 0 ) continue;
	} 
	/*--- on range un vecteur propre dans (bufxx, bufyy, bufzz)
	      et l'autre dans                 (bufyz, bufxy, bufxz) ---*/
	bufxx[z][y][x] = vec_propre[0];
	bufyy[z][y][x] = vec_propre[3];
	bufzz[z][y][x] = vec_propre[6];
	bufyz[z][y][x] = vec_propre[1];
	bufxy[z][y][x] = vec_propre[4];
	bufxz[z][y][x] = vec_propre[7];

	/*
	if ( _LOCAL_DEBUG_ == 1 )
	    fprintf( stderr, " M = (%2d,%2d,%2d), Valeurs = ( %f %f %f )\n", x,y,z,(float)val_propre[0],(float)val_propre[1],(float)val_propre[2] );
	*/

	/*--- si c'est un vaisseau blanc sur fond noir, on a :
	      les 2 plus petites VP sont negatives et grandes,
	      la troisieme est petite                          ---*/
	if ( (val_propre[0] >= 0.0) || (val_propre[1] >= 0.0) ) continue;
	if ( val_propre[1] >= (- 10.00) ) continue;

	/*--- calcul de la reponse du filtre gradient selon v en deux points :
	      M_{1,2} = M +/- sigma * v
	      on interpole les derivees selon X et Y en ces points et on
	      fait le produit scalaire avec v.                              ---*/
	xr = (double)(x) + sigma * vec_propre[0];   
	yr = (double)(y) + sigma * vec_propre[3];
	zr = (double)(z) + sigma * vec_propre[6];
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= dimz1) ) continue;
	xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	/*--- interpolation de la 1ere derivee ---*/
	coeff[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	coeff[0][1][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	coeff[1][0][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	coeff[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	coeff[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	coeff[0][1][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	coeff[1][0][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	coeff[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */

	gx  = coeff[0][0][0] * bufx[zi][yi][xi]     + coeff[0][1][0] * bufx[zi][yi][xi+1];
	gx += coeff[1][0][0] * bufx[zi][yi+1][xi]   + coeff[1][1][0] * bufx[zi][yi+1][xi+1];
	gx += coeff[0][0][1] * bufx[zi+1][yi][xi]   + coeff[0][1][1] * bufx[zi+1][yi][xi+1];
	gx += coeff[1][0][1] * bufx[zi+1][yi+1][xi] + coeff[1][1][1] * bufx[zi+1][yi+1][xi+1];

	gy  = coeff[0][0][0] * bufy[zi][yi][xi]     + coeff[0][1][0] * bufy[zi][yi][xi+1];
	gy += coeff[1][0][0] * bufy[zi][yi+1][xi]   + coeff[1][1][0] * bufy[zi][yi+1][xi+1];
	gy += coeff[0][0][1] * bufy[zi+1][yi][xi]   + coeff[0][1][1] * bufy[zi+1][yi][xi+1];
	gy += coeff[1][0][1] * bufy[zi+1][yi+1][xi] + coeff[1][1][1] * bufy[zi+1][yi+1][xi+1];

	gz  = coeff[0][0][0] * bufz[zi][yi][xi]     + coeff[0][1][0] * bufz[zi][yi][xi+1];
	gz += coeff[1][0][0] * bufz[zi][yi+1][xi]   + coeff[1][1][0] * bufz[zi][yi+1][xi+1];
	gz += coeff[0][0][1] * bufz[zi+1][yi][xi]   + coeff[0][1][1] * bufz[zi+1][yi][xi+1];
	gz += coeff[1][0][1] * bufz[zi+1][yi+1][xi] + coeff[1][1][1] * bufz[zi+1][yi+1][xi+1];

	r[0][0] = gx * vec_propre[0] + gy * vec_propre[3] + gz * vec_propre[6];

	/*--- 2eme derivee ---*/
	xr = (double)(x) - sigma * vec_propre[0];   
	yr = (double)(y) - sigma * vec_propre[3];
	zr = (double)(z) - sigma * vec_propre[6];
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= dimz1) ) continue;
	xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	/*--- interpolation de la 2eme derivee ---*/
	coeff[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	coeff[0][1][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	coeff[1][0][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	coeff[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	coeff[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	coeff[0][1][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	coeff[1][0][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	coeff[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */

	gx  = coeff[0][0][0] * bufx[zi][yi][xi]     + coeff[0][1][0] * bufx[zi][yi][xi+1];
	gx += coeff[1][0][0] * bufx[zi][yi+1][xi]   + coeff[1][1][0] * bufx[zi][yi+1][xi+1];
	gx += coeff[0][0][1] * bufx[zi+1][yi][xi]   + coeff[0][1][1] * bufx[zi+1][yi][xi+1];
	gx += coeff[1][0][1] * bufx[zi+1][yi+1][xi] + coeff[1][1][1] * bufx[zi+1][yi+1][xi+1];

	gy  = coeff[0][0][0] * bufy[zi][yi][xi]     + coeff[0][1][0] * bufy[zi][yi][xi+1];
	gy += coeff[1][0][0] * bufy[zi][yi+1][xi]   + coeff[1][1][0] * bufy[zi][yi+1][xi+1];
	gy += coeff[0][0][1] * bufy[zi+1][yi][xi]   + coeff[0][1][1] * bufy[zi+1][yi][xi+1];
	gy += coeff[1][0][1] * bufy[zi+1][yi+1][xi] + coeff[1][1][1] * bufy[zi+1][yi+1][xi+1];

	gz  = coeff[0][0][0] * bufz[zi][yi][xi]     + coeff[0][1][0] * bufz[zi][yi][xi+1];
	gz += coeff[1][0][0] * bufz[zi][yi+1][xi]   + coeff[1][1][0] * bufz[zi][yi+1][xi+1];
	gz += coeff[0][0][1] * bufz[zi+1][yi][xi]   + coeff[0][1][1] * bufz[zi+1][yi][xi+1];
	gz += coeff[1][0][1] * bufz[zi+1][yi+1][xi] + coeff[1][1][1] * bufz[zi+1][yi+1][xi+1];
	r[1][0] = gx * vec_propre[0] + gy * vec_propre[3] + gz * vec_propre[6];

	/*--- on verifie que les deux reponses ont bien des signes opposes ---*/
	if ( r[1][0] * r[0][0] >= 0.0 ) continue;

	/*--- 3eme derivee ---*/
	xr = (double)(x) + sigma * vec_propre[1];   
	yr = (double)(y) + sigma * vec_propre[4];
	zr = (double)(z) + sigma * vec_propre[7];
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= dimz1) ) continue;
	xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	/*--- interpolation de la 3eme derivee ---*/
	coeff[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	coeff[0][1][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	coeff[1][0][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	coeff[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	coeff[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	coeff[0][1][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	coeff[1][0][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	coeff[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */

	gx  = coeff[0][0][0] * bufx[zi][yi][xi]     + coeff[0][1][0] * bufx[zi][yi][xi+1];
	gx += coeff[1][0][0] * bufx[zi][yi+1][xi]   + coeff[1][1][0] * bufx[zi][yi+1][xi+1];
	gx += coeff[0][0][1] * bufx[zi+1][yi][xi]   + coeff[0][1][1] * bufx[zi+1][yi][xi+1];
	gx += coeff[1][0][1] * bufx[zi+1][yi+1][xi] + coeff[1][1][1] * bufx[zi+1][yi+1][xi+1];

	gy  = coeff[0][0][0] * bufy[zi][yi][xi]     + coeff[0][1][0] * bufy[zi][yi][xi+1];
	gy += coeff[1][0][0] * bufy[zi][yi+1][xi]   + coeff[1][1][0] * bufy[zi][yi+1][xi+1];
	gy += coeff[0][0][1] * bufy[zi+1][yi][xi]   + coeff[0][1][1] * bufy[zi+1][yi][xi+1];
	gy += coeff[1][0][1] * bufy[zi+1][yi+1][xi] + coeff[1][1][1] * bufy[zi+1][yi+1][xi+1];

	gz  = coeff[0][0][0] * bufz[zi][yi][xi]     + coeff[0][1][0] * bufz[zi][yi][xi+1];
	gz += coeff[1][0][0] * bufz[zi][yi+1][xi]   + coeff[1][1][0] * bufz[zi][yi+1][xi+1];
	gz += coeff[0][0][1] * bufz[zi+1][yi][xi]   + coeff[0][1][1] * bufz[zi+1][yi][xi+1];
	gz += coeff[1][0][1] * bufz[zi+1][yi+1][xi] + coeff[1][1][1] * bufz[zi+1][yi+1][xi+1];
	r[0][1] = gx * vec_propre[1] + gy * vec_propre[4] + gz * vec_propre[7];

	/*--- 4eme derivee ---*/
	xr = (double)(x) - sigma * vec_propre[1];   
	yr = (double)(y) - sigma * vec_propre[4];
	zr = (double)(z) - sigma * vec_propre[7];
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= dimz1) ) continue;
	xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	/*--- interpolation de la 4eme derivee ---*/
	coeff[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	coeff[0][1][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	coeff[1][0][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	coeff[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	coeff[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	coeff[0][1][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	coeff[1][0][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	coeff[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */

	gx  = coeff[0][0][0] * bufx[zi][yi][xi]     + coeff[0][1][0] * bufx[zi][yi][xi+1];
	gx += coeff[1][0][0] * bufx[zi][yi+1][xi]   + coeff[1][1][0] * bufx[zi][yi+1][xi+1];
	gx += coeff[0][0][1] * bufx[zi+1][yi][xi]   + coeff[0][1][1] * bufx[zi+1][yi][xi+1];
	gx += coeff[1][0][1] * bufx[zi+1][yi+1][xi] + coeff[1][1][1] * bufx[zi+1][yi+1][xi+1];

	gy  = coeff[0][0][0] * bufy[zi][yi][xi]     + coeff[0][1][0] * bufy[zi][yi][xi+1];
	gy += coeff[1][0][0] * bufy[zi][yi+1][xi]   + coeff[1][1][0] * bufy[zi][yi+1][xi+1];
	gy += coeff[0][0][1] * bufy[zi+1][yi][xi]   + coeff[0][1][1] * bufy[zi+1][yi][xi+1];
	gy += coeff[1][0][1] * bufy[zi+1][yi+1][xi] + coeff[1][1][1] * bufy[zi+1][yi+1][xi+1];

	gz  = coeff[0][0][0] * bufz[zi][yi][xi]     + coeff[0][1][0] * bufz[zi][yi][xi+1];
	gz += coeff[1][0][0] * bufz[zi][yi+1][xi]   + coeff[1][1][0] * bufz[zi][yi+1][xi+1];
	gz += coeff[0][0][1] * bufz[zi+1][yi][xi]   + coeff[0][1][1] * bufz[zi+1][yi][xi+1];
	gz += coeff[1][0][1] * bufz[zi+1][yi+1][xi] + coeff[1][1][1] * bufz[zi+1][yi+1][xi+1];
	r[1][1] = gx * vec_propre[1] + gy * vec_propre[4] + gz * vec_propre[7];

	/*--- on verifie que les deux reponses ont bien des signes oppposes ---*/
	if ( r[1][1] * r[0][1] >= 0.0 ) continue;

	/*--- on change les signes pour comparer ---*/
	if ( r[0][0] < 0.0 ) r[0][0] = (- r[0][0]);
	if ( r[1][0] < 0.0 ) r[1][0] = (- r[1][0]);
	if ( r[0][1] < 0.0 ) r[0][1] = (- r[0][1]);
	if ( r[1][1] < 0.0 ) r[1][1] = (- r[1][1]);

	/*
	if ( _LOCAL_DEBUG_ == 1 )
	    fprintf( stderr, " M = (%2d,%2d,%2d), R = ( %f %f %f %f )\n", x,y,z,(float)r[0][0],(float)r[1][0],(float)r[0][1],(float)r[1][1]);
        */

	if ( r[0][0] < r[1][0] ) r1 = r[0][0];
	else r1 = r[1][0];
	if ( r[0][1] < r[1][1] ) r2 = r[0][1];
	else r2 = r[1][1];

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
	xr = (double)x + (double)bufxx[z][y][x];
	zr = (double)x + (double)bufyy[z][y][x];
	yr = (double)x + (double)bufzz[z][y][x];
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= dimz1) ) continue;
	xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	/*--- interpolation de la 1ere reponse ---*/
	coeff[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	coeff[0][1][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	coeff[1][0][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	coeff[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	coeff[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	coeff[0][1][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	coeff[1][0][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	coeff[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */
	r[0][0]  = coeff[0][0][0] * bufr[zi][yi][xi]     + coeff[0][1][0] * bufr[zi][yi][xi+1];
	r[0][0] += coeff[1][0][0] * bufr[zi][yi+1][xi]   + coeff[1][1][0] * bufr[zi][yi+1][xi+1];
	r[0][0] += coeff[0][0][1] * bufr[zi+1][yi][xi]   + coeff[0][1][1] * bufr[zi+1][yi][xi+1];
	r[0][0] += coeff[1][0][1] * bufr[zi+1][yi+1][xi] + coeff[1][1][1] * bufr[zi+1][yi+1][xi+1];

	/* if ( r[0][0] > bufr[z][y][x] ) continue; */

	/*--- 2nde reponse ---*/
	xr = (double)x - (double)bufxx[z][y][x];
	zr = (double)x - (double)bufyy[z][y][x];
	yr = (double)x - (double)bufzz[z][y][x];
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= dimz1) ) continue;
	xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	/*--- interpolation de la 2nde reponse ---*/
	coeff[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	coeff[0][1][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	coeff[1][0][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	coeff[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	coeff[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	coeff[0][1][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	coeff[1][0][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	coeff[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */
	r[1][0]  = coeff[0][0][0] * bufr[zi][yi][xi]     + coeff[0][1][0] * bufr[zi][yi][xi+1];
	r[1][0] += coeff[1][0][0] * bufr[zi][yi+1][xi]   + coeff[1][1][0] * bufr[zi][yi+1][xi+1];
	r[1][0] += coeff[0][0][1] * bufr[zi+1][yi][xi]   + coeff[0][1][1] * bufr[zi+1][yi][xi+1];
	r[1][0] += coeff[1][0][1] * bufr[zi+1][yi+1][xi] + coeff[1][1][1] * bufr[zi+1][yi+1][xi+1];

	/* if ( r[1][0] > bufr[z][y][x] ) continue; */

	/*--- 3eme et 4eme reponse ---*/
	xr = (double)x + (double)bufyz[z][y][x];
	zr = (double)x + (double)bufxy[z][y][x];
	yr = (double)x + (double)bufxz[z][y][x];
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= dimz1) ) continue;
	xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	/*--- interpolation de la 1ere reponse ---*/
	coeff[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	coeff[0][1][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	coeff[1][0][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	coeff[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	coeff[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	coeff[0][1][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	coeff[1][0][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	coeff[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */
	r[0][1]  = coeff[0][0][0] * bufr[zi][yi][xi]     + coeff[0][1][0] * bufr[zi][yi][xi+1];
	r[0][1] += coeff[1][0][0] * bufr[zi][yi+1][xi]   + coeff[1][1][0] * bufr[zi][yi+1][xi+1];
	r[0][1] += coeff[0][0][1] * bufr[zi+1][yi][xi]   + coeff[0][1][1] * bufr[zi+1][yi][xi+1];
	r[0][1] += coeff[1][0][1] * bufr[zi+1][yi+1][xi] + coeff[1][1][1] * bufr[zi+1][yi+1][xi+1];

/* 	if ( r[0][1] > bufr[z][y][x] ) continue; */

	/*--- 2nde reponse ---*/
	xr = (double)x - (double)bufyz[z][y][x];
	zr = (double)x - (double)bufxy[z][y][x];
	yr = (double)x - (double)bufxz[z][y][x];
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= dimz1) ) continue;
	xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	/*--- interpolation de la 2nde reponse ---*/
	coeff[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	coeff[0][1][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	coeff[1][0][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	coeff[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	coeff[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	coeff[0][1][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	coeff[1][0][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	coeff[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */
	r[1][1]  = coeff[0][0][0] * bufr[zi][yi][xi]     + coeff[0][1][0] * bufr[zi][yi][xi+1];
	r[1][1] += coeff[1][0][0] * bufr[zi][yi+1][xi]   + coeff[1][1][0] * bufr[zi][yi+1][xi+1];
	r[1][1] += coeff[0][0][1] * bufr[zi+1][yi][xi]   + coeff[0][1][1] * bufr[zi+1][yi][xi+1];
	r[1][1] += coeff[1][0][1] * bufr[zi+1][yi+1][xi] + coeff[1][1][1] * bufr[zi+1][yi+1][xi+1];

/* 	if ( r[1][1] > bufr[z][y][x] ) continue; */

	if ( r[0][0] > bufr[z][y][x] ) continue;
	if ( r[1][0] > bufr[z][y][x] ) continue;
	if ( r[0][1] > bufr[z][y][x] ) continue;
	if ( r[1][1] > bufr[z][y][x] ) continue;
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

    /*--- derivee selon X ---*/
    local_par.derivative.x = VT_DERIVATIVE_1;
    local_par.derivative.y = VT_DERIVATIVE_0;
    local_par.derivative.z = VT_DERIVATIVE_0;
    if ( VT_RecFilterOnImage( im, &(ims->imx), &local_par ) != 1 ) {
	VT_Error( "unable to compute dI/dx", proc );
	return( -1 );
    }

    /*--- derivee selon Y ---*/
    local_par.derivative.x = VT_DERIVATIVE_0;
    local_par.derivative.y = VT_DERIVATIVE_1;
    local_par.derivative.z = VT_DERIVATIVE_0;
    if ( VT_RecFilterOnImage( im, &(ims->imy), &local_par ) != 1 ) {
	VT_Error( "unable to compute dI/dy", proc );
	return( -1 );
    }

    /*--- derivee selon Z ---*/
    local_par.derivative.x = VT_DERIVATIVE_0;
    local_par.derivative.y = VT_DERIVATIVE_0;
    local_par.derivative.z = VT_DERIVATIVE_1;
    if ( VT_RecFilterOnImage( im, &(ims->imz), &local_par ) != 1 ) {
	VT_Error( "unable to compute dI/dz", proc );
	return( -1 );
    }

    /*--- derivee selon XX ---*/
    local_par.derivative.x = VT_DERIVATIVE_2;
    local_par.derivative.y = VT_DERIVATIVE_0;
    local_par.derivative.z = VT_DERIVATIVE_0;
    if ( VT_RecFilterOnImage( im, &(ims->imxx), &local_par ) != 1 ) {
	VT_Error( "unable to compute d^2I/dx^2", proc );
	return( -1 );
    }

    /*--- derivee selon XY ---*/
    local_par.derivative.x = VT_DERIVATIVE_1;
    local_par.derivative.y = VT_DERIVATIVE_1;
    local_par.derivative.z = VT_DERIVATIVE_0;
    if ( VT_RecFilterOnImage( im, &(ims->imxy), &local_par ) != 1 ) {
	VT_Error( "unable to compute d^2I/dxdy", proc );
	return( -1 );
    }

    /*--- derivee selon XZ ---*/
    local_par.derivative.x = VT_DERIVATIVE_1;
    local_par.derivative.y = VT_DERIVATIVE_0;
    local_par.derivative.z = VT_DERIVATIVE_1;
    if ( VT_RecFilterOnImage( im, &(ims->imxz), &local_par ) != 1 ) {
	VT_Error( "unable to compute d^2I/dxdz", proc );
	return( -1 );
    }
    /*--- derivee selon YY ---*/
    local_par.derivative.x = VT_DERIVATIVE_0;
    local_par.derivative.y = VT_DERIVATIVE_2;
    local_par.derivative.z = VT_DERIVATIVE_0;
    if ( VT_RecFilterOnImage( im, &(ims->imyy), &local_par ) != 1 ) {
	VT_Error( "unable to compute d^2I/dy^2", proc );
	return( -1 );
    }

    /*--- derivee selon YZ ---*/
    local_par.derivative.x = VT_DERIVATIVE_0;
    local_par.derivative.y = VT_DERIVATIVE_1;
    local_par.derivative.z = VT_DERIVATIVE_1;
    if ( VT_RecFilterOnImage( im, &(ims->imyz), &local_par ) != 1 ) {
	VT_Error( "unable to compute d^2I/dydz", proc );
	return( -1 );
    }

    /*--- derivee selon ZZ ---*/
    local_par.derivative.x = VT_DERIVATIVE_0;
    local_par.derivative.y = VT_DERIVATIVE_0;
    local_par.derivative.z = VT_DERIVATIVE_2;
    if ( VT_RecFilterOnImage( im, &(ims->imzz), &local_par ) != 1 ) {
	VT_Error( "unable to compute d^2I/dz^2", proc );
	return( -1 );
    }
   
    return( 1 );
}
