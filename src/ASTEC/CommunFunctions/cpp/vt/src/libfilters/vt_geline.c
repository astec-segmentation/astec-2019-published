
#include <vt_geline.h>

#define COMPUTE2D_NEW(x,y) prod = (x) * v1[0] + (y) * v1[1];                \
	                   if ( prod > 0.0 )      f = prod * prod * irp;    \
                           else if ( prod < 0.0 ) f = prod * prod * irn;    \
	                   else                   f  = 0.0;                 \
	                   prod = (x) * v2[0] + (y) * v2[1];                \
	                   f += prod * prod * imin;                         \
	                   r = (double)255.0 - (double)128.0 * f;           \
	                   if ( r <= 0.0 ) new = (u8)0;          \
	                   else if ( r >= 255.0 ) new = (u8)255; \
	                   else new = (u8)( (int)(r+0.5) )





static int VT_AllocImages2D( vt_images *par, vt_image *im );





/* Initialisation des parametres
 */

void VT_Line( vt_line *par )
{
    /*--- multi-echelle ---*/
    par->first_coeff = 1.0;
    par->last_coeff = 1.0;
    par->nb_coeff = 1;
    /*--- type des vaisseaux ---*/
    par->type_structures = VT_WHITE;
    /*--- filtrage ---*/
    VT_RecFilters( &(par->par_filt) );
    par->par_filt.type_filter = VT_RECGAUSSIAN_MARTA;
    par->par_filt.derivative.x = VT_NODERIVATIVE;
    par->par_filt.derivative.y = VT_NODERIVATIVE;
    par->par_filt.derivative.z = VT_NODERIVATIVE;
    par->par_filt.length_continue.x = 15;
    par->par_filt.length_continue.y = 15;
    par->par_filt.length_continue.z = 15;
    par->par_filt.value_coefficient.x = 1.0;
    par->par_filt.value_coefficient.y = 1.0;
    par->par_filt.value_coefficient.z = 1.0;
}





/*============================================================*/

int VT_ExtractMaxima3D( vt_resline *res )
{
    int dimx, dimy, dimz, dimx1, dimy1, dimz1;
    register int x, y, z;
    int xi, yi, zi;
    r32 ***bufr, ***xbuf, ***ybuf, ***zbuf;
    r32 ***x2buf, ***y2buf, ***z2buf;
    u8 ***ebuf;
    double coeff[2][2][2], xr, yr, zr, dx, dy, dz, rep[2][2], r, n;

    /*--- tests (a faire)---*/

    /*--- initialisation ---*/
    dimx = res->imres.dim.x;   dimx1 = dimx - 1;
    dimy = res->imres.dim.y;   dimy1 = dimy - 1;
    dimz = res->imres.dim.z;   dimz1 = dimz - 1;

    /*--- buffers ---*/
    bufr = (r32 ***)(res->imres.array);
    xbuf = (r32 ***)(res->imdirx.array);
    ybuf = (r32 ***)(res->imdiry.array);
    zbuf = (r32 ***)(res->imdirz.array);
    x2buf = (r32 ***)(res->imdir2x.array);
    y2buf = (r32 ***)(res->imdir2y.array);
    z2buf = (r32 ***)(res->imdir2z.array);
    ebuf = (u8 ***)(res->imext.array);

    for ( z = 0; z < dimz; z ++ )
    for ( y = 0; y < dimy; y ++ )
    for ( x = 0; x < dimx; x ++ ) {
	/*--- initialisation ---*/
	ebuf[z][y][x] = (u8)0;
	r = bufr[z][y][x];
	if ( r <= 0.0 ) continue;
	
	/*--- vecteur : la norme doit etre egale a 1 ---*/
	xr = xbuf[z][y][x];   yr = ybuf[z][y][x];   zr = zbuf[z][y][x];
	n = xr * xr + yr * yr + zr * zr;
	if ( n < 0.5 ) continue;
	
	/*--- premier point ---*/
	xr += (double)(x);   yr += (double)(y);   zr += (double)(z);
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= dimz1) ) continue;
	xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	/*--- coefficients de l'interpolation ---*/
	coeff[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	coeff[0][1][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	coeff[1][0][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	coeff[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	coeff[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	coeff[0][1][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	coeff[1][0][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	coeff[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */
	rep[0][0]  = coeff[0][0][0] * bufr[zi][yi][xi]     + coeff[0][1][0] * bufr[zi][yi][xi+1];
	rep[0][0] += coeff[1][0][0] * bufr[zi][yi+1][xi]   + coeff[1][1][0] * bufr[zi][yi+1][xi+1];
	rep[0][0] += coeff[0][0][1] * bufr[zi+1][yi][xi]   + coeff[0][1][1] * bufr[zi+1][yi][xi+1];
	rep[0][0] += coeff[1][0][1] * bufr[zi+1][yi+1][xi] + coeff[1][1][1] * bufr[zi+1][yi+1][xi+1];

	/*--- deuxieme point ---*/
	xr = (double)(x) - (double)xbuf[z][y][x];
	yr = (double)(y) - (double)ybuf[z][y][x];
	zr = (double)(y) - (double)zbuf[z][y][x];
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= dimz1) ) continue;
	xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	/*--- coefficients de l'interpolation ---*/
	coeff[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	coeff[0][1][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	coeff[1][0][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	coeff[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	coeff[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	coeff[0][1][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	coeff[1][0][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	coeff[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */
	rep[1][0]  = coeff[0][0][0] * bufr[zi][yi][xi]     + coeff[0][1][0] * bufr[zi][yi][xi+1];
	rep[1][0] += coeff[1][0][0] * bufr[zi][yi+1][xi]   + coeff[1][1][0] * bufr[zi][yi+1][xi+1];
	rep[1][0] += coeff[0][0][1] * bufr[zi+1][yi][xi]   + coeff[0][1][1] * bufr[zi+1][yi][xi+1];
	rep[1][0] += coeff[1][0][1] * bufr[zi+1][yi+1][xi] + coeff[1][1][1] * bufr[zi+1][yi+1][xi+1];

	/*--- vecteur : la norme doit etre egale a 1 ---*/
	xr = x2buf[z][y][x];   yr = y2buf[z][y][x];   zr = z2buf[z][y][x];
	n = xr * xr + yr * yr + zr * zr;
	if ( n < 0.5 ) continue;

	/*--- troisieme point ---*/
	xr += (double)(x);   yr += (double)(y);   zr += (double)(z);
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= dimz1) ) continue;
	xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	/*--- coefficients de l'interpolation ---*/
	coeff[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	coeff[0][1][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	coeff[1][0][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	coeff[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	coeff[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	coeff[0][1][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	coeff[1][0][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	coeff[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */
	rep[0][1]  = coeff[0][0][0] * bufr[zi][yi][xi]     + coeff[0][1][0] * bufr[zi][yi][xi+1];
	rep[0][1] += coeff[1][0][0] * bufr[zi][yi+1][xi]   + coeff[1][1][0] * bufr[zi][yi+1][xi+1];
	rep[0][1] += coeff[0][0][1] * bufr[zi+1][yi][xi]   + coeff[0][1][1] * bufr[zi+1][yi][xi+1];
	rep[0][1] += coeff[1][0][1] * bufr[zi+1][yi+1][xi] + coeff[1][1][1] * bufr[zi+1][yi+1][xi+1];

	/*--- quatrieme point ---*/
	xr = (double)(x) - (double)x2buf[z][y][x];
	yr = (double)(y) - (double)y2buf[z][y][x];
	zr = (double)(y) - (double)z2buf[z][y][x];
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= dimz1) ) continue;
	xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	/*--- coefficients de l'interpolation ---*/
	coeff[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	coeff[0][1][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	coeff[1][0][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	coeff[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	coeff[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	coeff[0][1][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	coeff[1][0][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	coeff[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */
	rep[1][1]  = coeff[0][0][0] * bufr[zi][yi][xi]     + coeff[0][1][0] * bufr[zi][yi][xi+1];
	rep[1][1] += coeff[1][0][0] * bufr[zi][yi+1][xi]   + coeff[1][1][0] * bufr[zi][yi+1][xi+1];
	rep[1][1] += coeff[0][0][1] * bufr[zi+1][yi][xi]   + coeff[0][1][1] * bufr[zi+1][yi][xi+1];
	rep[1][1] += coeff[1][0][1] * bufr[zi+1][yi+1][xi] + coeff[1][1][1] * bufr[zi+1][yi+1][xi+1];

	if ( rep[0][0] > bufr[z][y][x] ) continue;
	if ( rep[1][0] > bufr[z][y][x] ) continue;
	if ( rep[0][1] > bufr[z][y][x] ) continue;
	if ( rep[1][1] > bufr[z][y][x] ) continue;

	/*--- c'est un extrema ---*/
	ebuf[z][y][x] = (u8)255;
    }
    return( 1 );
}





int VT_ExtractMaxRad3D( vt_resline *res, vt_images *ims, vt_line *par )
{
    int dimx, dimy, dimz, dimx1, dimy1, dimz1;
    register int x, y, z;
    int i, xi, yi, zi, s, end;
    r32 ***rbuf, ***xbuf, ***ybuf, ***zbuf;
    r32 ***x2buf, ***y2buf, ***z2buf;
    r32 ***gxbuf, ***gybuf, ***gzbuf, ***auxbuf=(r32 ***)NULL;
    r32 ***pbuf, ***nbuf, ***p2buf, ***n2buf;
    u8 ***ebuf, ***sbuf;
    double *sigma=(double*)NULL, step, dep;
    double coeff[2][2][2], xr, yr, zr, dx, dy, dz, rep[2][2], r, r1, r2;
    double vx, vy, vz, n, gx, gy, gz, mul=1.0;

    /*--- tests (a faire)---*/

    /*--- calcul des sigmas ---*/
    sigma = (double*)VT_Malloc( (unsigned int)( (par->nb_coeff + 1) * sizeof( double ) ) );
    sigma[0] = 0.0;
    sigma[1] = (double)(par->first_coeff);
    if ( par->nb_coeff > 1 ) {
	step = (double)(par->last_coeff - par->first_coeff) / (double)(par->nb_coeff - 1);
	for (i = 1; i < par->nb_coeff; i ++ )
	    sigma[i + 1] = (double)(par->first_coeff) + (double)(i) * step;
    }

    /*--- initialisation ---*/
    dimx = res->imres.dim.x;   dimx1 = dimx - 1;
    dimy = res->imres.dim.y;   dimy1 = dimy - 1;
    dimz = res->imres.dim.z;   dimz1 = dimz - 1;

    /*--- buffers ---*/
    rbuf = (r32 ***)(res->imres.array);
    xbuf = (r32 ***)(res->imdirx.array);
    ybuf = (r32 ***)(res->imdiry.array);
    zbuf = (r32 ***)(res->imdirz.array);
    x2buf = (r32 ***)(res->imdir2x.array);
    y2buf = (r32 ***)(res->imdir2y.array);
    z2buf = (r32 ***)(res->imdir2z.array);
    ebuf = (u8 ***)(res->imext.array);
    gxbuf = (r32 ***)(ims->imx.array);
    gybuf = (r32 ***)(ims->imy.array);
    gzbuf = (r32 ***)(ims->imy.array);
    sbuf = (u8 ***)(res->imscale.array);
    pbuf = (r32 ***)(res->imradp.array);
    nbuf = (r32 ***)(res->imradn.array);
    p2buf = (r32 ***)(res->imrad2p.array);
    n2buf = (r32 ***)(res->imrad2n.array);

    for ( z = 0; z < dimz; z ++ )
    for ( y = 0; y < dimy; y ++ )
    for ( x = 0; x < dimx; x ++ ) {
	/*--- initialisation ---*/
	r = rbuf[z][y][x];
	s = (int)sbuf[z][y][x];
	ebuf[z][y][x] = sbuf[z][y][x] = (u8)0;
	pbuf[z][y][x] = nbuf[z][y][x] = p2buf[z][y][x] = n2buf[z][y][x] = 0.0;

	if ( r <= 0.0 ) continue;
	
	/*--- vecteur : la norme doit etre egale a 1 ---*/
	vx = xbuf[z][y][x];   vy = ybuf[z][y][x];   vz = zbuf[z][y][x];
	n = vx * vx + vy * vy + vz * vz;
	if ( n < 0.5 ) continue;
	
	/*--- premier point ---*/
	xr = (double)(x) + vx;   yr = (double)(y) + vy;   zr = (double)(z) + vz;
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= dimz1) ) continue;
	xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	/*--- coefficients de l'interpolation ---*/
	coeff[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	coeff[0][1][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	coeff[1][0][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	coeff[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	coeff[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	coeff[0][1][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	coeff[1][0][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	coeff[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */
	rep[0][0]  = coeff[0][0][0] * rbuf[zi][yi][xi]     + coeff[0][1][0] * rbuf[zi][yi][xi+1];
	rep[0][0] += coeff[1][0][0] * rbuf[zi][yi+1][xi]   + coeff[1][1][0] * rbuf[zi][yi+1][xi+1];
	rep[0][0] += coeff[0][0][1] * rbuf[zi+1][yi][xi]   + coeff[0][1][1] * rbuf[zi+1][yi][xi+1];
	rep[0][0] += coeff[1][0][1] * rbuf[zi+1][yi+1][xi] + coeff[1][1][1] * rbuf[zi+1][yi+1][xi+1];

	if ( rep[0][0] > rbuf[z][y][x] ) continue;

	/*--- deuxieme point ---*/
	xr = (double)(x) - vx;   yr = (double)(y) - vy;   zr = (double)(z) - vz;
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= dimz1) ) continue;
	xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	/*--- coefficients de l'interpolation ---*/
	coeff[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	coeff[0][1][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	coeff[1][0][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	coeff[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	coeff[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	coeff[0][1][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	coeff[1][0][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	coeff[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */
	rep[1][0]  = coeff[0][0][0] * rbuf[zi][yi][xi]     + coeff[0][1][0] * rbuf[zi][yi][xi+1];
	rep[1][0] += coeff[1][0][0] * rbuf[zi][yi+1][xi]   + coeff[1][1][0] * rbuf[zi][yi+1][xi+1];
	rep[1][0] += coeff[0][0][1] * rbuf[zi+1][yi][xi]   + coeff[0][1][1] * rbuf[zi+1][yi][xi+1];
	rep[1][0] += coeff[1][0][1] * rbuf[zi+1][yi+1][xi] + coeff[1][1][1] * rbuf[zi+1][yi+1][xi+1];

	if ( rep[1][0] > rbuf[z][y][x] ) continue;

	/*--- vecteur : la norme doit etre egale a 1 ---*/
	vx = x2buf[z][y][x];   vy = y2buf[z][y][x];   vz = z2buf[z][y][x];
	n = xr * xr + yr * yr + zr * zr;
	if ( n < 0.5 ) continue;

	/*--- troisieme point ---*/
	xr = (double)(x) + vx;   yr = (double)(y) + vy;   zr = (double)(z) + vz;
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= dimz1) ) continue;
	xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	/*--- coefficients de l'interpolation ---*/
	coeff[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	coeff[0][1][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	coeff[1][0][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	coeff[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	coeff[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	coeff[0][1][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	coeff[1][0][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	coeff[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */
	rep[0][1]  = coeff[0][0][0] * rbuf[zi][yi][xi]     + coeff[0][1][0] * rbuf[zi][yi][xi+1];
	rep[0][1] += coeff[1][0][0] * rbuf[zi][yi+1][xi]   + coeff[1][1][0] * rbuf[zi][yi+1][xi+1];
	rep[0][1] += coeff[0][0][1] * rbuf[zi+1][yi][xi]   + coeff[0][1][1] * rbuf[zi+1][yi][xi+1];
	rep[0][1] += coeff[1][0][1] * rbuf[zi+1][yi+1][xi] + coeff[1][1][1] * rbuf[zi+1][yi+1][xi+1];

	if ( rep[0][1] > rbuf[z][y][x] ) continue;

	/*--- quatrieme point ---*/
	xr = (double)(x) - vx;   yr = (double)(y) - vy;   zr = (double)(z) - vz;
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= dimz1) ) continue;
	xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	/*--- coefficients de l'interpolation ---*/
	coeff[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	coeff[0][1][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	coeff[1][0][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	coeff[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	coeff[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	coeff[0][1][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	coeff[1][0][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	coeff[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */
	rep[1][1]  = coeff[0][0][0] * rbuf[zi][yi][xi]     + coeff[0][1][0] * rbuf[zi][yi][xi+1];
	rep[1][1] += coeff[1][0][0] * rbuf[zi][yi+1][xi]   + coeff[1][1][0] * rbuf[zi][yi+1][xi+1];
	rep[1][1] += coeff[0][0][1] * rbuf[zi+1][yi][xi]   + coeff[0][1][1] * rbuf[zi+1][yi][xi+1];
	rep[1][1] += coeff[1][0][1] * rbuf[zi+1][yi+1][xi] + coeff[1][1][1] * rbuf[zi+1][yi+1][xi+1];

	if ( rep[1][1] > rbuf[z][y][x] ) continue;

	/*--- c'est un extrema ---*/
	ebuf[z][y][x] = (u8)255;
	sbuf[z][y][x] = (u8)s;
	pbuf[z][y][x] = nbuf[z][y][x] = p2buf[z][y][x] = n2buf[z][y][x] = (float)sigma[ s ];

	/*--- calcul des rayons du vaisseau :
	      on cherche des structures blanches sur fond noir,
	      on cherche, dans les directions orthogonales du vaisseau
	      des points qui sont extrema locaux de la derivee premiere
	      dans cette direction. voir RQ en 2D.
	---*/
	for ( i = 0; i < 4; i ++ ) {
	  if ( i == 0 ) { 
	    mul = 1.0;   auxbuf = pbuf; 
	    vx = xbuf[z][y][x];   vy = ybuf[z][y][x];   vz = zbuf[z][y][x];
	  }
	  if ( i == 1 ) {
	    mul = (-1.0);   auxbuf = nbuf;
	    vx = xbuf[z][y][x];   vy = ybuf[z][y][x];   vz = zbuf[z][y][x];
	  }
	  if ( i == 2 ) { 
	    mul = 1.0;   auxbuf = p2buf; 
	    vx = x2buf[z][y][x];   vy = y2buf[z][y][x];   vz = z2buf[z][y][x];
	  }
	  if ( i == 3 ) { 
	    mul = (-1.0);   auxbuf = n2buf; 
	    vx = x2buf[z][y][x];   vy = y2buf[z][y][x];   vz = z2buf[z][y][x];
	  }
	  /*--- on interpole la derivee dans le sens du vecteur ---*/
	  /*-------------------------------------------------------*/
	  end = 0;
	  for ( dep = sigma[ s ]; (dep >= 0.0) && (end == 0); dep -= 1.0 ) {
	  /* for ( dep = 1.0; (dep <=  sigma[ s ]) && (end == 0); dep += 1.0 ) { */
	    /*--- au point + dep * vecteur ---*/
	    /*--- au point + dep * vecteur ---*/
	    xr = (double)(x) + mul * dep *  vx;
	    yr = (double)(y) + mul * dep *  vy;
	    zr = (double)(y) + mul * dep *  vz;
	    if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= dimz1) ) continue;
	    xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	    dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	    /*--- coefficients de l'interpolation ---*/
	    coeff[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	    coeff[0][1][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	    coeff[1][0][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	    coeff[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	    coeff[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	    coeff[0][1][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	    coeff[1][0][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	    coeff[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */
	    gx  = coeff[0][0][0] * gxbuf[zi][yi][xi]     + coeff[0][1][0] * gxbuf[zi][yi][xi+1];
	    gx += coeff[1][0][0] * gxbuf[zi][yi+1][xi]   + coeff[1][1][0] * gxbuf[zi][yi+1][xi+1];
	    gx += coeff[0][0][1] * gxbuf[zi+1][yi][xi]   + coeff[0][1][1] * gxbuf[zi+1][yi][xi+1];
	    gx += coeff[1][0][1] * gxbuf[zi+1][yi+1][xi] + coeff[1][1][1] * gxbuf[zi+1][yi+1][xi+1];
	    gy  = coeff[0][0][0] * gybuf[zi][yi][xi]     + coeff[0][1][0] * gybuf[zi][yi][xi+1];
	    gy += coeff[1][0][0] * gybuf[zi][yi+1][xi]   + coeff[1][1][0] * gybuf[zi][yi+1][xi+1];
	    gy += coeff[0][0][1] * gybuf[zi+1][yi][xi]   + coeff[0][1][1] * gybuf[zi+1][yi][xi+1];
	    gy += coeff[1][0][1] * gybuf[zi+1][yi+1][xi] + coeff[1][1][1] * gybuf[zi+1][yi+1][xi+1];
	    gz  = coeff[0][0][0] * gzbuf[zi][yi][xi]     + coeff[0][1][0] * gzbuf[zi][yi][xi+1];
	    gz += coeff[1][0][0] * gzbuf[zi][yi+1][xi]   + coeff[1][1][0] * gzbuf[zi][yi+1][xi+1];
	    gz += coeff[0][0][1] * gzbuf[zi+1][yi][xi]   + coeff[0][1][1] * gzbuf[zi+1][yi][xi+1];
	    gz += coeff[1][0][1] * gzbuf[zi+1][yi+1][xi] + coeff[1][1][1] * gzbuf[zi+1][yi+1][xi+1];
	    r = gx * vx + gy * vy + gz * vz;
	    
	    /*--- est-ce un point interessant a conserver ? 
	      les structures blanches sur fond noir ont un gradient negatif ---*/
	    if ( r > 0.0 ) continue;

	    /*--- au point + (dep-1) * vecteur ---*/
	    xr = (double)(x) + mul * (dep - 1.0) *  vx;
	    yr = (double)(y) + mul * (dep - 1.0) *  vy;
	    zr = (double)(z) + mul * (dep - 1.0) *  vz;
	    if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= dimz1) ) continue;
	    xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	    dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	    /*--- coefficients de l'interpolation ---*/
	    coeff[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	    coeff[0][1][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	    coeff[1][0][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	    coeff[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	    coeff[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	    coeff[0][1][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	    coeff[1][0][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	    coeff[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */
	    gx  = coeff[0][0][0] * gxbuf[zi][yi][xi]     + coeff[0][1][0] * gxbuf[zi][yi][xi+1];
	    gx += coeff[1][0][0] * gxbuf[zi][yi+1][xi]   + coeff[1][1][0] * gxbuf[zi][yi+1][xi+1];
	    gx += coeff[0][0][1] * gxbuf[zi+1][yi][xi]   + coeff[0][1][1] * gxbuf[zi+1][yi][xi+1];
	    gx += coeff[1][0][1] * gxbuf[zi+1][yi+1][xi] + coeff[1][1][1] * gxbuf[zi+1][yi+1][xi+1];
	    gy  = coeff[0][0][0] * gybuf[zi][yi][xi]     + coeff[0][1][0] * gybuf[zi][yi][xi+1];
	    gy += coeff[1][0][0] * gybuf[zi][yi+1][xi]   + coeff[1][1][0] * gybuf[zi][yi+1][xi+1];
	    gy += coeff[0][0][1] * gybuf[zi+1][yi][xi]   + coeff[0][1][1] * gybuf[zi+1][yi][xi+1];
	    gy += coeff[1][0][1] * gybuf[zi+1][yi+1][xi] + coeff[1][1][1] * gybuf[zi+1][yi+1][xi+1];
	    gz  = coeff[0][0][0] * gzbuf[zi][yi][xi]     + coeff[0][1][0] * gzbuf[zi][yi][xi+1];
	    gz += coeff[1][0][0] * gzbuf[zi][yi+1][xi]   + coeff[1][1][0] * gzbuf[zi][yi+1][xi+1];
	    gz += coeff[0][0][1] * gzbuf[zi+1][yi][xi]   + coeff[0][1][1] * gzbuf[zi+1][yi][xi+1];
	    gz += coeff[1][0][1] * gzbuf[zi+1][yi+1][xi] + coeff[1][1][1] * gzbuf[zi+1][yi+1][xi+1];
	    r1 = gx * vx + gy * vy + gz * vz;

	    /*--- on garde le point :
	      - structures blanches sur fond noir
	      r(M + d*v) < r(M + (d-1)*v) et r(M + d*v) <= r(M + (d+1)*v)
	      ---*/
	    if ( r >= r1 ) continue;

	    xr = (double)(x) + mul * (dep + 1.0) *  vx;
	    yr = (double)(y) + mul * (dep + 1.0) *  vy;
	    zr = (double)(z) + mul * (dep + 1.0) *  vz;
	    if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= dimz1) ) continue;
	    xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	    dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	    /*--- coefficients de l'interpolation ---*/
	    coeff[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	    coeff[0][1][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	    coeff[1][0][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	    coeff[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	    coeff[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	    coeff[0][1][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	    coeff[1][0][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	    coeff[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */
	    gx  = coeff[0][0][0] * gxbuf[zi][yi][xi]     + coeff[0][1][0] * gxbuf[zi][yi][xi+1];
	    gx += coeff[1][0][0] * gxbuf[zi][yi+1][xi]   + coeff[1][1][0] * gxbuf[zi][yi+1][xi+1];
	    gx += coeff[0][0][1] * gxbuf[zi+1][yi][xi]   + coeff[0][1][1] * gxbuf[zi+1][yi][xi+1];
	    gx += coeff[1][0][1] * gxbuf[zi+1][yi+1][xi] + coeff[1][1][1] * gxbuf[zi+1][yi+1][xi+1];
	    gy  = coeff[0][0][0] * gybuf[zi][yi][xi]     + coeff[0][1][0] * gybuf[zi][yi][xi+1];
	    gy += coeff[1][0][0] * gybuf[zi][yi+1][xi]   + coeff[1][1][0] * gybuf[zi][yi+1][xi+1];
	    gy += coeff[0][0][1] * gybuf[zi+1][yi][xi]   + coeff[0][1][1] * gybuf[zi+1][yi][xi+1];
	    gy += coeff[1][0][1] * gybuf[zi+1][yi+1][xi] + coeff[1][1][1] * gybuf[zi+1][yi+1][xi+1];
	    gz  = coeff[0][0][0] * gzbuf[zi][yi][xi]     + coeff[0][1][0] * gzbuf[zi][yi][xi+1];
	    gz += coeff[1][0][0] * gzbuf[zi][yi+1][xi]   + coeff[1][1][0] * gzbuf[zi][yi+1][xi+1];
	    gz += coeff[0][0][1] * gzbuf[zi+1][yi][xi]   + coeff[0][1][1] * gzbuf[zi+1][yi][xi+1];
	    gz += coeff[1][0][1] * gzbuf[zi+1][yi+1][xi] + coeff[1][1][1] * gzbuf[zi+1][yi+1][xi+1];
	    r2 = gx * vx + gy * vy + gz * vz;

	    /*--- on garde le point :
	      - structures blanches sur fond noir
	      r(M + d*v) < r(M + (d+1)*v) et r(M + d*v) <= r(M + (d+1)*v)
	      ---*/
	    if ( r <= r2 ) {
	      end = 1;
	      auxbuf[z][y][x] = dep + ( r2 - r1 ) / ( 2.0 * ( 2.0 * r - r1 - r2 ) );
	      if ( auxbuf[z][y][x] <= 0.0 ) auxbuf[z][y][x] = sigma[1];
	    }
	  }
	}
    }

    VT_Free( (void**)&sigma );
    return( 1 );
}





int VT_ExtractMaxima2D( vt_resline *res )
{
    int dimx, dimy, dimz, dimx1, dimy1;
    register int x, y, z;
    int xi, yi;
    r32 ***rbuf, ***xbuf, ***ybuf;
    u8 ***ebuf;
    double coeff[2][2], xr, yr, dx, dy, r1, r2, r, n;

    /*--- tests ---*/

    /*--- initialisation ---*/
    dimx = res->imres.dim.x;   dimx1 = dimx - 1;
    dimy = res->imres.dim.y;   dimy1 = dimy - 1;
    dimz = res->imres.dim.z;

    /*--- buffers ---*/
    rbuf = (r32 ***)(res->imres.array);
    xbuf = (r32 ***)(res->imdirx.array);
    ybuf = (r32 ***)(res->imdiry.array);
    ebuf = (u8 ***)(res->imext.array);

    for ( z = 0; z < dimz; z ++ )
    for ( y = 0; y < dimy; y ++ )
    for ( x = 0; x < dimx; x ++ ) {
	/*--- initialisation ---*/
	ebuf[z][y][x] = (u8)0;
	r = rbuf[z][y][x];
	if ( r <= 0.0 ) continue;
	
	/*--- vecteur : la norme doit etre egale a 1 ---*/
	xr = xbuf[z][y][x];   yr = ybuf[z][y][x];
	n = xr * xr + yr * yr;
	if ( n < 0.5 ) continue;
	
	/*--- premier point ---*/
	xr += (double)(x);   yr += (double)(y);
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) ) continue;
	xi = (int)xr;           yi = (int)yr;
	dx = xr - (double)xi;   dy = yr - (double)yi;
	/*--- coefficients de l'interpolation ---*/
	coeff[0][0] = (1.0 - dx) * (1.0 - dy); /* xi,   yi   */
	coeff[0][1] = dx         * (1.0 - dy); /* xi+1, yi   */
	coeff[1][0] = (1.0 - dx) * dy;         /* xi,   yi+1 */
	coeff[1][1] = dx         * dy;         /* xi+1, yi+1 */
	r1  = coeff[0][0] * rbuf[z][yi][xi]   + coeff[0][1] * rbuf[z][yi][xi+1];
	r1 += coeff[1][0] * rbuf[z][yi+1][xi] + coeff[1][1] * rbuf[z][yi+1][xi+1];
	if ( r1 > r ) continue;
	/*--- second point ---*/
	xr = (double)(x) - (double)xbuf[z][y][x];
	yr = (double)(y) - (double)ybuf[z][y][x];
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) ) continue;
	xi = (int)xr;           yi = (int)yr;
	dx = xr - (double)xi;   dy = yr - (double)yi;
	/*--- coefficients de l'interpolation ---*/
	coeff[0][0] = (1.0 - dx) * (1.0 - dy); /* xi,   yi   */
	coeff[0][1] = dx         * (1.0 - dy); /* xi+1, yi   */
	coeff[1][0] = (1.0 - dx) * dy;         /* xi,   yi+1 */
	coeff[1][1] = dx         * dy;         /* xi+1, yi+1 */
	r2  = coeff[0][0] * rbuf[z][yi][xi]   + coeff[0][1] * rbuf[z][yi][xi+1];
	r2 += coeff[1][0] * rbuf[z][yi+1][xi] + coeff[1][1] * rbuf[z][yi+1][xi+1];
	if ( r2 > r ) continue;
	
	/*--- c'est un extrema ---*/
	ebuf[z][y][x] = (u8)255;
    }
    return( 1 );
}





int VT_ExtractMaxRad2D( vt_resline *res, vt_images *ims, vt_line *par )
{
    int dimx, dimy, dimz, dimx1, dimy1;
    register int x, y, z;
    int i, xi, yi, s, end;
    r32 ***rbuf, ***xbuf, ***ybuf, ***pbuf, ***nbuf;
    r32 ***gxbuf, ***gybuf, ***auxbuf=(r32 ***)NULL;
    u8 ***ebuf, ***sbuf;
    double *sigma=(double*)NULL;
    double step, coeff[2][2], dep, xr, yr, dx, dy, r1, r2, r;
    double vx, vy, n, gx, gy, mul=1.0;

    /*--- tests ---*/

    /*--- calcul des sigmas ---*/
    sigma = (double*)VT_Malloc( (unsigned int)( (par->nb_coeff + 1) * sizeof( double ) ) );
    sigma[0] = 0.0;
    sigma[1] = (double)(par->first_coeff);
    if ( par->nb_coeff > 1 ) {
	step = (double)(par->last_coeff - par->first_coeff) / (double)(par->nb_coeff - 1);
	for (i = 1; i < par->nb_coeff; i ++ )
	    sigma[i + 1] = (double)(par->first_coeff) + (double)(i) * step;
    }

    /*--- initialisation ---*/
    dimx = res->imres.dim.x;   dimx1 = dimx - 1;
    dimy = res->imres.dim.y;   dimy1 = dimy - 1;
    dimz = res->imres.dim.z;

    /*--- buffers ---*/
    rbuf = (r32 ***)(res->imres.array);
    xbuf = (r32 ***)(res->imdirx.array);
    ybuf = (r32 ***)(res->imdiry.array);
    ebuf = (u8 ***)(res->imext.array);
    sbuf = (u8 ***)(res->imscale.array);
    pbuf = (r32 ***)(res->imradp.array);
    nbuf = (r32 ***)(res->imradn.array);
    gxbuf = (r32 ***)(ims->imx.array);
    gybuf = (r32 ***)(ims->imy.array);

    for ( z = 0; z < dimz; z ++ )
    for ( y = 0; y < dimy; y ++ )
    for ( x = 0; x < dimx; x ++ ) {
	/*--- initialisation ---*/
	r = rbuf[z][y][x];
	s = (int)sbuf[z][y][x];
	ebuf[z][y][x] = sbuf[z][y][x] = (u8)0;
	pbuf[z][y][x] = nbuf[z][y][x] = 0.0;
	
	if ( r <= 0.0 ) continue;
	
	/*--- vecteur : la norme doit etre egale a 1 ---*/
	vx = xbuf[z][y][x];   vy = ybuf[z][y][x];
	n = vx * vx + vy * vy;
	if ( n < 0.5 ) continue;

	/*--- premier point ---*/
	xr = (double)(x) + vx;   yr = (double)(y) + vy;
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) ) continue;
	xi = (int)xr;           yi = (int)yr;
	dx = xr - (double)xi;   dy = yr - (double)yi;
	/*--- coefficients de l'interpolation ---*/
	coeff[0][0] = (1.0 - dx) * (1.0 - dy); /* xi,   yi   */
	coeff[0][1] = dx         * (1.0 - dy); /* xi+1, yi   */
	coeff[1][0] = (1.0 - dx) * dy;         /* xi,   yi+1 */
	coeff[1][1] = dx         * dy;         /* xi+1, yi+1 */
	r1  = coeff[0][0] * rbuf[z][yi][xi]   + coeff[0][1] * rbuf[z][yi][xi+1];
	r1 += coeff[1][0] * rbuf[z][yi+1][xi] + coeff[1][1] * rbuf[z][yi+1][xi+1];
	if ( r1 > r ) continue;
	/*--- second point ---*/
	xr = (double)(x) - vx;
	yr = (double)(y) - vy;
	/*--- peut-on faire le calcul ? ---*/
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) ) continue;
	xi = (int)xr;           yi = (int)yr;
	dx = xr - (double)xi;   dy = yr - (double)yi;
	/*--- coefficients de l'interpolation ---*/
	coeff[0][0] = (1.0 - dx) * (1.0 - dy); /* xi,   yi   */
	coeff[0][1] = dx         * (1.0 - dy); /* xi+1, yi   */
	coeff[1][0] = (1.0 - dx) * dy;         /* xi,   yi+1 */
	coeff[1][1] = dx         * dy;         /* xi+1, yi+1 */
	r2  = coeff[0][0] * rbuf[z][yi][xi]   + coeff[0][1] * rbuf[z][yi][xi+1];
	r2 += coeff[1][0] * rbuf[z][yi+1][xi] + coeff[1][1] * rbuf[z][yi+1][xi+1];
	if ( r2 > r ) continue;
	
	/*--- c'est un extrema ---*/
	ebuf[z][y][x] = (u8)255;
	sbuf[z][y][x] = (u8)s;
	pbuf[z][y][x] = nbuf[z][y][x] = (float)sigma[ s ];
    
	/*--- si on ne connait pas le type de structures a chercher 
	      on passe au point suivant                             ---*/
	if ( (par->type_structures != VT_WHITE) && (par->type_structures != VT_BLACK) ) continue;

	/*--- calcul des rayons du vaisseau :
	      on cherche, dans la direction orthogonale du vaisseau
	      des points qui sont extrema locaux de la derivee premiere
	      dans cette direction.
	      RQ : il apparait que cela ne fonctionne pas tjs. Il faudrait
	           peut-etre calcule les extrema du gradient pour le plus
		   petit sigma et chercher un point qui soit extrema du 
		   gradient OU extrema de la derivee dans la direction
		   orthogonale.
	---*/

	for ( i = 0; i < 2; i ++ ) {
	  if ( i == 0 ) { mul = 1.0;      auxbuf = pbuf; }
	  if ( i == 1 ) { mul = (-1.0);   auxbuf = nbuf; }
	  /*--- on interpole la derivee dans le sens du vecteur ---*/
	  /*-------------------------------------------------------*/
	  end = 0;
	  for ( dep = sigma[ s ]; (dep >= 0.0) && (end == 0); dep -= 1.0 ) {
	  /* for ( dep = 1.0; (dep <=  sigma[ s ]) && (end == 0); dep += 1.0 ) { */
	    /*--- au point + dep * vecteur ---*/
	    xr = (double)(x) + mul * dep *  vx;
	    yr = (double)(y) + mul * dep *  vy;
	    if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) ) continue;
	    xi = (int)xr;           yi = (int)yr;
	    dx = xr - (double)xi;   dy = yr - (double)yi;
	    /*--- coefficients de l'interpolation ---*/
	    coeff[0][0] = (1.0 - dx) * (1.0 - dy); /* xi,   yi   */
	    coeff[0][1] = dx         * (1.0 - dy); /* xi+1, yi   */
	    coeff[1][0] = (1.0 - dx) * dy;         /* xi,   yi+1 */
	    coeff[1][1] = dx         * dy;         /* xi+1, yi+1 */
	    gx  = coeff[0][0] * gxbuf[z][yi][xi]   + coeff[0][1] * gxbuf[z][yi][xi+1];
	    gx += coeff[1][0] * gxbuf[z][yi+1][xi] + coeff[1][1] * gxbuf[z][yi+1][xi+1];
	    gy  = coeff[0][0] * gybuf[z][yi][xi]   + coeff[0][1] * gybuf[z][yi][xi+1];
	    gy += coeff[1][0] * gybuf[z][yi+1][xi] + coeff[1][1] * gybuf[z][yi+1][xi+1];
	    r = gx * vx + gy * vy;
	    
	    /*--- est-ce que la derivee a une valeur suffisante ? ---*/
	    /* if (r < rbuf[z][y][x]) continue; */
	    /*--- si je fais ca, le diametre augmente ---*/

	    /*--- est-ce un point interessant a conserver ? 
	      les structures blanches sur fond noir ont un gradient negatif
	      les structures noires sur fond blanc ont un gradient positif ---*/
	    if ( ( r > 0.0) && (par->type_structures == VT_WHITE) ) continue;
	    if ( ( r < 0.0) && (par->type_structures == VT_BLACK) ) continue;
	    /*--- au point + (dep-1) * vecteur ---*/
	    xr = (double)(x) + mul * (dep - 1.0) *  vx;
	    yr = (double)(y) + mul * (dep - 1.0) *  vy;
	    if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) ) continue;
	    xi = (int)xr;           yi = (int)yr;
	    dx = xr - (double)xi;   dy = yr - (double)yi;
	    /*--- coefficients de l'interpolation ---*/
	    coeff[0][0] = (1.0 - dx) * (1.0 - dy); /* xi,   yi   */
	    coeff[0][1] = dx         * (1.0 - dy); /* xi+1, yi   */
	    coeff[1][0] = (1.0 - dx) * dy;         /* xi,   yi+1 */
	    coeff[1][1] = dx         * dy;         /* xi+1, yi+1 */
	    gx  = coeff[0][0] * gxbuf[z][yi][xi]   + coeff[0][1] * gxbuf[z][yi][xi+1];
	    gx += coeff[1][0] * gxbuf[z][yi+1][xi] + coeff[1][1] * gxbuf[z][yi+1][xi+1];
	    gy  = coeff[0][0] * gybuf[z][yi][xi]   + coeff[0][1] * gybuf[z][yi][xi+1];
	    gy += coeff[1][0] * gybuf[z][yi+1][xi] + coeff[1][1] * gybuf[z][yi+1][xi+1];
	    r1 = gx * vx + gy * vy;
	    /*--- on garde le point :
	      - structures blanches sur fond noir
	      r(M + d*v) < r(M + (d-1)*v) et r(M + d*v) <= r(M + (d+1)*v)
	      - structures noires sur fond blanc : 
	      r(M + d*v) > r(M + (d-1)*v) et r(M + d*v) >= r(M + (d+1)*v)
	      ---*/
	    if ( ( r >= r1) && (par->type_structures == VT_WHITE) ) continue;
	    if ( ( r <= r1) && (par->type_structures == VT_BLACK) ) continue;
	    /*--- au point + (dep+1) * vecteur ---*/
	    xr = (double)(x) + mul * (dep + 1.0) *  vx;
	    yr = (double)(y) + mul * (dep + 1.0) *  vy;
	    if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) ) continue;
	    xi = (int)xr;           yi = (int)yr;
	    dx = xr - (double)xi;   dy = yr - (double)yi;
	    /*--- coefficients de l'interpolation ---*/
	    coeff[0][0] = (1.0 - dx) * (1.0 - dy); /* xi,   yi   */
	    coeff[0][1] = dx         * (1.0 - dy); /* xi+1, yi   */
	    coeff[1][0] = (1.0 - dx) * dy;         /* xi,   yi+1 */
	    coeff[1][1] = dx         * dy;         /* xi+1, yi+1 */
	    gx  = coeff[0][0] * gxbuf[z][yi][xi]   + coeff[0][1] * gxbuf[z][yi][xi+1];
	    gx += coeff[1][0] * gxbuf[z][yi+1][xi] + coeff[1][1] * gxbuf[z][yi+1][xi+1];
	    gy  = coeff[0][0] * gybuf[z][yi][xi]   + coeff[0][1] * gybuf[z][yi][xi+1];
	    gy += coeff[1][0] * gybuf[z][yi+1][xi] + coeff[1][1] * gybuf[z][yi+1][xi+1];
	    r2 = gx * vx + gy * vy;
	    /*--- on garde le point :
	      - structures blanches sur fond noir
	      r(M + d*v) < r(M + (d+1)*v) et r(M + d*v) <= r(M + (d+1)*v)
	      - structures noires sur fond blanc : 
	      r(M + d*v) > r(M + (d+1)*v) et r(M + d*v) >= r(M + (d+1)*v)
	      ---*/
	    if ( (( r <= r2) && (par->type_structures == VT_WHITE)) ||
		 (( r >= r2) && (par->type_structures == VT_BLACK)) ) {
	      end = 1;
	      auxbuf[z][y][x] = dep + ( r2 - r1 ) / ( 2.0 * ( 2.0 * r - r1 - r2 ) );
	      if ( auxbuf[z][y][x] <= 0.0 ) auxbuf[z][y][x] = sigma[1];
	    }
	  }
	}
    }

    VT_Free( (void**)&sigma );
    return( 1 );
}





int VT_Reconstruct2D( vt_image *ext, vt_image *scl, vt_line *par )
{
    char *proc="VT_Reconstruct2D";
    double *sigma=(double*)NULL;
    register int is, i, j, d;
    double step, s;
    register int x, y, z;
    int dimx, dimy, dimz;
    u8 ***ebuf, ***sbuf;

    /*--- tests ---*/
    if ( VT_Test2Image( ext, scl, proc ) == -1 ) return( -1 );
    if ( (ext->type != UCHAR) || (scl->type != UCHAR) ) {
	VT_Error( "bad type for images", proc );
	return( -1 );
    }
    
    /*--- calcul des sigmas ---*/
    sigma = (double*)VT_Malloc( (unsigned int)( (par->nb_coeff + 1) * sizeof( double ) ) );
    sigma[0] = 0.0;
    sigma[1] = (double)(par->first_coeff);
    if ( par->nb_coeff > 1 ) {
	step = (double)(par->last_coeff - par->first_coeff) / (double)(par->nb_coeff - 1);
	for (i = 1; i < par->nb_coeff; i ++ )
	    sigma[i + 1] = (double)(par->first_coeff) + (double)(i) * step;
    }

    /*--- buffers ---*/
    ebuf = (u8 ***)(ext->array);
    sbuf = (u8 ***)(scl->array);
    dimx = ext->dim.x;
    dimy = ext->dim.y;
    dimz = ext->dim.z;
    
    /*--- ----*/
    for ( z = 0; z < dimz; z ++ )
    for ( y = 0; y < dimy; y ++ )
    for ( x = 0; x < dimx; x ++ ) {
	if ( sbuf[z][y][x] == (u8)0 ) continue;
	s = sigma[ (int)sbuf[z][y][x] ] *  sigma[ (int)sbuf[z][y][x] ];
	is = (int)(s + 0.5);
	for ( j = 0; j <= is ; j ++ )
	for ( i = 0; i <= j ; i ++ ) {
	    d = i*i + j*j;
	    if ( d > s ) continue;
	    if ( x+i < dimx ) {
		if ( y+j < dimy ) ebuf[z][y+j][x+i] = (u8)255;
		if ( y-j >= 0   ) ebuf[z][y-j][x+i] = (u8)255;
	    }
	    if ( x-i >= 0 ) {
		if ( y+j < dimy ) ebuf[z][y+j][x-i] = (u8)255;
		if ( y-j >= 0   ) ebuf[z][y-j][x-i] = (u8)255;
	    }
	    if ( x+j < dimx ) {
		if ( y+i < dimy ) ebuf[z][y+i][x+j] = (u8)255;
		if ( y-i >= 0   ) ebuf[z][y-i][x+j] = (u8)255;
	    }
	    if ( x-j >= 0 ) {
		if ( y+i < dimy ) ebuf[z][y+i][x-j] = (u8)255;
		if ( y-i >= 0   ) ebuf[z][y-i][x-j] = (u8)255;
	    }
	}
    }

    VT_Free( (void**)&sigma );
    return( 1 );
}





int VT_GreyReconstruct2D( vt_resline *res )
{
    register int borne, i, j;
    register double xd, yd, imin, irp, irn, prod;
    double v1[2], v2[2], rp, rn, r, f;
    register int x, y, z;
    int dimx, dimy, dimz;
    u8 ***ebuf, ***sbuf, new;
    r32 ***xbuf, ***ybuf, ***pbuf, ***nbuf;
    /*--- tests ---*/
    
    /*--- buffers ---*/
    ebuf = (u8 ***)(res->imext.array);
    sbuf = (u8 ***)(res->imscale.array);
    xbuf = (r32 ***)(res->imdirx.array);
    ybuf = (r32 ***)(res->imdiry.array);
    pbuf = (r32 ***)(res->imradp.array);
    nbuf = (r32 ***)(res->imradn.array);
    dimx = res->imext.dim.x;
    dimy = res->imext.dim.y;
    dimz = res->imext.dim.z;
    
    /*--- initialisation ----*/
    for ( z = 0; z < dimz; z ++ )
    for ( y = 0; y < dimy; y ++ )
    for ( x = 0; x < dimx; x ++ ) 
      ebuf[z][y][x] = (u8)0;
    
    /*--- iso-surface a 127.0 ----*/
    for ( z = 0; z < dimz; z ++ )
    for ( y = 0; y < dimy; y ++ )
    for ( x = 0; x < dimx; x ++ ) {
	if ( sbuf[z][y][x] == (u8)0 ) continue;
	/*--- vecteurs ---*/
	v1[0] = xbuf[z][y][x];   v1[1] = ybuf[z][y][x];
	v2[0] = v1[1];           v2[1] = (- v1[0]);
	/*--- bornes : f vaut le max ---*/
	f = pbuf[z][y][x];
	if ( f < nbuf[z][y][x] ) f =  nbuf[z][y][x];
	borne = (int)(f + 0.5);
	if ( borne == 0 ) continue;

	/*--- 1 / rayons au carre : f vaut min positif ---*/
	rp = (double)pbuf[z][y][x];   rn = (double)nbuf[z][y][x];
	if ( (rp > 0.0) && (f > rp) ) f = rp;
	if ( (rn > 0.0) && (f > rn) ) f = rn;
	if ( rp <= 0.0 ) rp = f;
	if ( rn <= 0.0 ) rn = f;

	irp = (double)1.0 / ( rp * rp );
	irn = (double)1.0 / ( rn * rn );
	imin = (double)1.0 / ( f * f );

	for ( j = 0; j <= borne ; j ++ )
	for ( i = 0; i <= j ; i ++ ) {
	  /*--- coordonnees en double ---*/
	  xd = (double)i;
	  yd = (double)j;
	  
	  /*--- i,[j,-j] ---*/
	  if ( x+i < dimx ) {
	    if ( y+j < dimy ) {
	      COMPUTE2D_NEW( xd, yd );
	      if (new > ebuf[z][y+j][x+i]) ebuf[z][y+j][x+i] = new;
	    }
	    if ( y-j >= 0   ) {
	      COMPUTE2D_NEW( xd, (-yd) );
	      if (new > ebuf[z][y-j][x+i]) ebuf[z][y-j][x+i] = new;
	    }
	  }
	  /*--- -i,[j,-j] ---*/
	  if ( x-i >= 0 ) {
	    if ( y+j < dimy ) {
	      COMPUTE2D_NEW( (-xd), yd );
	      if (new > ebuf[z][y+j][x-i]) ebuf[z][y+j][x-i] = new;
	    }
	    if ( y-j >= 0   ) {
	      COMPUTE2D_NEW( (-xd), (-yd) );
	      if (new > ebuf[z][y-j][x-i]) ebuf[z][y-j][x-i] = new;
	    }
	  }
	  /*--- j,[i,-i] ---*/
	  if ( x+j < dimx ) {
	    if ( y+i < dimy ) {
	      COMPUTE2D_NEW( yd, xd );
	      if (new > ebuf[z][y+i][x+j]) ebuf[z][y+i][x+j] = new;
	    }
	    if ( y-i >= 0   ) {
	      COMPUTE2D_NEW( yd, (-xd) );
	      if (new > ebuf[z][y-i][x+j]) ebuf[z][y-i][x+j] = new;
	    }
	  }
	  /*--- -j,[i,-i] ---*/
	  if ( x-j >= 0 ) {
	    if ( y+i < dimy ) {
	      COMPUTE2D_NEW( (-yd), xd );
	      if (new > ebuf[z][y+i][x-j]) ebuf[z][y+i][x-j] = new;
	    }
	    if ( y-i >= 0   ) {
	      COMPUTE2D_NEW( (-yd), (-xd) );
	      if (new > ebuf[z][y-i][x-j]) ebuf[z][y-i][x-j] = new;
	    }
	  }
	}
    }

    return( 1 );
}





/*============================================================*/

void VT_InitRes( vt_resline *par )
{
    VT_Image( &(par->imres) );   /*--- valeur de la reponse maximale ---*/
    VT_Image( &(par->imext) );   /*--- extrema ---*/
    VT_Image( &(par->imscale) ); /*--- echelle ---*/
    VT_Image( &(par->imdirx) );  /*--- coordonnee X du vecteur orthogonal ---*/
    VT_Image( &(par->imdiry) );  /*--- coordonnee Y du vecteur orthogonal ---*/
    VT_Image( &(par->imradp) );  /*--- rayon dans le sens du vecteur ---*/
    VT_Image( &(par->imradn) );  /*--- rayon dans le sens oppose du vecteur ---*/

    VT_Image( &(par->imdirz) );  /*--- coordonnee Z du vecteur orthogonal ---*/
    VT_Image( &(par->imdir2x) );  /*--- coordonnee X du 2eme vecteur orthogonal ---*/
    VT_Image( &(par->imdir2y) );  /*--- coordonnee Y du 2eme vecteur orthogonal ---*/
    VT_Image( &(par->imdir2z) );  /*--- coordonnee Z du 2eme vecteur orthogonal ---*/
    VT_Image( &(par->imrad2p) );  /*--- rayon dans le sens du 2eme vecteur ---*/
    VT_Image( &(par->imrad2n) );  /*--- rayon dans le sens oppose du 2eme vecteur ---*/
}





int VT_AllocRes( vt_resline *par, vt_image *im, int dim )
{
    VT_InitFromImage( &(par->imres),   im, "geline.results.imres", (int)FLOAT );
    VT_InitFromImage( &(par->imext),   im, "geline.results.imext", (int)UCHAR );
    VT_InitFromImage( &(par->imscale), im, "geline.results.imscale", (int)UCHAR );
    VT_InitFromImage( &(par->imdirx),  im, "geline.results.imdirx", (int)FLOAT );
    VT_InitFromImage( &(par->imdiry),  im, "geline.results.imdiry", (int)FLOAT );
    VT_InitFromImage( &(par->imradp),  im, "geline.results.imradp", (int)FLOAT );
    VT_InitFromImage( &(par->imradn),  im, "geline.results.imradn", (int)FLOAT );

    VT_InitFromImage( &(par->imdirz),  im, "geline.results.imdirz", (int)FLOAT );
    VT_InitFromImage( &(par->imdir2x), im, "geline.results.imdir2x", (int)FLOAT );
    VT_InitFromImage( &(par->imdir2y), im, "geline.results.imdir2y", (int)FLOAT );
    VT_InitFromImage( &(par->imdir2z), im, "geline.results.imdir2z", (int)FLOAT );
    VT_InitFromImage( &(par->imrad2p),  im, "geline.results.imrad2p", (int)FLOAT );
    VT_InitFromImage( &(par->imrad2n),  im, "geline.results.imrad2n", (int)FLOAT );

    if ( VT_AllocImage( &(par->imres) ) != 1 ) {
	return( -1 );
    }
    if ( VT_AllocImage( &(par->imext) ) != 1 ) {
	VT_FreeImage( &(par->imres) );
	return( -1 );
    }
    if ( VT_AllocImage( &(par->imscale) ) != 1 ) {
	VT_FreeImage( &(par->imres) );
	VT_FreeImage( &(par->imext) );
	return( -1 );
    }
    if ( VT_AllocImage( &(par->imdirx) ) != 1 ) {
	VT_FreeImage( &(par->imres) );
	VT_FreeImage( &(par->imext) );
	VT_FreeImage( &(par->imscale) );
	return( -1 );
    }
    if ( VT_AllocImage( &(par->imdiry) ) != 1 ) {
	VT_FreeImage( &(par->imres) );
	VT_FreeImage( &(par->imext) );
	VT_FreeImage( &(par->imscale) );
	VT_FreeImage( &(par->imdirx) );
	return( -1 );
    }
     if ( VT_AllocImage( &(par->imradp) ) != 1 ) {
	VT_FreeImage( &(par->imres) );
	VT_FreeImage( &(par->imext) );
	VT_FreeImage( &(par->imscale) );
	VT_FreeImage( &(par->imdirx) );
	VT_FreeImage( &(par->imdiry) );
	return( -1 );
    }
     if ( VT_AllocImage( &(par->imradn) ) != 1 ) {
	VT_FreeImage( &(par->imres) );
	VT_FreeImage( &(par->imext) );
	VT_FreeImage( &(par->imscale) );
	VT_FreeImage( &(par->imdirx) );
	VT_FreeImage( &(par->imdiry) );
	VT_FreeImage( &(par->imradp) );
	return( -1 );
    }
   
    if ( dim == VT_2D ) return( 1 );

    if ( VT_AllocImage( &(par->imdirz) ) != 1 ) {
      VT_FreeRes( par, (int)VT_2D );
      return( -1 );
    }
    if ( VT_AllocImage( &(par->imdir2x) ) != 1 ) {
      VT_FreeRes( par, (int)VT_2D );
      VT_FreeImage( &(par->imdirz) );
      return( -1 );
    }
    if ( VT_AllocImage( &(par->imdir2y) ) != 1 ) {
      VT_FreeRes( par, (int)VT_2D );
      VT_FreeImage( &(par->imdirz) );
      VT_FreeImage( &(par->imdir2x) );
      return( -1 );
    }
    if ( VT_AllocImage( &(par->imdir2z) ) != 1 ) {
      VT_FreeRes( par, (int)VT_2D );
      VT_FreeImage( &(par->imdirz) );
      VT_FreeImage( &(par->imdir2x) );
      VT_FreeImage( &(par->imdir2y) );
      return( -1 );
    }
    if ( VT_AllocImage( &(par->imrad2p) ) != 1 ) {
      VT_FreeRes( par, (int)VT_2D );
      VT_FreeImage( &(par->imdirz) );
      VT_FreeImage( &(par->imdir2x) );
      VT_FreeImage( &(par->imdir2y) );
      VT_FreeImage( &(par->imdir2z) );
      return( -1 );
    }
    if ( VT_AllocImage( &(par->imrad2n) ) != 1 ) {
      VT_FreeRes( par, (int)VT_2D );
      VT_FreeImage( &(par->imdirz) );
      VT_FreeImage( &(par->imdir2x) );
      VT_FreeImage( &(par->imdir2y) );
      VT_FreeImage( &(par->imdir2z) );
      VT_FreeImage( &(par->imrad2n) );
      return( -1 );
    }
    return( 1 );
}





void VT_FreeRes( vt_resline *par, int dim )
{
    VT_FreeImage( &(par->imres) );
    VT_FreeImage( &(par->imext) );
    VT_FreeImage( &(par->imscale) );
    VT_FreeImage( &(par->imdirx) );
    VT_FreeImage( &(par->imdiry) );
    VT_FreeImage( &(par->imradp) );
    VT_FreeImage( &(par->imradn) );
    if ( dim != VT_2D ) {
      VT_FreeImage( &(par->imdirz) );
      VT_FreeImage( &(par->imdir2x) );
      VT_FreeImage( &(par->imdir2y) );
      VT_FreeImage( &(par->imdir2z) );
      VT_FreeImage( &(par->imrad2p) );
      VT_FreeImage( &(par->imrad2n) );
    }
    VT_InitRes( par );
}





/*============================================================*/

void VT_InitImages( vt_images *par )
{
  /*--- declarations pour 2D et 3D ---*/
  VT_Image( &(par->imx) );
  VT_Image( &(par->imy) );
  VT_Image( &(par->imxx) );
  VT_Image( &(par->imxy) );
  VT_Image( &(par->imyy) );
  VT_Image( &(par->imr) );
  VT_Image( &(par->ime) );
  /*--- declarations pour 3D seulement ----*/
  VT_Image( &(par->imz) );
  VT_Image( &(par->imxz) );
  VT_Image( &(par->imyz) );
  VT_Image( &(par->imzz) );
}





static int VT_AllocImages2D( vt_images *par, vt_image *im )
{
    VT_InitFromImage( &(par->imx),  im, "geline.images.imx", (int)FLOAT );
    VT_InitFromImage( &(par->imy),  im, "geline.images.imy", (int)FLOAT );
    VT_InitFromImage( &(par->imxx), im, "geline.images.imxx", (int)FLOAT );
    VT_InitFromImage( &(par->imxy), im, "geline.images.imxy", (int)FLOAT );
    VT_InitFromImage( &(par->imyy), im, "geline.images.imyy", (int)FLOAT );
    VT_InitFromImage( &(par->imr),  im, "geline.images.imr", (int)FLOAT );
    VT_InitFromImage( &(par->ime),  im, "geline.images.ime", (int)UCHAR );

    if ( VT_AllocImage( &(par->imx) ) != 1 ) {
	return( -1 );
    }
    if ( VT_AllocImage( &(par->imy) ) != 1 ) {
	VT_FreeImage( &(par->imx) );
	return( -1 );
    }
    if ( VT_AllocImage( &(par->imxx) ) != 1 ) {
	VT_FreeImage( &(par->imx) );
	VT_FreeImage( &(par->imy) );
	return( -1 );
    }
    if ( VT_AllocImage( &(par->imxy) ) != 1 ) {
	VT_FreeImage( &(par->imx) );
	VT_FreeImage( &(par->imy) );
	VT_FreeImage( &(par->imxx) );
	return( -1 );
    }
    if ( VT_AllocImage( &(par->imyy) ) != 1 ) {
	VT_FreeImage( &(par->imx) );
	VT_FreeImage( &(par->imy) );
	VT_FreeImage( &(par->imxx) );
	VT_FreeImage( &(par->imxy) );
	return( -1 );
    }
    if ( VT_AllocImage( &(par->imr) ) != 1 ) {
	VT_FreeImage( &(par->imx) );
	VT_FreeImage( &(par->imy) );
	VT_FreeImage( &(par->imxx) );
	VT_FreeImage( &(par->imxy) );
	VT_FreeImage( &(par->imyy) );
	return( -1 );
    }
    if ( VT_AllocImage( &(par->ime) ) != 1 ) {
	VT_FreeImage( &(par->imx) );
	VT_FreeImage( &(par->imy) );
	VT_FreeImage( &(par->imxx) );
	VT_FreeImage( &(par->imxy) );
	VT_FreeImage( &(par->imyy) );
	VT_FreeImage( &(par->imr) );
	return( -1 );
    }
    return( 1 );
}





int VT_AllocImages( vt_images *par, vt_image *im, int dim )
{
  if ( VT_AllocImages2D( par, im ) != 1 ) return( -1 );
  if ( dim == VT_2D ) return( 1 );

  VT_InitFromImage( &(par->imz),  im, "geline.images.imz", (int)FLOAT );
  VT_InitFromImage( &(par->imxz), im, "geline.images.imxz", (int)FLOAT );
  VT_InitFromImage( &(par->imyz), im, "geline.images.imyz", (int)FLOAT );
  VT_InitFromImage( &(par->imzz), im, "geline.images.imzz", (int)FLOAT );

  if ( VT_AllocImage( &(par->imz) ) != 1 ) {
    VT_FreeImages( par, (int)VT_2D );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imxz) ) != 1 ) {
    VT_FreeImages( par, (int)VT_2D );
    VT_FreeImage( &(par->imz) );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imyz) ) != 1 ) {
    VT_FreeImages( par, (int)VT_2D );
    VT_FreeImage( &(par->imz) );
    VT_FreeImage( &(par->imxz) );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imzz) ) != 1 ) {
    VT_FreeImages( par, (int)VT_2D );
    VT_FreeImage( &(par->imz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    return( -1 );
  }
  return( 1 );
}





void VT_FreeImages( vt_images *par, int dim )
{
    VT_FreeImage( &(par->imx) );
    VT_FreeImage( &(par->imy) );
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imr) );
    VT_FreeImage( &(par->ime) );
    if (dim != VT_2D ) {
      VT_FreeImage( &(par->imz) );
      VT_FreeImage( &(par->imxz) );
      VT_FreeImage( &(par->imyz) );
      VT_FreeImage( &(par->imzz) );
    }
    VT_InitImages( par );
}

