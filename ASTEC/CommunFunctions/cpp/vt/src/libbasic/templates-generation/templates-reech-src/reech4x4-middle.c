




/* Resampling procedure.

   Work for 3D images, not for vectorial ones.
   
   (double* mat) is the matrix which permits to get
   from resBuf into theBuf. 
   If one only have the matrix from theBuf into resBuf,
   it must be inverted first.

   Soit x le point transforme et ix=(int)x;
   nous allons distinguer les cas suivants :
    x < -0.5               => resultat = 0
    -0.5 <= x < 0.0        => ix=0, on n'interpole pas selon X
    0.0 < x && ix < dimx-1 => on interpole selon X
    x < dimx-0.5           => ix=dimx-1, on n'interpole pas selon X
    x >= dimx-0.5          => resultat = 0

*/

static void *_Reech3DTriLin4x4_TYPE ( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _LinearResamplingParam *p = (_LinearResamplingParam *)parameter;
  
  void* theBuf = p->theBuf;
  int *theDim = p->theDim;
  void* resBuf = p->resBuf;
  int *resDim = p->resDim;
  double* mat = p->mat; 

  size_t i, j, k;
  register int ix, iy, iz;
  register double x, y, z, dx, dy, dz, dxdy,dxdz,dydz,dxdydz;
  register double res;
  double v6, v5, v4;
  size_t rdimx=resDim[0], rdimy=resDim[1];
  int tdimx=theDim[0], tdimy=theDim[1], tdimz=theDim[2];
  int tdimxy=tdimx*tdimy;
  int toffset1=tdimxy+tdimx+1, toffset2=tdimxy-tdimx-1;
  register int t1dimx=tdimx-1, t1dimy=tdimy-1, t1dimz=tdimz-1;
  register double ddimx = (double)tdimx-0.5, ddimy = (double)tdimy-0.5;
  register double ddimz = (double)tdimz-0.5;
  register TYPE *tbuf = (TYPE*)theBuf;
  register TYPE *tpt;
  register TYPE *rbuf = (TYPE*)resBuf;

  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;
  size_t iend, jend;

  k = kfirst = first / (rdimx*rdimy);
  j = jfirst = (first - kfirst*(rdimx*rdimy)) / rdimx;
  i = ifirst = (first - kfirst*(rdimx*rdimy) - jfirst*rdimx);

  klast = last / (rdimx*rdimy);
  jlast = (last - klast*(rdimx*rdimy)) / rdimx;
  ilast = (last - klast*(rdimx*rdimy) - jlast*rdimx);

  rbuf += first;

  for ( ; k<=klast; k++, j=0 ) {
    if ( _verbose_ > 1 )
      fprintf( stderr, "Processing slice %lu\r", k );
    jend = (k==klast) ? jlast+1 : rdimy;
    for ( ; j<jend; j++, i=0 ) {
      iend = (j==jlast && k==klast) ? ilast+1 : rdimx;
      for ( ; i<iend; i++, rbuf++ ) {
        /* computation of the corresponding point coordinates in theBuf */
        x = mat[0] * i +  mat[1] * j + mat[2] * k + mat[3];
        if ((x <= -0.5) || ( x >= ddimx)) { *rbuf = 0; continue; }
        y = mat[4] * i +  mat[5] * j + mat[6] * k + mat[7];
        if ((y <= -0.5) || ( y >= ddimy)) { *rbuf = 0; continue; }
        z = mat[8] * i +  mat[9] * j + mat[10] * k + mat[11];
        if ((z <= -0.5) || ( z >= ddimz)) { *rbuf = 0; continue; }

        /* here, the point lies on the borders or completely inside
           the image */
        ix = (int)x;
        iy = (int)y;
        iz = (int)z;
        tpt = (TYPE *)tbuf;

        /* are we on the border or not ? */
        if ( (x > 0.0) && (ix < t1dimx) &&
             (y > 0.0) && (iy < t1dimy) &&
             (z > 0.0) && (iz < t1dimz) ) {
          /* the corresponding point is in the box defined
             by (ix[+1],iy[+1],iz[+1]) */
          dx = x - ix;
          dy = y - iy;
          dz = z - iz;
          dxdy = dx*dy;
          dxdz = dx*dz;
          dydz = dy*dz;
          dxdydz = dxdy*dz;

          /* we have
             v[7]=dxdydz;            coefficient of tbuf(ix+1,iy+1,iz+1)
             v[6]=dxdz-dxdydz;       coefficient of tbuf(ix+1,iy,  iz+1)
             v[5]=dxdy-dxdydz;       coefficient of tbuf(ix+1,iy+1,iz  )
             v[4]=dx-dxdy-v[6];      coefficient of tbuf(ix+1,iy  ,iz  )
             v[3]=dydz-dxdydz;       coefficient of tbuf(ix  ,iy+1,iz+1)
             v[2]=dz-dydz-v[6];      coefficient of tbuf(ix  ,iy  ,iz+1)
             v[1]=dy-dydz-v[5];      coefficient of tbuf(ix  ,iy+1,iz  )
             v[0]=1-dy-dz+dydz-v[4]; coefficient of tbuf(ix  ,iy  ,iz  )
          */
          tpt += (size_t)ix + (size_t)iy * (size_t)tdimx + (size_t)iz * (size_t)tdimxy + (size_t)toffset1;
          v6 = dxdz-dxdydz;
          v5 = dxdy-dxdydz;
          v4 = dx-dxdy-v6;

          res = 0;
          res += dxdydz * (*tpt);            /* tbuf(ix+1,iy+1,iz+1) */
          tpt --;
          res += (dydz-dxdydz) * (*tpt);     /* tbuf(ix  ,iy+1,iz+1) */
          tpt -= t1dimx;
          res += v6 * (*tpt);                /* tbuf(ix+1  ,iy,  iz+1) */
          tpt --;
          res += (dz-dydz-v6) * (*tpt);      /* tbuf(ix  ,iy  ,iz+1) */
          tpt -= toffset2;
          res += v5 * (*tpt);                /* tbuf(ix+1,iy+1,iz  ) */
          tpt --;
          res += (dy-dydz-v5) * (*tpt);      /* tbuf(ix  ,iy+1,iz  ) */
          tpt -= t1dimx;
          res += v4 * (*tpt);                /* tbuf(ix+1,iy  ,iz  ) */
          tpt --;
          res += (1-dy-dz+dydz-v4) * (*tpt); /* tbuf(ix  ,iy  ,iz  ) */
          *rbuf = (TYPE)_CONVERT_( res );
          continue;
        }
        /* here, we are sure we are on some border */
        tpt += (size_t)ix + (size_t)iy * (size_t)tdimx + (size_t)iz * (size_t)tdimxy;
        if ( (x < 0.0) || (ix == t1dimx) ) {
          if ( (y < 0.0) || (iy == t1dimy) ) {
            if ( (z < 0.0) || (iz == t1dimz) ) {
              *rbuf = *tpt;
              continue;
            }
            dz = z - iz;
            res  = (1-dz) * (*tpt); /* (1-dz)* tbuf(ix,iy,iz) */
            tpt += tdimxy;
            res += dz * (*tpt);     /* dz * tbuf(ix,iy,iz+1) */
            *rbuf = (TYPE)_CONVERT_( res );
            continue;
          }
          dy = y - iy;
          if ( (z < 0.0) || (iz == t1dimz) ) {
            res  = (1-dy) * (*tpt); /* (1-dy)* tbuf(ix,iy,iz) */
            tpt += tdimx;
            res += dy * (*tpt);     /* dy * tbuf(ix,iy+1,iz) */
            *rbuf = (TYPE)_CONVERT_( res );
            continue;
          }
          dz = z - iz;
          res = (1-dy)*(1-dz) * (*tpt); /* tbuf(ix,iy,iz) */
          tpt += tdimx;
          res += dy*(1-dz) * (*tpt);    /* tbuf(ix,iy+1,iz) */
          tpt += toffset2+1;
          res += (1-dy)*dz * (*tpt);    /* tbuf(ix,iy,iz+1) */
          tpt += tdimx;
          res += dy*dz * (*tpt);        /* tbuf(ix,iy+1,iz+1) */
          *rbuf = (TYPE)_CONVERT_( res );
          continue;
        }
        /* here we are sure that the border is either
           along the Y or the Z axis */
        dx = x - ix;
        if ( (y < 0.0) || (iy == t1dimy) ) {
          if ( (z < 0.0) || (iz == t1dimz) ) {
            res = (1-dx) * (*tpt); /* tbuf(ix,iy,iz) */
            tpt ++;
            res += dx * (*tpt);    /* tbuf(ix+1,iy,iz) */
            *rbuf = (TYPE)_CONVERT_( res );
            continue;
          }
          dz = z - iz;
          res = (1-dx)*(1-dz) * (*tpt); /* tbuf(ix,iy,iz) */
          tpt ++;
          res += dx*(1-dz) * (*tpt);    /* tbuf(ix+1,iy,iz) */
          tpt += tdimxy-1;
          res += (1-dx)*dz * (*tpt);    /* tbuf(ix,iy,iz+1) */
          tpt ++;
          res += dx*dz * (*tpt);        /* tbuf(ix+1,iy,iz+1) */
          *rbuf = (TYPE)_CONVERT_( res );
          continue;
        }
        /* here we are sure that the border is along the Z axis */
        dy = y - iy;
        res = (1-dx)*(1-dy) * (*tpt); /* tbuf(ix,iy,iz) */
        tpt ++;
        res += dx*(1-dy) * (*tpt);    /* tbuf(ix+1,iy,iz) */
        tpt += t1dimx;
        res += (1-dx)*dy * (*tpt);    /* tbuf(ix,iy+1,iz) */
        tpt ++;
        res += dx*dy * (*tpt);        /* tbuf(ix+1,iy+1,iz) */
        *rbuf = (TYPE)_CONVERT_( res );
      }
    }
  }
  chunk->ret = 1;
  return( (void*)NULL );
}





void Reech3DTriLin4x4_TYPE ( void* theBuf, /* buffer to be resampled */
                             int *theDim, /* dimensions of this buffer */
                             void* resBuf, /* result buffer */
                             int *resDim,  /* dimensions of this buffer */
                             double* mat   /* transformation matrix */
                             )
{
  char *proc = "Reech3DTriLin4x4_TYPE";
  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _LinearResamplingParam p;
  
  /* preparing parallelism
   */
  first = 0;
  last = (size_t)resDim[2] * (size_t)resDim[1] * (size_t)resDim[0] - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return;
  }
  
  p.theBuf = theBuf;
  p.theDim = theDim;
  p.resBuf = resBuf;
  p.resDim = resDim;
  p.mat = mat;
  p.gain = 1.0;
  p.bias = 0.0;
  
  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);
  
  /* processing
   */
  if ( processChunks( &_Reech3DTriLin4x4_TYPE, &chunks, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to resample image\n", proc );
    freeChunks( &chunks );
    return;
  }
  
  freeChunks( &chunks );
}





static void *_Reech3DTriLin4x4gb_TYPE ( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _LinearResamplingParam *p = (_LinearResamplingParam *)parameter;
  
  void* theBuf = p->theBuf;
  int *theDim = p->theDim;
  void* resBuf = p->resBuf;
  int *resDim = p->resDim;
  double* mat = p->mat; 
  float gain = p->gain;
  float bias = p->bias;
  
  size_t i, j, k;
  register int ix, iy, iz;
  register double x, y, z, dx, dy, dz, dxdy,dxdz,dydz,dxdydz;
  register double res;
  double v6, v5, v4;
  size_t rdimx=resDim[0], rdimy=resDim[1];
  int tdimx=theDim[0], tdimy=theDim[1], tdimz=theDim[2];
  int tdimxy=tdimx*tdimy;
  int toffset1=tdimxy+tdimx+1, toffset2=tdimxy-tdimx-1;
  register int t1dimx=tdimx-1, t1dimy=tdimy-1, t1dimz=tdimz-1;
  register double ddimx = (double)tdimx-0.5, ddimy = (double)tdimy-0.5;
  register double ddimz = (double)tdimz-0.5;
  register TYPE *tbuf = (TYPE*)theBuf;
  register TYPE *tpt;
  register TYPE *rbuf = (TYPE*)resBuf;
  register double b=bias;
  register double g=gain;
  
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;
  size_t iend, jend;
  
  k = kfirst = first / (rdimx*rdimy);
  j = jfirst = (first - kfirst*(rdimx*rdimy)) / rdimx;
  i = ifirst = (first - kfirst*(rdimx*rdimy) - jfirst*rdimx);
  
  klast = last / (rdimx*rdimy);
  jlast = (last - klast*(rdimx*rdimy)) / rdimx;
  ilast = (last - klast*(rdimx*rdimy) - jlast*rdimx);
  
  rbuf += first;
  
  for ( ; k<=klast; k++, j=0 ) {
    if ( _verbose_ > 1 )
      fprintf( stderr, "Processing slice %lu\r", k );
    jend = (k==klast) ? jlast+1 : rdimy;
    for ( ; j<jend; j++, i=0 ) {
      iend = (j==jlast && k==klast) ? ilast+1 : rdimx;
      for ( ; i<iend; i++, rbuf++ ) {
        /* computation of the corresponding point coordinates in theBuf */
        x = mat[0] * i +  mat[1] * j + mat[2] * k + mat[3];
        if ((x <= -0.5) || ( x >= ddimx)) { *rbuf = 0; continue; }
        y = mat[4] * i +  mat[5] * j + mat[6] * k + mat[7];
        if ((y <= -0.5) || ( y >= ddimy)) { *rbuf = 0; continue; }
        z = mat[8] * i +  mat[9] * j + mat[10] * k + mat[11];
        if ((z <= -0.5) || ( z >= ddimz)) { *rbuf = 0; continue; }

        /* here, the point lies on the borders or completely inside
           the image */
        ix = (int)x;
        iy = (int)y;
        iz = (int)z;
        tpt = (TYPE *)tbuf;

        /* are we on the border or not ? */
        if ( (x > 0.0) && (ix < t1dimx) &&
             (y > 0.0) && (iy < t1dimy) &&
             (z > 0.0) && (iz < t1dimz) ) {
          /* the corresponding point is in the box defined
             by (ix[+1],iy[+1],iz[+1]) */
          dx = x - ix;
          dy = y - iy;
          dz = z - iz;
          dxdy = dx*dy;
          dxdz = dx*dz;
          dydz = dy*dz;
          dxdydz = dxdy*dz;

          /* we have
             v[7]=dxdydz;            coefficient of tbuf(ix+1,iy+1,iz+1)
             v[6]=dxdz-dxdydz;       coefficient of tbuf(ix+1,iy,  iz+1)
             v[5]=dxdy-dxdydz;       coefficient of tbuf(ix+1,iy+1,iz  )
             v[4]=dx-dxdy-v[6];      coefficient of tbuf(ix+1,iy  ,iz  )
             v[3]=dydz-dxdydz;       coefficient of tbuf(ix  ,iy+1,iz+1)
             v[2]=dz-dydz-v[6];      coefficient of tbuf(ix  ,iy  ,iz+1)
             v[1]=dy-dydz-v[5];      coefficient of tbuf(ix  ,iy+1,iz  )
             v[0]=1-dy-dz+dydz-v[4]; coefficient of tbuf(ix  ,iy  ,iz  )
          */
          tpt += (size_t)ix + (size_t)iy * (size_t)tdimx + (size_t)iz * (size_t)tdimxy + (size_t)toffset1;
          v6 = dxdz-dxdydz;
          v5 = dxdy-dxdydz;
          v4 = dx-dxdy-v6;

          res = 0;
          res += dxdydz * (*tpt);            /* tbuf(ix+1,iy+1,iz+1) */
          tpt --;
          res += (dydz-dxdydz) * (*tpt);     /* tbuf(ix  ,iy+1,iz+1) */
          tpt -= t1dimx;
          res += v6 * (*tpt);                /* tbuf(ix+1  ,iy,  iz+1) */
          tpt --;
          res += (dz-dydz-v6) * (*tpt);      /* tbuf(ix  ,iy  ,iz+1) */
          tpt -= toffset2;
          res += v5 * (*tpt);                /* tbuf(ix+1,iy+1,iz  ) */
          tpt --;
          res += (dy-dydz-v5) * (*tpt);      /* tbuf(ix  ,iy+1,iz  ) */
          tpt -= t1dimx;
          res += v4 * (*tpt);                /* tbuf(ix+1,iy  ,iz  ) */
          tpt --;
          res += (1-dy-dz+dydz-v4) * (*tpt); /* tbuf(ix  ,iy  ,iz  ) */
          res = res * g + b;
          *rbuf = (TYPE)_CONVERT_( res );
          continue;
        }
        /* here, we are sure we are on some border */
        tpt += (size_t)ix + (size_t)iy * (size_t)tdimx + (size_t)iz * (size_t)tdimxy;
        if ( (x < 0.0) || (ix == t1dimx) ) {
          if ( (y < 0.0) || (iy == t1dimy) ) {
            if ( (z < 0.0) || (iz == t1dimz) ) {
              res = (double)(*tpt) * g + b;
              *rbuf = (TYPE)_CONVERT_( res );
              continue;
            }
            dz = z - iz;
            res  = (1-dz) * (*tpt); /* (1-dz)* tbuf(ix,iy,iz) */
            tpt += tdimxy;
            res += dz * (*tpt);     /* dz * tbuf(ix,iy,iz+1) */
            res = res * g + b;
            *rbuf = (TYPE)_CONVERT_( res );
            continue;
          }
          dy = y - iy;
          if ( (z < 0.0) || (iz == t1dimz) ) {
            res  = (1-dy) * (*tpt); /* (1-dy)* tbuf(ix,iy,iz) */
            tpt += tdimx;
            res += dy * (*tpt);     /* dy * tbuf(ix,iy+1,iz) */
            res = res * g + b;
            *rbuf = (TYPE)_CONVERT_( res );
            continue;
          }
          dz = z - iz;
          res = (1-dy)*(1-dz) * (*tpt); /* tbuf(ix,iy,iz) */
          tpt += tdimx;
          res += dy*(1-dz) * (*tpt);    /* tbuf(ix,iy+1,iz) */
          tpt += toffset2+1;
          res += (1-dy)*dz * (*tpt);    /* tbuf(ix,iy,iz+1) */
          tpt += tdimx;
          res += dy*dz * (*tpt);        /* tbuf(ix,iy+1,iz+1) */
          res = res * g + b;
          *rbuf = (TYPE)_CONVERT_( res );
          continue;
        }
        /* here we are sure that the border is either
           along the Y or the Z axis */
        dx = x - ix;
        if ( (y < 0.0) || (iy == t1dimy) ) {
          if ( (z < 0.0) || (iz == t1dimz) ) {
            res = (1-dx) * (*tpt); /* tbuf(ix,iy,iz) */
            tpt ++;
            res += dx * (*tpt);    /* tbuf(ix+1,iy,iz) */
            res = res * g + b;
            *rbuf = (TYPE)_CONVERT_( res );
            continue;
          }
          dz = z - iz;
          res = (1-dx)*(1-dz) * (*tpt); /* tbuf(ix,iy,iz) */
          tpt ++;
          res += dx*(1-dz) * (*tpt);    /* tbuf(ix+1,iy,iz) */
          tpt += tdimxy-1;
          res += (1-dx)*dz * (*tpt);    /* tbuf(ix,iy,iz+1) */
          tpt ++;
          res += dx*dz * (*tpt);        /* tbuf(ix+1,iy,iz+1) */
          res = res * g + b;
          *rbuf = (TYPE)_CONVERT_( res );
          continue;
        }
        /* here we are sure that the border is along the Z axis */
        dy = y - iy;
        res = (1-dx)*(1-dy) * (*tpt); /* tbuf(ix,iy,iz) */
        tpt ++;
        res += dx*(1-dy) * (*tpt);    /* tbuf(ix+1,iy,iz) */
        tpt += t1dimx;
        res += (1-dx)*dy * (*tpt);    /* tbuf(ix,iy+1,iz) */
        tpt ++;
        res += dx*dy * (*tpt);        /* tbuf(ix+1,iy+1,iz) */
        res = res * g + b;
        *rbuf = (TYPE)_CONVERT_( res );
      }
    }
  }
  chunk->ret = 1;
  return( (void*)NULL );
}





void Reech3DTriLin4x4gb_TYPE ( void* theBuf, /* buffer to be resampled */
                               int *theDim, /* dimensions of this buffer */
                               void* resBuf, /* result buffer */
                               int *resDim,  /* dimensions of this buffer */
                               double* mat,   /* transformation matrix */
                               float gain,
                               float bias )
{
  char *proc = "Reech3DTriLin4x4gb_TYPE";
  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _LinearResamplingParam p;
  
  /* preparing parallelism
   */
  first = 0;
  last = (size_t)resDim[2] * (size_t)resDim[1] * (size_t)resDim[0] - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return;
  }
  
  p.theBuf = theBuf;
  p.theDim = theDim;
  p.resBuf = resBuf;
  p.resDim = resDim;
  p.mat = mat;
  p.gain = gain;
  p.bias = bias;
  
  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);
  
  /* processing
   */
  if ( processChunks( &_Reech3DTriLin4x4gb_TYPE, &chunks, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to resample image\n", proc );
    freeChunks( &chunks );
    return;
  }
  
  freeChunks( &chunks );
}





static void *_Reech3DNearest4x4_TYPE ( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _LinearResamplingParam *p = (_LinearResamplingParam *)parameter;
  
  void* theBuf = p->theBuf;
  int *theDim = p->theDim;
  void* resBuf = p->resBuf;
  int *resDim = p->resDim;
  double* mat = p->mat; 
  
  size_t i, j, k;
  register int ix, iy, iz;
  register double x, y, z;
  size_t rdimx=resDim[0], rdimy=resDim[1];
  int tdimx=theDim[0], tdimy=theDim[1], tdimz=theDim[2];
  int tdimxy=tdimx*tdimy;
  register int t1dimx=tdimx-1, t1dimy=tdimy-1, t1dimz=tdimz-1;
  register TYPE *tbuf = (TYPE*)theBuf;
  register TYPE *rbuf = (TYPE*)resBuf;
  
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;
  size_t iend, jend;
  
  k = kfirst = first / (rdimx*rdimy);
  j = jfirst = (first - kfirst*(rdimx*rdimy)) / rdimx;
  i = ifirst = (first - kfirst*(rdimx*rdimy) - jfirst*rdimx);
  
  klast = last / (rdimx*rdimy);
  jlast = (last - klast*(rdimx*rdimy)) / rdimx;
  ilast = (last - klast*(rdimx*rdimy) - jlast*rdimx);
  
  rbuf += first;
  
  for ( ; k<=klast; k++, j=0 ) {
    if ( _verbose_ > 1 )
      fprintf( stderr, "Processing slice %lu\r", k );
    jend = (k==klast) ? jlast+1 : rdimy;
    for ( ; j<jend; j++, i=0 ) {
      iend = (j==jlast && k==klast) ? ilast+1 : rdimx;
      for ( ; i<iend; i++, rbuf++ ) {
        /* computation of the corresponding point coordinates in theBuf */
        x = mat[0] * i +  mat[1] * j + mat[2] * k + mat[3];
        ix = (int)(x+0.5);
        if (( x <= -0.5 ) || ( ix > t1dimx)) { *rbuf = 0; continue; }
        y = mat[4] * i +  mat[5] * j + mat[6] * k + mat[7];
        iy = (int)(y+0.5);
        if (( y <= -0.5 ) || ( iy > t1dimy)) { *rbuf = 0; continue; }
        z = mat[8] * i +  mat[9] * j + mat[10] * k + mat[11];
        iz = (int)(z+0.5);
        if (( z <= -0.5 ) || ( iz > t1dimz)) { *rbuf = 0; continue; }

        *rbuf = tbuf[ (size_t)ix + (size_t)iy * (size_t)tdimx + (size_t)iz * (size_t)tdimxy ];
      }
    }
  }
  chunk->ret = 1;
  return( (void*)NULL );
}





void Reech3DNearest4x4_TYPE ( void* theBuf, /* buffer to be resampled */
                              int *theDim,  /* dimensions of this buffer */
                              void* resBuf, /* result buffer */
                              int *resDim,  /* dimensions of this buffer */
                              double* mat   /* transformation matrix */
                              )
{
  char *proc = "Reech3DNearest4x4_TYPE";
  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _LinearResamplingParam p;
  
  /* preparing parallelism
   */
  first = 0;
  last = (size_t)resDim[2] * (size_t)resDim[1] * (size_t)resDim[0] - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return;
  }
  
  p.theBuf = theBuf;
  p.theDim = theDim;
  p.resBuf = resBuf;
  p.resDim = resDim;
  p.mat = mat;
  p.gain = 1.0;
  p.bias = 0.0;
  
  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);
  
  /* processing
   */
  if ( processChunks( &_Reech3DNearest4x4_TYPE, &chunks, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to resample image\n", proc );
    freeChunks( &chunks );
    return;
  }
  
  freeChunks( &chunks );
}





static void *_Reech2DTriLin4x4_TYPE ( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _LinearResamplingParam *p = (_LinearResamplingParam *)parameter;
  
  void* theBuf = p->theBuf;
  int *theDim = p->theDim;
  void* resBuf = p->resBuf;
  int *resDim = p->resDim;
  double* mat = p->mat; 
  
  size_t i, j, k;
  register int ix, iy;
  register double x, y, dx, dy, dxdy;
  register double res, v;
  size_t rdimx=resDim[0], rdimy=resDim[1];
  int tdimx=theDim[0], tdimy=theDim[1];
  int toffset=tdimx-1;
  register int t1dimx=tdimx-1, t1dimy=tdimy-1;
  register double ddimx = (double)tdimx-0.5, ddimy = (double)tdimy-0.5;
  register TYPE *tbuf = (TYPE*)theBuf;
  register TYPE *tpt;
  register TYPE *rbuf = (TYPE*)resBuf;
  
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;
  size_t iend, jend;
  
  k = kfirst = first / (rdimx*rdimy);
  j = jfirst = (first - kfirst*(rdimx*rdimy)) / rdimx;
  i = ifirst = (first - kfirst*(rdimx*rdimy) - jfirst*rdimx);
  
  klast = last / (rdimx*rdimy);
  jlast = (last - klast*(rdimx*rdimy)) / rdimx;
  ilast = (last - klast*(rdimx*rdimy) - jlast*rdimx);
  
  rbuf += first;
  
  for ( ; k<=klast; k++, j=0 ) {
    if ( _verbose_ > 1 )
      fprintf( stderr, "Processing slice %lu\r", k );
    /* tbuf represente le premier point du plan */
    tbuf  = (TYPE*)theBuf;
    tbuf += k*(tdimx * tdimy);
    jend = (k==klast) ? jlast+1 : rdimy;
    for ( ; j<jend; j++, i=0 ) {
      iend = (j==jlast && k==klast) ? ilast+1 : rdimx;
      for ( ; i<iend; i++, rbuf++ ) {
        /* computation of the corresponding point coordinates in theBuf */
        x = mat[0] * i +  mat[1] * j +              mat[3];
        if ((x <= -0.5) || ( x >= ddimx)) { *rbuf = 0; continue; }
        y = mat[4] * i +  mat[5] * j +              mat[7];
        if ((y <= -0.5) || ( y >= ddimy)) { *rbuf = 0; continue; }

        /* here, the point lies on the borders or completely inside
           the image */
        ix = (int)x;
        iy = (int)y;
        tpt = (TYPE *)tbuf;
        tpt += (size_t)ix + (size_t)iy * (size_t)tdimx;

        /* are we on the border or not ? */
        if ( (x > 0.0) && (ix < t1dimx) &&
             (y > 0.0) && (iy < t1dimy) ) {
          dx = x - ix;
          dy = y - iy;
          dxdy = dx*dy;
          /* we have
             v[5]=dxdy;         coefficient of tbuf(ix+1,iy+1)
             v[4]=dx-dxdy;      coefficient of tbuf(ix+1,iy  )
             v[1]=dy-dxdy;      coefficient of tbuf(ix  ,iy+1)
             v[0]=1-dx-dy+dxdy; coefficient of tbuf(ix  ,iy  )
          */
          v = dy-dxdy;
          res = 0;
          res += (1-dx-v) * (*tpt);  /* tbuf(ix  ,iy  ) */
          tpt ++;
          res += (dx-dxdy) * (*tpt); /* tbuf(ix+1,iy  ) */
          tpt += toffset;
          res += v * (*tpt);       /* tbuf(ix,iy+1  ) */
          tpt ++;
          res += dxdy * (*tpt);      /* tbuf(ix+1,iy+1) */
          *rbuf = (TYPE)_CONVERT_( res );
          continue;
        }

        /* here, we are sure we are on some border */
        if ( (x < 0.0) || (ix == t1dimx) ) {
          /* we just look at y */
          if ( (y < 0.0) || (iy == t1dimy) ) {
            *rbuf = *tpt;
            continue;
          }
          dy = y - iy;
          res  = (1-dy) * (*tpt); /* (1-dy)* tbuf(ix,iy) */
          tpt += tdimx;
          res += dy * (*tpt);     /* dy * tbuf(ix,iy+1) */
          *rbuf = (TYPE)_CONVERT_( res );
          continue;
        }
        dx = x - ix;
        res  = (1-dx) * (*tpt); /* (1-dx)* tbuf(ix,iy) */
        tpt ++;
        res += dx * (*tpt);     /* dx * tbuf(ix+1,iy) */
        *rbuf = (TYPE)_CONVERT_( res );
      }
    }
  }
  chunk->ret = 1;
  return( (void*)NULL );
}





void Reech2DTriLin4x4_TYPE ( void* theBuf, /* buffer to be resampled */
                             int *theDim, /* dimensions of this buffer */
                             void* resBuf, /* result buffer */
                             int *resDim,  /* dimensions of this buffer */
                             double* mat   /* transformation matrix */
                             )
{
  char *proc = "Reech2DTriLin4x4_TYPE";
  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _LinearResamplingParam p;
  
  /* preparing parallelism
   */
  first = 0;
  last = (size_t)resDim[2] * (size_t)resDim[1] * (size_t)resDim[0] - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return;
  }
  
  p.theBuf = theBuf;
  p.theDim = theDim;
  p.resBuf = resBuf;
  p.resDim = resDim;
  p.mat = mat;
  p.gain = 1.0;
  p.bias = 0.0;
  
  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);
  
  /* processing
   */
  if ( processChunks( &_Reech2DTriLin4x4_TYPE, &chunks, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to resample image\n", proc );
    freeChunks( &chunks );
    return;
  }
  
  freeChunks( &chunks );
}





static void *_Reech2DTriLin4x4gb_TYPE ( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _LinearResamplingParam *p = (_LinearResamplingParam *)parameter;
  
  void* theBuf = p->theBuf;
  int *theDim = p->theDim;
  void* resBuf = p->resBuf;
  int *resDim = p->resDim;
  double* mat = p->mat; 
  float gain = p->gain;
  float bias = p->bias;
  
  size_t i, j, k;
  register int ix, iy;
  register double x, y, dx, dy, dxdy;
  register double res, v;
  size_t rdimx=resDim[0], rdimy=resDim[1];
  int tdimx=theDim[0], tdimy=theDim[1];
  int toffset=tdimx-1;
  register int t1dimx=tdimx-1, t1dimy=tdimy-1;
  register double ddimx = (double)tdimx-0.5, ddimy = (double)tdimy-0.5;
  register TYPE *tbuf = (TYPE*)theBuf;
  register TYPE *tpt;
  register TYPE *rbuf = (TYPE*)resBuf;
  register double b=bias;
  register double g=gain;
  
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;
  size_t iend, jend;
  
  k = kfirst = first / (rdimx*rdimy);
  j = jfirst = (first - kfirst*(rdimx*rdimy)) / rdimx;
  i = ifirst = (first - kfirst*(rdimx*rdimy) - jfirst*rdimx);
  
  klast = last / (rdimx*rdimy);
  jlast = (last - klast*(rdimx*rdimy)) / rdimx;
  ilast = (last - klast*(rdimx*rdimy) - jlast*rdimx);
  
  rbuf += first;
  
  for ( ; k<=klast; k++, j=0 ) {
    if ( _verbose_ > 1 )
      fprintf( stderr, "Processing slice %lu\r", k );
    /* tbuf represente le premier point du plan */
    tbuf  = (TYPE*)theBuf;
    tbuf += k*(tdimx * tdimy);
    jend = (k==klast) ? jlast+1 : rdimy;
    for ( ; j<jend; j++, i=0 ) {
      iend = (j==jlast && k==klast) ? ilast+1 : rdimx;
      for ( ; i<iend; i++, rbuf++ ) {
        /* computation of the corresponding point coordinates in theBuf */
        x = mat[0] * i +  mat[1] * j +              mat[3];
        if ((x <= -0.5) || ( x >= ddimx)) { *rbuf = 0; continue; }
        y = mat[4] * i +  mat[5] * j +              mat[7];
        if ((y <= -0.5) || ( y >= ddimy)) { *rbuf = 0; continue; }

        /* here, the point lies on the borders or completely inside
           the image */
        ix = (int)x;
        iy = (int)y;
        tpt = (TYPE *)tbuf;
        tpt += (size_t)ix + (size_t)iy * (size_t)tdimx;

        /* are we on the border or not ? */
        if ( (x > 0.0) && (ix < t1dimx) &&
             (y > 0.0) && (iy < t1dimy) ) {
          dx = x - ix;
          dy = y - iy;
          dxdy = dx*dy;
          /* we have
             v[5]=dxdy;         coefficient of tbuf(ix+1,iy+1)
             v[4]=dx-dxdy;      coefficient of tbuf(ix+1,iy  )
             v[1]=dy-dxdy;      coefficient of tbuf(ix  ,iy+1)
             v[0]=1-dx-dy+dxdy; coefficient of tbuf(ix  ,iy  )
          */
          v = dy-dxdy;
          res = 0;
          res += (1-dx-v) * (*tpt);  /* tbuf(ix  ,iy  ) */
          tpt ++;
          res += (dx-dxdy) * (*tpt); /* tbuf(ix+1,iy  ) */
          tpt += toffset;
          res += v * (*tpt);       /* tbuf(ix,iy+1  ) */
          tpt ++;
          res += dxdy * (*tpt);      /* tbuf(ix+1,iy+1) */
          res = res * g + b;
          *rbuf = (TYPE)_CONVERT_( res );
          continue;
        }

        /* here, we are sure we are on some border */
        if ( (x < 0.0) || (ix == t1dimx) ) {
          /* we just look at y */
          if ( (y < 0.0) || (iy == t1dimy) ) {
            res = (double)(*tpt) * g + b;
            *rbuf = (TYPE)_CONVERT_( res );
            continue;
          }
          dy = y - iy;
          res  = (1-dy) * (*tpt); /* (1-dy)* tbuf(ix,iy) */
          tpt += tdimx;
          res += dy * (*tpt);     /* dy * tbuf(ix,iy+1) */
          res = (double)(*tpt) * g + b;
          *rbuf = (TYPE)_CONVERT_( res );
          continue;
        }
        dx = x - ix;
        res  = (1-dx) * (*tpt); /* (1-dx)* tbuf(ix,iy) */
        tpt ++;
        res += dx * (*tpt);     /* dx * tbuf(ix+1,iy) */
        res = res * g + b;
        *rbuf = (TYPE)_CONVERT_( res );
      }
    }
  }
  chunk->ret = 1;
  return( (void*)NULL );
}





void Reech2DTriLin4x4gb_TYPE ( void* theBuf, /* buffer to be resampled */
                               int *theDim, /* dimensions of this buffer */
                               void* resBuf, /* result buffer */
                               int *resDim,  /* dimensions of this buffer */
                               double* mat,   /* transformation matrix */
                               float gain,
                               float bias )
{
  char *proc = "Reech2DTriLin4x4gb_TYPE";
  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _LinearResamplingParam p;
  
  /* preparing parallelism
   */
  first = 0;
  last = (size_t)resDim[2] * (size_t)resDim[1] * (size_t)resDim[0] - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return;
  }
  
  p.theBuf = theBuf;
  p.theDim = theDim;
  p.resBuf = resBuf;
  p.resDim = resDim;
  p.mat = mat;
  p.gain = gain;
  p.bias = bias;
  
  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);
  
  /* processing
   */
  if ( processChunks( &_Reech2DTriLin4x4gb_TYPE, &chunks, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to resample image\n", proc );
    freeChunks( &chunks );
    return;
  }
  
  freeChunks( &chunks );
}





static void *_Reech2DNearest4x4_TYPE ( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _LinearResamplingParam *p = (_LinearResamplingParam *)parameter;
  
  void* theBuf = p->theBuf;
  int *theDim = p->theDim;
  void* resBuf = p->resBuf;
  int *resDim = p->resDim;
  double* mat = p->mat; 
  
  size_t i, j, k;
  register int ix, iy;
  register double x, y;
  size_t rdimx=resDim[0], rdimy=resDim[1];
  int tdimx=theDim[0], tdimy=theDim[1];
  register int t1dimx=tdimx-1, t1dimy=tdimy-1;
  register TYPE *tbuf = (TYPE*)theBuf;
  register TYPE *rbuf = (TYPE*)resBuf;
  
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;
  size_t iend, jend;
  
  k = kfirst = first / (rdimx*rdimy);
  j = jfirst = (first - kfirst*(rdimx*rdimy)) / rdimx;
  i = ifirst = (first - kfirst*(rdimx*rdimy) - jfirst*rdimx);
  
  klast = last / (rdimx*rdimy);
  jlast = (last - klast*(rdimx*rdimy)) / rdimx;
  ilast = (last - klast*(rdimx*rdimy) - jlast*rdimx);
  
  rbuf += first;
  
  for ( ; k<=klast; k++, j=0 ) {
    if ( _verbose_ > 1 )
      fprintf( stderr, "Processing slice %lu\r", k );
    /* tbuf represente le premier point du plan */
    tbuf  = (TYPE*)theBuf;
    tbuf += k*(tdimx * tdimy);
    jend = (k==klast) ? jlast+1 : rdimy;
    for ( ; j<jend; j++, i=0 ) {
      iend = (j==jlast && k==klast) ? ilast+1 : rdimx;
      for ( ; i<iend; i++, rbuf++ ) {
        /* computation of the corresponding point coordinates in theBuf */
        x = mat[0] * i +  mat[1] * j + mat[3];
        ix = (int)(x+0.5);
        if (( x <= -0.5 ) || ( ix > t1dimx)) { *rbuf = 0; continue; }
        y = mat[4] * i +  mat[5] * j + mat[7];
        iy = (int)(y+0.5);
        if (( y <= -0.5 ) || ( iy > t1dimy)) { *rbuf = 0; continue; }

        *rbuf = tbuf[ (size_t)ix + (size_t)iy * (size_t)tdimx ];
      }
    }
  }
  chunk->ret = 1;
  return( (void*)NULL );
}





void Reech2DNearest4x4_TYPE ( void* theBuf, /* buffer to be resampled */
                              int *theDim,  /* dimensions of this buffer */
                              void* resBuf, /* result buffer */
                              int *resDim,  /* dimensions of this buffer */
                              double* mat   /* transformation matrix */
                              )
{
  char *proc = "Reech2DNearest4x4_TYPE";
  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _LinearResamplingParam p;
  
  /* preparing parallelism
   */
  first = 0;
  last = (size_t)resDim[2] * (size_t)resDim[1] * (size_t)resDim[0] - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return;
  }
  
  p.theBuf = theBuf;
  p.theDim = theDim;
  p.resBuf = resBuf;
  p.resDim = resDim;
  p.mat = mat;
  p.gain = 1.0;
  p.bias = 0.0;
  
  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);
  
  /* processing
   */
  if ( processChunks( &_Reech2DNearest4x4_TYPE, &chunks, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to resample image\n", proc );
    freeChunks( &chunks );
    return;
  }
  
  freeChunks( &chunks );
}
