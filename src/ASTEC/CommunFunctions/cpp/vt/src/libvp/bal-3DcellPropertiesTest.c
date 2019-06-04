/*************************************************************************
 * bal_3DcellPropertiesTest.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Lun 26 nov 2018 11:42:23 CET
 *
 * ADDITIONS, CHANGES
 *
 */


/* random(), srandom(): Feature Test Macro Requirements for glibc
 * _SVID_SOURCE || _BSD_SOURCE || _XOPEN_SOURCE >= 500 || _XOPEN_SOURCE && _XOPEN_SOURCE_EXTENDED
 *
 * compilation with [gcc (GCC) 5.3.1 20151207 (Red Hat 5.3.1-2)] yields
 * "_BSD_SOURCE and _SVID_SOURCE are deprecated, use _DEFAULT_SOURCE"
 */
#define _XOPEN_SOURCE
#define _XOPEN_SOURCE_EXTENDED


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fcntl.h>
#include <sys/stat.h>

#ifndef WIN32
#include <unistd.h>
#endif

#include <math.h>

#include <vtmalloc.h>

#include <bal-transformation-cell.h>
#include <bal-transformation-tools.h>

#include <bal-3DcellPropertiesTest.h>


static int _verbose_ = 1;



/************************************************************
 *
 *
 *
 ************************************************************/


static char *_BaseName( char *p )
{
  int l;
  if ( p == (char*)NULL ) return( (char*)NULL );
  l = strlen( p ) - 1;
  while ( l >= 0 && p[l] != '/' ) l--;
  if ( l < 0 ) l = 0;
  if ( p[l] == '/' ) l++;
  return( &(p[l]) );
}





static int _ScilabOutputShapeTest( float *degree_angles, typeCellSequence *cellProperties, char *basename, float surface )
{
  char *proc = "_ScilabOutputShapeTest";
  float *angle = (float*)NULL;
  float *surf_error = (float*)NULL;
  float *surf_diff = (float*)NULL;
  float *volume = (float*)NULL;
  float surf_error_min, surf_error_max;
  float surf_diff_min = 0.0, surf_diff_max = 0.0;
  float surf_outer, surf_inner;
  int neighbor;
  float volume_min, volume_max;
  typeCell *cell, *c;
  int i, j, k, n;

  FILE *f;
  int fd;
  char filename[256];

  angle = (float*)vtmalloc( 8*cellProperties->n_data*sizeof(float), "x", proc );
  if ( angle == (float*)NULL ) {
      if ( _verbose_ )
          fprintf( stderr, "%s: allocation failed\n", proc );
      return( -1 );
  }

  /* erreur relative (pourcentage) sur les surfaces de contact
   * erreur = 100 * (mesure - (surface/8)) / (surface/8)
   */
  surf_error = (float*)vtmalloc( 8*cellProperties->n_data*sizeof(float), "surf_error", proc );
  if ( surf_error == (float*)NULL ) {
      vtfree( angle );
      if ( _verbose_ )
          fprintf( stderr, "%s: allocation failed\n", proc );
      return( -1 );
  }

  /*
   *
   */
  surf_diff = (float*)vtmalloc( 8*cellProperties->n_data*sizeof(float), "surf_diff", proc );
  if ( surf_diff == (float*)NULL ) {
      vtfree( surf_error );
      vtfree( angle );
      if ( _verbose_ )
          fprintf( stderr, "%s: allocation failed\n", proc );
      return( -1 );
  }

  /* volumes
   */
  volume = (float*)vtmalloc( cellProperties->n_data*sizeof(float), "volume", proc );
  if ( volume == (float*)NULL ) {
      vtfree( surf_diff );
      vtfree( surf_error );
      vtfree( angle );
      if ( _verbose_ )
          fprintf( stderr, "%s: allocation failed\n", proc );
      return( -1 );
  }

  surf_error_min = surf_error_max = 0.0;
  volume_min = volume_max = 0.0;

  /* surface errors
   */
  for ( n=0, i=0; i<cellProperties->n_data; i++ ) {
    if ( cellProperties->data[i].n_data != 11 ) {
        vtfree( volume );
        vtfree( surf_diff );
        vtfree( surf_error );
        vtfree( angle );
        if ( _verbose_ )
            fprintf( stderr, "%s: weird number of cells for test #%d\n", proc, i );
        return( -1 );
    }
    cell = &(cellProperties->data[i].data[2]);
    if ( cell->neighbors.n_data != 8 ) {
        vtfree( volume );
        vtfree( surf_diff );
        vtfree( surf_error );
        vtfree( angle );
        if ( _verbose_ )
            fprintf( stderr, "%s: weird number of neighbors for test #%d\n", proc, i );
        return( -1 );
    }

    /* volume
     */
    volume[i] = cell->npts;
    if ( i == 0 ) {
        volume_min = volume_max = volume[i];
    }
    else {
        if ( volume_min > volume[i] ) volume_min = volume[i];
        if ( volume_max < volume[i] ) volume_max = volume[i];
    }

    /* surfaces de contact
     */
    for ( j=0; j<8; j++, n++ ) {
        angle[n] = degree_angles[i];

        surf_error[n] = 100.0 * (cell->neighbors.data[j].surface - surface/8.0) / (surface/8.0);
        if ( n == 0 ) {
          surf_error_min = surf_error_max = surf_error[n];
        }
        else {
            if ( surf_error_min > surf_error[n] ) surf_error_min = surf_error[n];
            if ( surf_error_max < surf_error[n] ) surf_error_max = surf_error[n];
        }


        surf_outer = cell->neighbors.data[j].surface;
        neighbor = cell->neighbors.data[j].label;

        c = &(cellProperties->data[i].data[neighbor]);

        surf_inner = -1.0;
        for ( k=0; k<c->neighbors.n_data; k++ ) {
            if ( c->neighbors.data[k].label == 2 )
                surf_inner = c->neighbors.data[k].surface;
        }

        if ( surf_inner < 0.0 ) {
            surf_diff[n] = 0.0;
            if ( _verbose_ )
                fprintf( stderr, "%s: inner surface not found for contact (%d,%d)\n", proc, 2, neighbor );
        }
        else {
            surf_diff[n] = 100.0 * fabs( surf_outer - surf_inner ) / (surface/8.0);
        }
        if ( n == 0 ) {
          surf_diff_min = surf_diff_max = surf_diff[n];
        }
        else {
            if ( surf_diff_min > surf_diff[n] ) surf_diff_min = surf_diff[n];
            if ( surf_diff_max < surf_diff[n] ) surf_diff_max = surf_diff[n];
        }
    }

  }

  /* test
   */
  if ( 0 ) {
    for ( i=0; i<cellProperties->n_data; i++ ) {
      fprintf( stderr, "test #%3d: angle = %f, volume = %f\n", i, degree_angles[i], volume[i] );
      fprintf( stderr, "\t surf_angle =" );
      for ( j=8*i; j<8*(i+1); j++ )
          fprintf( stderr, " %f", angle[j] );
      fprintf( stderr, "\n" );
      fprintf( stderr, "\t surf_error =" );
      for ( j=8*i; j<8*(i+1); j++ )
          fprintf( stderr, " %f", surf_error[j] );
      fprintf( stderr, "\n" );
      fprintf( stderr, "\t surf_diff  =" );
      for ( j=8*i; j<8*(i+1); j++ )
          fprintf( stderr, " %f", surf_diff[j] );
      fprintf( stderr, "\n" );
    }
  }


  /* ecriture des fichiers
   */

  sprintf( filename, "%s.raw", basename );
  fd = open( filename, O_CREAT | O_TRUNC | O_WRONLY, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
  if ( fd == -1 ) {
    vtfree( volume );
    vtfree( surf_diff );
    vtfree( surf_error );
    vtfree( angle );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to open '%s' for writing\n", proc, filename );
    return( -1 );
  }

  sprintf( filename, "%s.sce", basename );
  f = fopen( filename, "w" );
  if ( f == (FILE*)NULL ) {
    close( fd );
    vtfree( volume );
    vtfree( surf_diff );
    vtfree( surf_error );
    vtfree( angle );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to open '%s' for writing\n", proc, filename );
    return( -1 );
  }

  /* scatter plot des errors absolues sur les surfaces
   */

  fprintf( f, "\n" );
  fprintf( f, "myfile=mopen('%s.raw','r');\n", _BaseName( basename ) );
  fprintf( f, "\n" );

  if ( write( fd, angle, 8*cellProperties->n_data*sizeof(float) ) == -1
       || write( fd, surf_error, 8*cellProperties->n_data*sizeof(float) ) == -1 ) {
      fclose( f );
      close( fd );
      vtfree( volume );
      vtfree( surf_diff );
      vtfree( surf_error );
      vtfree( angle );
      if ( _verbose_ )
          fprintf( stderr, "%s: error when writing raw data\n", proc );
      return( -1 );
  }

  fprintf( f, "\n" );
  fprintf( f, "////////////////////////////////////////////////////////\n" );
  fprintf( f, "// scatter plot\n" );

  fprintf( f, "\n" );
  fprintf( f, "// read data\n" );
  fprintf( f, "//\n" );
  fprintf( f, "surface_angles=mget( %d, 'f', myfile );\n", 8*cellProperties->n_data );
  fprintf( f, "surface_errors=mget( %d, 'f', myfile );\n", 8*cellProperties->n_data );

  fprintf( f, "\n" );
  fprintf( f, "// figure\n" );
  fprintf( f, "//\n" );

  fprintf( f, "figure;\n" );
  fprintf( f, "myfig=gcf();\n" );
  fprintf( f, "myfig.background = color(\"white\");\n" );

  fprintf( f, "// a=get(\"current_axes\");\n" );
  fprintf( f, "a=gca();\n" );
  fprintf( f, "set(a,\"auto_clear\",\"off\");\n" );
  fprintf( f, "// removing the trailing ';' allows to see all properties\n" );

  fprintf( f, "\n" );
  fprintf( f, "s = scatter(surface_angles, surface_errors);\n" );

  fprintf( f, "\n" );
  fprintf( f, "// a.data_bounds = [0,90;%f,%f];\n", surf_error_min, surf_error_max );
  fprintf( f, "a.font_size = 3;\n" );
  fprintf( f, "a.font_style = 8;\n" );
  fprintf( f, "\n" );
  fprintf( f, "a.title.text = \"Relative error (percentage) on contact surface\";\n" );
  fprintf( f, "a.x_label.text = \"angle (degrees)\";\n" );
  fprintf( f, "// a.y_label.text = \"Y Label\";\n" );
  fprintf( f, "// or \n" );
  fprintf( f, "// xtitle( \"Title\", \"X Label\", \"Y Label\" );\n" );
  fprintf( f, "a.title.font_size = 4;\n" );
  fprintf( f, "a.x_label.font_size = 3;\n" );
  fprintf( f, "// a.y_label.font_size = 3;\n" );
  fprintf( f, "\n" );
  fprintf( f, "// or \n" );
  fprintf( f, "// xlabel( 'X Label', 'fontsize', 4, 'fontname', 8 );\n" );
  fprintf( f, "// ylabel( 'Y Label', 'fontsize', 4, 'fontname', 8 );\n" );
  fprintf( f, "\n" );

  fprintf( f, "e = gce();\n" );
  fprintf( f, "// e.children(1).thickness = 3;\n" );
  fprintf( f, "e.children(1).foreground = 2;\n" );
  fprintf( f, "\n" );

  fprintf( f, "// xs2jpg(gcf(),'FIG%s_surferror.jpg');\n", _BaseName( basename )  );
  fprintf( f, "xs2png(gcf(),'FIG%s_surferror.png');\n", _BaseName( basename )  );
  fprintf( f, "\n" );

  fprintf( f, "// scatter plot\n" );
  fprintf( f, "////////////////////////////////////////////////////////\n" );
  fprintf( f, "\n" );

  /* scatter plot des differences des mesures absolues sur les surfaces
   */

  if ( write( fd, surf_diff, 8*cellProperties->n_data*sizeof(float) ) == -1 ) {
      fclose( f );
      close( fd );
      vtfree( volume );
      vtfree( surf_diff );
      vtfree( surf_error );
      vtfree( angle );
      if ( _verbose_ )
          fprintf( stderr, "%s: error when writing raw data\n", proc );
      return( -1 );
  }

  fprintf( f, "\n" );
  fprintf( f, "////////////////////////////////////////////////////////\n" );
  fprintf( f, "// scatter plot\n" );

  fprintf( f, "\n" );
  fprintf( f, "// read data\n" );
  fprintf( f, "//\n" );
  fprintf( f, "surface_diffs=mget( %d, 'f', myfile );\n", 8*cellProperties->n_data );

  fprintf( f, "\n" );
  fprintf( f, "// figure\n" );
  fprintf( f, "//\n" );

  fprintf( f, "figure;\n" );
  fprintf( f, "myfig=gcf();\n" );
  fprintf( f, "myfig.background = color(\"white\");\n" );

  fprintf( f, "// a=get(\"current_axes\");\n" );
  fprintf( f, "a=gca();\n" );
  fprintf( f, "set(a,\"auto_clear\",\"off\");\n" );
  fprintf( f, "// removing the trailing ';' allows to see all properties\n" );

  fprintf( f, "\n" );
  fprintf( f, "s = scatter(surface_angles, surface_diffs);\n" );

  fprintf( f, "\n" );
  fprintf( f, "// a.data_bounds = [0,90;%f,%f];\n", surf_diff_min, surf_diff_max );
  fprintf( f, "// a.data_bounds = [0,90;0,%f];\n", surf_diff_max );
  fprintf( f, "a.font_size = 3;\n" );
  fprintf( f, "a.font_style = 8;\n" );
  fprintf( f, "\n" );
  fprintf( f, "a.title.text = \"Relative error (percentage) on contact surface difference\";\n" );
  fprintf( f, "a.x_label.text = \"angle (degrees)\";\n" );
  fprintf( f, "// a.y_label.text = \"Y Label\";\n" );
  fprintf( f, "// or \n" );
  fprintf( f, "// xtitle( \"Title\", \"X Label\", \"Y Label\" );\n" );
  fprintf( f, "a.title.font_size = 4;\n" );
  fprintf( f, "a.x_label.font_size = 3;\n" );
  fprintf( f, "// a.y_label.font_size = 3;\n" );
  fprintf( f, "\n" );
  fprintf( f, "// or \n" );
  fprintf( f, "// xlabel( 'X Label', 'fontsize', 4, 'fontname', 8 );\n" );
  fprintf( f, "// ylabel( 'Y Label', 'fontsize', 4, 'fontname', 8 );\n" );
  fprintf( f, "\n" );

  fprintf( f, "e = gce();\n" );
  fprintf( f, "// e.children(1).thickness = 3;\n" );
  fprintf( f, "e.children(1).foreground = 2;\n" );
  fprintf( f, "\n" );

  fprintf( f, "// xs2jpg(gcf(),'FIG%s_surfdiff.jpg');\n", _BaseName( basename )  );
  fprintf( f, "xs2png(gcf(),'FIG%s_surfdiff.png');\n", _BaseName( basename )  );
  fprintf( f, "\n" );

  fprintf( f, "// scatter plot\n" );
  fprintf( f, "////////////////////////////////////////////////////////\n" );
  fprintf( f, "\n" );


  /* volumes
   */

  if ( write( fd, degree_angles, cellProperties->n_data*sizeof(float) ) == -1
       || write( fd, volume, cellProperties->n_data*sizeof(float) ) == -1 ) {
      fclose( f );
      close( fd );
      vtfree( volume );
      vtfree( surf_diff );
      vtfree( surf_error );
      vtfree( angle );
      if ( _verbose_ )
          fprintf( stderr, "%s: error when writing raw data\n", proc );
      return( -1 );
  }

  fprintf( f, "\n" );
  fprintf( f, "////////////////////////////////////////////////////////\n" );
  fprintf( f, "// volume\n" );

  fprintf( f, "\n" );
  fprintf( f, "// read data\n" );
  fprintf( f, "//\n" );
  fprintf( f, "angles=mget( %d, 'f', myfile );\n", cellProperties->n_data );
  fprintf( f, "volumes=mget( %d, 'f', myfile );\n", cellProperties->n_data );

  fprintf( f, "\n" );
  fprintf( f, "// figure\n" );
  fprintf( f, "//\n" );

  fprintf( f, "figure;\n" );
  fprintf( f, "myfig=gcf();\n" );
  fprintf( f, "myfig.background = color(\"white\");\n" );

  fprintf( f, "// a=get(\"current_axes\");\n" );
  fprintf( f, "a=gca();\n" );
  fprintf( f, "set(a,\"auto_clear\",\"off\");\n" );
  fprintf( f, "// removing the trailing ';' allows to see all properties\n" );

  fprintf( f, "\n" );
  fprintf( f, "s = scatter(angles, volumes);\n" );

  fprintf( f, "\n" );
  fprintf( f, "// a.data_bounds = [0,90;%f,%f];\n", volume_min, volume_max );
  fprintf( f, "a.font_size = 3;\n" );
  fprintf( f, "a.font_style = 8;\n" );
  fprintf( f, "\n" );
  fprintf( f, "a.title.text = \"Volume (voxels)\";\n" );
  fprintf( f, "a.x_label.text = \"angle (degrees)\";\n" );
  fprintf( f, "// a.y_label.text = \"Y Label\";\n" );
  fprintf( f, "// or \n" );
  fprintf( f, "// xtitle( \"Title\", \"X Label\", \"Y Label\" );\n" );
  fprintf( f, "a.title.font_size = 4;\n" );
  fprintf( f, "a.x_label.font_size = 3;\n" );
  fprintf( f, "// a.y_label.font_size = 3;\n" );
  fprintf( f, "\n" );
  fprintf( f, "// or \n" );
  fprintf( f, "// xlabel( 'X Label', 'fontsize', 4, 'fontname', 8 );\n" );
  fprintf( f, "// ylabel( 'Y Label', 'fontsize', 4, 'fontname', 8 );\n" );
  fprintf( f, "\n" );

  fprintf( f, "e = gce();\n" );
  fprintf( f, "// e.children(1).thickness = 3;\n" );
  fprintf( f, "e.children(1).foreground = 2;\n" );
  fprintf( f, "\n" );

  fprintf( f, "// xs2jpg(gcf(),'FIG%s_volume.jpg');\n", _BaseName( basename )  );
  fprintf( f, "xs2png(gcf(),'FIG%s_volume.png');\n", _BaseName( basename )  );
  fprintf( f, "\n" );

  fprintf( f, "// volumes\n" );
  fprintf( f, "////////////////////////////////////////////////////////\n" );
  fprintf( f, "\n" );


  fprintf( f, "\n" );
  fprintf( f, "mclose(myfile);\n" );
  fprintf( f, "\n" );



  /* close files,
   * release memory
   */

  fclose( f );
  close( fd );
  vtfree( volume );
  vtfree( surf_diff );
  vtfree( surf_error );
  vtfree( angle );
  return( 1 );
}





static int _compareneighbor ( const void * a, const void * b )
{
    typeSpatialNeighbor *da = (typeSpatialNeighbor*)a;
    typeSpatialNeighbor *db = (typeSpatialNeighbor*)b;

    if ( da->surface > db->surface ) return( -1 );
    if ( da->surface < db->surface ) return( 1 );
    return( 0 );
}





static int _ScilabOutputImageTest( float *sigma, typeCellSequence *cellProperties, char *basename )
{
  char *proc = "_ScilabOutputImageTest";

  int nTests = cellProperties->n_data;
  typeCell *firstCell, *cell;
  int label = -1;
  int c, n, i, j, k;
  float *surf = (float*)NULL;
  float *cumul = (float*)NULL;
  float *sum = (float*)NULL;
  float *scatter_sigma = (float*)NULL;
  float *scatter_diff = (float*)NULL;

  int nsurfaces;
  float *surf_percentage = (float*)NULL;
  float *surf_evolution = (float*)NULL;
  float csurf, csum, nsurf;

  FILE *f;
  int fd;
  char filename[256];

  /* ouverture des fichiers
   */

  sprintf( filename, "%s.raw", basename );
  fd = open( filename, O_CREAT | O_TRUNC | O_WRONLY, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
  if ( fd == -1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to open '%s' for writing\n", proc, filename );
    return( -1 );
  }

  sprintf( filename, "%s.sce", basename );
  f = fopen( filename, "w" );
  if ( f == (FILE*)NULL ) {
    close( fd );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to open '%s' for writing\n", proc, filename );
    return( -1 );
  }

  fprintf( f, "\n" );
  fprintf( f, "myfile=mopen('%s.raw','r');\n", _BaseName( basename ) );
  fprintf( f, "\n" );

  if ( write( fd, sigma, cellProperties->n_data*sizeof(float) ) == -1 ) {
    fclose( f );
    close( fd );
    if ( _verbose_ )
        fprintf( stderr, "%s: error when writing raw data\n", proc );
    return( -1 );
  }

  fprintf( f, "\n" );
  fprintf( f, "// read data\n" );
  fprintf( f, "//\n" );
  fprintf( f, "sigma=mget( %d, 'f', myfile );\n", cellProperties->n_data );
  fprintf( f, "\n" );


  /* loop on cells
   */

  for ( c=0; c<cellProperties->data[0].n_data; c++ ) {

    firstCell = &(cellProperties->data[0].data[c]);

    if ( firstCell->neighbors.n_data == 0 ) {
      fprintf( stderr, "%s: skip cell #%d, no neighbors\n", proc, c );
      continue;
    }

    qsort( firstCell->neighbors.data, firstCell->neighbors.n_data, sizeof(typeSpatialNeighbor), &_compareneighbor );

    surf = (float*)vtmalloc( nTests * (4 * firstCell->neighbors.n_data + 1) * sizeof(float), "surf", proc );

    if ( surf == (float*)NULL ) {
      fclose( f );
      close( fd );
      if ( _verbose_ )
        fprintf( stderr, "%s: allocation failed\n", proc );
      return( -1 );
    }

    cumul  = surf;
    cumul += nTests * firstCell->neighbors.n_data;
    sum    = cumul;
    sum   += nTests * firstCell->neighbors.n_data;
    scatter_sigma  = sum;
    scatter_sigma += nTests;
    scatter_diff   = scatter_sigma;
    scatter_diff  += nTests * firstCell->neighbors.n_data;

    /* process cell #c
     */

    for ( n=0; n<firstCell->neighbors.n_data; n++ ) {
      label = firstCell->neighbors.data[n].label;
      surf[ n * nTests ] = firstCell->neighbors.data[n].surface;

      for ( i=1; i<nTests; i++ ) {
        cell = &(cellProperties->data[i].data[c]);
        surf[ n * nTests + i ] = 0.0;
        for ( j=0; j<cell->neighbors.n_data; j++ ) {
          if ( label == cell->neighbors.data[j].label ) {
              surf[ n * nTests + i ] = cell->neighbors.data[j].surface;
          }
        }
      }
    }

    /* surface totale
     */
    for ( i=0; i<nTests; i++ ) {
      sum[i] = 0;
      for ( n=0; n<firstCell->neighbors.n_data; n++ )
        sum[i] += surf[ n * nTests + i ];
    }

    /* prorata surface totale
     */
    for ( i=0; i<nTests; i++ ) {
      for ( n=0; n<firstCell->neighbors.n_data; n++ ) {
        cumul[ n * nTests + i ] = 0.0;
        for ( j=0; j<=n; j++ )
          cumul[ n * nTests + i ] += surf[ j * nTests + i ];
        cumul[ n * nTests + i ] *= 100.0 / sum[i];
      }
    }

    for ( i=0; i<nTests; i++ ) {
      for ( n=0; n<firstCell->neighbors.n_data; n++ ) {
        scatter_sigma[n * nTests + i] = sigma[i];
        if ( i == 0 ) {
          scatter_diff[n * nTests + i] = 0.0;
        }
        else {
          scatter_diff[n * nTests + i] = 100.0 * fabs(surf[ n * nTests + i ] - surf[ n * nTests + (i-1) ]) / surf[ n * nTests + i-1 ];
        }
      }
    }


    if ( write( fd, surf, cellProperties->n_data*firstCell->neighbors.n_data*sizeof(float) ) == -1 ) {
      vtfree( surf );
      fclose( f );
      close( fd );
      if ( _verbose_ )
          fprintf( stderr, "%s: error when writing raw data\n", proc );
      return( -1 );
    }

    if ( write( fd, cumul, cellProperties->n_data*firstCell->neighbors.n_data*sizeof(float) ) == -1 ) {
      vtfree( surf );
      fclose( f );
      close( fd );
      if ( _verbose_ )
          fprintf( stderr, "%s: error when writing raw data\n", proc );
      return( -1 );
    }

    if ( write( fd, scatter_sigma, cellProperties->n_data*firstCell->neighbors.n_data*sizeof(float) ) == -1 ) {
      vtfree( surf );
      fclose( f );
      close( fd );
      if ( _verbose_ )
          fprintf( stderr, "%s: error when writing raw data\n", proc );
      return( -1 );
    }

    if ( write( fd, scatter_diff, cellProperties->n_data*firstCell->neighbors.n_data*sizeof(float) ) == -1 ) {
      vtfree( surf );
      fclose( f );
      close( fd );
      if ( _verbose_ )
          fprintf( stderr, "%s: error when writing raw data\n", proc );
      return( -1 );
    }


    vtfree( surf );


    fprintf( f, "\n" );
    fprintf( f, "////////////////////////////////////////////////////////\n" );
    fprintf( f, "// read data\n" );
    fprintf( f, "//\n" );
    for ( n=0; n<firstCell->neighbors.n_data; n++ ) {
      label = firstCell->neighbors.data[n].label;
      fprintf( f, "surface_contact_%03d_%03d=mget( %d, 'f', myfile );\n",
               c, label, cellProperties->n_data );
    }

    for ( n=0; n<firstCell->neighbors.n_data; n++ ) {
      label = firstCell->neighbors.data[n].label;
      fprintf( f, "surface_cumul_%03d_%03d=mget( %d, 'f', myfile );\n",
               c, label, cellProperties->n_data );
    }

    fprintf( f, "scatter_sigma_%03d=mget( %d, 'f', myfile );\n",
             c, cellProperties->n_data * firstCell->neighbors.n_data );
    fprintf( f, "scatter_diff_%03d=mget( %d, 'f', myfile );\n",
             c, cellProperties->n_data * firstCell->neighbors.n_data );





    fprintf( f, "\n" );
    fprintf( f, "// figure\n" );
    fprintf( f, "//\n" );

    fprintf( f, "figure;\n" );
    fprintf( f, "myfig=gcf();\n" );
    fprintf( f, "myfig.background = color(\"white\");\n" );

    fprintf( f, "// a=get(\"current_axes\");\n" );
    fprintf( f, "a=gca();\n" );
    fprintf( f, "set(a,\"auto_clear\",\"off\");\n" );
    fprintf( f, "// removing the trailing ';' allows to see all properties\n" );

    fprintf( f, "\n" );
    fprintf( f, "plot(sigma', [");
    for ( n=0; n<firstCell->neighbors.n_data; n++ ) {
      label = firstCell->neighbors.data[n].label;
      fprintf( f, "surface_contact_%03d_%03d'", c, label );
      if ( n<firstCell->neighbors.n_data-1 ) fprintf( f, ", ");
    }
    fprintf( f, "], \"thickness\", 2);\n" );

    fprintf( f, "\n" );
    fprintf( f, "// a.data_bounds = [%f,x;%f,x];\n", sigma[0], sigma[nTests-1] );
    fprintf( f, "a.font_size = 3;\n" );
    fprintf( f, "a.font_style = 8;\n" );
    fprintf( f, "\n" );
    fprintf( f, "a.title.text = \"Contact surfaces for cell #%d\";\n", c );
    fprintf( f, "a.x_label.text = \"sigma\";\n" );
    fprintf( f, "// a.y_label.text = \"Y Label\";\n" );
    fprintf( f, "// or \n" );
    fprintf( f, "// xtitle( \"Title\", \"X Label\", \"Y Label\" );\n" );
    fprintf( f, "a.title.font_size = 4;\n" );
    fprintf( f, "a.x_label.font_size = 3;\n" );
    fprintf( f, "// a.y_label.font_size = 3;\n" );
    fprintf( f, "\n" );
    fprintf( f, "// or \n" );
    fprintf( f, "// xlabel( 'X Label', 'fontsize', 4, 'fontname', 8 );\n" );
    fprintf( f, "// ylabel( 'Y Label', 'fontsize', 4, 'fontname', 8 );\n" );
    fprintf( f, "\n" );

    fprintf( f, "e = gce();\n" );
    fprintf( f, "// e.children(1).thickness = 3;\n" );
    fprintf( f, "e.children(1).foreground = 2;\n" );
    fprintf( f, "\n" );

    fprintf( f, "// xs2jpg(gcf(),'FIG%s_surface%03d.jpg');\n", _BaseName( basename ), c );
    fprintf( f, "// xs2png(gcf(),'FIG%s_surface%03d.png');\n", _BaseName( basename ), c );
    fprintf( f, "\n" );
    fprintf( f, "legend([");
    for ( n=0; n<firstCell->neighbors.n_data; n++ ) {
      label = firstCell->neighbors.data[n].label;
      fprintf( f, "'%d'", label );
      if ( n<firstCell->neighbors.n_data-1 ) fprintf( f, ";");
    }
    fprintf( f, "], -1);\n" );
    fprintf( f, "xs2png(gcf(),'FIG%s_surface%03d.png');\n", _BaseName( basename ), c );
    fprintf( f, "xdel( gcf().figure_id );\n" );




    fprintf( f, "\n" );
    fprintf( f, "// figure\n" );
    fprintf( f, "//\n" );

    fprintf( f, "figure;\n" );
    fprintf( f, "myfig=gcf();\n" );
    fprintf( f, "myfig.background = color(\"white\");\n" );

    fprintf( f, "// a=get(\"current_axes\");\n" );
    fprintf( f, "a=gca();\n" );
    fprintf( f, "set(a,\"auto_clear\",\"off\");\n" );
    fprintf( f, "// removing the trailing ';' allows to see all properties\n" );

    fprintf( f, "\n" );
    fprintf( f, "plot(sigma', [");
    for ( n=0; n<firstCell->neighbors.n_data; n++ ) {
      label = firstCell->neighbors.data[n].label;
      fprintf( f, "surface_cumul_%03d_%03d'", c, label );
      if ( n<firstCell->neighbors.n_data-1 ) fprintf( f, ", ");
    }
    fprintf( f, "], \"thickness\", 2);\n" );

    fprintf( f, "\n" );
    fprintf( f, "a.data_bounds = [%f,0;%f,100];\n", sigma[0], sigma[nTests-1] );
    fprintf( f, "a.font_size = 3;\n" );
    fprintf( f, "a.font_style = 8;\n" );
    fprintf( f, "\n" );
    fprintf( f, "a.title.text = \"Cumulative contact surfaces (percentage) for cell #%d\";\n", c );
    fprintf( f, "a.x_label.text = \"sigma\";\n" );
    fprintf( f, "// a.y_label.text = \"Y Label\";\n" );
    fprintf( f, "// or \n" );
    fprintf( f, "// xtitle( \"Title\", \"X Label\", \"Y Label\" );\n" );
    fprintf( f, "a.title.font_size = 4;\n" );
    fprintf( f, "a.x_label.font_size = 3;\n" );
    fprintf( f, "// a.y_label.font_size = 3;\n" );
    fprintf( f, "\n" );
    fprintf( f, "// or \n" );
    fprintf( f, "// xlabel( 'X Label', 'fontsize', 4, 'fontname', 8 );\n" );
    fprintf( f, "// ylabel( 'Y Label', 'fontsize', 4, 'fontname', 8 );\n" );
    fprintf( f, "\n" );

    fprintf( f, "e = gce();\n" );
    fprintf( f, "// e.children(1).thickness = 3;\n" );
    fprintf( f, "e.children(1).foreground = 2;\n" );
    fprintf( f, "\n" );

    fprintf( f, "// xs2jpg(gcf(),'FIG%s_cumulative%03d.jpg');\n", _BaseName( basename ), c );
    fprintf( f, "// xs2png(gcf(),'FIG%s_cumulative%03d.png');\n", _BaseName( basename ), c );
    fprintf( f, "\n" );
    fprintf( f, "legend([");
    for ( n=0; n<firstCell->neighbors.n_data; n++ ) {
      label = firstCell->neighbors.data[n].label;
      fprintf( f, "'%d'", label );
      if ( n<firstCell->neighbors.n_data-1 ) fprintf( f, ";");
    }
    fprintf( f, "], -1);\n" );
    fprintf( f, "xs2png(gcf(),'FIG%s_surfacecumul%03d.png');\n", _BaseName( basename ), c );
    fprintf( f, "xdel( gcf().figure_id );\n" );



    fprintf( f, "\n" );
    fprintf( f, "// figure\n" );
    fprintf( f, "//\n" );

    fprintf( f, "\n" );
    fprintf( f, "if 0 == 1 then;\n" );
    fprintf( f, "\n" );

    fprintf( f, "figure;\n" );
    fprintf( f, "myfig=gcf();\n" );
    fprintf( f, "myfig.background = color(\"white\");\n" );

    fprintf( f, "// a=get(\"current_axes\");\n" );
    fprintf( f, "a=gca();\n" );
    fprintf( f, "set(a,\"auto_clear\",\"off\");\n" );
    fprintf( f, "// removing the trailing ';' allows to see all properties\n" );

    fprintf( f, "\n" );
    fprintf( f, "s = scatter(scatter_sigma_%03d, scatter_diff_%03d);\n", c, c );

    fprintf( f, "\n" );
    fprintf( f, "// a.data_bounds = [0,90;x,x];\n" );
    fprintf( f, "a.font_size = 3;\n" );
    fprintf( f, "a.font_style = 8;\n" );
    fprintf( f, "\n" );
    fprintf( f, "a.title.text = \"Surface contact evolution (percentage)\";\n" );
    fprintf( f, "a.x_label.text = \"sigma\";\n" );
    fprintf( f, "// a.y_label.text = \"Y Label\";\n" );
    fprintf( f, "// or \n" );
    fprintf( f, "// xtitle( \"Title\", \"X Label\", \"Y Label\" );\n" );
    fprintf( f, "a.title.font_size = 4;\n" );
    fprintf( f, "a.x_label.font_size = 3;\n" );
    fprintf( f, "// a.y_label.font_size = 3;\n" );
    fprintf( f, "\n" );
    fprintf( f, "// or \n" );
    fprintf( f, "// xlabel( 'X Label', 'fontsize', 4, 'fontname', 8 );\n" );
    fprintf( f, "// ylabel( 'Y Label', 'fontsize', 4, 'fontname', 8 );\n" );
    fprintf( f, "\n" );

    fprintf( f, "e = gce();\n" );
    fprintf( f, "// e.children(1).thickness = 3;\n" );
    fprintf( f, "e.children(1).foreground = 2;\n" );
    fprintf( f, "\n" );

    fprintf( f, "// xs2jpg(gcf(),'FIG%s_surfacediff%03d.jpg');\n", _BaseName( basename ), c  );
    fprintf( f, "xs2png(gcf(),'FIG%s_surfacediff%03d.png');\n", _BaseName( basename ), c  );
    fprintf( f, "xdel( gcf().figure_id );\n" );

    fprintf( f, "\n" );
    fprintf( f, "end\n" );
    fprintf( f, "\n" );

    fprintf( f, "\n" );



  }

  fprintf( f, "\n" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );



  for ( i=0; i<nTests-1; i++ ) {

    nsurfaces = 0;
    for ( c=0; c<cellProperties->data[i].n_data; c++ ) {
      if ( cellProperties->data[i].data[c].neighbors.n_data > 0 )
        nsurfaces += cellProperties->data[i].data[c].neighbors.n_data;
    }

    surf_percentage = (float*)vtmalloc( 2 * nsurfaces * sizeof(float), "surf_percentage", proc );
    if ( surf_percentage == (float*)NULL ) {
      fclose( f );
      close( fd );
      if ( _verbose_ )
        fprintf( stderr, "%s: allocation failed\n", proc );
      return( -1 );
    }
    surf_evolution = surf_percentage;
    surf_evolution += nsurfaces;

    for ( k=0, c=0; c<cellProperties->data[i].n_data; c++ ) {

      firstCell = &(cellProperties->data[i].data[c]);
      for ( csum=0.0, n=0; n<firstCell->neighbors.n_data; n++ )
          csum += firstCell->neighbors.data[n].surface;

      for ( n=0; n<firstCell->neighbors.n_data; n++, k++ ) {
          csurf = firstCell->neighbors.data[n].surface;
          nsurf = 0.0;
          cell = &(cellProperties->data[i+1].data[c]);
          for ( j=0; j<cell->neighbors.n_data; j++ ) {
            if ( firstCell->neighbors.data[n].label == cell->neighbors.data[j].label )
              nsurf = cell->neighbors.data[j].surface;
          }
          surf_percentage[k] = 100.0 * csurf/csum;
          surf_evolution[k] = 100.0 * fabs(nsurf-csurf) / csurf;
      }
    }

    if ( write( fd, surf_percentage, nsurfaces*sizeof(float) ) == -1 ) {
      fclose( f );
      close( fd );
      if ( _verbose_ )
          fprintf( stderr, "%s: error when writing raw data\n", proc );
      return( -1 );
    }

    if ( write( fd, surf_evolution, nsurfaces*sizeof(float) ) == -1 ) {
      fclose( f );
      close( fd );
      if ( _verbose_ )
          fprintf( stderr, "%s: error when writing raw data\n", proc );
      return( -1 );
    }

    fprintf( f, "surf_percentage_%03d=mget( %d, 'f', myfile );\n",
             i, nsurfaces );
    fprintf( f, "surf_evolution_%03d=mget( %d, 'f', myfile );\n",
             i, nsurfaces );


    fprintf( f, "\n" );
    fprintf( f, "// figure\n" );
    fprintf( f, "//\n" );

    fprintf( f, "figure;\n" );
    fprintf( f, "myfig=gcf();\n" );
    fprintf( f, "myfig.background = color(\"white\");\n" );

    fprintf( f, "// a=get(\"current_axes\");\n" );
    fprintf( f, "a=gca();\n" );
    fprintf( f, "set(a,\"auto_clear\",\"off\");\n" );
    fprintf( f, "// removing the trailing ';' allows to see all properties\n" );

    fprintf( f, "\n" );
    fprintf( f, "s = scatter(surf_percentage_%03d, surf_evolution_%03d);\n", i, i );

    fprintf( f, "\n" );
    fprintf( f, "//a.data_bounds = [0,0;x,x];\n" );
    fprintf( f, "a.font_size = 3;\n" );
    fprintf( f, "a.font_style = 8;\n" );
    fprintf( f, "\n" );
    fprintf( f, "a.title.text = \"Surf. evolution from %f to %f / Surface contact (percentages)\";\n", sigma[i], sigma[i+1] );
    fprintf( f, "// a.x_label.text = \"sigma\";\n" );
    fprintf( f, "// a.y_label.text = \"Y Label\";\n" );
    fprintf( f, "// or \n" );
    fprintf( f, "// xtitle( \"Title\", \"X Label\", \"Y Label\" );\n" );
    fprintf( f, "a.title.font_size = 4;\n" );
    fprintf( f, "a.x_label.font_size = 3;\n" );
    fprintf( f, "// a.y_label.font_size = 3;\n" );
    fprintf( f, "\n" );
    fprintf( f, "// or \n" );
    fprintf( f, "// xlabel( 'X Label', 'fontsize', 4, 'fontname', 8 );\n" );
    fprintf( f, "// ylabel( 'Y Label', 'fontsize', 4, 'fontname', 8 );\n" );
    fprintf( f, "\n" );

    fprintf( f, "e = gce();\n" );
    fprintf( f, "// e.children(1).thickness = 3;\n" );
    fprintf( f, "e.children(1).foreground = 2;\n" );
    fprintf( f, "\n" );

    fprintf( f, "xs2png(gcf(),'FIG%s_scatterevolution%03d.png');\n", _BaseName( basename ), i+1  );
    fprintf( f, "xdel( gcf().figure_id );\n" );
    fprintf( f, "\n" );



    vtfree( surf_percentage );
  }




  fprintf( f, "\n" );
  fprintf( f, "mclose(myfile);\n" );
  fprintf( f, "\n" );



  /* close files,
   */

  fclose( f );
  close( fd );

  return( 1 );
}




/************************************************************
 *
 *
 *
 ************************************************************/

/*
 * int dl,                                      # cube de cote l, sphere de diametre l+1
 * int nTests,                                  #
 * float sigma,                                 # sigma > 0 => cell-based resampling
 * enumShape shape,
 * enumSurfaceEstimation surfaceEstimationType, # estimation de surface
 * char *basename
 */

int BAL_3DcellPropertiesShapeTest( typeCellSequence *cellProperties,
                              int dl,
                              float sigma,
                              enumShape shape,
                              enumSurfaceEstimation surfaceEstimationType,
                              char *basename )
{
  char *proc = "BAL_3DcellPropertiesShapeTest";

  bal_image theIm, trsfIm, *ptrIm;
  char imageName[256];
  u8 ***theBuf = (u8 ***)NULL;

  bal_transformation theTrsf, *ptrTrsf;
  char trsfName[256];

  int nTests = cellProperties->n_data;
  int l, x, y, z, a;
  float center[3];
  float dx, dy, dz, d, dmax;

  double radian_angle;
  float *degree_angles = (float*)NULL;
  float surface;

  /*-----------------------------------
   * image sans transformation
   */

  l = 2*dl;

  if ( BAL_AllocFullImage( &theIm, (char*)NULL, 2*l, 2*l, 2*l, 1, 1.0, 1.0, 1.0, UCHAR ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: template image allocation error\n", proc );
      return( -1 );
  }

  center[0] = (theIm.ncols-1.0)/2.0;
  center[1] = (theIm.nrows-1)/2.0;
  center[2] = (theIm.nplanes-1)/2.0;

  switch( shape ) {
  default :
      BAL_FreeImage( &theIm );
      if ( _verbose_ )
        fprintf( stderr, "%s: such shape not handled yet\n", proc );
      return( -1 );
  case _SQUARE_ :
      /* carre de cote l (en voxel)
       * surface = 6 * l * l
       */
      surface = 6 * l * l;
      theBuf = (u8 ***)theIm.array;
      for ( z=0; z<(int)theIm.nplanes; z++ )
      for ( y=0; y<(int)theIm.nrows; y++ )
      for ( x=0; x<(int)theIm.ncols; x++ ) {
        if ( z<dl || 3*dl<=z || y<dl || 3*dl<=y || x<dl || 3*dl<=x ) {
            if ( z < l ) {
              if ( y < l ) {
                if ( x < l )
                  theBuf[z][y][x] = 3;
                else
                  theBuf[z][y][x] = 4;
              }
              else {
                if ( x < l )
                  theBuf[z][y][x] = 5;
                else
                  theBuf[z][y][x] = 6;
              }
            }
            else {
              if ( y < l ) {
                if ( x < l )
                  theBuf[z][y][x] = 7;
                else
                  theBuf[z][y][x] = 8;
              }
              else {
                if ( x < l )
                  theBuf[z][y][x] = 9;
                else
                  theBuf[z][y][x] = 10;
              }
            }
        }
        else
            theBuf[z][y][x] = 2;
      }
      break;
  case _DISK_ :
      /* sphere de diametre l+1 (en voxel)
       * surface = 4 * pi * r * r = pi * (l+1) * (l+1)
       */
      surface = 3.141592653 * (l+1) * (l+1);
      theBuf = (u8 ***)theIm.array;
      dmax = 10;
      for ( z=0; z<(int)theIm.nplanes; z++ )
      for ( y=0; y<(int)theIm.nrows; y++ )
      for ( x=0; x<(int)theIm.ncols; x++ ) {
        dx = (float)x - center[0];
        dy = (float)y - center[1];
        dz = (float)z - center[2];
        d = sqrt( dx*dx + dy*dy + dz*dz ) * 1.0/(float)(l+1);
        if ( d > 0.5 ) {
            if ( z < l ) {
              if ( y < l ) {
                if ( x < l )
                  theBuf[z][y][x] = 3;
                else
                  theBuf[z][y][x] = 4;
              }
              else {
                if ( x < l )
                  theBuf[z][y][x] = 5;
                else
                  theBuf[z][y][x] = 6;
              }
            }
            else {
              if ( y < l ) {
                if ( x < l )
                  theBuf[z][y][x] = 7;
                else
                  theBuf[z][y][x] = 8;
              }
              else {
                if ( x < l )
                  theBuf[z][y][x] = 9;
                else
                  theBuf[z][y][x] = 10;
              }
            }
            if ( dmax > d ) dmax = d;
        }
        else
            theBuf[z][y][x] = 2;
      }
      break;
  }

  /* transformation
   */
  BAL_InitTransformation( &theTrsf );

  if ( nTests > 1 ) {
    if ( BAL_AllocImageFromImage( &trsfIm, (char*)NULL, &theIm, UCHAR ) != 1 ) {
        BAL_FreeImage( &theIm );
        if ( _verbose_ )
          fprintf( stderr, "%s: transformed image allocation error\n", proc );
        return( -1 );
    }
    if ( BAL_AllocTransformation( &theTrsf, RIGID_3D, (bal_image*)NULL ) != 1 ) {
      BAL_FreeImage( &trsfIm );
      BAL_FreeImage( &theIm );
      if ( _verbose_ )
        fprintf( stderr, "%s: transformation allocation error\n", proc );
      return( -1 );
    }
  }

  degree_angles = (float*)vtmalloc( nTests*sizeof(float), "degree_angles", proc );
  if ( degree_angles == (float*)NULL ) {
      if ( nTests > 1 ) BAL_FreeTransformation( &theTrsf );
      if ( nTests > 1 ) BAL_FreeImage( &trsfIm );
      BAL_FreeImage( &theIm );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate angle array\n", proc );
      return( -1 );
  }


  /*--------------------------------------
   * boucle sur les transformations/tests
   */

  for ( a=0; a<nTests; a++ ) {

    if ( _verbose_ >=1 )
      fprintf( stderr, "test #%3d\n", a );

    /* transformation
     */
    radian_angle = 0.0;
    if ( a == 0 ) {
        degree_angles[a] = 0.0;
        ptrTrsf = (bal_transformation*)NULL;
    }
    else {

        if ( BAL_SetTransformationToRandom( &theTrsf ) != 1 ) {
            vtfree( degree_angles );
            if ( nTests > 1 ) BAL_FreeTransformation( &theTrsf );
            if ( nTests > 1 ) BAL_FreeImage( &trsfIm );
            BAL_FreeImage( &theIm );
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to generate random transformation\n", proc );
            return( -1 );
        }

        if ( BAL_RotationAngle( &theTrsf, &radian_angle ) != 1 ) {
            vtfree( degree_angles );
            if ( nTests > 1 ) BAL_FreeTransformation( &theTrsf );
            if ( nTests > 1 ) BAL_FreeImage( &trsfIm );
            BAL_FreeImage( &theIm );
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to compute rotation angle\n", proc );
            return( -1 );
        }
        degree_angles[a] = radian_angle*180.0/3.141592653;

        theTrsf.mat.m[ 3] = center[0] - theTrsf.mat.m[0] * center[0] - theTrsf.mat.m[1] * center[1] - theTrsf.mat.m[ 2] * center[2];
        theTrsf.mat.m[ 7] = center[1] - theTrsf.mat.m[4] * center[0] - theTrsf.mat.m[5] * center[1] - theTrsf.mat.m[ 6] * center[2];
        theTrsf.mat.m[11] = center[2] - theTrsf.mat.m[8] * center[0] - theTrsf.mat.m[9] * center[1] - theTrsf.mat.m[10] * center[2];

        ptrTrsf = &theTrsf;
    }


    /* resampling
     */
    if ( sigma > 0 ) {
        if ( a == 0 ) {
            if ( _verbose_ >= 2 )
              fprintf( stderr, "%s: cell-based resampling\n", proc );
        }
        if ( BAL_ResampleCellImage( &theIm, &trsfIm, (bal_image*)NULL,
                                    ptrTrsf, sigma ) != 1 ) {
            vtfree( degree_angles );
            if ( nTests > 1 ) BAL_FreeTransformation( &theTrsf );
            if ( nTests > 1 ) BAL_FreeImage( &trsfIm );
            BAL_FreeImage( &theIm );
            if ( _verbose_ )
              fprintf( stderr, "%s: error when resampling cell-image '%s'\n", proc, imageName );
            return( -1 );
        }
        ptrIm = &trsfIm;
    }
    else {
        if ( a == 0 ) {
            if ( _verbose_ >= 2 )
              fprintf( stderr, "%s: nearest neighbor resampling\n", proc );
            ptrIm = &theIm;
        }
        else {
            if ( BAL_ResampleImage( &theIm, &trsfIm, ptrTrsf, NEAREST ) != 1 ) {
                vtfree( degree_angles );
              if ( nTests > 1 ) BAL_FreeTransformation( &theTrsf );
              if ( nTests > 1 ) BAL_FreeImage( &trsfIm );
              BAL_FreeImage( &theIm );
              if ( _verbose_ )
                fprintf( stderr, "%s: unable to compute resampling\n", proc );
              return( -1 );
            }
            ptrIm = &trsfIm;
        }
    }


    /* some writing
     */

    switch( shape ) {
    default :
        vtfree( degree_angles );
        if ( nTests > 1 ) BAL_FreeTransformation( &theTrsf );
        if ( nTests > 1 ) BAL_FreeImage( &trsfIm );
        BAL_FreeImage( &theIm );
        if ( _verbose_ )
          fprintf( stderr, "%s: such shape not handled yet\n", proc );
        return( -1 );
    case _SQUARE_ :
        if ( sigma > 0.0 )
            sprintf( imageName, "cube_cbr_%03d_l%06d_%05.1f.mha", a, l, degree_angles[a] );
        else
            sprintf( imageName, "cube_%03d_l%06d_%05.1f.mha", a, l, degree_angles[a] );
        break;
    case _DISK_ :
        if ( sigma > 0.0 )
            sprintf( imageName, "sphere_cbr_%03d_d%06d_%05.1f.mha", a, l+1, degree_angles[a] );
        else
            sprintf( imageName, "sphere_%03d_d%06d_%05.1f.mha", a, l+1, degree_angles[a] );
        break;
    }
    sprintf( trsfName, "trsf_%03d_%05.1f.trsf", a, degree_angles[a] );

    if ( 1 ) {
        if ( BAL_WriteImage( ptrIm, imageName ) != 1 ) {
          vtfree( degree_angles );
          if ( nTests > 1 ) BAL_FreeTransformation( &theTrsf );
          if ( nTests > 1 ) BAL_FreeImage( &trsfIm );
          BAL_FreeImage( &theIm );
          if ( _verbose_ )
            fprintf( stderr, "%s: error when writing image '%s'\n", proc, imageName );
          return( -1 );
        }
    }

    if ( 0 ) {
        if ( BAL_WriteTransformation( &theTrsf, trsfName ) != 1 ) {
            vtfree( degree_angles );
            if ( nTests > 1 ) BAL_FreeTransformation( &theTrsf );
            if ( nTests > 1 ) BAL_FreeImage( &trsfIm );
            BAL_FreeImage( &theIm );
            if ( _verbose_ )
              fprintf( stderr, "%s: error when writing transformation '%s'\n", proc, trsfName );
            return( -1 );
        }
    }

    /* computations
     */

    if ( BAL_CellImageProperties( ptrIm, &(cellProperties->data[a]),
                                  (typePropertyList*)NULL, surfaceEstimationType ) != 1 ) {
        vtfree( degree_angles );
        if ( nTests > 1 ) BAL_FreeTransformation( &theTrsf );
        if ( nTests > 1 ) BAL_FreeImage( &trsfIm );
        BAL_FreeImage( &theIm );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when computing contact for test #%d\n", proc, a );
        return( -1 );
    }

  }

  if ( nTests > 1 ) BAL_FreeTransformation( &theTrsf );
  if ( nTests > 1 ) BAL_FreeImage( &trsfIm );
  BAL_FreeImage( &theIm );

  /*
   * boucle sur les transformations/tests
   --------------------------------------*/

  if ( basename != (char*)NULL ) {
      if ( _ScilabOutputShapeTest( degree_angles, cellProperties, basename, surface ) != 1 ) {
          vtfree( degree_angles );
          if ( _verbose_ )
            fprintf( stderr, "%s: error when writing output '%s'\n", proc, basename );
          return( -1 );
      }
  }

  vtfree( degree_angles );

  return( 1 );
}





/************************************************************
 *
 *
 *
 ************************************************************/



int BAL_3DcellPropertiesImageTest( typeCellSequence *cellProperties,
                                   bal_image *theIm,
                                   float minSigma, float maxSigma,
                              enumSurfaceEstimation surfaceEstimationType,
                              char *basename )
{
  char *proc = "BAL_3DcellPropertiesImageTest";
  bal_image trsfIm;
  float *sigma;
  int nTests = cellProperties->n_data;
  int i;

  if ( BAL_AllocImageFromImage( &trsfIm, (char*)NULL, theIm, theIm->type ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: transformed image allocation error\n", proc );
      return( -1 );
  }

  sigma = (float*)vtmalloc( nTests*sizeof(float), "sigma", proc );
  if ( sigma == (float*)NULL ) {
      BAL_FreeImage( &trsfIm );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate angle array\n", proc );
      return( -1 );
  }

  for ( i=0; i<nTests; i++ ) {

    sigma[i] = minSigma + (float)i/(float)(nTests-1) * (maxSigma - minSigma);
    if ( _verbose_ )
      fprintf( stderr, " ... sigma = %f\n", sigma[i] );

    /* resampling
     */
    if ( sigma[i] > 0.0 ) {
        if ( BAL_ResampleCellImage( theIm, &trsfIm, (bal_image*)NULL,
                                    (bal_transformation*)NULL, sigma[i] ) != 1 ) {
            vtfree( sigma );
            BAL_FreeImage( &trsfIm );
            if ( _verbose_ )
              fprintf( stderr, "%s: error when resampling cell-image\n", proc );
            return( -1 );
        }
    }
    else {
        if ( BAL_ResampleImage( theIm, &trsfIm, (bal_transformation*)NULL, NEAREST ) != 1 ) {
          vtfree( sigma );
          BAL_FreeImage( &trsfIm );
          if ( _verbose_ )
            fprintf( stderr, "%s: error when resampling\n", proc );
          return( -1 );
        }
    }

    /* properties
     */
    if ( BAL_CellImageProperties( &trsfIm, &(cellProperties->data[i]),
                                  (typePropertyList*)NULL, surfaceEstimationType ) != 1 ) {
        vtfree( sigma );
        BAL_FreeImage( &trsfIm );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when computing contact for test #%d\n", proc, i );
        return( -1 );
    }
  }

  BAL_FreeImage( &trsfIm );

  /* scilab output
   */

  if ( basename != (char*)NULL ) {
      if ( _ScilabOutputImageTest( sigma, cellProperties, basename ) != 1 ) {
          vtfree( sigma );
          if ( _verbose_ )
            fprintf( stderr, "%s: error when writing output '%s'\n", proc, basename );
          return( -1 );
      }
  }

  vtfree( sigma );

  return( 1 );
}




