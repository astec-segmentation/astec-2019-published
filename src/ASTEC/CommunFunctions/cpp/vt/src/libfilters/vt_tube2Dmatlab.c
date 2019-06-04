/*************************************************************************
 * vt_tube2Dmatlab.c -
 *
 * $Id: vt_tube2Dmatlab.c,v 1.11 2006/05/16 09:33:34 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Tue Jul 11 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <unistd.h>

#include <vt_tube2D.h>
#include <vt_tube2Dmatlab.h>
#include <transfo.h>




















void VT_2DDrawImage( vt_image *theIm,
		     int slice,
		     int rawDataFileDesc,
		     FILE *commandFile )
{
  char *proc = "2DDrawImage";



  switch( theIm->type ) {
  default :
    fprintf( stderr, "%s: does not handle such image type\n", proc );
    return;
  case UCHAR :
  case USHORT :
  case FLOAT :
    fprintf( commandFile, "\n" );
    fprintf( commandFile, "\n" );
    fprintf( commandFile, "\n" );
    fprintf( commandFile, "%% image = %s\n", theIm->name );
  }

  
  switch( theIm->type ) {
  default :
    fprintf( stderr, "%s: does not handle such image type\n", proc );
    return;
  case UCHAR :
    if ( write( rawDataFileDesc, &((u8***)theIm->array)[slice][0][0], 
		theIm->dim.x * theIm->dim.y * sizeof( u8 ) ) == -1 ) {
      fprintf( stderr, "%s: error when writing\n", proc );
    }
    fprintf( commandFile, "[RAWIMAGE, DIMRAWIMAGE] = fread( fid, [%lu,%lu], 'uint%lu' );\n", 
	     theIm->dim.x, theIm->dim.y, 8*sizeof( u8 ) );
    break;
  case USHORT :
    if ( write( rawDataFileDesc, &((u16***)theIm->array)[slice][0][0], 
		theIm->dim.x * theIm->dim.y * sizeof( u16 ) ) == -1 ) {
      fprintf( stderr, "%s: error when writing\n", proc );
    }
    fprintf( commandFile, "[RAWIMAGE, DIMRAWIMAGE] = fread( fid, [%lu,%lu], 'uint%lu' );\n", 
	     theIm->dim.x, theIm->dim.y, 8*sizeof( u16 ) );
    break;
  case FLOAT :
    if ( write( rawDataFileDesc, &((r32***)theIm->array)[slice][0][0], 
		theIm->dim.x * theIm->dim.y * sizeof( r32 ) ) == -1 ) {
      fprintf( stderr, "%s: error when writing\n",proc  );
    }
    fprintf( commandFile, "[RAWIMAGE, DIMRAWIMAGE] = fread( fid, [%lu,%lu], 'float%lu' );\n", 
	     theIm->dim.x, theIm->dim.y, 8*sizeof( r32 ) );
    break;
  }

  /* matlab lit les colonnes d'abord,
     il faut donc transposer l'image
  */



  fprintf( commandFile, "imagesc( RAWIMAGE' );\n" );
  fprintf( commandFile, "colormap( gray );\n" );
  fprintf( commandFile, "axis([1 %lu 1 %lu]);\n",
	   theIm->dim.x, theIm->dim.y );
  fprintf( commandFile, "axis equal;\n" );
  fprintf( commandFile, "axis ij;\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );
}













void VT_3DDrawImage( vt_image *theIm,
		     int rawDataFileDesc,
		     FILE *commandFile )
{
  char *proc = "VT_3DDrawImage";


  switch( theIm->type ) {
  default :
    fprintf( stderr, "%s: does not handle such image type\n", proc );
    return;
  case UCHAR :
  case USHORT :
  case SSHORT :
  case FLOAT :
    fprintf( commandFile, "\n" );
    fprintf( commandFile, "\n" );
    fprintf( commandFile, "\n" );
    fprintf( commandFile, "%% image = %s\n", theIm->name );
  }

  switch( theIm->type ) {
  default :
    fprintf( stderr, "%s: does not handle such image type\n", proc );
    return;
  case UCHAR :
    fprintf( commandFile, "IMAGE = uint8( zeros([%lu %lu %lu]) );\n", 
	     theIm->dim.y, theIm->dim.x, theIm->dim.z );
    break;
  case USHORT :
    fprintf( commandFile, "IMAGE = uint16( zeros([%lu %lu %lu]) );\n", 
	     theIm->dim.y, theIm->dim.x, theIm->dim.z );
    break;
  case SSHORT :
    fprintf( commandFile, "IMAGE = int16( zeros([%lu %lu %lu]) );\n", 
	     theIm->dim.y, theIm->dim.x, theIm->dim.z );
    break;
  case FLOAT :
    fprintf( commandFile, "IMAGE = single( zeros([%lu %lu %lu]) );\n", 
	     theIm->dim.y, theIm->dim.x, theIm->dim.z );
    break;
  }

  
  
  switch( theIm->type ) {
  default :
    fprintf( stderr, "%s: does not handle such image type\n", proc );
    return;
  case UCHAR :
    if ( write( rawDataFileDesc, theIm->buf, theIm->dim.x * theIm->dim.y * theIm->dim.z * sizeof( u8 ) ) == -1 ) {
      fprintf( stderr, "%s: error when writing\n", proc );
    }
    fprintf( commandFile, "[RAWIMAGE, DIMRAWIMAGE] = fread( fid, [%lu,%lu], 'uint%lu' );\n", 
	     theIm->dim.x, theIm->dim.y * theIm->dim.z, 8*sizeof( u8 ) );
    break;
  case USHORT :
    if ( write( rawDataFileDesc, theIm->buf, theIm->dim.x * theIm->dim.y * theIm->dim.z * sizeof( u16 ) ) == -1 ) {
      fprintf( stderr, "%s: error when writing\n", proc );
    }
    fprintf( commandFile, "[RAWIMAGE, DIMRAWIMAGE] = fread( fid, [%lu,%lu], 'uint%lu' );\n", 
	     theIm->dim.x, theIm->dim.y * theIm->dim.z, 8*sizeof( u16 ) );
    break;
  case SSHORT :
    if ( write( rawDataFileDesc, theIm->buf, theIm->dim.x * theIm->dim.y * theIm->dim.z * sizeof( u16 ) ) == -1 ) {
      fprintf( stderr, "%s: error when writing\n", proc );
    }
    fprintf( commandFile, "[RAWIMAGE, DIMRAWIMAGE] = fread( fid, [%lu,%lu], 'uint%lu' );\n", 
	     theIm->dim.x, theIm->dim.y * theIm->dim.z, 8*sizeof( s16 ) );
    break;
  case FLOAT :
    if ( write( rawDataFileDesc, theIm->buf, theIm->dim.x * theIm->dim.y * theIm->dim.z * sizeof( r32 ) ) == -1 ) {
      fprintf( stderr, "%s: error when writing\n", proc );
    }
    fprintf( commandFile, "[RAWIMAGE, DIMRAWIMAGE] = fread( fid, [%lu,%lu], 'float%lu' );\n", 
	     theIm->dim.x, theIm->dim.y * theIm->dim.z, 8*sizeof( r32 ) );
    break;

  }
  
  fprintf( commandFile, "for z = 1:%lu\n", theIm->dim.z );
  fprintf( commandFile, "for y = 1:%lu\n", theIm->dim.y );
  fprintf( commandFile, "for x = 1:%lu\n", theIm->dim.x );
  fprintf( commandFile, "    IMAGE(y,x,z) = RAWIMAGE(x,(z-1)*%lu+y);\n", theIm->dim.y );
  fprintf( commandFile, "end\n" );
  fprintf( commandFile, "end\n" );
  fprintf( commandFile, "end\n" );
  fprintf( commandFile, "clear RAWIMAGE x y z;\n" );
  fprintf( commandFile, "\n" ); 
  fprintf( commandFile, "\n" ); 
  fprintf( commandFile, "\n" ); 
  fprintf( commandFile, "xlabel( 'X' );\n" );
  fprintf( commandFile, "ylabel( 'Y' );\n" );
  fprintf( commandFile, "zlabel( 'Z' );\n" );
  fprintf( commandFile, "%% slice en Z\n" );
  fprintf( commandFile, "hz = slice(double(IMAGE),[],[],[1] );\n" );
  fprintf( commandFile, "set( hz, 'FaceColor', 'interp', 'EdgeColor', 'none' );\n" );
  fprintf( commandFile, "%% slice en Y\n" );
  fprintf( commandFile, "hy = slice(double(IMAGE),[],[%lu],[] );\n", theIm->dim.y );
  fprintf( commandFile, "set( hy, 'FaceColor', 'interp', 'EdgeColor', 'none' );\n" );
  fprintf( commandFile, "%% slice en X\n" );
  fprintf( commandFile, "hx = slice(double(IMAGE),[%lu],[],[] );\n", theIm->dim.x );
  fprintf( commandFile, "set( hx, 'FaceColor', 'interp', 'EdgeColor', 'none' );\n" );
 
  fprintf( commandFile, "colormap( gray );\n" );

  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "shading interp;\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "%% p = patch( isosurface(IMAGE) );\n" );
  fprintf( commandFile, "%% p = patch( isosurface(IMAGE, 50) );\n" );
  fprintf( commandFile, "%% isonormals( IMAGE, p );\n" );
  fprintf( commandFile, "%% set(p, 'FaceColor', 'red', 'EdgeColor', 'none');\n" );
  fprintf( commandFile, "%% set(p, 'FaceColor', 'none', 'EdgeColor', 'red', 'EdgeLighting', 'phong');\n" );

  fprintf( commandFile, "daspect([1 1 1]);\n" );
  fprintf( commandFile, "view(3);\n" );
  fprintf( commandFile, "axis([1 %lu 1 %lu 1 %lu]);\n", theIm->dim.x, theIm->dim.y, theIm->dim.z );
  fprintf( commandFile, "%% camlight ;\n" );
  fprintf( commandFile, "%% lighting phong;\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );
}













void VT_2DDrawWeightedVectors( vt_image *imWeight,
			       vt_image *imTheta,
			       int rawDataFileDesc,
			       FILE *commandFile )
{
  char *proc = "VT_2DDrawWeightedVectors";
  int x, y;
  float ***theR = (float ***)imWeight->array;
  float ***theTHETA = (float ***)imTheta->array;
  
  float *cy, *cx = (float*)NULL;
  double max = 0;


  cx = (float*)malloc( 2*imWeight->dim.x * imWeight->dim.y * sizeof( float ) );
  cy = cx;
  cy += imWeight->dim.x * imWeight->dim.y;
  
  max = fabs( theR[0][0][0] );
  for ( y=0; y<(int)imWeight->dim.y; y++ )
  for ( x=0; x<(int)imWeight->dim.x; x++ ) {
    if ( max < fabs( theR[0][y][x] ) ) max = fabs( theR[0][y][x] );
  }
  

  for ( y=0; y<(int)imWeight->dim.y; y++ )
  for ( x=0; x<(int)imWeight->dim.x; x++ ) {
    cx[y*imWeight->dim.x+x] = fabs( theR[0][y][x] )/max * cos( theTHETA[0][y][x] );
    cy[y*imWeight->dim.x+x] = fabs( theR[0][y][x] )/max * sin( theTHETA[0][y][x] );
  }
  
  if ( write( rawDataFileDesc, cx, 2*imWeight->dim.x * imWeight->dim.y * sizeof( float ) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }
  
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "%% matlab transpose x et y par rapport a moi;\n" );
  fprintf( commandFile, "%% le maximum etait de %f\n", max );
  fprintf( commandFile, "[XORTHOEIGEN, NXORTHOEIGEN] = fread( fid, [%lu,%lu], 'float%lu' );\n", 
	   imWeight->dim.x, imWeight->dim.y, 8*sizeof( float ) );
  fprintf( commandFile, "[YORTHOEIGEN, NYORTHOEIGEN] = fread( fid, [%lu,%lu], 'float%lu' );\n", 
	   imWeight->dim.x, imWeight->dim.y, 8*sizeof( float ) );
  fprintf( commandFile, "quiver( XORTHOEIGEN', YORTHOEIGEN', 1.25, 'r' );\n" );
  fprintf( commandFile, "axis ij;\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );
  
}












void VT_Old3DDrawWeightedVectors( vt_image *imWeight,
			       vt_image *imTheta,
			       vt_image *imPhi,
			       int rawDataFileDesc,
			       FILE *commandFile )
{
  char *proc = "VT_Old3DDrawWeightedVectors";
  int x, y, z;
  float ***theR = (float ***)imWeight->array;
  float ***theTHETA = (float ***)imTheta->array;
  float ***thePHI = (float ***)imPhi->array;
  
  float *cz, *cy, *cx = (float*)NULL;
  double max = 0;
  double v[3];


  cx = (float*)malloc( 3*imWeight->dim.x * imWeight->dim.y * imWeight->dim.z * sizeof( float ) );
  cz = cy = cx;
  cy += imWeight->dim.x * imWeight->dim.y * imWeight->dim.z;
  cz += 2*imWeight->dim.x * imWeight->dim.y * imWeight->dim.z;
   
  max = fabs( theR[0][0][0] );
  for ( z=0; z<(int)imWeight->dim.z; z++ )
  for ( y=0; y<(int)imWeight->dim.y; y++ )
  for ( x=0; x<(int)imWeight->dim.x; x++ ) {
    if ( max < fabs( theR[z][y][x] ) ) max = fabs( theR[z][y][x] );
  }
  
  cx = (float*)malloc( 3*imWeight->dim.x * imWeight->dim.y * imWeight->dim.z * sizeof( float ) );
  cz = cy = cx;
  cy += imWeight->dim.x * imWeight->dim.y * imWeight->dim.z;
  cz += 2*imWeight->dim.x * imWeight->dim.y * imWeight->dim.z;


  for ( z=0; z<(int)imWeight->dim.z; z++ )
  for ( y=0; y<(int)imWeight->dim.y; y++ )
  for ( x=0; x<(int)imWeight->dim.x; x++ ) {
    SphericalAnglesToUnitVector( (double)theTHETA[z][y][x], (double)thePHI[z][y][x], v );
    cx[z*imWeight->dim.x*imWeight->dim.y+y*imWeight->dim.x+x] =
      fabs( theR[z][y][x] )/max * v[0];
    cy[z*imWeight->dim.x*imWeight->dim.y+y*imWeight->dim.x+x] =
      fabs( theR[z][y][x] )/max * v[1];
    cz[z*imWeight->dim.x*imWeight->dim.y+y*imWeight->dim.x+x] =
      fabs( theR[z][y][x] )/max * v[2];
  }
  
  if ( write( rawDataFileDesc, cx, 3*imWeight->dim.x*imWeight->dim.y*imWeight->dim.z * sizeof( float ) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }
  
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "[XRAWIMAGE, DIMXRAWIMAGE] = fread( fid, [%lu,%lu], 'float%lu' );\n", 
	   imWeight->dim.x, imWeight->dim.y*imWeight->dim.z, 8*sizeof( r32 ) );
  fprintf( commandFile, "[YRAWIMAGE, DIMYRAWIMAGE] = fread( fid, [%lu,%lu], 'float%lu' );\n", 
	   imWeight->dim.x, imWeight->dim.y*imWeight->dim.z, 8*sizeof( r32 ) );
  fprintf( commandFile, "[ZRAWIMAGE, DIMZRAWIMAGE] = fread( fid, [%lu,%lu], 'float%lu' );\n", 
	   imWeight->dim.x, imWeight->dim.y*imWeight->dim.z, 8*sizeof( r32 ) );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "for z = 1:%lu\n", imWeight->dim.z );
  fprintf( commandFile, "for y = 1:%lu\n", imWeight->dim.y );
  fprintf( commandFile, "for x = 1:%lu\n", imWeight->dim.x );
  fprintf( commandFile, "    XIMAGE(y,x,z) = XRAWIMAGE(x,(z-1)*%lu+y);\n", imWeight->dim.y );
  fprintf( commandFile, "end\n" );
  fprintf( commandFile, "end\n" );
  fprintf( commandFile, "end\n" );
  fprintf( commandFile, "clear XRAWIMAGE x y z;\n" );
  fprintf( commandFile, "for z = 1:%lu\n", imWeight->dim.z );
  fprintf( commandFile, "for y = 1:%lu\n", imWeight->dim.y );
  fprintf( commandFile, "for x = 1:%lu\n", imWeight->dim.x );
  fprintf( commandFile, "    YIMAGE(y,x,z) = YRAWIMAGE(x,(z-1)*%lu+y);\n", imWeight->dim.y );
  fprintf( commandFile, "end\n" );
  fprintf( commandFile, "end\n" );
  fprintf( commandFile, "end\n" );
  fprintf( commandFile, "clear YRAWIMAGE x y z;\n" );
  fprintf( commandFile, "for z = 1:%lu\n", imWeight->dim.z );
  fprintf( commandFile, "for y = 1:%lu\n", imWeight->dim.y );
  fprintf( commandFile, "for x = 1:%lu\n", imWeight->dim.x );
  fprintf( commandFile, "    ZIMAGE(y,x,z) = ZRAWIMAGE(x,(z-1)*%lu+y);\n", imWeight->dim.y );
  fprintf( commandFile, "end\n" );
  fprintf( commandFile, "end\n" );
  fprintf( commandFile, "end\n" );
  fprintf( commandFile, "clear ZRAWIMAGE x y z;\n" );
  fprintf( commandFile, "\n" ); 
  fprintf( commandFile, "\n" ); 
  fprintf( commandFile, "[X,Y,Z] = meshgrid(1:%lu,1:%lu,1:%lu);", imWeight->dim.x, 
	   imWeight->dim.y, imWeight->dim.z ); 
  fprintf( commandFile, "quiver3(X,Y,Z,XIMAGE, YIMAGE, ZIMAGE, 1.25, 'r' );\n" );
  fprintf( commandFile, "%% clear X Y Z;\n");
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );
  
}










void VT_3DDrawWeightedVectors( vt_image *imWeight,
			       vt_image *imTheta,
			       vt_image *imPhi,
			       int rawDataFileDesc,
			       FILE *commandFile,
			       double threshold )
{
  char *proc = "VT_3DDrawWeightedVectors";
  int x, y, z;
  float ***theR = (float ***)imWeight->array;
  float ***theTHETA = (float ***)imTheta->array;
  float ***thePHI = (float ***)imPhi->array;
  
  float *cz, *cy, *cx = (float*)NULL;
  int *iz, *iy, *ix = (int*)NULL;
  int n, nb = 0;
  double max = 0;
  double v[3];


  max = fabs( theR[0][0][0] );
  for ( z=0; z<(int)imWeight->dim.z; z++ )
  for ( y=0; y<(int)imWeight->dim.y; y++ )
  for ( x=0; x<(int)imWeight->dim.x; x++ ) {
    if ( fabs( theR[z][y][x] ) < threshold ) continue;
    nb ++;
    if ( max < fabs( theR[z][y][x] ) ) max = fabs( theR[z][y][x] );
  }
  
  cx = (float*)malloc( 3*nb* sizeof( float ) );
  cz = cy = cx;
  cy += nb;
  cz += 2*nb;
  
  ix = (int*)malloc( 3*nb* sizeof(int) );
  iz = iy = ix;
  iy += nb;
  iz += 2*nb;


  for ( n=0, z=0; z<(int)imWeight->dim.z; z++ )
  for ( y=0; y<(int)imWeight->dim.y; y++ )
  for ( x=0; x<(int)imWeight->dim.x; x++ ) {
    if ( fabs( theR[z][y][x] ) < threshold ) continue;
    ix[n] = x + 1;
    iy[n] = y + 1;
    iz[n] = z + 1;
    SphericalAnglesToUnitVector( (double)theTHETA[z][y][x], (double)thePHI[z][y][x], v );
    cx[n] = 1.25 * fabs( theR[z][y][x] )/max * v[0];
    cy[n] = 1.25 * fabs( theR[z][y][x] )/max * v[1];
    cz[n] = 1.25 * fabs( theR[z][y][x] )/max * v[2];
    n++;
  }
  
  if ( write( rawDataFileDesc, ix, 3*nb * sizeof( int ) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }

  if ( write( rawDataFileDesc, cx, 3*nb * sizeof( float ) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }
  
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "[IX, DIMIX] = fread( fid, %d, 'int%lu' );\n", nb, 8*sizeof( int ) );
  fprintf( commandFile, "[IY, DIMIY] = fread( fid, %d, 'int%lu' );\n", nb, 8*sizeof( int ) );
  fprintf( commandFile, "[IZ, DIMIZ] = fread( fid, %d, 'int%lu' );\n", nb, 8*sizeof( int ) );
  fprintf( commandFile, "[VX, DIMVX] = fread( fid, %d, 'float%lu' );\n", nb, 8*sizeof( float ) );
  fprintf( commandFile, "[VY, DIMVY] = fread( fid, %d, 'float%lu' );\n", nb, 8*sizeof( float ) );
  fprintf( commandFile, "[VZ, DIMVZ] = fread( fid, %d, 'float%lu' );\n", nb, 8*sizeof( float ) );
  fprintf( commandFile, "quiver3(IX,IY,IZ, VX, VY, VZ, 0, 'r' );\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );
  fprintf( commandFile, "\n" );

  free( cx );
  free( ix );
  
}


















