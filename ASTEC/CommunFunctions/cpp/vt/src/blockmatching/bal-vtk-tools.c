/*************************************************************************
 * bal-vtk-tools.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2018, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Lun 22 jan 2018 18:21:33 CET
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <vtmalloc.h>

#include <bal-vtk-tools.h>


static int _debug_ = 1;
static int _verbose_ = 1;



void BAL_SetVerboseInBalVtkTools( int v )
{
  _verbose_ = v;
}

void BAL_IncrementVerboseInBalVtkTools( )
{
  _verbose_ ++;
}

void BAL_DecrementVerboseInBalVtkTools( )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}

void BAL_SetDebugInBalVtkTools( int v )
{
  _debug_ = v;
}

void BAL_IncrementDebugInBalVtkTools( )
{
  _debug_ ++;
}

void BAL_DecrementDebugInBalVtkTools( )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}




/************************************************************
 *
 *
 *
 ************************************************************/




int BAL_TransformVtk( bal_vtkDataSet *theSet, bal_vtkDataSet *resSet, bal_transformation *theTrsf )
{
  char *proc = "BAL_TransformVtk";
  int i, p;
  int ipoints = -1;
  int inormals = -1;
  int opoints = -1;
  int onormals = -1;

  float *theBufPts = (float*)NULL;
  float *theBufNormals = (float*)NULL;
  float *resBufPts = (float*)NULL;
  float *resBufNormals = (float*)NULL;
  float n;

  bal_floatPoint thePt, resPt, theNo, resNo;


  if ( theSet->unit != theTrsf->transformation_unit ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: input mesh and transformation have different units\n", proc );
    return( -1 );
  }

  for ( i=0; i<theSet->n_data; i++ ) {
    switch( theSet->data[i].type ) {
    default :
      break;
    case _VTK_POINTS_ :
      ipoints = i; break;
    case _VTK_NORMALS_ :
      inormals = i; break;
    }
  }

  for ( i=0; i<resSet->n_data; i++ ) {
    switch( resSet->data[i].type ) {
    default :
      break;
    case _VTK_POINTS_ :
      opoints = i; break;
    case _VTK_NORMALS_ :
      onormals = i; break;
    }
  }

  if ( ipoints == -1 || opoints == -1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: no input or output points in vtk meshes?!\n", proc );
    return( -1 );
  }

  /* ...
   */
  if ( inormals >= 0 && onormals >= 0 ) {

    theBufPts = (float*)theSet->data[ipoints].data;
    theBufNormals = (float*)theSet->data[inormals].data;
    resBufPts = (float*)resSet->data[opoints].data;
    resBufNormals = (float*)resSet->data[onormals].data;

    for ( p=0; p<theSet->data[ipoints].n_primitives
          && p<theSet->data[inormals].n_primitives
          && p<resSet->data[opoints].n_primitives
          && p<resSet->data[onormals].n_primitives; p++ ) {
      thePt.x = theBufPts[3*p];
      thePt.y = theBufPts[3*p+1];
      thePt.z = theBufPts[3*p+2];
      theNo.x = thePt.x + theBufNormals[3*p];
      theNo.y = thePt.y + theBufNormals[3*p+1];
      theNo.z = thePt.z + theBufNormals[3*p+2];
      if ( BAL_TransformFloatPoint( &thePt, &resPt, theTrsf ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to transform point #1\n", proc );
        return( -1 );
      }
      if ( BAL_TransformFloatPoint( &theNo, &resNo, theTrsf ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to transform point #2\n", proc );
        return( -1 );
      }
      resNo.x -= resPt.x;
      resNo.y -= resPt.y;
      resNo.z -= resPt.z;
      n = sqrt( resNo.x*resNo.x + resNo.y*resNo.y + resNo.z*resNo.z );
      resNo.x /= n;
      resNo.y /= n;
      resNo.z /= n;
      resBufPts[3*p]   = resPt.x;
      resBufPts[3*p+1] = resPt.y;
      resBufPts[3*p+2] = resPt.z;
      resBufNormals[3*p]   = resNo.x;
      resBufNormals[3*p+1] = resNo.y;
      resBufNormals[3*p+2] = resNo.z;
    }

  }
  else {

    theBufPts = (float*)theSet->data[ipoints].data;
    resBufPts = (float*)resSet->data[opoints].data;

    for ( p=0; p<theSet->data[ipoints].n_primitives
          && p<resSet->data[opoints].n_primitives; p++ ) {
      thePt.x = theBufPts[3*p];
      thePt.y = theBufPts[3*p+1];
      thePt.z = theBufPts[3*p+2];
      if ( BAL_TransformFloatPoint( &thePt, &resPt, theTrsf ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to transform point #1\n", proc );
        return( -1 );
      }
      resBufPts[3*p]   = resPt.x;
      resBufPts[3*p+1] = resPt.y;
      resBufPts[3*p+2] = resPt.z;
    }

  }

  resSet->unit = theSet->unit;

  return( 1 );
}




/************************************************************
 *
 *
 *
 ************************************************************/

