/*************************************************************************
 * bal-transformation-averaging.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Dim 16 jui 2013 16:56:13 CEST
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
#include <math.h>

#include <bal-transformation-compose.h>
#include <bal-transformation-copy.h>
#include <bal-transformation-inversion.h>
#include <bal-transformation-averaging.h>



static int _verbose_ = 1;
static int _debug_ = 0;



void BAL_SetVerboseInBalTransformationAveraging( int v )
{
  _verbose_ = v;
}

void BAL_IncrementVerboseInBalTransformationAveraging(  )
{
  _verbose_ ++;
}

void BAL_DecrementVerboseInBalTransformationAveraging(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}

void BAL_SetDebugInBalTransformationAveraging( int d )
{
  _debug_ = d;
}

void BAL_IncrementDebugInBalTransformationAveraging(  )
{
  _debug_ ++;
}

void BAL_DecrementDebugInBalTransformationAveraging(  )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}














/***************************************************
 *
 * LS averaging atomic procedure
 *
 ***************************************************/



static int LS_Translation_2D_Average( bal_transformationList *listTrsf,
                                      bal_transformation *resTrsf )
{
  int i;
  
  BAL_SetTransformationToIdentity( resTrsf );

  /* should check that the transformation are all in REAL_UNIT ?
   */

  for ( i=0; i<listTrsf->n_selected_trsfs; i++ ) {
    resTrsf->mat.m[ 3] += listTrsf->pointer[i]->mat.m[ 3];
    resTrsf->mat.m[ 7] += listTrsf->pointer[i]->mat.m[ 7];
  }

  resTrsf->mat.m[ 3] /= (double)listTrsf->n_selected_trsfs;
  resTrsf->mat.m[ 7] /= (double)listTrsf->n_selected_trsfs;

  return( 1 );
}





static int LS_Translation_3D_Average( bal_transformationList *listTrsf,
                                      bal_transformation *resTrsf )
{
  int i;
  
  BAL_SetTransformationToIdentity( resTrsf );

  /* should check that the transformation are all in REAL_UNIT ?
   */

  for ( i=0; i<listTrsf->n_selected_trsfs; i++ ) {
    resTrsf->mat.m[ 3] += listTrsf->pointer[i]->mat.m[ 3];
    resTrsf->mat.m[ 7] += listTrsf->pointer[i]->mat.m[ 7];
    resTrsf->mat.m[11] += listTrsf->pointer[i]->mat.m[11];
  }

  resTrsf->mat.m[ 3] /= (double)listTrsf->n_selected_trsfs;
  resTrsf->mat.m[ 7] /= (double)listTrsf->n_selected_trsfs;
  resTrsf->mat.m[11] /= (double)listTrsf->n_selected_trsfs;

  return( 1 );
}





static int LS_Translation_Scaling_2D_Average( bal_transformationList *listTrsf __attribute__ ((unused)),
                                              bal_transformation *resTrsf )
{
  char *proc = "LS_Translation_Scaling_2D_Average";
  BAL_SetTransformationToIdentity( resTrsf );
  resTrsf->mat.m[0] = resTrsf->mat.m[5] = resTrsf->mat.m[10] = resTrsf->mat.m[15] = 0.0;

  if ( _verbose_ )
    fprintf( stderr, "%s: not implemented yet\n", proc );
  
  return( -1 );
}





static int LS_Translation_Scaling_3D_Average( bal_transformationList *listTrsf __attribute__ ((unused)),
                                              bal_transformation *resTrsf )
{
  char *proc = "LS_Translation_Scaling_2D_Average";
  BAL_SetTransformationToIdentity( resTrsf );
  resTrsf->mat.m[0] = resTrsf->mat.m[5] = resTrsf->mat.m[10] = resTrsf->mat.m[15] = 0.0;

  if ( _verbose_ )
    fprintf( stderr, "%s: not implemented yet\n", proc );
  
  return( -1 );
}





static int LS_Rigid_2D_Average( bal_transformationList *listTrsf __attribute__ ((unused)),
                                              bal_transformation *resTrsf )
{
  char *proc = "LS_Rigid_3D_Average";
  BAL_SetTransformationToIdentity( resTrsf );
  resTrsf->mat.m[0] = resTrsf->mat.m[5] = resTrsf->mat.m[10] = resTrsf->mat.m[15] = 0.0;

  if ( _verbose_ )
    fprintf( stderr, "%s: not implemented yet\n", proc );
  
  return( -1 );
}





static int LS_Rigid_3D_Average( bal_transformationList *listTrsf __attribute__ ((unused)),
                                              bal_transformation *resTrsf )
{
  char *proc = "LS_Rigid_3D_Average";
  BAL_SetTransformationToIdentity( resTrsf );
  resTrsf->mat.m[0] = resTrsf->mat.m[5] = resTrsf->mat.m[10] = resTrsf->mat.m[15] = 0.0;

  if ( _verbose_ )
    fprintf( stderr, "%s: not implemented yet\n", proc );
  
  return( -1 );
}





static int LS_Similitude_2D_Average( bal_transformationList *listTrsf __attribute__ ((unused)),
                                              bal_transformation *resTrsf )
{
  char *proc = "LS_Similitude_2D_Average";
  BAL_SetTransformationToIdentity( resTrsf );
  resTrsf->mat.m[0] = resTrsf->mat.m[5] = resTrsf->mat.m[10] = resTrsf->mat.m[15] = 0.0;

  if ( _verbose_ )
    fprintf( stderr, "%s: not implemented yet\n", proc );
  
  return( -1 );
}





static int LS_Similitude_3D_Average( bal_transformationList *listTrsf __attribute__ ((unused)),
                                              bal_transformation *resTrsf )
{
  char *proc = "LS_Similitude_3D_Average";
  BAL_SetTransformationToIdentity( resTrsf );
  resTrsf->mat.m[0] = resTrsf->mat.m[5] = resTrsf->mat.m[10] = resTrsf->mat.m[15] = 0.0;

  if ( _verbose_ )
    fprintf( stderr, "%s: not implemented yet\n", proc );
  
  return( -1 );
}





static int LS_Affine_2D_Average( bal_transformationList *listTrsf,
                                              bal_transformation *resTrsf )
{
  int i;

  BAL_SetTransformationToIdentity( resTrsf );
  resTrsf->mat.m[0] = resTrsf->mat.m[5] = 0.0;

  /* should check that the transformation are all in REAL_UNIT ?
   */

  for ( i=0; i<listTrsf->n_selected_trsfs; i++ ) {
    resTrsf->mat.m[ 0] += listTrsf->pointer[i]->mat.m[ 0];
    resTrsf->mat.m[ 1] += listTrsf->pointer[i]->mat.m[ 1];
    resTrsf->mat.m[ 3] += listTrsf->pointer[i]->mat.m[ 3];

    resTrsf->mat.m[ 4] += listTrsf->pointer[i]->mat.m[ 4];
    resTrsf->mat.m[ 5] += listTrsf->pointer[i]->mat.m[ 5];
    resTrsf->mat.m[ 7] += listTrsf->pointer[i]->mat.m[ 7];
  }
 
  resTrsf->mat.m[ 0] /= (double)listTrsf->n_selected_trsfs;
  resTrsf->mat.m[ 1] /= (double)listTrsf->n_selected_trsfs;
  resTrsf->mat.m[ 3] /= (double)listTrsf->n_selected_trsfs;

  resTrsf->mat.m[ 4] /= (double)listTrsf->n_selected_trsfs;
  resTrsf->mat.m[ 5] /= (double)listTrsf->n_selected_trsfs;
  resTrsf->mat.m[ 7] /= (double)listTrsf->n_selected_trsfs;

  return( 1 );
}





static int LS_Affine_3D_Average( bal_transformationList *listTrsf,
                                              bal_transformation *resTrsf )
{
  int i;

  BAL_SetTransformationToIdentity( resTrsf );
  resTrsf->mat.m[0] = resTrsf->mat.m[5] = resTrsf->mat.m[10] = 0.0;

  /* should check that the transformation are all in REAL_UNIT ?
   */

  for ( i=0; i<listTrsf->n_selected_trsfs; i++ ) {
    resTrsf->mat.m[ 0] += listTrsf->pointer[i]->mat.m[ 0];
    resTrsf->mat.m[ 1] += listTrsf->pointer[i]->mat.m[ 1];
    resTrsf->mat.m[ 2] += listTrsf->pointer[i]->mat.m[ 2];
    resTrsf->mat.m[ 3] += listTrsf->pointer[i]->mat.m[ 3];

    resTrsf->mat.m[ 4] += listTrsf->pointer[i]->mat.m[ 4];
    resTrsf->mat.m[ 5] += listTrsf->pointer[i]->mat.m[ 5];
    resTrsf->mat.m[ 6] += listTrsf->pointer[i]->mat.m[ 6];
    resTrsf->mat.m[ 7] += listTrsf->pointer[i]->mat.m[ 7];

    resTrsf->mat.m[ 8] += listTrsf->pointer[i]->mat.m[ 8];
    resTrsf->mat.m[ 9] += listTrsf->pointer[i]->mat.m[ 9];
    resTrsf->mat.m[10] += listTrsf->pointer[i]->mat.m[10];
    resTrsf->mat.m[11] += listTrsf->pointer[i]->mat.m[11];
  }
 
  resTrsf->mat.m[ 0] /= (double)listTrsf->n_selected_trsfs;
  resTrsf->mat.m[ 1] /= (double)listTrsf->n_selected_trsfs;
  resTrsf->mat.m[ 2] /= (double)listTrsf->n_selected_trsfs;
  resTrsf->mat.m[ 3] /= (double)listTrsf->n_selected_trsfs;

  resTrsf->mat.m[ 4] /= (double)listTrsf->n_selected_trsfs;
  resTrsf->mat.m[ 5] /= (double)listTrsf->n_selected_trsfs;
  resTrsf->mat.m[ 6] /= (double)listTrsf->n_selected_trsfs;
  resTrsf->mat.m[ 7] /= (double)listTrsf->n_selected_trsfs;

  resTrsf->mat.m[ 8] /= (double)listTrsf->n_selected_trsfs;
  resTrsf->mat.m[ 9] /= (double)listTrsf->n_selected_trsfs;
  resTrsf->mat.m[10] /= (double)listTrsf->n_selected_trsfs;
  resTrsf->mat.m[11] /= (double)listTrsf->n_selected_trsfs;

  return( 1 );
}





static int LS_LinearTrsf_Average(  bal_transformationList *listTrsf,
                                   bal_transformation *resTrsf )
{
  char *proc = "LS_LinearTrsf_Average";

  switch( resTrsf->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such transformation not handled yet\n", proc );
    return( -1 );
  case TRANSLATION_2D :
    return( LS_Translation_2D_Average( listTrsf, resTrsf ) );
  case TRANSLATION_3D :
    return( LS_Translation_3D_Average( listTrsf, resTrsf ) );
  case TRANSLATION_SCALING_2D :
    return( LS_Translation_Scaling_2D_Average( listTrsf, resTrsf ) );
  case TRANSLATION_SCALING_3D :
    return( LS_Translation_Scaling_3D_Average( listTrsf, resTrsf ) );
  case RIGID_2D :
    return( LS_Rigid_2D_Average( listTrsf, resTrsf ) );
  case RIGID_3D :
    return( LS_Rigid_3D_Average( listTrsf, resTrsf ) );
  case SIMILITUDE_2D :
    return( LS_Similitude_2D_Average( listTrsf, resTrsf ) );
  case SIMILITUDE_3D :
    return( LS_Similitude_3D_Average( listTrsf, resTrsf ) );
  case AFFINE_2D :
    return( LS_Affine_2D_Average( listTrsf, resTrsf ) );
  case AFFINE_3D :
    return( LS_Affine_3D_Average( listTrsf, resTrsf ) );
  }
  return( -1 );
}










/***************************************************
 *
 * WLS averaging atomic procedure
 *
 ***************************************************/



static int WLS_Translation_2D_Average( bal_transformationList *listTrsf,
                                      bal_transformation *resTrsf )
{
  int i;
  double sum = 0.0;

  BAL_SetTransformationToIdentity( resTrsf );

  /* should check that the transformation are all in REAL_UNIT ?
   */

  for ( i=0; i<listTrsf->n_selected_trsfs; i++ ) {
    resTrsf->mat.m[ 3] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 3];
    resTrsf->mat.m[ 7] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 7];
    sum += listTrsf->pointer[i]->weight;
  }

  resTrsf->mat.m[ 3] /= sum;
  resTrsf->mat.m[ 7] /= sum;

  return( 1 );
}





static int WLS_Translation_3D_Average( bal_transformationList *listTrsf,
                                      bal_transformation *resTrsf )
{
  int i;
  double sum = 0.0;
  
  BAL_SetTransformationToIdentity( resTrsf );

  /* should check that the transformation are all in REAL_UNIT ?
   */

  for ( i=0; i<listTrsf->n_selected_trsfs; i++ ) {
    resTrsf->mat.m[ 3] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 3];
    resTrsf->mat.m[ 7] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 7];
    resTrsf->mat.m[11] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[11];
    sum += listTrsf->pointer[i]->weight;
  }

  resTrsf->mat.m[ 3] /= sum;
  resTrsf->mat.m[ 7] /= sum;
  resTrsf->mat.m[11] /= sum;

  return( 1 );
}





static int WLS_Translation_Scaling_2D_Average( bal_transformationList *listTrsf __attribute__ ((unused)),
                                              bal_transformation *resTrsf )
{
  char *proc = "WLS_Translation_Scaling_2D_Average";
  BAL_SetTransformationToIdentity( resTrsf );

  if ( _verbose_ )
    fprintf( stderr, "%s: not implemented yet\n", proc );
  
  return( -1 );
}





static int WLS_Translation_Scaling_3D_Average( bal_transformationList *listTrsf __attribute__ ((unused)),
                                              bal_transformation *resTrsf )
{
  char *proc = "WLS_Translation_Scaling_2D_Average";
  BAL_SetTransformationToIdentity( resTrsf );

  if ( _verbose_ )
    fprintf( stderr, "%s: not implemented yet\n", proc );
  
  return( -1 );
}





static int WLS_Rigid_2D_Average( bal_transformationList *listTrsf __attribute__ ((unused)),
                                              bal_transformation *resTrsf )
{
  char *proc = "WLS_Rigid_3D_Average";
  BAL_SetTransformationToIdentity( resTrsf );

  if ( _verbose_ )
    fprintf( stderr, "%s: not implemented yet\n", proc );
  
  return( -1 );
}





static int WLS_Rigid_3D_Average( bal_transformationList *listTrsf __attribute__ ((unused)),
                                              bal_transformation *resTrsf )
{
  char *proc = "WLS_Rigid_3D_Average";
  BAL_SetTransformationToIdentity( resTrsf );

  if ( _verbose_ )
    fprintf( stderr, "%s: not implemented yet\n", proc );
  
  return( -1 );
}





static int WLS_Similitude_2D_Average( bal_transformationList *listTrsf __attribute__ ((unused)),
                                              bal_transformation *resTrsf )
{
  char *proc = "WLS_Similitude_2D_Average";
  BAL_SetTransformationToIdentity( resTrsf );

  if ( _verbose_ )
    fprintf( stderr, "%s: not implemented yet\n", proc );
  
  return( -1 );
}





static int WLS_Similitude_3D_Average( bal_transformationList *listTrsf __attribute__ ((unused)),
                                              bal_transformation *resTrsf )
{
  char *proc = "WLS_Similitude_3D_Average";
  BAL_SetTransformationToIdentity( resTrsf );

  if ( _verbose_ )
    fprintf( stderr, "%s: not implemented yet\n", proc );
  
  return( -1 );
}





static int WLS_Affine_2D_Average( bal_transformationList *listTrsf,
                                              bal_transformation *resTrsf )
{
  int i;
  double sum = 0.0;

  BAL_SetTransformationToIdentity( resTrsf );
  resTrsf->mat.m[0] = resTrsf->mat.m[5] = 0.0;

  /* should check that the transformation are all in REAL_UNIT ?
   */

  for ( i=0; i<listTrsf->n_selected_trsfs; i++ ) {
    resTrsf->mat.m[ 0] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 0];
    resTrsf->mat.m[ 1] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 1];
    resTrsf->mat.m[ 3] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 3];

    resTrsf->mat.m[ 4] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 4];
    resTrsf->mat.m[ 5] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 5];
    resTrsf->mat.m[ 7] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 7];

    sum += listTrsf->pointer[i]->weight;
  }

  resTrsf->mat.m[ 0] /= sum;
  resTrsf->mat.m[ 1] /= sum;
  resTrsf->mat.m[ 3] /= sum;

  resTrsf->mat.m[ 4] /= sum;
  resTrsf->mat.m[ 5] /= sum;
  resTrsf->mat.m[ 7] /= sum;

  return( 1 );
}





static int WLS_Affine_3D_Average( bal_transformationList *listTrsf,
                                  bal_transformation *resTrsf )
{
  int i;
  double sum = 0.0;

  BAL_SetTransformationToIdentity( resTrsf );
  resTrsf->mat.m[0] = resTrsf->mat.m[5] = resTrsf->mat.m[10] = 0.0;

  
  /* should check that the transformation are all in REAL_UNIT ?
   */

  for ( i=0; i<listTrsf->n_selected_trsfs; i++ ) {
    resTrsf->mat.m[ 0] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 0];
    resTrsf->mat.m[ 1] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 1];
    resTrsf->mat.m[ 2] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 2];
    resTrsf->mat.m[ 3] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 3];

    resTrsf->mat.m[ 4] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 4];
    resTrsf->mat.m[ 5] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 5];
    resTrsf->mat.m[ 6] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 6];
    resTrsf->mat.m[ 7] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 7];

    resTrsf->mat.m[ 8] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 8];
    resTrsf->mat.m[ 9] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[ 9];
    resTrsf->mat.m[10] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[10];
    resTrsf->mat.m[11] += listTrsf->pointer[i]->weight * listTrsf->pointer[i]->mat.m[11];
    
    sum += listTrsf->pointer[i]->weight;
  }

  resTrsf->mat.m[ 0] /= sum;
  resTrsf->mat.m[ 1] /= sum;
  resTrsf->mat.m[ 2] /= sum;
  resTrsf->mat.m[ 3] /= sum;

  resTrsf->mat.m[ 4] /= sum;
  resTrsf->mat.m[ 5] /= sum;
  resTrsf->mat.m[ 6] /= sum;
  resTrsf->mat.m[ 7] /= sum;

  resTrsf->mat.m[ 8] /= sum;
  resTrsf->mat.m[ 9] /= sum;
  resTrsf->mat.m[10] /= sum;
  resTrsf->mat.m[11] /= sum;

  return( 1 );
}




static int WLS_LinearTrsf_Average( bal_transformationList *listTrsf,
                                   bal_transformation *resTrsf )
{
  char *proc = "WLS_LinearTrsf_Average";
  
  if ( _debug_ ) {
    fprintf( stdout, "%s: called with %d/%d transformations\n", proc, 
             listTrsf->n_selected_trsfs, listTrsf->n_trsfs );
  }

  switch( resTrsf->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such transformation not handled yet\n", proc );
    return( -1 );
  case TRANSLATION_2D :
    return( WLS_Translation_2D_Average( listTrsf, resTrsf ) );
  case TRANSLATION_3D :
    return( WLS_Translation_3D_Average( listTrsf, resTrsf ) );
  case TRANSLATION_SCALING_2D :
    return( WLS_Translation_Scaling_2D_Average( listTrsf, resTrsf ) );
  case TRANSLATION_SCALING_3D :
    return( WLS_Translation_Scaling_3D_Average( listTrsf, resTrsf ) );
  case RIGID_2D :
    return( WLS_Rigid_2D_Average( listTrsf, resTrsf ) );
  case RIGID_3D :
    return( WLS_Rigid_3D_Average( listTrsf, resTrsf ) );
  case SIMILITUDE_2D :
    return( WLS_Similitude_2D_Average( listTrsf, resTrsf ) );
  case SIMILITUDE_3D :
    return( WLS_Similitude_3D_Average( listTrsf, resTrsf ) );
  case AFFINE_2D :
    return( WLS_Affine_2D_Average( listTrsf, resTrsf ) );
  case AFFINE_3D :
    return( WLS_Affine_3D_Average( listTrsf, resTrsf ) );
  }
  return( -1 );
}











/***************************************************
 *
 * General procedures on transformation list
 *
 ***************************************************/

static int _Transformation_Residual( bal_transformation *t,
                                     double *eps_r,
                                     double *eps_t,
                                     double *eps )
{
  char *proc = "_Transformation_Residual";

  *eps_r = *eps_t = *eps = 0.0;

  switch( t->type ) {
  default :
    if ( _verbose_ ) 
      fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
    return( -1 );
    
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_3D :
  case RIGID_3D :
  case SIMILITUDE_3D :
  case AFFINE_3D :
    *eps_t += t->mat.m[11] * t->mat.m[11];
    /* Falls through. */

  case TRANSLATION_2D :
  case TRANSLATION_SCALING_2D :
  case RIGID_2D :
  case SIMILITUDE_2D :
  case AFFINE_2D :
    *eps_t += t->mat.m[ 3] * t->mat.m[ 3];
    *eps_t += t->mat.m[ 7] * t->mat.m[ 7];
    break;
  }

  switch( t->type ) {
  case TRANSLATION_SCALING_3D :
  case SIMILITUDE_3D :
  case TRANSLATION_SCALING_2D :
  case SIMILITUDE_2D :

  default :
    if ( _verbose_ ) 
      fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
    return( -1 );

  case TRANSLATION_3D :
  case TRANSLATION_2D :
    break;
    
  case RIGID_3D :
  case RIGID_2D :
    /* from Norms_Tmatrix() in bal-matrix.c
     */
    *eps_r = 6.0 - 2.0 * (t->mat.m[ 0] + t->mat.m[ 5] + t->mat.m[10]);
    if ( *eps_r < 0.0 ) *eps_r = 0.0;
    break;
    
  case AFFINE_3D :
      *eps_r += t->mat.m[ 2] * t->mat.m[ 2];
      *eps_r += t->mat.m[ 6] * t->mat.m[ 6];
      *eps_r += t->mat.m[ 8] * t->mat.m[ 8];
      *eps_r += t->mat.m[ 9] * t->mat.m[ 9];
      *eps_r += (t->mat.m[10] - 1.0) * (t->mat.m[10] - 1.0);
      /* Falls through. */

    case AFFINE_2D :
      *eps_r += (t->mat.m[ 0] - 1.0) * (t->mat.m[ 0] - 1.0);
      *eps_r += t->mat.m[ 1] * t->mat.m[ 1];
      *eps_r += t->mat.m[ 4] * t->mat.m[ 4];
      *eps_r += (t->mat.m[ 5] - 1.0) * (t->mat.m[ 5] - 1.0);
      break;
      
  }

  *eps = *eps_r + *eps_t;
  return( 1 );
}





static int _LinearTrsf_Residuals( bal_transformationList *listTrsf,
                                  bal_transformation *resTrsf )
{
  char *proc = "_LinearTrsf_Residuals";
  bal_transformation invTrsf;
  bal_transformation tmpTrsf;
  int i;
  double eps_r, eps_t, eps;

  /* allocations
   */
  
  BAL_InitTransformation( &invTrsf );
  if ( BAL_AllocTransformation( &invTrsf, resTrsf->type, (bal_image *)NULL ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate inverse transformation\n", proc );
    return( -1 );
  }

  BAL_InitTransformation( &tmpTrsf );
  if ( BAL_AllocTransformation( &tmpTrsf, resTrsf->type, (bal_image *)NULL ) != 1 ) {
    BAL_FreeTransformation( &invTrsf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate auxiliary transformation\n", proc );
    return( -1 );
  }


  /* transformation inversion
   */

  if ( BAL_InverseTransformation( resTrsf, &invTrsf ) != 1 ) {
    BAL_FreeTransformation( &tmpTrsf );
    BAL_FreeTransformation( &invTrsf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to invert transformation\n", proc );
    return( -1 );
  }


  
  /* residual calculation
   */

  for ( i=0; i<listTrsf->n_trsfs; i++ ) {

    if ( BAL_TransformationComposition( &tmpTrsf, &invTrsf, listTrsf->pointer[i] ) != 1 ) {
      BAL_FreeTransformation( &tmpTrsf );
      BAL_FreeTransformation( &invTrsf );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compose with transformation #%d\n", proc, i );
      return( -1 );
    }

    if ( _Transformation_Residual( &tmpTrsf, &eps_r, &eps_t, &eps ) != 1 ) {
      BAL_FreeTransformation( &tmpTrsf );
      BAL_FreeTransformation( &invTrsf );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute residuals of transformation #%d\n", proc, i );
      return( -1 );
    }

    listTrsf->pointer[i]->error = sqrt( eps );

  }

  BAL_FreeTransformation( &tmpTrsf );
  BAL_FreeTransformation( &invTrsf );

  return( 1 );
}





/* inspired from BAL_SelectSmallestResiduals()
   in bal-field-tools.c
*/

static int _Select_Smallest_Residuals( bal_transformationList *listTrsf,
                                       bal_estimator *estimator )
{
  char *proc = "_Select_Smallest_Residuals";

  size_t selectedresiduals, h;
  
  size_t last, left, right; 
  size_t i, j, k;

  double m, s, t;

  bal_transformation *tmp;

  /* no selection at the beginning
   */
  selectedresiduals = listTrsf->n_trsfs;

  /* selection on error distribution
   */
  if ( estimator->standard_deviation_threshold > 0.0 ) {

    m = s = 0.0;
    for ( m=0.0, i=0; i<(size_t)listTrsf->n_trsfs; i++ )
      m += listTrsf->pointer[i]->error;
    m /= listTrsf->n_trsfs;
    for ( s=0.0, i=0; i<(size_t)listTrsf->n_trsfs; i++ )
      s += (listTrsf->pointer[i]->error - m) * (listTrsf->pointer[i]->error - m);
    s /= listTrsf->n_trsfs;
    s = sqrt( s );

    /* threshold
     */
    t = m + estimator->standard_deviation_threshold * s;

    if ( _debug_ ) {
      fprintf( stderr, "%s: selection on error distribution\n", proc );
      fprintf( stderr, "m=%f, s=%f, threshold=%f\n", m, s, t );
    }

    /* selection
     */
    for ( i=0, selectedresiduals=listTrsf->n_trsfs; i<selectedresiduals; i++ ) {
      if ( listTrsf->pointer[i]->error > t ) {
        /* swap with right */
        tmp = listTrsf->pointer[i];
        listTrsf->pointer[i] = listTrsf->pointer[selectedresiduals-1];
        listTrsf->pointer[selectedresiduals-1] = tmp;
        selectedresiduals --;
        i --;
      }
    }

    
  }



  /* selection on percentage
   */
  if ( estimator->retained_fraction  <= 0.0 ||  1.0 <= estimator->retained_fraction ) 
    return( selectedresiduals );

  h = (int) ( estimator->retained_fraction  * (double) listTrsf->n_trsfs );
  
  if ( h <= 0 || h >= selectedresiduals )
    return( selectedresiduals );
  
  if ( _debug_ ) {
    fprintf( stderr, "%s: selection percentage\n", proc );
    fprintf( stderr, "will retain %lu / %d points\n", h, listTrsf->n_trsfs );
  }

  if ( _debug_ ) 
    fprintf( stderr, "%s: selection percentage, use home made procedure\n", proc );
  left = 0; 
  right = selectedresiduals-1;
  do {
    /* swap left et (left+right)/2 */
    j = (left+right)/2;
    if ( _debug_ ) fprintf( stderr, "[%lu - %lu] : swap %lu <-> %lu (test <-> pivot)\n", left, right, left, j );
    tmp = listTrsf->pointer[left]; listTrsf->pointer[left] = listTrsf->pointer[j]; listTrsf->pointer[j] = tmp;
    last = left;
    for ( k = left+1; k <= right; k++ ) {
 
      if ( listTrsf->pointer[k]->error < listTrsf->pointer[left]->error ) {
        last ++;
        /* if ( k > last ) */ {
          if ( _debug_ ) fprintf( stderr, "[%lu - %lu] : swap %lu <-> %lu\n", left, right, k, last );
          if ( _debug_ && ( left == 235 && right == 267 ) )
            fprintf( stderr, "\t e[k] = %20.15f < %20.15f = e[left]\n", listTrsf->pointer[k]->error, listTrsf->pointer[left]->error );
          tmp = listTrsf->pointer[k]; listTrsf->pointer[k] = listTrsf->pointer[last]; listTrsf->pointer[last] = tmp;
        }
      }
    }
    if ( _debug_ ) fprintf( stderr, "[%lu - %lu] : swap %lu <-> %lu (test <-> last)\n", left, right, left, last );
    tmp = listTrsf->pointer[left]; listTrsf->pointer[left] = listTrsf->pointer[last]; listTrsf->pointer[last] = tmp;
    if ( last >  h ) right = last - 1;
    if ( last <  h ) left  = last + 1;
  } while ( last != h );                         
  
  if ( _debug_ ) {
    for ( i=0;i<(size_t)listTrsf->n_trsfs;i++ )
      fprintf( stderr, "residual[%lu] = %14.9f\n", i, listTrsf->pointer[i]->error );
  }

  return( h );

}





static int _LinearTrsf_Average( bal_transformationList *listTrsf,
                                bal_transformation *resTrsf,
                                bal_estimator *estimator )
{
  char *proc = "_LinearTrsf_Average";
  
  if ( _debug_ ) {
    fprintf( stdout, "%s: called with %d/%d transformations\n", proc, listTrsf->n_selected_trsfs, listTrsf->n_trsfs );
  }

  switch( estimator->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: estimator type not handled in switch\n", proc );
    return( -1 );
  case TYPE_LS :
  case TYPE_LTS :
    return( LS_LinearTrsf_Average( listTrsf, resTrsf ) );
  case TYPE_WLS :
  case TYPE_WLTS :
    return( WLS_LinearTrsf_Average( listTrsf, resTrsf ) );
  }
  return( -1 );
}





static int _LinearTrsf_Trimmed_Average( bal_transformationList *listTrsf,
                                        bal_transformation *resTrsf,
                                        bal_estimator *estimator )
{
  char *proc = "_LinearTrsf_Trimmed_Average";

  int nretainedtrsfs;

  int iter;
  double eps_r;
  double eps_t;
  double eps;
  double tol_r = 1e-4;
  double tol_t = 1e-2;

  bal_transformation T0, tmpTrsf, invTrsf;

  if ( _debug_ ) {
    fprintf( stdout, "%s: called with %d transformations\n", proc, listTrsf->n_trsfs );
  }


#ifdef _ORIGINAL_BALADIN_TOLERANCE_IN_LTS_
  switch( estimator->type ) {
  default :
    break;
  case TYPE_LS :
  case TYPE_LTS :
    tol_r=1e-4;
    tol_t=1e-2;
    break;
  case TYPE_WLS :
  case TYPE_WLTS :
    tol_r=1e-1;
    tol_t=1e-1;
  }
#endif

  
  /* set return matrix to identity
   */
  BAL_SetTransformationToIdentity( resTrsf );



  /* initial transformation estimation
   */
  listTrsf->n_selected_trsfs = listTrsf->n_trsfs;
  if ( _LinearTrsf_Average( listTrsf, resTrsf, estimator ) <= 0 ) {
    if ( _verbose_ )
      fprintf ( stderr, "%s: error when estimating initial transformation\n", proc );
    return( -1 );
  }



  /* allocations
   */
  BAL_InitTransformation( &T0 );
  if ( BAL_AllocTransformation( &T0, resTrsf->type, (bal_image *)NULL ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate matrix #1\n", proc );
    return( -1 );
  }
  
  BAL_InitTransformation( &invTrsf );
  if ( BAL_AllocTransformation( &invTrsf, resTrsf->type, (bal_image *)NULL ) != 1 ) {
    BAL_FreeTransformation( &T0 );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate matrix #2\n", proc );
    return( -1 );
  }

  BAL_InitTransformation( &tmpTrsf );
  if ( BAL_AllocTransformation( &tmpTrsf, resTrsf->type, (bal_image *)NULL ) != 1 ) {
    BAL_FreeTransformation( &invTrsf );
    BAL_FreeTransformation( &T0 );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate matrix #3\n", proc );
    return( -1 );
  }



  /* loop
   */
  for ( iter = 0, eps_r = 2.0 * tol_r, eps_t = 2.0 * tol_t;
        ((eps_r > tol_r) || (eps_t > tol_t)) && (iter < estimator->max_iterations);
        iter ++ ) {

    /* copy of current transformation
     */
    if ( BAL_CopyTransformation( resTrsf, &T0 ) != 1 ) {
      BAL_FreeTransformation( &tmpTrsf );
      BAL_FreeTransformation( &invTrsf );
      BAL_FreeTransformation( &T0 );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to copy transformation\n", proc );
      return( -1 );
    }

    /* residual computation
     */
    if ( _LinearTrsf_Residuals( listTrsf, resTrsf ) != 1 ) {
      BAL_FreeTransformation( &tmpTrsf );
      BAL_FreeTransformation( &invTrsf );
      BAL_FreeTransformation( &T0 );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute residuals at iteration #%d\n", proc, iter );
      return( -1 );
    }

    /* sort residuals
     */
    nretainedtrsfs = _Select_Smallest_Residuals( listTrsf, estimator );
    if ( nretainedtrsfs <= 0 ) {
      BAL_FreeTransformation( &tmpTrsf );
      BAL_FreeTransformation( &invTrsf );
      BAL_FreeTransformation( &T0 );
      if ( _verbose_ )
        fprintf (stderr, "%s: no retained residuals? Returned value was %d\n", proc, nretainedtrsfs );
      return( -1 );
    }
    listTrsf->n_selected_trsfs = nretainedtrsfs;

    /* transformation estimation
     */
    if ( _LinearTrsf_Average( listTrsf, resTrsf, estimator ) <= 0 ) {
      if ( BAL_CopyTransformation( &T0, resTrsf ) != 1 ) {
        BAL_FreeTransformation( &tmpTrsf );
        BAL_FreeTransformation( &invTrsf );
        BAL_FreeTransformation( &T0 );
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to copy transformation\n", proc );
        return( -1 );
      }
      BAL_FreeTransformation( &tmpTrsf );
      BAL_FreeTransformation( &invTrsf );
      BAL_FreeTransformation( &T0 );
      if ( _verbose_ )
        fprintf ( stderr, "%s: error when estimating transformation at iteration #%d\n", proc, iter );
      return( -1 );
    }
    
    /* --- Comparaison de T et T0 --- */
    if ( BAL_InverseTransformation( &T0, &invTrsf ) != 1 ) {
      BAL_FreeTransformation( &tmpTrsf );
      BAL_FreeTransformation( &invTrsf );
      BAL_FreeTransformation( &T0 );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to invert transformation at iteration #%d\n", proc, iter );
      return( -1 );
    }
    
    if ( BAL_TransformationComposition( &tmpTrsf, &invTrsf, resTrsf ) != 1 ) {
      BAL_FreeTransformation( &tmpTrsf );
      BAL_FreeTransformation( &invTrsf );
      BAL_FreeTransformation( &T0 );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compose with transformation at iteration #%d\n", proc, iter );
      return( -1 );
    }
    
    if ( _Transformation_Residual( &tmpTrsf, &eps_r, &eps_t, &eps ) != 1 ) {
      BAL_FreeTransformation( &tmpTrsf );
      BAL_FreeTransformation( &invTrsf );
      BAL_FreeTransformation( &T0 );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute residuals at iteration #%d\n", proc, iter );
      return( -1 );
    }
    eps_r = sqrt( eps_r );
    eps_t = sqrt( eps_t);
    
    if ( _verbose_ >= 3 ) {
      switch ( estimator->type ) {
      default: break;
      case TYPE_LS :
      case TYPE_LTS :
        fprintf( stderr, "      LTS: iteration #%2d ... \r", iter );
        break;
      case TYPE_WLS :
      case TYPE_WLTS :
        fprintf( stderr, "      WLTS: iteration #%2d ... \r", iter );
        break;
      }
    }

  }

  BAL_FreeTransformation( &tmpTrsf );
  BAL_FreeTransformation( &invTrsf );
  BAL_FreeTransformation( &T0 );

  return( 1 );
}





int BAL_ComputeAverageLinearTransformation( bal_transformationList *listTrsf,
                                            bal_transformation *resTrsf,
                                            bal_estimator *estimator )
{
  char *proc= "BAL_ComputeAverageLinearTransformation";
  int i;
  
  if ( _debug_ >= 3 ) {
    fprintf( stdout, "%s: called with %d transformations\n", proc ,listTrsf->n_trsfs );
  }

  listTrsf->n_selected_trsfs = listTrsf->n_trsfs;

  if ( listTrsf->n_selected_trsfs <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: no selected transformations (n_trsfs=%d, n_allocated_trsfs=%d)\n", proc, 
               listTrsf->n_trsfs, listTrsf->n_allocated_trsfs );
    return( -1 );
  }

  if ( BAL_IsTransformationLinear( resTrsf ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: result transformation is not of linear type\n", proc );
    return( -1 );
  }

  for ( i=0; i<listTrsf->n_selected_trsfs; i++ ) {
    if ( BAL_IsTransformationLinear( listTrsf->pointer[i] ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: transformation #%d/%d is not of linear type\n", proc, i, listTrsf->n_selected_trsfs );
      return( -1 );
    }
  }

  switch ( estimator->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such transformation estimation not handled yet\n", proc );
    return( -1 );
  case TYPE_LS :
  case TYPE_WLS :
    return( _LinearTrsf_Average( listTrsf, resTrsf, estimator ) );
  case TYPE_LTS :
  case TYPE_WLTS :
    return( _LinearTrsf_Trimmed_Average( listTrsf, resTrsf, estimator ) );
  }
}





int BAL_ComputeAverageVectorFieldTransformation( bal_transformationList *listTrsf,
                                                 bal_transformation *resTrsf,
                                                 bal_estimator *estimator )
{
  char *proc= "BAL_ComputeAverageVectorFieldTransformation";
  int i;

  listTrsf->n_selected_trsfs = listTrsf->n_trsfs;

  if ( listTrsf->n_selected_trsfs <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: no selected transformations (n_trsfs=%d, n_allocated_trsfs=%d)\n", proc, 
               listTrsf->n_trsfs, listTrsf->n_allocated_trsfs );
    return( -1 );
  }

  if ( BAL_IsTransformationVectorField( resTrsf ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: result transformation is not of vector field type\n", proc );
    return( -1 );
  }

  for ( i=0; i<listTrsf->n_selected_trsfs; i++ ) {
    if ( BAL_IsTransformationVectorField( listTrsf->pointer[i] ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: transformation #%d/%d is not of vector field type\n", proc, i, listTrsf->n_selected_trsfs );
      return( -1 );
    }
  }

  switch ( estimator->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such transformation estimation not handled yet\n", proc );
    return( -1 );
  case TYPE_LS :
  case TYPE_WLS :
    if ( _verbose_ )
      fprintf( stderr, "%s: not implemented yet\n", proc );
    return( -1 );
  case TYPE_LTS :
  case TYPE_WLTS :
    if ( _verbose_ )
      fprintf( stderr, "%s: not implemented yet\n", proc );
    return( -1 );
  }
}










/***************************************************
 *
 * General procedures on transformation list
 *
 ***************************************************/




int BAL_ComputeAverageTransformation( bal_transformationList *listTrsf,
                                      bal_transformation *resTrsf,
                                      bal_estimator *estimator )
{
  char *proc= "BAL_ComputeAverageTransformation";

  if ( _debug_ >= 3 ) {
    fprintf( stdout, "%s: called with %d transformations\n", proc, listTrsf->n_trsfs );
  }

  listTrsf->n_selected_trsfs = listTrsf->n_trsfs;

  switch( resTrsf->type ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
    return( -1 );

  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :

    if ( BAL_ComputeAverageLinearTransformation( listTrsf, resTrsf, estimator ) <= 0 ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: average linear transformation computation failed\n", proc );
      return( -1 );
    }

    break;
    
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :

    /* le champ est en voxel,
       avec le point de depart dans l'image de reference
    */
    
    if ( BAL_ComputeAverageVectorFieldTransformation( listTrsf, resTrsf, estimator ) <= 0 ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: average vector field transformation computation failed\n", proc );
      return( -1 );
    }

  }

  return( 1 );

}








/***************************************************
 *
 * General procedures on multiple transformations
 *
 ***************************************************/




int BAL_EstimateTransformationByAveraging( bal_transformationArray *theTrsfs,
                                        bal_transformationList *resTrsfs,
                                        bal_estimator *estimator,
                                        double sigma,
                                        int n_iterations,
                                        int useinverse )
{
  char *proc = "BAL_EstimateTransformationByAveraging";

  bal_transformationList resAuxTrsfs;
  bal_transformationList compAuxTrsfs;

  bal_transformationList *prevList;
  bal_transformationList *nextList;
  bal_transformationList *auxList;

  int iteration;
  int i, j, k;

  int halfinterval;
  bal_estimator smoothing;

  /* allocations of transformations,
   * initialization to identity
   */

  BAL_InitTransformationList( &resAuxTrsfs );
  BAL_InitTransformationList( &compAuxTrsfs );

  switch( resTrsfs->data[0].type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
    return( -1 );

  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :
    
    if ( BAL_FullAllocTransformationList( &resAuxTrsfs, 
                                          resTrsfs->n_allocated_trsfs,
                                          resTrsfs->data[0].type,
                                          (bal_image *)NULL ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when allocating first auxiliary transformation list\n", proc );
      return( -1 );
    }
    k = resTrsfs->n_allocated_trsfs;
    if ( useinverse )
      k = 2*resTrsfs->n_allocated_trsfs;
    if ( BAL_FullAllocTransformationList( &compAuxTrsfs, 
                                          k,
                                          resTrsfs->data[0].type,
                                          (bal_image *)NULL ) != 1 ) {
      BAL_FreeTransformationList( &resAuxTrsfs );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when allocating second auxiliary transformation list\n", proc );
      return( -1 );
    }
    break;
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
    if ( _verbose_ )
      fprintf( stderr, "%s: such transformation type not handled yet (vector field)\n", proc );
    return( -1 );
  }

  /* ... */
  prevList = resTrsfs;
  nextList = &resAuxTrsfs;



  for ( iteration=0; iteration<n_iterations; iteration++ ) {

    if ( _verbose_ >= 2 )
      fprintf( stdout, "%s: processing iteration #%2d\n", proc, iteration );



    /*  transformation estimation
       theTrsfs->array[j][i] is an estimation of T_{ i <- j }
       prevList->data[j] is the current estimation of T_{ j <- ref }
       
       the next estimation of  T_{ i <- ref } is the average
       of all T_{ i <- j } o T_{ j <- ref }
       should we add T_{ j <- i }^(-1) o T_{ j <- ref } ?
    */
    for ( i=0; i<prevList->n_trsfs; i++ ) {

      /* construction of transformation list
       */
      for ( j=0, k=0; j<prevList->n_trsfs; j++ ) {
        if ( BAL_DoesTransformationExist( &(theTrsfs->array[j][i]) ) == 1 ) {
          if ( BAL_TransformationComposition( compAuxTrsfs.pointer[k],
                                              &(theTrsfs->array[j][i]),
                                              &(prevList->data[j]) ) != 1 ) {
            BAL_FreeTransformationList( &compAuxTrsfs );
            BAL_FreeTransformationList( &resAuxTrsfs );
            if ( _verbose_ ) {
              fprintf( stderr, "%s: at iteration #%d\n", proc, iteration );
              fprintf( stderr, "\t unable to compose transformation #%d\n", j );
              fprintf( stderr, "\t with estimate T_{%3d <- %3d}\n", j, i );
              return( -1 );
            }
          }
          k++;
        }
      }

      compAuxTrsfs.n_trsfs = k;

      if ( useinverse ) {

        if ( 0 )
          fprintf( stdout, "%s: number of direct transformations: %d\n", proc, k );

        for ( j=0; j<prevList->n_trsfs; j++ ) {
          if ( BAL_DoesTransformationExist( &(theTrsfs->array[i][j]) ) == 1 ) {
            if ( BAL_InverseTransformation( &(theTrsfs->array[i][j]), compAuxTrsfs.pointer[k] ) != 1 ) {
              BAL_FreeTransformationList( &compAuxTrsfs );
              BAL_FreeTransformationList( &resAuxTrsfs );
              if ( _verbose_ ) {
                fprintf( stderr, "%s: at iteration #%d\n", proc, iteration );
                fprintf( stderr, "\t unable to invert transformation T_{%3d <- %3d}\n", i, j );
                return( -1 );
              }
            }
            if ( BAL_TransformationComposition( compAuxTrsfs.pointer[k],
                                                compAuxTrsfs.pointer[k],
                                                &(prevList->data[j]) ) != 1 ) {
              BAL_FreeTransformationList( &compAuxTrsfs );
              BAL_FreeTransformationList( &resAuxTrsfs );
              if ( _verbose_ ) {
                fprintf( stderr, "%s: at iteration #%d\n", proc, iteration );
                fprintf( stderr, "\t unable to compose transformation #%d\n", i );
                fprintf( stderr, "\t with inverted estimate T_{%3d <- %3d}\n", i, j );
                return( -1 );
              }
            }
            k++;
          }
        }

        if ( 0 )
          fprintf( stdout, "%s: number of direct+inverse transformations: %d\n", proc, k );


      }

      compAuxTrsfs.n_trsfs = k;

      /* the following line appears in the mean computation,
         thus is no more necessary here
         compAuxTrsfs.n_selected_trsfs = compAuxTrsfs.n_trsfs;
      */

      /* transformation estimation
       */
      if ( BAL_ComputeAverageTransformation( &compAuxTrsfs, &(nextList->data[i]), estimator ) != 1 ) {
        BAL_FreeTransformationList( &compAuxTrsfs );
        BAL_FreeTransformationList( &resAuxTrsfs );
        if ( _verbose_ ) {
          fprintf( stderr, "%s: at iteration #%d\n", proc, iteration );
          fprintf( stderr, "\t unable to compute average for transformation #%d\n", i );
          return( -1 );
        }
      }
      
      /* estimation of transformation #i is done
       */
      if ( _debug_ ) {
        fprintf( stdout, "\n" );
        fprintf( stdout, "at iteration %2d, estimation of transformation #%3d has been done with %3d transformations\n",
                 iteration, i, k );
        BAL_PrintTransformation( stdout, &(nextList->data[i]), (char *)NULL );
      }

    }



    /* transformation smoothing
     *
     * Warning. It is mandatory to use pointer array and *not* the data one
     * the auxiliary array compAuxTrsfs, since the ordered pointers should be correspond to 
     * the ordered data in case of robust estimation (sorting will modify the relative orders)
     *
     */
    if ( sigma > 0.0 ) {

      auxList = prevList;
      prevList = nextList;
      nextList = auxList;

      halfinterval = (int)(3.0 * sigma + 0.5);
      if ( halfinterval <= 0 ) halfinterval = 1;

      if ( 0 ) {
        fprintf( stdout, "\n\n\n" );
        fprintf( stdout, "at iteration %2d, smoothing with %d transformations\n",
                 iteration, 2*halfinterval+1 );
        fprintf( stdout, "\n" );
      }

      BAL_InitEstimator( &smoothing );
      smoothing.type = TYPE_WLS;

      for ( i=0; i<prevList->n_trsfs; i++ ) {

        for ( j=i-halfinterval, k=0; j<prevList->n_trsfs && j<=i+halfinterval; j++ ) {
          if ( j < 0 ) continue;
          if ( BAL_CopyTransformation( &(prevList->data[j]),
                                       compAuxTrsfs.pointer[k] ) != 1 ) {
            BAL_FreeTransformationList( &compAuxTrsfs );
            BAL_FreeTransformationList( &resAuxTrsfs );
            if ( _verbose_ ) {
              fprintf( stderr, "%s: at iteration #%d\n", proc, iteration );
              fprintf( stderr, "\t unable to copy transformation #%d for smoothing transformation #%d\n", j, i );
              return( -1 );
            }
          }
          compAuxTrsfs.pointer[k]->weight = exp( -(j-i)*(j-i)/(2*sigma*sigma) );
          k++;
        }
        compAuxTrsfs.n_trsfs = k;

        if ( BAL_ComputeAverageTransformation( &compAuxTrsfs, &(nextList->data[i]), &smoothing ) != 1 ) {
          BAL_FreeTransformationList( &compAuxTrsfs );
          BAL_FreeTransformationList( &resAuxTrsfs );
          if ( _verbose_ ) {
            fprintf( stderr, "%s: at iteration #%d\n", proc, iteration );
            fprintf( stderr, "\t unable to compute smoothed transformation #%d\n", i );
            return( -1 );
          }
        }

        /* smoothing of transformation #i is done
         */
        if ( 0 ) {
          fprintf( stdout, "\n" );
          fprintf( stdout, "at iteration %2d, smoothing of transformation #%3d has been done with %3d transformations\n",
                   iteration, i, k );
          for ( j=0; j<k; j++ ) {
            fprintf( stdout, "--- transformation #%2d weighted by %f, for smoothing transformation #%d\n",
                     j, compAuxTrsfs.pointer[j]->weight, i );
            BAL_PrintTransformation( stdout, compAuxTrsfs.pointer[j], (char *)NULL );
          }
          fprintf( stdout, "\n" );
          fprintf( stdout, "--- previous value of transformation #%d\n", i );
          BAL_PrintTransformation( stdout, &(prevList->data[i]), (char *)NULL );
          fprintf( stdout, "--- new value of transformation #%d\n", i );
          BAL_PrintTransformation( stdout, &(nextList->data[i]), (char *)NULL );
        }

      }
      
    }



    /* preparation of next iteration
     */

    auxList = prevList;
    prevList = nextList;
    nextList = auxList;
  }

  BAL_FreeTransformationList( &compAuxTrsfs );

  if ( BAL_CopyTransformationList( prevList, resTrsfs ) != 1 ) {
    BAL_FreeTransformationList( &resAuxTrsfs );
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to copy transformation list\n", proc );
      return( -1 );
    }
  }

  BAL_FreeTransformationList( &resAuxTrsfs );

  return( 1 );
}
