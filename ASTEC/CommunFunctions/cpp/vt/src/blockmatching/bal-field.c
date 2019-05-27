/*************************************************************************
 * bal-field.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2012, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mon Nov 19 17:45:00 CET 2012
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

#include <vtmalloc.h>

#include <bal-field.h>



static int _verbose_ = 1;
static int _debug_ = 0;


/*--------------------------------------------------*
 *
 * FIELD MANAGEMENT
 *
 *--------------------------------------------------*/



static void BAL_InitFieldGeometry( FIELD *field )
{
  field->vx = 1.0;
  field->vy = 1.0;
  field->vz = 1.0;

  field->geometry = _BAL_UNKNOWN_GEOMETRY_;

  /* matrices
   * - from voxel grid (integer) to real space
   * - from real space to voxel grid
   */
  _init_mat( &(field->to_real) );
  _init_mat( &(field->to_voxel) );
}



static void BAL_FreeFieldGeometry( FIELD *field )
{
  _free_mat( &(field->to_real) );
  _free_mat( &(field->to_voxel) );

  BAL_InitFieldGeometry( field );
}



int BAL_AllocFieldGeometry( FIELD *field )
{
  char *proc = "BAL_AllocFieldGeometry";

  /* desallocation of matrices if dimensions are not correct
   */
  if ( (field->to_real.l > 0 && field->to_real.l != 4)
       || (field->to_real.c > 0 && field->to_real.c != 4) )
    _free_mat( &(field->to_real) );

  if ( (field->to_voxel.l > 0 && field->to_voxel.l != 4)
       || (field->to_voxel.c > 0 && field->to_voxel.c != 4) )
    _free_mat( &(field->to_voxel) );

  /* allocations
   */
  if ( field->to_real.m == (double*)NULL ) {
    if ( _alloc_mat( &(field->to_real), 4, 4 ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: can not allocate 'to_real' matrix\n", proc );
        return( -1 );
    }
  }

  if ( field->to_voxel.m == (double*)NULL ) {
    if ( _alloc_mat( &(field->to_voxel), 4, 4 ) != 1 ) {
      _free_mat( &(field->to_real) );
        if ( _verbose_ )
          fprintf( stderr, "%s: can not allocate 'to_voxel' matrix\n", proc );
        return( -1 );
    }
  }

  /* set matrices to identity
   */
  field->vx = 1.0;
  field->vy = 1.0;
  field->vz = 1.0;
  _identity_mat( &(field->to_real) );
  _identity_mat( &(field->to_voxel) );

  return( 1 );
}





int BAL_SetFieldVoxelSizes( FIELD *field, typeVoxelSize vx, typeVoxelSize vy, typeVoxelSize vz )
{
  char *proc = "BAL_SetFieldVoxelSizes";

  switch ( field->geometry ) {
  case _BAL_UNKNOWN_GEOMETRY_ :
    /* geometry is set to identity
     */
    if ( BAL_AllocFieldGeometry( field ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate field geometry\n", proc );
      return( -1 );
    }
  case _BAL_HOMOTHETY_GEOMETRY_ :
    break;
  case _BAL_TRANSLATION_GEOMETRY_ :
  case _BAL_QFORM_GEOMETRY_ :
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: image orientation information may be inconsistent\n", proc );
  }

  /* recompute geometry information
   */
  field->vx = vx;
  field->vy = vy;
  field->vz = vz;

  _identity_mat( &(field->to_real) );
  field->to_real.m[ 0] = vx;
  field->to_real.m[ 5] = vy;
  field->to_real.m[10] = vz;
  field->to_real.m[15] = 1.0;

  _identity_mat( &(field->to_voxel) );
  field->to_voxel.m[ 0] = 1.0 / vx;
  field->to_voxel.m[ 5] = 1.0 / vy;
  field->to_voxel.m[10] = 1.0 / vz;
  field->to_voxel.m[15] = 1.0;

  field->geometry = _BAL_HOMOTHETY_GEOMETRY_;

  return( 1 );
}





static void BAL_InitField( FIELD *field )
{
  field->data = (typeScalarWeightedDisplacement*)NULL;
  field->pointer = (typeScalarWeightedDisplacement**)NULL;

  field->n_allocated_pairs = 0;

  field->n_computed_pairs = 0;

  field->n_selected_pairs = 0;

  field->unit =  VOXEL_UNIT;

  BAL_InitFieldGeometry( field );
}





static void Init_ScalarWeightedDisplacement( typeScalarWeightedDisplacement *d )
{
  d->origin.x = 0.0;
  d->origin.y = 0.0;
  d->origin.z = 0.0;
  d->vector.x = 0.0;
  d->vector.y = 0.0;
  d->vector.z = 0.0;
  d->valid = 0;
  d->rho = 0.0;
  d->error = 0.0;
}



/*---------------------------------------------------------
          Allocation d'un champ 3D de type double
----------------------------------------------------------*/

int BAL_AllocateField ( FIELD *field, size_t npairs )
{
  char *proc = "BAL_AllocateField";
  size_t i;

  BAL_InitField( field );

  if ( npairs <= 0 ) return( -1 );
  
  if ( BAL_AllocFieldGeometry( field ) != 1 ) {
    if ( _verbose_ )
        fprintf( stderr, "%s: geometry allocation failed\n", proc );
    return( -1 );
  }

  if ( _debug_ ) 
    fprintf( stderr, "%s: allocate %lu pairs\n", proc, npairs );

  field->data = (typeScalarWeightedDisplacement*)vtmalloc( npairs * sizeof(typeScalarWeightedDisplacement),
                                                           "field->data", proc );
  if ( field->data == (typeScalarWeightedDisplacement*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: data allocation failed\n", proc );
    BAL_FreeFieldGeometry( field );
    return( -1 );
  }

  field->pointer = (typeScalarWeightedDisplacement**)vtmalloc( npairs * sizeof(typeScalarWeightedDisplacement*),
                                                               "field->pointer", proc );
  if ( field->pointer == (typeScalarWeightedDisplacement**)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: pointer allocation failed\n", proc );
    BAL_FreeFieldGeometry( field );
    vtfree( field->data );
    return( -1 );
  }

  for ( i=0; i<npairs; i++ ) {
    field->pointer[i] = &( field->data[i] );
    Init_ScalarWeightedDisplacement( field->pointer[i] );
  }
  
  field->n_allocated_pairs = npairs;

  return( 1 );
}





/*----------------------------------------------------
         Desallocation d'un champ 3D de type double
-----------------------------------------------------*/

void BAL_FreeField ( FIELD * field )
{
  if ( field->data != NULL ) vtfree( field->data );
  if ( field->pointer != NULL ) vtfree( field->pointer );
  BAL_FreeFieldGeometry( field );
  BAL_InitField( field );
}





/*----------------------------------------------------
         output procedures
-----------------------------------------------------*/

void BAL_PrintField( FILE *f, FIELD *field )
{
  size_t i;
  
  fprintf( f, "Pairings: allocated = %lu, computed = %lu, selected = %lu\n",
           field->n_allocated_pairs, field->n_computed_pairs, field->n_selected_pairs );

  fprintf( f, "  - field geometry is: " );
  switch ( field->geometry ) {
  default :           fprintf( f, "not handled\n" ); break;
  case _BAL_UNKNOWN_GEOMETRY_ : fprintf( f, "_BAL_UNKNOWN_GEOMETRY_\n" ); break;
  case _BAL_HOMOTHETY_GEOMETRY_ : fprintf( f, "_BAL_HOMOTHETY_GEOMETRY_\n" ); break;
  case _BAL_TRANSLATION_GEOMETRY_ : fprintf( f, "_BAL_TRANSLATION_GEOMETRY_\n" ); break;
  case _BAL_QFORM_GEOMETRY_ :   fprintf( f, "_BAL_QFORM_GEOMETRY_\n" ); break;
  }
  fprintf( f, "  - conversion matrix voxel to real:\n" );
  _print_mat ( f, &(field->to_real), (char *)NULL );
  fprintf( f, "  - conversion matrix real to voxel:\n" );
  _print_mat ( f, &(field->to_voxel), (char *)NULL );


  for ( i=0; i<field->n_computed_pairs; i++ ) {
    fprintf( f, "#%2lu: (%f %f %f) <-> (%f %f %f) weight = %f",
             i,
             field->pointer[i]->origin.x, field->pointer[i]->origin.y, field->pointer[i]->origin.z,
             field->pointer[i]->vector.x, field->pointer[i]->vector.y, field->pointer[i]->vector.z,
             field->pointer[i]->rho );
    if ( 1 ) fprintf( f, " error = %20.15f", field->pointer[i]->error );
    fprintf( f, "\n" );
    
  }

}



void BAL_PrintSelectedPairsOfField( FILE *f, FIELD *field )
{
  size_t i;
  
  fprintf( f, "Pairings: allocated = %lu, computed = %lu, selected = %lu\n",
           field->n_allocated_pairs, field->n_computed_pairs, field->n_selected_pairs );
  for ( i=0; i<field->n_selected_pairs; i++ ) {
    fprintf( f, "#%2lu: (%f %f %f) <-> (%f %f %f) error = %f\n",
             i,
             field->pointer[i]->origin.x, field->pointer[i]->origin.y, field->pointer[i]->origin.z,
             field->pointer[i]->vector.x, field->pointer[i]->vector.y, field->pointer[i]->vector.z,
             field->pointer[i]->error );
  }

}



void BAL_PrintValidPairsOfField( FILE *f, FIELD *field )
{
  size_t i, j;
  
  fprintf( f, "Pairings: allocated = %lu, computed = %lu, selected = %lu\n",
           field->n_allocated_pairs, field->n_computed_pairs, field->n_selected_pairs );
  for ( i=0, j=0; i<field->n_computed_pairs; i++ ) {
    if ( field->pointer[i]->valid ) {
      fprintf( f, "#%2lu: (%f %f %f) <-> (%f %f %f) error = %f\n",
               j,
               field->pointer[i]->origin.x, field->pointer[i]->origin.y, field->pointer[i]->origin.z,
               field->pointer[i]->vector.x, field->pointer[i]->vector.y, field->pointer[i]->vector.z,
               field->pointer[i]->error );
      j++;
    }
  }

}





/*--------------------------------------------------*
 *
 * MISC
 *
 *--------------------------------------------------*/

void BAL_ChangeFieldToVoxelUnit(  FIELD *field )
{
  size_t i;
  bal_typeFieldPoint *pt;
  bal_typeFieldPoint res;
  double *m = (double*)NULL;
  typeVoxelSize vx = field->vx;
  typeVoxelSize vy = field->vy;
  typeVoxelSize vz = field->vz;

  /* nothing to do
   */
  if ( field->unit == VOXEL_UNIT ) return;


  /* is there a qmat matrix a la nifti
   */
  if ( field->to_voxel.l == 4 && field->to_voxel.c == 4 && field->to_voxel.m != (double*)NULL ) {
    m = field->to_voxel.m;

    switch( field->geometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
      for (i = 0; i < field->n_computed_pairs; i++) {
        pt = &(field->pointer[i]->origin);
        res.x = m[ 0] * pt->x                                ;
        res.y =                 m[ 5] * pt->y                ;
        res.z =                                 m[10] * pt->z;
        field->pointer[i]->origin.x = res.x;
        field->pointer[i]->origin.y = res.y;
        field->pointer[i]->origin.z = res.z;

        pt = &(field->pointer[i]->vector);
        res.x = m[ 0] * pt->x                                ;
        res.y =                 m[ 5] * pt->y                ;
        res.z =                                 m[10] * pt->z;
        field->pointer[i]->vector.x = res.x;
        field->pointer[i]->vector.y = res.y;
        field->pointer[i]->vector.z = res.z;
      }
      break;
    case _BAL_TRANSLATION_GEOMETRY_ :
      for (i = 0; i < field->n_computed_pairs; i++) {
        pt = &(field->pointer[i]->origin);
        res.x = m[ 0] * pt->x                                 + m[ 3];
        res.y =                 m[ 5] * pt->y                 + m[ 7];
        res.z =                                 m[10] * pt->z + m[11];
        field->pointer[i]->origin.x = res.x;
        field->pointer[i]->origin.y = res.y;
        field->pointer[i]->origin.z = res.z;

        pt = &(field->pointer[i]->vector);
        res.x = m[ 0] * pt->x                                ;
        res.y =                 m[ 5] * pt->y                ;
        res.z =                                 m[10] * pt->z;
        field->pointer[i]->vector.x = res.x;
        field->pointer[i]->vector.y = res.y;
        field->pointer[i]->vector.z = res.z;
      }
      break;
    case _BAL_QFORM_GEOMETRY_ :
      for (i = 0; i < field->n_computed_pairs; i++) {
        pt = &(field->pointer[i]->origin);
        res.x = m[ 0] * pt->x + m[ 1] * pt->y + m[ 2] * pt->z + m[ 3];
        res.y = m[ 4] * pt->x + m[ 5] * pt->y + m[ 6] * pt->z + m[ 7];
        res.z = m[ 8] * pt->x + m[ 9] * pt->y + m[10] * pt->z + m[11];
        field->pointer[i]->origin.x = res.x;
        field->pointer[i]->origin.y = res.y;
        field->pointer[i]->origin.z = res.z;

        pt = &(field->pointer[i]->vector);
        res.x = m[ 0] * pt->x + m[ 1] * pt->y + m[ 2] * pt->z;
        res.y = m[ 4] * pt->x + m[ 5] * pt->y + m[ 6] * pt->z;
        res.z = m[ 8] * pt->x + m[ 9] * pt->y + m[10] * pt->z;
        field->pointer[i]->vector.x = res.x;
        field->pointer[i]->vector.y = res.y;
        field->pointer[i]->vector.z = res.z;
      }
      break;
    }
  }
  else {
      for (i = 0; i < field->n_computed_pairs; i++) {
        pt = &(field->pointer[i]->origin);
        res.x = m[ 0] / vx                                ;
        res.y =                 m[ 5] / vy                ;
        res.z =                                 m[10] / vz;
        field->pointer[i]->origin.x = res.x;
        field->pointer[i]->origin.y = res.y;
        field->pointer[i]->origin.z = res.z;

        pt = &(field->pointer[i]->vector);
        res.x = m[ 0] / vx                                ;
        res.y =                 m[ 5] / vy                ;
        res.z =                                 m[10] / vz;
        field->pointer[i]->vector.x = res.x;
        field->pointer[i]->vector.y = res.y;
        field->pointer[i]->vector.z = res.z;
      }
  }


  
  field->unit = VOXEL_UNIT;
}





void BAL_ChangeFieldToRealUnit(  FIELD *field )
{
  size_t i;
  bal_typeFieldPoint *pt;
  bal_typeFieldPoint res;
  double *m = (double*)NULL;
  typeVoxelSize vx = field->vx;
  typeVoxelSize vy = field->vy;
  typeVoxelSize vz = field->vz;

  /* nothing to do
   */
  if ( field->unit == REAL_UNIT ) return;



  /* is there a qmat matrix a la nifti
   */
  if ( field->to_real.l == 4 && field->to_real.c == 4 && field->to_real.m != (double*)NULL ) {
    m = field->to_real.m;

    switch( field->geometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
      for (i = 0; i < field->n_computed_pairs; i++) {
          pt = &(field->pointer[i]->origin);
          res.x = m[ 0] * pt->x                                ;
          res.y =                 m[ 5] * pt->y                ;
          res.z =                                 m[10] * pt->z;
          field->pointer[i]->origin.x = res.x;
          field->pointer[i]->origin.y = res.y;
          field->pointer[i]->origin.z = res.z;

          pt = &(field->pointer[i]->vector);
          res.x = m[ 0] * pt->x                                ;
          res.y =                 m[ 5] * pt->y                ;
          res.z =                                 m[10] * pt->z;
          field->pointer[i]->vector.x = res.x;
          field->pointer[i]->vector.y = res.y;
          field->pointer[i]->vector.z = res.z;
      }
      break;
    case _BAL_TRANSLATION_GEOMETRY_ :
      for (i = 0; i < field->n_computed_pairs; i++) {
          pt = &(field->pointer[i]->origin);
          res.x = m[ 0] * pt->x                                 + m[ 3];
          res.y =                 m[ 5] * pt->y                 + m[ 7];
          res.z =                                 m[10] * pt->z + m[11];
          field->pointer[i]->origin.x = res.x;
          field->pointer[i]->origin.y = res.y;
          field->pointer[i]->origin.z = res.z;

          pt = &(field->pointer[i]->vector);
          res.x = m[ 0] * pt->x                                ;
          res.y =                 m[ 5] * pt->y                ;
          res.z =                                 m[10] * pt->z;
          field->pointer[i]->vector.x = res.x;
          field->pointer[i]->vector.y = res.y;
          field->pointer[i]->vector.z = res.z;
      }
      break;
    case _BAL_QFORM_GEOMETRY_ :
      for (i = 0; i < field->n_computed_pairs; i++) {
          pt = &(field->pointer[i]->origin);
          res.x = m[ 0] * pt->x + m[ 1] * pt->y + m[ 2] * pt->z + m[ 3];
          res.y = m[ 4] * pt->x + m[ 5] * pt->y + m[ 6] * pt->z + m[ 7];
          res.z = m[ 8] * pt->x + m[ 9] * pt->y + m[10] * pt->z + m[11];
          field->pointer[i]->origin.x = res.x;
          field->pointer[i]->origin.y = res.y;
          field->pointer[i]->origin.z = res.z;

          pt = &(field->pointer[i]->vector);
          res.x = m[ 0] * pt->x + m[ 1] * pt->y + m[ 2] * pt->z;
          res.y = m[ 4] * pt->x + m[ 5] * pt->y + m[ 6] * pt->z;
          res.z = m[ 8] * pt->x + m[ 9] * pt->y + m[10] * pt->z;
          field->pointer[i]->vector.x = res.x;
          field->pointer[i]->vector.y = res.y;
          field->pointer[i]->vector.z = res.z;
      }
      break;
    }
  }
  else {
      for (i = 0; i < field->n_computed_pairs; i++) {
          pt = &(field->pointer[i]->origin);
          res.x = m[ 0] * vx                                ;
          res.y =                 m[ 5] * vy                ;
          res.z =                                 m[10] * vz;
          field->pointer[i]->origin.x = res.x;
          field->pointer[i]->origin.y = res.y;
          field->pointer[i]->origin.z = res.z;

          pt = &(field->pointer[i]->vector);
          res.x = m[ 0] * vx                                ;
          res.y =                 m[ 5] * vy                ;
          res.z =                                 m[10] * vz;
          field->pointer[i]->vector.x = res.x;
          field->pointer[i]->vector.y = res.y;
          field->pointer[i]->vector.z = res.z;
      }
  }

  field->unit = REAL_UNIT;
}





/*--------------------------------------------------*
 *
 * MISC
 *
 *--------------------------------------------------*/

void CreateFileDef( FIELD *field,
                    char *nom_image, char *nom_champ )
{
  size_t i;
  FILE *def;
  
  def = fopen( nom_champ, "w");
  if ( def == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "CreateFileDef: unable to open '%s' for writing\n", nom_champ );
    return;
  }
  
  fprintf(def, "DEF=[\n");
  for(i = 0; i < field->n_computed_pairs; i++)
    fprintf(def, "%f %f %f %f %f %f %f\n", 
            field->pointer[i]->origin.x, field->pointer[i]->origin.y, field->pointer[i]->origin.z,
            field->pointer[i]->vector.x, field->pointer[i]->vector.y, field->pointer[i]->vector.z,
            field->pointer[i]->rho );
  fprintf(def, "];\n");
  fprintf(def, "X=DEF(:,1);\n");
  fprintf(def, "Y=DEF(:,2);\n");
  fprintf(def, "Z=DEF(:,3);\n");
  fprintf(def, "U=DEF(:,4);\n");
  fprintf(def, "V=DEF(:,5);\n");
  fprintf(def, "W=DEF(:,6);\n");
  fprintf(def, "RHO=DEF(:,7);\n");
  fprintf(def, "image_name = '%s';\n", nom_image );

  fclose(def);
}







