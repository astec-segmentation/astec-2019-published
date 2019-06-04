/*************************************************************************
 * bal-transformation-propagation.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2015, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mer 28 oct 2015 22:18:59 CET
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
#include <bal-transformation-propagation.h>



static int _verbose_ = 1;
static int _debug_ = 0;



void BAL_SetVerboseInBalTransformationPropagation( int v )
{
  _verbose_ = v;
}

void BAL_IncrementVerboseInBalTransformationPropagation(  )
{
  _verbose_ ++;
}

void BAL_DecrementVerboseInBalTransformationPropagation(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}

void BAL_SetDebugInBalTransformationPropagation( int d )
{
  _debug_ = d;
}

void BAL_IncrementDebugInBalTransformationPropagation(  )
{
  _debug_ ++;
}

void BAL_DecrementDebugInBalTransformationPropagation(  )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}




static int _EstimateTransformation( bal_transformationArray *theTrsfs,
                                    bal_transformation *resTrsf,
                                    int iref, int iprev, int icur )
{
    char *proc = "_EstimateTransformation";

    if ( _debug_ >= 2 ) {
        fprintf( stderr, "      %s( ..., ..., %d, %d, %d )\n", proc, iref, iprev, icur );
    }


    /* compose  T_{ icur <- iref } = T_{ icur <- iprev } o T_{ iprev <- iref }
     */

    /* check whether T_{ t <- ref } exists
     *               theTrsfs->array[ref][t]
     */
    if ( BAL_DoesTransformationExist( &(theTrsfs->array[iref][icur]) ) == 1 ) {
        if ( BAL_CopyTransformation( &(theTrsfs->array[iref][icur]), resTrsf ) != 1 ) {
            if ( _verbose_ )
                fprintf( stderr, "%s: unable to copy existing transformation T_{%d<-%d}\n",
                         proc, icur, iref );
            return( -1 );
        }
        return( 1 );
    }



    /* check whether T_{ iprev <- iref } exists
     *               theTrsfs->array[iref][iprev]
     * or create it from T_{ ief <- iprev }
     *               theTrsfs->array[iprev][iref]
     */
    if ( BAL_DoesTransformationExist( &(theTrsfs->array[iref][iprev]) ) != 1 ) {
        if ( BAL_DoesTransformationExist( &(theTrsfs->array[iprev][iref]) ) != 1 ) {
            if ( _verbose_ )
                fprintf( stderr, "%s: something weird with transformation #%d\n",
                         proc, iprev );
            return( -1 );
        }
        else {
            if ( _debug_ >= 2 ) {
                fprintf( stderr, "      invert transformation T_{%d<-%d}\n", iref, iprev );
            }
            if ( BAL_AllocTransformation( &(theTrsfs->array[iref][iprev]),
                                          resTrsf->type,
                                          &(resTrsf->vx) ) != 1 ) {
                if ( _verbose_ ) {
                    fprintf( stderr, "%s: unable to allocate transformation[%d][%d]\n",
                             proc, iref, iprev );
                }
                return( -1 );
            }
            if ( BAL_InverseTransformation( &(theTrsfs->array[iprev][iref]),
                                            &(theTrsfs->array[iref][iprev]) ) != 1 ) {
                if ( _verbose_ )
                    fprintf( stderr, "%s: unable to invert transformation T_{%d<-%d}\n",
                             proc, iref, iprev );
                return( -1 );
            }
        }
    }



    /* check whether T_{ icur <- iprev } exists
     *               theTrsfs->array[iprev][icur]
     * or create it from T_{ iprev <- icur }
     *               theTrsfs->array[icur][iprev]
     */
    if ( BAL_DoesTransformationExist( &(theTrsfs->array[iprev][icur]) ) != 1 ) {
        if ( BAL_DoesTransformationExist( &(theTrsfs->array[icur][iprev]) ) != 1 ) {
            if ( _verbose_ )
                fprintf( stderr, "%s: something weird with transformation #%d\n",
                         proc, icur );
            return( -1 );
        }
        else {
            if ( _debug_ >= 2 ) {
                fprintf( stderr, "      invert transformation T_{%d<-%d}\n", iprev, icur );
            }
            if ( BAL_AllocTransformation( &(theTrsfs->array[iprev][icur]),
                                          resTrsf->type,
                                          &(resTrsf->vx) ) != 1 ) {
                if ( _verbose_ ) {
                    fprintf( stderr, "%s: unable to allocate transformation[%d][%d]\n",
                             proc, iref, iprev );
                }
                return( -1 );
            }
            if ( BAL_InverseTransformation( &(theTrsfs->array[icur][iprev]),
                                            &(theTrsfs->array[iprev][icur]) ) != 1 ) {
                if ( _verbose_ )
                    fprintf( stderr, "%s: unable to invert transformation T_{%d<-%d}\n",
                             proc, iprev, icur );
                return( -1 );
            }
        }
    }


    
    /* allocate T_{ icur <- iref }
     *               theTrsfs->array[iref][icur]
     * compose T_{ icur <- iref } = T_{ icur <- iprev } o T_{ iprev <- iref }
     *               theTrsfs->array[iref][icur]
     *               theTrsfs->array[iprev][icur]
     *               theTrsfs->array[iref][iprev]
     *
     * Be careful, at the first call, iref and iprev can be equal, so
     * &(theTrsfs->array[iref][icur]) has been created/allocated at the
     * inversion step above
     *
     */
    if ( BAL_DoesTransformationExist( &(theTrsfs->array[iref][icur]) ) != 1 ) {
        if ( BAL_AllocTransformation( &(theTrsfs->array[iref][icur]),
                                      resTrsf->type,
                                      &(resTrsf->vx) ) != 1 ) {
            if ( _verbose_ ) {
                fprintf( stderr, "%s: unable to allocate transformation[%d][%d]\n",
                         proc, icur, iref );
            }
            return( -1 );
        }
    }

    if ( _debug_ >= 2 ) {
        fprintf( stderr, "      compose transformation T_{%d<-%d} = T_{%d<-%d} o T_{%d<-%d} \n",
                 icur, iref, icur, iprev, iprev, iref );
    }

    if ( BAL_TransformationComposition( &(theTrsfs->array[iref][icur]),
                                        &(theTrsfs->array[iprev][icur]),
                                        &(theTrsfs->array[iref][iprev])) != 1 ) {
        if ( _verbose_ ) {
            fprintf( stderr, "%s: unable to compose T_{%d<-%d} = T_{%d<-%d} o T_{%d<-%d}\n",
                     proc, icur, iref, icur, iprev, iprev, iref );
        }
        return( -1 );
    }
    if ( BAL_CopyTransformation( &(theTrsfs->array[iref][icur]), resTrsf ) != 1 ) {
        if ( _verbose_ )
            fprintf( stderr, "%s: unable to copy transformation T_{%d<-%d}\n",
                     proc, icur, iref );
        return( -1 );
    }

    return( 1 );
}




int BAL_EstimateTransformationByPropagation( bal_transformationArray *theTrsfs,
                                             bal_transformationList *resTrsfs,
                                             int referenceindex )
{
    char *proc = "BAL_EstimateTransformationByPropagation";
    int *doesTransformationExist = (int*)NULL;
    int i, j, p, t;



    doesTransformationExist = (int*)malloc( resTrsfs->n_trsfs * sizeof(int) );
    if ( doesTransformationExist == (int*)NULL ) {
        if ( _verbose_ )
            fprintf( stderr, "%s: unable to allocate auxiliary array\n", proc );
        return( -1 );
    }

    for ( t=0; t<resTrsfs->n_trsfs; t++ ) {
        doesTransformationExist[t] = 0;
        for ( j=0; j<theTrsfs->nrows && doesTransformationExist[t] == 0; j++ ) {
            if ( BAL_DoesTransformationExist( &(theTrsfs->array[j][t]) ) == 1 )
                doesTransformationExist[t] = 1;
        }
        for ( i=0; i<theTrsfs->ncols && doesTransformationExist[t] == 0; i++ ) {
            if ( BAL_DoesTransformationExist( &(theTrsfs->array[t][i]) ) == 1 )
                doesTransformationExist[t] = 1;
        }
        if ( _verbose_ >= 3 ) {
            if ( doesTransformationExist[t] == 0 )
                fprintf( stderr, "%s: #%3d does not exist\n", proc, t );
            else
                fprintf( stderr, "%s: #%3d does exist\n", proc, t );
        }
    }



    /* transformation estimation
     * theTrsfs->array[j][i] is an estimation of T_{ i <- j }
     * resTrsfs->data[j] is the current estimation of T_{ j <- ref }
     *
     * 1. resTrsfs->data[ref] is set to identity
     * 2. resTrsfs->data[i] = T_{ i <- j } o T_{ j <- ref }
     *    for valid i and j
     */

    if ( _debug_ >= 2 ) {
        fprintf( stderr, "  set resTrsfs[%d] to identity\n", referenceindex );
    }

    BAL_SetTransformationToIdentity( &(resTrsfs->data[referenceindex]) );
    if ( BAL_DoesTransformationExist( &(theTrsfs->array[referenceindex][referenceindex]) ) != 1 ) {
        if ( BAL_AllocTransformation( &(theTrsfs->array[referenceindex][referenceindex]),
                                      resTrsfs->data[referenceindex].type,
                                      &(resTrsfs->data[referenceindex].vx) ) != 1 ) {
            free( doesTransformationExist );
            if ( _verbose_ )
                fprintf( stderr, "%s: unable to allocate transformation[%d][%d]\n",
                         proc, referenceindex, referenceindex );
            return( -1 );
        }
    }
    BAL_SetTransformationToIdentity( &(theTrsfs->array[referenceindex][referenceindex]) );

    /* compose  T_{ t <- ref } = T_{ t <- p } o T_{ p <- ref }
     */

    if ( _debug_ >= 2 ) {
        fprintf( stderr, "  estimation backward\n" );
    }

    for ( p=referenceindex, t=referenceindex-1; t>=0; t-- ) {
        if ( doesTransformationExist[t] == 0 ) continue;

        if ( _debug_ >= 2 ) {
            fprintf( stderr, "    T_{ %d <- %d } = T_{ %d <- %d } o T_{ %d <- %d }\n", t, referenceindex, t, p, p, referenceindex );
        }

        if ( _EstimateTransformation( theTrsfs, &(resTrsfs->data[t]),
                                      referenceindex, p, t ) != 1 ) {
            free( doesTransformationExist );
            if ( _verbose_ )
                fprintf( stderr, "%s: unable to estimate transformation #%d\n", proc, t );
            return( -1 );
        }

        p = t;
    }

    if ( _debug_ >= 2 ) {
        fprintf( stderr, "  estimation forward\n" );
    }

    for ( p=referenceindex, t=referenceindex+1; t<resTrsfs->n_trsfs; t++ ) {
        if ( doesTransformationExist[t] == 0 ) continue;

        if ( _debug_ >= 2 ) {
            fprintf( stderr, "    T_{ %d <- %d } = T_{ %d <- %d } o T_{ %d <- %d }\n", t, referenceindex, t, p, p, referenceindex );
        }

        if ( _EstimateTransformation( theTrsfs, &(resTrsfs->data[t]),
                                      referenceindex, p, t ) != 1 ) {
            free( doesTransformationExist );
            if ( _verbose_ )
                fprintf( stderr, "%s: unable to estimate transformation #%d\n", proc, t );
            return( -1 );
        }

        p = t;
    }

    return( 1 );
}





