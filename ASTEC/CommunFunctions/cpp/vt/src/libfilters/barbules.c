/*************************************************************************
 * barbules.c -
 *
 * $Id$
 *
 * Copyrigh (c) INRIA 2000
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Mon Sep 11 18:20:24 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 */

#include <tube3Dutil.h>

static int _verbose_ = 0;


void _VerboseInTube3D()
{
  if ( _verbose_ <= 0 ) _verbose_ = 1;
  else _verbose_ ++;
}
void _NoVerboseInTube3D()
{
  _verbose_ = 0;
}


int remove_spurious_branches( unsigned char *
