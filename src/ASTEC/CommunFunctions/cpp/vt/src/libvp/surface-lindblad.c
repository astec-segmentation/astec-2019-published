/****************************************************
 * surface-lindblad.c -
 *
 * Copyright (c) INRIA 2018
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Mer  5 déc 2018 07:55:58 CET
 *
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */

#include <stdio.h>

#include <surface-lindblad.h>

#define A1  0.636
#define A2  0.669
#define A3  1.272
#define A4  1.272
#define A5  0.554
#define A6  1.305
#define A7  1.908
#define A8  0.927
#define A9  0.422
#define A10 1.338
#define A11 1.573
#define A12 1.190
#define A13 2.544

/* 0 1  -  4 5
 * 2 3  -  6 7
 *
 * grandes diagonales : 0-7 1-6 2-5 3-4
 *
 * 0 1 2 3 4 5 6 7
 */

float _surface_lindblad( int *binary )
{
  char * proc = "_surface_lindblad";

  if ( binary[0] > 0 ) {
    /* 1 . . . . . . . */
    if ( binary[1] > 0 ) {
      /* 1 1 . . . . . . */
      if ( binary[2] > 0 ) {
        /* 1 1 1 . . . . . */
        if ( binary[3] > 0 ) {
          /* 1 1 1 1 . . . . */
          if ( binary[4] > 0 ) {
            /* 1 1 1 1 1 . . . */
            if ( binary[5] > 0 ) {
              /* 1 1 1 1 1 1 . . */
              if ( binary[6] > 0) {
                /* 1 1 1 1 1 1 1 . */
                if ( binary[7] > 0) {
                  /* 1 1 1 1 1 1 1 1 */
                  return( 0.0 );
                }
                else {
                  return( A1 ); /* 1 1 1 1 1 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 1 1 1 1 1 0 . */
                if ( binary[7] > 0) {
                  return( A1 ); /* 1 1 1 1 1 1 0 1 */
                }
                else {
                  return( A2 ); /* 1 1 1 1 1 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 1 1 1 1 1 0 . . */
              if ( binary[6] > 0) {
                /* 1 1 1 1 1 0 1 . */
                if ( binary[7] > 0) {
                  return( A1 ); /* 1 1 1 1 1 0 1 1 */
                }
                else {
                  return( A2 ); /* 1 1 1 1 1 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 1 1 1 1 0 0 . */
                if ( binary[7] > 0) {
                  return( A3 ); /* 1 1 1 1 1 0 0 1 */
                }
                else {
                  return( A5 ); /* 1 1 1 1 1 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* binary[4] */
          else {
            /* 1 1 1 1 0 . . . */
            if ( binary[5] > 0 ) {
              /* 1 1 1 1 0 1 . . */
              if ( binary[6] > 0) {
                /* 1 1 1 1 0 1 1 . */
                if ( binary[7] > 0) {
                  return( A1 ); /* 1 1 1 1 0 1 1 1 */
                }
                else {
                  return( A3 ); /* 1 1 1 1 0 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 1 1 1 0 1 0 . */
                if ( binary[7] > 0) {
                  return( A2 ); /* 1 1 1 1 0 1 0 1 */
                }
                else {
                  return( A5 ); /* 1 1 1 1 0 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 1 1 1 1 0 0 . . */
              if ( binary[6] > 0) {
                /* 1 1 1 1 0 0 1 . */
                if ( binary[7] > 0) {
                  return( A2 ); /* 1 1 1 1 0 0 1 1 */
                }
                else {
                  return( A5 ); /* 1 1 1 1 0 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 1 1 1 0 0 0 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 1 1 1 1 0 0 0 1 */
                }
                else {
                  return( A8 ); /* 1 1 1 1 0 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* else binary[4] */
        } /* binary[3] */
        else {
          /* 1 1 1 0 . . . . */
          if ( binary[4] > 0 ) {
            /* 1 1 1 0 1 . . . */
            if ( binary[5] > 0 ) {
              /* 1 1 1 0 1 1 . . */
              if ( binary[6] > 0) {
                /* 1 1 1 0 1 1 1 . */
                if ( binary[7] > 0) {
                  return( A1 ); /* 1 1 1 0 1 1 1 1 */
                }
                else {
                  return( A2 ); /* 1 1 1 0 1 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 1 1 0 1 1 0 . */
                if ( binary[7] > 0) {
                  return( A3 ); /* 1 1 1 0 1 1 0 1 */
                }
                else {
                  return( A5 ); /* 1 1 1 0 1 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 1 1 1 0 1 0 . . */
              if ( binary[6] > 0) {
                /* 1 1 1 0 1 0 1 . */
                if ( binary[7] > 0) {
                  return( A3 ); /* 1 1 1 0 1 0 1 1 */
                }
                else {
                  return( A5 ); /* 1 1 1 0 1 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 1 1 0 1 0 0 . */
                if ( binary[7] > 0) {
                  return( A7 ); /* 1 1 1 0 1 0 0 1 */
                }
                else {
                  return( A9 ); /* 1 1 1 0 1 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* binary[4] */
          else {
            /* 1 1 1 0 0 . . . */
            if ( binary[5] > 0 ) {
              /* 1 1 1 0 0 1 . . */
              if ( binary[6] > 0) {
                /* 1 1 1 0 0 1 1 . */
                if ( binary[7] > 0) {
                  return( A4 ); /* 1 1 1 0 0 1 1 1 */
                }
                else {
                  return( A6 ); /* 1 1 1 0 0 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 1 1 0 0 1 0 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 1 1 1 0 0 1 0 1 */
                }
                else {
                  return( A11 ); /* 1 1 1 0 0 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 1 1 1 0 0 0 . . */
              if ( binary[6] > 0) {
                /* 1 1 1 0 0 0 1 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 1 1 1 0 0 0 1 1 */
                }
                else {
                  return( A11 ); /* 1 1 1 0 0 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 1 1 0 0 0 0 . */
                if ( binary[7] > 0) {
                  return( A12 ); /* 1 1 1 0 0 0 0 1 */
                }
                else {
                  return( A5 ); /* 1 1 1 0 0 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* else binary[4] */
        } /* else binary[3] */
      } /* binary[2] */
      else {
        /* 1 1 0 . . . . . */
        if ( binary[3] > 0 ) {
          /* 1 1 0 1 . . . . */
          if ( binary[4] > 0 ) {
            /* 1 1 0 1 1 . . . */
            if ( binary[5] > 0 ) {
              /* 1 1 0 1 1 1 . . */
              if ( binary[6] > 0) {
                /* 1 1 0 1 1 1 1 . */
                if ( binary[7] > 0) {
                  return( A1 ); /* 1 1 0 1 1 1 1 1 */
                }
                else {
                  return( A3 ); /* 1 1 0 1 1 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 1 0 1 1 1 0 . */
                if ( binary[7] > 0) {
                  return( A2 ); /* 1 1 0 1 1 1 0 1 */
                }
                else {
                  return( A5 ); /* 1 1 0 1 1 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 1 1 0 1 1 0 . . */
              if ( binary[6] > 0) {
                /* 1 1 0 1 1 0 1 . */
                if ( binary[7] > 0) {
                  return( A4 ); /* 1 1 0 1 1 0 1 1 */
                }
                else {
                  return( A6 ); /* 1 1 0 1 1 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 1 0 1 1 0 0 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 1 1 0 1 1 0 0 1 */
                }
                else {
                  return( A11 ); /* 1 1 0 1 1 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* binary[4] */
          else {
            /* 1 1 0 1 0 . . . */
            if ( binary[5] > 0 ) {
              /* 1 1 0 1 0 1 . . */
              if ( binary[6] > 0) {
                /* 1 1 0 1 0 1 1 . */
                if ( binary[7] > 0) {
                  return( A3 ); /* 1 1 0 1 0 1 1 1 */
                }
                else {
                  return( A7 ); /* 1 1 0 1 0 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 1 0 1 0 1 0 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 1 1 0 1 0 1 0 1 */
                }
                else {
                  return( A9 ); /* 1 1 0 1 0 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 1 1 0 1 0 0 . . */
              if ( binary[6] > 0) {
                /* 1 1 0 1 0 0 1 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 1 1 0 1 0 0 1 1 */
                }
                else {
                  return( A12 ); /* 1 1 0 1 0 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 1 0 1 0 0 0 . */
                if ( binary[7] > 0) {
                  return( A11 ); /* 1 1 0 1 0 0 0 1 */
                }
                else {
                  return( A5 ); /* 1 1 0 1 0 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* else binary[4] */
        } /* binary[3] */
        else {
          /* 1 1 0 0 . . . . */
          if ( binary[4] > 0 ) {
            /* 1 1 0 0 1 . . . */
            if ( binary[5] > 0 ) {
              /* 1 1 0 0 1 1 . . */
              if ( binary[6] > 0) {
                /* 1 1 0 0 1 1 1 . */
                if ( binary[7] > 0) {
                  return( A2 ); /* 1 1 0 0 1 1 1 1 */
                }
                else {
                  return( A5 ); /* 1 1 0 0 1 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 1 0 0 1 1 0 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 1 1 0 0 1 1 0 1 */
                }
                else {
                  return( A8 ); /* 1 1 0 0 1 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 1 1 0 0 1 0 . . */
              if ( binary[6] > 0) {
                /* 1 1 0 0 1 0 1 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 1 1 0 0 1 0 1 1 */
                }
                else {
                  return( A11 ); /* 1 1 0 0 1 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 1 0 0 1 0 0 . */
                if ( binary[7] > 0) {
                  return( A12 ); /* 1 1 0 0 1 0 0 1 */
                }
                else {
                  return( A5 ); /* 1 1 0 0 1 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* binary[4] */
          else {
            /* 1 1 0 0 0 . . . */
            if ( binary[5] > 0 ) {
              /* 1 1 0 0 0 1 . . */
              if ( binary[6] > 0) {
                /* 1 1 0 0 0 1 1 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 1 1 0 0 0 1 1 1 */
                }
                else {
                  return( A12 ); /* 1 1 0 0 0 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 1 0 0 0 1 0 . */
                if ( binary[7] > 0) {
                  return( A11 ); /* 1 1 0 0 0 1 0 1 */
                }
                else {
                  return( A5 ); /* 1 1 0 0 0 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 1 1 0 0 0 0 . . */
              if ( binary[6] > 0) {
                /* 1 1 0 0 0 0 1 . */
                if ( binary[7] > 0) {
                  return( A10 ); /* 1 1 0 0 0 0 1 1 */
                }
                else {
                  return( A6 ); /* 1 1 0 0 0 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 1 0 0 0 0 0 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 1 1 0 0 0 0 0 1 */
                }
                else {
                  return( A2 ); /* 1 1 0 0 0 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* else binary[4] */
        } /* else binary[3] */
      } /* else binary[2] */
    } /* binary[1] */
    else {
      /* 1 0 . . . . . . */
      if ( binary[2] > 0 ) {
        /* 1 0 1 . . . . . */
        if ( binary[3] > 0 ) {
          /* 1 0 1 1 . . . . */
          if ( binary[4] > 0 ) {
            /* 1 0 1 1 1 . . . */
            if ( binary[5] > 0 ) {
              /* 1 0 1 1 1 1 . . */
              if ( binary[6] > 0) {
                /* 1 0 1 1 1 1 1 . */
                if ( binary[7] > 0) {
                  return( A1 ); /* 1 0 1 1 1 1 1 1 */
                }
                else {
                  return( A3 ); /* 1 0 1 1 1 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 0 1 1 1 1 0 . */
                if ( binary[7] > 0) {
                  return( A4 ); /* 1 0 1 1 1 1 0 1 */
                }
                else {
                  return( A6 ); /* 1 0 1 1 1 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 1 0 1 1 1 0 . . */
              if ( binary[6] > 0) {
                /* 1 0 1 1 1 0 1 . */
                if ( binary[7] > 0) {
                  return( A2 ); /* 1 0 1 1 1 0 1 1 */
                }
                else {
                  return( A5 ); /* 1 0 1 1 1 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 0 1 1 1 0 0 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 1 0 1 1 1 0 0 1 */
                }
                else {
                  return( A11 ); /* 1 0 1 1 1 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* binary[4] */
          else {
            /* 1 0 1 1 0 . . . */
            if ( binary[5] > 0 ) {
              /* 1 0 1 1 0 1 . . */
              if ( binary[6] > 0) {
                /* 1 0 1 1 0 1 1 . */
                if ( binary[7] > 0) {
                  return( A3 ); /* 1 0 1 1 0 1 1 1 */
                }
                else {
                  return( A7 ); /* 1 0 1 1 0 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 0 1 1 0 1 0 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 1 0 1 1 0 1 0 1 */
                }
                else {
                  return( A12 ); /* 1 0 1 1 0 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 1 0 1 1 0 0 . . */
              if ( binary[6] > 0) {
                /* 1 0 1 1 0 0 1 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 1 0 1 1 0 0 1 1 */
                }
                else {
                  return( A9 ); /* 1 0 1 1 0 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 0 1 1 0 0 0 . */
                if ( binary[7] > 0) {
                  return( A11 ); /* 1 0 1 1 0 0 0 1 */
                }
                else {
                  return( A5 ); /* 1 0 1 1 0 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* else binary[4] */
        } /* binary[3] */
        else {
          /* 1 0 1 0 . . . . */
          if ( binary[4] > 0 ) {
            /* 1 0 1 0 1 . . . */
            if ( binary[5] > 0 ) {
              /* 1 0 1 0 1 1 . . */
              if ( binary[6] > 0) {
                /* 1 0 1 0 1 1 1 . */
                if ( binary[7] > 0) {
                  return( A2 ); /* 1 0 1 0 1 1 1 1 */
                }
                else {
                  return( A5 ); /* 1 0 1 0 1 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 0 1 0 1 1 0 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 1 0 1 0 1 1 0 1 */
                }
                else {
                  return( A11 ); /* 1 0 1 0 1 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 1 0 1 0 1 0 . . */
              if ( binary[6] > 0) {
                /* 1 0 1 0 1 0 1 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 1 0 1 0 1 0 1 1 */
                }
                else {
                  return( A8 ); /* 1 0 1 0 1 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 0 1 0 1 0 0 . */
                if ( binary[7] > 0) {
                  return( A12 ); /* 1 0 1 0 1 0 0 1 */
                }
                else {
                  return( A5 ); /* 1 0 1 0 1 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* binary[4] */
          else {
            /* 1 0 1 0 0 . . . */
            if ( binary[5] > 0 ) {
              /* 1 0 1 0 0 1 . . */
              if ( binary[6] > 0) {
                /* 1 0 1 0 0 1 1 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 1 0 1 0 0 1 1 1 */
                }
                else {
                  return( A12 ); /* 1 0 1 0 0 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 0 1 0 0 1 0 . */
                if ( binary[7] > 0) {
                  return( A10 ); /* 1 0 1 0 0 1 0 1 */
                }
                else {
                  return( A6 ); /* 1 0 1 0 0 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 1 0 1 0 0 0 . . */
              if ( binary[6] > 0) {
                /* 1 0 1 0 0 0 1 . */
                if ( binary[7] > 0) {
                  return( A11 ); /* 1 0 1 0 0 0 1 1 */
                }
                else {
                  return( A5 ); /* 1 0 1 0 0 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 0 1 0 0 0 0 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 1 0 1 0 0 0 0 1 */
                }
                else {
                  return( A2 ); /* 1 0 1 0 0 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* else binary[4] */
        } /* else binary[3] */
      } /* binary[2] */
      else {
        /* 1 0 0 . . . . . */
        if ( binary[3] > 0 ) {
          /* 1 0 0 1 . . . . */
          if ( binary[4] > 0 ) {
            /* 1 0 0 1 1 . . . */
            if ( binary[5] > 0 ) {
              /* 1 0 0 1 1 1 . . */
              if ( binary[6] > 0) {
                /* 1 0 0 1 1 1 1 . */
                if ( binary[7] > 0) {
                  return( A3 ); /* 1 0 0 1 1 1 1 1 */
                }
                else {
                  return( A7 ); /* 1 0 0 1 1 1 1 0 */

                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 0 0 1 1 1 0 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 1 0 0 1 1 1 0 1 */
                }
                else {
                  return( A12 ); /* 1 0 0 1 1 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 1 0 0 1 1 0 . . */
              if ( binary[6] > 0) {
                /* 1 0 0 1 1 0 1 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 1 0 0 1 1 0 1 1 */
                }
                else {
                  return( A12 ); /* 1 0 0 1 1 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 0 0 1 1 0 0 . */
                if ( binary[7] > 0) {
                  return( A10 ); /* 1 0 0 1 1 0 0 1 */
                }
                else {
                  return( A6 ); /* 1 0 0 1 1 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* binary[4] */
          else {
            /* 1 0 0 1 0 . . . */
            if ( binary[5] > 0 ) {
              /* 1 0 0 1 0 1 . . */
              if ( binary[6] > 0) {
                /* 1 0 0 1 0 1 1 . */
                if ( binary[7] > 0) {
                  return( A7 ); /* 1 0 0 1 0 1 1 1 */
                }
                else {
                  return( A13 ); /* 1 0 0 1 0 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 0 0 1 0 1 0 . */
                if ( binary[7] > 0) {
                  return( A12 ); /* 1 0 0 1 0 1 0 1 */
                }
                else {
                  return( A7 ); /* 1 0 0 1 0 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 1 0 0 1 0 0 . . */
              if ( binary[6] > 0) {
                /* 1 0 0 1 0 0 1 . */
                if ( binary[7] > 0) {
                  return( A12 ); /* 1 0 0 1 0 0 1 1 */
                }
                else {
                  return( A7 ); /* 1 0 0 1 0 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 0 0 1 0 0 0 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 1 0 0 1 0 0 0 1 */
                }
                else {
                  return( A3 ); /* 1 0 0 1 0 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* else binary[4] */
        } /* binary[3] */
        else {
          /* 1 0 0 0 . . . . */
          if ( binary[4] > 0 ) {
            /* 1 0 0 0 1 . . . */
            if ( binary[5] > 0 ) {
              /* 1 0 0 0 1 1 . . */
              if ( binary[6] > 0) {
                /* 1 0 0 0 1 1 1 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 1 0 0 0 1 1 1 1 */
                }
                else {
                  return( A9 ); /* 1 0 0 0 1 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 0 0 0 1 1 0 . */
                if ( binary[7] > 0) {
                  return( A11 ); /* 1 0 0 0 1 1 0 1 */
                }
                else {
                  return( A5 ); /* 1 0 0 0 1 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 1 0 0 0 1 0 . . */
              if ( binary[6] > 0) {
                /* 1 0 0 0 1 0 1 . */
                if ( binary[7] > 0) {
                  return( A11 ); /* 1 0 0 0 1 0 1 1 */
                }
                else {
                  return( A5 ); /* 1 0 0 0 1 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 0 0 0 1 0 0 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 1 0 0 0 1 0 0 1 */
                }
                else {
                  return( A2 ); /* 1 0 0 0 1 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* binary[4] */
          else {
            /* 1 0 0 0 0 . . . */
            if ( binary[5] > 0 ) {
              /* 1 0 0 0 0 1 . . */
              if ( binary[6] > 0) {
                /* 1 0 0 0 0 1 1 . */
                if ( binary[7] > 0) {
                  return( A12 ); /* 1 0 0 0 0 1 1 1 */
                }
                else {
                  return( A7 ); /* 1 0 0 0 0 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 0 0 0 0 1 0 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 1 0 0 0 0 1 0 1 */
                }
                else {
                  return( A3 ); /* 1 0 0 0 0 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 1 0 0 0 0 0 . . */
              if ( binary[6] > 0) {
                /* 1 0 0 0 0 0 1 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 1 0 0 0 0 0 1 1 */
                }
                else {
                  return( A3 ); /* 1 0 0 0 0 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 1 0 0 0 0 0 0 . */
                if ( binary[7] > 0) {
                  return( A4 ); /* 1 0 0 0 0 0 0 1 */
                }
                else {
                  return( A1 ); /* 1 0 0 0 0 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* else binary[4] */
        } /* else binary[3] */
      } /* else binary[2] */
    } /* else binary[1] */
  } /* binary[0] */
  else {
    /* 0 . . . . . . . */
    if ( binary[1] > 0 ) {
      /* 0 1 . . . . . . */
      if ( binary[2] > 0 ) {
        /* 0 1 1 . . . . . */
        if ( binary[3] > 0 ) {
          /* 0 1 1 1 . . . . */
          if ( binary[4] > 0 ) {
            /* 0 1 1 1 1 . . . */
            if ( binary[5] > 0 ) {
              /* 0 1 1 1 1 1 . . */
              if ( binary[6] > 0) {
                /* 0 1 1 1 1 1 1 . */
                if ( binary[7] > 0) {
                  return( A1 ); /* 0 1 1 1 1 1 1 1 */
                }
                else {
                  return( A4 ); /* 0 1 1 1 1 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 1 1 1 1 1 0 . */
                if ( binary[7] > 0) {
                  return( A3 ); /* 0 1 1 1 1 1 0 1 */
                }
                else {
                  return( A6 ); /* 0 1 1 1 1 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 0 1 1 1 1 0 . . */
              if ( binary[6] > 0) {
                /* 0 1 1 1 1 0 1 . */
                if ( binary[7] > 0) {
                  return( A3 ); /* 0 1 1 1 1 0 1 1 */
                }
                else {
                  return( A6 ); /* 0 1 1 1 1 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 1 1 1 1 0 0 . */
                if ( binary[7] > 0) {
                  return( A7 ); /* 0 1 1 1 1 0 0 1 */
                }
                else {
                  return( A12 ); /* 0 1 1 1 1 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* binary[4] */
          else {
            /* 0 1 1 1 0 . . . */
            if ( binary[5] > 0 ) {
              /* 0 1 1 1 0 1 . . */
              if ( binary[6] > 0) {
                /* 0 1 1 1 0 1 1 . */
                if ( binary[7] > 0) {
                  return( A2 ); /* 0 1 1 1 0 1 1 1 */
                }
                else {
                  return( A6 ); /* 0 1 1 1 0 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 1 1 1 0 1 0 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 0 1 1 1 0 1 0 1 */
                }
                else {
                  return( A11 ); /* 0 1 1 1 0 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 0 1 1 1 0 0 . . */
              if ( binary[6] > 0) {
                /* 0 1 1 1 0 0 1 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 0 1 1 1 0 0 1 1 */
                }
                else {
                  return( A11 ); /* 0 1 1 1 0 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 1 1 1 0 0 0 . */
                if ( binary[7] > 0) {
                  return( A9 ); /* 0 1 1 1 0 0 0 1 */
                }
                else {
                  return( A5 ); /* 0 1 1 1 0 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* else binary[4] */
        } /* binary[3] */
        else {
          /* 0 1 1 0 . . . . */
          if ( binary[4] > 0 ) {
            /* 0 1 1 0 1 . . . */
            if ( binary[5] > 0 ) {
              /* 0 1 1 0 1 1 . . */
              if ( binary[6] > 0) {
                /* 0 1 1 0 1 1 1 . */
                if ( binary[7] > 0) {
                  return( A3 ); /* 0 1 1 0 1 1 1 1 */
                }
                else {
                  return( A6 ); /* 0 1 1 0 1 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 1 1 0 1 1 0 . */
                if ( binary[7] > 0) {
                  return( A7 ); /* 0 1 1 0 1 1 0 1 */
                }
                else {
                  return( A12 ); /* 0 1 1 0 1 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 0 1 1 0 1 0 . . */
              if ( binary[6] > 0) {
                /* 0 1 1 0 1 0 1 . */
                if ( binary[7] > 0) {
                  return( A7 ); /* 0 1 1 0 1 0 1 1 */
                }
                else {
                  return( A12 ); /* 0 1 1 0 1 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 1 1 0 1 0 0 . */
                if ( binary[7] > 0) {
                  return( A13 ); /* 0 1 1 0 1 0 0 1 */
                }
                else {
                  return( A7 ); /* 0 1 1 0 1 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* binary[4] */
          else {
            /* 0 1 1 0 0 . . . */
            if ( binary[5] > 0 ) {
              /* 0 1 1 0 0 1 . . */
              if ( binary[6] > 0) {
                /* 0 1 1 0 0 1 1 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 0 1 1 0 0 1 1 1 */
                }
                else {
                  return( A10 ); /* 0 1 1 0 0 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 1 1 0 0 1 0 . */
                if ( binary[7] > 0) {
                  return( A12 ); /* 0 1 1 0 0 1 0 1 */
                }
                else {
                  return( A6 ); /* 0 1 1 0 0 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 0 1 1 0 0 0 . . */
              if ( binary[6] > 0) {
                /* 0 1 1 0 0 0 1 . */
                if ( binary[7] > 0) {
                  return( A12 ); /* 0 1 1 0 0 0 1 1 */
                }
                else {
                  return( A6 ); /* 0 1 1 0 0 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 1 1 0 0 0 0 . */
                if ( binary[7] > 0) {
                  return( A7 ); /* 0 1 1 0 0 0 0 1 */
                }
                else {
                  return( A3 ); /* 0 1 1 0 0 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* else binary[4] */
        } /* else binary[3] */
      } /* binary[2] */
      else {
        /* 0 1 0 . . . . . */
        if ( binary[3] > 0 ) {
          /* 0 1 0 1 . . . . */
          if ( binary[4] > 0 ) {
            /* 0 1 0 1 1 . . . */
            if ( binary[5] > 0 ) {
              /* 0 1 0 1 1 1 . . */
              if ( binary[6] > 0) {
                /* 0 1 0 1 1 1 1 . */
                if ( binary[7] > 0) {
                  return( A2 ); /* 0 1 0 1 1 1 1 1 */
                }
                else {
                  return( A6 ); /* 0 1 0 1 1 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 1 0 1 1 1 0 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 0 1 0 1 1 1 0 1 */
                }
                else {
                  return( A11 ); /* 0 1 0 1 1 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 0 1 0 1 1 0 . . */
              if ( binary[6] > 0) {
                /* 0 1 0 1 1 0 1 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 0 1 0 1 1 0 1 1 */
                }
                else {
                  return( A10 ); /* 0 1 0 1 1 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 1 0 1 1 0 0 . */
                if ( binary[7] > 0) {
                  return( A12 ); /* 0 1 0 1 1 0 0 1 */
                }
                else {
                  return( A6 ); /* 0 1 0 1 1 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* binary[4] */
          else {
            /* 0 1 0 1 0 . . . */
            if ( binary[5] > 0 ) {
              /* 0 1 0 1 0 1 . . */
              if ( binary[6] > 0) {
                /* 0 1 0 1 0 1 1 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 0 1 0 1 0 1 1 1 */
                }
                else {
                  return( A12 ); /* 0 1 0 1 0 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 1 0 1 0 1 0 . */
                if ( binary[7] > 0) {
                  return( A8 ); /* 0 1 0 1 0 1 0 1 */
                }
                else {
                  return( A5 ); /* 0 1 0 1 0 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 0 1 0 1 0 0 . . */
              if ( binary[6] > 0) {
                /* 0 1 0 1 0 0 1 . */
                if ( binary[7] > 0) {
                  return( A11 ); /* 0 1 0 1 0 0 1 1 */
                }
                else {
                  return( A6 ); /* 0 1 0 1 0 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 1 0 1 0 0 0 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 0 1 0 1 0 0 0 1 */
                }
                else {
                  return( A2 ); /* 0 1 0 1 0 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* else binary[4] */
        } /* binary[3] */
        else {
          /* 0 1 0 0 . . . . */
          if ( binary[4] > 0 ) {
            /* 0 1 0 0 1 . . . */
            if ( binary[5] > 0 ) {
              /* 0 1 0 0 1 1 . . */
              if ( binary[6] > 0) {
                /* 0 1 0 0 1 1 1 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 0 1 0 0 1 1 1 1 */
                }
                else {
                  return( A11 ); /* 0 1 0 0 1 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 1 0 0 1 1 0 . */
                if ( binary[7] > 0) {
                  return( A9 ); /* 0 1 0 0 1 1 0 1 */
                }
                else {
                  return( A5 ); /* 0 1 0 0 1 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 0 1 0 0 1 0 . . */
              if ( binary[6] > 0) {
                /* 0 1 0 0 1 0 1 . */
                if ( binary[7] > 0) {
                  return( A12 ); /* 0 1 0 0 1 0 1 1 */
                }
                else {
                  return( A6 ); /* 0 1 0 0 1 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 1 0 0 1 0 0 . */
                if ( binary[7] > 0) {
                  return( A7 ); /* 0 1 0 0 1 0 0 1 */
                }
                else {
                  return( A3 ); /* 0 1 0 0 1 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* binary[4] */
          else {
            /* 0 1 0 0 0 . . . */
            if ( binary[5] > 0 ) {
              /* 0 1 0 0 0 1 . . */
              if ( binary[6] > 0) {
                /* 0 1 0 0 0 1 1 . */
                if ( binary[7] > 0) {
                  return( A11 ); /* 0 1 0 0 0 1 1 1 */
                }
                else {
                  return( A6 ); /* 0 1 0 0 0 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 1 0 0 0 1 0 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 0 1 0 0 0 1 0 1 */
                }
                else {
                  return( A2 ); /* 0 1 0 0 0 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 0 1 0 0 0 0 . . */
              if ( binary[6] > 0) {
                /* 0 1 0 0 0 0 1 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 0 1 0 0 0 0 1 1 */
                }
                else {
                  return( A4 ); /* 0 1 0 0 0 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 1 0 0 0 0 0 . */
                if ( binary[7] > 0) {
                  return( A3 ); /* 0 1 0 0 0 0 0 1 */
                }
                else {
                  return( A1 ); /* 0 1 0 0 0 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* else binary[4] */
        } /* else binary[3] */
      } /* else binary[2] */
    } /* binary[1] */
    else {
      /* 0 0 . . . . . . */
      if ( binary[2] > 0 ) {
        /* 0 0 1 . . . . . */
        if ( binary[3] > 0 ) {
          /* 0 0 1 1 . . . . */
          if ( binary[4] > 0 ) {
            /* 0 0 1 1 1 . . . */
            if ( binary[5] > 0 ) {
              /* 0 0 1 1 1 1 . . */
              if ( binary[6] > 0) {
                /* 0 0 1 1 1 1 1 . */
                if ( binary[7] > 0) {
                  return( A2 ); /* 0 0 1 1 1 1 1 1 */
                }
                else {
                  return( A6 ); /* 0 0 1 1 1 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 0 1 1 1 1 0 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 0 0 1 1 1 1 0 1 */
                }
                else {
                  return( A10 ); /* 0 0 1 1 1 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 0 0 1 1 1 0 . . */
              if ( binary[6] > 0) {
                /* 0 0 1 1 1 0 1 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 0 0 1 1 1 0 1 1 */
                }
                else {
                  return( A11 ); /* 0 0 1 1 1 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 0 1 1 1 0 0 . */
                if ( binary[7] > 0) {
                  return( A12 ); /* 0 0 1 1 1 0 0 1 */
                }
                else {
                  return( A6 ); /* 0 0 1 1 1 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* binary[4] */
          else {
            /* 0 0 1 1 0 . . . */
            if ( binary[5] > 0 ) {
              /* 0 0 1 1 0 1 . . */
              if ( binary[6] > 0) {
                /* 0 0 1 1 0 1 1 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 0 0 1 1 0 1 1 1 */
                }
                else {
                  return( A12 ); /* 0 0 1 1 0 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 0 1 1 0 1 0 . */
                if ( binary[7] > 0) {
                  return( A11 ); /* 0 0 1 1 0 1 0 1 */
                }
                else {
                  return( A6 ); /* 0 0 1 1 0 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 0 0 1 1 0 0 . . */
              if ( binary[6] > 0) {
                /* 0 0 1 1 0 0 1 . */
                if ( binary[7] > 0) {
                  return( A8 ); /* 0 0 1 1 0 0 1 1 */
                }
                else {
                  return( A5 ); /* 0 0 1 1 0 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 0 1 1 0 0 0 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 0 0 1 1 0 0 0 1 */
                }
                else {
                  return( A2 ); /* 0 0 1 1 0 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* else binary[4] */
        } /* binary[3] */
        else {
          /* 0 0 1 0 . . . . */
          if ( binary[4] > 0 ) {
            /* 0 0 1 0 1 . . . */
            if ( binary[5] > 0 ) {
              /* 0 0 1 0 1 1 . . */
              if ( binary[6] > 0) {
                /* 0 0 1 0 1 1 1 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 0 0 1 0 1 1 1 1 */
                }
                else {
                  return( A11 ); /* 0 0 1 0 1 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 0 1 0 1 1 0 . */
                if ( binary[7] > 0) {
                  return( A12 ); /* 0 0 1 0 1 1 0 1 */
                }
                else {
                  return( A6 ); /* 0 0 1 0 1 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 0 0 1 0 1 0 . . */
              if ( binary[6] > 0) {
                /* 0 0 1 0 1 0 1 . */
                if ( binary[7] > 0) {
                  return( A9 ); /* 0 0 1 0 1 0 1 1 */
                }
                else {
                  return( A5 ); /* 0 0 1 0 1 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 0 1 0 1 0 0 . */
                if ( binary[7] > 0) {
                  return( A7 ); /* 0 0 1 0 1 0 0 1 */
                }
                else {
                  return( A3 ); /* 0 0 1 0 1 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* binary[4] */
          else {
            /* 0 0 1 0 0 . . . */
            if ( binary[5] > 0 ) {
              /* 0 0 1 0 0 1 . . */
              if ( binary[6] > 0) {
                /* 0 0 1 0 0 1 1 . */
                if ( binary[7] > 0) {
                  return( A11 ); /* 0 0 1 0 0 1 1 1 */
                }
                else {
                  return( A6 ); /* 0 0 1 0 0 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 0 1 0 0 1 0 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 0 0 1 0 0 1 0 1 */
                }
                else {
                  return( A4 ); /* 0 0 1 0 0 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 0 0 1 0 0 0 . . */
              if ( binary[6] > 0) {
                /* 0 0 1 0 0 0 1 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 0 0 1 0 0 0 1 1 */
                }
                else {
                  return( A2 ); /* 0 0 1 0 0 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 0 1 0 0 0 0 . */
                if ( binary[7] > 0) {
                  return( A3 ); /* 0 0 1 0 0 0 0 1 */
                }
                else {
                  return( A1 ); /* 0 0 1 0 0 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* else binary[4] */
        } /* else binary[3] */
      } /* binary[2] */
      else {
        /* 0 0 0 . . . . . */
        if ( binary[3] > 0 ) {
          /* 0 0 0 1 . . . . */
          if ( binary[4] > 0 ) {
            /* 0 0 0 1 1 . . . */
            if ( binary[5] > 0 ) {
              /* 0 0 0 1 1 1 . . */
              if ( binary[6] > 0) {
                /* 0 0 0 1 1 1 1 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 0 0 0 1 1 1 1 1 */
                }
                else {
                  return( A12 ); /* 0 0 0 1 1 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 0 0 1 1 1 0 . */
                if ( binary[7] > 0) {
                  return( A11 ); /* 0 0 0 1 1 1 0 1 */
                }
                else {
                  return( A6 ); /* 0 0 0 1 1 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 0 0 0 1 1 0 . . */
              if ( binary[6] > 0) {
                /* 0 0 0 1 1 0 1 . */
                if ( binary[7] > 0) {
                  return( A11 ); /* 0 0 0 1 1 0 1 1 */
                }
                else {
                  return( A6 ); /* 0 0 0 1 1 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 0 0 1 1 0 0 . */
                if ( binary[7] > 0) {
                  return( A6 ); /* 0 0 0 1 1 0 0 1 */
                }
                else {
                  return( A4 ); /* 0 0 0 1 1 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* binary[4] */
          else {
            /* 0 0 0 1 0 . . . */
            if ( binary[5] > 0 ) {
              /* 0 0 0 1 0 1 . . */
              if ( binary[6] > 0) {
                /* 0 0 0 1 0 1 1 . */
                if ( binary[7] > 0) {
                  return( A9 ); /* 0 0 0 1 0 1 1 1 */
                }
                else {
                  return( A7 ); /* 0 0 0 1 0 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 0 0 1 0 1 0 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 0 0 0 1 0 1 0 1 */
                }
                else {
                  return( A3 ); /* 0 0 0 1 0 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 0 0 0 1 0 0 . . */
              if ( binary[6] > 0) {
                /* 0 0 0 1 0 0 1 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 0 0 0 1 0 0 1 1 */
                }
                else {
                  return( A3 ); /* 0 0 0 1 0 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 0 0 1 0 0 0 . */
                if ( binary[7] > 0) {
                  return( A2 ); /* 0 0 0 1 0 0 0 1 */
                }
                else {
                  return( A1 ); /* 0 0 0 1 0 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* else binary[4] */
        } /* binary[3] */
        else {
          /* 0 0 0 0 . . . . */
          if ( binary[4] > 0 ) {
            /* 0 0 0 0 1 . . . */
            if ( binary[5] > 0 ) {
              /* 0 0 0 0 1 1 . . */
              if ( binary[6] > 0) {
                /* 0 0 0 0 1 1 1 . */
                if ( binary[7] > 0) {
                  return( A8 ); /* 0 0 0 0 1 1 1 1 */
                }
                else {
                  return( A5 ); /* 0 0 0 0 1 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 0 0 0 1 1 0 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 0 0 0 0 1 1 0 1 */
                }
                else {
                  return( A2 ); /* 0 0 0 0 1 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 0 0 0 0 1 0 . . */
              if ( binary[6] > 0) {
                /* 0 0 0 0 1 0 1 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 0 0 0 0 1 0 1 1 */
                }
                else {
                  return( A2 ); /* 0 0 0 0 1 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 0 0 0 1 0 0 . */
                if ( binary[7] > 0) {
                  return( A3 ); /* 0 0 0 0 1 0 0 1 */
                }
                else {
                  return( A1 ); /* 0 0 0 0 1 0 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* binary[4] */
          else {
            /* 0 0 0 0 0 . . . */
            if ( binary[5] > 0 ) {
              /* 0 0 0 0 0 1 . . */
              if ( binary[6] > 0) {
                /* 0 0 0 0 0 1 1 . */
                if ( binary[7] > 0) {
                  return( A5 ); /* 0 0 0 0 0 1 1 1 */
                }
                else {
                  return( A3 ); /* 0 0 0 0 0 1 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 0 0 0 0 1 0 . */
                if ( binary[7] > 0) {
                  return( A2 ); /* 0 0 0 0 0 1 0 1 */
                }
                else {
                  return( A1 ); /* 0 0 0 0 0 1 0 0 */
                } /* else binary[7] */
              } /* else binary[6] */
            } /* binary[5] */
            else {
              /* 0 0 0 0 0 0 . . */
              if ( binary[6] > 0) {
                /* 0 0 0 0 0 0 1 . */
                if ( binary[7] > 0) {
                  return( A2 ); /* 0 0 0 0 0 0 1 1 */
                }
                else {
                  return( A1 ); /* 0 0 0 0 0 0 1 0 */
                } /* else binary[7] */
              } /* binary[6] */
              else {
                /* 0 0 0 0 0 0 0 . */
                if ( binary[7] > 0) {
                  return( A1 ); /* 0 0 0 0 0 0 0 1 */
                }
                else {
                  /* 0 0 0 0 0 0 0 0 */
                  return( 0.0 );
                } /* else binary[7] */
              } /* else binary[6] */
            } /* else binary[5] */
          } /* else binary[4] */
        } /* else binary[3] */
      } /* else binary[2] */
    } /* else binary[1] */
  } /* else binary[0] */

  fprintf( stderr, "%s: this should not be reached\n", proc );
  return( -1.0 );
}


