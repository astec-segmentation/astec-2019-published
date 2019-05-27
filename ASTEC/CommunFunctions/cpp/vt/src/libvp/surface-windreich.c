/****************************************************
 * surface-windreich.c -
 *
 * Copyright (c) INRIA 2018
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mar 27 nov 2018 13:38:27 CET
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

#include <surface-windreich.h>



/* neighborhood :
 *  0  1  2      9 10 11     18 19 20
 *  3  4  5  -  12 13 14  -  21 22 23
 *  6  7  8     15 16 17     24 25 26
 *
 * 6-neighbors are 4 10 12 14 16 22
 *
 * W7:
 *  4, 22 : 1 0 0 0 0 1
 * 10, 16 : 0 1 0 0 1 0
 * 12, 14 : 0 0 1 1 0 0
 *
 * W2:
 * others configurations than W7 that have 2 '1'
 *
 * W3:
 * 1 1 1 0 0 0
 * 1 1 0 1 0 0
 * 1 0 1 0 1 0
 * 1 0 0 1 1 0
 * 0 1 1 0 0 1
 * 0 1 0 1 0 1
 * 0 0 1 0 1 1
 * 0 0 0 1 1 1
 *
 * W4:
 * others configurations than W7 that have 3 '1'
 *
 * W8:
 * 0 1 1 1 1 0
 * 1 0 1 1 0 1
 * 1 1 0 0 1 1
 */

#define GET_VAL(tab,index) tab[index]

#define W1 0.894
#define W2 1.3409
#define W3 1.5879
#define W4 2
#define W5 8.0/3.0
#define W6 10.0/3.0
#define W7 1.79
#define W8 2.68
#define W9 4.08



/* Voxel-based surface area estimation: from theory to practice,
 * G. Windreich, N. Kiryati, G. Lohmann,
 * Pattern Recognition, 36(11):2531-2541, (2003)
 */

float _surface_windreich( int *neighbors )
{
    char *proc = "_surface_windreich";
    int label;

    label = GET_VAL(neighbors,13);

    if ( GET_VAL(neighbors,4) == label ) {
        /*  6-neighbors :  0 . . . . . */
        if ( GET_VAL(neighbors,10) == label ) {
            /*  6-neighbors :  0 0 . . . . */
            if ( GET_VAL(neighbors,12) == label ) {
                /*  6-neighbors :  0 0 0 . . . */
                if ( GET_VAL(neighbors,14) == label ) {
                    /*  6-neighbors :  0 0 0 0 . . */
                    if ( GET_VAL(neighbors,16) == label ) {
                        /*  6-neighbors :  0 0 0 0 0 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  0 0 0 0 0 0 */
                            return( 0.0 );
                        } else {
                            /*  6-neighbors :  0 0 0 0 0 1 */
                            return( W1 );
                        }
                    } else {
                        /*  6-neighbors :  0 0 0 0 1 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  0 0 0 0 1 0 */
                            return( W1 );
                        } else {
                            /*  6-neighbors :  0 0 0 0 1 1 */
                            return( W2 );
                        }
                    }
                } else {
                    /*  6-neighbors :  0 0 0 1 . . */
                    if ( GET_VAL(neighbors,16) == label ) {
                        /*  6-neighbors :  0 0 0 1 0 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  0 0 0 1 0 0 */
                            return( W1 );
                        } else {
                            /*  6-neighbors :  0 0 0 1 0 1 */
                            return( W2 );
                        }
                    } else {
                        /*  6-neighbors :  0 0 0 1 1 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  0 0 0 1 1 0 */
                            return( W2 );
                        } else {
                            /*  6-neighbors :  0 0 0 1 1 1 */
                            return( W3 );
                        }
                    }
                }
            } else {
                /*  6-neighbors :  0 0 1 . . . */
                if ( GET_VAL(neighbors,14) == label ) {
                    /*  6-neighbors :  0 0 1 0 . . */
                    if ( GET_VAL(neighbors,16) == label ) {
                        /*  6-neighbors :  0 0 1 0 0 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  0 0 1 0 0 0 */
                            return( W1 );
                        } else {
                            /*  6-neighbors :  0 0 1 0 0 1 */
                            return( W2 );
                        }
                    } else {
                        /*  6-neighbors :  0 0 1 0 1 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  0 0 1 0 1 0 */
                            return( W2 );
                        } else {
                            /*  6-neighbors :  0 0 1 0 1 1 */
                            return( W3 );
                        }
                    }
                } else {
                    /*  6-neighbors :  0 0 1 1 . . */
                    if ( GET_VAL(neighbors,16) == label ) {
                        /*  6-neighbors :  0 0 1 1 0 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  0 0 1 1 0 0 */
                            return( W7 );
                        } else {
                            /*  6-neighbors :  0 0 1 1 0 1 */
                            return( W4 );
                        }
                    } else {
                        /*  6-neighbors :  0 0 1 1 1 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  0 0 1 1 1 0 */
                            return( W4 );
                        } else {
                            /*  6-neighbors :  0 0 1 1 1 1 */
                            return( W5 );
                        }
                    }
                }
            }
        } else {
            /*  6-neighbors :  0 1 . . . . */
            if ( GET_VAL(neighbors,12) == label ) {
                /*  6-neighbors :  0 1 0 . . . */
                if ( GET_VAL(neighbors,14) == label ) {
                    /*  6-neighbors :  0 1 0 0 . . */
                    if ( GET_VAL(neighbors,16) == label ) {
                        /*  6-neighbors :  0 1 0 0 0 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  0 1 0 0 0 0 */
                            return( W1 );
                        } else {
                            /*  6-neighbors :  0 1 0 0 0 1 */
                            return( W2 );
                        }
                    } else {
                        /*  6-neighbors :  0 1 0 0 1 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  0 1 0 0 1 0 */
                            return( W7 );
                        } else {
                            /*  6-neighbors :  0 1 0 0 1 1 */
                            return( W4 );
                        }
                    }
                } else {
                    /*  6-neighbors :  0 1 0 1 . . */
                    if ( GET_VAL(neighbors,16) == label ) {
                        /*  6-neighbors :  0 1 0 1 0 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  0 1 0 1 0 0 */
                            return( W2 );
                        } else {
                            /*  6-neighbors :  0 1 0 1 0 1 */
                            return( W3 );
                        }
                    } else {
                        /*  6-neighbors :  0 1 0 1 1 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  0 1 0 1 1 0 */
                            return( W4 );
                        } else {
                            /*  6-neighbors :  0 1 0 1 1 1 */
                            return( W5 );
                        }
                    }
                }
            } else {
                /*  6-neighbors :  0 1 1 . . . */
                if ( GET_VAL(neighbors,14) == label ) {
                    /*  6-neighbors :  0 1 1 0 . . */
                    if ( GET_VAL(neighbors,16) == label ) {
                        /*  6-neighbors :  0 1 1 0 0 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  0 1 1 0 0 0 */
                            return( W2 );
                        } else {
                            /*  6-neighbors :  0 1 1 0 0 1 */
                            return( W3 );
                        }
                    } else {
                        /*  6-neighbors :  0 1 1 0 1 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  0 1 1 0 1 0 */
                            return( W4 );
                        } else {
                            /*  6-neighbors :  0 1 1 0 1 1 */
                            return( W5 );
                        }
                    }
                } else {
                    /*  6-neighbors :  0 1 1 1 . . */
                    if ( GET_VAL(neighbors,16) == label ) {
                        /*  6-neighbors :  0 1 1 1 0 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  0 1 1 1 0 0 */
                            return( W4 );
                        } else {
                            /*  6-neighbors :  0 1 1 1 0 1 */
                            return( W5 );
                        }
                    } else {
                        /*  6-neighbors :  0 1 1 1 1 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  0 1 1 1 1 0 */
                            return( W8 );
                        } else {
                            /*  6-neighbors :  0 1 1 1 1 1 */
                            return( W6 );
                        }
                    }
                }
            }
        }
    } else {
        /*  6-neighbors :  1 . . . . . */
        if ( GET_VAL(neighbors,10) == label ) {
            /*  6-neighbors :  1 0 . . . . */
            if ( GET_VAL(neighbors,12) == label ) {
                /*  6-neighbors :  1 0 0 . . . */
                if ( GET_VAL(neighbors,14) == label ) {
                    /*  6-neighbors :  1 0 0 0 . . */
                    if ( GET_VAL(neighbors,16) == label ) {
                        /*  6-neighbors :  1 0 0 0 0 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  1 0 0 0 0 0 */
                            return( W1 );
                        } else {
                            /*  6-neighbors :  1 0 0 0 0 1 */
                            return( W7 );
                        }
                    } else {
                        /*  6-neighbors :  1 0 0 0 1 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  1 0 0 0 1 0 */
                            return( W2 );
                        } else {
                            /*  6-neighbors :  1 0 0 0 1 1 */
                            return( W4 );
                        }
                    }
                } else {
                    /*  6-neighbors :  1 0 0 1 . . */
                    if ( GET_VAL(neighbors,16) == label ) {
                        /*  6-neighbors :  1 0 0 1 0 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  1 0 0 1 0 0 */
                            return( W2 );
                        } else {
                            /*  6-neighbors :  1 0 0 1 0 1 */
                            return( W4 );
                        }
                    } else {
                        /*  6-neighbors :  1 0 0 1 1 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  1 0 0 1 1 0 */
                            return( W3 );
                        } else {
                            /*  6-neighbors :  1 0 0 1 1 1 */
                            return( W5 );
                        }
                    }
                }
            } else {
                /*  6-neighbors :  1 0 1 . . . */
                if ( GET_VAL(neighbors,14) == label ) {
                    /*  6-neighbors :  1 0 1 0 . . */
                    if ( GET_VAL(neighbors,16) == label ) {
                        /*  6-neighbors :  1 0 1 0 0 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  1 0 1 0 0 0 */
                            return( W2 );
                        } else {
                            /*  6-neighbors :  1 0 1 0 0 1 */
                            return( W4 );
                        }
                    } else {
                        /*  6-neighbors :  1 0 1 0 1 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  1 0 1 0 1 0 */
                            return( W3 );
                        } else {
                            /*  6-neighbors :  1 0 1 0 1 1 */
                            return( W5 );
                        }
                    }
                } else {
                    /*  6-neighbors :  1 0 1 1 . . */
                    if ( GET_VAL(neighbors,16) == label ) {
                        /*  6-neighbors :  1 0 1 1 0 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  1 0 1 1 0 0 */
                            return( W4 );
                        } else {
                            /*  6-neighbors :  1 0 1 1 0 1 */
                            return( W8 );
                        }
                    } else {
                        /*  6-neighbors :  1 0 1 1 1 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  1 0 1 1 1 0 */
                            return( W5 );
                        } else {
                            /*  6-neighbors :  1 0 1 1 1 1 */
                            return( W6 );
                        }
                    }
                }
            }
        } else {
            /*  6-neighbors :  1 1 . . . . */
            if ( GET_VAL(neighbors,12) == label ) {
                /*  6-neighbors :  1 1 0 . . . */
                if ( GET_VAL(neighbors,14) == label ) {
                    /*  6-neighbors :  1 1 0 0 . . */
                    if ( GET_VAL(neighbors,16) == label ) {
                        /*  6-neighbors :  1 1 0 0 0 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  1 1 0 0 0 0 */
                            return( W2 );
                        } else {
                            /*  6-neighbors :  1 1 0 0 0 1 */
                            return( W4 );
                        }
                    } else {
                        /*  6-neighbors :  1 1 0 0 1 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  1 1 0 0 1 0 */
                            return( W4 );
                        } else {
                            /*  6-neighbors :  1 1 0 0 1 1 */
                            return( W8 );
                        }
                    }
                } else {
                    /*  6-neighbors :  1 1 0 1 . . */
                    if ( GET_VAL(neighbors,16) == label ) {
                        /*  6-neighbors :  1 1 0 1 0 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  1 1 0 1 0 0 */
                            return( W3 );
                        } else {
                            /*  6-neighbors :  1 1 0 1 0 1 */
                            return( W5 );
                        }
                    } else {
                        /*  6-neighbors :  1 1 0 1 1 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  1 1 0 1 1 0 */
                            return( W5 );
                        } else {
                            /*  6-neighbors :  1 1 0 1 1 1 */
                            return( W6 );
                        }
                    }
                }
            } else {
                /*  6-neighbors :  1 1 1 . . . */
                if ( GET_VAL(neighbors,14) == label ) {
                    /*  6-neighbors :  1 1 1 0 . . */
                    if ( GET_VAL(neighbors,16) == label ) {
                        /*  6-neighbors :  1 1 1 0 0 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  1 1 1 0 0 0 */
                            return( W3 );
                        } else {
                            /*  6-neighbors :  1 1 1 0 0 1 */
                            return( W5 );
                        }
                    } else {
                        /*  6-neighbors :  1 1 1 0 1 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  1 1 1 0 1 0 */
                            return( W5 );
                        } else {
                            /*  6-neighbors :  1 1 1 0 1 1 */
                            return( W6 );
                        }
                    }
                } else {
                    /*  6-neighbors :  1 1 1 1 . . */
                    if ( GET_VAL(neighbors,16) == label ) {
                        /*  6-neighbors :  1 1 1 1 0 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  1 1 1 1 0 0 */
                            return( W5 );
                        } else {
                            /*  6-neighbors :  1 1 1 1 0 1 */
                            return( W6 );
                        }
                    } else {
                        /*  6-neighbors :  1 1 1 1 1 . */
                        if ( GET_VAL(neighbors,22) == label ) {
                            /*  6-neighbors :  1 1 1 1 1 0 */
                            return( W6 );
                        } else {
                            /*  6-neighbors :  1 1 1 1 1 1 */
                            return( W9 );
                        }
                    }
                }
            }
        }
    }

   fprintf( stderr, "%s: this should not be reached\n", proc );
   return( -1.0 );
}
