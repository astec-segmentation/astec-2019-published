
#include <math.h>

#include <chunks.h>

#include <vtmalloc.h>
#include <vt_unix.h>
#include <vt_copy.h>
#include <vt_inrimage.h>

#include <vt_removeLine.h>





static int _debug_ = 0;

static int _verbose_ = 1;

void VT_SetVerboseInVtRemoveLine( int v )
{
  _verbose_ = v;
}

void VT_IncrementVerboseInVtRemoveLine(  )
{
  _verbose_ ++;
}

void VT_DecrementVerboseInVtRemoveLine(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}





/************************************************************
 *
 *
 *
 ************************************************************/





void initRemoveLineParameter( typeRemoveLineParameter *p )
{
    p->method = _GLOBAL_;
    p->contrastSignificantFraction = 0.6;
    p->yRejectedFraction = 0.1;
    p->xzKeptFraction = 1.0 - 3.1415926 / 4.0;
    p->automatedChoices = 0;
}



void fprintfRemoveLineParameter( FILE *f, typeRemoveLineParameter *p )
{
    fprintf( f, "- method = " );
    switch( p->method ) {
    default : fprintf( f, "unknown method\n" ); break;
    case _LOCAL_ : fprintf( f, "_LOCAL_\n" ); break;
    case _REGIONAL_ : fprintf( f, "_REGIONAL_\n" ); break;
    case _GLOBAL_ : fprintf( f, "_GLOBAL_\n" ); break;
    }
    fprintf( f, "- contrastSignificantFraction = %f\n", p->contrastSignificantFraction );
    fprintf( f, "- yRejectedFraction = %f\n", p->yRejectedFraction );
    fprintf( f, "- xzKeptFraction = %f\n", p->xzKeptFraction );
    fprintf( f, "- automatedChoices = %d\n", p->automatedChoices );
}





/************************************************************
 *
 *
 *
 ************************************************************/





static void _initCorrection( typeCorrection *c )
{
    c->y = -1;
    c->a = 1.0;
    c->b = 0.0;
}



void VT_InitCorrectionList( typeCorrectionList *l )
{
    l->data = (typeCorrection*)NULL;
    l->n_data = 0;
    l->n_allocated_data = 0;
}



void VT_FreeCorrectionList( typeCorrectionList *l )
{
    if ( l->data != (typeCorrection*)NULL )
        vtfree( l->data );
    VT_InitCorrectionList( l );
}



static int _allocCorrectionList( typeCorrectionList *l, int size )
{
    char *proc = "_allocCorrectionList";
    int i;

    l->data = (typeCorrection*)vtmalloc( size*sizeof(typeCorrection), "l->data", proc );
    if ( l->data == (typeCorrection*)NULL ) {
        if ( _verbose_ )
            fprintf( stderr, "%s: allocation error\n", proc );
        return( -1 );
    }

    for ( i=0; i<size; i++ )
        _initCorrection( &(l->data[i]) );

    l->n_data = 0;
    l->n_allocated_data = size;

    return( 1 );
}



static int _corrections_to_be_allocated_ = 50;

static int _addCorrectionToCorrectionList( typeCorrectionList *l, typeCorrection *c )
{
  char *proc = "_addCorrectionToCorrectionList";
  int i, s =  l->n_allocated_data;
  typeCorrection *data;

  if ( l->n_data == l->n_allocated_data ) {
    s += _corrections_to_be_allocated_;
    data  = (typeCorrection*)vtmalloc( s*sizeof(typeCorrection), "data", proc );
    if ( data == (typeCorrection*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    }
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(typeCorrection) );
      vtfree( l->data );
    }
    for ( i=l->n_data; i<l->n_allocated_data; i++ )
        _initCorrection( &(data[i]) );
    l->n_allocated_data = s;
    l->data = data;
  }

  l->data[l->n_data] = *c;
  l->n_data ++;

  return( 1 );
}



int VT_ReadCorrectionList( char *name, typeCorrectionList *l )
{
    char *proc = "VT_ReadCorrectionList";
    FILE *f;
    typeCorrection c;
    int i, ret;
    int line = 0;

    f = fopen( name, "r" );
    if ( f == (FILE*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to open file '%s'\n", proc, name );
      return( -1 );
    }

    while ( (ret = fscanf( f, "%d %f %f", &(c.y), &(c.a), &(c.b) )) != EOF ) {
        line ++;
        if ( ret != 3 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: error at line %d in file '%s'\n", proc, line, name );
        }
        else {
            if ( _addCorrectionToCorrectionList( l, &c ) != 1 ) {
                fclose( f );
                if ( _verbose_ )
                    fprintf( stderr, "%s: error when adding line %d from file '%s'\n", proc, line, name );
                return( -1 );
            }
        }
    }

    fclose( f );

    if ( _verbose_ >= 2 ) {
      fprintf( stderr, "\n" );
      fprintf( stderr, "- %4d flagged XZ slices where changes occur\n", l->n_data );
      for ( i=0; i<l->n_data; i++ ) {
        fprintf( stderr, "#xz=%4d  -  a=%6.2f b=%6.2f\n",
                 l->data[i].y, l->data[i].a, l->data[i].b );
      }
      fprintf( stderr, "\n" );
    }

    return( 1 );
}



static int _fprintfCorrectionList( FILE *f, typeCorrectionList *l )
{
    char *proc = "_fprintfCorrectionList";
    int i;

    for ( i=0; i<l->n_data; i++ ) {
        if ( fprintf( f, "%d %f %f\n", l->data[i].y, l->data[i].a, l->data[i].b ) < 0 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: error when writing line %d\n", proc, i+1 );
            return( -1 );
        }
    }
    return( 1 );
}



int VT_WriteCorrectionList( char *name, typeCorrectionList *l )
{
    char *proc = "VT_WriteCorrectionList";
    FILE *f;

    f = fopen( name, "w" );
    if ( f == (FILE*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to open file '%s'\n", proc, name );
      return( -1 );
    }

    if ( _fprintfCorrectionList( f, l ) != 1 ) {
        fclose( f );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when writing file '%s'\n", proc, name );
        return( -1 );
    }

    fclose( f );

    return( 1 );
}





/************************************************************
 *
 * contrast measure between adjacent XZ sections
 *
 ************************************************************/





typedef struct _XZcontrast {
    int nSupToPrevious;
    int nSupToNext;
} _XZcontrast;



static void _initXZcontrast( _XZcontrast *m )
{
    m->nSupToPrevious = 0;
    m->nSupToNext = 0;
}



typedef struct _XZcontrastList {
  _XZcontrast *data;
  int n_allocated_data;
} _XZcontrastList;



static void _initXZcontrastList( _XZcontrastList *l )
{
    l->data = (_XZcontrast*)NULL;
    l->n_allocated_data = 0;
}



static void _freeXZcontrastList( _XZcontrastList *l )
{
    if ( l->data != (_XZcontrast*)NULL )
        vtfree( l->data );
    _initXZcontrastList( l );
}



static int _allocXZcontrastList( _XZcontrastList *l, int size )
{
    char * proc = "_allocZcontrastMeasureList";
    int i;

    l->data  = (_XZcontrast*)vtmalloc( size*sizeof(_XZcontrast), "l->data", proc );
    if ( l->data == (_XZcontrast*)NULL ) {
        if ( _verbose_ )
            fprintf( stderr, "%s: allocation error\n", proc );
        return( -1 );
    }

    for ( i=0; i<size; i++ )
        _initXZcontrast( &(l->data[i]) );

    l->n_allocated_data = size;

    return( 1 );
}


static void _fprintfXZcontrastList( FILE *f, _XZcontrastList *l ) __attribute__ ((unused));
static void _fprintfXZcontrastList( FILE *f, _XZcontrastList *l )
{
    int y;

    for ( y=0; y<l->n_allocated_data; y++ ) {
        fprintf( f, "#xz=%4d  -  #[i(y-1)<i(y)] = %7d - #[i(y-1)<i(y)] = %7d\n",
                 y, l->data[y].nSupToPrevious, l->data[y].nSupToNext );
    }
}



static void _fprintfXZcontrastEvolution( FILE *f, _XZcontrastList *l, int threshold, int nmeasures )
{
    int y;

    for ( y=0; y<l->n_allocated_data; y++ ) {
        if ( y%nmeasures == 0 ) {
           fprintf( f, "\n" );
           fprintf( f, "[%4d] ", y );
        }
        if ( y == 0 ) fprintf( f, "." );
        else {
          if ( l->data[y].nSupToPrevious > threshold ) {
              if ( l->data[y-1].nSupToNext > threshold )
                  fprintf( f, "*" );
              else
                  fprintf( f, "+" );
          }
          else if ( l->data[y-1].nSupToNext > threshold )
              fprintf( f, "-" );
          else
              fprintf( f, "." );
        }
    }
    fprintf( f, "\n" );
    fprintf( f, "\n" );
}



typedef struct _computeXZcontrastParam {
    vt_image *image;
    _XZcontrastList *list;
} _computeXZcontrastParam;



#define _COMPUTEXZCONTRAST( TYPE ) {                           \
    TYPE ***theArray = (TYPE***)image->array;                  \
    for ( y=first; y<=last; y++ ) {                            \
        m = &(list->data[y]);                                  \
        m->nSupToPrevious = 0;                                 \
        m->nSupToNext = 0;                                     \
        if ( y == 0 ) {                                        \
            for ( z=0; z<image->dim.z; z++ )                   \
            for ( x=0; x<image->dim.x; x++ ) {                 \
                if ( theArray[z][y][x] > theArray[z][y+1][x] ) \
                    m->nSupToNext ++;                          \
            }                                                  \
        }                                                      \
        else if ( y == image->dim.y-1 ) {                      \
            for ( z=0; z<image->dim.z; z++ )                   \
            for ( x=0; x<image->dim.x; x++ ) {                 \
                if ( theArray[z][y][x] > theArray[z][y-1][x] ) \
                    m->nSupToPrevious ++;                      \
            }                                                  \
        }                                                      \
        else {                                                 \
            for ( z=0; z<image->dim.z; z++ )                   \
            for ( x=0; x<image->dim.x; x++ ) {                 \
                if ( theArray[z][y][x] > theArray[z][y+1][x] ) \
                    m->nSupToNext ++;                          \
                if ( theArray[z][y][x] > theArray[z][y-1][x] ) \
                    m->nSupToPrevious ++;                      \
            }                                                  \
        }                                                      \
    }                                                          \
}



static void *_computeXZcontrastProcedure( void *par )
{
    char *proc = "_computeXZcontrastProcedure";

    typeChunk *chunk = (typeChunk *)par;
    void *parameter = chunk->parameters;
    size_t first = chunk->first;
    size_t last = chunk->last;

    _computeXZcontrastParam *p = (_computeXZcontrastParam*)parameter;
    vt_image *image = p->image;
    _XZcontrastList *list = p->list;

    _XZcontrast *m;
    size_t x, y, z;


    switch( image->type ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such image type not handled\n", proc );
      chunk->ret = 0;
      return( (void*)NULL );
      break;
    case SSHORT :
      _COMPUTEXZCONTRAST( s16 )
      break;
    case USHORT :
      _COMPUTEXZCONTRAST( u16 )
      break;
    }

    chunk->ret = 1;
    return( (void*)NULL );
}



static int _computeXZcontrast( vt_image *image, _XZcontrastList *list )
{
    char *proc = "_computeXZcontrast";
    typeChunks chunks;
    _computeXZcontrastParam parameters;
    int n;


    initChunks( &chunks );
    if ( buildChunks( &chunks, 0, image->dim.y-1, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when building chunks\n", proc );
      return( -1 );
    }

    parameters.image = image;
    parameters.list = list;
    for ( n=0; n<chunks.n_allocated_chunks; n++ )
      chunks.data[n].parameters = (void*)(&parameters);

    if ( processChunks( &_computeXZcontrastProcedure, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute contrast measures\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    freeChunks( &chunks );

    return( 1 );
}





/************************************************************
 *
 *
 *
 ************************************************************/




typedef struct _XZstatistics {
    float intensityMean;
    float intensityStddev;
    int flagged;

    float aToPrevious;
    float bToPrevious;
    float aToNext;
    float bToNext;

    int dToPrevious;
    int dToNext;

    float a;
    float b;
} _XZstatistics;



static void _initXZtatistics( _XZstatistics *s )
{
    s->intensityMean = 0.0;
    s->intensityStddev = 0.0;
    s->flagged = 0;

    s->aToPrevious = 1.0;
    s->bToPrevious = 0.0;
    s->aToNext = 1.0;
    s->bToNext = 0.0;

    s->dToPrevious = 0;
    s->dToNext = 0;

    s->a = 1.0;
    s->b = 0.0;
}



typedef struct _XZstatisticsList {
    _XZstatistics *data;
    int n_allocated_data;
} _XZstatisticsList;



static void _initXZstatisticsList( _XZstatisticsList *l )
{
    l->data = (_XZstatistics*)NULL;
    l->n_allocated_data = 0;
}



static void _freeXZstatisticsList( _XZstatisticsList *l )
{
    if ( l->data != (_XZstatistics*)NULL )
        vtfree( l->data );
    _initXZstatisticsList( l );
}



static int _allocXZstatisticsList( _XZstatisticsList *l, int size )
{
    char * proc = "_allocZstatisticsList";
    int i;

    l->data  = (_XZstatistics*)vtmalloc( size*sizeof(_XZstatistics), "l->data", proc );
    if ( l->data == (_XZstatistics*)NULL ) {
        if ( _verbose_ )
            fprintf( stderr, "%s: allocation error\n", proc );
        return( -1 );
    }

    for ( i=0; i<size; i++ )
        _initXZtatistics( &(l->data[i]) );

    l->n_allocated_data = size;

    return( 1 );
}



static void _fprintfXZstatisticsList( FILE *f, _XZstatisticsList *l ) __attribute__ ((unused));
static void _fprintfXZstatisticsList( FILE *f, _XZstatisticsList *l )
{
    int y;

    for ( y=0; y<l->n_allocated_data; y++ ) {
        fprintf( f, "#xz=%4d  -  mean = %f , std dev = %f , flag = %d\n",
                 y, l->data[y].intensityMean, l->data[y].intensityStddev, l->data[y].flagged );
    }
}



static int _copyXZstatisticsListToCorrectionList( _XZstatisticsList *theList,
                                                   typeCorrectionList *resList )
{
    char *proc = "_copyXZstatisticsListToCorrectionList";
    int n, i;
    typeCorrection c;

    for ( i=0, n=0; i<theList->n_allocated_data; i++ ) {
        if ( theList->data[i].flagged == 0 ) continue;
        n++;
    }

    if ( n == 0 ) return( 0 );

    if ( _allocCorrectionList( resList, n ) < 0 ) {
        if ( _verbose_ )
            fprintf( stderr, "%s: allocation error\n", proc );
        return( -1 );
    }

    for ( i=0, n=0; i<theList->n_allocated_data; i++ ) {
        if ( theList->data[i].flagged == 0 ) continue;
        c.y = i;
        c.a = theList->data[i].a;
        c.b = theList->data[i].b;
        if ( _addCorrectionToCorrectionList( resList, &c ) != 1 ) {
            VT_FreeCorrectionList( resList );
            if ( _verbose_ )
                fprintf( stderr, "%s: error when adding correction #%d\n", proc, n );
            return( -1 );
        }
        n++;
    }

    return( n );
}


typedef struct _computeXZstatisticsParam {
    vt_image *image;
    _XZstatisticsList *list;
} _computeXZstatisticsParam;



#define _COMPUTEXZSTATISTICS( TYPE ) {                         \
    TYPE ***theArray = (TYPE***)image->array;                  \
    for ( y=first; y<=last; y++ ) {                            \
        s = &(list->data[y]);                                  \
        if ( s->flagged == 0                                   \
             && ( y > 0 && list->data[y-1].flagged == 0 )      \
             && ( y < image->dim.y-1 && list->data[y+1].flagged == 0 ) ) \
        continue;                                              \
        for ( sum=0.0, z=0; z<image->dim.z; z++ )              \
        for ( x=0; x<image->dim.x; x++ ) {                     \
            sum += theArray[z][y][x];                          \
        }                                                      \
        m = sum / (double)(image->dim.x * image->dim.z);       \
        s->intensityMean = m;                                  \
        for ( sum=0.0, z=0; z<image->dim.z; z++ )              \
        for ( x=0; x<image->dim.x; x++ ) {                     \
            sum += (theArray[z][y][x] - m) * (theArray[z][y][x] - m); \
        }                                                      \
        sum /= (double)(image->dim.x * image->dim.z);          \
        s->intensityStddev = sqrt( sum );                      \
    }                                                          \
}



static void *_computeXZstatisticsProcedure( void *par )
{
    char *proc = "_computeXZstatisticsProcedure";

    typeChunk *chunk = (typeChunk *)par;
    void *parameter = chunk->parameters;
    size_t first = chunk->first;
    size_t last = chunk->last;

    _computeXZstatisticsParam *p = (_computeXZstatisticsParam*)parameter;
    vt_image *image = p->image;
    _XZstatisticsList *list = p->list;

    _XZstatistics *s;
    size_t x, y, z;
    double sum, m;


    switch( image->type ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such image type not handled\n", proc );
      chunk->ret = 0;
      return( (void*)NULL );
      break;
    case SSHORT :
      _COMPUTEXZSTATISTICS( s16 )
      break;
    case USHORT :
      _COMPUTEXZSTATISTICS( u16 )
      break;
    }

    chunk->ret = 1;
    return( (void*)NULL );
}



static int _computeXZstatistics( vt_image *image, _XZstatisticsList *list )
{
    char *proc = "_computeXZcontrast";
    typeChunks chunks;
    _computeXZstatisticsParam parameters;
    int n;


    initChunks( &chunks );
    if ( buildChunks( &chunks, 0, image->dim.y-1, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when building chunks\n", proc );
      return( -1 );
    }

    parameters.image = image;
    parameters.list = list;
    for ( n=0; n<chunks.n_allocated_chunks; n++ )
      chunks.data[n].parameters = (void*)(&parameters);

    if ( processChunks( &_computeXZstatisticsProcedure, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute contrast measures\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    freeChunks( &chunks );

    return( 1 );
}





/************************************************************
 *
 * slice selection based on contrast measure
 *
 ************************************************************/





static int _selectXZslicesLocally( _XZcontrastList *clist, _XZstatisticsList *slist,
                           int threshold )
{
    int i, nslices;

    for ( i=0; i<slist->n_allocated_data; i++ )
            slist->data[i].flagged = 0;

    for ( nslices=0, i=0; i<clist->n_allocated_data; i++ ) {
      if ( clist->data[i].nSupToPrevious > threshold ||
           clist->data[i].nSupToNext > threshold ) {
        slist->data[i].flagged = 1;
        nslices ++;
      }
    }

    return( nslices );
}



static int _selectXZsliceRegionally( _XZcontrastList *clist, _XZstatisticsList *slist,
                           int threshold, int s1, int s2 )
{
    int i;
    int increase, decrease, count, valid;
    int nslices;

    for ( i=0; i<slist->n_allocated_data; i++ )
        slist->data[i].flagged = 0;

    /* forward propagation
     */
    increase = 0;
    decrease = 0;
    count = 0;

    for ( i=1; i<clist->n_allocated_data-1; i++ ) {

        /* is there an increase of intensity in [i-1,i]?
         */
        if ( clist->data[i].nSupToPrevious > threshold ) {
            if ( increase == 0 ) {
              if ( decrease == 0 ) count ++;
              increase = 1;
              decrease = 0;
            }
        }

        /* is there an decrease of intensity in [i-1,i]?
         */
        else if ( clist->data[i-1].nSupToNext > threshold ) {
            if ( decrease == 0 ) {
              increase = 0;
              decrease = 1;
            }
        }

        /* neither increase nor decrease => plateau
         */
        else {
          if ( increase == 1 ) increase = 0;
          if ( decrease == 1 ) {
              decrease = 0;
              count --;
          }
        }
        slist->data[i].flagged += count;
    }


    /* backward progation
     */
    increase = 0;
    decrease = 0;
    count = 0;

    for ( i=clist->n_allocated_data-2; i>0; i-- ) {

        /* is there an increase of intensity in [i+1,i]?
         */
        if ( clist->data[i].nSupToNext > threshold ) {
            if ( increase == 0 ) {
              if ( decrease == 0 ) count ++;
              increase = 1;
              decrease = 0;
            }
        }
        /* is there an decrease of intensity in [i+1,i]?
         */
        else if ( clist->data[i+1].nSupToPrevious > threshold ) {
            if ( decrease == 0 ) {
              increase = 0;
              decrease = 1;
            }
        }
        else {
          if ( increase == 1 ) increase = 0;
          if ( decrease == 1 ) {
              decrease = 0;
              count --;
          }
        }
        slist->data[i].flagged += count;
    }

    /* marked points to be corrected and count them
     */
    nslices = 0;
    valid = 0;
    for ( i=0; i<slist->n_allocated_data; i++ ) {
        if ( slist->data[i].flagged >= 2 ) {
            slist->data[i].flagged = 1;
            nslices++;
            if ( i == s1 || i == s2 ) valid = 1;
        }
        else {
            slist->data[i].flagged = 0;
        }
    }

    if ( s1 >= 0 && s2 >= 0 ) {
      if ( valid ) return( nslices );
      return( 0 );
    }

    return( nslices );
}



static int _selectThresholdRegionally( _XZcontrastList *clist, _XZstatisticsList *slist,
                                            int minThreshold, int maxThreshold, int nThreshold )
{
  char *proc = "_selectThresholdRegionally";

  double dt;
  int threshold, i, nslices;
  int maxchgs, s1, s2;
  int nslicesMin = 0;
  int thresholdMin = 0;


  /* slices to be included
   * => slices of maximal changes
   */
  s1 = s2 = -1;
  maxchgs = 0;
  for ( i=1; i<slist->n_allocated_data; i++ ) {
      if ( maxchgs < clist->data[i].nSupToPrevious ) {
          maxchgs = clist->data[i].nSupToPrevious;
          s1 = i-1;
          s2 = i;
      }
      if ( maxchgs < clist->data[i].nSupToNext ) {
          maxchgs = clist->data[i].nSupToNext;
          s1 = i;
          s2 = i+1;
      }
  }

  if ( _verbose_ >= 2 )
    fprintf( stdout, "   %s: max changes detected in [%d,%d]\n", proc, s1, s2 );


  dt = (double)(maxThreshold - minThreshold) / ((double)(nThreshold-1));


  /* the best threshold
   * => minimal number of slices to be corrected
   * && a slice of maximal changes must be corrected
   */
  for ( i=0; i<nThreshold; i++ ) {
      if ( i == 0 )
          threshold = minThreshold;
      else if ( i == nThreshold-1 )
          threshold = maxThreshold;
      else
          threshold = (int)(minThreshold + i * dt + 0.5);

      nslices = _selectXZsliceRegionally( clist, slist, threshold, s1, s2 );

      if ( _verbose_ >= 2 )
        fprintf( stdout, "   %s: threshold %6d -> %4d changes\n", proc, threshold, nslices );

      if ( i == 0 ) {
        nslicesMin = nslices;
        thresholdMin = threshold;
      }
      else {
        if ( nslicesMin > nslices && nslices > 0 ) {
            nslicesMin = nslices;
            thresholdMin = threshold;
        }
      }
  }

  if ( _verbose_ >= 2 )
    fprintf( stdout, "   %s: best threshold %6d -> %4d changes\n", proc, thresholdMin, nslicesMin );

  return( thresholdMin );
}





/************************************************************
 *
 *
 *
 ************************************************************/





typedef struct _Ystatistics {
  int x;
  int z;

  float intensityMean;
  float intensityRobustMean;

  int nExcludedPoint;
  int *excludedPoint;
} _Ystatistics;



static void _initYstatistics( _Ystatistics *s )
{
    s->x = -1;
    s->z = -1;

    s->intensityMean = 0;
    s->intensityRobustMean = 0;

    s->nExcludedPoint = 0;
    s->excludedPoint = (int*)NULL;
}



static void _freeYstatistics( _Ystatistics *s )
{
    if ( s->excludedPoint != (int*)NULL )
        vtfree( s->excludedPoint );
    _initYstatistics( s );
}



static int _allocYstatistics( _Ystatistics *s, int n ) __attribute__ ((unused));
static int _allocYstatistics( _Ystatistics *s, int n )
{
    char *proc = "_allocYstatistics";
    int i;
    s->excludedPoint = (int*)vtmalloc( n*sizeof(int), "s->excludedPoint", proc );
    if ( s->excludedPoint == (int*)NULL ) {
       if ( _verbose_ )
           fprintf( stderr, "%s: allocation failed\n", proc );
       return( -1 );
    }
    for ( i=0; i<n; i++ ) s->excludedPoint[i] = -1;
    s->nExcludedPoint = n;
    return( 1 );
}



static int _compareYstatisticsMean( const void *x1, const void *x2 )
{
    _Ystatistics **v1 = (_Ystatistics**) x1;
    _Ystatistics **v2 = (_Ystatistics**) x2;
    /* to sort in increasing order
     */
    if ( (*v1)->intensityMean < (*v2)->intensityMean )
        return( -1 );
    else if ( (*v1)->intensityMean > (*v2)->intensityMean )
        return( 1 );
    return( 0 );
}



static int _compareYstatisticsRobustMean( const void *x1, const void *x2 ) __attribute__ ((unused));
static int _compareYstatisticsRobustMean( const void *x1, const void *x2 )
{
    _Ystatistics **v1 = (_Ystatistics**) x1;
    _Ystatistics **v2 = (_Ystatistics**) x2;
    /* to sort in increasing order
     */
    if ( (*v1)->intensityRobustMean < (*v2)->intensityRobustMean )
        return( -1 );
    else if ( (*v1)->intensityRobustMean > (*v2)->intensityRobustMean )
        return( 1 );
    return( 0 );
}



typedef struct _YstatisticsList {

  _Ystatistics *data;
  _Ystatistics **pointer;
  int n_allocated_data;
  int n_selected_data;

} _YstatisticsList;



static void _initYstatisticsList( _YstatisticsList *l )
{
    l->data = (_Ystatistics*)NULL;
    l->pointer = (_Ystatistics**)NULL;
    l->n_allocated_data = 0;
    l->n_selected_data = 0;
}



static void _freeYstatisticsList( _YstatisticsList *l )
{
    int i;

    if ( l->data != (_Ystatistics*)NULL ) {
        for ( i=0; i<l->n_allocated_data; i++ )
            _freeYstatistics( &(l->data[i]) );
        vtfree( l->data );
    }
    if ( l->pointer != (_Ystatistics**)NULL )
        vtfree( l->pointer );
    _initYstatisticsList( l );
}



static int _allocYstatisticsList( _YstatisticsList *l, int size )
{
    char *proc = "_allocYstatisticsList";
    int i;

    if ( size <= 0 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: negative or zero size\n",  proc );
      return( -1 );
    }

    l->data = (_Ystatistics *)vtmalloc( size * sizeof(_Ystatistics), "l->data", proc );
    if ( l->data == (_Ystatistics *)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate data\n",  proc );
      return( -1 );
    }

    l->pointer = (_Ystatistics **)vtmalloc( size * sizeof(_Ystatistics *), "l->pointer", proc );
    if ( l->pointer == (_Ystatistics **)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate pointer\n",  proc );
      vtfree( l->data );
      return( -1 );
    }

    for ( i=0; i<size; i++ ) {
      _initYstatistics( &(l->data[i]) );
      l->pointer[i] = &(l->data[i]);
    }

    l->n_allocated_data = size;

    return( 1 );
}



typedef struct _computeYmeanParam {
    vt_image *image;
    _YstatisticsList *list;
} _computeYmeanParam;



#define _COMPUTEYMEAN( TYPE ) {                           \
    TYPE ***theArray = (TYPE***)image->array;             \
    for ( i=first; z<=zlast; z++, x=0 ) {                 \
        xend = (z==zlast) ? xlast+1 : image->dim.x;       \
        for ( ; x<xend; x++, i++ ) {                      \
            s = &(list->data[i]);                         \
            s->x = x;                                     \
            s->z = z;                                     \
            for ( sum=0.0, y=0; y<image->dim.y; y++ )     \
                sum += theArray[z][y][x];                 \
            s->intensityMean = sum / (float)image->dim.y; \
        }                                                 \
    }                                                     \
}



static void *_computeYmeanProcedure( void *par )
{
    char *proc = "_computeYmeanProcedure";

    typeChunk *chunk = (typeChunk *)par;
    void *parameter = chunk->parameters;
    size_t first = chunk->first;
    size_t last = chunk->last;

    _computeYmeanParam *p = (_computeYmeanParam*)parameter;
    vt_image *image = p->image;
    _YstatisticsList *list = p->list;

    _Ystatistics *s;
    float sum;
    size_t i, x, y, z;
    size_t xfirst, zfirst;
    size_t xlast, zlast;
    size_t xend;

    z = zfirst = first / image->dim.x;
    x = xfirst = first - zfirst*image->dim.x;

    zlast = last / image->dim.x;
    xlast = last - zlast*image->dim.x;


    switch( image->type ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such image type not handled\n", proc );
      chunk->ret = 0;
      return( (void*)NULL );
      break;
    case SSHORT :
      _COMPUTEYMEAN( s16 )
      break;
    case USHORT :
      _COMPUTEYMEAN( u16 )
      break;
    }

    chunk->ret = 1;
    return( (void*)NULL );
}



static int _computeYmean( vt_image *image, _YstatisticsList *list )
{
    char *proc = "_computeYmean";
    typeChunks chunks;
    _computeYmeanParam parameters;
    int n;

    initChunks( &chunks );
    if ( buildChunks( &chunks, 0, image->dim.x*image->dim.z - 1, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when building chunks\n", proc );
      return( -1 );
    }

    parameters.image = image;
    parameters.list = list;
    for ( n=0; n<chunks.n_allocated_chunks; n++ )
      chunks.data[n].parameters = (void*)(&parameters);

    if ( processChunks( &_computeYmeanProcedure, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute Y mean measures\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    freeChunks( &chunks );

    return( 1 );
}



typedef struct _yvalue {
    int y;
    float value;
    float residual;
} _yvalue;



static int _compareYvalueResidual( const void *x1, const void *x2 )
{
    _yvalue *v1 = (_yvalue*) x1;
    _yvalue *v2 = (_yvalue*) x2;
    /* to sort in increasing order
     */
    if ( v1->residual < v2->residual )
        return( -1 );
    else if ( v1->residual > v2->residual )
        return( 1 );
    return( 0 );
}



typedef struct _computeYrobustMeanParam {
    vt_image *image;
    _YstatisticsList *list;
    float yRejectedFraction;
} _computeYrobustMeanParam;




#define _COMPUTEYROBUSTMEAN( TYPE ) {                                    \
    TYPE ***theArray = (TYPE***)image->array;                            \
    for ( i=first; i<=last; i++ ) {                                      \
        s = list->pointer[i];                                            \
        for ( y=0; y<dimy; y++ ) {                                       \
            line[y].value = theArray[s->z][y][s->x];                     \
            line[y].y = y;                                               \
            line[y].residual = fabs( line[y].value - s->intensityMean ); \
        }                                                                \
        qsort( line, dimy, sizeof(_yvalue), &_compareYvalueResidual );   \
        for ( n=0; n<niterations; n++ ) {                                \
            for ( sum=0.0, y=0; y<ltsy; y++ ) sum += line[y].value;      \
            sum /= (float)ltsy;                                          \
            for ( y=0; y<dimy; y++ ) line[y].residual = fabs( line[y].value - sum ); \
            qsort( line, dimy, sizeof(_yvalue), &_compareYvalueResidual ); \
        }                                                                \
        s->intensityRobustMean = sum;                                    \
        s->excludedPoint = (int*)vtmalloc( (dimy-ltsy)*sizeof(int), "s->excludedPoint", proc ); \
        if ( s->excludedPoint == (int*)NULL ) {                          \
            vtfree( line );                                              \
            if ( _verbose_ )                                             \
                fprintf( stderr, "%s: allocation error\n", proc );       \
            chunk->ret = 0;                                              \
            return( (void*)NULL );                                       \
        }                                                                \
        for ( j=0, y=ltsy; y<dimy; y++, j++ ) s->excludedPoint[j] = line[y].y; \
        s->nExcludedPoint = dimy-ltsy;                                   \
    }                                                                    \
}




static void *_computeYrobustMeanProcedure( void *par )
{
    char *proc = "_computeYrobustMeanProcedure";

    typeChunk *chunk = (typeChunk *)par;
    void *parameter = chunk->parameters;
    size_t first = chunk->first;
    size_t last = chunk->last;

    _computeYrobustMeanParam *p = (_computeYrobustMeanParam*)parameter;
    vt_image *image = p->image;
    _YstatisticsList *list = p->list;

    _Ystatistics *s;
    float sum = 0.0;
    size_t i, y;
    _yvalue *line;
    int j, n, niterations=0;
    size_t ltsy, dimy;

    dimy = image->dim.y;
    ltsy = (image->dim.y * (1.0 - p->yRejectedFraction) + 0.5);

    line = (_yvalue*)vtmalloc( image->dim.y*sizeof(_yvalue), "line", proc );
    if ( line == (_yvalue*)NULL ) {
        if ( _verbose_ )
            fprintf( stderr, "%s: allocation error\n", proc );
        chunk->ret = 0;
        return( (void*)NULL );
    }

    switch( image->type ) {
    default :
      vtfree( line );
      if ( _verbose_ )
        fprintf( stderr, "%s: such image type not handled\n", proc );
      chunk->ret = 0;
      return( (void*)NULL );
      break;
    case SSHORT :
      _COMPUTEYROBUSTMEAN( s16 )
      break;
    case USHORT :
      _COMPUTEYROBUSTMEAN( u16 )
      break;
    }

    vtfree( line );
    chunk->ret = 1;
    return( (void*)NULL );
}



static int _computeYrobustMean( vt_image *image, _YstatisticsList *list, float yRejectedFraction )
{
    char *proc = "_computeYrobustMean";
    typeChunks chunks;
    _computeYrobustMeanParam parameters;
    int n;

    initChunks( &chunks );
    if ( buildChunks( &chunks, 0, list->n_selected_data - 1, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when building chunks\n", proc );
      return( -1 );
    }

    parameters.image = image;
    parameters.list = list;
    parameters.yRejectedFraction = yRejectedFraction;

    for ( n=0; n<chunks.n_allocated_chunks; n++ )
      chunks.data[n].parameters = (void*)(&parameters);

    if ( processChunks( &_computeYrobustMeanProcedure, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute Y mean measures\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    freeChunks( &chunks );

    return( 1 );
}





/************************************************************
 *
 * slice selection based on statistics along Y line
 *
 ************************************************************/





static int _selectXZslicesGlobally( _YstatisticsList *clist, _XZstatisticsList *slist,
                           int threshold )
{
    int i, n, nslices;
    _Ystatistics *s;

    for ( i=0; i<slist->n_allocated_data; i++ )
            slist->data[i].flagged = 0;

    for ( i=0; i<clist->n_selected_data; i++ ) {
        s = clist->pointer[i];
        for ( n=0; n<s->nExcludedPoint; n++ )
            slist->data[ s->excludedPoint[n] ].flagged ++;
    }

    for ( nslices=0, i=0; i<slist->n_allocated_data; i++ ) {
        if ( slist->data[i].flagged > threshold )
            nslices ++;
        else
            slist->data[i].flagged = 0;
    }

    return( nslices );
}





/************************************************************
 *
 *
 *
 ************************************************************/





static int _imageSelectedYLine( vt_image *theIm, vt_image *resIm,
                                _YstatisticsList *list )
{
    char *proc = "_imageSelectedYLine";
    size_t x, z;
    int i;
    _Ystatistics *s;
    u8 ***theArray;

    VT_Image( resIm );
    VT_InitImage( resIm, (char*)NULL, theIm->dim.x, theIm->dim.z, 1, UCHAR );
    if ( VT_AllocImage( resIm ) != 1 ) {
        if ( _verbose_ )
            fprintf( stderr, "%s: allocation error\n", proc );
        return( -1 );
    }
    theArray = (u8***)resIm->array;

    for ( z=0; z<theIm->dim.z; z++ )
    for ( x=0; x<theIm->dim.x; x++) {
        theArray[0][z][x] = 0;
    }

    for ( i=0; i<list->n_selected_data; i++ ) {
        s = list->pointer[i];
        theArray[0][s->z][s->x] = 255;
    }

    return( 1 );
}



static int _imageOutliersYLine( vt_image *theIm, vt_image *resIm,
                                _YstatisticsList *list )
{
    char *proc = "_imageOutliersYLine";
    size_t x, y, z;
    int i, j;
    _Ystatistics *s;
    u8 ***theArray;

    VT_Image( resIm );
    VT_InitImage( resIm, (char*)NULL, theIm->dim.x, theIm->dim.y, theIm->dim.z, UCHAR );
    if ( VT_AllocImage( resIm ) != 1 ) {
        if ( _verbose_ )
            fprintf( stderr, "%s: allocation error\n", proc );
        return( -1 );
    }
    theArray = (u8***)resIm->array;

    for ( z=0; z<theIm->dim.z; z++ )
    for ( y=0; y<theIm->dim.y; y++ )
    for ( x=0; x<theIm->dim.x; x++) {
        theArray[z][y][x] = 0;
    }

    for ( i=0; i<list->n_selected_data; i++ ) {
        s = list->pointer[i];
        for ( j=0; j<s->nExcludedPoint; j++ )
            theArray[s->z][ s->excludedPoint[j] ][s->x] = 255;
    }

    return( 1 );
}





/************************************************************
 *
 * calcul des corrections lineaires a appliquer
 *
 ************************************************************/





static void _computeXZcorrections( _XZstatisticsList *list )
{
    int i;
    _XZstatistics *s, *p, *n;

    /* forward propagation
     * 1. compute the linear correction so that XZ-slice #i
     *    resembles to XZ-slice #(i-1)
     * 2. propagate correction to the first non-flagged XZ-slice
     */
    for ( i=1; i<list->n_allocated_data; i++ ) {
       if ( list->data[i].flagged == 0 ) continue;
       s = &(list->data[i]);
       p = &(list->data[i-1]);
       s->aToPrevious = p->intensityStddev / s->intensityStddev;
       s->bToPrevious = p->intensityMean - s->aToPrevious * s->intensityMean;
       s->dToPrevious = p->dToPrevious + 1;
       if ( p->flagged == 0 ) continue;
       s->bToPrevious += s->aToPrevious * p->bToPrevious;
       s->aToPrevious *= p->aToPrevious;
    }

    /* backward propagation
     * 1. compute the linear correction so that XZ-slice #i
     *    resembles to XZ-slice #(i-1)
     * 2. propagate correction to the first non-flagged XZ-slice
     */
    for ( i=list->n_allocated_data-1; i>0; i-- ) {
       if ( list->data[i].flagged == 0 ) continue;
       s = &(list->data[i]);
       n = &(list->data[i+1]);
       s->aToNext = n->intensityStddev / s->intensityStddev;
       s->bToNext = n->intensityMean - s->aToNext * s->intensityMean;
       s->dToNext = n->dToNext + 1;
       if ( n->flagged == 0 ) continue;
       s->bToNext += s->aToNext * n->bToNext;
       s->aToNext *= n->aToNext;
    }

    /* average corrections
     */
    if ( 0 ) {
        for ( i=0; i<list->n_allocated_data; i++ ) {
            if ( list->data[i].flagged == 0 ) continue;
            s = &(list->data[i]);
            s->a = (s->aToPrevious + s->aToNext) / 2.0;
            s->b = (s->bToPrevious + s->bToNext) / 2.0;
        }
    }

    /* weighted corrections
     */
    if ( 1 ) {
        for ( i=0; i<list->n_allocated_data; i++ ) {
            if ( list->data[i].flagged == 0 ) continue;
            s = &(list->data[i]);
            s->a = (s->aToPrevious * s->dToNext + s->aToNext * s->dToPrevious ) / (float)(s->dToPrevious + s->dToNext);
            s->b = (s->bToPrevious * s->dToNext + s->bToNext * s->dToPrevious ) / (float)(s->dToPrevious + s->dToNext);
        }
    }

}





/************************************************************
 *
 * application des corrections lineaires
 *
 ************************************************************/





typedef struct _correctionXZsliceParam {
    vt_image *theIm;
    vt_image *resIm;
    typeCorrectionList *list;
} _correctionXZsliceParam;



#define _CORRECTIONXZSLICE( TYPE, MIN, MAX ) { \
    TYPE ***theBuf = (TYPE***)theIm->array;    \
    TYPE ***resBuf = (TYPE***)resIm->array;    \
    for ( j=first; j<=last; j++ ) {            \
        c = &(list->data[j]);                  \
        y = c->y;                              \
        if ( y < 0 || y >= (int)theIm->dim.y ) continue; \
        for ( z=0; z<theIm->dim.z; z++)        \
        for ( x=0; x<theIm->dim.x; x++) {      \
          v = c ->a * (float)theBuf[z][y][x] + c->b; \
          i = ( v > 0.0 ) ? (int)(v+0.5) : (int)(v-0.5); \
          if ( i < MIN )                       \
              resBuf[z][y][x] = MIN;           \
          else if ( i > MAX )                  \
              resBuf[z][y][x] = MAX;           \
          else                                 \
              resBuf[z][y][x] = i;             \
        }                                      \
    }                                          \
}



static void *_correctionXZsliceProcedure( void *par )
{
    char *proc = "_correctionXZsliceProcedure";

    typeChunk *chunk = (typeChunk *)par;
    void *parameter = chunk->parameters;
    size_t first = chunk->first;
    size_t last = chunk->last;

    _correctionXZsliceParam *p = (_correctionXZsliceParam*)parameter;
    vt_image *theIm = p->theIm;
    vt_image *resIm = p->resIm;
    typeCorrectionList *list = p->list;
    typeCorrection *c;
    size_t x, z, j;
    float v;
    int y, i;

    switch( theIm->type ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such image type not handled\n", proc );
      chunk->ret = 0;
      return( (void*)NULL );
      break;
    case SSHORT :
      _CORRECTIONXZSLICE( s16, -32768, 32767 )
      break;
    case USHORT :
      _CORRECTIONXZSLICE( u16, 0, 65535 )
      break;
    }

    chunk->ret = 1;
    return( (void*)NULL );
}



static int _correctionXZslice( vt_image *theIm, vt_image *resIm, typeCorrectionList *list )
{
    char *proc = "_correctionXZslice";
    typeChunks chunks;
    _correctionXZsliceParam parameters;
    int n;


    initChunks( &chunks );
    if ( buildChunks( &chunks, 0, list->n_data-1, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when building chunks\n", proc );
      return( -1 );
    }

    parameters.theIm = theIm;
    parameters.resIm = resIm;
    parameters.list = list;
    for ( n=0; n<chunks.n_allocated_chunks; n++ )
      chunks.data[n].parameters = (void*)(&parameters);

    if ( processChunks( &_correctionXZsliceProcedure, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to correct XZ slices\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    freeChunks( &chunks );

    return( 1 );

}





/************************************************************
 *
 *
 *
 ************************************************************/





static int _removeLinesLocally( vt_image *theIm,
                                typeCorrectionList *l,
                                typeRemoveLineParameter *p )
{
  char *proc = "_removeLinesLocally";

  _XZcontrastList contrastList;
  _XZstatisticsList statisticsList;
  int nContrastThreshold;
  int i, nslices;


  if ( _debug_ )
    fprintf( stderr, "... entering %s\n", proc );


  _initXZcontrastList( &contrastList );
  _initXZstatisticsList( &statisticsList );


  /* calcul des contrastes
   */
  if ( _allocXZcontrastList( &contrastList, theIm->dim.y ) != 1 ) {
      if ( _verbose_ )
          fprintf( stderr, "%s: contrast allocation error\n", proc );
      return( -1 );
  }

  if ( _computeXZcontrast( theIm, &contrastList ) != 1 ) {
      _freeXZcontrastList( &contrastList );
      if ( _verbose_ )
          fprintf( stderr, "%s: contrast computation error\n", proc );
      return( -1 );
  }


  /* selection des sections XZ a corriger
   */
  if ( _allocXZstatisticsList( &statisticsList, theIm->dim.y ) != 1 ) {
      _freeXZcontrastList( &contrastList );
      if ( _verbose_ )
          fprintf( stderr, "%s: statistics allocation error\n", proc );
      return( -1 );
  }

  if ( p->contrastSignificantFraction <= 0.5 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: warning, weird (perhaps too small) value for contrast significant fraction (%f)\n",
               proc, p->contrastSignificantFraction );
  }

  /* selection
   * contrastList->data[i].nSupToPrevious > nContrastThreshold
   * -> la coupe #i est plus claire que la coupe #(i-1)
   *    contrastPos1 denombre les augmentations d'intensite
   * contrastList->data[i].nSupToNext > nContrastThreshold
   * -> la coupe #i est plus claire que la coupe #(i+1)
   */
  nContrastThreshold = (double)(theIm->dim.x * theIm->dim.z) * p->contrastSignificantFraction;

  nslices = _selectXZslicesLocally( &contrastList, &statisticsList, nContrastThreshold );


  /* print flagged slices
   */
  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "\n" );
    fprintf( stderr, "- %4d flagged XZ slices where changes occur\n", nslices );
    for ( i=0; i<contrastList.n_allocated_data; i++ ) {
      if ( statisticsList.data[i].flagged == 0 ) continue;
      fprintf( stderr, "#xz=%4d  -  #[i(y-1)<i(y)] = %7d - #[i(y-1)<i(y)] = %7d\n",
               i, contrastList.data[i].nSupToPrevious, contrastList.data[i].nSupToNext );
    }
    fprintf( stderr, "\n" );
  }

  _freeXZcontrastList( &contrastList );


  /* calcul des sections XZ a corriger
   */

  if ( _computeXZstatistics( theIm, &statisticsList ) != 1 ) {
      _freeXZstatisticsList( &statisticsList );
      if ( _verbose_ )
          fprintf( stderr, "%s: statistics computation error\n", proc );
      return( -1 );
  }


  /* calcul des corrections lineaires avec propagation
   */
  _computeXZcorrections( &statisticsList );

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "\n" );
    fprintf( stderr, "- %4d flagged XZ slices where changes occur\n", nslices );
    for ( i=0; i<statisticsList.n_allocated_data; i++ ) {
      if ( statisticsList.data[i].flagged == 0 ) continue;
      fprintf( stderr, "#xz=%4d  -  a=%6.2f b=%6.2f - a[-]=%6.2f b[-]=%6.2f - a[+]=%6.2f b[+]=%6.2f\n",
               i, statisticsList.data[i].a, statisticsList.data[i].b,
               statisticsList.data[i].aToPrevious, statisticsList.data[i].bToPrevious,
               statisticsList.data[i].aToNext, statisticsList.data[i].bToNext );
    }
    fprintf( stderr, "\n" );
  }

  if ( _copyXZstatisticsListToCorrectionList( &statisticsList, l ) <= 0 ) {
      _freeXZstatisticsList( &statisticsList );
      if ( _verbose_ )
          fprintf( stderr, "%s: error when copying list\n", proc );
      return( -1 );
  }

  _freeXZstatisticsList( &statisticsList );

  if ( _debug_ )
    fprintf( stderr, "... exiting %s\n", proc );

  return( 1 );
}







/************************************************************
 *
 *
 *
 ************************************************************/





static int _removeLinesRegionally( vt_image *theIm,
                                   typeCorrectionList *l,
                                   typeRemoveLineParameter *p )
{
  char *proc = "_removeLinesRegionally";

  _XZcontrastList contrastList;
  _XZstatisticsList statisticsList;
  int nContrastThreshold;
  int minThreshold, maxThreshold;
  float minCoefficient = 0.6;
  float maxCoefficient = 0.9;
  int nThreshold = 13;
  int i, nslices;


  if ( _debug_ )
    fprintf( stderr, "... entering %s\n", proc );


  _initXZcontrastList( &contrastList );
  _initXZstatisticsList( &statisticsList );


  /* calcul des contrastes
   */
  if ( _allocXZcontrastList( &contrastList, theIm->dim.y ) != 1 ) {
      if ( _verbose_ )
          fprintf( stderr, "%s: contrast allocation error\n", proc );
      return( -1 );
  }

  if ( _computeXZcontrast( theIm, &contrastList ) != 1 ) {
      _freeXZcontrastList( &contrastList );
      if ( _verbose_ )
          fprintf( stderr, "%s: contrast computation error\n", proc );
      return( -1 );
  }


  /* selection des sections XZ a corriger
   */
  if ( _allocXZstatisticsList( &statisticsList, theIm->dim.y ) != 1 ) {
      _freeXZcontrastList( &contrastList );
      if ( _verbose_ )
          fprintf( stderr, "%s: statistics allocation error\n", proc );
      return( -1 );
  }

  if ( p->contrastSignificantFraction <= 0.5 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: warning, weird (perhaps too small) value for contrast significant fraction (%f)\n",
               proc, p->contrastSignificantFraction );
  }

  /* selection
   * contrastList->data[i].nSupToPrevious > nContrastThreshold
   * -> la coupe #i est plus claire que la coupe #(i-1)
   *    contrastPos1 denombre les augmentations d'intensite
   * contrastList->data[i].nSupToNext > nContrastThreshold
   * -> la coupe #i est plus claire que la coupe #(i+1)
   */
  nContrastThreshold = ( (double)(theIm->dim.x * theIm->dim.z) * p->contrastSignificantFraction + 0.5 );

  if ( _verbose_ >= 3 ) {
      _fprintfXZcontrastEvolution( stderr, &contrastList, nContrastThreshold, 25 );
  }

  /* selection automatique ?
   */
  if ( p->contrastSignificantFraction > 0.0 && p->automatedChoices == 0 ) {
      nslices = _selectXZsliceRegionally( &contrastList, &statisticsList, nContrastThreshold, -1, -1 );
  }
  else {
      minThreshold = ( (double)(theIm->dim.x * theIm->dim.z) * minCoefficient + 0.5 );
      maxThreshold = ( (double)(theIm->dim.x * theIm->dim.z) * maxCoefficient + 0.5 );
      nContrastThreshold = _selectThresholdRegionally( &contrastList, &statisticsList, minThreshold, maxThreshold, nThreshold );
      if ( nContrastThreshold == 0 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to compute optimal coefficient\n", proc );
          return( -1 );
      }
      nslices = _selectXZsliceRegionally( &contrastList, &statisticsList, nContrastThreshold, -1, -1 );
  }

  /* print flagged slices
   */
  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "\n" );
    fprintf( stderr, "- %4d flagged XZ slices where changes occur\n", nslices );
    for ( i=0; i<contrastList.n_allocated_data; i++ ) {
      if ( statisticsList.data[i].flagged == 0 ) continue;
      fprintf( stderr, "#xz=%4d  -  #[i(y-1)<i(y)] = %7d - #[i(y-1)<i(y)] = %7d\n",
               i, contrastList.data[i].nSupToPrevious, contrastList.data[i].nSupToNext );
    }
    fprintf( stderr, "\n" );
  }

  _freeXZcontrastList( &contrastList );


  /* calcul des sections XZ a corriger
   */

  if ( _computeXZstatistics( theIm, &statisticsList ) != 1 ) {
      _freeXZstatisticsList( &statisticsList );
      if ( _verbose_ )
          fprintf( stderr, "%s: statistics computation error\n", proc );
      return( -1 );
  }


  /* calcul des corrections lineaires avec propagation
   */
  _computeXZcorrections( &statisticsList );

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "\n" );
    fprintf( stderr, "- %4d flagged XZ slices where changes occur\n", nslices );
    for ( i=0; i<statisticsList.n_allocated_data; i++ ) {
      if ( statisticsList.data[i].flagged == 0 ) continue;
      fprintf( stderr, "#xz=%4d  -  a=%6.2f b=%6.2f - a[-]=%6.2f b[-]=%6.2f - a[+]=%6.2f b[+]=%6.2f\n",
               i, statisticsList.data[i].a, statisticsList.data[i].b,
               statisticsList.data[i].aToPrevious, statisticsList.data[i].bToPrevious,
               statisticsList.data[i].aToNext, statisticsList.data[i].bToNext );
    }
    fprintf( stderr, "\n" );
  }

  if ( _copyXZstatisticsListToCorrectionList( &statisticsList, l ) <= 0 ) {
      _freeXZstatisticsList( &statisticsList );
      if ( _verbose_ )
          fprintf( stderr, "%s: error when copying list\n", proc );
      return( -1 );
  }

  _freeXZstatisticsList( &statisticsList );

  if ( _debug_ )
    fprintf( stderr, "... exiting %s\n", proc );

  return( 1 );
}





/************************************************************
 *
 *
 *
 ************************************************************/





int _removeLinesGlobally( vt_image *theIm,
                          typeCorrectionList *l,
                          typeRemoveLineParameter *p )
{
  char *proc = "_removeLinesGlobally";
  _YstatisticsList yStatisticsList;
  _XZstatisticsList xzStatisticsList;

  vt_image tmpIm;
  int nContrastThreshold;
  int nslices;
  int i;

  if ( _debug_ )
    fprintf( stderr, "... entering %s\n", proc );

  _initYstatisticsList( &yStatisticsList );
  if ( _allocYstatisticsList( &yStatisticsList, theIm->dim.x * theIm->dim.z ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }

  /* calcul des moyennes pour les lignes selon Y
   */
  if ( _computeYmean( theIm, &yStatisticsList ) != 1 ) {
      _freeYstatisticsList( &yStatisticsList );
      if ( _verbose_ )
        fprintf( stderr, "%s: mean computation error\n", proc );
      return( -1 );
  }

  /* sort wrt mean in ascending order
   * and Y line selection of smallest mean
   * => select Y lines that are the more likely in the bakground
   * (ie without embryo signal)
   */
  qsort( yStatisticsList.pointer, yStatisticsList.n_allocated_data,
         sizeof(_Ystatistics*), &_compareYstatisticsMean );

  yStatisticsList.n_selected_data = (theIm->dim.x * theIm->dim.z * p->xzKeptFraction + 0.5);

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "\n" );
    fprintf( stderr, "- %6d/%8lu Y lines selected\n", yStatisticsList.n_selected_data, theIm->dim.x * theIm->dim.z );
    if ( 0 ) {
      for ( i=0; i<10; i++ ) {
        fprintf( stderr, "#%6d: x=%4d z=%4d mean=%f robust mean=%f\n", i,
                 yStatisticsList.pointer[i]->x, yStatisticsList.pointer[i]->z,
                 yStatisticsList.pointer[i]->intensityMean,  yStatisticsList.pointer[i]->intensityRobustMean );
      }
      fprintf( stderr, "...\n" );
      for ( i=yStatisticsList.n_allocated_data-10; i<yStatisticsList.n_allocated_data; i++ ) {
        fprintf( stderr, "#%6d: x=%4d z=%4d mean=%f robust mean=%f\n", i,
                 yStatisticsList.pointer[i]->x, yStatisticsList.pointer[i]->z,
                 yStatisticsList.pointer[i]->intensityMean,  yStatisticsList.pointer[i]->intensityRobustMean );
      }
    }
  }

  if ( _debug_ ) {
      if ( _imageSelectedYLine( theIm, &tmpIm, &yStatisticsList ) != 1 ) {
          _freeYstatisticsList( &yStatisticsList );
          if ( _verbose_ )
            fprintf( stderr, "%s: image computation error\n", proc );
          return( -1 );
      }
      if ( VT_WriteInrimageWithName( &tmpIm, "selectedYlines.mha" ) != 1 ) {
          _freeYstatisticsList( &yStatisticsList );
          if ( _verbose_ )
            fprintf( stderr, "%s: image writing error\n", proc );
          return( -1 );
      }
      VT_FreeImage( &tmpIm );
  }


  /* calcul des moyennes robustes et des outliers
   * pour les lignes selectionnees
   *
   * chaque _Ystatistics contient
   * - intensityMean : la moyenne
   * - intensityRobustMean : la moyenne robuste
   * - nExcludedPoint : le nombre de points exclus pour le calcul de la moyenne robuste
   * - *excludedPoint : la liste des points exclus (le y) pour la moyenne robuste
   *
   * => permet de marquer les coordonnees Y qui ont ete exclues du calcul de la moyenne
   *    robuste
   */
  if ( _computeYrobustMean( theIm, &yStatisticsList, p->yRejectedFraction ) != 1 ) {
      _freeYstatisticsList( &yStatisticsList );
      if ( _verbose_ )
        fprintf( stderr, "%s: robuste mean computation error\n", proc );
      return( -1 );
  }


  /* selection des sections XZ a corriger
   */
  nContrastThreshold = ( (double)(yStatisticsList.n_selected_data) * p->contrastSignificantFraction + 0.5 );

  _initXZstatisticsList( &xzStatisticsList );

  if ( _allocXZstatisticsList( &xzStatisticsList, theIm->dim.y ) != 1 ) {
      _freeYstatisticsList( &yStatisticsList );
      if ( _verbose_ )
          fprintf( stderr, "%s: statistics allocation error\n", proc );
      return( -1 );
  }

  /* on compte le nombre de fois qu'une section XY (une valeur de Y) a ete exclue
   * du calcul de la moyenne robuste.
   * on a yStatisticsList.n_selected_data echantillons,
   * on estime que s'il y en a plus que (yStatisticsList.n_selected_data * p->contrastSignificantFraction)
   * alors la section doit etre corrigee
   */
  nslices = _selectXZslicesGlobally( &yStatisticsList, &xzStatisticsList, nContrastThreshold );
  if ( nslices == 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: no slices to be corrected were found\n", proc );
    return( 0 );
  }

  if ( _debug_ ) {
      if ( _imageOutliersYLine( theIm, &tmpIm, &yStatisticsList ) != 1 ) {
          _freeYstatisticsList( &yStatisticsList );
          if ( _verbose_ )
            fprintf( stderr, "%s: image computation error\n", proc );
          return( -1 );
      }
      if ( VT_WriteInrimageWithName( &tmpIm, "outliersYlines.mha" ) != 1 ) {
          _freeYstatisticsList( &yStatisticsList );
          if ( _verbose_ )
            fprintf( stderr, "%s: image writing error\n", proc );
          return( -1 );
      }
      VT_FreeImage( &tmpIm );
  }


  /* print flagged slices
   */
  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "\n" );
    fprintf( stderr, "- %4d flagged XZ slices where changes occur\n", nslices );
    for ( i=0; i<xzStatisticsList.n_allocated_data; i++ ) {
      if ( xzStatisticsList.data[i].flagged == 0 ) continue;
      fprintf( stderr, "#xz=%4d  -  #[outliers] = %7d\n",
               i, xzStatisticsList.data[i].flagged );
    }
    fprintf( stderr, "\n" );
  }

  _freeYstatisticsList( &yStatisticsList );


  /* calcul des sections XZ a corriger
   * pour chaque coupe, on calcule sa moyenne et son ecart-type
   */

  if ( _computeXZstatistics( theIm, &xzStatisticsList ) != 1 ) {
      _freeXZstatisticsList( &xzStatisticsList );
      if ( _verbose_ )
          fprintf( stderr, "%s: statistics computation error\n", proc );
      return( -1 );
  }


  /* calcul des corrections lineaires avec propagation
   * 1. une correction lineaire est calculee avec la coupe precedente
   *    et avec propagation
   * 2. une correction lineaire est calculee avec la coupe suivante
   *    et avec propagation
   * 3. la correction finale est une combinaison lineaire des 2 precedentes
   */
  _computeXZcorrections( &xzStatisticsList );

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "\n" );
    fprintf( stderr, "- %4d flagged XZ slices where changes occur\n", nslices );
    for ( i=0; i<xzStatisticsList.n_allocated_data; i++ ) {
      if ( xzStatisticsList.data[i].flagged == 0 ) continue;
      fprintf( stderr, "#xz=%4d  -  a=%6.2f b=%6.2f - a[-]=%6.2f b[-]=%6.2f - a[+]=%6.2f b[+]=%6.2f\n",
               i, xzStatisticsList.data[i].a, xzStatisticsList.data[i].b,
               xzStatisticsList.data[i].aToPrevious, xzStatisticsList.data[i].bToPrevious,
               xzStatisticsList.data[i].aToNext, xzStatisticsList.data[i].bToNext );
    }
    fprintf( stderr, "\n" );
  }

  if ( _copyXZstatisticsListToCorrectionList( &xzStatisticsList, l ) <= 0 ) {
      _freeXZstatisticsList( &xzStatisticsList );
      if ( _verbose_ )
          fprintf( stderr, "%s: error when copying list\n", proc );
      return( -1 );
  }

  _freeXZstatisticsList( &xzStatisticsList );

  if ( _debug_ )
    fprintf( stderr, "... exiting %s\n", proc );

  return( 1 );
}






/************************************************************
 *
 *
 *
 ************************************************************/





int VT_CorrectLines( vt_image *theIm, vt_image *resIm,
                     typeCorrectionList *l )
{
    char *proc = "VT_CorrectLines";

    if ( VT_CopyImage( theIm, resIm ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to copy image\n", proc );
      return( -1 );
    }

    if ( l->n_data > 0 ) {
      if ( _correctionXZslice( theIm, resIm, l ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to correct image\n", proc );
          return( -1 );
      }
    }

    return( 1 );
}




int VT_RemoveLines( vt_image *theIm, vt_image *resIm,
                    typeCorrectionList *l,
                    typeRemoveLineParameter *p )
{
  char *proc = "VT_RemoveLines";
  int ret;

  if ( _verbose_ >= 3 ) {
    fprintf( stderr, "----- parameters ----\n" );
    fprintfRemoveLineParameter( stderr, p );
    fprintf( stderr, "---------------------\n" );
  }

  switch( p->method ) {
  default :
      if ( _verbose_ )
          fprintf( stderr, "%s: such correction method not handled yet\n", proc );
      return( -1 );

  case _LOCAL_ :
      if ( _removeLinesLocally( theIm, l, p ) != 1 ) {
          if ( _verbose_ )
              fprintf( stderr, "%s: error when computing (local method)\n", proc );
          return( -1 );
      }
      break;

  case _REGIONAL_ :
      if ( _removeLinesRegionally( theIm, l, p ) != 1 ) {
          if ( _verbose_ )
              fprintf( stderr, "%s: error when computing (regional method)\n", proc );
          return( -1 );
      }
      break;

  case _GLOBAL_ :
    ret = _removeLinesGlobally( theIm, l, p );
    /* lines to be corrected were found
     */
    if ( ret == 1 ) {
      ;
    }
    /* no corrections to be done
     */
    else if ( ret == 0 ) {
      if ( VT_CopyImage( theIm, resIm ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to copy image\n", proc );
        return( -1 );
      }
      if ( _verbose_ )
        fprintf( stderr, "%s: no slices to be corrected were found\n", proc );
      return( 0 );
    }
    else if ( ret < 0 ) {
      if ( _verbose_ )
          fprintf( stderr, "%s: error when computing (global method)\n", proc );
      return( -1 );
    }
    else {
      if ( _verbose_ )
          fprintf( stderr, "%s: unexpected returned code\n", proc );
      return( -1 );
    }
    break;
  }

  if ( VT_CorrectLines( theIm, resIm, l ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to correct image\n", proc );
    return( -1 );
  }

  return( 0 );
}
