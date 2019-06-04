
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <ImageIO.h>
#include "litetiff.h"


static int _verbose_ = 1;
static int _debug_ = 0;


void setVerboseInLiteTiff( int v )
{
  _verbose_ = v;
}

void incrementVerboseLiteTiff( )
{
  _verbose_ ++;
}

void setDebugInLiteTiff( int d )
{
  _debug_ = d;
}

void incrementDebugInLiteTiff( )
{
  _debug_ ++;
}


/************************************************************
 *
 *
 *
 ************************************************************/

static int _consider_palette_as_grey_level_ = 0;

/************************************************************
 *
 *
 *
 ************************************************************/





static void _swap4bytes( unsigned char *pntr)
{
  unsigned char b0, b1, b2, b3;

  b0 = *pntr;
  b1 = *(pntr+1);
  b2 = *(pntr+2);
  b3 = *(pntr+3);

  *pntr = b3;
  *(pntr+1) = b2;
  *(pntr+2) = b1;
  *(pntr+3) = b0;
}



static void _swap2bytes( unsigned char *pntr )
{
    unsigned char b0, b1;

    b0 = *pntr;
    b1 = *(pntr+1);

    *pntr = b1;
    *(pntr+1) = b0;
}



typedef struct _tiffHeader {
   unsigned char data[8];
   unsigned short int *version;
   unsigned int *offset;
   int swap;
} _tiffHeader;



static void _initTiffHeader( _tiffHeader *h )
{
  int i;
  for ( i=0; i<8; i++ ) h->data[i] = 0;
  h->version = (unsigned short int*)&(h->data[2]);
  h->offset = (unsigned int*)&(h->data[4]);
  h->swap = 0;
}


static void _fprintfTiffHeader( FILE *f, _tiffHeader *h )
{
  fprintf( f, "Tiff header = %c%c version = %d offset = %u swap = %d\n",
           h->data[0], h->data[1], *(h->version), *(h->offset), h->swap );
}


/************************************************************
 *
 * _LsmHeader
 *
 ************************************************************/



typedef struct _LsmHeader {
  int MagicNumber;
  int StructureSize;

  int DimensionX;
  int DimensionY;
  int DimensionZ;

  /* number of channels
   */
  int DimensionChannels;

  /* Timestack size
   */
  int DimensionTime;

  int IntensityDataType;

  int ThumbnailX;
  int ThumbnailY;

  double VoxelSizeX;
  double VoxelSizeY;
  double VoxelSizeZ;

  double OriginX;
  double OriginY;
  double OriginZ;

  unsigned short int ScanType;
  unsigned short int SpectralScan;

  int DataType;
  int OffsetVectorOverlay;
  int OffsetInputLut;
  int OffsetOutputLut;
  int OffsetChannelColors;

  double TimeIntervall;

  int OffsetChannelDataTypes;

  int OffsetScanInformation;
  int OffsetKsData;
  int OffsetTimeStamps;
  int OffsetEventList;
  int OffsetRoi;
  int OffsetBleachRoi;
  int OffsetNextRecording;

  double DisplayAspectX;
  double DisplayAspectY;
  double DisplayAspectZ;
  double DisplayAspectTime;

  int OffsetMeanOfRoisOverlay;
  int OffsetTopoIsolineOverlay;
  int OffsetTopoProfileOverlay;
  int OffsetLinescanOverlay;

  int ToolbarFlags;
  int OffsetChannelWavelength;
  int OffsetChannelFactors;
  int ObjectiveSphereCorrection;
  int OffsetUnmixParameters;
} _LsmHeader;



static void _initLsmHeader( _LsmHeader *h )
{
  h->MagicNumber = 0;
  h->StructureSize = 0;

  h->DimensionX = 0;
  h->DimensionY = 0;
  h->DimensionZ = 0;

  /* number of channels
   */
  h->DimensionChannels = 0;

  /* Timestack size
   */
  h->DimensionTime = 0;

  h->IntensityDataType = 0;

  h->ThumbnailX = 0;
  h->ThumbnailY = 0;

  h->VoxelSizeX = 0.0;
  h->VoxelSizeY = 0.0;
  h->VoxelSizeZ = 0.0;

  h->OriginX = 0.0;
  h->OriginY = 0.0;
  h->OriginZ = 0.0;

  h->ScanType = 0;
  h->SpectralScan = 0;

  h->DataType = 0;
  h->OffsetVectorOverlay = 0;
  h->OffsetInputLut = 0;
  h->OffsetOutputLut = 0;
  h->OffsetChannelColors = 0;

  h->TimeIntervall = 0.0;

  h->OffsetChannelDataTypes = 0;

  h->OffsetScanInformation = 0;
  h->OffsetKsData = 0;
  h->OffsetTimeStamps = 0;
  h->OffsetEventList = 0;
  h->OffsetRoi = 0;
  h->OffsetBleachRoi = 0;
  h->OffsetNextRecording = 0;

  h->DisplayAspectX = 0.0;
  h->DisplayAspectY = 0.0;
  h->DisplayAspectZ = 0.0;
  h->DisplayAspectTime = 0.0;

  h->OffsetMeanOfRoisOverlay = 0;
  h->OffsetTopoIsolineOverlay = 0;
  h->OffsetTopoProfileOverlay = 0;
  h->OffsetLinescanOverlay = 0;

  h->ToolbarFlags = 0;
  h->OffsetChannelWavelength = 0;
  h->OffsetChannelFactors = 0;
  h->ObjectiveSphereCorrection = 0;
  h->OffsetUnmixParameters = 0;
}



static void _swapLsmHeader( _LsmHeader *h __attribute__ ((unused)) )
{
  char *proc = "_swapLsmHeader";
  if ( 1 )
    fprintf( stderr, "%s: not implemented yet, reading may be hazardous\n", proc );
}



static void _fprintfLsmHeader( FILE *f, _LsmHeader *h )
{
  fprintf( f, "   MagicNumber               = %d\n", h->MagicNumber );
  fprintf( f, "   StructureSize             = %d\n", h->StructureSize );

  fprintf( f, "   DimensionX                = %d\n", h->DimensionX );
  fprintf( f, "   DimensionY                = %d\n", h->DimensionY );
  fprintf( f, "   DimensionZ                = %d\n", h->DimensionZ );

  fprintf( f, "   DimensionChannels         = %d\n", h->DimensionChannels );

  fprintf( f, "   DimensionTime             = %d\n", h->DimensionTime );

  fprintf( f, "   IntensityDataType         = %d\n", h->IntensityDataType );

  fprintf( f, "   ThumbnailX                = %d\n", h->ThumbnailX );
  fprintf( f, "   ThumbnailY                = %d\n", h->ThumbnailY );

  fprintf( f, "   VoxelSizeX                = %g\n", h->VoxelSizeX );
  fprintf( f, "   VoxelSizeY                = %g\n", h->VoxelSizeY );
  fprintf( f, "   VoxelSizeZ                = %g\n", h->VoxelSizeZ );

  fprintf( f, "   OriginX                   = %g\n", h->OriginX );
  fprintf( f, "   OriginY                   = %g\n", h->OriginY );
  fprintf( f, "   OriginZ                   = %g\n", h->OriginZ );

  fprintf( f, "   ScanType                  = %d\n", h->ScanType );
  fprintf( f, "   SpectralScan              = %d\n", h->SpectralScan );

  fprintf( f, "   DataType                  = %d\n", h->DataType );
  fprintf( f, "   OffsetVectorOverlay       = %d\n", h->OffsetVectorOverlay );
  fprintf( f, "   OffsetInputLut            = %d\n", h->OffsetInputLut );
  fprintf( f, "   OffsetOutputLut           = %d\n", h->OffsetOutputLut );
  fprintf( f, "   OffsetChannelColors       = %d\n", h->OffsetChannelColors );

  fprintf( f, "   TimeIntervall             = %f\n", h->TimeIntervall );

  fprintf( f, "   OffsetChannelDataTypes    = %d\n", h->OffsetChannelDataTypes );

  fprintf( f, "   OffsetScanInformation     = %d\n", h->OffsetScanInformation );
  fprintf( f, "   OffsetKsData              = %d\n", h->OffsetKsData );
  fprintf( f, "   OffsetTimeStamps          = %d\n", h->OffsetTimeStamps );
  fprintf( f, "   OffsetEventList           = %d\n", h->OffsetEventList );
  fprintf( f, "   OffsetRoi                 = %d\n", h->OffsetRoi );
  fprintf( f, "   OffsetBleachRoi           = %d\n", h->OffsetBleachRoi );
  fprintf( f, "   OffsetNextRecording       = %d\n", h->OffsetNextRecording );

  fprintf( f, "   DisplayAspectX            = %f\n", h->DisplayAspectX );
  fprintf( f, "   DisplayAspectY            = %f\n", h->DisplayAspectY );
  fprintf( f, "   DisplayAspectZ            = %f\n", h->DisplayAspectZ );
  fprintf( f, "   DisplayAspectTime         = %f\n", h->DisplayAspectTime );

  fprintf( f, "   OffsetMeanOfRoisOverlay   = %d\n", h->OffsetMeanOfRoisOverlay );
  fprintf( f, "   OffsetTopoIsolineOverlay  = %d\n", h->OffsetTopoIsolineOverlay );
  fprintf( f, "   OffsetTopoProfileOverlay  = %d\n", h->OffsetTopoProfileOverlay );
  fprintf( f, "   OffsetLinescanOverlay     = %d\n", h->OffsetLinescanOverlay );

  fprintf( f, "   ToolbarFlags              = %d\n", h->ToolbarFlags );
  fprintf( f, "   OffsetChannelWavelength   = %d\n", h->OffsetChannelWavelength );
  fprintf( f, "   OffsetChannelFactors      = %d\n", h->OffsetChannelFactors );
  fprintf( f, "   ObjectiveSphereCorrection = %d\n", h->ObjectiveSphereCorrection );
  fprintf( f, "   OffsetUnmixParameters     = %d\n", h->OffsetUnmixParameters );
}

/************************************************************
 *
 * _tiffIfdEntry
 *
 ************************************************************/



typedef struct _tiffIfdEntry {
  unsigned char data[12];
  unsigned short int *tag;
  unsigned short int *fieldType;
  unsigned int *length;
  unsigned int *offset;
  char *values;
  unsigned int valuesLength;
} _tiffIfdEntry;



static void _setPointersTiffIfdEntry( _tiffIfdEntry *e )
{
    e->tag = (unsigned short int *)&(e->data[0]);
    e->fieldType = (unsigned short int *)&(e->data[2]);
    e->length = (unsigned int *)&(e->data[4]);
    e->offset = (unsigned int *)&(e->data[8]);
}



static void _initTiffIfdEntry( _tiffIfdEntry *e )
{
  int i;
  for ( i=0; i<12; i++ ) e->data[i] = 0;
  e->values = (char*)NULL;
  e->valuesLength = 0;
  _setPointersTiffIfdEntry( e );
}



static void _freeTiffIfdEntry( _tiffIfdEntry *e )
{
  if ( e->values != (char*)NULL )
    free( e->values );
  _initTiffIfdEntry( e );
}



static int _fieldTypeLength( _tiffIfdEntry *e )
{
  int l;
  switch( *(e->fieldType) ) {
  default : l = 0; break;
  case  1 : l = sizeof( unsigned char ); break;
  case  2 : l = sizeof( char ); break;
  case  3 : l = sizeof( unsigned short int ); break;
  case  4 : l = sizeof( unsigned int ); break;
  case  5 : l = 2 * sizeof( unsigned int ); break;
  case  6 : l = sizeof( char ); break;
  case  7 : l = sizeof( char ); break;
  case  8 : l = sizeof( short int ); break;
  case  9 : l = sizeof( int ); break;
  case 10 : l = 2 * sizeof( int ); break;
  case 11 : l = sizeof( float ); break;
  case 12 : l = sizeof( double ); break;
  }
  return( l );
}



static int _valueLength( _tiffIfdEntry *e )
{
  return( *(e->length) * _fieldTypeLength( e ) );
}



static void _swapValue( _tiffIfdEntry *e )
{
  char *proc = "_swapValue";
  int vl = _valueLength( e );

  if ( vl > 4 ) {
    /* swap offset */
    _swap4bytes( (unsigned char *)&(e->data[8]) );
  }
  else if ( vl == 4 ) {
    /* swap two 2-bytes values */
    if ( _fieldTypeLength( e ) == 2 ) {
      _swap2bytes( (unsigned char *)&(e->data[8]) );
      _swap2bytes( (unsigned char *)&(e->data[10]) );
    }
    else {
      /* swap one 4-bytes value */
     _swap4bytes( (unsigned char *)&(e->data[8]) );
    }
  }
  else if ( vl == 2 ) {
    _swap2bytes( (unsigned char *)&(e->data[8]) );
  }
  else if ( vl == 1 ) {
    ;
  }
  else {
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: weird case, value length is %d\n", proc, vl );
  }
}



/* 1  = entries are equal
 * 0  = entries are different
 * -1 = comparison can not be made
 */
static int _AreIfdEntriesEqual( _tiffIfdEntry *a, _tiffIfdEntry *b )
{
  char *proc = "_AreIfdEntriesEqual";

  if ( *(a->tag) != *(b->tag) ) {
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: entries have different tags\n", proc );
    return( 0 );
  }

  if ( *(a->fieldType) != *(b->fieldType) ) {
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: entries have different field types\n", proc );
    return( 0 );
  }

  if ( *(a->length) != *(b->length) ) {
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: entries have different lengths\n", proc );
    return( 0 );
  }

  if ( *(a->length) == (unsigned int)1 ) {
    switch( *(a->fieldType) ) {
    default :
      return( -1 );
    case  1 :
    case  2 :
    case  6 :
    case  7 :
      /* 1-byte */
      if ( *(a->offset) != *(b->offset) ) return( 0 );
      return( 1 );
    case  3 :
    case  8 :
      /* 2-bytes */
      if ( *(a->offset) != *(b->offset) ) return( 0 );
      return( 1 );
    case  4 :
    case  9 :
    case 11 :
      /* 4-bytes */
      if ( *(a->offset) != *(b->offset) ) return( 0 );
      return( 1 );
    case 5 :
    case 10 :
    case 12 :
      return( -1 );
    }
  }
  /* length > 1
   */
  else {
    return( -1 );
  }

  return( 1 );
}



/************************************************************
 *
 * _tiffIfd: list of _tiffIfdEntry
 *
 ************************************************************/



typedef struct _tiffIfd {
  _tiffIfdEntry *data;
  unsigned short int n_data;
  unsigned short int n_allocated_data;
  unsigned int offset;
} _tiffIfd;



static void _initTiffIfd( _tiffIfd *i )
{
  i->data = (_tiffIfdEntry*)NULL;
  i->n_data = 0;
  i->n_allocated_data = 0;
  i->offset = 0;
}



static void _freeTiffIfd( _tiffIfd *i )
{
  int j;
  if ( i->data != (_tiffIfdEntry*)NULL ) {
    for ( j=0; j<i->n_data; j++ )
      _freeTiffIfdEntry( &(i->data[j]) );
    free( i->data );
  }
  _initTiffIfd( i );
}



static int _size_to_be_allocated_ = 20;



static int _addTiffIfdEntryToIfd( _tiffIfd *l, _tiffIfdEntry *p )
{
    char *proc = "_addTiffIfdEntryToIfd";
    int i, s = l->n_allocated_data;
    _tiffIfdEntry *data;

    if ( l->n_data == l->n_allocated_data ) {
      s += _size_to_be_allocated_;
      data = (_tiffIfdEntry*)malloc( s * sizeof(_tiffIfdEntry) );
      if ( data == (_tiffIfdEntry*)NULL ) {
        if ( _verbose_ || _debug_ )
            fprintf( stderr, "%s: allocation error\n", proc );
        return( -1 );
      }
      if ( l->n_allocated_data > 0 ) {
        (void)memcpy( data, l->data, l->n_allocated_data*sizeof(_tiffIfdEntry) );
        free( l->data );
      }
      for ( i=l->n_allocated_data; i<s; i++ )
        _initTiffIfdEntry( &(data[i]) );
      l->n_allocated_data = s;
      l->data = data;
    }

    memcpy( &(l->data[l->n_data]), p, sizeof(_tiffIfdEntry) );
    _setPointersTiffIfdEntry( &(l->data[l->n_data]) );
    l->n_data ++;

    return( 1 );
}



static int _TestIfdCoherency( _tiffIfd *ifd0, _tiffIfd *ifd1, int verbose )
{
  char *proc = "_TestIfdCoherency";
  /* { 256, "ImageWidth" },
   * { 257, "ImageLength" },
   * { 262, "PhotometricInterpretation" }
   */
  unsigned short int tagtobetested[] = { 256, 257, 262, 0 };
  int e, e0, e1, tag=0;

  for ( tag=0; tagtobetested[tag] > 0; tag++ ) {

    /* search for tag in ifd #0
     */
    for ( e0=-1, e=0; e<ifd0->n_data; e++ ) {
      if ( *(ifd0->data[e].tag) == tagtobetested[tag] ) {
        e0 = e;
      }
    }

    if ( e0 == -1 ) {
      if ( verbose && _verbose_ )
        fprintf( stderr, "%s: tag #%d was not found in ifd #1\n",
                 proc, tagtobetested[tag] );
      return( -1 );
    }

    for ( e1=-1, e=0; e<ifd1->n_data; e++ ) {
      if ( *(ifd1->data[e].tag) == tagtobetested[tag] ) {
        e1 = e;
      }
    }

    if ( e1 == -1 ) {
      if ( verbose && _verbose_ )
        fprintf( stderr, "%s: tag #%d was not found in ifd #2\n",
                 proc, tagtobetested[tag] );
      return( -1 );
    }

    if ( e0 == -1 && e1 == -1 ) {
      continue;
    }
    else if ( e0 == -1 && e1 >= 0 ) {
      if ( verbose && _verbose_ )
        fprintf( stderr, "%s: tag #%d was not found in ifd #1 and found in ifd #2\n",
                 proc, tagtobetested[tag] );
      return( -1 );
    }
    else if ( e0 >= 0 && e1 == -1 ) {
      if ( verbose && _verbose_ )
        fprintf( stderr, "%s: tag #%d was found in ifd #1 and not found in ifd #2\n",
                 proc, tagtobetested[tag] );
      return( -1 );
    }
    else {
      if ( _AreIfdEntriesEqual( &(ifd0->data[e0]), &(ifd1->data[e1]) ) != 1 ) {
        if ( verbose && _verbose_ )
          fprintf( stderr, "%s: tag #%d values are different for ifd #1 and ifd #2\n",
                   proc, tagtobetested[tag] );
        return( -1 );
      }
    }

  }

  return( 1 );
}



/************************************************************
 *
 * _tiffIfdList: list of _tiffIfd
 * describe the whole image
 *
 ************************************************************/



typedef struct _tiffIfdList {
  _tiffIfd *data;
  int n_data;
  int n_allocated_data;
} _tiffIfdList;



static void _initTiffIfdList( _tiffIfdList *l )
{
  l->data = (_tiffIfd*)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;
}



static void _freeTiffIfdList( _tiffIfdList *l )
{
  int i;

  if ( l->data != (_tiffIfd*)NULL ) {
    for ( i=0;i<l->n_data; i++ )
      _freeTiffIfd( &(l->data[i]) );
    free ( l->data );
  }
  _initTiffIfdList( l );
}



static int _addTiffIfdToList( _tiffIfdList *l, _tiffIfd *p )
{
    char *proc = "_addTiffIfdToList";
    int i, s = l->n_allocated_data;
    _tiffIfd *data;

    if ( l->n_data == l->n_allocated_data ) {
      s += _size_to_be_allocated_;
      data = (_tiffIfd*)malloc( s * sizeof(_tiffIfd) );
      if ( data == (_tiffIfd*)NULL ) {
        if ( _verbose_ || _debug_ )
            fprintf( stderr, "%s: allocation error\n", proc );
        return( -1 );
      }
      if ( l->n_allocated_data > 0 ) {
        (void)memcpy( data, l->data, l->n_allocated_data*sizeof(_tiffIfd) );
        free( l->data );
      }
      for ( i=l->n_allocated_data; i<s; i++ )
        _initTiffIfd( &(data[i]) );
      l->n_allocated_data = s;
      l->data = data;
    }

    l->data[l->n_data] = *p;
    l->n_data ++;

    return( 1 );
}



static int _allocTiffIfdList( _tiffIfdList *l, int n )
{
    char *proc = "_allocTiffIfdList";
    int i;

    if ( n <= 0 ) return( 0 );
    l->data = (_tiffIfd*)malloc( n * sizeof(_tiffIfd) );
    if ( l->data == (_tiffIfd*)NULL ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: allocation failed\n", proc );
      return( -1 );
    }
    for ( i=0; i<n; i++ )
        _initTiffIfd( &(l->data[i]) );
    l->n_allocated_data = n;
    l->n_data = n;
    return( 1 );
}



/************************************************************
 *
 * _tiffIfdPtrList:
 *
 ************************************************************/



typedef struct _tiffIfdPtrList {
  _tiffIfd **data;
  int n_data;
  int n_allocated_data;
} _tiffIfdPtrList;



static void _initTiffIfdPtrList( _tiffIfdPtrList *l )
{
  l->data = (_tiffIfd**)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;
}



static void _freeTiffIfdPtrList( _tiffIfdPtrList *l )
{
  if ( l->data != (_tiffIfd**)NULL ) {
    free ( l->data );
  }
  _initTiffIfdPtrList( l );
}



static int _addTiffIfdToTiffIfdPtrList( _tiffIfdPtrList *l, _tiffIfd *p )
{
    char *proc = "_addTiffIfdToTiffIfdPtrList";
    int i, s = l->n_allocated_data;
    _tiffIfd **data;

    if ( l->n_data == l->n_allocated_data ) {
      s += _size_to_be_allocated_;
      data = (_tiffIfd**)malloc( s * sizeof(_tiffIfd*) );
      if ( data == (_tiffIfd**)NULL ) {
        if ( _verbose_ || _debug_ )
            fprintf( stderr, "%s: allocation error\n", proc );
        return( -1 );
      }
      if ( l->n_allocated_data > 0 ) {
        (void)memcpy( data, l->data, l->n_allocated_data*sizeof(_tiffIfd*) );
        free( l->data );
      }
      for ( i=l->n_allocated_data; i<s; i++ )
        data[i] = (_tiffIfd*)NULL;
      l->n_allocated_data = s;
      l->data = data;
    }

    l->data[l->n_data] = p;
    l->n_data ++;

    return( 1 );
}



typedef struct _tiffIfdPtrListList {
  _tiffIfdPtrList *data;
  int n_data;
  int n_allocated_data;
} _tiffIfdPtrListList;



static void _initTiffIfdPtrListList( _tiffIfdPtrListList *l )
{
  l->data = (_tiffIfdPtrList*)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;
}



static void _freeTiffIfdPtrListList( _tiffIfdPtrListList *l )
{
  int i;

  if ( l->data != (_tiffIfdPtrList*)NULL ) {
    for ( i=0;i<l->n_data; i++ )
      _freeTiffIfdPtrList( &(l->data[i]) );
    free ( l->data );
  }
  _initTiffIfdPtrListList( l );
}



static int _addTiffIfdToTiffIfdPtrListList( _tiffIfdPtrListList *l, _tiffIfd *p, int index )
{
  char *proc = "_addTiffIfdToTiffIfdPtrListList";
  int i, s = l->n_allocated_data;
  _tiffIfdPtrList *data;

  if ( index == -1 ) {

    if ( l->n_data == l->n_allocated_data ) {
      s += _size_to_be_allocated_;
      data = (_tiffIfdPtrList*)malloc( s * sizeof(_tiffIfdPtrList) );
      if ( data == (_tiffIfdPtrList*)NULL ) {
        if ( _verbose_ || _debug_ )
            fprintf( stderr, "%s: allocation error\n", proc );
        return( -1 );
      }
      if ( l->n_allocated_data > 0 ) {
        (void)memcpy( data, l->data, l->n_allocated_data*sizeof(_tiffIfdPtrList) );
        free( l->data );
      }
      for ( i=l->n_allocated_data; i<s; i++ )
        _initTiffIfdPtrList( &(data[i]) );
      l->n_allocated_data = s;
      l->data = data;
    }

    if ( _addTiffIfdToTiffIfdPtrList( &(l->data[ l->n_data ]),  p ) != 1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: unable to create list[%d]\n", proc, l->n_data );
      return( -1 );
    }

    l->n_data ++;

  }
  else if ( 0 <= index && index < l->n_data ) {
    if ( _addTiffIfdToTiffIfdPtrList( &(l->data[ index ]),  p ) != 1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: unable to add ifd to list[%d]\n", proc, index );
      return( -1 );
    }
  }
  else {
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: weird, this case should not occur\n", proc );
    return( -1 );
  }

  return( 1 );
}




static int _sortTiffIfdList( _tiffIfdPtrListList *resList, _tiffIfdList *theList )
{
  char *proc = "_sortiffIfdList";
  int i, j, listIndex;
  _tiffIfd *ifd0, *ifd;

  for ( i=0; i<theList->n_data; i++ ) {
    ifd = &(theList->data[i]);

    /* find a list with similar properties
     */
    for ( listIndex=-1, j=0; j<resList->n_data && listIndex==-1; j++ ) {
      if ( resList->data[j].n_data == 0 || resList->data[j].data == (_tiffIfd**)NULL )\
        continue;
      ifd0 = resList->data[j].data[0];
      if ( _TestIfdCoherency( ifd, ifd0, 0 ) == 1 ) {
        listIndex = j;
      }
    }

    /* create a new list or add to list
     */
    if ( _addTiffIfdToTiffIfdPtrListList( resList, ifd, listIndex ) != 1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: unable to add ifd #%d to pointer list\n", proc, i );
      return( -1 );
    }
  }
  return( 1 );
}



static int _selectIfdSubList( _tiffIfdPtrListList *theList )
{
  char *proc = "_selectIfdSubList";
  int i, e, xdim, ydim, xydim;
  int s=0;
  _tiffIfd *ifdptr;

  if ( theList->n_data <= 0 ) {
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: weird, this should not occur\n", proc );
    return( 0 );
  }
  if ( theList->n_data == 1 ) {
    return( 0 );
  }

  ifdptr = theList->data[0].data[0];
  for ( xdim=0, ydim=0, e=0; e<ifdptr->n_data; e++ ) {
    switch( *(ifdptr->data[e].tag) ) {
    default :
      break;
    case 256 :
      xdim = *(ifdptr->data[e].offset); break;
    case 257 :
      ydim = *(ifdptr->data[e].offset); break;
    }
  }
  xydim = xdim*ydim;

  for (i=1; i<theList->n_data; i++ ) {
    ifdptr = theList->data[i].data[0];
    for ( xdim=0, ydim=0, e=0; e<ifdptr->n_data; e++ ) {
      switch( *(ifdptr->data[e].tag) ) {
      default :
        break;
      case 256 :
        xdim = *(ifdptr->data[e].offset); break;
      case 257 :
        ydim = *(ifdptr->data[e].offset); break;
      }
    }
    if ( xydim < xdim*ydim ) {
      xydim = xdim*ydim;
      s = i;
    }
  }

  return( s );
}



/************************************************************
 *
 *
 *
 ************************************************************/

static off_t getFileSize( const char *name )
{
  char *proc = "getFileSize";
  struct stat stbuf;

  if ( name == (char*)NULL ) return( -1 );
  if ( stat( name, &stbuf ) != 0 ) {
    fprintf( stderr, "%s: unable to fill '%s' file info \n", proc, name );
    return( -1 );
  }
  return( stbuf.st_size );
}



/************************************************************
 *
 *
 *
 ************************************************************/



typedef struct _tiffTag {
  unsigned short int tag;
  char name[64];
} _tiffTag;


_tiffTag tagList[] =
{ { 254, "NewSubfileType" },
  { 255, "SubfileType" },
  { 256, "ImageWidth" },
  { 257, "ImageLength" },
  { 258, "BitsPerSample" },
  { 259, "Compression" },
  { 262, "PhotometricInterpretation" },
  { 263, "Thresholding" },
  { 264, "CellWidth" },
  { 265, "CellLength" },
  { 266, "FillOrder" },
  { 269, "DocumentName" },
  { 270, "ImageDescription" },
  { 271, "Make" },
  { 272, "Model" },
  { 273, "StripOffsets" },
  { 274, "Orientation" },
  { 277, "SamplesPerPixel" },
  { 278, "RowsPerStrip" },
  { 279, "StripByteCounts" },
  { 280, "MinSampleValue" },
  { 281, "MaxSampleValue" },
  { 282, "XResolution" },
  { 283, "YResolution" },
  { 284, "PlanarConfiguration" },
  { 285, "PageName" },
  { 286, "XPosition" },
  { 287, "YPosition" },
  { 288, "FreeOffsets" },
  { 289, "FreeByteCounts" },
  { 290, "GrayResponseUnit" },
  { 291, "GrayResponseCurve" },
  { 292, "T4Options" },
  { 293, "T6Options" },
  { 296, "ResolutionUnit" },
  { 297, "PageNumber" },
  { 301, "TransferFunction" },
  { 305, "Software" },
  { 306, "DateTime" },
  { 315, "Artist" },
  { 316, "HostComputer" },
  { 317, "Predictor" },
  { 318, "WhitePoint" },
  { 319, "PrimaryChromaticities" },
  { 320, "ColorMap" },
  { 321, "HalftoneHints" },
  { 322, "TileWidth" },
  { 323, "TileLength" },
  { 324, "TileOffsets" },
  { 325, "TileByteCounts" },
  { 332, "InkSet" },
  { 333, "InkNames" },
  { 334, "NumberOfInks" },
  { 336, "DotRange" },
  { 337, "TargetPrinter" },
  { 338, "ExtraSamples" },
  { 339, "SampleFormat" },
  { 340, "SMinSampleValue" },
  { 341, "SMaxSampleValue" },
  { 342, "TransferRange" },
  { 512, "JPEGProc" },
  { 513, "JPEGInterchangeFormat" },
  { 514, "JPEGInterchangeFormatLngth" },
  { 515, "JPEGRestartInterval" },
  { 517, "JPEGLosslessPredictors" },
  { 518, "JPEGPointTransforms" },
  { 519, "JPEGQTables" },
  { 520, "JPEGDCTables" },
  { 521, "JPEGACTables" },
  { 529, "YCbCrCoefficients" },
  { 530, "YCbCrSubSampling" },
  { 531, "YCbCrPositioning" },
  { 532, "ReferenceBlackWhite" },
  { 33432, "Copyright" },
  { 34412, "LSM header" },
  { 0, "unexpected tag" }
};

_tiffTag fieldTypeList[] = {
  { 1, "BYTE" },
  { 2, "ASCII" },
  { 3, "SHORT" },
  { 4, "LONG" },
  { 5, "RATIONAL" },
  { 6, "SBYTE" },
  { 7, "UNDEFINED" },
  { 8, "SSHORT" },
  { 9, "SLONG" },
  { 10, "SRATIONAL" },
  { 11, "FLOAT" },
  { 12, "DOUBLE" },
  { 0, "unexpected field type" }
};

/* The field types and their sizes are:
 *
1 = BYTE 8-bit unsigned integer.
2 = ASCII 8-bit byte that contains a 7-bit ASCII code; the last byte
must be NUL (binary zero).
3 = SHORT 16-bit (2-byte) unsigned integer.
4 = LONG 32-bit (4-byte) unsigned integer.
5 = RATIONAL Two LONGs: the first represents the numerator of a
fraction; the second, the denominator.

The value of the Count part of an ASCII field entry includes the NUL. If padding
is necessary, the Count does not include the pad byte. Note that there is no initial
“count byte” as in Pascal-style strings.

Any ASCII field can contain multiple strings, each terminated with a NUL. A
single string is preferred whenever possible. The Count for multi-string fields is
the number of bytes in all the strings in that field plus their terminating NUL
bytes. Only one NUL is allowed between strings, so that the strings following the
first string will often begin on an odd byte.

The reader must check the type to verify that it contains an expected value. TIFF
currently allows more than 1 valid type for some fields. For example, ImageWidth
and ImageLength are usually specified as having type SHORT. But images with
more than 64K rows or columns must use the LONG field type.

TIFF readers should accept BYTE, SHORT, or LONG values for any unsigned
integer field. This allows a single procedure to retrieve any integer value, makes
reading more robust, and saves disk space in some situations.

In TIFF 6.0, some new field types have been defined:

6 = SBYTE An 8-bit signed (twos-complement) integer.
7 = UNDEFINED An 8-bit byte that may contain anything, depending on
the definition of the field.
8 = SSHORT A 16-bit (2-byte) signed (twos-complement) integer.
9 = SLONG A 32-bit (4-byte) signed (twos-complement) integer.
10 = SRATIONAL Two SLONG’s: the first represents the numerator of a
fraction, the second the denominator.
11 = FLOAT Single precision (4-byte) IEEE format.
12 = DOUBLE Double precision (8-byte) IEEE format.

These new field types are also governed by the byte order (II or MM) in the TIFF
header.

Warning: It is possible that other TIFF field types will be added in the future.
Readers should skip over fields containing an unexpected field type.
*/



static void _fprintfTagDescription( FILE *f, int tag )
{
  int t;
  for ( t=0; tagList[t].tag != 0 && tagList[t].tag != tag; t++ )
      ;
  fprintf( f, "%26s ", tagList[t].name );
}



static void _fprintfFieldTypeDescription( FILE *f, int type )
{
  int t;
  for ( t=0; fieldTypeList[t].tag != 0 && fieldTypeList[t].tag != type; t++ )
      ;
  fprintf( f, "%9s ", fieldTypeList[t].name );
}



static void _fprintTiffIfdEntryValue( FILE *f, char *filebuffer,
                                      _tiffIfdEntry *e, int swap )
{
  char *proc = "_fprintTiffIfdEntryValue";
  unsigned int l;
  char *buf = (char*)NULL;
  unsigned char *u8;
  char *s8;
  unsigned short int *u16;
  unsigned int *u32;
  _LsmHeader lsmHeader;

  if ( *(e->tag) == 34412 ) {
    _initLsmHeader( &lsmHeader );
    memcpy( &lsmHeader, &(filebuffer[*(e->offset)]), sizeof(lsmHeader) );
    if ( swap ) _swapLsmHeader( &lsmHeader );
    fprintf( f, "\n" );
    _fprintfLsmHeader( f, &lsmHeader );
    return;
  }

  switch( *(e->fieldType) ) {
  default :
    fprintf( f, "type %d not handled yet", *(e->fieldType) ); break;
  /* BYTE
   */
  case 1 :
    if ( *(e->length) <= 4 ) {
      u8 = (unsigned char*)&(e->data[8]);
      for (l=0; l<*(e->length); l++ ) {
        fprintf( f, "%d ", u8[l] );
      }
    }
    else {
      fprintf( f, "not handled yet" );
    }
    break;
  /* ASCII
   */
  case 2 :
    if ( *(e->length) <= 4 ) {
      s8 = (char*)&(e->data[8]);
      for (l=0; l<*(e->length); l++ ) {
        fprintf( f, "%c", s8[l] );
      }
    }
    else {
      if ( filebuffer == (char*)NULL && e->values == (char*)NULL ) {
        fprintf( f, "NULL buffers" );
        break;
      }

      if ( filebuffer != (char*)NULL )
        s8 = (char*)&(filebuffer[*(e->offset)]);
      else
        s8 = e->values;

      for (l=0; l<*(e->length); l++ ) {
          fprintf( f, "%c", s8[l] );
      }
    }
    break;
  /* 3 = SHORT 16-bit (2-byte)
   */
  case 3 :
    if ( *(e->length) <= 2 ) {
      u16 = (unsigned short int*)&(e->data[8]);
      for (l=0; l<*(e->length); l++ ) {
        fprintf( f, "%d ", u16[l] );
      }
    }
    else {
      if ( filebuffer == (char*)NULL && e->values == (char*)NULL ) {
        fprintf( f, "NULL buffers" );
        break;
      }

      if ( filebuffer != (char*)NULL ) {
        buf = (char*)malloc( *(e->length) * 2 * sizeof(char) );
        if ( buf == (char*)NULL ) {
          if ( _verbose_ || _debug_ ) {
            fprintf( stderr, "%s: allocation failed for printing ", proc );
            _fprintfFieldTypeDescription( stderr, *(e->fieldType) );
            fprintf( stderr, "values\n" );
          }
          return;
        }
        memcpy( buf, &filebuffer[*(e->offset)], *(e->length) * 2 * sizeof(char) );
        if ( swap ) {
          for (l=0; l<*(e->length); l++ ) {
            _swap2bytes( (unsigned char *)&(buf[2*l]) );
          }
        }
        u16 = (unsigned short int*)buf;
      }
      else {
        u16 = (unsigned short int*)e->values;
      }

      if ( *(e->length) < 1000 ) {
        for (l=0; l<*(e->length); l++ ) {
          fprintf( f, "%d ", u16[l] );
        }
      }
      else {
        for (l=0; l<500; l++ ) {
          fprintf( f, "%d ", u16[l] );
        }
        fprintf( f, "...... " );
        for (l=*(e->length)-500; l<*(e->length); l++ ) {
          fprintf( f, "%d ", u16[l] );
        }
      }

      if ( filebuffer != (char*)NULL )
        free( buf );
    }
    break;
  /* 4 = LONG 32-bit (4-byte)
   */
  case 4 :
    if ( *(e->length) == 1 ) {
      fprintf( f, "%u", *(e->offset) );
    }
    else {
      if ( filebuffer == (char*)NULL && e->values == (char*)NULL ) {
        fprintf( f, "NULL buffers" );
        break;
      }

      if ( filebuffer != (char*)NULL ) {
        buf = (char*)malloc( *(e->length) * 4 * sizeof(char) );
        if ( buf == (char*)NULL ) {
          if ( _verbose_ || _debug_ ) {
            fprintf( stderr, "%s: allocation failed for printing ", proc );
            _fprintfFieldTypeDescription( stderr, *(e->fieldType) );
            fprintf( stderr, "values\n" );
          }
          return;
        }
        memcpy( buf, &filebuffer[*(e->offset)], *(e->length) * 4 * sizeof(char) );
        if ( swap ) {
          for (l=0; l<*(e->length); l++ ) {
            _swap4bytes( (unsigned char *)&(buf[4*l]) );
          }
        }
        u32 = (unsigned int *)buf;
      }
      else {
        u32 = (unsigned int*)e->values;
      }

      for (l=0; l<*(e->length); l++ ) {
        fprintf( f, "%u ", u32[l] );
      }

      if ( filebuffer != (char*)NULL )
        free( buf );
    }
    break;
  /* 5 = RATIONAL Two LONGs:
   */
  case 5 :
    if ( filebuffer == (char*)NULL && e->values == (char*)NULL ) {
      fprintf( f, "NULL buffers" );
      break;
    }

    if ( filebuffer != (char*)NULL && e->values == (char*)NULL ) {
      buf = (char*)malloc( 2 * 4 * sizeof(char) );
      if ( buf == (char*)NULL ) {
        if ( _verbose_ || _debug_ ) {
          fprintf( stderr, "%s: allocation failed for printing ", proc );
          _fprintfFieldTypeDescription( stderr, *(e->fieldType) );
          fprintf( stderr, "value\n" );
        }
        return;
      }
      memcpy( buf, &filebuffer[*(e->offset)], 8 );
      if ( swap ) {
        _swap4bytes( (unsigned char *)&(buf[0]) );
        _swap4bytes( (unsigned char *)&(buf[4]) );
      }
      u32 = (unsigned int *)buf;
    }
    else {
      u32 = (unsigned int*)e->values;
    }

    fprintf( f, "%u / %u", u32[0], u32[1] );

    if ( filebuffer != (char*)NULL )
        free( buf );
    break;
  }
}



static void _fprintTiffIfdEntry( FILE *f, char *filebuffer,
                                 _tiffIfdEntry *e, int swap )
{
  fprintf( f, "tag = %5d type = %1d length = %4u offset = %8u -- ",
           *(e->tag), *(e->fieldType), *(e->length), *(e->offset) );
  _fprintfTagDescription( f, *(e->tag) );
  _fprintfFieldTypeDescription( f, *(e->fieldType) );
  _fprintTiffIfdEntryValue( f, filebuffer, e, swap );
  fprintf( f, "\n" );
}



static void _fprintTiffIfd( FILE *f, char *filebuffer, _tiffIfd *ifd, int swap )
{
  int i;

  for ( i=0; i<ifd->n_data; i++ ) {
    fprintf( f, " %3d -- ", i );
    _fprintTiffIfdEntry( f, filebuffer, &(ifd->data[i]), swap );
  }
}



static void _fprintTiffIfdList( FILE *f, char *filebuffer, _tiffIfdList *l, int swap )
{
  int i;

  for ( i=0; i<l->n_data; i++ ) {
    fprintf( f, "--- IFD #%d : %d entries -------------------------\n", i, l->data[i].n_data );
    _fprintTiffIfd( f, filebuffer, &(l->data[i]), swap );
    fprintf( f, " --- offset of next IFD = %u ---------------------\n", l->data[i].offset );
  }
}

static void _fprintTiffIfdPtrListList( FILE *f, char *filebuffer, _tiffIfdPtrListList *l, int swap )
{
  int i;

  fprintf( f, "\n" );
  for ( i=0; i<l->n_data; i++ ) {
    fprintf( f, "=== IFD list #%d : %d entries =========================\n", i, l->data[i].n_data );
    fprintf( f, "   0 == \n" );
    _fprintTiffIfd( f, filebuffer, l->data[i].data[0], swap );
  }
  fprintf( f, "==================== ==================================\n" );
}



/************************************************************
 *
 *
 *
 ************************************************************/



int testLiteTiffHeader( char *magic __attribute__ ((unused)),
                   const char *name )
{
    /*
    if( !memcmp(magic,TIFF_LE_MAGIC,4) ||
        !memcmp(magic,TIFF_BE_MAGIC,4))
      return 0;
    */
    if ( (!strncmp(name+strlen(name)-4, ".tif", 4))
         || (!strncmp(name+strlen(name)-4, ".TIF", 4))
         || (!strncmp(name+strlen(name)-5, ".tiff", 5))
         || (!strncmp(name+strlen(name)-5, ".TIFF", 5))
         || (!strncmp(name+strlen(name)-4, ".lsm", 4)) )
      return 0;
    else
      return -1;
}






/************************************************************
 *
 *
 *
 ************************************************************/



/*
 * Image file header

A TIFF file begins with an 8-byte image file header, containing the
following information:

Bytes 0-1:	The first word of the file specifies the byte order
used within the file.  Legal values are:

  _II_	(hex 4949)
  _MM_	(hex 4D4D)

In the _II_  format, byte order is always from least
significant to most significant, for both 16-bit and 32-bit integers.
In the _MM_  format, byte order is always from most significant to
least significant, for both 16-bit and 32-bit integers.  In both
formats, character strings are stored into sequential byte locations.

All TIFF readers should support both byte orders_see Appendix G.

Bytes 2-3	The second word of the file is the TIFF _version
number._  This number, 42 (2A in hex), is not to be equated with the
current Revision of the TIFF specification.  In fact, the TIFF version
number (42) has never changed, and probably never will.  If it ever
does, it means that TIFF has changed in some way so radical that a
TIFF reader should give up immediately.  The number 42 was chosen for
its deep philosophical significance.  It can and should be used as
additional verification that this is indeed a TIFF file.

A TIFF file does not have a real version/revision number.
This was an explicit, conscious design decision.  In many file
formats, fields take on different meanings depending on a single
version number.  The problem is that as the file format _ages,_ it
becomes increasingly difficult to document which fields mean what in a
given version, and older software usually has to give up if it
encounters a file with a newer version number.  We wanted TIFF fields
to have a permanent and well-defined meaning, so that _older_ software
can usually read _newer_ TIFF files.  The bottom line is lower
software release costs and more reliable software.

Bytes 4-7	This long word contains the offset (in bytes) of the
first Image File Directory.  The directory may be at any location in
the file after the header but must begin on a word boundary.  In
particular, an Image File Directory may follow the image data it
describes.  Readers must simply follow the pointers, wherever they may
lead.

(The term _byte offset_ is always used in this document to
refer to a location with respect to the beginning of the file.  The
first byte of the file has an offset of 0.)
*/

static int _parseHeader( _tiffHeader *h )
{
  char *proc = "_parseHeader";
  ENDIANNESS hendianess = END_UNKNOWN;

  if ( h->data[0] == 73 && h->data[1] == 73 ) {
    hendianess = END_LITTLE;
  }
  else if ( h->data[0] == 77 && h->data[1] == 77 ) {
    hendianess = END_BIG;
  }
  else {
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: TIFF header not recognized (neither II or MM)\n", proc );
    return( -1 );
  }

  h->swap = ( _getEndianness() != hendianess );

  if ( h->swap ) {
    _swap2bytes( (unsigned char *)&(h->data[2]) );
    _swap4bytes( (unsigned char *)&(h->data[4]) );
  }

  if ( *(h->version) != 42 ) {
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: TIFF header not recognized (bad version mumber: %d)\n",
               proc, *(h->version) );
    return( -1 );
  }

  if ( _debug_ ) {
    switch( hendianess ) {
    default :
      fprintf( stderr, "%s: weird, this should not be reached\n", proc );
      break;
    case END_LITTLE :
      fprintf( stderr, "_II_ (little endian) tiff file" );
      break;
    case END_BIG :
      fprintf( stderr, "_MM_ (big endian) tiff file" );
      break;
    }
    if ( h->swap ) fprintf( stderr, ", swap is on\n" );
    else fprintf( stderr, ", swap is off\n" );
    fprintf( stderr, "     version = %d --- offset = %u\n",
             *(h->version), *(h->offset) );
  }

  return( 1 );
}



/*
Image file directory

An Image File Directory (IFD) consists of a 2-byte count of the number
of entries (i.e., the number of fields), followed by a sequence of
12-byte field entries, followed by a 4-byte offset of the next Image
File Directory (or 0 if none).  Do not forget to write the 4 bytes of
0 after the last IFD.

Each 12-byte IFD entry has the following format:

Bytes 0-1	contain the Tag for the field.
Bytes 2-3	contain the field Type.
Bytes 4-7	contain the Length (_Count_ might have been a better
term) of the field.
Bytes 8-11	contain the Value Offset, the file offset (in bytes)
of the Value for the field.  The Value is expected to begin on a word
boundary; the corresponding Value Offset will thus be an even number.
This file offset may point to anywhere in the file, including after
the image data.

The entries in an IFD must be sorted in ascending order by Tag.  Note
that this is not the order in which the fields are described in this
document.  For a numerically ordered list of tags, see Appendix E.
The Values to which directory entries point need not be in any
particular order in the file.

In order to save time and space, the Value Offset is interpreted to
contain the Value instead of pointing to the Value if the Value fits
into 4 bytes.  If the Value is less than 4 bytes, it is left-justified
within the 4-byte Value Offset, i.e., stored in the lower-numbered
bytes.  Whether or not the Value fits within 4 bytes is determined by
looking at the Type and Length of the field.

The Length is specified in terms of the data type, not the total
number of bytes.  A single 16-bit word (SHORT) has a Length of 1, not
2, for example.  The data types and their lengths are described below:

1 = BYTE	An 8-bit unsigned integer.
2 = ASCII	8-bit bytes that store ASCII codes; the last byte must be null.
3 = SHORT	A 16-bit (2-byte) unsigned integer.
4 = LONG	A 32-bit (4-byte) unsigned integer.
5 = RATIONAL	Two LONG_s:  the first represents the numerator of a
fraction, the second the denominator.

The value of the Length part of an ASCII field entry includes the
null.  If padding is necessary, the Length does not include the pad
byte.  Note that there is no _count byte,_ as there is in Pascal-type
strings.  The Length part of the field takes care of that.  The null
is not strictly necessary, but may make things slightly simpler for C
programmers.

The reader should check the type to ensure that it is what he expects.
TIFF currently allows more than 1 valid type for some fields.  For
example, ImageWidth and ImageLength were specified as having type
SHORT.  Very large images with more than 64K rows or columns are
possible with some devices even now.  Rather than add parallel LONG
tags for these fields, it is cleaner to allow both SHORT and LONG for
ImageWidth and similar fields.  See Appendix G for specific
recommendations.

Note that there may be more than one IFD.  Each IFD is said to define
a _subfile._   One potential use of subsequent subfiles is to describe
a _sub-image_  that is somehow related to the main image, such as a
reduced resolution version of the image.

If you have not already done so, you may wish to turn to Appendix G to
study the sample TIFF images.
*/

static int _readIfd( char *thebuf, _tiffIfd *ifd, _tiffHeader *h )
{
  char *proc = "_readIfd";
  char *buf = thebuf;
  int i;

  memcpy( &(ifd->n_data), buf, 2 );
  if ( h->swap )
    _swap2bytes( (unsigned char *)&(ifd->n_data) );
  buf += 2;

  if ( ifd->n_data <= 0 ) {
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: no ifd entries ?!\n", proc );
    return( -1 );
  }

  ifd->data = (_tiffIfdEntry*)malloc( ifd->n_data * sizeof(_tiffIfdEntry) );
  if ( ifd->data == (_tiffIfdEntry*)NULL ) {
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
    ifd->n_data = 0;
  }
  ifd->n_allocated_data = ifd->n_data;

  for ( i=0; i<ifd->n_data; i++, buf+=12 ) {
    _initTiffIfdEntry( &(ifd->data[i]) );
    memcpy( ifd->data[i].data, buf, 12 );
    if ( h->swap ) {
      _swap2bytes( (unsigned char *)&(ifd->data[i].data[0]) );
      _swap2bytes( (unsigned char *)&(ifd->data[i].data[2]) );
      _swap4bytes( (unsigned char *)&(ifd->data[i].data[4]) );
      _swapValue( &(ifd->data[i]) );
    }
  }

  memcpy( &(ifd->offset), buf, 4 );
  if ( h->swap )
    _swap4bytes( (unsigned char *)&(ifd->offset) );

  return( 1 );
}



static int _readUncompressedDataFromIfdPtrList( char *filebuffer, _tiffIfdPtrList *ifdList, _tiffHeader *header, _image *im )
{
  char *proc = "_readUncompressedDataFromIfdPtrList";

  _tiffIfd *ifdptr;
  int e, i, s, ns;
  unsigned int *u32StripOffsets;
  unsigned int *u32StripByteCounts;
  char *dataptr;
  int eStripOffsets;
  int eStripByteCounts;

  /* get data
   */
  dataptr = (char*)im->data;
  for ( i=0; i<ifdList->n_data; i++ ) {
    ifdptr = ifdList->data[i];

    for ( e=0, eStripOffsets = eStripByteCounts = -1; e<ifdptr->n_data && (eStripOffsets==-1 || eStripByteCounts==-1); e++ ) {
      switch( *(ifdptr->data[e].tag) ) {
      default  : break;
      case 273 : eStripOffsets = e;    break;
      case 279 : eStripByteCounts = e; break;
      }
    }

    if ( *(ifdptr->data[eStripOffsets].length) != *(ifdptr->data[eStripByteCounts].length) ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: different value length for StripOffsets (%d) and StripByteCounts (%d)\n",
                 proc, *(ifdptr->data[eStripOffsets].length), *(ifdptr->data[eStripByteCounts].length) );
      return( -1 );
    }

    ns = *(ifdptr->data[eStripOffsets].length);

    switch( *(ifdptr->data[eStripOffsets].fieldType) ) {
    default :
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: field type for StripOffsets (%d) not handled yet\n",
                 proc, *(ifdptr->data[eStripOffsets].fieldType) );
      return( -1 );
      /* end of *(ifdptr->data[eStripOffsets].fieldType) = default */
    case 3 :
      /* SHORT */
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: field type for StripOffsets (%d) not implemented yet\n",
                 proc, *(ifdptr->data[eStripOffsets].fieldType) );
      return( -1 );
      /* end of *(ifdptr->data[eStripOffsets].fieldType) = SHORT */
    case 4 :
      /* LONG */
      if ( ns == 1 ) {
        u32StripOffsets = ifdptr->data[eStripOffsets].offset;
        switch( *(ifdptr->data[eStripByteCounts].fieldType) ) {
        default :
          if ( _verbose_ || _debug_ )
            fprintf( stderr, "%s: field type for StripByteCounts (%d) not handled yet\n",
                     proc, *(ifdptr->data[eStripOffsets].fieldType) );
          return( -1 );
        case 4 :
          u32StripByteCounts = ifdptr->data[eStripByteCounts].offset;
          (void)memcpy( dataptr, &(filebuffer[*u32StripOffsets]), *u32StripByteCounts );
          dataptr += *u32StripByteCounts;
          break;
        }
      }
      /* ns > 1
       */
      else {
        u32StripOffsets = (unsigned int *)malloc( ns * sizeof(unsigned int) );
        if ( u32StripOffsets ==  (unsigned int *)NULL ) {
          if ( _verbose_ || _debug_ )
            fprintf( stderr, "%s: StripOffsets allocation failed\n", proc );
          return( -1 );
        }
        memcpy( u32StripOffsets, &(filebuffer[*(ifdptr->data[eStripOffsets].offset)]), ns*sizeof(unsigned int) );
        if ( header->swap ) {
          for ( s=0; s<ns; s++ )
            _swap4bytes( (unsigned char*)&(u32StripOffsets[s]) );
        }
        switch( *(ifdptr->data[eStripByteCounts].fieldType) ) {
        default :
          free( u32StripOffsets );
          if ( _verbose_ || _debug_ )
            fprintf( stderr, "%s: field type for StripByteCounts (%d) not handled yet\n",
                     proc, *(ifdptr->data[eStripOffsets].fieldType) );
          return( -1 );
        case 4 :
          u32StripByteCounts = (unsigned int *)malloc( ns * sizeof(unsigned int) );
          if ( u32StripByteCounts ==  (unsigned int *)NULL ) {
            free( u32StripOffsets );
            if ( _verbose_ || _debug_ )
              fprintf( stderr, "%s: StripByteCounts allocation failed\n", proc );
            return( -1 );
          }
          memcpy( u32StripByteCounts, &(filebuffer[*(ifdptr->data[eStripByteCounts].offset)]), ns*sizeof(unsigned int) );
          if ( header->swap ) {
            for ( s=0; s<ns; s++ )
              _swap4bytes( (unsigned char*)&(u32StripByteCounts[s]) );
          }
          for ( s=0; s<ns; s++ ) {
            (void)memcpy( dataptr, &(filebuffer[u32StripOffsets[s]]), u32StripByteCounts[s] );
            dataptr += u32StripByteCounts[s];
          }
          free( u32StripByteCounts );
          break;
        }
        free( u32StripOffsets );
      }
      break;
      /* end of *(ifdptr->data[eStripOffsets].fieldType) = LONG */
    }

  }
  /* end of getting data
   */

  return( 1 );
}



static unsigned short int *_extrapolatedColorMap( unsigned short int *u16ColorMap, unsigned int theColorLength,  unsigned int resColorLength )
{
  char *proc = "_extrapolatedColorMap";
  unsigned short int *newColorMap = (unsigned short int*)NULL;
  unsigned short int *u16RColor;
  unsigned short int *u16GColor;
  unsigned short int *u16BColor;
  unsigned short int *u16newRColor;
  unsigned short int *u16newGColor;
  unsigned short int *u16newBColor;
  int i;
  int ipc, inc;
  float c, dc;

  newColorMap = (unsigned short int *)malloc( resColorLength * 3 * sizeof(unsigned short int) );
  if ( newColorMap == (unsigned short int*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( (unsigned short int*)NULL );
  }

  u16RColor = u16ColorMap;
  u16GColor = u16RColor;
  u16BColor = u16RColor;
  u16GColor += theColorLength;
  u16BColor += theColorLength + theColorLength;

  u16newRColor = newColorMap;
  u16newGColor = u16newRColor;
  u16newBColor = u16newGColor;
  u16newGColor += resColorLength;
  u16newBColor += resColorLength + resColorLength;

  for ( i=0; i<resColorLength; i++ ) {
    c = (float)(theColorLength - 1) * (float)i / (float)(resColorLength-1);
    ipc = (int)c;
    if ( ipc <= 0 ) {
      u16newRColor[i] = u16RColor[0];
      u16newGColor[i] = u16GColor[0];
      u16newBColor[i] = u16BColor[0];
    }
    else if ( ipc >= theColorLength - 1 ) {
      u16newRColor[i] = u16RColor[theColorLength - 1];
      u16newGColor[i] = u16GColor[theColorLength - 1];
      u16newBColor[i] = u16BColor[theColorLength - 1];
    }
    else {
      dc = c - ipc;
      inc = ipc+1;
      u16newRColor[i] = (int)( (1.0-dc)*(float)u16RColor[ipc] + dc*(float)u16RColor[inc] + 0.5 );
      u16newGColor[i] = (int)( (1.0-dc)*(float)u16GColor[ipc] + dc*(float)u16GColor[inc] + 0.5 );
      u16newBColor[i] = (int)( (1.0-dc)*(float)u16BColor[ipc] + dc*(float)u16BColor[inc] + 0.5 );
    }
  }
  return( newColorMap );
}



static int _readUncompressedPaletteFromIfdPtrList( char *filebuffer, _tiffIfdPtrList *ifdList, _tiffHeader *header, _image *im )
{
  char *proc = "_readUncompressedPaletteFromIfdPtrList";
  _tiffIfd *ifdptr;
  int e, i, j, ns;
  unsigned int *u32StripOffsets;
  /* unsigned int *u32StripByteCounts; */
  unsigned short int *u16data = (unsigned short int*)im->data;

  int eStripOffsets = -1;
  int eStripByteCounts = -1;
  int eColorMap = -1;

  unsigned short int BitsPerSample[3] = { 1, 1, 1 };

  unsigned int colorLength;
  unsigned short int *u16ColorMap;
  unsigned short int *u16RColor;
  unsigned short int *u16GColor;
  unsigned short int *u16BColor;

  unsigned short int *u16TmpColorMap;

  int PlanarConfiguration = 1;
  size_t x, y;

  for ( i=0; i<ifdList->n_data; i++ ) {
    /* get some info
     */
    ifdptr = ifdList->data[i];
    for ( e=0; e<ifdptr->n_data; e++ ) {
      switch( *(ifdptr->data[e].tag) ) {
      default  : break;
      case 258 :
        BitsPerSample[0] = ((unsigned short int*)(ifdptr->data[e].offset))[0];
        break;
      case 273 : eStripOffsets = e;    break;
      case 279 : eStripByteCounts = e; break;
      case 284 :
        PlanarConfiguration = *(ifdptr->data[e].offset); break;
      case 320 : eColorMap = e;    break;
      }
    }

    /* some tests
     * StripOffsets field type is assumed to be LONG
     * StripByteCounts field type is assumed to be LONG
     * ColorMap field type is assumed to be SHORT
     *    Recall: image data are of the same type than the color map
     */
    if ( *(ifdptr->data[eStripOffsets].fieldType) != 4 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: field type for StripOffsets (%d) not handled yet\n",
                 proc, *(ifdptr->data[eStripOffsets].fieldType) );
      return( -1 );
    }
    if ( *(ifdptr->data[eStripByteCounts].fieldType) != 4 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: field type for StripByteCounts (%d) not handled yet\n",
                 proc, *(ifdptr->data[eStripByteCounts].fieldType) );
      return( -1 );
    }
    if ( *(ifdptr->data[eStripOffsets].length) != *(ifdptr->data[eStripByteCounts].length) ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: different value length for StripOffsets (%d) and StripByteCounts (%d)\n",
                 proc, *(ifdptr->data[eStripOffsets].length), *(ifdptr->data[eStripByteCounts].length) );
      return( -1 );
    }
    ns = *(ifdptr->data[eStripOffsets].length);
    if ( ns != 1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: case of several StripOffsets (%d) not implemented yet\n", proc, ns );
      return( -1 );
    }

    if ( *(ifdptr->data[eColorMap].fieldType) != 3 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: field type for ColorMap (%d) not handled yet\n",
                 proc, *(ifdptr->data[eColorMap].fieldType) );
      return( -1 );
    }
    if ( *(ifdptr->data[eColorMap].length) <= 0 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: no color map of negative or null length (%u)?!\n", proc, *(ifdptr->data[eColorMap].length) );
      return( -1 );
    }


    /* since the number of strips is 1,
     * get the values from the IFD offset
     * u32StripByteCounts = (unsigned int *)ifdptr->data[eStripByteCounts].offset;
     */
    u32StripOffsets = (unsigned int *)ifdptr->data[eStripOffsets].offset;


    /* color map
     */
    colorLength = *(ifdptr->data[eColorMap].length) / 3;
    if ( colorLength*3 != *(ifdptr->data[eColorMap].length) ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: weird length (%u): %u/3 * 3 = %d * 3 = %d != %u\n",
                 proc, *(ifdptr->data[eColorMap].length), *(ifdptr->data[eColorMap].length),
                 colorLength, colorLength*3, *(ifdptr->data[eColorMap].length) );
      return( -1 );
    }

    u16ColorMap = (unsigned short int*)malloc( colorLength*3 * sizeof(unsigned short int) );
    if ( u16ColorMap ==  (unsigned short int *)NULL ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: ColorMap allocation failed\n", proc );
      return( -1 );
    }
    memcpy( u16ColorMap, &(filebuffer[*(ifdptr->data[eColorMap].offset)]), colorLength*3 * sizeof(unsigned short int) );
    if ( 0 && header->swap ) {
      for ( j=0; j<(int)colorLength*3; j++ ) {
        _swap2bytes( (unsigned char*)&(u16ColorMap[j]) );
      }
    }
    u16RColor = u16ColorMap;
    u16GColor = u16RColor;
    u16BColor = u16RColor;
    u16GColor += colorLength;
    u16BColor += colorLength + colorLength;


    switch ( BitsPerSample[0] ) {
    default :
      free( u16ColorMap );
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: BitsPerSample[0] value (%d) not handled yet\n", proc, BitsPerSample[0] );
      return( -1 );
    case 8 :
      {
        unsigned char *ptrIndex = (unsigned char*)&(filebuffer[*u32StripOffsets]);
        unsigned int colorLengthThreshold = 256;
        if ( colorLength >= colorLengthThreshold ) {
          switch( PlanarConfiguration ) {
          default :
          case 1 :
            for ( j=0, y=0; y<im->ydim; y++ )
            for ( x=0; x<im->xdim; x++, j++ ) {
              u16data[(y*im->xdim+x)*im->vdim+0] = u16RColor[ptrIndex[j]];
              u16data[(y*im->xdim+x)*im->vdim+1] = u16GColor[ptrIndex[j]];
              u16data[(y*im->xdim+x)*im->vdim+2] = u16BColor[ptrIndex[j]];
            }
            break;
          case 2 :
            for ( j=0, y=0; y<im->ydim; y++ )
            for ( x=0; x<im->xdim; x++, j++ ) {
              u16data[(0*im->ydim+y)*im->xdim+x] = u16RColor[ptrIndex[j]];
              u16data[(1*im->ydim+y)*im->xdim+x] = u16GColor[ptrIndex[j]];
              u16data[(2*im->ydim+y)*im->xdim+x] = u16BColor[ptrIndex[j]];
            }
            break;
          }
        }
        else {
          free( u16ColorMap );
          if ( _verbose_ || _debug_ )
            fprintf( stderr, "%s: weird color length (%d)\n", proc, colorLength );
          return( -1 );
        }
      } /* case 8: */
      break;
    case 16 :
      {
        unsigned short int *ptrIndex = (unsigned short int *)&(filebuffer[*u32StripOffsets]);
        unsigned short int tmpIndex;
        unsigned int colorLengthThreshold = 65536;
        if ( colorLength >= colorLengthThreshold ) {
          switch( PlanarConfiguration ) {
          default :
          case 1 :
            for ( j=0, y=0; y<im->ydim; y++ )
            for ( x=0; x<im->xdim; x++, j++ ) {
              tmpIndex = ptrIndex[j];
              if ( header->swap ) _swap2bytes( (unsigned char*)&(tmpIndex) );
              u16data[(y*im->xdim+x)*im->vdim+0] = u16RColor[tmpIndex];
              u16data[(y*im->xdim+x)*im->vdim+1] = u16GColor[tmpIndex];
              u16data[(y*im->xdim+x)*im->vdim+2] = u16BColor[tmpIndex];
            }
            break;
          case 2 :
            for ( j=0, y=0; y<im->ydim; y++ )
            for ( x=0; x<im->xdim; x++, j++ ) {
              tmpIndex = ptrIndex[j];
              if ( header->swap ) _swap2bytes( (unsigned char*)&(tmpIndex) );
              u16data[(0*im->ydim+y)*im->xdim+x] = u16RColor[tmpIndex];
              u16data[(1*im->ydim+y)*im->xdim+x] = u16GColor[tmpIndex];
              u16data[(2*im->ydim+y)*im->xdim+x] = u16BColor[tmpIndex];
            }
            break;
          }
        }
        else {
          u16TmpColorMap = _extrapolatedColorMap( u16ColorMap, colorLength, colorLengthThreshold );
          if ( u16TmpColorMap == (unsigned short int*)NULL ) {
            free( u16ColorMap );
            if ( _verbose_ || _debug_ )
              fprintf( stderr, "%s: unable to extrapolate color map from %u to %u\n", proc, colorLength, colorLengthThreshold );
            return( -1 );
          }
          u16RColor = u16TmpColorMap;
          u16GColor = u16RColor;
          u16BColor = u16RColor;
          u16GColor += colorLengthThreshold;
          u16BColor += colorLengthThreshold + colorLengthThreshold;

          switch( PlanarConfiguration ) {
          default :
          case 1 :
            for ( j=0, y=0; y<im->ydim; y++ )
            for ( x=0; x<im->xdim; x++, j++ ) {
              tmpIndex = ptrIndex[j];
              if ( header->swap ) _swap2bytes( (unsigned char*)&(tmpIndex) );
              u16data[(y*im->xdim+x)*im->vdim+0] = u16RColor[tmpIndex];
              u16data[(y*im->xdim+x)*im->vdim+1] = u16GColor[tmpIndex];
              u16data[(y*im->xdim+x)*im->vdim+2] = u16BColor[tmpIndex];
            }
            break;
          case 2 :
            for ( j=0, y=0; y<im->ydim; y++ )
            for ( x=0; x<im->xdim; x++, j++ ) {
              tmpIndex = ptrIndex[j];
              if ( header->swap ) _swap2bytes( (unsigned char*)&(tmpIndex) );
              u16data[(0*im->ydim+y)*im->xdim+x] = u16RColor[tmpIndex];
              u16data[(1*im->ydim+y)*im->xdim+x] = u16GColor[tmpIndex];
              u16data[(2*im->ydim+y)*im->xdim+x] = u16BColor[tmpIndex];
            }
            break;
          }
          free( u16TmpColorMap );
        }
      } /* case 16 : */
      break;
    }

    u16data += im->xdim * im->ydim * 3;
    free( u16ColorMap );
  }

  return( 1 );
}



/* Loop until you get the number of unpacked bytes you are expecting:
 *   Read the next source byte into n.
 *   If n is between 0 and 127 inclusive, copy the next n+1 bytes literally.
 *   Else if n is between -127 and -1 inclusive, copy the next byte -n+1 times.
 *   Else if n is -128, noop.
 * Endloop
 */
static void _UnPackBitsStrip( char **buffer, size_t *bufferLength, char *strip, unsigned int stripLength )
{
  char *proc = "_UnPackBitsStrip";
  char *bufPtr = *buffer;
  char *stripPtr = strip;
  size_t bl = *bufferLength;
  unsigned int sl;
  int advance, n;

  for ( sl=stripLength; bl>0 && sl>0; ) {
    if ( 0 <= *stripPtr ) {
      advance = 1+(*stripPtr);
      stripPtr++;
      memcpy( bufPtr, stripPtr, advance );
      stripPtr += advance;
      bufPtr   += advance;
      sl -= (advance+1);
      bl -= advance;
    }
    else if ( -127 <= *stripPtr && *stripPtr <= -1 ) {
      advance = 1-(*stripPtr);
      stripPtr++;
      for ( n=0; n<advance; n++ ) {
        bufPtr[n] = *stripPtr;
      }
      stripPtr ++;
      bufPtr   += advance;
      sl -= 2;
      bl -= advance;
    }
    else if ( -128 == *stripPtr ) {
      stripPtr ++;
      sl --;
    }
    else {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: unexpected value\n", proc );
      stripPtr ++;
      sl --;
    }
    *bufferLength = bl;
    *buffer = bufPtr;
  }
}



static int _readPackBitsDataFromIfdPtrList( char *filebuffer, _tiffIfdPtrList *ifdList, _tiffHeader *header, _image *im )
{
  char *proc = "_readPackBitsDataFromIfdPtrList";

  _tiffIfd *ifdptr;
  int e, i, s, ns;
  unsigned int *u32StripOffsets;
  unsigned int *u32StripByteCounts;
  char *dataptr;
  size_t dataLength = im->xdim * im->ydim * im->zdim * (size_t)im->vdim * (size_t)im->wdim;
  int eStripOffsets;
  int eStripByteCounts;

  /* get data
   */
  dataptr = (char*)im->data;
  for ( i=0; i<ifdList->n_data; i++ ) {
    ifdptr = ifdList->data[i];
    eStripOffsets = eStripByteCounts = -1;

    for ( e=0; e<ifdptr->n_data; e++ ) {
      switch( *(ifdptr->data[e].tag) ) {
      default  : break;
      case 273 : eStripOffsets = e;    break;
      case 279 : eStripByteCounts = e; break;
      }
    }

    if ( *(ifdptr->data[eStripOffsets].length) != *(ifdptr->data[eStripByteCounts].length) ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: different value length for StripOffsets (%d) and StripByteCounts (%d)\n",
                 proc, *(ifdptr->data[eStripOffsets].length), *(ifdptr->data[eStripByteCounts].length) );
      return( -1 );
    }

    ns = *(ifdptr->data[eStripOffsets].length);

    switch( *(ifdptr->data[eStripOffsets].fieldType) ) {
    default :
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: field type for StripOffsets (%d) not handled yet\n",
                 proc, *(ifdptr->data[eStripOffsets].fieldType) );
      return( -1 );
      /* end of *(ifdptr->data[eStripOffsets].fieldType) = default */
    case 3 :
      /* SHORT */
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: field type for StripOffsets (%d) not implemented yet\n",
                 proc, *(ifdptr->data[eStripOffsets].fieldType) );
      return( -1 );
      /* end of *(ifdptr->data[eStripOffsets].fieldType) = SHORT */
    case 4 :
      /* LONG */
      if ( ns == 1 ) {
        u32StripOffsets = ifdptr->data[eStripOffsets].offset;
        switch( *(ifdptr->data[eStripByteCounts].fieldType) ) {
        default :
          if ( _verbose_ || _debug_ )
            fprintf( stderr, "%s: field type for StripByteCounts (%d) not handled yet\n",
                     proc, *(ifdptr->data[eStripOffsets].fieldType) );
          return( -1 );
        case 4 :
          u32StripByteCounts = ifdptr->data[eStripByteCounts].offset;
          _UnPackBitsStrip( &dataptr, &dataLength, &(filebuffer[*u32StripOffsets]), *u32StripByteCounts );
          break;
        }
      }
      /* ns > 1
       */
      else {
        u32StripOffsets = (unsigned int *)malloc( ns * sizeof(unsigned int) );
        if ( u32StripOffsets ==  (unsigned int *)NULL ) {
          if ( _verbose_ || _debug_ )
            fprintf( stderr, "%s: StripOffsets allocation failed\n", proc );
          return( -1 );
        }
        memcpy( u32StripOffsets, &(filebuffer[*(ifdptr->data[eStripOffsets].offset)]), ns*sizeof(unsigned int) );
        if ( header->swap ) {
          for ( s=0; s<ns; s++ )
            _swap4bytes( (unsigned char*)&(u32StripOffsets[s]) );
        }
        switch( *(ifdptr->data[eStripByteCounts].fieldType) ) {
        default :
          free( u32StripOffsets );
          if ( _verbose_ || _debug_ )
            fprintf( stderr, "%s: field type for StripByteCounts (%d) not handled yet\n",
                     proc, *(ifdptr->data[eStripOffsets].fieldType) );
          return( -1 );
        case 4 :
          u32StripByteCounts = (unsigned int *)malloc( ns * sizeof(unsigned int) );
          if ( u32StripByteCounts ==  (unsigned int *)NULL ) {
            free( u32StripOffsets );
            if ( _verbose_ || _debug_ )
              fprintf( stderr, "%s: StripByteCounts allocation failed\n", proc );
            return( -1 );
          }
          memcpy( u32StripByteCounts, &(filebuffer[*(ifdptr->data[eStripByteCounts].offset)]), ns*sizeof(unsigned int) );
          if ( header->swap ) {
            for ( s=0; s<ns; s++ )
              _swap4bytes( (unsigned char*)&(u32StripByteCounts[s]) );
          }
          for ( s=0; s<ns; s++ ) {
            _UnPackBitsStrip( &dataptr, &dataLength, &(filebuffer[u32StripOffsets[s]]), u32StripByteCounts[s] );
          }
          free( u32StripByteCounts );
          break;
        }
        free( u32StripOffsets );
      }
      break;
      /* end of *(ifdptr->data[eStripOffsets].fieldType) = LONG */
    }

  }
  /* end of getting data
   */

  return( 1 );
}



static int _getVoxelDimensionFromLiteTiff( char* buf, int length, _image *im )
{
  char *proc = "_getVoxelDimensionFromLiteTiff";
  int i;
  float s;

  if ( strncmp ( buf, "LiteTiff", 8 ) == 0 ) {
    for ( i=0; i<length-7; i++ ) {
      if ( buf[i] != 'V' ) continue;
      if ( strncmp( &(buf[i]), "VX=", 3 ) == 0 ) {
        i += 3;
        if ( sscanf( &(buf[i]) ,"%f", &s ) != 1 ) {
            if ( _verbose_ || _debug_ )
              fprintf( stderr, "%s: error when getting VX from ImageDescription\n", proc );
            return( -1 );
        }
        im->vx = s;
      }
      else if ( strncmp( &(buf[i]), "VY=", 3 ) == 0 ) {
        i += 3;
        if ( sscanf( &(buf[i]) ,"%f", &s ) != 1 ) {
            if ( _verbose_ || _debug_ )
              fprintf( stderr, "%s: error when getting VY from ImageDescription\n", proc );
            return( -1 );
        }
        im->vy = s;
      }
      else if ( strncmp( &(buf[i]), "VZ=", 3 ) == 0 ) {
        i += 3;
        if ( sscanf( &(buf[i]) ,"%f", &s ) != 1 ) {
            if ( _verbose_ || _debug_ )
              fprintf( stderr, "%s: error when getting VY from ImageDescription\n", proc );
            return( -1 );
        }
        im->vz = s;
      }
     }
    return( 1 );
  }
  return( 0 );
}



static int _getVoxelDimensionFromImageJ( char* buf, int length, _image *im )
{
  char *proc = "_getVoxelDimensionFromImageJ";
  int i;
  float s;

  if ( strncmp ( buf, "ImageJ=", 7 ) == 0 ) {
    for ( i=0; i<length-8; i++ ) {
      if ( buf[i] != 's' ) continue;
      if ( strncmp( &(buf[i]), "spacing=", 8 ) == 0 ) {
        i += 8;
        if ( sscanf( &(buf[i]) ,"%f", &s ) != 1 ) {
            if ( _verbose_ || _debug_ )
              fprintf( stderr, "%s: error when getting spacing from ImageDescription\n", proc );
            return( -1 );
        }
        im->vz = s;
        return( 1 );
      }
    }
  }
  return( 0 );
}



static int _getVoxelDimensionFromOmero( char* buf, int length, _image *im )
{
  char *proc = "_getVoxelDimensionFromOmero";
  int i, j;
  float s;

  for ( i=0; i<length-5; i++ ) {
    if ( buf[i] != '<' ) continue;
    if ( strncmp( &(buf[i]), "<OME ", 5) == 0 ) {
        i += 5;
        for ( j=i; j<length-16; j++ ) {
          if ( buf[j] != 'P' ) continue;
          if ( strncmp( &(buf[j]), "PhysicalSizeX=\"", 15 ) == 0 ) {
            j += 15;
            if ( sscanf( &(buf[j]) ,"%f", &s ) != 1 ) {
                if ( _verbose_ || _debug_ )
                  fprintf( stderr, "%s: error when getting PhysicalSizeX from ImageDescription\n", proc );
                return( -1 );
            }
            im->vx = s;
          }
          else if ( strncmp( &(buf[j]), "PhysicalSizeY=\"", 15 ) == 0 ) {
            j += 15;
            if ( sscanf( &(buf[j]) ,"%f", &s ) != 1 ) {
                if ( _verbose_ || _debug_ )
                  fprintf( stderr, "%s: error when getting PhysicalSizeY from ImageDescription\n", proc );
                return( -1 );
            }
            im->vy = s;
          }
          else if ( strncmp( &(buf[j]), "PhysicalSizeZ=\"", 15 ) == 0 ) {
            j += 15;
            if ( sscanf( &(buf[j]) ,"%f", &s ) != 1 ) {
                if ( _verbose_ || _debug_ )
                  fprintf( stderr, "%s: error when getting PhysicalSizeZ from ImageDescription\n", proc );
                return( -1 );
            }
            im->vz = s;
          }
        }
        return( 1 );
    }
  }
  return( 0 );
}



static void _getVoxelSizeFromIfdPtrList( _image *im, _tiffHeader *header, char *filebuffer, _tiffIfdPtrList *ifdPtrList )
{
  char * proc = "_getVoxelSizeFromIfdList";
  int e;
  int eImageDescription = -1;
  int eLsmHeader = -1;
  _tiffIfd *ifdptr = ifdPtrList->data[0];
  _LsmHeader lsmHeader;

  /* { 270, "ImageDescription" },
   * { 34412, "LSM header" },
   */
  for ( e=0; e<ifdptr->n_data; e++ ) {
    switch( *(ifdptr->data[e].tag) ) {
    default :
      break;
    case 270 :
      eImageDescription = e; break;
    case 34412 :
      eLsmHeader = e; break;
    }
  }


  if ( eLsmHeader != -1 ) {
    _initLsmHeader( &lsmHeader );
    memcpy( &lsmHeader, &(filebuffer[*(ifdptr->data[eLsmHeader].offset)]), sizeof(lsmHeader) );
    if ( header->swap ) _swapLsmHeader( &lsmHeader );
    im->vx = lsmHeader.VoxelSizeX * 1000000;
    im->vy = lsmHeader.VoxelSizeY * 1000000;
    im->vz = lsmHeader.VoxelSizeZ * 1000000;
    if ( _debug_ ) {
      fprintf( stderr, "%s: voxel sizes after ImageDescription = %lf %lf %lf from LSM\n", proc, im->vx, im->vy, im->vz );
    }
    return;
  }

  if ( eImageDescription == -1 ) return;

  if ( *(ifdptr->data[eImageDescription].length) <= 4 ) return;
  if ( _getVoxelDimensionFromLiteTiff( &(filebuffer[*(ifdptr->data[eImageDescription].offset)]),
                                      *(ifdptr->data[eImageDescription].length), im ) == 1 ) {
     if ( _debug_ ) {
       fprintf( stderr, "%s: voxel sizes after ImageDescription = %lf %lf %lf from LiteTiff\n", proc, im->vx, im->vy, im->vz );
     }
  }
  else if ( _getVoxelDimensionFromImageJ( &(filebuffer[*(ifdptr->data[eImageDescription].offset)]),
                                     *(ifdptr->data[eImageDescription].length), im ) == 1 ) {
    if ( _debug_ ) {
      fprintf( stderr, "%s: voxel sizes after ImageDescription = %lf %lf %lf from ImageJ\n", proc, im->vx, im->vy, im->vz );
    }
  }
  else if ( _getVoxelDimensionFromOmero( &(filebuffer[*(ifdptr->data[eImageDescription].offset)]),
                                      *(ifdptr->data[eImageDescription].length), im ) == 1 ) {
    if ( _debug_ ) {
      fprintf( stderr, "%s: voxel sizes after ImageDescription = %lf %lf %lf from Omero\n", proc, im->vx, im->vy, im->vz );
    }
  }
}



static int _planarToChunky( _image *im )
{
  char *proc = "_planarToChunky";
  size_t x, y, z, v;

#define _PLANARTOCHUNKY( TYPE ) {                                           \
  TYPE *tmpBuf, *imBuf;                                                     \
  imBuf = (TYPE*)im->data;                                                  \
  if ( im->wdim != sizeof(TYPE) ) {                                         \
    if ( _verbose_ || _debug_ ) {                                           \
      fprintf( stderr, "%s: bad type length\n", proc );                     \
    }                                                                       \
    return( -1 );                                                           \
  }                                                                         \
  tmpBuf = (TYPE*)malloc(im->xdim * im->ydim * im->vdim * sizeof(TYPE));    \
  if ( tmpBuf == (TYPE*)NULL ) {                                            \
    if ( _verbose_ || _debug_ )                                             \
      fprintf( stderr, "%s: allocation failed\n", proc );                   \
    return( -1 );                                                           \
  }                                                                         \
  for ( z=0; z<im->zdim; z++ ) {                                            \
    memcpy( tmpBuf, imBuf, im->xdim * im->ydim * im->vdim * sizeof(TYPE));  \
    for ( y=0; y<im->ydim; y++ )                                            \
    for ( x=0; x<im->xdim; x++ )                                            \
    for ( v=0; v<im->vdim; v++ )                                            \
      imBuf[(y*im->xdim+x)*im->vdim+v] = tmpBuf[(v*im->ydim+y)*im->xdim+x]; \
    imBuf += im->xdim * im->ydim * im->vdim;                                \
  }                                                                         \
  free( tmpBuf );                                                           \
}

  switch( im->wdim ) {
  default :
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: such word dimension (%d) not handled yet\n", proc, im->wdim);
    return( -1 );
  case 1 :
    _PLANARTOCHUNKY( unsigned char );
    break;
  case 2 :
    _PLANARTOCHUNKY( unsigned short int );
    break;
  }

  return( 1 );
}



static int _fillImageFromIfdPtrList( _image *im, _tiffHeader *header, char *filebuffer, _tiffIfdPtrList *ifdPtrList )
{
  char *proc = "_fillImageFromIfdPtrList";
  _tiffIfd *ifdptr;
  int e;
  unsigned short int *u16ptr;
  unsigned int *u32ptr;

  int eBitsPerSample = -1;
  int eColorMap = -1;
  unsigned short int BitsPerSample[3] = { 1, 1, 1 };
  int Compression = 1;
  int PhotometricInterpretation = -1;
  int PlanarConfiguration = 1;
  int ResolutionUnit = -1;
  int SamplesPerPixel = -1;
  unsigned int XResolution[2] = { 1, 1 };
  unsigned int YResolution[2] = { 1, 1 };
  int SampleFormat = 1;
  int TileLength = 0;
  int TileWidth = 0;


  /* fill image header
   * dimension
   * { 256, "ImageWidth" },
   * { 257, "ImageLength" },
   * { 258, "BitsPerSample" },
   * { 259, "Compression" },
   *   1 = No compression
   *   2 = CCITT Group 3 1-Dimensional Modified Huffman run-length encoding
   *   3 = T4-encoding: CCITT T.4 bi-level encoding
   *   4 = T6-encoding: CCITT T.6 bi-level encoding
   *   5 = LZW
   *   32773 = PackBits compression, a simple byte-oriented run-length scheme
   *   Default = 1.
   * { 262, "PhotometricInterpretation" },
   *   0 = WhiteIsZero
   *   1 = BlackIsZero
   *   2 = RGB
   *   3 = Palette color
   *   4 = Transparency Mask
   *   5 (Separated - usually CMYK).
   *   6 = YCbCr
   * { 277, "SamplesPerPixel" },
   * { 282, "XResolution" },
   * { 283, "YResolution" },
   * { 284, "PlanarConfiguration" },
   *   1 = Chunky format. The component values for each pixel are stored contiguously.
   *   2 = Planar format. The components are stored in separate “component planes.”
   * { 296, "ResolutionUnit" },
   *   1 = No absolute unit of measurement.
   *   2 = Inch.
   *   3 = Centimeter.
   *   Default is 2.
   * { 320, "ColorMap" },
   * { 322, "TileWidth" },
   * { 323, "TileLength" },
   * { 339, "SampleFormat" },
   *   1 = unsigned integer data
   *   2 = two’s complement signed integer data
   *   3 = IEEE floating point data [IEEE]
   *   4 = undefined data format
   *
   * Recall that ifd data have already been swapped (if needed), so that only data
   * given by pointer requires to be swapped (if needed)
   */
  ifdptr = ifdPtrList->data[0];

  for ( e=0; e<ifdptr->n_data; e++ ) {
    switch( *(ifdptr->data[e].tag) ) {
    default :
      break;
    case 256 :
      im->xdim = *(ifdptr->data[e].offset); break;
    case 257 :
      im->ydim = *(ifdptr->data[e].offset); break;
    case 258 :
      eBitsPerSample = e; break;
    case 259 :
      Compression = *(ifdptr->data[e].offset); break;
    case 262 :
      PhotometricInterpretation = *(ifdptr->data[e].offset); break;
    case 277 :
      SamplesPerPixel = *(ifdptr->data[e].offset); break;
    case 282 :
      u32ptr = (unsigned int*)&(filebuffer[*(ifdptr->data[e].offset)]);
      XResolution[0] = u32ptr[0];
      XResolution[1] = u32ptr[1];
      if ( header->swap ) {
        _swap4bytes( (unsigned char *)&(XResolution[0]) );
        _swap4bytes( (unsigned char *)&(XResolution[1]) );
      }
      break;
    case 283 :
      u32ptr = (unsigned int*)&(filebuffer[*(ifdptr->data[e].offset)]);
      YResolution[0] = u32ptr[0];
      YResolution[1] = u32ptr[1];
      if ( header->swap ) {
        _swap4bytes( (unsigned char *)&(YResolution[0]) );
        _swap4bytes( (unsigned char *)&(YResolution[1]) );
      }
      break;
    case 284 :
      PlanarConfiguration = *(ifdptr->data[e].offset); break;
    case 296 :
      ResolutionUnit = *(ifdptr->data[e].offset); break;
    case 320 :
      eColorMap = e; break;
    case 322 :
      TileWidth = *(ifdptr->data[e].offset); break;
    case 323 :
      TileLength = *(ifdptr->data[e].offset); break;
    case 339 :
      SampleFormat = *(ifdptr->data[e].offset); break;
    }
  }
  im->zdim = ifdPtrList->n_data;

  if ( eBitsPerSample == 1 ) {
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: no BitsPerSample IFD\n", proc );
    return( -1 );
  }

  /* { 296, "ResolutionUnit" },
   *   1 = No absolute unit of measurement.
   *   2 = Inch.
   *   3 = Centimeter.
   *   Default is 2.
   */
  switch( ResolutionUnit ) {
  default :
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: such ResolutionUnit (%d) not handled yet\n", proc, ResolutionUnit );
    return( -1 );
  case -1 :
    /* weird, no resolution unit
     */
    if ( XResolution[0] == 1 && XResolution[1] == 1 ) {
        im->vx = 1.0;
    }
    else {
      im->vx = (2.54 * (double)XResolution[1] ) / ( 10 * (double)XResolution[0] );
    }
    if ( YResolution[0] == 1 && YResolution[1] == 1 ) {
        im->vy = 1.0;
    }
    else {
      im->vy = (2.54 * (double)YResolution[1] ) / ( 10 * (double)YResolution[0] );
    }
    break;
  case 1 :
    im->vx = (double)XResolution[1] / (double)XResolution[0];
    im->vy = (double)YResolution[1] / (double)YResolution[0];
    break;
  case 2 :
    im->vx = (2.54 * (double)XResolution[1] ) / ( 10 * (double)XResolution[0] );
    im->vy = (2.54 * (double)YResolution[1] ) / ( 10 * (double)YResolution[0] );
    break;
  case 3 :
    im->vx = (double)XResolution[1] / ( 10 * (double)XResolution[0] );
    im->vy = (double)YResolution[1] / ( 10 * (double)YResolution[0] );
    break;
  }



  /* here, we can try to read and/or correct voxel size
   * from ImageDescription
   */
  _getVoxelSizeFromIfdPtrList( im, header, filebuffer, ifdPtrList );


  /* PhotometricInterpretation
   * 0 = WhiteIsZero
   * 1 = BlackIsZero
   * 2 = RGB
   * 3 = Palette color
   * 4 = Transparency Mask
   * 5 (Separated - usually CMYK).
   * 6 = YCbCr
   */
  switch( PhotometricInterpretation ) {
  default :
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: such photometric interpretation (%d) not handled yet\n",
               proc, PhotometricInterpretation );
    return( -1 );
  case 0 :
  case 1 :
    im->vdim = ( SamplesPerPixel == -1 ) ? 1 : SamplesPerPixel;
    /* { 258, "BitsPerSample" },
     */
    if ( eBitsPerSample != -1 ) {
      if ( *(ifdptr->data[eBitsPerSample].length) != im->vdim ) {
        if ( _verbose_ || _debug_ )
          fprintf( stderr, "%s: numbers (%d) of bits per sample is different from %d for PhotometricInterpretation=%d\n",
                     proc, *(ifdptr->data[eBitsPerSample].length), im->vdim, PhotometricInterpretation);
        return( -1 );
      }
      if ( im->vdim == 1 ) {
        BitsPerSample[0] = *(ifdptr->data[eBitsPerSample].offset);
      }
      else if ( im->vdim == 2 ) {
        u16ptr = (unsigned short int*)(ifdptr->data[eBitsPerSample].offset);
        BitsPerSample[0] = u16ptr[0];
        BitsPerSample[1] = u16ptr[1];
        if ( BitsPerSample[0] != BitsPerSample[1] ) {
          if ( _verbose_ || _debug_ )
            fprintf( stderr, "%s: numbers of bits per sample are different (%d,%d)\n",
                       proc, BitsPerSample[0], BitsPerSample[1] );
          return( -1 );
        }
      }
      else {
        if ( _verbose_ || _debug_ )
          fprintf( stderr, "%s: such vectorial dimension (%d) not handled yet\n", proc, im->vdim );
        return( -1 );
      }
    }
    break;
    /* end of case PhotometricInterpretation = 0, 1
     */
  case 2 :
    if ( SamplesPerPixel != 3 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: weird number (%d) of samples per pixel for PhotometricInterpretation=%d\n",
                   proc, SamplesPerPixel, PhotometricInterpretation);
      return( -1 );
    }
    im->vdim = 3;
    /* { 258, "BitsPerSample" },
     */
    if ( eBitsPerSample != -1 ) {
      if ( *(ifdptr->data[eBitsPerSample].length) != 3 ) {
        if ( _verbose_ || _debug_ )
          fprintf( stderr, "%s: weird numbers (%d) of bits per sample for PhotometricInterpretation=%d\n",
                     proc, *(ifdptr->data[eBitsPerSample].length), PhotometricInterpretation);
        return( -1 );
      }
      u16ptr = (unsigned short int*)&(filebuffer[*(ifdptr->data[eBitsPerSample].offset)]);
      BitsPerSample[0] = u16ptr[0];
      BitsPerSample[1] = u16ptr[1];
      BitsPerSample[2] = u16ptr[2];
      if ( header->swap ) {
        _swap2bytes( (unsigned char *)&(BitsPerSample[0]) );
        _swap2bytes( (unsigned char *)&(BitsPerSample[1]) );
        _swap2bytes( (unsigned char *)&(BitsPerSample[2]) );
      }
      if ( BitsPerSample[0] != BitsPerSample[1] || BitsPerSample[0] != BitsPerSample[2] ) {
        if ( _verbose_ || _debug_ )
          fprintf( stderr, "%s: numbers of bits per sample are different (%d,%d,%d)\n",
                     proc, BitsPerSample[0], BitsPerSample[1], BitsPerSample[2] );
        return( -1 );
      }
    }
    break;
  case 3 :
    /* Palette-color image
     * Here the "BitsPerSample" IFD refers to the index for the ColorMap
     * Translate the ColorMap fieldtype into BitsPerSample
     */

    if ( SamplesPerPixel != 1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: weird number (%d) of samples per pixel for PhotometricInterpretation=%d\n",
                   proc, SamplesPerPixel, PhotometricInterpretation);
      return( -1 );
    }
    if ( eColorMap == -1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: no ColorMap IFD for PhotometricInterpretation=%d\n",
                   proc, PhotometricInterpretation);
      return( -1 );
    }
    /* allows to fill image with color index
     */
    if ( _consider_palette_as_grey_level_ ) {
      im->vdim = 1;
      if ( eBitsPerSample != -1 ) {
        if ( *(ifdptr->data[eBitsPerSample].length) != 1 ) {
          if ( _verbose_ || _debug_ )
            fprintf( stderr, "%s: weird numbers (%d) of bits per sample for PhotometricInterpretation=%d\n",
                       proc, *(ifdptr->data[eBitsPerSample].length), PhotometricInterpretation);
          return( -1 );
        }
        u16ptr = (unsigned short int*)(ifdptr->data[eBitsPerSample].offset);
        BitsPerSample[0] = u16ptr[0];
      }
      else {
        if ( _verbose_ || _debug_ )
          fprintf( stderr, "%s: no BitsPerSample IFD for PhotometricInterpretation=%d\n",
                     proc, PhotometricInterpretation);
        return( -1 );
      }
    }
    else {
      im->vdim = 3;
      switch( *(ifdptr->data[eColorMap].fieldType) ) {
      default :
        if ( _verbose_ || _debug_ )
          fprintf( stderr, "%s: such color map field type (%d) not handled yet\n", proc, *(ifdptr->data[eColorMap].fieldType) );
        return( -1 );
      case 3 :
        /* USHORT type
         */
        im->wdim = 2;
        im->wordKind = WK_FIXED;
        im->sign = SGN_UNSIGNED;
        break;
      }
    }
    break;
  }

  /* { 339, "SampleFormat" },
   *   1 = unsigned integer data
   *   2 = two’s complement signed integer data
   *   3 = IEEE floating point data [IEEE]
   *   4 = undefined data format
   *
   * Do nothing for Palette Color Image
   */
  if ( PhotometricInterpretation != 3 || _consider_palette_as_grey_level_ ) {
    switch( SampleFormat ) {
    default :
    case 1 :
    case 4 :
      switch( BitsPerSample[0] ) {
      default :
      case 1 :
      case 4 :
        if ( _verbose_ || _debug_ )
          fprintf( stderr, "%s: such number of bits per sample (%d) not handled yet\n",
                    proc, BitsPerSample[0] );
        return( -1 );
      case 8 :
        /* unsigned char
         */
        im->wdim = 1;
        im->wordKind = WK_FIXED;
        im->sign = SGN_UNSIGNED;
        break;
      case 16 :
        /* unsigned short int
         */
        im->wdim = 2;
        im->wordKind = WK_FIXED;
        im->sign = SGN_UNSIGNED;
        break;
      }
      break;
    case 2 :
      switch( BitsPerSample[0] ) {
      default :
        if ( _verbose_ || _debug_ )
          fprintf( stderr, "%s: such number of bits per sample (%d) not handled yet\n",
                    proc, BitsPerSample[0] );
        return( -1 );
      case 8 :
        /* signed char
         */
        im->wdim = 1;
        im->wordKind = WK_FIXED;
        im->sign = SGN_SIGNED;
        break;
      case 16 :
        /* signed short int
         */
        im->wdim = 2;
        im->wordKind = WK_FIXED;
        im->sign = SGN_SIGNED;
        break;
      }
      break;
    case 3 :
      switch( BitsPerSample[0] ) {
      default :
        if ( _verbose_ || _debug_ )
          fprintf( stderr, "%s: such number of bits per sample (%d) not handled yet\n",
                    proc, BitsPerSample[0] );
        return( -1 );
      case 32 :
        /* float
         */
        im->wdim = 1;
        im->wordKind = WK_FLOAT;
        im->sign = SGN_SIGNED;
        break;
      }
      break;
    }
  }


  /* we are ready to fill the image
   */
  if ( TileWidth > 0 || TileLength > 0 ) {
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: tiles (%d x %d) reading has not been implemented yet\n",
                proc, TileWidth, TileLength );
    return( -1 );
  }


  /* read data
   * { 273, "StripOffsets" },
   *   SHORT or LONG.
   *   For each strip, the byte offset of that strip.
   *   N = StripsPerImage for PlanarConfiguration equal to 1.
   *     = SamplesPerPixel * StripsPerImage for PlanarConfiguration equal to 2
   * { 278, "RowsPerStrip" },
   *   SHORT or LONG.
   *   The number of rows in each strip (except possibly the last strip.)
   *   length of a row is ImageWidth
   *   StripsPerImage = floor ((ImageLength + RowsPerStrip - 1) / RowsPerStrip).
   *   The default is 2**32 - 1, which is effectively infinity.
   *   That is, the entire image is one strip.
   *   Recommendation: Set RowsPerStrip such that the size of each strip is about 8K bytes.
   * { 279, "StripByteCounts" },
   *   SHORT or LONG.
   *   For each strip, the number of bytes in that strip after any compression.
   *   N = StripsPerImage for PlanarConfiguration equal to 1.
   *     = SamplesPerPixel * StripsPerImage for PlanarConfiguration equal to 2
   */

  im->data = ImageIO_alloc( im->xdim * im->ydim * im->zdim * (size_t)im->vdim * (size_t)im->wdim );
  if ( im->data == (void*)NULL ) {
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: image buffer allocation failed\n", proc );
    return( -1 );
  }


  /* { 259, "Compression" },
   *   1 = No compression
   *   2 = CCITT Group 3 1-Dimensional Modified Huffman run-length encoding
   *   3 = T4-encoding: CCITT T.4 bi-level encoding
   *   4 = T6-encoding: CCITT T.6 bi-level encoding
   *   5 = LZW
   *   32773 = PackBits compression, a simple byte-oriented run-length scheme
   *   Default = 1.
   */
  switch( Compression ) {
  case 2 :
  case 3 :
  case 4 :
  case 5 :
  default :
     ImageIO_free(im->data);
     im->data = (void*)NULL;
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: such Compression (%d) not handled yet\n", proc, Compression );
    return( -1 );
  case 1 :
    switch( PhotometricInterpretation ) {
    default :
      ImageIO_free(im->data);
      im->data = (void*)NULL;
      fprintf( stderr, "%s: such photometric interpretation (%d) not handled yet for Compression=%d\n",
               proc, PhotometricInterpretation, Compression );
      return( -1 );
    case 0 :
    case 1 :
    case 2 :
      if ( _readUncompressedDataFromIfdPtrList( filebuffer, ifdPtrList, header, im ) != 1 ) {
        ImageIO_free(im->data);
        im->data = (void*)NULL;
        if ( _verbose_ || _debug_ )
          fprintf( stderr, "%s: unable to read data (Compression=%d, No compression, PhotometricInterpretation=%d)\n", proc, Compression, PhotometricInterpretation );
        return( -1 );
      }
      break;
    case 3 :
      if (_consider_palette_as_grey_level_ ) {
        if ( _readUncompressedDataFromIfdPtrList( filebuffer, ifdPtrList, header, im ) != 1 ) {
          ImageIO_free(im->data);
          im->data = (void*)NULL;
          if ( _verbose_ || _debug_ )
            fprintf( stderr, "%s: unable to read data (Compression=%d, No compression, PhotometricInterpretation=%d)\n", proc, Compression, PhotometricInterpretation );
          return( -1 );
        }
        break;
      }
      if ( _readUncompressedPaletteFromIfdPtrList( filebuffer, ifdPtrList, header, im ) != 1 ) {
        ImageIO_free(im->data);
        im->data = (void*)NULL;
        if ( _verbose_ || _debug_ )
          fprintf( stderr, "%s: unable to read data (Compression=%d, No compression, PhotometricInterpretation=%d)\n", proc, Compression, PhotometricInterpretation );
        return( -1 );
      }
      break;
    }
    break;
  case 32773 :
    switch( PhotometricInterpretation ) {
    default :
      ImageIO_free(im->data);
      im->data = (void*)NULL;
      fprintf( stderr, "%s: such photometric interpretation (%d) not handled yet for Compression=%d\n",
               proc, PhotometricInterpretation, Compression );
      return( -1 );
    case 0 :
    case 1 :
    case 2 :
      if ( _readPackBitsDataFromIfdPtrList( filebuffer, ifdPtrList, header, im ) != 1 ) {
          ImageIO_free(im->data);
          im->data = (void*)NULL;
         if ( _verbose_ || _debug_ )
           fprintf( stderr, "%s: unable to read data (Compression=%d, PhotometricInterpretation=%d)\n", proc, Compression, PhotometricInterpretation );
         return( -1 );
       }
      break;
    }
    break;
  }


  /* swapping data
   */
  if ( header->swap ) {
    switch ( im->wdim ) {
    default :
       ImageIO_free(im->data);
       im->data = (void*)NULL;
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: such word dim not handled\n", proc );
      return( -1 );
    case 1 :
        break;
    case 2 :
      {
        size_t l = im->xdim * im->ydim * im->zdim * (size_t)im->vdim;
        unsigned short int *theBuf = (unsigned short int *)im->data;
        unsigned short int *resBuf = (unsigned short int *)im->data;
        unsigned short int usi;
        for ( ; l>0; l--, theBuf++, resBuf++ ) {
          usi = *theBuf;
          *resBuf = ((usi >> 8) & 0xff) | (usi << 8);
        }
      }
      break;
    case 4 :
    {
      size_t l = im->xdim * im->ydim * im->zdim * (size_t)im->vdim;
      unsigned int *theBuf = (unsigned int *)im->data;
      unsigned int *resBuf = (unsigned int *)im->data;
      unsigned int ui;
      for ( ; l>0; l--, theBuf++, resBuf++ ) {
        ui = *theBuf;
        *resBuf = (ui << 24) | ((ui & 0xff00) << 8) | ((ui >> 8) & 0xff00) | ((ui >> 24) & 0xff);
      }
    }
    break;
    }
  }

  /* PhotometricInterpretation
   * 0 = WhiteIsZero
   * 1 = BlackIsZero
   * 2 = RGB
   * 3 = Palette color
   * 4 = Transparency Mask
   * 5 (Separated - usually CMYK).
   * 6 = YCbCr
   */
  switch( PhotometricInterpretation ) {
  default :
    break;
  case 0 :
    {
      size_t l = im->xdim * im->ydim * im->zdim * (size_t)im->vdim * (size_t)im->wdim;
      unsigned char *theBuf = (unsigned char *)im->data;
      unsigned char *resBuf = (unsigned char *)im->data;
      unsigned char uc;
      for ( ; l>0; l--, theBuf++, resBuf++ ) {
        uc = *theBuf;
        *resBuf = uc ^ 0xff;
      }
    }
    break;
  }

  /* { 284, "PlanarConfiguration" },
   *   1 = Chunky format. The component values for each pixel are stored contiguously.
   *   2 = Planar format. The components are stored in separate “component planes.”
   */

  switch( PlanarConfiguration ) {
  default :
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: unexpected value for PlanarConfiguration (%d)\n", proc, PlanarConfiguration );
    break;
  case 1 :
    break;
  case 2 :
    if ( im->vdim == 1 )
      break;
    if ( _debug_ ) {
      fprintf( stderr, "... %s: conversion from planar format to chunky format\n", proc );
    }
    if ( _planarToChunky( im ) != 1 ) {
      ImageIO_free(im->data);
      im->data = (void*)NULL;
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: planar to chunky conversion failed\n", proc );
      return( -1 );
    }
    break;
  }

  return( 1 );
}



int readLiteTiffImage( const char *name, _image *im )
{
  char *proc = "readLiteTiffImage";
  size_t filesize = 0;
  char *filebuffer = (char*)NULL;
  size_t nread;

  _tiffHeader header;
  unsigned int offset;
  _tiffIfd ifd;
  int nifd = 0;
  _tiffIfdList ifdList;

  _tiffIfdPtrListList ifdTrimmedList;
  int selectSubList = 0;

  /* reading the whole file
   * this should be modified to handle gzip'ed files
   */
  filesize = getFileSize( name );
  if ( filesize <= 0 ) {
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: size of file '%s' negative or null ?!", proc, name );
    return( -1 );
  }
  filebuffer = (char*)malloc( filesize );
  if ( filebuffer == (char*)NULL ) {
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: file buffer allocation failed\n", proc );
    return( -1 );
  }
  nread = ImageIO_read( im, filebuffer, filesize);
  if ( nread != filesize ) {
      free( filebuffer );
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: image buffer reading failed\n", proc );
      return( -1 );
  }



  /* the whole buffer has been stored in memory
   * parse it
   * 1. header
   * 2. IFDs
   */
  _initTiffHeader( &header );
  memcpy( header.data, filebuffer, 8 );
  if ( _parseHeader( &header ) != 1 ) {
    free( filebuffer );
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: error when parsing TIFF header\n", proc );
    return( -1 );
  }

  if ( _debug_ )
    fprintf( stderr, "%s: TIFF header recognized, swap = %d, filesize = %ld\n",
             proc, header.swap, filesize );

  /* reading IFDs
   */
  offset = *(header.offset);
  _initTiffIfdList( &ifdList );

  do {
    _initTiffIfd( &ifd );
    if ( _readIfd( &(filebuffer[offset]), &ifd, &header ) != 1 ) {
      _freeTiffIfdList( &ifdList );
      free( filebuffer );
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: error when reading ifd #%d\n", proc, nifd );
      return( -1 );
    }

    if ( _addTiffIfdToList( &ifdList, &ifd ) != 1 ) {
      _freeTiffIfdList( &ifdList );
      free( filebuffer );
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: unable to add ifd #%d to list\n", proc, nifd );
      return( -1 );
    }

    nifd ++;
    offset = ifd.offset;

  } while ( offset > 0 );

  if ( _debug_ >= 2 ) {
    _fprintfTiffHeader( stderr, &header );
    _fprintTiffIfdList( stderr, filebuffer, &ifdList, header.swap );
  }

  /* IFDs have been read
   * do some tests
   */

  if ( ifdList.n_data == 0 ) {
    _freeTiffIfdList( &ifdList );
    free( filebuffer );
    if ( _verbose_ || _debug_ ) {
      fprintf( stderr, "%s: no IFD was found in '%s'?!\n", proc, name );
    }
    return( -1 );
  }

  /* sort ifds:
   * there are perhaps several images
   */

  _initTiffIfdPtrListList( &ifdTrimmedList );
  if ( _sortTiffIfdList( &ifdTrimmedList, &ifdList ) != 1 ) {
    _freeTiffIfdList( &ifdList );
    free( filebuffer );
    if ( _verbose_ || _debug_ ) {
      fprintf( stderr, "%s: failed to sort ifds\n", proc );
    }
    return( -1 );
  }

  if ( ifdTrimmedList.n_data <= 0 ) {
    _freeTiffIfdPtrListList( &ifdTrimmedList );
    _freeTiffIfdList( &ifdList );
    free( filebuffer );
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: no image found after ifd sorting ?!\n", proc );
    return( -1 );
  }

  selectSubList = _selectIfdSubList( &ifdTrimmedList );

  if ( _verbose_ >= 2 || _debug_ ) {
    if ( ifdTrimmedList.n_data >= 2 ) {
      fprintf( stderr, "%s: select ifd sub list is #%d/%d\n", proc, selectSubList, ifdTrimmedList.n_data );
    }
    if ( _debug_ >= 1 )
      _fprintTiffIfdPtrListList( stderr, filebuffer, &ifdTrimmedList, header.swap );
  }

 if ( _fillImageFromIfdPtrList( im, &header, filebuffer, &(ifdTrimmedList.data[selectSubList]) ) != 1 ) {
   _freeTiffIfdPtrListList( &ifdTrimmedList );
   _freeTiffIfdList( &ifdList );
   free( filebuffer );
   if ( _verbose_ || _debug_ )
     fprintf( stderr, "%s: unable to fill image from selected ifds\n", proc );
   return( -1 );
 }


  _freeTiffIfdPtrListList( &ifdTrimmedList );
  _freeTiffIfdList( &ifdList );
  free( filebuffer );



  return( 1 );
}






/************************************************************
 *
 *
 *
 ************************************************************/


static void _computeResolution( unsigned int *rv,
                                double size )
{
  unsigned int best[2], n, d;
  double r, bestresidual;
  double epsilon = 1e-6;

  d = 1;
  n = (unsigned int)( size * d + 0.5 );
  bestresidual = size > (double)n/(double)d ? size - (double)n/(double)d : (double)n/(double)d - size;
  if ( bestresidual < epsilon ) {
    rv[0] = d;
    rv[1] = n;
    return;
  }
  best[0] = d;
  best[1] = n;

  for ( d=2; d<1000000; d++ ) {
    n = (unsigned int)( size * d + 0.5 );
    r = size > (double)n/(double)d ? size - (double)n/(double)d : (double)n/(double)d - size;
    if ( r < epsilon ) {
      rv[0] = d;
      rv[1] = n;
      return;
    }
    if ( r < bestresidual ) {
      best[0] = d;
      best[1] = n;
      bestresidual = r;
    }
  }
  rv[0] = best[0];
  rv[1] = best[1];
}



static int _fillTiffHeader( _tiffHeader *h )
{
  char *proc = "_fillTiffHeader";

  switch( _getEndianness() ) {
  default :
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: unknown endianess\n", proc );
    return( -1 );
  case END_LITTLE :
    h->data[0] = h->data[1] = 73;
    break;
  case END_BIG :
    h->data[0] = h->data[1] = 77;
    break;
  }

  *(h->version) = 42;
  /* ifds will be written
   * just after the header
   */
  *(h->offset) = 8;

  h->swap = 0;

  return( 1 );
}



static int _setVoxelDimension( _tiffIfdEntry *ImageDescription, _image *im )
{
  char *proc = "_setVoxelDimension";
  char str[256];
  int i, l;

  for ( i=0; i<256; i++ ) str[i] = '\0';
  sprintf( str, "LiteTiff\nVX=%f\nVY=%f\nVZ=%f\n", im->vx, im->vy, im->vz );
  l = strlen( str );

  _initTiffIfdEntry( ImageDescription );

  ImageDescription->values = (char*)malloc( l+1 );
  if ( ImageDescription->values == (char*)NULL ) {
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }
  memcpy( ImageDescription->values, str, l+1 );

  *(ImageDescription->tag)       = 270;
  *(ImageDescription->fieldType) = 2;
  *(ImageDescription->length)    = l+1;
  ImageDescription->valuesLength = l+1;
  return( 1 );
}



/* Required Fields for Grayscale Images
 *
 * ImageWidth                256  SHORT or LONG
 * ImageLength               257  SHORT or LONG
 * BitsPerSample             258  SHORT 4 or 8
 * Compression               259  SHORT 1 or 32773
 * PhotometricInterpretation 262  SHORT 0 or 1
 * StripOffsets              273  SHORT or LONG
 * RowsPerStrip              278  SHORT or LONG
 * StripByteCounts           279  LONG or SHORT
 * XResolution               282  RATIONAL
 * YResolution               283  RATIONAL
 * ResolutionUnit            296  SHORT 1 or 2 or 3

 * Required Fields for RGB Images
 * ImageWidth                256  SHORT or LONG
 * ImageLength               257  SHORT or LONG
 * BitsPerSample             258  SHORT 8,8,8
 * Compression               259  SHORT 1 or 32773
 * PhotometricInterpretation 262  SHORT 2
 * StripOffsets              273  SHORT or LONG
 * SamplesPerPixel           277  SHORT 3 or more
 * RowsPerStrip              278  SHORT or LONG
 * StripByteCounts           279  LONG or SHORT
 * XResolution               282  RATIONAL
 * YResolution               283  RATIONAL
 * ResolutionUnit            296  SHORT 1, 2 or 3
 *
 * 1 = BYTE	An 8-bit unsigned integer.
 * 2 = ASCII	8-bit bytes that store ASCII codes; the last byte must be null.
 * 3 = SHORT	A 16-bit (2-byte) unsigned integer.
 * 4 = LONG	A 32-bit (4-byte) unsigned integer.
 * 5 = RATIONAL	Two LONG_s:  the first represents the numerator of a fraction, the second the denominator.
 *
 * RowsPerStrip. SHORT or LONG. Readers must be able to handle any value
 * between 1 and 2**32-1. However, some readers may try to read an entire strip
 * into memory at one time. If the entire image is one strip, the application may run
 * out of memory. Recommendation: Set RowsPerStrip such that the size of each
 * strip is about 8K bytes. Do this even for uncompressed data because it is easy for
 * a writer and makes things simpler for readers. Note that extremely wide highresolution
 * images may have rows larger than 8K bytes; in this case, RowsPerStrip
 * should be 1, and the strip will be larger than 8K.
*/



int _addGeometryEntries( _tiffIfdList *list, _image *im )
{
  char *proc = "_addGeometryEntries";
  int i;
  unsigned int s;
  unsigned int v;
  _tiffIfdEntry ImageWidth;
  _tiffIfdEntry ImageLength;
  _tiffIfdEntry BitsPerSample;
  _tiffIfdEntry Compression;
  _tiffIfdEntry PhotometricInterpretation;
  _tiffIfdEntry ImageDescription;
  _tiffIfdEntry StripOffsets;
  _tiffIfdEntry SamplesPerPixel;
  _tiffIfdEntry RowsPerStrip;
  _tiffIfdEntry StripByteCounts;
  _tiffIfdEntry XResolution;
  _tiffIfdEntry YResolution;
  _tiffIfdEntry PlanarConfiguration;
  _tiffIfdEntry ResolutionUnit;
  _tiffIfdEntry SampleFormat;
  unsigned short int *BitsPerSampleValues;
  unsigned short int *SampleFormatValues;
  unsigned int XRationalValue[2];
  unsigned int YRationalValue[2];
  unsigned int RowsPerStripValues;
  unsigned int nStrip;
  unsigned int *StripOffsetsValues;
  unsigned int *StripByteCountsValues;



  _computeResolution( XRationalValue, im->vx );
  _computeResolution( YRationalValue, im->vy );

  _initTiffIfdEntry( &ImageWidth );
  *(ImageWidth.tag)       = 256;
  *(ImageWidth.fieldType) = 3;
  *(ImageWidth.length)    = 1;
  *(ImageWidth.offset)    = im->xdim;

  _initTiffIfdEntry( &ImageLength );
  *(ImageLength.tag)       = 257;
  *(ImageLength.fieldType) = 3;
  *(ImageLength.length)    = 1;
  *(ImageLength.offset)    = im->ydim;

  _initTiffIfdEntry( &BitsPerSample );
  *(BitsPerSample.tag)       = 258;
  *(BitsPerSample.fieldType) = 3;
  *(BitsPerSample.length)    = im->vdim;

  _initTiffIfdEntry( &Compression );
  *(Compression.tag)       = 259;
  *(Compression.fieldType) = 3;
  *(Compression.length)    = 1;
  *(Compression.offset)    = 1;

  /* RGB = unsigned char/short int with vdim = 3
   */
  _initTiffIfdEntry( &PhotometricInterpretation );
  *(PhotometricInterpretation.tag)       = 262;
  *(PhotometricInterpretation.fieldType) = 3;
  *(PhotometricInterpretation.length)    = 1;
  if ( im->vdim == 3 && (im->wdim == 1 || im->wdim == 2)
       && im->wordKind == WK_FIXED && im->sign == SGN_UNSIGNED )
    *(PhotometricInterpretation.offset)    = 2;
  else
    *(PhotometricInterpretation.offset)    = 1;

  RowsPerStripValues = (int)( 0.5 + (double)8192 / (double)(im->vdim * im->wdim * im->xdim) );
  if ( im->vdim * im->wdim * im->xdim * RowsPerStripValues < 8192 )
    RowsPerStripValues ++;
  nStrip = 0;
  while ( nStrip*RowsPerStripValues < im->ydim ) nStrip++;

  if ( _setVoxelDimension( &ImageDescription, im ) != 1 ) {
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: unable to fill entry 'ImageDescription'\n", proc );
    return( -1 );
  }

  _initTiffIfdEntry( &StripOffsets );
  *(StripOffsets.tag)       = 273;
  *(StripOffsets.fieldType) = 4;
  *(StripOffsets.length)    = nStrip;

  _initTiffIfdEntry( &SamplesPerPixel );
  *(SamplesPerPixel.tag)       = 277;
  *(SamplesPerPixel.fieldType) = 3;
  *(SamplesPerPixel.length)    = 1;
  *(SamplesPerPixel.offset)    = im->vdim;

  _initTiffIfdEntry( &RowsPerStrip );
  *(RowsPerStrip.tag)       = 278;
  *(RowsPerStrip.fieldType) = 3;
  *(RowsPerStrip.length)    = 1;
  *(RowsPerStrip.offset)    = RowsPerStripValues;

  _initTiffIfdEntry( &StripByteCounts );
  *(StripByteCounts.tag)       = 279;
  *(StripByteCounts.fieldType) = 4;
  *(StripByteCounts.length)    = nStrip;

  _initTiffIfdEntry( &XResolution );
  *(XResolution.tag)       = 282;
  *(XResolution.fieldType) = 5;
  *(XResolution.length)    = 1;
  XResolution.valuesLength = 8;

  _initTiffIfdEntry( &YResolution );
  *(YResolution.tag)       = 283;
  *(YResolution.fieldType) = 5;
  *(YResolution.length)    = 1;
  YResolution.valuesLength = 8;

  _initTiffIfdEntry( &PlanarConfiguration );
  *(PlanarConfiguration.tag)       = 284;
  *(PlanarConfiguration.fieldType) = 3;
  *(PlanarConfiguration.length)    = 1;
  *(PlanarConfiguration.offset)    = 1;

  _initTiffIfdEntry( &ResolutionUnit );
  *(ResolutionUnit.tag)       = 296;
  *(ResolutionUnit.fieldType) = 3;
  *(ResolutionUnit.length)    = 1;
  *(ResolutionUnit.offset)    = 1;

  _initTiffIfdEntry( &SampleFormat );
  *(SampleFormat.tag)       = 339;
  *(SampleFormat.fieldType) = 3;
  *(SampleFormat.length)    = im->vdim;


  /* fill IFDs
   */
  for ( i=0; i<list->n_data; i++ ) {

    if ( _addTiffIfdEntryToIfd( &(list->data[i]), &ImageWidth ) != 1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: unable to add entry 'ImageWidth' to ifd #%d\n", proc, i );
      return( -1 );
    }

    if ( _addTiffIfdEntryToIfd( &(list->data[i]), &ImageLength ) != 1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: unable to add entry 'ImageLength' to ifd #%d\n", proc, i );
      return( -1 );
    }

    if ( im->vdim <= 2 ) {
      BitsPerSampleValues = (unsigned short int*)BitsPerSample.offset;
      for ( v=0; v<im->vdim; v++ ) BitsPerSampleValues[v] = im->wdim * 8;
    }
    else {
      BitsPerSample.values = (char*)malloc( im->vdim * sizeof(unsigned short int) );
      if ( BitsPerSample.values == (char*)NULL ) {
        if ( _verbose_ || _debug_ )
          fprintf( stderr, "%s: allocation failed for 'BitsPerSample' for ifd #%d\n", proc, i );
        return( -1 );
      }
      BitsPerSampleValues = (unsigned short int*)BitsPerSample.values;
      BitsPerSample.valuesLength = im->vdim * sizeof(unsigned short int);
      for ( v=0; v<im->vdim; v++ ) BitsPerSampleValues[v] = im->wdim * 8;
    }
    if ( _addTiffIfdEntryToIfd( &(list->data[i]), &BitsPerSample ) != 1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: unable to add entry 'BitsPerSample' to ifd #%d\n", proc, i );
      return( -1 );
    }
    BitsPerSample.values = (char*)NULL;
    BitsPerSample.valuesLength = 0;

    if ( _addTiffIfdEntryToIfd( &(list->data[i]), &Compression ) != 1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: unable to add entry 'Compression' to ifd #%d\n", proc, i );
      return( -1 );
    }

    if ( _addTiffIfdEntryToIfd( &(list->data[i]), &PhotometricInterpretation ) != 1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: unable to add entry 'PhotometricInterpretation' to ifd #%d\n", proc, i );
      return( -1 );
    }

    if ( i == 0 ) {
      if ( _addTiffIfdEntryToIfd( &(list->data[i]), &ImageDescription ) != 1 ) {
        if ( _verbose_ || _debug_ )
          fprintf( stderr, "%s: unable to add entry 'ImageDescription' to ifd #%d\n", proc, i );
        return( -1 );
      }
    }

    if ( *(StripOffsets.length) > 1 ) {
      StripOffsets.values = (char*)malloc( *(StripOffsets.length) * sizeof (unsigned int) );
      if ( StripOffsets.values == (char*)NULL ) {
        if ( _verbose_ || _debug_ )
          fprintf( stderr, "%s: allocation failed for 'StripOffsets' for ifd #%d\n", proc, i );
        return( -1 );
      }
      StripOffsetsValues = (unsigned int*)StripOffsets.values;
      for ( s=0; s<nStrip; s++ ) StripOffsetsValues[s] = 0;
      StripOffsets.valuesLength = *(StripOffsets.length) * sizeof (unsigned int) ;
    }
    if ( _addTiffIfdEntryToIfd( &(list->data[i]), &StripOffsets ) != 1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: unable to add entry 'StripOffsets' to ifd #%d\n", proc, i );
      return( -1 );
    }
    StripOffsets.values = (char*)NULL;
    StripOffsets.valuesLength = 0;

    if ( _addTiffIfdEntryToIfd( &(list->data[i]), &SamplesPerPixel ) != 1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: unable to add entry 'SamplesPerPixel' to ifd #%d\n", proc, i );
      return( -1 );
    }

    if ( _addTiffIfdEntryToIfd( &(list->data[i]), &RowsPerStrip ) != 1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: unable to add entry 'RowsPerStrip' to ifd #%d\n", proc, i );
      return( -1 );
    }

    if ( *(StripByteCounts.length) > 1 ) {
      StripByteCounts.values = (char*)malloc( *(StripByteCounts.length) * sizeof (unsigned int) );
      if ( StripByteCounts.values == (char*)NULL ) {
        if ( _verbose_ || _debug_ )
          fprintf( stderr, "%s: allocation failed for 'StripByteCounts' for ifd #%d\n", proc, i );
        return( -1 );
      }
      StripByteCountsValues = (unsigned int*)StripByteCounts.values;
      for ( s=0; s<nStrip; s++ ) StripByteCountsValues[s] = RowsPerStripValues;
      if ( nStrip * RowsPerStripValues > im->ydim )
        StripByteCountsValues[nStrip-1] -= nStrip * RowsPerStripValues - im->ydim;
      for ( s=0; s<nStrip; s++ ) StripByteCountsValues[s] *= im->vdim * im->wdim * im->xdim;
      StripByteCounts.valuesLength = *(StripByteCounts.length) * sizeof (unsigned int) ;
    }
    else {
      *(StripByteCounts.offset) = (im->vdim * im->wdim * im->xdim * im->ydim);
    }
    if ( _addTiffIfdEntryToIfd( &(list->data[i]), &StripByteCounts ) != 1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: unable to add entry 'StripByteCounts' to ifd #%d\n", proc, i );
      return( -1 );
    }
    StripByteCounts.values = (char*)NULL;
    StripByteCounts.valuesLength = 0;

    XResolution.values = (char*)malloc( 2 * 4 );
    if ( XResolution.values == (char*)NULL ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: allocation failed for 'XResolution' for ifd #%d\n", proc, i );
    }
    memcpy( XResolution.values, XRationalValue, 8 );
    if ( _addTiffIfdEntryToIfd( &(list->data[i]), &XResolution ) != 1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: unable to add entry 'XResolution' to ifd #%d\n", proc, i );
      return( -1 );
    }
    XResolution.values = (char*)NULL;

    YResolution.values = (char*)malloc( 2 * 4 );
    if ( YResolution.values == (char*)NULL ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: allocation failed for 'YResolution' for ifd #%d\n", proc, i );
    }
    memcpy( YResolution.values, YRationalValue, 8 );
    if ( _addTiffIfdEntryToIfd( &(list->data[i]), &YResolution ) != 1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: unable to add entry 'YResolution' to ifd #%d\n", proc, i );
      return( -1 );
    }
    YResolution.values = (char*)NULL;


    if ( _addTiffIfdEntryToIfd( &(list->data[i]), &ResolutionUnit ) != 1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: unable to add entry 'ResolutionUnit' to ifd #%d\n", proc, i );
      return( -1 );
    }

    if ( _addTiffIfdEntryToIfd( &(list->data[i]), &PlanarConfiguration ) != 1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: unable to add entry 'PlanarConfiguration' to ifd #%d\n", proc, i );
      return( -1 );
    }

    if ( im->vdim <= 2 ) {
      SampleFormatValues = (unsigned short int*)SampleFormat.offset;
      for ( v=0; v<im->vdim; v++ ) {
          if ( im->wordKind == WK_FIXED )
            SampleFormatValues[v] = ( im->sign == SGN_UNSIGNED ) ? 1 : 2;
          else if ( im->wordKind == WK_FLOAT )
            SampleFormatValues[v] = 3;
          else
            SampleFormatValues[v] = 1;
      }
    }
    else {
      SampleFormat.values = (char*)malloc( im->vdim * sizeof(unsigned short int) );
      if ( SampleFormat.values == (char*)NULL ) {
        if ( _verbose_ || _debug_ )
          fprintf( stderr, "%s: allocation failed for 'SampleFormat' for ifd #%d\n", proc, i );
        return( -1 );
      }
      SampleFormatValues = (unsigned short int*)SampleFormat.values;
      SampleFormat.valuesLength = im->vdim * sizeof(unsigned short int);
      for ( v=0; v<im->vdim; v++ ) SampleFormatValues[v] = im->wdim * 8;
      for ( v=0; v<im->vdim; v++ ) {
          if ( im->wordKind == WK_FIXED )
            SampleFormatValues[v] = ( im->sign == SGN_UNSIGNED ) ? 1 : 2;
          else if ( im->wordKind == WK_FLOAT )
            SampleFormatValues[v] = 3;
          else
            SampleFormatValues[v] = 1;
      }
    }
    if ( _addTiffIfdEntryToIfd( &(list->data[i]), &SampleFormat ) != 1 ) {
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: unable to add entry 'SampleFormat' to ifd #%d\n", proc, i );
      return( -1 );
    }
    SampleFormat.values = (char*)NULL;
    SampleFormat.valuesLength = 0;

  }

  return( 1 );
}



int writeLiteTiffImage( char *name, _image *im )
{
  char *proc = "writeLiteTiffImage";
  _tiffHeader header;
  _tiffIfdList ifdList;
  unsigned int s;
  int i, e;
  size_t headerLength;
  size_t valuesLength;
  size_t offset;
  int eStripOffsets;
  int eStripByteCounts;
  unsigned int *StripOffsetsValues;
  unsigned int *StripByteCountsValues;



  _initTiffHeader( &header );
  if ( _fillTiffHeader( &header ) != 1 ) {
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: unable to fill header\n", proc );
    return( -1 );
  }

  _initTiffIfdList( &ifdList );
  if ( _allocTiffIfdList( &ifdList, im->zdim ) != 1 ) {
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: ifd list allocation failed\n", proc );
    return( -1 );
  }

  /* fill header
   */
  if ( _fillTiffHeader( &header ) != 1 ) {
    if ( _verbose_ || _debug_ ) {
      fprintf( stderr, "%s: unable to fill header\n", proc );
    }
    return( -1 );
  }

  /* fill IFDs
   */
  if ( _addGeometryEntries( &ifdList, im ) != 1 ) {
    _freeTiffIfdList( &ifdList );
    if ( _verbose_ || _debug_ ) {
      fprintf( stderr, "%s: unable to add geometry information to IFDs\n", proc );
    }
    return( -1 );
  }



  /* compute buffer lengths
   */

  /* header +IFDs length
   * header = 8;
   * IFD = 2 + 12 * entries + 4
   */
  headerLength = 8;
  for( i=0; i<ifdList.n_data; i++ ) {
    headerLength += 2 + 12 * ifdList.data[i].n_data + 4;
    ifdList.data[i].offset = headerLength;
  }
  ifdList.data[ ifdList.n_data-1 ].offset = 0;

  /* offsets for the values to be written
   */
  offset = headerLength;
  valuesLength = 0;
  for( i=0; i<ifdList.n_data; i++ ) {
    for ( e=0; e<ifdList.data[i].n_data; e++ ) {
      if ( ifdList.data[i].data[e].valuesLength == 0 )
        continue;
      valuesLength += ifdList.data[i].data[e].valuesLength;
      *(ifdList.data[i].data[e].offset) = offset;
      offset += ifdList.data[i].data[e].valuesLength;
    }
  }

  /* offset for the data to be written
   */
  for( i=0; i<ifdList.n_data; i++ ) {

    eStripOffsets = -1;
    eStripByteCounts = -1;
    for ( e=0; e<ifdList.data[i].n_data; e++ ) {
      switch( *(ifdList.data[i].data[e].tag) ) {
      default : break;
      case 273 : eStripOffsets = e; break;
      case 279 : eStripByteCounts = e; break;
      }
    }

    if ( *(ifdList.data[i].data[eStripOffsets].length) == 1 ) {
      *(ifdList.data[i].data[eStripOffsets].offset) = offset;
      offset += *(ifdList.data[i].data[eStripByteCounts].offset);
    }
    else {
      StripOffsetsValues = (unsigned int*)(ifdList.data[i].data[eStripOffsets].values);
      StripByteCountsValues = (unsigned int*)(ifdList.data[i].data[eStripByteCounts].values);
      for ( s=0; s<*(ifdList.data[i].data[eStripOffsets].length); s++ ) {
        StripOffsetsValues[s] = offset;
        offset += StripByteCountsValues[s];
      }
    }

  }


  if ( _debug_ >= 2 ) {
    fprintf( stderr, "--------------------------------------------------\n" );
    fprintf( stderr, "--- IFDs   length = %lu\n", headerLength );
    fprintf( stderr, "--- values length = %lu\n", valuesLength );
    _fprintTiffIfdList( stderr, (char*)NULL, &ifdList, header.swap );
  }


  /* write image
   */
  _openWriteImage( im, name );

  if ( !im->fd ) {
    _freeTiffIfdList( &ifdList );
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: error: unable to open file \'%s\'\n", proc, name );
    return( ImageIO_OPENING );
  }

  if ( ImageIO_write( im, (void*)header.data, 8 ) != 8 ) {
    ImageIO_close( im );
    im->openMode = OM_CLOSE;
    _freeTiffIfdList( &ifdList );
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: error: unable to write header\n", proc  );
    return( -1 );
  }

  for( i=0; i<ifdList.n_data; i++ ) {
    if ( ImageIO_write( im, (void*)&(ifdList.data[i].n_data), 2 ) != 2 ) {
      ImageIO_close( im );
      im->openMode = OM_CLOSE;
      _freeTiffIfdList( &ifdList );
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: error: unable to write number of entries of IFD #%d\n", proc, i );
      return( -1 );
    }
    for ( e=0; e<ifdList.data[i].n_data; e++ ) {
      if ( ImageIO_write( im, (void*)ifdList.data[i].data[e].data, 12 ) != 12 ) {
        ImageIO_close( im );
        im->openMode = OM_CLOSE;
        _freeTiffIfdList( &ifdList );
        if ( _verbose_ || _debug_ )
          fprintf( stderr, "%s: error: unable to write entry #%d of IFD #%d\n", proc, e, i );
        return( -1 );
      }
    }
    if ( ImageIO_write( im, (void*)&(ifdList.data[i].offset), 4 ) != 4 ) {
      ImageIO_close( im );
      im->openMode = OM_CLOSE;
      _freeTiffIfdList( &ifdList );
      if ( _verbose_ || _debug_ )
        fprintf( stderr, "%s: error: unable to write next offset of IFD #%d\n", proc, i );
      return( -1 );
    }

  }

  for( i=0; i<ifdList.n_data; i++ ) {
    for ( e=0; e<ifdList.data[i].n_data; e++ ) {
      if ( ifdList.data[i].data[e].valuesLength == 0 )
        continue;
      if ( ImageIO_write( im, (void*)ifdList.data[i].data[e].values,
                          ifdList.data[i].data[e].valuesLength ) != ifdList.data[i].data[e].valuesLength ) {
        ImageIO_close( im );
        im->openMode = OM_CLOSE;
        _freeTiffIfdList( &ifdList );
        if ( _verbose_ || _debug_ )
          fprintf( stderr, "%s: error: unable to write values of entry #%d of IFD #%d\n", proc, e, i );
        return( -1 );
      }
    }
  }

  if ( ImageIO_write( im, im->data, im->vdim*im->xdim*im->ydim*im->zdim*im->wdim ) != im->vdim*im->xdim*im->ydim*im->zdim*im->wdim ) {
    ImageIO_close( im );
    im->openMode = OM_CLOSE;
    _freeTiffIfdList( &ifdList );
    if ( _verbose_ || _debug_ )
      fprintf( stderr, "%s: error: unable to write image data\n", proc );
    return( -1 );
  }

  ImageIO_close( im );
  im->openMode = OM_CLOSE;
  _freeTiffIfdList( &ifdList );

  return( 1 );
}



/************************************************************
 *
 *
 *
 ************************************************************/



PTRIMAGE_FORMAT createLiteTiffFormat()
{
  PTRIMAGE_FORMAT f=(PTRIMAGE_FORMAT) ImageIO_alloc(sizeof(IMAGE_FORMAT));
  _initImageFormat( f );

  f->testImageFormat = &testLiteTiffHeader;
  f->readImageHeader = &readLiteTiffImage;
  f->readImageData = NULL;
  /* f->writeImageHeader = &writeLiteTiffImage; */
  f->writeImage = &writeLiteTiffImage;

  strcpy(f->readingFileExtension,".tif,.TIF,.tiff,.TIFF,.lsm");
  strcpy(f->writingFileExtension,".tif,.TIF,.tiff,.TIFF,.tif.gz,.TIF.gz,.tiff.gz,.TIFF.gz");
  strcpy(f->realName,"HomeMadeLightTiff");
  return f;
}
