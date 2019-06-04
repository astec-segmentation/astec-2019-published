/*************************************************************************
 * bal-vtk.c -
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

#include <vtmalloc.h>

#include <bal-vtk.h>


static int _debug_ = 0;
static int _verbose_ = 1;



void BAL_SetVerboseInBalVtk( int v )
{
  _verbose_ = v;
}

void BAL_IncrementVerboseInBalVtk( )
{
  _verbose_ ++;
}

void BAL_DecrementVerboseInBalVtk( )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}


void BAL_SetDebugInBalVtk( int v )
{
  _debug_ = v;
}

void BAL_IncrementDebugInBalVtk( )
{
  _debug_ ++;
}

void BAL_DecrementDebugInBalVtk( )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}





/************************************************************
 *
 *
 *
 ************************************************************/

static void BAL_InitVtkData( bal_vtkData *d )
{
  d->attribute = _VTK_NO_ATTRIBUTE_;
  d->type = _VTK_UNKNOWN_;
  d->n_primitives = 0;
  d->n_data = 0;
  d->data = (void*)NULL;
  d->name = (char*)NULL;
  d->dataLength = 0;
  d->iodone = 0;
}





static void BAL_FreeVtkData( bal_vtkData *d )
{
  if ( d->data != (void*)NULL)
      vtfree( d->data );
  if ( d->name != (char*)NULL)
      vtfree( d->name );
  BAL_InitVtkData( d );
}





void BAL_InitVtkDataSet( bal_vtkDataSet *l )
{
  l->n_data = 0;
  l->n_allocated_data = 0;
  l->data = (bal_vtkData*)NULL;

  l->unit = REAL_UNIT;
  l->vx = 1.0;
  l->vy = 1.0;
  l->vz = 1.0;
}





void BAL_FreeVtkDataSet( bal_vtkDataSet *l )
{
  int d;
  for ( d=0; d<l->n_data; d++ )
    BAL_FreeVtkData( &(l->data[d]) );
  if ( l->data != (bal_vtkData*)NULL )
    vtfree( l->data );
  BAL_InitVtkDataSet( l );
}





static int _nextVTkData( bal_vtkDataSet *l )
{
  char *proc = "_nextVTkData";
  bal_vtkData *data;
  int i;
  int n = l->n_allocated_data;

  if ( l->n_data == l->n_allocated_data ) {
    n += 10;
    data = (bal_vtkData*)vtmalloc( n * sizeof(bal_vtkData), "data", proc );
    if ( data == (bal_vtkData*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: allocation failed\n", proc );
      return( -1 );
    }
    for ( i=l->n_data; i<n; i++ )
      BAL_InitVtkData( &(data[i]) );
    if ( l->n_data > 0 ) {
      memcpy( data, l->data, l->n_data*sizeof(bal_vtkData) );
      vtfree( l->data );
    }
    l->data = data;
    l->n_allocated_data = n;
  }
  l->n_data ++;
  return( l->n_data-1 );
}





int BAL_CopyVtkDataSet( bal_vtkDataSet *theSet, bal_vtkDataSet *resSet )
{
  char *proc = "BAL_CopyVtkDataSet";
  int p;

  if ( theSet->n_data <= 0 )
    return( 1 );

  resSet->data = (bal_vtkData*)vtmalloc( theSet->n_data*sizeof(bal_vtkData), "resSet->data", proc );
  if ( resSet->data == (bal_vtkData*)NULL ) {
    if ( _verbose_ )
      fprintf(stderr, "%s: allocation error\n", proc );
    return( -1 );
  }

  resSet->n_data = resSet->n_allocated_data = theSet->n_data;

  for ( p=0; p<theSet->n_data; p++ ) {

    resSet->data[p] = theSet->data[p];

    if ( theSet->data[p].name != (char*)NULL ) {
      resSet->data[p].name = (char*)vtmalloc( strlen(theSet->data[p].name)+1, "res->data[p].name", proc );
      if ( resSet->data[p].name == (char*)NULL ) {
        BAL_FreeVtkDataSet( resSet );
        if ( _verbose_ )
          fprintf(stderr, "%s: allocation error (name)\n", proc );
        return( -1 );
      }
      memcpy( resSet->data[p].name, theSet->data[p].name, strlen(theSet->data[p].name)+1 );
    }

    if ( theSet->data[p].data != (void*)NULL ) {
      resSet->data[p].data = (void*)vtmalloc( theSet->data[p].dataLength, "res->data[p].data", proc );
      if ( resSet->data[p].data == (char*)NULL ) {
        BAL_FreeVtkDataSet( resSet );
        if ( _verbose_ )
          fprintf(stderr, "%s: allocation error (data)\n", proc );
        return( -1 );
      }
      memcpy( resSet->data[p].data, theSet->data[p].data, theSet->data[p].dataLength );
    }

  }

  resSet->unit = theSet->unit;
  resSet->vx = theSet->vx;
  resSet->vy = theSet->vy;
  resSet->vz = theSet->vz;

  return( 1 );
}





/************************************************************
 *
 *
 *
 ************************************************************/



static void _fprintfString( FILE *f, char *str )
{
  char *s = str;
  while ( *s != '\0' && *s != '\n' ) {
    fprintf( f, "%c", *s );
    s++;
  }
}



static void _fprintfEnumVtkDataEncoding( FILE *f, enumVtkDataEncoding type )
{
    switch( type ){
    default : break;
    case _VTK_BINARY_DATA_ : fprintf( f, "BINARY" ); break;
    case _VTK_ASCII_DATA_ : fprintf( f, "ASCII" ); break;
    }
}



static void _fprintfEnumVtkAttributeType( FILE *f, enumVtkAttributeType type )
{
    switch( type ){
    default : break;
    case _VTK_CELL_ATTRIBUTE_ : fprintf( f, "CELL_DATA " ); break;
    case _VTK_POINT_ATTRIBUTE_ : fprintf( f, "POINT_DATA " ); break;
    }
}



static void _fprintfEnumVtkData( FILE *f, enumVtkData type )
{
    switch( type ){
    default : fprintf( f, "default" ); break;
    case _VTK_UNKNOWN_ : fprintf( f, "unknown" ); break;
    case _VTK_POINTS_ : fprintf( f, "POINTS" ); break;
    case _VTK_LINES_ : fprintf( f, "LINES" ); break;
    case _VTK_POLYGONS_ : fprintf( f, "POLYGONS" ); break;
    case _VTK_TRIANGLE_STRIPS_ : fprintf( f, "TRIANGLE STRIPS" ); break;
    case _VTK_VERTICES_ : fprintf( f, "VERTICES" ); break;
    case _VTK_SCALARS_ : fprintf( f, "SCALARS" ); break;
    case _VTK_COLOR_SCALARS_ : fprintf( f, "COLOR_SCALARS_" ); break;
    case _VTK_LOOKUP_TABLE_ : fprintf( f, "LOOKUP_TABLE" ); break;
    case _VTK_VECTORS_ : fprintf( f, "VECTORS" ); break;
    case _VTK_NORMALS_ : fprintf( f, "NORMALS" ); break;
    case _VTK_TEXTURE_COORDINATES_ : fprintf( f, "TEXTURE_COORDINATES" ); break;
    case _VTK_TENSORS_ : fprintf( f, "TENSORS_" ); break;
    case _VTK_FIELD_ : fprintf( f, "FIELD" ); break;
    }
}





static void _fprintfVtkData( FILE *f, bal_vtkData *d )
{
  _fprintfEnumVtkAttributeType( f, d->attribute );
  _fprintfEnumVtkData( f, d->type );
  fprintf( f," n_primitives=%d", d->n_primitives );
  fprintf( f," n_data=%d", d->n_data );
  fprintf( f," (%p)", d->data );
  if ( d->name != (char*)NULL ) {
      fprintf( f, " '%s'", d->name );
      fprintf( f," (%p)", d->name );
  }
  fprintf( f," dataLength=%d", d->dataLength );
  fprintf( f," iodone=%d", d->iodone );
  fprintf( f, "\n" );
}





void BAL_FprintfVtkDataSet( FILE *f, bal_vtkDataSet *s )
{
  int n;
  for ( n=0; n<s->n_data; n++ ) {
    fprintf( stderr, "--- data #%d\n", n );
    fprintf( stderr, "    " );
    _fprintfVtkData( f, &(s->data[n]) );
  }
}





/************************************************************
 *
 *
 *
 ************************************************************/



static void _swap4octets( bal_vtkData *dataSet )
{
  int i;
  unsigned int u, *in, *out;
  union {
    unsigned char uc[2];
    unsigned short us;
  } twobytes;
  twobytes.us = 255;

  if ( twobytes.uc[1] == 0 ) {
    in = out = (unsigned int *)dataSet->data;
    for ( i=0; i<dataSet->n_data; i++ ) {
        u = *in++;
        *out++ =  (u << 24) | ((u & 0xff00) << 8) | ((u >> 8) & 0xff00) | ((u >> 24) & 0xff);
    }
  }
}





static int _readPointsData( FILE *f, bal_vtkData *dataSet,
                            enumVtkDataEncoding dataEncoding )
{
  char *proc = "_readPointsData";
  int i;
  float *pts = (float*)dataSet->data;

  switch( dataEncoding ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such mode not handled yet\n", proc );
    return( -1 );
  case _VTK_ASCII_DATA_ :
    for ( i=0; i<3*dataSet->n_primitives; i++ ) {
      if ( fscanf( f, "%f", &(pts[i]) ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to read point #%d\n", proc, i );
        return( -1 );
      }
    }
    break;
  case _VTK_BINARY_DATA_ :
    if ( fread( dataSet->data, dataSet->dataLength, 1, f ) < 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to read points\n", proc );
     return( -1 );
    }
    _swap4octets( dataSet );
    break;
  }

  dataSet->iodone = 1;

  return( 1 );
}





static int _readIntData( FILE *f, bal_vtkData *dataSet,
                            enumVtkDataEncoding dataEncoding )
{
  char *proc = "_readIntData";
  int i;
  int *pts = (int*)dataSet->data;

  switch( dataEncoding ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such mode not handled yet\n", proc );
    return( -1 );
  case _VTK_ASCII_DATA_ :
    for ( i=0; i<dataSet->n_data; i++ ) {
      if ( fscanf( f, "%d", &(pts[i]) ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to read value #%d\n", proc, i );
        return( -1 );
      }
    }
    break;
  case _VTK_BINARY_DATA_ :
    if ( fread( dataSet->data, dataSet->dataLength, 1, f ) < 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to read primitive\n", proc );
     return( -1 );
    }
    _swap4octets( dataSet );
    break;
  }

  dataSet->iodone = 1;

  return( 1 );
}





static int _readTensorData( FILE *f, bal_vtkData *dataSet,
                            enumVtkDataEncoding dataEncoding )
{
  char *proc = "_readTensorData";
  int i;
  float *pts = (float*)dataSet->data;

  switch( dataEncoding ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such mode not handled yet\n", proc );
    return( -1 );
  case _VTK_ASCII_DATA_ :
    for ( i=0; i<dataSet->n_primitives; i++ ) {
      if ( fscanf( f, "%f %f %f", &(pts[9*i]), &(pts[9*i+1]), &(pts[9*i+2]) ) != 3 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to read tensor[0-2] #%d\n", proc, i );
        return( -1 );
      }
      if ( fscanf( f, "%f %f %f", &(pts[9*i+3]), &(pts[9*i+4]), &(pts[9*i+5]) ) != 3 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to read tensor[3-5] #%d\n", proc, i );
        return( -1 );
      }
      if ( fscanf( f, "%f %f %f", &(pts[9*i+6]), &(pts[9*i+7]), &(pts[9*i+8]) ) != 3 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to read tensor[6-8] #%d\n", proc, i );
        return( -1 );
      }
    }
    break;
  case _VTK_BINARY_DATA_ :
    if ( fread( dataSet->data, dataSet->dataLength, 1, f ) < 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to read points\n", proc );
     return( -1 );
    }
    _swap4octets( dataSet );
    break;
  }

  dataSet->iodone = 1;

  return( 1 );
}





static int _readVtkPolyData( FILE *f, bal_vtkDataSet *dataSet,
                             enumVtkDataEncoding dataEncoding )
{
  char *proc = "_readVtkPolyData";
  char str[512];
  long int pos;
  int i, nprimitive, ndata;
  int idata;
  int ipoint = -1;
  bal_vtkData *ptrData;
  int hasreaddata;
  int nameLength;
  enumVtkData primitive;
  enumVtkAttributeType dataAttribute = _VTK_NO_ATTRIBUTE_;


  do {

    /* get position
     */
    pos = ftell( f );

    if ( _debug_ >= 2 ) {
      fprintf( stderr, "%s: analysing from %ld\n", proc, pos );
    }

    if ( fgets( str, 512, f ) == NULL ) {
      /* assume end of file
       */
      return( 1 );
    }

    /* skip empty line
     */
    i = 0;
    while ( str[i] == ' ' || str[i] == '\t' )
      i++;
    if ( str[i] == '\n' ) {
      continue;
    }

    if ( 1 && _debug_ >= 2 ) {
      fprintf( stderr, "\t line='" );
      _fprintfString( stderr, str );
      fprintf( stderr, "'\n" );
    }



    /* analysing line
     */
    if ( strncmp( str, "POINTS ", 7 ) == 0 ) {

      dataAttribute = _VTK_NO_ATTRIBUTE_;

      if ( sscanf( &(str[7]), "%d", &nprimitive ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to read POINTS number\n", proc );
        return( -1 );
      }
      for ( i=7; ('0' <= str[i] && str[i] <= '9'); i++ )
        ;
      for ( ; str[i] == ' '; i++ )
        ;
      if ( strncmp( &(str[i]), "float", 5 ) != 0 ) {
        if ( _verbose_ ) {
          fprintf( stderr, "%s: such POINT type not handled yet\n", proc );
          fprintf( stderr, "\t %s\n", &(str[i]) );
        }
        return( -1 );
      }
      idata = _nextVTkData( dataSet );
      if ( idata < 0 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to get next data id\n", proc );
        return( -1 );
      }

      ptrData = &(dataSet->data[idata]);
      ptrData->dataLength = 3 * nprimitive * sizeof( float );
      ptrData->data = (void*)vtmalloc( ptrData->dataLength, "ptrData->dataLength", proc );
      if ( ptrData->data == (void*)NULL ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: allocation error\n", proc );
        return( -1 );
      }
      ptrData->type = _VTK_POINTS_;
      ptrData->n_primitives = nprimitive;
      ptrData->n_data = 3 * nprimitive;
      ptrData->iodone = 0;

    }

    else if ( strncmp( str, "LINES ", 6 ) == 0
              || strncmp( str, "POLYGONS ", 9 ) == 0
              || strncmp( str, "TRIANGLE_STRIPS ", 16 ) == 0
              || strncmp( str, "VERTICES ", 9 ) == 0 ) {

      dataAttribute = _VTK_NO_ATTRIBUTE_;

      if ( strncmp( str, "LINES ", 6 ) == 0 ) {
        primitive = _VTK_LINES_;
        nameLength = 6;
      }
      else if ( strncmp( str, "POLYGONS ", 9 ) == 0 ) {
        primitive = _VTK_POLYGONS_;
        nameLength = 9;
      }
      else if ( strncmp( str, "TRIANGLE_STRIPS ", 16 ) == 0 ) {
        primitive = _VTK_TRIANGLE_STRIPS_;
        nameLength = 16;
      }
      else if ( strncmp( str, "VERTICES ", 9 ) == 0 ) {
        primitive = _VTK_VERTICES_;
        nameLength = 9;
      }
      else {
        primitive = _VTK_UNKNOWN_;
        nameLength = 0;
      }

      if ( sscanf( &(str[nameLength]), "%d", &nprimitive ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to read primitive number\n", proc );
        return( -1 );
      }
      for ( i=nameLength; ('0' <= str[i] && str[i] <= '9'); i++ )
        ;
      for ( ; str[i] == ' '; i++ )
        ;
      if ( sscanf( &(str[i]), "%d", &ndata ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to read primitive data length\n", proc );
        return( -1 );
      }
      idata = _nextVTkData( dataSet );
      if ( idata < 0 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to get next data id\n", proc );
        return( -1 );
      }

      ptrData = &(dataSet->data[idata]);
      ptrData->dataLength = ndata * sizeof(int);
      ptrData->data = (void*)vtmalloc( ptrData->dataLength, "ptrData->dataLength", proc );
      if ( ptrData->data == (void*)NULL ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: allocation error\n", proc );
        return( -1 );
      }
      ptrData->type = primitive;
      ptrData->n_primitives = nprimitive;
      ptrData->n_data = ndata;
      ptrData->iodone = 0;

    }

    else if ( strncmp( str, "SCALARS ", 8 ) == 0
              || strncmp( str, "COLOR_SCALARS ", 14 ) == 0
              || strncmp( str, "LOOKUP_TABLE ", 13 ) == 0
              || strncmp( str, "VECTORS ", 8 ) == 0
              || strncmp( str, "NORMALS ", 8 ) == 0
              || strncmp( str, "TEXTURE_COORDINATES ", 20 ) == 0
              || strncmp( str, "TENSORS ", 8 ) == 0
              || strncmp( str, "FIELD ", 6 ) == 0 ) {

      if ( strncmp( str, "SCALARS ", 8 ) == 0 ) {
        primitive = _VTK_SCALARS_;
        nameLength = 8;
      }
      else if ( strncmp( str, "COLOR_SCALARS ", 14 ) == 0 ) {
          primitive = _VTK_COLOR_SCALARS_;
          nameLength = 14;
       }
      else if ( strncmp( str, "LOOKUP_TABLE ", 13 ) == 0 ) {
          primitive = _VTK_LOOKUP_TABLE_;
          nameLength = 13;
      }
      else if ( strncmp( str, "VECTORS ", 8 ) == 0 ) {
          primitive = _VTK_VECTORS_;
          nameLength = 8;
       }
      else if ( strncmp( str, "NORMALS ", 8 ) == 0 ) {
          primitive = _VTK_NORMALS_;
          nameLength = 8;
      }
      else if ( strncmp( str, "TEXTURE_COORDINATES ", 20 ) == 0 ) {
          primitive = _VTK_TEXTURE_COORDINATES_;
          nameLength = 20;
      }
      else if ( strncmp( str, "TENSORS ", 8 ) == 0 ) {
          primitive = _VTK_TENSORS_;
          nameLength = 8;
       }
      else if ( strncmp( str, "FIELD ", 6 ) == 0 ) {
          primitive = _VTK_FIELD_;
          nameLength = 6;
      }
      else {
        primitive = _VTK_UNKNOWN_;
        nameLength = 0;
      }

      switch( primitive ) {
      default :
        if ( _verbose_ ) {
          fprintf( stderr, "%s: such attribute '", proc );
          _fprintfEnumVtkData( stderr, primitive );
          fprintf( stderr, "' not handled yet\n" );
        }
        return( -1 );
      case _VTK_VECTORS_ :
      case _VTK_NORMALS_ :
      case _VTK_TENSORS_ :

        idata = _nextVTkData( dataSet );
        if ( idata < 0 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to get next data id\n", proc );
          return( -1 );
        }
        ptrData = &(dataSet->data[idata]);

        for ( i=nameLength; str[i]!=' '; i++ )
           ;
        ptrData->name = (char*)vtmalloc( i+1-nameLength, "ptrData->name", proc );
        if ( ptrData->name == (char*)NULL ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: allocation error\n", proc );
          return( -1 );
        }
        for ( i=nameLength; str[i]!=' '; i++ ) {
          ptrData->name[i-nameLength] = str[i];
        }
        ptrData->name[i-nameLength] = '\0';
        for ( ; str[i] == ' '; i++ )
          ;
        if ( strncmp( &(str[i]), "float", 5 ) != 0 ) {
          if ( _verbose_ ) {
            fprintf( stderr, "%s: such ATTRIBUTE type not handled yet\n", proc );
            fprintf( stderr, "\t %s\n", &(str[i]) );
          }
          return( -1 );
        }

        ptrData->attribute = dataAttribute;
        ptrData->type = primitive;

        switch( dataAttribute ) {
        default :
        case _VTK_CELL_ATTRIBUTE_ :
          if ( _verbose_ )
            fprintf( stderr, "%s: such data ATTRIBUTE not handled yet\n", proc );
          break;
        case _VTK_POINT_ATTRIBUTE_ :
          ptrData->n_primitives = dataSet->data[ipoint].n_primitives;
          switch( primitive ) {
          default :
            if ( _verbose_ )
              fprintf( stderr, "%s: such attribute not handled yet\n", proc );
            return( -1 );
          case _VTK_VECTORS_ :
          case _VTK_NORMALS_ :
            ptrData->n_data = 3 * ptrData->n_primitives; break;
          case _VTK_TENSORS_ :
            ptrData->n_data = 9 * ptrData->n_primitives; break;
          }
          ptrData->dataLength = ptrData->n_data * sizeof(float);
          ptrData->data = (void*)vtmalloc( ptrData->dataLength, "ptrData->dataLength", proc );
          if ( ptrData->data == (void*)NULL ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: allocation error\n", proc );
            return( -1 );
          }
          ptrData->iodone = 0;
          break;
        }


        break;
      }

    }

    else if ( strncmp( str, "CELL_DATA ", 10 ) == 0 ) {
      dataAttribute = _VTK_CELL_ATTRIBUTE_;
      if ( _verbose_ )
        fprintf( stderr, "%s: such ATTRIBUTES (CELL_DATA) not handled yet\n", proc );
      return( -1 );
    }

    else if ( strncmp( str, "POINT_DATA ", 11 ) == 0 ) {
      dataAttribute = _VTK_POINT_ATTRIBUTE_;
      if ( sscanf( &(str[11]), "%d", &nprimitive ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to read POINTS number\n", proc );
        return( -1 );
      }
      for ( ipoint=-1, i=0; i<dataSet->n_data && ipoint == -1; i++ ) {
        if ( dataSet->data[i].type == _VTK_POINTS_ )
          ipoint = i;
      }
      if ( dataSet->data[ipoint].n_primitives != nprimitive ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: POINT_DATA has a different count than POINTS\n", proc );
        return( -1 );
      }
    }

    /* Reading data
     */
    else {
      if ( _debug_ >= 2 ) {
        fprintf( stderr, "%s: will read data\n", proc );
      }
      fseek( f, pos, SEEK_SET );
      hasreaddata = 0;
      for ( i=0; i<dataSet->n_data; i++ ) {
        if ( dataSet->data[i].iodone == 1 ) continue;
        if ( dataSet->data[i].dataLength == 0 ) continue;
        switch( dataSet->data[i].type ) {
        default :
          if ( _verbose_ )
            fprintf( stderr, "%s: such type not handled yet\n", proc );
          return( -1 );
        case _VTK_POINTS_ :
        case _VTK_VECTORS_ :
        case _VTK_NORMALS_ :
          if ( _readPointsData( f, &(dataSet->data[i]), dataEncoding ) != 1 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: error when reading points\n", proc );
            return( -1 );
          }
          hasreaddata = 1;
          break;
        case _VTK_LINES_ :
        case _VTK_POLYGONS_ :
        case _VTK_TRIANGLE_STRIPS_ :
        case _VTK_VERTICES_ :
          if ( _readIntData( f, &(dataSet->data[i]), dataEncoding ) != 1 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: error when reading primitive\n", proc );
            return( -1 );
          }
          hasreaddata = 1;
          break;
        case _VTK_TENSORS_ :
          if ( _readTensorData( f, &(dataSet->data[i]), dataEncoding ) != 1 ) {
            if ( _verbose_ )
                fprintf( stderr, "%s: error when reading primitive\n", proc );
              return( -1 );
          }
          hasreaddata = 1;
          break;
        }
      }
      if ( hasreaddata == 0 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: no data were read?!\n", proc );
        return( -1 );
      }
    }


  } while( 1 );

  /* points
   */
  if ( fgets( str, 512, f ) == NULL ) {
  }

  return( 1 );
}





int BAL_ReadVtkDataSet( bal_vtkDataSet *dataset, char *name )
{
  char *proc = "BAL_ReadVtkDataSet";
  FILE *f;
  char str[512];
  enumVtkDataEncoding dataEncoding;

  f = fopen( name, "r" );
  if ( f == (FILE*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to open file '%s'\n", proc, name );
    return( -1 );
  }

  /* read header
   */
  if ( fgets( str, 512, f ) == NULL ) {
    fclose( f );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to read line #1 of file '%s'\n", proc, name );
    return( -1 );
  }
  if  ( strncmp( str, "# vtk DataFile Version", 22 ) != 0 ) {
    fclose( f );
    if ( _verbose_ ) {
      fprintf( stderr, "%s: file '%s' seems not be vtk file\n", proc, name );
      fprintf( stderr, "\t line #1 is '%s'\n", str );
    }
    return( -1 );
  }

  if ( fgets( str, 512, f ) == NULL ) {
    fclose( f );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to read line #2 of file '%s'\n", proc, name );
    return( -1 );
  }

  if ( fgets( str, 512, f ) == NULL ) {
    fclose( f );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to read line #3 of file '%s'\n", proc, name );
    return( -1 );
  }
  if ( strncmp( str, "BINARY", 6 ) == 0 ) {
    dataEncoding = _VTK_BINARY_DATA_;
  }
  else if ( strncmp( str, "ASCII", 5 ) == 0 ) {
    dataEncoding = _VTK_ASCII_DATA_;
  }
  else {
    fclose( f );
    if ( _verbose_ )
      fprintf( stderr, "%s: unknown data type '%s'\n", proc, str );
    return( -1 );
  }



  /*
   */
  if ( fgets( str, 512, f ) == NULL ) {
    fclose( f );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to read line\n", proc );
    return( -1 );
  }


  if ( strncmp( str, "DATASET ", 8 ) != 0 ) {
      fclose( f );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to recognize DATASET format\n", proc );
      return( -1 );
  }

  if ( strncmp( &(str[8]), "POLYDATA", 8 ) == 0 ) {
    if ( _readVtkPolyData( f, dataset, dataEncoding ) != 1 ) {
      fclose( f );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when reading PolyData\n", proc );
      return( -1 );
    }
  }
  else {
      fclose( f );
      if ( _verbose_ ) {
        fprintf( stderr, "%s: DATASET format not handled yet\n", proc );
        fprintf( stderr, "\t %s\n", str );
      }
      return( -1 );
  }

  fclose( f );

  return( 1 );
}





/************************************************************
 *
 *
 *
 ************************************************************/





static int _writePointsData( FILE *f, bal_vtkData *dataSet,
                            enumVtkDataEncoding dataEncoding )
{
  char *proc = "_writePointsData";
  int p;
  float *pts = (float*)dataSet->data;

  switch( dataEncoding ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such data encoding not handled yet\n", proc );
    return( -1 );
  case _VTK_ASCII_DATA_ :
    for ( p=0; p<dataSet->n_primitives; p++ ) {
      if ( fprintf( f, "%f %f %f\n", pts[3*p], pts[3*p+1], pts[3*p+2] ) < 0 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: error when writing data\n", proc );
        return( -1 );
      }
    }
    fprintf( f, "\n" );
    break;
  case _VTK_BINARY_DATA_ :
    _swap4octets( dataSet );
    if ( fwrite( dataSet->data, dataSet->dataLength, 1, f ) < 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: error when writing data\n", proc );
        return( -1 );
    }
    _swap4octets( dataSet );
    fprintf( f, "\n" );
    break;
  }

  dataSet->iodone = 1;

  return( 1 );
}





static int _writeIntData( FILE *f, bal_vtkData *dataSet,
                            enumVtkDataEncoding dataEncoding )
{
  char *proc = "_writeIntData";
  int i, j, n;
  int *data = (int*)dataSet->data;

  switch( dataEncoding ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such data encoding not handled yet\n", proc );
    return( -1 );
  case _VTK_ASCII_DATA_ :
    for ( i=0; i<dataSet->n_data; ) {
      n = data[i];
      fprintf( f, "%d ", n );
      i++;
      for ( j=0; j<n; j++, i++ ) {
        fprintf( f, "%d", data[i] );
        if ( j < n-1 ) fprintf( f, " " );
      }
      fprintf( f, "\n" );
    }
    fprintf( f, "\n" );
    break;
  case _VTK_BINARY_DATA_ :
    _swap4octets( dataSet );
    if ( fwrite( dataSet->data, dataSet->dataLength, 1, f ) < 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: error when writing data\n", proc );
        return( -1 );
    }
    _swap4octets( dataSet );
    fprintf( f, "\n" );
    break;
  }

  dataSet->iodone = 1;

  return( 1 );
}





int BAL_WriteVtkDataSet( bal_vtkDataSet *dataSet, char *name, enumVtkDataEncoding dataEncoding )
{
  char *proc = "BAL_WriteVtkDataSet";
  FILE *f;
  int i, j;
  int ipoint = -1;

  for ( i=0; i<dataSet->n_data; i++ )
      dataSet->data[i].iodone = 0;

  f = fopen( name, "w" );
  if ( f == (FILE*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to open file '%s'\n", proc, name );
    return( -1 );
  }

  /* write header
   */
  fprintf( f, "# vtk DataFile Version 2.0\n" );
  fprintf( f, "%s\n", proc );
  _fprintfEnumVtkDataEncoding( f, dataEncoding );
  fprintf( f, "\n" );
  fprintf( f, "DATASET POLYDATA\n" );


  /* write points
   */
  for ( i=0; i<dataSet->n_data; i++ ) {
    switch( dataSet->data[i].type ) {
    default :
      break;
    case _VTK_POINTS_ :
      fprintf( f, "POINTS %d float\n", dataSet->data[i].n_primitives );
      if ( _writePointsData( f, &(dataSet->data[i]), dataEncoding ) != 1 ) {
        fclose( f );
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to write POINTS\n", proc );
        return( -1 );
      }
      break;
    }
  }

  /* write structures
   */

  for ( i=0; i<dataSet->n_data; i++ ) {
    switch( dataSet->data[i].type ) {
    default :
      break;
    case _VTK_LINES_ :
    case _VTK_POLYGONS_ :
      if ( dataSet->data[i].type == _VTK_LINES_ )
        fprintf( f, "LINES %d %d\n", dataSet->data[i].n_primitives, dataSet->data[i].n_data );
      else
        fprintf( f, "POLYGONS %d %d\n", dataSet->data[i].n_primitives, dataSet->data[i].n_data );
      if ( _writeIntData( f, &(dataSet->data[i]), dataEncoding ) != 1 ) {
        fclose( f );
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to write primitives\n", proc );
        return( -1 );
      }
      break;
    }
  }

  /* write data attributes
   */

  for ( i=0; i<dataSet->n_data; i++ ) {
    switch( dataSet->data[i].attribute ) {
    default :
      break;
    case _VTK_POINT_ATTRIBUTE_ :
      for ( ipoint=-1, j=0; j<dataSet->n_data && ipoint == -1; j++ ) {
        if ( dataSet->data[j].type == _VTK_POINTS_ )
          ipoint = j;
      }
      fprintf( f, "POINT_DATA %d\n", dataSet->data[ipoint].n_primitives );
      switch( dataSet->data[i].type ) {
      default :
        break;
      case _VTK_SCALARS_ :
      case _VTK_COLOR_SCALARS_ :
      case _VTK_LOOKUP_TABLE_ :
      case _VTK_VECTORS_ :
      case _VTK_TEXTURE_COORDINATES_ :
      case _VTK_TENSORS_ :
      case _VTK_FIELD_ :
          break;
      case _VTK_NORMALS_ :
        fprintf( f, "NORMALS Normals float\n" );
        if ( _writePointsData( f, &(dataSet->data[i]), dataEncoding ) != 1 ) {
          fclose( f );
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to write NORMALS\n", proc );
          return( -1 );
        }
      }

      break;
    }
  }

  fclose( f );
  return( 1 );
}





/************************************************************
 *
 *
 *
 ************************************************************/










/************************************************************
 *
 *
 *
 ************************************************************/


