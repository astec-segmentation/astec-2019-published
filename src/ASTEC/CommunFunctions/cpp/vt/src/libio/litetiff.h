
#ifndef LITETIFF_H
#define LITETIFF_H

#ifdef __cplusplus
extern "C" {
#endif

#include <ImageIO.h>

extern void setVerboseInLiteTiff( int v );
extern void incrementVerboseLiteTiff( );
extern void setDebugInLiteTiff( int d );
extern void incrementDebugInLiteTiff( );

extern PTRIMAGE_FORMAT createLiteTiffFormat();

#ifdef __cplusplus
}
#endif

#endif
