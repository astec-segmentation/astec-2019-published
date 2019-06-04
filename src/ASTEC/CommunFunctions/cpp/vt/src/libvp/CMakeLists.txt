####################################################################
#
# $Id$
#
# Copyright (c) INRIA 2013
#
# AUTHOR:
# Etienne Delclaux (etienne.delclaux@inria.fr)
# From Gregoire Malandain (gregoire.malandain@inria.fr)
# 
#

# Name
project(libvp)
set(LIB_NAME vp)

## #################################################################
## Input Files 
## #################################################################
# Inludes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})
include_directories(${VT_INCLUDE_DIRS})
include_directories(${DAVID_INCLUDE_DIRS})
include_directories(${BLOCKMATCHING_INCLUDE_DIRS})


SET( SRC_FILES
        api-adhocFuse.c
	api-cellProperties.c
        api-cellfilter.c
	api-linearCombination.c
	api-maskSeeds.c
	bal-2DcellPropertiesTest.c
	bal-3DcellPropertiesTest.c
	bal-cellProperties.c
	surface-lindblad.c
	surface-windreich.c
	vt_adhocFuse.c
        vt_cellfilter.c
	vt_maskSeeds.c
        vt_temporal-matching.c
        vt_removeLine.c
)

SET( EXE_NAMES 
        cellfilter
	mc-adhocFuse
	mc-cellProperties
	mc-extractMIPembryo
	mc-linearCombination
	mc-maskSeeds
        mc-removeLine
	temporal-matching
	test-2DcellProperties
	test-3DcellProperties
)
  
## #################################################################
## Build
## #################################################################

add_library(${LIB_NAME} ${SRC_FILES})
target_link_libraries( ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt david
blockmatchingLIB )

foreach(E ${EXE_NAMES})
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt
  david blockmatchingLIB)
endforeach(E)
