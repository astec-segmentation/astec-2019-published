############################################################
#
# $Id$
#
# Copyright (c) INRIA 2012
#
# AUTHOR:
# Etienne Delclaux (etienne.delclaux@inria.fr)
# From Gregoire Malandain (gregoire.malandain@inria.fr)
# 
#

# Name
project(libfilters)
SET(LIB_NAME tubes)


## #################################################################
## Input Files 
## #################################################################
# Inludes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})
include_directories(${VT_INCLUDE_DIRS})
# include(${VTK_USE_FILE})

# SRC files
SET(SRC_FILES tube3D.c
	tube3Dutil.c
	vt_tube3D.c
	vt_tube3Delf.c
	vt_tube2D.c	
	vt_tube2Dmatlab.c
	vt_contoursMatlab.c
	vt_tubeutils.c
	vt_link.c
	vt_contours.c
	vt_isocontours.c
)


SET(EXE_NAMES vessel
	road
	reechSlices
	test-contours
	test-hessien3D
	reechSlices2
	test-contours-3D
	drawCenters
	)


if (  VTK_FOUND ) 
  SET( SRC_FILES ${SRC_FILES}
          D_graphe.cpp
	  G_classif.cpp
  )
  SET( EXECPP_NAMES extraction_arbre
	  classification
  ) 
else( VTK_FOUND )
  message( WARNING, " vtk was either not found or too old" )
  message( WARNING, " library ' libfilters' building will NOT be complete" )  
  message( WARNING, " executables ' extraction_arbre' and 'classification' will NOT be built" )
endif( VTK_FOUND )

	
## #################################################################
## Build 
## #################################################################

# Lib
add_library(${LIB_NAME} ${SRC_FILES})
if (  VTK_FOUND ) 
  target_link_libraries( ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt ${VTK_LIBRARIES})
else( VTK_FOUND )
  target_link_libraries( ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt)
endif( VTK_FOUND )

# binaries
foreach(E ${EXE_NAMES})
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt)
endforeach(E)

if (  VTK_FOUND ) 
  foreach(E ${EXECPP_NAMES})
    add_executable(${E} ${E}.cpp)
    target_link_libraries(${E} ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt ${VTK_LIBRARIES})
  endforeach(E)
endif( VTK_FOUND )



## #################################################################
## TEST 
## #################################################################
#if(${FILTERS_BUILD_TEST})
#  
#  SET(TEST_NAME test-tube3D)
#
#  add_executable(${TEST_NAME} ${TEST_NAME}.c)
#  target_link_libraries(${TEST_NAME} ${LIB_NAME} basic io vt)
#  add_test(${NAME}_${T} ${NAME}_${TEST_NAME})
#endif(${FILTERS_BUILD_TEST})


## #################################################################
## History. 
## #################################################################

# SET(INSIGHT_NAME test-tube3D-insight)
# SET(PURIFY_NAME test-tube3D-purify)
# SET(QUANTIFY_NAME test-tube3D-quantify)
