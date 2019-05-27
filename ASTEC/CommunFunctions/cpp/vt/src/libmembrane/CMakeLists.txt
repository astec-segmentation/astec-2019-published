####################################################################
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
project(membrane)
set(LIB_NAME membraneLIB)


## #################################################################
## Input Files 
## #################################################################
# Inludes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})
include_directories(${VT_INCLUDE_DIRS})
include_directories(${FILTERS_INCLUDE_DIRS})
include_directories(${BLOCKMATCHING_INCLUDE_DIRS})

# Src files
SET(SRC_FILES  api-membrane.c
               api-TVmembrane.c
               mt_membrane2D.c
	       mt_membrane3D.c
	       mt_dice.c
	       mt_anisotropicHist.c
	       mt_mergeSegmentations.c
	       vt_symmetryPlane.c
	       vt_planeRegistration.c
               vt_tree.c
	       )

# Exe names
SET(EXE_NAMES membrane
        test-fields
        test-TVsample
        TVmembrane
        dice
        anisotropicHist
        anisotropicThreshold
        directionHistogram
        directionHistogramMaxima
        symmetryPlane
        mergeSegmentations
        artefacts_acylYFP
        addPlane
        diceMaximisation
        planeRegistration
        pointCloudRegistration
        progressiveThreshold
        test-tri
        test-kd-tree
        ICP
        )
  
## #################################################################
## Build
## #################################################################

# build library and link
add_library(${LIB_NAME} ${LIB_STATUS} ${SRC_FILES})
target_link_libraries( ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt tubes blockmatchingLIB)

# build execs and link
foreach(E ${EXE_NAMES})
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt tubes blockmatchingLIB)
endforeach(E)

