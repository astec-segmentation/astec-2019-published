####################################################################
#
# $Id$
#
# Copyright (c) INRIA 2013
#
# AUTHOR:
# Gael Michelin (gael.michelin@inria.fr)
# From Gregoire Malandain (gregoire.malandain@inria.fr)
# 
#

# Name
project(libtools)
set(LIB_NAME tools)


## #################################################################
## Input Files 
## #################################################################
# Inludes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})
include_directories(${VT_INCLUDE_DIRS})
include_directories(${NOESIS_INCLUDE_DIRS})
include_directories(${BLOCKMATCHING_INCLUDE_DIRS})

# Src files
SET(SRC_FILES  vt_bound.c
               vt_overlapTools.c
               vt_tree.c
  )

# Exe names
SET(EXE_NAMES
  boundingboxes
  deleteBorders
  bordersMask
  nlabels
  labelBorders
  mergeLabels
  fuseLabels
  associateLabels
  extractLabels
  segmentationOverlapping
  overlapPruning
  overlapAnalysis
  overlapInterpretation
  pointSetDistances
  thicknessMeasure
  interfacesRegistrationEvaluation
  pointSetsNNSubsets
  multiplyPointSetCoordinates
  pointSetToImage
  pointSetsDistanceAdjusting
  nonZerosImage
  labelBarycenters
  )

# For exe name with suffixes needed: presents in another project
#SET(EXE_NAMES_DUPLICATED 
#  skiz # vt-exec
#)
#message(WARNING "Bin names duplicated : skiz (vt-exec)")

## #################################################################
## Build
## #################################################################

add_library(${LIB_NAME}  ${SRC_FILES})
target_link_libraries( ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt noesis )

foreach(E ${EXE_NAMES_DUPLICATED})
  add_executable(libelf_${E} ${E}.c)
  target_link_libraries(libelf_${E} ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt noesis)
endforeach(E)

foreach(E ${EXE_NAMES})
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt noesis)
endforeach(E)



