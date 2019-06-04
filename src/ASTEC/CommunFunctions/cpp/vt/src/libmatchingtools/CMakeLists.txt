####################################################################
#
# $Id$
#
# Copyright (c) INRIA 2013
#
# AUTHOR:
# Gregoire Malandain (gregoire.malandain@inria.fr)
# From Etienne Delclaux (etienne.delclaux@inria.fr)
#

# Name
project(libmatchingtools)
set(LIB_NAME matchingtools)


## #################################################################
## Input Files 
## #################################################################
# Includes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})
include_directories(${BLOCKMATCHING_INCLUDE_DIRS})

# Src files
SET(SRC_FILES bal-transformation-averaging.c
	      bal-transformation-list-tools.c
              bal-transformation-propagation.c
)

# Exe names
SET(EXE_NAMES changeMultipleTrsfs
	      meanTrsfs
	      multipleTrsfs
	      statsMultipleTrsfs
	      trsf2xxxlab
)
  
## #################################################################
## Build
## #################################################################

# build library and link
add_library(${LIB_NAME} ${LIB_STATUS} ${SRC_FILES})
target_link_libraries( ${LIB_NAME} blockmatchingLIB  ${ZLIB_LIBRARIES} basic io)

# build execs and link
foreach(E ${EXE_NAMES})
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} ${LIB_NAME} blockmatchingLIB basic io)
endforeach(E)

