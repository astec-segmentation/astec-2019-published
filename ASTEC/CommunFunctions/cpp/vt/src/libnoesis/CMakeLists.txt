############################################################
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
project(libnoesis)
SET(LIB_NAME noesis)

## #################################################################
## Input Files 
## #################################################################
# Inludes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})
include_directories(${VT_INCLUDE_DIRS})

# Src files
SET(SRC_FILES chamfer.c
#	convert.c
	ccparameters.c
)

## #################################################################
## Build Directories 
## #################################################################
# Build libs and link
add_library(${LIB_NAME} ${SRC_FILES})
target_link_libraries( ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt )





