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
project(libdavid)
SET(LIB_NAME david)

## #################################################################
## Input Files 
## #################################################################
# Inludes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})
include_directories(${VT_INCLUDE_DIRS})

# Src files
SET(SRC_FILE ccparameters.c)

# Exe names
SET(EXE_NAMES 
	ext-param-2
	ext-param
)

## #################################################################
## Build Directories 
## #################################################################

# Build lib and link
add_library(${LIB_NAME} ${SRC_FILE})   
target_link_libraries( ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt )

# Build binaries and link
foreach(E ${EXE_NAMES})
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt)
endforeach(E)






