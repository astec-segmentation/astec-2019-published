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
SET( BASIC_INCLUDE_DIRS   ${CMAKE_CURRENT_SOURCE_DIR}/libbasic )
SET( IO_INCLUDE_DIRS      ${CMAKE_CURRENT_SOURCE_DIR}/libio )
SET( VT_INCLUDE_DIRS      ${CMAKE_CURRENT_SOURCE_DIR}/libvt ) 
SET( FILTERS_INCLUDE_DIRS      ${CMAKE_CURRENT_SOURCE_DIR}/libfilters ) 
SET( BLOCKMATCHING_INCLUDE_DIRS      ${CMAKE_CURRENT_SOURCE_DIR}/blockmatching ) 
SET( NOESIS_INCLUDE_DIRS  ${CMAKE_CURRENT_SOURCE_DIR}/libnoesis )
SET( DAVID_INCLUDE_DIRS   ${CMAKE_CURRENT_SOURCE_DIR}/libdavid )

# Independants
add_subdirectory(libio)
add_subdirectory(libbasic)

# Depend on BASIC and IO
add_subdirectory(blockmatching)
add_subdirectory(libvt)

# Depend on BASIC, IO and BLOCKMATCHING
add_subdirectory(libmatchingtools)

# Depend on BASIC, IO, and VT
add_subdirectory(libtools)
add_subdirectory(libdavid)
add_subdirectory(libmembrane)
add_subdirectory(libnoesis)
add_subdirectory(libfilters)
add_subdirectory(vt-exec)


# Depend on BASIC, IO, VT, and DAVID
add_subdirectory(libvp)







