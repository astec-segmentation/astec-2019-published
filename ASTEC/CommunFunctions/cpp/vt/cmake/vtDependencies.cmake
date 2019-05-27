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

## #################################################################
## ZLIB
## #################################################################
# Look for Z lib. Automatically defines following values : 
# ZLIB_FOUND         : Bool : true if find, else false
# ZLIB_INCLUDE_DIRS  : ZLIB HEADERS
# ZLIB_LIBRARIES     : FILE FOR LINKING WITH ZLIB
find_package( ZLIB REQUIRED )
include_directories( ${ZLIB_INCLUDE_DIRS} )

if(vt_USE_OPENMP)
  ## #################################################################
  # Look for OpenMP. Automatically definecs following values : 
  # OPENMP_FOUND     : Bool : true if find, else false
  # OpenMP_C_FLAGS   : OpenMP flags for C compiler
  # OpenMP_CXX_FLAGS : OpenMP flags for CXX compiler
  find_package( OpenMP )
  if (OPENMP_FOUND)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

## #################################################################
## M Library
## #################################################################
# Unix specific : link the m library (automated on windows)
if(UNIX)
    link_libraries(m)
    link_libraries(pthread)
endif(UNIX)

## #################################################################
## LEMON
## #################################################################

# find_package()
# REQUIRED option stops processing with an error message if the package cannot be found.
# EXACT option requests that the version be matched exactly
#

#
# message: Display a message to the user.
#
#  message([STATUS|WARNING|AUTHOR_WARNING|FATAL_ERROR|SEND_ERROR]
#          "message to display" ...)
#
# The optional keyword determines the type of message:
#
#  (none)         = Important information
#  STATUS         = Incidental information
#  WARNING        = CMake Warning, continue processing
#  AUTHOR_WARNING = CMake Warning (dev), continue processing
#  SEND_ERROR     = CMake Error, continue but skip generation
#  FATAL_ERROR    = CMake Error, stop all processing
#

find_package (LEMON QUIET)
if ( LEMON_FOUND )
  message(STATUS "lemon was found")
  include_directories(${LEMON_INCLUDE_DIRS})
  link_directories(${LEMON_LIBRARY_DIR})
#  message( "" )
#  message( "LEMON_LIBRARIES=${LEMON_LIBRARIES}" )
#  message( "LEMON_LIBRARY_DIR=${LEMON_LIBRARY_DIR}" )
#  message( "LEMON_INCLUDE_DIRS=${LEMON_INCLUDE_DIRS}" )
#  message( "LEMON_USE_FILE=${LEMON_USE_FILE}" )
#  message( "" )
else( LEMON_FOUND )
  message(STATUS "lemon was NOT found")
endif( LEMON_FOUND )


## #################################################################
## VTK
## #################################################################

find_package (VTK QUIET)
if ( VTK_FOUND )
  if ( ${VTK_MAJOR_VERSION} LESS 6 )
    message(WARNING "\n *** VTK version must be 6 or greater *** \n")
    set( VTK_FOUND 0 )
  else( ${VTK_MAJOR_VERSION} LESS 6 )
    message(STATUS "vtk was found")
    include(${VTK_USE_FILE}) 
    # include_directories(${VTK_INCLUDE_DIRS})
    # link_directories(${VTK_LIBRARY_DIR})
#    message( "" )
#    message( "VTK_LIBRARIES=${VTK_LIBRARIES}" )
#    message( "VTK_LIBRARY_DIRS=${VTK_LIBRARY_DIRS}" )
#    message( "VTK_INCLUDE_DIRS=${VTK_INCLUDE_DIRS}" )
#    message( "VTK_USE_FILE=${VTK_USE_FILE}" )
#    message( "" )
  endif( ${VTK_MAJOR_VERSION} LESS 6 )
else( VTK_FOUND )
  message(STATUS "vtk was NOT found")
endif( VTK_FOUND )

## #################################################################
## NIFTI
## #################################################################

find_package (NIFTI QUIET)
if ( NIFTI_FOUND )
  message(STATUS "nifti was found")
  include(${NIFTI_USE_FILE}) 
  include_directories(${NIFTI_INCLUDE_DIRS})
  link_directories(${NIFTI_LIBRARY_DIRS})
  link_libraries(niftiio znz)
#  message( "" )
#  message( "NIFTI_LIBRARIES=${NIFTI_LIBRARIES}" )
#  message( "NIFTI_LIBRARY_DIRS=${NIFTI_LIBRARY_DIRS}" )
#  message( "NIFTI_INCLUDE_DIRS=${NIFTI_INCLUDE_DIRS}" )
#  message( "NIFTI_USE_FILE=${NIFTI_USE_FILE}" )
#  message( "" )
  SET( CPPFLAGS "${CPPFLAGS} -D_NIFTICLIB_" )
  add_definitions( ${CPPFLAGS} )
#  MESSAGE( "CPPFLAGS=${CPPFLAGS}" )
else( NIFTI_FOUND )
  message(STATUS "nifti was NOT found")
endif( NIFTI_FOUND )



## #################################################################
## KLB
## #################################################################

find_package (KLB QUIET)
if ( KLB_FOUND )
  message(STATUS "klb was found")
  include(${KLB_USE_FILE}) 
  include_directories(${KLB_INCLUDE_DIRS})
  link_directories(${KLB_LIBRARY_DIRS})
  link_libraries(klb)
#  link_libraries(klb z bzip2)
#  message( "" )
#  message( "KLB_LIBRARIES=${KLB_LIBRARIES}" )
#  message( "KLB_LIBRARY_DIRS=${KLB_LIBRARY_DIRS}" )
#  message( "KLB_INCLUDE_DIRS=${KLB_INCLUDE_DIRS}" )
#  message( "KLB_USE_FILE=${KLB_USE_FILE}" )
#  message( "" )
  SET( CPPFLAGS "${CPPFLAGS} -D_KLB_" )
  add_definitions( ${CPPFLAGS} )
#  MESSAGE( "CPPFLAGS=${CPPFLAGS}" )
else( KLB_FOUND )
  message(STATUS "klb was NOT found")
endif( KLB_FOUND )



## #################################################################
## LIBTIFF
## #################################################################

find_package (LIBTIFF QUIET)
if ( LIBTIFF_FOUND )
  message(STATUS "libtiff was found")
  include(${LIBTIFF_USE_FILE}) 
  include_directories(${LIBTIFF_INCLUDE_DIRS})
  link_directories(${LIBTIFF_LIBRARY_DIRS})
  link_libraries(tiff)
#  message( "" )
#  message( "LIBTIFF_LIBRARIES=${LIBTIFF_LIBRARIES}" )
#  message( "LIBTIFF_LIBRARY_DIRS=${LIBTIFF_LIBRARY_DIRS}" )
#  message( "LIBTIFF_INCLUDE_DIRS=${LIBTIFF_INCLUDE_DIRS}" )
#  message( "LIBTIFF_USE_FILE=${LIBTIFF_USE_FILE}" )
#  message( "" )
  SET( CPPFLAGS "${CPPFLAGS} -D_LIBLIBTIFF_" )
  add_definitions( ${CPPFLAGS} )
#  MESSAGE( "CPPFLAGS=${CPPFLAGS}" )
else( LIBTIFF_FOUND )
  message(STATUS "libtiff was NOT found (from cmake/vtDependencies)")
endif( LIBTIFF_FOUND )









