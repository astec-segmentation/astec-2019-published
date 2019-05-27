## #################################################################
## Configure cmake variables
## #################################################################

set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type")

set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

mark_as_advanced(CMAKE_BACKWARDS_COMPATIBILITY)
mark_as_advanced(CMAKE_BUILD_TYPE)
mark_as_advanced(CMAKE_INSTALL_PREFIX)
if(APPLE)
  mark_as_advanced(CMAKE_OSX_ARCHITECTURES)
  mark_as_advanced(CMAKE_OSX_DEPLOYMENT_TARGET)
  mark_as_advanced(CMAKE_OSX_SYSROOT)
endif()
mark_as_advanced(BUILD_TESTING)

## #################################################################
## Install prefix
## #################################################################

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install" CACHE PATH "${PROJECT_NAME} install prefix" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

## #################################################################
## Configure arch
## #################################################################

if(CMAKE_SIZEOF_VOID_P MATCHES 8)
    set(DTK_PLATFORM_64 1)
else()
    set(DTK_PLATFORM_32 1)
endif()

## #################################################################
## Configure path
## #################################################################

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_BINARY_DIR})

set(${PROJECT_NAME}_ARCHIVE_OUTPUT_DIRECTORY lib)
set(${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY bin)
set(${PROJECT_NAME}_LIBRARY_OUTPUT_DIRECTORY lib)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${${PROJECT_NAME}_LIBRARY_OUTPUT_DIRECTORY})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${${PROJECT_NAME}_ARCHIVE_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY})

set(LIBRARY_INSTALL_OUTPUT_PATH ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_LIBRARY_OUTPUT_DIRECTORY})
set(ARCHIVE_INSTALL_OUTPUT_PATH ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_ARCHIVE_OUTPUT_DIRECTORY})
set(RUNTIME_INSTALL_OUTPUT_PATH ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY})
set(EXECUTABLE_INSTALL_OUTPUT_PATH ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY})

set(${PROJECT_NAME}_INCLUDE_DIR ${${PROJECT_NAME}_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}) 
set(${PROJECT_NAME}_LIBRARY_DIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(${PROJECT_NAME}_RUNTIME_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set(${PROJECT_NAME}_ARCHIVE_DIR ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
set(${PROJECT_NAME}_CMAKE_DIR ${PROJECT_SOURCE_DIR}/cmake)
set(${PROJECT_NAME}_USE_FILE ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Use.cmake)

set(${PROJECT_NAME}_INSTALL_INCLUDE_DIR ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR} ${CMAKE_INSTALL_PREFIX}/include ${CMAKE_INSTALL_PREFIX}/include/dtk)
set(${PROJECT_NAME}_INSTALL_LIBRARY_DIR ${LIBRARY_INSTALL_OUTPUT_PATH})
set(${PROJECT_NAME}_INSTALL_RUNTIME_DIR ${RUNTIME_INSTALL_OUTPUT_PATH})
set(${PROJECT_NAME}_INSTALL_CMAKE_DIR ${CMAKE_INSTALL_PREFIX}/cmake)
set(${PROJECT_NAME}_INSTALL_USE_FILE ${CMAKE_INSTALL_PREFIX}/cmake/${PROJECT_NAME}Use.cmake)

## #################################################################
## Setup configration files
## #################################################################

if(EXISTS ${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in)
configure_file( ## Build tree configure file
  ${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
  ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  @ONLY IMMEDIATE)
endif()

if(EXISTS ${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.install.cmake.in)
configure_file( ## Install tree configure file
  ${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.install.cmake.in
  ${PROJECT_BINARY_DIR}/install/${PROJECT_NAME}Config.cmake
  @ONLY IMMEDIATE)
endif()

## #################################################################
## Setup find file
## #################################################################

string(TOUPPER ${PROJECT_NAME} FND_NAME)

if(EXISTS ${PROJECT_SOURCE_DIR}/cmake/Find${FND_NAME}.cmake.in)
configure_file( ## Build tree configure file
  ${PROJECT_SOURCE_DIR}/cmake/Find${FND_NAME}.cmake.in
  ${PROJECT_BINARY_DIR}/CMakeModules/Find${FND_NAME}.cmake
  @ONLY IMMEDIATE)
endif()

## #################################################################
## Setup use file
## #################################################################

if(EXISTS ${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Use.cmake.in)
configure_file( ## Common use file
  ${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Use.cmake.in
  ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Use.cmake
  @ONLY IMMEDIATE)
endif()

## #################################################################
## Path relocation
## #################################################################

if (APPLE)
  set(CMAKE_INSTALL_NAME_DIR "${LIBRARY_INSTALL_OUTPUT_PATH}")
endif ()

if (UNIX AND NOT APPLE)
  set(CMAKE_INSTALL_RPATH "${LIBRARY_INSTALL_OUTPUT_PATH}")
endif ()

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

## #################################################################
## Uninstall target
## #################################################################

if(EXISTS ${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Uninstall.cmake.in)

configure_file("${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Uninstall.cmake.in"
               "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Uninstall.cmake"
               IMMEDIATE @ONLY)

add_custom_target(${PROJECT_NAME}-uninstall
  "${CMAKE_COMMAND}" -P
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Uninstall.cmake")

endif()

## #################################################################
## Install cmake files
## #################################################################

if( EXISTS ${${PROJECT_NAME}_SOURCE_DIR}/cmake/${PROJECT_NAME}Dependencies.cmake
AND EXISTS ${${PROJECT_NAME}_SOURCE_DIR}/cmake/${PROJECT_NAME}Build.cmake
AND EXISTS ${${PROJECT_NAME}_SOURCE_DIR}/cmake/${PROJECT_NAME}Pack.cmake
AND EXISTS ${${PROJECT_NAME}_SOURCE_DIR}/cmake/${PROJECT_NAME}Reporting.cmake
AND EXISTS ${${PROJECT_NAME}_BINARY_DIR}/install/${PROJECT_NAME}Config.cmake
AND EXISTS ${${PROJECT_NAME}_BINARY_DIR}/${PROJECT_NAME}Use.cmake
AND EXISTS ${${PROJECT_NAME}_BINARY_DIR}/${PROJECT_NAME}Uninstall.cmake)

install(FILES
  ${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Dependencies.cmake
  ${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Build.cmake
  ${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Pack.cmake
  ${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Reporting.cmake
  ${PROJECT_BINARY_DIR}/install/${PROJECT_NAME}Config.cmake
  ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Use.cmake
  ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Uninstall.cmake
  DESTINATION
  cmake)

endif()

if(EXISTS ${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Plugin.cmake)
install(FILES
  ${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Plugin.cmake
  DESTINATION
  cmake)
endif()
