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
project(vt-exec)
set(LIB_NAME exec)

## #################################################################
## Input Files 
## #################################################################
# Includes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})
include_directories(${VT_INCLUDE_DIRS})

# Src files
SET(SRC_FILES 
       api-connexe.c
#       api-execTemplate.c
       api-linearFilter.c
       api-meanImages.c
       api-morpho.c
       api-regionalext.c
       api-tree.c
       api-watershed.c
)

# Exe names
SET(EXE_NAMES 
        Arit  
	Convert  
	Create  
	Ctable  
	Ext  
	Ifmt  
	Logic  
	Tran  
	amincir  
	axeMedianOfComponents  
	bbox  
	binary_noise  
	caracterisation  
	chamfer  
	changevals  
	connexe  
	contours  
	convolution  
	copy  
	createImage  
	detrend  
	distance  
	drawShapes
 #     execTemplate
        extImage  
	filters  
	getVoxelValue
	gradient_unitaire  
	grey2color  
	histo  
	histoCopy
	histoJoint
	histoMerge
	histoStats
	laplacien  
	linearFilter  
	localFilter
	meanImages  
	melImages
	mergeSlices
	mescorr  
#	minimum  
	mip_project  
	morpho  
	newTran  
	new_skiz  
	normalisation  
	norme_gradient  
	parcellisation  
        patchLogic
	recfilters  
	reech3d  
	regionalext  
#	regionalmax  
	relabel  
	sbinter  
	setVoxelSize  
	setVoxelValue  
	seuillage    
	slow_connexe  
	smooth  
	splitImage
	surface06  
	tgrowing
	thinning  
        tree
#        tsingle-axon
	watershed  
	zpar
)

# For exe name with suffixes needed: presents in another project
SET(EXE_NAMES_DUPLICATED 
	skeleton # libtopo
	skiz # libelf
)

message(STATUS "Bin names duplicated : skeleton (libtopo), skiz (libelf)") 
# message(WARNING "Bin names duplicated : skeleton (libtopo), skiz (libelf)") 


## #################################################################
## Build Directories 
## #################################################################

add_library(${LIB_NAME} ${SRC_FILES})
target_link_libraries( ${LIB_NAME} basic io vt)

foreach(E ${EXE_NAMES})
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} ${LIB_NAME} basic io vt)
endforeach(E)

foreach(E ${EXE_NAMES_DUPLICATED})
  add_executable(vt-exec_${E} ${E}.c)
  target_link_libraries(vt-exec_${E} ${LIB_NAME} basic io vt)
endforeach(E)
