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
project(libvt)
set(LIB_NAME vt)

## #################################################################
## Input Files 
## #################################################################
# Inludes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})

# Src files
SET(FILES vt_3m.c
	vt_ambdd.c
	vt_amincir.c
	vt_amliste.c
	vt_amseuil.c
	vt_average.c
	vt_bbox.c
	vt_bdd.c
	vt_bdd2.c
	vt_bdd3.c
	vt_bdd_amincir.c
	vt_boolean.c
	vt_bytes.c
	vt_caracterisation.c
	vt_cc.c
	vt_connexe.c
	vt_contours.c
	vt_copy.c
	vt_daneucmapsc.c
	vt_daneucmapss.c
	vt_distance.c
	vt_distfield.c
	vt_dst_amincir.c
	vt_error.c
	vt_etoilebarre.c
	vt_eucmap.c
	vt_eucmapsc.c
	vt_eucmapss.c
	vt_extract.c
	vt_filters.c
	vt_gausssmooth.c
	vt_gb.c
	vt_getval.c
	vt_greynumbers.c
	vt_histo-tools.c
	vt_histo.c
	vt_image.c
	vt_inrimage.c
	vt_localsmooth.c
	vt_maxgrad.c
	vt_median.c
	vt_mescorr.c
	vt_mip.c
	vt_names.c
	vt_old_amincir.c
	vt_old_t26.c
	vt_old_t26noopt.c
	vt_recfilters.c
	vt_recline.c
	vt_seuil.c
	vt_skiz.c
	vt_slconnexe.c
	vt_t04t08.c
	vt_t06.c
	vt_t06noopt.c
	vt_t06t26.c
	vt_t06t26_simple.c
	vt_t26.c
	vt_t26noopt.c
	vt_unix.c
	vt_zerocross.c
	vt_morpho.c
	vt_reech.c
#	vt_histo-seuil.c
#	vt_stats_bdd2.c
#	vt_stats_t06.c 
#	vt_stats_t06noopt.c 
#	vt_stats_t06t26.c 
#	vt_stats_t06t26_simple.c 
#	vt_stats_t26.c
#	vt_stats_t26noopt.c
  )


## #################################################################
## Build 
## #################################################################
# Build library
add_library(${LIB_NAME} ${FILES})
target_link_libraries(${LIB_NAME} ${ZLIB_LIBRARIES} basic io)


