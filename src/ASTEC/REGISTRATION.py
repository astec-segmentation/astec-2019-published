import numpy as np
import sys
import os
from scipy import ndimage as nd
import pickle as pkl

sys.path.append(os.path.join(os.path.dirname(__file__), "CommunFunctions"))
from ImageHandling import imread, imsave, SpatialImage
from cpp_wrapping import *

# move into commonTools
def readLUT(file):
    '''
	Return a dictionnary of integer key-to-key correspondances
	'''
    lut = {}
    if os.path.exists(file):
        f = open(file)
        for line in f:
            li = line.strip()
            if not li.startswith('#'):
                info = li.split()
                if len(info) == 2:
                    if not lut.has_key(int(info[0])):
                        lut[int(info[0])] = None
                    lut[int(info[0])] = int(info[1])
    return lut


def symmetry_plane(image_bin, image_seg=None, image_direction_histogram="tmp_sphere.inr",
                   equation_output=None, trsf_output=None, plane_output=None,
                   frac=None, maxima=None,
                   d=10, dmin=2, delta=2, p=0.1, sigma=None, iterations=None, realSize=True, keep_all=False,
                   verbose=False):
    """


	maxima: pour specifier le ou les maxima (0=plus gros, 1=deuxieme plus gros, etc.) que l'on souhaite utiliser pour l'extraction du plan de symetrie.
			Si non specifie, on prend par defaut la liste des maxima de hauteur >= hauteur max * frac
	frac: pour specifier la tolerance de hauteur que l'on se donne pour l'extraction des maxima en fonction de la hauteur maximale de l'histogramme des directions.
			Si non specifie, le programme par defaut prend la valeur 0.5 pour ce parametre.

	"""
    ###### Example of pipeline of shell commands #######
    # symmetryPlane bin_t0${i}_on_t099.inr -sphere HISTO/bin_t0${i}_on_t099_R15A32.inr -sigma 0.1 -weq SYM/planeEq_t0${i}_on_t099_max_0.txt -d 0 -dmin 1
    #
    # diceMaximisation WAT/OUT_t0${i}_on_t099-wat.inr SYM/planeEq_t0${i}_on_t099_max_0_dmax.txt -n `cat SYM/planeEq_t0${i}_on_t099_max_0.txt | grep -v "#" | awk -F':' '{print $2}'` -delta 10
    #
	# symmetryPlane BIN/bin_t0${i}_on_t099.inr -n `cat SYM/planeEq_t0${i}_on_t099_max_0_dmax.txt | grep new | awk -F':' '{print $2}'` -weq SYM/planeEq_t0${i}_on_t099_max_0_dmax_a.txt -d 10 -dmin 10 -p 0.1

	# assert(os.path.exists(flo_file_bin) and os.path.exists(flo_file_hist))
    # symmetryPlane(flo_file_bin, flo_file_hist, equation_output=flo_file_sym_eq, trsf_output=flo_file_sym_alignment_trsf, plane_output=flo_file_sym_plane,
    #			  maximum=maximum_ref, d=d_ref, dmin=dmin_ref, delta=delta_ref, p=p_ref, sigma=sigma_ref, realSize=True, verbose=verbose)

    # Histogramme directionnel
    if not os.path.exists(image_direction_histogram):
        assert (os.path.exists(image_bin))
        if not os.path.dirname(image_direction_histogram):
            try:
                os.makedirs(os.path.dirname(image_direction_histogram))
            except:
                print "Unable to create experiences path. Check the complete path..."
                raise
        directionHistogram(image_bin, image_direction_histogram, verbose=verbose)

    # Extraction de la liste des directions candidates

    if image_seg:
        # Histogramme maxima extraction
        candidates = directionHistogramMaxima(image_direction_histogram, file_out=None, maxima=maxima, frac=frac,
                                              verbose=False)

        n = candidates.shape[0]
        means_max = []
        i_max = 0
        tmp_eq_files = []
        tmp_plane_files = []
        tmp_distribution_files = []
        tmp_trsf_files = []
        # boucle sur les candidats
        for i in range(n):
            tmp_eq_files[i] = 'tmp_symmetryPlane_' + str(i) + '.eq'
            tmp_plane_files[i] = None
            tmp_distribution_files[i] = None
            tmp_trsf_files[i] = None
            if plane_output:
                tmp_plane_files[i] = 'tmp_symmetryPlane_' + str(i) + '.inr'
            if distribution_output:
                tmp_distribution_files[i] = 'tmp_symmetryPlane_' + str(i) + '.distribution'
            if trsf_output:
                tmp_trsf_files[i] = 'tmp_symmetryPlane_' + str(i) + '.trsf'
            symmetryPlane(image_bin, path_input_sphere=None, normal=candidates[i], equation_output=tmp_eq_files[i],
                          plane_output=tmp_plane_files[i], distribution_output=tmp_distribution_files[i],
                          trsf_output=tmp_trsf_files[i],
                          d=d, dmin=dmin, delta=delta, p=p, sigma=sigma, realSize=realSize, iterations=iterations,
                          verbose=verbose, lazy=True)
            d, m, s = dice(image_seg, symmetry=tmp_eq_files[i])
            if means_max == [] or m > max(means_max):
                i_max = i
            means_max[i] = m
        # Copie des fichiers correspondants au plan de symetrie du candidat optimal
        cmd = "cp " + tmp_eq_files[i_max] + " " + 'tmp_symmetryPlane.eq'
        if verbose:
            print cmd
        os.system(cmd)
        if equation_output:
            cmd = "cp " + tmp_eq_files[i_max] + " " + equation_output
            if verbose:
                print cmd
            os.system(cmd)
        if plane_output:
            copy(tmp_plane_files[i_max], plane_output, verbose=verbose)
        if distribution_output:
            cmd = "cp " + tmp_distribution_files[i_max] + " " + distribution_output
            if verbose:
                print cmd
            os.system(cmd)
        if trsf_output:
            cmd = "cp " + tmp_trsf_files[i_max] + " " + trsf_output
            if verbose:
                print cmd
            os.system(cmd)
        # Effacement des donnees temporaires
        if not keep_all:
            for i in range(n):
                cmd = "rm " + tmp_eq_files[i]
                if tmp_plane_files[i]:
                    cmd += " " + tmp_plane_files[i]
				if tmp_distribution_files[i]:
					cmd += " " + tmp_distribution_files[i]
                if tmp_trsf_files[i]:
                    cmd += " " + tmp_trsf_files[i]
            if verbose:
                print cmd
            os.system(cmd)
    else:
        maximum = None
        if len(maxima) == 1 and type(maxima[0]) == int:
            maximum = maxima[0]
        symmetryPlane(image_bin, path_input_sphere=image_direction_histogram, normal=None,
                      equation_output='tmp_symmetryPlane.eq', plane_output=plane_output,
                      distribution_output=distribution_output, trsf_output=trsf_output,
                      maximum=maximum, d=d, dmin=dmin, delta=delta, p=p, sigma=sigma, realSize=realSize,
                      iterations=iterations, verbose=verbose, lazy=True)
        # Copie du fichier correspondants a l'equation du plan de symetrie
        if equation_output:
            cmd = "cp tmp_symmetryPlane.eq " + equation_output
            if verbose:
                print cmd
            os.system(cmd)
    # Lecture de l'equation du plan de symetrie
    f = open("tmp_symmetryPlane.eq")
    n = []
    for line in f:
        li = line.strip()
        # print li
        if not li.startswith('#'):
            # print li
            l = li.split()
            for value in l[1:5]:
                n.append(float(value))
    # Effacement du fichier d'equation temporaire
    cmd = "rm tmp_symmetryPlane.eq"
    if verbose:
        print cmd
    os.system(cmd)
    return n


def spatial_registration(ref_seg_post_file, flo_seg_post_file, ref_fused_file=None, flo_fused_file=None,  # Input images
                         path_trsf_flo_ref=None, path_pairs_ref_flo=None, path_dices_ref_flo=None,
                         path_residuals_ref_flo=None,  # Outputs registration
                         embryoKey_ref=None, embryoKey_flo=None, folder_tmp='WORKSPACE/',  # Temporary files and paths
                         init_ref=0.9, init_flo=0.9, realScale_ref=True, realScale_flo=True,
                         # arguments for membrane_renforcement ref and flo
                         sensitivity_ref=0.99, sensitivity_flo=0.99,  # arguments for anisotropicHist ref and flo
                         maximum_ref=None, d_ref=None, dmin_ref=None, delta_ref=None, p_ref=None, sigma_ref=None,
                         # arguments for symmetryPlane ref
                         maximum_flo=None, d_flo=None, dmin_flo=None, delta_flo=None, p_flo=None, sigma_flo=None,
                         # arguments for symmetryPlane flo
                         trsf_type='rigid', estimator='lts', lts_fraction=0.9,  # arguments for planeRegistration
                         background_ref=1, background_flo=1,  # Background labels for planeRegistration
                         keep_mem=False,
                         keep_bin=False,
                         keep_hist=False,
                         keep_sym=False,
                         keep_inter=False,
                         verbose=False):
    '''
	# Temporary folders:

	folder_mem = folder_tmp + 'mem/'
	folder_bin = folder_tmp + 'bin/'
	folder_hist = folder_tmp + 'hist/'
	folder_sym = folder_tmp + 'sym/'
	folder_inter = folder_tmp + 'inter/'

	# Keep temporary folders ?

	keep_mem = False,
	keep_bin = False,
	keep_hist = False,
	keep_sym = False,
	keep_inter = False,

	# Parameters:

	ref_seg_post_file, flo_seg_post_file, ref_fused_file=None, flo_fused_file=None,  				# Input images
	path_trsf_flo_ref=None, path_pairs_ref_flo=None, path_dices_ref_flo=None, path_residuals_ref_flo=None, 		# Outputs registration
	embryoKey_ref=None, embryoKey_flo=None, folder_tmp='WORKSPACE/', 			# Temporary files and paths
	init_ref=0.9, init_flo=0.9, realScale_ref=True, realScale_flo=True, 		# arguments for membrane_renforcement ref and flo
	sensitivity_ref=0.99, sensitivity_flo=0.99, 		# arguments for anisotropicHist ref and flo
	maximum_ref=None, d_ref=None, dmin_ref=None, delta_ref=None, p_ref=None, sigma_ref=None, 		# arguments for symmetryPlane ref
	maximum_flo=None, d_flo=None, dmin_flo=None, delta_flo=None, p_flo=None, sigma_flo=None, 		# arguments for symmetryPlane flo
	trsf_type='rigid', estimator='lts', lts_fraction=0.9,		 				# arguments for planeRegistration
	background_ref = 1, background_flo = 1, 			# Background labels for planeRegistration

	# Verbosity:

	verbose=False
	'''

    delete_workspace = not (keep_mem or keep_bin or keep_hist or keep_sym or keep_inter)

    ### Ref ###

    if not embryoKey_ref:
        embryoKey_ref = ref_seg_post_file.split(os.path.sep)[-1]
        embryoKey_ref = embryoKey_ref.split('-')[1] + '_' + embryoKey_ref.split('_')[-1].split('.')[0]

    ### Flo ###

    if not embryoKey_flo:
        embryoKey_flo = flo_seg_post_file.split(os.path.sep)[-1]
        embryoKey_flo = embryoKey_flo.split('-')[1] + '_' + embryoKey_flo.split('_')[-1].split('.')[0]

    ### TEMPORARY PATH VARIABLES ###

    folder_mem = folder_tmp + 'mem/'
    folder_bin = folder_tmp + 'bin/'
    folder_hist = folder_tmp + 'hist/'
    folder_sym = folder_tmp + 'sym/'
    folder_inter = folder_tmp + 'inter/'

    ### Flo ###

    flo_file_mem_prefix = folder_mem + embryoKey_flo + "_mem"
    flo_file_bin_prefix = folder_bin + embryoKey_flo + "_bin"
    flo_file_bin = flo_file_bin_prefix + '.inr'
    flo_file_theta = flo_file_bin_prefix + '.theta.inr'
    flo_file_phi = flo_file_bin_prefix + '.phi.inr'
    flo_file_hist = folder_hist + embryoKey_flo + "_directionHistogram.inr"
    flo_file_sym_prefix = folder_sym + embryoKey_flo + "_sym"
    flo_file_sym_eq = flo_file_sym_prefix + '.eq'
    flo_file_sym_plane = flo_file_sym_prefix + '.plane.inr'
    flo_file_sym_alignment_trsf = flo_file_sym_prefix + '.trsf'

    ### Ref ###

    ref_file_mem_prefix = folder_mem + embryoKey_ref + "_mem"
    ref_file_bin_prefix = folder_bin + embryoKey_ref + "_bin"
    ref_file_bin = ref_file_bin_prefix + '.inr'
    ref_file_theta = ref_file_bin_prefix + '.theta.inr'
    ref_file_phi = ref_file_bin_prefix + '.phi.inr'
    ref_file_hist = folder_hist + embryoKey_ref + "_directionHistogram.inr"
    ref_file_sym_prefix = folder_sym + embryoKey_ref + "_sym"
    ref_file_sym_eq = ref_file_sym_prefix + '.eq'
    ref_file_sym_plane = ref_file_sym_prefix + '.plane.inr'
    ref_file_sym_alignment_trsf = ref_file_sym_prefix + '.trsf'

    ### Registration ###

    ref_flo_file_trsf = folder_inter + embryoKey_flo + '_' + embryoKey_ref + '.trsf'
    ref_flo_file_pairs = folder_inter + embryoKey_ref + '_' + embryoKey_flo + '.pairs'
    ref_flo_file_res = folder_inter + embryoKey_ref + '_' + embryoKey_flo + '_residuals.m'
    ref_flo_file_dices = folder_inter + embryoKey_ref + '_' + embryoKey_flo + '.dices'

    ### STUFF ###

    if not os.path.isdir(folder_tmp):
        try:
            os.makedirs(folder_tmp)
        except:
            print "Unable to create experiences path. Check the complete path..."
            raise

    # Fonction de rehaussement de membranes :
    if not (os.path.exists(ref_file_bin) or os.path.exists(ref_file_mem_prefix + '.ext.inr')):
        assert (os.path.exists(ref_fused_file))
        if not os.path.isdir(folder_mem):
            try:
                os.makedirs(folder_mem)
            except:
                print "Unable to create experiences path. Check the complete path..."
                raise
        membrane_renforcement(ref_fused_file, prefix_output=ref_file_mem_prefix, init=init_ref, realScale=realScale_ref,
                              verbose=verbose)

    if not (os.path.exists(flo_file_bin) or os.path.exists(flo_file_mem_prefix + '.ext.inr')):
        assert (os.path.exists(flo_fused_file))
        if not os.path.isdir(folder_mem):
            try:
                os.makedirs(folder_mem)
            except:
                print "Unable to create experiences path. Check the complete path..."
                raise
        membrane_renforcement(flo_fused_file, prefix_output=flo_file_mem_prefix, init=init_flo, realScale=realScale_flo,
                              verbose=verbose)

    # Binarisation de membranes :
    if not os.path.exists(ref_file_bin):
        assert (os.path.exists(ref_file_mem_prefix + '.ext.inr'))
        if not os.path.isdir(folder_bin):
            try:
                os.makedirs(folder_bin)
            except:
                print "Unable to create experiences path. Check the complete path..."
                raise
        anisotropicHist(ref_file_mem_prefix + '.ext.inr', path_output=ref_file_bin, sensitivity=sensitivity_ref,
                        verbose=verbose)

    if not os.path.exists(flo_file_bin):
        assert (os.path.exists(flo_file_mem_prefix + '.ext.inr'))
        if not os.path.isdir(folder_bin):
            try:
                os.makedirs(folder_bin)
            except:
                print "Unable to create experiences path. Check the complete path..."
                raise
        anisotropicHist(flo_file_mem_prefix + '.ext.inr', path_output=flo_file_bin, sensitivity=sensitivity_flo,
                        verbose=verbose)

    if not keep_mem:
        cmd = 'rm -rf ' + folder_mem
        if verbose:
            print cmd
        os.system(cmd)

    # Calcul de l'histogramme directionnel
    if not os.path.exists(ref_file_hist):
        assert (os.path.exists(ref_file_bin))
        if not os.path.isdir(folder_hist):
            try:
                os.makedirs(folder_hist)
            except:
                print "Unable to create experiences path. Check the complete path..."
                raise
        directionHistogram(ref_file_bin, ref_file_hist, verbose=verbose)
    if not os.path.exists(flo_file_hist):
        assert (os.path.exists(flo_file_bin))
        if not os.path.isdir(folder_hist):
            try:
                os.makedirs(folder_hist)
            except:
                print "Unable to create experiences path. Check the complete path..."
                raise
        directionHistogram(flo_file_bin, flo_file_hist, verbose=verbose)

    # Calcul de l'equation du plan de symetrie
    assert (os.path.exists(ref_file_bin))
    assert (os.path.exists(flo_file_bin))
    assert (os.path.exists(ref_file_hist))
    assert (os.path.exists(flo_file_hist))
    if not os.path.isdir(folder_sym):
        try:
            os.makedirs(folder_sym)
        except:
            print "Unable to create experiences path. Check the complete path..."
            raise
    if not os.path.exists(flo_file_sym_eq):
        assert (os.path.exists(flo_file_bin) and os.path.exists(flo_file_hist))
        symmetryPlane(flo_file_bin, flo_file_hist, equation_output=flo_file_sym_eq,
                      trsf_output=flo_file_sym_alignment_trsf, plane_output=flo_file_sym_plane,
                      maximum=maximum_ref, d=d_ref, dmin=dmin_ref, delta=delta_ref, p=p_ref, sigma=sigma_ref,
                      realSize=True, verbose=verbose)
    if not os.path.exists(ref_file_sym_eq):
        assert (os.path.exists(ref_file_bin) and os.path.exists(ref_file_hist))
        symmetryPlane(ref_file_bin, ref_file_hist, equation_output=ref_file_sym_eq,
                      trsf_output=ref_file_sym_alignment_trsf, plane_output=ref_file_sym_plane,
                      maximum=maximum_flo, d=d_flo, dmin=dmin_flo, delta=delta_flo, p=p_flo, sigma=sigma_flo,
                      realSize=True, verbose=verbose)

    if not keep_bin:
        cmd = 'rm -rf ' + folder_bin
        if verbose:
            print cmd
        os.system(cmd)

    if not keep_hist:
        cmd = 'rm -rf ' + folder_hist
        if verbose:
            print cmd
        os.system(cmd)

    # Registration
    assert (os.path.exists(ref_seg_post_file))
    assert (os.path.exists(flo_seg_post_file))
    assert (os.path.exists(ref_file_sym_eq))
    assert (os.path.exists(flo_file_sym_eq))
    if not os.path.isdir(folder_inter):
        try:
            os.makedirs(folder_inter)
        except:
            print "Unable to create experiences path. Check the complete path..."
            raise
    planeRegistration(ref_seg_post_file, flo_seg_post_file, ref_file_sym_eq, flo_file_sym_eq,
                      path_trsf_ref_to_flo=ref_flo_file_trsf, path_residuals=ref_flo_file_res,
                      path_dices=ref_flo_file_dices, path_pairs_ref_flo=ref_flo_file_pairs,
                      background_ref=background_ref, background_flo=background_flo,
                      trsf_type=trsf_type, estimator=estimator, lts_fraction=lts_fraction, verbose=verbose, lazy=True)

    if not keep_sym:
        cmd = 'rm -rf ' + folder_sym
        if verbose:
            print cmd
        os.system(cmd)

    if not path_trsf_flo_ref:
        path_trsf_flo_ref = os.path.basename(ref_flo_file_trsf)
    if path_trsf_flo_ref:
        cmd = "cp " + ref_flo_file_trsf + ' ' + path_trsf_flo_ref
        if verbose:
            print cmd
        os.system(cmd)

    if not path_pairs_ref_flo:
        path_pairs_ref_flo = os.path.basename(ref_flo_file_pairs)
    if path_pairs_ref_flo:
        cmd = "cp " + ref_flo_file_pairs + ' ' + path_pairs_ref_flo
        if verbose:
            print cmd
        os.system(cmd)

    if not path_dices_ref_flo:
        path_dices_ref_flo = os.path.basename(ref_flo_file_dices)
    if path_dices_ref_flo:
        cmd = "cp " + ref_flo_file_dices + ' ' + path_dices_ref_flo
        if verbose:
            print cmd
        os.system(cmd)

    if not path_residuals_ref_flo:
        path_residuals_ref_flo = os.path.basename(ref_flo_file_res)
    if path_residuals_ref_flo:
        cmd = "cp " + ref_flo_file_res + ' ' + path_residuals_ref_flo
        if verbose:
            print cmd
        os.system(cmd)

    if not keep_inter:
        cmd = 'rm -rf ' + folder_inter
        if verbose:
            print cmd
        os.system(cmd)

    if delete_workspace and folder_tmp:
        cmd = 'rm -rf ' + folder_tmp
        if verbose:
            print cmd
        os.system(cmd)

    return path_trsf_flo_ref, path_pairs_ref_flo, path_dices_ref_flo, path_residuals_ref_flo


def to_256_bits_key(n, background=[0, 1], key_no_correspondences='no_correspondences', lut=None):
    '''
	Returns :
		0 if n in background
		255 if n is equal to key_no_correspondences or not lut.has_key(n) or not n in lut
		((n - d )% 255 ) + 1 if key_no_correspondences is None and lut is None
		((n - d )% 254 ) + 1 otherwise (default)
		where n <- lut[n] if lut is of type dict (may be of type list instead if one does not expect a lut but a set of cells with existing correspondences)
		(default = None), otherwise original parameter value is kept.
	This method is appropriated for Fiji LUT "glasbeyBGD"
	'''
    modulo = 254
    if key_no_correspondences == None and not lut:
        modulo = 255
    else:
        if type(n) == type(key_no_correspondences) and n == key_no_correspondences:
            return 255

    assert (type(n) == int)
    if n in background:
        return 0
    if lut:
        if type(lut) == dict:
            if not lut.has_key(n):
                return 255
            n = lut[n]
        else:
            if not n in lut:
                return 255

    d = len([i for i in background if i < n])  # dividende
    l = ((n - d) % modulo) + 1
    return l


def segmentationRelabellingAtTime(path_segmentation, path_output, lineage, time_point, reverse_lut=None, lut=None,
                                  trsf=None, template=None, voxelsize=None, iso=None, dimensions=None,
                                  backgroundLabels=[0, 1], visu=True, verbose=False):
    '''
	Builds the segmented image path_output encoded in 8 bits from input segmentation path_segmentation, with a relabellization and, if visu=True, a thin erosion of the cells for vizualization
	(the user may choose the 'glasbeyBGD' LUT in Fiji for path_output vizualisation). The relabelling is processed with respect to the lineage relabelling, to lut (or reverse_lut), and to time_point.

	Parameters:

		lut: if specified, do not specify reverse_lut.
		     The lut option can be of two types:
		     	- list (or set or tuple): in that case, associates value 255 if the relabelled value l of path_segmentation is not in lut. Otherwise, associate the "to_256_bits_key" conversion of l.
		     	- dict: in that case, associates value 255 if the relabelled value l of path_segmentation is not in lut.keys(). Otherwise, associate the "to_256_bits_key" conversion of lut[l].
		     This parameter may be useful when the user wants to synchronize a reference image considering a floating sequence registration onto the reference and when the reference-to-floating correspondences dictionary is known.
		     In that case of use, the lut parameter should be set with the list given by correspondences.keys(), so that cells with no correspondences will be automatically set to output value 255.

		reverse_lut: if specified, do not specify lut. Must be a dict type and is equivalent than to set lut option with the value {v:k for k,v in reverse_luts.items()} (ie the reverse version of reverse_lut).
					 May be a useful parameter when the user wants to synchronise a floating image onto a reference sequence and when the reference-to-floating correspondences dictionary is known.
					 In that case of use, the reverse_lut parameter should be set with the correspondences dictionary.

	Optional parameters:

		trsf (default=None):  if None or False, path_output will have the same field of view as path_segmentation.
							  If trsf is a path for a transformation or a transformation matrix 4*4 (np.array or list) that transforms spatially a target field into path_segmentation (with convention trsf = T_segmentation<-target)
							  then the function applies this transformation to path_segmentation (with optional dimensions, voxelsize, iso and template parameters).

		voxelsize (default=None): forces path_output voxelsizes to the given values (must be a list or tuple of length 3) (may be used independently to trsf parameter).

		iso (default=None): forces path_output to have isotropic voxel size of given value (float type) (may be used independently to trsf parameter).

		dimensions (default=None): forces the path_output dimensions (must be an integer list or tuple of length 3) (may be used independently to trsf parameter).

		template (default=None): path to a template image that forces the path_output dimensions and voxelsize to correspond to the tmeplate image (must be a str type). If used, do not use iso, voxelsize neither dimensions parameters.

		backgroundLabels (default=[0,1]): labels of the input image to be set to 0 in the output image

		visu (default=True): option for vizualisation. If True, then a thin erosion of original cells borders is processed to see well cell borders in the output image.

		verbose (defaut=False): option for verbosity of the function.
	'''
    assert os.path.exists(path_segmentation)
    if template:
        assert os.path.exists(template)

    if not os.path.dirname(path_output):
        try:
            os.makedirs(os.path.dirname(path_output))
        except:
            print "Unable to create output path. Check the complete path..."
            raise

    if type(trsf) != np.ndarray and trsf == None:
        trsf = False

    out_original_labels = {}
    rev_lut = lineage.reverse_luts[time_point]

    if reverse_lut:
        assert not lut
        lut = {v: k for k, v in reverse_lut.items()}

    if not lut:
        for k, v in rev_lut.items():
			k_256_bits = to_256_bits_key(k, backgroundLabels)
            if not out_original_labels.has_key(k_256_bits):
                out_original_labels[k_256_bits] = set()
            out_original_labels[k_256_bits].add(v)
    else:
        for k, v in rev_lut.items():
            k_256_bits = to_256_bits_key(k, backgroundLabels, lut=lut)
            if not out_original_labels.has_key(k_256_bits):
                out_original_labels[k_256_bits] = set()
            out_original_labels[k_256_bits].add(v)

    # We add the background correspondence...
    if not out_original_labels.has_key(0):
        out_original_labels[0] = set()
    for v in backgroundLabels:
        out_original_labels[0].add(v)
        out_original_labels[0].add(v)

    txt = ""
    for k, v in out_original_labels.items():
        for val in v:
            txt += "%d %d\n" % (val, k)

    path_tmp_lut = "tmp_segmentationRelabellingAtTime_t%03d.lut" % time_point
    if verbose:
        print "Writing temporary file %s" % path_tmp_lut
    f = open(path_tmp_lut, "w")
    f.write(txt)
    f.close()

    tmp = path_output + ".tmp_segmentationRelabellingAtTime_t%03d.inr" % time_point

    trsf_tmp_file = 'tmp_segmentationRelabellingAtTime_t%03d.trsf' % time_point
    # Check if need to compute the trsf (->set cpt_trsf) or if it is a filename (-> copy to temporary file)
    if type(trsf) != bool:
        if type(trsf) == str:
            assert os.path.exists(
                trsf), "Error: unexpected value '%s' for trsf parameter (file not found). See help." % trsf
            cmd = "cp %s %s" % (trsf, trsf_tmp_file)
            if verbose:
                print cmd
            os.system(cmd)
    # If trsf is a tabular, then write the temporary trsf file
    if type(trsf) != bool:
        if type(trsf) == list:
            trsf = np.array(trsf)
        if type(trsf) == np.ndarray:
            f = open(trsf_tmp_file, "w")
            f.write(str(trsf).replace('[', '').replace(']', ''))
            f.close()

    if template:
        assert (not dimensions and not voxelsize and not iso)
    assert type(trsf) != bool or not trsf

    if voxelsize or dimensions or iso or template or type(trsf) != bool:
        if type(trsf) != bool:
            assert os.path.exists(trsf_tmp_file), "Unexpected problem while setting the transformation T_flo<-ref."
            apply_trsf(path_segmentation, path_trsf=trsf_tmp_file, path_output=path_output, nearest=True,
                       template=template, voxelsize=voxelsize, dimensions=dimensions, iso=iso, verbose=verbose)
        else:
            apply_trsf(path_segmentation, path_output=path_output, nearest=True, template=template, voxelsize=voxelsize,
                       dimensions=dimensions, iso=iso, verbose=verbose)

    if visu:
        if voxelsize or dimensions or iso or type(trsf) != bool:
            labelBorders(path_output, tmp, verbose)
        else:
            labelBorders(path_segmentation, tmp, verbose)
        Logic(tmp, tmp, Mode='inv', verbose=verbose)

    if voxelsize or dimensions or iso or type(trsf) != bool:
        fuseLabelsWithLUT(path_output, path_output, path_tmp_lut, U8=True, verbose=verbose)
    else:
        fuseLabelsWithLUT(path_segmentation, path_output, path_tmp_lut, U8=True, verbose=verbose)

    if os.path.exists(path_tmp_lut):
        cmd = 'rm ' + path_tmp_lut
        if verbose:
            print cmd
        os.system(cmd)

    if visu:
        Logic(tmp, path_output, path_output, Mode='mask', verbose=verbose)
        cmd = "rm " + tmp
        if verbose:
            print cmd
        os.system(cmd)

    # RM TEMPORARY FILE
    if type(trsf) != bool:
        if os.path.exists(trsf_tmp_file):
            cmd = "rm %s" % trsf_tmp_file
            if verbose:
                print cmd
            os.system(cmd)


def associateSegmentationsAtTimesFromLineageCorrespondences(path_ref, path_flo, path_ref_out, path_flo_out,
                                                            lineage_ref, lineage_flo, correspondences, time_ref,
                                                            time_flo,
                                                            trsf=False, estimator='lts', lts_fraction=1.0,
                                                            # template=None,
                                                            voxelsize=None, iso=None, dimensions=None,
                                                            backgroundLabels=[0, 1], visu=True, verbose=False):
    '''
	Builds the segmented images path_ref_out and path_flo_out from input segmentations path_ref and path_flo, with a relabellization and thin erosion of the cells for vizualization
	(the user may choose the 'glasbeyBGD' LUT in Fiji). The relabelling is processed with respect to the lineage_ref, lineage_flo, correspondences, time_ref and time_flo parameters.
	Optional parameters:
		trsf (default=False): if False, path_flo_out will have the same dimensions as path_flo.
							  If True or 'rigid', the function will compute automatically a rigid transformation that registers lineage_ref with lineage_flo at given time_ref and time_flo,
							  under the constraint of given correspondences, and apply the transformation to path_flo_out so that dimension(path_flo_out)=dimension(path_ref_out).
							  If 'affine', the same as if True, but with an affine transformation instead of a rigid.
							  In these two cases, a transformation will be automatically computed with the parameters 'estimator' and 'lts_fraction' that are related to the pointCloudRegistration function.
							  If trsf is a path for a transformation or a transformation matrix 4*4 (np.array or list) that registers spatially path_flo into path_ref (with convention trsf = T_flo<-ref)
							  then the function applies this transformation to path_flo so that dimension(path_flo_out)=dimension(path_ref_out).
		voxelsize (default=None): forces path_ref_out and path_flo_out voxelsizes to the given values (must be a list or tuple of length 3).
		iso (default=None): forces path_ref_out and path_flo_out to have isotropic voxel size of given value (float type).
		dimensions (default=None): forces the path_ref_out and path_flo_out output dimensions (must be an integer list or tuple of length 3).
		backgroundLabels (default=[0,1]): labels of the input images to be set to 0 in the output images
		visu (default=True): option for vizualisation. If True, then a thin erosion of original cells borders is processed to see well cell borders in the output images.
		verbose (defaut=False): option for verbosity of the function.
	'''

    assert (os.path.exists(path_ref))
    assert (os.path.exists(path_flo))
    key_no_correspondences = 'no_correspondences'

    if type(trsf) != np.ndarray and trsf == None:
        trsf = False

    def get_correspondence(lineage_ref, lineage_flo, time_ref, time_flo, correspondences,
                           key_no_correspondences='no_correspondences'):
        '''
		Returns a dict structure with each item as follows:
		keys, values:
			- key_no_correspondences, (set_ref, set_flo) :
				key for labels that belong to lineages at given times and that do not exist in correspondences dictionary
				set_ref for relabellings that belong to lineages_ref at time_ref and not belong to correspondences (ie correspondences.has_key(l) returns False for each element l of set_ref)
				set_flo for relabellings that belong to lineages_flo at time_flo and not belong to correspondences (ie l in correspondences.values() returns False for each element l of set_flo)
			- label : (set_ref, set_flo) :
				set_ref for set of relabellings from lineage_ref that is or are associated to the set of relabellings set_flo from lineage_flo (all elements of set_ref are existing cells of lineage_ref at time_ref)
				set_flo for set of relabellings from lineage_flo that is or are associated to the set of relabellings set_ref from lineage_ref (all elements of set_flo are existing cells of lineage_flo at time_flo)
				A property that must be respected is that at least one of the two sets is a singleton (ie has one and only one element). The other set can have 0, 1 or more elements.
				The key label is a label which value is the "greatest common ancestor" of elements of set_ref so that all the elements of set_flo are associated to this ancestor or its descendents.

		'''
        ref_correspondences = correspondences.keys()
        flo_correspondences = correspondences.values()

        relabels_ref = lineage_ref.reverse_luts[time_ref].keys()
        relabels_flo = lineage_flo.reverse_luts[time_flo].keys()

        not_corresponding_ref = set(relabels_ref).difference(ref_correspondences).intersection(relabels_ref)
        not_corresponding_flo = set(relabels_flo).difference(flo_correspondences).intersection(relabels_flo)

        out = {key_no_correspondences: (not_corresponding_ref,
                                        not_corresponding_flo)}  # ref and flo sets of cells that have to be visualized as cells with no correspondences

        relabels_ref = list(set(relabels_ref).difference(
            not_corresponding_ref))  # we remove from relabels_ref the cells with no correspondences
        # relabels_ref.sort()
        relabels_flo = list(set(relabels_flo).difference(
            not_corresponding_flo))  # we remove from relabels_flo the cells with no correspondences
        # relabels_flo.sort()

        relabels_flo_as_ancestors = set()
        scanned_flo_set = set()

        for relabel_ref in relabels_ref:
            assert relabel_ref in ref_correspondences
            # Cas ou relabel_ref existe dans le dict 'correspondences'
            corres_flo = correspondences[relabel_ref]
            if corres_flo in relabels_flo:
                # Cas simple : correspondance 1 a 1 trouvee
                assert not out.has_key(relabel_ref)
                out[relabel_ref] = ({relabel_ref}, {corres_flo})
                scanned_flo_set.add(corres_flo)
            else:
                # Cas possibles :
                # 1. relabels_flo contient des descendants de corres_flo -> on associe relabel_ref a tous les descendants de corres_flo qui appartiennent au dictionnaire 'correspondences'
                # 2. relabels_flo contient un ancetre de corres_flo -> on traite le cas dans le sens "inverse" en cherchant tous les correspondants de l'ancetre dans relabels_ref qui appartiennent au dictionnaire 'correspondences'
                #	 -> on procede a la recherche "inverse" afin de determiner les correspondances de l'ancetre dans relabels_ref (boucle finale sur relabels_flo_loop)
                # 3. relabels_flo ne contient aucune cellule de la lignee de corres_flo -> out[relabel_ref]=(relabel_ref, set())

                corres_flo_descendents = lineage_flo.relabelled_descendents(corres_flo)
                corres_flo_ancestors = lineage_flo.relabelled_ancestors(corres_flo)

                # Test cas 1:
                intersection = set(relabels_flo).intersection(corres_flo_descendents)
                if intersection:
                    # Cas 1:
                    assert (not set(relabels_flo).intersection(
                        corres_flo_ancestors))  # cohabitation logiquement impossible entre descendants et ancetres d'une meme cellule a un instant donne du developpement d'un embryon

                    corresponding = intersection.intersection(
                        correspondences.values())  # extraction des descendants du label en correspondance avec relabel_ref qui appartiennent au dictionnaire 'correspondences' -> out key=relabel_ref
                    not_corresponding = intersection.difference(correspondences.values()).intersection(
                        intersection)  # extraction des descendants du label en correspondance avec relabel_ref qui n'appartiennent pas au dictionnaire 'correspondences' -> out key=key_no_correspondences
                    assert not not_corresponding
                    assert not out.has_key(relabel_ref)
                    out[relabel_ref] = ({relabel_ref}, corresponding)
                    scanned_flo_set = scanned_flo_set.union(corresponding)
                else:
                    # Cas 2 ou 3:
                    # Test cas 2:
                    intersection = set(relabels_flo).intersection(corres_flo_ancestors)
                    assert (len(
                        intersection) <= 1)  # on ne devrait pas avoir une cohabitation de plusieurs ancetres d'un unique label a un pas de temps donne
                    if intersection:
                        # Cas 2:
                        intersection_value = intersection.pop()
                        assert intersection_value in flo_correspondences
                        relabels_flo_as_ancestors.add(intersection_value)
                    else:
                        # Cas 3:
                        assert not out.has_key(relabel_ref)
                        out[relabel_ref] = ({relabel_ref}, set())

        rev_correspondences = {v: k for k, v in correspondences.items()}
        for relabel_flo in relabels_flo_as_ancestors:
            # Cette boucle concerne tous les labels de relabels_flo qui s'associent a un ancetre de cellules existentes dans relabels_ref...
            assert rev_correspondences.has_key(relabel_flo)
            corres_ref = rev_correspondences[relabel_flo]
			corres_ref_descendents = lineage_ref.relabelled_descendents(corres_ref)
			intersection = set(relabels_ref).intersection(corres_ref_descendents)
			assert len(intersection) > 0
            corresponding = intersection.intersection(
                rev_correspondences.values())  # extraction des descendants du label en correspondance avec relabel_flo qui appartiennent au dictionnaire 'rev_correspondences' -> out key=corres_ref
            not_corresponding = intersection.difference(rev_correspondences.values()).intersection(
                intersection)  # extraction des descendants du label en correspondance avec relabel_flo qui n'appartiennent pas au dictionnaire 'rev_correspondences' -> out key=key_no_correspondences (theoriquement, deja ajoute a ce stade)
            assert len(not_corresponding) == len(not_corresponding.intersection(out[key_no_correspondences][0]))
            assert not out.has_key(corres_ref)
			out[corres_ref] = (corresponding, {relabel_flo})
            scanned_flo_set.add(relabel_flo)

        # Etape finale : aucun element de relabels_ref ni de relabels_flo ne doit avoir ete omis
        last_loop_flo = scanned_flo_set.difference(relabels_flo)
        for relabel_flo in last_loop_flo:  # labels de relabels_flo qui n'ont pas encore ete scannes, i.e. qui n'ont soit pas de correspondants dans rev_correspondences (->key_no_correspondences), soit pas de correspondants dont la lignee appartient a relabels_ref (->key=rev_correspondences[relabel_flo])
            assert relabel_flo in flo_correspondences
            assert not out.has_key(rev_correspondences[relabel_flo])
            out[rev_correspondences[relabel_flo]] = (set(), {relabel_flo})

        return out

    # extraction des correspondences entre cellules de reference et flottantes avec nouvel etiquetage
    out_relabelled = get_correspondence(lineage_ref, lineage_flo, time_ref, time_flo, correspondences,
                                        key_no_correspondences)

    out_original_labels = {}
    rev_lut_ref = lineage_ref.reverse_luts[time_ref]
    rev_lut_flo = lineage_flo.reverse_luts[time_flo]

    for k, (v_ref, v_flo) in out_relabelled.items():
        k_256_bits = to_256_bits_key(k, backgroundLabels, key_no_correspondences)
        if not out_original_labels.has_key(k_256_bits):
            out_original_labels[k_256_bits] = (set(), set())
        for relabel_ref in v_ref:
            assert (rev_lut_ref.has_key(relabel_ref))
            out_original_labels[k_256_bits][0].add(rev_lut_ref[relabel_ref])
        for relabel_flo in v_flo:
            assert (rev_lut_flo.has_key(relabel_flo))
            out_original_labels[k_256_bits][1].add(rev_lut_flo[relabel_flo])

    # We add the background correspondence...
    if not out_original_labels.has_key(0):
        out_original_labels[0] = (set(), set())
    for v in backgroundLabels:
        out_original_labels[0][0].add(v)
        out_original_labels[0][1].add(v)

    txt_ref = ""
    txt_flo = ""
    for k, (v_ref, v_flo) in out_original_labels.items():
        for val in v_ref:
            txt_ref += "%d %d\n" % (val, k)
        for val in v_flo:
            txt_flo += "%d %d\n" % (val, k)

    path_tmp_lut_ref = "tmp_associateSegmentationsAtTimesFromLineageCorrespondences_t%03d_t%03d_ref.lut" % (
    time_ref, time_flo)
    path_tmp_lut_flo = 'tmp_associateSegmentationsAtTimesFromLineageCorrespondences_t%03d_t%03d_flo.lut' % (
    time_ref, time_flo)
    if verbose:
        print "Wrinting temporary file %s" % path_tmp_lut_ref
    f = open(path_tmp_lut_ref, "w")
    f.write(txt_ref)
    f.close()
    if verbose:
        print "Wrinting temporary file %s" % path_tmp_lut_flo
    g = open(path_tmp_lut_flo, "w")
    g.write(txt_flo)
    g.close()

    tmp_ref = path_ref_out + ".tmp_t%03d_t%03d.inr" % (time_ref, time_flo)
    tmp_flo = path_flo_out + ".tmp_t%03d_t%03d.inr" % (time_ref, time_flo)

    cpt_trsf = None
    trsf_tmp_file = 'tmp_associateSegmentationsAtTimesFromLineageCorrespondences_flo_ref.trsf'
    # Check if need to compute the trsf (->set cpt_trsf) or if it is a filename (-> copy to temporary file)
    if type(trsf) == bool and trsf:
        cpt_trsf = 'rigid'
    if type(trsf) != bool:
        if type(trsf) == str:
            if trsf == 'rigid':
                cpt_trsf = 'rigid'
            if trsf == 'affine':
                cpt_trsf = 'affine'
            else:
                assert os.path.exists(
                    trsf), "Error: unexpected value '%s' for trsf parameter (file not found). See help." % trsf
                cmd = "cp %s %s" % (trsf, trsf_tmp_file)
                if verbose:
                    print cmd
                os.system(cmd)
    # If need to compute trsf (->set trsf as a np.array)
    if cpt_trsf:
        point_cloud_ref = lineage_ref.relabelled_barycenters_at_time(time_ref)
        point_cloud_flo = lineage_flo.relabelled_barycenters_at_time(time_flo)
        if type(backgroundLabels) == list or type(backgroundLabels) == tuple or type(backgroundLabels) == set:
            for v in backgroundLabels:
                point_cloud_ref.pop(v, None)
        if type(backgroundLabels) == int:
            point_cloud_ref.pop(backgroundLabels, None)

        out = pointCloudRegistration(point_cloud_ref, point_cloud_flo, correspondences, skip_not_found=True,
									 trsf_type=cpt_trsf, estimator=estimator, lts_fraction=lts_fraction,
									 lazy=False, verbose=verbose)
        trsf = np.array(out['trsf'])
    # If trsf is a tabular, then write the temporary trsf file
    if type(trsf) != bool:
        if type(trsf) == list:
            trsf = np.array(trsf)
        if type(trsf) == np.ndarray:
            f = open(trsf_tmp_file, "w")
            f.write(str(trsf).replace('[', '').replace(']', ''))
            f.close()

    if voxelsize or dimensions or iso:
        assert (not dimensions and not voxelsize) or (not dimensions and not iso) or (
                    not voxelsize and not iso)  # not the two options at the same time
        apply_trsf(path_ref, path_output=path_ref_out, nearest=True, voxelsize=voxelsize, dimensions=dimensions,
                   iso=iso, verbose=verbose)
    # apply_trsf(path_flo, path_output=path_flo_out, nearest=True, voxelsize=voxelsize, dimensions=dimensions, iso=iso, verbose=verbose)

    # Checks the existence of trsf_tmp_file if trsf was mentionned
    assert type(trsf) != bool or not trsf
    if type(trsf) != bool:
        assert os.path.exists(trsf_tmp_file), "Unexpected problem while setting the transformation T_flo<-ref."
        if voxelsize or dimensions or iso:
            apply_trsf(path_flo, path_trsf=trsf_tmp_file, path_output=path_flo_out, template=path_ref_out, nearest=True,
                       verbose=verbose)
        else:
            apply_trsf(path_flo, path_trsf=trsf_tmp_file, path_output=path_flo_out, template=path_ref, nearest=True,
                       verbose=verbose)
    else:
        if voxelsize or dimensions or iso:
            apply_trsf(path_flo, path_output=path_flo_out, nearest=True, voxelsize=voxelsize, dimensions=dimensions,
                       iso=iso, verbose=verbose)

    if visu:
        if voxelsize or dimensions or iso:
            labelBorders(path_ref_out, tmp_ref, verbose)
        # cmd=path_labelborders + " " + path_ref_out + ' ' + tmp_ref
        else:
            labelBorders(path_ref, tmp_ref, verbose)
        if type(trsf) != bool or voxelsize or dimensions or iso:
            labelBorders(path_flo_out, tmp_flo, verbose)
        else:
            labelBorders(path_flo, tmp_flo, verbose)
        Logic(tmp_ref, tmp_ref, Mode='inv', verbose=verbose)
        Logic(tmp_flo, tmp_flo, Mode='inv', verbose=verbose)

    if voxelsize or dimensions or iso:
        fuseLabelsWithLUT(path_ref_out, path_ref_out, path_tmp_lut_ref, U8=True, verbose=verbose)
    else:
        fuseLabelsWithLUT(path_ref, path_ref_out, path_tmp_lut_ref, U8=True, verbose=verbose)
    if type(trsf) != bool:
        fuseLabelsWithLUT(path_flo_out, path_flo_out, path_tmp_lut_flo, U8=True, verbose=verbose)
    else:
        fuseLabelsWithLUT(path_flo, path_flo_out, path_tmp_lut_flo, U8=True, verbose=verbose)

    if os.path.exists(path_tmp_lut_ref):
        cmd = 'rm ' + path_tmp_lut_ref
        if verbose:
            print cmd
        os.system(cmd)
    if os.path.exists(path_tmp_lut_flo):
        cmd = 'rm ' + path_tmp_lut_flo
        if verbose:
            print cmd
        os.system(cmd)

    if visu:
        Logic(tmp_ref, path_ref_out, path_ref_out, Mode='mask', verbose=verbose)
        Logic(tmp_flo, path_flo_out, path_flo_out, Mode='mask', verbose=verbose)
        cmd = "rm " + tmp_ref + ' ' + tmp_flo
        if verbose:
            print cmd
        os.system(cmd)

    # RM TEMPORARY FILE
    if type(trsf) != bool:
        if os.path.exists(trsf_tmp_file):
            cmd = "rm %s" % trsf_tmp_file
            if verbose:
                print cmd
            os.system(cmd)


def associateRefFloLabels(path_ref, path_flo, path_pairs_ref_flo, path_ref_out, path_flo_out, path_labels_out=None,
                          path_trsf_flo_ref=None, zeros=True, backgroundLabels=[0, 1], ref=False, force=False,
                          visu=False, verbose=False):
    '''
   path_ref : correspond au parametre path_label_ref de la methode planeRegistration
   path_flo : correspond au parametre path_label_flo de la methode planeRegistration
   path_pairs_ref_flo : fichier de correspondances label-ref / label-flo (tel qu'en sortie de planeRegistration / pointCloudRegistration path_pairs_ref_flo %s) (accepte aussi les dictionnaires avec les items de correspondences sous forme de key=ref_label: value=flo_label)
   path_trsf_flo_ref : fichier de transformation flo<-ref (tel qu'en sortie de planeRegistration path_trsf_flo_ref %s) permettant de calculer un repositionnement de l'image flottante dans le referentiel de l'image de reference (accepte aussi les transformations sous forme de numpy.ndarray ou de "list of lists")
   zeros : les zeros en entree sont conserves en sortie si True
   backgroundLabels : si zeros a True, met a zero les labels stipules dans la liste passee en argument
   force : force la valeur de sortie des labels (tableau d'entree n*3) si True
   ref : utilise les labels de l'image in comme labels de reference si True
   visu : (defaut: False) mise a zero des voxels a la frontiere entre deux labels differents, utile pour la visualisation des images
   verbose : option de verbosite
  '''
    flag_remove_path_correspondences = False
    flag_remove_path_trsf = False
    assert (os.path.exists(path_ref))
    assert (os.path.exists(path_flo))
    if type(path_pairs_ref_flo) == dict:
        from morpheme_lineage import write_correspondences
        correspondences = path_pairs_ref_flo
        path_pairs_ref_flo = "tmp_ref_flo_associateRefFloLabels.pairs"
        write_correspondences(correspondences, path_pairs_ref_flo)
        flag_remove_path_correspondences = True
    assert (type(path_pairs_ref_flo) == str and os.path.exists(path_pairs_ref_flo))

    tmp_ref = path_ref_out + ".tmp.inr"
	tmp_flo = path_flo_out + ".tmp.inr"

    if (type(path_trsf_flo_ref) == list or type(path_trsf_flo_ref) == np.ndarray):
        trsf = np.array(path_trsf_flo_ref)
        path_trsf_flo_ref = "tmp_ref_flo_associateRefFloLabels.trsf"
        f = open(path_trsf_flo_ref, "w")
        f.write(str(trsf).replace('[', '').replace(']', ''))
        f.close()
        flag_remove_path_trsf = True

    if path_trsf_flo_ref:
        assert (os.path.exists(path_trsf_flo_ref))
        apply_trsf(path_flo, path_trsf=path_trsf_flo_ref, path_output=path_flo_out, template=path_ref, nearest=True,
                   verbose=verbose)

    if zeros and len(backgroundLabels):
        resetLabels(path_ref, path_ref_out, backgroundLabels, verbose=verbose)
        if path_trsf_flo_ref:
            resetLabels(path_flo_out, path_flo_out, backgroundLabels, verbose=verbose)
        else:
            resetLabels(path_flo, path_flo_out, backgroundLabels, verbose=verbose)

    if visu:
        if zeros and len(backgroundLabels):
            labelBorders(path_ref_out, tmp_ref, verbose)
            # cmd=path_labelborders + " " + path_ref_out + ' ' + tmp_ref
        else:
            labelBorders(path_ref_out, tmp_ref, verbose)
		if path_trsf_flo_ref or (zeros and len(backgroundLabels)):
			labelBorders(path_flo_out, tmp_flo, verbose)
		else:
			labelBorders(path_flo, tmp_flo, verbose)
		Logic(tmp_ref, tmp_ref, Mode='inv', verbose=verbose)
		Logic(tmp_flo, tmp_flo, Mode='inv', verbose=verbose)

    if zeros and len(backgroundLabels):
        associateLabels(path_ref_out, path_flo_out, path_pairs_ref_flo, path_ref_out, path_flo_out,
                        path_labels_out=path_labels_out, zeros=zeros, ref=ref, force=force, verbose=verbose)
    else:
        if path_trsf_flo_ref:
            associateLabels(path_ref, path_flo_out, path_pairs_ref_flo, path_ref_out, path_flo_out,
                            path_labels_out=path_labels_out, zeros=zeros, ref=ref, force=force, verbose=verbose)
        else:
            associateLabels(path_ref, path_flo, path_pairs_ref_flo, path_ref_out, path_flo_out,
                            path_labels_out=path_labels_out, zeros=zeros, ref=ref, force=force, verbose=verbose)

    if flag_remove_path_correspondences:
        cmd = 'rm ' + path_pairs_ref_flo
        if verbose:
            print cmd
        os.system(cmd)

    if flag_remove_path_trsf:
        cmd = 'rm ' + path_trsf_flo_ref
        if verbose:
            print cmd
        os.system(cmd)

    if visu:
        Logic(tmp_ref, path_ref_out, path_ref_out, Mode='mask', verbose=verbose)
        Logic(tmp_flo, path_flo_out, path_flo_out, Mode='mask', verbose=verbose)
        cmd = "rm " + tmp_ref + ' ' + tmp_flo
        if verbose:
            print cmd
        os.system(cmd)


def sisters_association(point_cloud_ref, point_cloud_flo, correspondences, sisters_ref, sisters_flo,
                        threshold_angle=90.0,
                        background_ref=1, background_flo=1, trsf_type='affine', estimator='lts', lts_fraction=1.0,
                        bash_options=None,
                        verbose=False):
    '''
	Function that determines correspondences between reference and floating pairs of labels, given a reference
	and a floating point-clouds (pairs of labels to be associated must belong to these point-clouds) and given a set of correspondences.
	Inputs:
		point_cloud_ref : dict
		point_cloud_flo : dict
		correspondences : dict
		sisters_ref     : list or tuple of length 2
		sisters_flo     : list or tuple of length 2
	Optional association parameter:
		threshold_angle : positive angle in degrees (default value is 90).
						  Sisters association is kept for the result iif sisters axes form an absolute angle that is smaller than the
						  given threshold

	Optional parameters (for pointCloudRegistration cpp_wrapping imported function):
		background_ref  : label de fond a ignorer dans l'image ref correspondante
		background_flo  : label de fond a ignorer dans l'image flo correspondante
		trsf_type       : computes 'affine' (set as default) or 'rigid' transformation
		estimator       : transformation estimator:
   							'wlts': weighted least trimmed squares
   							'lts': least trimmed squares (default)
   							'wls': weighted least squares
							'ls': least squares
		lts_fraction    : for trimmed estimations, fraction of pairs that are kept (default: 1.0)
		bash_options    : parameter of type str (None by default) to set iif the user wants to add specific bash options for the
						  execution of the program pointCloudRegistration from library vt (morpheme_privat)
	Outputs:
		new_correspondences : dict that is a copy of the input correspondences with the new keys corresponding to sisters_ref labels
							  and values corresponding to associated sisters_flo labels so that it minimizes the squares distances
							  between associated labels with respect to the T_flo<-ref transformation computed in function of the input
							  correspondences
		(angle, volume_ratio_ref, volume_ratio_flo) : tuple of length 3:
							angle: value of the angle between division axes of sisters_ref and sisters_flo
							volume_ratio_ref: volume ratio between reference sisters under constraint that ratio <=1 (if the volumes
											  are not known, returns 1)
							volume_ratio_flo: volume ratio between floating sisters following the association determined by this function
											  with the reference sisters (means that it may be > 1) (if the volumes are not known, returns 1)
	'''
    from math import pi
    if len(sisters_ref) != 2 or len(sisters_flo) != 2:
        if verbose:
            print "sisters_ref or sisters_flo does not have the expected number of elements (2)."
        return {}, (None, None, None)

    out = pointCloudRegistration(point_cloud_ref, point_cloud_flo, correspondences, skip_not_found=True,
                                 background_ref=background_ref, background_flo=background_flo,
                                 trsf_type=trsf_type, estimator=estimator, lts_fraction=lts_fraction,
                                 lazy=False, verbose=verbose)
    trsf = np.array(out['trsf'])
    daughter_ref_0 = np.array(point_cloud_ref[sisters_ref[0]][0:3] + (1,))
    daughter_ref_1 = np.array(point_cloud_ref[sisters_ref[1]][0:3] + (1,))
    daughter_flo_0 = np.array(point_cloud_flo[sisters_flo[0]][0:3] + (1,))
    daughter_flo_1 = np.array(point_cloud_flo[sisters_flo[1]][0:3] + (1,))
    trsf_daughter_ref_0 = np.dot(trsf, daughter_ref_0)
    trsf_daughter_ref_1 = np.dot(trsf, daughter_ref_1)
    trsf_vector_ref = trsf_daughter_ref_1 - trsf_daughter_ref_0
    vector_flo = daughter_flo_1 - daughter_flo_0
    ps = np.vdot(trsf_vector_ref, vector_flo)
    psn = ps / (np.linalg.norm(trsf_vector_ref) * np.linalg.norm(
        vector_flo))  # psn est le produit scalaire normalise entre les deux vecteurs definis par les centres de gravite des cellules filles juste apres la division
    angle_degree = np.arccos(abs(psn)) * 180 / pi
    # Pour l'appariement de cellules, la minimisation au sens des moindres carres de distance equivaut a associer daughter_ref_0 a daughter_flo_0 et daughter_ref_1 a daughter_ref_1 ssi le psn > 0, et sinon on inverse les associations
    volume_ratio_ref = 1
    volume_ratio_flo = 1
    new_correspondences = {}
    # new_correspondences = correspondences.copy()
    if psn > 0:
        if angle_degree < threshold_angle:
			new_correspondences[sisters_ref[0]] = sisters_flo[0]
            new_correspondences[sisters_ref[1]] = sisters_flo[1]
        if len(point_cloud_ref[sisters_ref[0]]) == len(point_cloud_ref[sisters_ref[1]]) == len(
                point_cloud_flo[sisters_flo[0]]) == len(point_cloud_flo[sisters_flo[1]]) == 4:
            volume_ratio_ref = point_cloud_ref[sisters_ref[1]][3] / point_cloud_ref[sisters_ref[0]][3]
            volume_ratio_flo = point_cloud_flo[sisters_flo[1]][3] / point_cloud_flo[sisters_flo[0]][3]
    else:
        if angle_degree < threshold_angle:
            new_correspondences[sisters_ref[0]] = sisters_flo[1]
            new_correspondences[sisters_ref[1]] = sisters_flo[0]
        if len(point_cloud_ref[sisters_ref[0]]) == len(point_cloud_ref[sisters_ref[1]]) == len(
                point_cloud_flo[sisters_flo[0]]) == len(point_cloud_flo[sisters_flo[1]]) == 4:
            volume_ratio_ref = point_cloud_ref[sisters_ref[1]][3] / point_cloud_ref[sisters_ref[0]][3]
            volume_ratio_flo = point_cloud_flo[sisters_flo[0]][3] / point_cloud_flo[sisters_flo[1]][3]
    if volume_ratio_ref > 1:
        volume_ratio_ref = 1 / volume_ratio_ref
        volume_ratio_flo = 1 / volume_ratio_flo
    return new_correspondences, (angle_degree, volume_ratio_ref, volume_ratio_flo)


# label_mere_ref=lineage.
# volume_naissance_mere_ref=?
#
# return new_correspondences, (angle_degree, volume_naissance_mere_ref, volume_naissance_fille_1_ref, volume_naissance_fille_2_ref, volume_naissance_mere_flo, volume_naissance_fille_1_flo, volume_naissance_fille_2_flo )


def propagate_all_correspondences(lineage_ref, lineage_flo, correspondences, starting_time_ref=None,
                                  threshold_angle=90.0, stop_when_blocked=False, verbose=False):
    """
	Function for propagation of all the cell-to-cell correspondences between given lineage_ref and lineage_flo (as Morpheme_lineage instances).
	Propagation will start at specified starting_time_ref (by default, it starts at the first time-point of lineage_ref).
	The parameter 'correspondences' must be a dictionary with keys (resp. values) as relabellings of reference (resp. floating) embryo cells.
	This parameter can also be a filename at readable format for the "readLUT" function that contains the relabelled cell-to-cell
	correspondences between ref and flo embryos.
	Optional parameter for sisters_association:
			threshold_angle : positive angle in degrees (default value is 90).
						      Sisters association is kept for the result iif sisters axes form an absolute angle that is smaller than the
						      given threshold


	Returns:
		new_correspondences    : dict structure with all the cell-to-cell correspondences built by this function.
		unpropagated_cells_ref : list of relabelled cells of lineage_ref that have not been propagated due to unrespected threshold condition.
		scalars                : dict structure with keys as relabelled cells that have been propagated (or tried to be propagated) with
								 their associated scalar measures computed by the function "sisters_association".
	"""
    import morpheme_lineage as ml
    if type(correspondences) == str:
        correspondences = readLUT(correspondences)
    assert type(correspondences) == dict
    new_correspondences = correspondences.copy()
    if starting_time_ref == None:
        starting_time_ref = lineage_ref.timepoints()[0]

    scalars = {}
    unpropagated_cells_ref = []
    next_time_ref = starting_time_ref
    while next_time_ref < lineage_ref.timepoints()[-1]:
        if verbose:
            print "Time-point ref t%03d" % next_time_ref
        next_labels_ref, next_time_ref = lineage_ref.relabelled_next_death(next_time_ref)
        # Propagation of correspondences for all cells of lineage_ref that die at time next_time_ref
        for next_label_ref in next_labels_ref:
            if next_label_ref in new_correspondences.keys():
                if verbose:
                    print "relabel ref %04d" % next_label_ref
                # correspondences = new_correspondences
                point_cloud_ref, point_cloud_flo, daughters_ref, daughters_flo = ml.lineages_relabelled_correspondences_propagation(
                    lineage_ref, lineage_flo,
					new_correspondences,
                    next_label_ref)
                if daughters_ref and daughters_flo:
                    if len(daughters_ref) == 2 and (
                            new_correspondences.has_key(daughters_ref[0]) or new_correspondences.has_key(
                            daughters_ref[1])):
                        if verbose:
                            print "-> daughters already in the correspondence dictionary"
                    else:
                        daughters_correspondences, scalars[next_label_ref] = sisters_association(point_cloud_ref,
																								 point_cloud_flo,
																								 new_correspondences,
																								 daughters_ref,
																								 daughters_flo,
																								 threshold_angle=threshold_angle)
                        new_correspondences.update(daughters_correspondences)
                        if not daughters_correspondences:
                            unpropagated_cells_ref.append(next_label_ref)
                            if verbose:
                                print "-> daughters division axes form an angle %.1f which is higher than the tolerance" % \
                                      scalars[next_label_ref][0]
                            if stop_when_blocked:
                                return new_correspondences, unpropagated_cells_ref, scalars

                else:
                    if verbose:
                        print "-> end of lineage correspondences"
            else:
                if verbose:
                    print "-> skipping relabel ref %04d (not in the correspondence dictionary)" % next_label_ref
        next_time_ref = next_time_ref + 1

    return new_correspondences, unpropagated_cells_ref, scalars


def temporal_affine_registration(lineage_ref, lineage_flo, correspondences, time_ref_min=None, time_ref_max=None,
                                 time_flo_min=None, time_flo_max=None, weighting_param=None):
    """
	Affine temporal registration of a floating lineage lineage_flo onto a reference lineage lineage_ref
	with respect to a cell-to-cell correspondence dictionary correspondences.
	The registration is processed with a basic linear regression of data given by corresponding cell death time-points.
	Inputs:
		lineage_ref
		lineage_flo
		correspondences
	Optional:
		time_ref_min
		time_ref_max
		time_flo_min
		time_flo_max
			-> in order to specify the delimiters
		weighting_param
			-> in order to specify weighing parameters (those for the function <weighting_fun_ncells>) ; default = None
	Outputs:
		a, b : providing the relation t_flo = a * t_ref + b
	"""
    if time_ref_min == None:
        time_ref_min = lineage_ref.timepoints()[0]
    if time_flo_min == None:
        time_flo_min = lineage_flo.timepoints()[0]
    if time_ref_max == None:
        time_ref_max = lineage_ref.timepoints()[-1]
    if time_flo_max == None:
        time_flo_max = lineage_flo.timepoints()[-1]

    times_death_ref = []
    times_death_flo = []
    ncells_at_death_ref = []
    ncells_at_death_flo = []

    for k, v in correspondences.items():
        if not lineage_ref.exists(k):
            print "Warning : relabel %d not found in lineage_ref." % k
        else:
            if not lineage_flo.exists(v):
                print "Warning : relabel %d not found in lineage_flo." % v
            else:
                t_ref = lineage_ref.relabelled_death(k)
                t_flo = lineage_flo.relabelled_death(v)
                if (t_ref < time_ref_max and t_flo < time_flo_max) and (
                        t_ref >= time_ref_min and t_flo >= time_flo_min):
                    times_death_ref.append(t_ref)
					times_death_flo.append(t_flo)
                    ncells_at_death_ref.append(lineage_ref.ncells(t_ref))
                    ncells_at_death_flo.append(lineage_flo.ncells(t_flo))

    X = times_death_ref
    Y = times_death_flo

    # Linear regression on (X,Y)
    # -> find (a,b) that minimize square error given by sum(|y-(ax+b)|^2) for

    if weighting_param == None or weighting_param == False:
        c = np.cov(X, Y, bias=False)
        a = c[0, 1] / c[0, 0];
        b = np.mean(Y) - a * np.mean(X)
        return a, b
    W = []
    if type(weighting_param) is not tuple:
        W = [weighting_fun_ncells((ncells_at_death_ref[i], ncells_at_death_flo[i]))[-1] for i in
             range(len(ncells_at_death_ref))]
    else:
        # print "weighting_param = %s"% str(weighting_param)
        W = [weighting_fun_ncells((ncells_at_death_ref[i], ncells_at_death_flo[i], weighting_param))[-1] for i in
             range(len(ncells_at_death_ref))]
    M = np.zeros((2, 2))
    N = np.zeros((2, 1))
    for i in range(len(X)):
        M[0, 0] = M[0, 0] + W[i] * X[i] ** 2
		M[0, 1] = M[0, 1] + W[i] * X[i]
        N[0] = N[0] + W[i] * X[i] * Y[i]
        N[1] = N[1] + W[i] * Y[i]
    M[1, 0] = M[0, 1]
    M[1, 1] = sum(W)

    D = np.dot(np.linalg.inv(M), N)
    return D[0, 0], D[1, 0]


def temporal_affine_registration_robust(lineage_ref, lineage_flo, correspondences, init_ref_min=20, init_flo_min=10):
    """
	NOT TO BE USED
	"""
    time_ref_first = lineage_ref.timepoints()[0]
    time_ref_last = lineage_ref.timepoints()[-1]
    time_flo_first = lineage_flo.timepoints()[0]
    time_flo_last = lineage_flo.timepoints()[-1]

    a = {}
    b = {}
    for t_ref in lineage_ref.timepoints()[init_ref_min:]:
        for t_flo in lineage_flo.timepoints()[init_flo_min:]:
            a[(t_ref, t_flo)], b[(t_ref, t_flo)] = temporal_affine_registration(lineage_ref, lineage_flo,
                                                                                correspondences, time_ref_max=t_ref,
                                                                                time_flo_max=t_flo)
    return a, b


def volumes_variation_study_for_pairings_validity(lineage_ref, lineage_flo, correspondences, label_ref,
                                                  field_volume='real_volume'):
    '''

	'''
    assert correspondences.has_key(label_ref)
    label_flo = correspondences[label_ref]
    daughters_ref = lineage_ref.relabelled_daughters(label_ref)
    if len(daughters_ref) != 2:
        print "Not the expected number of daughters for reference relabelled cell %d." % label_ref
        return None
    daughter_1_ref = daughters_ref[0]
    daughter_2_ref = daughters_ref[1]
    assert correspondences.has_key(daughter_1_ref)
    assert correspondences.has_key(daughter_2_ref)
    daughter_1_flo = correspondences[daughter_1_ref]
    daughter_2_flo = correspondences[daughter_2_ref]

    volumes_label_ref = lineage_ref.relabelled_volumes(label_ref, field_volume=field_volume)
    volumes_label_flo = lineage_flo.relabelled_volumes(label_flo, field_volume=field_volume)
    ratio_death_birth_label_ref = lineage_ref.relabelled_growth_ratio(label_ref,
                                                                      field_volume=field_volume)  # Must be close to 1
    ratio_death_birth_label_flo = lineage_flo.relabelled_growth_ratio(label_flo,
                                                                      field_volume=field_volume)  # Must be close to 1

    volume_daughter_1_ref = lineage_ref.relabelled_volume_at_time(daughter_1_ref,
                                                                  lineage_ref.relabelled_birth(daughter_1_ref),
                                                                  field_volume=field_volume)
    volume_daughter_2_ref = lineage_ref.relabelled_volume_at_time(daughter_2_ref,
                                                                  lineage_ref.relabelled_birth(daughter_2_ref),
                                                                  field_volume=field_volume)
    volume_daughter_1_flo = lineage_flo.relabelled_volume_at_time(daughter_1_flo,
                                                                  lineage_flo.relabelled_birth(daughter_1_flo),
                                                                  field_volume=field_volume)
    volume_daughter_2_flo = lineage_flo.relabelled_volume_at_time(daughter_2_flo,
                                                                  lineage_flo.relabelled_birth(daughter_2_flo),
                                                                  field_volume=field_volume)

    ratio_daughters_mother_ref = (volume_daughter_1_ref + volume_daughter_2_ref) / volumes_label_ref[
        lineage_ref.relabelled_death(label_ref)]  # Must be close to 1
    ratio_daughters_mother_flo = (volume_daughter_1_flo + volume_daughter_2_flo) / volumes_label_flo[
        lineage_flo.relabelled_death(label_flo)]  # Must be close to 1

    ratio_d_1_daughters_ref = volume_daughter_1_ref / (volume_daughter_1_ref + volume_daughter_2_ref)
    ratio_d_1_daughters_flo = volume_daughter_1_flo / (volume_daughter_1_flo + volume_daughter_2_flo)
    # Must have ratio_d_1_daughters_ref ~ ratio_d_1_daughters_flo and (1-ratio_d_1_daughters_ref) close to (1-ratio_d_1_daughters_flo)
    squared_ratio_difference = (ratio_d_1_daughters_ref - ratio_d_1_daughters_flo) ** 2  # Must be close to 0

    D = {
        'volumes_label_ref': volumes_label_ref,
        'volumes_label_flo': volumes_label_flo,
        'ratio_death_birth_mother_ref': ratio_death_birth_label_ref,  # Must be close to 1
        'ratio_death_birth_mother_flo': ratio_death_birth_label_flo,  # Must be close to 1
        'volume_daughter_1_ref': volume_daughter_1_ref,
        'volume_daughter_2_ref': volume_daughter_2_ref,
        'volume_daughter_1_flo': volume_daughter_1_flo,
        'volume_daughter_2_flo': volume_daughter_2_flo,
        'ratio_daughters_mother_ref': ratio_daughters_mother_ref,  # Must be close to 1
        'ratio_daughters_mother_flo': ratio_daughters_mother_flo,  # Must be close to 1
        'ratio_d_1_daughters_ref': ratio_d_1_daughters_ref,
        'ratio_d_1_daughters_flo': ratio_d_1_daughters_flo,
        'squared_ratio_difference': squared_ratio_difference  # Should be close to 0
        # 'grand_mother':
    }
    return D


def temporal_registration_distance_to_regression_for_pairings_validity(lineage_ref, lineage_flo, correspondences,
																	   label_ref, a=None, b=None,
																	   include_end_of_time=False, predictive=False):
	"""

    """
	assert correspondences.has_key(label_ref)
    label_flo = correspondences[label_ref]

    if a == None or b == None:
        a, b = temporal_affine_registration(lineage_ref, lineage_flo, correspondences,
                                            include_end_of_time=include_end_of_time)

    # Alignement temporel : t_flo = a * t_ref + b

    d_r = lineage_ref.relabelled_death(label_ref)  # death ref
    # b_r=lineage_ref.relabelled_birth(label_ref) # birth ref
    d_f = lineage_flo.relabelled_death(label_flo)  # death flo
    # b_f=lineage_flo.relabelled_birth(label_flo) # birth flo

    # d=0
    if predictive:
        # return b_f - a*b_r - b, d_f - a*d_r - b
        return d_f - a * d_r - b
    # denom=np.sqrt(a**2+1)
    # return np.abs(a*b_r-b_f+b)/denom, np.abs(a*d_r-d_f+b)/denom
    return np.abs(a * d_r - d_f + b) / np.sqrt(a ** 2 + 1)


def temporal_registration_weighting(weighting_mode=None, parameters=None):
    """
	weighting mode: function to be applied for the given pair of labels. If None, returns 1.
	opt: options for the weighting_mode function (tuple)

	"""
    if weighting_mode == None:
        return 1
    if type(weighting_mode) == str:
        # TODO
        return 1
    else:
        # weighting_mode is a function
        return weighting_mode(parameters)


def weighting_fun_ncells(parameters):
    """
	Function (a,b,c,d) = weighting_fun_ncells(parameters)
	    parameters = ( n, m )
	    parameters = ( lineage_ref, lineage_flo, label_ref, [ label_flo | correspondences ] )
	    parameters = ( ... , ('alpha', alpha))
	    parameters = ( ... , ('func_alpha', func_alpha))
	    parameters = ( ... , ('func_alpha_normalized', func_alpha_normalized))
	    parameters = ( ... , ('beta', beta))
	    parameters = ( ... , ('func_beta', func_beta))
	    parameters = ( ... , ('func', func))
	Note: the mandatory parameters can be given using tuples of length 2 with first value being the
	name of the parameter as well as the optional parameters.
	E.g. ('label_ref', 50 ) is equivalent to label_ref = 50
	Input parameters description:
		n: int > 0
		m: int > 0
		lineage_ref: morpheme_lineage.Morpheme_lineage instance
		lineage_flo: morpheme_lineage.Morpheme_lineage instance
		label_ref: int
		label_flo: int (if "correspondences" not specified)
		correspondences: dict (if "label_flo" not specified)
		alpha: float (=1 by default)
		beta: float (=1 by default)
		func_alpha_normalized: boolean (default=True)
		func_alpha, func_beta, func: <function>
			'exponential': lambda x, coef=1: numpy.exp(-float(x)*coef), # (default for func_alpha)
			'inverse':lambda x,coef=1: 1/(float(x)*coef),
			'inverse_sqrt':lambda x,coef=1: 1/numpy.sqrt(float(x)*coef),# (default for func_beta)
			'inverse_square':lambda x,coef=1: 1/(float(x)*coef)**2,
			'inverse_1':lambda x,coef=1: 1/(1+float(x)*coef),
			'identity':lambda x,coef=1: x*coef),
			'one':lambda x,coef=1: 1
			(where coef = alpha | beta respectively for func_alpha | func_beta)
			Or any custom function taking one mandatory and one optional (default=1) parameters
		Note: ('func', func) is a shortcut for ('func_alpha', func), ('func_beta',func).
	Output values description:
		a = func_alpha( abs(n - m) ) if not func_alpha_normalized or func_alpha( 2. * abs(n - m) / (n + m) ) if func_alpha_normalized
		b = func_beta( n )
		c = func_beta( m )
		d = a*b*c
	E.g. of use:
		a,b,c,d=weighting_fun_ncells([lineage_ref, lineage_flo, 2, correspondences, ('alpha', 10), ('beta', 1./64), ('func_alpha','exponential'), ('func_alpha_normalized', True)])
	"""
    program = "weighting_fun_ncells"
    import numpy
    default_lambdas = {
		'exponential': lambda x, coef=1: numpy.exp(-float(x) * coef),
        'inverse': lambda x, coef=1: 1 / (float(x) * coef),
        'inverse_sqrt': lambda x, coef=1: 1 / numpy.sqrt(float(x) * coef),
        'inverse_square': lambda x, coef=1: 1 / (float(x) * coef) ** 2,
        'inverse_1': lambda x, coef=1: 1 / (1 + float(x) * coef),
        'identity': lambda x, coef=1: float(x) * coef,
        'one': lambda x, coef=1: 1
    }

    # headers=dict(parameters)
    # Parameters parsing
    assert len(parameters) >= 1, "%s: Unexpected number of parameters. See function help." % program
    # Mandatory parameters
    ncell_ref = None
    ncell_flo = None
    lineage_ref = None
    lineage_flo = None
    label_ref = None
    label_flo = None
    correspondences = None
    for i in range(min(len(parameters), 4)):
        if type(parameters[i]) == tuple:
            break
        if i == 0:
            ncell_ref = parameters[i]
        if i == 1:
            ncell_flo = parameters[i]
        if i == 2:
            lineage_ref = ncell_ref
            lineage_flo = ncell_flo
            label_ref = parameters[i]
            # the two first arguments correspond to the lineages
			# assert not type(label)
			ncell_ref = None
			ncell_flo = None
		if i == 3:
            if type(parameters[i]) == dict:
                correspondences = parameters[i]
            else:
                assert type(label_ref) == type(
					parameters[i]), "%s: Unexpected type for the %d-th mandatory parameter." % (program, i + 1)
                label_flo = parameters[i]
            i = i + 1

    assert len(parameters) == 2 or i == 2 or i == 4, "%s: Unexpected number of parameters. See function help." % program

    headers = {}
    if i == 2 or i == 4:
        headers = dict(parameters[i:])
    # print "headers=%s"%str(headers)

    # Mandatory parameters
    if headers.has_key('n'):
        assert ncell_ref == None, "%s: Multiple definition of 'n'." % program
        assert lineage_ref == None and lineage_flo == None and label_ref == None and label_flo == None and correspondences == None, "%s: The {'n', 'm'} and {'lineage_ref', 'lineage_flo', 'label_ref', 'label_flo' | 'correspondences'}  parameters cannot be defined simutaneously." % program
        ncell_ref = int(headers['n'])
        del headers['n']
    if headers.has_key('m'):
        assert ncell_flo == None, "%s: Multiple definition of 'm'." % program
        assert lineage_ref == None and lineage_flo == None and label_ref == None and label_flo == None and correspondences == None, "%s: The {'n', 'm'} and {'lineage_ref', 'lineage_flo', 'label_ref', 'label_flo' | 'correspondences'}  parameters cannot be defined simutaneously." % program
        ncell_flo = int(headers['m'])
        del headers['m']
    if headers.has_key('lineage_ref'):
        assert lineage_ref == None, "%s: Multiple definition of 'lineage_ref'." % program
        assert i == 0, "%s: Multiple definition of 'lineage_ref' (see help)." % program
        lineage_ref = headers['lineage_ref']
        del headers['lineage_ref']
    if headers.has_key('lineage_flo'):
        assert lineage_flo == None, "%s: Multiple definition of 'lineage_flo'." % program
        assert i < 2, "%s: Multiple definition of 'lineage_flo' (see help)." % program
        lineage_flo = headers['lineage_flo']
        del headers['lineage_flo']
    if headers.has_key('label_ref'):
        assert label_ref == None, "%s: Multiple definition of 'label_ref'." % program
        if i > 0:
            lineage_ref = ncell_ref
            ncell_ref = None
            if i == 2:
                lineage_flo = ncell_flo
                ncell_flo = None
        label_ref = int(headers['label_ref'])
        del headers['label_ref']
    if headers.has_key('label_flo'):
        assert label_flo == None, "%s: Multiple definition of 'label_flo'." % program
        label_flo = int(headers['label_flo'])
        del headers['label_flo']
    if headers.has_key('correspondences'):
        assert correspondences == None, "%s: Multiple definition of 'correspondences'." % program
        correspondences = headers['correspondences']
        del headers['correspondences']

    assert ncell_ref != None or lineage_ref != None
    assert ncell_flo != None or lineage_flo != None
    assert ncell_ref != None or label_ref != None
    if correspondences != None:
        assert ncell_flo != None or label_flo == None, "%s: The 'correspondences' and 'label_flo' parameters cannot be defined simultaneously." % program
    if label_flo == None:
        assert ncell_flo != None or correspondences != None, "%s: Missing 'label_flo' or 'correspondences' parameter (see documentation)." % program
        if ncell_flo == None:
            assert type(
                correspondences) == dict, "%s: 'correspondences' parameter of wrong type '%s' (expected type 'dict') ." % (
            program, str(type(correspondences)))
            assert correspondences.has_key(
                label_ref), "%s: 'correspondences' does not have the key corresponding to label_ref %d." % (
            program, label_ref)
            label_flo = correspondences[label_ref]
            assert label_flo != None

    # Optional parameters default values
    func_alpha_normalized = True

    # Optional parameters
    if headers.has_key('func_alpha_normalized'):
        func_alpha_normalized = bool(headers['func_alpha_normalized'])
        del headers['func_alpha_normalized']

    # Optional parameters default values
    alpha = 1
    beta = 1

    # Optional parameters
    if headers.has_key('alpha'):
        alpha = float(headers['alpha'])
        del headers['alpha']
    if headers.has_key('beta'):
        beta = float(headers['beta'])
        del headers['beta']

    # Optional parameters default values
    func_alpha = default_lambdas['exponential']
    func_beta = default_lambdas['inverse_sqrt']

    # Optional parameters
    if headers.has_key('func'):
        assert not headers.has_key(
            'func_alpha'), "%s: The 'func' parameter cannot be defined simultaneously with 'func_alpha' or 'func_beta'." % program
        assert not headers.has_key(
            'func_beta'), "%s: The 'func' parameter cannot be defined simultaneously with 'func_alpha' or 'func_beta'." % program
        if default_lambdas.has_key(headers['func']):
            func_alpha = default_lambdas[headers['func']]
            func_beta = default_lambdas[headers['func']]
        else:
            assert type(headers['func']) != str, "%s: Unknown function name %s for 'func' option." % (
            program, headers['func'])
            func_alpha = headers['func']
            func_beta = headers['func']
        del headers['func']

    if headers.has_key('func_alpha'):
        if default_lambdas.has_key(headers['func_alpha']):
            func_alpha = default_lambdas[headers['func_alpha']]
        else:
            assert type(headers['func_alpha']) != str, "%s: Unknown function name %s for 'func_alpha' option." % (
            program, headers['func_alpha'])
            func_alpha = headers['func_alpha']
        del headers['func_alpha']

    if headers.has_key('func_beta'):
        if default_lambdas.has_key(headers['func_beta']):
            func_beta = default_lambdas[headers['func_beta']]
        else:
            assert type(headers['func_beta']) != str, "%s: Unknown function name %s for 'func_beta' option." % (
            program, headers['func_beta'])
            func_beta = headers['func_beta']
        del headers['func_beta']

    assert len(headers) == 0, "%s: Unexpected optional parameters (unknown or duplicated ones): %s" % (
    program, str(headers))

    # Stuff
    if ncell_ref == None:
        assert lineage_ref.exists(label_ref), "%s: Specified label_ref '%d' not found in lineage_ref." % (
        program, int(label_ref))
        t_ref = lineage_ref.death(label_ref)
        ncell_ref = lineage_ref.ncells(t_ref)
    if ncell_flo == None:
        assert lineage_flo.exists(label_flo), "%s: Specified label_flo '%d' not found in lineage_flo." % (
        program, int(label_flo))
        t_flo = lineage_flo.death(label_flo)
        ncell_flo = lineage_flo.ncells(t_flo)

    # print "ncell_ref = %d ; ncell_flo = %d" % (ncell_ref, ncell_flo)
    a = 0
    if func_alpha_normalized:
        a = func_alpha(2. * abs(ncell_ref - ncell_flo) / (ncell_ref + ncell_flo), alpha)
    else:
        a = func_alpha(abs(ncell_ref - ncell_flo), alpha)
    b = func_beta(ncell_ref, beta)
    c = func_beta(ncell_flo, beta)

    return (a, b, c, a * b * c)


###################################
### VISU TEMPORAL REGISTRATION ####
###################################

def generate_registered_sequences(lin_ref, lin_flo, correspondences, func_temporal_reg,
                                  format_images_in_ref, format_images_in_flo,
                                  format_images_out_ref, format_images_out_flo, format_trsfs_ref=None,
                                  template_image_ref=None, iso=None, dimensions=None, voxelsize=None,
                                  backgroundLabels=[0, 1],
                                  visu=True, verbose=False):
    '''


	'''
    from lineage import timeNamed, timesNamed

    path_images_ref_in = os.path.sep.join(format_images_in_ref.split(os.path.sep)[:-1])
    path_images_flo_in = os.path.sep.join(format_images_in_flo.split(os.path.sep)[:-1])
    path_images_ref_out = os.path.sep.join(format_images_out_ref.split(os.path.sep)[:-1])
    path_images_flo_out = os.path.sep.join(format_images_out_flo.split(os.path.sep)[:-1])

    if path_images_ref_in == '':
        path_images_ref_in = os.path.curdir
    if path_images_flo_in == '':
        path_images_flo_in = os.path.curdir
    if path_images_ref_out == '':
        path_images_ref_out = os.path.curdir
    if path_images_flo_out == '':
        path_images_flo_out = os.path.curdir

    assert format_images_out_ref.count('$TIME') == 1 or (
                format_images_out_ref.count('$TIMEREF') == 1 and format_images_out_ref.count('$TIMEFLO') == 1)
    assert format_images_out_flo.count('$TIME') == 1 or (
			format_images_out_flo.count('$TIMEREF') == 1 and format_images_out_flo.count('$TIMEFLO') == 1)

	assert os.path.isdir(path_images_ref_in)

	assert os.path.isdir(path_images_flo_in)

	if not os.path.isdir(path_images_ref_out):
		os.mkdir(path_images_ref_out)
    assert os.path.isdir(path_images_ref_out)

    if not os.path.isdir(path_images_flo_out):
        os.mkdir(path_images_flo_out)
    assert os.path.isdir(path_images_flo_out)

    if format_trsfs_ref:
        path_trsfs_ref = os.path.sep.join(format_trsfs_ref.split(os.path.sep)[:-1])
        if path_trsfs_ref == '':
            path_trsfs_ref = os.path.curdir
        assert os.path.isdir(path_trsfs_ref)

    if func_temporal_reg is not None:
        if type(func_temporal_reg) is list:
            assert len(func_temporal_reg) == 2
            # func_temporal_reg = [a, b] where a and b are the values of the coefficients for the affine function t_flo = a * t_ref + b
            func_temporal_reg = lambda x, a=func_temporal_reg[0], b=func_temporal_reg[1]: a * x + b
    # else:
    #	# Here we assume that func_temporal_reg is a function that takes one numeric argument and returns a numeric value
    #	x=range(t_beg_ref,t_end_ref+1)
    #	y=[func_temporal_reg(i) for i in x]

    for time_ref in lin_ref.timepoints():
        time_flo = int(func_temporal_reg(time_ref))
		if time_flo < lin_flo.timepoints()[0]:
			# time_flo=lin_flo.timepoints()[0]
			continue
		if time_flo > lin_flo.timepoints()[-1]:
			# time_flo=lin_flo.timepoints()[-1]
			continue

        path_ref_in_at_time = timeNamed(format_images_in_ref, time_ref)
        path_flo_in_at_time = timeNamed(format_images_in_flo, time_flo)
        assert os.path.exists(path_ref_in_at_time), "Input image file '%s' not found." % path_ref_in_at_time
		assert os.path.exists(path_flo_in_at_time), "Input image file '%s' not found." % path_flo_in_at_time
        path_trsf_at_time = None
        if format_trsfs_ref:
            path_trsf_at_time = timeNamed(format_trsfs_ref, time_ref)
            assert os.path.exists(path_trsf_at_time), "Input transformation file '%s' not found." % path_trsf_at_time
        if verbose:
            print "TIME_REF = t%03d, TIME_FLO = t%03d" % (time_ref, time_flo)

        # path_ref_out=timeNamed(format_images_out_ref,time_ref)
        path_ref_out = ''
        path_flo_out = ''
        if format_images_out_ref.count('$TIME') == 1:
            path_ref_out = timeNamed(format_images_out_ref, time_ref)
        else:
            path_ref_out = timesNamed(format_images_out_ref, '$TIMEREF', time_ref, '$TIMEFLO', time_flo)
        if format_images_out_flo.count('$TIME') == 1:
            path_flo_out = timeNamed(format_images_out_flo, time_ref)
        else:
            path_flo_out = timesNamed(format_images_out_flo, '$TIMEREF', time_ref, '$TIMEFLO', time_flo)

        segmentationRelabellingAtTime(path_ref_in_at_time, path_ref_out, lin_ref, time_ref, lut=correspondences.keys(),
                                      trsf=path_trsf_at_time, template=template_image_ref,
                                      voxelsize=voxelsize, iso=iso, dimensions=dimensions,
                                      backgroundLabels=backgroundLabels, visu=visu, verbose=verbose)

        assert os.path.exists(path_ref_out), "Failed to generate output file '%s'" % path_ref_out

        # Calcul de la transformation
        out = pointCloudRegistration(lin_ref.relabelled_barycenters_at_time(time_ref),
									 lin_flo.relabelled_barycenters_at_time(time_flo),
									 correspondences, trsf_type='rigid', skip_not_found=True, lazy=False,
									 verbose=verbose)
		trsf_f_r = out['trsf']
        # path_trsf_r_template = "/home/gmicheli/DIG-EM/Data/140317-Patrick-St8/FUSE/REG/140317-Patrick-St8_reg_compose_t%03d.trsf"%time_ref
        path_trsf_tmp = path_images_flo_out + os.path.sep + "generate_registered_sequences_tmp_path.trsf"
        f = open(path_trsf_tmp, "w")
        np.savetxt(f, trsf_f_r, '%f')
        # f.write(str(np.array(trsf_f_r)).replace('[','').replace(']',''))
        f.close()
        # Composition de la transformation avec la transformation de la reference vers le template
        if path_trsf_at_time:
            compose_trsf(path_trsf_tmp, path_trsf_at_time, path_trsf_tmp, verbose=verbose)

        segmentationRelabellingAtTime(path_flo_in_at_time, path_flo_out, lin_flo, time_flo, reverse_lut=correspondences,
                                      trsf=path_trsf_tmp, template=path_ref_out,
                                      voxelsize=None, iso=None, dimensions=None, backgroundLabels=backgroundLabels,
                                      visu=visu, verbose=verbose)

        assert os.path.exists(path_flo_out), "Failed to generate output file '%s'" % path_flo_out

        if os.path.exists(path_trsf_tmp):
            cmd = 'rm -f %s' % path_trsf_tmp
            if verbose:
                print cmd
            os.system(cmd)


def plot_temporal_correspondences(lin_ref, lin_flo, correspondences, func_temporal_reg=None, marker_cloud=None,
                                  label_cloud=None, marker_reg=None, label_reg=None, block=0, loc=4, xlabel=None,
                                  ylabel=None, title=None):
    '''
	Function to plot a cloud of temporal correspondences between division time-points of corresponding cells from two developing ascidian embryos.

	Inputs:
		marker_cloud="r+" or whatever instead
		label_cloud='Mapped cell division times' or whatever instead
		xlabel='Reference sample timespan' or whatever instead
		ylabel='Floating sample timespan' or whatever instead
		title='Embryo temporal alignment' or whatever instead

		func_temporal_reg = [a, b] where a and b are the values of the coefficients for the affine function t_flo = a * t_ref + b
		marker_reg='b--'

		loc=4 -> plt.legend(loc=loc)
		block=0 -> plt.show(block=block)

	'''
    import matplotlib.pyplot as plt

    t_beg_ref = lin_ref.timepoints()[0]
    t_beg_flo = lin_flo.timepoints()[0]

    t_end_ref = lin_ref.timepoints()[-1]
    t_end_flo = lin_flo.timepoints()[-1]

    X_tot = [lin_ref.relabelled_death(r) for r in correspondences.keys()]
    Y_tot = [lin_flo.relabelled_death(r) for r in correspondences.values()]

    X = [X_tot[i] for i in range(len(X_tot)) if (X_tot[i] != t_end_ref and Y_tot[i] != t_end_flo)]
    Y = [Y_tot[i] for i in range(len(Y_tot)) if (X_tot[i] != t_end_ref and Y_tot[i] != t_end_flo)]

    x = []
    y = []
    if func_temporal_reg is not None:
        if type(func_temporal_reg) is list:
            assert len(func_temporal_reg) == 2
            # func_temporal_reg = [a, b] where a and b are the values of the coefficients for the affine function t_flo = a * t_ref + b
            x = [t_beg_ref, t_end_ref]
            y = [func_temporal_reg[0] * i + func_temporal_reg[1] for i in x]

        else:
            # Here we assume that func_temporal_reg is a function that takes one numeric argument and returns a numeric value
            x = range(t_beg_ref, t_end_ref + 1)
            y = [func_temporal_reg(i) for i in x]

    if not marker_cloud:
        marker_cloud = "r+"
    if not label_cloud:
        label_cloud = 'Mapped cell division times'

    plt.plot(X, Y, marker_cloud, label=label_cloud)
    if x:
        if not marker_reg:
            marker_reg = 'b--'
        if not label_reg:
            label_reg = "Linear regression"
        plt.plot(x, y, marker_reg, label=label_reg)
    if xlabel is None:
        xlabel = 'Reference sample timespan'
    if ylabel is None:
        ylabel = 'Floating sample timespan'
    if title is None:
        title = 'Embryo temporal alignment'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.legend(loc=loc)
    plt.axis([t_beg_ref, t_end_ref, t_beg_flo, t_end_flo])
    plt.show(block=block)


def plot_descent_lineage_with_temporal_correspondences(lin_ref, lin_flo, correspondences, relabel_ref,
                                                       func_temporal_reg=None, weights=None, s=70, cmap=None,
                                                       flt_format=None, loc=4, xlabel=None, ylabel=None, title=None,
                                                       block=0):
    '''
	weights : None or dict containing the keys one wants to plot with values which are the weights the user wants to plot.
			If None, the the weights will be automatically set to the cell level in the lineage to be plotted (1, 2, 3, ...).
	'''
    import matplotlib.pyplot as plt

    if xlabel is None:
        xlabel = 'Reference sample timespan'
    if ylabel is None:
        ylabel = 'Floating sample timespan'
    if title is None:
        title = 'Lineage divisions alignments from cell %d' % relabel_ref

    t_beg_ref = lin_ref.timepoints()[0]
    t_beg_flo = lin_flo.timepoints()[0]

    t_end_ref = lin_ref.timepoints()[-1]
    t_end_flo = lin_flo.timepoints()[-1]

    assert lin_ref.exists(relabel_ref)
    assert correspondences.has_key(relabel_ref)
    relabels_ref = [relabel_ref] + lin_ref.relabelled_descendents(relabel_ref)

    X = [lin_ref.relabelled_death(r) for r in relabels_ref if r in correspondences.keys()]
    Y = [lin_flo.relabelled_death(correspondences[r]) for r in relabels_ref if r in correspondences.keys()]

    c = []

    for r in relabels_ref:
        if not r in correspondences.keys():
            continue
        if weights:
            assert weights.has_key(r), "No entry for reference relabel %d" % r
            c = c + [weights[r]]
        else:
            r0 = r
            i = 1
            while r0 is not relabel_ref:
                i = i + 1
                r0 = lin_ref.relabelled_mother(r0)
                assert r0 is not None, "Unexpectedly unable to find the descent level of label %d in the sub-lineage of label %d " % (
                r, relabel_ref)
            c = c + [i]

    assert len(c) == len(X)

    the_Xs = []
    the_Ys = []
    for r in relabels_ref:
        if r is relabel_ref:
            continue
        if not r in correspondences.keys():
            continue
        if lin_ref.relabelled_death(r) == lin_ref.timepoints()[-1] and lin_flo.relabelled_death(correspondences[r]) == \
                lin_flo.timepoints()[-1]:
            continue
        the_Xs.append([lin_ref.relabelled_death(r)])
        the_Ys.append([lin_flo.relabelled_death(correspondences[r])])
        r0 = lin_ref.relabelled_mother(r)
        assert r0 is not None, "Unexpectedly unable to find the mother of label %d in the sub-lineage of label %d " % (
        r, relabel_ref)
        while not r0 in correspondences.keys():
            r1 = lin_ref.relabelled_mother(r0)
            assert r1 is not None, "Unexpectedly unable to find the mother of label %d in the sub-lineage of label %d " % (
            r0, relabel_ref)
            r0 = r1
        the_Xs[-1].append(lin_ref.relabelled_death(r0))
        the_Ys[-1].append(lin_flo.relabelled_death(correspondences[r0]))

    # print the_Xs
    # print the_Ys

    x = []
    y = []
    if func_temporal_reg is not None:
        if type(func_temporal_reg) is list:
            assert len(func_temporal_reg) == 2
            # func_temporal_reg = [a, b] where a and b are the values of the coefficients for the affine function t_flo = a * t_ref + b
            x = [t_beg_ref, t_end_ref]
            y = [func_temporal_reg[0] * i + func_temporal_reg[1] for i in x]

        else:
            # Here we assume that func_temporal_reg is a function that takes one numeric argument and returns a numeric value
            x = range(t_beg_ref, t_end_ref + 1)
            y = [func_temporal_reg(i) for i in x]

    if cmap is None:
        cmap = plt.cm.gnuplot2

    for i in range(len(the_Xs)):
        the_x = the_Xs[i]
        the_y = the_Ys[i]
        plt.plot(the_x, the_y, 'k')

    cax = plt.scatter(X, Y, c=c, cmap=cmap, s=s, marker='o')
    if x:
        plt.plot(x, y, marker_reg, label=label_reg)
    plt.axis([t_beg_ref, t_end_ref, t_beg_flo, t_end_flo])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc=loc)
    plt.title(title)

    if weights:
        if angle:
            cbar = plt.colorbar(cax, ticks=[1, 10, 20, 30, 40, 50, 60, 70, 80, 87.5], orientation='vertical')
            cbar.set_ticklabels([str(angle) + ' degrees'.decode("utf8") for angle in
                                 [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]])  # horizontal colorbar
        else:
            import numpy as np
            cbar = plt.colorbar(cax, ticks=np.linspace(min(c), max(c), 5), orientation='vertical')
            if not flt_format:
                flt_format = "%f"
            cbar.set_ticklabels([flt_format % v for v in np.linspace(min(c), max(c), 5)])  # horizontal colorbar
    else:
        cbar = plt.colorbar(cax, ticks=range(1, max(c) + 1), orientation='vertical')
        cbar.set_ticklabels(["Level " + str(v) for v in range(1, max(c) + 1)])  # horizontal colorbar
    plt.show(block=block)


def plot_angles_histogram(angles, xlabel=None, ylabel=None, title=None, bins=None, normed=True, block=0):
    """

	"""
    import matplotlib.pyplot as plt

    # Histogrammes des angles
    if bins is None:
        bins = [0.5 + _x for _x in range(0, 90, 2)]
    plt.hist(angles, bins=bins, normed=normed)
    if xlabel is None:
        xlabel = 'Angle (degrees)'
    if ylabel is None:
        ylabel = 'Proportion of cell divisions'
    if title is None:
        title = 'Histogram of angles between corresponding cells division axes (Ref/Flo)'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show(block=block)


def scatter_temporal_correspondences(lin_ref, lin_flo, correspondences, func_temporal_reg=None, weights=None,
                                     marker_cloud=None, label_cloud=None, marker_reg=None,
                                     s=70, label_reg=None, block=0, loc=4, bins=None, xlabel=None, ylabel=None,
                                     title=None, cmap=None, c=None):
    '''

	'''
    import matplotlib.pyplot as plt

    t_beg_ref = lin_ref.timepoints()[0]
    t_beg_flo = lin_flo.timepoints()[0]

    t_end_ref = lin_ref.timepoints()[-1]
    t_end_flo = lin_flo.timepoints()[-1]

    X_tot = [lin_ref.relabelled_death(r) for r in correspondences.keys()]
    Y_tot = [lin_flo.relabelled_death(r) for r in correspondences.values()]

    # X90=[X90_tot[i] for i in range(len(X90_tot)) if (X90_tot[i] != t_end_ref and Y90_tot[i] != t_end_flo)]
    # Y90=[Y90_tot[i] for i in range(len(Y90_tot)) if (X90_tot[i] != t_end_ref and Y90_tot[i] != t_end_flo)]
    X = [X_tot[i] for i in range(len(X_tot)) if (X_tot[i] != t_end_ref and Y_tot[i] != t_end_flo)]
    Y = [Y_tot[i] for i in range(len(Y_tot)) if (X_tot[i] != t_end_ref and Y_tot[i] != t_end_flo)]

    x = []
    y = []
    if func_temporal_reg is not None:
        if type(func_temporal_reg) is list:
			assert len(func_temporal_reg) == 2
            # func_temporal_reg = [a, b] where a and b are the values of the coefficients for the affine function t_flo = a * t_ref + b
            x = [t_beg_ref, t_end_ref]
            y = [func_temporal_reg[0] * i + func_temporal_reg[1] for i in x]

        else:
            # Here we assume that func_temporal_reg is a function that takes one numeric argument and returns a numeric value
            x = range(t_beg_ref, t_end_ref + 1)
            y = [func_temporal_reg(i) for i in x]

    if not marker_reg:
        marker_reg = 'b--'
    if not label_reg:
        label_reg = "Linear regression"

    if xlabel is None:
        xlabel = 'Reference sample timespan'
    if ylabel is None:
        ylabel = 'Floating sample timespan'
    if title is None:
        title = 'Embryos temporal alignment'
    if cmap is None:
        cmap = plt.cm.gnuplot2

    if bins is None:
        bins = [0.5 + val for val in range(0, 90)]
    # colors = cmap([b/90 for b in bins])
    # hist, bin_edges = np.histogram(angles, bins)
    # ax1.bar(bin_edges[:-1], hist, color=colors, alpha=0.8) # HISTOGRAMME NORMAL
    # Hist=[sum(hist[:i]) for i in range(1,len(hist)+1)]
    # ax2.set_colorbar()
    if c:
        cax = plt.scatter(X, Y, c=c, cmap=cmap, s=s, marker='o')
    else:
        cax = plt.scatter(X, Y, s=s, marker='o')
    if x:
        plt.plot(x, y, marker_reg, label=label_reg)
    # ax2.plot(x,y45,'y--', label="Linear regression (angle max=45)")
    # plt.p
t(x,yinit,'g--', label="Linear regression (initial mapping)")
    plt.axis([t_beg_ref, t_end_ref, t_beg_flo, t_end_flo])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc=loc)
    plt.title(title)

    if c:
        cbar = plt.colorbar(cax, ticks=[1, 10, 20, 30, 40, 50, 60, 70, 80, 87.5], orientation='vertical')
        cbar.set_ticklabels([str(angle) + ' degrees'.decode("utf8") for angle in
                             [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]])  # horizontal colorbar

    plt.show(block=block)


###################################################
### INTRA-SEQUENCE-RELATED REGISTRATION METHODS ###
###################################################

def compute_intra_sequence_two_by_two_registration(fused_files_format, trsf_files_format, begin, end, delta=1,
                                                   verbose=False):
    """
	Step-by-step registration using sequence of fused images

	Inputs:
		fused_files_format: fused images names for blockmatching rigid registration (eg. '/fakepath/171106-Leonard-St8_fuse_t$TIME.inr')
		begin: index for sequence beginning
		end: index for sequence ending
		delta: delta between two indices in the sequence
		verbose: level of verbosity
	Output:
		trsf_files_format: intra registration two-by-two (step-by-step) trsf file names (eg. '/fakepath/171106-Leonard-St8_reg_t$TIMEFLO_t$TIMEREF.trsf')
					ie each transformation is a transformation T_flo<-ref enabling to resample the floating image onto reference image referential
	"""
    from lineage import timeNamed, timesNamed

    trsf_files_path = os.path.sep.join(trsf_files_format.split(os.path.sep)[:-1])
    if not os.path.isdir(trsf_files_path):
        os.mkdir(trsf_files_path)
    assert os.path.isdir(trsf_files_path)

    for t in range(begin, end):
        time_ref = t + delta  # Time point of Segmentation
        time_flo = t
        if verbose:
            print 'Starting the step-by-step registration between ' + str(time_flo) + " and " + str(time_ref)
        fused_file_ref = timeNamed(fused_files_format, time_ref)  # Reference image file
        fused_file_flo = timeNamed(fused_files_format, time_flo)  # Floating image file
        trsf_file = timesNamed(trsf_files_format, '$TIMEFLO', time_flo, '$TIMEREF', time_ref)  # trsf file
        assert os.path.exists(fused_file_ref), "Input image file %s not found." % fused_file_ref
        assert os.path.exists(fused_file_flo), "Input image file %s not found." % fused_file_flo
        rigid_registration(fused_file_ref, fused_file_flo, None, "/dev/null", trsf_file, verbose=verbose)


def compute_intra_sequence_respatialization_trsfs(trsf_files_two_by_two_format, trsf_files_composed_format, begin, end,
                                                  time_point_ref=None, index_format='%03d', verbose=False):
    '''
	Recomputing the trsfs with a unique reference

	Inputs:
		trsf_files_two_by_two_format: intra registration two-by-two (step-by-step) trsf file names (eg. '/fakepath/171106-Leonard-St8_reg_t$TIMEFLO_t$TIMEREF.trsf')
					ie each transformation is a transformation T_flo<-ref enabling to resample the floating image onto reference image referential
		begin: index for sequence beginning
		end: index for sequence ending
		time_point_ref: time-point of the sequence which will be used as reference for output transformations (ie for this time-point, the resulting trsf will be identity)
					Default is time_point_ref=begin
		index_format: format of the indices (default is '%03d')
		verbose: level of verbosity
	Output:
		trsf_files_composed_format: intra registration recomputed trsf file names so that each image can be put in a unique referential.
					Two admitted formats: '/fakepath/171106-Leonard-St8_reg_compose_t$TIME_t$TIME.trsf'
										  '/fakepath/171106-Leonard-St8_reg_compose_t$TIME.trsf'
					Output trsfs enable the resampling of each image onto a unique referential (T_flo<-ref)
	Return value:
		time_point_ref: %d
	'''

    if time_point_ref == None:
        time_point_ref = begin

    input_trsfs_path = os.path.sep.join(trsf_files_two_by_two_format.split(os.path.sep)[:-1])
    assert os.path.isdir(input_trsfs_path)
    trsf_files_path = os.path.sep.join(trsf_files_composed_format.split(os.path.sep)[:-1])
    if not os.path.isdir(trsf_files_path):
        os.mkdir(trsf_files_path)
    assert os.path.isdir(trsf_files_path)

    intrareg_step_format = trsf_files_two_by_two_format.replace('$TIMEREF', index_format).replace('$TIMEFLO',
                                                                                                  index_format)
    intrareg_comp_format = trsf_files_composed_format.split('$TIME')
    assert len(intrareg_comp_format) == 2 or len(
        intrareg_comp_format) == 3, "Naming convention not respected for the 2nd parameter of function compute_intra_sequence_respatialization_trsfs. Refer to function 'help'."
    if len(intrareg_comp_format) == 3:
        intrareg_comp_format = intrareg_comp_format[0] + index_format + intrareg_comp_format[1] + (
                    index_format % time_point_ref) + intrareg_comp_format[2]
    else:
        intrareg_comp_format = trsf_files_composed_format.replace('$TIME', index_format)
    multiple_trsfs(intrareg_step_format, intrareg_comp_format, begin, end, time_point_ref, verbose=verbose)
    return time_point_ref


def compute_optimized_intra_sequence_respatialization_trsfs(image_files_format, trsf_files_composed_format,
                                                            optimized_template, trsf_files_optimized_format, begin, end,
                                                            threshold=2, iso=1.0, margin=10, verbose=False):
    '''
	Usage of function "changeMultipleTrsfs" for
		- computing the minimal window for building a 3D+t subsequence of images of the full sequence
		- building a template image into which the original images will be transformed

	Inputs:
		image_files_format: format for the stack of images on which the user wants to be based on for the optimal box computation (ideally a sequence of segmented images)
				Accepted format example: '/fakepath/171106-Leonard-St8_fuse_seg_post_t$TIME.inr'
		trsf_files_composed_format: intra registration trsf file names so that each image can be put in a unique spatial referential
				Accepted format example: '/fakepath/171106-Leonard-St8_reg_compose_t$TIME_t000.trsf' (/!\ only one occurence of '$TIME' is accepted)
		begin: index for sequence beginning
		end: index for sequence ending
		threshold: threshold on input templates/images to compute the useful bounding box (else it is the entire image)
		iso: make voxels isotropic for the output template (%f, default is 1.0) (if no value is given, uses the smallest voxel size from the template(s))
		margin: add a margin (in voxels)
		verbose: level of verbosity
	Outputs:
		optimized_template: output template image corresponding to the output transformations
		trsf_files_optimized_format: format 'a la printf' for output transformations will allow to resample input image into the resulting
             	template (thus still of the form T_{i<-ref}) reference is changed if '-index-reference' is used
				Accepted format example: '/fakepath/171106-Leonard-St8_reg_compose_t$TIME.trsf'
	'''

    assert image_files_format.count('$TIME') == 1
    assert trsf_files_composed_format.count('$TIME') == 1
    assert trsf_files_optimized_format.count('$TIME') == 1

    input_trsfs_path = os.path.sep.join(trsf_files_composed_format.split(os.path.sep)[:-1])
    assert os.path.isdir(input_trsfs_path)
    input_images_path = os.path.sep.join(image_files_format.split(os.path.sep)[:-1])
    assert os.path.isdir(input_images_path)
    trsf_files_path = os.path.sep.join(trsf_files_optimized_format.split(os.path.sep)[:-1])
    if not os.path.isdir(trsf_files_path):
        os.mkdir(trsf_files_path)
    assert os.path.isdir(trsf_files_path)
    template_path = os.path.sep.join(optimized_template.split(os.path.sep)[:-1])
    if not os.path.isdir(template_path):
        os.mkdir(template_path)
    assert os.path.isdir(template_path)

    intrareg_comp_format = trsf_files_composed_format.replace('$TIME', '%03d')
    intrareg_change_trsf_format = trsf_files_optimized_format.replace('$TIME', '%03d')
    intrareg_change_images_format = image_files_format.replace('$TIME', '%03d')
    change_multiple_trsfs(intrareg_comp_format, intrareg_change_trsf_format, intrareg_change_images_format,
                          optimized_template,
                          begin, end,
                          reference_index=None, trsf_type='rigid',
                          threshold=threshold, iso=iso, margin=margin, verbose=verbose)


def compose_transformation_stack_with_a_transformation(format_input_trsfs, file_trsf, format_output_trsfs, begin, end,
                                                       verbose=False):
    '''
	Function for stack of transformations composition with an outsider transformation.
	Useful for example in order to force a sequence realignment under a given orientation.

	Input data:
		format_input_trsfs: str containing transformations format with index appearing as '$TIME' (eg. '/fakepath/171106-Leonard-St8_reg_compose_t$TIME.trsf')
		file_trsf: str containing path to the file of transformation which will be used for the composition
		begin: index for sequence beginning (MANDATORY)
		end: index for sequence ending (MANDATORY)
		verbose; level of verbosity

	Output data:
		format_outputs: str containing output transformations format with index appearing as '$TIME' (eg. '/fakepath/COMPOSE_GERMINAL/171106-Leonard-St8_reg_compose_germinal_t$TIME.trsf')
	'''
    input_trsfs_path = os.path.sep.join(format_input_trsfs.split(os.path.sep)[:-1])
    assert os.path.isdir(input_trsfs_path)
    assert os.path.exists(file_trsf), "Input transformation file '%s' not found." % file_trsf
    outputs_path = os.path.sep.join(format_output_trsfs.split(os.path.sep)[:-1])
    if not os.path.isdir(outputs_path):
        os.mkdir(outputs_path)
    assert os.path.isdir(outputs_path)

    assert format_input_trsfs.count('$TIME') == 1
    assert format_output_trsfs.count('$TIME') == 1

    for t in range(begin, end + 1):
        path_trsf = timeNamed(format_input_trsfs, t)
        assert os.path.exists(path_trsf), "Input transformation file '%s' not found." % path_trsf
        path_out = timeNamed(format_output_trsfs, t)
        compose_trsf(path_trsf, file_trsf, path_output=path_out, lazy=True, verbose=verbose)


def intra_sequence_realignment(format_images, format_outputs, format_trsfs, template_image=None, begin=None, end=None,
                               delta=1, nearest=True, iso=None, threshold=None, margin=None, visu=False, verbose=False):
    '''
	Function for embryo image sequence

	Input data:
		format_images: str containing input images format with index appearing as '$TIME' (eg. '/fakepath/171106-Leonard-St8_glas_seg_post_t$TIME.inr')
					   (which will be also used for the optimal box computation if 'iso' parameter is specified)
		format_trsfs: str containing rigid transformations format with index appearing as '$TIME' (eg. '/fakepath/171106-Leonard-St8_reg_compose_t$TIME.trsf')
		begin: index for sequence beginning (MANDATORY)
		end: index for sequence ending (MANDATORY)
		delta: delta between two indices in the sequence (default=1)
		nearest : do not interpolate (take the nearest value) if True, to use when applying on label images (default = True). Use nearest = True for segmented images.
		visu:   option enabling labels thin erosion. Useful for segmented images 3D visualisation for example.

		template_image: template image corresponding to the transformations to be applied on the images
			/!\ THIS PARAMETER MUST NOT BE USED WITH 'iso' PARAMETER AT THE SAME TIME
		iso: make voxels isotropic for the output template (%f, default is 1.0) (if no value is given, uses the smallest voxel size from the template(s))
			/!\ THIS PARAMETER MUST NOT BE USED WITH 'template_image' PARAMETER AT THE SAME TIME
		/!\ EITHER 'template_image' OR 'iso' PARAMETER MUST BE PROVIDED (BUT NOT BOTH)
		# Options to be used in complement with 'iso' parameter
		threshold: threshold on input templates/images to compute the useful bounding box (else it is the entire image)
		margin: add a margin (in voxels)
		verbose: level of verbosity

	Output data:
		format_outputs: str containing output images format with index appearing as '$TIME' (eg. '/fakepath/res/171106-Leonard-St8_glas_seg_post_realignment_t$TIME.mha')
	'''
    assert not begin == None
    assert not end == None
    assert format_trsfs, "The trsfs option is not optional"
    assert format_images.count('$TIME') == 1
    assert format_trsfs.count('$TIME') == 1
    assert format_outputs.count('$TIME') == 1

    input_trsfs_path = os.path.sep.join(format_trsfs.split(os.path.sep)[:-1])
    assert os.path.isdir(input_trsfs_path)
    input_images_path = os.path.sep.join(format_images.split(os.path.sep)[:-1])
    assert os.path.isdir(input_images_path)
    outputs_path = os.path.sep.join(format_outputs.split(os.path.sep)[:-1])
    if not os.path.isdir(outputs_path):
        os.mkdir(outputs_path)
    assert os.path.isdir(outputs_path)

    flag_rm = False

    format_optimized_trsfs = ''

    if template_image == None:
        assert not iso == None, "Either 'template_image' or 'iso' parameter must be provided (but not both). Exiting."
        flag_rm = True
        template_image = outputs_path + os.path.sep + "intra_sequence_realignment_template.inr"
        format_optimized_trsfs = outputs_path + os.path.sep + "intra_sequence_realignment_tmp_t$TIME.trsf"
        compute_optimized_intra_sequence_respatialization_trsfs(format_images, format_trsfs, template_image,
                                                                format_optimized_trsfs, begin, end, threshold=threshold,
                                                                iso=iso, margin=margin, verbose=verbose)
    else:
        assert iso == None, "Either 'template_image' or 'iso' parameter must be provided (but not both). Exiting."
        format_optimized_trsfs = format_trsfs

    for t in range(begin, end + 1):
        if verbose:
            print 'Realignment at t=%d ' % t
        image_in = timeNamed(format_images, t)  # Original image file
        image_out = timeNamed(format_outputs, t)  # Output image file
        trsf_file = timeNamed(format_optimized_trsfs, t)  # Transformation to be applied
        assert os.path.exists(image_in), "Input image file '%s' not found." % image_in
        assert os.path.exists(trsf_file), "Input transformation file '%s' not found." % trsf_file
        apply_trsf(image_in, path_trsf=trsf_file, path_output=image_out,
                   template=template_image, nearest=nearest, lazy=True, verbose=verbose)
        if flag_rm:
            cmd = 'rm -f ' + trsf_file
            if verbose:
                print cmd
            os.system(cmd)
        if visu:
            erodeLabels(image_out, image_out, lazy=True, verbose=verbose)

    if flag_rm:
        cmd = 'rm -f ' + template_image
        if verbose:
            print cmd
        os.system(cmd)

