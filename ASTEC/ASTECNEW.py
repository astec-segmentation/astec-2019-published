
import os
import imp
import sys
import time
import shutil
import multiprocessing
import numpy as np
from scipy import ndimage as nd
import copy

import ACE
import MARS
import commonTools
import nomenclature
import lineage
import reconstruction
import CommunFunctions.cpp_wrapping as cpp_wrapping

from CommunFunctions.ImageHandling import imread, imsave, SpatialImage

#
#
#
#
#

monitoring = commonTools.Monitoring()


########################################################################################
#
# classes
# - computation environment
# - computation parameters
#
########################################################################################


class AstecEnvironment(object):

    def __init__(self):

        #
        # fusion data paths
        #
        self.path_fuse_exp = None
        self.path_fuse_exp_files = None

        # segmentation data paths
        #
        self.path_seg_exp = None
        self.path_mars_exp_files = None
        self.path_seg_exp_files = None
        self.path_seg_exp_lineage = None

        #
        #
        #
        self.path_reconstruction = None
        self.temporary_path = None

        #
        #
        #
        self.path_logdir = None
        self.path_history_file = None
        self.path_log_file = None

    def update_from_file(self, parameter_file, start_time):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print ("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = imp.load_source('*', parameter_file)

        self.path_fuse_exp = nomenclature.replaceFlags(nomenclature.path_fuse_exp, parameters)
        self.path_fuse_exp_files = nomenclature.replaceFlags(nomenclature.path_fuse_exp_files, parameters)

        self.path_seg_exp = nomenclature.replaceFlags(nomenclature.path_seg_exp, parameters)
        self.path_mars_exp_files = nomenclature.replaceFlags(nomenclature.path_mars_exp_files, parameters)
        self.path_seg_exp_files = nomenclature.replaceFlags(nomenclature.path_seg_exp_files, parameters)
        self.path_seg_exp_lineage = nomenclature.replaceFlags(nomenclature.path_seg_exp_lineage, parameters)

        self.path_logdir = nomenclature.replaceFlags(nomenclature.path_seg_logdir, parameters)
        self.path_history_file = nomenclature.replaceFlags(nomenclature.path_seg_historyfile, parameters)
        self.path_log_file = nomenclature.replaceFlags(nomenclature.path_seg_logfile, parameters, start_time)

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('AstecEnvironment\n')

            logfile.write('- path_fuse_exp = ' + str(self.path_fuse_exp) + '\n')
            logfile.write('- path_fuse_exp_files = ' + str(self.path_fuse_exp_files) + '\n')

            logfile.write('- path_seg_exp = ' + str(self.path_seg_exp) + '\n')
            logfile.write('- path_mars_exp_files = ' + str(self.path_mars_exp_files) + '\n')
            logfile.write('- path_seg_exp_files = ' + str(self.path_seg_exp_files) + '\n')
            logfile.write('- path_seg_exp_lineage = ' + str(self.path_seg_exp_lineage) + '\n')

            logfile.write('- path_reconstruction = ' + str(self.path_reconstruction) + '\n')
            logfile.write('- temporary_path = ' + str(self.temporary_path) + '\n')

            logfile.write('- path_logdir = ' + str(self.path_logdir) + '\n')
            logfile.write('- path_history_file = ' + str(self.path_history_file)+'\n')
            logfile.write('- path_log_file = ' + str(self.path_log_file)+'\n')
            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('AstecEnvironment')

        print('- path_fuse_exp = ' + str(self.path_fuse_exp))
        print('- path_fuse_exp_files = ' + str(self.path_fuse_exp_files))

        print('- path_seg_exp = ' + str(self.path_seg_exp))
        print('- path_mars_exp_files = ' + str(self.path_mars_exp_files))
        print('- path_seg_exp_files = ' + str(self.path_seg_exp_files))
        print('- path_seg_exp_lineage = ' + str(self.path_seg_exp_lineage))

        print('- path_reconstruction = ' + str(self.path_reconstruction))
        print('- temporary_path = ' + str(self.temporary_path))

        print('- path_logdir = ' + str(self.path_logdir))
        print('- path_history_file = ' + str(self.path_history_file))
        print('- path_log_file = ' + str(self.path_log_file))
        print("")


#
#
#
#
#


class AstecParameters(object):

    def __init__(self):
        #
        #
        #
        self.intensity_transformation = 'Identity'
        self.intensity_enhancement = None
        self.cell_normalization_min_method = 'cellinterior'
        self.cell_normalization_max_method = 'cellborder'

        #
        # membrane enhancement parameters
        #
        self.ace = ACE.AceParameters()

        #
        #
        #
        self.keep_reconstruction = True

        #
        # watershed
        #
        self.watershed_seed_sigma = 0.6
        self.watershed_membrane_sigma = 0.0
        self.watershed_seed_hmin_min_value = 4
        self.watershed_seed_hmin_max_value = 18
        self.watershed_seed_hmin_delta_value = 2

        #
        # to decide whether there will be division
        #
        self.seed_selection_tau = 25

        #
        # threshold
        # cells deformed from previous timepoint that does not have any seed
        # and whose volume (in voxels) is below this threshold are discarded
        # they will correspond to dead-end in the lineage
        #
        self.minimum_volume_unseeded_cell = 100

        #
        # magic values for the volume checking
        # - volume_minimal_value is in volxel units
        #
        self.volume_ratio_tolerance = 0.1
        self.volume_ratio_threshold = 0.5
        self.volume_minimal_value = 1000

        #
        # images suffixes/formats
        #
        self.result_image_suffix = 'inr'
        self.default_image_suffix = 'inr'

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('AstecParameters\n')

            logfile.write('- intensity_transformation = ' + str(self.intensity_transformation) + '\n')
            logfile.write('- intensity_enhancement = ' + str(self.intensity_enhancement) + '\n')
            logfile.write('- cell_normalization_min_method = ' + str(self.cell_normalization_min_method) + '\n')
            logfile.write('- cell_normalization_max_method = ' + str(self.cell_normalization_max_method) + '\n')

            self.ace.write_parameters(log_file_name)

            logfile.write('- keep_reconstruction = ' + str(self.keep_reconstruction) + '\n')

            logfile.write('- watershed_seed_sigma = ' + str(self.watershed_seed_sigma) + '\n')
            logfile.write('- watershed_membrane_sigma = ' + str(self.watershed_membrane_sigma) + '\n')
            logfile.write('- watershed_seed_hmin_min_value = ' + str(self.watershed_seed_hmin_min_value) + '\n')
            logfile.write('- watershed_seed_hmin_max_value = ' + str(self.watershed_seed_hmin_max_value) + '\n')
            logfile.write('- watershed_seed_hmin_delta_value = ' + str(self.watershed_seed_hmin_delta_value) + '\n')

            logfile.write('- seed_selection_tau = ' + str(self.seed_selection_tau) + '\n')

            logfile.write('- minimum_volume_unseeded_cell = ' + str(self.minimum_volume_unseeded_cell) + '\n')

            logfile.write('- volume_ratio_tolerance = ' + str(self.volume_ratio_tolerance) + '\n')
            logfile.write('- volume_ratio_threshold = ' + str(self.volume_ratio_threshold) + '\n')
            logfile.write('- volume_minimal_value = ' + str(self.volume_minimal_value) + '\n')

            logfile.write('- result_image_suffix = ' + str(self.result_image_suffix) + '\n')
            logfile.write('- default_image_suffix = '+str(self.default_image_suffix) + '\n')

            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('AstecParameters')

        print('- intensity_transformation = ' + str(self.intensity_transformation))
        print('- intensity_enhancement = ' + str(self.intensity_enhancement))
        print('- cell_normalization_min_method = ' + str(self.cell_normalization_min_method))
        print('- cell_normalization_max_method = ' + str(self.cell_normalization_max_method))

        self.ace.print_parameters()

        print('- keep_reconstruction = ' + str(self.keep_reconstruction))

        print('- watershed_seed_sigma = ' + str(self.watershed_seed_sigma))
        print('- watershed_membrane_sigma = ' + str(self.watershed_membrane_sigma))
        print('- watershed_seed_hmin_min_value = ' + str(self.watershed_seed_hmin_min_value))
        print('- watershed_seed_hmin_max_value = ' + str(self.watershed_seed_hmin_max_value))
        print('- watershed_seed_hmin_delta_value = ' + str(self.watershed_seed_hmin_delta_value))

        print('- seed_selection_tau = ' + str(self.seed_selection_tau))

        print('- minimum_volume_unseeded_cell = ' + str(self.minimum_volume_unseeded_cell))

        print('- volume_ratio_tolerance = ' + str(self.volume_ratio_tolerance))
        print('- volume_ratio_threshold = ' + str(self.volume_ratio_threshold))
        print('- volume_minimal_value = ' + str(self.volume_minimal_value))

        print('- result_image_suffix = ' + str(self.result_image_suffix))
        print('- default_image_suffix = ' + str(self.default_image_suffix))

        print("")

    def update_from_file(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = imp.load_source('*', parameter_file)

        #
        # reconstruction method
        #

        if hasattr(parameters, 'intensity_transformation'):
            self.intensity_transformation = parameters.intensity_transformation
        if hasattr(parameters, 'astec_intensity_transformation'):
            self.intensity_transformation = parameters.astec_intensity_transformation

        if hasattr(parameters, 'intensity_enhancement'):
            self.intensity_enhancement = parameters.intensity_enhancement
        if hasattr(parameters, 'astec_intensity_enhancement'):
            self.intensity_enhancement = parameters.astec_intensity_enhancement

        if hasattr(parameters, 'cell_normalization_min_method'):
            self.cell_normalization_min_method = parameters.cell_normalization_min_method
        if hasattr(parameters, 'astec_cell_normalization_min_method'):
            self.cell_normalization_min_method = parameters.astec_cell_normalization_min_method

        if hasattr(parameters, 'cell_normalization_max_method'):
            self.cell_normalization_max_method = parameters.cell_normalization_max_method
        if hasattr(parameters, 'astec_cell_normalization_max_method'):
            self.cell_normalization_max_method = parameters.astec_cell_normalization_max_method

        #
        #
        #
        self.ace.update_from_file(parameter_file)

        #
        #
        #
        if hasattr(parameters, 'keep_reconstruction'):
            if parameters.keep_reconstruction is not None:
                self.keep_reconstruction = parameters.keep_reconstruction
        if hasattr(parameters, 'astec_keep_reconstruction'):
            if parameters.astec_keep_reconstruction is not None:
                self.keep_reconstruction = parameters.astec_keep_reconstruction

        #
        # watershed
        #
        if hasattr(parameters, 'astec_watershed_seed_sigma'):
            if parameters.astec_watershed_seed_sigma is not None:
                self.watershed_seed_sigma = parameters.astec_watershed_seed_sigma
        if hasattr(parameters, 'watershed_seed_sigma'):
            if parameters.watershed_seed_sigma is not None:
                self.watershed_seed_sigma = parameters.watershed_seed_sigma
        if hasattr(parameters, 'astec_sigma1'):
            if parameters.astec_sigma1 is not None:
                self.watershed_seed_sigma = parameters.astec_sigma1

        if hasattr(parameters, 'astec_watershed_membrane_sigma'):
            if parameters.astec_watershed_membrane_sigma is not None:
                self.watershed_membrane_sigma = parameters.astec_watershed_membrane_sigma
        if hasattr(parameters, 'watershed_membrane_sigma'):
            if parameters.watershed_membrane_sigma is not None:
                self.watershed_membrane_sigma = parameters.watershed_membrane_sigma
        if hasattr(parameters, 'astec_sigma2'):
            if parameters.astec_sigma2 is not None:
                self.watershed_seed_sigma = parameters.astec_sigma2

        if hasattr(parameters, 'watershed_seed_hmin_min_value'):
            if parameters.watershed_seed_hmin_min_value is not None:
                self.watershed_seed_hmin_min_value = parameters.watershed_seed_hmin_min_value
        if hasattr(parameters, 'astec_h_min_min'):
            if parameters.astec_h_min_min is not None:
                self.watershed_seed_hmin_min_value = parameters.astec_h_min_min

        if hasattr(parameters, 'watershed_seed_hmin_max_value'):
            if parameters.watershed_seed_hmin_max_value is not None:
                self.watershed_seed_hmin_max_value = parameters.watershed_seed_hmin_max_value
        if hasattr(parameters, 'astec_h_min_max'):
            if parameters.astec_h_min_max is not None:
                self.watershed_seed_hmin_max_value = parameters.astec_h_min_max

        if hasattr(parameters, 'watershed_seed_hmin_delta_value'):
            if parameters.watershed_seed_hmin_delta_value is not None:
                self.watershed_seed_hmin_delta_value = parameters.watershed_seed_hmin_delta_value

        #
        # images suffixes/formats
        #
        if hasattr(parameters, 'RESULT_IMAGE_SUFFIX_FUSE'):
            if parameters.result_image_suffix is not None:
                self.result_image_suffix = parameters.RESULT_IMAGE_SUFFIX_FUSE
        if hasattr(parameters, 'result_image_suffix'):
            if parameters.result_image_suffix is not None:
                self.result_image_suffix = parameters.result_image_suffix

        if hasattr(parameters, 'default_image_suffix'):
            if parameters.default_image_suffix is not None:
                self.default_image_suffix = parameters.default_image_suffix
                if not hasattr(parameters, 'result_image_suffix') \
                        and not hasattr(parameters, 'RESULT_IMAGE_SUFFIX_FUSE'):
                    self.result_image_suffix = parameters.default_image_suffix


########################################################################################
#
# some internal procedures
#
########################################################################################

#
# create seeds from previous segmentation
# cells are eroded either with a maximum number of iterations (10 for 'true' cells,
# 25 for the background) or with less iterations if the object to be eroded
# disappears
# Note (GM 15/07/2018): it should be more efficient to use distance maps,
# and to threshold them
#

def _erode_cell(parameters):
    """

    :param parameters:
    :return:
    """
    #
    # Erodes the label i in the label image
    # tmp : binary SpatialImage
    # max_size_cell : size max allow for a cell (here put at np.inf)
    # size_cell : size of the cell to erode
    # iterations : maximum number of iterations for normal cells
    # out_iterations : maximum number of iterations for exterior
    # bb : bounding box if tmp in the global image (necessary when computing in parallel)
    # i : label of the cell to erode
    #

    proc = '_erode_cell'

    tmp, max_size_cell, size_cell, iterations, out_iterations, bb, i = parameters

    #
    # background
    #
    if i == 1:
        nb_iter = out_iterations
    #
    # regular cell
    #
    else:
        nb_iter = iterations

    opened = nd.binary_erosion(tmp, iterations=nb_iter)
    while len(nd.find_objects(opened)) != 1 and nb_iter >= 0:
        nb_iter -= 1
        opened = nd.binary_erosion(tmp, iterations=nb_iter)

    if max_size_cell < size_cell:
        monitoring.to_log_and_console('    ' + proc + ': should not reach this, quite weird', 2)
        num = 1
    else:
        num = i
    return opened, num, i, bb


def _build_seeds_from_previous_segmentation(label_image, output_image, max_size_cell=np.inf, min_size_cell=1000,
                                            cell_iterations=10, background_iterations=25, nprocessors=26):
    """
    Erodes all the labels in the segmented image seg
    :param label_image: image whose cells are to be eroded
    :param output_image:
    :param max_size_cell: size maximum of a cell in number of voxels (infinity ?)
    :param min_size_cell: size minimum of a cell in number of voxels
           cells below this size are not eroded
    :param cell_iterations: maximum number of erosion iteration for cells
    :param background_iterations: maximum number of erosion iteration for background
    :param nprocessors: number maximum of processors allowed to be used
    :return:
    """

    proc = '_build_seeds_from_previous_segmentation'
    seg = imread(label_image)

    bboxes = nd.find_objects(seg)
    a = np.unique(seg)

    # print str(a)
    # print str(bboxes)

    pool = multiprocessing.Pool(processes=nprocessors)
    mapping = []

    seeds = np.zeros_like(seg)

    for i in a:
        tmp = seg[bboxes[i - 1]] == i
        size_cell = np.sum(tmp)
        if size_cell > min_size_cell:
            mapping.append((tmp, max_size_cell, size_cell, cell_iterations, background_iterations, bboxes[i - 1], i))
        else:
            monitoring.to_log_and_console('     .. skip cell ' + str(i) + ', size (' + str(size_cell) + ') <= '
                                          + str(min_size_cell), 2)

    outputs = pool.map(_erode_cell, mapping)
    pool.close()
    pool.terminate()

    for seed, num, i, bb in outputs:
        if num != i:
            monitoring.to_log_and_console('    ' + proc + ': should not reach this, quite weird', 2)
        seeds[bb][seed] = num

    seeds._set_resolution(seg._get_resolution())
    imsave(output_image, seeds)

    return


########################################################################################
#
#
#
########################################################################################

#
# compute the seeds for a range of 'h' values
#

def _extract_seeds_in_cell(parameters):
    """
    Return the seeds in seeds_sub_image stricly included in cell c in cell_segmentation
    """
    #
    # cell_segmentation is a sub-image (extracted from the propagated segmentation at t-1) with 'c' for the cell
    # and 1 for the background
    # seeds_sub_image is a sub-image of the extracted seeds
    # c is the cell label
    #
    cell_segmentation, seeds_sub_image, c = parameters

    #
    # check whether the cell_segmentation has only two labels
    #
    if len(np.unique(cell_segmentation)) != 2:
        monitoring.to_log_and_console('    sub-image of cell ' + str(c) + ' contains '
                                      + str(len(np.unique(cell_segmentation))) + ' labels', 2)
        return

    #
    # get the seeds that intersect the cell 'c'
    #
    labels = list(np.unique(seeds_sub_image[cell_segmentation == c]))

    #
    # remove 0 label (correspond to non-minima regions)
    # Note: check whether 0 is inside labels list (?)
    #
    labels.remove(0)

    nb = len(labels)

    return nb, labels, c


def _cell_based_h_minima(first_segmentation, cells, bounding_boxes, membrane_image, environment, parameters,
                         nprocessors=26):
    """

    :param first_segmentation:
    :param cells:
    :param bounding_boxes:
    :param membrane_image:
    :param environment:
    :param parameters:
    :param nprocessors:
    :return:
    """

    proc = '_cell_based_h_minima'

    MARS.monitoring.copy(monitoring)
    #
    # h-minima extraction with h = max value
    # the difference image is kept for further computation
    #
    h_min = parameters.watershed_seed_hmin_max_value

    input_image = membrane_image
    seed_image = commonTools.add_suffix(membrane_image, "_seed_h" + str('{:03d}'.format(h_min)),
                                        new_dirname=environment.temporary_path,
                                        new_extension=parameters.default_image_suffix)
    difference_image = commonTools.add_suffix(membrane_image, "_seed_diff_h" + str('{:03d}'.format(h_min)),
                                              new_dirname=environment.temporary_path,
                                              new_extension=parameters.default_image_suffix)
    sigma = parameters.watershed_seed_sigma

    if not os.path.isfile(seed_image) or not os.path.isfile(difference_image) \
            or monitoring.forceResultsToBeBuilt is True:
        #
        # computation of regional minima
        # -> keeping the 'difference' image allows to speep up the further computation
        #    for smaller values of h
        #
        MARS.build_seeds(input_image, difference_image, seed_image, sigma, h_min, environment, parameters)
        #
        # select only the 'seeds' that are totally included in cells
        #
        cpp_wrapping.only_keep_seeds_in_cell(seed_image, first_segmentation, seed_image)

    #
    # collect the number of seeds found for each cell
    #
    #

    n_seeds = {}
    parameter_seeds = {}

    checking = True

    while checking:

        #
        # for each cell,
        # 2. build a sub-image (corresponding to the bounding box) from the propagated segmentation from t-1
        #    with the cell labeled at 'c' and the rest at '1'
        # 3. build a sub-image from the seeds extracted at h
        #

        im_segmentation = imread(first_segmentation)
        im_seed = imread(seed_image)

        mapping = []

        for c in cells:
            cell_segmentation = np.ones_like(im_segmentation[bounding_boxes[c]])
            cell_segmentation[im_segmentation[bounding_boxes[c]] == c] = c
            mapping.append((cell_segmentation, im_seed[bounding_boxes[c]], c))

        del im_seed
        del im_segmentation

        pool = multiprocessing.Pool(processes=nprocessors)
        outputs = pool.map(_extract_seeds_in_cell, mapping)
        pool.close()
        pool.terminate()

        #
        # outputs are
        # - nb: the number of labels/seeds that are totally inside cell 'c'
        # - labels: the list of these labels
        # - c: the id of the cell
        #
        returned_n_seeds = []
        for nb, labels, c in outputs:
            returned_n_seeds.append(nb)
            n_seeds.setdefault(c, []).append(nb)
            parameter_seeds.setdefault(c, []).append([h_min, parameters.watershed_seed_sigma])

        #
        # next h value
        #
        h_min -= parameters.watershed_seed_hmin_delta_value

        #
        # still compute while
        # h has not reach the minimum value
        #  and
        # there is at least one cell with a number of seeds in [1, 2]
        # it stops then if all cells have more than 2 seeds
        #
        # Note: I did not catch the utility of 'or returned_n_seeds == []'
        #
        checking = h_min >= parameters.watershed_seed_hmin_min_value \
                   and (((np.array(returned_n_seeds) <= 2) & (np.array(returned_n_seeds) != 0)).any()
                        or returned_n_seeds == [])

        if checking:

            #
            # compute seeds fot this new value of h
            # seeds are computed on the previous 'difference' image
            # therefore smoothing has already be done (to get the first difference image)
            # and is no more required -> sigma = 0.0
            #
            input_image = difference_image
            seed_image = commonTools.add_suffix(membrane_image, "_seed_h" + str('{:03d}'.format(h_min)),
                                                new_dirname=environment.temporary_path,
                                                new_extension=parameters.default_image_suffix)
            difference_image = commonTools.add_suffix(membrane_image, "_seed_diff_h" + str('{:03d}'.format(h_min)),
                                                      new_dirname=environment.temporary_path,
                                                      new_extension=parameters.default_image_suffix)
            sigma = 0.0

            if not os.path.isfile(seed_image) or not os.path.isfile(difference_image) \
                    or monitoring.forceResultsToBeBuilt is True:
                MARS.build_seeds(input_image, difference_image, seed_image, sigma, h_min, environment, parameters,
                                 operation_type='max')
                cpp_wrapping.only_keep_seeds_in_cell(seed_image, first_segmentation, seed_image)

            if not os.path.isfile(seed_image) or not os.path.isfile(difference_image):
                monitoring.to_log_and_console("       " + proc + ": computation failed at h = " + str(h_min), 2)
                monitoring.to_log_and_console("\t Exiting.")
                sys.exit(1)

    return n_seeds, parameter_seeds


########################################################################################
#
#
#
########################################################################################

#
#
#

def _select_seed_parameters(n_seeds, parameter_seeds, tau=25):
    """
    Return the correct h-minima value for each cell
    :param n_seeds: { cell: [#seeds, ] }: dict, key: cell, values: list of #seeds
    :param parameter_seeds: { cell: [[h_min, sigma], ]}: dict matching nb_cells, key: cell, values: list of parameters
    :param tau: magic threshold (see page 72 of L. Guignard PhD thesis)
    :return:
    """

    selected_parameter_seeds = {}
    unseeded_cells = []

    #
    # the selection whether a cell should divide or not is based on the length
    # of the plateau of h values that yield a division (see section 2.3.3.5, pages 70-72
    # of L. Guignard PhD thesis)
    #
    # it can also divided into 2 if there is no h value that gives one seed
    #

    for c, s in n_seeds.iteritems():
        nb_2 = np.sum(np.array(s) == 2)
        nb_3 = np.sum(np.array(s) >= 2)
        score = nb_2*nb_3
        if (s.count(1) or s.count(2)) != 0:
            if score >= tau:
                h, sigma = parameter_seeds[c][np.where(np.array(s) == 2)[0][0]]
                nb_final = 2
            elif s.count(1) != 0:
                h, sigma = parameter_seeds[c][np.where(np.array(s) == 1)[0][0]]
                nb_final = 1
            else:
                h, sigma = parameter_seeds[c][np.where(np.array(s) == 2)[0][0]]
                nb_final = 2
            selected_parameter_seeds[c] = [h, sigma, nb_final]
        elif s.count(3) != 0:
            h, sigma = parameter_seeds[c][s.index(3)]
            selected_parameter_seeds[c] = [h, sigma, 3]
        else:
            unseeded_cells.append(c)
            selected_parameter_seeds[c] = [0, 0, 0]
    return selected_parameter_seeds, unseeded_cells


########################################################################################
#
#
#
########################################################################################

#
# this one is similar to _extract_seeds_in_cell()
#

def _extract_seeds(c, cell_segmentation, cell_seeds=None, bb=None, accept_3_seeds=False):
    """
    Return the seeds from cell_seeds stricly included in cell c from cell_segmentation
    (the labels of the seeds go from 1 to 3)
    :param c: cell label
    :param cell_segmentation: sub-image with 'c' for the cell and '1' for the background
    :param cell_seeds: (sub-)image of labeled h-minima
    :param bb: dilated bounding box of the cell
    :param accept_3_seeds: True if 3 seeds can be accepted as a possible choice
    :return:
    """

    proc = "_extract_seeds"

    #
    # sub-image containing the seeds
    #
    if type(cell_seeds) != SpatialImage:
        seeds_image = imread(cell_seeds)
        if bb is not None:
            seeds = seeds_image[bb]
        else:
            seeds = copy.deepcopy(seeds_image)
        del seeds_image
    else:
        if bb is not None:
            seeds = cell_seeds[bb]
        else:
            seeds = copy.deepcopy(cell_seeds)

    #
    # seeds that intersects the cell
    #
    labels = list(np.unique(seeds[cell_segmentation == c]))
    labels.remove(0)

    #
    # returns
    #
    if len(labels) == 1:
        return 1, (seeds == labels[0]).astype(np.uint8)
    elif len(labels) == 2:
        return 2, ((seeds == labels[0]) + 2 * (seeds == labels[1])).astype(np.uint8)
    elif len(labels) == 3 and not accept_3_seeds:
        #
        # weird, return 3 seeds but label two of them
        #
        monitoring.to_log_and_console("       " + proc + ": weird case, there are 3 seeds but only two are labeled", 2)
        return 3, ((seeds == labels[0]) + 2 * (seeds == labels[1])).astype(np.uint8)
    elif len(labels) == 3 and accept_3_seeds:
        return 3, ((seeds == labels[0]) + 2 * (seeds == labels[1]) + 3 * (seeds == labels[2])).astype(np.uint8)
    else:
        return 0, None


#
#
#


def _build_seeds_from_selected_parameters(selected_parameter_seeds,
                                          segmentation_from_previous, seeds_from_previous, selected_seeds,
                                          cells, unseeded_cells, bounding_boxes, membrane_image,
                                          environment, parameters):
    """

    :param selected_parameter_seeds:
    :param segmentation_from_previous:
    :param seeds_from_previous:
    :param selected_seeds:
    :param cells:
    :param unseeded_cells:
    :param bounding_boxes:
    :param membrane_image:
    :param environment:
    :param parameters:
    :return:
    """

    proc = '_build_seeds_from_selected_parameters'

    first_segmentation = imread(segmentation_from_previous)

    #
    # temporary array of spatial images
    # (to avoid multiple readings of the same image)
    #
    seed_image_list = {}

    #
    #
    #
    new_seed_image = np.zeros_like(first_segmentation, dtype=np.uint16)

    #
    # correspondances: correspondances between cells of previous segmentations and new seed labels
    # divided_cells: siblings
    #
    label_max = 2
    correspondances = {}
    divided_cells = []

    #
    # if one want to keep these informations
    #
    # h_min_information = {}
    # sigma_information = {}

    for c in cells:

        #
        # cells for which no seeds were found for all value of h
        #
        if c in unseeded_cells:
            continue

        #
        # selected_parameter_seeds[c][0] : h_min
        # selected_parameter_seeds[c][1] : sigma
        # # selected_parameter_seeds[c][2] : number of cells (2 or 3 means division)
        #
        h_min = selected_parameter_seeds[c][0]

        #
        # add the seed image to list if required
        #

        if h_min not in seed_image_list:
            seed_image = commonTools.add_suffix(membrane_image, "_seed_h" + str('{:03d}'.format(h_min)),
                                                new_dirname=environment.temporary_path,
                                                new_extension=parameters.default_image_suffix)
            if not os.path.isfile(seed_image):
                monitoring.to_log_and_console("       " + proc + ": '" + str(seed_image).split(os.path.sep)[-1]
                                              + "' was not found", 2)
                monitoring.to_log_and_console("\t Exiting.")
                sys.exit(1)
            seed_image_list[h_min] = imread(seed_image)

        #
        # get the seeds totally included in the cell
        # that was already done in _cell_based_h_minima()
        #
        # cell_segmentation is a sub-image with 'c' for the cell and 1 for the background
        # cell_seeds is a sub-image (same dimensions) of the h-minima image
        #

        cell_segmentation = np.ones_like(first_segmentation[bounding_boxes[c]])
        cell_segmentation[first_segmentation[bounding_boxes[c]] == c] = c
        cell_seeds = seed_image_list[h_min][bounding_boxes[c]]

        #
        # n_seeds: number of seeds totally included in the cell
        # labeled_seeds: sub-image with seeds numbered from 1
        #
        n_seeds, labeled_seeds = _extract_seeds(c, cell_segmentation, cell_seeds)

        #
        # 1 seed
        #
        if n_seeds == 1:
            monitoring.to_log_and_console('      .. process cell ' + str(c) + ' -> ' + str(label_max), 3)
            correspondances[c] = [label_max]
            #
            # if one want to keep h_min and sigma information
            # t designs the previous time, thus t+delta_t is the current time of the image to be segmented
            # h_min_information[(t + delta_t) * 10 ** 4 + label_max] = right_parameters[c][0]
            # sigma_information[(t + delta_t) * 10 ** 4 + label_max] = right_parameters[c][1]
            #
            # here labeled_seeds has only 0 and 1's
            #
            new_seed_image[bounding_boxes[c]][labeled_seeds == 1] = label_max
            label_max += 1
        elif n_seeds == 2 or n_seeds == 3:
            monitoring.to_log_and_console('      .. process cell ' + str(c) + ' -> ' + str(label_max) + ', '
                                          + str(label_max + 1), 3)
            #
            # case n_seeds == 3
            # since _extract_seeds() has been called with 'accept_3_seeds=False'
            # => there are only the two first labeled seeds in 'labeled_seeds'
            #
            correspondances[c] = [label_max, label_max+1]
            divided_cells.append((label_max, label_max+1))
            new_seed_image[bounding_boxes[c]][labeled_seeds == 1] = label_max
            # h_min_information[(t + delta_t) * 10 ** 4 + label_max] = right_parameters[c][0]
            # sigma_information[(t + delta_t) * 10 ** 4 + label_max] = right_parameters[c][1]
            label_max += 1
            new_seed_image[bounding_boxes[c]][labeled_seeds == 2] = label_max
            # h_min_information[(t + delta_t) * 10 ** 4 + label_max] = right_parameters[c][0]
            # sigma_information[(t + delta_t) * 10 ** 4 + label_max] = right_parameters[c][1]
            label_max += 1
        else:
            monitoring.to_log_and_console("       " + proc + ": weird, there were " + str(n_seeds)
                                          + " seeds found for cell " + str(c), 2)

    #
    # create background seed
    # 1. create a background cell
    # 2. get the seeds from the read h-minima image with the smallest h
    # 3. add all the seeds
    #

    background_cell = np.ones_like(first_segmentation)
    background_cell[first_segmentation != 1] = 0

    h_min = min(seed_image_list.keys())
    n_seeds, labeled_seeds = _extract_seeds(1, background_cell, seed_image_list[h_min])
    new_seed_image[labeled_seeds > 0] = 1
    correspondances[1] = [1]

    #
    # create seeds for cell for no seed found
    #
    if len(unseeded_cells) > 0:
        first_seeds = imread(seeds_from_previous)
        for c in unseeded_cells:
            vol = np.sum(first_segmentation == c)
            if vol <= parameters.minimum_volume_unseeded_cell:
                monitoring.to_log_and_console("       " + proc + ": cell " + str(c)
                                              + " from previous time point will have no lineage", 2)
                monitoring.to_log_and_console("       .. volume = " + str(vol), 2)
            else:
                monitoring.to_log_and_console('      .. process cell ' + str(c) + ' -> ' + str(label_max), 3)
                correspondances[c] = [label_max]
                new_seed_image[first_seeds == c] = label_max
                label_max += 1
        del first_seeds
    #
    #
    #

    imsave(selected_seeds, new_seed_image)

    #
    #
    #
    del new_seed_image

    for i in seed_image_list.keys():
        del seed_image_list[i]

    del first_segmentation

    #
    #
    #
    return label_max, correspondances, divided_cells


########################################################################################
#
#
#
########################################################################################


def _compute_volumes(im):
    """

    :param im:
    :return:
    """
    labels = np.unique(im)
    volume = nd.sum(np.ones_like(im), im, index=np.int16(labels))
    return dict(zip(labels, volume))


def _volume_checking(previous_segmentation, segmentation_from_selection,
                     correspondances, n_seeds, parameter_seeds, environment, parameters):
    """

    :param previous_segmentation:
    :param segmentation_from_selection:
    :param correspondances: for each parent cell, give the list of children cells
    :param n_seeds: for each parent cell, give the number of seeds for each couple of parameters [h-min, sigma]
    :param parameter_seeds: for each parent cell, give the list of used parameters [h-min, sigma]
    :param environment:
    :param parameters:
    :return:
    """

    proc = "_volume_checking"

    #
    # compute volumes
    #

    prev_seg = imread(previous_segmentation)
    curr_seg = imread(segmentation_from_selection)

    prev_volumes = _compute_volumes(prev_seg)
    curr_volumes = _compute_volumes(curr_seg)

    #
    # process whole embryo
    #

    prev_embryo_volume = prev_seg.size - prev_volumes[1]
    curr_embryo_volume = curr_seg.size - curr_volumes[1]

    #
    # kept from Leo, very weird formula
    # volume_ratio is the opposite of the fraction of volume lose for the
    # volume from previous time compared to current time
    #
    # volume_ratio < 0 => previous volume > current volume
    # volume_ratio > 0 => previous volume < current volume
    #

    volume_ratio = 1.0 - prev_embryo_volume/curr_embryo_volume

    #
    # default value of parameters.volume_ratio_tolerance is 0.1
    #

    if -parameters.volume_ratio_tolerance <= volume_ratio <= parameters.volume_ratio_tolerance:
        pass
    else:
        if volume_ratio < 0:
            monitoring.to_log_and_console('    .. embryo volume has strongly diminished', 2)
        else:
            monitoring.to_log_and_console('    .. embryo volume has strongly increased', 2)

    #
    # lists of (parent (cell at t), children (cell(s) at t+dt))
    #
    # large_volume_ratio          : volume(mother)   <  SUM volume(childrens)
    # small_volume_ratio          : volume(mother)   >  SUM volume(childrens)
    # abnormal_large_volume_ratio : volume(mother)   << SUM volume(childrens)
    # abnormal_small_volume_ratio : volume(mother)   >> SUM volume(childrens)
    # small_volume_daughter       : volume(children) <  threshold
    #
    large_volume_ratio = []
    small_volume_ratio = []
    abnormal_large_volume_ratio = []
    abnormal_small_volume_ratio = []
    small_volume_daughter = []

    for mother_c, sisters_c in correspondances.iteritems():
        #
        # skip background
        #
        if mother_c == 1:
            continue

        #
        # check whether the volumes exist
        #
        if prev_volumes.has_key(mother_c) is False:
            monitoring.to_log_and_console('    ' + proc + ': no volume for cell ' + str(mother_c)
                                          + ' in previous segmentation', 2)
        for s in sisters_c:
            if curr_volumes.has_key(s) is False:
                monitoring.to_log_and_console('    ' + proc + ': no volume for cell ' + str(s)
                                              + ' in current segmentation', 2)

        print str(mother_c) + " " + str(sisters_c) + " " + str(len(sisters_c))

        #
        # compute ratios
        #
        # volume_ratio = 0  <=> volume(mother) = SUM volume(childrens)
        # volume_ratio < 0  <=> volume(mother) > SUM volume(childrens)
        # volume_ratio > 0  <=> volume(mother) < SUM volume(childrens)
        #
        volume_ratio = 1.0 - prev_volumes[mother_c] / np.sum([curr_volumes.get(s, 1) for s in sisters_c])

        #
        # admissible ratio, check whether the daughter cell(s) are large enough
        # default value of parameters.volume_ratio_tolerance is 0.1
        #
        if -parameters.volume_ratio_tolerance <= volume_ratio <= parameters.volume_ratio_tolerance:
            for s in sisters_c:
                if curr_volumes[s] < parameters.volume_minimal_value:
                    small_volume_daughter.append((mother_c, s))
        else:
        #
        # non-admissible ratios
        # default value of parameters.volume_ratio_threshold is 0.5
        #
            if volume_ratio > 0:
                large_volume_ratio.append((mother_c, sisters_c))
                if volume_ratio > parameters.volume_ratio_threshold:
                    abnormal_large_volume_ratio.append(mother_c)
            elif volume_ratio < 0:
                small_volume_ratio.append((mother_c, sisters_c))
                if volume_ratio < -parameters.volume_ratio_threshold:
                    abnormal_small_volume_ratio.append(mother_c)
            else:
                monitoring.to_log_and_console('    ' + proc + ': should not reach this point', 2)
                monitoring.to_log_and_console('    mother cell was ' + str(mother_c), 2)
                monitoring.to_log_and_console('    daughter cell(s) was(ere) ' + str(sisters_c), 2)

    #
    # here we look at cells that experiment a large decrease of volume
    # ie vol(mother) >> vol(daughter(s))
    # this is the step (1) of section 2.3.3.6 of L. Guignard thesis
    # [corresponds to the list to_look_at in historical astec code]
    #
    for cell in abnormal_small_volume_ratio:
        monitoring.to_log_and_console('      process cell with large decrease of volume', 2)
        #
        # this is similar to _select_seed_parameters()
        # it has already been done ?!
        #
        s = n_seeds[cell]

        #
        # np.sum(np.array(s) == 2) is equivalent to s.count(2)
        #
        nb_2 = np.sum(np.array(s) == 2)
        nb_3 = np.sum(np.array(s) >= 2)
        score = nb_2 * nb_3

        if s.count(1) > 0 or s.count(2) > 0:

            if score >= parameters.seed_selection_tau:
                nb_final = 2
            elif s.count(1) != 0:
                nb_final = 1
            else:
                nb_final = 2

            #
            # w
            #
            if nb_final == 1 and s.count(2) != 0:
                #
                # get the h-min image corresponding to the first case (seeds == 2)
                #
                h_min = parameter_seeds[c][s.index(2)]
                seed_image_name = commonTools.add_suffix(membrane_image, "_seed_h" + str('{:03d}'.format(h_min)),
                                                         new_dirname=environment.temporary_path,
                                                         new_extension=parameters.default_image_suffix)
                if not os.path.isfile(seed_image_name):
                    monitoring.to_log_and_console("       " + proc + ": '" + str(seed_image_name).split(os.path.sep)[-1]
                                                  + "' was not found", 2)
                    monitoring.to_log_and_console("\t Exiting.")
                    sys.exit(1)
                seed_image = imread(seed_image_name)

                pass
                # ... #
            #
            #
            #
        elif s.count(3) != 0:
            pass
            #
            #
            #
        else:
            monitoring.to_log_and_console('    ' + proc + ': should not reach this point', 2)
            monitoring.to_log_and_console('    mother cell was ' + str(c), 2)

    return



########################################################################################
#
#
#
########################################################################################

#
#
#
#
#


def astec_process(previous_time, current_time, lineage_tree_information, experiment, environment, parameters):
    """

    :param previous_time:
    :param current_time:
    :param lineage_tree_information:
    :param experiment:
    :param environment:
    :param parameters:
    :return:
    """

    proc = "astec_process"
    default_width = 3

    acquisition_time = str('{:0{width}d}'.format(current_time, width=default_width))

    #
    # nothing to do if the segmentation image exists
    #

    astec_name = nomenclature.replaceTIME(environment.path_seg_exp_files, current_time)
    astec_image = commonTools.find_file(environment.path_seg_exp, astec_name, monitoring=None, verbose=False)

    if astec_image is not None:
        if monitoring.forceResultsToBeBuilt is False:
            monitoring.to_log_and_console('    astec image already existing', 2)
            return
        else:
            monitoring.to_log_and_console('    astec image already existing, but forced', 2)

    astec_image = os.path.join(environment.path_seg_exp, astec_name + '.' + parameters.result_image_suffix)

    #
    # build or read the membrane image to be segmented
    #

    reconstruction.monitoring.copy(monitoring)

    membrane_image = reconstruction.build_membrane_image(current_time, environment, parameters,
                                                         previous_time=previous_time)
    if membrane_image is None:
        monitoring.to_log_and_console("    .. " + proc + ": no membrane image was found/built for time "
                                      + str(current_time), 2)
        return False

    #
    # build seeds by eroding previous segmentation and deforming it
    #
    # erosion iterations are set by default in voxel units
    # there is also a volume defined in voxel units
    #
    # it may be worth trying to deform first the previous segmentation and then
    # extract the seeds
    #
    monitoring.to_log_and_console('    build seeds from previous segmentation', 2)

    previous_segmentation = reconstruction.get_segmentation_image(previous_time, environment)
    if previous_segmentation is None:
        monitoring.to_log_and_console("    .. " + proc + ": no segmentation image was found for time "
                                      + str(previous_time), 2)
        return False

    undeformed_seeds = commonTools.add_suffix(previous_segmentation, '_undeformed_seeds_from_previous',
                                              new_dirname=environment.temporary_path,
                                              new_extension=parameters.default_image_suffix)
    if not os.path.isfile(undeformed_seeds):
        _build_seeds_from_previous_segmentation(previous_segmentation, undeformed_seeds)

    deformed_seeds = commonTools.add_suffix(previous_segmentation, '_deformed_seeds_from_previous',
                                            new_dirname=environment.temporary_path,
                                            new_extension=parameters.default_image_suffix)

    if not os.path.isfile(deformed_seeds):
        deformation = reconstruction.get_deformation_from_current_to_previous(current_time, environment,
                                                                              parameters, previous_time)
        cpp_wrapping.apply_transformation(undeformed_seeds, deformed_seeds, deformation,
                                          interpolation_mode='nearest', monitoring=monitoring)

    #
    # watershed segmentation with seeds extracted from previous segmentation
    #
    #
    monitoring.to_log_and_console('    watershed from previous segmentation', 2)
    segmentation_from_previous = commonTools.add_suffix(membrane_image, '_watershed_from_previous',
                                                        new_dirname=environment.temporary_path,
                                                        new_extension=parameters.default_image_suffix)
    if not os.path.isfile(segmentation_from_previous):
        MARS.watershed(deformed_seeds, membrane_image, segmentation_from_previous, environment, parameters, sigma=0.0)

    #
    # bounding_boxes: bounding boxes for each cell from the watershed segmentation
    # cells: list of cell labels
    #
    first_segmentation = imread(segmentation_from_previous)
    cells = list(np.unique(first_segmentation))
    cells.remove(1)
    bounding_boxes = dict(zip(range(1, max(cells) + 1), nd.find_objects(first_segmentation)))
    del first_segmentation

    #
    # for each cell and a collection of h values,
    # get the number of seeds
    #
    monitoring.to_log_and_console('    estimation of h-minima for h in ['
                                  + str(parameters.watershed_seed_hmin_min_value) + ','
                                  + str(parameters.watershed_seed_hmin_max_value) + ']', 2)
    n_seeds, parameter_seeds = _cell_based_h_minima(segmentation_from_previous, cells, bounding_boxes, membrane_image,
                                                    environment, parameters)

    #
    # from above results (ie, the number of seeds for every value of h),
    # select the parameter (ie h value)
    # Note: sigma (smoothing parameter to extract seeds) is also kept here, meaning that
    #       it also could be used for selection
    #
    monitoring.to_log_and_console('    parameter selection', 2)
    selected_parameter_seeds, unseeded_cells = _select_seed_parameters(n_seeds, parameter_seeds,
                                                                       tau=parameters.seed_selection_tau)

    #
    # print out the list of cells without seeds and the list of cells that may divide
    #
    if len(unseeded_cells) > 0:
        string = ""
        for i in range(len(unseeded_cells)):
            string = str(unseeded_cells[i])
            if i < len(unseeded_cells)-1:
                string += ', '
        monitoring.to_log_and_console('    .. cells with no seeds: ' + string, 2)
    string = ""
    for c in cells:
        if c in unseeded_cells:
            continue
        if selected_parameter_seeds[c][2] > 1:
            if string == "":
                string = str(c)
            else:
                string += ", " + str(c)
    if string != "":
        monitoring.to_log_and_console('    .. cells that will divide: ' + string, 2)


    #
    # build an image of seeds with selected parameters
    #
    monitoring.to_log_and_console('    build seed image from all h-minima images', 2)

    selected_seeds = commonTools.add_suffix(membrane_image, '_seeds_from_selection',
                                                        new_dirname=environment.temporary_path,
                                                        new_extension=parameters.default_image_suffix)

    output = _build_seeds_from_selected_parameters(selected_parameter_seeds,
                                                    segmentation_from_previous, deformed_seeds,
                                                    selected_seeds,
                                                    cells, unseeded_cells, bounding_boxes, membrane_image,
                                                    environment, parameters)

    label_max, correspondances, divided_cells = output

    #
    # watershed segmentation with the selected seeds
    #

    monitoring.to_log_and_console('    watershed from selection of seeds', 2)
    segmentation_from_selection = commonTools.add_suffix(membrane_image, '_watershed_from_selection',
                                                        new_dirname=environment.temporary_path,
                                                        new_extension=parameters.default_image_suffix)
    if not os.path.isfile(segmentation_from_selection):
        MARS.watershed(selected_seeds, membrane_image, segmentation_from_selection, environment, parameters, sigma=0.0)

    #
    # volume checking
    #
    monitoring.to_log_and_console('    volume checking', 2)
    _volume_checking(previous_segmentation, segmentation_from_selection, correspondances, n_seeds, parameter_seeds, environment, parameters)

    sys.exit(1)

    lineage_tree = lineage_tree_information.get('lin_tree', {})


    tmp = lineage_tree_information.get('volumes_information', {})
    volumes_t_1 = {k % 10 ** 4: v for k, v in tmp.iteritems() if k / 10 ** 4 == t}
    h_min_information = {}

    treated = []



    return


#
# check whether a lineage file exists
# loops over the time points
#


def astec_control(experiment, environment, parameters):
    """

    :param experiment:
    :param environment:
    :param parameters:
    :return:
    """

    proc = "astec_control"
    default_width = 3

    #
    # make sure that the result directory exists
    #

    if not os.path.isdir(environment.path_seg_exp):
        monitoring.to_log_and_console(proc + ": weird, '" + str(environment.path_seg_exp) + "' does not exists", 1)
        monitoring.to_log_and_console("\t Exiting")
        sys.exit(1)

    if not os.path.isdir(environment.path_logdir):
        os.makedirs(environment.path_logdir)

    monitoring.to_log_and_console('', 1)

    #
    # re-read the lineage file, if any
    # and check whether any time point should be re-computed
    #

    first_time_point = experiment.first_time_point + experiment.delay_time_point
    last_time_point = experiment.last_time_point + experiment.delay_time_point

    lineage_tree_information = lineage.read_lineage_tree(environment.path_seg_exp_lineage)
    # print environment.path_seg_exp_lineage
    # print lineage_tree_information

    if len(lineage_tree_information) > 0 and 'lin_tree' in lineage_tree_information:
        monitoring.to_log_and_console("    .. test '" + str(environment.path_seg_exp_lineage) + "'", 1)
        cellat = {}
        for y in lineage_tree_information['lin_tree']:
            t = y/10**4
            if t not in cellat:
                cellat[t] = 1
            else:
                cellat[t] += 1

        restart = -1
        t = first_time_point
        while restart == -1 and t <= last_time_point:
            #
            # possible time point of segmentation, test if ok
            #
            time_value = t + experiment.delta_time_point
            acquisition_time = str('{:0{width}d}'.format(time_value, width=default_width))
            segmentation_file = environment.path_seg_exp_files.replace(nomenclature.FLAG_TIME, acquisition_time)
            if not os.path.isfile(os.path.join(environment.path_seg_exp, segmentation_file)):
                monitoring.to_log_and_console("       image '" + segmentation_file + "' not found", 1)
                restart = t
            else:
                if cellat[t] == 0:
                    monitoring.to_log_and_console("       lineage of image '" + segmentation_file + "' not found", 1)
                    restart = t
                else:
                    try:
                        segmentation_image = imread(segmentation_file)
                    except IOError:
                        monitoring.to_log_and_console("       error in image '" + segmentation_file + "'", 1)
                        restart = t
            #
            #
            #

            if restart == -1:
                monitoring.to_log_and_console("       time '" + str(t) + "' seems ok", 1)
            t += 1
        first_time_point = restart
        monitoring.to_log_and_console("    .. " + proc + ": restart computation at time '"
                                      + str(first_time_point) + "'", 1)
    else:
        monitoring.to_log_and_console("    .. " + proc + ": start computation at time '"
                                      + str(first_time_point) + "'", 1)

    #
    #
    #

    for current_time in range(first_time_point + experiment.delay_time_point + experiment.delta_time_point,
                              last_time_point + experiment.delay_time_point + 1, experiment.delta_time_point):

        acquisition_time = str('{:0{width}d}'.format(current_time, width=default_width))
        previous_time = current_time - experiment.delta_time_point

        #
        # start processing
        #

        monitoring.to_log_and_console('... astec processing of time ' + acquisition_time, 1)
        start_time = time.time()

        #
        # set and make temporary directory
        #

        environment.temporary_path = os.path.join(environment.path_seg_exp, "TEMP_$TIME")
        environment.temporary_path = environment.temporary_path.replace(nomenclature.FLAG_TIME, acquisition_time)

        if not os.path.isdir(environment.temporary_path):
            os.makedirs(environment.temporary_path)

        if parameters.keep_reconstruction is True:
            environment.path_reconstruction = os.path.join(environment.path_seg_exp, "RECONSTRUCTION")
            if not os.path.isdir(environment.path_reconstruction):
                os.makedirs(environment.path_reconstruction)
        else:
            environment.path_reconstruction = environment.temporary_path

        #
        # process
        #

        ret = astec_process(previous_time, current_time, lineage_tree_information, experiment, environment, parameters)
        if ret is False:
            monitoring.to_log_and_console('    an error occurs when processing time ' + acquisition_time, 1)
            return

        #
        # cleaning
        #

        if monitoring.keepTemporaryFiles is False:
            shutil.rmtree(environment.temporary_path)

        #
        # end processing for a time point
        #
        end_time = time.time()

        monitoring.to_log_and_console('    computation time = ' + str(end_time - start_time) + ' s', 1)
        monitoring.to_log_and_console('', 1)

    return
