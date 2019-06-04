
import os
import imp
import sys
import time
import shutil
import numpy as np
from scipy import ndimage as nd

import ACE
import commonTools
import nomenclature
import reconstruction

import CommunFunctions.cpp_wrapping as cpp_wrapping
from CommunFunctions.ImageHandling import SpatialImage, imread, imsave

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


class MarsEnvironment(object):

    def __init__(self):

        #
        # fused data paths
        #
        self.path_fuse_exp = None
        self.path_fuse_exp_files = None

        #
        # mars data paths
        #
        self.path_seg_exp = None
        self.path_mars_exp_files = None

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

        self.path_logdir = nomenclature.replaceFlags(nomenclature.path_seg_logdir, parameters)
        self.path_history_file = nomenclature.replaceFlags(nomenclature.path_seg_historyfile, parameters)
        self.path_log_file = nomenclature.replaceFlags(nomenclature.path_seg_logfile, parameters, start_time)

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('MarsEnvironment\n')

            logfile.write('- path_fuse_exp = ' + str(self.path_fuse_exp) + '\n')
            logfile.write('- path_fuse_exp_files = ' + str(self.path_fuse_exp_files) + '\n')

            logfile.write('- path_seg_exp = ' + str(self.path_seg_exp) + '\n')
            logfile.write('- path_mars_exp_files = ' + str(self.path_mars_exp_files) + '\n')

            logfile.write('- path_reconstruction = ' + str(self.path_reconstruction) + '\n')
            logfile.write('- temporary_path = ' + str(self.temporary_path) + '\n')

            logfile.write('- path_logdir = ' + str(self.path_logdir) + '\n')
            logfile.write('- path_history_file = ' + str(self.path_history_file)+'\n')
            logfile.write('- path_log_file = ' + str(self.path_log_file)+'\n')
            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('MarsEnvironment')

        print('- path_fuse_exp = ' + str(self.path_fuse_exp))
        print('- path_fuse_exp_files = ' + str(self.path_fuse_exp_files))

        print('- path_seg_exp = ' + str(self.path_seg_exp))
        print('- path_mars_exp_files = ' + str(self.path_mars_exp_files))

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


class MarsParameters(object):

    def __init__(self):
        #
        #
        #
        self.first_time_point = -1
        self.last_time_point = -1

        #
        #
        #
        self.intensity_transformation = 'Identity'
        self.intensity_enhancement = None

        #
        # membrane enhancement parameters
        #
        self.ace = ACE.AceParameters()

        #
        #
        #
        self.keep_reconstruction = True

        #
        # watershed parameters
        #
        self.watershed_seed_sigma = 0.6
        self.watershed_membrane_sigma = 0.15
        self.watershed_seed_hmin = 4

        #
        # images suffixes/formats
        #
        self.result_image_suffix = 'inr'
        self.default_image_suffix = 'inr'

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('MarsParameters\n')

            logfile.write('- first_time_point = ' + str(self.first_time_point) + '\n')
            logfile.write('- last_time_point = ' + str(self.last_time_point) + '\n')

            logfile.write('- intensity_transformation = ' + str(self.intensity_transformation) + '\n')
            logfile.write('- intensity_enhancement = ' + str(self.intensity_enhancement) + '\n')

            self.ace.write_parameters(log_file_name)

            logfile.write('- keep_reconstruction = ' + str(self.keep_reconstruction) + '\n')

            logfile.write('- watershed_seed_sigma = ' + str(self.watershed_seed_sigma) + '\n')
            logfile.write('- watershed_membrane_sigma = ' + str(self.watershed_membrane_sigma) + '\n')
            logfile.write('- watershed_seed_hmin = ' + str(self.watershed_seed_hmin) + '\n')

            logfile.write('- result_image_suffix = ' + str(self.result_image_suffix) + '\n')
            logfile.write('- default_image_suffix = ' + str(self.default_image_suffix) + '\n')

            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('MarsParameters')

        print('- first_time_point = ' + str(self.first_time_point))
        print('- last_time_point = ' + str(self.last_time_point))

        print('- intensity_transformation = ' + str(self.intensity_transformation))
        print('- intensity_enhancement = ' + str(self.intensity_enhancement))

        self.ace.print_parameters()

        print('- keep_reconstruction = ' + str(self.keep_reconstruction))

        print('- watershed_seed_sigma = ' + str(self.watershed_seed_sigma))
        print('- watershed_membrane_sigma = ' + str(self.watershed_membrane_sigma))
        print('- watershed_seed_hmin = ' + str(self.watershed_seed_hmin))

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
        #
        #
        if hasattr(parameters, 'mars_begin'):
            self.first_time_point = parameters.mars_begin
        if hasattr(parameters, 'mars_end'):
            self.last_time_point = parameters.mars_end

        #
        # reconstruction methods
        #

        if hasattr(parameters, 'mars_method'):
            if parameters.mars_method == 1:
                self.intensity_transformation = 'Identity'
                self.intensity_enhancement = None
            elif parameters.mars_method == 2:
                self.intensity_transformation = None
                self.intensity_enhancement = 'GACE'

        if hasattr(parameters, 'intensity_transformation'):
            self.intensity_transformation = parameters.intensity_transformation
        if hasattr(parameters, 'mars_intensity_transformation'):
            self.intensity_transformation = parameters.mars_intensity_transformation

        if hasattr(parameters, 'intensity_enhancement'):
            self.intensity_enhancement = parameters.intensity_enhancement
        if hasattr(parameters, 'mars_intensity_enhancement'):
            self.intensity_enhancement = parameters.mars_intensity_enhancement

        #
        #
        #
        self.ace.update_from_file(parameter_file)

        #
        #
        #
        if hasattr(parameters, 'mars_keep_reconstruction'):
            if parameters.mars_keep_reconstruction is not None:
                self.keep_reconstruction = parameters.mars_keep_reconstruction
        if hasattr(parameters, 'keep_reconstruction'):
            if parameters.keep_reconstruction is not None:
                self.keep_reconstruction = parameters.keep_reconstruction

        #
        # watershed parameters
        #
        if hasattr(parameters, 'mars_sigma1'):
            if parameters.mars_sigma1 is not None:
                self.watershed_seed_sigma = parameters.mars_sigma1
        if hasattr(parameters, 'mars_watershed_seed_sigma'):
            if parameters.mars_watershed_seed_sigma is not None:
                self.watershed_seed_sigma = parameters.mars_watershed_seed_sigma
        if hasattr(parameters, 'watershed_seed_sigma'):
            if parameters.watershed_seed_sigma is not None:
                self.watershed_seed_sigma = parameters.watershed_seed_sigma

        if hasattr(parameters, 'mars_sigma2'):
            if parameters.mars_sigma2 is not None:
                self.watershed_membrane_sigma = parameters.mars_sigma2
        if hasattr(parameters, 'mars_watershed_membrane_sigma'):
            if parameters.mars_watershed_membrane_sigma is not None:
                self.watershed_membrane_sigma = parameters.mars_watershed_membrane_sigma
        if hasattr(parameters, 'watershed_membrane_sigma'):
            if parameters.watershed_membrane_sigma is not None:
                self.watershed_membrane_sigma = parameters.watershed_membrane_sigma

        if hasattr(parameters, 'mars_h_min'):
            if parameters.mars_h_min is not None:
                self.watershed_seed_hmin = parameters.mars_h_min
        if hasattr(parameters, 'mars_watershed_seed_hmin'):
            if parameters.mars_watershed_seed_hmin is not None:
                self.watershed_seed_hmin = parameters.mars_watershed_seed_hmin
        if hasattr(parameters, 'watershed_seed_hmin'):
            if parameters.watershed_seed_hmin is not None:
                self.watershed_seed_hmin = parameters.watershed_seed_hmin

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

def build_seeds(input_image, output_difference_image, output_seed_image, sigma, h, environment, parameters,
                operation_type='min'):
    """

    :param input_image:
    :param output_difference_image:
    :param output_seed_image:
    :param sigma:
    :param h:
    :param environment:
    :param parameters:
    :param operation_type:
    :return:
    """

    proc = "build_seeds"
    #
    # output_difference_image is the ouput of 'regional_minima'. It is the difference between
    # the reconstructed image after addition of height h and the original image. This way, valley are transformed
    # into hills. Note that the resulting image has its values in [0, h].
    # Such an image can be used as an input for further regional minima computation with *smaller* h
    #

    if os.path.isfile(output_seed_image) and monitoring.forceResultsToBeBuilt is False:
        return

    monitoring.to_log_and_console("    .. seed extraction '" + str(input_image).split(os.path.sep)[-1]
                                  + "' with h = " + str(h), 2)

    if sigma > 0.0:
        monitoring.to_log_and_console("       smoothing '" + str(input_image).split(os.path.sep)[-1]
                                      + "' with sigma = " + str(sigma), 2)
        seed_preimage = commonTools.add_suffix(input_image, "_smoothing_for_seeds",
                                               new_dirname=environment.temporary_path,
                                               new_extension=parameters.default_image_suffix)
        if not os.path.isfile(seed_preimage) or monitoring.forceResultsToBeBuilt is True:
            cpp_wrapping.linear_smoothing(input_image, seed_preimage, sigma,
                                          real_scale=True, other_options="-o 2", monitoring=monitoring)
    else:
        seed_preimage = input_image

    if not os.path.isfile(seed_preimage):
        monitoring.to_log_and_console("       '" + str(seed_preimage).split(os.path.sep)[-1] + "' does not exist", 2)
        monitoring.to_log_and_console("\t Exiting.")
        sys.exit(1)

    #
    # le resultat de regional minima est la difference entre l'image originale et
    # l'image reconstruite apres soustraction d'une hauteur h
    #

    if output_difference_image is None:
        local_difference_image = commonTools.add_suffix(input_image, "_seed_diff_h" + str('{:03d}'.format(h)),
                                                        new_dirname=environment.temporary_path,
                                                        new_extension=parameters.default_image_suffix)
    else:
        local_difference_image = output_difference_image

    if operation_type.lower() == 'min':
        monitoring.to_log_and_console("       extract regional minima '"
                                      + str(seed_preimage).split(os.path.sep)[-1] + "' with h = " + str(h), 2)

        if not os.path.isfile(local_difference_image) or monitoring.forceResultsToBeBuilt is True:
            cpp_wrapping.regional_minima(seed_preimage, local_difference_image, h=h,
                                         monitoring=monitoring)
    elif operation_type.lower() == 'max':
        monitoring.to_log_and_console("       extract regional maxima '"
                                      + str(seed_preimage).split(os.path.sep)[-1] + "' with h = " + str(h), 2)

        if not os.path.isfile(local_difference_image) or monitoring.forceResultsToBeBuilt is True:
            cpp_wrapping.regional_maxima(seed_preimage, local_difference_image, h=h,
                                         monitoring=monitoring)
    else:
        monitoring.to_log_and_console("       " + proc + ": unknown operation type '" + str(operation_type) + "'", 2)
        monitoring.to_log_and_console("\t Exiting.")
        sys.exit(1)

    if not os.path.isfile(local_difference_image):
        monitoring.to_log_and_console("       '" + str(local_difference_image).split(os.path.sep)[-1]
                                      + "' does not exist", 2)
        monitoring.to_log_and_console("\t Exiting.")
        sys.exit(1)

    #
    # seuillage par hysteresis avec un seuil haut a h
    #
    if operation_type.lower() == 'min':
        monitoring.to_log_and_console("       label regional minima '"
                                      + str(local_difference_image).split(os.path.sep)[-1] + "'", 2)
    elif operation_type.lower() == 'max':
        monitoring.to_log_and_console("       label regional maxima '"
                                      + str(local_difference_image).split(os.path.sep)[-1] + "'", 2)

    if not os.path.isfile(output_seed_image) or monitoring.forceResultsToBeBuilt is True:
        cpp_wrapping.connected_components(local_difference_image, output_seed_image, high_threshold=h,
                                          monitoring=monitoring)

    if output_difference_image is None:
        os.remove(local_difference_image)

    return


def watershed(seed_image, membrane_image, result_image, environment, parameters, sigma=0.0):
    """

    :param seed_image:
    :param membrane_image:
    :param result_image:
    :param environment:
    :param parameters:
    :param sigma:
    :return:
    """

    if sigma > 0.0:
        monitoring.to_log_and_console("    .. smoothing '" + str(membrane_image).split(os.path.sep)[-1] + "'", 2)
        height_image = commonTools.add_suffix(membrane_image, "_smoothing_for_watershed",
                                              new_dirname=environment.temporary_path,
                                              new_extension=parameters.default_image_suffix)
        if not os.path.isfile(height_image) or monitoring.forceResultsToBeBuilt is True:
            cpp_wrapping.linear_smoothing(membrane_image, height_image, parameters.watershed_membrane_sigma,
                                          real_scale=True, other_options="-o 2", monitoring=monitoring)
    else:
        height_image = membrane_image

        #
        # watershed
        #
    if not os.path.isfile(seed_image):
        monitoring.to_log_and_console("       '" + str(seed_image).split(os.path.sep)[-1] + "' does not exist", 2)
        monitoring.to_log_and_console("\t Exiting.")
        sys.exit(1)
    if not os.path.isfile(height_image):
        monitoring.to_log_and_console("       '" + str(height_image).split(os.path.sep)[-1] + "' does not exist", 2)
        monitoring.to_log_and_console("\t Exiting.")
        sys.exit(1)

    monitoring.to_log_and_console("    .. watershed '" + str(height_image).split(os.path.sep)[-1] + "'", 2)

    if not os.path.isfile(result_image) or monitoring.forceResultsToBeBuilt is True:
        cpp_wrapping.watershed(seed_image, height_image, result_image, monitoring=monitoring)

    return


def _mars_watershed(template_image, membrane_image, mars_image, environment, parameters):
    """

    :param template_image:
    :param membrane_image:
    :param mars_image:
    :param environment:
    :param parameters:
    :return:
    """
    #
    # computation of seed image
    # - [smoothing]
    # - regional minima
    # - hysteresis thresholding
    #

    seed_image = commonTools.add_suffix(template_image, "_seed_h"
                                        + str('{:03d}'.format(parameters.watershed_seed_hmin)),
                                        new_dirname=environment.temporary_path,
                                        new_extension=parameters.default_image_suffix)
    if not os.path.isfile(seed_image) or monitoring.forceResultsToBeBuilt is True:
        build_seeds(membrane_image, None, seed_image, parameters.watershed_seed_sigma, parameters.watershed_seed_hmin,
                    environment, parameters)

    #
    #
    #
    watershed(seed_image, membrane_image, mars_image, environment, parameters, parameters.watershed_membrane_sigma)

    return


def _volume_diagnosis(mars_image, ncells=10):
    """

    :param mars_image:
    :param ncells:
    :return:
    """

    proc = "_volume_diagnosis"
    if not os.path.isfile(mars_image):
        monitoring.to_log_and_console("    "+proc+": error, '"+str(mars_image)+"' was not found", 2)
        return

    image = imread(mars_image)
    labels = np.unique(image)
    volumes = nd.sum(np.ones_like(image), image, index=np.int16(labels))
    list_for_sort = list()
    for i in range(len(labels)):
        list_for_sort.append([volumes[i],labels[i]])

    #
    # statistics without the background
    #
    m = np.mean(volumes[1:])
    s = np.std(volumes[1:])

    list_for_sort.sort()

    monitoring.to_log_and_console("    .. diagnosis on cell volumes, smallest cells to be looked at", 1)
    monitoring.to_log_and_console("       mean cell volume = " + str(m) + ", standard deviation = " + str(s), 1)
    for i in range(len(labels)):
        if i <= ncells or list_for_sort[i][0] <= m - 2*s:
            monitoring.to_log_and_console('       cell #'+'{:3d}'.format(list_for_sort[i][1])+' volume='
                                          +'{:10.1f}'.format(list_for_sort[i][0]), 1)

    del image
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

def mars_process(current_time, environment, parameters):
    """

    :param current_time:
    :param environment:
    :param parameters:
    :return:
    """

    #
    # nothing to do if the segmentation image exists
    #

    mars_name = nomenclature.replaceTIME(environment.path_mars_exp_files, current_time)
    mars_image = commonTools.find_file(environment.path_seg_exp, mars_name, monitoring=None, verbose=False)

    if mars_image is not None:
        if monitoring.forceResultsToBeBuilt is False:
            monitoring.to_log_and_console('    mars image already existing', 2)
            #
            # compute diagnosis anyway
            #
            mars_image = os.path.join(environment.path_seg_exp, mars_name + '.' + parameters.result_image_suffix)
            _volume_diagnosis(mars_image)
            return
        else:
            monitoring.to_log_and_console('    mars image already existing, but forced', 2)

    mars_image = os.path.join(environment.path_seg_exp, mars_name + '.' + parameters.result_image_suffix)

    #
    #
    #

    input_name = nomenclature.replaceTIME(environment.path_fuse_exp_files, current_time)
    input_image = commonTools.find_file(environment.path_fuse_exp, input_name, monitoring)
    if input_image is None:
        monitoring.to_log_and_console("    .. image '" + input_name + "' not found in '"
                                      + str(environment.path_fuse_exp) + "'", 2)
        monitoring.to_log_and_console("       skip time " + str(current_time), 2)
        return

    input_image = os.path.join(environment.path_fuse_exp, input_image)

    #
    # build the membrane image to be segmented
    #

    reconstruction.monitoring.copy(monitoring)
    membrane_image = reconstruction.build_membrane_image(current_time, environment, parameters)

    if membrane_image is None or not os.path.isfile(membrane_image):
        monitoring.to_log_and_console("       '" + str(membrane_image).split(os.path.sep)[-1]
                                      + "' does not exist", 2)
        monitoring.to_log_and_console("\t Exiting.")
        sys.exit(1)

    #
    # compute the seeded watershed
    #

    _mars_watershed(input_image, membrane_image, mars_image, environment, parameters)

    #
    #
    #
    _volume_diagnosis(mars_image)

    return


#
#
#
#
#


def mars_control(experiment, environment, parameters):
    """

    :param experiment:
    :param environment:
    :param parameters:
    :return:
    """

    default_width = 3

    #
    # make sure that the result directory exists
    #

    if not os.path.isdir(environment.path_seg_exp):
        os.makedirs(environment.path_seg_exp)

    if not os.path.isdir(environment.path_logdir):
        os.makedirs(environment.path_logdir)

    monitoring.to_log_and_console('', 1)

    #
    #
    #

    if parameters.first_time_point < 0 or parameters.last_time_point < 0:
        monitoring.to_log_and_console("... time interval does not seem to be defined in the parameter file")
        monitoring.to_log_and_console("    set parameters 'begin' and 'end'")
        monitoring.to_log_and_console("\t Exiting")
        sys.exit(1)

    if parameters.first_time_point > parameters.last_time_point:
        monitoring.to_log_and_console("... weird time interval: 'begin' = " + str(parameters.first_time_point)
                                      + ", 'end' = " + str(parameters.last_time_point))

    for time_value in range(parameters.first_time_point + experiment.delay_time_point,
                            parameters.last_time_point + experiment.delay_time_point + 1, experiment.delta_time_point):

        acquisition_time = str('{:0{width}d}'.format(time_value, width=default_width))

        #
        # start processing
        #
        monitoring.to_log_and_console('... mars processing of time ' + acquisition_time, 1)
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
        # processing
        #

        mars_process(time_value, environment, parameters)

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
