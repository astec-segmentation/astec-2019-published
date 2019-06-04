
import os
import imp
import shutil
import sys
import time

import commonTools
import nomenclature
from CommunFunctions.ImageHandling import imread
import CommunFunctions.cpp_wrapping as cpp_wrapping

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


class IntraRegEnvironment(object):

    def __init__(self):

        #
        # fusion data paths
        #
        self.path_fuse_exp = None
        self.path_fuse_exp_files = None

        #
        # segmentation data paths
        #
        self.path_seg_exp = None
        self.path_mars_exp_files = None
        self.path_seg_exp_files = None

        self.path_post_exp = None
        self.path_post_exp_files = None

        #
        # registration data path
        #
        self.path_intrareg_exp = None

        self.path_intrareg_cotrsf = None
        self.path_intrareg_cotrsf_files = None

        self.path_intrareg_trsf = None
        self.path_intrareg_trsf_files = None

        #
        # resampled data
        #
        self.path_intrareg_fuse = None
        self.path_intrareg_fuse_files = None
        self.path_intrareg_seg = None
        self.path_intrareg_seg_files = None
        self.path_intrareg_post = None
        self.path_intrareg_post_files = None

        #
        # movies
        #
        self.path_intrareg_movies = None

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
            print("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = imp.load_source('*', parameter_file)

        self.path_fuse_exp = nomenclature.replaceFlags(nomenclature.path_fuse_exp, parameters)
        self.path_fuse_exp_files = nomenclature.replaceFlags(nomenclature.path_fuse_exp_files, parameters)

        self.path_seg_exp = nomenclature.replaceFlags(nomenclature.path_seg_exp, parameters)
        self.path_mars_exp_files = nomenclature.replaceFlags(nomenclature.path_mars_exp_files, parameters)
        self.path_seg_exp_files = nomenclature.replaceFlags(nomenclature.path_seg_exp_files, parameters)

        self.path_post_exp = nomenclature.replaceFlags(nomenclature.path_post_exp, parameters)
        self.path_post_exp_files = nomenclature.replaceFlags(nomenclature.path_post_exp_files, parameters)

        self.path_intrareg_exp = nomenclature.replaceFlags(nomenclature.path_intrareg_exp, parameters)

        self.path_intrareg_cotrsf = nomenclature.replaceFlags(nomenclature.path_intrareg_cotrsf, parameters)
        self.path_intrareg_cotrsf_files = nomenclature.replaceFlags(nomenclature.path_intrareg_cotrsf_files, parameters)
        self.path_intrareg_trsf = nomenclature.replaceFlags(nomenclature.path_intrareg_trsf, parameters)
        self.path_intrareg_trsf_files = nomenclature.replaceFlags(nomenclature.path_intrareg_trsf_files, parameters)

        self.path_intrareg_fuse = nomenclature.replaceFlags(nomenclature.path_intrareg_fuse, parameters)
        self.path_intrareg_fuse_files = nomenclature.replaceFlags(nomenclature.path_intrareg_fuse_files, parameters)
        self.path_intrareg_seg = nomenclature.replaceFlags(nomenclature.path_intrareg_seg, parameters)
        self.path_intrareg_seg_files = nomenclature.replaceFlags(nomenclature.path_intrareg_seg_files, parameters)
        self.path_intrareg_post = nomenclature.replaceFlags(nomenclature.path_intrareg_post, parameters)
        self.path_intrareg_post_files = nomenclature.replaceFlags(nomenclature.path_intrareg_post_files, parameters)

        self.path_intrareg_movies = nomenclature.replaceFlags(nomenclature.path_intrareg_movies, parameters)

        self.path_logdir = nomenclature.replaceFlags(nomenclature.path_intrareg_logdir, parameters)
        self.path_history_file = nomenclature.replaceFlags(nomenclature.path_intrareg_historyfile, parameters)
        self.path_log_file = nomenclature.replaceFlags(nomenclature.path_intrareg_logfile, parameters, start_time)

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('IntraRegEnvironment\n')

            logfile.write('- path_fuse_exp = ' + str(self.path_fuse_exp) + '\n')
            logfile.write('- path_fuse_exp_files = ' + str(self.path_fuse_exp_files) + '\n')

            logfile.write('- path_seg_exp = ' + str(self.path_seg_exp) + '\n')
            logfile.write('- path_mars_exp_files = ' + str(self.path_mars_exp_files) + '\n')
            logfile.write('- path_seg_exp_files = ' + str(self.path_seg_exp_files) + '\n')

            logfile.write('- path_post_exp = ' + str(self.path_post_exp) + '\n')
            logfile.write('- path_post_exp_files = ' + str(self.path_post_exp_files) + '\n')

            logfile.write('- path_intrareg_exp = ' + str(self.path_intrareg_exp) + '\n')

            logfile.write('- path_intrareg_cotrsf = ' + str(self.path_intrareg_cotrsf) + '\n')
            logfile.write('- path_intrareg_cotrsf_files = ' + str(self.path_intrareg_cotrsf_files) + '\n')
            logfile.write('- path_intrareg_trsf = ' + str(self.path_intrareg_trsf) + '\n')
            logfile.write('- path_intrareg_trsf_files = ' + str(self.path_intrareg_trsf_files) + '\n')

            logfile.write('- path_intrareg_fuse = ' + str(self.path_intrareg_fuse) + '\n')
            logfile.write('- path_intrareg_fuse_files = ' + str(self.path_intrareg_fuse_files) + '\n')
            logfile.write('- path_intrareg_seg = ' + str(self.path_intrareg_seg) + '\n')
            logfile.write('- path_intrareg_seg_files = ' + str(self.path_intrareg_seg_files) + '\n')
            logfile.write('- path_intrareg_post = ' + str(self.path_intrareg_post) + '\n')
            logfile.write('- path_intrareg_post_files = ' + str(self.path_intrareg_post_files) + '\n')

            logfile.write('- path_intrareg_movies = ' + str(self.path_intrareg_movies) + '\n')

            logfile.write('- path_logdir = ' + str(self.path_logdir) + '\n')
            logfile.write('- path_history_file = ' + str(self.path_history_file)+'\n')
            logfile.write('- path_log_file = ' + str(self.path_log_file)+'\n')
            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('IntraRegEnvironment')

        print('- path_fuse_exp = ' + str(self.path_fuse_exp))
        print('- path_fuse_exp_files = ' + str(self.path_fuse_exp_files))

        print('- path_seg_exp = ' + str(self.path_seg_exp))
        print('- path_mars_exp_files = ' + str(self.path_mars_exp_files))
        print('- path_seg_exp_files = ' + str(self.path_seg_exp_files))

        print('- path_post_exp = ' + str(self.path_post_exp))
        print('- path_post_exp_files = ' + str(self.path_post_exp_files))

        print('- path_intrareg_exp = ' + str(self.path_intrareg_exp))

        print('- path_intrareg_cotrsf = ' + str(self.path_intrareg_cotrsf))
        print('- path_intrareg_cotrsf_files = ' + str(self.path_intrareg_cotrsf_files))
        print('- path_intrareg_trsf = ' + str(self.path_intrareg_trsf))
        print('- path_intrareg_trsf_files = ' + str(self.path_intrareg_trsf_files))

        print('- path_intrareg_fuse = ' + str(self.path_intrareg_fuse))
        print('- path_intrareg_fuse_files = ' + str(self.path_intrareg_fuse_files))
        print('- path_intrareg_seg = ' + str(self.path_intrareg_seg))
        print('- path_intrareg_seg_files = ' + str(self.path_intrareg_seg_files))
        print('- path_intrareg_post = ' + str(self.path_intrareg_post))
        print('- path_intrareg_post_files = ' + str(self.path_intrareg_post_files))

        print('- path_intrareg_movies = ' + str(self.path_intrareg_movies))

        print('- path_logdir = ' + str(self.path_logdir))
        print('- path_history_file = ' + str(self.path_history_file))
        print('- path_log_file = ' + str(self.path_log_file))
        print("")


#
#
#
#
#


class IntraRegParameters(object):

    def __init__(self):

        #
        # Co-registration parameters
        #

        self.registration = commonTools.RegistrationParameters()
        self.registration.transformation_type = 'rigid'
        self.registration.prefix = 'intra_registration_'

        #
        # intra-sequence transformation parameters
        #

        self.reference_index = None
        self.template_type = "FUSION"
        self.template_threshold = None
        self.resolution = 0.6
        self.margin = None

        #
        # resampling parameters
        #
        self.sigma_segmentation_images = 1.0

        self.resample_fusion_images = True
        self.resample_segmentation_images = False
        self.resample_post_segmentation_images = False

        self.movie_fusion_images = True
        self.movie_segmentation_images = False
        self.movie_post_segmentation_images = False

        self.xy_movie_fusion_images = []
        self.xz_movie_fusion_images = []
        self.yz_movie_fusion_images = []

        self.xy_movie_segmentation_images = []
        self.xz_movie_segmentation_images = []
        self.yz_movie_segmentation_images = []

        self.xy_movie_post_segmentation_images = []
        self.xz_movie_post_segmentation_images = []
        self.yz_movie_post_segmentation_images = []

        #
        # images suffixes/formats
        #

        self.result_image_suffix = 'inr'
        self.default_image_suffix = 'inr'

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('IntraRegParameters\n')

            #
            # Co-registration parameters
            #

            self.registration.write_parameters(log_file_name)

            #
            # intra-sequence transformation parameters
            #

            logfile.write('- reference_index = ' + str(self.reference_index) + '\n')
            logfile.write('- template_type = ' + str(self.template_type) + '\n')
            logfile.write('- template_threshold = ' + str(self.template_threshold) + '\n')
            logfile.write('- resolution = ' + str(self.resolution) + '\n')
            logfile.write('- margin = ' + str(self.margin) + '\n')

            #
            # resampling parameters
            #

            logfile.write('- sigma_segmentation_images = ' + str(self.sigma_segmentation_images) + '\n')
            logfile.write('- resample_fusion_images = ' + str(self.resample_fusion_images) + '\n')
            logfile.write('- resample_segmentation_images = ' + str(self.resample_segmentation_images) + '\n')
            logfile.write('- resample_post_segmentation_images = ' + str(self.resample_post_segmentation_images) + '\n')

            #
            # movie parameters
            #

            logfile.write('- movie_fusion_images = ' + str(self.movie_fusion_images) + '\n')
            logfile.write('- movie_segmentation_images = ' + str(self.movie_segmentation_images) + '\n')
            logfile.write('- movie_post_segmentation_images = ' + str(self.movie_post_segmentation_images) + '\n')

            logfile.write('- xy_movie_fusion_images = ' + str(self.xy_movie_fusion_images) + '\n')
            logfile.write('- xz_movie_fusion_images = ' + str(self.xz_movie_fusion_images) + '\n')
            logfile.write('- yz_movie_fusion_images = ' + str(self.yz_movie_fusion_images) + '\n')

            logfile.write('- xy_movie_segmentation_images = ' + str(self.xy_movie_segmentation_images) + '\n')
            logfile.write('- xz_movie_segmentation_images = ' + str(self.xz_movie_segmentation_images) + '\n')
            logfile.write('- yz_movie_segmentation_images = ' + str(self.yz_movie_segmentation_images) + '\n')

            logfile.write('- xy_movie_post_segmentation_images = ' + str(self.xy_movie_post_segmentation_images) + '\n')
            logfile.write('- xz_movie_post_segmentation_images = ' + str(self.xz_movie_post_segmentation_images) + '\n')
            logfile.write('- yz_movie_post_segmentation_images = ' + str(self.yz_movie_post_segmentation_images) + '\n')

            #
            # images suffixes/formats
            #

            logfile.write('- result_image_suffix = ' + str(self.result_image_suffix) + '\n')
            logfile.write('- default_image_suffix = '+str(self.default_image_suffix) + '\n')

            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('IntraRegParameters')

        #
        # Co-registration parameters
        #

        self.registration.print_parameters()

        #
        # intra-sequence transformation parameters
        #

        print('- reference_index = ' + str(self.reference_index))
        print('- template_type = ' + str(self.template_type))
        print('- template_threshold = ' + str(self.template_threshold))
        print('- resolution = ' + str(self.resolution))
        print('- margin = ' + str(self.margin))

        #
        # resampling parameters
        #

        print('- sigma_segmentation_images = ' + str(self.sigma_segmentation_images))
        print('- resample_fusion_images = ' + str(self.resample_fusion_images))
        print('- resample_segmentation_images = ' + str(self.resample_segmentation_images))
        print('- resample_post_segmentation_images = ' + str(self.resample_post_segmentation_images))

        #
        # movie parameters
        #

        print('- movie_fusion_images = ' + str(self.movie_fusion_images))
        print('- movie_segmentation_images = ' + str(self.movie_segmentation_images))
        print('- movie_post_segmentation_images = ' + str(self.movie_post_segmentation_images))

        print('- xy_movie_fusion_images = ' + str(self.xy_movie_fusion_images))
        print('- xz_movie_fusion_images = ' + str(self.xz_movie_fusion_images))
        print('- yz_movie_fusion_images = ' + str(self.yz_movie_fusion_images))

        print('- xy_movie_segmentation_images = ' + str(self.xy_movie_segmentation_images))
        print('- xz_movie_segmentation_images = ' + str(self.xz_movie_segmentation_images))
        print('- yz_movie_segmentation_images = ' + str(self.yz_movie_segmentation_images))

        print('- xy_movie_post_segmentation_images = ' + str(self.xy_movie_post_segmentation_images))
        print('- xz_movie_post_segmentation_images = ' + str(self.xz_movie_post_segmentation_images))
        print('- yz_movie_post_segmentation_images = ' + str(self.yz_movie_post_segmentation_images))

        #
        # images suffixes/formats
        #

        print('- result_image_suffix = ' + str(self.result_image_suffix))
        print('- default_image_suffix = ' + str(self.default_image_suffix))

        print("")

    def update_from_file(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        margin_is_updated = False
        parameters = imp.load_source('*', parameter_file)

        #
        # co-registration parameters
        #

        self.registration.update_from_file(parameter_file)

        #
        # intra-sequence transformation parameters
        #

        if hasattr(parameters, 'intra_registration_reference_index'):
            if parameters.intra_registration_reference_index is not None:
                self.reference_index = parameters.intra_registration_reference_index
        if hasattr(parameters, 'intra_registration_template_type'):
            if parameters.intra_registration_template_type is not None:
                self.template_type = parameters.intra_registration_template_type
        if hasattr(parameters, 'intra_registration_template_threshold'):
            if parameters.intra_registration_template_threshold is not None:
                self.template_threshold = parameters.intra_registration_template_threshold
        if hasattr(parameters, 'intra_registration_resolution'):
            if parameters.intra_registration_resolution is not None:
                self.resolution = parameters.intra_registration_resolution
        if hasattr(parameters, 'intra_registration_margin'):
            if parameters.intra_registration_margin is not None:
                self.margin = parameters.intra_registration_margin
                margin_is_updated = True

        #
        # resampling parameters
        #

        if hasattr(parameters, 'intra_registration_sigma_segmentation_images'):
            if parameters.intra_registration_sigma_segmentation_images is not None:
                self.sigma_segmentation_images = parameters.intra_registration_sigma_segmentation_images

        if hasattr(parameters, 'intra_registration_resample_fusion_images'):
            if parameters.intra_registration_resample_fusion_images is not None:
                self.resample_fusion_images = parameters.intra_registration_resample_fusion_images
        if hasattr(parameters, 'intra_registration_resample_segmentation_images'):
            if parameters.intra_registration_resample_segmentation_images is not None:
                self.resample_segmentation_images = parameters.intra_registration_resample_segmentation_images
        if hasattr(parameters, 'intra_registration_resample_post_segmentation_images'):
            if parameters.intra_registration_resample_post_segmentation_images is not None:
                self.resample_post_segmentation_images = parameters.intra_registration_resample_post_segmentation_images

        if hasattr(parameters, 'intra_registration_movie_fusion_images'):
            if parameters.intra_registration_movie_fusion_images is not None:
                self.movie_fusion_images = parameters.intra_registration_movie_fusion_images
        if hasattr(parameters, 'intra_registration_movie_segmentation_images'):
            if parameters.intra_registration_movie_segmentation_images is not None:
                self.movie_segmentation_images = parameters.intra_registration_movie_segmentation_images
        if hasattr(parameters, 'intra_registration_movie_post_segmentation_images'):
            if parameters.intra_registration_movie_post_segmentation_images is not None:
                self.movie_post_segmentation_images = parameters.intra_registration_movie_post_segmentation_images

        if hasattr(parameters, 'intra_registration_xy_movie_fusion_images'):
            if parameters.intra_registration_xy_movie_fusion_images is not None:
                self.xy_movie_fusion_images = parameters.intra_registration_xy_movie_fusion_images
        if hasattr(parameters, 'intra_registration_xz_movie_fusion_images'):
            if parameters.intra_registration_xz_movie_fusion_images is not None:
                self.xz_movie_fusion_images = parameters.intra_registration_xz_movie_fusion_images
        if hasattr(parameters, 'intra_registration_yz_movie_fusion_images'):
            if parameters.intra_registration_yz_movie_fusion_images is not None:
                self.yz_movie_fusion_images = parameters.intra_registration_yz_movie_fusion_images

        if hasattr(parameters, 'intra_registration_xy_movie_segmentation_images'):
            if parameters.intra_registration_xy_movie_segmentation_images is not None:
                self.xy_movie_segmentation_images = parameters.intra_registration_xy_movie_segmentation_images
        if hasattr(parameters, 'intra_registration_xz_movie_segmentation_images'):
            if parameters.intra_registration_xz_movie_segmentation_images is not None:
                self.xz_movie_segmentation_images = parameters.intra_registration_xz_movie_segmentation_images
        if hasattr(parameters, 'intra_registration_yz_movie_segmentation_images'):
            if parameters.intra_registration_yz_movie_segmentation_images is not None:
                self.yz_movie_segmentation_images = parameters.intra_registration_yz_movie_segmentation_images

        if hasattr(parameters, 'intra_registration_xy_movie_post_segmentation_images'):
            if parameters.intra_registration_xy_movie_post_segmentation_images is not None:
                self.xy_movie_post_segmentation_images = parameters.intra_registration_xy_movie_post_segmentation_images
        if hasattr(parameters, 'intra_registration_xz_movie_post_segmentation_images'):
            if parameters.intra_registration_xz_movie_post_segmentation_images is not None:
                self.xz_movie_post_segmentation_images = parameters.intra_registration_xz_movie_post_segmentation_images
        if hasattr(parameters, 'intra_registration_yz_movie_post_segmentation_images'):
            if parameters.intra_registration_yz_movie_post_segmentation_images is not None:
                self.yz_movie_post_segmentation_images = parameters.intra_registration_yz_movie_post_segmentation_images

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

        #
        # default for margin is none (ie no margin)
        # which is convenient for fusion images
        # however, when dealing with segmentation images as template, a little margin will be great
        # thus, define the default as 10 (if no specification from the user)
        #
        if parameters.intra_registration_template_type.lower() == 'segmentation' \
                or parameters.intra_registration_template_type.lower() == 'seg' \
                or parameters.intra_registration_template_type.lower() == 'post-segmentation' \
                or parameters.intra_registration_template_type.lower() == 'post_segmentation' \
                or parameters.intra_registration_template_type.lower() == 'post':
            if margin_is_updated is False:
                parameters.intra_registration_margin = 10


########################################################################################
#
# some internal procedures
#
########################################################################################

def _check_data(experiment, data_path, file_format, suffix=None):
    """
    Check whether all the images (from the first time point to the last one) exist
    :param experiment:
    :param data_path:
    :param file_format:
    :param suffix:
    :return:
    """

    first_time_point = experiment.first_time_point + experiment.delay_time_point
    last_time_point = experiment.last_time_point + experiment.delay_time_point

    for current_time in range(first_time_point + experiment.delay_time_point + experiment.delta_time_point,
                              last_time_point + experiment.delay_time_point + 1, experiment.delta_time_point):

        input_name = nomenclature.replaceTIME(file_format, current_time)

        if suffix is None:
            input_image = commonTools.find_file(data_path, input_name, monitoring)

            if input_image is None:
                monitoring.to_log_and_console("    .. image '" + input_name + "' not found in '"
                                              + str(data_path) + "'", 2)
                return False

        else:
            input_image = input_name + '.' + suffix
            if not os.path.isfile(os.path.join(data_path, input_image)):
                monitoring.to_log_and_console("    .. image '" + input_image + "' not found in '"
                                              + str(data_path) + "'", 2)
                return False

    return True


def _get_file_suffix(experiment, data_path, file_format):
    """

    :param experiment:
    :param data_path:
    :param file_format:
    :return:
    """

    proc = "_get_file_suffix"

    first_time_point = experiment.first_time_point + experiment.delay_time_point
    last_time_point = experiment.last_time_point + experiment.delay_time_point

    suffixes = {}
    nimages = 0
    nfiles = 0

    #
    # get and count suffixes for images
    #
    for current_time in range(first_time_point + experiment.delay_time_point + experiment.delta_time_point,
                              last_time_point + experiment.delay_time_point + 1, experiment.delta_time_point):

        file_prefix = nomenclature.replaceTIME(file_format, current_time)

        for f in os.listdir(data_path):
            if len(f) <= len(file_prefix):
                pass
            if f[0:len(file_prefix)] == file_prefix and f[len(file_prefix)] == '.':
                suffix = f[len(file_prefix) + 1:len(f)]
                suffixes[suffix] = suffixes.get(suffix, 0) + 1
                nfiles += 1

        nimages += 1

    for s, n in suffixes.items():
        if n == nimages:
            return s

    if nfiles < nimages:
        monitoring.to_log_and_console(proc + ": weird, not enough images '" + str(file_format)
                                      + "' were found in '" + str(data_path) + "'", 0)
        monitoring.to_log_and_console("\t Exiting.", 0)
        exit(1)

    monitoring.to_log_and_console(proc + ": no common suffix for '" + str(file_format)
                                  + "' was found in '" + str(data_path) + "'", 2)
    monitoring.to_log_and_console("\t time point range was ["+str(first_time_point)+", "+str(last_time_point)+"]")
    return None


########################################################################################
#
#
#
########################################################################################


def _coregistration_control(experiment, environment, parameters):
    """
    Perform the co-registration of any couple of two successive images
    Resulting transformations are computed with fused image at t as floating image
    and used image at t+delta_t as reference image

    :param experiment:
    :param environment:
    :param parameters:
    :return:
    """

    proc = "_coregistration_control"

#    if _check_data(experiment, environment.path_fuse_exp, environment.path_fuse_exp_files) is False:
#        monitoring.to_log_and_console(proc + ": error, some fused data are missing", 1)
#        monitoring.to_log_and_console("\t Exiting")
#        sys.exit(1)

    if not os.path.isdir(environment.path_intrareg_cotrsf):
        os.makedirs(environment.path_intrareg_cotrsf)

    first_time_point = experiment.first_time_point + experiment.delay_time_point
    last_time_point = experiment.last_time_point + experiment.delay_time_point

    for reference_time in range(first_time_point + experiment.delay_time_point + experiment.delta_time_point,
                                last_time_point + experiment.delay_time_point + 1, experiment.delta_time_point):

        floating_time = reference_time - experiment.delta_time_point

        trsf_name = nomenclature.replaceTIME(environment.path_intrareg_cotrsf_files, floating_time,
                                             nomenclature.FLAG_TIMEFLO)
        trsf_name = os.path.join(environment.path_intrareg_cotrsf, nomenclature.replaceTIME(trsf_name, reference_time,
                                                                                            nomenclature.FLAG_TIMEREF))
        if not os.path.isfile(trsf_name) or monitoring.forceResultsToBeBuilt is True:

            floating_prefix = nomenclature.replaceTIME(environment.path_fuse_exp_files, floating_time)
            floating_name = commonTools.find_file(environment.path_fuse_exp, floating_prefix, monitoring)
            if floating_name is None:
                monitoring.to_log_and_console(proc + ": error, image '" + str(floating_prefix) + "' was not found", 1)
                monitoring.to_log_and_console("\t Exiting")
                sys.exit(1)

            reference_prefix = nomenclature.replaceTIME(environment.path_fuse_exp_files, reference_time)
            reference_name = commonTools.find_file(environment.path_fuse_exp, reference_prefix, monitoring)
            if reference_name is None:
                monitoring.to_log_and_console(proc + ": error, image '" + str(reference_prefix) + "' was not found", 1)
                monitoring.to_log_and_console("\t Exiting")
                sys.exit(1)

            floating_image = os.path.join(environment.path_fuse_exp, floating_name)
            reference_image = os.path.join(environment.path_fuse_exp, reference_name)
            monitoring.to_log_and_console("       co-registering '" + floating_name + "'", 2)
            cpp_wrapping.linear_registration(reference_image, floating_image, None, trsf_name, None,
                                             py_hl=parameters.pyramid_highest_level,
                                             py_ll=parameters.pyramid_lowest_level,
                                             transformation_type=parameters.transformation_type,
                                             transformation_estimator=parameters.transformation_estimation_type,
                                             lts_fraction=parameters.lts_fraction,
                                             normalization=parameters.normalization,
                                             monitoring=monitoring)

    return


def _transformations_from_reference(experiment, environment, parameters, temporary_dir):
    """
    Combine the transformations issued from the co-registration of pairs of successive images
    to get transformations from one given image
    :param experiment:
    :param environment:
    :param parameters:
    :return:
    """

    if not os.path.isdir(temporary_dir):
        os.makedirs(temporary_dir)

    first_time_point = experiment.first_time_point + experiment.delay_time_point
    last_time_point = experiment.last_time_point + experiment.delay_time_point

    if parameters.reference_index is None:
        reference_index = first_time_point
    else:
        reference_index = parameters.reference_index

    name_input = environment.path_intrareg_cotrsf_files.replace(nomenclature.FLAG_TIMEFLO, '%03d')
    name_input = name_input.replace(nomenclature.FLAG_TIMEREF, '%03d')
    format_input = os.path.join(environment.path_intrareg_cotrsf, name_input)

    name_output = environment.path_intrareg_trsf_files.replace(nomenclature.FLAG_TIME, '%03d')
    format_output = os.path.join(temporary_dir, name_output)

    #
    #
    #

    cpp_wrapping.multiple_trsfs(format_input, format_output, first_time_point, last_time_point,
                                reference_index, trsf_type=parameters.registration.transformation_type,
                                monitoring=monitoring)

    return


def _transformations_and_template(experiment, environment, parameters, temporary_dir, result_dir, result_template):
    """
    From transformations from one given image, compute the template to resample all images.
    :param experiment:
    :param environment:
    :param parameters:
    :param temporary_dir:
    :param result_dir:
    :param result_template:
    :return:
    """

    proc = "_transformations_and_template"

    if os.path.isfile(result_template) and monitoring.forceResultsToBeBuilt is not True:
        return

    monitoring.to_log_and_console("       warning: this stage may be long", 2)

    if not os.path.isdir(result_dir):
        os.makedirs(result_dir)

    #
    # which format to use?
    # if requested, try to use segmentation image (threshold should be set to 2)
    # else use fusion images
    #
    template_format = None

    if parameters.template_type.lower() == 'segmentation' or parameters.template_type.lower() == 'seg':
        suffix = _get_file_suffix(experiment, environment.path_seg_exp, environment.path_seg_exp_files)

        if suffix is None:
            monitoring.to_log_and_console(proc + ": no consistent naming was found in '"
                                          + str(environment.path_seg_exp) + "'", 1)
            monitoring.to_log_and_console("\t switch to fused images as templates", 1)

        else:
            monitoring.to_log_and_console("       ... build template from segmentation images of '"
                                          + str(environment.path_seg_exp) + "'", 2)
            template_name = environment.path_seg_exp_files.replace(nomenclature.FLAG_TIME, '%03d')
            template_name += '.' + suffix
            template_format = os.path.join(environment.path_seg_exp, template_name)

    elif parameters.template_type.lower() == 'post-segmentation' \
            or parameters.template_type.lower() == 'post_segmentation' or parameters.template_type.lower() == 'post':
        suffix = _get_file_suffix(experiment, environment.path_post_exp, environment.path_post_exp_files)

        if suffix is None:
            monitoring.to_log_and_console(proc + ": no consistent naming was found in '"
                                          + str(environment.path_post_exp) + "'", 1)
            monitoring.to_log_and_console("\t switch to fused images as templates", 1)

        else:
            monitoring.to_log_and_console("       ... build template from post-segmentation images of '"
                                          + str(environment.path_post_exp) + "'", 2)
            template_name = environment.path_post_exp_files.replace(nomenclature.FLAG_TIME, '%03d')
            template_name += '.' + suffix
            template_format = os.path.join(environment.path_post_exp, template_name)

    if template_format is None:
        suffix = _get_file_suffix(experiment, environment.path_fuse_exp, environment.path_fuse_exp_files)
        if suffix is None:
            monitoring.to_log_and_console(proc + ": no consistent naming was found in '"
                                          + str(environment.path_fuse_exp) + "'", 1)
        else:
            monitoring.to_log_and_console("       ... build template from fusion images of '"
                                          + str(environment.path_fuse_exp) + "'", 2)
            template_name = environment.path_fuse_exp_files.replace(nomenclature.FLAG_TIME, '%03d')
            template_name += '.' + suffix
            template_format = os.path.join(environment.path_fuse_exp, template_name)

    #
    # other parameters
    #

    first_time_point = experiment.first_time_point + experiment.delay_time_point
    last_time_point = experiment.last_time_point + experiment.delay_time_point

    name_trsf = environment.path_intrareg_trsf_files.replace(nomenclature.FLAG_TIME, '%03d')
    format_input = os.path.join(temporary_dir, name_trsf)
    format_output = os.path.join(result_dir, name_trsf)

    if parameters.reference_index is None:
        reference_index = first_time_point
    else:
        reference_index = parameters.reference_index

    #
    #
    #

    cpp_wrapping.change_multiple_trsfs(format_input, format_output, first_time_point, last_time_point, reference_index,
                                       result_template, trsf_type=parameters.registration.transformation_type,
                                       resolution=parameters.resolution, threshold=parameters.template_threshold,
                                       margin=parameters.margin, format_template=template_format, monitoring=monitoring)

    return


def _resample_images(experiment, environment, parameters, dir_input, format_input, dir_output, format_output,
                     trsf_dir, template_image, interpolation_mode='linear'):
    """
    resample all images given a set of transformations and a template
    :param experiment:
    :param environment:
    :param parameters:
    :param dir_input:
    :param format_input:
    :param dir_output:
    :param format_output:
    :param template_image:
    :param interpolation_mode:
    :return:
    """

    proc = "_resample_images"

    if not os.path.isdir(dir_output):
        os.makedirs(dir_output)

    #
    # in case the template has been gziped, or copied into an other format
    #

    b = os.path.basename(template_image)
    d = os.path.dirname(template_image)
    local_template_name = commonTools.find_file(d, b, monitoring)
    if local_template_name is None:
        monitoring.to_log_and_console(proc + ": template '" + str(b) + "' was not found in '" + str(d) + "'", 1)
        monitoring.to_log_and_console("\t resampling will not be done")
        return
    local_template_image = os.path.join(d, local_template_name)

    #
    #
    #

    first_time_point = experiment.first_time_point + experiment.delay_time_point
    last_time_point = experiment.last_time_point + experiment.delay_time_point

    for t in range(first_time_point + experiment.delay_time_point, last_time_point + experiment.delay_time_point + 1,
                   experiment.delta_time_point):

        output_name = nomenclature.replaceTIME(format_output, t)
        output_image = commonTools.find_file(dir_output, output_name, monitoring=None, verbose=False)

        if output_image is None or monitoring.forceResultsToBeBuilt is True:
            input_name = nomenclature.replaceTIME(format_input, t)
            input_image = commonTools.find_file(dir_input, input_name, monitoring)
            output_image = os.path.join(dir_output, output_name + '.' + str(parameters.result_image_suffix))

            if input_image is None:
                monitoring.to_log_and_console(proc + ": image '" + str(input_name) + "' was not found in '"
                                              + str(dir_input) + "'", 1)
            else:
                monitoring.to_log_and_console("       resampling '" + str(input_image) + "'", 2)
                trsf_name = os.path.join(trsf_dir, nomenclature.replaceTIME(environment.path_intrareg_trsf_files, t))
                input_image = os.path.join(dir_input, input_image)
                cpp_wrapping.apply_transformation(input_image, output_image, trsf_name, local_template_image,
                                                  interpolation_mode=interpolation_mode,
                                                  cell_based_sigma=parameters.sigma_segmentation_images,
                                                  monitoring=monitoring)

    return


def _make_movies(experiment, parameters, dir_input, format_input, dir_output, xylist, xzlist, yzlist):
    """

    :param experiment:
    :param parameters:
    :param dir_input:
    :param format_input:
    :param dir_output:
    :param xylist:
    :param xzlist:
    :param yzlist:
    :return:
    """

    if not os.path.isdir(dir_output):
        os.makedirs(dir_output)

    firstindex = experiment.first_time_point + experiment.delay_time_point
    lastindex = experiment.last_time_point + experiment.delay_time_point

    suffix = _get_file_suffix(experiment, dir_input, format_input)

    format_resample = format_input.replace(nomenclature.FLAG_TIME, '%03d')
    format_resample = os.path.join(dir_input, format_resample + '.' + str(suffix))

    xy = []
    xz = []
    yz = []

    #
    # should we switch to default behavior?
    #

    if len(xylist) == 0 and len(xzlist) == 0 and len(yzlist) == 0:
        #
        # read the first image and set XY slice in the middle
        #
        first_prefix = nomenclature.replaceTIME(format_input, firstindex)
        first_name = commonTools.find_file(dir_input, first_prefix, monitoring=None, verbose=False)
        first_image = imread(os.path.join(dir_input, first_name))
        xy.append(int(first_image.shape[2]/2))
        del first_image
    else:
        xy = xylist
        xz = xzlist
        yz = yzlist

    #
    # processing
    #

    if len(xy) > 0:
        for s in xy:
            name_output = format_input.replace(nomenclature.FLAG_TIME, str(firstindex) + '-' + str(lastindex))
            name_output += '_xy' + str(s)
            name_output = os.path.join(dir_output, name_output + '.' + str(parameters.result_image_suffix))
            if os.path.isfile(name_output) is False or monitoring.forceResultsToBeBuilt is True:
                monitoring.to_log_and_console("       process xy=" + str(s), 2)
                cpp_wrapping.crop_sequence(format_resample, name_output, firstindex, lastindex, 'xy', s,
                                           monitoring=monitoring)

    if len(xz) > 0:
        for s in xz:
            name_output = format_input.replace(nomenclature.FLAG_TIME, str(firstindex) + '-' + str(lastindex))
            name_output += 'xz' + str(s)
            name_output = os.path.join(dir_output, name_output + '.' + str(parameters.result_image_suffix))
            if os.path.isfile(name_output) is False or monitoring.forceResultsToBeBuilt is True:
                monitoring.to_log_and_console("       process xz=" + str(s), 2)
                cpp_wrapping.crop_sequence(format_resample, name_output, firstindex, lastindex, 'xz', s,
                                           monitoring=monitoring)

    if len(yz) > 0:
        for s in yz:
            name_output = format_input.replace(nomenclature.FLAG_TIME, str(firstindex) + '-' + str(lastindex))
            name_output += '_yz' + str(s)
            name_output = os.path.join(dir_output, name_output + '.' + str(parameters.result_image_suffix))
            if os.path.isfile(name_output) is False or monitoring.forceResultsToBeBuilt is True:
                monitoring.to_log_and_console("       process yz=" + str(s), 2)
                cpp_wrapping.crop_sequence(format_resample, name_output, firstindex, lastindex, 'yz', s,
                                           monitoring=monitoring)

    return


########################################################################################
#
#
#
########################################################################################


def intraregistration_control(experiment, environment, parameters):
    """

    :param experiment:
    :param environment:
    :param parameters:
    :return:
    """

    proc = "intraregistration_control"

#
    #
    # start processing
    #
    start_time = time.time()

    #
    # if template does not exists, compute transformations
    #
    firstindex = experiment.first_time_point + experiment.delay_time_point
    lastindex = experiment.last_time_point + experiment.delay_time_point

    result_dir = environment.path_intrareg_trsf + "_t" + str(firstindex) + "-" + str(lastindex)
    result_template = "template" + "_t" + str(firstindex) + "-" + str(lastindex) + "." + parameters.result_image_suffix
    result_template = os.path.join(result_dir, result_template)

    if not os.path.isdir(result_dir) or os.path.isfile(result_template) or monitoring.forceResultsToBeBuilt is True:

        if experiment.delta_time_point > 1:
            monitoring.to_log_and_console(proc + ": warning, delta_time=" + str(experiment.delta_time_point)
                                          + ", this step may be fragile", 1)
        #
        # co-registration of any 2 successive images
        #

        monitoring.to_log_and_console("    .. co-registrations", 2)
        _coregistration_control(experiment, environment, parameters.registration)

        #
        # composition of transformations by propagation
        # results in a set of transformation from the reference image (given by the reference_index)
        # towards each image
        #

        monitoring.to_log_and_console("    .. transformation composition", 2)
        temporary_dir = os.path.join(environment.path_intrareg_cotrsf, 'TEMP')
        if not os.path.isdir(temporary_dir):
            os.makedirs(temporary_dir)

        _transformations_from_reference(experiment, environment, parameters, temporary_dir)

        _get_file_suffix(experiment, environment.path_fuse_exp, environment.path_fuse_exp_files)

        #
        # re-composition of transformations
        # template creation for further resampling operations
        #

        monitoring.to_log_and_console("    .. transformation recomposition and template generation", 2)

        _transformations_and_template(experiment, environment, parameters, temporary_dir, result_dir, result_template)

        if monitoring.keepTemporaryFiles is False:
            shutil.rmtree(temporary_dir)

    #
    # template image and resampling transformations have been computed
    # resample images if required
    #
    # NOTE: il y a un probleme de coherence, puisque la premiere image de segmentation peut etre issue
    # de mars et donc etre nommee differemment. Pour le reechantillonage, on pourrait utiliser
    # reconstruction.get_segmentation_image().
    #
    if parameters.resample_fusion_images is True or parameters.movie_fusion_images is True \
            or len(parameters.xy_movie_fusion_images) > 0 or len(parameters.xz_movie_fusion_images) > 0 \
            or len(parameters.yz_movie_fusion_images) > 0:
        monitoring.to_log_and_console("    .. resampling fusion images", 2)
        _resample_images(experiment, environment, parameters, environment.path_fuse_exp,
                         environment.path_fuse_exp_files, environment.path_intrareg_fuse,
                         environment.path_intrareg_fuse_files, result_dir, result_template)

    if parameters.resample_segmentation_images is True or parameters.movie_segmentation_images is True \
            or len(parameters.xy_movie_segmentation_images) > 0 or len(parameters.xz_movie_segmentation_images) > 0 \
            or len(parameters.yz_movie_segmentation_images) > 0:
        monitoring.to_log_and_console("    .. resampling segmentation images", 2)
        _resample_images(experiment, environment, parameters, environment.path_seg_exp, environment.path_seg_exp_files,
                         environment.path_intrareg_seg, environment.path_intrareg_seg_files, result_dir,
                         result_template, interpolation_mode='nearest')

    if parameters.resample_post_segmentation_images is True or parameters.movie_post_segmentation_images is True \
            or len(parameters.xy_movie_post_segmentation_images) > 0 \
            or len(parameters.xz_movie_post_segmentation_images) > 0 \
            or len(parameters.yz_movie_post_segmentation_images) > 0:
        monitoring.to_log_and_console("    .. resampling post-segmentation images", 2)
        _resample_images(experiment, environment, parameters, environment.path_post_exp,
                         environment.path_post_exp_files,
                         environment.path_intrareg_post, environment.path_intrareg_post_files, result_dir,
                         result_template, interpolation_mode='nearest')

    #
    # make 3D=2D+t images = movie of evolving slices with respect to time
    #
    if parameters.movie_fusion_images is True or len(parameters.xy_movie_fusion_images) > 0 \
            or len(parameters.xz_movie_fusion_images) > 0 or len(parameters.yz_movie_fusion_images) > 0:
        monitoring.to_log_and_console("    .. movies from fusion images", 2)
        _make_movies(experiment, parameters, environment.path_intrareg_fuse, environment.path_intrareg_fuse_files,
                     environment.path_intrareg_movies, parameters.xy_movie_fusion_images,
                     parameters.xz_movie_fusion_images, parameters.yz_movie_fusion_images)

    if parameters.movie_segmentation_images is True or len(parameters.xy_movie_segmentation_images) > 0 \
            or len(parameters.xz_movie_segmentation_images) > 0 or len(parameters.yz_movie_segmentation_images) > 0:
        monitoring.to_log_and_console("    .. movies from segmentation images", 2)
        _make_movies(experiment, parameters, environment.path_intrareg_seg, environment.path_intrareg_seg_files,
                     environment.path_intrareg_movies, parameters.xy_movie_segmentation_images,
                     parameters.xz_movie_segmentation_images, parameters.yz_movie_segmentation_images)

    if parameters.movie_post_segmentation_images is True or len(parameters.xy_movie_post_segmentation_images) > 0 \
            or len(parameters.xz_movie_post_segmentation_images) > 0 \
            or len(parameters.yz_movie_post_segmentation_images) > 0:
        monitoring.to_log_and_console("    .. movies from post-segmentation images", 2)
        _make_movies(experiment, parameters, environment.path_intrareg_post, environment.path_intrareg_post_files,
                     environment.path_intrareg_movies, parameters.xy_movie_post_segmentation_images,
                     parameters.xz_movie_post_segmentation_images, parameters.yz_movie_post_segmentation_images)

    #
    # end processing
    #

    end_time = time.time()

    monitoring.to_log_and_console('    computation time = ' + str(end_time - start_time) + ' s', 1)
    monitoring.to_log_and_console('', 1)

    return
