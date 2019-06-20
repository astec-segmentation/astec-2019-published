
import os
import imp
import sys
import time
import math
import platform
import shutil
import subprocess
import numpy as np
from scipy import ndimage as nd

import commonTools
import nomenclature
from CommunFunctions.ImageHandling import SpatialImage, imread, imsave
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
# - channel environment
# - computation environment
# - computation parameters
#
########################################################################################


class FusionChannel(object):

    def __init__(self):
        #
        # raw data directories
        #
        self.path_angle1 = None
        self.path_angle2 = None
        self.path_angle3 = None
        self.path_angle4 = None

        #
        # fused data paths
        #
        self.path_fuse_exp = None

        #
        # temporary_paths
        #
        self.temporary_paths = list()

    def update_main_channel_from_file(self, parameter_file):
        proc = 'FusionChannel.update_main_channel_from_file'
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print ("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = imp.load_source('*', parameter_file)

        self.path_angle1 = nomenclature.replaceFlags(nomenclature.path_rawdata_angle1, parameters).replace("//", "/")
        self.path_angle2 = nomenclature.replaceFlags(nomenclature.path_rawdata_angle2, parameters).replace("//", "/")
        self.path_angle3 = nomenclature.replaceFlags(nomenclature.path_rawdata_angle3, parameters).replace("//", "/")
        self.path_angle4 = nomenclature.replaceFlags(nomenclature.path_rawdata_angle4, parameters).replace("//", "/")

        if not os.path.isdir(self.path_angle1) or not os.path.isdir(self.path_angle2) \
                or not os.path.isdir(self.path_angle3) or not os.path.isdir(self.path_angle4):
            monitoring.to_log_and_console(proc + ": at least one raw data directory for main channel does not exist")
            monitoring.to_log_and_console("- " + str(self.path_angle1))
            monitoring.to_log_and_console("- " + str(self.path_angle2))
            monitoring.to_log_and_console("- " + str(self.path_angle3))
            monitoring.to_log_and_console("- " + str(self.path_angle4))
            monitoring.to_log_and_console("\t Exiting")
            sys.exit(1)

        self.path_fuse_exp = nomenclature.replaceFlags(nomenclature.path_fuse_exp, parameters)
        return

    def update_channel_x_from_file(self, channel_id, parameter_file):
        proc = 'FusionChannel.update_channel_x_from_file'
        if parameter_file is None:
            return False
        if not os.path.isfile(parameter_file):
            print ("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = imp.load_source('*', parameter_file)

        #
        # build paths to raw data
        #

        if hasattr(parameters, 'DIR_RAWDATA_CHANNEL_' + str(channel_id)):
            path_rawdata = os.path.join(nomenclature.FLAG_PATH_EMBRYO,
                                        getattr(parameters, 'DIR_RAWDATA_CHANNEL_' + str(channel_id)))
        elif hasattr(parameters, 'DIR_RAWDATA_CHANNEL' + str(channel_id)):
            path_rawdata = os.path.join(nomenclature.FLAG_PATH_EMBRYO,
                                        getattr(parameters, 'DIR_RAWDATA_CHANNEL' + str(channel_id)))
        elif hasattr(parameters, 'DIR_RAWDATA'):
            path_rawdata = os.path.join(nomenclature.FLAG_PATH_EMBRYO, parameters.DIR_RAWDATA)
        else:
            path_rawdata = os.path.join(nomenclature.FLAG_PATH_EMBRYO, nomenclature.FLAG_DIR_RAWDATA)

        path_rawdata = nomenclature.replaceFlags(path_rawdata, parameters)

        if not os.path.isdir(path_rawdata):
            return False

        if hasattr(parameters, 'DIR_LEFTCAM_STACKZERO_CHANNEL_' + str(channel_id)):
            self.path_angle1 = os.path.join(path_rawdata,
                                            getattr(parameters, 'DIR_LEFTCAM_STACKZERO_CHANNEL_' + str(channel_id)))
        elif hasattr(parameters, 'DIR_LEFTCAM_STACKZERO_CHANNEL' + str(channel_id)):
            self.path_angle1 = os.path.join(path_rawdata,
                                            getattr(parameters, 'DIR_LEFTCAM_STACKZERO_CHANNEL' + str(channel_id)))
        elif hasattr(parameters, 'DIR_LEFTCAM_STACKZERO'):
            self.path_angle1 = os.path.join(path_rawdata, parameters.DIR_LEFTCAM_STACKZERO)
        else:
            self.path_angle1 = os.path.join(path_rawdata, nomenclature.FLAG_DIR_LEFTCAM_STACKZERO)

        if hasattr(parameters, 'DIR_RIGHTCAM_STACKZERO_CHANNEL_' + str(channel_id)):
            self.path_angle2 = os.path.join(path_rawdata,
                                            getattr(parameters, 'DIR_RIGHTCAM_STACKZERO_CHANNEL_' + str(channel_id)))
        elif hasattr(parameters, 'DIR_RIGHTCAM_STACKZERO_CHANNEL' + str(channel_id)):
            self.path_angle2 = os.path.join(path_rawdata,
                                            getattr(parameters, 'DIR_RIGHTCAM_STACKZERO_CHANNEL' + str(channel_id)))
        elif hasattr(parameters, 'DIR_RIGHTCAM_STACKZERO'):
            self.path_angle2 = os.path.join(path_rawdata, parameters.DIR_RIGHTCAM_STACKZERO)
        else:
            self.path_angle2 = os.path.join(path_rawdata, nomenclature.FLAG_DIR_RIGHTCAM_STACKZERO)

        if hasattr(parameters, 'DIR_LEFTCAM_STACKONE_CHANNEL_' + str(channel_id)):
            self.path_angle3 = os.path.join(path_rawdata,
                                            getattr(parameters, 'DIR_LEFTCAM_STACKONE_CHANNEL_' + str(channel_id)))
        elif hasattr(parameters, 'DIR_LEFTCAM_STACKONE_CHANNEL' + str(channel_id)):
            self.path_angle3 = os.path.join(path_rawdata,
                                            getattr(parameters, 'DIR_LEFTCAM_STACKONE_CHANNEL' + str(channel_id)))
        elif hasattr(parameters, 'DIR_LEFTCAM_STACKONE'):
            self.path_angle3 = os.path.join(path_rawdata, parameters.DIR_LEFTCAM_STACKONE)
        else:
            self.path_angle3 = os.path.join(path_rawdata, nomenclature.FLAG_DIR_LEFTCAM_STACKONE)

        if hasattr(parameters, 'DIR_RIGHTCAM_STACKONE_CHANNEL_' + str(channel_id)):
            self.path_angle4 = os.path.join(path_rawdata,
                                            getattr(parameters, 'DIR_RIGHTCAM_STACKONE_CHANNEL_' + str(channel_id)))
        elif hasattr(parameters, 'DIR_RIGHTCAM_STACKONE_CHANNEL' + str(channel_id)):
            self.path_angle4 = os.path.join(path_rawdata,
                                            getattr(parameters, 'DIR_RIGHTCAM_STACKONE_CHANNEL' + str(channel_id)))
        elif hasattr(parameters, 'DIR_RIGHTCAM_STACKONE'):
            self.path_angle4 = os.path.join(path_rawdata, parameters.DIR_RIGHTCAM_STACKONE)
        else:
            self.path_angle4 = os.path.join(path_rawdata, nomenclature.FLAG_DIR_RIGHTCAM_STACKONE)

        self.path_angle1 = nomenclature.replaceFlags(self.path_angle1, parameters)
        self.path_angle2 = nomenclature.replaceFlags(self.path_angle2, parameters)
        self.path_angle3 = nomenclature.replaceFlags(self.path_angle3, parameters)
        self.path_angle4 = nomenclature.replaceFlags(self.path_angle4, parameters)

        if os.path.isdir(self.path_angle1) or os.path.isdir(self.path_angle2) \
                or os.path.isdir(self.path_angle3) or os.path.isdir(self.path_angle4):
            if not os.path.isdir(self.path_angle1) or not os.path.isdir(self.path_angle2) \
                    or not os.path.isdir(self.path_angle3) or not os.path.isdir(self.path_angle4):
                monitoring.to_log_and_console(proc + ": at least one raw data directory for channel '" + str(channel_id)
                                              + "' does not exist")
                monitoring.to_log_and_console("- " + str(self.path_angle1))
                monitoring.to_log_and_console("- " + str(self.path_angle2))
                monitoring.to_log_and_console("- " + str(self.path_angle3))
                monitoring.to_log_and_console("- " + str(self.path_angle4))
                monitoring.to_log_and_console("\t Exiting")
                sys.exit(1)
        else:
            return False

        self.path_fuse_exp = None
        if hasattr(parameters, 'EXP_FUSE_CHANNEL_' + str(channel_id)):
            self.path_fuse_exp = os.path.join(nomenclature.path_fuse, nomenclature.DIR_STAGE_FUSE + '_'
                                              + getattr(parameters, 'EXP_FUSE_CHANNEL_' + str(channel_id)))
        elif hasattr(parameters, 'EXP_FUSE_CHANNEL' + str(channel_id)):
            self.path_fuse_exp = os.path.join(nomenclature.path_fuse, nomenclature.DIR_STAGE_FUSE + '_'
                                              + getattr(parameters, 'EXP_FUSE_CHANNEL' + str(channel_id)))

        return True

    def has_same_raw_data_dirs(self, c):
        if self.path_angle1 != c.path_angle1 and self.path_angle2 != c.path_angle2 \
                and self.path_angle3 != c.path_angle3 and self.path_angle4 != c.path_angle4:
            return False
        return True

    def update_path_fuse_exp(self, c, suffixe):
        self.path_fuse_exp = c.path_fuse_exp + suffixe

    def write_parameters(self, log_file_name, desc=None):
        with open(log_file_name, 'a') as logfile:

            logfile.write("\n")
            if desc is not None:
                logfile.write("- " + str(desc) + " =\n")

            logfile.write('  FusionChannel\n')
            logfile.write('  - path_angle1 = ' + str(self.path_angle1)+'\n')
            logfile.write('  - path_angle2 = ' + str(self.path_angle2)+'\n')
            logfile.write('  - path_angle3 = ' + str(self.path_angle3)+'\n')
            logfile.write('  - path_angle4 = ' + str(self.path_angle4)+'\n')

            logfile.write('  - path_fuse_exp = ' + str(self.path_fuse_exp)+'\n')

            for j in range(0, len(self.temporary_paths)):
                logfile.write('  - temporary path #' + str(j) + ' = ' + str(self.temporary_paths[j]) + '\n')

            logfile.write("\n")
        return

    def print_parameters(self, desc=None):

        print("")
        if desc is not None:
            print("- " + str(desc) + " =")

        print('  FusionChannel')
        print('  - path_angle1 = ' + str(self.path_angle1))
        print('  - path_angle2 = ' + str(self.path_angle2))
        print('  - path_angle3 = ' + str(self.path_angle3))
        print('  - path_angle4 = ' + str(self.path_angle4))

        print('  - path_fuse_exp = ' + str(self.path_fuse_exp))

        for j in range(0, len(self.temporary_paths)):
            print('  - temporary path #' + str(j) + ' = ' + str(self.temporary_paths[j]))

        print("")


class FusionEnvironment(object):

    def __init__(self):

        #
        # Channels
        #
        self.channel = list()

        #
        # raw data file names
        # assumed to be the same for all channels
        #
        self.path_angle1_files = None
        self.path_angle2_files = None
        self.path_angle3_files = None
        self.path_angle4_files = None

        #
        # fused data paths
        #
        self.path_fuse = None
        self.path_fuse_exp_files = None

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

        self.channel[0].update_main_channel_from_file(parameter_file)

        #
        # build other channels
        #
        channel2 = FusionChannel()
        if channel2.update_channel_x_from_file('2', parameter_file) is True:
            if channel2.has_same_raw_data_dirs(self.channel[0]) is False:
                if channel2.path_fuse_exp is None:
                    channel2.update_path_fuse_exp(self.channel[0], '_CHANNEL_2')
                channel2.path_fuse_exp = nomenclature.replaceFlags(channel2.path_fuse_exp, parameters)
                if channel2.path_fuse_exp == self.channel[0].path_fuse_exp:
                    print ("Error: channel #2 result directory is the same than channel #1. Exiting.")
                    sys.exit(1)
                self.channel.append(channel2)

        channel3 = FusionChannel()
        if channel3.update_channel_x_from_file('3', parameter_file) is True:
            if channel3.has_same_raw_data_dirs(self.channel[0]) is False:
                if channel3.path_fuse_exp is None:
                    channel3.update_path_fuse_exp(self.channel[0], '_CHANNEL_3')
                channel3.path_fuse_exp = nomenclature.replaceFlags(channel3.path_fuse_exp, parameters)
                if channel3.path_fuse_exp == self.channel[0].path_fuse_exp:
                    print ("Error: channel #3 result directory is the same than channel #1. Exiting.")
                    sys.exit(1)
                elif channel3.path_fuse_exp == self.channel[1].path_fuse_exp:
                    print ("Error: channel #3 result directory is the same than channel #2. Exiting.")
                    sys.exit(1)
                self.channel.append(channel3)

        self.path_angle1_files = nomenclature.replaceFlags(nomenclature.path_rawdata_angle1_files, parameters)
        self.path_angle2_files = nomenclature.replaceFlags(nomenclature.path_rawdata_angle2_files, parameters)
        self.path_angle3_files = nomenclature.replaceFlags(nomenclature.path_rawdata_angle3_files, parameters)
        self.path_angle4_files = nomenclature.replaceFlags(nomenclature.path_rawdata_angle4_files, parameters)

        self.path_fuse = nomenclature.replaceFlags(nomenclature.path_fuse, parameters)

        self.path_fuse_exp_files = nomenclature.replaceFlags(nomenclature.path_fuse_exp_files, parameters)

        self.path_logdir = nomenclature.replaceFlags(nomenclature.path_fuse_logdir, parameters)
        self.path_history_file = nomenclature.replaceFlags(nomenclature.path_fuse_historyfile, parameters)
        self.path_log_file = nomenclature.replaceFlags(nomenclature.path_fuse_logfile, parameters, start_time)

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('FusionEnvironment\n')

            for c in range(0, len(self.channel)):
                self.channel[c].write_parameters(log_file_name, 'channel #' + str(c))

            logfile.write('- path_angle1_files = ' + str(self.path_angle1_files)+'\n')
            logfile.write('- path_angle2_files = ' + str(self.path_angle2_files)+'\n')
            logfile.write('- path_angle3_files = ' + str(self.path_angle3_files)+'\n')
            logfile.write('- path_angle4_files = ' + str(self.path_angle4_files)+'\n')

            logfile.write('- path_fuse = ' + str(self.path_fuse)+'\n')
            logfile.write('- path_fuse_exp_files = ' + str(self.path_fuse_exp_files)+'\n')

            logfile.write('- path_logdir = ' + str(self.path_logdir) + '\n')
            logfile.write('- path_history_file = ' + str(self.path_history_file)+'\n')
            logfile.write('- path_log_file = ' + str(self.path_log_file)+'\n')
            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('FusionEnvironment')

        for c in range(0, len(self.channel)):
            self.channel[c].print_parameters('channel #' + str(c))

        print('- path_angle1_files = ' + str(self.path_angle1_files))
        print('- path_angle2_files = ' + str(self.path_angle2_files))
        print('- path_angle3_files = ' + str(self.path_angle3_files))
        print('- path_angle4_files = ' + str(self.path_angle4_files))

        print('- path_fuse = ' + str(self.path_fuse))
        print('- path_fuse_exp_files = ' + str(self.path_fuse_exp_files))

        print('- path_logdir = ' + str(self.path_logdir))
        print('- path_history_file = ' + str(self.path_history_file))
        print('- path_log_file = ' + str(self.path_log_file))
        print("")


#
#
#
#
#


class FusionParameters(object):

    def __init__(self):
        #
        # acquisition parameters
        #
        self.acquisition_orientation = 'left'
        self.acquisition_mirrors = False
        self.acquisition_resolution = None

        #
        # Correction of slit lines
        #
        self.acquisition_slit_line_correction = False

        #
        # fused image parameters
        #
        self.target_resolution = (0.3, 0.3, 0.3)

        #
        # Cropping of acquisition images (before fusion)
        #
        self.acquisition_cropping = True
        self.acquisition_cropping_margin_x_0 = 40
        self.acquisition_cropping_margin_x_1 = 40
        self.acquisition_cropping_margin_y_0 = 40
        self.acquisition_cropping_margin_y_1 = 40

        #
        # Registration parameters
        #
        self.pre_registration = commonTools.RegistrationParameters()
        self.pre_registration.prefix = 'fusion_preregistration_'
        self.pre_registration.compute_registration = False
        self.pre_registration.transformation_type = 'translation'

        self.registration = commonTools.RegistrationParameters()
        self.registration.prefix = 'fusion_registration_'

        #
        # Cropping of fused image (after fusion)
        #
        self.fusion_cropping = True
        self.fusion_cropping_margin_x_0 = 40
        self.fusion_cropping_margin_x_1 = 40
        self.fusion_cropping_margin_y_0 = 40
        self.fusion_cropping_margin_y_1 = 40

        #
        # images suffixes/formats
        #
        self.result_image_suffix = 'inr'
        self.default_image_suffix = 'inr'

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('FusionParameters\n')

            logfile.write('- acquisition_orientation = '+str(self.acquisition_orientation)+'\n')
            logfile.write('- acquisition_mirrors     = '+str(self.acquisition_mirrors)+'\n')
            logfile.write('- acquisition_resolution  = '+str(self.acquisition_resolution)+'\n')

            logfile.write('- acquisition_slit_line_correction = '+str(self.acquisition_slit_line_correction)+'\n')

            logfile.write('- target_resolution  = '+str(self.target_resolution)+'\n')

            logfile.write('- acquisition_cropping = '+str(self.acquisition_cropping)+'\n')
            logfile.write('- acquisition_cropping_margin_x_0 = '+str(self.acquisition_cropping_margin_x_0)+'\n')
            logfile.write('- acquisition_cropping_margin_x_1 = '+str(self.acquisition_cropping_margin_x_1)+'\n')
            logfile.write('- acquisition_cropping_margin_y_0 = '+str(self.acquisition_cropping_margin_y_0)+'\n')
            logfile.write('- acquisition_cropping_margin_y_1 = '+str(self.acquisition_cropping_margin_y_1)+'\n')

            self.pre_registration.write_parameters(log_file_name)
            self.registration.write_parameters(log_file_name)

            logfile.write('- fusion_cropping = '+str(self.fusion_cropping)+'\n')
            logfile.write('- fusion_cropping_margin_x_0 = '+str(self.fusion_cropping_margin_x_0)+'\n')
            logfile.write('- fusion_cropping_margin_x_1 = '+str(self.fusion_cropping_margin_x_1)+'\n')
            logfile.write('- fusion_cropping_margin_y_0 = '+str(self.fusion_cropping_margin_y_0)+'\n')
            logfile.write('- fusion_cropping_margin_y_1 = '+str(self.fusion_cropping_margin_y_1)+'\n')

            logfile.write('- result_image_suffix = '+str(self.result_image_suffix) + '\n')
            logfile.write('- default_image_suffix = '+str(self.default_image_suffix) + '\n')

            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('FusionParameters')

        print('- acquisition_orientation = '+str(self.acquisition_orientation))
        print('- acquisition_mirrors     = '+str(self.acquisition_mirrors))
        print('- acquisition_resolution  = '+str(self.acquisition_resolution))

        print('- acquisition_slit_line_correction = '+str(self.acquisition_slit_line_correction))

        print('- target_resolution  = '+str(self.target_resolution))

        print('- acquisition_cropping = '+str(self.acquisition_cropping))
        print('- acquisition_cropping_margin_x_0 = '+str(self.acquisition_cropping_margin_x_0))
        print('- acquisition_cropping_margin_x_1 = '+str(self.acquisition_cropping_margin_x_1))
        print('- acquisition_cropping_margin_y_0 = '+str(self.acquisition_cropping_margin_y_0))
        print('- acquisition_cropping_margin_y_1 = '+str(self.acquisition_cropping_margin_y_1))

        self.pre_registration.print_parameters()
        self.registration.print_parameters()

        print('- fusion_cropping = '+str(self.fusion_cropping))
        print('- fusion_cropping_margin_x_0 = '+str(self.fusion_cropping_margin_x_0))
        print('- fusion_cropping_margin_x_1 = '+str(self.fusion_cropping_margin_x_1))
        print('- fusion_cropping_margin_y_0 = '+str(self.fusion_cropping_margin_y_0))
        print('- fusion_cropping_margin_y_1 = '+str(self.fusion_cropping_margin_y_1))

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
        # acquisition parameters
        #
        if hasattr(parameters, 'raw_ori'):
            if parameters.raw_ori is not None:
                self.acquisition_orientation = parameters.raw_ori
        elif hasattr(parameters, 'acquisition_orientation'):
            if parameters.acquisition_orientation is not None:
                self.acquisition_orientation = parameters.acquisition_orientation

        if hasattr(parameters, 'raw_mirrors'):
            if parameters.raw_mirrors is not None:
                self.acquisition_mirrors = parameters.raw_mirrors
        elif hasattr(parameters, 'acquisition_mirrors'):
            if parameters.acquisition_mirrors is not None:
                self.acquisition_mirrors = parameters.acquisition_mirrors

        if hasattr(parameters, 'raw_resolution'):
            if parameters.raw_resolution is not None:
                if type(parameters.raw_resolution) is tuple or type(parameters.raw_resolution) is list:
                    if len(parameters.raw_resolution) == 3:
                        self.acquisition_resolution = parameters.raw_resolution
                    else:
                        print("Error in'" + parameter_file + "'")
                        print("\t 'raw_resolution' has length " + str(len(parameters.raw_resolution))
                              + " instead of 3.")
                        print("\t Exiting.")
                        sys.exit(1)
                else:
                    print("Error in'" + parameter_file + "'")
                    print("\t type of 'raw_resolution' (" + str(type(parameters.raw_resolution))
                          + ") is not handled")
                    print("\t Exiting.")
                    sys.exit(1)
        elif hasattr(parameters, 'acquisition_resolution'):
            if parameters.acquisition_resolution is not None:
                if type(parameters.acquisition_resolution) is tuple or type(parameters.acquisition_resolution) is list:
                    if len(parameters.acquisition_resolution) == 3:
                        self.acquisition_resolution = parameters.acquisition_resolution
                    else:
                        print("Error in'" + parameter_file + "'")
                        print("\t 'acquisition_resolution' has length " + str(len(parameters.acquisition_resolution))
                              + " instead of 3.")
                        print("\t Exiting.")
                        sys.exit(1)
                else:
                    print("Error in'" + parameter_file + "'")
                    print("\t type of 'acquisition_resolution' (" + str(type(parameters.acquisition_resolution))
                          + ") is not handled")
                    print("\t Exiting.")
                    sys.exit(1)

        #
        # correction of slit lines
        #
        if hasattr(parameters, 'acquisition_slit_line_correction'):
            if parameters.acquisition_slit_line_correction is not None:
                self.acquisition_slit_line_correction = parameters.acquisition_slit_line_correction

        #
        # fused image parameters
        #
        if hasattr(parameters, 'target_resolution'):
            if parameters.target_resolution is not None:
                self.target_resolution = parameters.target_resolution

        #
        # Cropping of acquisition images (before fusion)
        #
        if hasattr(parameters, 'raw_crop'):
            if parameters.raw_crop is not None:
                self.acquisition_cropping = parameters.raw_crop
        if hasattr(parameters, 'raw_margin_x_0'):
            if parameters.raw_margin_x_0 is not None:
                self.acquisition_cropping_margin_x_0 = parameters.raw_margin_x_0
        if hasattr(parameters, 'raw_margin_x_1'):
            if parameters.raw_margin_x_1 is not None:
                self.acquisition_cropping_margin_x_1 = parameters.raw_margin_x_1
        if hasattr(parameters, 'raw_margin_y_0'):
            if parameters.raw_margin_y_0 is not None:
                self.acquisition_cropping_margin_y_0 = parameters.raw_margin_y_0
        if hasattr(parameters, 'raw_margin_y_1'):
            if parameters.raw_margin_y_1 is not None:
                self.acquisition_cropping_margin_y_1 = parameters.raw_margin_y_1

        #
        # registration parameters
        self.pre_registration.update_from_file(parameter_file)
        self.registration.update_from_file(parameter_file)

        #
        # Cropping of fused image (after fusion)
        #
        if hasattr(parameters, 'fusion_crop'):
            if parameters.fusion_crop is not None:
                self.fusion_cropping = parameters.fusion_crop
        if hasattr(parameters, 'fusion_margin_x_0'):
            if parameters.fusion_margin_x_0 is not None:
                self.fusion_cropping_margin_x_0 = parameters.fusion_margin_x_0
        if hasattr(parameters, 'fusion_margin_x_1'):
            if parameters.fusion_margin_x_1 is not None:
                self.fusion_cropping_margin_x_1 = parameters.fusion_margin_x_1
        if hasattr(parameters, 'fusion_margin_y_0'):
            if parameters.fusion_margin_y_0 is not None:
                self.fusion_cropping_margin_y_0 = parameters.fusion_margin_y_0
        if hasattr(parameters, 'fusion_margin_y_1'):
            if parameters.fusion_margin_y_1 is not None:
                self.fusion_cropping_margin_y_1 = parameters.fusion_margin_y_1

        #
        # images suffixes/formats
        #
        if hasattr(parameters, 'RESULT_IMAGE_SUFFIX_FUSE'):
            if parameters.RESULT_IMAGE_SUFFIX_FUSE is not None:
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


__extension_to_be_converted__ = ['.h5', '.tif', '.tiff', '.TIF', '.TIFF']
__extension_with_resolution__ = ['.inr', '.inr.gz', '.mha', '.mha.gz']


def _read_image_name(data_path, temporary_path, file_name, resolution, default_extension='inr'):
    """
    Read an image. Eventually, unzip a compressed file, and convert the image
    to a format known by executables
    :param data_path: path to data directory
    :param temporary_path: directory for temporary file
    :param file_name: image file
    :param resolution: resolution of the result image
            required to write the output image with the right resolution
    :return:
    """

    proc = "_read_image_name"

    #
    # test whether the extension is zip
    #
    f = file_name
    full_name = os.path.join(data_path, f)

    if f[len(f)-4:len(f)] == '.zip':

        prefix = f[0:len(f)-4]

        #
        # maybe the file has already be unzipped
        #
        file_names = []
        for f in os.listdir(temporary_path):
            if len(f) <= len(prefix):
                pass
            if f[0:len(prefix)] == prefix:
                if f[len(prefix):len(f)] in commonTools.recognized_extensions:
                    file_names.append(f)

        if len(file_names) > 1:
            monitoring.to_log_and_console(proc + ": already several images with name '"
                                          + str(prefix) + "' were found in '" + str(temporary_path) + "'")
            monitoring.to_log_and_console("\t " + str(file_names))
            monitoring.to_log_and_console("\t Exiting")
            sys.exit(1)

        elif len(file_names) == 0:
            #
            # unzipping
            #
            monitoring.to_log_and_console("    .. unzipping '" + str(f) + "'", 2)
            #
            # there are issues with unzip
            # seems to occur when files are zipped with zip 3.0
            #
            if platform.system() == 'Linux':
                command_line = 'unzip ' + os.path.join(data_path, f) + ' -d ' + str(temporary_path)
            elif platform.system() == 'Darwin':
                command_line = 'tar xvf ' + os.path.join(data_path, f) + ' -C ' + str(temporary_path)
            else:
                command_line = 'unzip ' + os.path.join(data_path, f) + ' -d ' + str(temporary_path)
            if monitoring.verbose >= 3 or monitoring.debug > 0:
                monitoring.to_log("* Launch: " + command_line)
                with open(monitoring.logfile, 'a') as logfile:
                    subprocess.call(command_line, shell=True, stdout=logfile, stderr=subprocess.STDOUT)
            else:
                subprocess.call(command_line, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

            #
            # parsing again the temporay directory
            #
            file_names = []
            for f in os.listdir(temporary_path):
                if len(f) <= len(prefix):
                    pass
                if f[0:len(prefix)] == prefix:
                    file_names.append(f)

        if len(file_names) == 0:
            monitoring.to_log_and_console(proc + ": no image with name '" + str(prefix)
                                          + "' was found in '" + str(temporary_path) + "'")
            monitoring.to_log_and_console("\t Exiting")
            sys.exit(1)

        if len(file_names) > 1:
            monitoring.to_log_and_console(proc + ": several images with name '"
                                          + str(prefix) + "' were found in '" + str(temporary_path) + "'")
            monitoring.to_log_and_console("\t " + str(file_names))
            monitoring.to_log_and_console("\t Exiting")
            sys.exit(1)
        #
        #
        #
        f = file_names[0]
        full_name = os.path.join(temporary_path, f)

    #
    # test whether the file has to be converted into a more 'readable' format
    # if yes, set the resolution if required
    #

    file_has_been_converted = False
    for extension in __extension_to_be_converted__:
        if f[len(f)-len(extension):len(f)] == extension:
            prefix = f[0:len(f) - len(extension)]

            #
            # new file name
            # check whether it has already been converted
            #
            new_full_name = os.path.join(temporary_path, prefix) + '.' + str(default_extension)
            if not os.path.isfile(new_full_name):
                monitoring.to_log_and_console("    .. converting '" + str(f) + "'", 2)
                image = imread(full_name)
                if type(resolution) is tuple and len(resolution) == 3:
                    image.resolution = resolution
                    monitoring.to_log("    * resolution of '" + full_name + "' has been set to "
                                      + str(image.resolution))
                elif type(resolution) is list and len(resolution) == 3:
                    image.resolution = (resolution(0), resolution(1), resolution(2))
                    monitoring.to_log("    * resolution of '" + full_name + "' has been set to "
                                      + str(image.resolution))
                else:
                    monitoring.to_log("    * resolution of '" + full_name + "' is " + str(image.resolution)
                                      + "(default/read values)")
                #
                # remove unzipped file to avoid having two files in the directory
                # verify that it is not the input file!
                #
                if os.path.dirname(full_name) == temporary_path:
                    os.remove(full_name)

                #
                # save converted file
                #
                imsave(new_full_name, image)
                file_has_been_converted = True

            full_name = new_full_name
            break

    #
    # test whether the input format is supposed to have the resolution set
    #

    if file_has_been_converted is False:
        for extension in __extension_with_resolution__:
            if f[len(f) - len(extension):len(f)] == extension:
                file_has_been_converted = True
                break

    #
    # if no conversion occurs, the resolution has not been set yet
    #

    if file_has_been_converted is False:
        if (type(resolution) is tuple or type(resolution) is list) and len(resolution) == 3:
            monitoring.to_log_and_console("    .. changing resolution '" + str(f) + "'", 2)
            image = imread(full_name)
            if type(resolution) is tuple and len(resolution) == 3:
                image.resolution = resolution
                monitoring.to_log("    * resolution of '" + full_name + "' has been set to " + str(image.resolution))
            elif type(resolution) is list and len(resolution) == 3:
                image.resolution = (resolution(0), resolution(1), resolution(2))
                monitoring.to_log("    * resolution of '" + full_name + "' has been set to " + str(image.resolution))
            imsave(full_name, image)
        else:
            monitoring.to_log("    * resolution of '" + full_name + "' is left unchanged")

    return full_name


def _analyze_data_directory(data_dir):
    """
    Parse a directory containing images
    :param data_dir:
    :return:
    1. the common prefix of image file names
    2. the number of characters used to encoded the variable part
       (time points)
    3. the list of the variable parts
    4. the common suffix of image file names
       may be longer than just the file extension
    """

    proc = "_analyze_data_directory"
    images = []
    extensions = []
    #
    # recognize images and extensions
    #
    for f in os.listdir(data_dir):
        for e in commonTools.recognized_extensions:
            if f[len(f)-len(e):len(f)] == e:
                if e not in extensions:
                    extensions.append(e)
                    if len(extensions) > 1:
                        print proc + ": several image extensions were found in '" + data_dir + "'"
                        print "\t -> " + str(extensions)
                        print "\t Exiting."
                        sys.exit(1)
                images.append(f)

    if len(images) == 0:
        print proc + ": no images were found in '" + data_dir + "'"
        print "\t Exiting."
        sys.exit(1)

    #
    # one image case
    #

    if len(images) == 1:
        suffix = extensions[0]
        time_length = 0
        im = images[0]
        length = len(im) - 1 - len(suffix)
        for i in range(0, 3):
            if '0' <= im[length-i] <= '9':
                time_length += 1
            else:
                break
        prefix = im[0:len(im) - time_length - len(suffix)]
        time_points = im[len(im) - time_length - len(suffix):len(im) - len(suffix)]
        return prefix, time_length, time_points, suffix

    #
    # several images
    # 1. check that image names are of the same length
    # 2. get prefix = common part at beginning
    # 3. get suffix = common part at end
    # 4. get length for time point encoding
    # 5. get list of time points
    #

    for i in range(1, len(images)):
        if len(images[0]) != len(images[i]):
            print proc + ": image names are of different lengths in '" + data_dir + "'"
            print "\t -> " + images[0] + ", " + images[i]
            print "\t Exiting."
            sys.exit(1)

    prefix = ''
    for j in range(0, len(images[0])):
        ok = True
        for i in range(1, len(images)):
            if images[0][j] != images[i][j]:
                ok = False
                break
        if ok is True:
            prefix += images[0][j]
        else:
            break

    suffix = ''
    for j in range(len(images[0]) - 1, -1, -1):
        ok = True
        for i in range(1, len(images)):
            if images[0][j] != images[i][j]:
                ok = False
                break
        if ok is True:
            suffix += images[0][j]
        else:
            break
    suffix = suffix[::-1]

    time_length = len(images[0]) - len(prefix) - len(suffix)

    time_points = []
    for i in range(0, len(images)):
        time_points.append(images[i][len(images[i]) - time_length - len(suffix):len(images[i]) - len(suffix)])

    return prefix, time_length, time_points, suffix


########################################################################################
#
# cropping
#
########################################################################################


def _crop_bounding_box(the_image):
    """
    Compute a bounding box to crop an image (ie extract a subimage)
    :param the_image:
    :return:
    """

    #
    # build a 2D binary image from the MIP projection
    #

    the_selection = commonTools.add_suffix(the_image, "_cropselection")
    cpp_wrapping.mip_projection_for_crop(the_image, the_selection, None, monitoring)

    #
    # read input image
    #
    selection = imread(the_selection)

    #
    # the get the connected component (4-connectivity)
    # there should have only two components: the background and the selected component
    #
    cc_image, cc_n = nd.label(selection)

    #
    # compute the volumes of each connected component
    # and create a dictionary of tuples (label, volume)
    #
    labels = np.unique(cc_image)
    volumes = nd.sum(np.ones_like(cc_image), cc_image, index=np.int16(labels))
    dict_volumes = dict(zip(labels, volumes))

    #
    # remove the background
    # then get the label associated to the largest connected component
    dict_volumes.pop(0)
    max_label = dict_volumes.keys()[np.argmax(dict_volumes.values())]

    #
    # get the bounding boxes for all objects
    # it is not necessary to searched for all labels
    # seems that there is no bounding box computed for label #0
    #
    # boundingBoxes = nd.find_objects(ccImage, max_label=maxLabel)
    # maxBox = boundingBoxes[int(maxLabel)-1]
    #
    max_box = nd.find_objects(cc_image, max_label=max_label)[int(max_label)-1]

    del cc_image

    return max_box


def _crop_disk_image(the_image, res_image, the_max_box=None,
                     margin_x_0=40, margin_x_1=40, margin_y_0=40, margin_y_1=40):
    """
    Crop an image on disk in XY plane
    :param the_image:
    :param res_image:
    :param the_max_box:
    :param margin_x_0:
    :param margin_x_1:
    :param margin_y_0:
    :param margin_y_1:
    :return:
    """

    #
    # 2D bounding box
    #
    if the_max_box is None:
        max_box = _crop_bounding_box(the_image)
    else:
        max_box = the_max_box

    #
    # 2D bounding box + margin
    #
    image = imread(the_image)

    xmin = max(max_box[0].start - margin_x_0, 0)
    xmax = min(image.shape[0], max_box[0].stop + margin_x_1)
    ymin = max(max_box[1].start - margin_y_0, 0)
    ymax = min(image.shape[1], max_box[1].stop + margin_y_1)

    new_box = (slice(xmin, xmax, None),
               slice(ymin, ymax, None),
               slice(0, image.shape[2]))

    new_image = SpatialImage(image[new_box])
    new_image._set_resolution(image._get_resolution())

    imsave(res_image, new_image)

    monitoring.to_log_and_console("       crop from [0," + str(image.shape[0]) + "]x[0,"
                                  + str(image.shape[1]) + "] to [" + str(xmin) + ","
                                  + str(xmax) + "]x[" + str(ymin) + "," + str(ymax) + "]", 2)

    return


########################################################################################
#
# computation of a rotation matrix
#
########################################################################################


def _axis_rotation_matrix(axis, angle, min_space=None, max_space=None):
    """ Return the transformation matrix from the axis and angle necessary
    this is a rigid transformation (rotation) that preserves the center of
    the field of view
    axis : axis of rotation ("X", "Y" or "Z")
    angle : angle of rotation (in degree)
    min_space : coordinates of the bottom point (usually (0, 0, 0))
    max_space : coordinates of the top point (usually im shape)
    """
    i = np.linalg.inv
    d = np.dot
    if axis not in ["X", "Y", "Z"]:
        raise Exception("Unknown axis : " + str(axis))
    rads = math.radians(angle)
    s = math.sin(rads)
    c = math.cos(rads)

    centering = np.identity(4)
    if min_space is None and max_space is not None:
        min_space = np.array([0., 0., 0.])

    if max_space is not None:
        space_center = (max_space-min_space)/2.
        offset = -1.*space_center
        centering[:3, 3] = offset

    rot = np.identity(4)
    if axis == "X":
        rot = np.array([[1., 0., 0., 0.],
                        [0., c, -s, 0.],
                        [0., s, c, 0.],
                        [0., 0., 0., 1.]])
    elif axis == "Y":
        rot = np.array([[c,   0., s,  0.],
                        [0., 1., 0., 0.],
                        [-s, 0., c, 0.],
                        [0., 0., 0., 1.]])

    elif axis == "Z":
        rot = np.array([[c, -s,  0., 0.],
                        [s, c, 0., 0.],
                        [0., 0., 1., 0.],
                        [0., 0., 0., 1.]])

    return d(i(centering), d(rot, centering))


def _init_rotation_matrix(axis, angle, ref_center=None, flo_center=None):

    if axis not in ["X", "Y", "Z"]:
        raise Exception("Unknown axis : " + str(axis))
    rads = math.radians(angle)
    s = math.sin(rads)
    c = math.cos(rads)

    rot = np.identity(3)
    if axis == "X":
        rot = np.array([[1., 0., 0.],
                        [0., c, -s],
                        [0., s, c]])
    elif axis == "Y":
        rot = np.array([[c, 0., s],
                        [0., 1., 0.],
                        [-s, 0., c]])
    elif axis == "Z":
        rot = np.array([[c, -s, 0.],
                        [s, c, 0.],
                        [0., 0., 1.]])

    if ref_center is not None:
        if flo_center is not None:
            trs = flo_center - np.dot(rot, ref_center)
        else:
            trs = ref_center - np.dot(rot, ref_center)
    else:
        if flo_center is not None:
            trs = flo_center - np.dot(rot, flo_center)
        else:
            trs = np.array[0., 0., 0.]

    mat = np.identity(4)
    mat[0:3, 0:3] = rot
    (mat.T)[3:4, 0:3] = trs

    return mat


########################################################################################
#
# function for the ad hoc computation of weights
# for the linear combination of images of the 4 cameras
#
########################################################################################


def _histogram(image, nbins=256):
    """Return histogram of image.

        Unlike `np.histogram`, this function returns the centers of bins and
        does not rebin integer arrays. For integer arrays, each integer value has
        its own bin, which improves speed and intensity-resolution.

        Parameters
        ----------
        image : array
        Input image.
        nbins : int
        Number of bins used to calculate histogram. This value is ignored for
        integer arrays.

        Returns
        -------
        hist : array
        The values of the histogram.
        bin_centers : array
        The values at the center of the bins.
        """

    proc = "_histogram"
    if not isinstance(image, SpatialImage):
        print proc + ": argument image is not an ndarray"
        return

    # For integer types, histogramming with bincount is more efficient.
    if np.issubdtype(image.dtype, np.integer):
        offset = 0
        if np.min(image) < 0:
            offset = np.min(image)
        hist = np.bincount(image.ravel() - offset)
        bin_centers = np.arange(len(hist)) + offset

        # clip histogram to start with a non-zero bin
        idx = np.nonzero(hist)[0][0]
        return hist[idx:], bin_centers[idx:]
    else:
        hist, bin_edges = np.histogram(image.flat, nbins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
        return hist, bin_centers


def _threshold_otsu(image, nbins=256):
    """Return threshold value based on Otsu's method.

        Parameters
        ----------
        image : array
        Input image.
        nbins : int
        Number of bins used to calculate histogram. This value is ignored for
        integer arrays.

        Returns
        -------
        threshold : float
        Threshold value.

        References
        ----------
        .. [1] Wikipedia, http://en.wikipedia.org/wiki/Otsu's_Method
        """

    proc = "_threshold_otsu"
    if not isinstance(image, SpatialImage):
        print proc + ": argument image is not an ndarray"
        return

    hist, bin_centers = _histogram(image, nbins)
    hist = hist.astype(float)

    # class probabilities for all possible thresholds
    weight1 = np.cumsum(hist)
    weight2 = np.cumsum(hist[::-1])[::-1]
    # class means for all possible thresholds
    mean1 = np.cumsum(hist * bin_centers) / weight1
    mean2 = (np.cumsum((hist * bin_centers)[::-1]) / weight2[::-1])[::-1]

    # Clip ends to align class 1 and class 2 variables:
    # The last value of `weight1`/`mean1` should pair with zero values in
    # `weight2`/`mean2`, which do not exist.
    variance12 = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2

    idx = np.argmax(variance12)
    threshold = bin_centers[:-1][idx]
    return threshold


def _exp_func(x, length=500, speed=5):
    """ Decay function used to take into account the remotness to the camera
    x : value to compute
    length : lenght of the function
    speed : speed of the function
    """

    return .1+np.exp(-((np.float32(x)*speed)/length))


def _build_mask(image, direction):
    """Return the mask on a given image from the decay function
    im : intensity image (SpatialImage)
    direction : if True the camera is in the side of the first slices in Z
    """

    proc = "_build_mask"
    if not isinstance(image, SpatialImage):
        print proc + ": argument image is not an ndarray"
        return

    th = _threshold_otsu(image)
    im_th = np.zeros_like(image)
    im_th[image > th] = 1
    if direction is False:
        im_th = im_th[:, :, -1::-1]
    im_th_sum = np.cumsum(im_th, axis=2)
    if direction is False:
        im_th_sum = im_th_sum[:, :, -1::-1]
    mask = _exp_func(im_th_sum, np.max(im_th_sum))
    return mask


########################################################################################
#
#
#
########################################################################################

def _linear_registration(path_ref, path_flo, path_output, path_output_trsf, path_init_trsf=None,
                         parameters=None):
    if parameters is not None:
        cpp_wrapping.linear_registration(path_ref, path_flo, path_output, path_output_trsf, path_init_trsf,
                                         py_hl=parameters.pyramid_highest_level,
                                         py_ll=parameters.pyramid_lowest_level,
                                         transformation_type=parameters.transformation_type,
                                         transformation_estimator=parameters.transformation_estimation_type,
                                         lts_fraction=parameters.lts_fraction,
                                         normalization=parameters.normalization,
                                         monitoring=monitoring)
    else:
        cpp_wrapping.linear_registration(path_ref, path_flo, path_output, path_output_trsf, path_init_trsf,
                                         monitoring=monitoring)


########################################################################################
#
#
#
########################################################################################

def fusion_process(input_image_list, fused_image, channel, parameters):
    """
    
    :param input_image_list: a list of list of images to be fused. One list per channel
    :param fused_image: a generic name for the fusion result
    :param channel:
    :param parameters:
    :return:
    """

    proc = 'fusion_process'

    if monitoring.debug > 1:
        print ""
        print proc + " was called with:"
        print "- input_image_list = " + str(input_image_list)
        print "- fused_image = " + str(fused_image)
        for c in range(0, len(channel)):
            channel[c].print_parameters('channel #' + str(c))
        print ""

    #
    # nothing to do if the fused image exists
    #

    do_something = False
    for c in range(0, len(channel)):
        if os.path.isfile(os.path.join(channel[c].path_fuse_exp, fused_image)):
            if monitoring.forceResultsToBeBuilt is False:
                monitoring.to_log_and_console('    fused channel #' + str(c) + ' already existing', 2)
            else:
                monitoring.to_log_and_console('    fused channel #' + str(c) + ' already existing, but forced', 2)
                do_something = True
        else:
            do_something = True

    if do_something is False:
        return

    #
    # how to copy a list:
    # NB: 'res_images = inputImages' acts like pointers
    #

    res_image_list = input_image_list[:]

    #
    # slit line correction
    # these corrections are done on original data (without resampling) on channel[0]
    # the compute corrections are then used for the other channels
    # Crop could be done beforehand to reduce the computational burden
    #

    if parameters.acquisition_slit_line_correction is True:

        the_image_list = res_image_list[:]
        res_image_list = list()
        corrections = list()

        #
        # build the file names
        #

        for c in range(0, len(channel)):

            the_images = the_image_list[c]
            res_images = []

            for i in range(0, len(the_images)):
                res_images.append(commonTools.add_suffix(input_image_list[c][i], "_line_corrected",
                                                         new_dirname=channel[c].temporary_paths[i],
                                                         new_extension=parameters.default_image_suffix))
                if c == 0:
                    corrections.append(commonTools.add_suffix(input_image_list[c][i], "_line_corrected",
                                                              new_dirname=channel[c].temporary_paths[i],
                                                              new_extension='.txt'))
            res_image_list.append(res_images)

        #
        # is there something to do ?
        # check whether one corrected line image is missing
        # for channel #0, check also whether the correction file is missing
        #

        do_something = [False] * len(input_image_list[0])

        for c in range(0, len(channel)):

            the_images = the_image_list[c]
            res_images = res_image_list[c]

            for i in range(0, len(the_images)):
                if os.path.isfile(res_images[i]):
                    if monitoring.forceResultsToBeBuilt is True:
                        do_something[i] = True
                else:
                    do_something[i] = True
                if c == 0:
                    if os.path.isfile(corrections[i]):
                        if monitoring.forceResultsToBeBuilt is True:
                            do_something[i] = True
                    else:
                        do_something[i] = True

        #
        # process
        # corrections are computed on channel #0
        # and are used for other channels
        #

        for c in range(0, len(channel)):

            the_images = the_image_list[c]
            res_images = res_image_list[c]

            for i in range(0, len(the_images)):

                monitoring.to_log_and_console("    .. correcting slit lines of '"
                                              + the_images[i].split(os.path.sep)[-1] + "'", 2)

                if do_something[i] is False:
                    monitoring.to_log_and_console("       nothing to do", 2)
                    continue

                if c == 0:
                    cpp_wrapping.slitline_correction(the_images[i], res_images[i],
                                                     output_corrections=corrections[i],
                                                     monitoring=monitoring)
                else:
                    cpp_wrapping.slitline_correction(the_images[i], res_images[i],
                                                     input_corrections=corrections[i],
                                                     monitoring=monitoring)

    #
    # to do: linear filtering to compensate for resolution change
    # for a change of voxel size from x0 to x1
    # smooth with a Gaussian of sigma = \sqrt(2)^(ln(x0/x1) / ln(2))
    #

    #
    # first change of resolution
    # - for X and Y: target resolution (supposed to be larger than original)
    # - for Z: original resolution (supposed to be larger than target)
    #

    the_image_list = res_image_list[:]
    res_image_list = list()

    for c in range(0, len(channel)):

        the_images = the_image_list[c]
        res_images = []

        #
        # build the file names
        #

        for i in range(0, len(the_images)):
            res_images.append(commonTools.add_suffix(input_image_list[c][i], "_resample",
                                                     new_dirname=channel[c].temporary_paths[i],
                                                     new_extension=parameters.default_image_suffix))
        res_image_list.append(res_images)

        #
        # process
        #

        for i in range(0, len(the_images)):

            im = imread(the_images[i])
            if type(parameters.target_resolution) == int or type(parameters.target_resolution) == float:
                resampling_resolution = [parameters.target_resolution, parameters.target_resolution, im.voxelsize[2]]
            elif (type(parameters.target_resolution) == list or type(parameters.target_resolution) == tuple) \
                    and len(parameters.target_resolution) == 3:
                resampling_resolution = [parameters.target_resolution[0],
                                         parameters.target_resolution[1], im.voxelsize[2]]
            else:
                monitoring.to_log_and_console(proc+': unable to set target resolution for first resampling', 0)
                monitoring.to_log_and_console("\t target resolution was '"+str(parameters.target_resolution)+"'", 0)
                monitoring.to_log_and_console("\t image resolution was '" + str(im.voxelsize) + "'", 0)
                monitoring.to_log_and_console("Exiting.", 0)
                sys.exit(1)
            del im

            monitoring.to_log_and_console("    .. resampling '" + the_images[i].split(os.path.sep)[-1]
                                          + "' at " + str(resampling_resolution), 2)
            if not os.path.isfile(res_images[i]) or monitoring.forceResultsToBeBuilt is True:
                cpp_wrapping.apply_transformation(the_images[i], res_images[i], the_transformation=None,
                                                  template_image=None,
                                                  voxel_size=resampling_resolution, interpolation_mode='linear',
                                                  monitoring=monitoring)
            else:
                monitoring.to_log_and_console("       already existing", 2)

    #
    # 2D crop of resampled acquisition images
    #

    if parameters.acquisition_cropping is True:

        the_image_list = res_image_list[:]
        res_image_list = list()

        #
        # build the file names
        #

        for c in range(0, len(channel)):

            the_images = the_image_list[c]
            res_images = []

            for i in range(0, len(the_images)):
                res_images.append(commonTools.add_suffix(input_image_list[c][i], "_crop",
                                                         new_dirname=channel[c].temporary_paths[i],
                                                         new_extension=parameters.default_image_suffix))
            res_image_list.append(res_images)

        #
        # is there something to do ?
        # check whether one cropped image is missing
        #

        do_something = [False] * len(input_image_list[0])

        for c in range(0, len(channel)):

            the_images = the_image_list[c]
            res_images = res_image_list[c]

            for i in range(0, len(the_images)):
                if os.path.isfile(res_images[i]):
                    if monitoring.forceResultsToBeBuilt is True:
                        do_something[i] = True
                else:
                    do_something[i] = True

        #
        # process
        # bounding box are computed on channel #0
        # and are used for other channels
        #

        box_list = list()

        for c in range(0, len(channel)):

            the_images = the_image_list[c]
            res_images = res_image_list[c]

            for i in range(0, len(the_images)):

                monitoring.to_log_and_console("    .. cropping '"
                                              + the_images[i].split(os.path.sep)[-1] + "'", 2)

                if do_something[i] is False:
                    monitoring.to_log_and_console("       nothing to do", 2)
                    continue

                if c == 0:
                    box = _crop_bounding_box(the_images[i])
                    box_list.append(box)
                else:
                    box = box_list[i]

                if not os.path.isfile(res_images[i]) or monitoring.forceResultsToBeBuilt is True:
                    _crop_disk_image(the_images[i], res_images[i], box,
                                     parameters.acquisition_cropping_margin_x_0,
                                     parameters.acquisition_cropping_margin_x_1,
                                     parameters.acquisition_cropping_margin_y_0,
                                     parameters.acquisition_cropping_margin_y_1)
                else:
                    monitoring.to_log_and_console("       already existing", 2)

    #
    # Mirroring of 'right' images if required
    #

    if parameters.acquisition_mirrors is False:

        the_image_list = res_image_list[:]
        res_image_list = list()

        for c in range(0, len(channel)):

            the_images = the_image_list[c]
            res_images = []

            #
            # build the file names
            #

            for i in range(0, len(the_images)):
                if i == 0 or i == 2:
                    res_images.append(the_images[i])
                else:
                    res_images.append(commonTools.add_suffix(input_image_list[c][i], "_mirror",
                                                             new_dirname=channel[c].temporary_paths[i],
                                                             new_extension=parameters.default_image_suffix))
            res_image_list.append(res_images)

            #
            # process
            #

            for i in range(0, len(the_images)):

                if i == 0 or i == 2:
                    continue
                monitoring.to_log_and_console("    .. mirroring  #" + str(i) + " '"
                                              + the_images[i].split(os.path.sep)[-1], 2)
                if not os.path.isfile(res_images[i]) or monitoring.forceResultsToBeBuilt is True:
                    the_im = imread(the_images[i])
                    res_im = SpatialImage(the_im.copy())[-1::-1, :, :]
                    res_im._set_resolution(the_im._get_resolution())
                    imsave(res_images[i], res_im)
                    del the_im
                    del res_im
                else:
                    monitoring.to_log_and_console("       already existing", 2)

    #
    # 1. Putting all images in a common reference
    # - resampling of first image in an isotropic grid = reference image
    # - co-registration of other images
    # 2. Compute weights with an ad-hoc method
    #

    the_image_list = res_image_list[:]
    res_image_list = list()
    init_trsfs = []
    prereg_trsfs = []
    res_trsfs = []
    unreg_weight_images = []
    weight_images = []

    #
    # build the file names
    #

    for c in range(0, len(channel)):

        the_images = the_image_list[c]
        res_images = []

        for i in range(0, len(the_images)):
            res_images.append(commonTools.add_suffix(input_image_list[c][i], "_reg",
                                                     new_dirname=channel[c].temporary_paths[i],
                                                     new_extension=parameters.default_image_suffix))
            if c == 0:
                init_trsfs.append(
                    commonTools.add_suffix(input_image_list[c][i], "_init",
                                           new_dirname=channel[c].temporary_paths[i], new_extension="trsf"))
                prereg_trsfs.append(commonTools.add_suffix(input_image_list[c][i], "_prereg",
                                                           new_dirname=channel[c].temporary_paths[i],
                                                           new_extension="trsf"))
                res_trsfs.append(commonTools.add_suffix(input_image_list[c][i], "_reg",
                                                        new_dirname=channel[c].temporary_paths[i],
                                                        new_extension="trsf"))
                unreg_weight_images.append(commonTools.add_suffix(input_image_list[c][i], "_init_weight",
                                                                  new_dirname=channel[c].temporary_paths[i],
                                                                  new_extension=parameters.default_image_suffix))
                weight_images.append(commonTools.add_suffix(input_image_list[c][i], "_weight",
                                                            new_dirname=channel[c].temporary_paths[i],
                                                            new_extension=parameters.default_image_suffix))

        res_image_list.append(res_images)

    #
    # is there something to do ?
    # check whether fused images are missing
    # if one image is missing, re-process channel #0 to get weights and sum of weights
    #

    do_something = [False] * len(input_image_list[0])

    for c in range(0, len(channel)):

        the_images = the_image_list[c]
        res_images = res_image_list[c]

        if not os.path.isfile(os.path.join(channel[c].path_fuse_exp, fused_image)):
            do_something = [True] * len(input_image_list[0])
            break

        for i in range(0, len(the_images)):
            if os.path.isfile(res_images[i]):
                if monitoring.forceResultsToBeBuilt is True:
                    do_something[i] = True
            else:
                do_something[i] = True

    for i in range(1, len(input_image_list[0])):
        if do_something[i] is True:
            do_something[0] = True

    #
    # default angle for initial rotation matrix
    #

    if parameters.acquisition_orientation.lower() == 'left':
        default_angle = 270.0
    elif parameters.acquisition_orientation.lower() == 'right':
        default_angle = 90.0
    else:
        monitoring.to_log_and_console(proc + ": unknown acquisition orientation '"
                                      + str(parameters.acquisition_orientation) + "'", 0)
        monitoring.to_log_and_console("Exiting.", 0)
        sys.exit(1)

    #
    # process
    # transformations and weights are computed on channel #0
    # and are used for other channels
    # - first image is just resampled to the destination resolution
    # - other images are co-registered with the first image
    #

    fusion_box = None

    for c in range(0, len(channel)):

        the_images = the_image_list[c]
        res_images = res_image_list[c]
        full_image = None

        for i in range(0, len(the_images)):

            monitoring.to_log_and_console("    .. process '"
                                          + the_images[i].split(os.path.sep)[-1] + "' for fusion", 2)

            if do_something[i] is False:
                monitoring.to_log_and_console("       nothing to do", 2)
                continue

            #
            # first image: resampling only
            #
            if i == 0:
                #
                # image center
                #
                im = imread(the_images[i])
                ref_center = np.multiply(im.shape[:3], im.resolution) / 2.0
                del im

                #
                # resampling first image
                #
                monitoring.to_log_and_console("       resampling '" + the_images[i].split(os.path.sep)[-1]
                                              + "' at " + str(parameters.target_resolution), 2)
                if not os.path.isfile(res_images[i]) or monitoring.forceResultsToBeBuilt is True:
                    cpp_wrapping.apply_transformation(the_images[i], res_images[i], the_transformation=None,
                                                      template_image=None,
                                                      voxel_size=parameters.target_resolution,
                                                      interpolation_mode='linear',
                                                      monitoring=monitoring)
                else:
                    monitoring.to_log_and_console("       already existing", 2)

            #
            # other images:
            # - channel #0: co-registration
            # - other channels: resampling with transformation of channel #0
            #
            else:

                if c == 0:

                    monitoring.to_log_and_console("       co-registering '" + the_images[i].split(os.path.sep)[-1]
                                                  + "'", 2)

                    #
                    # initial transformation
                    #
                    if not os.path.isfile(init_trsfs[i]) or monitoring.forceResultsToBeBuilt is True:
                        if i == 1:
                            angle = 0.0
                        else:
                            angle = default_angle
                        monitoring.to_log_and_console("       angle used for '" + init_trsfs[i].split(os.path.sep)[-1]
                                                      + "' is " + str(angle), 2)
                        im = imread(the_images[i])
                        flo_center = np.multiply(im.shape[:3], im.resolution) / 2.0
                        del im

                        #
                        # the initial transformation was computed with _axis_rotation_matrix(). To compute the
                        # translation, it preserves the center of the field of view of the floating image.
                        # However it seems more coherent to compute a translation that put the center of FOV of the
                        # floating image onto he FOV of the reference one.
                        #
                        # the call to _axis_rotation_matrix() was
                        # rotation_matrix = _axis_rotation_matrix(axis="Y", angle=angle, min_space=(0, 0, 0),
                        #                                         max_space=np.multiply(im.shape[:3], im.resolution))
                        # Note: it requires that 'im' is deleted after the call
                        #
                        # it can be mimicked by
                        # _ init_rotation_matrix(axis="Y", angle=angle, ref_center=flo_center, flo_center=flo_center)
                        #

                        rotation_matrix = _init_rotation_matrix(axis="Y", angle=angle, ref_center=ref_center,
                                                                flo_center=flo_center)

                        np.savetxt(init_trsfs[i], rotation_matrix)
                        del rotation_matrix

                    #
                    # a two-fold registration, translation then affine, could be preferable
                    #
                    if not os.path.isfile(res_images[i]) or not os.path.isfile(res_trsfs[i]) \
                            or monitoring.forceResultsToBeBuilt is True:
                        if parameters.pre_registration.compute_registration is True:
                            _linear_registration(res_images[0], the_images[i], res_images[i],
                                                 prereg_trsfs[i], init_trsfs[i], parameters.pre_registration)
                            _linear_registration(res_images[0], the_images[i], res_images[i],
                                                 res_trsfs[i], prereg_trsfs[i], parameters.registration)
                        else:
                            _linear_registration(res_images[0], the_images[i], res_images[i],
                                                 res_trsfs[i], init_trsfs[i], parameters.registration)
                    else:
                        monitoring.to_log_and_console("       already existing", 2)
                #
                # other channels
                #
                else:

                    monitoring.to_log_and_console("       resampling '" + the_images[i].split(os.path.sep)[-1] + "'", 2)
                    if not os.path.isfile(res_images[i]) or monitoring.forceResultsToBeBuilt is True:
                        cpp_wrapping.apply_transformation(the_images[i], res_images[i],
                                                          the_transformation=res_trsfs[i], template_image=res_images[0],
                                                          voxel_size=None, interpolation_mode='linear',
                                                          monitoring=monitoring)
                    else:
                        monitoring.to_log_and_console("       already existing", 2)

            #
            # compute weighting masks on channel #0
            # - mask is computed on an untransformed image
            #   however, resolution may have changes, or it can be cropped
            #   or it can be mirrored
            # - mask are then transformed with the computed transformation
            #
            if c == 0:

                monitoring.to_log_and_console("       computing weights for fusion", 2)

                if i % 2 == 1:
                    direction = False
                else:
                    direction = True

                if not os.path.isfile(unreg_weight_images[i]) or monitoring.forceResultsToBeBuilt is True:
                    im = imread(the_images[i])
                    unreg_weight = _build_mask(im, direction)
                    unreg_weight._set_resolution(im._get_resolution())
                    imsave(unreg_weight_images[i], unreg_weight)
                    del im
                    del unreg_weight
                else:
                    monitoring.to_log_and_console("       already existing", 2)

                monitoring.to_log_and_console("       resampling '" + unreg_weight_images[i].split(os.path.sep)[-1]
                                              + "'", 2)
                if i == 0:
                    if not os.path.isfile(weight_images[i]) or monitoring.forceResultsToBeBuilt is True:
                        cpp_wrapping.apply_transformation(unreg_weight_images[i], weight_images[i],
                                                          the_transformation=None, template_image=None,
                                                          voxel_size=parameters.target_resolution,
                                                          interpolation_mode='linear', monitoring=monitoring)
                    else:
                        monitoring.to_log_and_console("       already existing", 2)
                else:
                    if not os.path.isfile(weight_images[i]) or monitoring.forceResultsToBeBuilt is True:
                        cpp_wrapping.apply_transformation(unreg_weight_images[i], weight_images[i],
                                                          the_transformation=res_trsfs[i], template_image=res_images[0],
                                                          voxel_size=None, interpolation_mode='linear',
                                                          monitoring=monitoring)
                    else:
                        monitoring.to_log_and_console("       already existing", 2)

        #
        # compute fused image as a linear combination of co-registered images
        # the sun of weights have been precomputed to mimic historical behavior
        #
        # do not forget to cast the result on 16 bits
        #

        monitoring.to_log_and_console("    .. combining images", 2)

        if parameters.fusion_cropping is True:
            tmp_fused_image = commonTools.add_suffix(fused_image, "_uncropped_fusion",
                                                     new_dirname=channel[c].temporary_paths[4],
                                                     new_extension=parameters.default_image_suffix)
        else:
            tmp_fused_image = os.path.join(channel[c].path_fuse_exp, fused_image)

        cpp_wrapping.linear_combination(weight_images, res_images, tmp_fused_image, monitoring=monitoring)

        #
        # save image if fusion is required
        #
        if parameters.fusion_cropping is True:

            if c == 0:
                fusion_box = _crop_bounding_box(tmp_fused_image)

            monitoring.to_log_and_console("    .. cropping '" + fused_image.split(os.path.sep)[-1], 2)
            _crop_disk_image(tmp_fused_image, os.path.join(channel[c].path_fuse_exp, fused_image), fusion_box,
                             parameters.fusion_cropping_margin_x_0,
                             parameters.fusion_cropping_margin_x_1,
                             parameters.fusion_cropping_margin_y_0,
                             parameters.fusion_cropping_margin_y_1)

    return


#
#
# read the raw data
#
#


def fusion_preprocess(input_images, fused_image, time_point, environment, parameters):
    """

    :param input_images:
    :param fused_image:
    :param time_point:
    :param environment:
    :param parameters:
    :return:
    """

    proc = "fusion_preprocess"

    if monitoring.debug > 1:

        print proc + " was called with:"
        print "- input_images = " + str(input_images)
        print "- fused_image = " + str(fused_image)
        print "- time_point = " + str(time_point)
        print ""

    monitoring.to_log_and_console('... fusion of time ' + time_point, 1)
    if len(environment.channel) > 1:
        monitoring.to_log_and_console('    there are ' + str(len(environment.channel)) + ' channels to be fused', 1)

    #
    # check whether there exists some unfused channel
    #

    do_something = False
    for c in range(0, len(environment.channel)):
        if os.path.isfile(os.path.join(environment.channel[c].path_fuse_exp, fused_image)):
            if not monitoring.forceResultsToBeBuilt:
                monitoring.to_log_and_console('    channel #' + str(c) + ' already existing', 2)
            else:
                monitoring.to_log_and_console('    channel #' + str(c) + ' already existing, but forced', 2)
                do_something = True
        else:
            do_something = True

    if do_something is False:
        monitoring.to_log_and_console('    nothing to do', 2)
        return

    #
    # start processing
    #

    start_time = time.time()

    #
    # directory for auxiliary files
    #
    # ANGLE_0: LC/Stack0000
    # ANGLE_1: RC/Stack0000
    # ANGLE_2: LC/Stack0001
    # ANGLE_3: RC/Stack0001
    #

    for c in range(0, len(environment.channel)):
        environment.channel[c].temporary_paths = list()
        environment.channel[c].temporary_paths.append(os.path.join(environment.channel[c].path_fuse_exp,
                                                                   "TEMP_$TIME", "ANGLE_0"))
        environment.channel[c].temporary_paths.append(os.path.join(environment.channel[c].path_fuse_exp,
                                                                   "TEMP_$TIME", "ANGLE_1"))
        environment.channel[c].temporary_paths.append(os.path.join(environment.channel[c].path_fuse_exp,
                                                                   "TEMP_$TIME", "ANGLE_2"))
        environment.channel[c].temporary_paths.append(os.path.join(environment.channel[c].path_fuse_exp,
                                                                   "TEMP_$TIME", "ANGLE_3"))
        environment.channel[c].temporary_paths.append(os.path.join(environment.channel[c].path_fuse_exp,
                                                                   "TEMP_$TIME"))

    #
    # recall that time_point is a string here
    # nomenclature.replaceTIME() can not be used
    #
    for c in range(0, len(environment.channel)):
        for i in range(0, len(environment.channel[c].temporary_paths)):
            environment.channel[c].temporary_paths[i] = \
                environment.channel[c].temporary_paths[i].replace(nomenclature.FLAG_TIME, time_point)
            if not os.path.isdir(environment.channel[c].temporary_paths[i]):
                os.makedirs(environment.channel[c].temporary_paths[i])

    #
    # get image file names
    # - may involve unzipping and conversion
    #
    monitoring.to_log_and_console('    get original images', 2)

    image_list = list()

    for c in range(0, len(environment.channel)):
        images = list()
        images.append(_read_image_name(environment.channel[c].path_angle1,
                                       environment.channel[c].temporary_paths[0],
                                       input_images[0],
                                       parameters.acquisition_resolution, parameters.default_image_suffix))
        images.append(_read_image_name(environment.channel[c].path_angle2,
                                       environment.channel[c].temporary_paths[1],
                                       input_images[1],
                                       parameters.acquisition_resolution, parameters.default_image_suffix))
        images.append(_read_image_name(environment.channel[c].path_angle3,
                                       environment.channel[c].temporary_paths[2],
                                       input_images[2],
                                       parameters.acquisition_resolution, parameters.default_image_suffix))
        images.append(_read_image_name(environment.channel[c].path_angle4,
                                       environment.channel[c].temporary_paths[3],
                                       input_images[3],
                                       parameters.acquisition_resolution, parameters.default_image_suffix))
        image_list.append(images)

    #
    #
    #
    monitoring.to_log_and_console('    fuse images', 2)

    fusion_process(image_list, fused_image, environment.channel, parameters)

    #
    # remove temporary files if required
    #

    if monitoring.keepTemporaryFiles is False:
        for c in range(0, len(environment.channel)):
            shutil.rmtree(environment.channel[c].temporary_paths[4])

    #
    # end processing
    #

    end_time = time.time()
    monitoring.to_log_and_console('    computation time = ' + str(end_time - start_time) + ' s', 1)
    monitoring.to_log_and_console('', 1)

    return


#
#
# Parse the raw data directories and identify data to be fused for each time point
# - for each time point fusion_preprocess() is then called
#
#


def fusion_control(experiment, environment, parameters):
    """

    :param experiment:
    :param environment:
    :param parameters:
    :return:
    """

    proc = 'fusion_control'
    default_width = 3

    #
    # make sure that the result directory exists
    # (although it should have been verified in 1-fuse.py)
    #

    for c in range(0, len(environment.channel)):
        if not os.path.isdir(environment.channel[c].path_fuse_exp):
            os.makedirs(environment.channel[c].path_fuse_exp)

    if not os.path.isdir(environment.path_logdir):
        os.makedirs(environment.path_logdir)

    monitoring.to_log_and_console('', 1)

    #
    # if data directories of the main channel are different, parse them
    #

    if environment.channel[0].path_angle1 != environment.channel[0].path_angle2 \
            and environment.channel[0].path_angle1 != environment.channel[0].path_angle3 \
            and environment.channel[0].path_angle1 != environment.channel[0].path_angle4 \
            and environment.channel[0].path_angle2 != environment.channel[0].path_angle3 \
            and environment.channel[0].path_angle2 != environment.channel[0].path_angle4 \
            and environment.channel[0].path_angle3 != environment.channel[0].path_angle4:

        prefix1, time_length1, time_points1, suffix1 = _analyze_data_directory(environment.channel[0].path_angle1)
        prefix2, time_length2, time_points2, suffix2 = _analyze_data_directory(environment.channel[0].path_angle2)
        prefix3, time_length3, time_points3, suffix3 = _analyze_data_directory(environment.channel[0].path_angle3)
        prefix4, time_length4, time_points4, suffix4 = _analyze_data_directory(environment.channel[0].path_angle4)

        if monitoring.debug > 0:
            print ""
            print "analysis of '" + str(environment.channel[0].path_angle1) + "'"
            print "   -> " + prefix1
            print "   -> " + str(time_length1)
            print "   -> " + str(time_points1)
            print "   -> " + suffix1
            print "analysis of '" + str(environment.channel[0].path_angle2) + "'"
            print "   -> " + prefix2
            print "   -> " + str(time_length2)
            print "   -> " + str(time_points2)
            print "   -> " + suffix2
            print "analysis of '" + str(environment.channel[0].path_angle3) + "'"
            print "   -> " + prefix3
            print "   -> " + str(time_length3)
            print "   -> " + str(time_points3)
            print "   -> " + suffix3
            print "analysis of '" + str(environment.channel[0].path_angle4) + "'"
            print "   -> " + prefix4
            print "   -> " + str(time_length4)
            print "   -> " + str(time_points4)
            print "   -> " + suffix4
            print ""

        #
        # loop over acquisitions
        # 1. case where all acquisition have to be processed
        #    begin < 0 or end < 0 or begin > end or delta < 0
        # 2. only a few acquisitions have to be processed
        #

        extra_zeros = ''
        if time_length1 < default_width:
            extra_zeros = (default_width - time_length1) * '0'

        if experiment.first_time_point < 0 or experiment.last_time_point < 0 or experiment.delta_time_point < 0 \
                or experiment.first_time_point > experiment.last_time_point:

            for time_point in time_points1:

                #
                # fused image name
                #

                fused_image = environment.path_fuse_exp_files.replace(nomenclature.FLAG_TIME, extra_zeros + time_point)\
                              + '.' + parameters.result_image_suffix

                #
                # input image names
                #

                images = list()

                images.append(prefix1 + time_point + suffix1)
                im = prefix2 + time_point + suffix2
                if time_point not in time_points2:
                    print proc + ": image '" + im + "' not found in '" + environment.path_angle2 + "'"
                    print "\t Exiting."
                    sys.exit(1)
                else:
                    images.append(im)
                im = prefix3 + time_point + suffix3
                if time_point not in time_points3:
                    print proc + ": image '" + im + "' not found in '" + environment.path_angle3 + "'"
                    print "\t Exiting."
                    sys.exit(1)
                else:
                    images.append(im)
                im = prefix4 + time_point + suffix4
                if time_point not in time_points4:
                    print proc + ": image '" + im + "' not found in '" + environment.path_angle4 + "'"
                    print "\t Exiting."
                    sys.exit(1)
                else:
                    images.append(im)

                #
                # process
                #

                fusion_preprocess(images, fused_image, extra_zeros + time_point, environment, parameters)

        else:

            if experiment.first_time_point < 0 or experiment.last_time_point < 0:
                monitoring.to_log_and_console("... time interval does not seem to be defined in the parameter file")
                monitoring.to_log_and_console("    set parameters 'begin' and 'end'")
                monitoring.to_log_and_console("\t Exiting")
                sys.exit(1)

            for time_value in range(experiment.first_time_point, experiment.last_time_point + 1,
                                    experiment.delta_time_point):

                acquisition_time = str('{:0{width}d}'.format(time_value, width=time_length1))
                fused_time = str('{:0{width}d}'.format(time_value + experiment.delay_time_point, width=time_length1))

                #
                # fused image name
                #

                fused_image = environment.path_fuse_exp_files.replace(nomenclature.FLAG_TIME, extra_zeros+fused_time) \
                              + '.' + parameters.result_image_suffix

                #
                # input image names
                #

                images = list()

                im = prefix1 + acquisition_time + suffix1
                if acquisition_time not in time_points1:
                    monitoring.to_log_and_console("    .. image '" + im + "' not found in '"
                                                  + environment.channel[0].path_angle1 + "'", 2)
                    monitoring.to_log_and_console("       skip time " + str(acquisition_time), 2)
                    continue
                else:
                    images.append(im)
                im = prefix2 + acquisition_time + suffix2
                if acquisition_time not in time_points2:
                    monitoring.to_log_and_console("    .. image '" + im + "' not found in '"
                                                  + environment.channel[0].path_angle2 + "'", 2)
                    monitoring.to_log_and_console("       skip time " + str(acquisition_time), 2)
                    continue
                else:
                    images.append(im)
                im = prefix3 + acquisition_time + suffix3
                if acquisition_time not in time_points3:
                    monitoring.to_log_and_console("    .. image '" + im + "' not found in '"
                                                  + environment.channel[0].path_angle3 + "'", 2)
                    monitoring.to_log_and_console("       skip time " + str(acquisition_time), 2)
                    continue
                else:
                    images.append(im)
                im = prefix4 + acquisition_time + suffix4
                if acquisition_time not in time_points4:
                    monitoring.to_log_and_console("    .. image '" + im + "' not found in '"
                                                  + environment.channel[0].path_angle4 + "'", 2)
                    monitoring.to_log_and_console("       skip time " + str(acquisition_time), 2)
                    continue
                else:
                    images.append(im)

                #
                # process
                #

                fusion_preprocess(images, fused_image, extra_zeros + acquisition_time, environment, parameters)

    #
    # here data directories are not different, we have to rely on built names
    #

    else:

        if experiment.first_time_point < 0 or experiment.last_time_point < 0:
            monitoring.to_log_and_console("... time interval does not seem to be defined in the parameter file")
            monitoring.to_log_and_console("    set parameters 'begin' and 'end'")
            monitoring.to_log_and_console("\t Exiting")
            sys.exit(1)

        for time_value in range(experiment.first_time_point, experiment.last_time_point+1, experiment.delta_time_point):

            acquisition_time = str('{:0{width}d}'.format(time_value, width=default_width))

            #
            # fused image name
            #

            fused_image = nomenclature.replaceTIME(environment.path_fuse_exp_files,
                                                   time_value+experiment.delay_time_point) \
                          + '.' + parameters.result_image_suffix

            #
            # input image names
            #

            images = list()

            name = nomenclature.replaceTIME(environment.path_angle1_files, time_value)
            im = commonTools.find_file(environment.channel[0].path_angle1, name, monitoring)
            if im is None:
                monitoring.to_log_and_console("    .. image '" + name + "' not found in '"
                                              + environment.channel[0].path_angle1 + "'", 2)
                monitoring.to_log_and_console("       skip time " + str(acquisition_time), 2)
                continue
            else:
                images.append(im)

            name = nomenclature.replaceTIME(environment.path_angle2_files, time_value)
            im = commonTools.find_file(environment.channel[0].path_angle2, name, monitoring)
            if im is None:
                monitoring.to_log_and_console("    .. image '" + name + "' not found in '"
                                              + environment.channel[0].path_angle2 + "'", 2)
                monitoring.to_log_and_console("       skip time " + str(acquisition_time), 2)
                continue
            else:
                images.append(im)

            name = nomenclature.replaceTIME(environment.path_angle3_files, time_value)
            im = commonTools.find_file(environment.channel[0].path_angle3, name, monitoring)
            if im is None:
                monitoring.to_log_and_console("    .. image '" + name + "' not found in '"
                                              + environment.channel[0].path_angle3 + "'", 2)
                monitoring.to_log_and_console("       skip time " + str(acquisition_time), 2)
                continue
            else:
                images.append(im)

            name = nomenclature.replaceTIME(environment.path_angle4_files, time_value)
            im = commonTools.find_file(environment.channel[0].path_angle4, name, monitoring)
            if im is None:
                monitoring.to_log_and_console("    .. image '" + name + "' not found in '"
                                              + environment.channel[0].path_angle4 + "'", 2)
                monitoring.to_log_and_console("       skip time " + str(acquisition_time), 2)
                continue
            else:
                images.append(im)

            #
            # process
            #

            fusion_preprocess(images, fused_image, acquisition_time, environment, parameters)

    return
