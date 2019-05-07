
import os
import imp
import sys
import time
import numpy as np
from scipy import ndimage as nd

import commonTools
import nomenclature

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


class ManualCorrectionEnvironment(object):

    def __init__(self):

        #
        # mars data paths
        #
        self.path_mars_exp_files = None

        #
        # segmentation data paths
        #
        self.path_seg_exp = None
        self.path_seg_exp_files = None

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

        self.path_mars_exp_files = nomenclature.replaceFlags(nomenclature.path_mars_exp_files, parameters)

        self.path_seg_exp = nomenclature.replaceFlags(nomenclature.path_seg_exp, parameters)
        self.path_seg_exp_files = nomenclature.replaceFlags(nomenclature.path_seg_exp_files, parameters)

        self.path_logdir = nomenclature.replaceFlags(nomenclature.path_seg_logdir, parameters)
        self.path_history_file = nomenclature.replaceFlags(nomenclature.path_seg_historyfile, parameters)
        self.path_log_file = nomenclature.replaceFlags(nomenclature.path_seg_logfile, parameters, start_time)

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('ManualCorrectionEnvironment\n')

            logfile.write('- path_mars_exp_files = ' + str(self.path_mars_exp_files) + '\n')

            logfile.write('- path_seg_exp = ' + str(self.path_seg_exp) + '\n')
            logfile.write('- path_seg_exp_files = ' + str(self.path_seg_exp_files) + '\n')

            logfile.write('- path_logdir = ' + str(self.path_logdir) + '\n')
            logfile.write('- path_history_file = ' + str(self.path_history_file)+'\n')
            logfile.write('- path_log_file = ' + str(self.path_log_file)+'\n')
            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('ManualCorrectionEnvironment')

        print('- path_mars_exp_files = ' + str(self.path_mars_exp_files))

        print('- path_seg_exp = ' + str(self.path_seg_exp))
        print('- path_seg_exp_files = ' + str(self.path_seg_exp_files))

        print('- path_logdir = ' + str(self.path_logdir))
        print('- path_history_file = ' + str(self.path_history_file))
        print('- path_log_file = ' + str(self.path_log_file))
        print("")


#
#
#
#
#


class ManualCorrectionParameters(object):

    def __init__(self):

        #
        #
        #
        self.first_time_point = -1
        self.last_time_point = -1

        #
        #
        #
        self.input_image = None
        self.output_image = None
        self.mapping_file = None

        #
        #
        #
        self.smallest_cells = 8
        self.largest_cells = 8

        #
        # images suffixes/formats
        #
        self.result_image_suffix = 'inr'
        self.default_image_suffix = 'inr'

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('ManualCorrectionParameters\n')

            logfile.write('- first_time_point = ' + str(self.first_time_point) + '\n')
            logfile.write('- last_time_point = ' + str(self.last_time_point) + '\n')

            logfile.write('- input_image = ' + str(self.input_image) + '\n')
            logfile.write('- output_image = ' + str(self.output_image) + '\n')
            logfile.write('- mapping_file = ' + str(self.mapping_file) + '\n')

            logfile.write('- result_image_suffix = ' + str(self.result_image_suffix) + '\n')
            logfile.write('- default_image_suffix = '+str(self.default_image_suffix) + '\n')

            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('ManualCorrectionParameters')

        print('- first_time_point = ' + str(self.first_time_point))
        print('- last_time_point = ' + str(self.last_time_point))

        print('- input_image = ' + str(self.input_image))
        print('- output_image = ' + str(self.output_image))
        print('- mapping_file = ' + str(self.mapping_file))

        print('- result_image_suffix = ' + str(self.result_image_suffix))
        print('- default_image_suffix = ' + str(self.default_image_suffix))

        print("")

    def update_from_args(self, args):
        self.input_image = args.input_image
        self.output_image = args.output_image
        self.mapping_file = args.mapping_file
        if int(args.smallest_cells) >= 0:
            self.smallest_cells = args.smallest_cells
        if int(args.largest_cells) >= 0:
            self.largest_cells = args.largest_cells

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
        #
        #
        if hasattr(parameters, 'mancor_input_seg_file'):
            if parameters.mancor_input_seg_file is not None and len(str(parameters.mancor_input_seg_file)) > 0:
                self.input_image = parameters.mancor_input_seg_file
        #
        # for back-compatibility
        #
        if hasattr(parameters, 'mancor_seg_file'):
            if parameters.mancor_seg_file is not None and len(str(parameters.mancor_seg_file)) > 0:
                self.input_image = parameters.mancor_seg_file

        if hasattr(parameters, 'mancor_output_seg_file'):
            if parameters.mancor_output_seg_file is not None and len(str(parameters.mancor_output_seg_file)) > 0:
                self.output_image = parameters.mancor_output_seg_file
        if hasattr(parameters, 'mancor_mapping_file'):
            if parameters.mancor_mapping_file is not None and len(str(parameters.mancor_mapping_file)) > 0:
                self.mapping_file = parameters.mancor_mapping_file

        #
        # images suffixes/formats
        #
        if hasattr(parameters, 'result_image_suffix'):
            if parameters.result_image_suffix is not None:
                self.result_image_suffix = parameters.result_image_suffix
        if hasattr(parameters, 'default_image_suffix'):
            if parameters.default_image_suffix is not None:
                self.default_image_suffix = parameters.default_image_suffix


########################################################################################
#
# some internal procedures
#
########################################################################################


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

def correction_process(input_image, output_image, parameters):
    """

    :param input_image:
    :param output_image:
    :param parameters:
    :return:
    """

    proc = "correction_process"

    if monitoring.debug > 2:
        print ""
        print proc + " was called with:"
        print "- input_image = " + str(input_image)
        print "- output_image = " + str(output_image)
        print ""

    #
    # nothing to do if the corrected image exists
    #

    if os.path.isfile(output_image):
        if monitoring.forceResultsToBeBuilt is False:
            monitoring.to_log_and_console('    corrected image already existing', 2)
            return
        else:
            monitoring.to_log_and_console('    corrected image already existing, but forced', 2)

    #
    #
    #

    im = imread(input_image)
    voxelsize = im._get_resolution()
    vol = voxelsize[0] * voxelsize[1] * voxelsize[2]

    #
    #
    #
    if parameters.mapping_file is not None and len(str(parameters.mapping_file)) > 0 \
        and os.path.isfile(parameters.mapping_file):

        mapping = np.arange(np.max(im)+1)
        # print str(mapping)

        #
        # corrections to be done
        #
        corrected_mapping = commonTools.read_lut(parameters.mapping_file)

        if len(corrected_mapping) > 0:
            for k, v in corrected_mapping.iteritems():
                mapping[k] = v
            im = mapping[im]
        else:
            monitoring.to_log_and_console('    no corrections to be done (no valid mapping file)', 2)

    else:
        monitoring.to_log_and_console('    no corrections to be done (no valid mapping file)', 2)

    #
    # get
    # - the list of cell label
    # - the list of cell volume
    #
    # build a dictionary and sort it (increasing order) wrt the volume

    cell_label = np.unique(im)
    cell_volume = nd.sum(np.ones_like(im), im, index=np.int16(cell_label))

    d = dict(zip(cell_label, cell_volume))
    s = sorted(d, key=d.__getitem__)

    monitoring.to_log_and_console('    Number of cells: ' + str(len(cell_label)), 0)
    monitoring.to_log_and_console('    Maximal label: ' + str(np.max(im)), 0)
    # monitoring.to_log_and_console('    Cell ids: ' + str(cell_label), 0)

    monitoring.to_log_and_console('    Sorted cell volumes: ', 0)
    monitoring.to_log_and_console('      Id :    voxels          (um^3)', 0)

    if (int(parameters.smallest_cells) <= 0 and int(parameters.largest_cells) <= 0) \
            or (int(parameters.smallest_cells) + int(parameters.largest_cells) >= len(s)):
        for l in s:
            monitoring.to_log_and_console('    {:>4d} : {:>9d} {:>15s}'.format(l, int(d[l]),
                                                                               '({:.2f})'.format(d[l]*vol)), 0)
    else:
        if int(parameters.smallest_cells) > 0:
            for l in s[:int(parameters.smallest_cells)]:
                monitoring.to_log_and_console('    {:>4d} : {:>9d} {:>15s}'.format(l, int(d[l]),
                                                                                   '({:.2f})'.format(d[l]*vol)), 0)
        monitoring.to_log_and_console('       ...', 0)
        if int(parameters.largest_cells) > 0:
            for l in s[-int(parameters.largest_cells):]:
                monitoring.to_log_and_console('    {:>4d} : {:>9d} {:>15s}'.format(l, int(d[l]),
                                                                                   '({:.2f})'.format(d[l]*vol)), 0)

    imsave(output_image, SpatialImage(im, voxelsize=voxelsize).astype(np.uint16))

    return

#
#
#
#
#


def correction_control(experiment, environment, parameters):
    """

    :param experiment:
    :param environment:
    :param parameters:
    :return:
    """

    monitoring.to_log_and_console('', 1)

    #
    # case #1:
    # - input and/or output images are given by the user
    # case #2
    # - input image is the mars image and output image is the corresponding segmentation image
    #

    if (parameters.input_image is not None and len(str(parameters.input_image)) > 0) \
        or (parameters.output_image is not None and len(str(parameters.output_image)) > 0):

        #
        # input image
        #
        if parameters.input_image is not None and len(str(parameters.input_image)) > 0:
            if not os.path.isfile(parameters.input_image):
                monitoring.to_log_and_console("       '" + str(parameters.input_image).split(os.path.sep)[-1]
                                              + "' does not exist", 2)
                monitoring.to_log_and_console("\t Exiting.")
                sys.exit(1)
            input_image = parameters.input_image
        else:
            input_name = nomenclature.replaceTIME(environment.path_mars_exp_files,
                                                  experiment.first_time_point + experiment.delay_time_point)
            input_name = commonTools.find_file(environment.path_mars_exp, input_name, monitoring=monitoring)
            input_image = os.path.join(environment.path_mars_exp, input_name)

        #
        # output image
        #
        if parameters.output_image is not None and len(str(parameters.output_image)) > 0:
            output_image = parameters.output_image
        else:
            output_name = nomenclature.replaceTIME(environment.path_seg_exp_files,
                                                   experiment.first_time_point + experiment.delay_time_point) \
                          + '.' + parameters.result_image_suffix
            output_image = os.path.join(environment.path_seg_exp, output_name)

        #
        # start processing
        #
        monitoring.to_log_and_console("... correction of '" + str(input_image).split(os.path.sep)[-1] + "'", 1)
        start_time = time.time()

        correction_process(input_image, output_image, parameters)

        #
        # end processing for a time point
        #
        end_time = time.time()
        monitoring.to_log_and_console('    computation time = ' + str(end_time - start_time) + ' s', 1)
        monitoring.to_log_and_console('', 1)

    else:

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

            input_mars_name = nomenclature.replaceTIME(environment.path_mars_exp_files, time_value)
            input_name = commonTools.find_file(environment.path_seg_exp, input_mars_name, monitoring=monitoring)
            if input_name is None:
                monitoring.to_log_and_console("    mars image '" + str(input_mars_name) + "' not found: skip time "
                                              + str(time_value), 1)
                continue

            input_image = os.path.join(environment.path_seg_exp, input_name)
            output_name = nomenclature.replaceTIME(environment.path_seg_exp_files, time_value) \
                          + '.' + parameters.result_image_suffix
            output_image = os.path.join(environment.path_seg_exp, output_name)

            #
            # start processing
            #
            monitoring.to_log_and_console("... correction of '" + str(input_image).split(os.path.sep)[-1] + "'", 1)
            start_time = time.time()

            correction_process(input_image, output_image, parameters)

            #
            # end processing for a time point
            #
            end_time = time.time()
            monitoring.to_log_and_console('    computation time = ' + str(end_time - start_time) + ' s', 1)
            monitoring.to_log_and_console('', 1)
