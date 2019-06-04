
import os
import sys
import subprocess

from ImageHandling import imread, imsave, SpatialImage

#
# path_to_bins = '/user/gmicheli/home/DIG-EM/Codes/Packages/ASTEC-170210/ASTEC/CommunFunctions/cpp/build/bin/'
#
#
path_to_bins = os.path.join(os.path.dirname(__file__), 'cpp')+os.path.sep
#
# path_filters = path_to_bins + 'recfilters'
# path_linearfilters = path_to_bins + 'linearFilter'
# path_reech3d = path_to_bins + 'reech3d'
# path_apply_trsf = path_to_bins + "applyTrsf"
# path_block = path_to_bins + 'blockmatching'
# path_morpho = path_to_bins + "morpho"
# path_regional_max = path_to_bins + 'regionalmax'
# path_regional_ext = path_to_bins + 'regionalext'
# path_connexe = path_to_bins + 'connexe'
# path_watershed = path_to_bins + 'watershed'
# path_gradient_norm = path_to_bins + 'norme_gradient'
#
#
# path_membrane = path_to_bins + 'membrane'
# path_anisotropicHist = path_to_bins + 'anisotropicHist'
# path_TVmembrane = path_to_bins + 'TVmembrane'
# path_seuillage = path_to_bins + 'seuillage'
# path_Arit = path_to_bins + 'Arit'
# path_Logic = path_to_bins + 'Logic'
# path_create_image = path_to_bins + 'createImage'
# path_linearFilter = path_to_bins + 'linearFilter'
# path_copy = path_to_bins + 'copy'
# path_symmetryPlane = path_to_bins + 'symmetryPlane'
# path_dice = path_to_bins + 'dice'
# path_diceMaximisation = path_to_bins + 'diceMaximisation'
# path_directionHistogram = path_to_bins + 'directionHistogram'
# path_directionHistogramMaxima = path_to_bins + 'directionHistogramMaxima'
# path_planeRegistration = path_to_bins + 'planeRegistration'
# path_pointCloudRegistration = path_to_bins + 'pointCloudRegistration'
# path_setvoxelsize = path_to_bins + 'setVoxelSize'
# path_compose_trsf = path_to_bins + 'composeTrsf'
# path_multiple_trsfs = path_to_bins + 'multipleTrsfs'
# path_change_multiple_trsfs = path_to_bins + 'changeMultipleTrsfs'
# path_non_zeros_image = path_to_bins + 'nonZerosImage'
# path_setvoxelvalue = path_to_bins + 'setVoxelValue'
# path_fuselabels = path_to_bins + "fuseLabels"
# path_labelborders = path_to_bins + "labelBorders"
# path_associateLabels = path_to_bins + "associateLabels"
#
# path_boundingboxes = path_to_bins + 'boundingboxes'
# path_cropImage = path_to_bins + 'cropImage'
# path_patchLogic = path_to_bins + 'patchLogic'
# path_mc_adhocfuse = path_to_bins + 'mc-adhocFuse'
#

############################################################
#
#
#
############################################################


def _write_error_msg(text, monitoring):
    """
    Write an error message
    :param text:
    :param monitoring:
    :return:
    """
    if monitoring is not None:
        monitoring.to_log_and_console(text, 0)
    else:
        print(text)


def _find_exec(executable_file, monitoring=None):
    """
    Try to find the executable file 'executable_file'
    :param executable_file:
    :param monitoring:
    :return:
    """
    cmd = 'which' + ' ' + str(executable_file)
    path_to_exec = ""
    try:
        which_exec = subprocess.check_output(cmd, shell=True)
        path_to_exec = which_exec.split('\n')[0]
    except subprocess.CalledProcessError:
        try_file = os.path.join(os.path.dirname(__file__), 'cpp', 'vt', 'build', 'bin', str(executable_file))
        if os.path.isfile(try_file):
            return try_file
        try_file = os.path.join(os.path.dirname(__file__), 'cpp', str(executable_file))
        if os.path.isfile(try_file):
            return try_file

        _write_error_msg("findExec: can not find executable '" + str(executable_file) + "'", monitoring)
        _write_error_msg("... Exiting", monitoring)

        sys.exit(1)

    return path_to_exec


def path_to_vt():
    """
    Try to find the executable file 'executable_file'
    :return:
    """
    cmd = 'which blockmatching'
    path_to_exec = ""
    try:
        which_exec = subprocess.check_output(cmd, shell=True)
        path_to_exec = which_exec.split('\n')[0]
    except subprocess.CalledProcessError:
        return None

    return os.path.dirname(path_to_exec)


############################################################
#
# function for the fusion steps
#
############################################################


def _launch_inline_cmd(command_line, monitoring=None):
    """

    :param command_line:
    :param monitoring:
    :return:
    """

    if monitoring is not None and (monitoring.verbose >= 3 or monitoring.debug > 0):
        monitoring.to_log("* Launch: " + command_line)
        with open(monitoring.logfile, 'a') as logfile:
            subprocess.call(command_line, shell=True, stdout=logfile, stderr=subprocess.STDOUT)
    else:
        subprocess.call(command_line, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    return


############################################################
#
# function for the fusion steps
#
############################################################


def apply_transformation(the_image, res_image, the_transformation=None,
                         template_image=None,
                         voxel_size=None,
                         dimensions=None,
                         interpolation_mode='linear',
                         cell_based_sigma=0.0,
                         monitoring=None,
                         return_image=False):
    """

    :param the_image: path to the image to be resampled
    :param res_image: path to the resampled image (ie result)
    :param the_transformation: path to the transformation to be used
    :param template_image: path to the template image (used to specify result image geometry)
    :param voxel_size: to specify the output voxel size(s); may be used to change image
                resolution by resampling
    :param dimensions: dimensions of the result image; may be used to change image
                resolution by resampling
    :param nearest: do not interpolate (take the nearest value)
             to be used when applying on label images (default = False)
    :param cell_based_sigma: smoothing parameter when resampling with nearest=True.
           cell_based_sigma=0.0 means no smoothing.
    :param monitoring: control structure (for verboseness and log informations)
    :param return_image: if True, return the result image as an spatial image
                  (default = False; nothing is returned)
    :return: no returned value if return_image = False
            if return_image = True, return the result image as an spatial image
    """

    proc = "apply_transformation"

    path_to_exec = _find_exec('applyTrsf')

    command_line = path_to_exec + " " + the_image + " " + res_image

    if the_transformation is not None:
        command_line += " -trsf " + the_transformation

    if template_image is not None:
        command_line += " -template " + template_image

    if type(interpolation_mode) == str:
        if interpolation_mode.lower() == 'linear':
            command_line += " -linear"
        elif interpolation_mode.lower() == 'nearest':
            if (type(cell_based_sigma) == int and cell_based_sigma > 0) \
                    or (type(cell_based_sigma) == float and cell_based_sigma > 0.0):
                command_line += " -cellbased"
                command_line += " -cell-based-sigma " + str(cell_based_sigma)
            else:
                command_line += " -nearest"
        else:
            # default
            pass
    else:
        # default
        pass

    if voxel_size is not None:
        if type(voxel_size) == int or type(voxel_size) == float:
            command_line += " -iso "+str(voxel_size)
        elif type(voxel_size) == tuple or type(voxel_size) == list:
            if len(voxel_size) == 3:
                command_line += " -vs "+str(voxel_size[0]) + " " + str(voxel_size[1]) + " " + str(voxel_size[2])
            else:
                _write_error_msg(proc + ": unhandled length for voxel_size '" + str(len(voxel_size)) + "'", monitoring)
                _write_error_msg("\t Exiting", monitoring)
                sys.exit(1)
        else:
            _write_error_msg(proc + ": unhandled type for voxel_size '" + str(type(voxel_size)) + "'", monitoring)
            _write_error_msg("\t Exiting", monitoring)
            sys.exit(1)

    if dimensions is not None:
        if type(dimensions) == tuple or type(dimensions) == list:
            if len(dimensions) == 3:
                command_line += " -vs " + str(dimensions[0]) + " " + str(dimensions[1]) + " " + str(dimensions[2])
            else:
                _write_error_msg(proc + ": unhandled length for dimensions '" + str(len(dimensions)) + "'", monitoring)
                _write_error_msg("\t Exiting", monitoring)
                sys.exit(1)
        else:
            _write_error_msg(proc + ": unhandled type for dimensions '" + str(type(dimensions)) + "'", monitoring)
            _write_error_msg("\t Exiting", monitoring)
            sys.exit(1)

    _launch_inline_cmd(command_line, monitoring=monitoring)

    if return_image is True:
        out = imread(res_image)
        return out

    return


def slitline_correction(the_image, res_image,
                        output_corrections=None, input_corrections=None,
                        other_options=None, monitoring=None):
    """

    :param the_image:
    :param res_image:
    :param output_corrections
    :param input_corrections
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('mc-removeLine')

    #
    # global method
    # -xz 0.5: 50% of y lines are used to compute slit line to be corrected
    #         the ones with the smallest intensity average (assumed to be background)
    # -y 0.1: robust mean (with 10% outliers) is computed along selected y lines
    # -c 0.2: lines that contains more than 20% of outliers are subject to correction
    #

    command_line = path_to_exec + " " + the_image + " " + res_image
    if input_corrections is None:
        command_line += " -method g -xz 0.5 -y 0.1 -c 0.2"
        if output_corrections is not None:
            command_line += " -output-corrections " + output_corrections
    else:
        command_line += " -input-corrections " + input_corrections

    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def mip_projection_for_crop(the_image, res_image, other_options=None, monitoring=None):
    """

    :param the_image:
    :param res_image:
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('mc-extractMIPembryo')

    command_line = path_to_exec + " " + the_image + " " + res_image

    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def linear_registration(path_ref, path_flo, path_output,
                        path_output_trsf, path_init_trsf=None,
                        py_hl=6, py_ll=3,
                        transformation_type='affine',
                        transformation_estimator='wlts',
                        lts_fraction=0.55,
                        normalization=True,
                        other_options=None,
                        monitoring=None):
    """
    Compute the transformation that register the floating image onto the reference image
    :param path_ref: path to the reference image
    :param path_flo: path to the floating image
    :param path_output: path to the floating image after registration and resampling
    :param path_output_trsf: path to the computed transformation

    :param py_ll: pyramid lowes
        :param path_init_trsf: path to the initial registration (default=None)
    :param py_hl: pyramid highest level (default = 6)

    t level (default = 3)
    :param transformation_type: type of transformation to be computed (default is 'affine')
    :param transformation_estimator: transformation estimator (default is 'wlts')
    :param lts_fraction: least trimmed squares fraction (default = 0.55)
    :param normalization:
    :param other_options: other options to be passed to 'blockmatching'
           see blockmatching options for details
    :param monitoring: control structure (for verboseness and log informations)
    :return: no returned value
    """

    path_to_exec = _find_exec('blockmatching')

    command_line = path_to_exec + " -ref " + path_ref + " -flo " + path_flo
    if path_output is not None:
        command_line += " -res " + path_output
    if path_init_trsf is not None:
        command_line += " -init-res-trsf " + path_init_trsf + " -composition-with-initial"
    command_line += " -res-trsf " + path_output_trsf

    command_line += " -pyramid-highest-level " + str(py_hl) + " -pyramid-lowest-level " + str(py_ll)
    command_line += " -trsf-type " + transformation_type
    command_line += " -estimator " + transformation_estimator
    command_line += " -lts-fraction " + str(lts_fraction)
    if normalization is False:
        # monitoring.to_log_and_console("       non-normalized registration", 2)
        command_line += " -no-normalisation"

    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def linear_combination(the_weights, the_images, res_image, other_options=None, monitoring=None):
    """

    :param the_weights:
    :param the_images:
    :param res_image:
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('mc-linearCombination')

    command_line = path_to_exec + " -weights"
    for i in range(len(the_weights)):
        command_line += " " + str(the_weights[i])

    command_line += " -images"
    for i in range(len(the_images)):
        command_line += " " + str(the_images[i])

    command_line += " -res " + str(res_image)

    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


############################################################
#
# functions for intra-registration
#
############################################################


def multiple_trsfs(format_input, format_output, first, last, reference, trsf_type='rigid', other_options=None,
                   monitoring=None):
    """

    :param format_input:
    :param format_output:
    :param first:
    :param last:
    :param reference:
    :param trsf_type:
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('multipleTrsfs')

    command_line = path_to_exec + " " + format_input + " -res " + format_output
    command_line += " -method propagation "
    command_line += " -first " + str(first) + " -last " + str(last)
    command_line += " -reference " + str(reference)
    command_line += " -trsf-type " + trsf_type

    #
    #
    #
    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def change_multiple_trsfs(format_input, format_output, first, last, reference, result_template, trsf_type='rigid',
                          resolution=None, threshold=None, margin=None, format_template=None, other_options=None,
                          monitoring=None):
    """

    :param format_input:
    :param format_output:
    :param first:
    :param last:
    :param reference:
    :param result_template:
    :param trsf_type:
    :param resolution:
    :param threshold:
    :param margin:
    :param format_template:
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('changeMultipleTrsfs')

    command_line = path_to_exec + " -format " + format_input + " -res-format " + format_output
    command_line += " -first " + str(first) + " -last " + str(last)
    command_line += " -reference " + str(reference)
    command_line += " -result-template " + str(result_template)

    command_line += " -trsf-type " + trsf_type

    if resolution is not None:
        if (type(resolution) is tuple or type(resolution) is list) and len(resolution) == 3:
            command_line += " -result-voxel-size "
            command_line += str(resolution[0]) + " " + str(resolution[1]) + " " + str(resolution[2])
        elif type(resolution) is int or type(resolution) is float:
            command_line += " -result-isotropic " + str(resolution)

    if threshold is not None:
        command_line += " -threshold " + str(threshold)

    if margin is not None:
        command_line += " -margin " + str(margin)

    if format_template is not None:
        command_line += " -template-format " + str(format_template)

    #
    #
    #
    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def crop_sequence(format_input, path_output, firstindex, lastindex, orientation, sliceindex, monitoring=None):
    """

    :param format_input:
    :param path_output:
    :param firstindex:
    :param lastindex:
    :param orientation:
    :param sliceindex:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('cropImage')
    command_line = path_to_exec + " -format " + format_input + " " + path_output
    command_line += " -first " + str(firstindex) + " -last " + str(lastindex)
    if orientation.lower() == 'xy':
        command_line += " -xy " + str(sliceindex)
    elif orientation.lower() == 'xz':
        command_line += " -xz " + str(sliceindex)
    elif orientation.lower() == 'yz':
        command_line += " -yz " + str(sliceindex)
    else:
        return

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


############################################################
#
# functions for the MARS step
#
############################################################


def linear_smoothing(path_input, path_output, filter_value=1.0, real_scale=False, filter_type='deriche',
                     other_options=None, monitoring=None):
    """

    :param path_input: path to the image to filter
    :param path_output: path to the output image
    :param filter_value: sigma of the gaussian filter for each axis (default is 1.0)
    :param real_scale: scale values are in 'real' units (will be divided by the voxel size to get 'voxel' values)
           if this option is at True (default=False)
    :param filter_type: gaussian type, can be ['deriche'|'fidrich'|'young-1995'|'young-2002'|...
           ...|'gabor-young-2002'|'convolution'] or None (default is 'deriche')
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('linearFilter')

    command_line = path_to_exec + " " + path_input + " " + path_output

    #
    # filter parameter value
    #
    command_line += " -sigma " + str(filter_value)
    if real_scale is True:
        command_line += " -unit real"
    else:
        command_line += " -unit voxel"

    #
    # filter type
    #
    if type is not None:
        command_line += " -gaussian-type " + str(filter_type)
    command_line += " -x 0 -y 0 -z 0"

    #
    # add points at borders
    #
    command_line += " -cont 10"

    #
    #
    #
    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def regional_maxima(path_input, path_output, h=1, other_options=None, monitoring=None):
    """

    :param path_input: path to the input image
    :param path_output: path to the output image
    :param h: h-maxima parameter value
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('regionalext')

    #
    #
    #
    command_line = path_to_exec + " " + path_input + " -diff " + path_output
    command_line += " -max -h " + str(h)

    #
    #
    #
    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def regional_minima(path_input, path_output, h=1, other_options=None, monitoring=None):
    """

    :param path_input: path to the input image
    :param path_output: path to the output image
    :param h: h-minima parameter value
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('regionalext')

    #
    #
    #
    command_line = path_to_exec + " " + path_input + " -diff " + path_output
    command_line += " -min -h " + str(h)

    #
    #
    #
    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def connected_components(path_input, path_output, low_threshold=1, high_threshold=None,
                         other_options=None, monitoring=None):
    """
    Label connected components
    :param path_input:
    :param path_output:
    :param low_threshold:
    :param high_threshold: high threshold for hysteresis thresholding
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('connexe')

    #
    #
    #
    command_line = path_to_exec + " " + path_input + " " + path_output
    command_line += " -lt " + str(low_threshold)
    if high_threshold is not None:
        command_line += " -ht " + str(high_threshold)

    #
    # force output type
    #
    command_line += " -labels"

    #
    # force output image type
    #
    command_line += " -o 2"

    #
    #
    #
    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def watershed(path_seeds, path_gradient, path_output, other_options=None, monitoring=None):
    """
    Perform the watershed operation
    :param path_seeds: path to the seeds image
    :param path_gradient: path to the intensity/gradient image
    :param path_output: path to the output image
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('watershed')

    #
    #
    #
    command_line = path_to_exec + " -seeds " + path_seeds + " -gradient " + path_gradient + " " + path_output

    #
    #
    #
    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def global_normalization_to_u8(path_input, path_output, min_percentile=0.01, max_percentile=0.99,
                               other_options=None, monitoring=None):
    """

    :param path_input:
    :param path_output:
    :param min_percentile:
    :param max_percentile:
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('mc-adhocFuse')

    #
    #
    #
    command_line = path_to_exec + " -intensity-image " + path_input + " -result-intensity-image " + path_output
    command_line += " -min-method global -max-method global"
    command_line += " -min-percentile " + str(min_percentile) + " -max-percentile " + str(max_percentile)
    #
    #
    #
    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def cell_normalization_to_u8(path_input, path_segmentation, path_output, min_percentile=0.01, max_percentile=0.99,
                             cell_normalization_min_method='cellinterior',
                             cell_normalization_max_method='cellborder',
                             other_options=None, monitoring=None):
    """

    :param path_input:
    :param path_segmentation:
    :param path_output:
    :param min_percentile:
    :param max_percentile:
    :param cell_normalization_min_method:
    :param cell_normalization_max_method:
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('mc-adhocFuse')

    #
    #
    #
    command_line = path_to_exec + " -intensity-image " + path_input + " -segmentation-image " + path_segmentation
    command_line += " -result-intensity-image " + path_output
    command_line += " -min-method " + cell_normalization_min_method
    command_line += " -max-method " + cell_normalization_max_method
    command_line += " -min-percentile " + str(min_percentile) + " -max-percentile " + str(max_percentile)
    command_line += " -sigma 5.0"

    #
    #
    #
    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


############################################################
#
# functions for the membrane detection and/or enhancement
#
############################################################


def membrane_extraction(path_input, prefix_output='tmp_membrane',
                        path_mask=None, scale=0.9, real_scale=True, other_options=None, monitoring=None):
    """
    Membrane (plane-like) structures enhancement and oriented centerplanes extraction
    (using the method detailed in [Michelin et al. 2014] inspired by [Krissian 2000])
    :param path_input: path to the image to filter
    :param prefix_output: Write 3 files named <prefix_output>.ext.inr,<prefix_output>.theta.inr,<prefix_output>.phi.inr
    :param path_mask: binary image (u8 or u16) such that the response function is only computed for non-null voxels
           from this mask
    :param scale: detection scale parameter (should be set as the semi-thickness of the membrane)
    :param real_scale: set as True if the scale parameter is given in real coordinates system (default),
           set as False if given in voxel coordinates
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('membrane')

    #
    #
    #
    command_line = path_to_exec + " " + path_input + " " + prefix_output
    if path_mask is not None and os.path.isfile(path_mask):
        command_line += " -mask " + path_mask
    command_line += " -single -init " + str(scale)
    if real_scale is True:
        command_line += " -real"

    #
    #
    #
    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def seuillage(path_input, path_output, low_threshold=1, high_threshold=None, other_options=None, monitoring=None):
    """

    :param path_input:
    :param path_output:
    :param low_threshold: low threshold (default: 1)
    :param high_threshold:
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('seuillage')

    #
    #
    #
    command_line = path_to_exec + " " + path_input + " " + path_output
    command_line += " -sb " + str(low_threshold)
    if high_threshold is not None:
        command_line += " -sh " + str(high_threshold)

    #
    #
    #
    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def anisotropic_histogram(path_input_prefix, path_output, path_mask=None,
                          manual=False, manual_sigma=7, sensitivity=0.98, other_options=None, monitoring=None):
    """
    Centerplanes image binarisation using an adaptative anisotropic threshold method detailed in [Michelin 2016]
    :param path_input_prefix: generic prefix to input data
    :param path_output: output binary image
    :param path_mask: binary image (u8 or u16) such that the thresholding is only computed for non-null voxels
           from this mask (8 bits image of same size as input image).
    :param manual: if True, enables manual initialisation of sigma value for histograms fitting (default: False)
    :param manual_sigma: the sigma value for histogram fitting in case of manual mode (default: 20)
    :param sensitivity: computes the anisotropic thresholds following a sensitivity criterion (true positive rate):
           threshold = #(membrane class >= threshold) / #(membrane class)
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('anisotropicHist')

    #
    #
    #
    command_line = path_to_exec + " " + path_input_prefix + ".ext.inr " + path_input_prefix + ".hist.txt"
    command_line += " -bin-out " + path_output
    if path_mask is not None and os.path.isfile(path_mask):
        command_line += " -mask " + path_mask
    command_line += " -sensitivity " + str(sensitivity)
    command_line += " -v"

    if manual is True:
        from math import exp
        amplitude = 1.0 / (manual_sigma * exp(-0.5))
        lmin = manual_sigma / 3.0
        lmax = manual_sigma * 5.0
        command_line += ' -rayleighcentered ' + str(amplitude) + ' ' + str(manual_sigma)
        command_line += ' -lmin ' + str(lmin) + ' -lmax ' + str(lmax)
    else:
        command_line += " -auto"

    #
    #
    #
    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def arithmetic_operation(path_first_input, path_second_input, path_output, other_options=None, monitoring=None):
    """

    :param path_first_input:
    :param path_second_input:
    :param path_output:
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('Arit')

    #
    #
    #
    command_line = path_to_exec + " " + path_first_input + " " + path_second_input + " " + path_output
    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def tensor_voting_membrane(path_input, prefix_input, path_output, path_mask=None,
                           scale_tensor_voting=3.6, sigma_smoothing=0.9, real_scale=True, sample=0.2, monitoring=None):
    """
    Grey-level image reconstruction from an image of binarised membranes associated to images of orientations.
    :param path_input: path to input image which is contains binarised planar structures,
           with associated input angle images <path_prefix>.theta.inr and <path_prefix>.phi.inr in the same folder.
    :param prefix_input: prefix name for temporary images
    :param path_output: path to the output image
    :param path_mask: mask image applied to the input image to restrict the domain of tensor voting tokens
           (must be same dimensions with input image)
    :param scale_tensor_voting: scaling for tensor voting method (only on isotropic images) (default = 3.6 um)
    :param sigma_smoothing: scale for gaussien smooting after tensor voting (default = 0.9 um)
    :param real_scale: True (default) if scale in real coordinates. False if scale in voxel coordinates.
    :param sample: multiplying parameter for decrease the time-cost of the function by diminishing the
           number of voting token; 0 < sample <= 1, default = 0.2
    :param monitoring:
    :return:
    """

    proc = "tensor_voting_membrane"

    #
    # check whether the input image is empty
    #
    path_to_exec = _find_exec('nonZerosImage')
    command_line = path_to_exec + " " + path_input
    isnonzero = subprocess.call(command_line, shell=True)
    if isnonzero == 0:
        if monitoring is not None:
            monitoring.to_log(proc + ": '" + str(path_input).split(os.path.sep)[-1] + "' only contains 0.")
            monitoring.to_log("\t Exiting.")
        else:
            print(proc + ": '" + str(path_input).split(os.path.sep)[-1] + "' only contains 0.")
            print("\t Exiting.")
        sys.exit(1)

    #
    # tensor voting
    #

    path_to_exec = _find_exec('TVmembrane')
    command_line = path_to_exec + " " + path_input + " -output-eigenvalues " + prefix_input
    if path_mask is not None and os.path.isfile(path_mask):
        command_line += " -mask " + path_mask
    command_line += " -scale " + str(scale_tensor_voting) + " -hessian"
    command_line += " -sample " + str(sample)

    _launch_inline_cmd(command_line, monitoring=monitoring)

    #
    # eigenvalues substraction
    #
    arithmetic_operation(prefix_input + ".imvp3.inr", prefix_input + ".imvp1.inr", prefix_input + ".tv.inr",
                         other_options='-sub')

    #
    # smoothing
    #
    #
    input_image = prefix_input + ".tv.inr"
    if sigma_smoothing > 0.0:
        linear_smoothing(input_image, prefix_input + ".lf.inr", filter_value=sigma_smoothing,
                         real_scale=real_scale, monitoring=monitoring)
        input_image = prefix_input + ".lf.inr"

    #
    # copy into 1-byte image
    #
    path_to_exec = _find_exec('copy')
    command_line = path_to_exec + " -norma -o 1 " + input_image + " " + path_output

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def bounding_boxes(image_labels, path_bboxes=None, monitoring=None):
    """
    Calcul des bounding-boxes de chaque label de l'image d'entree.
    Si path_bboxes est renseigne, le resultat est sauvegarde dans ce fichier.
    Output : dictionnaire D dont les cles sont les labels d'image et les
    valeurs sont les listes contenant dans l'ordre les informations de
    [volume, xmin, ymin, zmin, xmax, ymax, zmax] correspondant a ces labels
    avec volume > 0 et label > 0.
    :param image_labels:
    :param path_bboxes:
    :param monitoring:
    :return:
    """

    if path_bboxes is None:
        file_boxes = 'tmp_bounding_boxes.txt'
    else:
        file_boxes = path_bboxes

    #
    #
    #

    path_to_exec = _find_exec('boundingboxes')
    command_line = path_to_exec + " " + image_labels + " " + file_boxes
    _launch_inline_cmd(command_line, monitoring=monitoring)

    #
    #
    #

    f = open(file_boxes, 'r')
    lines = f.readlines()
    f.close()

    boxes = {}

    for line in lines:
        if not line.lstrip().startswith('#'):
            li = line.split()
            if int(li[1]):
                boxes[int(li[0])] = map(int, li[1:])

    if path_bboxes is None:
        os.remove(file_boxes)

    return boxes


def crop_image(path_input, path_output, bbox, monitoring=None):
    """
    crop an image on disk
    :param path_input:
    :param path_output:
    :param bbox: [volume, xmin, ymin, zmin, xmax, ymax, zmax]
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('cropImage')
    command_line = path_to_exec + " " + path_input + " " + path_output
    command_line += " -origin " + str(bbox[1]) + " " + str(bbox[2]) + " " + str(bbox[3])
    command_line += " -dim " + str(bbox[4]-bbox[1]+1) + " " + str(bbox[5]-bbox[2]+1) + " " + str(bbox[6]-bbox[3]+1)

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def mathematical_morphology(path_input, path_output, other_options=None, monitoring=None):
    """

    :param path_input:
    :param path_output:
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('morpho')
    command_line = path_to_exec + " " + path_input + " " + path_output
    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def logical_operation(path_first_input, path_second_input, path_output, other_options=None, monitoring=None):
    """

    :param path_first_input:
    :param path_second_input:
    :param path_output:
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('Logic')

    #
    #
    #
    command_line = path_to_exec
    if other_options is not None:
        command_line += " " + other_options

    command_line += " " + path_first_input

    if path_second_input is not None:
        command_line += " " + path_second_input

    command_line += " " + path_output

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def patch_logical_operation(path_first_input, path_second_input, path_output, bbox,
                            other_options=None, monitoring=None):
    """

    :param path_first_input:
    :param path_second_input:
    :param path_output:
    :param bbox: [volume, xmin, ymin, zmin, xmax, ymax, zmax]
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('patchLogic')

    #
    #
    #
    command_line = path_to_exec
    if other_options is not None:
        command_line += " " + other_options

    command_line += " " + path_first_input

    if path_second_input is not None:
        command_line += " " + path_second_input

    command_line += " " + path_output

    command_line += " -origin " + str(bbox[1]) + " " + str(bbox[2]) + " " + str(bbox[3])

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def create_image(image_out, image_template, other_options=None, monitoring=None):
    """

    :param image_out:
    :param image_template:
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('createImage')

    #
    #
    #
    command_line = path_to_exec + " " + image_out

    if image_template is not None:
        command_line += " -template " + image_template

    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


############################################################
#
#
#
############################################################


def non_linear_registration(path_ref, path_flo, path_affine, path_vector, affine_trsf, vector_trsf,
                            py_hl=5, py_ll=3,
                            transformation_estimator='wlts',
                            lts_fraction=0.55,
                            other_options=None,
                            monitoring=None):
    """
    Compute the non-linear transformation that register the floating image onto the reference image
    :param path_ref: path to the reference image
    :param path_flo: path to the floating image
    :param path_affine: path to the floating image after affine registration
    :param path_vector: path to the floating image after affine o non-linear registration
    :param affine_trsf: path to the affine transformation
    :param vector_trsf: path to the non-linear registration (affine o non-linear)
    :param py_hl:
    :param py_ll:
    :param transformation_estimator:
    :param lts_fraction:
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('blockmatching')

    #
    # affine registration
    #

    command_line = path_to_exec + " -ref " + path_ref + " -flo " + path_flo + " -res " + path_affine
    command_line += " -res-trsf " + affine_trsf

    command_line += " -pyramid-highest-level " + str(py_hl) + " -pyramid-lowest-level " + str(py_ll)

    command_line += " -trsf-type affine"

    command_line += " -estimator " + transformation_estimator
    command_line += " -lts-fraction " + str(lts_fraction)
    command_line += " -py-gf"

    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    #
    # non-linear registration
    #

    command_line = path_to_exec + " -ref " + path_ref + " -flo " + path_flo + " -res " + path_vector
    command_line += " -init-trsf " + affine_trsf
    command_line += " -res-trsf " + vector_trsf
    command_line += " -composition-with-initial"

    command_line += " -pyramid-highest-level " + str(py_hl) + " -pyramid-lowest-level " + str(py_ll)

    command_line += " -trsf-type vectorfield"

    command_line += " -estimator " + transformation_estimator
    #
    # was not set in previous version ...
    #
    # command_line += " -lts-fraction " + str(lts_fraction)
    #
    command_line += " -py-gf"
    command_line += " -elastic-sigma 2.0 2.0 2.0"
    command_line += " -fluid-sigma 2.0 2.0 2.0"

    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return


def only_keep_seeds_in_cell(seed_image, cell_image, seed_result,
                            other_options=None, monitoring=None):
    """

    :param seed_image:
    :param cell_image:
    :param seed_result:
    :param other_options:
    :param monitoring:
    :return:
    """

    path_to_exec = _find_exec('mc-maskSeeds')

    #
    #
    #
    command_line = path_to_exec

    command_line += " -seed-image " + seed_image
    command_line += " -cell-image " + cell_image
    command_line += " -output-image " + seed_result

    if other_options is not None:
        command_line += " " + other_options

    _launch_inline_cmd(command_line, monitoring=monitoring)

    return

############################################################
#
#
#
############################################################

def obsolete_seuillage(path_input, path_output='tmp_threshold.inr',sb=1, sh=None,grey=False,lazy=True, verbose=False):
    ''' Manual image threshold
    path_input : path to the image to threshold
    path_output: image threshold
    sb : low threshold (seuil bas) (defalut : 1)
    sh : high threshold (seuil haut) (optional)
    grey : set to True if the user wants to keep original voxel values for the thresholded ones.
           if set to False (default), the output image will be binary
           (out[i][j][k]=255 iif sb<=in[i][j][k](<=sh), out[i][j][k]=0 otherwise)
    lazy : do not return the output image if True
    '''
    options=""
    if sh:
      options+=" -sh " + str(sh)
    path_seuillage = _find_exec('seuillage')
    cmd=path_seuillage + ' ' + path_input + ' ' + path_output +\
              ' -sb '+str(sb) + options
    if verbose:
      print cmd
    os.system(cmd)
    if not lazy:
        out = imread(path_output)
        os.system('rm ' + path_output)
        return out

def obsolete_membrane_renforcement(path_input, prefix_output='tmp_membrane', path_mask=None, init=0.9, realScale=True, lazy=True, verbose=False):
    '''
    Membrane (plane-like) structures enhancement and oriented centerplanes extraction
    (using the method detailed in [Michelin et al. 2014] inspired by [Krissian 2000])
    path_input : path to the image to filter
    prefix_output : Write 3 files call <prefix_output>.ext.inr,<prefix_output>.theta.inr,<prefix_output>.phi.inr
    path_mask : binary image (u8 or u16) such that the response function is only computed for non-null voxels from this mask
    init : enhancement scale parameter (should be set as the semi-thickness of the membrane)
    realScale : set as True if the scale parameter is given in real coordinates system (defaut), set as False if given in voxel coordinates
    lazy : do not return the output image if True
    '''
    options=''
    if path_mask:
      assert(os.path.exists(path_mask))
      options += ' -mask ' + str(path_mask)
    if realScale:
      options += ' -real'

    path_membrane = _find_exec('membrane')
    cmd=path_membrane + ' ' + path_input + ' ' + prefix_output +\
              ' -single -init '+str(init) + options
    if verbose:
      print cmd
    os.system(cmd)
    if os.path.exists(prefix_output+'.rep.inr'):
      cmd='rm ' + prefix_output+'.rep.inr'
      if verbose:
        print cmd
      os.system(cmd)

    if not lazy:
        out_ext = imread(prefix_output+'.ext.inr')
        out_theta = imread(prefix_output+'.theta.inr')
        out_phi = imread(prefix_output+'.phi.inr')
        os.system('rm ' + prefix_output+'.*')
        return out_ext, out_theta, out_phi

def _recfilter(path_input, path_output='tmp.inr', filter_value=2, rad_min=1, lazy=False):
    ''' Perform a gaussian filtering on an intensity image
    path_input : path to the image to filter
    path_output : path to the temporary output image
    filter_value : sigma of the gaussian filter
    rad_min : TO REMOVE, NOT USED
    lazy : do not return the output image if True
    '''
    print "recfilter: WARNING: This function is obsolete. The user should replace its use by function linearfilter."
    os.system(path_filters + ' ' + path_input +\
              ' ' + path_output +\
              ' -cont 10 -sigma ' + str(filter_value) +\
              ' -x 0 -y 0 -z 0 -o 2')
    if not lazy:
        out = imread(path_output)
        os.system('rm ' + path_output)
        return out  


def obsolete_linearfilter(path_input, path_output='tmp.inr', filter_value=2, rad_min=1, realScale=False, type='deriche', verbose=False, lazy=False):
    ''' Perform a gaussian filtering on an intensity image
    path_input : path to the image to filter
    path_output : path to the temporary output image
    filter_value : sigma of the gaussian filter for each axis (default is 1.0)
    rad_min : TO REMOVE, NOT USED
    realScale : scale values are in 'real' units (will be divided by the voxel size to get 'voxel' values) if this option is at True (default=False)
    type : gaussian type, which can be ['deriche'|'fidrich'|'young-1995'|'young-2002'|'gabor-young-2002'|'convolution'] or None (default is 'deriche')
    lazy : do not return the output image if True
    '''
    opt=""
    path_linearfilters = _find_exec('linearFilter')
    if realScale:
      opt += " -unit real"
    else:
      opt += ' -unit voxel'
    if type:
      opt += " -gaussian-type " + str(type)
    cmd=path_linearfilters + ' ' + path_input +\
              ' ' + path_output +\
              ' -cont 10 -sigma ' + str(filter_value) + opt +\
              ' -x 0 -y 0 -z 0 -o 2'
    if verbose:
      print cmd
    os.system(cmd)
    if not lazy:
        out = imread(path_output)
        os.system('rm ' + path_output)
        return out

def _regionalmax(path_input, path_output, h_min):
    ''' Perform the h-minima operation on a given image
    path_input : path to the input image
    path_output : path to the output image
    h_min : h-minima parameter value
    '''
    temp_out=path_output.replace('.inr','_out_regionalmax.inr')
    if os.path.exists(path_regional_max):
      os.system(path_regional_max + ' ' + path_input +\
              ' -diff ' + path_output +' '+temp_out+\
              ' -h ' + str(h_min) +\
              ' -inv')
    else:
      if os.path.exists(path_regional_ext):
        print "Warning : using " + path_regional_ext + " instead of " + path_regional_max + " (not found)."
        os.system(path_regional_ext + ' ' + path_input +\
              ' -diff ' + path_output +' '+temp_out+\
              ' -h ' + str(h_min) +\
              ' -min')
      else:
        print "Error : did not found " + path_regional_max + " neither "+ path_regional_ext + " binary functions. Exiting."
        return
    
    os.system('rm ' + temp_out)



def _connexe(path_input, path_output, high_th):
    ''' Perform the connected componant operation
    path_input : path to the input image
    path_output : path to the output image
    high_th : high threshold for an hysteresis filtering on the componants
    '''
    os.system(path_connexe + ' ' + path_input +\
              ' ' + path_output +\
              ' -lt 1 -ht ' + str(high_th) +\
              ' -labels -o 2')


def obsolete_watershed(path_seeds, path_int, path_output=None, lazy=True, temporary_folder='', verbose=False):
    ''' Perform the watershed operation
    path_seeds : path to the seeds image
    path_int : path to the intensity image
    path_output : path to the output image
    lazy : do not return the output image if True
    temporary_folder : folder for temporary files 
                       (by default, temporary files are written in the workspace)
    '''
    path_watershed = _find_exec('watershed')
    cmd=""   
    if type(path_seeds)!=str:
        imsave(os.path.join(temporary_folder,"seeds.inr"), path_seeds)
        path_seeds = os.path.join(temporary_folder,"seeds.inr")
        cmd+=" "+path_seeds
    if type(path_int)!=str:
        imsave(os.path.join(temporary_folder,"intensity.inr"), path_int)
        path_int = os.path.join(temporary_folder,"intensity.inr")
        cmd+=" "+path_int
    if path_output is None:
        lazy = False
        path_output = os.path.join(temporary_folder,"seg.inr")
        cmd+=" "+path_output
 
    command=path_watershed + ' ' + path_seeds +\
              ' ' + path_int +\
              ' ' + path_output

    if verbose:
      print command
    os.system(command)

    if not lazy:
        out=imread(path_output)
        if cmd:
            cmd='rm '+cmd
            if verbose:
                print cmd
            os.system(cmd)
        return out


def _reech3d(im_path, output_shape):
    ''' Perform a resampling operation
    im_path : path to the image to resample
    output_shape : desired output shape
    '''
    tmp_file=im_path.replace('.inr','_temp.inr')
    os.system(path_reech3d + ' '+im_path+'  '+tmp_file+
              ' -x ' + str(int(output_shape[0])) +
              ' -y ' + str(int(output_shape[1])) +
              ' -z ' + str(int(output_shape[2])))
    out = imread(tmp_file)
    os.system('rm '+tmp_file)
    return  out

def _reech(path_flo, path_output, voxelsize):
    ''' Perform a resampling operation
    path_flo : path to the image to resample
    path_output : path to the output image
    voxelsize : size of the voxels in $\mu m**3$ (x, y, z)
    '''
    os.system(path_reech3d +
              " " + path_flo + 
              " " + path_output +
              " -linear"
              " -iso " + str(voxelsize))





def obsolete_non_linear_registration(image_file_flo,image_file_ref, affine_image, affine_trsf,vectorfield_image,vectorfield_trsf, verbose=False):
    ''' Compute the non-linear transformation that register the floating image onto the reference image
    image_file_flo : path to the floating image
    image_file_ref : path to the reference image
    affine_image : path to the floating image after affine registration
    affine_trsf : path to the affine transformation
    vectorfield_image : path to the floating image after affine o non-linear registration
    vectorfield_trsf : path to the non-linear registration (affine o non-linear)
    '''
    path_block = _find_exec('blockmatching')
    cmd=path_block +\
              " -ref " + image_file_ref+\
              " -flo " + image_file_flo+\
              " -res " + affine_image +\
              " -res-trsf " +affine_trsf+\
              " -trsf-type affine" +\
              " -estimator wlts" +\
              " -py-gf" +\
              " -pyramid-highest-level 5" +\
              " -pyramid-lowest-level 3" +\
              " -lts-fraction 0.55"
    if verbose:
      print cmd
    os.system(cmd)
    
    cmd=path_block +\
              " -ref " + image_file_ref+\
              " -flo " + image_file_flo+\
              " -res " + vectorfield_image+\
              " -res-trsf " + vectorfield_trsf+\
              " -init-trsf " +affine_trsf+\
              " -trsf-type vectorfield" +\
              " -estimator wlts" +\
              " -py-gf" +\
              " -pyramid-highest-level 5" +\
              " -pyramid-lowest-level 3" +\
              " -elastic-sigma 2.0 2.0 2.0" +\
              " -fluid-sigma 2.0 2.0 2.0" +\
              " -composition-with-initial"
    if verbose:
      print cmd
    os.system(cmd)
    

def _linear_registration(path_ref, path_flo, path_trsf, path_output, path_output_trsf, py_hl=6, py_ll=3, lts_fraction=0.55):
    ''' Compute the linear transformation that register the floating image onto the reference image
    path_flo : path to the floating image
    path_ref : path to the reference image
    path_output : path to the floating image after affine registration
    path_trsf : path to the initial registration
    path_output_trsf : path to the affine transformation
    py_hl : pyramid highest level (default = 6)
    py_ll : pyramid lowest level (default = 3)
    lts_fraction : least trimmed squares fraction (default = 0.55)
    '''
    os.system(path_block +
              " -ref " + path_ref +
              " -flo " + path_flo +
              " -res " + path_output +
              " -res-trsf " + path_output_trsf +
              " -init-res-trsf " + path_trsf +
              " -trsf-type affine" +
              " -composition-with-initial "
              " -estimator wlts" +
              " -pyramid-highest-level " + str(py_hl) +
              " -pyramid-lowest-level " + str(py_ll) +
              " -lts-fraction " + str(lts_fraction))

def rigid_registration(path_ref, path_flo, path_trsf, path_output, path_output_trsf, py_hl=6, py_ll=3, lts_fraction=0.55, verbose=False):
    ''' Compute the rigid transformation that register the floating image onto the reference image
    path_flo : path to the floating image
    path_ref : path to the reference image
    path_output : path to the floating image after affine registration
    path_trsf : path to the initial registration (can be set to None or "/dev/null" if one does not want to initialize the trsf)
    path_output_trsf : path to the affine transformation
    py_hl : pyramid highest level (default = 6)
    py_ll : pyramid lowest level (default = 3)
    lts_fraction : least trimmed squares fraction (default = 0.55)
    '''
    cmd=path_block + \
              " -ref " + path_ref + \
              " -flo " + path_flo + \
              " -res " + path_output + \
              " -res-trsf " + path_output_trsf + \
              " -trsf-type rigid" + \
              " -estimator wlts" + \
              " -pyramid-highest-level " + str(py_hl) + \
              " -pyramid-lowest-level " + str(py_ll) + \
              " -lts-fraction " + str(lts_fraction)
    if path_trsf :
      cmd = cmd + " -init-trsf " + path_trsf
    if verbose:
      print cmd
    os.system(cmd)

def apply_trsf(path_flo, path_trsf=None, path_output="tmp_seeds.inr",
               template=None, nearest=True, voxelsize=None, iso=None, dimensions=None, lazy=True, verbose=False):
    ''' Apply a transformation to a given image
    path_flo : path to the floating image
    path_trsf : path to the transformation
    path_output : path to the output image
    template : path to the template image
    voxelsize : to specify the output voxelsize 
    iso : iso=iso <=> voxelsize=(iso, iso, iso)
    dimensions : dimensions of the result image
    nearest : do not interpolate (take the nearest value) if True, to use when applying on label images (default = True)
    '''
    #assert path_trsf != None 

    path_to_exec = _find_exec('applyTrsf')
    # command_line = path_apply_trsf + " " + path_flo + " " + path_output
    command_line = path_to_exec + " " + path_flo + " " + path_output
    if path_trsf:
      command_line += " -trsf " + path_trsf
    if not template is None:
        command_line += " -template " + template
    if not nearest is None:
        command_line += " -nearest"
    if voxelsize:
      assert type(voxelsize)==tuple or type(voxelsize)==list 
      assert len(voxelsize)==3
      command_line += " -vs %f %f %f"%(voxelsize[0], voxelsize[1], voxelsize[2])
    if iso:
      command_line += " -iso %f"%float(iso)
    if dimensions:
      assert type(dimensions)==tuple or type(dimensions)==list 
      assert len(dimensions)==3
      command_line += " -dim %d %d %d"%(int(dimensions[0]), int(dimensions[1]), int(dimensions[2]))
    if verbose:
      print command_line
    os.system(command_line)
    if not lazy:
        out=imread(path_output)
        if path_output=='tmp_seeds.inr':
            os.system('rm -f tmp_seeds.inr')
        return out


def obsolete_find_local_minima(path_out, path_ref, h_min, mask=None, sigma=0.6, verbose=False):
    ''' Find local minima in an intensity image
    path_out : path to the output seeds image
    path_ref : path to the reference intensity image
    h_min : value of the h-minima operator value
    mask : mask on the intensity image
    sigma : value of the gaussian filter in voxels
    '''
    from os import path
    path_mask_out=path_out.replace('.inr','_mask_'+str(h_min)+'.inr')
    tmp_min=path_out.replace('.inr','_local_minima_out.inr')
    tmp_filt=path_out.replace('.inr','_local_minima_filter'+str(sigma)+'.inr') 
    if not path.exists(tmp_filt) and mask==None:
        #recfilter(path_ref, tmp_filt, filter_value=sigma, lazy=True)
        obsolete_linearfilter(path_ref, tmp_filt, filter_value=sigma, realScale=True, type='deriche', lazy=True, verbose=verbose)

    path_regional_ext = _find_exec('regionalext')
    if mask==None:
        cmd=path_regional_ext + ' ' + tmp_filt + ' ' +\
              ' -diff ' + path_mask_out + ' ' +\
              tmp_min + ' ' +\
              '-h ' + str(h_min) + ' ' +\
              '-min'
        if verbose:
            print cmd
        os.system(cmd)
    else:
        cmd=path_regional_ext + ' ' + mask + ' ' + \
            '-diff ' + path_mask_out + ' ' + \
            tmp_min + ' ' + \
            '-h ' + str(h_min) + ' -max'
        if verbose:
            print cmd
        os.system(cmd)

    path_connexe = _find_exec('connexe')
    cmd=path_connexe + ' ' + path_mask_out + ' ' +\
              path_out + ' ' +\
              '-sb 1 -sh ' + str(h_min) +\
              ' -labels -o 2'
    if verbose:
      print cmd
    os.system(cmd)
    try:
        im=imread(path_out.replace('\\', ''))
    except:
        im=None
    cmd='rm -f '+tmp_filt+' '+tmp_min
    if verbose:
      print cmd
    os.system(cmd);
    return im, path_mask_out

def obsolete_morpho(image_input,image_output,paramstre,verbose=False):
  ''' Morphological operation
  '''
  path_morpho = _find_exec('morpho')
  cmd=path_morpho+' '+image_input+' '+' '+image_output+' '+ paramstre
  if verbose:
    print cmd
  os.system(cmd)

  
  
  

def outer_detection(im_ref_tmp, radius, seg_ref_tmp):
    ''' Compute the detection of the outer of the embryo
    im_ref_tmp : intensity image for the outer detection (SpatialImage)
    radius : radius of the grey closing to perform
    seg_ref_tmp : segmented reference image (SpatialImage)
    '''
    path_morpho = _find_exec('morpho')
    from copy import deepcopy
    if radius!='0':
        imsave("tmp_bounds.inr", im_ref_tmp)
        os.system(path_morpho + " tmp_bounds.inr closed.inr -clo -R " + radius)
        im=imread("closed.inr")
    else:
        im=deepcopy(im_ref_tmp)
    imax=np.max(im)
    h=np.histogram(im, range=(0, imax), bins=imax)
    cumhist=np.cumsum(h[0][::-1]),h[1][::-1]
    vol=np.sum(seg_ref_tmp!=1)#*1.10
    low=np.max(cumhist[1][cumhist[0]>vol])
    im_th=np.zeros_like(im)
    #im=imread("closed.inr")
    im_th[im>=low]=1 # Cytoplasm
    if radius!='0':
        imsave("tmp.inr", SpatialImage(im_th))
        os.system((path_morpho + " tmp.inr closing.inr -clo -R " + radius))
    else:
        imsave("closing.inr", SpatialImage(im_th))
    os.system((path_morpho + " closing.inr erode.inr -ero -R 5"))
    imE=imread("closing.inr")
    imE=nd.binary_fill_holes(imE)
    mask=np.uint8(imE)
    bounds=nd.binary_dilation(mask, structure=nd.generate_binary_structure(3, 1))-mask
    im_refB=im_ref_tmp.copy()
    im_refB[bounds.astype(np.bool)]=np.max(im_ref_tmp)
    imsave('tmp.inr', SpatialImage(im_refB))
    os.system(path_filters + " tmp.inr out_bounds.inr -x 0 -y 0 -z 0 -sigma 1 -o 2")
    return imread('out_bounds.inr'), bounds.astype(np.bool)

def obsolete_gradient_norm(image_input,gradient_output):
    ''' Perform the gradient norm of an intensity image
    im_input : input image (SpatialImage)
    path_input : path to the input image
    path_output : path to the output image
    '''
    path_gradient_norm = _find_exec('norme_gradient')
    os.system(path_gradient_norm + ' ' + image_input + ' ' + gradient_output + ' -sigma 1')

def readMatrixFile(file,t=int,comments='#'):
    '''
    Reads a matrix file.
    Possibility to specify the matrix type (int <default>, float, ...) with option 't'.
    Possibility to specify the line comment characters with option 'comments' (initially set to '#').
    '''
    M=[]
    if os.path.exists(file):
        f=open(file)
        for line in f:
            li=line.strip()
            if not li.startswith(str(comments)):
                info=li.split()
                v=[]
                for val in info:
                    v.append(t(val))
                if len(v):
                    M.append(v)
        f.close()
    return M






def obsolete_anisotropicHist(path_input="temp_membrane.ext.inr", path_output='tmp_membrane.bin.inr', path_mask=None, manual=False,manual_sigma=7,sensitivity=0.98, keepAll=False, lazy=True, verbose=False):
    ''' binarisation des membranes par seuillage anisotropic adaptatif
    Centerplanes image binarisation using an adaptative anisotropic threshold method detailed in [Michelin 2016]
    Generate temp_membrane.bin.inr and path_output
    path_input : path to the image to filter
    path_output : path to the temporary output image
    path_mask : binary image (u8 or u16) such that the thresholding is only computed for non-null voxels from this mask 
                (8 bits image of same size as input image). 
    manual : if True, enables manual initialisation of sigma value for histograms fitting (default: False)
    manual_sigma : the sigma value for histogram fitting in case of manual mode (default: 20)
    sensitivity : computes the anisotropic thresholds following a sensitivity criterion (true positive rate) : threshold = #(membrane class >= threshold) / #(membrane class) 
    lazy : do not return the output image if True
    '''
    mask_option=''
    if path_mask != None:
      assert(os.path.exists(path_mask))
      mask_option = ' -mask ' + str(path_mask)

    assert(os.path.dirname(path_output) == '' or os.path.exists(os.path.dirname(path_output)))

    path_theta = ".".join(path_input.rstrip('.gz').split('.')[0:-2])+".theta.inr"
    path_phi = ".".join(path_input.rstrip('.gz').split('.')[0:-2])+".phi.inr"

    assert(os.path.exists(path_input) and os.path.exists(path_theta) and os.path.exists(path_phi))

    path_hist = ".".join(path_input.rstrip('.gz').split('.')[0:-2])+".hist.txt"
    path_anisotropicHist = _find_exec('anisotropicHist')
    if not manual:
        #For Autoparametrization
        cmd=path_anisotropicHist + ' ' + path_input + ' ' + path_hist + ' -bin-out ' + str(path_output) +' -v -auto -sensitivity '+str(sensitivity) + mask_option
        if verbose:
          print cmd
        os.system(cmd)
    else:
        #Manual
        from math import exp
        amplitude=1.0/(manual_sigma*exp(-0.5))
        lmin=manual_sigma/3.0
        lmax=manual_sigma*5.0
        cmd=path_anisotropicHist + ' ' + path_input + ' ' + path_hist + ' -bin-out ' + str(path_output) + ' -v -rayleighcentered '+str(amplitude)+' '+str(manual_sigma)+' -lmin '+str(lmin)+' -lmax '+str(lmax)+' -sensitivity '+str(sensitivity) + mask_option
        if verbose:
          print cmd
        os.system(cmd)

    
    if True:

      i=1
      path_output_theta = ".".join(path_output.rstrip('.gz').split('.')[0:i])+".theta.inr"
      path_output_phi = ".".join(path_output.rstrip('.gz').split('.')[0:i])+".phi.inr"
      if (not os.path.exists(path_output_theta)) or (not os.path.exists(path_output_phi)):
        if len(path_output.split(os.path.sep)[-1].rstrip('.gz').split('.')) > 2 :
          i=2
          path_output_theta = ".".join(path_output.rstrip('.gz').split('.')[0:i])+".theta.inr"
          path_output_phi = ".".join(path_output.rstrip('.gz').split('.')[0:i])+".phi.inr"

      if path_theta != path_output_theta:
        cmd=path_copy + " " + path_theta + " " + path_output_theta
        if verbose:
          print cmd
        os.system(cmd)
      if path_phi != path_output_phi:
        cmd=path_copy + " " + path_phi + " " + path_output_phi
        if verbose:
          print cmd
        os.system(cmd)

      if not keepAll:
        if os.path.exists(path_hist):
          cmd='rm ' + path_hist
          if verbose:
            print cmd
          os.system(cmd)

    if not lazy:
        out = imread(path_output)
        os.system('rm ' + path_output)
        return out        

def obsolete_nonZerosImage(image_in, verbose=False):
  """
  Returns True if image_in contains only zeros, False otherwise.  
  """
  path_non_zeros_image = _find_exec('nonZerosImage')
  cmd=path_non_zeros_image + ' ' + str(image_in) 

  if verbose:
    cmd = cmd + ' -v'
    print cmd

  return bool(os.system(cmd))

def obsolete_TVmembrane(path_input="temp_membrane.bin.inr", path_output='tmp_TVmembrane.VP31.inr', path_mask=None, scale=3.6, sigma_LF=0.9, realScale=True, sample=0.2, keepAll=False, lazy=True, verbose=False):
    '''
    Grey-level image reconstruction from an image of binarised membranes associated to images of orientations.
    path_input : path to input image which is contains binarised planar structures, 
                 with associated input angle images <path_prefix>.theta.inr and <path_prefix>.phi.inr in the same folder.
    path_output : path to the output images prefix (final reconstructed image is path_output+'.u8.inr')
    path_mask : mask image applied to the input image to restrict the domain of tensor voting tokens (must be same dimensions with input image)
    scale : scaling for tensor voting method (only on isotropic images) (default = 3.6 um)
    realScale : True (default) if scale in real coordinates. Falseif scale in voxel coordinates.
    sample : multiplying parameter for decrease the time-cost of the function by diminishing the number of voting token;
             0 < sample <= 1, default = 0.2
    sigma_LF : scale for gaussien smooting after tensor voting (default = 0.9 um)
    keepAll : option to be set to True to keep all the intermediary images (imvp*, tv, lf, u8...) (default = False)

    lazy : do not return the output image if True

    '''

    prefixe_output = ".".join(path_output.rstrip('.gz').split('.')[0:-1])
    mask_option=''
    if path_mask != None:
      assert(os.path.exists(path_mask))
      mask_option = ' -mask ' + str(path_mask)

    assert(os.path.exists(path_input))
    assert(os.path.exists(os.path.dirname(path_output)))

    if obsolete_nonZerosImage(path_input):
      # Tensor Voting
      tvOpt=""
      if realScale:
        tvOpt += ' -real'
      path_TVmembrane = _find_exec('TVmembrane')
      cmd=path_TVmembrane + ' ' + path_input +' -output-eigenvalues ' + prefixe_output + mask_option + ' -scale '+str(scale) +' -hessian -sample '+str(sample) + tvOpt
      if verbose:
        print cmd
      os.system(cmd)
      #Substract eigenvalues
      path_Arit = _find_exec('Arit')
      cmd=path_Arit + ' ' + prefixe_output+'.imvp3.inr' +' ' + prefixe_output+'.imvp1.inr ' +prefixe_output+'.tv.inr  -sub '
      if verbose:
        print cmd
      os.system(cmd)
      #Gaussian
      lfOpt=""
      if realScale: 
        lfOpt += " -unit real"
      else:
        lfOpt += " -unit voxel"
      path_linearFilter = _find_exec('linearFilter')
      cmd=path_linearFilter + ' ' + prefixe_output+'.tv.inr' +' ' +prefixe_output+'.lf.inr  -x 0 -y 0 -z 0 -sigma ' + str(sigma_LF) + lfOpt 
      if verbose:
        print cmd
      os.system(cmd)
      #u16 conversion  
      #os.system(path_copy + ' ' + path_output+'.lf.inr' +' ' +path_output+'.u16.inr  -norma -o 2')
      #u8 conversion
      path_copy = _find_exec('copy')
      cmd=path_copy + ' ' + prefixe_output+'.lf.inr' +' ' +path_output + ' -norma -o 1'
      if verbose:
        print cmd
      os.system(cmd)

      if not keepAll:
        #ext_to_erase = ['', '.lf.inr', '.tv.inr', '.imxx.inr', '.imxy.inr', '.imxz.inr', '.imyy.inr', '.imyz.inr', '.imzz.inr', 
        #                '.iszero.inr', '.imvp1.inr', '.imvp2.inr', '.imvp3.inr', '.imtheta1.inr', '.imtheta2.inr', '.imtheta3.inr', '.imphi1.inr', '.imphi2.inr', '.imphi3.inr']
        ext_to_erase = ['', '.lf.inr', '.tv.inr', '.imvp1.inr', '.imvp2.inr', '.imvp3.inr']
        files_to_erase=(' '+prefixe_output).join(ext_to_erase)
        files_to_erase.replace(' '+path_output.lstrip().rstrip(),'') # in case path_output would correspond to one of the file names automatically generated by TVmembrane...
        cmd='rm ' + files_to_erase
        if verbose:
          print cmd
        os.system(cmd)
    else:
      # The input image is only compounded of zeros : the output image will only contain zeros too
      if verbose:
        print "WARNING : TVmembrane INPUT IMAGE " + path_input + " SEEMS TO CONTAIN ONLY VOXELS WITH NULL VALUE /!\\"
      path_copy = _find_exec('copy')
      cmd=path_copy + ' ' + path_input + ' ' + path_output + ' -o 1'
      if verbose:
        print cmd
      os.system(cmd)

    if not lazy:
        out = imread(path_output)
        os.system('rm ' + path_output)
        return out  

def obsolete_copy(path_input, path_output, normalize=False, lazy=True, verbose=False):
    ''' copy of the input image
    path_input : path for input image
    path_output : path for output image
    normalize : if True, automatic normalization of input image before its copy
          (linear transformation (ax+b) of intensity

    lazy : do not return the output image if True

    '''
    path_copy = _find_exec('copy')
    if normalize:
      cmd=path_copy + ' ' + path_input + ' ' + path_output + '-norma -o 1'
      if verbose:
        print cmd
      os.system(cmd)
    else:
      cmd=path_copy + ' ' + path_input + ' ' + path_output
      if verbose:
        print cmd
      os.system(cmd)

    if not lazy:
        out = imread(path_output)
        os.system('rm ' + path_output)
        return out  


def directionHistogram(path_input,path_output,rayon=None,sigma=None,verbose=0,lazy=True):
    '''
    Compute the orientations histogram from an image of oriented binarized membranes
    
    path_input : path for input image
    path_output : path for output histogram which will be a 32-bits image

    optional parameters:
    rayon : ray parameter for output histogram (related to histogram granulosity)
    sigma :  standard deviation for voting accumulation (kernel density estimation) (defaut : PI/32)
    verbose : verbosity level

    lazy : do not return the output image if True
    '''
    options=''
    if rayon:
        options += ' -rayon ' + str(rayon)
    if sigma:
        options += ' -sigma ' + str(sigma)

    commande_shell=path_directionHistogram + ' ' + path_input + ' ' + path_output + ' ' + options

    if verbose:
        print commande_shell
    os.system(commande_shell)

    if not lazy:
        out=imread(path_output)
        os.system('rm ' + path_output)
        return out

def directionHistogramMaxima(path_input, file_out=None, maxima=None, frac=None, verbose=False):
  """
  Extract the maxima (modes) of the orientations histogram provided by the function 'directionHistogram'.

  path_input : input image containing the orientations histogram 
  file_out : output file ; if not specified, the function returns an array

  Optional parameters:

  maxima : list of integers where each integer %d corresponds to the %d-th highest direction maximum to be extracted (index starting at 0)
           (not specified by default)
  frac : fraction %lf between 0 and 1 so that local maxima are extracted iif their value is greater or equal to the global maximum multipied by the fraction %lf
         (default is not specified)
  Note for users : the options 'maxima' and 'frac' must not be used simultaneously. If none is specified, the function behaves like if frac=0.5

  verbose : level of verbosity

  """
  from numpy import loadtxt
  keep_file=True
  options=""
  if not file_out:
    keep_file=False
    file_out="tmp_directionHistogramMaxima.txt"
  if maxima:
    options += ' -max ' + str(maxima).lstrip('[').rstrip(']').replace(',','')
  if frac:
    options += ' -frac ' + str(frac)
  if verbose:
    options += ' -v'

  commande_shell=path_directionHistogramMaxima+" "+path_input+" "+file_out+" "+options
  if verbose:
    print commande_shell
  os.system(commande_shell)

  normals = loadtxt(file_out, delimiters=" ")
  if not keep_file:
    commande_shell="rm "+file_out
    if verbose:
      print commande_shell
    os.system(commande_shell)
  return normals

def dice(path_image, file_out=None, symmetry=None, binary=False, verbose=False, lazy=True):
  '''
  Genere un fichier texte synthetisant les dice label a label par rapport au plan de symetrie de l'equation specifiee 
  (argument 'symmetry' valant soit un str deja formate, soit une liste de 4 valeurs, soit le path du fichier contenant l'equation).
  Selon le mode, cela peut etre realise de label a label ou en binarisant les labels > 1.
  La fonction retourne trois elements :
   - un dictionnaire contenant tous les labels l de l'image en tant que keys, et pour chacune de ces keys, la value associee
     est un dictonnaire avec les keys suivantes :
      - 'dice_right_left' -> valeur = liste des valeurs de dice du label l reflechi de l'hemisphere droite vers l'hemisphere 
        gauche ([] si l est exclusivement dans l'hemisphere gauche)
      - 'dice_left_right' -> valeur = liste des valeurs de dice du label l reflechi de l'hemisphere gauche vers l'hemisphere 
        droite ([] si l est exclusivement dans l'hemisphere droite)
      - 'histo_right' -> valeur = histogramme de ce label l dans l'hemisphere droite (0 si l est exclusivement dans l'hemisphere gauche)
      - 'histo_left' -> valeur = histogramme de ce label l dans l'hemisphere gauche (0 si l est exclusivement dans l'hemisphere droite)
      - 'dice_max' -> valeur = maximum des dices de ce label (reflexion gauche-droite et droite-gauche confondues)
   - la moyenne des dice_max mesures sur l'ensemble des labels de l'image
   - l'ecart-type des dice_max mesures sur l'ensemble des labels de l'image
  Usage : dice [image-in] [image-ext] [filename-out] [-symmetry|s|n %f [%f %f %f]]
   [-graines|-seed|-intersection] [-bin] [-wi] [-v] [-D] [-help]
   si 'image-in' est '-', on prendra stdin
   si 'image-ext' est absent, on prendra stdin
   si 'filename-out' est absent, on prendra stdout
   si les trois sont absents, on prendra stdin et stdout
   -[graines|seed] : ecrit le tableau de confusion a la place des Dice
   -intersection : ecrit le tableau des intersections label/label
   -symmetry %d : calcule le dice entre les volumes a gauche et a droite du 
                  plan d'equation x = %d
   -n %f %f %f %f : calcule le dice entre les volumes de part et d'autre du 
                  plan d'equation n[0]*x+n[1]*y+n[2]*z+n[3] = 0
   -bin : calcul du dice sur les images binaires (sans les labels)
   -wi : ecrit toutes les images intermediaires
   -v : mode verbose
   -D : mode debug
  '''
  keep_file=True
  options=""
  if not file_out:
    keep_file=False
    file_out="tmp_dice.txt"
  if symmetry:
    if type(symmetry)==str and os.path.exists(symmetry): # means user gave a file-name of the equation of symmetry plane
      f=open(symmetry)
      n=[]
      for line in f:
        li=line.strip()
        #print li
        if not li.startswith('#'):
          #print li
          l=li.split()
          for value in l[1:5]:
            n.append(float(value))
    else:
      n=symmetry
    options += ' -n ' + str(n).replace('n:','').replace(',','').replace('[','').replace(']','').replace("'","")
  if binary:
    options += " -bin "

  commande_shell = path_dice + " " + path_image + ' ' + file_out + options 
  if verbose:
    print commande_shell
  os.system(commande_shell)

  if verbose:
    print "Read dices file "+file_out+" ..."
  dices, mean_dice_max, std_dev_dice_max=read_dices(file_out)

  if not keep_file:
    commande_shell="rm "+file_out
    if verbose:
      print commande_shell
    os.system(commande_shell)

  return dices, mean_dice_max, std_dev_dice_max

def read_dices(file):
  """
  La fonction retourne trois elements :
   - un dictionnaire contenant tous les labels l de l'image en tant que keys, et pour chacune de ces keys, la value associee
     est un dictonnaire avec les keys suivantes :
      - 'dice_right_left' -> valeur = liste des valeurs de dice du label l reflechi de l'hemisphere droite vers l'hemisphere 
        gauche ([] si l est exclusivement dans l'hemisphere gauche)
      - 'dice_left_right' -> valeur = liste des valeurs de dice du label l reflechi de l'hemisphere gauche vers l'hemisphere 
        droite ([] si l est exclusivement dans l'hemisphere droite)
      - 'histo_right' -> valeur = histogramme de ce label l dans l'hemisphere droite (0 si l est exclusivement dans l'hemisphere gauche)
      - 'histo_left' -> valeur = histogramme de ce label l dans l'hemisphere gauche (0 si l est exclusivement dans l'hemisphere droite)
      - 'dice_max' -> valeur = maximum des dices de ce label (reflexion gauche-droite et droite-gauche confondues)
   - la moyenne des dice_max mesures sur l'ensemble des labels de l'image
   - l'ecart-type des dice_max mesures sur l'ensemble des labels de l'image
  """
  from numpy import sqrt
  f = open(file, 'r')
  x = f.readlines()
  y=[line.strip('\n').split() for line in x if line.strip('\n').split() != []]
  assert(y[0][0] == 'Dice')
  labels_right=[int(val) for val in y[0][1:]]
  labels_left=[int(suby[0]) for suby in y[1:] if len(suby)>1]
  tab=[]
  histo_left=[]
  histo_right=[]
  flag_g=0
  flag_d=0
  for suby in y[1:]:
    if len(suby)==1:
      if suby[0]=="HistGauche":
        flag_g=1
        flag_d=0
      else:
        if suby[0]=="HistDroite":
          flag_g=0
          flag_d=1
        else:
          if flag_g:
            histo_left.append(int(suby[0]))
          else:
            if flag_d:
              histo_right.append(int(suby[0]))
            else:
              print "Error ? At line : " + suby
    else:
      assert(len(suby[1:])==len(labels_right))
      tab.append(suby[1:])

  assert(len(labels_right)==len(histo_right))
  assert(len(labels_left)==len(histo_left))
  assert(len(tab)==len(labels_left))
  dico={}
  for i in range(len(labels_right)):
    if not dico.has_key(labels_right[i]):
      dico[labels_right[i]]={}
    dico[labels_right[i]]['histo_right']=histo_right[i]
    dico[labels_right[i]]['dice_right_left']=[float(subt[i]) for subt in tab] 
    dico[labels_right[i]]['histo_left']=[]
    dico[labels_right[i]]['dice_left_right']=[]
  for i in range(len(labels_left)):
    if not dico.has_key(labels_left[i]):
      dico[labels_left[i]]={}
      dico[labels_left[i]]['histo_right']=[]
      dico[labels_left[i]]['dice_right_left']=[] 
    dico[labels_left[i]]['histo_left']=histo_left[i]
    dico[labels_left[i]]['dice_left_right']=[float(val) for val in tab[i] ]
  mean_dice_max=0
  std_dev_dice_max=0
  for k in dico.keys():
    dico[k]['dice_max']=max(dico[k]['dice_left_right']+dico[k]['dice_right_left'])
    mean_dice_max = mean_dice_max + dico[k]['dice_max']
  mean_dice_max = mean_dice_max / len(dico)
  for k in dico.keys():
    std_dev_dice_max = std_dev_dice_max + (dico[k]['dice_max'] - mean_dice_max)**2
  std_dev_dice_max = sqrt(std_dev_dice_max/len(dico))

  return dico, mean_dice_max, std_dev_dice_max


def diceMaximisation(path_input_seg, path_fileout, normal, delta=10, voxel=False, path_plane=None, verbose=False):
  '''
  Ajuste (par translation) la position d'un plan de symetrie initialise par le parametre "normal"
  si 'image-in' est '-', on prendra stdin
  si 'image-ext' est absent, on prendra stdin
  si 'filename-out' est absent, on prendra stdout
  si les trois sont absents, on prendra stdin et stdout
  -symmetry %f : calcule le dice entre les volumes a gauche et a droite du 
                 plan d'equation x = %f  -> not implemented here
  -n %f %f %f %f : calcule le dice entre les volumes de part et 
                 d'autre du plan d'equation %f*x+%f*y+%f*z+%f = 0   -> normal
  -voxel : equation du plan en coordonnees voxelliques  -> voxel
  -delta %d : fenetre de calcul du maximum de dice  -> delta
  -plane %s : sauve dans %s l'image du plan  -> path_plane
  -v : mode verbose
  -D : mode debug


  '''


def symmetryPlane(path_input_image, path_input_sphere=None, normal=None, equation_output='tmp_symmetryPlane', plane_output=None, distribution_output=None, trsf_output=None, 
                  maximum=None, d=None, dmin=None, delta=None, p=None, sigma=None, realSize=True, iterations=None, verbose=0, lazy=True):
    '''
    Calcule l'equation du plan de symetrie d'une image de membranes binarisees et orientees 
    path_input : path de l'image binarisee (implique l'existence des images d'angles associees selon la convention de nomenclature habituelle)
    path_output :

    # Exemple d'utilisation de l'executable dans un shell :
    symmetryPlane bin/bin_t004.inr -weq sym/sym_eq_t004.txt -sphere sym/directionHistogram_t004.inr -max 0 -plane sym/plane_t004.inr.gz -distribution sym/distribution_t004.inr 

    ### Aide de la commande shell 'symmetryPlane' ###

   -trsf %s : calcule la transformation permettant d'aligner le plan avec les
              axes de l'image --> trsf_output
   -angles %lf %lf : angles de la direction consideree (theta + phi) --> NOT HANDLED YET
   -n %lf %lf %lf : composantes de la direction consideree --> normal
   -n %lf %lf %lf %lf : fixe le plan de depart pour l'algorithme least square --> normal
   -sigma %lf : ecart-type pour la ponderation selon l'ecart angulaire (defaut : PI/64) --> sigma
   -sphere %s : direction deduite du fichier histogramme directionnel %s --> path_input_sphere
   -max %d : direction extraite la %d-eme plus grande (0=plus grand: default) --> maximum
   -bin : aucune ponderation de la distribution --> NOT HANDLED YET
   -plane %s : calcule la position du plan de symetrie --> plane_output
   -weq %s : ecrit dans %s l'equation du plan [-voxel : equation en coor. reelles] --> equation_output
   -iter %d : nombre d'iterations maxi --> iterations
   -d %lf : option de ponderation des points par une gaussienne d'e-t %lf 
            par rapport a la distance au plan de l'iteration precedante --> d
   -delta %lf : pas de diminution de l'e-t des distances (default: 2.0) --> delta
   -dmin %d : e-t minimal pour la distance au plan pour le least square (default: 5.0) --> dmin
   -real-d : distances in real coordinates (default) --> realSize = True (default)
   -vox[el]-d : distances in voxel coordinates --> realSize=False
   -p %lf : option de ponderation des points par une gaussienne d'e-t %lf 
            par rapport a l'ecart angulaire au plan de l'iteration precedante (default: PI/64) --> p
   -distribution %s : sauve l'image des distributions de l'angle principal --> distribution_output
   -v : mode verbose --> verbose

    '''


    ###### Exemple de chaine d'execution de commandes shell #######
    # symmetryPlane ../BIN/bin_t0${i}_on_t099.inr -sphere ../HISTO/bin_t0${i}_on_t099_R15A32.inr -sigma 0.1 -weq ../SYM/planeEq_t0${i}_on_t099_max_0.txt -d 0 -dmin 1 
    #
    # diceMaximisation ../WAT/OUT_t0${i}_on_t099-wat.inr ../SYM/planeEq_t0${i}_on_t099_max_0_dmax.txt -n `cat ../SYM/planeEq_t0${i}_on_t099_max_0.txt | grep -v "#" | awk -F':' '{print $2}'` -delta 10 
    #
    # symmetryPlane ../BIN/bin_t0${i}_on_t099.inr -n `cat ../SYM/planeEq_t0${i}_on_t099_max_0_dmax.txt | grep new | awk -F':' '{print $2}'` -weq ../SYM/planeEq_t0${i}_on_t099_max_0_dmax_a.txt -d 10 -dmin 10 -p 0.1 
    
    options = ''
    if equation_output:
        options += ' -weq ' + equation_output 
    if normal:
        options += ' -n ' + str(normal).replace('n:','').replace(',','').replace('[','').replace(']','').replace("'","") # normal peut etre le 'string' deja formate ou bien une list ou array contenant les 4 flottants attendus
    if path_input_sphere:
        options += ' -sphere ' + path_input_sphere
    if plane_output:
        options += ' -plane ' + plane_output
    if distribution_output:
        options += ' -distribution ' + distribution_output
    if trsf_output:
        options += ' -trsf ' + trsf_output
    if maximum:
        options += ' -max ' + str(maximum)
    if d:
        options += " -d " + str(d)
    if dmin:
        options += ' -dmin ' + str(dmin)
    if delta:
        options += " -delta " + str(delta)
    if p:
        options += ' -p ' + str(p)
    if sigma:
        options += ' -sigma ' + str(sigma)
    if realSize:
        options += ' -real-d'
    else:
        options += ' -voxel-d'
    if iterations:
        options += ' -iter ' + str(iterations)
    if verbose:
        options += ' -v'

    # Calcul du plan de symetrie
    commande_shell=path_symmetryPlane + ' ' + path_input_image + ' '  + options
    if verbose:
        print commande_shell
    os.system(commande_shell)

    if not lazy:
        fout=open(equation_output)
        out=[]
        for line in fout:
            li=line.strip()
            #print li
            if not li.startswith('#'):
                print li
                l=li.split()
                for value in l[1:5]:
                   print value
                   #print str(float(value))
                   out.append(float(value))
        fout.close()
        os.system('rm ' + equation_output)
        return out

def pointCloudRegistration(path_label_ref, path_label_flo, labels_pairing, 
                           skip_not_found=False, background_ref = 1, background_flo = 1,
                           path_trsf_ref_to_flo='', path_residuals='', path_dices='',  path_pairs_ref_flo='', 
                           trsf_type='affine', estimator='lts', lts_fraction=1.0, verbose=0, bash_options=None, lazy=True):
    """
  Realise le recalage de deux nuages de points labellises astreintes a la contrainte d'une initialisation des appariements de points.
    trsf_type = ['affine' | 'rigid']
    estimator = ['lts']
    Les options "path_label_ref" et "path_label_flo" doivent contenir une variable de type str (nom de fichier contenant les appariements de barycentres de labels au format text avec un barycentre par ligne)
    ou bien de type dict (retour de la fonction Morpheme_lineage.extract_label_barycenters_at_time)    
    L'option "labels_pairing" doit contenir une variable de type str (nom de fichier contenant les appariements de labels au format texte avec un appariement par ligne)
    ou bien de type dict (keys = labels de la reference, values = labels associes du flottant)
  # pointCloudRegistration -label-ref ../../../Data/140317-Patrick-St8/FUSE/SEG/POST/140317-Patrick-St8_fuse_seg_post_t030.inr -label-flo
  #   ../../../Data/160707-Ralph-St8/FUSE/SEG/POST/160707-Ralph-St8_fuse_seg_post_t020.inr -pairs rigid_Patrick_t030_Ralph_t020.pairs -affine  -estimator lts -lts-fraction 0.9 -pairs-out foo2.pairs -trsf foo2.trsf -background-ref 1 -background-flo 1
  # pointCloudRegistration -help
  Usage : pointCloudRegistration [-label-[ref|in] %s] [-label-[flo|ext] %s] [-background-[ref|in] %d] [-background-[flo|ext] %d]
   [-pairs|pair|p] %d %d [...] | %s] [-skip[-not-found]]
   [-rigid | -affine] 
   [-estimator-type|-estimator|-es-type wlts|lts|wls|ls]
   [-lts-fraction %lf] [-lts-deviation %f] [-lts-iterations %d]
   [-fluid-sigma|-lts-sigma[-ll|-hl] %lf %lf %lf]
   [-pairs-out %s] [-pairs-in %s] [-trsf %s] [-dices %s] [-residuals %s]
   [-inv] [-swap] [-v] [-D] [-help]
   -label-[ref|flo] %s : nom des fichiers image ou texte de labels ref|flo
      N.B.: le cas echeant, le fichier texte attendu doit respecter le format suivant : 
      1ere ligne (optionnel) : en-tete au format "Voxelsize %f %f %f" stipulant le voxelsize d'origine de l'image de laquelle proviennent les coordonnees barycentriques des labels
      liste des barycentres  : une ligne par barycentre ; chaque ligne respecte le format "%d %f %f %f" ou "%d %f %f %f %f" correspondant au label, coordonnee en x, en y, en z et (le cas echeant) au volume (ou poids) du label 
      lignes de commentaires : commencent avec le caractere "#"
   -background-[ref|flo] %d : label de fond a ignorer dans l'image [ref|flo] correspondante 
   -pairs %d %d [...] | %s : associe les paires de labels specifiees (ref-flo [ref-flo [...]])
   -skip-not-found : poursuit le calcul du recalage malgre d'eventuelles associations de labels inexistantes dans les labels d'entree
   -trsf %s : fichier dans lequel est enregistree la transformation T_flo<-ref calculee
   -pairs-in %s : fichier dans lequel sont enregistres les appariements utilises pour le calcul
   -pairs-out %s : fichier dans lequel sont enregistres les appariements trouves apres calcul
   ou seulement sur l'iteration donnee
   -residuals %s : enregistre les valeurs de residus pour les paires de labels mis en correspondance
   -dices %s : calcule l'indice de similarite de Dice entre chaque paire de labels en correspondance ;
               OPTION VALABLE UNIQUEMENT AVEC DES PARAMETRES -label-[ref|flo] DE TYPE IMAGE 
   ### transformation type ###
   -rigid : computes rigid transformation (set as default)
   -affine : computes affine transformation 
   ### transformation estimation ###
   [-estimator-type|-estimator|-es-type %s] # transformation estimator
   wlts: weighted least trimmed squares
   lts: least trimmed squares
   wls: weighted least squares
   ls: least squares
   [-lts-fraction %lf] # for trimmed estimations, fraction of pairs that are kept
   [-lts-deviation %lf] # for trimmed estimations, defines the threshold to discard
   pairings, ie 'average + this_value * standard_deviation'
   [-lts-iterations %d] # for trimmed estimations, the maximal number of iterations
   [-fluid-sigma|-lts-sigma] %lf %lf %lf] # sigma for fluid regularization,
   ie field interpolation and regularization for pairings (only for vector field)
   -inv : inverse 'image-in'
   -swap : swap 'image-in' (si elle est codee sur 2 octets)
   -v : mode verbose
   -D : mode debug

    """
    assert(path_label_ref and path_label_flo and labels_pairing)
    tmp_ref='tmp_labels_ref.tmp'
    tmp_flo='tmp_labels_flo.tmp'
    tmp_file='tmp_labels_pairing.tmp'
    tmp_lazy_trsf='tmp_trsf.tmp'
    tmp_lazy_pairs='tmp_pairs.tmp'
    tmp_lazy_dices=''
    tmp_lazy_residuals=''
    if not lazy:
        path_trsf_ref_to_flo=tmp_lazy_trsf
        path_dices=tmp_lazy_dices
        path_residuals=tmp_lazy_residuals
        path_pairs_ref_flo=tmp_lazy_pairs
    options = ''
    if path_label_ref:
        if type(path_label_ref)==dict:
            f = open(tmp_ref, 'w')
            for k,v in path_label_ref.iteritems():
                print >>f, k, v[0], v[1], v[2]
            f.close()
            options += " -label-ref " + str(tmp_ref)
        else:
            options += ' -label-ref ' + str(path_label_ref)
    if path_label_flo:
        if type(path_label_flo)==dict:
            f = open(tmp_flo, 'w')
            for k,v in path_label_flo.iteritems():
                print >>f, k, v[0], v[1], v[2]
            f.close()
            options += " -label-flo " + str(tmp_flo)
        else:
            options += ' -label-flo ' + str(path_label_flo)
    if background_ref:
        options += ' -background-ref ' + str(background_ref)
    if background_flo:
        options += ' -background-flo ' + str(background_flo)
    if labels_pairing:
        if type(labels_pairing)==dict:
            f = open(tmp_file, 'w')
            for k,v in labels_pairing.iteritems():
                print >>f, k, v
            f.close()
            options += " -pairs " + str(tmp_file)
        else:
            options += " -pairs " + str(labels_pairing)
    if skip_not_found:
        options += ' -skip-not-found'
    if trsf_type:
        options += ' -'+(str(trsf_type).strip().lstrip('-'))
    if estimator:
        options += ' -estimator ' + str(estimator)
    if lts_fraction:
        options += ' -lts-fraction ' + str(lts_fraction)
    if path_trsf_ref_to_flo:
        options += ' -trsf ' + str(path_trsf_ref_to_flo)
    if path_residuals:
        options += ' -residuals ' + str(path_residuals)
    if path_dices:
        options += ' -dices ' + str(path_dices)
    if path_pairs_ref_flo:
        options += ' -pairs-out ' + str(path_pairs_ref_flo)

    if verbose:
        options += ' -v'
    if bash_options:
        options += ' ' + bash_options

    commande_shell = path_pointCloudRegistration + ' ' + options
    if verbose:
        print commande_shell
    os.system(commande_shell)
    if type(path_label_ref)==dict:
        cmd='rm ' + tmp_ref
        if verbose:
            print cmd
        os.system(cmd)
    if type(path_label_flo)==dict:
        cmd='rm ' + tmp_flo
        if verbose:
            print cmd
        os.system(cmd)
    if type(labels_pairing)==dict:
        cmd='rm ' + tmp_file
        if verbose:
            print cmd
        os.system(cmd)
    if not lazy:
        out={}
        lut={}
        p=[]
        residuals=[]
        if os.path.exists(str(path_pairs_ref_flo)):
            p=readMatrixFile(str(path_pairs_ref_flo))
            os.system('rm ' + path_pairs_ref_flo)
            out['pairs']=p
        #if os.path.exists(path_residuals):
        #    r=readMatrixFile(str(path_residuals),float)
        #    os.system('rm ' + path_pairs_ref_flo)
        #    out['residuals']=r
        #if os.path.exists(str(path_dices)):
        #    d=readMatrixFile(str(path_dices))
        #    os.system('rm ' + path_dices)
        #    out['dices']=d
        if os.path.exists(str(path_trsf_ref_to_flo)):
            t=readMatrixFile(str(path_trsf_ref_to_flo), float)
            os.system('rm ' + path_trsf_ref_to_flo)
            out['trsf']=t

        return out

def planeRegistration(path_label_ref, path_label_flo, plane_eq_ref, plane_eq_flo, path_trsf_ref_to_flo='', path_residuals='', path_dices='', path_pairs_ref_flo='tmp_planeRegistration_pairs', 
                      background_ref = 1, background_flo = 1,
                      trsf_type='affine', estimator='lts', lts_fraction=1.0, verbose=0, bash_options=None, lazy=True):
    '''
    Realise le recalage de deux images labellisees (le label de fond doit valoir 0) astreintes a la contrainte d'une superposition initiale de "plans de symetrie" prealablement calcules

    trsf_type = ['affine' | 'rigid']
    estimator = ['lts']
    '''
    #   planeRegistration -label-ext ../WAT/OUT_t0${i}_on_t099-wat.inr -label-in ../../../WAT/wat_t005_on_t000_noborder.hdr -p-in `cat ../../../SYM/PlaneEq_t005_on_t000_max_0_a.txt` -p-ext `cat ../SYM/planeEq_t0${i}_on_t099_max_0_dmax_a.txt | grep -v "#" | ...
    #   ... awk -F':' '{print $2}'` -affine -estimator lts -lts-fraction 0.8 -pairs-out paires_T005-t0${i}_on_t099.txt -trsf trsf_T005-t0${i}_on_t099.txt -residuals residus_T005-t0${i}_on_t099.m 

    # planeRegistration -label-ref POST/140317-Patrick-St8_fuse_seg_post_t001.inr -label-flo POST/160707-Ralph-St8_fuse_seg_post_t004.inr -background-ref 1 -background-flo 1 \
    #      -p-ref -0.913401 0.074868 0.400117 26.907440 -p-flo -0.914614 -0.353682 -0.195935 122.856049 -affine -es-type lts -lts-fraction 0.9 \
    #      -trsf WORK/inter/Ralph_t004_Patrick_t001.trsf -pairs-out WORK/inter/Ralph_t004_Patrick_t001.pairs \
    #      -residuals WORK/inter/Ralph_t004_Patrick_t001_residuals.m -dices WORK/inter/Ralph_t004_Patrick_t001_dices.txt -v

    tmp_ref='tmp_labels_ref.tmp'
    tmp_flo='tmp_labels_flo.tmp'
    tmp_file='tmp_labels_pairing.tmp'
    tmp_lazy_trsf='tmp_trsf.tmp'
    tmp_lazy_pairs='tmp_pairs.tmp'
    tmp_lazy_dices=''
    tmp_lazy_residuals=''
    if not lazy:
        path_trsf_ref_to_flo=tmp_lazy_trsf
        path_dices=tmp_lazy_dices
        path_residuals=tmp_lazy_residuals
        path_pairs_ref_flo=tmp_lazy_pairs

    options = ''
    if path_label_ref:
        options += ' -label-ref ' + str(path_label_ref)
    if path_label_flo:
        options += ' -label-flo ' + str(path_label_flo)
    if plane_eq_ref:
        n=''
        if os.path.exists(str(plane_eq_ref)):
            file=open(str(plane_eq_ref))
            for line in file:
                li=line.strip()
                if not li.startswith('#'):
                    n=li
            file.close()
        else:
            n=plane_eq_ref
        options += ' -p-ref ' + str(n).replace('n:','').replace(',','').replace('[','').replace(']','').replace("'","") # normal peut etre le 'string' deja formate ou bien un vecteur contenant les 4 flottants attendus
    if plane_eq_flo:
        if os.path.exists(str(plane_eq_flo)):
            file=open(str(plane_eq_flo))
            for line in file:
                li=line.strip()
                if not li.startswith('#'):
                    n=li
            file.close()
        else:
            n=plane_eq_flo
        options += ' -p-flo ' + str(n).replace('n:','').replace(',','').replace('[','').replace(']','').replace("'","") # normal peut etre le 'string' deja formate ou bien un vecteur contenant les 4 flottants attendus
    if trsf_type:
        options += ' -'+(str(trsf_type).strip().lstrip('-'))
    if estimator:
        options += ' -estimator ' + str(estimator)
    if lts_fraction:
        options += ' -lts-fraction ' + str(lts_fraction)
    if path_trsf_ref_to_flo:
        options += ' -trsf ' + str(path_trsf_ref_to_flo)
    if path_residuals:
        options += ' -residuals ' + str(path_residuals)
    if path_dices:
        options += ' -dices ' + str(path_dices)
    if path_pairs_ref_flo:
        options += ' -pairs-out ' + str(path_pairs_ref_flo)
    if background_ref:
        options += ' -background-ref ' + str(background_ref)
    if background_flo:
        options += ' -background-flo ' + str(background_flo)

    if verbose:
        options += ' -v'
    if bash_options:
        options += ' ' + bash_options

    commande_shell = path_planeRegistration + ' ' + options
    if verbose:
        print commande_shell
    os.system(commande_shell)

    if not lazy:
        out={}
        lut={}
        p=[]
        residuals=[]
        if os.path.exists(str(path_pairs_ref_flo)):
            p=readMatrixFile(str(path_pairs_ref_flo))
            os.system('rm ' + path_pairs_ref_flo)
            out['pairs']=p
        #if os.path.exists(path_residuals):
        #    r=readMatrixFile(str(path_residuals),float, comments="%")
        #    os.system('rm ' + path_pairs_ref_flo)
        #    out['residuals']=r
        #if os.path.exists(str(path_dices)):
        #    d=readMatrixFile(str(path_dices))
        #    os.system('rm ' + path_dices)
        #    out['dices']=d
        if os.path.exists(str(path_trsf_ref_to_flo)):
            t=readMatrixFile(str(path_trsf_ref_to_flo), float)
            os.system('rm ' + path_trsf_ref_to_flo)
            out['trsf']=t

        return out

def associateLabels(path_ref, path_flo, path_pairs_ref_flo, path_ref_out, path_flo_out, path_labels_out=None, zeros=True, ref=False, force=False, visu=False, verbose=False):    
  '''
   path_ref : correspond au parametre path_label_ref de la methode planeRegistration
   path_flo : correspond au parametre path_label_flo de la methode planeRegistration
   path_pairs_ref_flo : fichier de correspondances label-ref / label-flo (tel qu'en sortie de planeRegistration path_pairs_ref_flo %s)
   zeros : les zeros en entree sont conserves en sortie si True
   force : force la valeur de sortie des labels (tableau d'entree n*3) si True
   ref : utilise les labels de l'image in comme labels de reference si True
   verbose : option de verbosite
  '''
  assert(os.path.exists(path_ref))
  assert(os.path.exists(path_flo))
  assert(os.path.exists(path_pairs_ref_flo))
  if not path_labels_out:
    path_labels_out="/dev/null"

  options=''
  if ref:
    options += ' -ref'
  if force:
    options += ' -f'
  if zeros:
    options += ' -z'

  cmd=path_associateLabels + ' ' + path_ref + ' ' + path_flo  + ' ' + path_pairs_ref_flo + ' ' + path_ref_out + ' ' + path_flo_out + ' ' + path_labels_out + ' ' + options
  if verbose:
    print cmd
  os.system(cmd)

def fuseLabelsWithLUT(image_in, image_out, path_lut, U8=False, lazy=True, verbose=False):
  """
  transmet en parametre un fichier de look-up-table dans lequel chaque ligne respecte le format : <ancien_label> <nouveau_label> (/!\ contraintes : uniquement des entiers positifs ; uniquement des images entree/sortie encodees sur le meme nombre de bits /!\).
  U8 option : forces U8 conversion of image_out if True (default=False)
  """
  cmd=path_fuselabels + " " + str(image_in) + " " + str(image_out) + " -lut " + str(path_lut) 
  if verbose:
    print cmd
  os.system(cmd)
  if U8:
    cmd=path_copy + " " + str(image_out) + " " + str(image_out) + " -o 1"
    if verbose:
      print cmd
    os.system(cmd)

  if not lazy:
    out=imread(str(image_out))
    if image_out != image_in:
      cmd="rm "+image_out 
      os.system(cmd)
    return out


def setVoxelSize(path_input, vx, vy, vz):
    ''' Set the voxel resolution of specified image
    '''
    os.system(path_setvoxelsize + ' ' + path_input + ' ' + str(vx) +\
              ' ' +str(vy) + ' ' + str(vz))
    
def connexe_with_options(path_input, path_output='tmp_threshold', path_seeds=None, sb=None, sh=None, tcc=None, ncc=None, connectivity=None, binary=False, label=False, size=False, sort=False, twoD=False, lazy=True):
    ''' Manual Connected Component with threshold
    path_input : path to the image to threshold
    path_output: path to image thresholded
    path_seeds : path to seeds image for seeded connected component extraction
                 implies binary output (unless changed afterwards)

    sb  : minimal threshold
    sh  : maximal threshold
    tcc : minimal size of connected components
    ncc : maximal number  of connected components
    connectivity : connectivity (4, 8, 6, 10, 18 or 26 (default))
    
    binary: if True, all valid connected components have the same value
    label : if True, one label per connected component

    size : the value attributed to points of a connected component is
          the size of the connected components (allows to select w.r.t to size afterwards)
    sort : the labels ordered the connected components with decreasing size
          implies -label-output (unless changed afterwards)

    twoD : slice by slice computation; each XY slice of a 3D volume is processed as an
    independent 2D image (eg labels restarted at 1 for each slice)

    lazy: do not return the output image if True
    '''

    options = ''
    if sb != None:
      option += ' -sb ' + str(sb)
    if sh != None:
      option += ' -sh ' + str(sh)
    if tcc != None:
      option += ' -tcc ' + str(tcc)
    if ncc != None:
      option += ' -ncc ' + str(ncc)
    if connectivity != None:
      option += ' -con ' + str(connectivity)
    if binary:
      option += ' -bin'
    if label:
      option += ' -label'
    if path_seeds:
      option += ' -seeds ' + str(path_seeds)
    if size:
      option += ' -size'
    if sort:
      option += ' -sort'
    if TwoD:
      option += ' -2D'

    os.system(path_connexe + ' ' + path_input + ' ' + path_output +\
              option)
    if not lazy:
        out = imread(path_output)
        os.system('rm ' + path_output)
        return out  

def compose_trsf(path_trsf_1, path_trsf_2, path_output="tmp_compose.inr", lazy=True, verbose=False):
  '''
  Transformations composition 
  '''
  assert(os.path.exists(path_trsf_1) and os.path.exists(path_trsf_2))
  cmd=path_compose_trsf + ' -res ' + path_output + ' -trsfs ' + path_trsf_1 + ' ' + path_trsf_2
  if verbose:
    print cmd
  os.system(cmd)
  if not lazy:
    out=imread(path_output)
    cmd='rm '+path_output
    os.system(cmd)
    return out

def _multiple_trsfs(format_in, format_out, first_index, last_index, reference_index, trsf_type='rigid', method='propagation', nfloatingbefore=None, nfloatingafter=None, verbose=False):
  '''
  Given the transformations between couple of images (typically successive
  images in a series), compute the transformations for every images with
  respect to the same reference image

  format_in    # format 'a la printf' of transformations to be processed
               # must contain two '%d'
               # the first one if for the floating image
               # the second one if for the reference image
  format_out # format 'a la printf' for output transformations
  first_index %d     # first value of the index in the format
  last_index %d      # last value of the index in the format
  reference_index %d # index of the reference image for result transformations
   # if 'before' or 'after' are non-null, the interval of transformation
   # for a given reference 'ref' is [ref-before, ref-1] U [ref+1, ref+after]
  trsf_type %s # transformation type
    translation2D, translation3D, translation-scaling2D, translation-scaling3D,
  rigid2D, rigid3D, rigid, similitude2D, similitude3D, similitude,
  affine2D, affine3D, affine, vectorfield2D, vectorfield3D, vectorfield, vector
  method %s # 
    average: compute transformation estimation by averaging composition
           of transformation
    propagation: compose transformations from the reference one
  nfloatingbefore %d # relative left half-interval for floating images
  nfloatingafter %d  # relative right half-interval for floating images
  '''
  cmd=path_multiple_trsfs+ \
      ' ' + format_in + ' -res ' + format_out + \
      ' -trsf-type ' + trsf_type + \
      ' -f ' + str(first_index) + ' -l ' + str(last_index) + \
      ' -ref ' +  str(reference_index) + ' -method ' + method 
  if nfloatingbefore != None:
    cmd = cmd + ' -nfloatingbefore ' + str(nfloatingbefore)
  if nfloatingafter != None:
    cmd = cmd + ' -nfloatingafter ' + str(nfloatingafter)

  if verbose:
    print cmd
  os.system(cmd)

def _change_multiple_trsfs(format_trsf_in,
                          format_trsf_out, 
                          format_image_in, template_image_out, 
                          first_index, last_index, 
                          reference_index=None, trsf_type='rigid', 
                          threshold=2, iso=None, margin=10, verbose=False):
  '''
  Given a list of transformation towards a reference,
  compute a new template that contains all transformed images as well
  as the new transformations

format_trsf_in # format 'a la printf' of transformation files
             # to be processed. It must contain one '%d'
             # depicts transformations of the form T_{i<-ref}
             # (ie allows to resample image I_i in geometry of I_ref)
format_trsf_out # format 'a la printf' for output transformations
             # will allow to resample input image into the resulting
             # template (thus still of the form T_{i<-ref})
             # reference is changed if '-index-reference' is used
format_image_in # format 'a la printf' of image corresponding to the input images
             # one template for all
template_image_out # output template image corresponding to the output transformations
first_index  # first value of the index in the format (%d)
last_index   # last value of the index in the format (%d)
reference_index # index of the reference transformation
             # the corresponding image will only be translated
             # if none, only translations will change in transformations
             # (ie reference is not changed) (%d)
trsf_type    # transformation type, which can be :
             # [translation2D, translation3D, translation-scaling2D, translation-scaling3D,
             # rigid2D, rigid3D, rigid, similitude2D, similitude3D, similitude,
             # affine2D, affine3D, affine, vectorfield2D, vectorfield3D, vectorfield, vector]
threshold    # threshold on input templates/images to compute the
             # useful bounding box (else it is the entire image)
iso          # make voxels isotropic for the output template (%f)
             # (if no value is given, uses the smallest voxel size from the template(s))
margin       # add a margin (in voxels)
verbose      # vorbose mode
  '''
  cmd=path_change_multiple_trsfs+ \
      ' -trsf-format ' + format_trsf_in + ' -res-trsf-format ' + format_trsf_out +\
      ' -trsf-type ' + trsf_type + \
      ' -f ' + str(first_index) + ' -l ' + str(last_index) + \
      ' -image-format ' +  str(format_image_in) + ' -result-template ' + str(template_image_out)
  if reference_index != None:
    cmd = cmd + ' -index-reference ' + str(reference_index)
  if threshold != None:
    cmd = cmd + ' -threshold ' + str(threshold)
  if iso != None:
    cmd = cmd + ' -res-iso ' + str(iso)
  if margin != None:
    cmd = cmd + ' -margin ' + str(margin)

  if verbose:
    print cmd
  os.system(cmd)

def obsolete_Arit(image_in, image_ext_or_out, image_out=None, Mode=None, Type='', lazy=True, verbose=False):
  """
  Addition, soustraction, multiplication, division, minimum, maximum d'images
  accepted modes : {'add', 'sub', 'mul', 'div', 'min', 'max'}
  accepted types :   '-o 1'    : unsigned char
                     '-o 2'    : unsigned short int
                     '-o 2 -s' : short int
                     '-o 4 -s' : int
                     '-r'      : float
                     ''      : default (image_in type)
  """
  mode=Mode.lstrip('-')
  existing_modes=['add', 'sub', 'mul', 'div', 'min', 'max']
  assert(existing_modes.count(mode)), "Mode option is wrong or missing, see help."
  existing_types=['-o 1', '-o 2', '-o 2 -s', '-o 4 -s', '-r', '']
  assert(existing_types.count(Type)), "Type option is wrong, see help"
  path_Arit = _find_exec('Arit')
  cmd=path_Arit + ' ' + str(image_in) + ' ' + str(image_ext_or_out)
  if image_out:
    cmd += ' ' + str(image_out)
  cmd += ' -' + mode + ' ' + Type

  if verbose:
    print cmd
  os.system(cmd)
  if not lazy:
    out=None
    if image_out:
      out=imread(str(image_out))
      cmd='rm '+image_out
    else:
      out=imread(str(image_ext_or_out))
      cmd='rm '+image_ext_or_out
    os.system(cmd)
    return out

def obsolete_Logic(image_in, image_ext_or_out, image_out=None, Mode=None, lazy=True, verbose=False):
  """
  Usage : Logic [-inv | [-et|-and] | [-ou|-or] | [-xou|-xor] | -mask]
  accepted modes : {'inv', 'et', 'and', 'ou', 'or', 'xou', 'xor', 'mask'}
  """
  mode=Mode.lstrip('-')
  existing_modes=['inv', 'et', 'and', 'ou', 'or', 'xou', 'xor', 'mask']
  assert(existing_modes.count(mode)), "Mode option is wrong or missing, see help."
  path_Logic = _find_exec('Logic')
  cmd=path_Logic + ' ' + str(image_in) + ' ' + str(image_ext_or_out)
  if image_out:
    cmd += ' ' + str(image_out)
  cmd += ' -' + mode 

  if verbose:
    print cmd
  os.system(cmd)
  if not lazy:
    out=None
    if image_out:
      out=imread(str(image_out))
      cmd='rm '+image_out
    else:
      out=imread(str(image_ext_or_out))
      cmd='rm '+image_ext_or_out
    os.system(cmd)
    return out

def obsolete_createImage(image_out, image_template, options='', lazy=True, verbose=False):
  """
  Creation d'image vide a partir de template
  accepted types :   '-o 1'    : unsigned char
                     '-o 2'    : unsigned short int
                     '-o 2 -s' : short int
                     '-o 4 -s' : int
                     '-r'      : float
                     ''      : default (image_template type)
  other options :    [-dim %d %d [%d] | [-x %d] [-y %d] [-z %d]] [-v %d] [-template %s]
                     [-voxel | -pixel | -vs %f %f [%f] | [-vx %f] [-vy %f] [-vz %f] ]
                     [-value %lf]
  """
  path_create_image = _find_exec('createImage')
  cmd=path_create_image + ' ' + str(image_out) + ' -template ' + str(image_template) + ' ' + options

  if verbose:
    print cmd
  os.system(cmd)
  if not lazy:
    out=None
    if image_out:
      out=imread(str(image_out))
      cmd='rm '+image_out
    else:
      out=imread(str(image_ext_or_out))
      cmd='rm '+image_ext_or_out
    os.system(cmd)
    return out

def setVoxelValue(image_in, image_out, x, y, z, value, lazy=True, verbose=False):
  """
  Fonction permettant d'attribuer au voxel de coordonnees (x,y,z) d'une image une valeur specifique.
  """
  cmd=path_setvoxelvalue + " " + str(image_in) + " " + str(image_out) + " -x " + str(x) + " -y " + str(y) + " -z " + str(z) + " -i " + str(value)

  if verbose:
    print cmd
  os.system(cmd)
  if not lazy:
    out=None
    out=imread(str(image_out))
    if image_out != image_in:
      cmd="rm "+image_out 
      os.system(cmd)
    return out

def setLabelValue(image_in, image_out, label_in, label_out, lazy=True, verbose=False):
  """
  Fonction permettant d'attribuer a un label d'image label_in une nouvelle valeur label_out.
  """
  cmd=path_fuselabels + " " + str(image_in) + " " + str(image_out) + " -p " + str(label_out) + " " + str(label_in)
  if verbose:
    print cmd
  os.system(cmd)
  if not lazy:
    out=imread(str(image_out))
    if image_out != image_in:
      cmd="rm "+image_out 
      os.system(cmd)
    return out

def resetLabels(image_in, image_out, labels, lazy=True, verbose=False):
  """
  Mise a zero des labels specifies.
  L'argument labels peut etre de type int ou list 
  """
  cmd=path_fuselabels + " " + str(image_in) + " " + str(image_out) + " -s " + str(labels).strip('[]').replace(',',' ')
  if verbose:
    print cmd
  os.system(cmd)
  if not lazy:
    out=imread(str(image_out))
    if image_out != image_in:
      cmd="rm "+image_out 
      os.system(cmd)
    return out

def erodeLabels(image_in, image_out, r=None, lazy=True, verbose=False):
  """
  Fonction permettant d'eroder tous les labels non nuls d'une image de labels.
  Possibilite de specifier un rayon d'erosion r. Sinon, l'erosion se limite a la mise a zero des voxels en bordure de labels.  
  """
  image_tmp=image_out
  flag_rm=False
  if image_tmp==image_in:
    if image_tmp.count(os.path.sep):
      image_tmp=os.path.sep.join(image_tmp.split(os.path.sep)[:-1])+os.path.sep+"erodeLabels_tmp.inr"
    else:
      image_tmp="erodeLabels_tmp.inr"
    flag_rm=True
  assert(image_in != image_tmp)
  cmd=path_labelborders + " " + str(image_in) + ' ' + str(image_tmp)
  if verbose:
    print cmd
  os.system(cmd)
  if r:
    cmd=path_morpho + ' -dil -r ' + str(r) + ' ' + str(image_tmp) + ' ' + str(image_tmp)
    if verbose:
      print cmd
    os.system(cmd)
  cmd = path_Logic + ' -inv ' + str(image_tmp) + ' ' + str(image_tmp)
  if verbose:
    print cmd
  os.system(cmd)
  cmd=path_Logic + " -mask " + str(image_tmp) + ' ' + str(image_in) + ' ' + str(image_out)
  if verbose:
    print cmd
  os.system(cmd)
  if flag_rm:
    cmd='rm -f '+image_tmp
    if verbose:
      print cmd
    os.system(cmd)
  if not lazy:
    out=imread(str(image_out))
    if image_out != image_in:
      cmd="rm "+image_out 
      os.system(cmd)
    return out

def labelBorders(image_in, image_out, lazy=True, verbose=False):
  """
  Fonction permettant d'eroder tous les labels non nuls d'une image de labels.
  Possibilite de specifier un rayon d'erosion r. Sinon, l'erosion se limite a la mise a zero des voxels en bordure de labels.  
  """
  assert(os.path.exists(image_in))
  cmd=path_labelborders + " " + str(image_in) + ' ' + str(image_out)
  if verbose:
    print cmd
  os.system(cmd)
  if not lazy:
    out=imread(str(image_out))
    if image_out != image_in:
      cmd="rm "+image_out 
      os.system(cmd)
    return out

#####
def obsolete_boudingboxes(image_labels, file_out=None, verbose=False):
  """
  Calcul des bounding-boxes de chaque label de l'image d'entree.
  Si file_out est renseigne, le resultat est sauvegarde dans ce fichier.
  Output : dictionnaire D dont les cles sont les labels d'image et les 
  valeurs sont les listes contenant dans l'ordre les informations de 
  [volume, xmin, ymin, zmin, xmax, ymax, zmax] correspondant a ces labels
  avec volume > 0 et label > 0.
  """
  assert(os.path.exists(image_labels))
  keep_file=True
  if not file_out:
    keep_file=False
    file_out='boundingboxes.txt'
  path_boundingboxes = _find_exec('boundingboxes')
  cmd=path_boundingboxes+' '+image_labels+' '+file_out
  if verbose:
    print cmd
  os.system(cmd)
  D={}
  f=open(file_out, 'r')
  lines=f.readlines()
  f.close()

  D={}
  for line in lines:
    if not line.lstrip().startswith('#'):
      l=line.split()
      if int(l[1]):
        D[int(l[0])]=map(int,l[1:])

  if not keep_file:
    cmd='rm '+file_out
    if verbose:
      print cmd
    os.system(cmd)

  return D

def obsolete_cropImage(image_in, image_out, bbox, verbose=False):
  """
  Croping an image.
  Inputs:
    - image_in : path to image in
    - bbox : [xmin, ymin, zmin, xmax, ymax, zmax]
  Outputs:
    - image_out : path to cropped image
  """
  assert(os.path.exists(image_in))
  bbox=bbox[-6:]
  path_cropImage = _find_exec('cropImage')
  cmd=path_cropImage + ' ' + image_in + " " + image_out + " -origin " + str(bbox[0]) + ' ' + str(bbox[1]) + ' ' + str(bbox[2]) + ' -dim ' + str(bbox[3]-bbox[0]+1) + " " + str(bbox[4]-bbox[1]+1) + " " + str(bbox[5]-bbox[2]+1)
  if verbose:
    print cmd
  os.system(cmd)

def obsolete_patchLogic(image_patch, image_ext, image_out, origin, Mode=None, verbose=False):
  '''
  Usage : Logic [[-et|-and] | [-ou|-or] | [-xou|-xor]]
  accepted modes : {'et', 'and', 'ou', 'or', 'xou', 'xor'}
  bbox : tuple of 3 elements (origin), 6 elements (bounding box with origin in the 3 first elements), or 7 elements (volume first, then bounding box) of the patch.
  '''
  mode=Mode.lstrip('-')
  existing_modes=['et', 'and', 'ou', 'or', 'xou', 'xor']
  assert(existing_modes.count(mode)), "Mode option is wrong or missing, see help."

  assert(os.path.exists(image_patch) and os.path.exists(image_ext))
  if len(origin)==7:
    origin=origin[-6:]
  path_patchLogic = _find_exec('patchLogic')
  cmd = path_patchLogic + ' ' + image_patch + " " + image_ext + " " + image_out + " -origin " + str(origin[0]) + ' ' + str(origin[1]) + ' ' + str(origin[2]) + ' -' + mode
  if verbose:
    print cmd
  os.system(cmd)

def obsolete_mc_adhocFuse(image_fuse, image_seg, image_out, min_percentile=0.01, max_percentile=0.99, min_method='cellinterior', max_method='cellborder', sigma=5.0, verbose=False):
  '''
  Function for fused images enhancement knowing a segmentation propagation
  Usage: mc-adhocFuse -intensity-image|-ii %s
   -reconstructed-image|-ri %s
   [-minimum-mask-image|-min-mi %s]
   [-maximum-mask-image|-max-mi %s]
   [-segmentation-image|-si %s]
   [-result-minimum-image|-rmini %s]
   [-result-maximum-image|-rmaxi %s]
   [-result-intensity-image|-rii %s]
   [-result-fused-image|-rfi %s]
   [-min-percentile|-min-p %f]
   [-max-percentile|-max-p %f]
   [-min-method|-method-min global|cell|cellborder|cellinterior|voxel]
   [-max-method|-method-max global|cell|cellborder|cellinterior|voxel]
   [-sigma %f]
   [-parallel|-no-parallel] [-max-chunks %d]
   [-parallelism-type|-parallel-type default|none|openmp|omp|pthread|thread]
   [-omp-scheduling|-omps default|static|dynamic-one|dynamic|guided]
   [-inv] [-swap] [output-image-type | -type s8|u8|s16|u16...]
   [-verbose|-v] [-no-verbose|-noverbose|-nv]
   [-debug|-D] [-no-debug|-nodebug]
   [-allow-pipe|-pipe] [-no-allow-pipe|-no-pipe|-nopipe]
   [-print-parameters|-param]
   [-print-time|-time] [-no-time|-notime]
   [-trace-memory|-memory] [-no-memory|-nomemory]
   [-help|-h]
  --------------------------------------------------
    -intensity-image|-ii %s
    -reconstructed-image|-ri %s
    -minimum-mask-image|-min-mi %s
    -maximum-mask-image|-max-mi %s
    -segmentation-image|-si %s
    -result-minimum-image|-rmini %s
    -result-maximum-image|-rmaxi %s
    -result-intensity-image|-rii %s
    -result-fused-image|-rfi %s
  # ...
    -min-percentile|-min-p %f
    -max-percentile|-max-p %f
    -min-method|-method-min global|cell|cellborder|cellinterior|voxel
    -max-method|-method-max global|cell|cellborder|cellinterior|voxel
    -sigma %f
  # parallelism parameters
   -parallel|-no-parallel:
   -max-chunks %d:
   -parallelism-type|-parallel-type default|none|openmp|omp|pthread|thread:
   -omp-scheduling|-omps default|static|dynamic-one|dynamic|guided:
  # general image related parameters
    -inv: inverse 'image-in'
    -swap: swap 'image-in' (if encoded on 2 or 4 bytes)
     output-image-type: -o 1    : unsigned char
                        -o 2    : unsigned short int
                        -o 2 -s : short int
                        -o 4 -s : int
                        -r      : float
    -type s8|u8|s16|u16|... 
     default is type of input image
  # general parameters 
    -verbose|-v: increase verboseness
      parameters being read several time, use '-nv -v -v ...'
      to set the verboseness level
    -no-verbose|-noverbose|-nv: no verboseness at all
    -debug|-D: increase debug level
    -no-debug|-nodebug: no debug indication
    -allow-pipe|-pipe: allow the use of stdin/stdout (with '-')
    -no-allow-pipe|-no-pipe|-nopipe: do not allow the use of stdin/stdout
    -print-parameters|-param:
    -print-time|-time:
    -no-time|-notime:
    -trace-memory|-memory:
    -no-memory|-nomemory:
    -h: print option list
    -help: print option list + details
  --------------------------------------------------
  '''
  assert(os.path.exists(image_fuse))
  assert(os.path.exists(image_seg))
  path_mc_adhocfuse = _find_exec('mc-adhocFuse')
  arguments=[path_mc_adhocfuse, '-intensity-image ', image_fuse, '-segmentation-image', image_seg, \
  '-min-percentile', str(min_percentile), '-max-percentile', str(max_percentile), \
  '-min-method', min_method, '-max-method', max_method, '-sigma', str(sigma), '-result-intensity-image',\
   image_out]
  cmd=' '.join(arguments)
  if verbose:
    print cmd
  os.system(cmd)
