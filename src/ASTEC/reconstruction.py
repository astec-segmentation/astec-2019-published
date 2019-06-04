
import os

import ACE
import commonTools
import nomenclature
import CommunFunctions.cpp_wrapping as cpp_wrapping

monitoring = commonTools.Monitoring()


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

def get_deformation_from_current_to_previous(current_time, environment, parameters, previous_time):
    """

    :param current_time:
    :param environment:
    :param parameters:
    :param previous_time:
    :return:
    """

    proc = 'get_deformation_from_current_to_previous'

    #
    # image to be registered
    # floating image = image at previous time
    # reference image = image at current time
    #

    current_name = nomenclature.replaceTIME(environment.path_fuse_exp_files, current_time)
    current_image = commonTools.find_file(environment.path_fuse_exp, current_name, monitoring=None, verbose=False)

    if current_image is None:
        monitoring.to_log_and_console("    .. " + proc + " no fused image was found for time "
                                      + str(current_time), 2)
        return None

    current_image = os.path.join(environment.path_fuse_exp, current_image)

    previous_name = nomenclature.replaceTIME(environment.path_fuse_exp_files, previous_time)
    previous_image = commonTools.find_file(environment.path_fuse_exp, previous_name, monitoring=None, verbose=False)

    if previous_image is None:
        monitoring.to_log_and_console("    .. " + proc + " no fused image was found for time "
                                      + str(previous_time), 2)
        return None

    previous_image = os.path.join(environment.path_fuse_exp, previous_image)

    #
    # temporary files for registration
    #

    affine_image = commonTools.add_suffix(previous_image, "_affine",
                                          new_dirname=environment.temporary_path,
                                          new_extension=parameters.default_image_suffix)

    vector_image = commonTools.add_suffix(previous_image, "_vector", new_dirname=environment.temporary_path,
                                          new_extension=parameters.default_image_suffix)

    affine_trsf = commonTools.add_suffix(previous_image, "_affine",
                                         new_dirname=environment.temporary_path,
                                         new_extension="trsf")

    vector_trsf = commonTools.add_suffix(previous_image, "_vector",
                                         new_dirname=environment.temporary_path,
                                         new_extension="trsf")

    #
    #
    #

    if os.path.isfile(vector_trsf):
        return vector_trsf

    #
    # compute non-linear transformation
    #
    monitoring.to_log_and_console("    .. non-linear registration of time point " + str(previous_time) + " on "
                                  + str(current_time), 2)
    cpp_wrapping.non_linear_registration(current_image, previous_image, affine_image, vector_image,
                                         affine_trsf, vector_trsf, monitoring=monitoring)

    if os.path.isfile(vector_trsf):
        return vector_trsf

    return None


#
#
#
#
#


def get_segmentation_image(current_time, environment):
    """

    :param current_time:
    :param environment:
    :return:
    """

    proc = 'get_segmentation_image'

    if current_time is None:
        monitoring.to_log_and_console("    .. " + proc + ": no valid time: "
                                      + str(current_time), 2)
        return None

    seg_name = nomenclature.replaceTIME(environment.path_seg_exp_files, current_time)
    seg_image = commonTools.find_file(environment.path_seg_exp, seg_name, monitoring=None, verbose=False)

    if seg_image is None:
        seg_name = nomenclature.replaceTIME(environment.path_mars_exp_files, current_time)
        seg_image = commonTools.find_file(environment.path_seg_exp, seg_name, monitoring=None, verbose=False)

    if seg_image is None:
        monitoring.to_log_and_console("    .. " + proc + ": no segmentation image was found for time "
                                      + str(current_time), 2)
        monitoring.to_log_and_console("       \t  " + "was looking for file '" + str(seg_name) + "'", 2)
        monitoring.to_log_and_console("       \t  " + "in directory '" + str(environment.path_seg_exp) + "'", 2)
        return None

    return os.path.join(environment.path_seg_exp, seg_image)


#
#
#
#
#


def get_previous_deformed_segmentation(current_time, environment, parameters, previous_time=None):
    """

    :param current_time:
    :param environment:
    :param parameters:
    :param previous_time:
    :return:
    """

    #
    #
    #

    proc = 'get_previous_deformed_segmentation'

    #
    # it will generate (in the temporary_path directory)
    #
    # files related to the deformation computation
    # - $EN_fuse_t$TIME_affine.inr
    # - $EN_fuse_t$TIME_affine.trsf
    # - $EN_fuse_t$TIME_vector.inr
    # - $EN_fuse_t$TIME_vector.trsf
    #
    # the deformed segmentation
    # - $EN_[mars,seg]_t$TIME_deformed.inr
    #

    prev_segimage = get_segmentation_image(previous_time, environment)
    if prev_segimage is None:
        monitoring.to_log_and_console("    .. " + proc + ": no segmentation image was found for time "
                                      + str(previous_time), 2)
        return None

    prev_def_segimage = commonTools.add_suffix(prev_segimage, "_deformed", new_dirname=environment.temporary_path,
                                               new_extension=parameters.default_image_suffix)

    if os.path.isfile(prev_def_segimage):
        return prev_def_segimage

    deformation = get_deformation_from_current_to_previous(current_time, environment, parameters, previous_time)

    if deformation is None:
        monitoring.to_log_and_console("    .. " + proc + ": no deformation was found for time "
                                      + str(current_time) + " towards " + str(previous_time), 2)
        return None

    monitoring.to_log_and_console("    .. resampling of '" + str(prev_segimage).split(os.path.sep)[-1] + "'", 2)

    cpp_wrapping.apply_transformation(prev_segimage, prev_def_segimage, deformation,
                                      interpolation_mode='nearest', monitoring=monitoring)

    return prev_def_segimage


########################################################################################
#
#
#
########################################################################################


def build_membrane_image(current_time, environment, parameters, previous_time=None):
    """

    :param current_time:
    :param environment:
    :param parameters:
    :param previous_time:
    :return:
    """
    #
    # build input image(s) for segmentation
    # 0. build names
    # 1. do image enhancement if required
    # 2. transform input image if required
    # 3. mix the two results
    #
    # If any transformation is performed on the input image, the final result image (input of the watershed)
    # will be named (in the 'temporary_path' directory)
    # - 'input_image'_membrane
    # if fusion is required
    # - 'input_image'_enhanced
    #   is the membrance-enhanced image (by tensor voting method)
    # - 'input_image'_intensity
    #   is the intensity-transformed image (by normalization)
    #

    input_name = nomenclature.replaceTIME(environment.path_fuse_exp_files, current_time)
    input_image = commonTools.find_file(environment.path_fuse_exp, input_name, monitoring)

    if input_image is None:
        monitoring.to_log_and_console("    .. image '" + input_name + "' not found in '"
                                      + str(environment.path_fuse_exp) + "'", 2)
        return None
    input_image = os.path.join(environment.path_fuse_exp, input_image)

    enhanced_image = None
    intensity_image = None
    membrane_image = None

    if parameters.intensity_enhancement is None or parameters.intensity_enhancement.lower() == 'none':

        #
        # no membrane enhancement image
        # only set intensity_image
        #
        if parameters.intensity_transformation is None or parameters.intensity_transformation.lower() == 'none' \
                or parameters.intensity_transformation.lower() == 'identity':
            #
            # nothing to do
            #
            intensity_image = input_image
        elif parameters.intensity_transformation.lower() == 'normalization_to_u8' \
                or parameters.intensity_transformation.lower() == 'global_normalization_to_u8' \
                or parameters.intensity_transformation.lower() == 'cell_normalization_to_u8':
            intensity_image = commonTools.add_suffix(input_image, "_membrane",
                                                     new_dirname=environment.path_reconstruction,
                                                     new_extension=parameters.default_image_suffix)
        else:
            monitoring.to_log_and_console("    unknown intensity transformation method: '"
                                          + str(parameters.intensity_transformation) + "'", 2)

    elif parameters.intensity_enhancement.lower() == 'gace' or parameters.intensity_enhancement.lower() == 'glace':
        #
        # only set enhanced_image
        # or set the 3 names: intensity_image, enhanced_image, and membrane_image
        #
        if parameters.intensity_transformation is None or parameters.intensity_transformation.lower() == 'none':
            enhanced_image = commonTools.add_suffix(input_image, "_membrane",
                                                    new_dirname=environment.path_reconstruction,
                                                    new_extension=parameters.default_image_suffix)
        elif parameters.intensity_transformation.lower() == 'identity':
            intensity_image = input_image
            enhanced_image = commonTools.add_suffix(input_image, "_enhanced", new_dirname=environment.temporary_path,
                                                    new_extension=parameters.default_image_suffix)
            membrane_image = commonTools.add_suffix(input_image, "_membrane",
                                                    new_dirname=environment.path_reconstruction,
                                                    new_extension=parameters.default_image_suffix)
        elif parameters.intensity_transformation.lower() == 'normalization_to_u8' \
                or parameters.intensity_transformation.lower() == 'global_normalization_to_u8' \
                or parameters.intensity_transformation.lower() == 'cell_normalization_to_u8':
            intensity_image = commonTools.add_suffix(input_image, "_intensity", new_dirname=environment.temporary_path,
                                                     new_extension=parameters.default_image_suffix)
            enhanced_image = commonTools.add_suffix(input_image, "_enhanced", new_dirname=environment.temporary_path,
                                                    new_extension=parameters.default_image_suffix)
            membrane_image = commonTools.add_suffix(input_image, "_membrane",
                                                    new_dirname=environment.path_reconstruction,
                                                    new_extension=parameters.default_image_suffix)
        else:
            monitoring.to_log_and_console("    unknown intensity transformation method: '"
                                          + str(parameters.intensity_transformation) + "'", 2)
            return None
    else:
        monitoring.to_log_and_console("    unknown membrane enhancement method: '"
                                      + str(parameters.intensity_enhancement) + "'", 2)
        return None

    #
    # computation of membrane enhancement image, if required
    #

    if parameters.intensity_enhancement is None or parameters.intensity_enhancement.lower() == 'none':
        pass
    elif parameters.intensity_enhancement.lower() == 'gace':
        if (not os.path.isfile(enhanced_image) and (membrane_image is None or not os.path.isfile(membrane_image))) \
                or monitoring.forceResultsToBeBuilt is True:
            monitoring.to_log_and_console("    .. membrane global enhancement of '"
                                          + str(input_image).split(os.path.sep)[-1] + "'", 2)
            ACE.monitoring.copy(monitoring)
            ACE.global_membrane_enhancement(input_image, enhanced_image, temporary_path=environment.temporary_path,
                                            parameters=parameters.ace)
    elif parameters.intensity_enhancement.lower() == 'glace':
        if (not os.path.isfile(enhanced_image) and (membrane_image is None or not os.path.isfile(membrane_image))) \
                or monitoring.forceResultsToBeBuilt is True:
            monitoring.to_log_and_console("    .. membrane cell enhancement of '"
                                          + str(input_image).split(os.path.sep)[-1] + "'", 2)
            ACE.monitoring.copy(monitoring)
            previous_deformed_segmentation = get_previous_deformed_segmentation(current_time, environment, parameters,
                                                                                previous_time)
            ACE.cell_membrane_enhancement(input_image, previous_deformed_segmentation, enhanced_image,
                                          temporary_path=environment.temporary_path,
                                          parameters=parameters.ace)
    else:
        monitoring.to_log_and_console("    unknown membrane enhancement method: '"
                                      + str(parameters.intensity_enhancement) + "'", 2)
        return None

    #
    # computation of transformed intensity image if required
    #

    arit_options = "-o 2"

    if parameters.intensity_transformation is None or parameters.intensity_transformation.lower() == 'none':
        pass
    elif parameters.intensity_transformation.lower() == 'identity':
        pass
    elif parameters.intensity_transformation.lower() == 'normalization_to_u8' \
            or parameters.intensity_transformation.lower() == 'global_normalization_to_u8':
        if (not os.path.isfile(intensity_image) and (membrane_image is None or not os.path.isfile(membrane_image))) \
                or monitoring.forceResultsToBeBuilt is True:
            monitoring.to_log_and_console("    .. intensity global normalization of '"
                                          + str(input_image).split(os.path.sep)[-1] + "'", 2)
            cpp_wrapping.global_normalization_to_u8(input_image, intensity_image,
                                                    min_percentile=0.01, max_percentile=0.99,
                                                    other_options=None, monitoring=monitoring)
        arit_options = "-o 1"
    elif parameters.intensity_transformation.lower() == 'cell_normalization_to_u8':
        if (not os.path.isfile(intensity_image) and (membrane_image is None or not os.path.isfile(membrane_image))) \
                or monitoring.forceResultsToBeBuilt is True:
            monitoring.to_log_and_console("    .. intensity cell-based normalization of '"
                                          + str(input_image).split(os.path.sep)[-1] + "'", 2)
            if previous_time is None:
                monitoring.to_log_and_console("       previous time point was not given", 2)
                return None
            previous_deformed_segmentation = get_previous_deformed_segmentation(current_time, environment, parameters,
                                                                                previous_time)
            cpp_wrapping.cell_normalization_to_u8(input_image, previous_deformed_segmentation, intensity_image,
                                                  min_percentile=0.01, max_percentile=0.99,
                                                  cell_normalization_min_method=parameters.cell_normalization_min_method,
                                                  cell_normalization_max_method=parameters.cell_normalization_max_method,
                                                  other_options=None, monitoring=monitoring)
        arit_options = "-o 1"
    else:
        monitoring.to_log_and_console("    unknown intensity transformation method: '"
                                      + str(parameters.intensity_transformation) + "'", 2)
        return None

    #
    # fusion of the two images (if fusion is required)
    #

    if enhanced_image is None:
        return intensity_image
    else:
        if intensity_image is None:
            return enhanced_image
        else:
            #
            # mix images
            #
            arit_options += " -max"
            if not os.path.isfile(membrane_image) or monitoring.forceResultsToBeBuilt is True:
                monitoring.to_log_and_console("       fusion of intensity and enhancement", 2)
                cpp_wrapping.arithmetic_operation(intensity_image, enhanced_image, membrane_image,
                                                  other_options=arit_options)
            return membrane_image

    #
    # should not reach this point
    #
