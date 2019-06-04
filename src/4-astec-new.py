#!/usr/bin/python2.7

import os
import time
from argparse import ArgumentParser

#
# local imports
# add ASTEC subdirectory
#


import ASTEC.commonTools as commonTools
import ASTEC.ASTECNEW as ASTEC
import ASTEC.nomenclature as nomenclature
from ASTEC.CommunFunctions.cpp_wrapping import path_to_vt


#
#
#
#
#


def _set_options(my_parser):
    """

    :param my_parser:
    :return:
    """
    proc = "_set_options"
    if not isinstance(my_parser, ArgumentParser):
        print proc + ": argument is not of type ArgumentParser"
        return
    #
    # common parameters
    #

    my_parser.add_argument('-p', '--parameters',
                           action='store', dest='parameterFile', const=None,
                           help='python file containing parameters definition')
    my_parser.add_argument('-e', '--embryo-rep',
                           action='store', dest='embryo_path', const=None,
                           help='path to the embryo data')

    #
    # control parameters
    #

    my_parser.add_argument('-k', '--keep-temporary-files',
                           action='store_const', dest='keepTemporaryFiles',
                           default=False, const=True,
                           help='keep temporary files')

    my_parser.add_argument('-f', '--force',
                           action='store_const', dest='forceResultsToBeBuilt',
                           default=False, const=True,
                           help='force building of results')

    my_parser.add_argument('-v', '--verbose',
                           action='count', dest='verbose', default=2,
                           help='incrementation of verboseness')
    my_parser.add_argument('-nv', '--no-verbose',
                           action='store_const', dest='verbose', const=0,
                           help='no verbose at all')
    my_parser.add_argument('-d', '--debug',
                           action='count', dest='debug', default=0,
                           help='incrementation of debug level')
    my_parser.add_argument('-nd', '--no-debug',
                           action='store_const', dest='debug', const=0,
                           help='no debug information')

    return


#
#
# main function
#
#


def main():
    #
    # initialization
    #
    start_time = time.localtime()
    monitoring = commonTools.Monitoring()
    experiment = commonTools.Experiment()
    parameters = ASTEC.AstecParameters()
    environment = ASTEC.AstecEnvironment()

    #
    # reading command line arguments
    #
    parser = ArgumentParser(description='Astec')
    _set_options(parser)
    args = parser.parse_args()

    monitoring.update_from_args(args)
    experiment.update_from_args(args)

    #
    # reading parameter files
    # and updating parameters
    #
    parameterFile = commonTools.get_parameter_file(args.parameterFile)
    environment.update_from_file(parameterFile, start_time)
    environment.path_history_file = nomenclature.replaceEXECUTABLE(environment.path_history_file, __file__)
    environment.path_log_file = nomenclature.replaceEXECUTABLE(environment.path_log_file, __file__)

    if not os.path.isdir(environment.path_logdir):
        os.makedirs(environment.path_logdir)

    experiment.update_from_file(parameterFile)
    parameters.update_from_file(parameterFile)

    #
    # make fusion directory and subdirectory if required
    # => allows to write log and history files
    #    and to copy parameter file
    #
    # for i in range(0, len(environment.channel)):
    #    if not os.path.isdir(environment.channel[i].path_fuse_exp):
    #        os.makedirs(environment.channel[i].path_fuse_exp)

    #
    # write history information in history file
    #
    commonTools.write_history_information(environment.path_history_file,
                                          experiment,
                                          parameterFile,
                                          start_time,
                                          os.path.dirname(__file__),
                                          path_to_vt())

    #
    # define log file
    # and write some information
    #
    monitoring.logfile = environment.path_log_file
    ASTEC.monitoring.copy(monitoring)

    experiment.write_parameters(monitoring.logfile)
    environment.write_parameters(monitoring.logfile)
    parameters.write_parameters(monitoring.logfile)

    #
    # copy parameter file
    #
    commonTools.copy_date_stamped_file(parameterFile, environment.path_logdir, start_time)

    #
    # processing
    #
    ASTEC.astec_control(experiment, environment, parameters)

    #
    # end of execution
    # write execution time in both log and history file
    #
    endtime = time.localtime()
    with open(environment.path_log_file, 'a') as logfile:
        logfile.write("\n")
        logfile.write('Total execution time = '+str(time.mktime(endtime)-time.mktime(start_time))+' sec\n')
        logfile.write("\n")

    with open(environment.path_history_file, 'a') as logfile:
        logfile.write('# Total execution time = '+str(time.mktime(endtime)-time.mktime(start_time))+' sec\n')
        logfile.write("\n\n")


#
#
# main call
#
#


if __name__ == '__main__':
    main()