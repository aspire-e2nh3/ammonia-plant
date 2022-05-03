'''
This file will handle all interaction with the options configuration file.
'''
from genericpath import exists
import os
import configparser

from numpy import empty

class SSConfig:
    """
    A config class to store and handle all options for the steady-state ammonia-plant
    """

    # DEFAULTS: to be driven by a static .ini file stored in utils

    VALID_SECTIONS = ['plant', 'n2compressor', 'h2compressor', 'precooler', 'recompressor', 'reactor',
                      'heater', 'heat_exchanger', 'chiller', 'condenser','param']

    #VALID = {'MINIMIZERS': VALID_MINIMIZERS,
    #         'FITTING': VALID_FITTING,
    #         'JACOBIAN': VALID_JACOBIAN,
    #         'PLOTTING': VALID_PLOTTING,
    #         'LOGGING': VALID_LOGGING}
    # At a later date we could add valid ranges


    def __init__(self, file_name=None):
        """
        Initialise the options from a given file.
        """
        # stores the file name to be used to reset options for multiple
        # problem groups
        self.stored_file_name = file_name
        self.error_message = []
        self._results_dir = ''
        config = configparser.ConfigParser(inline_comment_prefixes="#")

        for section in self.VALID_SECTIONS:
            config.add_section(section)

        if file_name is not None:
            config.read(file_name)

        plant = config["plant"]
        self.plant_convergence      = plant.getfloat('convergence')
        self.plant_pressure         = plant.getfloat('pressure')
        self.plant_h2               = plant.getfloat('h2')
        self.plant_ratio_n          = plant.getfloat('ratio_n')
        self.plant_n2               = self.plant_h2/self.plant_ratio_n

        plant = config["plant"]
        self.precooler_water_mfr    = config.getfloat('precooler', 'water_mfr')
        self.precooler_T_cold_in    = config.getfloat('precooler', 'T_cold_in')

        self.n2compressor_dT        = config.getfloat('n2compressor', 'dT')

        self.h2compressor_dT        = config.getfloat('h2compressor', 'dT')

        self.recompressor_eta       = config.getfloat('recompressor', 'eta')

        self.reactor_T_IN           = config.getfloat('reactor', 'T_IN')
        self.reactor_T_1c           = config.getfloat('reactor', 'T_1c')
        self.reactor_length         = config.getfloat('reactor', 'length')
        self.reactor_diameter       = config.getfloat('reactor', 'diameter')
        self.reactor_minimum_step   = config.getfloat('reactor', 'minimum_step')

        self.max_mol                = 5
        self.max_iter               = 200


class OutOps:
    """A configuration class to store output options."""

    VALID_SECTIONS = ['PLOTTING', 'WRITING', 'TERMINAL']

    def __init__(self, file_name=None):
        """
        Initialise the options from a given file.
        """
        # stores the file name to be used to reset options for multiple
        # problem groups
        self.stored_file_name = file_name
        self.error_message = []
        self._results_dir = ''
        config = configparser.ConfigParser(inline_comment_prefixes="#")

        for section in self.VALID_SECTIONS:
            config.add_section(section)

        if file_name is not None:
            config.read(file_name)

        self.REACTOR_BED = config.getboolean('PLOTTING', 'REACTOR_BED')
        self.POWER_PIE_CHART = config.getboolean('PLOTTING', 'POWER_PIE_CHART')

        self.ITERATION_HISTORY = config.getboolean('WRITING', 'ITERATION_HISTORY')
        self.STEADY_STATE_SOLUTION = config.getboolean('WRITING', 'STEADY_STATE_SOLUTION')

        self.TERMINAL_LOG = config.getboolean('TERMINAL', 'LOG')
        self.TERMINAL_LOG_SHORT = config.getboolean('TERMINAL', 'LOG_SHORT')
        self.TERMINAL_END_LOG_DETAIL = config.getboolean('TERMINAL', 'END_LOG_DETAIL')


def read_list_of_floats(s):
    """
    Utility function to allow ranges to be read by the config parser

    :param s: string to convert to a list
    :type s: string

    :return: list of floats
    :rtype: list
    """
    if s[0] == '[' and s[-1] == ']':
        lst = [float(item) for item in s[1:-1].split(",")]
    else:
        raise ValueError
    return lst
