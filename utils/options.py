'''
This file will handle all interaction with the options configuration file.
'''
import os
import configparser

class SSConfig:
    """
    A config class to store and handle all options for the steady-state ammonia-plant
    """

    # DEFAULTS: to be driven by a static .ini file stored in utils

    VALID_SECTIONS = ['plant', 'n2compressor', 'h2compressor', 'precooler', 'recompressor', 'reactor',
                      'heater', 'heat_exchanger', 'chiller', 'condenser']

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

        self.plant_convergence     = config.getfloat('plant', 'convergence')
        self.plant_pressure        = config.getfloat('plant', 'pressure')
        self.plant_n2              = config.getfloat('plant', 'n2')
        self.plant_h2              = config.getfloat('plant', 'h2')

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

        self.heat_exchanger_eff     = config.getfloat('heat_exchanger', 'eff')

        self.chiller_water_mfr      = config.getfloat('chiller', 'water_mfr')
        self.chiller_eff            = config.getfloat('chiller', 'eff')
        self.chiller_T_cold_in      = config.getfloat('chiller', 'T_cold_in')

        self.condenser_water_mfr    = config.getfloat('condenser', 'water_mfr')
        self.condenser_eff          = config.getfloat('condenser', 'eff')
        self.condenser_T_cold_in    = config.getfloat('condenser', 'T_cold_in')

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
        self.TERMINAL_END_LOG_DETAIL = config.getboolean('TERMINAL', 'END_LOG_DETAIL')
