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

    VALID_SECTIONS = ['reactor', 'n2compressor', 'h2compressor',
                      'recompressor', 'heater', 'mixer', 'heat_exchanger_counter',
                      'heat_exchanger_water2gas', 'condensor']

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

        self.reactor_T_0_0 = config.getfloat('reactor','T_0_0')
        self.reactor_T_1_0 = config.getfloat('reactor','T_1_0')
        self.reactor_P_0 = config.getfloat('reactor','P_0')
        self.reactor_length = config.getfloat('reactor','length')
        self.reactor_diameter = config.getfloat('reactor','diameter')
        self.reactor_mini = config.getfloat('reactor','mini')
        self.reactor_convergence = config.getfloat('reactor','convergence')
        self.n2 = config.getfloat('reactor','n2')
        self.h2 = config.getfloat('reactor','h2')
        self.reactor_dT_res = config.getfloat('reactor','dT_res')
        self.recycle = config.getfloat('reactor','recycle')


        self.n2compressor_dT = config.getfloat('n2compressor','dT')

        self.h2compressor_dT = config.getfloat('h2compressor','dT')

        self.recompressor_eta = config.getfloat('recompressor','eta')

        self.he_counter_eff = config.getfloat('he_counter','eff')

        self.he_water2gas_eff = config.getfloat('he_water2gas','eff')

        self.condensor_water_mfr = config.getfloat('condensor','water_mfr')
        self.condensor_eff = config.getfloat('condensor','eff')

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
