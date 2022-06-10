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

    VALID_SECTIONS = ['plant', 'n2 compressor', 'h2 compressor', 'precooler', 'recompressor', 'reactor',
                      'heat exchanger', 'condenser']

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
        self.plant_max_mol                = plant.getfloat('max_mol')
        self.plant_max_iter               = plant.getfloat('max_iter')

        precooler = config["precooler"]
        self.precooler_water_mfr    = precooler.getfloat('water_mfr')
        self.precooler_T_cold_in    = precooler.getfloat('T_cold_in')

        self.n2compressor_dT        = config.getfloat('n2 compressor', 'dT')

        self.h2compressor_dT        = config.getfloat('h2 compressor', 'dT')

        self.recompressor_eta       = config.getfloat('recompressor', 'eta')

        reactor = config["reactor"]
        self.reactor_T_IN           = reactor.getfloat('T_IN')
        self.reactor_T_1c           = reactor.getfloat('T_1c')
        self.reactor_length         = reactor.getfloat('length')
        self.reactor_diameter       = reactor.getfloat('diameter')
        self.reactor_minimum_step   = reactor.getfloat('minimum_step')

        
        he = config["heat exchanger"]
        self.he_r1 = he.getfloat('r1')
        self.he_r2 = he.getfloat('r2')
        self.he_r3 = he.getfloat('r3')
        self.he_r4 = he.getfloat('r4')
        self.he_r5 = he.getfloat('r5')
        self.he_hyd = 2 * (self.he_r3 - self.he_r2)  # hydraulic diameter of cooling channel [m]
        self.he_length = he.getfloat('length')  # length of heat exchanger [m]
        self.he_numb = he.getint('numb')  # number of counterflow heat exchangers

        self.he_ix = he.getint('elements')  # number of elements along heat exchanger
        self.he_max_relax = he.getfloat('max relax')
        self.he_eval = he.getfloat('convergence')
        self.he_jints = he.getint('max iter')  # max number of iterations

        self.he_kval = he.getfloat('kval')  # thermal conductivity of pipe [W/mK]
        self.he_htc_ext = he.getfloat('htc ext')
        self.he_kval_insul = he.getfloat('kval insul')
        self.he_T_ext = he.getfloat('T ext')
'''
        c = config["condenser"]

        self.c_r1 = 0.002  # inner rad of inner ammonia pipe [m]
        self.c_r2 = 0.003  # outer rad of inner pipe [m]
        self.c_r3 = 0.02  # inner rad of outer coolant pipe [m]
        self.c_hyd = 2 * (self.c_r3 - self.c_r2)  # hydraulic diameter of cooling channel [m]
        self.c_length = 5  # length of heat exchanger [m]
        self.c_numb = 100  # number of counterflow heat exchangers
        
        self.c_ix = 500  # number of elements along heat exchanger
        self.c_jints = 100  # max number of iterations
        self.c_kval = 16  # thermal conductivity of pipe [W/mK]
        self.kints = 10000
        self.eval = 1e-5

        self.dhvap = 23370  # average value for enthalpy of condensation [J/mol]
        self.mcool = 0.3  # mass flow of coolant [kg/s]
        self.cpcool = 4186  # heat capacity of coolant [kg/s]
        self.Tcoolin = 273  # coolant inlet temp [K]
        self.rcool = 1000  # density of coolant [kg/m3]
        self.kcool = 0.6  # thermal conductivity of coolant [W/mK]
        self.vcool = 0.001  # dynamic viscosity of coolant [Pas]
        self.velcool = self.mcool / (np.pi * (self.r3 ** 2 - self.r2 ** 2)) / self.numb
        self.reycool = self.rcool * self.velcool * self.hyd / self.vcool
        self.prcool = self.vcool * self.cpcool / self.kcool
        self.nuscool = 0.023 * self.reycool ** 0.8 * self.prcool ** 0.4
        self.htc2 = self.nuscool * self.kcool / self.hyd
        used:?
        self.abar = 239.69  # constant for sat curve of nh3
        self.bbar = 0.0964  # constant for sat curve of nh3

'''

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
