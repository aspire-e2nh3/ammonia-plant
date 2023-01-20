'''
This file will handle all interaction with the options configuration file.
'''
from genericpath import exists
import os
import configparser

import numpy as np

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
        self.plant_max_mol                = plant.getfloat('max mol')
        self.plant_max_iter               = plant.getfloat('max iter')
        self.plant_recycle_estimate       = plant.getfloat('recycle estimate')


        precooler = config["precooler"]
        self.precooler_water_mfr    = precooler.getfloat('water mass flow rate')
        self.precooler_T_cold_in    = precooler.getfloat('T water in')
        self.precooler_T_outlet      = precooler.getfloat('T precooler outlet')

        self.n2compressor_dT        = config.getfloat('n2 compressor', 'dT')

        self.h2compressor_dT        = config.getfloat('h2 compressor', 'dT')

        self.recompressor_eta       = config.getfloat('recompressor', 'eta')

        reactor = config["reactor"]

        self.reactor_T_1c           = reactor.getfloat('T_1c')
        self.reactor_length         = reactor.getfloat('length')
        self.reactor_diameter       = reactor.getfloat('tube diameter')
        self.reactor_num_tubes      = reactor.getint('number of tubes')|1


        self.reactor_r              = self.reactor_diameter/2
        self.reactor_cs_area        = np.pi * self.reactor_r ** 2  # m^2
        self.reactor_in_circum      = np.pi * self.reactor_r * 2

        min = reactor.getfloat('minimum_step')
        log = reactor.getint('log divs')
        split = reactor.getfloat('split point')*self.reactor_length
        lin = reactor.getint('lin divs')
        # generate vect

        self.reactor_vect = [0] + np.geomspace(min, split, log, endpoint=False).tolist() + \
                    np.linspace(split, self.reactor_length, lin-1).tolist()

        self.reactor_vectlen = len(self.reactor_vect)

        self.reactor_strength       = reactor.getfloat('wall strength') #MPa
        self.reactor_sf             = reactor.getint('safety factor')
        self.reactor_shell_density  = reactor.getfloat('shell density')
        self.reactor_cat_blk_density = reactor.getfloat('catalyst bulk density')
        self.reactor_void_frac      = reactor.getfloat('catalyst void fraction')
        self.reactor_cat_density    = reactor.getfloat('catalyst bulk density')*(1-self.reactor_void_frac)
        self.reactor_cat_size       = reactor.getfloat('catalyst particle size')/1000 # convert to m
        self.reactor_stress         = self.reactor_strength * 10 ** 6 / self.reactor_sf
        self.reactor_thickness      = self.plant_pressure * 10 ** 5 * self.reactor_r / self.reactor_stress
        self.reactor_OD             = self.reactor_diameter + self.reactor_thickness * 2
        self.reactor_out_circum     = np.pi * self.reactor_OD
        self.reactor_shell_volume   = self.reactor_thickness * np.pi * (self.reactor_length * 2 * self.reactor_r + self.reactor_r ** 2 * 4 / 3)
        self.reactor_cat_volume     = self.reactor_r ** 2 * self.reactor_length * np.pi
        self.reactor_shell_mass     = self.reactor_shell_volume * self.reactor_shell_density
        self.reactor_cat_mass       = self.reactor_cat_volume * self.reactor_cat_density


        self.reactor_shell_kval = reactor.getfloat('shell thermal conductivity')
        self.reactor_insul_kval = reactor.getfloat('insulation thermal conductivity')
        self.reactor_insul_thickness = reactor.getfloat('insulation thickness')
        self.reactor_external_htc = reactor.getfloat('external heat transfer coefficient')
        self.reactor_ext_T = reactor.getfloat('surrounding T')

        self.reactor_coolant_mfr = reactor.getfloat('coolant mass flow rate')
        self.reactor_coolant_channel = reactor.getfloat('coolant channel width')
        self.reactor_coolant_rho = reactor.getfloat('coolant density')
        self.reactor_coolant_mu = reactor.getfloat('coolant viscosity')
        self.reactor_coolant_cp = reactor.getfloat('coolant heat capacity')
        self.reactor_coolant_k = reactor.getfloat('coolant thermal conductivity')
        self.reactor_coolant_T = reactor.getfloat('coolant inlet T')


        
        he = config["heat exchanger"]
        self.he_concentric = 'concentric' in he.get('type')
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

        c = config["condenser"]
        self.c_shell = 'shell' in c.get('type')
        self.c_r1 = c.getfloat('r1')  # inner rad of inner ammonia pipe [m]
        self.c_r2 = c.getfloat('r2')  # outer rad of inner pipe [m]
        self.c_r3 = c.getfloat('r3')  # inner rad of outer coolant pipe [m]
        self.c_r4 = c.getfloat('r4')
        self.c_length = c.getfloat('length')
        self.c_numb = c.getint('numb')
        if self.c_shell:
            self.c_hyd = 2 * (self.c_r3**2 - self.c_numb * self.c_r2**2) / (self.c_r3 + self.c_numb * self.c_r2)
            if (self.c_r3**2) < (self.c_numb * self.c_r2**2):
                print('ERROR: Condenser flow area negative. Resize condenser')
            elif (self.c_r3**2) < (self.c_numb * self.c_r2**2)/0.907:
                print('ERROR: Condenser packing better than perfect hexagonal. Resize condenser')
            elif (self.c_r3**2) < (self.c_numb * self.c_r2**2)/0.82:
                print('CAUTION: Condenser packing better than random packing. Ensure condenser sizing is correct')
        else:
            self.c_hyd = 2 * (self.c_r3 - self.c_r2)  # hydraulic diameter of cooling channel [m]
        self.c_ix = c.getint('ix')  # number of elements along heat exchanger
        self.c_jints = c.getint('jints')  # max number of iterations
        self.c_kval = c.getfloat('kval')  # thermal conductivity of pipe [W/mK]
        self.c_kints = c.getint('kints')
        self.c_eval = c.getfloat('eval')

        self.c_dhvap = c.getfloat('dhvap')  # average value for enthalpy of condensation [J/mol]
        self.c_mcool = c.getfloat('mcool')  # mass flow of coolant [kg/s]
        self.c_cpcool = c.getfloat('cpcool')   # heat capacity of coolant [kg/s]
        self.c_Tcoolin = c.getfloat('Tcoolin')  # coolant inlet temp [K]
        self.c_rcool = c.getfloat('rcool')  # density of coolant [kg/m3]
        self.c_kcool = c.getfloat('kcool')  # thermal conductivity of coolant [W/mK]
        self.c_vcool = c.getfloat('vcool')  # dynamic viscosity of coolant [Pas]
        if self.c_shell:
            self.c_velcool = self.c_mcool /self.c_rcool/ (np.pi * (self.c_r3 ** 2 - self.c_numb * self.c_r2 ** 2))
            self.c_reycool = self.c_rcool * self.c_velcool * self.c_hyd / self.c_vcool
            self.c_prcool = self.c_vcool * self.c_cpcool / self.c_kcool
            self.c_nuscool = 0.023 * self.c_reycool ** 0.8 * self.c_prcool ** 0.4
            self.c_htc2 = self.c_nuscool * self.c_kcool / self.c_hyd
        else:
            self.c_velcool = self.c_mcool/self.c_rcool / (np.pi * (self.c_r3 ** 2 - self.c_r2 ** 2)) / self.c_numb
            self.c_reycool = self.c_rcool * self.c_velcool * self.c_hyd / self.c_vcool
            self.c_prcool = self.c_vcool * self.c_cpcool / self.c_kcool
            self.c_nuscool = 0.023 * self.c_reycool ** 0.8 * self.c_prcool ** 0.4
            self.c_htc2 = self.c_nuscool * self.c_kcool / self.c_hyd
        if self.c_velcool > 100:
            print('CAUTION: Condenser coolant velocity over 100m/s')



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
