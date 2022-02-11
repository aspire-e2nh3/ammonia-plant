'''
This file will handle all interaction with the options configuration file.
'''

import os

class Options:
    """
    An options class to store and handle all options for fitbenchmarking
    """

    DEFAULTS = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                            'default_options.ini'))
    VALID_SECTIONS = ['MINIMIZERS', 'FITTING', 'JACOBIAN',
                      'PLOTTING', 'LOGGING']
    VALID_MINIMIZERS = \
        {'bumps': ['amoeba', 'lm-bumps', 'newton', 'de', 'mp'],
         'dfo': ['dfogn', 'dfols'],
         'gradient_free': ['HillClimbingOptimizer',
                           'RepulsingHillClimbingOptimizer',
                           'SimulatedAnnealingOptimizer',
                           'RandomSearchOptimizer',
                           'RandomRestartHillClimbingOptimizer',
                           'RandomAnnealingOptimizer',
                           'ParallelTemperingOptimizer',
                           'ParticleSwarmOptimizer',
                           'EvolutionStrategyOptimizer',
                           'BayesianOptimizer',
                           'TreeStructuredParzenEstimators',
                           'DecisionTreeOptimizer'],
         'gsl': ['lmsder', 'lmder', 'nmsimplex', 'nmsimplex2',
                 'conjugate_pr', 'conjugate_fr', 'vector_bfgs',
                 'vector_bfgs2', 'steepest_descent'],
         'levmar': ['levmar'],
         'mantid': ['BFGS',
                    'Conjugate gradient (Fletcher-Reeves imp.)',
                    'Conjugate gradient (Polak-Ribiere imp.)',
                    'Damped GaussNewton', 'Levenberg-Marquardt',
                    'Levenberg-MarquardtMD', 'Simplex',
                    'SteepestDescent', 'Trust Region'],
         'matlab': ['Nelder-Mead Simplex'],
         'matlab_curve': ['Levenberg-Marquardt', 'Trust-Region'],
         'matlab_opt': ['levenberg-marquardt', 'trust-region-reflective'],
         'matlab_stats': ['Levenberg-Marquardt'],
         'minuit': ['minuit'],
         'ralfit': ['gn', 'gn_reg', 'hybrid', 'hybrid_reg'],
         'scipy': ['Nelder-Mead', 'Powell', 'CG', 'BFGS',
                   'Newton-CG', 'L-BFGS-B', 'TNC', 'SLSQP'],
         'scipy_ls': ['lm-scipy', 'trf', 'dogbox'],
         'scipy_go': ['differential_evolution', 'shgo', 'dual_annealing']}
    VALID_FITTING = \
        {'algorithm_type': ['all', 'ls', 'deriv_free', 'general', 'simplex',
                            'trust_region', 'levenberg-marquardt',
                            'gauss_newton', 'bfgs', 'conjugate_gradient',
                            'steepest_descent', 'global_optimization'],
         'software': ['bumps', 'dfo', 'gradient_free', 'gsl', 'levmar',
                      'mantid', 'matlab', 'matlab_curve', 'matlab_opt',
                      'matlab_stats', 'minuit', 'ralfit', 'scipy',
                      'scipy_ls', 'scipy_go'],
         'jac_method': ['scipy', 'analytic', 'default', 'numdifftools'],
         'cost_func_type': ['nlls', 'weighted_nlls', 'hellinger_nlls',
                            'poisson']}
    VALID_JACOBIAN = \
        {'scipy': ['2-point', '3-point', 'cs'],
         'analytic': ['cutest'],
         'default': ['default'],
         'numdifftools': ['central',
                          'complex', 'multicomplex',
                          'forward', 'backward']}
    VALID_PLOTTING = \
        {'make_plots': [True, False],
         'comparison_mode': ['abs', 'rel', 'both'],
         'table_type': ['acc', 'runtime', 'compare', 'local_min'],
         'colour_map': plt.colormaps()}
    VALID_LOGGING = \
        {'level': ['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'ERROR',
                   'CRITICAL'],
         'append': [True, False],
         'external_output': ['debug', 'display', 'log_only']}

    VALID = {'MINIMIZERS': VALID_MINIMIZERS,
             'FITTING': VALID_FITTING,
             'JACOBIAN': VALID_JACOBIAN,
             'PLOTTING': VALID_PLOTTING,
             'LOGGING': VALID_LOGGING}

    DEFAULT_MINIMZERS = \
        {'bumps': ['amoeba', 'lm-bumps', 'newton', 'mp'],
         'dfo': ['dfogn', 'dfols'],
         'gradient_free': ['HillClimbingOptimizer',
                           'RepulsingHillClimbingOptimizer',
                           'SimulatedAnnealingOptimizer',
                           'RandomSearchOptimizer',
                           'RandomRestartHillClimbingOptimizer',
                           'RandomAnnealingOptimizer',
                           'ParallelTemperingOptimizer',
                           'ParticleSwarmOptimizer',
                           'EvolutionStrategyOptimizer'],
         'gsl': ['lmsder', 'lmder', 'nmsimplex', 'nmsimplex2',
                 'conjugate_pr', 'conjugate_fr', 'vector_bfgs',
                 'vector_bfgs2', 'steepest_descent'],
         'levmar': ['levmar'],
         'mantid': ['BFGS',
                    'Conjugate gradient (Fletcher-Reeves imp.)',
                    'Conjugate gradient (Polak-Ribiere imp.)',
                    'Damped GaussNewton', 'Levenberg-Marquardt',
                    'Levenberg-MarquardtMD', 'Simplex',
                    'SteepestDescent', 'Trust Region'],
         'matlab': ['Nelder-Mead Simplex'],
         'matlab_curve': ['Levenberg-Marquardt', 'Trust-Region'],
         'matlab_opt': ['levenberg-marquardt', 'trust-region-reflective'],
         'matlab_stats': ['Levenberg-Marquardt'],
         'minuit': ['minuit'],
         'ralfit': ['gn', 'gn_reg', 'hybrid', 'hybrid_reg'],
         'scipy': ['Nelder-Mead', 'Powell', 'CG', 'BFGS',
                   'Newton-CG', 'L-BFGS-B', 'TNC', 'SLSQP'],
         'scipy_ls': ['lm-scipy', 'trf', 'dogbox'],
         'scipy_go': ['differential_evolution', 'shgo', 'dual_annealing']}
    DEFAULT_FITTING = \
        {'num_runs': 5,
         'algorithm_type': ['all'],
         'software': ['scipy', 'scipy_ls'],
         'jac_method': ['scipy'],
         'cost_func_type': 'weighted_nlls'}
    DEFAULT_JACOBIAN = \
        {'analytic': ['cutest'],
         'scipy': ['2-point'],
         'default': ['default'],
         'numdifftools': ['central']}
    DEFAULT_PLOTTING = \
        {'make_plots': True,
         'colour_map': 'magma_r',
         'colour_ulim': 100,
         'cmap_range': [0.2, 0.8],
         'comparison_mode': 'both',
         'table_type': ['acc', 'runtime', 'compare', 'local_min'],
         'results_dir': 'fitbenchmarking_results'}
    DEFAULT_LOGGING = \
        {'file_name': 'fitbenchmarking.log',
         'append': False,
         'level': 'INFO',
         'external_output': 'log_only'}
    DEFAULTS = {'MINIMIZERS': DEFAULT_MINIMZERS,
                'FITTING': DEFAULT_FITTING,
                'JACOBIAN': DEFAULT_JACOBIAN,
                'PLOTTING': DEFAULT_PLOTTING,
                'LOGGING': DEFAULT_LOGGING}

    def __init__(self, file_name=None):
        """
        Initialise the options from a file if file is given.
        Priority is values in the file, failing that, values are taken from
        DEFAULTS (stored in ./default_options.ini)

        :param file_name: The options file to load
        :type file_name: str
        """
        # stores the file name to be used to reset options for multiple
        # problem groups
        self.stored_file_name = file_name
        self.error_message = []
        self._results_dir = ''
        config = configparser.ConfigParser(converters={'list': read_list,
                                                       'str': str,
                                                       'rng': read_range},
                                           allow_no_value=True)

        for section in self.VALID_SECTIONS:
            config.add_section(section)

        if file_name is not None:
            config.read(file_name)

            # Checks that the user defined sections are valid
            if config.sections() != self.VALID_SECTIONS:
                raise OptionsError(
                    "Invalid options sections set, {0}, the valid sections "
                    "are {1}".format(config.sections(), self.VALID_SECTIONS))
            config.sections()

            # Checks that the options within the sections are valid
            for key in self.VALID_SECTIONS:
                default_options_list = list(self.DEFAULTS[key].keys())
                user_options_list = [option[0] for option in config.items(key)]
                if not set(user_options_list) <= set(default_options_list):
                    raise OptionsError(
                        "Invalid options key set in the {2} Section: \n{0}, \n"
                        " the valid keys are: \n{1}".format(
                            user_options_list,
                            default_options_list,
                            key))

        minimizers = config['MINIMIZERS']
        self._minimizers = {}
        self.minimizer_alg_type = {}
        for key in self.VALID_FITTING["software"]:
            self._minimizers[key] = self.read_value(minimizers.getlist,
                                                    key)

        fitting = config['FITTING']
        self.num_runs = self.read_value(fitting.getint, 'num_runs')
        self.algorithm_type = self.read_value(
            fitting.getlist, 'algorithm_type')
        self.software = self.read_value(fitting.getlist, 'software')
        self.jac_method = self.read_value(fitting.getlist, 'jac_method')
        self.cost_func_type = self.read_value(fitting.getstr, 'cost_func_type')

        jacobian = config['JACOBIAN']
        self.num_method = {}
        for key in self.VALID_FITTING["jac_method"]:
            self.num_method[key] = self.read_value(jacobian.getlist,
                                                   key)

        plotting = config['PLOTTING']
        self.make_plots = self.read_value(plotting.getboolean, 'make_plots')
        self.colour_map = self.read_value(plotting.getstr, 'colour_map')
        self.colour_ulim = self.read_value(plotting.getfloat, 'colour_ulim')
        self.cmap_range = self.read_value(plotting.getrng, 'cmap_range')

        self.comparison_mode = self.read_value(plotting.getstr,
                                               'comparison_mode')
        self.table_type = self.read_value(plotting.getlist, 'table_type')
        self.results_dir = self.read_value(plotting.getstr, 'results_dir')

        logging = config['LOGGING']
        self.log_append = self.read_value(logging.getboolean, 'append')
        self.log_file = self.read_value(logging.getstr, 'file_name')
        self.log_level = self.read_value(logging.getstr, 'level')

        self.external_output = self.read_value(logging.getstr,
                                               'external_output')

        if self.error_message != []:
            raise OptionsError('\n'.join(self.error_message))