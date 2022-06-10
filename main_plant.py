import os

import pandas
import pandas as pd
from plant_functions import *
import matplotlib.pyplot as plt
import math
import numpy as np
import copy
import pyromat as pm
import argparse
from utils.options import SSConfig, OutOps
from tristan_condenser import *

parser = argparse.ArgumentParser()
parser.add_argument(
    "-c",
    "--configuration",
    help="specify the '*.ini' file path that configures the model")
parser.add_argument(
    "-o",
    "--output",
    help="specify the '*.ini' file path that is to define the outputs from the model")
args = parser.parse_args()


#  future would be to add output options file too

def get_configs(args):
    """Retrive config objects."""
    default_steady_state = 'default_ssconfig.ini'
    default_outputs = 'default_output_ops.ini'

    if args.configuration is None:
        fcfg = os.path.join(os.path.dirname(__file__),
                            'utils',
                            default_steady_state)
    else:
        fcfg = args.configuration  # pathname from input argument
    cfg = SSConfig(fcfg)  # initialise

    if args.output is None:
        fops = os.path.join(os.path.dirname(__file__),
                            'utils',
                            default_outputs)
    else:
        fops = args.output  # pathname from input argument
    ops = OutOps(fops)  # initialise

    return cfg, ops


def evaluate_loop(cfg, ops, id_run):
    """Solve a single loop.

    :param cfg: config parser object with steady state configuration
    :param ops: config parser object with output options configuration

    :return: stream_data - pandas dataframe of states throughout loop
             power_data - pandas dataframe of power consumption
    """

    if ops.TERMINAL_LOG:
        print("Log for run_%d" % id_run)

    # initial inputs - define bed structure, quench ratio
    bed1 = Bed(cfg.reactor_length,
               cfg.reactor_diameter,
               cfg.reactor_minimum_step,
               log_divs=19, lin_divs=100)
    shell_mass, cat_mass = bed1.shell_mass, bed1.cat_mass


    # set reactor inlet temp
    reactor_in_temp = cfg.reactor_T_1c

    # initialise iteration limits
    inlet_temp = cfg.reactor_T_IN
    HTHE_P = 27000 * cfg.plant_h2/0.3

    recycle_estimate = 6

    # initialise power consumption dictionary
    power_consumption = {}

    # make n2 and h2 input lines from psa and electrolysis functions
    [Pipe_N2_LP, power_consumption["PSA"]] = psa_estimate(cfg.plant_n2)  # MAYBE ADD AIR_IN?
    [Pipe_H2_LP, power_consumption["electrolysis"]] = electrolysis(cfg.plant_h2)

    # compress N2 and H2 lines

    [Pipe_N2_HP, power_consumption["N2_comp"], heatlost1] = compressor(Pipe_N2_LP,
                                                                             cfg.plant_pressure)
    [Pipe_H2_HP, power_consumption["H2_comp"], heatlost2] = compressor(Pipe_H2_LP,
                                                                             cfg.plant_pressure)

    # mix H2 and N2 lines
    Pipe_IN_uncooled = mixer(Pipe_H2_HP, Pipe_N2_HP)

    '''
    if ops.TERMINAL_LOG:
        print('reactor weight = %3.3f' % shell_mass)
        print('catalyst weight = %3.3f' % cat_mass)
        print('N2 Compressor Polytropic Index = %3.3f' % n1)
        print('H2 Compressor Polytropic Index = %3.3f' % n2)
        print('H2 LP Feed = %3.1f' % Pipe_H2_LP.T, 'K')
        print('N2 LP Feed = %3.1f' % Pipe_N2_LP.T, 'K')
        print('H2 HP Feed = %3.1f' % Pipe_H2_HP.T, 'K')
        print('N2 HP Feed = %3.1f' % Pipe_N2_HP.T, 'K')
        print('Mixed Feed pre-cooling = %3.1f' % Pipe_IN_uncooled.T, 'K')
    '''
    [Pipe_IN, power_consumption["inflow cooling"], effectiveness_precooler] = heat_exchanger_water2gas(Pipe_IN_uncooled,
                                                                                                       T_end=cfg.reactor_T_IN,
                                                                                                       T_cin=cfg.precooler_water_mfr)

    Pipe_RE = State(recycle_estimate * cfg.plant_h2,
                    recycle_estimate * cfg.plant_n2,
                    0.5 * cfg.plant_n2,
                    cfg.reactor_T_IN,
                    cfg.plant_pressure - 2)
    count = 0
    stop = 0
    he_det = Heat_Exchanger_Details()
    Pipe_purged = State(0,0,0,273,200)
    note = ''
    while (stop == 0):
        count += 1
        if count >= cfg.max_iter:
            stop = 1
            note = "did not converge"

        #mix in recycle stream
        Pipe_1a = mixer(Pipe_IN, Pipe_RE)

        # recompress recycled stream
        [Pipe_1b, power_consumption["recompressor"], _] = compressor(Pipe_1a,
                                                                     p_out=cfg.plant_pressure,
                                                                     eta=cfg.recompressor_eta)

        # Add heat from Low Temp Heat Exchanger to Pipe 1b to make Pipe 1c
        Pipe_1c = heater_power(Pipe_1b,HTHE_P)

        # heat/cool to 673K
        [Pipe_1d, power_consumption["heater"]] = heater(Pipe_1c, cfg.reactor_T_1c, no_cooling=False)

        #Pipe_1d = copy.copy(Pipe_1c)
        #power_consumption["heater"] = 0
        #print(Pipe_1d.T)
        # Initialise bed data


        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ BED 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        [Pipe_2a,exotherm_reac,heatloss_reac,Bed_data] = reactor(Pipe_1d,bed1)
        #print(Pipe_2a.T)
        Pipe_2a.p -= 5*Pipe_1d.mol_tot/cfg.max_mol

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 1 (Pipe 6 to Pipe 4) ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # estimate heat exchange variant

        [Pipe_2b, Pipe_1c_fake, HTHE_P_new, heat_ex_heatloss, he_det] = tristan_heat_exchanger(Pipe_2a, Pipe_1b, he_det)
        effectiveness_heatex = he_det.last_run_eff


        #[Pipe_2b, Pipe_1c_fake, HTHE_DelT_new, effectiveness_heatex] = heat_exchanger_counter(Pipe_2a, Pipe_1b, T2out=cfg.reactor_T_1c)
        HTHE_DelT_resid = abs(HTHE_P - HTHE_P_new)
        HTHE_P = HTHE_P_new

        #print(f'HE eff = {effectiveness_heatex:0.4f}')
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 2 (outlet to water) ~~~~~~~~~~~~~~~~~~~~~~~~~~~


        # print(effectiveness_chiller)
        # print('\n')
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Condenser ~~~~~~~~~~~~~~~~~~~~~~~~~~~


        #[Pipe_RE, power_consumption["Condenser"], ammonia_removed, condenser_water_out_temp] = condenser_crude(Pipe_2c, water_mass_flow=cfg.condenser_water_mfr, T_cin=cfg.condenser_T_cold_in)

        [Pipe_2c, power_consumption["Condenser"], condenser_water_out_temp] = tristan_condenser(Pipe_2b, Condenser_Details())

        ammonia_produced = Pipe_2b.NH3 - Pipe_2c.NH3
        ammonia_removed = ammonia_produced/Pipe_2b.NH3

        if (Pipe_2c.mol_tot + Pipe_IN.mol_tot) > cfg.max_mol:
            bool_purged = True
            print(f'Purged stream = {Pipe_2c.mol_tot:3.3f}')
            k = (cfg.max_mol - Pipe_IN.mol_tot) / Pipe_2c.mol_tot
            [Pipe_RE, Pipe_purged] = Pipe_2c.split(k)
            note = "purged"
        else:
            Pipe_RE = copy.copy(Pipe_2c)





        #set up purge stream if non-converging




        if ops.TERMINAL_LOG:
            print(' i: %i, convergence: %1.6f' %(count,(2*Pipe_IN.N2 - ammonia_produced)/(2*Pipe_IN.N2)))
        if abs(2*(Pipe_IN.N2 - Pipe_purged.N2) - (ammonia_produced + Pipe_purged.NH3))/(2*Pipe_IN.N2) < cfg.plant_convergence:
            stop = 1

    if ops.TERMINAL_LOG_SHORT:
        print('run_%i iter: %i' %(id_run, count))

    power_consumption["n2_m_s"] = cfg.plant_n2
    power_consumption["h2_m_s"] = cfg.plant_h2
    power_consumption["ammonia_produced"] = ammonia_produced

    power_consumption['reactor_mol_in'] = Pipe_1d.mol_tot
    power_consumption["reactor_exotherm"] = exotherm_reac
    power_consumption["reactor_heatloss"] = heatloss_reac

    power_consumption["conversion"] = (Pipe_2a.NH3 - Pipe_1d.NH3)/(2*Pipe_1d.N2)
    power_consumption["reactor_in_temp"] = Bed_data[0][3]
    power_consumption["reactor_max_temp"] = max(i[3] for i in Bed_data)
    power_consumption["reactor_out_temp"] = Bed_data[-1][3]


    power_consumption["heat_exchanger_power_exchange"] = HTHE_P
    power_consumption['heat_exchanger_effectiveness'] = effectiveness_heatex
    power_consumption["heat_exchanger_heatloss"] = heat_ex_heatloss

    vel_in_reactor = Pipe_1d.volume_fr/bed1.cs_area
    recycle_ratio_mol = ((Pipe_RE.mol_tot+Pipe_IN.mol_tot) / Pipe_IN.mol_tot)
    recycle_ratio_mass = ((Pipe_RE.mass_tot+Pipe_IN.mass_tot) / Pipe_IN.mass_tot)
    power_consumption["ammonia_%_removed"] = ammonia_removed * 100


    power_consumption["recycle_ratio_mol"] = recycle_ratio_mol
    power_consumption["recycle_ratio_mass"] = recycle_ratio_mass




    power_consumption["Purged"] = Pipe_purged.mol_tot/Pipe_2c.mol_tot
    power_consumption["note"] = note

    if ops.TERMINAL_LOG or ops.TERMINAL_LOG_SHORT:
        if note != '':
            print("note = " + note)
        '''
        print('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')
        # Mix recycle stream in
        print('New Feed (cooled) = %3.1f' % Pipe_IN.T + ' K')
        print('    Prechill eff = %1.3f' % effectiveness_precooler)
        print('Recycle = %3.1f' % Pipe_RE.T + ' K')
        print('    Mixed with Recycle = %3.1f' % Pipe_1a.T + ' K')
        # recompress recycled stream
        print('Recompressed = %3.1f' % Pipe_1b.T + ' K')
        # Add heat from Low Temp Heat Exchanger to Pipe 1b to make Pipe 1c
        print('Post heat exchanger stream = %3.1f' % Pipe_1c.T + ' K')
        print('Reheated stream = %3.1f' % Pipe_1d.T + ' K')
        print('Bed 1 length = %2.2fm, conversion = %2.2f' % (bed1.vect[-1], (Pipe_2a.NH3 - Pipe_1c.NH3) / (
                2 * Pipe_1c.N2) * 100) + '%' + ', T = %3.1f' % Pipe_2a.T + 'K')
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 1 (Pipe 6 to Pipe 4) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('exotherm minus heatloss = %3.3f W ' % exotherm_minus_heatloss)
        print('Immediately post reactor = %3.1f' % Pipe_2a.T + ' K')
        print('    Velocity of Flow in reac = %3.1f'% vel_in_reactor)
        print('Single cooled post reactor = %3.1f' % Pipe_2b.T + ' K')
        print('    HTHE P = %3.1f' % HTHE_P)
        print('    heat exchanger eff = %1.3f' % effectiveness_heatex)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Condensor ~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('ammonia removed = %2.2f' % ((ammonia_removed) * 100) + ' %')
        print('    initial ammonia molar = %1.3f' % Pipe_2b.yNH3)
        print('    final ammonia molar = %1.3f' % Pipe_2c.yNH3)
        print('    condenser water out temp = %3.1f' % condenser_water_out_temp)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Recycle ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # reheat recycle stream? loosing too much energy rn
        print('Ammonia produced = %2.4f g/s' % (ammonia_produced*17.03))


        print('recycle ratio mass = %2.3f' % recycle_ratio_mass)
        print('recycle ratio mol = %2.3f' % recycle_ratio_mol)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 1 (pt 2?) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # print('power consumption = ', power_consumption)
    '''
    if ops.REACTOR_BED:

        Bed_data_T = np.array(Bed_data).T.tolist()
        x_plot_data = bed1.vect
        y_plot_data = Bed_data_T[3]
        plt.plot(x_plot_data, y_plot_data)
        plt.title("Reactor Bed Temeprature Profile")
        plt.xlabel("Position Along Length (m)")
        plt.ylabel("Temperature (K)")

    # print(bedvect)

    stream_list = [Pipe_IN.store(),
                   Pipe_1a.store(),
                   Pipe_1b.store(),
                   Pipe_1c.store(),
                   Pipe_1d.store(),
                   Pipe_2a.store(),
                   Pipe_2b.store(),
                   Pipe_2c.store(),
                   Pipe_RE.store(),
                   Pipe_purged.store()]

    pipe_locs = ["IN", "1a", "1b", "1c", "1d", "2a", "2b", "2c", "RE", 'Purge']

    stream_data = pd.DataFrame(stream_list,
                               columns=['h2_mol_s',
                                        'n2_mol_s',
                                        'nh3_mol_s',
                                        'temperature',
                                        'pressure'
                                        ],
                               index=pipe_locs
                               )

    power_data = pd.DataFrame.from_dict(power_consumption, orient='index',
                                        columns=["run_%d" % id_run])

    if ops.TERMINAL_END_LOG_DETAIL:
        print('\nStream data for run_%d' % id_run)
        print(stream_data)
        print('\nPower data for run_%d' % id_run)
        print(power_data)

    return stream_data, power_data


def multi_run(cfg, ops, param=None, vals=None):
    """ """
    # preparing for multiple runs
    if vals is None or param is None:
        n_runs = 1
        rewrite_config = False
        cfg.plant_n2 = cfg.plant_h2 / cfg.plant_ratio_n  # just incase selected parameter is h2; to maintain ratio
    else:
        n_runs = len(vals)
        rewrite_config = True

    run_headers = []
    power_lst = []
    n2_lst = []
    h2_lst = []
    nh3_lst = []
    temperature_lst = []
    pressure_lst = []

    for i in range(n_runs):
        if rewrite_config:
            exec('cfg.' + param + ' = %0.5f' % vals[i])  # this is bad practise I know
            cfg.plant_n2 = cfg.plant_h2 / cfg.plant_ratio_n  # just incase selected parameter is h2; to maintain ratio
        run_headers.append("run_%d" % i)
        stream_temp, power_temp = evaluate_loop(cfg, ops, i)
        power_lst.append(power_temp.iloc[:, 0].tolist())
        n2_lst.append(stream_temp['h2_mol_s'].tolist())
        h2_lst.append(stream_temp['n2_mol_s'].tolist())
        nh3_lst.append(stream_temp['nh3_mol_s'].tolist())
        temperature_lst.append(stream_temp['temperature'].tolist())
        pressure_lst.append(stream_temp['pressure'].tolist())


    stream_indices = stream_temp.index
    power_indices = power_temp.index.values.tolist()
    power_data = pd.DataFrame(tps_lst(power_lst), index=power_indices, columns=run_headers)
    n2_data  = pd.DataFrame(tps_lst(n2_lst), index=stream_indices, columns=run_headers)
    h2_data = pd.DataFrame(tps_lst(h2_lst), index=stream_indices, columns=run_headers)
    nh3_data = pd.DataFrame(tps_lst(nh3_lst), index=stream_indices, columns=run_headers)
    temperature_data = pd.DataFrame(tps_lst(temperature_lst), index=stream_indices, columns=run_headers)
    pressure_data = pd.DataFrame(tps_lst(pressure_lst), index=stream_indices, columns=run_headers)


    return power_data, n2_data, h2_data, nh3_data, temperature_data, pressure_data


def params():
    # hardcoding param_sweep for now, will eventually be improved
    chosen_param = 'plant_h2'
    rng = np.linspace(0.05, 0.3, 6)  # 0.1,0.4,13 for more detail
    return chosen_param, rng


def main():
    """ Main to run the ammonia plant"""


    [chosen_param, rng] = params()

    cfg, ops = get_configs(args)

    power, n2, h2, nh3, temperature, pressure = multi_run(cfg, ops, param=chosen_param, vals=rng)

    createFolder('./outputs/')
    power.to_csv('outputs/power.csv')
    n2.to_csv('outputs/n2.csv')
    h2.to_csv('outputs/h2.csv')
    nh3.to_csv('outputs/nh3.csv')
    temperature.to_csv('outputs/temperature.csv')
    pressure.to_csv('outputs/pressure.csv')


    # power has a unique set of indices, see the following print out as an example:
    print('Output for power:')
    print(power)  # look at the index column to see what keys are valid (PSA, etc...)

    # n2, h2, nh3, temperature, pressure are separate DataFrames with the same indices.
    # See the following print out as an example:
    print('Output for a state variable:')
    print(temperature)  # look at the index column to see what keys are valid (IN, 1a, 1b, etc...)

    #chosen_solution = "ammonia_produced"
    #plant_figure(chosen_solution, power)


def plant_figure(chosen_solution, power):
    [chosen_param, rng] = params()
    plt.figure
    plt.plot(rng, power.loc[chosen_solution])
    plt.xlabel(chosen_param)
    plt.ylabel(chosen_solution)
    plt.show()

    # decide how to pass in the variables for multiple runs


def single_run():
    cfg, ops = get_configs(args)
    ops.TERMINAL_LOG = True
    streamtemp,powertemp = evaluate_loop(cfg,ops,1)
    print(streamtemp)
    print(powertemp)
    plt.show()


def read_and_plot():
    '''
    'PSA', 'electrolysis', 'N2_comp', 'H2_comp', 'inflow cooling',
       'recompressor', 'heater', 'Condenser', 'reactor_mol',
       'recycle_ratio_mol', 'recycle_ratio_mass', 'ammonia_produced',
       'ammonia_removed', 'n2_m_s', 'h2_m_s', 'exotherm_minus_heatloss', 'Purged'
    '''
    chosen_solution = "condensor"
    power = pandas.read_csv('outputs/power.csv',index_col=0,header=0)
    print(power.index)
    plant_figure(chosen_solution, power)


if __name__ == "__main__":
    single_run()
