import os
import pandas as pd
from isobaric_reactor import *
import matplotlib.pyplot as plt
import math
import numpy as np
import copy
import pyromat as pm
import argparse
from utils.options import SSConfig, OutOps

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

def main():
    """ Main to run the ammonia plant"""
    
    if args.configuration is None:
        fcfg = os.path.join(os.path.dirname(__file__),'utils','default_ssconfig.ini')
    else:
        fcfg = args.configuration # pathname from input argument
    cfg = SSConfig(fcfg) # initialise

    if args.output is None:
        fops = os.path.join(os.path.dirname(__file__),'utils','default_output_ops.ini')
    else:
        fops = args.output # pathname from input argument
    ops = OutOps(fops) # initialise 

    # initial inputs - define bed structure, quench ratio
    bed1 = Bed(cfg.reactor_length,
               cfg.reactor_diameter,
               cfg.reactor_minimum_step,
               False)
    shell_mass, cat_mass = bed1.mass()

    # set reactor inlet temp
    reactor_in_temp = cfg.reactor_T_1c

    # initialise iteration limits
    inlet_temp = cfg.reactor_T_IN
    HTHE_DelT = reactor_in_temp - inlet_temp
    HTHE_DelT_resid = 1
    recycle_estimate = 6

    #initialise power consumption dictionary
    power_consumption = {}


    # make n2 and h2 input lines from psa and electrolysis functions
    [Pipe_N2_LP, power_consumption["PSA"]] = psa_estimate(cfg.plant_n2)                    # MAYBE ADD AIR_IN?
    [Pipe_H2_LP, power_consumption["electrolysis"],Water_in] = electrolysis(cfg.plant_h2)

    # compress N2 and H2 lines

    [Pipe_N2_HP, power_consumption["N2_comp"],heatlost1, n1] = ptcompressor(Pipe_N2_LP,
                                                                            cfg.plant_pressure,
                                                                            t_out=cfg.reactor_T_IN+cfg.n2compressor_dT)
    [Pipe_H2_HP, power_consumption["H2_comp"],heatlost2, n2] = ptcompressor(Pipe_H2_LP,
                                                                            cfg.plant_pressure,
                                                                            t_out=cfg.reactor_T_IN+cfg.h2compressor_dT)

    # mix H2 and N2 lines
    Pipe_IN_uncooled = mixer(Pipe_H2_HP,Pipe_N2_HP)

    if ops.TERMINAL_LOG:
        print('reactor weight = %3.3f' % shell_mass)
        print('catalyst weight = %3.3f' % cat_mass)
        print('N2 Compressor Polytropic Index = %3.3f' % n1)
        print('H2 Compressor Polytropic Index = %3.3f' % n2)
        print('H2 LP Feed = %3.1f' % Pipe_H2_LP.T, 'K\n')
        print('N2 LP Feed = %3.1f' % Pipe_N2_LP.T, 'K\n')
        print('H2 HP Feed = %3.1f' % Pipe_H2_HP.T, 'K\n')
        print('N2 HP Feed = %3.1f' % Pipe_N2_HP.T, 'K\n')
        print('Mixed Feed pre-cooling = %3.1f' % Pipe_IN_uncooled.T, 'K\n')

    [Pipe_IN,power_consumption["inflow cooling"],effectiveness_precooler] = heat_exchanger_water2gas(Pipe_IN_uncooled,
                                                                                                    T_end=cfg.reactor_T_IN,
                                                                                                    T_cin=cfg.precooler_water_mfr)

    Pipe_RE = State(recycle_estimate * cfg.plant_h2,
                    recycle_estimate * cfg.plant_n2,
                    0.5 * cfg.plant_n2,
                    cfg.reactor_T_IN,
                    cfg.plant_pressure - 2)
    count = 0
    stop = 0
    while (stop == 0):
        count += 1
        # Mix recycle stream in
        Pipe_1a = mixer(Pipe_IN, Pipe_RE)


        # recompress recycled stream
        [Pipe_1b,power_consumption["recompressor"],_] = compressor(Pipe_1a,
                                                                   p_out=cfg.plant_pressure,
                                                                   eta=cfg.recompressor_eta)


        # Add heat from Low Temp Heat Exchanger to Pipe 1b to make Pipe 1c
        Pipe_1c = copy.copy(Pipe_1b)
        Pipe_1c.T += HTHE_DelT

        [Pipe_1d,power_consumption["heater"]] = heater(Pipe_1c, cfg.reactor_T_1c)


        # Initialise bed data
        Bed_data = []
        Bed_iterator = copy.copy(Pipe_1d)
        Bed_data.append(Bed_iterator.store())

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ BED 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        for X in range(bed1.vectlen - 1):
            # setup new step
            dX = bed1.vect[X + 1] - bed1.vect[X]

            # run reactor step
            Bed_iterator = reactorStep(Bed_iterator, dX, bed1.area)

            # store data
            Bed_data.append(Bed_iterator.store())

        Pipe_2a = copy.copy(Bed_iterator)
        Pipe_2a.p -= 2

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 1 (Pipe 6 to Pipe 4) ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # estimate heat exchange variant
        [Pipe_2b, Pipe_1c_fake, HTHE_DelT_new, effectiveness_heatex] = heat_exchanger_counter(Pipe_2a,
                                                                                            Pipe_1b,
                                                                                            T2out=cfg.reactor_T_1c)
        #print(effectiveness_heex)
        HTHE_DelT_resid = abs(HTHE_DelT - HTHE_DelT_new)

        HTHE_DelT = HTHE_DelT_new

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 2 (outlet to water) ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        [Pipe_2c, power_consumption["Chiller"], effectiveness_chiller] = heat_exchanger_water2gas(Pipe_2b,
                                                                                                  T_cin=cfg.chiller_water_mfr,
                                                                                                  cool_to_sat_point=True)
        #print(effectiveness_cool)
        #print('\n')
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Condenser ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        [Pipe_RE, power_consumption["Condenser"], ammonia_removed, condenser_water_out_temp] = condenser_crude(Pipe_2c,
                                                                    water_mass_flow=cfg.condenser_water_mfr,
                                                                    T_cin=cfg.condenser_T_cold_in)

        if ops.TERMINAL_LOG:
            print('%i, %1.4f' % (count, HTHE_DelT_resid))
        if HTHE_DelT_resid < cfg.plant_convergence:
            stop = 1

    Bed_data_T = np.array(Bed_data).T.tolist()

    if ops.TERMINAL_LOG:
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
        print('Bed 1 length = %2.2fm, conversion = %2.2f' % (bed1.vect[-1], (Bed_iterator.NH3 - Pipe_1c.NH3) / (
                    2 * Pipe_1c.N2) * 100) + '%' + ', T = %3.1f' % Bed_iterator.T + 'K')
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 1 (Pipe 6 to Pipe 4) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('Immediately post reactor = %3.1f' % Pipe_2a.T + ' K')
        print('Single cooled post reactor = %3.1f' % Pipe_2b.T + ' K')
        print('    HTHE del T = %3.1f' % HTHE_DelT_new)
        print('    heat exchanger eff = %1.3f' % effectiveness_heatex)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 2 (outlet to water) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('Double cooled chiller outlet T = %3.1f' % Pipe_2c.T)
        print('    chiller eff = %1.3f' % effectiveness_chiller)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Condensor ~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('ammonia removed = %2.2f' % ((ammonia_removed) * 100) + ' %')
        print('    initial ammonia molar = %1.3f' % Pipe_2c.yNH3)
        print('    recycle ammonia molar = %1.3f' % Pipe_RE.yNH3)
        print('    condenser water out temp = %3.1f' % condenser_water_out_temp)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Recycle ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # reheat recycle stream? loosing too much energy rn
        print('Ammonia produced = %2.4f g/s' % float((Pipe_2c.mNH3 - Pipe_RE.mNH3) * 1e3))

        recycle_ratio_mol = (Pipe_RE.mol_tot / Pipe_IN.mol_tot)
        recycle_ratio_mass = (Pipe_RE.mass_tot / Pipe_IN.mass_tot)
        print('recycle ratio mass = %2.3f' % recycle_ratio_mass, ', recycle ratio mol = %2.3f' % recycle_ratio_mol)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 1 (pt 2?) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('power consumption = ', power_consumption)
    
    if ops.REACTOR_BED:
            x_plot_data = bed1.vect
            y_plot_data = Bed_data_T[3]
            plt.plot(x_plot_data,y_plot_data)
            plt.title("Reactor Bed Temeprature Profile")
            plt.xlabel("Position Along Length (m)")
            plt.ylabel("Temperature (K)")
            plt.show()
    # print(bedvect)

    stream_list = [Pipe_IN.store(),
                   Pipe_1a.store(),
                   Pipe_1b.store(),
                   Pipe_1c.store(),
                   Pipe_2a.store(),
                   Pipe_2b.store(),
                   Pipe_2c.store(),
                   Pipe_RE.store()]

    stream_data = pd.DataFrame(stream_list,
                               columns = ['n2_mol_s',
                                          'h2_mol_s',
                                          'nh3_mol_s',
                                          'temperature',
                                          'pressure'
                                          ]
                               )

    pipe_locs = ["IN", "1a", "1b", "1c", "2a", "2b", "2c", "RE"]
    # stream_locs = ["Pipe_" + s for s in pipe_lables]
    stream_data.insert(loc=0, column="pipe_locs", value=pipe_locs)

    # preparing for multiple runs
    n_runs = 1
    run_headers = []
    for i in range(n_runs):
        run_headers.append(["run_%0.0f" % i])
    print(run_headers)

    power_data = pd.DataFrame.from_dict(power_consumption, orient='index',
                                        columns = run_headers)

    if ops.TERMINAL_END_LOG_DETAIL:
        print('\n Stream data :')
        print(stream_data)
        print('\n Power data :')
        print(power_data)


if __name__ == "__main__":
    main()
