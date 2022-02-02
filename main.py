from ammonia_plant.isobaric_reactor import *
import matplotlib.pyplot as plt
import seaborn as sns
import math
import numpy as np
import copy
import pyromat as pm

# initial inputs

bed1 = Bed(1, 0.05, 0.001, False)
quench_ratio = 0.75 # amount into first bed/total amount
bed2 = Bed(3, 0.05, 0.001, False)

bed3bool = 0
q2fact = 1
bed3 = Bed(7, 0.1, 0.001, False)

Criterion = 0.1



total_mol_N2 = 0.15  # mol
total_mol_H2 = total_mol_N2 * 3 # mol


Reactor_Pressure = 200  # bar
Inlet_Temperature = 673  # K



power_consumption = {}

[Pipe_N2_LP, power_consumption["PSA"]] = psa_estimate(total_mol_N2) #              MAYBE ADD AIR IN?
[Pipe_H2_LP, power_consumption["electrolysis"],Water_in] = electrolysis(total_mol_H2)

print('H2 LP Feed = %3.1f' % Pipe_H2_LP.T, 'K\n')
print('N2 LP Feed = %3.1f' % Pipe_N2_LP.T, 'K\n')

[Pipe_N2_HP, power_consumption["N2_comp"]] = compressor(Pipe_N2_LP, Reactor_Pressure)
[Pipe_H2_HP, power_consumption["H2_comp"]] = compressor(Pipe_H2_LP, Reactor_Pressure)

print('H2 HP Feed = %3.1f' % Pipe_H2_HP.T, 'K\n')
print('N2 HP Feed = %3.1f' % Pipe_N2_HP.T, 'K\n')

Pipe_IN = mixer(Pipe_H2_HP,Pipe_N2_HP)

print('Mixed Feed pre-cooling = %3.1f' % Pipe_IN.T, 'K\n')

[Pipe_IN_cooled,power_consumption["inflow cooling"]] = heat_exchanger_water2gas(Pipe_IN, 1, 298)

Pipe_recycle = State(4*total_mol_H2,4*total_mol_N2, 0.05*total_mol_N2, 298, Reactor_Pressure-1)

HTHE_DelT = 200  #
LTHE_DelT = 150
HTHE_DelT_resid = 50


print(power_consumption)

stop = 0
while (stop == 0):

    print('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')
    # Mix recycle stream in
    print('New Feed (cooled) = %3.1f' % Pipe_IN_cooled.T, 'K\n')
    print('Recycle = %3.1f' % Pipe_recycle.T, 'K\n')
    Pipe_1a = mixer(Pipe_IN_cooled, Pipe_recycle)
    print('Mixed with Recycle = %3.1f' % Pipe_1a.T, 'K\n')

    # recompress recycled stream
    [Pipe_1b,power_consumption["recompressor"]] = compressor(Pipe_1a, 200, 0.7)
    print('Recompressed = %3.1f' % Pipe_1b.T, 'K\n')

    # Add heat from Low Temp Heat Exchanger to Pipe 1b to make Pipe 1c
    Pipe_1c = copy.deepcopy(Pipe_1b)
    Pipe_1c.T += LTHE_DelT
    print('Reheated stream = %3.1f' % Pipe_1c.T, 'K\n')

    # split pipe for quench stream:
    [Pipe_2a, Pipe_3] = Pipe_1c.split(quench_ratio)

    # Add heat from High Temp Heat Exchanger to Pipe 2a to make Pipe 2b
    Pipe_2b = copy.deepcopy(Pipe_2a)
    Pipe_2b.T += HTHE_DelT
    print('Double Reheated stream = %3.1f' % Pipe_2b.T, 'K\n')
    # Initialise bed data
    Bed_data = []
    Bed_iterator = copy.deepcopy(Pipe_2b)
    Bed_data.append(Bed_iterator.store())

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ BED 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for X in range(bed1.vectlen - 1):
        # setup new step
        dX = bed1.vect[X + 1] - bed1.vect[X]

        # run reactor step
        Bed_iterator = reactorStep(Bed_iterator, dX, bed1.area)

        # store data
        Bed_data.append(Bed_iterator.store())


    print('Bed 1 length = %2.1fm, conversion = %2.2f' % (bed1.vect[-1], Bed_iterator.NH3/(2*Pipe_2a.N2)*100) + '%' + ', T = %3.1f' % Bed_iterator.T + 'K')
    print('Ammonia produced = %2.2f kg/h' % (Bed_iterator.NH3 * 3600 * 17.03 / 1000))

    # ~~~~~~~~~~~~~~~~~~~~~~~~ QUENCH 1-2 ~~~~~~~~~~~~~~~~~~~~~~~~~
    # quench pipe by mixing with unused unheated stream Pipe_5
    Bed_iterator = mixer(Bed_iterator, Pipe_3)
    Bed_data.append(Bed_iterator.store())
    print('Midbed post-quench T = %3.1f' % Bed_iterator.T + 'K')
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ BED 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for X in range(bed2.vectlen - 1):
        # setup new step
        dX = bed1.vect[X + 1] - bed1.vect[X]

        # run reactor step
        Bed_iterator = reactorStep(Bed_iterator, dX, bed1.area)

        # store data
        Bed_data.append(Bed_iterator.store())

    print('Bed 2 length = %2.1fm, conversion = %2.2f' % (bed2.vect[-1], Bed_iterator.NH3 / (2 * Pipe_1c.N2) * 100) + '%' + ', T = %3.1f' % Bed_iterator.T + 'K')
    print('Ammonia produced = %2.2f kg/h' % (Bed_iterator.NH3 * 3600 * 17.03 / 1000))

    Pipe_4a = copy.deepcopy(Bed_iterator)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 1 (Pipe 6 to Pipe 4) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('Immediately post reactor = %3.1f' % Pipe_4a.T, 'K\n')
    [Pipe_4b, Pipe_2b_fake, HTHE_DelT_new] = heat_exchanger_hotgas2coldgas(Pipe_4a, Pipe_2a)
    HTHE_DelT_resid = abs(HTHE_DelT - HTHE_DelT_new)
    print('Single cooled post reactor = %3.1f' % Pipe_4b.T,
          'K\n Post HE Inlet into Reactor = %3.1f' % Pipe_2b_fake.T,
          'K\n HTHE del T = %3.1f' % HTHE_DelT_new, ' Resid = %3.1f' % HTHE_DelT_resid)
    HTHE_DelT = HTHE_DelT_new

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 2 (Pipe 7 to inlet) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    [Pipe_4c, Pipe_1c_fake, LTHE_DelT_new] = heat_exchanger_hotgas2coldgas(Pipe_4b, Pipe_1b)
    LTHE_DelT_resid = abs(LTHE_DelT - LTHE_DelT_new)
    print('Double cooled post reactor = %3.1f' % Pipe_4c.T,
          'K\nPost HE quench stream = %3.1f' % Pipe_1c.T,
          'K\n LTHE del T = %3.1f' % LTHE_DelT_new, ' Resid = %3.1f' % LTHE_DelT_resid)
    LTHE_DelT = LTHE_DelT_new


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 2 (outlet to water) ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    [Pipe_4d,power_consumption["Chiller"]] = heat_exchanger_water2gas(Pipe_4c, 1)
    print('Triple cooled chiller outlet T = %3.1f' % Pipe_4d.T)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Condensor ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    [Pipe_recycle, power_consumption["Condensor"]] = condensor(Pipe_4d, 0.8)
    print('recycle stream composition = ', Pipe_recycle.store())


    if HTHE_DelT_resid < Criterion:
        if LTHE_DelT_resid < Criterion:
            stop = 1

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Recycle ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reheat recycle stream? loosing too much energy rn


recycle_ratio_mol = (Pipe_recycle.mol_tot / Pipe_IN.mol_tot)
recycle_ratio_mass = (Pipe_recycle.mass_tot / Pipe_IN.mass_tot)
print('recycle ratio mass = %2.3f' % recycle_ratio_mass, ', recycle ratio mol = %2.3f' % recycle_ratio_mol)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 1 (pt 2?) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('power consumption', power_consumption)

Bed_data_T = np.array(Bed_data).T.tolist()

plot = 1
if plot:
    x_plot_data = bed1.vect + [bed1.length + x for x in bed2.vect]
    if bed3bool:
        x_plot_data.append([bed1.length + bed2.length + x for x in bed3.vect])
    # y_plot_data = [x/(2*total_mol_N2) * 100 for x in Bed_data_T[2]]
    y_plot_data = Bed_data_T[3]

    sns.set(style='whitegrid')
    sns.scatterplot(x=x_plot_data, y=y_plot_data)
    plt.show()

# print(bedvect)
