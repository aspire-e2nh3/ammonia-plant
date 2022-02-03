from isobaric_reactor import *
import matplotlib.pyplot as plt
import seaborn as sns
import math
import numpy as np
import copy
import pyromat as pm

# initial inputs - define bed structure, quench ratio
bed1 = Bed(2.0, 0.05, 0.001, False)
'''
quench_ratio = 0.6 # amount into first bed/total amount
bed2 = Bed(3, 0.05, 0.001, False)
'''

# criterion for stopping heat exchanger integration
Criterion = 0.001

# set reactor inlet temp
reactor_in_temp = 663

# input mols of n2 and h2
total_mol_N2 = 0.1     # mol
total_mol_H2 = total_mol_N2 * 3 # mol

# initialise iteration limits
inlet_temp = 398
HTHE_DelT = reactor_in_temp - inlet_temp
HTHE_DelT_resid = 1
recycle_estimate = 6.2



#default pressure
Reactor_Pressure = 200  # bar

#initialise power consumption dictionary
power_consumption = {}


# make n2 and h2 input lines from psa and electrolysis functions
[Pipe_N2_LP, power_consumption["PSA"]] = psa_estimate(total_mol_N2) #              MAYBE ADD AIR IN?
[Pipe_H2_LP, power_consumption["electrolysis"],Water_in] = electrolysis(total_mol_H2)

print('H2 LP Feed = %3.1f' % Pipe_H2_LP.T, 'K\n')
print('N2 LP Feed = %3.1f' % Pipe_N2_LP.T, 'K\n')

# compress N2 and H2 lines
[Pipe_N2_HP, power_consumption["N2_comp"]] = compressor(Pipe_N2_LP, Reactor_Pressure)
[Pipe_H2_HP, power_consumption["H2_comp"]] = compressor(Pipe_H2_LP, Reactor_Pressure)

print('H2 HP Feed = %3.1f' % Pipe_H2_HP.T, 'K\n')
print('N2 HP Feed = %3.1f' % Pipe_N2_HP.T, 'K\n')

Pipe_IN = mixer(Pipe_H2_HP,Pipe_N2_HP)

print('Mixed Feed pre-cooling = %3.1f' % Pipe_IN.T, 'K\n')

[Pipe_IN_cooled,power_consumption["inflow cooling"]] = heat_exchanger_water2gas(Pipe_IN, 10, inlet_temp)

Pipe_RE = State(recycle_estimate*total_mol_H2,recycle_estimate*total_mol_N2, 0.5*total_mol_N2, inlet_temp, Reactor_Pressure-2)



print(power_consumption)

stop = 0
while (stop == 0):

    print('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')
    # Mix recycle stream in
    print('New Feed (cooled) = %3.1f' % Pipe_IN_cooled.T, 'K\n')
    print('Recycle = %3.1f' % Pipe_RE.T, 'K\n')
    Pipe_1a = mixer(Pipe_IN_cooled, Pipe_RE)
    print('Mixed with Recycle = %3.1f' % Pipe_1a.T, 'K\n')

    # recompress recycled stream
    [Pipe_1b,power_consumption["recompressor"]] = compressor(Pipe_1a, 200, 0.7)
    print('Recompressed = %3.1f' % Pipe_1b.T, 'K\n')

    # Add heat from Low Temp Heat Exchanger to Pipe 1b to make Pipe 1c
    Pipe_1c = copy.deepcopy(Pipe_1b)
    Pipe_1c.T += HTHE_DelT
    print('Reheated stream = %3.1f' % Pipe_1c.T, 'K\n')
    #Pipe_1c.T = 673
    #print('Reheated stream set to 673K.')


    # Initialise bed data
    Bed_data = []
    Bed_iterator = copy.deepcopy(Pipe_1c)
    Bed_data.append(Bed_iterator.store())

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ BED 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for X in range(bed1.vectlen - 1):
        # setup new step
        dX = bed1.vect[X + 1] - bed1.vect[X]

        # run reactor step
        Bed_iterator = reactorStep(Bed_iterator, dX, bed1.area)

        # store data
        Bed_data.append(Bed_iterator.store())


    print('Bed 1 length = %2.2fm, conversion = %2.2f' % (bed1.vect[-1], Bed_iterator.NH3/(2*Pipe_1c.N2)*100) + '%' + ', T = %3.1f' % Bed_iterator.T + 'K')


    '''
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
    '''
    Pipe_2a = copy.deepcopy(Bed_iterator)
    Pipe_2a.p -= 2

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 1 (Pipe 6 to Pipe 4) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('Immediately post reactor = %3.1f' % Pipe_2a.T, 'K')
    [Pipe_2b, Pipe_1c_fake, HTHE_DelT_new] = heat_exchanger_parallel(Pipe_2a, Pipe_1b,T2out=reactor_in_temp)
    '''
    Pipe_2b = copy.deepcopy(Pipe_2a)
    Pipe_1c_fake = copy.deepcopy(Pipe_1b)
    HTHE_DelT_new = max(673 - Pipe_1c_fake.T,0)
    Pipe_2b.T += - HTHE_DelT_new*Pipe_1c_fake.cp*Pipe_1c_fake.mass_tot/(Pipe_2b.mass_tot * Pipe_2b.cp)
    Pipe_1c_fake.T = 673
    '''

    HTHE_DelT_resid = abs(HTHE_DelT - HTHE_DelT_new)
    print('Single cooled post reactor = %3.1f' % Pipe_2b.T,
          'K\n Post HE Inlet into Reactor = %3.1f' % Pipe_1c_fake.T,
          'K\n HTHE del T = %3.1f' % HTHE_DelT_new, ' Resid = %3.1f' % HTHE_DelT_resid)
    HTHE_DelT = HTHE_DelT_new
    '''
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 2 (Pipe 7 to inlet) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #[Pipe_4c, Pipe_1c_fake, LTHE_DelT_new] = heat_exchanger_hotgas2coldgas(Pipe_4b, Pipe_1b)
    Pipe_4c = copy.deepcopy(Pipe_4a)
    Pipe_1c_fake = copy.deepcopy(Pipe_1b)
    LTHE_DelT_new = max(573 - Pipe_1c_fake.T,0)
    Pipe_4c.T += - HTHE_DelT_new*Pipe_1c_fake.cp*Pipe_1c_fake.mass_tot/(Pipe_4c.mass_tot * Pipe_4c.cp)
    Pipe_1c_fake.T = 573
    LTHE_DelT_resid = abs(LTHE_DelT - LTHE_DelT_new)
    print('Double cooled post reactor = %3.1f' % Pipe_4c.T,
          'K\nPost HE quench stream = %3.1f' % Pipe_1c.T,
          'K\n LTHE del T = %3.1f' % LTHE_DelT_new, ' Resid = %3.1f' % LTHE_DelT_resid)
    LTHE_DelT = LTHE_DelT_new
    '''

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 2 (outlet to water) ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    [Pipe_2c,power_consumption["Chiller"]] = heat_exchanger_water2gas(Pipe_2b, 10)
    print('Double cooled chiller outlet T = %3.1f' % Pipe_2c.T)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Condensor ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    [Pipe_RE, power_consumption["Condensor"]] = condensor(Pipe_2c, e2=0.8, water_mass_flow=10)
    print('Recycle stream composition = ', Pipe_RE.store())

    if HTHE_DelT_resid < Criterion:
        stop = 1

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Recycle ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reheat recycle stream? loosing too much energy rn
print('Ammonia produced = %2.4f g/s' % float((Pipe_2c.mNH3-Pipe_RE.mNH3)*1e3))


recycle_ratio_mol = (Pipe_RE.mol_tot / Pipe_IN.mol_tot)
recycle_ratio_mass = (Pipe_RE.mass_tot / Pipe_IN.mass_tot)
print('recycle ratio mass = %2.3f' % recycle_ratio_mass, ', recycle ratio mol = %2.3f' % recycle_ratio_mol)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 1 (pt 2?) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('power consumption = ', power_consumption)

Bed_data_T = np.array(Bed_data).T.tolist()

print(Pipe_IN_cooled.store())
print(Pipe_1a.store())
print(Pipe_1b.store())
print(Pipe_1c.store())
print(Pipe_2a.store())
print(Pipe_2b.store())
print(Pipe_2c.store())
print(Pipe_RE.store())
plot = 1
if plot:
    x_plot_data = bed1.vect

    y_plot_data = Bed_data_T[3]

    sns.set(style='whitegrid')
    sns.scatterplot(x=x_plot_data, y=y_plot_data)
    plt.show()

# print(bedvect)
