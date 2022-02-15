from isobaric_reactor import *
import matplotlib.pyplot as plt
import math
import numpy as np
import copy
import pyromat as pm

# initial inputs - define bed structure, quench ratio
bed1 = Bed(2.0, 0.05, 0.001, True)
bed1.mass()
'''
quench_ratio = 0.6 # amount into first bed/total amount
bed2 = Bed(3, 0.05, 0.001, False)
'''

# criterion for stopping heat exchanger integration
Criterion = 0.01

# set reactor inlet temp
reactor_in_temp = 663

# input mols of n2 and h2
total_mol_N2 = 0.1   # mol
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
[Pipe_N2_LP, power_consumption["PSA"]] = psa_estimate(total_mol_N2)                    # MAYBE ADD AIR_IN?
[Pipe_H2_LP, power_consumption["electrolysis"],Water_in] = electrolysis(total_mol_H2)

print('H2 LP Feed = %3.1f' % Pipe_H2_LP.T, 'K\n')
print('N2 LP Feed = %3.1f' % Pipe_N2_LP.T, 'K\n')

# compress N2 and H2 lines
[Pipe_N2_HP, power_consumption["N2_comp"],heatlost1] = ptcompressor(Pipe_N2_LP, Reactor_Pressure, t_out=inlet_temp+200)
[Pipe_H2_HP, power_consumption["H2_comp"],heatlost2] = ptcompressor(Pipe_H2_LP, Reactor_Pressure, t_out=inlet_temp+100)

print('H2 HP Feed = %3.1f' % Pipe_H2_HP.T, 'K\n')
print('N2 HP Feed = %3.1f' % Pipe_N2_HP.T, 'K\n')

# mix H2 and N2 lines
Pipe_IN = mixer(Pipe_H2_HP,Pipe_N2_HP)

print('Mixed Feed pre-cooling = %3.1f' % Pipe_IN.T, 'K\n')

[Pipe_IN_cooled,power_consumption["inflow cooling"],_] = heat_exchanger_water2gas(Pipe_IN, 10, inlet_temp)

Pipe_RE = State(recycle_estimate*total_mol_H2,recycle_estimate*total_mol_N2, 0.5*total_mol_N2, inlet_temp-80, Reactor_Pressure-2)

print(power_consumption)

stop = 0
while (stop == 0):

    # Mix recycle stream in
    Pipe_1a = mixer(Pipe_IN_cooled, Pipe_RE)


    # recompress recycled stream
    [Pipe_1b,power_consumption["recompressor"],_] = compressor(Pipe_1a, p_out=200, eta=0.7)


    # Add heat from Low Temp Heat Exchanger to Pipe 1b to make Pipe 1c
    Pipe_1c = copy.copy(Pipe_1b)
    Pipe_1c.T += HTHE_DelT

    [Pipe_1d,power_consumption["heater"]] = heater(Pipe_1c, reactor_in_temp)

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

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 1 (Pipe 1b to Pipe 2a) ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # estimate heat exchange variant
    [Pipe_2b, Pipe_1c_fake, HTHE_DelT_new, effectiveness_heex] = heat_exchanger_counter(Pipe_2a, Pipe_1b)
    #print(effectiveness_heex)
    HTHE_DelT_resid = abs(HTHE_DelT - HTHE_DelT_new)

    HTHE_DelT = HTHE_DelT_new

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 2 (2b to water) ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    [Pipe_2c,power_consumption["Chiller"],effectiveness_cool] = heat_exchanger_water2gas(Pipe_2b, e=0.8)
    #print(effectiveness_cool)
    #print('\n')
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Condenser ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    [Pipe_RE, power_consumption["Condenser"]] = condensor(Pipe_2c, e2=0.8, water_mass_flow=10)

    print(Pipe_RE.yH2)
    print(HTHE_DelT_resid)
    if HTHE_DelT_resid < Criterion:
        stop = 1


print('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')
# Mix recycle stream in
print('New Feed (cooled) = %3.1f' % Pipe_IN_cooled.T, 'K')
print('Recycle = %3.1f' % Pipe_RE.T, 'K')

print('Mixed with Recycle = %3.1f' % Pipe_1a.T, 'K')

# recompress recycled stream

print('Recompressed = %3.1f' % Pipe_1b.T, 'K')

# Add heat from Low Temp Heat Exchanger to Pipe 1b to make Pipe 1c

print('Reheated stream = %3.1f' % Pipe_1c.T, 'K')

print('heated to 663K stream = %3.1f' % Pipe_1d.T, 'K')

print('Bed 1 length = %2.2fm, conversion = %2.2f' % (bed1.vect[-1], Bed_iterator.NH3/(2*Pipe_1c.N2)*100) + '%' + ', T = %3.1f' % Bed_iterator.T + 'K')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 1 (Pipe 6 to Pipe 4) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('Immediately post reactor = %3.1f' % Pipe_2a.T, 'K')

print('Single cooled post reactor = %3.1f' % Pipe_2b.T,
      'K\n Post HE Inlet into Reactor = %3.1f' % Pipe_1c_fake.T,
      'K\n HTHE del T = %3.1f' % HTHE_DelT_new, ' Resid = %3.1f' % HTHE_DelT_resid)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 2 (outlet to water) ~~~~~~~~~~~~~~~~~~~~~~~~~~~


print('Double cooled chiller outlet T = %3.1f' % Pipe_2c.T)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Condensor ~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Recycle ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reheat recycle stream? loosing too much energy rn
print('Ammonia produced = %2.4f g/s' % float((Pipe_2c.mNH3-Pipe_RE.mNH3)*1e3))


recycle_ratio_mol = (Pipe_RE.mol_tot / Pipe_IN.mol_tot)
recycle_ratio_mass = (Pipe_RE.mass_tot / Pipe_IN.mass_tot)
print('recycle ratio mass = %2.3f' % recycle_ratio_mass, ', recycle ratio mol = %2.3f' % recycle_ratio_mol)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 1 (pt 2?) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('power consumption = ', power_consumption)

print('\n Stream data :')
Bed_data_T = np.array(Bed_data).T.tolist()

print(bed1.vect)

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


    plt.plot(x_plot_data,y_plot_data)
    plt.show()


