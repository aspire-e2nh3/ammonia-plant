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
reactor_in_temp = 673
cool_water_temp = 273
precool = 1

# input mols of n2 and h2
total_mol_N2 = 0.1   # mol
total_mol_H2 = total_mol_N2 * 3 # mol

# initialise iteration limits
inlet_temp = 373
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
[Pipe_N2_HP, power_consumption["N2_comp"],heatlost1,n1] = ptcompressor(Pipe_N2_LP, Reactor_Pressure, t_out=inlet_temp+350)
[Pipe_H2_HP, power_consumption["H2_comp"],heatlost2,n2] = ptcompressor(Pipe_H2_LP, Reactor_Pressure, t_out=inlet_temp+250)

print('H2 HP Feed = %3.1f' % Pipe_H2_HP.T, 'K\n')
print('N2 HP Feed = %3.1f' % Pipe_N2_HP.T, 'K\n')

print('N2 Compressor Polytropic Index = %3.3f' % n1)
print('H2 Compressor Polytropic Index = %3.3f' % n2)

# mix H2 and N2 lines
Pipe_IN = mixer(Pipe_H2_HP,Pipe_N2_HP)

print('Mixed Feed pre-cooling = %3.1f' % Pipe_IN.T, 'K\n')

if precool == 1:
    [Pipe_IN_cooled,power_consumption["inflow cooling"],effectiveness_prechill] = heat_exchanger_water2gas(Pipe_IN,
                                                                                                           T_end=inlet_temp,
                                                                                                           effectiveness=0.8,
                                                                                                           cool_to_sat_point=False)
else:
    Pipe_IN_cooled = Pipe_IN
    effectiveness_prechill = 0

Pipe_RE = State(recycle_estimate*total_mol_H2,recycle_estimate*total_mol_N2, 0.5*total_mol_N2, inlet_temp-80, Reactor_Pressure-2)

print(power_consumption)



stop = 0
count = 0
while (stop == 0):
    count +=1
    if count >1000:
        stop = 1
        print('Failed to converge')
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
    [Pipe_2b, Pipe_1c_fake, HTHE_DelT_new, effectiveness_heatex] = heat_exchanger_counter(Pipe_2a, Pipe_1b, T2out=reactor_in_temp)

    HTHE_DelT_resid = abs(HTHE_DelT - HTHE_DelT_new)

    HTHE_DelT = HTHE_DelT_new

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 2 (2b to water) ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    [Pipe_2c,power_consumption["Chiller"],effectiveness_chiller] = heat_exchanger_water2gas(Pipe_2b,
                                                                                            T_cin=cool_water_temp,
                                                                                            cool_to_sat_point=True)
    #print('\n')
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Condenser ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    [Pipe_RE, power_consumption["Condenser"], ammonia_removed, condenser_water_out_temp] = condenser_crude(Pipe_2c, water_mass_flow=10, T_cin=cool_water_temp)

    print('%i, %1.4f' % (count, HTHE_DelT_resid))
    print()
    if HTHE_DelT_resid < Criterion:
        stop = 1


print('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')
# Mix recycle stream in
print('New Feed (cooled) = %3.1f' % Pipe_IN_cooled.T + ' K')
print('    Prechill eff = %1.3f' % effectiveness_prechill)
print('Recycle = %3.1f' % Pipe_RE.T + ' K')
print('    Mixed with Recycle = %3.1f' % Pipe_1a.T + ' K')
# recompress recycled stream
print('Recompressed = %3.1f' % Pipe_1b.T + ' K')
# Add heat from Low Temp Heat Exchanger to Pipe 1b to make Pipe 1c
print('Post heat exchanger stream = %3.1f' % Pipe_1c.T + ' K')
print('Reheated stream = %3.1f' % Pipe_1d.T + ' K')
print('Bed 1 length = %2.2fm, conversion = %2.2f' % (bed1.vect[-1], (Bed_iterator.NH3 - Pipe_1c.NH3)/(2*Pipe_1c.N2)*100) + '%' + ', T = %3.1f' % Bed_iterator.T + 'K')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 1 (Pipe 6 to Pipe 4) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('Immediately post reactor = %3.1f' % Pipe_2a.T + ' K')
print('Single cooled post reactor = %3.1f' % Pipe_2b.T + ' K')
print('    HTHE del T = %3.1f' % HTHE_DelT_new)
print('    heat exchanger eff = %1.3f' % effectiveness_heatex)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 2 (outlet to water) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('Double cooled chiller outlet T = %3.1f' % Pipe_2c.T)
print('    chiller eff = %1.3f' % effectiveness_chiller)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Condensor ~~~~~~~~~~~~~~~~~~~~~~~~~~
print('ammonia removed = %2.2f' % ((ammonia_removed)*100) + ' %')
print('    initial ammonia molar = %1.3f' % Pipe_2c.yNH3)
print('    recycle ammonia molar = %1.3f' % Pipe_RE.yNH3)
print('    condenser water out temp = %3.1f' % condenser_water_out_temp)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Recycle ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reheat recycle stream? loosing too much energy rn
print('Ammonia produced = %2.4f g/s' % float((Pipe_2c.mNH3-Pipe_RE.mNH3)*1e3))


recycle_ratio_mol = (Pipe_RE.mol_tot / Pipe_IN.mol_tot)
recycle_ratio_mass = (Pipe_RE.mass_tot / Pipe_IN.mass_tot)
print('recycle ratio mass = %2.3f' % recycle_ratio_mass, ', recycle ratio mol = %2.3f' % recycle_ratio_mol)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 1 (pt 2?) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('power consumption = ', power_consumption)

print_stream_data = 0
if print_stream_data == 1:
    print('\n Stream data :')


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
    Bed_data_T = np.array(Bed_data).T.tolist()

    x_plot_data = bed1.vect

    y_plot_data = Bed_data_T[3]


    plt.plot(x_plot_data,y_plot_data)

    plt.title("Reactor Bed Temeprature Profile")
    plt.xlabel("Position Along Length (m)")
    plt.ylabel("Temperature (K)")

    plt.show()


