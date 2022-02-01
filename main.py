from ammonia_plant.isobaric_reactor import *
import matplotlib.pyplot as plt
import seaborn as sns
import math
import numpy as np
import pyromat as pm

# initial inputs

bed1 = Bed(0.6, 0.05, 0.001, False)
q1fact = 0.3
bed2 = Bed(3, 0.05, 0.001, False)

bed3bool = 0
q2fact = 1
bed3 = Bed(7, 0.1, 0.001, False)

Criterion = 0.1



total_mol_N2 = 0.15  # mol
total_mol_H2 = total_mol_N2 * 3

initial_mol_N2 = total_mol_N2/(1+q1fact+bed3bool*q2fact)
initial_mol_H2 = initial_mol_N2 *3
Reactor_Pressure = 200  # bar
Inlet_Temperature = 298  # K
Inlet_NH3_mol = 0  # mol

air_input = [total_mol_N2,298,1]
[N2_LP_Pipe,PSA_power] = PSA_estimate(air_input)
[H2_LP_Pipe,electrolysis_power,Water_in] = electrolysis(total_mol_H2)

power_consumption = {}
power_consumption["PSA"] = PSA_power
power_consumption["electrolysis"] = electrolysis_power

N2_HP_Pipe = compressor(N2_LP_Pipe, Reactor_Pressure)
H2_HP_Pipe = compressor(H2_LP_Pipe, Reactor_Pressure)

Pipe_1 = State(3 * initial_mol_N2, initial_mol_N2, Inlet_NH3_mol, Inlet_Temperature, Reactor_Pressure)



Pipe_9 = State(6 * initial_mol_N2,2*initial_mol_N2, 0.07*initial_mol_N2, 313, Reactor_Pressure)

HE_1_Del_T = 400  # ?????????????
HE_1_Del_T_resid = 50




stop = 0
while (stop == 0):

    print('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')
    Bed1input = mixer(Binlet,B_recycle)
    Bed1input[3] += HE_1_Del_T
    [Bed1input,P_recomp] = compressor(Bed1input,200,0.7)
    power_consumption_loop.append(P_recomp)
    Bj = Bed1input

    Bed_data = []
    Bed_data.append(Bj)

    print('Post HE Inlet into Reactor = %3.1f' % Bed1input[3], 'K\n')

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ BED 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for X in range(bed1.vectlen - 1):
        # setup new step
        dX = bed1.vect[X + 1] - bed1.vect[X]
        Bi = Bj + [dX, bed1.area]

        # run reactor step
        Bj = reactorStep(Bi)

        # print if necs
        # if bedvect[X] > lastPrintX:
        #    print('Reactor length = %2.1f, conversion = %2.2f' % (bedvect[X], Bi[2] / (2*initial_mol_N2) * 100) + '%' + ', T = %3.1f' %Bi[3])
        #    lastPrintX = np.ceil(bedvect[X]*10)/10

        # store data
        Bed_data.append(Bj)

    Bi = Bj

    print('Bed 1 length = %2.1fm, conversion = %2.2f' % (
    bed1.vect[-1], Bi[2] / (2 * Bed1input[1]) * 100) + '%' + ', T = %3.1f' % Bi[3] + 'K')
    print('Ammonia produced = %2.2f kg/h' % (Bi[2] * 3600 * 17.03 / 1000))

    # ~~~~~~~~~~~~~~~~~~~~~~~~ QUENCH 1-2 ~~~~~~~~~~~~~~~~~~~~~~~~~

    Bq1 = [q1fact * 3 * initial_mol_N2, q1fact * initial_mol_N2, q1fact * Inlet_NH3_mol, Inlet_Temperature,
           Reactor_Pressure]
    Bed2input = mixer(Bi, Bq1)
    Bed_data.append(Bed2input)
    Bj = Bed2input
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ BED 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for X in range(bed2.vectlen - 1):
        # setup new step
        dX = bed2.vect[X + 1] - bed2.vect[X]
        Bi = Bj + [dX, bed2.area]

        # run reactor step
        Bj = reactorStep(Bi)

        # print if necs
        # if bedvect[X] > lastPrintX:
        #    print('Reactor length = %2.1f, conversion = %2.2f' % (bedvect[X], Bi[2] / (2*initial_mol_N2) * 100) + '%' + ', T = %3.1f' %Bi[3])
        #    lastPrintX = np.ceil(bedvect[X]*10)/10

        # store data
        Bed_data.append(Bj)

    Bi = Bj
    print('Bed 2 length = %2.1fm, conversion = %2.2f' % (bed2.vect[-1], Bi[2] / (2 * Bed2input[1]) * 100) +
          '%' + ', T = %3.1f' % Bi[3] + 'K')
    print('Ammonia produced = %2.2f kg/h' % (Bi[2] * 3600 * 17.03 / 1000))

    # ~~~~~~~~~~~~~~~~~~~~~~~~ QUENCH 2-3 ~~~~~~~~~~~~~~~~~~~~~~~~~
    if bed3bool:
        Bq2 = [q2fact * 3 * initial_mol_N2, q2fact * initial_mol_N2, q2fact * Inlet_NH3_mol, Inlet_Temperature,
               Reactor_Pressure]
        Bj = mixer(Bi, Bq2)
        Bed_data.append(Bj)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ BED 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        for X in range(bed3.vectlen - 1):
            # setup new step
            dX = bed3.vect[X + 1] - bed3.vect[X]
            Bi = Bj + [dX, bed3.area]

            # run reactor step
            Bj = reactorStep(Bi)

            # print if necs
            # if bedvect[X] > lastPrintX:
            #    print('Reactor length = %2.1f, conversion = %2.2f' % (bedvect[X], Bi[2] / (2*initial_mol_N2) * 100) + '%' + ', T = %3.1f' %Bi[3])
            #    lastPrintX = np.ceil(bedvect[X]*10)/10

            # store data
            Bed_data.append(Bj)

        Bi = Bj
        print('Bed 3 length = %2.1fm, conversion = %2.2f' % (bed3.vect[-1], Bi[2] / (2 * initial_mol_N2 * (1 + q1fact + q2fact))
                                                             * 100) + '%' + ', T = %3.1f' % Bi[3] + 'K')

        print('Ammonia produced = %2.2f kg/h' % (Bi[2] * 3600 * 17.03 / 1000))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 1 (outlet to inlet) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('Pre HE Outlet from Reactor = ', Bi)
    [B_Sep_input, Bed1input, HE_1_Del_T_new] = heat_exchanger_hotgas2coldgas(Bi, Binlet)
    HE_1_Del_T_resid = abs(HE_1_Del_T - HE_1_Del_T_new)
    print('Post HE Outlet to condensor = %3.1f' % B_Sep_input[3], 'K\nPost HE Inlet into Reactor = %3.1f' % Bed1input[3], 'K\n\n HE1 del T = %3.1f' % HE_1_Del_T_new, '\n Resid = %3.1f' % HE_1_Del_T_resid)
    HE_1_Del_T = HE_1_Del_T_new
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 2 (outlet to water) ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    [Bj, Q] = heat_exchanger_water2gas(B_Sep_input,1)
    print('HE2 outlet T = %3.1f' % Bj[3])

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Condensor ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    [B_recycle, DelHact] = condensor(Bj,0.8)
    print('recycle stream composition = ', B_recycle, '\n Q coolant = %5.1f + %5.1f = %5.1f' % (Q, DelHact,Q+DelHact), ' [W] ')

    if HE_1_Del_T_resid > 100:
        stop = 1
    elif HE_1_Del_T < 200:
        stop = 1
    elif HE_1_Del_T_resid < Criterion:
        stop = 1

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Recycle ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reheat recycle stream? loosing too much energy rn

mass_inlet = (1 + q1fact + q2fact*bed3bool) * (Binlet[0] * 2 + Binlet[1] * 28 + Binlet[2] * 17)
mass_recycle = B_recycle[0] * 2 + B_recycle[1] * 28 + B_recycle[2] * 17
mol_inlet = (1 + q1fact + q2fact*bed3bool) * (Binlet[0] + Binlet[1] + Binlet[2])
mol_recycle = B_recycle[0] + B_recycle[1] + B_recycle[2]

recycle_ratio_mol = (mol_recycle / mol_inlet)
recycle_ratio_mass = (mass_recycle / mass_inlet)
print('recycle ratio mass = %2.3f' % recycle_ratio_mass, ', recycle ratio mol = %2.3f' % recycle_ratio_mol)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Heat exchanger 1 (pt 2?) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('power consumption', power_consumption, 'Q',  Q, 'delHact', DelHact)

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
