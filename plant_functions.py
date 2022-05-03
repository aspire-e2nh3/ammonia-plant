# this is a basic setup for an isobaric reactor.

# input should be a timeseries of N2, H2 flows?
# no power use yet - this will come?
# assumptions
# - isobaric - constant pressure in reactor
import copy
import math
import pyromat as pm
import numpy as np
import os

pm.config['unit_energy'] = 'J'  # default for pyromat is kJ


def sat_point_lookup(f="utils/sat_point_data.txt", T_sat=0, p_sat=0):
    """
    looks up saturation point from data stored in sat_point_data.txt or elsewhere
    :param f: saturation data in form T,p as columns
    :param T_sat: temperature to request partial pressure at
    :param p_sat: pressure to request temperature at
    :return:
    """
    with open(f,"r") as file:
        data = [[float(y) for y in x.strip("\n").split(",")] for x in file.readlines()]
    transposed_data = np.transpose(data)
    if T_sat != 0:
        return float(np.interp(T_sat, transposed_data[0], transposed_data[1]))
    elif p_sat != 0:
        return float(np.interp(p_sat, transposed_data[1], transposed_data[0]))
    else:
        return "Error: no interp value requested"


def heat_exchanger_parallel(s1, s2, effectiveness=0.8):
    '''
    Heat exchanger: e-NTU form.

    :param  s1: state of hot stream input
    :param  s2: state of cold stream input
    :param  effectiveness: efficiency?

    :return s1_out:
    '''
    s1.update_fast()
    s2.update_fast()

    s1_out = copy.copy(s1)
    s2_out = copy.copy(s2)

    Cmin = min(s1_out.cp * s1_out.mass_tot, s2_out.cp * s2_out.mass_tot)
    # Cmax = max(s1.cp*s1.mass_tot, s2.cp*s2.mass_tot)
    # Cr = Cmin / Cmax

    Q = effectiveness * (Cmin * (s1_out.T - s2_out.T))

    s1_out.T += - Q / (s1_out.cp * s1_out.mass_tot)
    s2_out.T += Q / (s2_out.cp * s2_out.mass_tot)


    s1_out.update_fast()
    s2_out.update_fast()
    return s1_out, s2_out, Q / (s2_out.cp * s2_out.mass_tot)


def heat_exchanger_counter(s1, s2, effectiveness=0.75, set_cold_out=False, T2out=273):  # mol,mol,mol,K,bar,K,m/s,mm check units!!
    """
    Heat exchanger: e-NTU form. defaults to eff =0.75 with no input
    :param s1: state of hot stream input
    :param s2: state of cold stream input
    :param T2out: desired temp of output stream
    :param effectiveness: effectiveness (if set directly)
    :return:
    """
    s1.update_fast()
    s2.update_fast()

    s1_out = copy.copy(s1)
    s2_out = copy.copy(s2)

    if set_cold_out:
        s2_out.T = effectiveness * s1.T + (1 - effectiveness) * s2.T
        s1_out.T = effectiveness * s2.T + (1 - effectiveness) * s1.T
    else:
        effectiveness = (T2out - s2.T) / (s1.T - s2.T)
        s2_out.T = T2out
        s1_out.T = effectiveness * s2.T + (1 - effectiveness) * s1.T

    Q = s2.cp * s2.mass_tot * (s2.T - s2_out.T)

    s1_out.update_fast()
    s2_out.update_fast()
    return s1_out, s2_out, -(s2.T - s2_out.T), effectiveness


def heat_exchanger_water2gas(s, T_end=0, cool_to_sat_point=False, effectiveness=0, water_mfr=10,  T_cin=273, Vmax=5, D=0.006):  # -,kg,K,m/s,m
    """
    saturation
    """
    s.update_fast()
    s_out = copy.copy(s)

    if cool_to_sat_point:
        T_end = sat_point_lookup(p_sat=s_out.yNH3 * s_out.p)
    elif T_end:
        pass
    else:
        return "Error: unknown required cooling temp"

    C_mix = s.cp * s.mass_tot  # J/kg/K * kg/s

    Q = C_mix * (s.T - T_end)

    C_cool = water_mfr * 4180  # 1 kg/s * 4.18 J/kg/K

    Cmin = min(C_cool, C_mix)
    Cmax = max(C_cool, C_mix)
    Cr = Cmin / Cmax

    if effectiveness == 0:
        effectiveness = Q / (Cmin * (s_out.T - T_cin))


    NTU = (1/(Cr-1)) * np.log((effectiveness - 1) / (Cr * effectiveness - 1))  #

    water_temp_out = T_cin + Q * effectiveness / (C_cool)

    U1 = 300  # W/m^2/K

    A1 = NTU * Cmin / U1

    Q_flow = s_out.mass_tot / s_out.rho  # m^3/s?
    num_p = 4 * Q_flow / (Vmax * np.pi * D ** 2)
    l1 = A1 / (num_p * np.pi * D)

    Re = s_out.rho * Vmax * D / s_out.mu
    F_fact = 0.316 * Re ** -0.25
    Del_P = F_fact * s_out.rho * Vmax ** 2 * l1 / 2 / D

    s_out.T = T_end
    s_out.p += - Del_P * 10 ** -5
    s_out.update_fast()
    return s_out, Q, effectiveness


def tristan_heat_exchanger(sh_in,sc_in,he):
    sh = copy.copy(sh_in)
    sc = copy.copy(sc_in)
    dq_hot2cold = np.full(he.ix, 0, float)
    dq_cold2ext = np.full(he.ix, 0, float)

    dqsumsaved = 0
    if he.last_run_bool:
        Th = np.linspace(sh.T, sh.T - he.last_run_eff*(sh.T - sc.T), he.ix+1)
        Tc = np.linspace(sc.T + he.last_run_eff*(sh.T - sc.T), sc.T, he.ix+1)
    else:
        Th = np.full(he.ix + 1, sh.T, float)
        Tc = np.full(he.ix + 1, sc.T, float)

    pie = 3.141592658
    dx = he.Length / he.ix  # length of finite segment [m]



    rstart = 0.1  # RELAX FACTOR =0.1 seems to work well
    rstep = rstart
    miniter = (1 - rstart) / rstep

    if he.last_run_bool:
        miniter = 0
        rstart = 1

    #MAIN COMPUTATIONAL LOOP
    j = 0
    stop = False
    while stop == False:
        j += 1
        if j > he.jints:
            stop=True





        #define relaxation factor for smooth convergence

        relax = min(rstart + j * rstep, he.max_relax)


        #FORWARD PASS ON MIXTURE SIDE OF HEAT EXCHANGER

        dq_hot2cold_sum = 0 #reset heat integration variable
        dq_cold2ext_sum = 0


        for ii in range(0, he.ix):

            sh.T = Th[ii]
            sh.update_fast()


            # Determine if mixture temp low enough for condensation to take place


            #Use law of mixtures to calculate average fluid properties and overall heat transfer coefficient for an element


            vel_h = sh.volume_fr / (he.numb * 2 * pie * he.r1 ** 2)
            rey_h = sh.rho * 2 * he.r1 * vel_h / sh.mu
            pr_h = sh.mu * sh.cp / sh.k
            nus_h = 0.023 * rey_h ** 0.8 * pr_h ** 0.4
            htc_h = nus_h * sh.k / (2 * he.r1)

            vel_c = sc.volume_fr / (he.numb * pie * (he.r3 ** 2 - he.r2 ** 2))
            rey_c = sc.rho * 2 * he.hyd * vel_c / sc.mu
            pr_c = sc.mu * sc.cp / sc.k
            nus_c = 0.023 * rey_c ** 0.8 * pr_c ** 0.4
            htc_c = nus_c * sc.k / he.hyd

            Uval = 1 / (1 / (htc_h * he.r1 * 2 * pie * dx) + np.log(he.r2 / he.r1) / (2 * pie * he.kval * dx) + 1 / (htc_c * he.r2 * 2 * pie * dx))

            Uval_ext = 2 * np.pi * dx/ (1 / (htc_c * he.r3)
                                        + np.log(he.r4 / he.r3) / he.kval
                                        + np.log(he.r5 / he.r4) / he.kval_insul
                                        + 1 / (he.htc_ext * he.r5))

            #determine heat flow from mixture to coolant in that element and heat lost to enviroment
            dq_hot2cold[ii] = Uval * (Th[ii] - Tc[ii]) * he.numb * relax
            dq_cold2ext[ii] = Uval_ext * (Tc[ii] - he.T_ext) * he.numb * relax

            dq_hot2cold_sum += dq_hot2cold[ii] #integrate total heat lost to coolant
            dq_cold2ext_sum += dq_cold2ext[ii]  # integrate total heat lost from coolant to external

            delT = dq_hot2cold[ii]/(sh.cp*sh.mass_tot)
            #reduce temperature of mixture for next element based on heat lost from this element

            Th[ii + 1] = Th[ii] - delT

        #COUNTER PASS ON COOLANT CHANNEL
        Tc[he.ix] = sc_in.T
        for ii in range(0, he.ix): #loop to update coolant temperatures
            Icon = he.ix - ii
            sc.T = Tc[Icon]
            sc.update_fast()
            Tc[Icon-1] = Tc[Icon] + (dq_hot2cold[Icon-1] + dq_cold2ext[Icon-1]) / (sc.cp*sc.mass_tot) #increase temp of coolant in line with heat flow from mixture

        #convergence check - has exit mole fraction changed from last main loop iteration
        if (abs(dq_hot2cold_sum - dqsumsaved) < he.eval) & (j > miniter): #if loop to check for convergence
            stop = True
        dqsumsaved = dq_hot2cold_sum
        #print(j, end=' ')
        #if j % 100 == 0:
        #    print('\n', end=' ')

    #print('\n', end=' ')
    he_eff = (Tc[0] - Tc[he.ix])/(Th[0]-Tc[he.ix])
    he.last_run_eff = he_eff


    #return(sh,sc,Th,Tc)
    return sh, sc, dq_hot2cold_sum, dq_cold2ext_sum, he


def tristan_condenser(s_in, c):


    #DECLARE ARRAYS
    Tmix = np.full(c.ix+1,s_in.T,float)
    Tcool = np.full(c.ix+1,c.Tcoolin,float)
    dq = np.full(c.ix,0,float)
    savxnh3 = np.full(c.ix,0,float)
    savnnh3 = np.full(c.ix,0,float)
    pie = 3.141592658


     #INITIALISING
    dx = c.Length / c.ix #length of finite segment [m]
    NH3saved = 0 #initialise convergence variable for ammonia concentration


    rstart = 0.1  # RELAX FACTOR =0.1 seems to work well
    rstep = rstart
    minint = (1 - rstart) / rstep

    #MAIN COMPUTATIONAL LOOP
    j = 0
    stop = False
    while stop == False:
        j += 1
        if j > c.jints:
            stop=True

        s = copy.copy(s_in)

        #define relaxation factor for smooth convergence

        relax = min(rstart + j * rstep, 1)


        #FORWARD PASS ON MIXTURE SIDE OF HEAT EXCHANGER

        dqsum = 0 #reset heat integration variable
        qflowsum = 0

        for ii in range(0, c.ix):

            s.T = Tmix[ii]
            s.update_fast()


            # Determine if mixture temp low enough for condensation to take place
            pnh3 = s.p * s.yNH3 #partial pressure of ammonia
            Tsat = c.abar * pnh3 ** c.bbar #determine saturation temperature of ammonia
            if Tmix[ii] > Tsat:
                conact = 0 #mixture temp above sat temp so no condensation in this element
            elif Tmix[ii] <= Tsat:
                conact = 1 #mixture temp below or equal to sat temp so condensation occuring in this element

            #Use law of mixtures to calculate average fluid properties and overall heat transfer coefficient for an element

            velmix = s.volume_fr / (c.numb * pie * c.r1 ** 2)
            reymix = s.rho * 2 * c.r1 * velmix / s.mu
            prmix = s.mu * s.cp / s.k
            nusmix = 0.023 * reymix ** 0.8 * prmix ** 0.4
            htc1 = nusmix * s.k / (2 * c.r1)

            Uval = 2 * np.pi * dx / (1 / (htc1 * c.r1) + np.log(c.r2 / c.r1) / c.kval + 1 / (c.htc2 * c.r2))

            #determine heat flow from mixture to coolant in that element
            dq[ii] = Uval * (Tmix[ii] - Tcool[ii]) * c.numb * relax
            dqsum = dqsum + dq[ii] #integrate total heat lost to coolant

            #determine molar flow of nh3 lost due to condensation in each element
            nlost = 0

            losint = s.NH3 / c.kints
            #increment condensation molar flow until sum of sensible heat reduction in gas + enthalpy of condensation equates to heat transferred out of element

            qout = 0
            # need to change to stop statement.

            dTdP = c.abar * c.bbar * (s.p * s.yNH3) ** (c.bbar - 1)
            delT = 0
            while qout<dq[ii]:
                nlost = nlost + losint
                dP = s.p * (s.yNH3 - (s.NH3 - nlost) / (s.mol_tot - nlost))
                delT = dTdP * dP
                qflow = (s.cp*s.mass_tot) * delT
                qout = qflow + conact * nlost * c.dhvap


            qflowsum += qflow
            s.NH3 -= conact * nlost  # reduce molar flow of gaseous ammonia as ammonia is condensed
            s.update()

            savnnh3[ii] = s.NH3
            savxnh3[ii] = s.yNH3
            #reduce temperature of mixture for next element based on heat lost from this element

            Tmix[ii + 1] = Tmix[ii] - delT
            #print(ii, dq[ii], qout, losint, nlost, qflow, delT)

        #COUNTER PASS ON COOLANT CHANNEL
        Tcool[c.ix] = c.Tcoolin
        for ii in range(0, c.ix): #loop to update coolant temperatures
            Icon = c.ix - ii
            Tcool[Icon-1] = Tcool[Icon] + dq[Icon-1] / (c.mcool * c.cpcool) #increase temp of coolant in line with heat flow from mixture


        #convergence check - has exit mole fraction changed from last main loop iteration
        if (abs(s.NH3 - NH3saved) < c.eval) & (j > minint): #if loop to check for convergence
            stop = True
        NH3saved = s.NH3

        #print(j, end=' ')
        #if j % 100 == 0:
        #    print('\n', end=' ')

    #print('\n', end=' ')
    #CALCULATE PRESSURE DROP IN EACH CHANNEL
    fricmix = 0.3164 * reymix ** -0.25
    delpmix = fricmix * s.rho * velmix ** 2 * c.Length / (2 * 2 * c.r1)
    s.T = Tmix[c.ix]
    s.p -= delpmix/1e5
    s.update_fast()


    friccool = 0.3164 * c.reycool ** -0.25
    delpcool = friccool * c.rcool * c.velcool ** 2 * c.Length / (2 * 2 * c.r1)


    #return Tmix,Tcool,savnnh3,s,delpmix,delpcool,dqsum
    return s,dqsum,Tcool[0]


def condenser(s, effectiveness=0.8, water_mfr=1, T_cin=10+273, Vmax=5, D=0.006):  ### check units!!

    s.update_fast()
    s_out = copy.deepcopy(s)

    Del_H_evap = 22.7 * 1000 * s_out.NH3  # J/mol * mol

    C_mix = s_out.cp * s_out.mass_tot

    C_cool = water_mfr * 4180  # 1 kg/s * 4180 J/kg/K

    Cmin = min(C_cool, C_mix)
    Cmax = max(C_cool, C_mix)
    Cr = Cmin / Cmax

    NTU = -np.log(1 - effectiveness)

    T_cout = s.T - (s.T - T_cin)*np.exp(-NTU)

    X = min(1,C_cool * (T_cout - T_cin)/Del_H_evap)

    #pressure drop calcs
    U2 = 750  # W/m^2/K
    A2 = NTU * Cmin / U2
    Q_flow = s_out.mass_tot / s_out.rho  # m^3/s?
    num_p = 4 * Q_flow / (Vmax * np.pi * D ** 2)
    l1 = A2 / (num_p * np.pi * D)
    Re = s_out.rho * Vmax * D / s_out.mu
    F_fact = 0.316 * Re ** -0.25
    Del_P = F_fact * s_out.rho * Vmax ** 2 * l1 / 2 / D


    s_out.p += - Del_P * 10 ** -5
    s_out.NH3 = s_out.NH3 * (1.00001-X)
    s_out.update_fast()
    return s_out, X, Del_H_evap*X,T_cout


def condenser_crude(s, water_mass_flow=1, T_cin=10+273):
    initial_pp = s.p * s.yNH3
    C_cool = water_mass_flow * 4180  # 1 kg/s * 4180 J/kg/K
    T_cout = T_cin

    s_out = copy.copy(s)

    criterion = 0.01
    stop = 0
    count = 0
    while stop == 0:
        final_pp = sat_point_lookup(T_sat=T_cout)
        ammonia_removed = (1-final_pp/initial_pp)


        #condensor guess

        s_out.NH3 = s.NH3 * (1-ammonia_removed)
        Q_latent = s.NH3 * ammonia_removed * 22.7 * 1000
        s_out.T = T_cout
        s_out.update()
        T_cout_new = T_cin + (Q_latent - (s_out.cp*s.mass_tot*s_out.T) + (s.cp*s.mass_tot*s.T))/C_cool

        if abs(T_cout_new - T_cout) < criterion:
            stop = 1
        T_cout = T_cout_new


    power = C_cool * (T_cout - T_cin)
    s_out.update_fast()
    return s_out, power, ammonia_removed, T_cout


def reactor(s_in, b):  # mol/s, K, Pa
    """
    An iterative function to determine the change in state variables and reactants over the length of a reactor Bed step

    :param  s: state class. Should include at least some NH3 or reaction will shoot up
            dX: X step.
            area: area of reactor bed.

    :return s: state class
    """
    s = copy.copy(s_in)
    exotherm_tot = 0
    heat_loss_tot = 0
    bed_data = [s.store()]

    for X in range(b.vectlen - 1):
        # setup new step
        dX = b.vect[X + 1] - b.vect[X]
        # [s_bed_temp, exotherm] = reactorStep(s_bed_temp, dX, b.cs_area)

        # initial concentrations
        T = s.T
        P = s.p

        # activity coefficients for all species
        N2fuga = 0.93431737 + 0.3101804 * 10 ** -3 * T + 0.295895 * 10 ** -3 * P - 0.270729 * 10 ** -6 * T ** 2 + \
                 0.4775207 * 10 ** -6 * P ** 2
        H2fuga = math.exp(math.exp(-3.8402 * T ** 0.125 + 0.541) * P - math.exp(-0.1263 * T ** 0.5 - 15.980) * P ** 2 +
                          300 * math.exp(-0.011901 * T - 5.941) * (math.exp(-P / 300) - 1))
        NH3fuga = 0.1438996 + 0.2028538 * 10 ** -2 * T - 0.4487672 * 10 ** -3 * P - 0.1142945 * 10 ** -5 * T ** 2 + \
                  0.2761216 * 10 ** -6 * P ** 2
        a_N2 = s.yN2 * N2fuga * P  # bar
        a_H2 = s.yH2 * H2fuga * P  # bar
        a_NH3 = s.yNH3 * NH3fuga * P  # bar

        # reaction rate constant
        Ea = 1.7056 * 10 ** 5  # J/mol
        R = 8.31446261815324  # J/K/mol
        k0 = 8.8490 * 10 ** 17/3600  # mol/m^3/s
        k_r = k0 * math.exp(-Ea / (R * T))  # mol/m^3/s ?

        # equilibrium constant
        K_eq = math.pow(10, -2.691122 * math.log(T, 10) - 5.519265 * 10 ** -5 * T + 1.848863 * 10 ** -7 * T ** 2
                        + 2001.6 / T + 2.67899)

        # Reaction Rate per unit volume_fr
        alpha = 0.5
        RR_NH3_plus = 2 * k_r * (K_eq ** 2 * a_N2 * (a_H2 ** 3 / (a_NH3 ** 2)) ** alpha)
        RR_NH3_minus = 2 * k_r * (- (a_NH3 ** 2 / (a_H2 ** 3)) ** (1 - alpha))
        RR_NH3 = RR_NH3_plus + RR_NH3_minus

        RR_N2 = -1 / 2 * RR_NH3
        RR_H2 = -3 / 2 * RR_NH3

        # Heat of Reaction
        Del_H = 4.184 * (-(0.54526 + 846.609 / T + 459.734 * 10 ** 6 * T ** -3) * P - 5.34685 * T
                         - 0.2525 * 10 ** -3 * T ** 2 + 1.69197 * 10 ** -6 * T ** 3 - 9157.09)  # J/mol

        # change in molars
        s.N2 += dX * RR_N2 * b.cs_area
        s.H2 += dX * RR_H2 * b.cs_area
        s.NH3 += dX * RR_NH3 * b.cs_area

        # change in temp
        exotherm = -dX * b.cs_area * Del_H * RR_NH3

        vel = s.volume_fr / b.cs_area
        rey = s.rho * 2 * b.r * vel / s.mu
        pr = s.mu * s.cp / s.k
        nus = 0.023 * rey ** 0.8 * pr ** 0.4
        htc = nus * s.k / b.r

        Uval = 2 * np.pi * dX / (1 / (htc * b.r)
                                 + np.log((b.r + b.thickness) / b.r) / b.shell_kval
                                 + np.log((b.r + b.thickness + b.insul_thickness)/(b.r + b.thickness))/b.insul_kval
                                 + 1 / (b.external_htc * (b.r + b.thickness)))

        heat_loss = Uval * (s.T - b.coolant_T)

        s.T += (exotherm - heat_loss) / (s.cp * s.mass_tot)

        # update state
        s.update_fast()

        bed_data.append(s.store())
        exotherm_tot += exotherm
        heat_loss_tot += heat_loss

    return s, exotherm_tot, heat_loss_tot, bed_data


def reactor_coolant(s_in, b):  # mol/s, K, Pa
    """
    An iterative function to determine the change in state variables and reactants over the length of a reactor Bed step

    :param  s: state class. Should include at least some NH3 or reaction will shoot up
            dX: X step.
            area: area of reactor bed.

    :return s: state class
    """
    s = copy.copy(s_in)
    exotherm_tot = 0
    heat_loss_tot = 0
    bed_data = [s.store()]

    for X in range(b.vectlen - 1):
        # setup new step
        dX = b.vect[X + 1] - b.vect[X]
        # [s_bed_temp, exotherm] = reactorStep(s_bed_temp, dX, b.cs_area)

        # initial concentrations
        T = s.T
        P = s.p

        # activity coefficients for all species
        N2fuga = 0.93431737 + 0.3101804 * 10 ** -3 * T + 0.295895 * 10 ** -3 * P - 0.270729 * 10 ** -6 * T ** 2 + \
                 0.4775207 * 10 ** -6 * P ** 2
        H2fuga = math.exp(math.exp(-3.8402 * T ** 0.125 + 0.541) * P - math.exp(-0.1263 * T ** 0.5 - 15.980) * P ** 2 +
                          300 * math.exp(-0.011901 * T - 5.941) * (math.exp(-P / 300) - 1))
        NH3fuga = 0.1438996 + 0.2028538 * 10 ** -2 * T - 0.4487672 * 10 ** -3 * P - 0.1142945 * 10 ** -5 * T ** 2 + \
                  0.2761216 * 10 ** -6 * P ** 2
        a_N2 = s.yN2 * N2fuga * P  # bar
        a_H2 = s.yH2 * H2fuga * P  # bar
        a_NH3 = s.yNH3 * NH3fuga * P  # bar

        # reaction rate constant
        Ea = 1.7056 * 10 ** 5  # J/mol
        R = 8.31446261815324  # J/K/mol
        k0 = 8.8490 * 10 ** 17  # mol/m^3
        k_r = k0 * math.exp(-Ea / (R * T))  # mol/m^3/h

        # equilibrium constant
        K_eq = math.pow(10, -2.691122 * math.log(T, 10) - 5.519265 * 10 ** -5 * T + 1.848863 * 10 ** -7 * T ** 2
                        + 2001.6 / T + 2.67899)

        # Reaction Rate per unit volume_fr
        alpha = 0.5
        RR_NH3 = k_r * (K_eq ** 2 * a_N2 * (a_H2 ** 3 / (a_NH3 ** 2)) ** alpha
                        - (a_NH3 ** 2 / (a_H2 ** 3)) ** (1 - alpha)) / 3600
        RR_N2 = -1 / 2 * RR_NH3
        RR_H2 = -3 / 2 * RR_NH3

        # Heat of Reaction
        Del_H = 4.184 * (-(0.54526 + 846.609 / T + 459.734 * 10 ** 6 * T ** -3) * P - 5.34685 * T
                         - 0.2525 * 10 ** -3 * T ** 2 + 1.69197 * 10 ** -6 * T ** 3 - 9157.09)  # J/mol

        # change in molars
        s.N2 += dX * RR_N2 * b.cs_area
        s.H2 += dX * RR_H2 * b.cs_area
        s.NH3 += dX * RR_NH3 * b.cs_area

        # change in temp
        exotherm = -dX * b.cs_area * Del_H * RR_NH3

        vel = s.volume_fr / b.cs_area
        rey = s.rho * 2 * b.r * vel / s.mu
        pr = s.mu * s.cp / s.k
        nus = 0.023 * rey ** 0.8 * pr ** 0.4
        htc = nus * s.k / b.r

        vel_c = (b.coolant_mfr / b.coolant_rho) / (2 * b.coolant_channel * b.length)
        rey_c = b.coolant_rho * np.pi * (b.r + b.thickness) * vel_c / b.coolant_mu
        pr_c = b.coolant_mu * b.coolant_cp / b.coolant_k
        nus_c = 0.023 * rey_c ** 0.8 * pr_c ** 0.4
        htc_c = nus_c * b.coolant_k / (b.r + b.thickness)

        Uval = 2 * np.pi * dX / (1 / (htc * b.in_circum)
                                 + np.log((b.r + b.thickness) / b.r) / b.shell_kval
                                 + 1 / (htc_c * (b.r + b.thickness)))

        heat_loss = Uval * (s.T - b.coolant_T)

        s.T += (exotherm - heat_loss) / (s.cp * s.mass_tot)

        # update state
        s.update_fast()

        bed_data.append(s.store())
        exotherm_tot += exotherm
        heat_loss_tot += heat_loss

    return s, exotherm_tot, heat_loss_tot, bed_data


def reactor_selfcooled(s_in, b):  # mol/s, K, Pa
    """
    An iterative function to determine the change in state variables and reactants over the length of a reactor Bed step

    :param  s: state class. Should include at least some NH3 or reaction will shoot up
            dX: X step.
            area: area of reactor bed.

    :return s: state class
    """
    s = copy.copy(s_in)

    ext_bed = np.full(273,0,float)
    int_bed = np.full(273,0,float)


    exotherm_tot = 0
    heat_loss_tot = 0
    bed_data = [s.store()]

    for X in range(b.vectlen - 1):
        # setup new step
        dX = b.vect[X + 1] - b.vect[X]
        # [s_bed_temp, exotherm] = reactorStep(s_bed_temp, dX, b.cs_area)

        # initial concentrations
        T = s.T
        P = s.p

        # activity coefficients for all species
        N2fuga = 0.93431737 + 0.3101804 * 10 ** -3 * T + 0.295895 * 10 ** -3 * P - 0.270729 * 10 ** -6 * T ** 2 + \
                 0.4775207 * 10 ** -6 * P ** 2
        H2fuga = math.exp(math.exp(-3.8402 * T ** 0.125 + 0.541) * P - math.exp(-0.1263 * T ** 0.5 - 15.980) * P ** 2 +
                          300 * math.exp(-0.011901 * T - 5.941) * (math.exp(-P / 300) - 1))
        NH3fuga = 0.1438996 + 0.2028538 * 10 ** -2 * T - 0.4487672 * 10 ** -3 * P - 0.1142945 * 10 ** -5 * T ** 2 + \
                  0.2761216 * 10 ** -6 * P ** 2
        a_N2 = s.yN2 * N2fuga * P  # bar
        a_H2 = s.yH2 * H2fuga * P  # bar
        a_NH3 = s.yNH3 * NH3fuga * P  # bar

        # reaction rate constant
        Ea = 1.7056 * 10 ** 5  # J/mol
        R = 8.31446261815324  # J/K/mol
        k0 = 8.8490 * 10 ** 17  # mol/m^3
        k_r = k0 * math.exp(-Ea / (R * T))  # mol/m^3/h

        # equilibrium constant
        K_eq = math.pow(10, -2.691122 * math.log(T, 10) - 5.519265 * 10 ** -5 * T + 1.848863 * 10 ** -7 * T ** 2
                        + 2001.6 / T + 2.67899)

        # Reaction Rate per unit volume_fr
        alpha = 0.5
        RR_NH3 = k_r * (K_eq ** 2 * a_N2 * (a_H2 ** 3 / (a_NH3 ** 2)) ** alpha
                        - (a_NH3 ** 2 / (a_H2 ** 3)) ** (1 - alpha)) / 3600
        RR_N2 = -1 / 2 * RR_NH3
        RR_H2 = -3 / 2 * RR_NH3

        # Heat of Reaction
        Del_H = 4.184 * (-(0.54526 + 846.609 / T + 459.734 * 10 ** 6 * T ** -3) * P - 5.34685 * T
                         - 0.2525 * 10 ** -3 * T ** 2 + 1.69197 * 10 ** -6 * T ** 3 - 9157.09)  # J/mol

        # change in molars
        s.N2 += dX * RR_N2 * b.cs_area
        s.H2 += dX * RR_H2 * b.cs_area
        s.NH3 += dX * RR_NH3 * b.cs_area

        # change in temp
        exotherm = -dX * b.cs_area * Del_H * RR_NH3

        vel = s.volume_fr / b.cs_area
        rey = s.rho * 2 * b.r * vel / s.mu
        pr = s.mu * s.cp / s.k
        nus = 0.023 * rey ** 0.8 * pr ** 0.4
        htc = nus * s.k / b.r

        vel_c = (b.coolant_mfr / b.coolant_rho) / (2 * b.coolant_channel * b.length)
        rey_c = b.coolant_rho * np.pi * (b.r + b.thickness) * vel_c / b.coolant_mu
        pr_c = b.coolant_mu * b.coolant_cp / b.coolant_k
        nus_c = 0.023 * rey_c ** 0.8 * pr_c ** 0.4
        htc_c = nus_c * b.coolant_k / (b.r + b.thickness)

        Uval = 2 * np.pi * dX / (1 / (htc * b.in_circum)
                                 + np.log((b.r + b.thickness) / b.r) / b.shell_kval
                                 + 1 / (htc_c * (b.r + b.thickness)))

        heat_loss = Uval * (s.T - b.coolant_T)

        s.T += (exotherm - heat_loss) / (s.cp * s.mass_tot)

        # update state
        s.update_fast()

        bed_data.append(s.store())
        exotherm_tot += exotherm
        heat_loss_tot += heat_loss

    return s, exotherm_tot, heat_loss_tot, bed_data


def mixer(s1, s2):
    """
    Function to mix together two pipe streams into one.
    :param  s1: input state 1
            s2: input state 2
    :return s_out: output state, mixed, at lowest pressure


    """
    # Tav = (s1.cp * s1.T * s1.mass_tot + s2.cp * s2.T * s2.mass_tot) \
    #       / (((s1.cp + s2.cp) / 2) * (s1.mass_tot + s2.mass_tot))
    Tav = (s1.T*s1.mass_tot + s2.T*s2.mass_tot)/(s1.mass_tot + s2.mass_tot)


    s_out = State(s1.H2 + s2.H2, s1.N2 + s2.N2, s1.NH3 + s2.NH3, Tav, min(s1.p, s2.p))

    s_out.T = (s1.cp * s1.T * s1.mass_tot + s2.cp * s2.T * s2.mass_tot) / (s_out.cp * s_out.mass_tot)

    s_out.update_fast()
    return s_out


def compressor(s, p_out, eta=0.7):
    """
    Function to return rate of work for a compressor based on a
    target pressure.

    Inputs: INPUT - 5x1 list of [float], standard input of mol H2, mol N2, mol NH3, T[K] and p[bar]
            eta - float, efficiency of compressor

    Outputs: OUTPUT - 5x1 list of [float], standard output of mol H2, mol N2, mol NH3, T[K] and p[bar]
             w - float, output
    """
    s.update_fast()
    s_out = copy.deepcopy(s)

    y = s_out.gamma
    r_p = p_out / s.p
    a = (y - 1) / y
    power = s.cp * s.mass_tot*s.T / eta * (r_p ** a - 1)
    s_out.T = s.T * (1 + r_p ** a / eta - 1 / eta)
    s_out.p = p_out
    s_out.update_fast()
    return s_out, power,0


def ptcompressor(s, p_out, t_out=0, eta=0.7):
    """
    Function to return rate of work for a polytropic compressor based on a
    target pressure and outlet temperature.

    Inputs: INPUT - 5x1 list of [float], standard input of mol H2, mol N2, mol NH3, T[K] and p[bar]
            eta - float, efficiency of compressor

    Outputs: OUTPUT - 5x1 list of [float], standard output of mol H2, mol N2, mol NH3, T[K] and p[bar]
             w - float, output
    """
    n = 1/(1-(math.log(t_out/s.T)/math.log(p_out/s.p)))
    a = n/(n-1)


    s.update_fast()
    s_out = copy.deepcopy(s)

    y = s_out.gamma
    R = (y-1)*s_out.cp
    power = s_out.mass_tot*R*a*(t_out - s_out.T)
    Q_out = s_out.mass_tot*R*a*(t_out - s_out.T) - s_out.mass_tot*s_out.cp*(t_out - s_out.T)
    s_out.T = t_out
    s_out.p = p_out
    s_out.update_fast()
    return s_out, power, Q_out, n


def psa_estimate(N2_mol, p_out=10, eta=0.7):
    """
    Function to return power for PSA to separate given mol of N2.

    Inputs:     N2_mol - mol flow rate of N2
                Pout - float, target pressure
                eta - float, efficiency of compressor

    Outputs:    OUTPUT - 5x1 list of [float], standard output of mol H2, mol N2, mol NH3, T[K] and p[bar]
                w - float, work rate required (W) for PSA
    """
    [N2in, Tin, Pin] = [N2_mol,298,1]
    s_out = State(0, N2in, 0, Tin, Pin)
    [s_out,power,_,_] = ptcompressor(s_out,10,373)
    return s_out, power


def electrolysis(H2mol, eta=0.65):  # mol/s to W
    """
        Function to return power for PEM electrolysis to separate given mol of H2.

        Inputs:     H2 mol input + output

        Outputs:    OUTPUT - 5x1 list of [float], standard output of mol H2, mol N2, mol NH3, T[K] and p[bar]
                    w - float, work rate required (W) for PSA
                    H20_in - required water to produce hydrogen
        """
    W = 241.83 / eta * H2mol * 1000
    s_out = State(H2mol, 0, 0, 298, 10)
    H20_in = H2mol * 18.015 / 2.016
    return s_out, W


def heater(s, T_end, no_cooling = True):
    if (s.T > T_end) & no_cooling:
        return s, 0
    s_out = copy.copy(s)
    power = s_out.cp * s_out.mass_tot * (T_end - s_out.T)
    s_out.T = T_end
    s_out.update_fast()
    return s_out, power


def heater_power(s,P):
    s_out = copy.copy(s)
    s_out.T += P/(s_out.mass_tot * s_out.cp)
    s_out.update_fast()
    return s_out


class Bed(object):
    '''
    A class to store properties of a 1D Bed model for a reactor.
    also generates a vector to describe end
    '''
    def __init__(self, length, diam, mini, log_divs=19,lin_divs=100):
        self.length = length
        self.diam = diam
        self.r = diam/2
        self.cs_area = math.pi * self.r ** 2  # m^2
        self.in_circum = math.pi * self.r * 2
        self.surrounding_T = 298
        # generate vect

        self.vect = [0] + np.geomspace(mini, length / lin_divs, log_divs, endpoint=False).tolist() + \
                    np.linspace(length / lin_divs, length, lin_divs-1).tolist()

        self.vectlen = len(self.vect)


        self.strength = 148
        self.sf = 2
        self.pressure = 200
        self.shell_density = 7700
        self.cat_density = 7900/2
        self.stress = self.strength * 10 ** 6 / self.sf
        self.thickness = self.pressure * 10 ** 5 * self.r / self.stress
        self.OD = diam + self.thickness * 2
        self.out_circum = math.pi * self.OD
        self.shell_volume = self.thickness * math.pi * (self.length * 2 * self.r + self.r ** 2 * 4 / 3)
        self.cat_volume = self.r ** 2 * self.length * math.pi
        self.shell_mass = self.shell_volume * self.shell_density
        self.cat_mass = self.cat_volume * self.cat_density
        self.shell_kval = 16
        self.insul_kval = 0.05
        self.insul_thickness = 0.05
        self.external_htc = 12

        self.coolant_mfr = 0.03
        self.coolant_rho = 680.2 # kg/m3
        self.coolant_channel = 0.005
        self.coolant_mu = 0.13*0.001 # Pa/s
        self.coolant_cp = 2701 # J/kg/K
        self.coolant_k = 0.0779 # W/m/K
        self.coolant_T = 273 # K


class State(object):
    '''
    Class to store variables about mass and molar flow rate, temperature, pressure etc of data.

    :param H2: molar flow rate [mol/s] of H2 gas.
    :param N2: molar flow rate [mol/s] of N2 gas.
    :param NH3: molar flow rate [mol/s] of NH3 gas.
    :param T: Temperature of flow [K].
    :param P: Pressure of flow [bar].
    '''

    def __init__(self, H2, N2, NH3, T, p):
        self.H2 = H2
        self.N2 = N2
        self.NH3 = NH3
        self.T = T
        self.p = p
        if H2 + N2 + NH3 > 0:
            self.update_fast()
        else:
            self.mol_tot = 0


    def update(self):
        # masses of components in [kg]
        self.mH2 = self.H2 * 2.016 / 1000
        self.mN2 = self.N2 * 28.0134 / 1000
        self.mNH3 = self.NH3 * 17.03 / 1000

        # tot mass
        self.mass_tot = self.mH2 + self.mN2 + self.mNH3

        # total mol
        self.mol_tot = self.H2 + self.N2 + self.NH3

        # mol
        self.yH2 = self.H2 / self.mol_tot
        self.yN2 = self.N2 / self.mol_tot
        self.yNH3 = self.NH3 / self.mol_tot

        #mass frac
        self.myH2 = self.mH2 / self.mass_tot
        self.myN2 = self.mN2 / self.mass_tot
        self.myNH3 = self.mNH3 / self.mass_tot

    def update_fast(self):
        self.update()

        # cp [J/kg/K]
        self.cp = self.myH2*14615 + self.myN2*1113 + self.myNH3*2600
        # mu [Ns/m^2]
        #self.mu = 0.000009 * self.myH2 + 0.000018 * self.myN2 + 0.00001 * self.myNH3

        # mu from sutherland's formula [Ns/m^2]
        mu_N2 = 1.663 * 10 ** -5 * (self.T / 273) ** (3 / 2) * (273 + 107) / (self.T + 107)
        mu_H2 = 8.411 * 10 ** -5 * (self.T / 273) ** (3 / 2) * (273 + 97) / (self.T + 97)
        mu_NH3 = 0.919 * 10 ** -5 * (self.T / 273) ** (3 / 2) * (273 + 370) / (self.T + 370)
        self.mu = self.myN2 * mu_N2 + self.myH2 * mu_H2 + self.myNH3 * mu_NH3

        # k [W/K/m]
        #self.k = 0.2 * self.yH2 + 0.028 * self.yN2 + 0.026 * self.yNH3

        # laminar thermal conductivity from data from engineering toolbox [W/K/m]
        self.k = (0.17 + (self.T - 250) * (0.24 - 0.17) / (400 - 250)) * self.myH2 + \
                 (0.0311 + (self.T - 300) * (0.047 - 0.0311) / (600 - 300)) * self.myN2 + \
                 (0.0334 + (self.T - 325) * (0.047 - 0.0311) / (725 - 325)) * self.myNH3

        #rho [kg/m^3]
        #self.rho = float(1/(self.yH2*(4126*self.T)/(self.p*1e5) + self.yN2*(297*self.T)/(self.p*1e5) + self.yNH3/pm.get('ig.NH3').d(self.T, self.p)))

        self.rho = self.p*1e5 / ( self.myH2 *(4126 * self.T) +
                               self.myN2 * (297 * self.T) +
                               self.myNH3 * (488 * self.T))

        #gamma [-]
        self.gamma = (1 + 1 / (self.yN2 / 0.4 + self.yH2 / 0.4 + self.yNH3 / 0.31))

        # volume_fr [m^3
        self.volume_fr = self.mass_tot / self.rho


    def update_slow(self):
        self.update()

        self.cp = float(pm.get('ig.N2').cp(self.T, self.p) * self.myN2 + pm.get('ig.H2').cp(self.T, self.p) * self.myH2 +
                         pm.get('ig.NH3').cp(self.T, self.p) * self.myNH3)
        # rho [kg/m^3]
        self.rho = float( 1 / (self.myN2 / pm.get('ig.N2').d(self.T, self.p) +
                               self.myH2 / pm.get('ig.H2').d(self.T, self.p) +
                               self.myNH3 / pm.get('ig.NH3').d(self.T, self.p)))
        # volume_fr [m^3
        self.volume_fr = self.mass_tot / self.rho

        # k from data from engineering toolbox [W/K/m]
        self.k = (0.17 + (self.T-250)*(0.24-0.17)/(400-250)) * self.myH2 + \
                 (0.0311 + (self.T-300)*(0.047-0.0311)/(600-300)) * self.myN2 + \
                 (0.0334 + (self.T-325)*(0.047-0.0311)/(725-325)) * self.myNH3

        # gamma from addition of specific heats(1+ 1/sum of mol%/gamma)[-]
        self.gamma = float((1 + 1 / ((self.yN2 / 0.4) + (self.yH2 / 0.4) + self.yNH3 / (pm.get('ig.NH3').gam(self.T, self.p) - 1))))

        mu_N2 = 1.663 * 10 ** -5 * (self.T / 273) ** (3 / 2) * (273 + 107) / (self.T + 107)
        mu_H2 = 8.411 * 10 ** -5 * (self.T / 273) ** (3 / 2) * (273 + 97) / (self.T + 97)
        mu_NH3 = 0.919 * 10 ** -5 * (self.T / 273) ** (3 / 2) * (273 + 370) / (self.T + 370)
        self.mu = self.myN2 * mu_N2 + self.myH2 * mu_H2 + self.myNH3 * mu_NH3

    def store(self):
        return [self.H2, self.N2, self.NH3, self.T, self.p]

    def split(self,f):
        s1 = State(self.H2 * f, self.N2 * f, self.NH3 * f, self.T, self.p)
        s2 = State(self.H2 * (1 - f), self.N2 * (1 - f), self.NH3 * (1 - f), self.T, self.p)
        return s1, s2

    def __copy__(self):
        return State(self.H2,self.N2,self.NH3,self.T,self.p)


class Condenser_Details(object):

    def __init__(self,mass_flow=0.3):
        pie = 3.141592658
        self.r1 = 0.002 #inner rad of inner ammonia pipe [m]
        self.r2 = 0.003 #outer rad of inner pipe [m]
        self.r3 = 0.02 #inner rad of outer coolant pipe [m]
        self.hyd = 2 * (self.r3 - self.r2) #hydraulic diameter of cooling channel [m]
        self.Length = 5 #length of heat exchanger [m]
        self.numb = 100 #number of counterflow heat exchangers
        self.ix = 500 #number of elements along heat exchanger
        self.jints = 100 #max number of iterations
        self.kval = 16 #thermal conductivity of pipe [W/mK]

        self.dhvap = 23370 #average value for enthalpy of condensation [J/mol]
        self.mcool = mass_flow #mass flow of coolant [kg/s]
        self.cpcool = 4186 #heat capacity of coolant [kg/s]
        self.Tcoolin = 273 #coolant inlet temp [K]
        self.rcool = 1000 #density of coolant [kg/m3]
        self.kcool = 0.6 #thermal conductivity of coolant [W/mK]
        self.vcool = 0.001 #dynamic viscosity of coolant [Pas]
        self.velcool = self.mcool / (pie * (self.r3 ** 2 - self.r2 ** 2)) / self.numb
        self.reycool = self.rcool * self.velcool * self.hyd / self.vcool
        self.prcool = self.vcool * self.cpcool / self.kcool
        self.nuscool = 0.023 * self.reycool ** 0.8 * self.prcool ** 0.4
        self.htc2 = self.nuscool * self.kcool / self.hyd
        self.abar = 239.69 #constant for sat curve of nh3
        self.bbar = 0.0964 #constant for sat curve of nh3
        self.kints = 10000
        self.eval = 1e-5


class Heat_Exchanger_Details(object):

    def __init__(self):
        pie = 3.141592658
        self.r1 = 0.002 #inner rad of inner ammonia pipe [m]
        self.r2 = 0.003 #outer rad of inner pipe [m]
        self.r3 = 0.0035 #inner rad of outer coolant pipe [m]
        self.hyd = 2 * (self.r3 - self.r2) #hydraulic diameter of cooling channel [m]
        self.Length = 3 #length of heat exchanger [m]
        self.numb = 6 #number of counterflow heat exchangers
        self.ix = 500 #number of elements along heat exchanger
        self.jints = 100 #max number of iterations
        self.kval = 16 #thermal conductivity of pipe [W/mK]
        self.last_run_bool = True
        self.last_run_eff = 0.8
        self.max_relax = 1
        self.eval = 1e-5
        self.r4 = 0.0045
        self.r5 = 0.01
        self.htc_ext = 12
        self.kval_insul = 0.05
        self.T_ext = 298



def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)

def tps_lst(list_of_lists):
    """"""
    array = np.array(list_of_lists)
    transpose = array.T
    transpose_list = transpose.tolist()

    return transpose_list