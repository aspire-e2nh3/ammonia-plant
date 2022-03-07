# this is a basic setup for an isobaric reactor.

# input should be a timeseries of N2, H2 flows?
# no power use yet - this will come?
# assumptions
# - isobaric - constant pressure in reactor
import copy
import math
import pyromat as pm
import numpy as np

pm.config['unit_energy'] = 'J'  # default for pyromat is kJ


def sat_point_lookup(f="sat_point_data.txt", T_sat=0, p_sat=0):
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

def BedBlock(s,b,log=False):
    s_bed_temp = copy.copy(s)
    exotherm_tot = 0
    heat_loss_tot = 0
    if log:
        bed_data = [s_bed_temp.store()]
    for X in range(b.vectlen - 1):
        # setup new step
        dX = b.vect[X + 1] - b.vect[X]
        [s_bed_temp, exotherm] = reactorStep(s_bed_temp, dX, b.cs_area)
        if log:
            bed_data.append(s_bed_temp.store())
        exotherm_tot += exotherm
        heat_loss_tot += b.htc * (s_bed_temp.T - b.surrounding_T) * dX * b.circum

    if log:
        return s_bed_temp,exotherm_tot,heat_loss_tot, bed_data
    else:
        return s_bed_temp,exotherm_tot,heat_loss_tot,False

def reactorStep(s, dX, area):  # mol/s, K, Pa
    """
    An iterative function to determine the change in state variables and reactants over the length of a reactor Bed step

    :param  s: state class. Should include at least some NH3 or reaction will shoot up
            dX: X step.
            area: area of reactor bed.

    :return s: state class
    """

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

    # Reaction Rate per unit volume
    alpha = 0.5
    RR_NH3 = k_r * (K_eq ** 2 * a_N2 * (a_H2 ** 3 / (a_NH3 ** 2)) ** alpha
                    - (a_NH3 ** 2 / (a_H2 ** 3)) ** (1 - alpha)) / 3600
    RR_N2 = -1 / 2 * RR_NH3
    RR_H2 = -3 / 2 * RR_NH3

    # Heat of Reaction
    Del_H = 4.184 * (-(0.54526 + 846.609 / T + 459.734 * 10 ** 6 * T ** -3) * P - 5.34685 * T
                     - 0.2525 * 10 ** -3 * T ** 2 + 1.69197 * 10 ** -6 * T ** 3 - 9157.09)  # J/mol

    # change in molars
    s.N2 += dX * RR_N2 * area
    s.H2 += dX * RR_H2 * area
    s.NH3 += dX * RR_NH3 * area

    # change in temp
    exotherm = -dX * area * Del_H * RR_NH3
    s.T += exotherm / (s.cp * s.mass_tot)

    # update state
    s.update_fast()
    return s, exotherm


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


def compressor(s, p_out, t_out=0, eta=0.7):
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
    return s_out, power, 0

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
    return s_out, W, H20_in


def heater(s, T_end):
    if s.T>T_end:
        return s, 0
    s_out = copy.copy(s)
    power = s_out.cp * s_out.mass_tot * (T_end - s_out.T)
    s_out.T = T_end
    s_out.update_fast()
    return s_out, power


class Bed(object):
    '''
    A class to store properties of a 1D Bed model for a reactor.
    also generates a vector to describe end
    '''
    def __init__(self, length, diam, mini, newvectmethod=True, divs=19):
        self.length = length
        self.diam = diam
        self.r = diam/2
        self.cs_area = math.pi * self.r ** 2  # m^2
        self.circum = math.pi * self.r * 2
        self.htc = 1
        self.surrounding_T = 298
        # generate vect
        if not newvectmethod:
            self.vect = [0] + np.geomspace(mini, length / 10, divs, endpoint=False).tolist() + \
                        np.linspace(length / 10, length, 10).tolist()
        else:
            self.vect = [0] + np.geomspace(mini, length / 100, divs, endpoint=False).tolist() + \
                        np.linspace(length / 100, length, 100).tolist()

        self.vectlen = len(self.vect)

    def mass(self,pressure=200, stress=148, sf=2, shell_density=7700, cat_density=7900):
        '''
        calculates weight from steel with stresss at 450C = 148MPa, sf=2, density = 7700 by default. cracking not modeled but could be in future
        '''
        stress = stress*10**6/sf
        thickness = pressure*10**5*self.r/stress
        shell_volume = thickness * math.pi * (self.length * 2 * self.r + self.r ** 2 * 4 / 3)
        cat_volume = self.r**2*self.length*math.pi
        shell_mass = shell_volume * shell_density
        cat_mass = cat_volume * cat_density
        return shell_mass,cat_mass


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
        self.update_fast()


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
                               self.myNH3 * (297 * self.T))

        #gamma [-]
        self.gamma = (1 + 1 / (self.yN2 / 0.4 + self.yH2 / 0.4 + self.yNH3 / 0.31))

        # volume [m^3
        self.volume = self.mass_tot / self.rho


    def update_slow(self):
        self.update()

        self.cp = float(pm.get('ig.N2').cp(self.T, self.p) * self.myN2 + pm.get('ig.H2').cp(self.T, self.p) * self.myH2 +
                         pm.get('ig.NH3').cp(self.T, self.p) * self.myNH3)
        # rho [kg/m^3]
        self.rho = float( 1 / (self.myN2 / pm.get('ig.N2').d(self.T, self.p) +
                               self.myH2 / pm.get('ig.H2').d(self.T, self.p) +
                               self.myNH3 / pm.get('ig.NH3').d(self.T, self.p)))
        # volume [m^3
        self.volume = self.mass_tot/self.rho

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

