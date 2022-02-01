# this is a basic setup for an isobaric reactor.

# input should be a timeseries of N2, H2 flows?
# no power use yet - this will come?
# assumptions
# - isobaric - constant pressure in reactor
import math
import pyromat as pm
import numpy as np
pm.config['unit_energy'] = 'J'  # default for pyromat is kJ


def sat_point_lookup(P_sat):
    # ammonia saturation point T(K), p(bar)
    data = [[2.3968597e+02, 2.4048963e+02, 2.4129599e+02, 2.4210506e+02, 2.4291683e+02, 2.4373133e+02, 2.4454856e+02,
             2.4536853e+02, 2.4619125e+02, 2.4701673e+02, 2.4784497e+02, 2.4867599e+02, 2.4950980e+02, 2.5034641e+02,
             2.5118582e+02, 2.5202804e+02, 2.5287309e+02, 2.5372097e+02, 2.5457170e+02, 2.5542527e+02, 2.5628171e+02,
             2.5714102e+02, 2.5800322e+02, 2.5886830e+02, 2.5973628e+02, 2.6060718e+02, 2.6148099e+02, 2.6235773e+02,
             2.6323742e+02, 2.6412005e+02, 2.6500564e+02, 2.6589421e+02, 2.6678575e+02, 2.6768028e+02, 2.6857781e+02,
             2.6947835e+02, 2.7038191e+02, 2.7128850e+02, 2.7219812e+02, 2.7311080e+02, 2.7402654e+02, 2.7494535e+02,
             2.7586724e+02, 2.7679222e+02, 2.7772030e+02, 2.7865150e+02, 2.7958581e+02, 2.8052326e+02, 2.8146386e+02,
             2.8240760e+02, 2.8335451e+02, 2.8430460e+02, 2.8525787e+02, 2.8621434e+02, 2.8717401e+02, 2.8813691e+02,
             2.8910303e+02, 2.9007239e+02, 2.9104500e+02, 2.9202087e+02, 2.9300001e+02, 2.9398244e+02, 2.9496816e+02,
             2.9595719e+02, 2.9694953e+02, 2.9794520e+02, 2.9894421e+02, 2.9994657e+02, 3.0095228e+02, 3.0196138e+02,
             3.0297385e+02, 3.0398972e+02, 3.0500899e+02, 3.0603169e+02, 3.0705781e+02, 3.0808737e+02, 3.0912039e+02,
             3.1015687e+02, 3.1119682e+02, 3.1224026e+02, 3.1328720e+02, 3.1433765e+02, 3.1539162e+02, 3.1644913e+02,
             3.1751018e+02, 3.1857479e+02, 3.1964297e+02, 3.2071473e+02, 3.2179008e+02, 3.2286904e+02, 3.2395162e+02,
             3.2503783e+02, 3.2612768e+02, 3.2722118e+02, 3.2831835e+02, 3.2941920e+02, 3.3052374e+02, 3.3163198e+02,
             3.3274394e+02, 3.3385963e+02, 3.3497906e+02, 3.3610224e+02, 3.3722919e+02, 3.3835992e+02, 3.3949444e+02,
             3.4063276e+02, 3.4177490e+02, 3.4292087e+02, 3.4407068e+02, 3.4522435e+02, 3.4638188e+02, 3.4754330e+02,
             3.4870861e+02, 3.4987782e+02, 3.5105096e+02, 3.5222803e+02, 3.5340905e+02, 3.5459403e+02, 3.5578298e+02,
             3.5697592e+02, 3.5817286e+02, 3.5937381e+02, 3.6057879e+02, 3.6178780e+02, 3.6300088e+02, 3.6421802e+02,
             3.6543924e+02, 3.6666455e+02, 3.6789397e+02, 3.6912752e+02, 3.7036520e+02, 3.7160703e+02, 3.7285303e+02,
             3.7410320e+02, 3.7535757e+02, 3.7661614e+02, 3.7787893e+02, 3.7914596e+02, 3.8041723e+02, 3.8169277e+02,
             3.8297258e+02, 3.8425669e+02, 3.8554510e+02, 3.8683783e+02, 3.8813489e+02, 3.8943630e+02, 3.9074208e+02,
             3.9205224e+02, 3.9336678e+02, 3.9468574e+02, 3.9600912e+02, 3.9733693e+02, 3.9866920e+02, 4.0000594e+02,
             4.0134715e+02, 4.0269287e+02, 4.0367440e+02],
            [1.0167000e-01, 1.0466000e-01, 1.0836000e-01, 1.1227000e-01, 1.1617000e-01, 1.2043000e-01, 1.2462000e-01,
             1.2919000e-01, 1.3359000e-01, 1.3858000e-01, 1.4340000e-01, 1.4857000e-01, 1.5372000e-01, 1.5927000e-01,
             1.6490000e-01, 1.7085000e-01, 1.7689000e-01, 1.8315000e-01, 1.8975000e-01, 1.9634000e-01, 2.0355000e-01,
             2.1048000e-01, 2.1835000e-01, 2.2578000e-01, 2.3423000e-01, 2.4220000e-01, 2.5109000e-01, 2.5981000e-01,
             2.6918000e-01, 2.7870000e-01, 2.8856000e-01, 2.9897000e-01, 3.0955000e-01, 3.2071000e-01, 3.3184000e-01,
             3.4402000e-01, 3.5574000e-01, 3.6904000e-01, 3.8160000e-01, 3.9562000e-01, 4.0909000e-01, 4.2438000e-01,
             4.3912000e-01, 4.5495000e-01, 4.7074000e-01, 4.8802000e-01, 5.0497000e-01, 5.2317000e-01, 5.4168000e-01,
             5.6085000e-01, 5.8107000e-01, 6.0124000e-01, 6.2332000e-01, 6.4496000e-01, 6.6864000e-01, 6.9141000e-01,
             7.1726000e-01, 7.4168000e-01, 7.6892000e-01, 7.9561000e-01, 8.2429000e-01, 8.5346000e-01, 8.8366000e-01,
             9.1552000e-01, 9.4791000e-01, 9.8208000e-01, 1.0161800e+00, 1.0534900e+00, 1.0900600e+00, 1.1228400e+00,
             1.1573400e+00, 1.1913600e+00, 1.2263900e+00, 1.2624500e+00, 1.3004000e+00, 1.3386400e+00, 1.3779900e+00,
             1.4185100e+00, 1.4602100e+00, 1.5031400e+00, 1.5473400e+00, 1.5938600e+00, 1.6407200e+00, 1.6889600e+00,
             1.7386100e+00, 1.7897300e+00, 1.8423500e+00, 1.8965100e+00, 1.9535300e+00, 2.0122600e+00, 2.0700900e+00,
             2.1254700e+00, 2.1879600e+00, 2.2493900e+00, 2.3095700e+00, 2.3744100e+00, 2.4379300e+00, 2.5063800e+00,
             2.5767500e+00, 2.6456900e+00, 2.7199700e+00, 2.7927300e+00, 2.8711400e+00, 2.9498500e+00, 3.0385300e+00,
             3.1098000e+00, 3.1971100e+00, 3.2805200e+00, 3.3661200e+00, 3.4539500e+00, 3.5440700e+00, 3.6412200e+00,
             3.7362300e+00, 3.8361900e+00, 3.9388100e+00, 4.0441900e+00, 4.1363700e+00, 4.2552400e+00, 4.3662700e+00,
             4.4744300e+00, 4.5852600e+00, 4.7049000e+00, 4.8338800e+00, 4.9504300e+00, 5.0763300e+00, 5.1987300e+00,
             5.3378100e+00, 5.4735600e+00, 5.6055300e+00, 5.7406900e+00, 5.8942700e+00, 6.0325000e+00, 6.1859200e+00,
             6.3350700e+00, 6.5003600e+00, 6.6571000e+00, 6.8132200e+00, 6.9909900e+00, 7.1503400e+00, 7.3274600e+00,
             7.5089700e+00, 7.6949800e+00, 7.8754400e+00, 8.0445800e+00, 8.2704400e+00, 8.4644100e+00, 8.6573400e+00,
             8.8889500e+00, 9.0857000e+00, 9.2987800e+00, 9.5414100e+00, 9.7588900e+00, 9.9877600e+00, 1.0241760e+01,
             1.0481950e+01, 1.0727780e+01, 1.0953950e+01]]
    return float(np.interp(P_sat, [10*x for x in data[1]], data[0]))


def heat_exchanger_hotgas2coldgas(INPUT1, INPUT2, e1=0.8):  # mol,mol,mol,K,bar,K,m/s,mm check units!!
    '''

    :param INPUT1:
    :param INPUT2:
    :param e1:
    :return:
    '''
    [H2in1, N2in1, NH3in1, Tin1, Pin1] = INPUT1  # hot stream from reactor
    [H2in2, N2in2, NH3in2, Tin2, Pin2] = INPUT2  # cold stream into reactor

    mN2in1 = N2in1 * 28.0134 / 1000
    mH2in1 = H2in1 * 2.016 / 1000
    mNH3in1 = NH3in1 * 17.03 / 1000
    mtot1 = mN2in1 + mH2in1 + mNH3in1

    mN2in2 = N2in2 * 28.0134 / 1000
    mH2in2 = H2in2 * 2.016 / 1000
    mNH3in2 = NH3in2 * 17.03 / 1000
    mtot2 = mN2in2 + mH2in2 + mNH3in2

    N2data = pm.get('ig.N2')
    H2data = pm.get('ig.H2')
    NH3data = pm.get('ig.NH3')

    C_1 = float(mN2in1 * N2data.cp(Tin1, Pin1) + mH2in1 * H2data.cp(Tin1, Pin1) + mNH3in1 * NH3data.cp(Tin1, Pin1))  # J/kg/K
    C_2 = float(mN2in2 * N2data.cp(Tin2, Pin2) + mH2in2 * H2data.cp(Tin2, Pin2) + mNH3in2 * NH3data.cp(Tin2, Pin2))  # J/kg/K

    Cmin = min(C_1, C_2)  # should be C1
    Cmax = max(C_1, C_2)  # should be C2
    Cr = Cmin/Cmax

    Q = e1*(Cmin*(Tin1-Tin2))

    Tout1 = Tin1 - Q / C_1
    Tout2 = Tin2 + Q / C_2

    # NTU = np.log(1-e1*(1+Cr))/(1+Cr)
    #
    # U1 = 300 # W/m^2/K
    #
    # A1 = NTU * Cmin / U1
    #
    # rho_mix = yN2 * N2data.d(T_hot_in, Pin) + yH2 * H2data.d(T_hot_in, Pin) + yNH3 * NH3data.d(T_hot_in, Pin)
    #
    # Q_flow = m_mix / rho_mix
    # num_p = 4 * Q_flow / (Vmax * np.pi * D ** 2)
    # l1 = A1 / (num_p * np.pi * D)
    #
    # mu_N2 = 1.663 * 10 ** -5 * (T_hot_in / 273) ** (3 / 2) * (273 + 107) / (T_hot_in + 107)
    # mu_H2 = 8.411 * 10 ** -5 * (T_hot_in / 273) ** (3 / 2) * (273 + 97) / (T_hot_in + 97)
    # mu_NH3 = 0.919 * 10 ** -5 * (T_hot_in / 273) ** (3 / 2) * (273 + 370) / (T_hot_in + 370)
    #
    # mu_mix = yN2 * mu_N2 + yH2 * mu_H2 + yNH3 * mu_NH3
    #
    # Re = rho_mix*Vmax*D/mu_mix
    # F_fact = 0.316*Re**-0.25
    # Del_P = F_fact*rho_mix*Vmax^2*l1/2/D

    OUTPUT1 = [H2in1,N2in1,NH3in1,Tout1,Pin1] # out to condensor
    OUTPUT2 = [H2in2,N2in2,NH3in2,Tout2,Pin2] # out to reactor

    return OUTPUT1,OUTPUT2, Tout2-Tin2


def heat_exchanger_water2gas(INPUT,water_mass_flow,T_cold_in=10+273,Vmax=5,D=0.006): ### mol,mol,mol,K,bar,K,m/s,mm check units!!

    [H2in, N2in, NH3in, T_hot_in, Pin] = INPUT
    tot_in_mol = N2in + H2in + NH3in
    yN2 = N2in / tot_in_mol
    yH2 = H2in / tot_in_mol
    yNH3 = NH3in / tot_in_mol

    mN2in = N2in * 28.0134/1000
    mH2in = H2in * 2.016/1000
    mNH3in = NH3in * 17.03/1000
    m_mix = mN2in + mH2in + mNH3in

    N2data = pm.get('ig.N2')
    H2data = pm.get('ig.H2')
    NH3data = pm.get('ig.NH3')

    P_NH3 = yNH3*Pin

    T_sat = sat_point_lookup(P_NH3)

    C_mix = float(mN2in * N2data.cp(T_hot_in, Pin) + mH2in * H2data.cp(T_hot_in, Pin) + mNH3in * NH3data.cp(T_hot_in, Pin)) # J/kg/K

    Q = C_mix * (T_hot_in - T_sat)

    C_cool = water_mass_flow * 4180# 1 kg/s * 4.18 J/kg/K

    Cmin = min(C_cool, C_mix)
    Cmax = max(C_cool, C_mix)
    Cr = Cmin/Cmax


    e1 = Q/(Cmin*(T_hot_in-T_cold_in))

    NTU = - np.log(1-e1*(1+Cr))/(1+Cr) #


    U1 = 300 # W/m^2/K

    A1 = NTU * Cmin / U1

    rho_mix = float(yN2 * N2data.d(T_hot_in, Pin) + yH2 * H2data.d(T_hot_in, Pin) + yNH3 * NH3data.d(T_hot_in, Pin)) # kg/m^3

    Q_flow = m_mix / rho_mix # m^3/s?
    num_p = 4 * Q_flow / (Vmax * np.pi * D ** 2)
    l1 = A1 / (num_p * np.pi * D)

    mu_N2 = 1.663 * 10 ** -5 * (T_hot_in / 273) ** (3 / 2) * (273 + 107) / (T_hot_in + 107)
    mu_H2 = 8.411 * 10 ** -5 * (T_hot_in / 273) ** (3 / 2) * (273 + 97) / (T_hot_in + 97)
    mu_NH3 = 0.919 * 10 ** -5 * (T_hot_in / 273) ** (3 / 2) * (273 + 370) / (T_hot_in + 370)

    mu_mix = yN2 * mu_N2 + yH2 * mu_H2 + yNH3 * mu_NH3

    Re = rho_mix*Vmax*D/mu_mix
    F_fact = 0.316*Re**-0.25
    Del_P = F_fact*rho_mix*Vmax**2*l1/2/D

    OUTPUT = [H2in,N2in,NH3in,T_sat,Pin-Del_P*10**-5]
    return OUTPUT, Q

def condensor(INPUT,e2=0.8,T_cold_in=10+273,Vmax=5,D=6*0.001): ### check units!!

    [H2in, N2in, NH3in, T_hot_in, Pin] = INPUT
    tot_in_mol = N2in + H2in + NH3in
    yN2 = N2in / tot_in_mol
    yH2 = H2in / tot_in_mol
    yNH3 = NH3in / tot_in_mol

    mN2in = N2in * 28.0134/1000
    mH2in = H2in * 2.016/1000
    mNH3in = NH3in * 17.03/1000
    m_mix = mN2in + mH2in + mNH3in

    N2data = pm.get('ig.N2')
    H2data = pm.get('ig.H2')
    NH3data = pm.get('ig.NH3')

    Del_H_c = 22.7*NH3in # kJ/mol * mol

    C_mix = float(mN2in * N2data.cp(T_hot_in, Pin) + mH2in * H2data.cp(T_hot_in, Pin) + mNH3in * NH3data.cp(T_hot_in, Pin)) # J/kg/K

    C_cool = 1 * 4180# 1 kg/s * 4180 J/kg/K

    Cmin = min(C_cool, C_mix)
    Cmax = max(C_cool, C_mix)
    Cr = Cmin/Cmax



    #e2 = Del_H_c/(Cmin*(T_hot_in-T_cold_in))

    Del_H_act = e2*Del_H_c #J

    NTU = -np.log(1-e2)


    U2 = 750 # W/m^2/K

    A1 = NTU * Cmin / U2
    rho_mix = float(yN2 * N2data.d(T_hot_in, Pin) + yH2 * H2data.d(T_hot_in, Pin) + yNH3 * NH3data.d(T_hot_in, Pin))
    Q_flow = m_mix/rho_mix
    num_p = 4 * Q_flow / (Vmax * np.pi * D ** 2)
    l1 = A1 / (num_p * np.pi * D)

    mu_N2 = 1.663 * 10 ** -5 * (T_hot_in / 273) ** (3 / 2) * (273 + 107) / (T_hot_in + 107)
    mu_H2 = 8.411 * 10 ** -5 * (T_hot_in / 273) ** (3 / 2) * (273 + 97) / (T_hot_in + 97)
    mu_NH3 = 0.919 * 10 ** -5 * (T_hot_in / 273) ** (3 / 2) * (273 + 370) / (T_hot_in + 370)

    mu_mix = yN2 * mu_N2 + yH2 * mu_H2 + yNH3 * mu_NH3

    Re = rho_mix*Vmax*D/mu_mix
    F_fact = 0.316*Re**-0.25
    Del_P = F_fact*rho_mix*Vmax**2*l1/2/D

    OUTPUT = [H2in,N2in,NH3in*(1-e2),T_hot_in,Pin-Del_P*10**-5]
    return OUTPUT, Del_H_act

def reactorStep(INPUT): #mol/s, K, Pa
    '''
    An iterative function to determine the change in state variables and reactants over the length of a reactor Bed step

    :param INPUT: molar flow rates of H2,N2,NH3, and Temp and Pressure of inlet to Bed step

    :return: molar flow rates of H2,N2,NH3, and Temp and Pressure of outlet from Bed step
    '''
    [H2in, N2in, NH3in, Tin, Pin, DelX, Areac] = INPUT

    # initial concentrations
    tot_in_mol = N2in + H2in + NH3in
    yN2 = N2in / tot_in_mol
    yH2 = H2in / tot_in_mol
    yNH3 = NH3in / tot_in_mol

    mN2in = N2in * 28.0134/1000
    mH2in = H2in * 2.016/1000
    mNH3in = NH3in * 17.03/1000
    mdot = mN2in + mH2in + mNH3in



    # activity coefficients for all species
    N2fuga = 0.93431737 + 0.3101804*10**-3 * Tin + 0.295895*10**-3 * Pin - 0.270729*10**-6 * Tin**2 + \
             0.4775207*10**-6 * Pin**2
    H2fuga = math.exp(math.exp(-3.8402*Tin**0.125+0.541)*Pin - math.exp(-0.1263*Tin**0.5-15.980)*Pin**2 +
                      300*math.exp(-0.011901*Tin-5.941)*(math.exp(-Pin/300)-1))
    NH3fuga = 0.1438996 + 0.2028538*10**-2*Tin - 0.4487672*10**-3*Pin - 0.1142945*10**-5*Tin**2 + \
              0.2761216*10**-6*Pin**2
    a_N2 = yN2 * N2fuga * Pin #bar
    a_H2 = yH2 * H2fuga * Pin #bar
    a_NH3 = yNH3 * NH3fuga * Pin #bar

    # reaction rate constant
    Ea = 1.7056*10**5 #J/mol
    R = 8.31446261815324 #J/K/mol
    k0 = 8.8490*10**17 #mol/m^3
    k_r = k0*math.exp(-Ea/(R*Tin)) #mol/m^3/h

    # equilibrium constant
    K_eq = math.pow(10,-2.691122*math.log(Tin,10) - 5.519265*10**-5*Tin + 1.848863*10**-7*Tin**2 + 2001.6/Tin + 2.67899)

    # Reaction Rate
    alpha = 0.5
    RR_NH3 = k_r*(K_eq**2*a_N2*(a_H2**3/(a_NH3**2))**alpha - (a_NH3**2/(a_H2**3))**(1-alpha))/3600
    RR_N2 = -1/2*RR_NH3
    RR_H2 = -3/2*RR_NH3

    # Cp
    N2 = pm.get('ig.N2')
    H2 = pm.get('ig.H2')
    NH3 = pm.get('ig.NH3')
    m_Cp_tot = float(mN2in*N2.cp(Tin,Pin) + mH2in*H2.cp(Tin,Pin) + mNH3in*NH3.cp(Tin,Pin)) # J/kg/K

    # Heat of Reaction
    Del_H = 4.184*(-(0.54526 + 846.609/Tin + 459.734*10**6*Tin**-3)*Pin - 5.34685*Tin -0.2525*10**-3*Tin**2 + 1.69197*10**-6*Tin**3 -9157.09) #J/mol


    # production
    N2out = N2in + DelX * RR_N2 * Areac  #### check these - conversion percentages are not moles ------- need to fix
    H2out = H2in + DelX * RR_H2 * Areac
    NH3out = NH3in + DelX * RR_NH3 * Areac

    Tout = Tin - DelX * Areac * Del_H * RR_NH3 / (m_Cp_tot)

    #print((2*N2out+NH3out),(2*N2in+NH3in))
    #print((2*H2out + 3*NH3out), (2*H2in+3*NH3in))
    OUTPUT = [H2out, N2out, NH3out, Tout, Pin]
    return OUTPUT

def mixer(INPUT1,INPUT2):
    [H2in1, N2in1, NH3in1, Tin1, Pin1] = INPUT1
    [H2in2, N2in2, NH3in2, Tin2, Pin2] = INPUT2

    mN2in1 = N2in1 * 28.0134 / 1000
    mH2in1 = H2in1 * 2.016 / 1000
    mNH3in1 = NH3in1 * 17.03 / 1000
    mtot1 = mN2in1 + mH2in1+mNH3in1

    mN2in2 = N2in2 * 28.0134 / 1000
    mH2in2 = H2in2 * 2.016 / 1000
    mNH3in2 = NH3in2 * 17.03 / 1000
    mtot2 = mN2in2 + mH2in2 + mNH3in2


    N2 = pm.get('ig.N2')
    H2 = pm.get('ig.H2')
    NH3 = pm.get('ig.NH3')

    Pout = min(Pin1,Pin2)
    Tav = (Tin1*mtot1 + Tin2*mtot2)/((mtot1 + mtot2))

    Qin1 = (mN2in1*N2.cp(Tin1, Pin1) + mH2in1*H2.cp(Tin1, Pin1) + mNH3in1*NH3.cp(Tin1, Pin1)) * Tin1

    Qin2 = (mN2in2*N2.cp(Tin2, Pin2) + mH2in2*H2.cp(Tin2, Pin2) + mNH3in2*NH3.cp(Tin2, Pin2)) * Tin2

    QoutperT = ((mN2in1 + mN2in2)*N2.cp(Tav, Pout) + mH2in1*H2.cp(Tav, Pout) + mNH3in1*NH3.cp(Tav, Pout) + mN2in2*N2.cp(Tav, Pout) + mH2in2*H2.cp(Tav, Pout) + mNH3in2*NH3.cp(Tav, Pout))

    Tout = float((Qin1 + Qin2)/QoutperT)

    OUTPUT = [H2in1 + H2in2,N2in1 + N2in2,NH3in1 + NH3in2,Tout,Pout]
    return OUTPUT

def compressor(INPUT, Pout, eta=0.7):
    """
    Function to return rate of work for a compressor based on a
    target pressure.

    Inputs: INPUT - 5x1 list of [float], standard input of mol H2, mol N2, mol NH3, T[K] and p[bar]
            eta - float, efficiency of compressor

    Outputs: OUTPUT - 5x1 list of [float], standard output of mol H2, mol N2, mol NH3, T[K] and p[bar]
             w - float, output
    """
    [H2in, N2in, NH3in, Tin, Pin] = INPUT
    tot_in_mol = N2in + H2in + NH3in
    yN2 = N2in / tot_in_mol
    yH2 = H2in / tot_in_mol
    yNH3 = NH3in / tot_in_mol

    mN2in = N2in * 28.0134/1000
    mH2in = H2in * 2.016/1000
    mNH3in = NH3in * 17.03/1000
    m_mix = mN2in + mH2in + mNH3in

    N2data = pm.get('ig.N2')
    H2data = pm.get('ig.H2')
    NH3data = pm.get('ig.NH3')


    C_mix = float(mN2in * N2data.cp(Tin, Pin) + mH2in * H2data.cp(Tin, Pin) + mNH3in * NH3data.cp(Tin, Pin)) # J/kg/K

    y = float(1 + 1/((yN2/0.4) + (yH2/0.4) + yNH3/(NH3data.gam(Tin,Pin)-1)))
    r_p = Pout/Pin
    a = (y-1)/y
    w = C_mix*Tin/eta*(r_p**a-1)
    Tout = Tin*(1+r_p**a/eta-1/eta)
    OUTPUT = [H2in, N2in, NH3in, Tout, Pout]
    return OUTPUT, w

def compressor_power(INPUT, Pout, eta=0.7):
    """
    Function to return rate of work for a compressor based on a
    target pressure.

    Inputs: INPUT - 5x1 list of [float], standard input of mol H2, mol N2, mol NH3, T[K] and p[bar]
            Pout - float, target pressure
            eta - float, efficiency of compressor

    Outputs: w - float, work rate required (W)
    """
    [H2in, N2in, NH3in, Tin, Pin] = INPUT
    tot_in_mol = N2in + H2in + NH3in
    yN2 = N2in / tot_in_mol
    yH2 = H2in / tot_in_mol
    yNH3 = NH3in / tot_in_mol

    mN2in = N2in * 28.0134/1000
    mH2in = H2in * 2.016/1000
    mNH3in = NH3in * 17.03/1000
    m_mix = mN2in + mH2in + mNH3in

    N2data = pm.get('ig.N2')
    H2data = pm.get('ig.H2')
    NH3data = pm.get('ig.NH3')


    C_mix = float(mN2in * N2data.cp(Tin, Pin) + mH2in * H2data.cp(Tin, Pin) + mNH3in * NH3data.cp(Tin, Pin)) # J/kg/K

    y = float(1 + 1/((yN2/0.4) + (yH2/0.4) + yNH3/(NH3data.gam(Tin,Pin)-1)))
    r_p = Pout/Pin
    a = (y-1)/y
    W = C_mix*Tin/eta*(r_p**a-1)
    Tout = Tin*(1+r_p**a/eta-1/eta)
    return W

def PSA_estimate(INPUT, Pout=10, eta=0.7):
    """
    Function to return power for PSA to separate given mol of N2.

    Inputs:     INPUT - 3x1 list of [float],mol N2, T[K] and p[bar]
                Pout - float, target pressure
                eta - float, efficiency of compressor

    Outputs:    OUTPUT - 5x1 list of [float], standard output of mol H2, mol N2, mol NH3, T[K] and p[bar]
                w - float, work rate required (W) for PSA
    """
    [N2in, Tin, Pin] = INPUT
    mN2in = N2in * 28.0134/1000

    N2data = pm.get('ig.N2')

    C_mix = float(mN2in * N2data.cp(Tin, Pin))/0.79 # J/kg/K

    y = 1.4
    r_p = Pout/Pin
    a = (y-1)/y
    W = C_mix*Tin/eta*(r_p**a-1)
    Tout = Tin*(1+r_p**a/eta-1/eta)
    output_state = State(0, N2in, 0, Tout, Pout)
    return output_state, W


def electrolysis(H2mol,eta = 0.7): #mol/s to kW
    """
        Function to return power for PEM electrolysis to separate given mol of H2.

        Inputs:     H2 mol input + output

        Outputs:    OUTPUT - 5x1 list of [float], standard output of mol H2, mol N2, mol NH3, T[K] and p[bar]
                    w - float, work rate required (W) for PSA
                    H20_in - required water to produce hydrogen
        """
    W = 241.83 / eta * H2mol
    OUTPUT = [H2mol, 0, 0, 298, 1]
    H20_in = H2mol*18/2
    return OUTPUT, W,  H20_in

def heater(State,T_end):
    power = State.cp*State.mass_tot*(T_end - State.T)
    State.T = T_end
    return State,power

class Bed(object):
    '''
    A class to store properties of a 1D Bed model for a reactor.
    also generates a vector to describe end
    '''

    def __init__(self, length, diam, mini, newvectmethod = True):
        self.length = length
        self.diam = diam
        self.area = math.pi * (diam / 2) ** 2 # m^2
        # generate vect
        if not newvectmethod:
            self.vect = [0] + np.geomspace(mini, length / 10, 19, endpoint=False).tolist() + np.linspace(length / 10, length, 10).tolist()
        else:
            maxi = 0
            scale_max = np.floor(np.log10(mini))+1
            scaler = 0
            OUTPUT = [0]
            while maxi < self.length:
                maxi = maxi + mini*10**scaler
                OUTPUT.append(maxi)
                if np.log10(maxi) + 0.001 >= scale_max:
                    scaler += 1
                    scale_max += 1
            self.vect = OUTPUT

        self.vectlen = len(self.vect)

class State(object):
    '''
    Class to store variables about mass and molar flow rate, temperature, pressure etc of data.
    Species: molar flow rate [mol/s]
    ySpecies: molar concentration [mol/mol]
    mSpecies: mass flow rate [kg/s
    :param H2: molar flow rate [mol/s] of H2 gas.
    :param N2: molar flow rate [mol/s] of N2 gas.
    :param NH3: molar flow rate [mol/s] of NH3 gas.
    :param T: Temperature of flow [K].
    :param P: Pressure of flow [bar].
    '''
    def __init__(self,H2,N2,NH3,T,P):
        self.H2 = H2
        self.N2 = N2
        self.NH3 = NH3
        self.T = T
        self.P = P
        self.update()


    def update(self):
        self.mH2 = self.H2 * 2.016 / 1000
        self.mN2 = self.N2 * 28.0134 / 1000
        self.mNH3 = self.NH3 * 17.03 / 1000
        self.mol_tot = self.H2 + self.N2 + self.NH3
        self.mass_tot = self.H2 * 2.016 + self.N2 * 28.0134 + self.NH3 * 17.03
        self.yH2 = self.H2 / self.mol_tot
        self.yN2 = self.N2 / self.mol_tot
        self.yNH3 = self.NH3 / self.mol_tot

        self.cp = float((pm.get('ig.N2').cp(self.T, self.P) * self.mN2 + pm.get('ig.H2').cp(self.T, self.P) * self.mH2 +
                   pm.get('ig.NH3').cp(self.T, self.P) * self.mNH3) / self.mass_tot)

        self.rho = float(pm.get('ig.N2').d(self.T, self.P) * self.yN2 + pm.get('ig.H2').d(self.T, self.P) * self.yH2 + \
                   pm.get('ig.NH3').d(self.T, self.P) * self.yNH3)

        mu_N2 = 1.663 * 10 ** -5 * (self.T / 273) ** (3 / 2) * (273 + 107) / (self.T + 107)
        mu_H2 = 8.411 * 10 ** -5 * (self.T / 273) ** (3 / 2) * (273 + 97) / (self.T + 97)
        mu_NH3 = 0.919 * 10 ** -5 * (self.T / 273) ** (3 / 2) * (273 + 370) / (self.T + 370)
        self.mu = self.yN2 * mu_N2 + self.yH2 * mu_H2 + self.yNH3 * mu_NH3