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


def sat_point_temp_lookup(P_sat):
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
             1.0481950e+01, 1.0727780e+01, 1.0953950e+01]]  # data in K, MPa
    return float(np.interp(P_sat, [10 * x for x in data[1]], data[0]))

def sat_point_pressure_lookup(T_sat):
    # ammonia saturation point T(K), p(bar)
    data = [
        [2.3968597e+02, 2.4048963e+02, 2.4129599e+02, 2.4210506e+02, 2.4291683e+02, 2.4373133e+02, 2.4454856e+02,
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
         1.0481950e+01, 1.0727780e+01, 1.0953950e+01]]  # data in K, MPa
    return float(np.interp(T_sat, data[0],[10 * x for x in data[1]]))


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
    s1.update_special()
    s2.update_special()

    s1_out = copy.copy(s1)
    s2_out = copy.copy(s2)

    Cmin = min(s1_out.cp * s1_out.mass_tot, s2_out.cp * s2_out.mass_tot)
    # Cmax = max(s1.cp*s1.mass_tot, s2.cp*s2.mass_tot)
    # Cr = Cmin / Cmax

    Q = effectiveness * (Cmin * (s1_out.T - s2_out.T))

    s1_out.T += - Q / (s1_out.cp * s1_out.mass_tot)
    s2_out.T += Q / (s2_out.cp * s2_out.mass_tot)


    s1_out.update()
    s2_out.update()
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
    s1.update_special()
    s2.update_special()

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

    s1_out.update()
    s2_out.update()
    return s1_out, s2_out, -(s2.T - s2_out.T), effectiveness


def heat_exchanger_water2gas(s, T_end=0, cool_to_sat_point=False, effectiveness=0, water_mfr=10,  T_cin=273, Vmax=5, D=0.006):  # -,kg,K,m/s,m
    """
    saturation
    """
    s.update_special()
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
    s_out.update()
    return s_out, Q, effectiveness


def condenser(s, effectiveness=0.8, water_mfr=1, T_cin=10+273, Vmax=5, D=0.006):  ### check units!!

    s.update_special()
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
    s_out.update()
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
    s_out.update()
    return s_out, power, ammonia_removed, T_cout


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

    # Reaction Rate
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
    s.T += -dX * area * Del_H * RR_NH3 / (s.cp * s.mass_tot)

    # update state
    s.update()
    return s


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

    s_out.update()
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
    s.update_special()
    s_out = copy.deepcopy(s)

    y = s_out.gamma
    r_p = p_out / s.p
    a = (y - 1) / y
    power = s.cp * s.mass_tot*s.T / eta * (r_p ** a - 1)
    s_out.T = s.T * (1 + r_p ** a / eta - 1 / eta)
    s_out.p = p_out
    s_out.update()
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


    s.update_special()
    s_out = copy.deepcopy(s)

    y = s_out.gamma
    R = (y-1)*s_out.cp
    power = s_out.mass_tot*R*a*(t_out - s_out.T)
    Q_out = s_out.mass_tot*R*a*(t_out - s_out.T) - s_out.mass_tot*s_out.cp*(t_out - s_out.T)
    s_out.T = t_out
    s_out.p = p_out
    s_out.update()
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
    s_out.update()
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
        self.area = math.pi * self.r ** 2  # m^2
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
        self.update_special()

    def update(self):
        # masses of components in [kg]
        self.mH2 = self.H2 * 2.016 / 1000
        self.mN2 = self.N2 * 28.0134 / 1000
        self.mNH3 = self.NH3 * 17.03 / 1000
        self.mass_tot = self.mH2 + self.mN2 + self.mNH3

        # total mol
        self.mol_tot = self.H2 + self.N2 + self.NH3
        # mol
        self.yH2 = self.H2 / self.mol_tot
        self.yN2 = self.N2 / self.mol_tot
        self.yNH3 = self.NH3 / self.mol_tot

        # cp [J/kg/K]
        self.cp = float((pm.get('ig.N2').cp(self.T, self.p) * self.mN2 + pm.get('ig.H2').cp(self.T, self.p) * self.mH2 +
                         pm.get('ig.NH3').cp(self.T, self.p) * self.mNH3) / self.mass_tot)

    def update_special(self):
        self.update()
        # rho [kg/m^3]
        self.rho = float(pm.get('ig.N2').d(self.T, self.p) * self.yN2 + pm.get('ig.H2').d(self.T, self.p) * self.yH2 + \
                         pm.get('ig.NH3').d(self.T, self.p) * self.yNH3)

        # gamma from addition of specific heats(1+ 1/sum of mol%/gamma)[-]
        self.gamma = float((1 + 1 / ((self.yN2 / 0.4) + (self.yH2 / 0.4) + self.yNH3 / (pm.get('ig.NH3').gam(self.T, self.p) - 1))))

        # mu from sutherland's formula [Ns/m^2]
        mu_N2 = 1.663 * 10 ** -5 * (self.T / 273) ** (3 / 2) * (273 + 107) / (self.T + 107)
        mu_H2 = 8.411 * 10 ** -5 * (self.T / 273) ** (3 / 2) * (273 + 97) / (self.T + 97)
        mu_NH3 = 0.919 * 10 ** -5 * (self.T / 273) ** (3 / 2) * (273 + 370) / (self.T + 370)
        self.mu = self.yN2 * mu_N2 + self.yH2 * mu_H2 + self.yNH3 * mu_NH3

    def store(self):
        return [self.H2, self.N2, self.NH3, self.T, self.p]

    def split(self,f):
        s1 = State(self.H2 * f, self.N2 * f, self.NH3 * f, self.T, self.p)
        s2 = State(self.H2 * (1 - f), self.N2 * (1 - f), self.NH3 * (1 - f), self.T, self.p)
        return s1, s2

    def __copy__(self):
        return State(self.H2,self.N2,self.NH3,self.T,self.p)

