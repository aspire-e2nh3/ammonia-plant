import pyromat as pm
import numpy as np


pm.config['unit_energy'] = 'J'  # default for pyromat is kJ


def rp2w(gname, m, st_0, r_p, eta=0.7):
    """Function to return rate of work for a compressor based on a
    target pressure ratio.

    Inputs: gname - str, name of substance compatible with PYroMat
            m - mass flow rate (kg/s)
            st_0 - tuple of numpy array objects [T(K), p(bar)],
                   input state compatible with PYroMat
            eta - float, efficiency of compressor
            r_p - float, pressure ratio

    Outputs: st_1 - tuple of numpy array objects [T, p],
                    output state compatible with PYroMat
             w - numpy array object, work rate required (W)
    """
    T_0, p_0 = st_0      # unpack input state tuple
    gas = pm.get(gname)  # retrieve substance object
    c_p = gas.cp(T_0, p_0)   # retrieve thermodynamic properties
    y = gas.gam(T_0, p_0)
    a = (y-1)/y
    w = m*c_p*T_0/eta*(r_p**a-1)
    T_1 = T_0*(1+r_p**a/eta-1/eta)
    p_1 = p_0*r_p
    st_1 = [T_1, p_1]
    return w, st_1


def p12w(gname, m, st_0, p_1, eta=0.7):
    """Function to return rate of work for a compressor based on a
    target pressure.

    Inputs: gname - str, name of substance compatible with PYroMat
            m - mass flow rate (kg/s)
            st_0 - tuple of numpy array objects [T(K), p(bar)],
                   input state compatible with PYroMat
            p_1 - numpy array, target pressure (bar)
            eta - float, efficiency of compressor

    Outputs: st_1 - tuple of numpy array objects [T, p],
                    output state compatible with PYroMat
             w - numpy array object, work rate required (W)
    """
    T_0, p_0 = st_0      # unpack input state tuple
    gas = pm.get(gname)  # retrieve substance object
    c_p = gas.cp(T_0, p_0)   # retrieve thermodynamic properties
    y = gas.gam(T_0, p_0)
    r_p = p_1/p_0
    a = (y-1)/y
    w = m*c_p*T_0/eta*(r_p**a-1)
    T_1 = T_0*(1+r_p**a/eta-1/eta)
    st_1 = [T_1, p_1]
    return w, st_1


def w2rp(gname, m, st_0, w, eta=0.7):
    """Function to obtain the pressure ratio delivered for a compressor based
    on supplied power.

    Inputs: gname - str, name of substance compatible with PYroMat
            m - mass flow rate (kg/s)
            st_0 - tuple of numpy array objects [T(K), p(bar)],
                   input state compatible with PYroMat
            w - numpy array, power in (W)
            eta - float, efficiency of compressor (default is 0.7)

    Outputs: st_1 - tuple of numpy array objects [T, p],
                    output state compatible with PYroMat
             w - numpy array object, work rate required (W)
    """

    T_0, p_0 = st_0      # unpack input state tuple
    gas = pm.get(gname)  # retrieve substance object
    c_p = gas.cp(T_0, p_0)   # retrieve thermodynamic properties
    y = gas.gam(T_0, p_0)
    a = (y-1)/y

    r_p = ((w*eta/m/c_p/T_0)+1)**(y/(y-1))
    T_1 = T_0*(1+r_p**a/eta-1/eta)

    p_1 = p_0*r_p
    st_1 = [T_1, p_1]

    return r_p, st_1


def w2p1(gname, m, st_0, w, eta=0.7):
    """Function to obtain the pressure delivered for a compressor based
    on supplied power.

    Inputs: gname - str, name of substance compatible with PYroMat
            m - mass flow rate (kg/s)
            st_0 - tuple of numpy array objects [T(K), p(bar)],
                   input state compatible with PYroMat
            w - numpy array, power in (W)
            eta - float, efficiency of compressor (default is 0.7)

    Outputs: st_1 - tuple of numpy array objects [T, p],
                    output state compatible with PYroMat
             w - numpy array object, work rate required (W)
    """

    T_0, p_0 = st_0      # unpack input state tuple
    gas = pm.get(gname)  # retrieve substance object
    c_p = gas.cp(T_0, p_0)   # retrieve thermodynamic properties
    y = gas.gam(T_0, p_0)
    a = (y-1)/y

    r_p = ((w*eta/m/c_p/T_0)+1)**(1/a)
    T_1 = T_0*(1+r_p**a/eta-1/eta)

    p_1 = p_0*r_p
    st_1 = [T_1, p_1]

    return st_1


def main():
    """ Runs some example problems for H2, N2 and air compression
    given a certain available power.
    """
    gname = 'ig.N2'  # chosen gas
    m = 1  # mass flow rate (kg/s)
    eta = 0.7  # isentropic efficiency of compressor
    T_0 = np.array(300)  # input temperature (K)
    p_0 = np.array(1)  # input pressure (bar)
    st_0 = [T_0, p_0]  # converting to state tuple

    # Example 1
    w = 1e5  # available power (W)
    st_1 = w2p1(gname, m, st_0, w, eta)
    T_1, p_1 = st_1

    print('Example 1:\n'
          'A theoretical compressor with:\n'
          '\tEfficiency = %.2f\n'
          '\tPower supplied = %.2f W\n'
          'Inputs:\n'
          '\tSubstance = %s\n'
          '\tMass flow rate = %.2f kg/s\n'
          '\tTemperature = %.2f K\n'
          '\tPressure = %.2f bar\n'
          'Using "w2p1..."\n'
          'Outputs:\n'
          '\tTemperature = %.2f K\n'
          '\tPressure = %.2f bar\n' %
          (eta, w, gname, m, T_0, p_0, T_1, p_1))

    # Example 2
    p_1 = 200  # output pressure
    w, st_1 = p12w(gname, m, st_0, p_1, eta)
    T_1 = st_1[0]

    print('Example 2:\n'
          'A theoretical compressor with:\n'
          '\tEfficiency = %.2f\n'
          'Inputs:\n'
          '\tSubstance = %s\n'
          '\tMass flow rate = %.2f kg/s\n'
          '\tInput Temperature = %.2f K\n'
          '\tInput Pressure = %.2f bar\n'
          '\tOutput Pressure = %.2f bar\n'
          'Using "p12w..."\n'
          'Outputs:\n'
          '\tOutput Temperature = %.2f K\n'
          '\tPower needed = %.2f W\n' %
          (eta, gname, m, T_0, p_0, p_1, T_1, w))


if __name__ == "__main__":
    main()
