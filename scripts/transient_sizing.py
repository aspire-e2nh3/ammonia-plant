import os
import argparse
import pandas as pd
from scipy.integrate import trapezoid
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from sklearn.metrics import balanced_accuracy_score

parser = argparse.ArgumentParser()
parser.add_argument(
    "-sc",
    "--scale",
    help="specify the capacity of the transient energy source in W")
args = parser.parse_args()


def read_wind_data(pfolder, fname):
    """Read experimental data from a tab delimited file.

    :param pfolder: string of folder path
    :param fname: string of filename of sasmodel data

    :returns: 
    """
    
    df = pd.read_csv(os.path.join(pfolder,fname),
                             header = 0)
    df.columns = ['timestamp', 'wind']
    
    return df


def raw2trans(raw, scale=1e5):
    """Converts raw data from timestamp and power in MW to
    ellapsed time in seconds and normalised power."""

    norm = scale*((raw.wind-raw.wind.min())/(raw.wind.max()-raw.wind.min()))
    secs = []
    t_0 = datetime.strptime(raw.timestamp[0], ' %Y-%m-%d %H:%M:%S')
    for t in raw.timestamp:
        dt = datetime.strptime(t, ' %Y-%m-%d %H:%M:%S')-t_0
        secs.append(dt.total_seconds())

    df_trns = pd.DataFrame(None,columns=['time','power'])
    df_trns.time = secs
    df_trns.power = norm
    return df_trns

def find_baseload(trns):
    """Find the baseload power from the transient signal."""

    E = trapezoid(trns.power, trns.time)
    E_dot = E/trns.time.iat[-1]

    df_bsld = pd.DataFrame(None, columns=['time','power'])
    df_bsld.time = trns.time
    df_bsld.power = E_dot

    return df_bsld

def find_balance(trns, bsld):
    """Find the balance signal from the transient signal 
    and corresponding baseload."""

    df_bal = pd.DataFrame(None, columns=['time','power'])
    df_bal.time = trns.time
    df_bal.power = trns.power - bsld.power

    return df_bal

def find_charges(bal):
    """Identify all charge and discharge instances."""
    bal_p = pd.DataFrame(None, columns=['time','power'])
    bal_p.time = bal.time
    bal_n = pd.DataFrame(None, columns=['time','power'])
    bal_n.time = bal.time

    bal_p.power = bal.power
    bal_p.power[bal_p.power<0] = 0

    bal_n.power = bal.power
    bal_n.power[bal_n.power>0] = 0
    cnt_p = 0
    cnt_n = 0
    pow_p_mem = -1
    pow_n_mem = -1
    charge_loc_p = []
    charge_loc_n = []
    for pow_p, id_p, pow_n, id_n in zip(bal_p.power, bal_p.index, bal_n.power, bal_p.index):

        if (pow_p == 0) and (pow_p_mem != 0):
            charge_loc_p.append(id_p)
            cnt_p = cnt_p + 1

        if (pow_n == 0) and (pow_n_mem != 0):
            charge_loc_n.append(id_n)
            cnt_n = cnt_n + 1

        pow_p_mem = pow_p
        pow_n_mem = pow_n
    
    #for i in range(5):
    #    plt.figure
    #    plt.plot(bal_p.time.iloc[charge_loc_p[i]:charge_loc_p[i+1]],
    #            bal_p.power.iloc[charge_loc_p[i]:charge_loc_p[i+1]],
    #            bal_n.time.iloc[charge_loc_p[i]:charge_loc_p[i+1]],
    #            bal_n.power.iloc[charge_loc_p[i]:charge_loc_p[i+1]])
    #    plt.xlabel('Time (s)')
    #    plt.ylabel('Power (W)')
    #    plt.show()

    ch_loc = pd.DataFrame(None, columns=["time","power"])
    ch_loc.time = bal.time.iloc[charge_loc_p]
    ch_loc.power = 0

    ints_p = []
    ints_n = []
    start_locs = charge_loc_p[0:-2]
    end_locs = charge_loc_p[1:-1]
    t_0 = np.array(bal.time.iloc[start_locs].tolist())
    t_1 = np.array(bal.time.iloc[end_locs].tolist())
    for loc_0, loc_1 in zip(start_locs, end_locs):
        ints_p.append(trapezoid(bal_p.power.iloc[loc_0:loc_1],
                                bal_p.time.iloc[loc_0:loc_1]))
        ints_n.append(trapezoid(bal_n.power.iloc[loc_0:loc_1],
                                bal_n.time.iloc[loc_0:loc_1]))
    
    cycles = pd.DataFrame(None, columns=["id_0","id_1","t_0","t_1","dt","int_p","int_n"])
    cycles.id_0 = start_locs
    cycles.id_1 = end_locs
    cycles.t_0 = t_0
    cycles.t_1 = t_1
    cycles.dt = t_1 - t_0
    cycles.int_p = ints_p
    cycles.int_n = ints_n

    return bal_p, bal_n, ch_loc, cycles


def main():
    """Main function to run transient sizing method."""

    if args.scale is None:
        scale = 1e5
    else:
        scale = float(args.scale)

    pfolder = os.path.join(os.path.dirname(__file__),"data")
    fname = "UK_gridwatch_wind_2020.csv"
    data = read_wind_data(pfolder, fname)
    trns = raw2trans(data, scale)
    bsld = find_baseload(trns)
    bal = find_balance(trns,bsld)
    bal_p, bal_n, ch_loc, cycles = find_charges(bal)

    plt.figure
    plt.plot(bal_p.time,
             bal_p.power/1e3,
             bal_n.time,
             bal_n.power/1e3,
             ch_loc.time,
             ch_loc.power/1e3,'+')
    plt.legend(['Charging','Discharging','Cycle Bounds'])
    plt.xlabel('Time (s)')
    plt.ylabel('Power (kW)')
    plt.title("Source Capacity: %3.1f kW \n"
              "Largest Charging Storage: %3.1f kWh \n"
              "Largest Discharging Storage: %3.1f kWh \n"
              % (scale/1e3,
              cycles.int_p.max()/3.6e6,
              cycles.int_n.min()/3.6e6
                )
             )
    plt.show()

    plt.figure
    plt.plot(cycles.dt, cycles.int_p/3.6e6, '+',
             cycles.dt, cycles.int_n/3.6e6, '+')
    plt.legend(['Charging','Discharging'])
    plt.xlabel('Duration (s)')
    plt.ylabel('Energy Stored (kWh)')
    plt.title("Source Capacity: %3.1f kW \n"
              "Largest Charging Cycle: %3.1f kWh \n"
              "Largest Discharging Cycle: %3.1f kWh \n"
              % (scale/1e3,
              cycles.int_p.max()/3.6e6,
              cycles.int_n.min()/3.6e6
                )
             )
    plt.show()

    
    

if __name__ == "__main__":
    main()
