import numpy as np
from plant_functions import *
import matplotlib.pyplot as plt
#INPUT VARIABLES

def condensor_test():
    cndsr_det = Condenser_Details()
    state_det = State(1.5, 0.5, 0.2, 350, 200)
    [Tmix,Tcool,NH3,s_out,delpmix,delpcool,dqsum] = tristan_condenser(state_det,cndsr_det)
    print(delpmix, delpcool, dqsum,s_out.yNH3,s_out.NH3,state_det.NH3-s_out.NH3)
    plt.plot(np.linspace(0, 1,cndsr_det.ix+1), Tmix-250)
    plt.plot(np.linspace(0, 1,cndsr_det.ix+1), Tcool-250)
    plt.plot(np.linspace(0, 1,cndsr_det.ix), NH3*120/0.25)
    plt.xlabel('Distance along exchanger')
    plt.ylabel('mix temp')
    plt.show()

def heat_exchanger_test():
    he_det = Heat_Exchanger_Details()
    sh_det = State(1.6, 0.5, 0.25, 770, 200)
    sc_det = State(1.9, 0.6, 0.05, 290, 200)
    [sh_out,sc_out,sh_data,sc_data] = tristan_heat_exchanger(sh_det, sc_det, he_det, last_HE_run=1)
    plt.plot(np.linspace(0, 1,he_det.ix+1), sh_data-273)
    plt.plot(np.linspace(0, 1,he_det.ix+1), sc_data-273)
    plt.plot()
    plt.xlabel('Distance along exchanger')
    plt.ylabel('Temp /K')
    plt.show()

#heat_exchanger_test()
#condensor_test()