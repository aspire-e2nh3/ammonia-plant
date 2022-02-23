import numpy as np
from isobaric_reactor import *
#INPUT VARIABLES

class Condenser_Details(object):

    def __init__(self):
        pie = 3.141592658
        self.r1 = 0.002 #inner rad of inner ammonia pipe [m]
        self.r2 = 0.003 #outer rad of inner pipe [m]
        self.r3 = 0.02 #inner rad of outer coolant pipe [m]
        self.hyd = 4 * (self.r3 ** 2 - self.r2 ** 2) * pie / (2 * pie * (self.r2 + self.r3)) #hydraulic diameter of cooling channel [m]
        self.Length = 5 #length of heat exchanger [m]
        self.numb = 100 #number of counterflow heat exchangers
        self.ix = 500 #number of elements along heat exchanger
        self.jints = 100 #max number of iterations
        self.kval = 16 #thermal conductivity of pipe [W/mK]
        '''
        self.mwnh3 = 0.017 #molar weight of ammonia [kg/mol]
        self.mwn2 = 0.028 #molar weight of nitrogen [kg/mol]
        self.mwh2 = 0.02 #molar weight of hydrogen [kg/mol]
        '''
        self.dhvap = 23370 #average value for enthalpy of condensation [J/mol]
        self.mcool = 0.3 #mass flow of coolant [kg/s]
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
        '''
        self.mnh3 = 0.00347 #mass flow of ammonia into condenser [kg/s]
        self.mh2 = 0.0198 #mass flow of hydrogen into condenser [kg/s]
        self.mn2 = 0.092 #mass flow of nitrogen into condenser [kg/s]
        self.nnh3 = mnh3 / mwnh3 #molar flow rate of ammonia [moles/s]
        self.nh2 = mh2 / mwh2 #molar flow rate of hydrogen [moles/s]
        self.nn2 = mn2 / mwn2 #molar flow rate of nitrogen [moles/s]
        '''
        self.knh3 = 0.026 #thermal conductivity of ammonia [W/mK]
        self.vnh3 = 0.00001 #dynamic viscosity of ammonia [Pas]
        self.cpnh3 = 2600 #specific heat capacity of ammonia [J/kgK]
        self.rnh3 = 8.2 #density of ammonia [kg/m3]
        self.kn2 = 0.028 #thermal conductivity of nitrogen [W/mK]
        self.vn2 = 0.000018 #dynamic viscosity of nitrogen [Pas]
        self.cpn2 = 1113 #specific heat capacity of nitrogen [J/kgK]
        self.rn2 = 51 #density of nitrogen [kg/m3]
        self.kh2 = 0.2 #thermal conductivity of hydrogen [W/mK]
        self.vh2 = 0.000009 #dynamic viscosity of hydrogen [Pas]
        self.cph2 = 14615 #specific heat capacity of hydrogen [J/kgK]
        self.rh2 = 10 #density of hydrogen [kg/m3]
        '''
        self.ptot = 200 #total condenser pressure [bar]
        
        self.xnh3 = nnh3 / (nnh3 + nh2 + nn2) #determine molar fraction of ammonia
        self.pnh3 = ptot * xnh3 #determine partial pressure of ammonia
        
        self.Tmixin = 350 #set mixture inlet temperature to cooler/condenser [K]
        '''
        self.abar = 239.69 #constant for sat curve of nh3
        self.bbar = 0.0964 #constant for sat curve of nh3





def tristan_condenser(s, c):

    #DECLARE ARRAYS
    Tmix = np.full(c.ix,s.T,float)
    Tcool = np.full(c.ix,c.Tcoolin,float)
    dq = np.full(c.ix,0,float)
    savxnh3 = np.full(c.ix,0,float)
    savnnh3 = np.full(c.ix,0,float)
    pie = 3.141592658


     #INITIALISING
    dx = c.Length / c.ix #length of finite segment [m]
    nnh3saved = 0 #initialise convergence variable for ammonia concentration
    Evalue = 0.0000001 #acceptable error in ammonia mole fraction



    #MAIN COMPUTATIONAL LOOP
    j = 0
    stop = False
    while stop == False:
        j += 1
        if j > c.jints:
            stop=True

        #define relaxation factor for smooth convergence
        rstart = 0.1 #RELAX FACTOR =0.1 seems to work well
        rstep = rstart
        minint = (1 - rstart) / rstep
        relax = min(rstart + j * rstep,1)
        s_calc = copy.copy(s)

        #FORWARD PASS ON MIXTURE SIDE OF HEAT EXCHANGER
        nnh3 = s.NH3 #reset molar flow rate of ammonia at inlet [moles/s]
        xnh3 = s.yNH3
        dqsum = 0 #reset heat integration variable


        for ii in range(0, c.ix):
            # Determine if mixture temp low enough for condensation to take place
            pnh3 = s.p * s.yNH3 #partial pressure of ammonia
            Tsat = 239.69 * pnh3 ** 0.0964 #determine saturation temperature of ammonia
            if Tmix[ii] > Tsat:
                conact = 0 #mixture temp above sat temp so no condensation in this element
            elif Tmix[ii] <= Tsat:
                conact = 1 #mixture temp below or equal to sat temp so condensation occuring in this element


            #Use law of mixtures to calculate average fluid properties and overall heat transfer coefficient for an element
            s.NH3 = nnh3
            s.T = Tmix[ii]
            s.update_fast()

            velmix = s.volume/ (c.numb * 2 * pie * c.r1 ** 2)
            reymix = s.rho * 2 * c.r1 * velmix / s.mu
            prmix = s.mu * s.cp / s.k
            nusmix = 0.023 * reymix ** 0.8 * prmix ** 0.4
            htc1 = nusmix * s.k / (2 * c.r1)


            Uval = 1 / (1 / (htc1 * c.r1 * 2 * pie * dx) + np.log(c.r2 / c.r1) / (2 * pie * c.kval * dx) + 1 / (c.htc2 * c.r2 * 2 * pie * dx))

            #determine heat flow from mixture to coolant in that element
            dq[ii] = Uval * (Tmix[ii] - Tcool[ii]) * c.numb * relax
            dqsum = dqsum + dq[ii] #integrate total heat lost to coolant

            #determine molar flow of nh3 lost due to condensation in each element
            nlost = 0
            kints = 100000
            losint = nnh3 / kints
            #increment condensation molar flow until sum of sensible heat reduction in gas + enthalpy of condensation equates to heat transferred out of element

            qout = 0
            # need to change to stop statement.
            while qout<dq[ii]:
                nlost = nlost + losint
                dTdP = c.abar * c.bbar * (s.p * xnh3) ** (c.bbar - 1)
                dP = s.p * (xnh3 - (nnh3 - nlost) / (nnh3 - nlost + s.N2 + s.H2))
                delT = dTdP * dP
                qflow = (s.cp*s.mass_tot) * delT
                qout = qflow + conact * nlost * c.dhvap

            nnh3 = nnh3 - conact * nlost  # reduce molar flow of gaseous ammonia as ammonia is condensed
            s.update_fast() #determine new mole fraction of ammonia after condensation
            savnnh3[ii] = nnh3
            savxnh3[ii] = xnh3
            #reduce temperature of mixture for next element based on heat lost from this element

            Tmix[ii + 1] = Tmix[ii] - delT

        #COUNTER PASS ON COOLANT CHANNEL
        Tcool[c.ix] = c.Tcoolin
        for ii in range(0, c.ix): #loop to update coolant temperatures
            Icon = c.ix - ii
            Tcool[Icon] = Tcool[Icon + 1] + dq[Icon] / (c.mcool * c.cpcool) #increase temp of coolant in line with heat flow from mixture


        #convergence check - has exit mole fraction changed from last main loop iteration
        if (abs(nnh3 - nnh3saved) < Evalue) & (j > 2 * minint): #if loop to check for convergence
            stop = True
        nnh3saved = nnh3

    s.NH3 = nnh3
    s.T = Tmix[c.ix]
    s.update_fast()
    #CALCULATE PRESSURE DROP IN EACH CHANNEL
    fricmix = 0.3164 * reymix ** -0.25
    delpmix = fricmix * s.rho * velmix ** 2 * c.Length / (2 * 2 * c.r1)
    friccool = 0.3164 * c.reycool ** -0.25
    delpcool = friccool * c.rcool * c.velcool ** 2 * c.Length / (2 * 2 * c.r1)

    #PRINT SELECTED RESULTS TO SPREADSHEET
    '''
    For i = 1 To ix
    ActiveSheet.Cells(i, 8) = i * dx
    ActiveSheet.Cells(i, 9) = Tmix(i)
    ActiveSheet.Cells(i, 10) = Tcool(i)
    ActiveSheet.Cells(i, 11) = savnnh3(i)
    ActiveSheet.Cells(i, 12) = savxnh3(i)
    ActiveSheet.Cells(2, 1) = " mixture velocity"
    ActiveSheet.Cells(3, 1) = " coolant velocity"
    ActiveSheet.Cells(2, 2) = velmix
    ActiveSheet.Cells(3, 2) = velcool
    ActiveSheet.Cells(2, 3) = "m/s"
    ActiveSheet.Cells(3, 3) = "m/s"
    ActiveSheet.Cells(4, 1) = " mixture pressure drop"
    ActiveSheet.Cells(5, 1) = " coolant pressure drop"
    ActiveSheet.Cells(4, 2) = delpmix
    ActiveSheet.Cells(5, 2) = delpcool
    ActiveSheet.Cells(4, 3) = "Pa"
    ActiveSheet.Cells(5, 3) = "Pa"
    ActiveSheet.Cells(6, 1) = "heat out of mix"
    ActiveSheet.Cells(7, 1) = "heat in to coolant"
    ActiveSheet.Cells(7, 2) = mcool * cpcool * (Tcool(1) - Tcool(ix))
    ActiveSheet.Cells(6, 2) = dqsum
    ActiveSheet.Cells(6, 3) = "W"
    ActiveSheet.Cells(7, 3) = "W"
    '''
    return dqsum,delpmix

def main():
    cndsr_det = Condenser_Details()
    state_det = State(1.5, 0.5, 0.2, 325, 200)
    [dqsum,delpmix] = tristan_condenser(state_det,cndsr_det)
    print(dqsum,delpmix)


if __name__ == "__main__":
    main()
