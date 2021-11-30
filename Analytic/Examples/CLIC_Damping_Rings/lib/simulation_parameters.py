import numpy as np
from scipy.constants import c
from scipy.constants import physical_constants

r0 = physical_constants["classical electron radius"][0]                     #Electron Classical Radius [m]
E0 = physical_constants["electron mass energy equivalent in MeV"][0]*1e-3   #Electron Rest Mass [Gev]


h   = 2852       #Harmonic Number
Z   = 1.0        #Number of Charges
Nb  = 4.4e9     #Bunch Population
a_p = 13e-5      #Momentum Compaction Factor
En  = 2.86       #Electron Total Energy [GeV]

Ek  = En - E0    #Kinetic Energy  [GeV]
P0  = np.sqrt(En**2 - E0**2) / c #Reference momentum [GeV/c]

U0  = 0. #0.0038        #Energy Loss / Turn [Gev]
#p_increment = U0 * 1e+9 * 1.602 * 1e-19 / c    # kg*m/s //turn

gamma_rel = En / E0
beta_rel  = np.sqrt(1 - 1/gamma_rel**2)
gamma_tt  = 1 / np.sqrt(a_p)  #Gamma Transition
SlipF     = (1/gamma_tt**2) - (1 / gamma_rel**2)  #Slip Factor

V0 = 0.0045     #RF-Voltage [GV]
#-----------------------------------

Nemit_x0 = 5.6644e-7 
Nemit_y0 = 3.7033e-9 
emit_x0  = Nemit_x0 / beta_rel / gamma_rel
emit_y0  = Nemit_y0 / beta_rel / gamma_rel
#Sigma_E0 = 1.086625e-03
bl0 = 0.00158
