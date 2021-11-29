import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.constants import c, hbar, m_e
from scipy.constants import physical_constants

from lib.loadlattice import *
from lib.IBS_Piwinski import *
from lib.IBS_Nagaitsev import *
from lib.Equilibrium_emit import *
from lib.General_functions import *
from lib.simulation_parameters import *

# --------------------------------------------------------- #
# 1)"loadlattice.py" reads the CLIC DRs lattice             #
# --------------------------------------------------------- #
# 2)"IBS_Piwinski.py" includes the evaluation of the IBS    #
#    growth rates using Piwinski's model.                   #
# --------------------------------------------------------- #
# 3)"IBS_Nagaitsev.py" includes the evaluation of the IBS   #
#    growth rates using Nagaitsev's method.                 #
# --------------------------------------------------------- #
# 4)"Equilibrium_emit.py" includes all damping times and    #
#    equilibrium emittances, considering Synchrotron        #
#    radiation.                                             #
# --------------------------------------------------------- #
# 5)"General_functions.py" include functions to evaluate:   #
#    Emittance, sigma, bunchlength and Energy spread.       #
# --------------------------------------------------------- #
# 6)"simulation_parameters.py" includes all beam parameters #
#    like: Energy, emittances, Voltage, Harmonic number,etc #
# --------------------------------------------------------- #

print '\n!-------------------------------------------------------!'
print '!~~~~~~~~~~~~~~~~~~~ Lattice Analysis ~~~~~~~~~~~~~~~~~~!'
print '!-------------------------------------------------------!'
Posi, dels, betx, alfx, etax, etadx, bety, alfy, etay, etady, DimT, Circumference = loadfulllattice()

Phi_rf = 0.
Phi_s  = 0.
frev    = beta_rel * c  / Circumference
omega_0 = 2 * np.pi * frev #[Turns/s]
W_rf    = h * omega_0
T0      = 1 / frev
W_sync  = omega_0 * np.sqrt(h * Z * V0 * abs(SlipF * np.cos(0.))/(2 * np.pi * beta_rel**2 * En))
Qs      = W_sync / omega_0
omega_s = Qs * omega_0
sync_T  = 1 / Qs
lag = 180.0
f_RF = W_rf / (2 * np.pi)

print '\n!-------------------------------------------------------!'
print '!~~~~~~~~~~~~~~~~~~ General Parameters ~~~~~~~~~~~~~~~~~!'
print '!-------------------------------------------------------!'
print 'E0 = {0:3.2e} [GeV] ||  Ek = {1:1.4f} [GeV]  ||  En = {2:3.2f} [GeV]'.format(E0, Ek, En)
print 'V0  = {0:1.1e} [GV]  ||  RCIR = {1:2.2f} [m]  ||  h  = {2:1.0f}'.format(V0, Circumference, h)
print "-------------------------------------------------------"
print 'Np = {1:1.0e}       ||  Z  = {0:2.0f}'.format(Z,Nb)
print "-------------------------------------------------------"
print 'gamma_rel = {0:1.5f}  ||  beta_rel = {1:1.5f}'.format(gamma_rel, beta_rel)
print 'gamma_tt  = {0:1.5f}    ||  SlipF    = {1:1.5f}'.format(gamma_tt, SlipF)
print "-------------------------------------------------------"
print 't_rev = {0:1.4e} [s] ||  f0 = {1:1.1f} [Turns/s] || Ts = {2:1.0f}  '.format(T0, frev, sync_T)

print '\n!-------------------------------------------------------!'
print '!~~~~~~~~~~~~~~~~~~~ Beam parameters ~~~~~~~~~~~~~~~~~~~!'
print '!-------------------------------------------------------!'
#bl0 = BunchLength(Circumference, h, En, SlipF, Sigma_E0, beta_rel, V0, U0, Z)
Sigma_E0 = EnergySpread(Circumference, h, En, SlipF, bl0, beta_rel, V0, U0, Z)
Sigma_M0 = Sigma_E0 / beta_rel**2
emit_s0  = Sigma_M0 * bl0 * En * 1e9
print 'Phys_emit_X = {0:4.6e}  ||  Phys_emit_Y = {1:4.6e}'.format(emit_x0,emit_y0)
print 'Bunch Length  = {0:4.8f}  ||  Longitud_emit = {1:6.2e}'.format(bl0, emit_s0)
print 'Energy Spread = {0:1.5e} ||  Moment Spread = {1:4.6e}'.format(Sigma_E0, Sigma_M0)
'''
print '\n!-------------------------------------------------------!'
print '!~~~~~~~~~~~~~~~~~ Equilibrium Values ~~~~~~~~~~~~~~~~~~!'
print '!-------------------------------------------------------!'
Jx0 = 1.03166127
Jy0 = 1.00076635
Js0 = 1.96992250
tx_s = 2 * En * Circumference / (Jx0 * U0 * c)   # !!! in [s] !!!
ty_s = 2 * En * Circumference / (Jy0 * U0 * c)   # !!! in [s] !!!
ts_s = 2 * En * Circumference / (Js0 * U0 * c)   # !!! in [s] !!!

eqSigE = EnergySpread(Circumference, h, En, SlipF, eqBL, beta_rel, V0, U0, Z)
eqSigM = eqSigE / beta_rel**2

print 'Phys_EQ_emit_X = {0:4.6e}  ||  Phys_EQ_emit_Y = {1:4.6e}'.format(eqEx, eqEy)
print 'Phys_EQ_Sig_M  = {0:4.6e}  ||  Phys_EQ_BunchL = {1:4.6e}'.format(eqSigE, eqBL)
print '------------------------------------------------'
print 'Damping times: tau_x = {0:4.3e} [s]  ||  tau_y = {1:4.3e} [s]  ||  tau_s = {2:4.3e} [s]'.format(tx_s ,ty_s, ts_s)
print 'Damping times: tau_x = {0:1.0f}  [Turns]  ||  tau_y = {1:1.0f}  [Turns]  ||  tau_s = {2:1.0f}   [Turns]'.format(tx_s*frev, ty_s*frev, ts_s*frev)
'''
print '\n!-------------------------------------------------------!'
print '!~~~~~~~~~~~~~~~~~~~ Tracking Module ~~~~~~~~~~~~~~~~~~~!'
print '!-------------------------------------------------------!'
IBS = NagaitsevIBS()
#IBS = PiwinskiIBS()
IBS.set_beam_parameters(Z, Nb, En)
IBS.set_optic_functions(Posi, dels, betx, alfx, etax, etadx, bety, alfy, etay, etady)

Exi, Eyi, bli = emit_x0, emit_y0, bl0
Sig_Ei, Sig_Mi = Sigma_E0, Sigma_M0


beam_evol = []

turns = 16000
dt = 1 / frev
IBS_step = 50

for i in xrange(turns):
    print 'Turn = {}'.format(i)
    beam_evol.append([i, Exi, Eyi, Sig_Mi, bli])

    #!- To include Energy or Intensity variations include the -!#
    #!- command below with the corresponding variations       -!#
    #IBS.set_beam_parameters(Z, Nb, En)

    if (i % IBS_step == 0):
        IBS.calculate_integrals(Exi, Eyi, Sig_Mi, bli)

    Exi, Eyi, Sig_Mi = IBS.emit_evol(Exi, Eyi, Sig_Mi, bli, dt)
    #!- To include Synchrotron radiation change to:           -!#
    #Exi, Eyi, Sig_Mi = emit_evol_with_SR(Exi, Eyi, Sig_Mi, bli, eqEx, eqEy, eqSigM, tx_s, ty_s, ts_s, dt)
    
    Sig_Ei = Sig_Mi * beta_rel**2
    bli = BunchLength(Circumference, h, En, SlipF, Sig_Ei, beta_rel, V0, U0, Z)

np.savetxt(r'./output/Nagaitsev.txt', np.array(beam_evol).T.tolist(), fmt='%e')
#np.savetxt(r'./output/Piwinski.txt', np.array(beam_evol).T.tolist(), fmt='%e')

