import pylab as plt
import numpy as np

SMALL_SIZE = 12
MEDIUM_SIZE = 16
BIGGER_SIZE = 18
plt.rcParams["font.family"] = "serif"
plt.rc('font', size=BIGGER_SIZE)         # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)   # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)   # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)   # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

NagaiData = np.loadtxt('./Nagaitsev.txt')
nag_trn, nag_epsx, nag_epsy, nag_delt, nag_bunl = np.array(NagaiData[0]), np.array(NagaiData[1]), np.array(NagaiData[2]), np.array(NagaiData[3]), np.array(NagaiData[4])

PiwinData = np.loadtxt('./Piwinski.txt')
piw_trn, piw_epsx, piw_epsy, piw_delt, piw_bunl = np.array(PiwinData[0]), np.array(PiwinData[1]), np.array(PiwinData[2]), np.array(PiwinData[3]), np.array(PiwinData[4])

#madxData = np.loadtxt('./MAD-X.txt')
#mad_trn, mad_epsx, mad_epsy, mad_delt, mad_bunl = np.array(madxData[0]), np.array(madxData[1]), np.array(madxData[2]), np.array(madxData[3]), np.array(madxData[4])

fig, (ax1,ax2,ax3) = plt.subplots(1, 3, figsize=(20,6))

ax1.plot(piw_trn, piw_epsx, 'b', label='Piwinski')
ax1.plot(nag_trn, nag_epsx, 'r', label='Nagaitsev')

ax2.plot(piw_trn, piw_epsy, 'b', label='Piwinski')
ax2.plot(nag_trn, nag_epsy, 'r', label='Nagaitsev')

ax3.plot(piw_trn, piw_delt * 1e3, 'b', label='Piwinski')
ax3.plot(nag_trn, nag_delt * 1e3, 'r', label='Nagaitsev')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ax1.set_ylabel('$\epsilon_x$ [m]')
ax1.set_xlabel('Turns')
ax1.legend()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ax2.set_ylabel('$\epsilon_y$ [m]')
ax2.set_xlabel( 'Turns')
ax2.legend()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ax3.set_ylabel('$\sigma_\delta \;[10^{-3}]$')
ax3.set_xlabel('Turns')
ax3.legend()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plt.tight_layout()
plt.show()
plt.savefig('results_comparison.png')
