# CLIC DRs Example
This is a full example of IBS growth rates evaluation for the case of the CLIC Damping rings.

The "IBS_evolution.py" is the main script that reads and prints all the parameters and then does the iterative process of IBS.

The scripts that are called from the main file are:
* "loadlattice.py" reads the CLIC DRs lattice
* "IBS_Piwinski.py" includes the evaluation of the IBS growth rates using Piwinski's model
* "IBS_Nagaitsev.py" includes the evaluation of the IBS growth rates using Nagaitsev's method
* "Equilibrium_emit.py" includes all damping times and equilibrium emittances, considering Synchrotron radiation
* "General_functions.py" include functions to evaluate: Emittance, sigma, bunchlength and Energy spread
* "simulation_parameters.py" includes all beam parameters like: Energy, emittances, Voltage, Harmonic number,etc

