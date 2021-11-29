import numpy as np

def mean2(numb):
    return np.mean( (numb - np.mean(numb))**2 )

def mean3(numbx , numbpx): 
    return np.mean( (numbx - np.mean(numbx)) * (numbpx - np.mean(numbpx)) )

def emittance(x, px):
    return np.sqrt(mean2(x) * mean2(px) - mean3(x,px)**2)

def Sigma(beta, emit, eta, Sigma_M):
    return np.sqrt( beta * emit + (eta * Sigma_M)**2 )

def BunchLength(Circumferance, Harmonic_Num, Energy_total, SlipF, Sigma_E, beta_rel, RF_Voltage, Energy_loss, Z):
    '~~~ from Wiedermanns book ~~~'
    return Sigma_E * Circumferance * np.sqrt(SlipF * Energy_total / (2 * np.pi * beta_rel * Harmonic_Num * np.sqrt(Z**2 * RF_Voltage**2 - Energy_loss**2)))

def EnergySpread(Circumferance, Harmonic_Num, Energy_total, SlipF, BL, beta_rel, RF_Voltage, Energy_loss, Z):
    '~~~ from Wiedermanns book ~~~'
    return BL / (Circumferance * np.sqrt(SlipF * Energy_total / (2 * np.pi * beta_rel * Harmonic_Num * np.sqrt(Z**2 * RF_Voltage**2 - Energy_loss**2))))
