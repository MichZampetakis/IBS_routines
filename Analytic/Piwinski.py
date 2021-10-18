import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d
from scipy.constants import c, epsilon_0, e, hbar, m_p
from scipy.constants import physical_constants
from general_functions import *

class PiwinskiIBS():

    def __init__(self, *args, **kwargs):
        pass

    def set_beam_parameters(self, Ncharges, Npart, Total_Energy):
        E0p = physical_constants["proton mass energy equivalent in MeV"][0]*1e-3
        self.Npart  = Npart
        self.Ncharg = Ncharges
        self.EnTot  = Total_Energy
        self.E_rest = 193.6999
        self.gammar = self.EnTot / self.E_rest
        self.betar  = np.sqrt(1 - 1/self.gammar**2)
        self.mass_i = (self.E_rest * m_p) / E0p   #Ions Mass
        self.c_rad  = (self.Ncharg * e)**2 / (4 * np.pi * epsilon_0 * c**2 * self.mass_i)

    def set_optic_functions(self, position, dels, bet_x, bet_y, alf_x, alf_y, eta_x, eta_dx, eta_y, eta_dy):
        self.posit  = position
        self.Circu  = position[len(position)-1]
        self.bet_x  = bet_x
        self.bet_y  = bet_y
        self.alf_x  = alf_x
        self.alf_y  = alf_y
        self.eta_x  = eta_x
        self.eta_dx = eta_dx
        self.eta_y  = eta_y
        self.eta_dy = eta_dy
        self.alf_x  = alf_x
        self.H_x    = self._Hi(bet_x, alf_x, eta_x, eta_dx)
        self.H_y    = self._Hi(bet_y, alf_y, eta_y, eta_dy)
        self.dels   = dels
        self.frev   = self.betar * c / position[len(position)-1]
        self.x_set  = np.linspace(0,1,50)
        self.dimt   = len(position)

    def _Hi(self, beta , alpha , eta , eta_d):
        return (1 / beta) * ( eta**2 + (beta * eta_d + alpha * eta)**2 )

    def Sigma(self, beta, emit, eta, Sigma_M):
        return np.sqrt( beta * emit + (eta * Sigma_M)**2 )

    def SigmaY(self, Ey, Sigma_M):
        return np.sqrt(self.bet_y * Ey + (self.eta_y * Sigma_M)**2)

    def Aconst(self, Ex, Ey, Sig_M, bl):
        return self.c_rad**2 * c * self.Npart / (64 * np.pi**2 * self.gammar**4 * self.betar**3 * Ex * Ey * bl * Sig_M)

    def sigmaH_inv(self, Ex, Ey, Sig_M):
        self.sigH_inv = 1 / np.sqrt((1/Sig_M**2) + (self.H_x / Ex) + (self.H_y / Ey))

    def a_b_c(self, Ex, Ey, Sig_M):
        self.sigmaH_inv(Ex, Ey, Sig_M)
        sigY = self.SigmaY(Ey, Sig_M)

        self.a1 = self.sigH_inv * np.sqrt(self.bet_x / Ex) / self.gammar
        self.b1 = self.sigH_inv * np.sqrt(self.bet_y / Ey) / self.gammar
        self.c1 = self.sigH_inv * np.sqrt(2. * sigY / self.c_rad) * self.betar

    def PQ(self, ab):
        PQ_list = []
        for x in self.x_set:
            PQ_list.append(np.sqrt(ab**2 + (1 - ab**2) * x**2))
        return PQ_list

    def Ffunc(self, aa, bb, cc):
        FF = []
        P = self.PQ(aa)
        Q = self.PQ(bb)
        for i in xrange(len(self.x_set)):
            FF.append(8 * np.pi * (1 - 3 * self.x_set[i]**2) * (2 * np.log(cc * (1/P[i] + 1/Q[i]) / 2.) - np.euler_gamma) / (P[i] * Q[i]))
        return FF

    def calculate_Ffunctions(self):
        F1bc = self.Ffunc(1/self.a1, self.b1/self.a1, self.c1/self.a1)
        F1ac = self.Ffunc(1/self.b1, self.a1/self.b1, self.c1/self.b1)
        Fabc = self.Ffunc(self.a1, self.b1, self.c1)

        invF1bc = map(list, zip(*F1bc))
        invF1ac = map(list, zip(*F1ac))
        invFabc = map(list, zip(*Fabc))

        self.IF1bc = integrate.simps(invF1bc, self.x_set)
        self.IF1ac = integrate.simps(invF1ac, self.x_set)
        self.IFabc = integrate.simps(invFabc, self.x_set)

    def interp_integr(self, gx, gy, gs):
        int_gx = interp1d(self.posit, gx)
        int_gy = interp1d(self.posit, gy)
        int_gs = interp1d(self.posit, gs)

        TTx = integrate.quad(int_gx, self.posit[0], self.posit[-1])[0] / self.Circu
        TTy = integrate.quad(int_gy, self.posit[0], self.posit[-1])[0] / self.Circu
        TTs = integrate.quad(int_gs, self.posit[0], self.posit[-1])[0] / self.Circu

        return TTx, TTy, TTs

    def dels_integr(self, gx, gy, gs):
        TTx = np.sum(gx * self.dels) / self.Circu
        TTy = np.sum(gy * self.dels) / self.Circu
        TTs = np.sum(gs * self.dels) / self.Circu

        return TTx, TTy, TTs

    def calculate_integrals(self, Ex, Ey, Sig_M, bl):
        Acon = self.Aconst(Ex, Ey, Sig_M, bl)

        gx = Acon * (self.IF1bc + self.H_x * self.sigH_inv**2 * self.IFabc / Ex)
        gy = Acon * (self.IF1ac + self.H_y * self.sigH_inv**2 * self.IFabc / Ey)
        gs = Acon * self.sigH_inv**2 * self.IFabc / Sig_M**2

        self.Tx, self.Ty, self.Ts = self.interp_integr(gx, gy, gs)
#         self.Tx, self.Ty, self.Ts = self.dels_integr(gx, gy, gs)

    def calculate_growth_times(self, Ex, Ey, Sig_M, bl):
        self.a_b_c(Ex, Ey, Sig_M)
        self.calculate_Ffunctions()
        self.calculate_integrals(Ex, Ey, Sig_M, bl)


    def emit_evol(self, Ex, Ey, Sig_M, bl, dt):
        #self.calculate_growth_times(Ex, Ey, Sig_M, bl)
        #print self.Tx, self.Ty, self.Ts
        Evolemx = Ex * np.exp(2. * dt * float(self.Tx))
        Evolemy = Ey * np.exp(2. * dt * float(self.Ty))
        EvolsiM = Sig_M  * np.exp(1. * dt * float(self.Ts))

        return Evolemx, Evolemy, EvolsiM

    # Run if you want to evaluate the emittance evolution using Nagaitsev's Integrals, including Synchrotron Radiation.
    # Give damping times in [s], not turns!!!
    def emit_evol_with_SR(self, Emit_x, Emit_y, Sig_M, BunchL, EQemitX, EQemitY, EQsigmM, tau_x, tau_y, tau_s, dt):
        #Ixx, Iyy, Ipp = self.Nagaitsev_Integrals(Emit_x, Emit_y, Sig_M, BunchL)

        Evolemx = (- EQemitX + np.exp(dt * 2 * (float(self.Ixx) - 1.0 / tau_x)) * (EQemitX + Emit_x * (float(self.Ixx) * tau_x - 1.0))) / (float(self.Ixx) * tau_x - 1.0)
        Evolemy = (- EQemitY + np.exp(dt * 2 * (float(self.Iyy) - 1.0 / tau_y)) * (EQemitY + Emit_y * (float(self.Iyy) * tau_y - 1.0))) / (float(self.Iyy) * tau_y - 1.0)
        EvolsiM = np.sqrt((- EQsigmM**2 + np.exp(dt * 2 * (float(self.Ipp) - 1.0 / tau_s)) * (EQsigmM**2 + Sig_M**2 * (float(self.Ipp) * tau_s - 1.0))) / (float(self.Ipp) * tau_s - 1.0))

        return Evolemx, Evolemy, EvolsiM
