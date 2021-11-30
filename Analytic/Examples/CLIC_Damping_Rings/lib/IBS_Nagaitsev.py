import numpy as np
import scipy.integrate as integrate
from scipy.constants import c, hbar, m_e
from scipy.constants import physical_constants
from General_functions import *

class NagaitsevIBS():

    def __init__(self, *args, **kwargs):
        pass

    def _Hi(self, beta , alpha , eta , eta_d):
        return (1 / beta) * ( eta**2 + (beta * eta_d + alpha * eta)**2 )      

    def _Phi(self, beta , alpha , eta , eta_d):
        return eta_d + alpha * eta / beta

    def set_beam_parameters(self, Ncharges, Npart, Total_Energy):
        self.Npart  = Npart
        self.Ncharg = Ncharges
        self.EnTot  = Total_Energy
        self.E_rest = physical_constants["electron mass energy equivalent in MeV"][0]*1e-3
        self.gammar = self.EnTot / self.E_rest
        self.betar  = np.sqrt(1 - 1/self.gammar**2)
        self.c_rad  = physical_constants["classical electron radius"][0]

    def set_optic_functions(self, position, dels, bet_x, alf_x, eta_x, eta_dx, bet_y, alf_y, eta_y, eta_dy):
        self.posit  = position
        self.Circu  = position[len(position)-1]
        self.bet_x  = bet_x
        self.bet_y  = bet_y
        self.eta_x  = eta_x
        self.eta_dx = eta_dx
        self.eta_y  = eta_y
        self.alf_x  = alf_x
        self.phi_x  = self._Phi(bet_x, alf_x, eta_x, eta_dx)
        self.dels   = dels
        self.bx_bar = sum(bet_x * dels) / self.Circu
        self.by_bar = sum(bet_y * dels) / self.Circu
        self.dx_bar = sum(eta_x * dels) / self.Circu
        self.dy_bar = sum(eta_y * dels) / self.Circu
        self.frev   = self.betar * c / position[len(position)-1]

    # Evaluate Integrals' constant and MEAN Coulomb's logarithm
    def CoulogConst(self, Emit_x, Emit_y, Sig_M, BunchL):
        Etrans = 5e8 * (self.gammar * self.EnTot - self.E_rest) * (Emit_x / self.bx_bar)
        TempeV = 2.0 * Etrans
        sigxcm = 100 * np.sqrt(Emit_x * self.bx_bar + (self.dx_bar * Sig_M)**2)
        sigycm = 100 * np.sqrt(Emit_y * self.by_bar + (self.dy_bar * Sig_M)**2)
        sigtcm = 100 * BunchL
        volume = 8.0 * np.sqrt(np.pi**3) * sigxcm * sigycm * sigtcm
        densty = self.Npart  / volume
        debyul = 743.4 * np.sqrt(TempeV / densty) / self.Ncharg
        rmincl = 1.44e-7 * self.Ncharg**2 / TempeV
        rminqm = hbar * c * 1e5 / (2.0 * np.sqrt(2e-3 * Etrans * self.E_rest))
        rmin   = max(rmincl,rminqm)
        rmax   = min(sigxcm,debyul)
        coulog = np.log(rmax/rmin)
        Ncon   = self.Npart * self.c_rad**2 * c / (12 * np.pi * self.betar**3 * self.gammar**5 * BunchL)
        return Ncon * coulog

    # Solves Carlson's complete elliptic integrals of
    # the 2nd kind using the duplication theorem.
    def RDiter(self, x, y, z):
        R = []
        for i, j, k in zip(x, y, z):
            x0 = i
            y0 = j
            z0 = k
            if (x0 < 0) and (y0 <= 0) and (z0 <= 0): 
                print 'Elliptic Integral Calculation Failed. Wrong input values!' 
                return
            x = x0
            y = y0
            z = [z0]
            li = []
            Sn = []
            differ = 10e-4
            for n in range(0,1000):
                xi = x
                yi = y 
                li.append( np.sqrt(xi * yi) + np.sqrt(xi * z[n]) + np.sqrt(yi * z[n]) )
                x = (xi + li[n]) / 4.
                y = (yi + li[n]) / 4.
                z.append((z[n] + li[n]) / 4.)
                if ((abs(x - xi)/x0 < differ) and (abs(y - yi)/y0 < differ) and (abs(z[n] - z[n+1])/z0 < differ)): break    
            lim = n 
            mi = (xi + yi + 3 * z[lim]) / 5.
            Cx = 1 - (xi / mi)
            Cy = 1 - (yi / mi)
            Cz = 1 - (z[n] / mi)
            En = max(Cx,Cy,Cz)
            if En >= 1 :
                print 'Something went wrong with En'
                return
            summ = 0    
            for m in range(2,6): Sn.append( (Cx**m + Cy**m + 3 * Cz**m) / (2 * m) )  
            for m in range(0,lim): summ += 1 / (np.sqrt(z[m]) * (z[m] + li[m]) * 4**m)
        
            Ern = 3 * En**6 / (1 - En)**(3/2.)   
            rn = - Sn[2-2]**3 / 10. + 3 * Sn[3-2]**2 / 10. + 3 * Sn[2-2] * Sn[4-2] / 5.
            R.append(3 * summ + (1 + 3 * Sn[2-2] / 7. + Sn[3-2] / 3. + 3 * Sn[2-2]**2 / 22. + 3 * Sn[4-2] / 11. + 3 * Sn[2-2] * Sn[3-2] / 13. + 3 * Sn[5-2] / 13. + rn) / (4**lim * mi**(3/2.))) 
        return R

    # Evaluates the growth rates by solving Carlson's complete elliptic
    # integrals of the 2nd kind using the duplication theorem.
    def _Nagaitsev_Integrals(self, Emit_x, Emit_y, Sig_M, BunchL):
        const = self.CoulogConst(Emit_x, Emit_y, Sig_M, BunchL)
        sigx  = np.sqrt(self.bet_x * Emit_x + (self.eta_x * Sig_M)**2)
        sigy  = np.sqrt(self.bet_y * Emit_y + (self.eta_y * Sig_M)**2)
        ax = self.bet_x / Emit_x
        ay = self.bet_y / Emit_y
        a_s = ax * (self.eta_x**2 / self.bet_x**2 + self.phi_x**2) + 1 / Sig_M**2
        a1  = (ax + self.gammar**2 * a_s) / 2.
        a2  = (ax - self.gammar**2 * a_s) / 2.
        denom = np.sqrt(a2**2 + self.gammar**2 * ax**2 * self.phi_x**2)
        #--------------------------------------------------------------------------------
        l1 = ay
        l2 = a1 + denom    
        l3 = a1 - denom
        #--------------------------------------------------------------------------------
        R1 = self.RDiter(1/l2, 1/l3, 1/l1) / l1
        R2 = self.RDiter(1/l3, 1/l1, 1/l2) / l2
        R3 = 3 * np.sqrt(l1 * l2 / l3) - l1 * R1 / l3 - l2 * R2 / l3
        #--------------------------------------------------------------------------------
        Nagai_Sp  = ( 2 * R1 - R2 * (1 - 3 * a2 / denom) - R3 * (1 + 3 * a2 / denom) ) * 0.5 * self.gammar**2
        Nagai_Sx  = ( 2 * R1 - R2 * (1 + 3 * a2 / denom) - R3 * (1 - 3 * a2 / denom) ) * 0.5
        Nagai_Sxp = 3 * self.gammar**2 * self.phi_x**2 * ax * (R3 - R2) / denom
        #--------------------------------------------------------------------------------
        Ixi = self.bet_x / (self.Circu * sigx * sigy) * (Nagai_Sx + Nagai_Sp * (self.eta_x**2 / self.bet_x**2 + self.phi_x**2) + Nagai_Sxp) 
        Iyi = self.bet_y / (self.Circu * sigx * sigy) * (R2 + R3 - 2 * R1)
        Ipi = Nagai_Sp   / (self.Circu * sigx * sigy)
        #--------------------------------------------------------------------------------
        Ix = integrate.simps(Ixi, self.posit) * const / Emit_x
        Iy = integrate.simps(Iyi, self.posit) * const / Emit_y
        Ip = integrate.simps(Ipi, self.posit) * const / Sig_M**2

        return Ix, Iy, Ip

    # Call to evaluate the growth rates
    def calculate_integrals(self, Emit_x, Emit_y, Sig_M, BunchL):
        self.Ixx, self.Iyy, self.Ipp = self._Nagaitsev_Integrals(Emit_x, Emit_y, Sig_M, BunchL)

    # Evaluate the emittance evolution using Nagaitsev's Integrals.
    def emit_evol(self, Emit_x, Emit_y, Sig_M, BunchL, dt):
        Evolemx = Emit_x * np.exp(dt * float(self.Ixx))
        Evolemy = Emit_y * np.exp(dt * float(self.Iyy))
        EvolsiM = Sig_M  * np.exp(dt * float(0.5 * self.Ipp))

        return Evolemx, Evolemy, EvolsiM

    # Evaluate the emittance evolution using Nagaitsev's Integrals, including Synchrotron Radiation!!!
    # Give damping times in [s], not turns!!!
    def emit_evol_with_SR(self, Emit_x, Emit_y, Sig_M, BunchL, EQemitX, EQemitY, EQsigmM, tau_x, tau_y, tau_s, dt):
        Evolemx = (- EQemitX + np.exp(dt * 2 * (float(self.Ixx / 2.) - 1.0 / tau_x)) 
                   * (EQemitX + Emit_x * (float(self.Ixx / 2.) * tau_x - 1.0))) / (float(self.Ixx / 2.) * tau_x - 1.0)
        Evolemy = (- EQemitY + np.exp(dt * 2 * (float(self.Iyy / 2.) - 1.0 / tau_y)) 
                   * (EQemitY + Emit_y * (float(self.Iyy / 2.) * tau_y - 1.0))) / (float(self.Iyy / 2.) * tau_y - 1.0)
        EvolsiM = np.sqrt((- EQsigmM**2 + np.exp(dt * 2 * (float(self.Ipp / 2.) - 1.0 / tau_s)) 
                           * (EQsigmM**2 + Sig_M**2 * (float(self.Ipp / 2.) * tau_s - 1.0))) / (float(self.Ipp / 2.) * tau_s - 1.0))

        return Evolemx, Evolemy, EvolsiM
