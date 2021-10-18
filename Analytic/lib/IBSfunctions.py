import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d
from scipy.constants import c, epsilon_0, e, hbar, m_p
from scipy.constants import physical_constants
from general_functions import *

class NagaitsevIBS():

    def __init__(self, *args, **kwargs):
        pass

    def _Phi(self, beta , alpha , eta , eta_d):
        return eta_d + alpha * eta / beta

    def set_beam_parameters(self, Ncharges, Npart, Total_Energy):
        E0p = physical_constants["proton mass energy equivalent in MeV"][0]*1e-3
        
        self.Npart  = Npart
        self.Ncharg = Ncharges
        self.EnTot  = Total_Energy
        #self.E_rest = physical_constants["electron mass energy equivalent in MeV"][0]*1e-3
        self.E_rest = 193.6999
        self.gammar = self.EnTot / self.E_rest
        self.betar  = np.sqrt(1 - 1/self.gammar**2)
        self.mass_i = (self.E_rest * m_p) / E0p   #Ions Mass
        self.c_rad  = (self.Ncharg * e)**2 / (4 * np.pi * epsilon_0 * c**2 * self.mass_i)

    def set_optic_functions(self, position, dels, bet_x, bet_y, eta_x, eta_dx, eta_y, alf_x):
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

    def Nagaitsev_Integrals(self, Emit_x, Emit_y, Sig_M, BunchL):
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


    def calculate_integrals(self, Emit_x, Emit_y, Sig_M, BunchL):
        self.Ixx, self.Iyy, self.Ipp = self.Nagaitsev_Integrals(Emit_x, Emit_y, Sig_M, BunchL)
        #print self.Ixx, self.Iyy, self.Ipp


    # Run if you want to evaluate the emittance evolution using Nagaitsev's Integrals.
    def emit_evol(self, Emit_x, Emit_y, Sig_M, BunchL, dt):
        #Ixx, Iyy, Ipp = self.Nagaitsev_Integrals(Emit_x, Emit_y, Sig_M, BunchL)

        Evolemx = Emit_x * np.exp(dt * float(self.Ixx))
        Evolemy = Emit_y * np.exp(dt * float(self.Iyy))
        EvolsiM = Sig_M  * np.exp(dt * float(0.5 * self.Ipp))

        return Evolemx, Evolemy, EvolsiM


    # Run if you want to evaluate the emittance evolution using Nagaitsev's Integrals, including Synchrotron Radiation.
    # Give damping times in [s], not turns!!!
    def emit_evol_with_SR(self, Emit_x, Emit_y, Sig_M, BunchL, EQemitX, EQemitY, EQsigmM, tau_x, tau_y, tau_s, dt):
        #Ixx, Iyy, Ipp = self.Nagaitsev_Integrals(Emit_x, Emit_y, Sig_M, BunchL)

        Evolemx = (- EQemitX + np.exp(dt * 2 * (float(self.Ixx / 2.) - 1.0 / tau_x)) * (EQemitX + Emit_x * (float(self.Ixx / 2.) * tau_x - 1.0))) / (float(self.Ixx / 2.) * tau_x - 1.0)
        Evolemy = (- EQemitY + np.exp(dt * 2 * (float(self.Iyy / 2.) - 1.0 / tau_y)) * (EQemitY + Emit_y * (float(self.Iyy / 2.) * tau_y - 1.0))) / (float(self.Iyy / 2.) * tau_y - 1.0)
        EvolsiM = np.sqrt((- EQsigmM**2 + np.exp(dt * 2 * (float(self.Ipp / 2.) - 1.0 / tau_s)) * (EQsigmM**2 + Sig_M **2 * (float(self.Ipp / 2.) * tau_s - 1.0))) / (float(self.Ipp / 2.) * tau_s - 1.0))

        return Evolemx, Evolemy, EvolsiM

    def emit_evol_kicks(self, part_coord_dict, BunchL):
        Emit_x = emittance(part_coord_dict['x'],part_coord_dict['px'])
        Emit_y = emittance(part_coord_dict['y'],part_coord_dict['py'])
        Sig_M  = np.std(part_coord_dict['pz'])

        p_std_x = np.std(part_coord_dict['px'])
        p_std_y = np.std(part_coord_dict['py'])

        Ixx, Iyy, Ipp = self.Nagaitsev_Integrals(Emit_x, Emit_y, Sig_M, BunchL)

        #rho = self.line_density(40, part_coord_dict)

        DSx = p_std_x * np.sqrt(2 * Ixx / self.frev)
        DSy = p_std_y * np.sqrt(2 * Iyy / self.frev)
        DSz = Sig_M   * np.sqrt(2 * Ipp / self.frev) * self.betar**2

        return DSx, DSy, DSz


    def line_density(self, n_slices, beam_coords):
        z_cut_head = np.max(beam_coords['z'])
        z_cut_tail = np.min(beam_coords['z'])
        slice_width = (z_cut_head - z_cut_tail) / float(n_slices)

        bin_edges = np.linspace(z_cut_tail - 1e-7 * slice_width,  z_cut_head + 1e-7 * slice_width, num=n_slices+1, dtype=np.float64)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.

        bunch_length_rms = np.std(beam_coords['z'])
        factor_distribution = bunch_length_rms * 2 * np.sqrt(np.pi)

        #counts, bin_edges = np.histogram(parameters['z'], bin_edges)
        #Rho = np.interp(beam_coords['z'], bin_centers, counts/slice_width * factor_distribution/Npart)
        #kick_factor = np.mean(Rho)

        counts_normed, bin_edges = np.histogram(beam_coords['z'], bin_edges, density=True)
        Rho_normed = np.interp(beam_coords['z'], bin_centers, counts_normed * factor_distribution)
        kick_factor_normed = np.mean(Rho_normed)

        return Rho_normed

    def calculate_emittance_change(self, part_coord_dict, BunchL):
        self.DSx, self.DSy, self.DSz = self.emit_evol_kicks(part_coord_dict, BunchL)

    def track(self, part_coord_dict):
        rho = self.line_density(40, part_coord_dict)
        Dkick_x = np.random.normal(loc = 0, scale = self.DSx, size = part_coord_dict['px'].shape) * np.sqrt(rho)
        Dkick_y = np.random.normal(loc = 0, scale = self.DSy, size = part_coord_dict['py'].shape) * np.sqrt(rho)
        Dkick_p = np.random.normal(loc = 0, scale = self.DSz, size = part_coord_dict['pz'].shape) * np.sqrt(rho)

        part_coord_dict['px'] += Dkick_x
        part_coord_dict['py'] += Dkick_y
        part_coord_dict['pz'] += Dkick_p

        
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

    def emit_evol_kicks(self, part_coord_dict, BunchL):
        Emit_x = emittance(part_coord_dict['x'],part_coord_dict['px'])
        Emit_y = emittance(part_coord_dict['y'],part_coord_dict['py'])
        Sig_M  = np.std(part_coord_dict['pz'])

        p_std_x = np.std(part_coord_dict['px'])
        p_std_y = np.std(part_coord_dict['py'])

        self.calculate_growth_times(Emit_x, Emit_y, Sig_M, BunchL)
        #Ixx, Iyy, Ipp = self.Nagaitsev_Integrals(Emit_x, Emit_y, Sig_M, BunchL)

        #rho = self.line_density(40, part_coord_dict)

        DSx = p_std_x * np.sqrt(2 * 2 * self.Tx / self.frev)
        DSy = p_std_y * np.sqrt(2 * 2 * self.Ty / self.frev)
        DSz = Sig_M   * np.sqrt(2 * 2 * self.Ts / self.frev) * self.betar**2

        return DSx, DSy, DSz


    def line_density(self, n_slices, beam_coords):
        z_cut_head = np.max(beam_coords['z'])
        z_cut_tail = np.min(beam_coords['z'])
        slice_width = (z_cut_head - z_cut_tail) / float(n_slices)

        bin_edges = np.linspace(z_cut_tail - 1e-7 * slice_width,  z_cut_head + 1e-7 * slice_width, num=n_slices+1, dtype=np.float64)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.

        bunch_length_rms = np.std(beam_coords['z'])
        factor_distribution = bunch_length_rms * 2 * np.sqrt(np.pi)

        #counts, bin_edges = np.histogram(parameters['z'], bin_edges)
        #Rho = np.interp(beam_coords['z'], bin_centers, counts/slice_width * factor_distribution/Npart)
        #kick_factor = np.mean(Rho)

        counts_normed, bin_edges = np.histogram(beam_coords['z'], bin_edges, density=True)
        Rho_normed = np.interp(beam_coords['z'], bin_centers, counts_normed * factor_distribution)
        kick_factor_normed = np.mean(Rho_normed)

        return Rho_normed

    def calculate_emittance_change(self, part_coord_dict, BunchL):
        self.DSx, self.DSy, self.DSz = self.emit_evol_kicks(part_coord_dict, BunchL)

    def track(self, part_coord_dict):
        rho = self.line_density(40, part_coord_dict)
        Dkick_x = np.random.normal(loc = 0, scale = self.DSx, size = part_coord_dict['px'].shape) * np.sqrt(rho)
        Dkick_y = np.random.normal(loc = 0, scale = self.DSy, size = part_coord_dict['py'].shape) * np.sqrt(rho)
        Dkick_p = np.random.normal(loc = 0, scale = self.DSz, size = part_coord_dict['pz'].shape) * np.sqrt(rho)

        part_coord_dict['px'] += Dkick_x
        part_coord_dict['py'] += Dkick_y
        part_coord_dict['pz'] += Dkick_p