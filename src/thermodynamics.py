#! /usr/bin/env python
# -*- coding: utf-8 -*-

#############################################################################
# NiceProp - Interactively learning NICFD
# Authors: ir. A. Giuffre', Dr. ir. M. Pini
# Content: Computation of thermodynamic properties
# 2021 - TU Delft - All rights reserved
#############################################################################


import CoolProp
from plot import *
from tqdm import tqdm

class ThermodynamicModel:
    """ Class to define the thermodynamic properties of the working fluid under analysis """
    def __init__(self, settings, MachEdge=1.0, WriteTecplotData=False):
        """
        Arguments:
        :param settings: dictionary collecting the information read from configuration file

        Optional Arguments:
        :param MachEdge: Mach number at the edge of the BL, used only to compute Gruneisen number
        :param WriteTecplotData: boolean to enable generation of a data file for Tecplot

         Methods:
        - TsContour: compute non-ideality factors in the reduced T-s thermodynamic plane
        - PTContour: compute non-ideality factors in the reduced P-T thermodynamic plane
        """
        self.fluid = settings['fluid']
        self.lib = settings['equation of state']
        self.samples = settings['samples']
        self.process = settings['process']
        self.alpha = settings['alpha']
        self.beta = settings['beta']
        self.labels = settings['labels']
        self.Tr_vec = settings['Tr vec']
        self.Pr_vec = settings['Pr vec']
        self.sr_vec = settings['sr vec']
        self.EoS = CoolProp.AbstractState(self.lib, self.fluid)
        self.Tc = self.EoS.T_critical()
        self.Pc = self.EoS.p_critical()
        self.EoS.update(CoolProp.PT_INPUTS, self.Pc, self.Tc)
        self.sc = self.EoS.smass()
        self.MM = self.EoS.molar_mass()
        self.R = CoolProp.CoolProp.PropsSI('gas_constant', self.lib + '::' + self.fluid) / self.MM
        self.cp0 = CoolProp.CoolProp.PropsSI('CP0MASS', 'T', self.Tc, 'P', self.Pc, self.lib + '::' + self.fluid)
        self.cv0 = self.cp0 - self.R  # Mayer's relation
        self.N = 2 * self.cv0 / self.R
        self.gamma_id = self.cp0 / self.cv0
        self.MachEdge = MachEdge

        # initialize thermodynamic and themo-physical properties of interest
        self.FundDerGamma = np.zeros((self.samples, self.samples))
        self.Z = np.zeros((self.samples, self.samples))
        self.gamma = np.zeros((self.samples, self.samples))
        self.gamma_PT = np.zeros((self.samples, self.samples))
        self.gamma_Pv = np.zeros((self.samples, self.samples))
        self.gamma_Tv = np.zeros((self.samples, self.samples))
        self.Eckert = np.zeros((self.samples, self.samples))
        self.Gruneisen = np.zeros((self.samples, self.samples))
        self.mu = np.zeros((self.samples, self.samples))
        self.k = np.zeros((self.samples, self.samples))
        self.Cv = np.zeros((self.samples, self.samples))
        self.Cp = np.zeros((self.samples, self.samples))

        # compute inlet state and initialize outlet state
        self.Pt_in = np.zeros(len(settings['inlet input 1']))
        self.Tt_in = np.zeros(len(settings['inlet input 1']))
        self.ht_in = np.zeros(len(settings['inlet input 1']))
        self.s_in = np.zeros(len(settings['inlet input 1']))
        self.Dt_in = np.zeros(len(settings['inlet input 1']))
        self.P_out = np.zeros(len(settings['inlet input 1']))
        self.T_out = np.zeros(len(settings['inlet input 1']))
        self.D_out = np.zeros(len(settings['inlet input 1']))
        for ii in range(len(settings['inlet input 1'])):
            if settings['inlet thermodynamic plane'][ii] == 'PT':
                self.Pt_in[ii] = settings['inlet input 1'][ii] * self.Pc
                self.Tt_in[ii] = settings['inlet input 2'][ii] * self.Tc
                self.EoS.update(CoolProp.PT_INPUTS, self.Pt_in[ii], self.Tt_in[ii])
                self.ht_in[ii] = self.EoS.hmass()
            elif settings['inlet thermodynamic plane'][ii] == 'Ts':
                self.Tt_in[ii] = settings['inlet input 1'][ii] * self.Tc
                s_in = settings['inlet input 2'][ii] * self.sc
                self.EoS.update(CoolProp.SmassT_INPUTS, s_in, self.Tt_in[ii])
                self.Pt_in[ii] = self.EoS.p()
                self.ht_in[ii] = self.EoS.hmass()
            else:
                raise ValueError("The inlet thermodynamic state must be specified in 'PT' or 'Ts'")
            self.EoS.update(CoolProp.PT_INPUTS, self.Pt_in[ii], self.Tt_in[ii])
            self.Dt_in[ii] = self.EoS.rhomass()
            self.s_in[ii] = self.EoS.smass()
        self.Tecplot = WriteTecplotData

    def TsContour(self):
        """
        Compute non-ideality factors and normalized transport properties in the reduced T-s thermodynamic plane.

        Notes:
        - self.process: 'compression' or 'expansion'
        - self.alpha: total-to-static volumetric flow ratio defining the outlet state of the transformation.
        It must be equal to zero if using beta.
        - self.beta: total-to-static pressure ratio defining the outlet state of the transformation.
        It must be equal to zero if using alpha.
        - self.labels: list of labels assigned to each thermodynamic transformation
        - self.Tr_vec: array defining the range of reduced T over which the non-ideality factors are computed
        - self.sr_vec: array defining the range of reduced s over which the non-ideality factors are computed
        """
        print("\n Computing non-ideality factors and normalized transport properties "
              "in the reduced T-s thermodynamic plane for: %8s" % self.fluid)
        T_vec = self.Tr_vec * self.Tc
        s_vec = self.sr_vec * self.sc

        Tc_isobaric = np.zeros(self.samples)
        T_isobaric = np.zeros((len(self.Tt_in), self.samples))
        s_sat_v = []
        s_sat_l = []
        T_Widom_line = []
        s_Widom_line = []

        # compute outlet state
        if np.all(self.alpha == 0.0):        # transformation defined by beta
            if self.process == 'expansion':
                plot_process = True
                for ii in range(len(self.Tt_in)):
                    self.P_out[ii] = self.Pt_in[ii] / self.beta[ii]
                    self.EoS.update(CoolProp.PSmass_INPUTS, self.P_out[ii], self.s_in[ii])
                    self.T_out[ii] = self.EoS.T()
            elif self.process == 'compression':
                plot_process = True
                for ii in range(len(self.Tt_in)):
                    self.P_out[ii] = self.Pt_in[ii] * self.beta[ii]
                    self.EoS.update(CoolProp.PSmass_INPUTS, self.P_out[ii], self.s_in[ii])
                    self.T_out[ii] = self.EoS.T()
            else:
                plot_process = False

        elif np.all(self.beta == 0.0):       # transformation defined by alpha
            if self.process == 'expansion':
                plot_process = True
                for ii in range(len(self.Tt_in)):
                    self.D_out[ii] = self.Dt_in[ii] / self.alpha[ii]
                    self.EoS.update(CoolProp.DmassSmass_INPUTS, self.D_out[ii], self.s_in[ii])
                    self.T_out[ii] = self.EoS.T()
                    self.P_out[ii] = self.EoS.p()
            elif self.process == 'compression':
                plot_process = True
                for ii in range(len(self.Tt_in)):
                    self.D_out[ii] = self.Dt_in[ii] * self.alpha[ii]
                    self.EoS.update(CoolProp.DmassSmass_INPUTS, self.D_out[ii], self.s_in[ii])
                    self.T_out[ii] = self.EoS.T()
                    self.P_out[ii] = self.EoS.p()
            else:
                plot_process = False
        else:
            raise ValueError("The thermodynamic transformation must be defined in terms of alpha or beta "
                             "(all entries in the other option must be null)")

        # compute the saturation curve
        for ii in range(self.samples):
            if T_vec[ii] <= self.Tc:
                try:
                    self.EoS.update(CoolProp.QT_INPUTS, 1.0, T_vec[ii])
                    s_sat_v.append(self.EoS.smass())
                except:
                    s_sat_v.append(np.nan)
                try:
                    self.EoS.update(CoolProp.QT_INPUTS, 0.0, T_vec[ii])
                    s_sat_l.append(self.EoS.smass())
                except:
                    s_sat_l.append(np.nan)
        s_sat_v = np.asarray(s_sat_v)
        s_sat_l = np.asarray(s_sat_l)
        s_sat = np.concatenate((s_sat_l, s_sat_v[::-1]))
        T_sat = np.linspace(min(T_vec), self.Tc, len(s_sat_v))
        T_sat = np.concatenate((T_sat, T_sat[::-1]))

        if ((self.process == 'expansion') or (self.process == 'compression')) and (self.lib == 'REFPROP'):
            # compute isobaric lines for plotting purposes
            for ii in range(len(self.Pt_in)):
                for jj in range(self.samples):
                    self.EoS.update(CoolProp.PSmass_INPUTS, self.Pt_in[ii], s_vec[jj])
                    T_isobaric[ii, jj] = self.EoS.T()

        # compute the thermodynamic properties of interest (liquid, vapor and supercritical region)
        s_matrix = np.linspace(min(s_vec), max(s_vec), self.samples)

        for ii in tqdm(range(self.samples)):
            if T_vec[ii] <= self.Tc:
                for jj in range(self.samples):
                    # compute thermo-physical properties in the vapor region
                    if s_matrix[jj] >= s_sat_v[ii]:
                        try:
                            self.EoS.update(CoolProp.SmassT_INPUTS, s_matrix[jj], T_vec[ii])
                            rho = self.EoS.rhomass()
                            P = self.EoS.p()
                            self.FundDerGamma[ii, jj] = self.EoS.fundamental_derivative_of_gas_dynamics()
                            self.Z[ii, jj] = P / (rho * self.R * T_vec[ii])
                            cp = self.EoS.cpmass()
                            cv = self.EoS.cvmass()
                            self.Cp[ii, jj] = cp
                            self.Cv[ii, jj] = cv
                            dP_dT_v = self.EoS.first_partial_deriv(CoolProp.iP, CoolProp.iT, CoolProp.iDmass)
                            dP_dv_T = (- 1 / (rho ** 2) *
                                       self.EoS.first_partial_deriv(CoolProp.iDmass, CoolProp.iP, CoolProp.iT)) ** (-1)
                            dv_dT_P = - 1 / (rho ** 2) * \
                                      self.EoS.first_partial_deriv(CoolProp.iDmass, CoolProp.iT, CoolProp.iP)
                            self.gamma[ii, jj] = cp / cv
                            self.gamma_Tv[ii, jj] = 1 + 1 / (rho * cv) * dP_dT_v
                            self.gamma_Pv[ii, jj] = - 1 / (P * rho) * cp / cv * dP_dv_T
                            self.gamma_PT[ii, jj] = 1 / (1 - P / cp * dv_dT_P)
                            self.Eckert[ii, jj] = self.EoS.speed_sound() ** 2 / (cp * T_vec[ii])
                            self.Gruneisen[ii, jj] = np.sqrt(self.Eckert[ii, jj] * (self.FundDerGamma[ii, jj] - 1) /
                                                             self.MachEdge ** 2)
                            self.mu[ii, jj] = self.EoS.viscosity()
                            self.k[ii, jj] = self.EoS.conductivity()
                        except:
                            self.FundDerGamma[ii, jj] = np.nan
                            self.Z[ii, jj] = np.nan
                            self.gamma[ii, jj] = np.nan
                            self.gamma_PT[ii, jj] = np.nan
                            self.gamma_Pv[ii, jj] = np.nan
                            self.gamma_Tv[ii, jj] = np.nan
                            self.Eckert[ii, jj] = np.nan
                            self.Gruneisen[ii, jj] = np.nan
                            self.mu[ii, jj] = np.nan
                            self.k[ii, jj] = np.nan
                            self.Cp[ii, jj] = np.nan
                            self.Cv[ii, jj] = np.nan
                    # compute thermo-physical properties in the liquid region
                    elif s_matrix[jj] <= s_sat_l[ii]:
                        try:
                            self.EoS.update(CoolProp.SmassT_INPUTS, s_matrix[jj], T_vec[ii])
                            rho = self.EoS.rhomass()
                            P = self.EoS.p()
                            self.FundDerGamma[ii, jj] = self.EoS.fundamental_derivative_of_gas_dynamics()
                            self.Z[ii, jj] = P / (rho * self.R * T_vec[ii])
                            cp = self.EoS.cpmass()
                            cv = self.EoS.cvmass()
                            self.Cp[ii, jj] = cp
                            self.Cv[ii, jj] = cv
                            dP_dT_v = self.EoS.first_partial_deriv(CoolProp.iP, CoolProp.iT, CoolProp.iDmass)
                            dP_dv_T = (- 1 / (rho ** 2) *
                                       self.EoS.first_partial_deriv(CoolProp.iDmass, CoolProp.iP, CoolProp.iT)) ** (-1)
                            dv_dT_P = - 1 / (rho ** 2) * \
                                      self.EoS.first_partial_deriv(CoolProp.iDmass, CoolProp.iT, CoolProp.iP)
                            self.gamma[ii, jj] = cp / cv
                            self.gamma_Tv[ii, jj] = 1 + 1 / (rho * cv) * dP_dT_v
                            self.gamma_Pv[ii, jj] = - 1 / (P * rho) * cp / cv * dP_dv_T
                            self.gamma_PT[ii, jj] = 1 / (1 - P / cp * dv_dT_P)
                            self.Eckert[ii, jj] = self.EoS.speed_sound() ** 2 / (cp * T_vec[ii])
                            self.Gruneisen[ii, jj] = np.sqrt(self.Eckert[ii, jj] * (self.FundDerGamma[ii, jj] - 1) /
                                                             self.MachEdge ** 2)
                            self.mu[ii, jj] = self.EoS.viscosity()
                            self.k[ii, jj] = self.EoS.conductivity()
                        except:
                            self.FundDerGamma[ii, jj] = np.nan
                            self.Z[ii, jj] = np.nan
                            self.gamma[ii, jj] = np.nan
                            self.gamma_PT[ii, jj] = np.nan
                            self.gamma_Pv[ii, jj] = np.nan
                            self.gamma_Tv[ii, jj] = np.nan
                            self.Eckert[ii, jj] = np.nan
                            self.Gruneisen[ii, jj] = np.nan
                            self.mu[ii, jj] = np.nan
                            self.k[ii, jj] = np.nan
                            self.Cp[ii, jj] = np.nan
                            self.Cv[ii, jj] = np.nan
                    else:
                        # two-phase region, no calculation of properties
                        self.FundDerGamma[ii, jj] = np.nan
                        self.Z[ii, jj] = np.nan
                        self.gamma[ii, jj] = np.nan
                        self.gamma_PT[ii, jj] = np.nan
                        self.gamma_Pv[ii, jj] = np.nan
                        self.gamma_Tv[ii, jj] = np.nan
                        self.Eckert[ii, jj] = np.nan
                        self.Gruneisen[ii, jj] = np.nan
                        self.mu[ii, jj] = np.nan
                        self.k[ii, jj] = np.nan
                        self.Cp[ii, jj] = np.nan
                        self.Cv[ii, jj] = np.nan
            # compute thermo-physical properties in the supercritical region
            else:
                for jj in range(self.samples):
                    # if s_matrix[jj] >= self.sc:
                    try:
                        self.EoS.update(CoolProp.SmassT_INPUTS, s_matrix[jj], T_vec[ii])
                        rho = self.EoS.rhomass()
                        P = self.EoS.p()
                        cp = self.EoS.cpmass()
                        cv = self.EoS.cvmass()
                        self.Cp[ii, jj] = cp
                        self.Cv[ii, jj] = cv
                        self.FundDerGamma[ii, jj] = self.EoS.fundamental_derivative_of_gas_dynamics()
                        self.Z[ii, jj] = P / (rho * self.R * T_vec[ii])
                        dP_dT_v = self.EoS.first_partial_deriv(CoolProp.iP, CoolProp.iT, CoolProp.iDmass)
                        dP_dv_T = (- 1 / (rho ** 2) *
                                   self.EoS.first_partial_deriv(CoolProp.iDmass, CoolProp.iP, CoolProp.iT)) ** (-1)
                        dv_dT_P = - 1 / (rho ** 2) * \
                                  self.EoS.first_partial_deriv(CoolProp.iDmass, CoolProp.iT, CoolProp.iP)
                        self.gamma[ii, jj] = cp / cv
                        self.gamma_Tv[ii, jj] = 1 + 1 / (rho * cv) * dP_dT_v
                        self.gamma_Pv[ii, jj] = - 1 / (P * rho) * cp / cv * dP_dv_T
                        self.gamma_PT[ii, jj] = 1 / (1 - P / cp * dv_dT_P)
                        self.Eckert[ii, jj] = self.EoS.speed_sound() ** 2 / (cp * T_vec[ii])
                        self.Gruneisen[ii, jj] = np.sqrt(self.Eckert[ii, jj] * (self.FundDerGamma[ii, jj] - 1) /
                                                         self.MachEdge ** 2)
                        self.mu[ii, jj] = self.EoS.viscosity()
                        self.k[ii, jj] = self.EoS.conductivity()
                    except:
                        self.FundDerGamma[ii, jj] = np.nan
                        self.Z[ii, jj] = np.nan
                        self.gamma[ii, jj] = np.nan
                        self.gamma_PT[ii, jj] = np.nan
                        self.gamma_Pv[ii, jj] = np.nan
                        self.gamma_Tv[ii, jj] = np.nan
                        self.Eckert[ii, jj] = np.nan
                        self.Gruneisen[ii, jj] = np.nan
                        self.mu[ii, jj] = np.nan
                        self.k[ii, jj] = np.nan
                        self.Cp[ii, jj] = np.nan
                        self.Cv[ii, jj] = np.nan

        # compute critical isobaric line
        for jj in range(self.samples):
            self.EoS.update(CoolProp.PSmass_INPUTS, self.Pc, s_vec[jj])
            Tc_isobaric[jj] = self.EoS.T()

        # compute Widom line if supercritical region is accounted for
        T_Widom_line, s_Widom_line = self.ComputeWidomLine( s_matrix )

        # plotting
        print("\n Plotting results ...")
        self.plotClass = Plot(s_matrix, T_vec, self.sc, self.Tc, s_sat, T_sat, self.s_in,
                              self.s_in + 1e-10, self.Tt_in, self.T_out, s_vec, T_isobaric, self.labels, plot_process,
                              xc_iso=s_vec, yc_iso=Tc_isobaric, x_Widom=s_Widom_line, y_Widom=T_Widom_line)
        self.plotClass('output/' + self.fluid + '/Ts')

        print("  Creating contour of fundamental derivative...")
        self.plotClass.PlotTsContour(self.FundDerGamma, 'FundamentalDerivative', r'$\Gamma$ [-]')

        print("  Creating contour of compressibility factor...")
        self.plotClass.PlotTsContour(self.Z, 'Z', '$Z$ [-]')

        print("  Creating contour of gamma...")
        self.plotClass.PlotTsContour(self.gamma, 'gamma', r'$\gamma$ [-]', contour_bounds=(0, 2), powerNorm=True)

        print("  Creating contour of gamma PT...")
        self.plotClass.PlotTsContour(self.gamma_PT, 'gamma_PT', r'$\gamma_{PT}$ [-]')

        print("  Creating contour of gamma Pv...")
        self.plotClass.PlotTsContour(self.gamma_Pv, 'gamma_Pv', r'$\gamma_{Pv}$ [-]', powerNorm=True)

        print("  Creating contour of gamma Tv...")
        self.plotClass.PlotTsContour(self.gamma_Tv, 'gamma_Tv', r'$\gamma_{Tv}$ [-]')

        print("  Creating contour of Eckert number...")
        self.plotClass.PlotTsContour(self.Eckert, 'Eckert', '$Eck$ [-]')

        print("  Creating contour of Gruneisen parameter...")
        self.plotClass.PlotTsContour(self.Gruneisen, 'Gruneisen', '$Gr$ [-]')

        print("  Creating contour of normalized dynamic viscosity...")
        self.plotClass.PlotTsContour(self.mu/1.802/10**-5, 'Viscosity', '$\mu$ [-]', powerNorm=True)

        print("  Creating contour of normalized thermal conductivity...")
        self.plotClass.PlotTsContour(self.k/24.76/10**-3, 'Conductivity', '$\kappa$ [-]', powerNorm=True)

        if self.Tecplot:
            print("  Generating data file for Tecplot...")
            self.plotClass.WriteTsThermoData(self.Z)

        return

    def PTContour(self):
        """
        Compute non-ideality factors and normalized transport properties in the reduced P-T thermodynamic plane.

        Notes:
        - self.process: 'compression' or 'expansion'
        - self.alpha: total-to-static volumetric flow ratio defining the outlet state of the transformation.
        It must be equal to zero if using beta.
        - self.beta: total-to-static pressure ratio defining the outlet state of the transformation.
        It must be equal to zero if using alpha.
        - self.labels: list of labels assigned to each thermodynamic transformation
        - self.Tr_vec: array defining the range of reduced T over which the non-ideality factors are computed
        - self.Pr_vec: array defining the range of reduced P over which the non-ideality factors are computed
        """
        print("\n Computing non-ideality factors and normalized transport properties "
              "in the reduced P-T thermodynamic plane for: %8s" % self.fluid)
        T_grid = np.zeros((self.samples, self.samples))
        P_grid = np.zeros((self.samples, self.samples))
        T_sat = np.linspace(self.Tc * self.Tr_vec[0], self.Tc, int(self.samples / 2))
        T_isentropic = np.zeros((len(self.Tt_in), self.samples))
        P_isentropic = np.zeros((len(self.Tt_in), self.samples))
        P_sat = np.zeros(int(self.samples / 2))
        Pc_isobaric = np.zeros(self.samples)
        Tc_isobaric = np.zeros(self.samples)
        P_Widom_line = []
        T_Widom_line = []

        # compute outlet state
        if np.all(self.alpha == 0.0):  # transformation defined by beta
            if self.process == 'expansion':
                plot_process = True
                for ii in range(len(self.Tt_in)):
                    self.P_out[ii] = self.Pt_in[ii] / self.beta[ii]
                    self.EoS.update(CoolProp.PSmass_INPUTS, self.P_out[ii], self.s_in[ii])
                    self.T_out[ii] = self.EoS.T()
            elif self.process == 'compression':
                plot_process = True
                for ii in range(len(self.Tt_in)):
                    self.P_out[ii] = self.Pt_in[ii] * self.beta[ii]
                    self.EoS.update(CoolProp.PSmass_INPUTS, self.P_out[ii], self.s_in[ii])
                    self.T_out[ii] = self.EoS.T()
            else:
                plot_process = False

        elif np.all(self.beta == 0.0):  # transformation defined by alpha
            if self.process == 'expansion':
                plot_process = True
                for ii in range(len(self.Tt_in)):
                    self.D_out[ii] = self.Dt_in[ii] / self.alpha[ii]
                    self.EoS.update(CoolProp.DmassSmass_INPUTS, self.D_out[ii], self.s_in[ii])
                    self.T_out[ii] = self.EoS.T()
                    self.P_out[ii] = self.EoS.p()
            elif self.process == 'compression':
                plot_process = True
                for ii in range(len(self.Tt_in)):
                    self.D_out[ii] = self.Dt_in[ii] * self.alpha[ii]
                    self.EoS.update(CoolProp.DmassSmass_INPUTS, self.D_out[ii], self.s_in[ii])
                    self.T_out[ii] = self.EoS.T()
                    self.P_out[ii] = self.EoS.p()
            else:
                plot_process = False
        else:
            raise ValueError("The thermodynamic transformation must be defined in terms of alpha or beta "
                             "(all entries in the other option must be null)")

        for jj in range(len(T_sat)):
            self.EoS.update(CoolProp.QT_INPUTS, 1.0, T_sat[jj])
            P_sat[jj] = self.EoS.p()
        P_grad = np.linspace(self.Pc, self.Pc * self.Pr_vec[-1], int(self.samples / 2) + 1)

        for jj in range(self.samples):
            # vapor region
            if jj < (self.samples / 2):
                T_grid[jj, jj:self.samples] = np.linspace(T_sat[jj] * 1.0001, self.Tr_vec[-1] * self.Tc,
                                                          self.samples - jj)  # 1.0001
                P_grid[jj, jj:self.samples] = P_sat[jj]
            # supercritical region
            else:
                T_grid[jj, int(self.samples / 2):self.samples] = np.linspace(self.Tc, self.Tr_vec[-1] * self.Tc,
                                                                             int(self.samples / 2))
                P_grid[jj, int(self.samples / 2):self.samples] = P_grad[jj - int(self.samples / 2) + 1]

        if (self.process == 'expansion') or (self.process == 'compression'):
            # compute isentropic lines for plotting purposes
            for ii in range(len(self.Pt_in)):
                T_isentropic[ii, :] = np.linspace(self.T_out[ii], self.Tt_in[ii], self.samples)
                for jj in range(self.samples):
                    self.EoS.update(CoolProp.SmassT_INPUTS, self.s_in[ii], T_isentropic[ii, jj])
                    P_isentropic[ii, jj] = self.EoS.p()

        # compute the thermodynamic properties of interest
        for ii in tqdm(range(self.samples)):
            for jj in range(self.samples):
                if P_grid[ii, jj] > 0.0:
                    self.EoS.update(CoolProp.PT_INPUTS, P_grid[ii, jj], T_grid[ii, jj])
                    rho = self.EoS.rhomass()
                    cp = self.EoS.cpmass()
                    cv = self.EoS.cvmass()
                    self.Cp[ii, jj] = cp
                    self.Cv[ii, jj] = cv
                    self.FundDerGamma[ii, jj] = self.EoS.fundamental_derivative_of_gas_dynamics()
                    self.Z[ii, jj] = P_grid[ii, jj] / (rho * self.R * T_grid[ii, jj])
                    dP_dT_v = self.EoS.first_partial_deriv(CoolProp.iP, CoolProp.iT, CoolProp.iDmass)
                    dP_dv_T = (- 1 / (rho ** 2) *
                               self.EoS.first_partial_deriv(CoolProp.iDmass, CoolProp.iP, CoolProp.iT)) ** (-1)
                    dv_dT_P = - 1 / (rho ** 2) * self.EoS.first_partial_deriv(CoolProp.iDmass, CoolProp.iT, CoolProp.iP)
                    self.gamma[ii, jj] = cp / cv
                    self.gamma_Tv[ii, jj] = 1 + 1 / (rho * cv) * dP_dT_v
                    self.gamma_Pv[ii, jj] = - 1 / (P_grid[ii, jj] * rho) * cp / cv * dP_dv_T
                    self.gamma_PT[ii, jj] = 1 / (1 - P_grid[ii, jj] / cp * dv_dT_P)
                    self.Eckert[ii, jj] = self.EoS.speed_sound() ** 2 / (cp * T_grid[ii, jj])
                    self.Gruneisen[ii, jj] = np.sqrt(self.Eckert[ii, jj] * (self.FundDerGamma[ii, jj] - 1) /
                                                     self.MachEdge ** 2)
                    self.mu[ii, jj] = self.EoS.viscosity()
                    self.k[ii, jj] = self.EoS.conductivity()
                else:
                    self.FundDerGamma[ii, jj] = np.nan
                    self.Z[ii, jj] = np.nan
                    self.gamma[ii, jj] = np.nan
                    self.gamma_PT[ii, jj] = np.nan
                    self.gamma_Pv[ii, jj] = np.nan
                    self.gamma_Tv[ii, jj] = np.nan
                    self.Eckert[ii, jj] = np.nan
                    self.Gruneisen[ii, jj] = np.nan
                    self.mu[ii, jj] = np.nan
                    self.k[ii, jj] = np.nan
                    self.Cp[ii, jj] = np.nan
                    self.Cv[ii, jj] = np.nan

        # compute critical isobaric line
        Pc_isobaric[:] = self.Pc
        Tc_isobaric = np.linspace(self.Tc, self.Tr_vec[-1] * self.Tc, self.samples)

        # compute Widom line if supercritical region is accounted for
        P_Widom_line, T_Widom_line = self.ComputePTWidomLine(T_grid, P_grid)

        # plotting
        print("\n Plotting results ...")
        self.plotClass = Plot(T_grid, P_grid, self.Tc, self.Pc, T_sat, P_sat, self.Tt_in, self.T_out,
                              self.Pt_in, self.P_out, T_isentropic, P_isentropic, self.labels, plot_process,
                              xc_iso=Tc_isobaric, yc_iso=Pc_isobaric, x_Widom=T_Widom_line, y_Widom=P_Widom_line)

        self.plotClass('output/' + self.fluid + '/PT')

        print("  Creating contour of fundamental derivative...")
        self.plotClass.PlotPTContour(self.FundDerGamma, 'FundamentalDerivative', r'$\Gamma$ [-]')

        print("  Creating contour of compressibility factor...")
        self.plotClass.PlotPTContour(self.Z, 'Z', '$Z$ [-]')

        print("  Creating contour of gamma...")
        self.plotClass.PlotPTContour(self.gamma, 'gamma', r'$\gamma$ [-]')

        print("  Creating contour of gamma_PT...")
        self.plotClass.PlotPTContour(self.gamma_PT, 'gamma PT', r'$\gamma_{PT}$ [-]')

        print("  Creating contour of gamma_Pv...")
        self.plotClass.PlotPTContour(self.gamma_Pv, 'gamma Pv', r'$\gamma_{Pv}$ [-]')

        print("  Creating contour of gamma_Tv...")
        self.plotClass.PlotPTContour(self.gamma_Tv, 'gamma Tv', r'$\gamma_{Tv}$ [-]')

        print("  Creating contour of Eckert number...")
        self.plotClass.PlotPTContour(self.Eckert, 'Eckert', '$Eck$ [-]')

        print("  Creating contour of Gruneisen parameter...")
        self.plotClass.PlotPTContour(self.Gruneisen, 'Gruneisen', '$Gr$ [-]')

        print("  Creating contour of normalized dynamic viscosity...")
        self.plotClass.PlotPTContour(self.mu/1.802/10**-5, 'Viscosity', '$\mu$ [-]', powerNorm=True)

        print("  Creating contour of normalized thermal conductivity...")
        self.plotClass.PlotPTContour(self.k/24.76/10**-3, 'Conductivity', '$\kappa$ [-]', powerNorm=True)

        if self.Tecplot:
            print("  Generating data file for Tecplot...")
            self.plotClass.WritePTThermoData(self.Z)

        return

    def ComputeWidomLine(self, s_matrix ):
        """
        Compute Widom line in the reduced T-s thermodynamic plane.

        Notes:
        """
        print("\n Computing Widom line in the T-s thermodynamic plane for: %8s" % self.fluid)

        T_vec = self.Tr_vec * self.Tc

        Cp_Widom_line = []
        T_Widom_line = []
        s_Widom_line = []

        for ii in tqdm(range(self.samples)):
            if T_vec[ii] > self.Tc:
                Cp_Widom_line.append( np.nanmax(self.Cp[ii, :]) )
                an_array = np.array( self.Cp[ii, :] )
                max_index = np.nanargmax(an_array)
                T_Widom_line.append( T_vec[ii] )
                s_Widom_line.append( s_matrix[max_index] )

        Cp_Widom_line = np.asarray(Cp_Widom_line)
        T_Widom_line = np.asarray(T_Widom_line)
        s_Widom_line = np.asarray(s_Widom_line)

        # for ii in range(len(Cp_Widom_line)):
        #     print(" Widom line (T,s,Cp): %8f " % T_Widom_line[ii], s_Widom_line[ii], Cp_Widom_line[ii])

        return T_Widom_line, s_Widom_line

    def ComputePTWidomLine(self, T_grid, P_grid ):
        """
        Compute Widom line in the reduced P-T thermodynamic plane.

        Notes:
        """
        print("\n Computing Widom line in the P-T thermodynamic plane for: %8s" % self.fluid)

        Cp_Widom_line = []
        P_Widom_line = []
        T_Widom_line = []

        for ii in tqdm(range(self.samples)):
            if P_grid[ii,-1] > self.Pc:
                Cp_Widom_line.append( np.nanmax(self.Cp[ii, :]) )
                an_array = np.array( self.Cp[ii, :] )
                max_index = np.nanargmax(an_array)
                P_Widom_line.append( P_grid[ii,-1] )
                T_Widom_line.append( T_grid[ii,max_index] )

        Cp_Widom_line = np.asarray(Cp_Widom_line)
        P_Widom_line = np.asarray(P_Widom_line)
        T_Widom_line = np.asarray(T_Widom_line)

        # for ii in range(len(Cp_Widom_line)):
        #     print(" Widom line (P,T,Cp): %8f " % P_Widom_line[ii], T_Widom_line[ii], Cp_Widom_line[ii])

        return P_Widom_line, T_Widom_line
