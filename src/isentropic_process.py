#! /usr/bin/env python
# -*- coding: utf-8 -*-

#############################################################################
# NiceProp - Interactively learning NICFD
# Authors: ir. A. Giuffre', Dr. ir. M. Pini
# Content: Isentropic expansion(s)/compression(s) in stationary flow devices
# 2021 - TU Delft - All rights reserved
#############################################################################


import CoolProp
from src.IO import *
from src.plot import *
import scipy.optimize as opt
from scipy.interpolate import UnivariateSpline


class IsentropicFlowModel:
    """ Class to compute isentropic expansion(s)/compression(s) in nozzle(s)/diffuser(s) """
    def __init__(self, thermo, settings):
        """
        Arguments:
        :param thermo: instance of the ThermodynamicModel class
        :param settings: dictionary collecting the information read from configuration file

        Methods:
        - IdealProcess
        - NozzleExpansion
        - DiffuserCompression
        - FindSonicState
        - NozzleGeometry
        - NozzleEquations
        - ConicalDiffuserEquations
        - RadialDiffuserEquations
        """
        self.thermo = thermo
        self.massflow = settings['mass flow rate']
        self.velocity = settings['velocity']
        self.geometry = settings['geometry']

        # initialize arrays
        self.P_vec = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.h_vec = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.ht_vec = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.Pt_vec = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.T_vec = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.D_vec = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.V_vec = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.M_vec = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.A_vec = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.R_vec = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.x_norm = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.R_norm = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.gamma = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.Z_vec = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.Z_mean = np.zeros(len(self.thermo.Tt_in))
        self.gamma_Pv_vec = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.gamma_Pv_mean = np.zeros(len(self.thermo.Tt_in))
        self.thermo.Pt_in_diff = np.zeros(len(self.thermo.Tt_in))
        self.Cp_diff = np.zeros(len(self.thermo.Tt_in))
        self.P_diff = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.T_diff = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.D_diff = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.V_diff = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.alpha_diff = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.M_diff = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.A_diff = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.L_diff = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))
        self.R_diff = np.zeros((len(self.thermo.Tt_in), self.thermo.samples))

    def IdealProcess(self):
        """
        Given the inlet and outlet states, compute flow quantities along an ideal expansion/compression.
        """
        print("\n Computing ideal process(es) ...")

        # iterate over the selected inlet states
        for ii in range(len(self.thermo.Tt_in)):
            if self.thermo.alpha != 0:
                self.D_vec[ii, :] = np.linspace(self.thermo.Dt_in[ii], self.thermo.D_out[ii], self.thermo.samples)
            elif self.thermo.beta != 0:
                self.P_vec[ii, :] = np.linspace(self.thermo.Pt_in[ii], self.thermo.P_out[ii], self.thermo.samples)

            # compute thermodynamic properties along the isentropic expansion
            for jj in range(self.thermo.samples):
                if self.thermo.alpha != 0:
                    self.thermo.EoS.update(CoolProp.DmassSmass_INPUTS, self.D_vec[ii, jj], self.thermo.s_in[ii])
                    self.P_vec[ii, jj] = self.thermo.EoS.p()
                elif self.thermo.beta != 0:
                    self.thermo.EoS.update(CoolProp.PSmass_INPUTS, self.P_vec[ii, jj], self.thermo.s_in[ii])
                    self.D_vec[ii, jj] = self.thermo.EoS.rhomass()

                self.h_vec[ii, jj] = self.thermo.EoS.hmass()
                self.T_vec[ii, jj] = self.thermo.EoS.T()
                self.Z_vec[ii, jj] = self.P_vec[ii, jj] / (self.D_vec[ii, jj] * self.thermo.R * self.T_vec[ii, jj])
                self.Z_mean[ii] = np.mean(self.Z_vec[ii, :])
                self.gamma[ii, jj] = self.thermo.EoS.cpmass() / self.thermo.EoS.cvmass()
                dP_dv_T = (- 1 / (self.D_vec[ii, jj] ** 2) *
                           self.thermo.EoS.first_partial_deriv(CoolProp.iDmass, CoolProp.iP, CoolProp.iT)) ** (-1)
                self.gamma_Pv_vec[ii, jj] = - 1 / (self.P_vec[ii, jj] * self.D_vec[ii, jj]) * \
                                            self.thermo.EoS.cpmass() / self.thermo.EoS.cvmass() * dP_dv_T
                self.gamma_Pv_mean[ii] = np.log(self.P_vec[ii, 0] / self.P_vec[ii, -1]) / \
                                         np.log(self.D_vec[ii, 0] / self.D_vec[ii, -1])

                if self.thermo.process == 'expansion':
                    if jj == 0:
                        self.V_vec[ii, jj] = 0.
                    else:
                        self.V_vec[ii, jj] = np.sqrt(2 * (self.thermo.ht_in[ii] - self.thermo.EoS.hmass()))

                    self.M_vec[ii, jj] = self.V_vec[ii, jj] / self.thermo.EoS.speed_sound()
                    self.ht_vec[ii, jj] = self.h_vec[ii, jj] + self.V_vec[ii, jj] ** 2 / 2
                    self.thermo.EoS.update(CoolProp.HmassSmass_INPUTS, self.ht_vec[ii, jj], self.thermo.s_in[ii])
                    self.Pt_vec[ii, jj] = self.thermo.EoS.p()

            # plotting
            if self.thermo.process == 'expansion':
                self.thermo.plotClass('output/' + self.thermo.fluid + '/expansion')
                self.thermo.plotClass.PlotExpansion(self.thermo.Pt_in, self.thermo.Dt_in, self.P_vec, self.D_vec,
                                                    self.Z_vec, self.gamma_Pv_vec, self.M_vec)
            elif self.thermo.process == 'compression':
                self.thermo.plotClass('output/' + self.thermo.fluid + '/compression')
                self.thermo.plotClass.PlotCompression(self.thermo.Pt_in, self.thermo.Dt_in, self.P_vec, self.D_vec,
                                                      self.Z_vec, self.gamma_Pv_vec)

        return

    def NozzleExpansion(self, switch, k_in=6.0, k_out=3.0, 
                        x_norm=np.array([]), R_norm=np.array([]), R_throat_vec=np.array([])):
        """
        Given the inlet and outlet states, design a convergent-divergent nozzle and compute the isentropic
        flow along it.

        Arguments:
        :param switch: 'design' to let NiceProp design the nozzle for the given expansion process or 'off-design'
        to use a user-defined nozzle geometry, see optional arguments x_norm, R_norm, A_vec.

        Optional arguments:
        :param k_in: inlet-to-throat nozzle length / throat radius
        :param k_out: outlet-to-throat nozzle length / throat radius
        Adjust the previous parameters to get the desired shape of the nozzle.
        Note that some combinations of values may lead to unfeasible nozzle geometry.
        :param x_norm, R_norm, R_throat_vec: optional arguments to compute the isentropic flow along a
        user-defined nozzle geometry. The three vectors must have the same shape.

        Notes:
        - self.geometry: 'rectangular' or 'circular'
        - self.massflow: used to compute the area and the radius variation of the nozzle (it must be >= 0).
        Leave 1.0 if you're interested only in the normalized radius variation.
        - self.velocity: flow velocity at the inlet of the nozzle (it must be >= 0).
        Increase if the inlet nozzle area is too large compared to the sections downstream.
        """
        print("\n Computing isentropic expansion(s) in nozzle(s) ...")

        # iterate over the selected inlet states
        for ii in range(len(self.thermo.Tt_in)):
            # check input values
            if (self.massflow[ii] <= 0) or (self.velocity[ii] <= 0):
                raise ValueError(" The values assigned to 'mass flow rate' and 'velocity' must be larger than zero")

            # compute inlet and throat pressure
            h_in = self.thermo.ht_in[ii] - self.velocity[ii] ** 2 / 2
            self.thermo.EoS.update(CoolProp.HmassSmass_INPUTS, h_in, self.thermo.s_in[ii])
            P_in = self.thermo.EoS.p()
            P_throat = opt.fsolve(self.FindSonicState, (self.thermo.P_out[ii]), args=(ii), full_output=False, xtol=1e-12)

            if switch == 'design':    # design mode: compute nozzle geometry
                if P_throat < self.thermo.P_out[ii]:
                    raise ValueError(" The selected P out doesn't lead to supersonic conditions, please increase the "
                                     "expansion/volumetric flow ratio")

                # compute the inlet radius
                self.thermo.EoS.update(CoolProp.PSmass_INPUTS, P_in, self.thermo.s_in[ii])
                D_in = self.thermo.EoS.rhomass()
                self.velocity[ii] = np.sqrt(2 * (self.thermo.ht_in[ii] - self.thermo.EoS.hmass()))
                A_in = self.massflow[ii] / (D_in * self.velocity[ii])

                # compute the throat radius
                self.thermo.EoS.update(CoolProp.PSmass_INPUTS, P_throat, self.thermo.s_in[ii])
                D_throat = self.thermo.EoS.rhomass()
                V_throat = np.sqrt(2 * (self.thermo.ht_in[ii] - self.thermo.EoS.hmass()))
                A_throat = self.massflow[ii] / (D_throat * V_throat)

                # compute the outlet radius
                self.thermo.EoS.update(CoolProp.PSmass_INPUTS, self.thermo.P_out[ii], self.thermo.s_in[ii])
                D_out = self.thermo.EoS.rhomass()
                V_out = np.sqrt(2 * (self.thermo.ht_in[ii] - self.thermo.EoS.hmass()))
                A_out = self.massflow[ii] / (D_out * V_out)

                # define nozzle geometry
                R_norm = np.zeros(self.thermo.samples)
                if self.geometry == 'rectangular':
                    R_in = A_in / 2
                    R_throat = A_throat / 2
                    R_out = A_out / 2
                elif self.geometry == 'circular':
                    R_in = np.sqrt(A_in / np.pi)
                    R_throat = np.sqrt(A_throat / np.pi)
                    R_out = np.sqrt(A_out / np.pi)
                else:
                    raise ValueError(" The available choices for 'geometry' are 'rectangular' or 'circular'")
                R_throat_vec = np.append(R_throat_vec, R_throat)

                while True:
                    # inlet to throat
                    data = (k_in, 1, R_in / R_throat, 0.0075)
                    a, b, c, d = opt.fsolve(self.NozzleGeometry, (0.1, 0.1, 0.1, 0.1),
                                            args=data, full_output=False, xtol=1.0e-04)
                    x_in = np.linspace(0, k_in, int(self.thermo.samples / 2))
                    R_norm[0:int(self.thermo.samples / 2)] = a + b * np.tanh(c * x_in[::-1] - d)

                    # throat to outlet
                    data = (k_out, 1, R_out / R_throat, 0.0075)
                    a, b, c, d = opt.fsolve(self.NozzleGeometry, (0.1, 0.1, 0.1, 0.1),
                                            args=data, full_output=False, xtol=1.0e-04)
                    x_out = np.linspace(0, k_out, int(self.thermo.samples / 2) + 1)
                    R_norm[int(self.thermo.samples / 2):] = a + b * np.tanh(c * x_out[1:] - d)

                    # inlet to outlet
                    self.x_norm[ii, 0:int(self.thermo.samples / 2)] = x_in
                    self.x_norm[ii, int(self.thermo.samples / 2):] = x_out[1:] + x_in[-1]
                    spl = UnivariateSpline(self.x_norm[ii, :], R_norm, s=0)
                    self.R_norm[ii, :] = spl(self.x_norm[ii, :])
                    fig, ax = plt.subplots()
                    ax.plot(self.x_norm[ii, :], self.R_norm[ii, :], lw=2, color='black')
                    ax.plot(self.x_norm[ii, :], -self.R_norm[ii, :], lw=2, color='black')
                    ax.grid(1)
                    plt.show()
                    flag = input("\n Is the nozzle shape ok? (y/n) ")
                    if flag == "y" or flag == 'Y':
                        break
                    else:
                        k_in, k_out = [float(x) for x in input(" Stop and modify 'velocity' and/or 'mass flow rate' "
                                                               "or input new values of k_in, k_out: ").split()]
            elif switch == 'off-design':       # off-design mode: given nozzle geometry
                self.x_norm[ii, :] = x_norm[ii, :]
                self.R_norm[ii, :] = R_norm[ii, :]
            else:
                raise ValueError(" The available choices for 'switch' are 'design' or 'off-design'")

            # set nozzle area variation
            self.R_vec[ii, :] = self.R_norm[ii, :] * R_throat_vec[ii]
            if self.geometry == 'rectangular':
                self.A_vec[ii, :] = self.R_vec[ii, :] * 2
            elif self.geometry == 'circular':
                self.A_vec[ii, :] = np.pi * self.R_vec[ii, :] ** 2
            else:
                raise ValueError(" The available choices for 'geometry' are 'rectangular' or 'circular'")

            # forward propagation: compute P from inlet to outlet
            for jj in range(self.thermo.samples):
                if jj == 0:
                    self.P_vec[ii, jj] = P_in
                elif jj == int(self.thermo.samples / 2):
                    self.P_vec[ii, jj] = P_throat * 0.999  # impose supersonic flow solution
                else:
                    data = (ii, self.massflow[ii], self.A_vec[ii, jj])
                    self.P_vec[ii, jj] = opt.fsolve(self.NozzleEquations, (self.P_vec[ii, jj - 1]),
                                                    args=data, full_output=False, xtol=1.0e-08)

            # smoothing of P in the neighbourhood of the throat
            base_unit = int(self.thermo.samples / 10)
            throat_idx = int(self.thermo.samples / 2 - 1)
            k_in = self.x_norm[ii, throat_idx] - self.x_norm[ii, 0]
            k_out = self.x_norm[ii, -1] - self.x_norm[ii, throat_idx]
            inlet_idx = throat_idx - int(base_unit * k_out / k_in)
            outlet_idx = throat_idx + base_unit
            P_vec = np.append(self.P_vec[ii, 0:inlet_idx], self.P_vec[ii, throat_idx])
            P_vec = np.append(P_vec, self.P_vec[ii, outlet_idx:])
            x_tmp = np.append(self.x_norm[ii, 0:inlet_idx], self.x_norm[ii, throat_idx])
            x_tmp = np.append(x_tmp, self.x_norm[ii, outlet_idx:])
            spl = UnivariateSpline(x_tmp, P_vec, s=0)
            self.P_vec[ii, :] = spl(self.x_norm[ii, :])

            # compute thermodynamic properties along the isentropic expansion
            for jj in range(self.thermo.samples):
                self.thermo.EoS.update(CoolProp.PSmass_INPUTS, self.P_vec[ii, jj], self.thermo.s_in[ii])
                self.D_vec[ii, jj] = self.thermo.EoS.rhomass()
                self.V_vec[ii, jj] = np.sqrt(2 * (self.thermo.ht_in[ii] - self.thermo.EoS.hmass()))
                self.T_vec[ii, jj] = self.thermo.EoS.T()
                self.M_vec[ii, jj] = self.V_vec[ii, jj] / self.thermo.EoS.speed_sound()
                self.h_vec[ii, jj] = self.thermo.EoS.hmass()
                self.ht_vec[ii, jj] = self.h_vec[ii, jj] + self.V_vec[ii, jj] ** 2 / 2
                self.Z_vec[ii, jj] = self.P_vec[ii, jj] / (self.D_vec[ii, jj] * self.thermo.R * self.T_vec[ii, jj])
                self.Z_mean[ii] = np.mean(self.Z_vec[ii, :])
                self.gamma[ii, jj] = self.thermo.EoS.cpmass() / self.thermo.EoS.cvmass()
                dP_dv_T = (- 1 / (self.D_vec[ii, jj] ** 2) *
                           self.thermo.EoS.first_partial_deriv(CoolProp.iDmass, CoolProp.iP, CoolProp.iT)) ** (-1)
                self.gamma_Pv_vec[ii, jj] = - 1 / (self.P_vec[ii, jj] * self.D_vec[ii, jj]) * \
                                            self.thermo.EoS.cpmass() / self.thermo.EoS.cvmass() * dP_dv_T
                self.thermo.EoS.update(CoolProp.HmassSmass_INPUTS, self.ht_vec[ii, jj], self.thermo.s_in[ii])
                self.Pt_vec[ii, jj] = self.thermo.EoS.p()
            self.gamma_Pv_mean[ii] = np.log(self.P_vec[ii, 0] / self.P_vec[ii, -1]) / \
                                     np.log(self.D_vec[ii, 0] / self.D_vec[ii, -1])

        # plotting
        self.thermo.plotClass('output/' + self.thermo.fluid + '/expansion')
        self.thermo.plotClass.PlotNozzle(self.thermo.Pt_in, self.thermo.Dt_in, self.P_vec, self.D_vec,
                                         self.Z_vec, self.gamma_Pv_vec, self.gamma_Pv_mean, self.thermo.gamma_id,
                                         self.M_vec, self.A_vec, self.x_norm, self.R_norm)
        writeNozzleCoordinates(self.thermo.plotClass.results_dir, self.thermo.plotClass.labels,
                               self.x_norm, self.R_norm, R_throat_vec, self.massflow, self.velocity, self.geometry)
        writeNozzleFlow(self.thermo.plotClass.results_dir, self.thermo.plotClass.labels,
                        self.x_norm, self.P_vec, self.D_vec, self.M_vec)

        return

    def DiffuserCompression(self, beta_diff=1.2, eps=5, H_Rin=0.1, alpha_in=65):
        """
        Given the inlet and outlet conditions, compute the isentropic states along a compression process.
        Additionally, design a conical/radial diffuser for the given target pressure recovery and compute the
        isentropic flow along it.
        The first part of the compression process resembles the transformation occurring in an impeller,
        whereas the second corresponds to the one occurring in a diffuser.

        Optional arguments:
        :param beta_diff: pressure ratio of the diffuser (it must be > 1), used to define diffuser inlet state.
        Decrease the default value if the inlet Mach number is too large (it must be < 1)
        :param eps: semi-diffusion angle (only for conical diffuser), common values = 4-5 deg
        :param H_Rin: ratio among diffuser height and inlet radius (only for radial diffuser), common values = 0.05-0.4
        :param alpha_in: diffuser inlet absolute flow angle (only for radial diffuser), common values = 60-80 deg

        Notes:
        - self.geometry: 'conical' (for axial machines) or 'radial' (for radial machines)
        - self.massflow: used to compute the area and the radius variation of the diffuser (it must be >= 0).
        Use 1.0 if you're interested only in the normalized radius variation.
        - self.velocity: flow velocity at the outlet of the diffuser (it must be >= 0).
        Decrease if the outlet Mach number is too large (it must be < 1).
        """
        # check values of optional arguments
        if beta_diff <= 1:
            raise ValueError("The value assigned to 'beta_diff' must be larger than one")

        # iterate over the selected inlet states
        for ii in range(len(self.thermo.Tt_in)):
            if (self.massflow[ii] <= 0) or (eps <= 0) or (self.velocity[ii] <= 0):
                raise ValueError("The values assigned to 'mass flow rate', 'velocity' and 'eps' "
                                 "must be larger than zero")

            # define the inlet and outlet states of the diffuser
            self.thermo.EoS.update(CoolProp.PSmass_INPUTS, self.thermo.P_out[ii], self.thermo.s_in[ii])
            h_out = self.thermo.EoS.hmass()
            D_out = self.thermo.EoS.rhomass()
            ht_out = h_out + self.velocity[ii] ** 2 / 2
            ht_in = ht_out
            self.thermo.EoS.update(CoolProp.HmassSmass_INPUTS, ht_in, self.thermo.s_in[ii])
            self.thermo.Pt_in_diff[ii] = self.thermo.EoS.p()
            self.P_diff[ii, 0] = self.thermo.P_out[ii] / beta_diff
            self.Cp_diff[ii] = (self.thermo.P_out[ii] - self.P_diff[ii, 0]) / \
                               (self.thermo.Pt_in_diff[ii] - self.P_diff[ii, 0])
            self.thermo.EoS.update(CoolProp.PSmass_INPUTS, self.P_diff[ii, 0], self.thermo.s_in[ii])
            h_in = self.thermo.EoS.hmass()
            self.D_diff[ii, 0] = self.thermo.EoS.rhomass()
            self.T_diff[ii, 0] = self.thermo.EoS.T()
            self.V_diff[ii, 0] = np.sqrt(2 * (ht_in - h_in))
            self.M_diff[ii, 0] = self.V_diff[ii, 0] / self.thermo.EoS.speed_sound()
            if self.M_diff[ii, 0] >= 1:
                raise ValueError("Supersonic flow at the inlet of the diffuser; reduce 'beta_diff'")

            # define the radius distribution along the diffuser
            if self.geometry == 'conical':
                print("\n Computing isentropic compression(s) in conical diffuser(s) ...")
                A_in = self.massflow[ii] / (self.D_diff[ii, 0] * self.V_diff[ii, 0])
                A_out = self.massflow[ii] / (D_out * self.velocity[ii])
                self.R_diff[ii, 0] = np.sqrt(A_in / np.pi)
                self.R_diff[ii, -1] = np.sqrt(A_out / np.pi)
                self.L_diff[ii, :] = np.linspace(0, (self.R_diff[ii, -1] - self.R_diff[ii, 0]) /
                                                 np.tan(np.deg2rad(eps)), self.thermo.samples)
                self.R_diff[ii, :] = self.R_diff[ii, 0] + self.L_diff[ii, :] * np.tan(np.deg2rad(eps))
                self.A_diff[ii, :] = np.pi * self.R_diff[ii, :] ** 2

                # compute the isentropic states along the diffuser
                for jj in range(self.thermo.samples - 1):
                    self.V_diff[ii, jj + 1] = opt.fsolve(self.ConicalDiffuserEquations, (self.V_diff[ii, jj]),
                                                         args=(self.massflow[ii], self.A_diff[ii, jj + 1], ht_in,
                                                               self.thermo.s_in[ii]),
                                                         full_output=False, xtol=1.0e-06)
                    self.D_diff[ii, jj + 1] = self.massflow[ii] / (self.A_diff[ii, jj + 1] * self.V_diff[ii, jj + 1])
                    self.thermo.EoS.update(CoolProp.DmassSmass_INPUTS, self.D_diff[ii, jj + 1], self.thermo.s_in[ii])
                    self.P_diff[ii, jj + 1] = self.thermo.EoS.p()
                    self.T_diff[ii, jj + 1] = self.thermo.EoS.T()
                    self.M_diff[ii, jj + 1] = self.V_diff[ii, jj + 1] / self.thermo.EoS.speed_sound()
            elif self.geometry == 'radial':
                print("\n Computing isentropic compression(s) in radial diffuser(s) ...")
                self.alpha_diff[ii, 0] = np.deg2rad(alpha_in)
                Vm = self.V_diff[ii, 0] * np.cos(self.alpha_diff[ii, 0])
                Vt = self.V_diff[ii, 0] * np.sin(self.alpha_diff[ii, 0])
                A_in = self.massflow[ii] / (self.D_diff[ii, 0] * Vm)
                self.R_diff[ii, 0] = np.sqrt(A_in / (2 * np.pi * H_Rin))
                H = H_Rin * self.R_diff[ii, 0]
                mass = self.D_diff[ii, 0] * Vm * self.R_diff[ii, 0] * H
                momentum = self.R_diff[ii, 0] * Vt
                energy = ht_in
                self.R_diff[ii, -1] = np.sqrt((mass ** 2 + (D_out * H * momentum) ** 2) /
                                              ((D_out * H * self.velocity[ii]) ** 2))
                self.R_diff[ii, :] = np.linspace(self.R_diff[ii, 0], self.R_diff[ii, -1], self.thermo.samples)
                self.A_diff[ii, :] = 2 * np.pi * self.R_diff[ii, :] * H

                # compute the isentropic states along the diffuser
                for jj in range(self.thermo.samples - 1):
                    self.P_diff[ii, jj + 1] = opt.fsolve(self.RadialDiffuserEquations, (self.P_diff[ii, jj]),
                                                         args=(mass, momentum, energy, self.R_diff[ii, jj + 1],
                                                               self.thermo.s_in[ii], H), 
                                                         full_output=False, xtol=1.0e-06)
                    self.thermo.EoS.update(CoolProp.PSmass_INPUTS, self.P_diff[ii, jj + 1], self.thermo.s_in[ii])
                    h = self.thermo.EoS.hmass()
                    self.D_diff[ii, jj + 1] = self.thermo.EoS.rhomass()
                    self.T_diff[ii, jj + 1] = self.thermo.EoS.T()
                    Vt = momentum / self.R_diff[ii, jj + 1]
                    self.V_diff[ii, jj + 1] = np.sqrt(2 * (energy - h))
                    Vm = np.sqrt(self.V_diff[ii, jj + 1] ** 2 - Vt ** 2)
                    self.alpha_diff[ii, jj + 1] = np.arctan(Vt / Vm)
                    self.M_diff[ii, jj + 1] = self.V_diff[ii, jj + 1] / self.thermo.EoS.speed_sound()

            # impose linear pressure variation among inlet and outlet state
            self.P_vec[ii, :] = np.linspace(self.thermo.Pt_in[ii], self.thermo.P_out[ii], self.thermo.samples)
            for jj in range(self.thermo.samples):
                self.thermo.EoS.update(CoolProp.PSmass_INPUTS, self.P_vec[ii, jj], self.thermo.s_in[ii])
                self.D_vec[ii, jj] = self.thermo.EoS.rhomass()
                self.T_vec[ii, jj] = self.thermo.EoS.T()
                self.Z_vec[ii, jj] = self.P_vec[ii, jj] / (self.D_vec[ii, jj] * self.thermo.R * self.T_vec[ii, jj])
                self.Z_mean[ii] = np.mean(self.Z_vec[ii, :])
                dP_dv_T = (- 1 / (self.D_vec[ii, jj] ** 2) *
                           self.thermo.EoS.first_partial_deriv(CoolProp.iDmass, CoolProp.iP, CoolProp.iT)) ** (-1)
                self.gamma_Pv_vec[ii, jj] = - 1 / (self.P_vec[ii, jj] * self.D_vec[ii, jj]) * \
                                            self.thermo.EoS.cpmass() / self.thermo.EoS.cvmass() * dP_dv_T
            self.gamma_Pv_mean[ii] = np.log(self.P_vec[ii, 0] / self.P_vec[ii, -1]) / \
                                     np.log(self.D_vec[ii, 0] / self.D_vec[ii, -1])

        # plotting
        self.thermo.plotClass('output/' + self.thermo.fluid + '/compression')
        if self.geometry == 'conical':
            self.thermo.plotClass.PlotConicalDiffuser(beta_diff, self.Z_vec, self.gamma_Pv_vec, self.gamma_Pv_mean,
                                                      self.thermo.gamma_id, self.P_vec, self.P_diff, self.M_diff,
                                                      self.R_diff)
        elif self.geometry == 'radial':
            self.thermo.plotClass.PlotRadialDiffuser(beta_diff, self.Z_vec, self.gamma_Pv_vec, self.gamma_Pv_mean,
                                                     self.thermo.gamma_id, self.P_vec, self.P_diff, self.M_diff,
                                                     self.alpha_diff, self.R_diff)
        return

    def FindSonicState(self, p, *data):
        """
        Non-linear system of eqs. defining the pressure at which sonic state is established throughout an expansion.
        Conservation of energy.
        """
        loc = data
        P = p

        try:
            self.thermo.EoS.update(CoolProp.PSmass_INPUTS, P, self.thermo.s_in[loc])
            V = np.sqrt(2 * (self.thermo.ht_in[loc] - self.thermo.EoS.hmass()))
            M = V / self.thermo.EoS.speed_sound()
            res = M - 1.0

        except ValueError:
            res = 10

        return res

    def NozzleGeometry(self, p, *data):
        """
        Non-linear system of eqs. determining the nozzle profile from inlet to throat or from throat to outlet.
        Fit tanh function with null derivative at both boundaries.
        """
        x_out, R_in, R_out, toll = data
        a, b, c, d = p

        # R(0) = R_throat
        dx1 = (a + b * np.tanh(-d)) - R_in
        # R(x_out) = R_out
        dx2 = (a + b * np.tanh(c * x_out - d)) - R_out
        # dR/dx|0 --> 0
        dx3 = b * c / (np.cosh(d) ** 2) - toll
        # dR/dx|x_out --> 0
        dx4 = b * c / (np.cosh(d - c * x_out) ** 2) - toll

        return dx1, dx2, dx3, dx4

    def NozzleEquations(self, p, *data):
        """
        Non-linear system of eqs. defining the isentropic compressible flow in a nozzle, given its geometry.
        Conservation of mass + conservation of energy.
        """
        loc, massFlow, A = data
        P = p

        try:
            self.thermo.EoS.update(CoolProp.PSmass_INPUTS, P, self.thermo.s_in[loc])
            D = self.thermo.EoS.rhomass()
            V = np.sqrt(2 * (self.thermo.ht_in[loc] - self.thermo.EoS.hmass()))
            A_new = massFlow / (D * V)
            res = (A_new - A) / A

        except ValueError:
            res = 10

        return res

    def ConicalDiffuserEquations(self, p, *data):
        """
        Non-linear system of eqs. defining the isentropic compressible flow in a conical diffuser.
        Conservation of mass + conservation of energy.
        """
        massFlow, A, ht, s = data
        V = p

        try:
            D = massFlow / (A * V)
            self.thermo.EoS.update(CoolProp.DmassSmass_INPUTS, D, s)
            h1 = self.thermo.EoS.hmass()
            h2 = ht - V ** 2 / 2
            res = 2 * (h1 - h2) / (h1 + h2)

        except ValueError:
            res = 10

        return res

    def RadialDiffuserEquations(self, p, *data):
        """
        Non-linear system of eqs. defining the isentropic compressible flow in a radial diffuser.
        Conservation of mass, tangential momentum and energy.
        """
        mass, momentum, energy, R, s, H = data
        P = p

        try:
            self.thermo.EoS.update(CoolProp.PSmass_INPUTS, P, s)
            h = self.thermo.EoS.hmass()
            D = self.thermo.EoS.rhomass()
            Vt = momentum / R
            V = np.sqrt(2 * (energy - h))
            Vm = np.sqrt(V ** 2 - Vt ** 2)
            mass_new = R * D * Vm * H
            res = (mass_new - mass) / mass

        except ValueError:
            res = 10

        return res

