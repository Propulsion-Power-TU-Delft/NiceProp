#! /usr/bin/env python
# -*- coding: utf-8 -*-

#############################################################################
# NiceProp - Interactively learning NICFD
# Authors: ir. A. Giuffre', Dr. ir. M. Pini
# Content: Results rendering and plotting
# 2021 - TU Delft - All rights reserved
#############################################################################


import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os


class Plot:
    """ Class to plot results computed with the methods of classes ThermodynamicModel and IsentropicFlowModel """
    def __init__(self, x, y, xc, yc, x_sat, y_sat, x_in, x_out, y_in, y_out, x_iso, y_iso, labels,
                 plot_process, sc_iso=np.array([]), Tc_iso=np.array([]), cmap='viridis'):
        """
        Arguments:
        :param x: 2D array defining the x-axis of the grid used for contour plots
        :param y: 2D array defining the y-axis of the grid used for contour plots
        :param xc: x-coordinate of the critical state
        :param yc: y-coordinate of the critical state
        :param x_sat: 1D array defining the x-coordinates of the saturation curve
        :param y_sat: 1D array defining the y-coordinates of the saturation curve
        :param x_in: x-coordinate of the inlet state
        :param x_out: x-coordinate of the outlet state
        :param y_in: y-coordinate of the inlet state
        :param y_out: y-coordinate of the outlet state
        :param x_iso: 1D array defining the x-coordinates of the iso curve
        :param y_iso: 1D array defining the y-coordinates of the iso curve
        :param labels: labels associated to the thermodynamic transformations
               If there are more than one inlet states, then x_sat, y_sat, x_iso, y_iso are 2D arrays;
               x_in, x_out, y_in, y_out are 1D arrays.
        :param plot_process: flag to activate/deactivate plot of isentropic process(es) over the contour plots

        Optional arguments:
        :param sc_iso: 1D array defining the x-coordinate of the critical line in the T-s plane
        :param Tc_iso: 1D array defining the y-coordinate of the critical line in the T-s plane
        :param cmap: colormap used for all the plots

        Methods:
        - FindNearest
        - PlotTsContour
        - PlotPTContour
        - PlotExpansion
        - PlotCompression
        - PlotNozzle
        - PlotConicalDiffuser
        - PlotRadialDiffuser
        """
        self.x = x
        self.y = y
        self.xc = xc
        self.yc = yc
        self.x_sat = x_sat
        self.y_sat = y_sat
        self.x_in = x_in
        self.x_out = x_out
        self.y_in = y_in
        self.y_out = y_out
        self.x_iso = x_iso
        self.y_iso = y_iso
        self.sc_iso = sc_iso
        self.Tc_iso = Tc_iso
        self.labels = labels
        self.cmap = cmap
        self.plot_process = plot_process

        config = {  # setup matplotlib to use latex for output
            "pgf.texsystem": "pdflatex",
            "text.usetex": True,
            "font.family": "DejaVu Sans",
            "axes.titlesize": 30,
            "axes.labelsize": 30,
            "font.size": 30,
            "legend.fontsize": 20,
            "legend.frameon": True,
            "xtick.labelsize": 20,
            "ytick.labelsize": 20,
            "xtick.major.pad": 10,
            "ytick.major.pad": 10,
            "figure.autolayout": True,
            "figure.figsize": (9.6, 7.2)}

        mpl.rcParams.update(config)

    def __call__(self, folder):
        home_dir, _ = os.path.split(os.path.dirname(__file__))
        self.results_dir = os.path.join(home_dir, folder)
        self.jpeg_dir = os.path.join(home_dir, folder, 'jpeg')
        self.tiff_dir = os.path.join(home_dir, folder, 'tiff')
        if not os.path.isdir(self.results_dir):
            os.makedirs(self.results_dir)
        if not os.path.isdir(self.jpeg_dir):
            os.makedirs(self.jpeg_dir)
        if not os.path.isdir(self.tiff_dir):
            os.makedirs(self.tiff_dir)

    def FindNearest(self, x, vec):
        """
        Given a 1D array and a scalar, it returns the index of the element of the array that is closest to x
        :param x: scalar value
        :param vec: 1D array
        :return idx: index of the element of the array that is closest to the prescribed scalar value
        """
        vec = np.asarray(vec)
        idx = (np.abs(vec - x)).argmin()

        return idx

    def PlotTsContour(self, z, title, cbar_label, contour_bounds=np.array([0, 0]),
                      levels=50, powerNorm=False, powerNormCoeff=0.4):
        """ Create contour plot of the thermodynamic quantity z in the reduced T-s plane """
        fig, ax = plt.subplots()

        if self.plot_process:
            # plot isobaric lines corresponding to the selected inlet states
            for ii in range(len(self.x_in)):
                ax.plot(self.x_iso / self.xc, self.y_iso[ii, :] / self.yc, color='xkcd:light grey')

            # plot isentropic transformations(s)
            markers = ['s', 'o', '^', 'D', 'P', 'X', '*']
            for ii in range(len(self.x_in)):
                ax.plot(self.x_in[ii] / self.xc, self.y_in[ii] / self.yc,
                        marker=markers[ii], ms=10, color='black', label=self.labels[ii])
                ax.plot([self.x_in[ii] / self.xc, self.x_out[ii] / self.xc],
                        [self.y_in[ii] / self.yc, self.y_out[ii] / self.yc], color='black')

        # plot saturation curve and divide the thermodynamic plane into liquid, vapor and supercritical regions
        ax.plot(self.x_sat / self.xc, self.y_sat / self.yc, color='black')
        ax.plot(self.sc_iso / self.xc, self.Tc_iso / self.yc, linestyle='dashed', color='red')

        # create contour plot of the selected thermodynamic quantity z
        if powerNorm:
            if contour_bounds[0] == 0 and contour_bounds[1] == 0:
                CS = ax.contourf(self.x / self.xc, self.y / self.yc, z, levels,
                                 cmap=self.cmap, origin='upper', norm=mpl.colors.PowerNorm(gamma=powerNormCoeff))
            else:
                CS = ax.contourf(self.x / self.xc, self.y / self.yc, z,
                                 np.linspace(contour_bounds[0], contour_bounds[1], levels),
                                 cmap=self.cmap, origin='upper', norm=mpl.colors.PowerNorm(gamma=powerNormCoeff))
        else:
            if contour_bounds[0] == 0 and contour_bounds[1] == 0:
                CS = ax.contourf(self.x / self.xc, self.y / self.yc, z, levels, cmap=self.cmap, origin='upper')
            else:
                CS = ax.contourf(self.x / self.xc, self.y / self.yc, z,
                                 np.linspace(contour_bounds[0], contour_bounds[1], levels),
                                 cmap=self.cmap, origin='upper')

        # plot settings
        cbar = fig.colorbar(CS, shrink=1.0, format='%.2f')
        cbar.locator = mpl.ticker.MaxNLocator(nbins=8)
        cbar.update_ticks()
        cbar.ax.set_xlabel(cbar_label)
        cbar.ax.xaxis.set_label_coords(1.75, -0.04)
        ax.set_xlabel('$s_\mathrm{r}$ [-]')
        ax.set_ylabel('$T_\mathrm{r}$ [-]')
        ax.set_xlim(min(self.x_iso) / self.xc, max(self.x_iso) / self.xc)
        ax.set_ylim(min(self.y) / self.yc, max(self.y) / self.yc)

        if self.plot_process:
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles, labels, loc='upper left')

        fig.savefig(self.jpeg_dir + '/' + title + '.jpeg')
        fig.savefig(self.tiff_dir + '/' + title + '.tiff')
        plt.close(fig)

        return

    def PlotPTContour(self, z, title, cbar_label, levels=50, powerNorm=False, powerNormCoeff=0.4):
        """ Create contour plot of the thermodynamic quantity z in the reduced P-T plane """
        fig, ax = plt.subplots()

        # plot saturation curve and divide the thermodynamic plane into liquid, vapor and supercritical regions
        ax.plot(self.x_sat / self.xc, self.y_sat / self.yc, color='black')
        ax.hlines(1.0, 1.0, np.max(self.x) / self.xc, linestyle='dashed', color='red')

        if self.plot_process:
            # plot isentropic transformations(s)
            markers = ['s', 'o', '^', 'D', 'P', 'X', '*']
            for ii in range(len(self.x_in)):
                ax.plot(self.x_in[ii] / self.xc, self.y_in[ii] / self.yc,
                        marker=markers[ii], ms=10, color='black', label=self.labels[ii])
                ax.plot(self.x_iso[ii, :] / self.xc, self.y_iso[ii, :] / self.yc, color='black')

        # create contour plot of the selected thermodynamic quantity z
        if powerNorm:
            CS = ax.contourf(self.x / self.xc, self.y / self.yc, z, levels,
                             cmap=self.cmap, origin='upper', norm=mpl.colors.PowerNorm(gamma=powerNormCoeff))
        else:
            CS = ax.contourf(self.x / self.xc, self.y / self.yc, z, levels, cmap=self.cmap, origin='upper')

        # plot settings
        cbar = fig.colorbar(CS, shrink=1.0, format='%.2f')
        cbar.locator = mpl.ticker.MaxNLocator(nbins=8)
        cbar.update_ticks()
        cbar.ax.set_xlabel(cbar_label)
        cbar.ax.xaxis.set_label_coords(1.75, -0.04)
        ax.set_xlabel('$T_\mathrm{r}$ [-]')
        ax.set_ylabel('$P_\mathrm{r}$ [-]')
        ax.set_xlim(left=np.min(self.x_sat) / self.xc)
        ax.set_ylim(bottom=np.min(self.y_sat) / self.yc)

        if self.plot_process:
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles, labels, loc='upper left')

        fig.savefig(self.jpeg_dir + '/' + title + '.jpeg')
        fig.savefig(self.tiff_dir + '/' + title + '.tiff')
        plt.close(fig)

        return

    def PlotExpansion(self, Pt_in, Dt_in, P_vec, D_vec, Z_vec, gamma_Pv_vec, M_vec):
        """ Plot trend of gamma Pv, Mach, Z along the isentropic expansion """

        # iterate over inlet states
        for ii in range(len(self.labels)):
            fig1, ax1 = plt.subplots()  # gamma_Pv, Z trend vs beta
            fig2, ax2 = plt.subplots()  # gamma_Pv, Z trend vs alpha
            ax1t1 = ax1.twinx()
            ax1t2 = ax1.twinx()
            ax2t1 = ax2.twinx()
            ax2t2 = ax2.twinx()

            ax1.plot(Pt_in[ii] / P_vec[ii, :], M_vec[ii, :], lw=2, color='xkcd:green', label='$M$ [-]')
            ax1t1.plot(Pt_in[ii] / P_vec[ii, :], Z_vec[ii, :], lw=2, color='xkcd:black',
                       label='$Z$', linestyle='dashed')
            ax1t2.plot(Pt_in[ii] / P_vec[ii, :], gamma_Pv_vec[ii, :], lw=2, color='xkcd:grey',
                       label='$\gamma_{Pv}$', linestyle='dashed')

            ax2.plot(Dt_in[ii] / D_vec[ii, :], M_vec[ii, :], lw=2, color='xkcd:green',
                     label='$M$ [-]')
            ax2t1.plot(Dt_in[ii] / D_vec[ii, :], Z_vec[ii, :], lw=2, color='xkcd:black',
                       label='$Z$', linestyle='dashed')
            ax2t2.plot(Dt_in[ii] / D_vec[ii, :], gamma_Pv_vec[ii, :], lw=2, color='xkcd:grey',
                       label='$\gamma_{Pv}$', linestyle='dashed')

            ax1.grid(1)
            ax2.grid(1)

            ax1.set_xlabel(r'$\beta_\mathrm{ts}$ [-]')
            ax1.set_ylabel('$M$ [-]')
            ax1t1.set_ylabel('$Z$ [-]')
            ax1t2.set_ylabel(r'$\gamma_{Pv}$ [-]')
            ax1t2.spines['right'].set_position(('outward', 75))
            ax1.yaxis.label.set_color(color='xkcd:green')
            ax1t2.yaxis.label.set_color(color='xkcd:grey')

            ax2.set_xlabel(r'$\alpha_\mathrm{ts}$ [-]')
            ax2.set_ylabel('$M$ [-]')
            ax2t1.set_ylabel('$Z$ [-]')
            ax2t2.set_ylabel(r'$\gamma_{Pv}$ [-]')
            ax2t2.spines['right'].set_position(('outward', 75))
            ax2.yaxis.label.set_color(color='xkcd:green')
            ax2t2.yaxis.label.set_color(color='xkcd:grey')

            fig1.savefig(self.jpeg_dir + '/NICFD_vs_beta_' + self.labels[ii] + '.jpeg')
            fig2.savefig(self.jpeg_dir + '/NICFD_vs_alpha_' + self.labels[ii] + '.jpeg')

            fig1.savefig(self.tiff_dir + '/NICFD_vs_beta_' + self.labels[ii] + '.tiff')
            fig2.savefig(self.tiff_dir + '/NICFD_vs_alpha_' + self.labels[ii] + '.tiff')

            plt.close(fig1)
            plt.close(fig2)

        return

    def PlotCompression(self, Pt_in, Dt_in, P_vec, D_vec, Z_vec, gamma_Pv_vec):
        """ Plot trend of gamma Pv, Mach, Z along the isentropic expansion """

        # iterate over inlet states
        for ii in range(len(self.labels)):
            fig1, ax1 = plt.subplots()  # gamma_Pv, Z trend vs beta
            fig2, ax2 = plt.subplots()  # gamma_Pv, Z trend vs alpha
            ax1t1 = ax1.twinx()
            ax2t1 = ax2.twinx()

            ax1.plot(P_vec[ii, :] / Pt_in[ii], Z_vec[ii, :], lw=2, color='xkcd:green',
                     label='$Z$ [-]')
            ax1t1.plot(P_vec[ii, :] / Pt_in[ii], gamma_Pv_vec[ii, :], lw=2, color='xkcd:grey',
                       label='$\gamma_{Pv}$', linestyle='dashed')

            ax2.plot(D_vec[ii, :] / Dt_in[ii], Z_vec[ii, :], lw=2, color='xkcd:green',
                     label='$Z$ [-]', linestyle='dashed')
            ax2t1.plot(D_vec[ii, :] / Dt_in[ii], gamma_Pv_vec[ii, :], lw=2, color='xkcd:grey',
                       label='$\gamma_{Pv}$', linestyle='dashed')

            ax1.grid(1)
            ax2.grid(1)

            ax1.set_xlabel(r'$\beta_\mathrm{ts}$ [-]')
            ax1.set_ylabel('$Z$ [-]')
            ax1t1.set_ylabel(r'$\gamma_{Pv}$ [-]')
            ax1.yaxis.label.set_color(color='xkcd:green')
            ax1t1.yaxis.label.set_color(color='xkcd:grey')

            ax2.set_xlabel(r'$\alpha_\mathrm{ts}$ [-]')
            ax2.set_ylabel('$Z$ [-]')
            ax2t1.set_ylabel(r'$\gamma_{Pv}$ [-]')
            ax2.yaxis.label.set_color(color='xkcd:green')
            ax2t1.yaxis.label.set_color(color='xkcd:grey')

            fig1.savefig(self.jpeg_dir + '/NICFD_vs_beta_' + self.labels[ii] + '.jpeg')
            fig2.savefig(self.jpeg_dir + '/NICFD_vs_alpha_' + self.labels[ii] + '.jpeg')

            fig1.savefig(self.tiff_dir + '/NICFD_vs_beta_' + self.labels[ii] + '.tiff')
            fig2.savefig(self.tiff_dir + '/NICFD_vs_alpha_' + self.labels[ii] + '.tiff')

            plt.close(fig1)
            plt.close(fig2)

        return

    def PlotNozzle(self, Pt_in, Dt_in, P_vec, D_vec, Z_vec, gamma_Pv_vec, gamma_Pv_mean, gamma_id, M_vec, A_vec,
                   x_norm, R_norm):
        """ Plot nozzle shape plus trends of gamma Pv, Mach, Area ratio along the isentropic expansion """
        new_colors = [plt.get_cmap(self.cmap)(1. * i / len(self.labels)) for i in range(len(self.labels))]

        fig1, ax1 = plt.subplots()      # gamma_Pv trend
        fig2, ax2 = plt.subplots()      # Mach trend
        fig3, ax3 = plt.subplots()      # non-dimensional nozzle shape
        fig4, ax4 = plt.subplots()      # normalized A trend
        fig5, ax5 = plt.subplots()      # Z trend
        fig6, ax6 = plt.subplots()      # beta and alpha trend

        # iterate over inlet states
        for ii in range(len(self.labels)):

            ax1.plot(x_norm[ii, :], gamma_Pv_vec[ii, :], lw=2, color=new_colors[ii], label=self.labels[ii])
            ax1.plot(x_norm[ii, :], gamma_Pv_mean[ii] * np.ones(len(x_norm[ii, :])),
                     linestyle='dashed', color=new_colors[ii], label=r'$\overline{\gamma}_{Pv}$')
            ax2.plot(x_norm[ii, :], M_vec[ii, :], lw=2, color=new_colors[ii], label=self.labels[ii])

            ax3.plot(x_norm[ii, :], R_norm[ii, :], lw=2, color=new_colors[ii], label=self.labels[ii])
            ax3.plot(x_norm[ii, :], -R_norm[ii, :], lw=2, color=new_colors[ii])

            ax4.plot(x_norm[ii, :], A_vec[ii, :] / A_vec[ii, int(len(A_vec[ii, :]) / 2 - 1)],
                     lw=2, color=new_colors[ii], label=self.labels[ii])
            ax5.plot(x_norm[ii, :], Z_vec[ii, :], lw=2, color=new_colors[ii], label=self.labels[ii])

            ax6.plot(x_norm[ii, :], Pt_in[ii] / P_vec[ii, :], lw=2, color=new_colors[ii],
                     label=self.labels[ii] + ': ' + r'$\beta_\mathrm{ts}$')
            ax6.plot(x_norm[ii, :], Dt_in[ii] / D_vec[ii, :], lw=2, color=new_colors[ii],
                     label=self.labels[ii] + ': ' + r'$\alpha_\mathrm{ts}$', linestyle='dashed')

        ax1.plot(x_norm[0, :], gamma_id * np.ones(len(x_norm[0, :])),
                 linestyle='dashed', color='xkcd:grey', label=r'$\gamma$')
        ax2.hlines(1.0, 0.0, x_norm[0, -1], linestyle='dashed', color='xkcd:grey', label='$M=1$')
        ax3.vlines(x_norm[0, np.argmin(R_norm[0, :])], -1.0, 1.0, linestyle='dashed', color='xkcd:grey', label='throat')
        ax5.plot(x_norm[0, :], np.ones(len(x_norm[0, :])), linestyle='dashed', color='xkcd:grey', label='$Z=1$')

        ax1.grid(1)
        ax2.grid(1)
        ax3.grid(1)
        ax4.grid(1)
        ax5.grid(1)
        ax6.grid(1)
        
        ax1.set_xlabel(r'$x_\mathrm{norm}$ [-]')
        ax2.set_xlabel(r'$x_\mathrm{norm}$ [-]')
        ax3.set_xlabel(r'$x_\mathrm{norm}$ [-]')
        ax4.set_xlabel(r'$x_\mathrm{norm}$ [-]')
        ax5.set_xlabel(r'$x_\mathrm{norm}$ [-]')
        ax6.set_xlabel(r'$x_\mathrm{norm}$ [-]')
        ax1.set_ylabel(r'$\gamma_{Pv}$ [-]')
        ax2.set_ylabel('$M$ [-]')
        ax3.set_ylabel(r'$y_\mathrm{norm}$ [-]')
        ax4.set_ylabel(r'$A_\mathrm{norm}$ [-]')
        ax5.set_ylabel('$Z$ [-]')

        handles1, labels1 = ax1.get_legend_handles_labels()
        handles2, labels2 = ax2.get_legend_handles_labels()
        handles3, labels3 = ax3.get_legend_handles_labels()
        handles4, labels4 = ax4.get_legend_handles_labels()
        handles5, labels5 = ax5.get_legend_handles_labels()
        handles6, labels6 = ax6.get_legend_handles_labels()

        ax1.legend(handles1, labels1)
        ax2.legend(handles2, labels2)
        ax3.legend(handles3, labels3, loc='center right')
        ax4.legend(handles4, labels4)
        ax5.legend(handles5, labels5)
        ax6.legend(handles6, labels6)

        fig1.savefig(self.jpeg_dir + '/gamma_Pv.jpeg')
        fig2.savefig(self.jpeg_dir + '/Mach.jpeg')
        fig3.savefig(self.jpeg_dir + '/nozzle.jpeg')
        fig4.savefig(self.jpeg_dir + '/A_norm.jpeg')
        fig5.savefig(self.jpeg_dir + '/Z.jpeg')
        fig6.savefig(self.jpeg_dir + '/beta_alpha.jpeg')
        
        fig1.savefig(self.tiff_dir + '/gamma_Pv.tiff')
        fig2.savefig(self.tiff_dir + '/Mach.tiff')
        fig3.savefig(self.tiff_dir + '/nozzle.tiff')
        fig4.savefig(self.tiff_dir + '/A_norm.tiff')
        fig5.savefig(self.tiff_dir + '/Z.tiff')
        fig6.savefig(self.tiff_dir + '/beta_alpha.tiff')

        plt.close(fig1)
        plt.close(fig2)
        plt.close(fig3)
        plt.close(fig4)
        plt.close(fig5)
        plt.close(fig6)

        return

    def PlotConicalDiffuser(self, beta_diff, Z_vec, gamma_Pv_vec, gamma_Pv_mean, gamma_id,
                            P_vec, P_diff, M_diff, R_diff):
        """ Plot diffuser shape plus trends of gamma Pv, Mach, compression ratio along the isentropic compression """
        new_colors = [plt.get_cmap(self.cmap)(1. * i / len(self.labels)) for i in range(len(self.labels))]
        beta_in_diff = np.zeros(len(self.labels))
        gamma_Pv_in_diff = np.zeros(len(self.labels))

        fig1, ax1 = plt.subplots()      # gamma_Pv trend
        fig2, ax2 = plt.subplots()      # Z trend
        fig3, ax3 = plt.subplots()      # Mach trend (diffuser)
        fig4, ax4 = plt.subplots()      # diffuser shape
        fig5, ax5 = plt.subplots()      # normalized P trend (diffuser)
        # iterate over inlet states
        for ii in range(len(self.labels)):
            R_norm = R_diff[ii, :] / R_diff[ii, 0]
            unit_vec = np.linspace(0, 1, len(R_norm))
            beta = P_vec[ii, :] / P_vec[ii, 0]
            beta_in_diff[ii] = beta[-1] / beta_diff
            gamma_Pv_in_diff[ii] = gamma_Pv_vec[ii, self.FindNearest(beta_in_diff[ii], beta)]
            ax1.plot(beta, gamma_Pv_vec[ii, :], lw=2, color=new_colors[ii], label=self.labels[ii])
            ax1.plot(beta, gamma_Pv_mean[ii] * np.ones(len(beta)),
                     linestyle='dashed', color=new_colors[ii], label=r'$\overline{\gamma}_{Pv}$')
            ax2.plot(beta, Z_vec[ii, :], lw=2, color=new_colors[ii], label=self.labels[ii])
            ax3.plot(unit_vec, M_diff[ii, :], lw=2, color=new_colors[ii], label=self.labels[ii])
            ax4.plot(unit_vec, R_norm, lw=2, color=new_colors[ii], label=self.labels[ii])
            ax4.plot(unit_vec, -R_norm, lw=2, color=new_colors[ii])
            ax5.plot(unit_vec, P_diff[ii, :] / P_diff[ii, 0], lw=2, color=new_colors[ii], label=self.labels[ii])
        ax1.plot(beta, gamma_id * np.ones(len(beta)), linestyle='dashed', color='xkcd:grey', label=r'$\gamma$')
        ax1.scatter(beta_in_diff, gamma_Pv_in_diff, marker='o', s=50, color='black', label='diffuser inlet')
        ax2.plot(beta, np.ones(len(beta)), linestyle='dashed', color='xkcd:grey', label='$Z=1$')

        ax1.grid(1)
        ax2.grid(1)
        ax3.grid(1)
        ax4.grid(1)
        ax5.grid(1)
        
        ax1.set_xlabel(r'$\beta$ [-]')
        ax2.set_xlabel(r'$\beta$ [-]')
        ax3.set_xlabel(r'$x_\mathrm{norm}$ [-]')
        ax4.set_xlabel(r'$x_\mathrm{norm}$ [-]')
        ax5.set_xlabel(r'$x_\mathrm{norm}$ [-]')
        ax1.set_ylabel(r'$\gamma_{Pv}$ [-]')
        ax2.set_ylabel('$Z$ [-]')
        ax3.set_ylabel('$M$ [-]')
        ax4.set_ylabel(r'$R_\mathrm{norm}$ [-]')
        ax5.set_ylabel(r'$\beta$ [-]')
        
        handles1, labels1 = ax1.get_legend_handles_labels()
        handles2, labels2 = ax2.get_legend_handles_labels()
        handles3, labels3 = ax3.get_legend_handles_labels()
        handles4, labels4 = ax4.get_legend_handles_labels()
        handles5, labels5 = ax5.get_legend_handles_labels()
        
        ax1.legend(handles1, labels1)
        ax2.legend(handles2, labels2)
        ax3.legend(handles3, labels3)
        ax4.legend(handles4, labels4, loc='center right')
        ax5.legend(handles5, labels5)
        
        fig1.savefig(self.jpeg_dir + '/gamma_Pv.jpeg')
        fig2.savefig(self.jpeg_dir + '/Z.jpeg')
        fig3.savefig(self.jpeg_dir + '/Mach_diff_conical.jpeg')
        fig4.savefig(self.jpeg_dir + '/conical_diffuser.jpeg')
        fig5.savefig(self.jpeg_dir + '/beta_diff_conical.jpeg')
        
        fig1.savefig(self.tiff_dir + '/gamma_Pv.tiff')
        fig2.savefig(self.tiff_dir + '/Z.tiff')
        fig3.savefig(self.tiff_dir + '/Mach_diff_conical.tiff')
        fig4.savefig(self.tiff_dir + '/conical_diffuser.tiff')
        fig5.savefig(self.tiff_dir + '/beta_diff_conical.tiff')
        
        plt.close(fig1)
        plt.close(fig2)
        plt.close(fig3)
        plt.close(fig4)
        plt.close(fig5)

        return

    def PlotRadialDiffuser(self, beta_diff, Z_vec, gamma_Pv_vec, gamma_Pv_mean, gamma_id,
                           P_vec, P_diff, M_diff, alpha_diff, R_diff):
        """
        Plot diffuser shape plus trends of gamma Pv, Mach, compression ratio, flow angle along the isentropic compression
        """
        new_colors = [plt.get_cmap(self.cmap)(1. * i / len(self.labels)) for i in range(len(self.labels))]
        beta_in_diff = np.zeros(len(self.labels))
        gamma_Pv_in_diff = np.zeros(len(self.labels))
        x_vec = np.zeros((R_diff.shape[0], R_diff.shape[1], 4))
        y_vec = np.zeros((R_diff.shape[0], R_diff.shape[1], 4))
        theta = np.linspace(0, 2 * np.pi, 1000)

        fig1, ax1 = plt.subplots()      # gamma_Pv trend
        fig2, ax2 = plt.subplots()      # Z trend
        fig3, ax3 = plt.subplots()      # Mach trend (diffuser)
        fig4, ax4 = plt.subplots()      # diffuser shape
        fig5, ax5 = plt.subplots()      # normalized P trend (diffuser)
        fig6, ax6 = plt.subplots()      # flow angle trend (diffuser)
        # iterate over inlet states
        for ii in range(len(self.labels)):
            R_norm = R_diff[ii, :] / R_diff[ii, 0]
            x_in = R_norm[0] * np.cos(theta)
            y_in = R_norm[0] * np.sin(theta)
            x_out = R_norm[-1] * np.cos(theta)
            y_out = R_norm[-1] * np.sin(theta)
            unit_vec = np.linspace(0, 1, len(R_norm))
            beta = P_vec[ii, :] / P_vec[ii, 0]
            beta_in_diff[ii] = beta[-1] / beta_diff
            csi_init = [np.pi / 2, 0.0, - np.pi / 2, - np.pi]
            x_vec[ii, 0, 0] = x_in[0]
            y_vec[ii, 0, 0] = y_in[0]
            x_vec[ii, 0, 1] = x_in[249]
            y_vec[ii, 0, 1] = y_in[249]
            x_vec[ii, 0, 2] = x_in[499]
            y_vec[ii, 0, 2] = y_in[499]
            x_vec[ii, 0, 3] = x_in[749]
            y_vec[ii, 0, 3] = y_in[749]
            # compute flow path through the diffuser
            for kk in range(4):
                csi = csi_init[kk]
                for jj in range(len(R_norm) - 1):
                    gamma = np.arcsin(R_norm[jj] / R_norm[jj + 1] * np.sin(np.pi - alpha_diff[ii, jj]))
                    csi_new = alpha_diff[ii, jj] - gamma + csi
                    x_vec[ii, jj + 1, kk] = R_norm[jj + 1] * np.sin(csi_new)
                    y_vec[ii, jj + 1, kk] = R_norm[jj + 1] * np.cos(csi_new)
                    csi = csi_new
            gamma_Pv_in_diff[ii] = gamma_Pv_vec[ii, self.FindNearest(beta_in_diff[ii], beta)]
            ax1.plot(beta, gamma_Pv_vec[ii, :], lw=2, color=new_colors[ii], label=self.labels[ii])
            ax1.plot(beta, gamma_Pv_mean[ii] * np.ones(len(beta)),
                     linestyle='dashed', color=new_colors[ii], label=r'$\overline{\gamma}_{Pv}$')
            ax2.plot(beta, Z_vec[ii, :], lw=2, color=new_colors[ii], label=self.labels[ii])
            ax3.plot(unit_vec, M_diff[ii, :], lw=2, color=new_colors[ii], label=self.labels[ii])
            ax4.plot(x_in, y_in, lw=2, color='xkcd:grey')
            ax4.plot(x_out, y_out, lw=2, color=new_colors[ii], label=self.labels[ii])
            ax4.plot(x_vec[ii, :, 0], y_vec[ii, :, 0], linestyle='dashed', color=new_colors[ii])
            ax4.plot(x_vec[ii, :, 1], y_vec[ii, :, 1], linestyle='dashed', color=new_colors[ii])
            ax4.plot(x_vec[ii, :, 2], y_vec[ii, :, 2], linestyle='dashed', color=new_colors[ii])
            ax4.plot(x_vec[ii, :, 3], y_vec[ii, :, 3], linestyle='dashed', color=new_colors[ii])
            ax5.plot(unit_vec, P_diff[ii, :] / P_diff[ii, 0], lw=2, color=new_colors[ii], label=self.labels[ii])
            ax6.plot(unit_vec, np.rad2deg(alpha_diff[ii, :]), lw=2, color=new_colors[ii], label=self.labels[ii])
        ax1.plot(beta, gamma_id * np.ones(len(beta)), linestyle='dashed', color='xkcd:grey', label=r'$\gamma$')
        ax1.scatter(beta_in_diff, gamma_Pv_in_diff, marker='o', s=50, color='black', label='diffuser inlet')
        ax2.plot(beta, np.ones(len(beta)), linestyle='dashed', color='xkcd:grey', label='$Z=1$')

        ax1.grid(1)
        ax2.grid(1)
        ax3.grid(1)
        ax4.grid(1)
        ax5.grid(1)
        ax6.grid(1)
        
        ax1.set_xlabel(r'$\beta$ [-]')
        ax2.set_xlabel(r'$\beta$ [-]')
        ax3.set_xlabel(r'$R_\mathrm{norm}$ [-]')
        ax4.set_xlabel(r'$x_\mathrm{norm}$ [-]')
        ax5.set_xlabel(r'$R_\mathrm{norm}$ [-]')
        ax6.set_xlabel(r'$R_\mathrm{norm}$ [-]')
        
        ax1.set_ylabel(r'$\gamma_{Pv}$ [-]')
        ax2.set_ylabel('$Z$ [-]')
        ax3.set_ylabel('$M$ [-]')
        ax4.set_ylabel(r'$y_\mathrm{norm}$ [-]')
        ax5.set_ylabel(r'$\beta$ [-]')
        ax6.set_ylabel(r'$\alpha$ [deg]')
        ax4.axis('square')
        
        handles1, labels1 = ax1.get_legend_handles_labels()
        handles2, labels2 = ax2.get_legend_handles_labels()
        handles3, labels3 = ax3.get_legend_handles_labels()
        handles4, labels4 = ax4.get_legend_handles_labels()
        handles5, labels5 = ax5.get_legend_handles_labels()
        handles6, labels6 = ax6.get_legend_handles_labels()
        
        ax1.legend(handles1, labels1)
        ax2.legend(handles2, labels2)
        ax3.legend(handles3, labels3)
        ax4.legend(handles4, labels4, loc='upper right')
        ax5.legend(handles5, labels5)
        ax6.legend(handles6, labels6)
        
        fig1.savefig(self.jpeg_dir + '/gamma_Pv.jpeg')
        fig2.savefig(self.jpeg_dir + '/Z.jpeg')
        fig3.savefig(self.jpeg_dir + '/Mach_diff_radial.jpeg')
        fig4.savefig(self.jpeg_dir + '/radial_diffuser.jpeg')
        fig5.savefig(self.jpeg_dir + '/beta_diff_radial.jpeg')
        fig6.savefig(self.jpeg_dir + '/angle_diff_radial.jpeg')
        
        fig1.savefig(self.tiff_dir + '/gamma_Pv.tiff')
        fig2.savefig(self.tiff_dir + '/Z.tiff')
        fig3.savefig(self.tiff_dir + '/Mach_diff_radial.tiff')
        fig4.savefig(self.tiff_dir + '/radial_diffuser.tiff')
        fig5.savefig(self.tiff_dir + '/beta_diff_radial.tiff')
        fig6.savefig(self.tiff_dir + '/angle_diff_radial.tiff')
        
        plt.close(fig1)
        plt.close(fig2)
        plt.close(fig3)
        plt.close(fig4)
        plt.close(fig5)
        plt.close(fig6)

        return
