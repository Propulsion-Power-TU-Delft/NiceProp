#! /usr/bin/env python
# -*- coding: utf-8 -*-

#############################################################################
# NiceProp - Interactively learning NICFD
# Authors: ir. A. Giuffre', Dr. ir. M. Pini
# Content: GUI and main program
# 2021 - TU Delft - All rights reserved
#############################################################################


import time
import warnings
import pyfiglet
import tkinter as tk
from src.IO import *
import src.thermodynamics as thermo
import src.isentropic_process as isos


class App(tk.Frame):
    """ Class to take user input from a GUI and run the main program """
    def __init__(self, master=None, **kw):
        tk.Frame.__init__(self, master=master, **kw)
        tk.Label(self, text=" \n ", justify='left', font=('calibre', 12)).grid(row=0, column=0)
        tk.Label(self, text="Configuration file", justify='left', font=('calibre', 12)).grid(row=1, column=0)
        self.q1 = tk.Entry(self, font=('calibre', 12))
        self.q1.grid(row=1, column=1)
        tk.Label(self, text="Input name without extension", justify='left', font=('calibre', 12)).grid(row=2, column=0)
        tk.Label(self, text="   The file must be located in input directory",
                 justify='left', font=('calibre', 12)).grid(row=3, column=0)
        tk.Button(self, text="Run", command=self.run, justify='left', font=('calibre', 12)).grid(row=4, column=1)

    def run(self):
        """ Run the main program and retrieve the execution time """
        start_time = time.time()

        # Read configuration file
        settings = readConfigFile(self.q1.get() + '.txt')
        root.destroy()

        # Run
        warnings.filterwarnings("ignore")
        banner = pyfiglet.figlet_format("NiceProp", justify='auto')
        print("\n\n******************************************************")
        print(banner)
        print(" Interactively learning NICFD")
        print(" Authors: ir. A. Giuffre', Dr. ir. M. Pini")
        print(" Delft University of Technology - All rights reserved")
        print("******************************************************\n")

        thermodynamics = thermo.ThermodynamicModel(settings)
        if settings['thermodynamic plane'] == 'Ts':
            thermodynamics.TsContour()
        elif settings['thermodynamic plane'] == 'PT':
            thermodynamics.PTContour()
        else:
            raise ValueError("The available choices for 'Thermodynamic plane' are 'PT' or 'Ts'")

        if (settings['process'] == 'expansion') or (settings['process'] == 'compression'):
            flow = isos.IsentropicFlowModel(thermodynamics, settings)
            if settings['design'] == 'y' or settings['design'] == 'Y':
                if settings['process'] == 'expansion':
                    if settings['geometry flag'] == 'y' or settings['geometry flag'] == 'Y':
                        massflow, velocity, geometry, x_norm, R_norm, R_throat = \
                            readNozzleCoordinates(settings['labels'], settings['samples'])
                        flow.NozzleExpansion('off-design', x_norm=x_norm, R_norm=R_norm, R_throat_vec=R_throat)
                    else:
                        flow.NozzleExpansion('design')
                elif settings['process'] == 'compression':
                    flow.DiffuserCompression()
                else:
                    raise ValueError("The available choices for 'Thermodynamic process' are "
                                     "'expansion' or 'compression'")
            else:
                flow.IdealProcess()
            print("  Average value(s) of compressibility factor: " + str(flow.Z_mean))
            print("  Average value(s) of isentropic pressure-volume exponent: " + str(flow.gamma_Pv_mean))

        print("\n Elapsed time: %10.1f seconds" % (time.time() - start_time))


if __name__ == '__main__':
    root = tk.Tk()
    root.geometry("500x180")
    root.title('Welcome to NiceProp!')
    App(root).grid()
    root.mainloop()
