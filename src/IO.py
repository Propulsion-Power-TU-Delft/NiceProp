#! /usr/bin/env python
# -*- coding: utf-8 -*-

#############################################################################
# NiceProp - Interactively learning NICFD
# Authors: ir. A. Giuffre', Dr. ir. M. Pini
# Content: Functions for input/output
# 2021 - TU Delft - All rights reserved
#############################################################################


import numpy as np


def readConfigFile(config_file):
    """
    Read data from configuration file.

    :param config_file: name of the configuration file (including .txt); by default, it is assumed to be in input dir
    :return settings: dictionary collecting the information read from configuration file
    """
    settings = {}
    file = open('../input/' + config_file, 'r')
    file.readline()
    file.readline()
    settings['fluid'] = file.readline().rstrip('\n')
    file.readline()
    settings['equation of state'] = file.readline().rstrip('\n')
    file.readline()
    settings['thermodynamic plane'] = file.readline().rstrip('\n')
    file.readline()
    file.readline()
    file.readline()
    settings['process'] = file.readline().rstrip('\n')
    file.readline()
    settings['design'] = file.readline().split()[0]
    file.readline()
    settings['alpha'] = np.asarray(file.readline().split(), dtype=float)
    file.readline()
    settings['beta'] = np.asarray(file.readline().split(), dtype=float)
    file.readline()
    settings['inlet thermodynamic plane'] = file.readline().split()
    file.readline()
    settings['inlet input 1'] = np.asarray(file.readline().split(), dtype=float)
    file.readline()
    settings['inlet input 2'] = np.asarray(file.readline().split(), dtype=float)
    file.readline()
    settings['labels'] = file.readline().split()
    file.readline()
    settings['geometry flag'] = file.readline().split()[0]
    file.readline()
    file.readline()
    file.readline()
    settings['mass flow rate'] = np.asarray(file.readline().split(), dtype=float)
    file.readline()
    settings['velocity'] = np.asarray(file.readline().split(), dtype=float)
    file.readline()
    settings['geometry'] = file.readline().rstrip('\n')
    file.readline()
    file.readline()
    file.readline()
    settings['samples'] = int(file.readline())
    file.readline()
    Tr = file.readline().split()
    settings['Tr vec'] = np.linspace(float(Tr[0]), float(Tr[1]), settings['samples'])
    file.readline()
    sr = file.readline().split()
    settings['sr vec'] = np.linspace(float(sr[0]), float(sr[1]), settings['samples'])
    file.readline()
    Pr = file.readline().split()
    settings['Pr vec'] = np.linspace(float(Pr[0]), float(Pr[1]), settings['samples'])

    return settings


def isSubstring(s1, s2):
    """ Returns true if s1 is substring of s2 """
    m = len(s1)
    n = len(s2)
    jj = 0

    for ii in range(n - m + 1):
        for jj in range(m):
            if s2[ii + jj] != s1[jj]:
                break
        if (jj + 1) == m:
            return ii

    return -1


def readNozzleCoordinates(labels, samples):
    """ Read nozzle coordinates from a txt file, in terms of normalized length, radius and actual area variation """
    x_norm = np.ndarray((len(labels), samples))
    R_norm = np.ndarray((len(labels), samples))
    R_throat = np.zeros(len(labels))
    massFlow = np.zeros(len(labels))
    V_in = np.zeros(len(labels))
    geometry = [''] * len(labels)
    ii = -1
    f = open('input/nozzleCoordinates.txt', "r")
    lines = f.readlines()

    for line in lines:
        if isSubstring('NOZZLE GEOMETRY', line) != -1:
            ii += 1
            jj = 0
            pass
        elif isSubstring('label', line) != -1:
            pass
        elif isSubstring('R throat', line) != -1:
            R_throat[ii] = float(line.split()[2])
        elif isSubstring('mass flow', line) != -1:
            massFlow[ii] = float(line.split()[2])
        elif isSubstring('inlet velocity', line) != -1:
            V_in[ii] = float(line.split()[2])
        elif isSubstring('geometry', line) != -1:
            geometry[ii] = line.split()[1]
        elif isSubstring('x norm', line) != -1:
            pass
        else:
            x_norm[ii, jj] = float(line.split()[0])
            R_norm[ii, jj] = float(line.split()[1])
            jj += 1

    f.close()

    return massFlow, V_in, geometry, x_norm, R_norm, R_throat


def writeNozzleCoordinates(results_dir, labels, x_norm, R_norm, R_throat, massFlow, V_in, geometry):
    """
    Write nozzle coordinates in a txt file, in terms of normalized length, radius and actual throat radius,
    together with design mass flow rate, inlet velocity and type of nozzle (rectangular or circular).
    These information can be used as input for new runs of NiceProp or for further analysis, e.g. CFD simulations.
    """
    f = open(results_dir + '/nozzleCoordinates.txt', "w+")
    for ii in range(x_norm.shape[0]):
        f.write("********** NEW NOZZLE GEOMETRY **********" + "\n")
        f.write("label:          %s" % labels[ii] + "\n")
        f.write("R throat:       %2.8f" % R_throat[ii] + " m\n")
        f.write("mass flow:      %2.4f" % massFlow[ii] + " kg/s\n")
        f.write("inlet velocity: %3.2f" % V_in[ii] + " m/s\n")
        f.write("geometry:       %s" % geometry[ii] + "\n")
        f.write("x norm         R norm" + "\n")
        for jj in range(x_norm.shape[1]):
            f.write("%2.8f     %2.8f" % (x_norm[ii, jj], R_norm[ii, jj]) + "\n")
    f.close()

    return


def writeNozzleFlow(results_dir, labels, x_norm, P_vec, D_vec, M_vec):
    """
    Write nozzle flow properties in a txt file, in terms of normalized length, pressure, density and Mach
    """
    f = open(results_dir + '/nozzleFlow.txt', "w+")
    for ii in range(x_norm.shape[0]):
        f.write("********** NEW NOZZLE GEOMETRY **********" + "\n")
        f.write("label: %s" % labels[ii] + "\n")
        f.write("x norm [-]     P [Pa]              rho [kg/m3]     Mach [-]" + "\n")
        for jj in range(x_norm.shape[1]):
            f.write("%2.8f     %2.8f     %2.8f     %2.8f"
                    % (x_norm[ii, jj], P_vec[ii, jj], D_vec[ii, jj], M_vec[ii, jj]) + "\n")
    f.close()

    return

