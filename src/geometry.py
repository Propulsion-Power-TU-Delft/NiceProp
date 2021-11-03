#! /usr/bin/env python
# -*- coding: utf-8 -*-

#############################################################################
# NiceProp - Interactively learning NICFD
# Authors: ir. A. Giuffre', Dr. ir. M. Pini
# Content: Geometrical definition of stationary flow devices
# 2021 - TU Delft - All rights reserved
#############################################################################


import numpy as np
import scipy.optimize as opt
from scipy.special import comb


class NozzleGeometry:
    """ Class to create the 2D profile of a converging or a converging-diverging nozzle """
    def __init__(self, inlet, throat, outlet, n_points):
        """
        Arguments:
        :param inlet: list collecting the (x,y) coordinates of the nozzle inlet
        :param throat: list collecting the (x,y) coordinates of the nozzle throat
        :param outlet: list collecting the (x,y) coordinates of the nozzle outlet
        :param n_points: number of points to be used to describe the nozzle shape

        Methods:
        - bernsteinPoly
        - bezierCurve
        - setBezierProfile
        - setTanhProfile
        - tanhNozzleProfile
        """
        # define inlet, throat, and outlet points, and horizontal lines passing through them
        self._inlet = Point(inlet[0], inlet[1])
        self._throat = Point(throat[0], throat[1])
        self._outlet = Point(outlet[0], outlet[1])
        self._angle_inlet = np.deg2rad(0.0)
        self._angle_throat = np.deg2rad(0.0)
        self._angle_outlet = np.deg2rad(0.0)
        self._angle_vertical = np.deg2rad(90.0)
        self._n_points = n_points

        self._inlet_line = Line(self._inlet, np.cos(self._angle_inlet), np.sin(self._angle_inlet))
        self._throat_line = Line(self._throat, np.cos(self._angle_throat), np.sin(self._angle_throat))
        self._outlet_line = Line(self._outlet, np.cos(self._angle_outlet), np.sin(self._angle_outlet))

        # define inlet-throat midpoint and vertical line passing through it
        self._in_throat = Point((inlet[0] + throat[0]) / 2, (inlet[1] + throat[1]) / 2)
        self._in_throat_line = Line(self._in_throat, np.cos(self._angle_vertical), np.sin(self._angle_vertical))

        # define throat-outlet midpoint and vertical line passing through it
        self._throat_out = Point((throat[0] + outlet[0]) / 2, (throat[1] + outlet[1]) / 2)
        self._throat_out_line = Line(self._throat_out, np.cos(self._angle_vertical), np.sin(self._angle_vertical))

        # define two additional control points for converging section
        self._converging1 = self._inlet_line.intersect(self._in_throat_line)
        self._converging2 = self._throat_line.intersect(self._in_throat_line)

        # define two additional control points for diverging section
        self._diverging1 = self._outlet_line.intersect(self._throat_out_line)
        self._diverging2 = self._throat_line.intersect(self._throat_out_line)

        # define control points for converging section
        self.cp_converging = np.array([self._inlet, self._converging1, self._converging2, self._throat])
        self.x_cp_converging = np.zeros(self.cp_converging.shape[0])
        self.y_cp_converging = np.zeros(self.cp_converging.shape[0])
        for ii in range(self.cp_converging.shape[0]):
            self.x_cp_converging[ii] = self.cp_converging[ii].x
            self.y_cp_converging[ii] = self.cp_converging[ii].y

        # define control points for diverging section
        self.cp_diverging = np.array([self._throat, self._diverging2, self._diverging1, self._outlet])
        self.x_cp_diverging = np.zeros(self.cp_diverging.shape[0])
        self.y_cp_diverging = np.zeros(self.cp_diverging.shape[0])
        for ii in range(self.cp_diverging.shape[0]):
            self.x_cp_diverging[ii] = self.cp_diverging[ii].x
            self.y_cp_diverging[ii] = self.cp_diverging[ii].y

        # initialize vector of coordinates for the entire nozzle
        self.x = np.zeros(self._n_points)
        self.y = np.zeros(self._n_points)
        self.x_norm = np.zeros(self._n_points)
        self.y_norm = np.zeros(self._n_points)

    def bernsteinPoly(self, i, n, t):
        """ Bernstein polynomial of n, i as a function of t """
        return comb(n, i) * (t ** (n - i)) * (1 - t) ** i

    def bezierCurve(self, cp, x_cp, y_cp):
        """
        Compute Bezier curve defined by a set of control points.
        The final curve is defined by n_points.
        Check http://processingjs.nihongoresources.com/bezierinfo/ for further information.
        """
        t = np.linspace(0.0, 1.0, self._n_points)
        polynomial_array = np.array([self.bernsteinPoly(i, cp.shape[0] - 1, t) for i in range(0, cp.shape[0])])

        x = np.dot(x_cp, polynomial_array)
        y = np.dot(y_cp, polynomial_array)

        return x[::-1], y[::-1]

    def setBezierProfile(self, profile='converging-diverging'):
        """
        Define 2D nozzle profile using 3rd order bezier curves.
        The bezier curves are controlled by the inlet, throat, and outlet coordinates,
        and by the angles at the same sections.
        """
        if profile == 'converging':
            self.x, self.y = self.bezierCurve(self.cp_converging, self.x_cp_converging, self.y_cp_converging)

        elif profile == 'converging-diverging':
            self._n_points = int(self._n_points / 2)
            x_converging, y_converging = self.bezierCurve(self.cp_converging, self.x_cp_converging, self.y_cp_converging)
            self._n_points += 1
            x_diverging, y_diverging = self.bezierCurve(self.cp_diverging, self.x_cp_diverging, self.y_cp_diverging)
            self.x = np.concatenate((x_converging, x_diverging[1:]))
            self.y = np.concatenate((y_converging, y_diverging[1:]))

        else:
            raise ValueError("The available choices for nozzle profile are 'converging' or 'converging-diverging'")

        self.x_norm = self.x / self._throat.y
        self.y_norm = self.y / self._throat.y

    def setTanhProfile(self):
        """
        Define converging-diverging 2D nozzle profile using hyperbolic tangent shapes.
        ***** This method is deprecated *****
        """
        # inlet to throat
        data = ((self._throat.x - self._inlet.x) / self._throat.y, 1, self._inlet.y / self._throat.y, 0.0075)
        a, b, c, d = opt.fsolve(self.tanhNozzleProfile, (0.1, 0.1, 0.1, 0.1),
                                args=data, full_output=False, xtol=1.0e-04)
        x_in = np.linspace(0, (self._throat.x - self._inlet.x) / self._throat.y, int(self._n_points / 2))
        self.x_norm[0:int(self._n_points / 2)] = x_in
        self.y_norm[0:int(self._n_points / 2)] = a + b * np.tanh(c * x_in[::-1] - d)

        # throat to outlet
        data = ((self._outlet.x - self._throat.x) / self._throat.y, 1, self._outlet.y / self._throat.y, 0.0075)
        a, b, c, d = opt.fsolve(self.tanhNozzleProfile, (0.1, 0.1, 0.1, 0.1),
                                args=data, full_output=False, xtol=1.0e-04)
        x_out = np.linspace(0, (self._outlet.x - self._throat.x) / self._throat.y, int(self._n_points / 2) + 1)
        self.x_norm[int(self._n_points / 2):] = x_out[1:] + x_in[-1]
        self.y_norm[int(self._n_points / 2):] = a + b * np.tanh(c * x_out[1:] - d)
        
    def tanhNozzleProfile(self, p, *data):
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


class Point:
    """
    Class to create 4D point used to generalize the calculation of rational curves.
    See documentation about Bspline and Bezier for further information.
    """
    def __init__(self, x=0.0, y=0.0, z=0.0, w=1.0):
        self.x = x
        self.y = y
        self.z = z
        self.w = w

    def __call__(self):
        return self.x, self.y, self.z

    def towards(self, target, t):
        return Point((1.0 - t) * self.x + t * target.x, (1.0 - t) * self.y + t * target.y,
                     (1.0 - t) * self.z + t * target.z)


class Line:
    """
    Class to create a parametric line.
    P -> Point; T -> Direction
    """
    def __init__(self, P, Tx, Ty, Tz=0.0):
        self.P = P
        self.T = [Tx, Ty, Tz]

    def __call__(self, u):
        return Point(self.P.x + self.T[0] * u, self.P.y + self.T[1] * u, self.P.z + self.T[2] * u)

    def intersect(self, other):
        P1, T1, P2, T2 = self.P, self.T, other.P, other.T
        u = (T2[1] * (P2.x - P1.x) - T2[0] * (P2.y - P1.y)) / (T1[0] * T2[1] - T1[1] * T2[0])
        return self(u)

