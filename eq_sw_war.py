# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 20:23:59 2019

@author: Usuario
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import newton
from Funciones_salto_estacionario import phi
from equation import Equation
from eq_sw import SWEquation
from nosteadyexc import NoSteadyError

class SWEquationWAR(SWEquation):
    """ 1D shallow water equation, vars [h,q=hu] 

        h_t +  q_x               = 0
        q_t + (q^2/h + gh^2/2)_x = gh H_x

    """
    def steady(self, H,x):
        U0 = np.zeros((self.dim(), len(H)))
        """ Returns an arbitrary steady state for the equation.
            Input: 
                x: spatial coordinates
            Output:
                (nvars, len(x)) numpy array with the values
                If nvars = 1, this must still be a (1,len(x)) matrix;
                a len(x) array will not work!
        """
        # # Water at rest solution: [h, q](x) = [eta + H(x), 0]
        # arbitrary_eta = 1
        # U0 = np.zeros((2, len(x)))
        # U0[0,:] = arbitrary_eta + self.H(x)
        # return U0
        eta = 2.
        U0[0]= eta + H
        return U0
    

    def steady_constraint(self, HConstr, uConstr, H,x, U0):
        """ Returns a steady state solution of the equation, u*, constrained
            to u*(xConstr) = uConstr
            Input:
                xConstr: double, x to fix the constraint
                uConstr: double or (dims,1) np array, u*(x) = uConstr for u* we look for
                x: values of x at which to evaluate u*
            Output:
                (nvars, len(x)) numpy array with the values
                If nvars = 1, this must still be a (1,len(x)) matrix;
                a len(x) array will not work! """
        Ustar = np.zeros((self.dim(), len(H)))
 
        (hi, qi, ui) = (uConstr[0], uConstr[1], uConstr[1]/uConstr[0])


        for j in range(len(H)):
            Ustar[0,j] = uConstr[0] - HConstr + H[j]
 
  
        return Ustar