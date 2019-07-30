# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d # UnivariateSpline
from scipy.optimize import newton
from Funciones_salto_estacionario import phi, phiu

class Equation:


    # Identifiers for equation
    LINEAR = 0   # eq_linear.py
    BURGERS = 1  # eq_burgers.py
    SW = 2       # eq_sw.py

    SEED = 11235813 # seed for reproducibility

    ###########################################################################
    # Following functions MUST be overridden by subclasses for anything to work:
    ###########################################################################

    def dim(self):
        """ Returns dimension of the problem: 1 for scalars """
        raise NotImplementedError

    def F(self, U):
        """ Flux function """
        raise NotImplementedError

    def dF(self, U):
        """ Derivative of flux function """
        raise NotImplementedError

    ###########################################################################
    # Following functions MUST be overridden for some purposes:
    ###########################################################################

    # def __init__(self):
        # """ Must be overridden if any initialization is required
            # In particular, if the class needs a random number, their __init__
            # needs to call np.random.seed(self.SEED) """
        # pass

    # def H(self, x):
        # """ Return H(x).

            # Required for SWE, other non-implemented equations """
        # raise NotImplementedError

    # def Hx(self, x):
        # """ Return H_x(x). 

            # Required for non-wb solver """
        # raise NotImplementedError

    def S(self, U):
        """ Return S(U). 

            Required for non-wb solver """
        raise NotImplementedError

    def upw_criterion(self, uStencil):
        """ Returns a pair (l,r) with the velocity for upwind criterion at
            left and right intercells, for cell at center of uStencil.

            Required for upwind """
        raise NotImplementedError


    def steady(self, x):
        """ Returns an arbitrary steady state for the equation.
            Input: 
                x: spatial coordinates
            Output:
                (nvars, len(x)) numpy array with the values
                If nvars = 1, this must still be a (1,len(x)) matrix;
                a len(x) array will not work!

            Required for InitCond.STEADY
        """
        raise NotImplementedError

    def steady_constraint(self, xConstr, uConstr, x):
        """ Returns a steady state solution of the equation u*, constrained
            to u*(xConstr) = uConstr
            Input:
                xConstr: double, x to fix the constraint
                uConstr: double or (dims,1) np array, u*(x) = uConstr for u* we look for
                x: values of x at which to evaluate u*
            Output:
                (nvars, len(x)) numpy array with the values
                If nvars = 1, this must still be a (1,len(x)) matrix;
                a len(x) array will not work!

            Required for well-balanced solver """
        raise NotImplementedError


    def Pi(self, V):
        """ Return the projection of vector V into the space of non-conservative
            subsystems. E.g. for SWE [h, q], only the first equation is
            conservative; therefore Pi([4,2]) returns [0,2].
            Defaults to assuming no conservative subsystems are present 
            (ie Pi(V)=V); but this makes the scheme lose conservativeness.

            Required for well-balanced solver to be conservative """
        print "[WARNING] Using default Pi; this will lose conservativeness"
        return V



    ###########################################################################
    # Following functions can be used in subclasses, but should be overridden:
    ###########################################################################

    def eig_of_dF(self, U):
        """ Returns eigenvalues of dF(U), as
            numpy array shaped like U, with
            eig[0,i] < eig[1,i] < ... for all i

            Strongly advised for Rusanov flux """
        print "[WARNING] Using default (very slow) eigenvalue computation"
        eig = np.zeros(U.shape)
        for i in range(U.shape[1]):
            # np.newaxis forces a shape (2,) numpy array to (2,1):
            X = np.linalg.eig(self.dF(U[:,i, np.newaxis]))[0]
            X.sort() # force eigenvalue order
            eig[:,i] = X
        return eig

    def prepare_plot(self,x,u,t):
        """ Plot x and u in whichever way is appropriate for the equation.
            This function will be called by io_manager.
            This function should produce a finished plot, including title,
            legend, labels and so on; io_manager will do plt.show() or savefig() 
            as required 

            Suggested for nicer plotting """
        print "[WARNING] Using default (unlabelled) plotting"
        for var in u.shape[0]:    
            plt.plot(x,u[var])
            plt.legend()
            plt.title(t)


    ###########################################################################
    # Following functions are convenience aliases and can be left alone:
    ###########################################################################
    def SHx(self, x, U):
        """ Return S(U) H_x(x) """
        return self.S(U)*self.Hx(x)

    def max_vel(self, U):
        """ Returns maximum velocity for CFL computation """
        return np.amax(np.abs(self.eig_of_dF(U)))