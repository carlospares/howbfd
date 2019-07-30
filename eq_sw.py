# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import newton
from Funciones_salto_estacionario import phi
from equation import Equation

class SWEquation(Equation):
    """ 1D shallow water equation, vars [h,q=hu] 

        h_t +  q_x               = 0
        q_t + (q^2/h + gh^2/2)_x = gh H_x

    """
    
    # other function parameters
    g = 9.8

    def F(self, U):
        """ Flux function """
        ret = np.empty(U.shape)
        h = U[0,:]
        q = U[1,:]
        ret[0,:] = q
        ret[1,:] = q*q/h + 0.5*self.g*h*h
        return ret

    def dF(self, U):
        """ Derivative of flux function """
        DF = np.zeros((U.shape[1], self.dim(), self.dim()))
        for i in range(U.shape[1]):
            h = U[0,i]
            q = U[1,i]
            DF[i] = np.array([[0,1],[-q*q/h*h + self.g*h,2*q/h]])
        return DF

    def eig_of_dF(self, U):
        """ Returns eigenvalues of dF(U), as
            numpy array shaped like U, with
            eig[0,i] < eig[1,i] < ... for all i """
        eig = np.zeros(U.shape)
        eig[0,:] = U[1,:]/U[0,:] - np.sqrt(self.g*U[0,:])
        eig[1,:] = U[1,:]/U[0,:] + np.sqrt(self.g*U[0,:])
        return eig


    def S(self, U):
        """ Return S(U) """
        return np.array([ 0, self.g*U[0] ])

    def upw_criterion(self, uStencil):
        """ Returns a pair (l,r) with the velocity for upwind criterion at
            left and right intercells, for cell at center of uStencil
        """
        print "[ERROR] Upwind only implemented for scalar equations!"
        raise NotImplementedError

    def dim(self):
        """ Returns dimension of the problem: 1 for scalars """
        return 2

    def steady(self, H):
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
        HConst = 1.
        uConst = [1., 5.]
        return self.steady_constraint(HConst, uConst, H)

    def Froude(self, u, h):
        return abs(u) / np.sqrt(self.g*h)

    def Pi(self, V):
        """ Return projection of V into the space of non-conservative subsys.
            (i.e. Pi([4,2]) = [0,2]) """
        PV = np.zeros_like(V)
        PV[1] = V[1]
        return PV

    def steady_constraint(self, HConstr, uConstr, H):
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
        # polyNewton = self.solve_steady_poly_newton(qi, hi, Hj)
        # print self.find_steady_poly_roots(qi, hi, Hj)

        Fr_i = self.Froude(ui, hi)
        for j in range(len(H)):
            (hsuperc, hsubc) = phi(hi, ui, HConstr, H[j])
            # hstar = polyNewton[j] # Halley's method
            Ustar[0,j] = hsuperc if Fr_i > 1 else hsubc
            Ustar[1,j] = uConstr[1]
        return Ustar

    # # TO DO: very slow. Delete me?
    # def find_steady_poly_roots(self, qi, hi, Hj):
    #     i = (len(Hj)-1)/2 # central point is the constraint
    #     Ci = 0.5*qi*qi/hi/hi + self.g*(hi - Hj[i])
    #     # hstar = np.empty_like(Hj)
    #     # for j in range(len(Hj)):
    #     #     hstar[j] = max(np.roots([self.g, -(Ci+self.g*Hj[j]), 0, 0.5*qi*qi]))
    #     # return hstar
    #     hstars = np.zeros((len(Hj), 3))
    #     for j in range(len(Hj)):
    #         hstars[j] = np.roots([self.g, -(Ci+self.g*Hj[j]), 0, 0.5*qi*qi])
    #     return hstars
    
    # def solve_steady_poly_newton(self, qi, hi, Hj):
    #     """ find h*(x) at x_j so that (h*, qi) is a steady state.
    #         Input:
    #             qi: double, qi = q(x_i) = q*(x_i) = q(x_j) forall j in stencil
    #             hi: double, hi = h(x_i) = h*(x_i)
    #             Hj: numpy array with [H(x_j)] forall j in stencil
    #         Output:
    #             numpy array with [h*(x_j)] forall j in stencil
    #     """
    #     i = (len(Hj)-1)/2 # central point is the constraint
    #     Ci = 0.5*qi*qi/hi/hi + self.g*(hi - Hj[i])
    #     hstar = np.empty_like(Hj)
    #     for j in range(len(Hj)):
    #         hstar[j] = newton(steady_poly, hi, args=(qi, Ci, Hj[j], self.g), 
    #                           fprime=steady_poly_prime, 
    #                           fprime2=steady_poly_second) # for Halley's method
    #     return hstar

    def prepare_plot(self,x,u,H,t):
        """ Plot x and u in whichever way is appropriate for the equation.
            This function will be called by io_manager.
            This function should produce a finished plot, including title,
            legend, labels and so on; io_manager will do plt.show() or savefig() 
            as required """
        plt.subplot(211)
        plt.title(t)
        plt.plot(x, -H, 'b', label='-H')
        plt.plot(x, u[0]-H, 'g', label='$\eta$')
#        plt.plot(x, u[0], 'r', label='h')
        plt.legend()
        plt.subplot(212)
        plt.plot(x, u[1]/u[0], label='u')
        plt.legend()


# def steady_poly(h, q, Ci, Hj, g):
#     """ Polynomial in h, such that P(h) = 0 iff h is a value at x_j """
#     return g*h*h*h - (Ci+g*Hj)*h*h + 0.5*q*q
# def steady_poly_prime(h, q, Ci, Hj, g):
#     """ Derivative in h of steady_poly """
#     return 3*g*h*h - 2*(Ci+g*Hj)*h
# def steady_poly_second(h, q, Ci, Hj, g):
#     """ Second derivative in h of steady_poly """
#     return 6*g*h - 2*(Ci+g*Hj)