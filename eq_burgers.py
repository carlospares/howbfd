# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import numpy as np
import matplotlib.pyplot as plt
from equation import Equation
from functionH import FunH
from initcond import InitCond
from boundary import BoundaryCond

class BurgersEquation(Equation):
    """ 1D scalar Burgers' equation with mass term
    
        u_t + u u_x = S(u)Hx

    """

    def F(self, U):
        """ Flux function """
        return U*U/2

    def dF(self, U):
        """ Derivative of flux function """
        return U

    def eig_of_dF(self, U):
        """ Returns eigenvalues of dF(U), as
            numpy array shaped like U, with
            eig[0,i] < eig[1,i] < ... for all i """
        return U

    def S(self, U):
        """ Return S(U) """
        #return U*U # std burger's case
        #return U # std burger's case
        return (U - 1.0)# MMSburg case
        
    def discH_jumpF(self, ui, uip1, i, dH, x, t):
        # depends on S
        delta = 0.5*ui*ui*( np.exp( 2.0*dH ) - 1. )
        return delta

    def Piplus(self,ui, uip1):
        uip12 = .5*(ui + uip1)
        return 1.*(uip12 > 0 )
    
    def Piminus(self,ui, uip1):
        uip12 = .5*(ui  + uip1)
        return 1.*(uip12 < 0 )

        
    def dim(self):
        """ Returns dimension of the problem: 1 for scalars """
        return 1

    def Pi(self, V):
        """ Return projection of V into the space of non-conservative subsys.
            (i.e. everywhere) """
        return V

    def steady(self, H,x):
        """ Returns an arbitrary steady state for the equation.
            Input: 
                x: spatial coordinates
            Output:
                (nvars, len(x)) numpy array with the values
                If nvars = 1, this must still be a (1,len(x)) matrix;
                a len(x) array will not work!
        """
        HConstr = 0.
        uConstr = 1.
        return self.steady_constraint(HConstr, uConstr, H,x)

    def steady_constraint(self, HConstr, uConstr, H,x, U = 0):
        """ Returns a steady state solution of the equation u*, constrained
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
        Ustar[0] = uConstr*np.exp(H-HConstr)
        return Ustar

    def prepare_plot(self,x,u,H,t):
        """ Plot x and u in whichever way is appropriate for the equation.
            This function will be called by io_manager.
            This function should produce a finished plot, including title,
            legend, labels and so on; io_manager will do plt.show() or savefig() 
            as required """
        plt.plot(x,u[0], label='u')
        plt.legend()
        plt.title(t)
        
    def exact(self, x, t, H, params):
        U = np.zeros((self.dim(), len(x)))
        if params.funh == FunH.FLAT and params.init == InitCond.SHOCK and params.perturb_init == InitCond.PERT_NONE and params.boundary == BoundaryCond.IN_OUT:
            U[0] = 1. + 1.*(x < 0.5 + 1.5*t)
            return U
            
        if params.funh == FunH.FLAT and params.init == InitCond.RAREFACTION and params.perturb_init == InitCond.PERT_NONE and params.boundary == BoundaryCond.IN_OUT:
            U[0] = 1.*(x <= t) + (x/t)*(x>t)*(x<2*t) + 2.*(x >= 2*t)
            return U
            
        if params.funh == FunH.IDENT and params.init == InitCond.SHOCK and params.perturb_init == InitCond.PERT_NONE and params.boundary == BoundaryCond.IN_OUT:
            U[0] = 1. + 1.*(x < 0.5 + 1.5*t)
            return U
        if params.funh == FunH.IDENT and params.init == InitCond.STEADY and params.perturb_init == InitCond.PERT_NONE and (params.boundary in [BoundaryCond.FORCE_STEADY, BoundaryCond.FORCE_STEADY_ARBITRARY, BoundaryCond.FORCE_STEADY_INIT]):
            U[0] = np.exp(x)
            return U
        if params.funh == FunH.MMSburg and params.init == InitCond.MMSburg and params.perturb_init == InitCond.PERT_NONE and (params.boundary in [BoundaryCond.FORCE_STEADY, BoundaryCond.FORCE_STEADY_ARBITRARY, BoundaryCond.FORCE_STEADY_INIT]):
            U[0] = np.exp(-(x-5.0-t)*(x-5.0-t))
            return U
            
        else:
            return Equation.exact(self, x, t, H,params) # return default exact method from class Equation
