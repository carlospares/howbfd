# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import numpy as np
import matplotlib.pyplot as plt
from equation import Equation
from initcond import InitCond
from boundary import BoundaryCond
from functionH import FunH

class LinearEquation(Equation):
    """ 1D scalar linear transport equation with mass term
    
        u_t + alpha u_x = u

    """

    alpha = 1.0

    def F(self, U):
        """ Flux function """
        return self.alpha*U

    def dF(self, U):
        """ Derivative of flux function """
        return self.alpha*np.ones(np.size(U))

    def eig_of_dF(self, U):
        """ Returns eigenvalues of dF(U), as
            numpy array shaped like U, with
            eig[0,i] < eig[1,i] < ... for all i """
        return self.alpha*np.ones(np.size(U))


    def S(self, U):
        """ Return S(U) """
        return U

    
    def Piplus(self, ui, uip1):
        return 1.*(self.alpha >  0)
    
    def Piminus(self, ui, uip1):
        return 1.*(self.alpha < 0)

    def dim(self):
        """ Returns dimension of the problem: 1 for scalars """
        return 1

    def Pi(self, V):
        """ Return projection of V into the space of non-conservative subsys.
            (i.e. everywhere) """
        return V

    def steady(self, H, x):
        """ Returns an arbitrary steady state for the equation.
            Input: 
                x: spatial coordinates
            Output:
                (nvars, len(x)) numpy array with the values
                If nvars = 1, this must still be a (1,len(x)) matrix;
                a len(x) array will not work!
        """
        return self.steady_constraint(0, 1, H, x)

    def steady_constraint(self, HConstr, uConstr, H, x, U = 0):
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
        Ustar[0] = uConstr*np.exp((H - HConstr)/self.alpha)
        return Ustar
        
    def exact(self, x, t, H, params):
        alpha = self.alpha
        x0 = x - alpha*t
        U = np.zeros((self.dim(), len(x)))
        if params.funh == FunH.IDENT and params.init == InitCond.STEADY and params.perturb_init == InitCond.PERT_PATCH and params.boundary == BoundaryCond.FORCE_STEADY:
            U[0] = (x <= (-0.5 + alpha*t))*np.exp(x/alpha) + (x >= (-0.3+alpha*t))*np.exp(x/alpha) +\
                    (x > (-0.5 + alpha*t))*(x < (-0.3+alpha*t))*(np.exp(t) + np.exp(x/alpha))
            return U
                    
        elif params.funh == FunH.IDENT and params.init == InitCond.STEADY and params.perturb_init == InitCond.PERT_MGAUSS and params.boundary == BoundaryCond.FORCE_STEADY:
            U[0] = (-0.3*np.exp(-200*x0*x0) + np.exp(x0/alpha))*np.exp(t)
            return U
            
        elif params.funh == FunH.IDENT and params.init == InitCond.SIN and params.perturb_init == InitCond.PERT_NONE and params.boundary == BoundaryCond.IN_OUT:
            U[0] = np.exp((x + 1)/alpha)*(x<= -1 + alpha*t) + (x > -1 + alpha*t)*(1+np.sin(2*np.pi*x0))*np.exp(t)
            return U
            
        elif params.funh == FunH.FLAT and params.init == InitCond.SIN and params.perturb_init == InitCond.PERT_NONE and params.boundary == BoundaryCond.IN_OUT:
            U[0] = 1*(x<= -1 + alpha*t) + (x > -1 + alpha*t)*(1+np.sin(2*np.pi*x0))
            return U
        elif params.init==InitCond.ORDER_TEST:
            U[0] = np.exp(1)*(((x-1)**6*(1 - 6*(x-2) + 21*(x-2)**2 - 56*(x-2)**3 + 126*(x-2)**4-252*(x-2)**5))*( x >1)*(x < 2) + 1.*(x >= 2))
#            U[0] = np.exp(1)*(((x-1)**5*(1 - 5*(x-2) + 15*(x-2)**2 - 35*(x-2)**3 + 70*(x-2)**4))*( x >1)*(x < 2) + 1.*(x >= 2))
#            U[0] =  np.exp(1)*(((x-1)**3 -3*(x-1)**3*(x-2) + 6*(x-1)**3*(x-2)**2)*(x > 1)*(x<2)+ 1.*(x>=2.))
#            U[0] = np.exp(1)*((-2*(x-1)**3 + 3.*(x-1)**2)*(x > 1)*(x<2)+ 1.*(x>=2.))
#            U[0] = np.exp(3.)*np.exp(-1./(1 -(x-3.)**2))*(x < 4)*(x > 2)
#            U[0] = np.exp(3.)*(1 -(x-3.)**2)*(x < 4)*(x > 2)
            return U
            
        else:
            return Equation.exact(self, x, t, H, params) # return default exact method from class Equation
        
            

    def prepare_plot(self,x,u,H,t):
        """ Plot x and u in whichever way is appropriate for the equation.
            This function will be called by io_manager.
            This function should produce a finished plot, including title,
            legend, labels and so on; io_manager will do plt.show() or savefig() 
            as required """
    
        plt.plot(x,u[0], label='u')
        plt.legend()
        plt.title(t)
        
        # alpha = self.alpha
        # exact = (x <= (-0.5 + alpha*t))*np.exp(x) + (x >= (-0.3+alpha*t))*np.exp(x) + (x > (-0.5 + alpha*t))*(x < (-0.3+alpha*t))*(np.exp(alpha*t) + np.exp(x))
        # x0 = x - alpha*t
        # exact = (-0.3*np.exp(-200*x0*x0) + np.exp(x0))*np.exp(0.05*t)
        # plt.plot(x, exact-u[0], 'r')