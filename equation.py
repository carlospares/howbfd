import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d # UnivariateSpline
from scipy.optimize import newton
from Funciones_salto_estacionario import phi, phiu

class Equation:
    # Identifiers for equation
    LINEAR = 0
    BURGERS = 1
    SW = 2 # 1D shallow water equation, vars [h,q=hu]

    SEED = 11235813 # for reproducibility

    def F(self, U):
        """ Flux function """
        raise NotImplementedError

    def dF(self, U):
        """ Derivative of flux function """
        raise NotImplementedError

    def eig_of_dF(self, U):
        """ Returns eigenvalues of dF(U), as
            numpy array shaped like U, with
            eig[0,i] < eig[1,i] < ... for all i """
        print "[WARNING] Using default (very slow) eigenvalue computation"
        eig = np.zeros(U.shape)
        for i in range(U.shape[1]):
            # np.newaxis forces a shape (2,) numpy array to (2,1):
            X = np.linalg.eig(self.dF(U[:,i, np.newaxis]))[0]
            X.sort() # force eigenvalue order
            eig[:,i] = X
        return eig

    def SHx(self, x, U):
        """ Return S(U) H_x(x) """
        return self.S(U)*self.Hx(x)

    def H(self, x):
        """ Return H(x) """
        raise NotImplementedError

    def Hx(self, x):
        """ Return H_x(x) """
        raise NotImplementedError

    def S(self, U):
        """ Return S(U) """
        raise NotImplementedError

    def upw_criterion(self, uStencil):
        """ Returns a pair (l,r) with the velocity for upwind criterion at
            left and right intercells, for cell at center of uStencil
        """
        raise NotImplementedError

    def max_vel(self, U):
        """ Returns maximum velocity for CFL computation """
        return np.amax(np.abs(self.eig_of_dF(U)))

    def dim(self):
        """ Returns dimension of the problem: 1 for scalars """
        raise NotImplementedError

    def steady(self, x):
        """ Returns an arbitrary steady state for the equation.
            Input: 
                x: spatial coordinates
            Output:
                (nvars, len(x)) numpy array with the values
                If nvars = 1, this must still be a (1,len(x)) matrix;
                a len(x) array will not work!
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
                a len(x) array will not work! """
        raise NotImplementedError

    def prepare_plot(self,x,u,t):
        """ Plot x and u in whichever way is appropriate for the equation.
            This function will be called by io_manager.
            This function should produce a finished plot, including title,
            legend, labels and so on; io_manager will do plt.show() or savefig() 
            as required """
        for var in u.shape[0]:    
            plt.plot(x,u[var])
            plt.legend()
            plt.title(t)