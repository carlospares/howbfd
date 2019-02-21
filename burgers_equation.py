import numpy as np
import matplotlib.pyplot as plt
from equation import Equation

class BurgersEquation(Equation):

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

    def H(self, x):
        """ Return H(x) """
        return x

    def Hx(self, x):
        """ Return H_x(x) """
        return np.ones_like(x)

    def S(self, U):
        """ Return S(U) """
        return U*U

    def upw_criterion(self, uStencil):
        """ Returns a pair (l,r) with the velocity for upwind criterion at
            left and right intercells, for cell at center of uStencil
        """
        mid = int((len(uStencil)-1)/2) # approx. u at intercell
        return ((uStencil[0,mid]+uStencil[0,mid-1])/2, 
                (uStencil[0,mid]+uStencil[0,mid+1])/2)

    def dim(self):
        """ Returns dimension of the problem: 1 for scalars """
        return 1

    def steady(self, x):
        """ Returns an arbitrary steady state for the equation.
            Input: 
                x: spatial coordinates
            Output:
                (nvars, len(x)) numpy array with the values
                If nvars = 1, this must still be a (1,len(x)) matrix;
                a len(x) array will not work!
        """
        U0 = np.zeros((1, len(x)))
        U0[0] = np.exp(x)
        return U0

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
        Ustar = np.zeros((self.dim(), len(x)))
        Ustar[0] = uConstr*np.exp(x - xConstr)
        return Ustar

    def prepare_plot(self,x,u,t):
        """ Plot x and u in whichever way is appropriate for the equation.
            This function will be called by io_manager.
            This function should produce a finished plot, including title,
            legend, labels and so on; io_manager will do plt.show() or savefig() 
            as required """
        plt.plot(x,u[0], label='u')
        plt.legend()
        plt.title(t)