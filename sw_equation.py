import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d # UnivariateSpline
from scipy.optimize import newton
from Funciones_salto_estacionario import phi, phiu
from equation import Equation

class SWEquation(Equation):
    # Identifiers for SW topography
    H_FLAT = 0
    H_PWPOLY = 2

    # other function parameters
    g = 9.8

    def __init__(self, x, H=H_FLAT, noise_amplit=0):
        self.H_type = H
        self.noise_amplit = noise_amplit
        if noise_amplit != 0:
            np.random.seed(self.SEED) # so we get consistent results
            Hnoise = noise_amplit*np.random.rand(len(x))
            # self.Hnoiseinterp = UnivariateSpline(x, Hnoise, k=1) # probably faster than interp1d! To do: why does it break?
            self.Hnoiseinterp = interp1d(x, Hnoise, kind="nearest") # allows to evaluate as a function
            # to do (or think?): this is convenient but might introduce machine-error

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

    def H(self, x):
        """ Return H(x) """
        if self.H_type == self.H_FLAT:
            H = 0.1*np.ones_like(x)
        elif self.H_type == self.H_PWPOLY:
            H = (0.13+0.05*(x-10)*(x-10))*(x>=8)*(x<=12)+0.33*((x<8)+(x>12))

        if self.noise_amplit != 0:
            H += self.Hnoiseinterp(x)
        return H

    def Hx(self, x):
        """ Return H_x(x) """
        if self.noise_amplit != 0:
            print "[ERROR] Tried to compute H_x but H has noise on, not smooth"
            raise NotImplementedError

        if self.H_type == self.H_FLAT:
            return np.zeros_like(x)
        elif self.H_type == self.H_PWPOLY:
            return 0.1*(x-10)*(x>=8)*(x<=12)

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

    def steady(self, x):
        """ Returns an arbitrary steady state for the equation.
            Input: 
                x: spatial coordinates
            Output:
                (nvars, len(x)) numpy array with the values
                If nvars = 1, this must still be a (1,len(x)) matrix;
                a len(x) array will not work!
        """
        # Water at rest solution: [h, q](x) = [eta + H(x), 0]
        arbitrary_eta = 1
        U0 = np.zeros((2, len(x)))
        U0[0,:] = arbitrary_eta + self.H(x)
        return U0

    def steady_constraint(self, xConstr, uConstr, x):
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
        Ustar = np.zeros((self.dim(), len(x)))
        # if uConstr[1] != 0:
        #     Ustar[1] = uConstr[1]*np.ones_like(x)
        #     Ustar[0] = self.solve_steady_poly(uConstr[1], uConstr[0], 
        #                                           self.H(x))
        #     print "[TO DO] Debugging because of {}".format(uConstr[1])
        # else:
        #     Ustar[0,:] = uConstr[0] - self.H(xConstr) + self.H(x)
        Hj = self.H(x)
        for j in range(len(x)):
            (hsuperc, hsubc) = phi(uConstr[0], uConstr[1]/uConstr[0], self.H(xConstr), Hj[j])
            usuperc = phiu(uConstr[0], uConstr[1]/uConstr[0], hsuperc)
            usubc = phiu(uConstr[0], uConstr[1]/uConstr[0], hsubc)
            # Ustar[0,j] = hsuperc
            # Ustar[1,j] = usuperc
            Ustar[0,j] = hsubc
            Ustar[1,j] = usubc
        return Ustar

    # def solve_steady_poly(self, qi, hi, Hj):
    #     """ find h*(x) at x_j so that (h*, qi) is a steady state.
    #         Input:
    #             qi: double, qi = q(x_i) = q*(x_i) = q(x_j) forall j in stencil
    #             hi: double, hi = h(x_i) = h*(x_i)
    #             Hj: numpy array with [H(x_j)] forall j in stencil
    #         Output:
    #             numpy array with [h*(x_j)] forall j in stencil
    #     """
    #     i = (len(Hj)-1)/2 # central point is the constraint
    #     Ci = 0.5*qi*qi/hi/hi + self.g*hi - self.g*Hj[i]
    #     hstar = np.empty_like(Hj)
    #     for j in range(len(Hj)):
    #         hstar[j] = newton(steady_poly, hi, args=(qi, Ci, Hj[j], self.g), 
    #                           fprime=steady_poly_prime, 
    #                           fprime2=steady_poly_second) # for Halley's method
    #     return hstar

    def prepare_plot(self,x,u,t):
        """ Plot x and u in whichever way is appropriate for the equation.
            This function will be called by io_manager.
            This function should produce a finished plot, including title,
            legend, labels and so on; io_manager will do plt.show() or savefig() 
            as required """
        plt.subplot(211)
        plt.title(t)
        H = self.H(x)
        plt.plot(x, -H, 'b', label='-H')
        plt.plot(x, u[0]-H, 'g', label='$\eta$')
        plt.plot(x, u[0], 'r', label='h')
        plt.legend()
        plt.subplot(212)
        plt.plot(x, u[1]/u[0], label='u')
        plt.legend()

# def steady_poly(h, q, Ci, Hj, g):
#     """ Polynomial in h, such that P(h) = 0 iff h is a value at x_j """
#     return h*h*h - (Hj+ Ci/g)*h*h + 0.5*q*q/g
# def steady_poly_prime(h, q, Ci, Hj, g):
#     """ Derivative in h of steady_poly """
#     return 3*h*h - 2*(Hj + Ci/g)*h
# def steady_poly_second(h, q, Ci, Hj, g):
#     """ Second derivative in h of steady_poly """
#     return 6*h - 2*(Hj + Ci/g)