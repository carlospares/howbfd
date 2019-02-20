import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline, interp1d

class Equation:
    # Identifiers for equation
    LINEAR = 0
    BURGERS = 1
    SWE_REST = 2 # 1D shallow water equation, vars [h,q=hu]

    SWE_H_FLAT = 0
    SWE_H_PWPOLY = 2

    linear_alpha = 0.05 # advection velocity for linear
    swe_g = 9.8
    swe_eta = 1

    SEED = 11235813 # for reproducibility

    def __init__(self, eqn, x, swe_H=SWE_H_FLAT, swe_noise_amplit=0):
        self.eq = eqn
        self.swe_H = swe_H
        self.swe_noise_amplit = swe_noise_amplit
        if swe_noise_amplit != 0:
            np.random.seed(self.SEED) # so we get consistent results
            Hnoise = swe_noise_amplit*np.random.rand(len(x))
            # self.Hnoiseinterp = UnivariateSpline(x, Hnoise, k=1) # allows to evaluate as a function
            self.Hnoiseinterp = interp1d(x, Hnoise, kind="nearest")
            # to do (or think?): this is convenient but might introduce machine-error

    def F(self, U):
        """ Flux function """
        if self.eq==Equation.LINEAR:
            return self.linear_alpha*U
        elif self.eq==Equation.BURGERS:
            return U*U/2
        elif self.eq==Equation.SWE_REST:
            ret = np.empty(U.shape)
            h = U[0,:]
            q = U[1,:]
            ret[0,:] = q
            ret[1,:] = q*q/h + 0.5*self.swe_g*h*h
            return ret

    def dF(self, U):
        """ Derivative of flux function """
        if self.eq==Equation.LINEAR:
            return self.linear_alpha*np.ones(np.size(U))
        elif self.eq==Equation.BURGERS:
            return U
        elif self.eq==Equation.SWE_REST or True:
            DF = np.zeros((U.shape[1], self.dim(), self.dim()))
            for i in range(U.shape[1]):
                h = U[0,i]
                q = U[1,i]
                DF[i] = np.array([[0,1],[-q*q/h*h + self.swe_g*h,2*q/h]])
            return DF

    def eig_of_dF(self, U):
        """ Returns eigenvalues of dF(U), as
            numpy array shaped like U, with
            eig[0,i] < eig[1,i] < ... for all i """
        if self.eq==Equation.LINEAR:
            return self.linear_alpha*np.ones(np.size(U))
        elif self.eq==Equation.BURGERS:
            return U
        elif self.eq==Equation.SWE_REST:
            eig = np.zeros(U.shape)
            eig[0,:] = U[1,:]/U[0,:] - np.sqrt(self.swe_g*U[0,:])
            eig[1,:] = U[1,:]/U[0,:] + np.sqrt(self.swe_g*U[0,:])
            return eig
        else: # default: slow, should be avoided!
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
        if self.eq == Equation.LINEAR:
            return self.linear_alpha*x
        elif self. eq == Equation.BURGERS:
            return x
        elif self.eq == Equation.SWE_REST:
            return self.swe_H_eval(x) # see later

    def Hx(self, x):
        """ Return H_x(x) """
        if self.eq == Equation.LINEAR:
            return self.linear_alpha*np.ones_like(x)
        elif self. eq == Equation.BURGERS:
            return np.ones_like(x)
        elif self.eq == Equation.SWE_REST:
            self.swe_Hx_eval(x) # see later

    def S(self, U):
        """ Return S(U) """
        if self.eq == Equation.LINEAR:
            return U
        elif self.eq == Equation.BURGERS:
            return U*U
        elif self.eq == Equation.SWE_REST:
            return np.array([ 0, self.swe_g*U[0] ])

    def upw_criterion(self, uStencil):
        """ Returns a pair (l,r) with the velocity for upwind criterion at
            left and right intercells, for cell at center of uStencil
        """
        if self.eq==Equation.LINEAR:
            return (np.sign(self.linear_alpha),np.sign(self.linear_alpha))
        elif self.eq==Equation.BURGERS:
            mid = int((len(uStencil)-1)/2) # approx. u at intercell
            return ((uStencil[0,mid]+uStencil[0,mid-1])/2, 
                    (uStencil[0,mid]+uStencil[0,mid+1])/2)
        elif self.eq==Equation.SWE_REST:
            print "[ERROR] Upwind only implemented for scalar equations!"
            sys.exit()

    def max_vel(self, U):
        """ Returns maximum velocity for CFL computation """
        return np.amax(np.abs(self.eig_of_dF(U)))

    def dim(self):
        """ Returns dimension of the problem: 1 for scalars """
        if self.eq != Equation.SWE_REST:
            return 1
        else:
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
        if self.eq in [Equation.LINEAR, Equation.BURGERS]:
            U0 = np.zeros((1, len(x)))
            U0[0] = np.exp(x)
        elif self.eq == Equation.SWE_REST:
            U0 = np.zeros((2, len(x)))
            U0[0,:] = self.swe_eta + self.swe_H_eval(x)
        return U0

    def steady_constraint(self, xConstr, uConstr, x):
        """ Returns a steady state solution of the equation u*, constrained
            to u*(xConstr) = uConstr
            Input:
                xConstr: x to fix the constraint
                uConstr: u*(x) = uConstr for u* we look for
                x: values of x at which to evaluate u*
            Output:
                (nvars, len(x)) numpy array with the values
                If nvars = 1, this must still be a (1,len(x)) matrix;
                a len(x) array will not work! """
        Ustar = np.zeros((self.dim(), len(x)))
        if self.eq in [Equation.LINEAR, Equation.BURGERS]:
            Ustar[0] = uConstr*np.exp(x - xConstr)
        elif self.eq == Equation.SWE_REST:
            if uConstr[1] != 0:
                print "[TO DO] Not yet implemented. Ignoring q for steady..."
                # sys.exit()
            # else:
            Ustar[0,:] = uConstr[0] - self.swe_H_eval(xConstr) + self.swe_H_eval(x)
        return Ustar

    def swe_H_eval(self, x):
        """ H(x) for SWE """
        if self.swe_H == Equation.SWE_H_FLAT:
            H = 0.1*np.ones_like(x)
        elif self.swe_H == Equation.SWE_H_PWPOLY:
            H = (0.13+0.05*(x-10)*(x-10))*(x>=8)*(x<=12)+0.33*((x<8)+(x>12))

        if self.swe_noise_amplit != 0:
            H += self.Hnoiseinterp(x)
        return H


    def swe_Hx_eval(self, x):
        """ H_x(x) for SWE """
        if self.swe_noise_amplit != 0:
            print "[ERROR] Tried to compute H_x but H has noise on, not smooth"
            sys.exit()

        if self.swe_H == Equation.SWE_H_FLAT:
            return np.zeros_like(x)
        elif self.swe_H == Equation.SWE_H_PWPOLY:
            return 0.1*(x-10)*(x>=8)*(x<=12)

    def prepare_plot(self,x,u,t):
        """ Plot x and u in whichever way is appropriate for the equation.
            This function will be called by io_manager.
            This function should produce a finished plot, including title,
            legend, labels and so on; io_manager will do plt.show() or savefig() 
            as required """
        if self.eq in [Equation.LINEAR, Equation.BURGERS]: # 1D
            plt.plot(x,u[0], label='u')
            plt.legend()
        elif self.eq == Equation.SWE_REST:
            plt.subplot(211)
            plt.title(t)
            H = self.swe_H_eval(x)
            plt.plot(x, -H, 'b', label='-H')
            plt.plot(x, u[0]-H, 'g', label='$\eta$')
            plt.plot(x, u[0], 'r', label='h')
            plt.legend()
            plt.subplot(212)
            plt.plot(x, u[1]/u[0], label='u')
            plt.legend()