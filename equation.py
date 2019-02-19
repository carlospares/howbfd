import numpy as np
import sys

class Equation:
    # Identifiers for equation
    LINEAR = 0
    BURGERS = 1
    SWE_REST = 2 # 1D shallow water equation, vars [h,q=hu]

    SWE_H_FLAT = 0
    SWE_H_NOISE = 1

    linear_alpha = 0.05 # advection velocity for linear
    swe_g = 9.8
    swe_eta = 1

    def __init__(self, eqn, swe_H=SWE_H_FLAT):
        self.eq = eqn
        self.swe_H = swe_H

    def g(self,ui,uj,xi,xj):
        """ Input:
                ui: u[i] for i the center point of the stencil
                uj: u[j] for a single j, or array of values for all j in the stencil
                xi: x[i] for the center point of the stencil
                xj: x[j] (like uj). If array, it must be length(xj)==length(uj)
            Output:
                If uj is a number, returns g_i(x_j) 
                If uj is an array, returns [g_i(x_j) for every j in the stencil],
                                            which can be unpacked with * 
        """
        if self.eq==Equation.LINEAR:
            return self.linear_alpha*(uj - ui*np.exp(xj-xi))
        elif self.eq==Equation.BURGERS:
            return uj*uj*0.5 - 0.5*ui*ui*np.exp(2*(xj-xi))
        elif self.eq==Equation.SWE_REST:
            G = np.zeros(uj.shape)
            G[1,:] = self.swe_g*0.5*(self.swe_eta + self.swe_H_eval(xj))**2
            return self.F(uj) - G


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
            print U.shape
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
            print "Using default (very slow) eigenvalue computation"
            eig = np.zeros(U.shape)
            for i in range(U.shape[1]):
                # np.newaxis forces a shape (2,) numpy array to (2,1):
                X = np.linalg.eig(self.dF(U[:,i, np.newaxis]))[0]
                X.sort() # force eigenvalue order
                eig[:,i] = X
            return eig

    def SHx(self, x, U):
        return self.S(U)*self.Hx(x)

    def H(self, x):
        if self.eq == Equation.LINEAR:
            return self.linear_alpha*x
        elif self. eq == Equation.BURGERS:
            return x
        elif self.eq == Equation.SWE_REST:
            return self.swe_H_eval(x)

    def Hx(self, x):
        if self.eq == Equation.LINEAR:
            return self.linear_alpha*np.ones_like(x)
        elif self. eq == Equation.BURGERS:
            return np.ones_like(x)
        elif self.eq == Equation.SWE_REST:
            if self.swe_H == Equation.SWE_H_FLAT:
                return self.swe_H_eval(x)

    def S(self, U):
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
            mid = int((len(uStencil)-1)/2) # approx. u at intercell
            uL = (uStencil[1,mid]+uStencil[1,mid-1]) / \
                    (uStencil[0,mid]+uStencil[0,mid-1])
            uR = (uStencil[1,mid]+uStencil[1,mid+1]) / \
                    (uStencil[0,mid]+uStencil[0,mid+1])
            return (uL, uR)

    def max_vel(self, u):
        """ Returns maximum velocity for CFL computation """
        return np.amax(np.abs(self.eig_of_dF(U)))
        # if self.eq==Equation.LINEAR:
        #     return abs(self.linear_alpha)
        # elif self.eq==Equation.BURGERS:
        #     return np.amax(np.abs(u))
        # elif self.eq==Equation.SWE_REST:
        #     # TO DO 
        #     return 0.1 #np.amax(np.abs(u[1,:]/u[0,:])) # assume h != 0

    def dim(self):
        """ Returns dimension of the problem: 1 for scalars """
        if self.eq != Equation.SWE_REST:
            return 1
        else:
            return 2

    def steady(self, x):
        """ Returns a steady state for the equation.
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

    def swe_H_eval(self, x):
        if self.swe_H == Equation.SWE_H_FLAT:
            return 0.1*np.ones_like(x)
        elif self.swe_H == Equation.SWE_H_NOISE:
            np.random.seed(len(x)) # so we get consistent results
            return np.random.rand(len(x))


    def swe_Hx_eval(self, x):
        if self.swe_H == Equation.SWE_H_FLAT:
            return np.zeros_like(x)
        if self.swe_H == Equation.SWE_H_NOISE:
            print "Derivative of noise? Not happening, sorry"
            sys.exit() # if done, remove sys import