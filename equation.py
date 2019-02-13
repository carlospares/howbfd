import numpy as np
import sys

class Equation:
    # Identifiers for equation
    LINEAR = 0
    BURGERS = 1
    SWE_REST = 2 # 1D shallow water equation, vars [h,q=hu]

    linear_alpha = 0.05 # advection velocity for linear
    swe_g = 9.8
    swe_eta = 1

    def __init__(self, eqn):
        self.eq = eqn

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
            G[1,:] = swe_g*0.5*(uj[0,:] )
            return self.F(uj) - G


    def F(self, U):
        """ Flux function """
        if self.eq==Equation.LINEAR:
            return self.linear_alpha*U
        elif self.eq==Equation.BURGERS:
            return U*U/2
        elif self.eq==Equation.SWE_REST:
            ret = np.empty(U.shape)
            q = U[0,:]
            h = U[1,:]
            ret[0,:] = q
            ret[1,:] = q*q/h + 0.5*g*h*h
            return ret

    def dF(self, U):
        """ Derivative of flux function """
        if self.eq==Equation.LINEAR:
            return self.linear_alpha*np.ones(np.size(U))
        elif self.eq==Equation.BURGERS:
            return U
        elif self.eq==Equation.SWE_REST:
            print "Rusanov for SWE at rest not implemented yet. Bye!"
            sys.exit() # remove import sys once done

    def SHx(self, U):
        if self.eq==Equation.LINEAR:
            return self.linear_alpha*U
        elif self.eq==Equation.BURGERS:
            return U*U
        elif self.eq==Equation.SWE_REST:
            print "no-wb for SWE at rest not implemented yet. Bye!"
            sys.exit() # remove import sys once done

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
        if self.eq==Equation.LINEAR:
            return abs(self.linear_alpha)
        elif self.eq==Equation.BURGERS:
            return np.amax(np.abs(u))
        elif self.eq==Equation.SWE_REST:
            return np.amax(np.abs(u[1,:]/u[0,:])) # assume h != 0

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
            h = (0.1 + 0.1*x)*(x < 0)  + (1+x)*(x >= 0)
            u = np.exp(0.1*x)*(x < 0) + (np.exp(0.9+x))*(x>=0)
            U0[0,:] = h
            U0[1,:] = h*u
        return U0