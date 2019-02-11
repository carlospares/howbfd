import numpy as np

class Equation:
    # Identifiers for equation
    LINEAR = 0
    BURGERS = 1
    alpha = 0.05 # advection velocity for linear

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
            return self.alpha*(uj - ui*np.exp(xj-xi))
        elif self.eq==Equation.BURGERS:
            return uj*uj*0.5 - 0.5*ui*ui*np.exp(2*(xj-xi)) 


    def F(self, U):
        """ Flux function """
        if self.eq==Equation.LINEAR:
            return self.alpha*U
        elif self.eq==Equation.BURGERS:
            return U*U/2

    def dF(self, U):
        """ Derivative of flux function """
        if self.eq==Equation.LINEAR:
            return self.alpha
        elif self.eq==Equation.BURGERS:
            return U

    def SHx(self, U):
        if self.eq==Equation.LINEAR:
            return self.alpha*U
        elif self.eq==Equation.BURGERS:
            return U*U

    def upw_criterion(self, uLeft, uRight):
        """ Returns the velocity for upwind criterion at intercell
            between values uLefr and uRight
        """
        if self.eq==Equation.LINEAR:
            return np.sign(self.alpha)
        elif self.eq==Equation.BURGERS:
            return 0.5*(uLeft + uRight)

    def max_vel(self, u):
        """ Returns maximum velocity for CFL computation """
        if self.eq==Equation.LINEAR:
            return abs(self.alpha)
        elif self.eq==Equation.BURGERS:
            return abs(max(u))
