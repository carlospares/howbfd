# from generic_equation import GenEquation
import wenorec as wr
import numpy as np
import sys

class Flux:
    UPWIND = 0
    RUSANOV = 1

    def __init__(self, flux, order, is_wb):
        self.numflux = flux
        self.order = order
        self.wb = is_wb

    def flux(self, u, x, eqn, dt=0.1):
        """ Computes the numerical flux
          Input:
              u: values of u in stencil of 2gw-1 cells centered in u_i
              x: spatial values corresponding to cells u
              eqn: object of clas Equation
        """
        if self.numflux == self.UPWIND:
            return self.upwind(u, x, eqn)
        elif self.numflux == self.RUSANOV:
            return self.rusanov(u, x, eqn)

    def upwind(self, u, x, eqn):
        nvars = eqn.dim()
        if nvars != 1:
            print "[ERROR] Upwind only implemented for scalar equations!"
            sys.exit()
        Gl = np.zeros(1)
        Gr = np.zeros(1)
        i = (u.shape[1]-1)/2
        (critL, critR) = eqn.upw_criterion(u)
        if self.wb:
            ustar = eqn.steady_constraint(x[i], u[:,i], x)
            phi = eqn.F(u) - eqn.F(ustar)
            # phi = eqn.g(u[:,i], u, x[i], x) # g evaluated at the stencil
        else:
            phi = phi = eqn.F(u)

        Grm = wr.wenorec(self.order, phi[0,1:-1]) # at i+1/2^-
        Grp = wr.wenorec(self.order, phi[0,-1:1:-1]) # at i+1/2^+
        Glm = wr.wenorec(self.order, phi[0,0:-2]) # at i-1/2^-
        Glp = wr.wenorec(self.order, phi[0,-2:0:-1]) # at i-1/2^+
        Gr[0] = (critR >= 0)*Grm + (critR < 0)*Grp
        Gl[0] = (critL >= 0)*Glm + (critL < 0)*Glp
        return (Gl, Gr)

    def rusanov(self, u, x, eqn):
        nvars = eqn.dim()
        i = (u.shape[1]-1)/2
        alpha = np.amax(np.abs(eqn.eig_of_dF(u)))
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
        if self.wb:
            ustar = eqn.steady_constraint(x[i], u[:,i], x)
            phip = eqn.F(u) - eqn.F(ustar) + alpha*(u - ustar)
            phim = eqn.F(u) - eqn.F(ustar) - alpha*(u - ustar)
        else:
            phip = eqn.F(u) + alpha*u  # phi plus
            phim = eqn.F(u) - alpha*u  # phi minus

        for var in range(nvars):
            Grm = wr.wenorec(self.order, phip[var,1:-1]) # phip at i+1/2^-
            Glm = wr.wenorec(self.order, phip[var,0:-2]) # phip at i-1/2^-

            Grp = wr.wenorec(self.order, phim[var,-1:1:-1]) # phim at i+1/2^+
            Glp = wr.wenorec(self.order, phim[var,-2:0:-1]) # phim at i-1/2^+
            Gl[var] = 0.5*(Glm + Glp)
            Gr[var] = 0.5*(Grm + Grp)
        return (Gl, Gr)