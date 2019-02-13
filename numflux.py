from equation import Equation
import wenorec as wr
import numpy as np

class Flux:
    UPWIND = 0
    RUSANOV = 1

    def __init__(self, flux, order, is_wb):
        self.numflux = flux
        self.order = order
        self.wb = is_wb

    def flux(self, u, x, eqn):
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
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
        i = (len(u)-1)/2
        (critL, critR) = eqn.upw_criterion(u)
        if self.wb:
            phi = eqn.g(u[:,i], u, x[i], x) # g evaluated at the stencil
        else:
            phi = phi = eqn.F(u)

        for var in range(nvars):
            Grm = wr.wenorec(self.order, phi[var,1:-1]) # at i+1/2^-
            Grp = wr.wenorec(self.order, phi[var,-1:1:-1]) # at i+1/2^+
            Glm = wr.wenorec(self.order, phi[var,0:-2]) # at i-1/2^-
            Glp = wr.wenorec(self.order, phi[var,-2:0:-1]) # at i-1/2^+
            Gr[var] = (critR >= 0)*Grm + (critR < 0)*Grp
            Gl[var] = (critL >= 0)*Glm + (critL < 0)*Glp
        return (Gl, Gr)

    def rusanov(self, u, x, eqn):
        nvars = eqn.dim()
        i = (len(u)-1)/2
        alpha = max(abs(eqn.dF(u)))
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
        if self.wb:
            phip = eqn.g(u[:,i], u, x[i], x) + alpha*u  # phi plus
            phim = eqn.g(u[:,i], u, x[i], x) - alpha*u  # phi minus
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
