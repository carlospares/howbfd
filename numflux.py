from equation import Equation
import wenorec as wr

class Flux:
    UPWIND = 0
    RUSANOV = 1

    def __init__(self, flux, order, is_wb):
        self.numflux = flux
        self.order = order
        self.wb = is_wb

    # def flux(self, fluxL, fluxR, eq, uL, uR):
    #   """ Computes the numerical flux
    #       Input:
    #           fluxL: flux reconst. at left of interface
    #           flurR: same at right of interface
    #           eq: object of class Equation
    #           uL: current value at cell to left of interface
    #           uR: same at right of interface
    #   """
    #   if self.numflux == self.UPWIND:
    #       crit = eq.upw_criterion(uL, uR)
    #       return self.upwind(fluxL, fluxR, crit)
    #   elif self.numflux == self.RUSANOV:
    #       return self.rusanov(fluxL, fluxR, uL, uR, eq)

    def upwind(self, u, x, eqn):
        i = (len(u)-1)/2
        if self.wb:
            phi = eqn.g(u[i], u, x[i], x) # g evaluated at the stencil
        else:
            phi = phi = eqn.F(u)
        Grm = wr.wenorec(self.order, phi[1:-1]) # at i+1/2^-
        Grp = wr.wenorec(self.order, phi[-1:1:-1]) # at i+1/2^+
        Glm = wr.wenorec(self.order, phi[0:-2]) # at i-1/2^-
        Glp = wr.wenorec(self.order, phi[-2:0:-1]) # at i-1/2^+
        (critL, critR) = eqn.upw_criterion(u)
        Gr = (critR >= 0)*Grm + (critR < 0)*Grp
        Gl = (critL >= 0)*Glm + (critL < 0)*Glp
        return (Gl, Gr)

    def rusanov(self, u, x, eqn):
        i = (len(u)-1)/2
        alpha = max(abs(eqn.dF(u)))
        if self.wb:
            phip = eqn.g(u[i], u, x[i], x) + alpha*u  # phi plus
            phim = eqn.g(u[i], u, x[i], x) - alpha*u  # phi minus
        else:
            phip = eqn.F(u) + alpha*u  # phi plus
            phim = eqn.F(u) - alpha*u  # phi minus

        Grm = wr.wenorec(self.order, phip[1:-1]) # phip at i+1/2^-
        Glm = wr.wenorec(self.order, phip[0:-2]) # phip at i-1/2^-

        Grp = wr.wenorec(self.order, phim[-1:1:-1]) # phim at i+1/2^+
        Glp = wr.wenorec(self.order, phim[-2:0:-1]) # phim at i-1/2^+
        return (0.5*(Glm + Glp), 0.5*(Grm + Grp))


    def flux(self, u, x, eqn):
        if self.numflux == self.UPWIND:
            return self.upwind(u, x, eqn)
        elif self.numflux == self.RUSANOV:
            return self.rusanov(u, x, eqn)
