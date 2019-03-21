# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

# from generic_equation import GenEquation
import wenorec as wr
import numpy as np

class Flux:
    UPWIND = 0
    RUSANOV = 1

    def __init__(self, flux, order, is_wb=True, is_conservative=True):
        self.numflux = flux
        self.order = order
        self.wb = is_wb
        self.conservative = is_conservative # forced to True in main

    def flux(self, u, x, eqn, dt=0.1):
        """ Computes a numerical flux corresponding to F  s.t.
            u_t + (1/dx)(Fr - Fl) = 0
          Input:
              u: values of u in stencil of 2gw-1 cells centered in u_i
              x: spatial values corresponding to cells u
              eqn: object of class Equation
        """
        if self.wb: # is well balanced
            if self.conservative: # new conservative version
                if self.numflux == self.UPWIND:
                    return self.upwind(u, x, eqn)
                elif self.numflux == self.RUSANOV:
                    return self.rusanov(u, x, eqn)
            else: # old, non conservative version (kept for completeness)
                if self.numflux == self.UPWIND:
                    return self.upwind_nonconservative(u, x, eqn)
                elif self.numflux == self.RUSANOV:
                    return self.rusanov_nonconservative(u, x, eqn)
        else: # basic WENO
            if self.numflux == self.UPWIND:
                return self.upwind_nowb(u, x, eqn)
            elif self.numflux == self.RUSANOV:
                return self.rusanov_nowb(u, x, eqn)

    def upwind(self, u, x, eqn): # conservative
        nvars = eqn.dim()
        if nvars != 1:
            print "[ERROR] Upwind only implemented for scalar equations!"
            raise NotImplementedError
        Fl = np.zeros(1)
        Fr = np.zeros(1)
        Gl = np.zeros(1)
        Gr = np.zeros(1)
        i = (u.shape[1]-1)/2
        (critL, critR) = eqn.upw_criterion(u)
        ustar = eqn.steady_constraint(x[i], u[:,i], x)
        phi = eqn.F(ustar)
        
        Grm = wr.wenorec(self.order, phi[0,1:-1]) # at i+1/2^-
        Grp = wr.wenorec(self.order, phi[0,-1:1:-1]) # at i+1/2^+
        Glm = wr.wenorec(self.order, phi[0,0:-2]) # at i-1/2^-
        Glp = wr.wenorec(self.order, phi[0,-2:0:-1]) # at i-1/2^+
        Gr[0] = (critR >= 0)*Grm + (critR < 0)*Grp
        Gl[0] = (critL >= 0)*Glm + (critL < 0)*Glp

        phi = eqn.F(u)
        Frm = wr.wenorec(self.order, phi[0,1:-1]) # at i+1/2^-
        Frp = wr.wenorec(self.order, phi[0,-1:1:-1]) # at i+1/2^+
        Flm = wr.wenorec(self.order, phi[0,0:-2]) # at i-1/2^-
        Flp = wr.wenorec(self.order, phi[0,-2:0:-1]) # at i-1/2^+
        Fr[0] = (critR >= 0)*Frm + (critR < 0)*Frp
        Fl[0] = (critL >= 0)*Flm + (critL < 0)*Flp

        return (Fl-Gl, Fr-Gr)

    def rusanov(self, u, x, eqn): # conservative
        nvars = eqn.dim()
        i = (u.shape[1]-1)/2
        alpha = np.amax(np.abs(eqn.eig_of_dF(u)))
        Fl = np.zeros(nvars)
        Fr = np.zeros(nvars)
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)

        ustar = eqn.steady_constraint(x[i], u[:,i], x)
        phip = eqn.F(ustar) + alpha*ustar
        phim = eqn.F(ustar) - alpha*ustar
        Fp = eqn.F(u) + alpha*u
        Fm = eqn.F(u) - alpha*u

        for var in range(nvars):
            Grm = wr.wenorec(self.order, phip[var,1:-1]) # phip at i+1/2^-
            Glm = wr.wenorec(self.order, phip[var,0:-2]) # phip at i-1/2^-
            Grp = wr.wenorec(self.order, phim[var,-1:1:-1]) # phim at i+1/2^+
            Glp = wr.wenorec(self.order, phim[var,-2:0:-1]) # phim at i-1/2^+

            Frm = wr.wenorec(self.order, Fp[var,1:-1]) # Fp at i+1/2^-
            Flm = wr.wenorec(self.order, Fp[var,0:-2]) # Fp at i-1/2^-
            Frp = wr.wenorec(self.order, Fm[var,-1:1:-1]) # Fm at i+1/2^+
            Flp = wr.wenorec(self.order, Fm[var,-2:0:-1]) # Fm at i-1/2^+

            Gl[var] = 0.5*(Flm + Flp) - 0.5*(Glm + Glp)
            Gr[var] = 0.5*(Frm + Frp) - 0.5*(Grm + Grp)
        return (Gl, Gr)


    def upwind_nonconservative(self, u, x, eqn):
        nvars = eqn.dim()
        if nvars != 1:
            print "[ERROR] Upwind (nc) only implemented for scalar equations!"
            raise NotImplementedError
        Gl = np.zeros(1)
        Gr = np.zeros(1)
        i = (u.shape[1]-1)/2
        (critL, critR) = eqn.upw_criterion(u)
        ustar = eqn.steady_constraint(x[i], u[:,i], x)
        phi = eqn.F(u) - eqn.F(ustar)

        Grm = wr.wenorec(self.order, phi[0,1:-1]) # at i+1/2^-
        Grp = wr.wenorec(self.order, phi[0,-1:1:-1]) # at i+1/2^+
        Glm = wr.wenorec(self.order, phi[0,0:-2]) # at i-1/2^-
        Glp = wr.wenorec(self.order, phi[0,-2:0:-1]) # at i-1/2^+
        Gr[0] = (critR >= 0)*Grm + (critR < 0)*Grp
        Gl[0] = (critL >= 0)*Glm + (critL < 0)*Glp
        return (Gl, Gr)

    def rusanov_nonconservative(self, u, x, eqn):
        nvars = eqn.dim()
        i = (u.shape[1]-1)/2
        alpha = np.amax(np.abs(eqn.eig_of_dF(u)))
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
        ustar = eqn.steady_constraint(x[i], u[:,i], x)
        phip = eqn.F(u) - eqn.F(ustar) + alpha*(u - ustar)
        phim = eqn.F(u) - eqn.F(ustar) - alpha*(u - ustar)

        for var in range(nvars):
            Grm = wr.wenorec(self.order, phip[var,1:-1]) # phip at i+1/2^-
            Glm = wr.wenorec(self.order, phip[var,0:-2]) # phip at i-1/2^-

            Grp = wr.wenorec(self.order, phim[var,-1:1:-1]) # phim at i+1/2^+
            Glp = wr.wenorec(self.order, phim[var,-2:0:-1]) # phim at i-1/2^+
            Gl[var] = 0.5*(Glm + Glp)
            Gr[var] = 0.5*(Grm + Grp)
        return (Gl, Gr)


    def upwind_nowb(self, u, x, eqn):
        nvars = eqn.dim()
        if nvars != 1:
            print "[ERROR] Upwind (nc) only implemented for scalar equations!"
            raise NotImplementedError
        Gl = np.zeros(1)
        Gr = np.zeros(1)
        i = (u.shape[1]-1)/2
        (critL, critR) = eqn.upw_criterion(u)
        phi = eqn.F(u)

        Grm = wr.wenorec(self.order, phi[0,1:-1]) # at i+1/2^-
        Grp = wr.wenorec(self.order, phi[0,-1:1:-1]) # at i+1/2^+
        Glm = wr.wenorec(self.order, phi[0,0:-2]) # at i-1/2^-
        Glp = wr.wenorec(self.order, phi[0,-2:0:-1]) # at i-1/2^+
        Gr[0] = (critR >= 0)*Grm + (critR < 0)*Grp
        Gl[0] = (critL >= 0)*Glm + (critL < 0)*Glp
        return (Gl, Gr)

    def rusanov_nowb(self, u, x, eqn):
        nvars = eqn.dim()
        i = (u.shape[1]-1)/2
        alpha = np.amax(np.abs(eqn.eig_of_dF(u)))
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
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