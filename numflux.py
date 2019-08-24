# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

# from generic_equation import GenEquation
import wenorec as wr
import numpy as np
from nosteadyexc import NoSteadyError

class Flux:
    UPWIND = 100
    RUSANOV = 101
    
    def __init__(self, cf):
        self.numflux = cf.numflux
        self.order = cf.order
        self.wb = cf.well_balanced
        self.conservative = cf.is_conservative

    def flux(self, u, x, H, eqn, dt=0.1):
        """ Computes a numerical flux corresponding to F  s.t.
            u_t + (1/dx)(Fr - Fl) = 0
          Input:
              u: values of u in stencil of 2gw-1 cells centered in u_i
              x: spatial values corresponding to cells u
              eqn: object of class Equation
        """
        if self.wb: # is well balanced
            if self.conservative: # new version, conservative where RHS is 0
                if self.numflux == self.UPWIND:
                    return self.upwind(u, x, H,eqn)
                elif self.numflux == self.RUSANOV:
                    return self.rusanov(u, x, H,eqn)
            else: # old, non conservative version (kept for completeness)
                if self.numflux == self.UPWIND:
                    return self.upwind_nonconservative(u, x, H,eqn)
                elif self.numflux == self.RUSANOV:
                    return self.rusanov_nonconservative(u, x, H,eqn)
        else: # basic WENO
            if self.numflux == self.UPWIND:
                return self.upwind_nowb(u, x, eqn)
            elif self.numflux == self.RUSANOV:
                return self.rusanov_nowb(u, x, eqn)

    def upwind(self, u, x, H,eqn):
        nvars = eqn.dim()
        if nvars != 1:
            print "[ERROR] Upwind only implemented for scalar equations!"
            raise NotImplementedError
        Gl = np.zeros(1)
        Gr = np.zeros(1)
        i = (u.shape[1]-1)/2
        (critL, critR) = eqn.upw_criterion(u)
        noSteady = 0
        try: 
            ustar = eqn.steady_constraint(H[i], u[:,i], H)
            phi = eqn.F(u) - eqn.Pi(eqn.F(ustar))
        except NoSteadyError, e: # no steady state exists! Default to basic WENO
            print "NoSteadyError triggered: {}".format(str(e))
            phi = eqn.F(u)
            noSteady = 1

        Grm = wr.wenorec(self.order, phi[0,1:-1]) # at i+1/2^-
        Grp = wr.wenorec(self.order, phi[0,-1:1:-1]) # at i+1/2^+
        Glm = wr.wenorec(self.order, phi[0,0:-2]) # at i-1/2^-
        Glp = wr.wenorec(self.order, phi[0,-2:0:-1]) # at i-1/2^+
        Gr[0] = (critR >= 0)*Grm + (critR < 0)*Grp
        Gl[0] = (critL >= 0)*Glm + (critL < 0)*Glp
        return (Gl, Gr, noSteady)

    def rusanov(self, u, x, H, eqn):
        nvars = eqn.dim()
        i = (u.shape[1]-1)/2
        alpha = np.amax(np.abs(eqn.eig_of_dF(u)))
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
        noSteady = 0
        try:
            ustar = eqn.steady_constraint(H[i], u[:,i], H)
            phip = eqn.F(u) - eqn.Pi(eqn.F(ustar)) + alpha*(u - ustar)
            phim = eqn.F(u) - eqn.Pi(eqn.F(ustar)) - alpha*(u - ustar)
        except NoSteadyError, e: # no steady state exists! Default to basic WENO
            #print "NoSteadyError triggered: {}".format(str(e))
            phip = eqn.F(u) + alpha*u
            phim = eqn.F(u) - alpha*u
            noSteady = 1

        for var in range(nvars):
            Grm = wr.wenorec(self.order, phip[var,1:-1]) # phip at i+1/2^-
            Glm = wr.wenorec(self.order, phip[var,0:-2]) # phip at i-1/2^-

            Grp = wr.wenorec(self.order, phim[var,-1:1:-1]) # phim at i+1/2^+
            Glp = wr.wenorec(self.order, phim[var,-2:0:-1]) # phim at i-1/2^+
            Gl[var] = 0.5*(Glm + Glp)
            Gr[var] = 0.5*(Grm + Grp)
        return (Gl, Gr, noSteady)


    def upwind_nonconservative(self, u, x, H, eqn):
        nvars = eqn.dim()
        if nvars != 1:
            print "[ERROR] Upwind (nc) only implemented for scalar equations!"
            raise NotImplementedError
        Gl = np.zeros(1)
        Gr = np.zeros(1)
        i = (u.shape[1]-1)/2
        (critL, critR) = eqn.upw_criterion(u)
        ustar = eqn.steady_constraint(H[i], u[:,i], H)
        phi = eqn.F(u) - eqn.F(ustar)

        Grm = wr.wenorec(self.order, phi[0,1:-1]) # at i+1/2^-
        Grp = wr.wenorec(self.order, phi[0,-1:1:-1]) # at i+1/2^+
        Glm = wr.wenorec(self.order, phi[0,0:-2]) # at i-1/2^-
        Glp = wr.wenorec(self.order, phi[0,-2:0:-1]) # at i-1/2^+
        Gr[0] = (critR >= 0)*Grm + (critR < 0)*Grp
        Gl[0] = (critL >= 0)*Glm + (critL < 0)*Glp
        return (Gl, Gr)

    def rusanov_nonconservative(self, u, x, H, eqn):
        nvars = eqn.dim()
        i = (u.shape[1]-1)/2
        alpha = np.amax(np.abs(eqn.eig_of_dF(u)))
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
        ustar = eqn.steady_constraint(H[i], u[:,i], H)
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