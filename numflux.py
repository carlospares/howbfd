# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

# from generic_equation import GenEquation
import wenorec as wr
import numpy as np
from nosteadyexc import NoSteadyError

class Flux:
    UPWIND = 100
    RUSANOV = 101
    RUSANOVG = 102#Rusanov with global alpha
    SW_SPLIT = 103
    
    def __init__(self, cf):
        self.numflux = cf.numflux
        self.order = cf.order
        self.wb = cf.well_balanced
        self.conservative = cf.is_conservative

    def flux(self, u, x, H, alpha, eqn, dt=0.1):
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
                elif self.numflux == self.RUSANOVG:
                    return self.rusanovg(u, x, H,alpha, eqn)
                elif self.numflux == self.SW_SPLIT:
                    return self.sw_split(u,x,H,eqn)
            else: # old, non conservative version (kept for completeness)
                if self.numflux == self.UPWIND:
                    return self.upwind_nonconservative(u, x, H,eqn)
                elif self.numflux == self.RUSANOV:
                    return self.rusanov_nonconservative(u, x, H,eqn)
                elif self.numflux == self.RUSANOVG:
                    return self.rusanovg_nonconservative(u, x, H,alpha, eqn)
                elif self.numflux == self.SW_SPLIT:
                    return self.sw_split_nonconservative(u,x,H,eqn)
        else: # basic WENO
            if self.numflux == self.UPWIND:
                return self.upwind_nowb(u, x, eqn)
            elif self.numflux == self.RUSANOV:
                return self.rusanov_nowb(u, x, eqn)
            elif self.numflux == self.RUSANOVG:
                return self.rusanovg_nowb(u, x, alpha, eqn)
            elif self.numflux == self.SW_SPLIT:
                return self.sw_split_nowb(u,x,eqn)

    def upwind(self, u, x, H, eqn):
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
            ustar = eqn.steady_constraint(H[i], u[:,i], H,x)
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

#    def rusanov(self, u, x, H, eqn):
#        nvars = eqn.dim()
#        i = (u.shape[1]-1)/2
#        alpha = np.amax(np.abs(eqn.eig_of_dF(u)))
#        Gl = np.zeros(nvars)
#        Gr = np.zeros(nvars)
#        noSteady = 0
#    
#        try:
#            ustar = eqn.steady_constraint(H[i], u[:,i], H,x, u)
#            phip = eqn.F(u) - eqn.Pi(eqn.F(ustar)) + alpha*(u - eqn.Pi(ustar))  # phi plus
#            phim = eqn.F(u) - eqn.Pi(eqn.F(ustar)) - alpha*(u - eqn.Pi(ustar)) # phi minus
#
#
#        except NoSteadyError, e: # no steady state exists! Default to basic WENO
#            #print "NoSteadyError triggered: {}".format(str(e))
#            phip = eqn.F(u) + alpha*u  # phi plus
#            phim = eqn.F(u) - alpha*u # phi minus
#            noSteady = 1
#
#        for var in range(nvars):
#            Grm = wr.wenorec(self.order, phip[var,1:-1]) # phip at i+1/2^-
#            Glm = wr.wenorec(self.order, phip[var,0:-2]) # phip at i-1/2^-
#
#            Grp = wr.wenorec(self.order, phim[var,-1:1:-1]) # phim at i+1/2^+
#            Glp = wr.wenorec(self.order, phim[var,-2:0:-1]) # phim at i-1/2^+
##            print 'wb'
##            print phip[var,1:-1],  phip[var,0:-2],phim[var,-1:1:-1], phim[var,-2:0:-1]
#            Gl[var] = 0.5*(Glm + Glp)
#            Gr[var] = 0.5*(Grm + Grp)
#
#        return (Gl, Gr, noSteady)
    
    def rusanov(self, u, x, H, eqn):
        tol = 1.e-12
        nvars = eqn.dim()
        i = (u.shape[1]-1)/2
        alpha = np.amax(np.abs(eqn.eig_of_dF(u)))
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
        noSteady = 0
        I = nvars - np.sum(eqn.Pi(np.ones(nvars)))
        F = eqn.F(u)

    
        try:
            ustar = eqn.steady_constraint(H[i], u[:,i], H,x, u)
            phip = eqn.F(u) - eqn.Pi(eqn.F(ustar)) + alpha*(u - eqn.Pi(ustar))  # phi plus
            phim = eqn.F(u) - eqn.Pi(eqn.F(ustar)) - alpha*(u - eqn.Pi(ustar)) # phi minus


        except NoSteadyError, e: # no steady state exists! Default to basic WENO
            #print "NoSteadyError triggered: {}".format(str(e))
            phip = eqn.F(u) + alpha*u  # phi plus
            phim = eqn.F(u) - alpha*u # phi minus
            noSteady = 1

        for var in range(nvars):
            if var < I and max(abs(F[var,0:-1] - F[var,i]))< tol:
                Glm = F[var,i]
                Grm = F[var,i]
            else:        
                Grm = wr.wenorec(self.order, phip[var,1:-1]) # phip at i+1/2^-
                Glm = wr.wenorec(self.order, phip[var,0:-2]) # phip at i-1/2^-

            if  var < I and max(abs(F[var,1:] - F[var,i]))< tol:  
                Grp = F[var,i]
                Glp = F[var,i]
            else:
                Grp = wr.wenorec(self.order, phim[var,-1:1:-1]) # phim at i+1/2^+
                Glp = wr.wenorec(self.order, phim[var,-2:0:-1]) # phim at i-1/2^+
#          print 'wb'
#            print phip[var,1:-1],  phip[var,0:-2],phim[var,-1:1:-1], phim[var,-2:0:-1]
            Gl[var] = 0.5*(Glm + Glp)
            Gr[var] = 0.5*(Grm + Grp)

        return (Gl, Gr, noSteady)

    def rusanovg(self, u, x, H, alpha, eqn):
        tol = 1.e-12
        nvars = eqn.dim()
        i = (u.shape[1]-1)/2
        I = nvars - np.sum(eqn.Pi(np.ones(nvars)))
#        alpha = np.amax(np.abs(eqn.eig_of_dF(u)))
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
        noSteady = 0
        F = eqn.F(u)
        
        try:
            ustar = eqn.steady_constraint(H[i], u[:,i], H,x, u)
            phip = F - eqn.Pi(eqn.F(ustar)) + alpha*(u - eqn.Pi(ustar))  # phi plus
            phim = F - eqn.Pi(eqn.F(ustar)) - alpha*(u - eqn.Pi(ustar)) # phi minus


        except NoSteadyError, e: # no steady state exists! Default to basic WENO
            #print "NoSteadyError triggered: {}".format(str(e))
            phip = F + alpha*u  # phi plus
            phim = F - alpha*u # phi minus
            noSteady = 1

        for var in range(nvars):
            if var < I and max(abs(F[var,0:-1] - F[var,i]))< tol:
                Glm = F[var,i]
                Grm = F[var,i]
            else:        
                Grm = wr.wenorec(self.order, phip[var,1:-1]) # phip at i+1/2^-
                Glm = wr.wenorec(self.order, phip[var,0:-2]) # phip at i-1/2^-

            if  var < I and max(abs(F[var,1:] - F[var,i]))< tol:  
                Grp = F[var,i]
                Glp = F[var,i]
            else:
                Grp = wr.wenorec(self.order, phim[var,-1:1:-1]) # phim at i+1/2^+
                Glp = wr.wenorec(self.order, phim[var,-2:0:-1]) # phim at i-1/2^+
#            print phip[var,1:-1],  phip[var,0:-2],phim[var,-1:1:-1], phim[var,-2:0:-1]
            Gl[var] = 0.5*(Glm + Glp)
            Gr[var] = 0.5*(Grm + Grp)

        return (Gl, Gr, noSteady)

    def sw_split(self, u, x, H, eqn):
        tol = 1.e-12
        nvars = eqn.dim()
        i = (u.shape[1]-1)/2
#        alpha = np.amax(np.abs(eqn.eig_of_dF(u)))
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
        noSteady = 0
        I = nvars - np.sum(eqn.Pi(np.ones(nvars)))
        F = eqn.F(u)

        try:
            ustar = eqn.steady_constraint(H[i], u[:,i], H,x, u)
#            phip = eqn.Fp(u) - eqn.Pi(eqn.Fp(ustar))   # phi plus
#            phim = eqn.Fm(u) - eqn.Pi(eqn.Fm(ustar)) # phi minus
            phip = eqn.Fp(u) - eqn.Pi(eqn.Fp(ustar))   # phi plus
            phim = eqn.Fm(u) - eqn.Pi(eqn.Fm(ustar)) # phi minus

        except NoSteadyError, e: # no steady state exists! Default to basic WENO
            #print "NoSteadyError triggered: {}".format(str(e))
            phip = eqn.Fp(u)  # phi plus
            phim = eqn.Fm(u) # phi minus
            noSteady = 1

        for var in range(nvars):
            if var < I and max(abs(F[var,0:-1] - F[var,i]))< tol:
                Glm = F[var,i]
                Grm = F[var,i]
            else:        
                Grm = wr.wenorec(self.order, phip[var,1:-1]) # phip at i+1/2^-
                Glm = wr.wenorec(self.order, phip[var,0:-2]) # phip at i-1/2^-

            if  var < I and max(abs(F[var,1:] - F[var,i]))< tol:  
                Grp = F[var,i]
                Glp = F[var,i]
            else:
                Grp = wr.wenorec(self.order, phim[var,-1:1:-1]) # phim at i+1/2^+
                Glp = wr.wenorec(self.order, phim[var,-2:0:-1]) # phim at i-1/2^+
#            print phip[var,1:-1],  phip[var,0:-2],phim[var,-1:1:-1], phim[var,-2:0:-1]
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
        ustar = eqn.steady_constraint(H[i], u[:,i], H,x,u)
        phi = eqn.F(u) - eqn.F(ustar)

        Grm = wr.wenorec(self.order, phi[0,1:-1]) # at i+1/2^-
        Grp = wr.wenorec(self.order, phi[0,-1:1:-1]) # at i+1/2^+
        Glm = wr.wenorec(self.order, phi[0,0:-2]) # at i-1/2^-
        Glp = wr.wenorec(self.order, phi[0,-2:0:-1]) # at i-1/2^+
        Gr[0] = (critR >= 0)*Grm + (critR < 0)*Grp
        Gl[0] = (critL >= 0)*Glm + (critL < 0)*Glp
        return (Gl, Gr, 0)
    
    def rusanov_nonconservative(self, u, x, H, eqn):
        nvars = eqn.dim()
        i = (u.shape[1]-1)/2
        alpha = np.amax(np.abs(eqn.eig_of_dF(u)))
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
        ustar = eqn.steady_constraint(H[i], u[:,i], H, x, u)
        phip = eqn.F(u) - eqn.F(ustar) + alpha*(u - ustar)
        phim = eqn.F(u) - eqn.F(ustar) - alpha*(u - ustar)

        for var in range(nvars):
            Grm = wr.wenorec(self.order, phip[var,1:-1]) # phip at i+1/2^-
            Glm = wr.wenorec(self.order, phip[var,0:-2]) # phip at i-1/2^-
            Grp = wr.wenorec(self.order, phim[var,-1:1:-1]) # phim at i+1/2^+
            Glp = wr.wenorec(self.order, phim[var,-2:0:-1]) # phim at i-1/2^+
            Gl[var] = 0.5*(Glm + Glp)
            Gr[var] = 0.5*(Grm + Grp)
        return (Gl, Gr, 0)

    def rusanovg_nonconservative(self, u, x, H, alpha,eqn):
        nvars = eqn.dim()
        i = (u.shape[1]-1)/2
#        alpha = np.amax(np.abs(eqn.eig_of_dF(u)))
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
        ustar = eqn.steady_constraint(H[i], u[:,i], H,x,u)
        phip = eqn.F(u) - eqn.F(ustar) + alpha*(u - ustar)
        phim = eqn.F(u) - eqn.F(ustar) - alpha*(u - ustar)

        for var in range(nvars):
            Grm = wr.wenorec(self.order, phip[var,1:-1]) # phip at i+1/2^-
            Glm = wr.wenorec(self.order, phip[var,0:-2]) # phip at i-1/2^-

            Grp = wr.wenorec(self.order, phim[var,-1:1:-1]) # phim at i+1/2^+
            Glp = wr.wenorec(self.order, phim[var,-2:0:-1]) # phim at i-1/2^+
            Gl[var] = 0.5*(Glm + Glp)
            Gr[var] = 0.5*(Grm + Grp)
        return (Gl, Gr, 0)
    
    def sw_split_nonconservative(self, u, x, H, eqn):
        nvars = eqn.dim()
        i = (u.shape[1]-1)/2
#        alpha = np.amax(np.abs(eqn.eig_of_dF(u)))
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
        noSteady = 0
        try:
            ustar = eqn.steady_constraint(H[i], u[:,i], H,x, u)
#            phip = eqn.Fp(u) - eqn.Pi(eqn.Fp(ustar))   # phi plus
#            phim = eqn.Fm(u) - eqn.Pi(eqn.Fm(ustar)) # phi minus
            phip = eqn.Fp(u) - eqn.Fp(ustar)   # phi plus
            phim = eqn.Fm(u) - eqn.Fm(ustar) # phi minus

        except NoSteadyError, e: # no steady state exists! Default to basic WENO
            #print "NoSteadyError triggered: {}".format(str(e))
            phip = eqn.Fp(u)  # phi plus
            phim = eqn.Fm(u) # phi minus
            noSteady = 1

        for var in range(nvars):
            Grm = wr.wenorec(self.order, phip[var,1:-1]) # phip at i+1/2^-
            Glm = wr.wenorec(self.order, phip[var,0:-2]) # phip at i-1/2^-

            Grp = wr.wenorec(self.order, phim[var,-1:1:-1]) # phim at i+1/2^+
            Glp = wr.wenorec(self.order, phim[var,-2:0:-1]) # phim at i-1/2^+
#            print 'wb'
#            print phip[var,1:-1],  phip[var,0:-2],phim[var,-1:1:-1], phim[var,-2:0:-1]
#            print Grm, Glm, Grp, Glp
            Gl[var] = 0.5*(Glm + Glp)
            Gr[var] = 0.5*(Grm + Grp)
        return (Gl, Gr, noSteady)


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
        return (Gl, Gr, 0)
    
    def rusanovg_nowb(self, u, x,alpha, eqn):
        
        nvars = eqn.dim()
        i = (u.shape[1]-1)/2
  
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
        phip = eqn.F(u) + alpha*u
        phim = eqn.F(u) - alpha*u
#        print x, phip, phim
        for var in range(nvars):
            Gl[var] = .5*wr.wenorec(self.order, phim[var,::-1]) # phim at i-1/2^+
#            print phim[var,::-1], Gl[var]
            Gr[var] = .5*wr.wenorec(self.order, phip[var,:]) # phip at i+1/2^-
#            print phip[var,:], Gr[var]
        return (Gl, Gr, 0)

    def rusanov_nowb(self, u, x, eqn):
        nvars = eqn.dim()
        i = (u.shape[1]-1)/2
        alpha = np.amax(np.abs(eqn.eig_of_dF(u)))
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
        phip = eqn.F(u) + alpha*u 
        phim = eqn.F(u) - alpha*u
#        print x, phip,  phim
        for var in range(nvars):
#            print x
            Grm = wr.wenorec(self.order, phip[var,1:-1]) # phip at i+1/2^-
#            print phip[var,1:-1], Grm
            Glm = wr.wenorec(self.order, phip[var,0:-2]) # phip at i-1/2^-
#            print phip[var,0:-2], Glm

            Grp = wr.wenorec(self.order, phim[var,-1:1:-1]) # phim at i+1/2^+
#            print phim[var,-1:1:-1], Grp
            Glp = wr.wenorec(self.order, phim[var,-2:0:-1]) # phim at i-1/2^+
#            print phim[var,-2:0:-1], Glp
 
            Gl[var] = 0.5*(Glm + Glp)
            Gr[var] = 0.5*(Grm + Grp)
            
        return (Gl, Gr, 0)
    
    def sw_split_nowb(self, u, x,eqn):
        
        nvars = eqn.dim()
        i = (u.shape[1]-1)/2
  
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
        phip = eqn.Fp(u)
        phim = eqn.Fm(u)
#        print x, phip, phim
        for var in range(nvars):
            Gl[var] = .5*wr.wenorec(self.order, phim[var,::-1]) # phim at i-1/2^+
#            print phim[var,::-1], Gl[var]
            Gr[var] = .5*wr.wenorec(self.order, phip[var,:]) # phip at i+1/2^-
#            print phip[var,:], Gr[var]
        return (Gl, Gr, 0)

    
#    def sw_split_nowb(self, u, x,  eqn):
#        nvars = eqn.dim()
#        i = (u.shape[1]-1)/2
##        alpha = np.amax(np.abs(eqn.eig_of_dF(u)))
#        Gl = np.zeros(nvars)
#        Gr = np.zeros(nvars)
#        phip = eqn.Fp(u)  # phi plus
#        phim = eqn.Fm(u) # phi minus
#        for var in range(nvars):
#            Grm = wr.wenorec(self.order, phip[var,1:-1]) # phip at i+1/2^-
#            Glm = wr.wenorec(self.order, phip[var,0:-2]) # phip at i-1/2^-
#
#            Grp = wr.wenorec(self.order, phim[var,-1:1:-1]) # phim at i+1/2^+
#            Glp = wr.wenorec(self.order, phim[var,-2:0:-1]) # phim at i-1/2^+
##            print 'wb'
##            print phip[var,1:-1],  phip[var,0:-2],phim[var,-1:1:-1], phim[var,-2:0:-1]
##            print Grm, Glm, Grp, Glp
#            Gl[var] = 0.5*(Glm + Glp)
#            Gr[var] = 0.5*(Grm + Grp)
#        return (Gl, Gr, noSteady)
