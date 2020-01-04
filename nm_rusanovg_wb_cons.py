# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import wenorec as wr
import numpy as np
from nosteadyexc import NoSteadyError
from nummeth import NumericalMethod
from equation import Equation
from functionH import FunH

class RusanovGWBCons(NumericalMethod):
    """ 1D scalar linear transport equation with mass term
    
        u_t + alpha u_x = u

    """
    def __init__(self, cf):
        self.order = cf.order

    def tend(self, x, u, nm, bdry, funH, initCond, eqn, gw, dx, dt, cf):
        nvars = eqn.dim()
        N = len(x)
        xGhost = np.zeros(N+2*gw)
        bdry.x_expand_with_bcs(xGhost, x, gw) 
        uGhost = np.zeros((nvars, N+2*gw)) 
        bdry.expand_with_bcs(uGhost, u, gw, eqn, initCond,funH, xGhost)  # apply BC to u
        tend = np.zeros((nvars,N))
        alpha = np.amax(np.abs(eqn.eig_of_dF(uGhost)))

        fails = 0
        for i in range(N):
            iOff = i+gw # i with offset for {u,x}Ghost
            u_st = uGhost[:,iOff-gw:iOff+gw+1] # u at the stencil for ui, size 2gw+1
            x_st = xGhost[  iOff-gw:iOff+gw+1] # x at the stencil for ui
            (Gl, Gr, fail) = self.flux(u_st, x_st, funH.H(x_st), alpha, eqn)
            fails += fail
            tend[:,i] = -(Gr - Gl)/dx
            if fail==1:
                print 'fails at ', x[i]
                tend[:,i] += eqn.S(u[:,i])*funH.Hx(x[i])
        if fails>0:
            print "{}/{} stencils failed to find a steady state solution this timestep".format(fails, N)
        return tend
    
    def flux(self, u, x, H,alpha, eqn):
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
        return (Gl, Gr,  noSteady)