# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import wenorec as wr
import numpy as np
from nosteadyexc import NoSteadyError
from nummeth import NumericalMethod
from equation import Equation
from functionH import FunH

class UpwindWBCons(NumericalMethod):
    """ 1D scalar linear transport equation with mass term
    
        u_t + alpha u_x = u

    """
    def __init__(self, cf):
        self.order = cf.order

    def tend(self, x, u, nm, bdry, funH, initCond, eqn, gw, dx, dt, cf, tloc):
        nvars = eqn.dim()
        N = len(x)
        xGhost = np.zeros(N+2*gw)
        bdry.x_expand_with_bcs(xGhost, x, gw) 
        uGhost = np.zeros((nvars, N+2*gw)) 
        bdry.expand_with_bcs(uGhost, u, gw, eqn, initCond,funH, xGhost, tloc)  # apply BC to u
        tend = np.zeros((nvars,N))

        fails = 0
        for i in range(N):
            iOff = i+gw # i with offset for {u,x}Ghost
            u_st = uGhost[:,iOff-gw:iOff+gw+1] # u at the stencil for ui, size 2gw+1
            x_st = xGhost[  iOff-gw:iOff+gw+1] # x at the stencil for ui
            (Gl, Gr, fail) = self.flux(u_st, x_st, funH.H(x_st, tloc), eqn)
            fails += fail
            tend[:,i] = -(Gr - Gl)/dx
            if fail==1:
                print ('fails at ', x[i])
                tend[:,i] += eqn.S(u[:,i])*funH.Hx(x[i], tloc)
        if fails>0:
            print ("{}/{} stencils failed to find a steady state solution this timestep".format(fails, N))
        return tend
    
    def flux(self, u, x, H, eqn):
        nvars = eqn.dim()
        Grm = np.zeros(nvars)
        Grp = np.zeros(nvars)
        Glm = np.zeros(nvars)
        Glp = np.zeros(nvars)
        i = (u.shape[1]-1)/2
        noSteady = 0
        try: 
            ustar = eqn.steady_constraint(H[i], u[:,i], H,x,u)
            phi = eqn.F(u) - eqn.Pi(eqn.F(ustar))
        except (NoSteadyError, e): # no steady state exists! Default to basic WENO
            print ("NoSteadyError triggered: {}".format(str(e)))
            phi = eqn.F(u)
            noSteady = 1

        for var in range(nvars):
            Grm[var] = wr.wenorec(self.order, phi[var,1:-1]) # at i+1/2^-
            Grp[var] = wr.wenorec(self.order, phi[var,-1:1:-1]) # at i+1/2^+
            Glm[var] = wr.wenorec(self.order, phi[var,0:-2]) # at i-1/2^-
            Glp[var] = wr.wenorec(self.order, phi[var,-2:0:-1]) # at i-1/2^+
            
        Gr = np.dot(eqn.Piplus(u[:,i], u[:,i+1]),Grm) + np.dot(eqn.Piminus(u[:,i], u[:,i+1]),Grp)
        Gl = np.dot(eqn.Piplus(u[:,i-1], u[:,i]),Glm) + np.dot(eqn.Piminus(u[:,i-1], u[:,i]),Glp)

 
        return (Gl, Gr, noSteady)
    
