# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import wenorec as wr
import numpy as np
from nosteadyexc import NoSteadyError
from nummeth import NumericalMethod
from equation import Equation
from functionH import FunH

class ScalarUpwindWB(NumericalMethod):
    """ 1D scalar linear transport equation with mass term
    
        u_t + alpha u_x = u

    """
    def __init__(self, cf):
        self.order = cf.order

    def tend(self, x, u, nm, bdry, funH, initCond, eqn, gw, dx, dt, cf, tloc):
        nvars = eqn.dim()
        if nvars != 1:
            print "[ERROR] Upwind (nc) only implemented for scalar equations!"
            raise NotImplementedError
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
            (Gl, Gr, fail) = self.flux(u_st, x_st, funH.H(x_st), eqn)
            fails += fail
            tend[:,i] = -(Gr - Gl)/dx
            if fail==1:
                print 'fails at ', x[i]
                tend[:,i] += eqn.S(u[:,i])*funH.Hx(x[i], tloc)
        if fails>0:
            print "{}/{} stencils failed to find a steady state solution this timestep".format(fails, N)
        return tend
    
    def flux(self, u, x, H, eqn):
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
            phi = eqn.F(u) - eqn.F(ustar)
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
