# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import wenorec as wr
import numpy as np
from nosteadyexc import NoSteadyError
from nummeth import NumericalMethod
from equation import Equation
from functionH import FunH

class SWUpwindWB(NumericalMethod):
    """ 1D scalar linear transport equation with mass term
    
        u_t + alpha u_x = u

    """
    def __init__(self, cf):
        self.order = cf.order

    def tend(self, x, u, nm, bdry, funH, initCond, eqn, gw, dx, dt, cf):
        nvars = eqn.dim()
        N = len(x)
        tend = np.zeros((nvars,N))
        xGhost = np.zeros(N+2*gw)
        alpha = np.zeros(N+2)
        bdry.x_expand_with_bcs(xGhost, x, gw) 
        uGhost = np.zeros((nvars, N+2*gw)) 

        bdry.expand_with_bcs(uGhost, u, gw, eqn, initCond,funH, xGhost)  # apply BC to u
        
        alpha = np.amax(np.abs(eqn.eig_of_dF(uGhost)))
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
                tend[:,i] += eqn.S(u[:,i])*funH.Hx(x[i])
        if fails>0:
            print "{}/{} stencils failed to find a steady state solution this timestep".format(fails, N)
        return tend
    
    def flux(self, u, x, H, eqn):
        g=9.812
        nvars = eqn.dim()
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
        i = (u.shape[1]-1)/2
        noSteady = 0
        try: 
            ustar = eqn.steady_constraint(H[i], u[:,i], H,x, u)
            phi = eqn.F(u) - eqn.Pi(eqn.F(ustar))
        except NoSteadyError, e: # no steady state exists! Default to basic WENO
            print "NoSteadyError triggered: {}".format(str(e))
            phi = eqn.F(u)
            noSteady = 1
        
        for var in range(nvars):
            Grm = wr.wenorec(self.order, phi[var,1:-1]) # at i+1/2^-
            Grp = wr.wenorec(self.order, phi[var,-1:1:-1]) # at i+1/2^+
            Glm = wr.wenorec(self.order, phi[var,0:-2]) # at i-1/2^-
            Glp = wr.wenorec(self.order, phi[var,-2:0:-1]) # at i-1/2^+
            
        um = .5*(u[:,i-1] + u[:,i])
        up = .5*(u[:,i] + u[:,i+1])
        lambdam1 = um[1]/um[0] - np.sqrt(g*um[0])
        lambdam2 = um[1]/um[0] + np.sqrt(g*um[0])
        lambdap1 = up[1]/up[0] - np.sqrt(g*up[0])
        lambdap2 = up[1]/up[0] + np.sqrt(g*up[0])


        if lambdap2 < 0:
            Gr[0] = Grp
        elif lambdap1 > 0:
            Gr[0] = Grm
        else:
            Gr[0] = .5*(Grp + Grm)
        if lambdam2 < 0:
            Gl[0] = Glp
        elif lambdam1 > 0:
            Gl[0] = Glm
        else:
            Gl[0] = .5*(Glp + Glm)
  
        return (Gl, Gr, noSteady)