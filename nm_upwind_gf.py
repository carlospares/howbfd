# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import wenorec as wr
import ode_integrators as odi
import numpy as np
from nosteadyexc import NoSteadyError
from nummeth import NumericalMethod
from equation import Equation
from functionH import FunH

class UpwindGF(NumericalMethod):
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
        fstar = self.gf(uGhost, xGhost, funH.Hx, eqn, gw, dx) #it returns the integral of the source term in the extended mesh


        fails = 0
        fail = 0
        for i in range(N):
            iOff = i+gw # i with offset for {u,x}Ghost
            u_st = uGhost[:,iOff-gw:iOff+gw+1] # u at the stencil for ui, size 2gw+1
            fstar_st = fstar[:,iOff-gw:iOff+gw+1]
            #if i==2:
            #    print i,iOff-gw, iOff+gw+1,u_st,fstar_st
            x_st = xGhost[  iOff-gw:iOff+gw+1] # x at the stencil for ui
            (Gl, Gr) = self.flux(u_st, x_st, funH.H(x_st), fstar_st, eqn)
            #fails += fail
            tend[:,i] = -(Gr - Gl)/dx
            if fail==1:
                print 'fails at ', x[i]
                tend[:,i] += eqn.S(u[:,i])*funH.Hx(x[i])
        if fails>0:
            print "{}/{} stencils failed to find a steady state solution this timestep".format(fails, N)
        return tend
    
    def gf(self, u, x, Hx, eqn, gw, dx):
        nsteps = 4
        nvars = eqn.dim()
        N = len(x)-2*gw
        fstar = np.zeros((nvars,max(N+2*gw,N+nsteps)))

        fstar[:,0:nsteps] =  eqn.F(u[:,0]) ### initatilization of the multistep method

        for i in range(N+min(2*gw-nsteps,0)):
            iOff = nsteps + i #+max(gw,nsteps) # i with offset for {fstar}Ghost

            sumSHx=odi.odeint(nsteps,'AB', eqn, Hx, u, x, iOff)
            fstar[:,iOff] = fstar[:,iOff-1] + dx*sumSHx

        if nsteps< 2*gw :
            fstar[:,N+nsteps:N+nsteps+(2*gw-nsteps)] = fstar[:,N+nsteps-1]
    
        
        return fstar

    def flux(self, u, x, H, fstar, eqn):
        nvars = eqn.dim()
        Grm = np.zeros(nvars)
        Grp = np.zeros(nvars)
        Glm = np.zeros(nvars)
        Glp = np.zeros(nvars)
        i = (u.shape[1]-1)/2
        phi = eqn.F(u) - fstar
  
        for var in range(nvars):
            Grm[var] = wr.wenorec(self.order, phi[var,1:-1]) # at i+1/2^-
            Grp[var] = wr.wenorec(self.order, phi[var,-1:1:-1]) # at i+1/2^+
            Glm[var] = wr.wenorec(self.order, phi[var,0:-2]) # at i-1/2^-
            Glp[var] = wr.wenorec(self.order, phi[var,-2:0:-1]) # at i-1/2^+
            
        Gr = np.dot(eqn.Piplus(u[:,i], u[:,i+1]),Grm) + np.dot(eqn.Piminus(u[:,i], u[:,i+1]),Grp)
        Gl = np.dot(eqn.Piplus(u[:,i-1], u[:,i]),Glm) + np.dot(eqn.Piminus(u[:,i-1], u[:,i]),Glp)

 
        return (Gl, Gr)
    
