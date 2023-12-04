# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import wenorec as wr
import numpy as np
from nosteadyexc import NoSteadyError
from nummeth import NumericalMethod
from equation import Equation
from functionH import FunH

class Rusanov(NumericalMethod):
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
        for i in range(N):
            iOff = i+gw # i with offset for {u,x}Ghost
            u_st = uGhost[:,iOff-gw:iOff+gw+1] # u at the stencil for ui, size 2gw+1
            x_st = xGhost[  iOff-gw:iOff+gw+1] # x at the stencil for ui
            (Gl, Gr) = self.flux(u_st, x_st, eqn)
            tend[:,i] = -(Gr - Gl)/dx + eqn.S(u[:,i])*funH.Hx(x[i], tloc)
        
        if cf.funh==FunH.DISC:
            ind = np.where(x>=0)[0][0]
            Ssing=  eqn.S(u[:,ind-1])*(funH.H(x[ind], tloc)- funH.H(x[ind-1], tloc))/dx
            Ssingp = Ssing
            Ssingm = 0
            tend[:,ind-1] += Ssingm
            tend[:,ind] += Ssingp
        return tend
    
    def flux(self, u, x, eqn):
        nvars = eqn.dim()
        i = (u.shape[1]-1)/2
        alpha = np.amax(np.abs(eqn.eig_of_dF(u)))
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
        phip = eqn.F(u) + alpha*u 
        phim = eqn.F(u) - alpha*u
        for var in range(nvars):
            Grm = wr.wenorec(self.order, phip[var,1:-1]) # phip at i+1/2^-
            Glm = wr.wenorec(self.order, phip[var,0:-2]) # phip at i-1/2^-
            Grp = wr.wenorec(self.order, phim[var,-1:1:-1]) # phim at i+1/2^+
            Glp = wr.wenorec(self.order, phim[var,-2:0:-1]) # phim at i-1/2^+
            Gl[var] = 0.5*(Glm + Glp)
            Gr[var] = 0.5*(Grm + Grp)
            
        return (Gl, Gr)

 
