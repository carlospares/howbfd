# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import wenorec as wr
import numpy as np
from nosteadyexc import NoSteadyError
from nummeth import NumericalMethod
from equation import Equation
from functionH import FunH

class FluxSplit(NumericalMethod):
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
        fl = np.zeros((nvars,N+3))
        alpha = np.amax(np.abs(eqn.eig_of_dF(uGhost)))

        for i in range(N+2):
 
            u_st = uGhost[:,i:i+2*gw-1] # u at the stencil for ui, size 2gw+1
            x_st = xGhost[  i:i+2*gw-1] # x at the stencil for ui
            (Gl, Gr) = self.flux(u_st, x_st, eqn)
            fl[:,i] +=  Gl
            fl[:,i+1] += Gr  

        
        for i in range(N):
            tend[:,i] = -(fl[:,i+2] - fl[:,i+1])/dx+ eqn.S(u[:,i])*funH.Hx(x[i])
            
        if cf.funh==FunH.DISC:
            ind = np.where(x>=0)[0][0]
#            Ssing=  .5*(eqn.S(u[:,ind-1]) + eqn.S(u[:,ind]))*(funH.H(x[ind])- funH.H(x[ind-1]))/dx
            Ssing=  eqn.S(u[:,ind-1])*(funH.H(x[ind])- funH.H(x[ind-1]))/dx
            Ssingp = Ssing
            Ssingm = 0
            tend[:,ind-1] += Ssingm
            tend[:,ind] += Ssingp

        return tend
    
    def flux(self, u, x,eqn):       
        nvars = eqn.dim()
        i = (u.shape[1]-1)/2
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
        phip = eqn.Fp(u) 
        phim = eqn.Fm(u) 
        for var in range(nvars):
            Gl[var] = .5*wr.wenorec(self.order, phim[var,::-1]) # phim at i-1/2^+
            Gr[var] = .5*wr.wenorec(self.order, phip[var,:]) # phip at i+1/2^-
        return (Gl, Gr)
    