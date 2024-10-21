# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import wenorec as wr
import numpy as np
from nosteadyexc import NoSteadyError
from nummeth import NumericalMethod
from equation import Equation
from functionH import FunH

class Upwind(NumericalMethod):
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
        fl = np.zeros((nvars,N+3))

        for i in range(N+2):
 
            u_st = uGhost[:,i:i+2*gw-1] # u at the stencil for ui, size 2gw+1
            x_st = xGhost[  i:i+2*gw-1] # x at the stencil for ui
            (Gl, Gr) = self.flux(u_st, x_st, eqn)
            mid = i + gw -1
            fl[:,i] +=  np.dot(eqn.Piminus(uGhost[:,mid-1], uGhost[:,mid]), Gl)
            fl[:,i+1] += np.dot(eqn.Piplus(uGhost[:,mid], uGhost[:,mid+1]), Gr)
        
        for i in range(N):
            tend[:,i] = -(fl[:,i+2] - fl[:,i+1])/dx+ eqn.S(u[:,i])*funH.Hx(x[i], tloc) - eqn.sigma(u[:,i])
        
            
        if cf.funh==FunH.DISC:
            ind = np.where(x>=0)[0][0]
            Ssing=  eqn.S(0.5*(u[:,ind-1] + u[:,ind]))*(funH.H(x[ind], tloc)- funH.H(x[ind-1], tloc))/dx
#            Ssing=  .5*(eqn.S(u[:,ind-1]) + eqn.S(u[:,ind]))*(funH.H(x[ind])- funH.H(x[ind-1]))/dx
#            Ssing=  eqn.S(u[:,ind-1])*(funH.H(x[ind])- funH.H(x[ind-1]))/dx
            Ssingp = np.dot(eqn.Piplus(u[:,ind-1], u[:,ind]), Ssing)
            Ssingm = np.dot(eqn.Piminus(u[:,ind-1], u[:,ind]), Ssing)

            tend[:,ind-1] += Ssingm
            tend[:,ind] += Ssingp

        return tend
    
    def flux(self, u, x, eqn):
        nvars = eqn.dim()
        i = (u.shape[1]-1)/2
        alpha = np.amax(np.abs(eqn.eig_of_dF(u)))
        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars) 
        phi = eqn.F(u)
        for var in range(nvars):
            Gl[var] = wr.wenorec(self.order, phi[var,::-1]) # phim at i-1/2^+
            Gr[var] = wr.wenorec(self.order, phi[var,:]) # phip at i+1/2^-
        return (Gl, Gr)
    
