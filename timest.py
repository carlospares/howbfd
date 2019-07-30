# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 12:10:50 2019

@author: Usuario
"""
import numpy as np
from eq_factory import equation_factory
from boundary import BoundaryCond
from numflux import Flux
from functionH import FunH



class TimeStepping:
    EULER = 0
    TVDRK2 = 1
    TVDRK3 = 2
    def __init__(self, timest):
        self.timest = timest
    
    def update(self, x,u, flux, bdry,funH,eqn,wb, N, gw,nvars, dx, dt):
        if self.timest == self.EULER:
            return self.euler(x,u, flux, bdry,funH, eqn,wb,  N, gw,nvars, dx,dt)
        if self.timest == self.TVDRK2:
            return self.tvdrk2(x,u, flux, bdry,funH, eqn,wb,  N, gw,nvars, dx,dt)
        if self.timest == self.TVDRK3:
            return self.tvdrk3(x,u, flux, bdry,funH, eqn,wb,  N, gw,nvars, dx,dt)
        
    def euler(self, x,u, flux, bdry,funH, eqn, wb, N, gw, nvars, dx,dt):
        tend = np.zeros((nvars,N)) 
        xGhost = np.zeros(N+2*gw)
        bdry.x_expand_with_bcs(xGhost, x, gw) 
        uGhost = np.zeros((nvars, N+2*gw)) 
        bdry.expand_with_bcs(uGhost, u, gw, eqn, funH,xGhost=xGhost)  # apply BC to u
        for i in range(N):
            iOff = i+gw # i with offset for {u,x}Ghost
            u_st = uGhost[:,iOff-gw:iOff+gw+1] # u at the stencil for ui, size 2gw+1
            x_st = xGhost[iOff-gw:iOff+gw+1] # x at the stencil for ui
            (Gl, Gr) = flux.flux(u_st, x_st, funH.H(x_st), eqn, dt)
            tend[:,i] = -(Gr - Gl)/dx
            if not wb:
                tend[:,i] += eqn.S(u[:,i])*funH.Hx(x[i])
        return u + dt*tend
    
    def tvdrk2(self, x,u, flux, bdry,funH,eqn, wb, N, gw, nvars, dx,dt):
        u1 = self.euler(x,u, flux, bdry,funH, eqn, wb, N, gw, nvars, dx,dt)
        u2 = self.euler(x,u1, flux, bdry,funH, eqn, wb, N, gw, nvars, dx,dt)
        return .5*(u + u2)
        
    def tvdrk3(self, x,u, flux, bdry,funH, eqn, wb, N, gw, nvars, dx,dt):
        u1 = self.euler(x,u, flux, bdry,funH, eqn, wb, N, gw, nvars, dx,dt)
        u2 = self.euler(x,u1, flux, bdry,funH, eqn, wb, N, gw, nvars, dx,dt)
        uint = 3/4.*u + 1/4.*u2
        u3 = self.euler(x,uint, flux, bdry,funH,eqn, wb, N, gw, nvars, dx,dt)
        return 1/3.*u + 2/3.*u3
        
    