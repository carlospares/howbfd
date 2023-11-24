# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 12:10:50 2019

@author: Carlos Parés Madroñal
"""
import numpy as np
from eq_factory import equation_factory
from boundary import BoundaryCond
#from numflux import Flux
#from functionH import FunH



class TimeStepping:
    EULER = 1
    TVDRK2 = 2
    TVDRK3 = 3
    def __init__(self, cf):
        self.timest = cf.timest
        
    def order(self):
        if self.timest == self.EULER:
            return 1.
        elif self.timest == self.TVDRK2:
            return 2.
        elif self.timest == self.TVDRK3:
            return 3.
        else:
            return None
    
    def update(self, x,u, nm, bdry,funH,initCond,eqn, gw, dx, dt, cf, t):
        if self.timest == self.EULER:
            return self.euler(x, u, nm, bdry,funH,initCond, eqn, gw, dx, dt, cf, t)
        if self.timest == self.TVDRK2:
            return self.tvdrk2(x, u, nm, bdry,funH,initCond, eqn, gw, dx, dt, cf, t)
        if self.timest == self.TVDRK3:
            return self.tvdrk3(x, u, nm, bdry,funH,initCond, eqn, gw, dx, dt, cf, t)
             

        
    def euler(self, x, u, nm, bdry, funH, initCond, eqn, gw, dx, dt, cf, t):
        return u + dt*nm.tend(x, u, nm, bdry, funH, initCond, eqn, gw, dx, dt, cf, t)

    
    def tvdrk2(self, x, u, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf, t):
        tloc = t
        u1 = self.euler(x, u, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf, tloc)
        tloc = t + 0.5*dt
        u2 = self.euler(x, u1, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf, tloc)
        return .5*(u + u2)
        
    def tvdrk3(self, x, u, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf, t):
        tloc = t
        u1 = self.euler(x, u, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf, tloc)
        tloc = t + dt
        u2 = self.euler(x, u1, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf, tloc)
        uint = 3/4.*u + 1/4.*u2
        tloc = t + 0.5*dt
        u3 = self.euler(x, uint, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf, tloc)
        return 1/3.*u + 2/3.*u3
    
