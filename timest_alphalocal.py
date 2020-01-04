# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 12:10:50 2019

@author: Carlos Parés Madroñal
"""
import numpy as np
from eq_factory import equation_factory
from boundary import BoundaryCond
from numflux import Flux
from functionH import FunH



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
    
    def update(self, x,u, flux, bdry,funH,initCond,eqn, gw, dx, dt, cf):

        if self.timest == self.EULER:
            if cf.well_balanced:
                return self.euler(x, u, flux, bdry,funH,initCond, eqn, gw, dx, dt, cf)
            else:
                return self.eulernwb(x, u, flux, bdry,funH,initCond, eqn, gw, dx, dt, cf)
        if self.timest == self.TVDRK2:
            if cf.well_balanced:
                return self.tvdrk2(x, u, flux, bdry,funH,initCond, eqn, gw, dx, dt, cf)
            else:
                return self.tvdrk2nwb(x, u, flux, bdry,funH,initCond, eqn, gw, dx, dt, cf)
        if self.timest == self.TVDRK3:
            if cf.well_balanced:
                return self.tvdrk3(x, u, flux, bdry,funH,initCond, eqn, gw, dx, dt, cf)
            else:
                return self.tvdrk3nwb(x, u, flux, bdry,funH,initCond, eqn, gw, dx, dt, cf)

 
        
    def euler(self, x, u, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf):
        nvars = eqn.dim()
        N = len(x)
        tend = np.zeros((nvars,N))
        xGhost = np.zeros(N+2*gw)
        alpha = np.zeros(N+2)
        bdry.x_expand_with_bcs(xGhost, x, gw) 
        uGhost = np.zeros((nvars, N+2*gw)) 

        bdry.expand_with_bcs(uGhost, u, gw, eqn, initCond,funH, xGhost)  # apply BC to u
        for i in range(N+2):
            ust = uGhost[:,i:i+gw+1]
            alpha[i] = np.amax(np.abs(eqn.eig_of_dF(ust)))
        # import matplotlib.pyplot as plt
        # print uGhost
        # plt.title("uGhost at timest")
        # plt.plot(uGhost[0], '-*')
        # plt.show()
        fails = 0
        
        for i in range(N):
            iOff = i+gw # i with offset for {u,x}Ghost
            u_st = uGhost[:,iOff-gw:iOff+gw+1] # u at the stencil for ui, size 2gw+1
            x_st = xGhost[  iOff-gw:iOff+gw+1] # x at the stencil for ui
            alpha_st = alpha[i:i+gw+1]
            (Gl, Gr, fail) = flux.flux(u_st, x_st, funH.H(x_st), alpha_st,eqn, dt)
            fails += fail
            tend[:,i] = -(Gr - Gl)/dx

            
            if not cf.well_balanced:
                tend[:,i] += eqn.S(u[:,i])*funH.Hx(x[i])
        
        if (cf.funh==FunH.DISC and not cf.well_balanced):
            ind = np.where(x>=0)[0][0]
#            Ssing=  .5*(eqn.S(u[:,ind-1]) + eqn.S(u[:,ind]))*(funH.H(x[ind])- funH.H(x[ind-1]))/dx
            Ssing=  eqn.S(u[:,ind-1])*(funH.H(x[ind])- funH.H(x[ind-1]))/dx
            Ssingp = Ssing
            Ssingm = 0
#            Ssingp = .5*Ssing
#            Ssingm = .5*Ssing
            tend[:,ind-1] += Ssingm
            tend[:,ind] += Ssingp
#            print 'Euler'
#            print ind, x[ind-1], x[ind], funH.H(x[ind-1]), funH.H(x[ind])

        if fails>0:
            print "{}/{} stencils failed to find a steady state solution this timestep".format(fails, N)
        return u + dt*tend

    def eulernwb(self, x, u, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf):
    
        nvars = eqn.dim()
        N = len(x)
        fl = np.zeros((nvars,N+3))
        tend = np.zeros((nvars,N)) 
        xGhost = np.zeros(N+2*gw)
        bdry.x_expand_with_bcs(xGhost, x, gw) 
        uGhost = np.zeros((nvars, N+2*gw)) 
        bdry.expand_with_bcs(uGhost, u, gw, eqn, initCond,funH, xGhost)  # apply BC to u
        # import matplotlib.pyplot as plt
        # print uGhost
        # plt.title("uGhost at timest")
        # plt.plot(uGhost[0], '-*')
        # plt.show()
        fails = 0
        for i in range(N+2):
 
            u_st = uGhost[:,i:i+2*gw-1] # u at the stencil for ui, size 2gw+1
            x_st = xGhost[  i:i+2*gw-1] # x at the stencil for ui

#            print 'euler', x_st
            (Gl, Gr, fail) = flux.flux(u_st, x_st, funH.H(x_st), 0., eqn, dt)
#            print 'nwb'
#            print x_st, u_st, Gl, Gr
            fails += fail
            fl[:,i] +=  Gl
            fl[:,i+1] += Gr  
        
        
#        print fl
        for i in range(N):
            tend[:,i] = -(fl[:,i+2] - fl[:,i+1])/dx+ eqn.S(u[:,i])*funH.Hx(x[i])
            

               
        if (cf.funh==FunH.DISC and not cf.well_balanced):
            ind = np.where(x>=0)[0][0]
#            Ssing=  .5*(eqn.S(u[:,ind-1]) + eqn.S(u[:,ind]))*(funH.H(x[ind])- funH.H(x[ind-1]))/dx
            Ssing=  eqn.S(u[:,ind-1])*(funH.H(x[ind])- funH.H(x[ind-1]))/dx
            Ssingp = Ssing
            Ssingm = 0
#            Ssingp = .5*Ssing
#            Ssingm = .5*Ssing
            tend[:,ind-1] += Ssingm
            tend[:,ind] += Ssingp
#            print 'Euler'
#            print ind, x[ind-1], x[ind], funH.H(x[ind-1]), funH.H(x[ind])

        if fails>0:
            print "{}/{} stencils failed to find a steady state solution this timestep".format(fails, N)
        return u + dt*tend

    
    def tvdrk2(self, x, u, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf):
        u1 = self.euler(x, u, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf)
        u2 = self.euler(x, u1, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf)
        return .5*(u + u2)
        
    def tvdrk3(self, x, u, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf):
        u1 = self.euler(x, u, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf)
        u2 = self.euler(x, u1, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf)
        uint = 3/4.*u + 1/4.*u2
        u3 = self.euler(x, uint, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf)
        return 1/3.*u + 2/3.*u3
    def tvdrk2nwb(self, x, u, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf):
        u1 = self.eulernwb(x, u, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf)
        u2 = self.eulernwb(x, u1, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf)
        return .5*(u + u2)
        
    def tvdrk3nwb(self, x, u, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf):
        u1 = self.eulernwb(x, u, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf)
        u2 = self.eulernwb(x, u1, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf)
        uint = 3/4.*u + 1/4.*u2
        u3 = self.eulernwb(x, uint, flux, bdry, funH, initCond, eqn, gw, dx, dt, cf)
        return 1/3.*u + 2/3.*u3
        
    