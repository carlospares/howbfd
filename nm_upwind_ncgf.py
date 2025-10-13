# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import wenorec as wr
import ode_integrators as odi
import numpy as np
from nosteadyexc import NoSteadyError
from nummeth import NumericalMethod
from equation import Equation
from functionH import FunH
from howbfd_io import IoManager, parse_command_line

### Get config file from command line, or load default:
config = parse_command_line() # from howbdf_io, defaults to howbdf_config

nsteps = config.steps
multmeth = config.ode
compute_source=config.compute_source

class UpwindNCGF(NumericalMethod):
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
        #print('In tend ', xGhost)
        uGhost = np.zeros((nvars, N+2*gw)) 
        bdry.expand_with_bcs(uGhost, u, gw, eqn, initCond,funH, xGhost, tloc)  # apply BC to u
        tend = np.zeros((nvars,N))
#        fstar, bstar = self.gf(uGhost, xGhost, funH.Hx, funH.H, eqn, initCond,funH, gw, dx, tloc) #it returns the integral of the source term in the extended mesh

        #return

        fails = 0
        fail = 0
        for i in range(N):
            iOff = i+gw # i with offset for {u,x}Ghost
            iOff2 = i+nsteps
            u_st = uGhost[:,iOff-gw:iOff+gw+1] # u at the stencil for ui, size 2gw+1
#            fstar_st = fstar[:,iOff-gw:iOff+gw+1]
#            bstar_st = bstar[iOff2-gw:iOff2+gw+1]
            x_st = xGhost[  iOff-gw:iOff+gw+1] # x at the stencil for ui

#            (Gl, Gr) = self.flux(u_st, x_st, funH.H(x_st, tloc), fstar_st, bstar_st, eqn)
            (Gl, Gr) = self.flux(u_st, x_st, funH.H(x_st, tloc), eqn)
            #fails += fail
            tend[:,i] = -(Gr - Gl)/dx
            #tend[:,i] = -(Gr + Gl)/dx #written like the non conservative systems
           #print ('fails at ', tend[:,i])
            
            if fail==1:
                print ('fails at ', x[i])
                tend[:,i] += eqn.S(u[:,i])*funH.Hx(x[i], tloc)
        #exit()
        if fails>0:
            print ("{}/{} stencils failed to find a steady state solution this timestep".format(fails, N))
        return tend
    

    def flux(self, u, x, H, eqn):
        nvars = eqn.dim()
        Glm = np.zeros(nvars)
        Grp = np.zeros(nvars)
        Grm = np.zeros(nvars)
        Glp = np.zeros(nvars)

        Gl = np.zeros(nvars)
        Gr = np.zeros(nvars)
        phi = np.zeros((nvars, u.shape[-1]))
        t=0.0

        i = (u.shape[1]-1)/2
        i = int(i)
        
#        for var in range(nvars):
        dH = H[:]-H[i]
        S=np.zeros(u.shape)

        v=0
        S[nvars-1,:]=eqn.discH_jumpF(u, u[:,i],dH,v)
        phi[:,:] = eqn.F(u[:, :]) - eqn.F(u[:, [i]]) - S

        for var in range(nvars):
            Grm[var] = wr.wenorec(self.order, phi[var,1:-1]) # at i+1/2^-
            Grp[var] = wr.wenorec(self.order, phi[var,-1:1:-1]) # at i+1/2^+

        v=1
        S[nvars-1,:]=eqn.discH_jumpF(u, u[:,i],dH,v)
        phi[:,:] = eqn.F(u[:, :]) - eqn.F(u[:, [i]]) - S

        for var in range(nvars):
            Glm[var] = wr.wenorec(self.order, phi[var,0:-2]) # at i-1/2^-
            Glp[var] = wr.wenorec(self.order, phi[var,-2:0:-1]) # at i-1/2^+
#            Glm[var] = wr.wenorec(self.order, -phi[var,0:-2]) # at i-1/2^-
#            Glp[var] = wr.wenorec(self.order, -phi[var,-2:0:-1]) # at i-1/2^+


        Gr = np.dot(eqn.Piplus(u[:,i], u[:,i+1]),Grm) + np.dot(eqn.Piminus(u[:,i], u[:,i+1]),Grp)
        Gl = np.dot(eqn.Piplus(u[:,i-1], u[:,i]),Glm) + np.dot(eqn.Piminus(u[:,i-1], u[:,i]),Glp)
            
        return (Gl, Gr)
    
