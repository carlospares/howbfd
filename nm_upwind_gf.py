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

class UpwindGF(NumericalMethod):
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
        fstar, bstar = self.gf(uGhost, xGhost, funH.Hx, funH.H, eqn, initCond,funH, gw, dx, tloc) #it returns the integral of the source term in the extended mesh

        #return

        fails = 0
        fail = 0
        for i in range(N):
            iOff = i+gw # i with offset for {u,x}Ghost
            iOff2 = i+nsteps
            u_st = uGhost[:,iOff-gw:iOff+gw+1] # u at the stencil for ui, size 2gw+1
            fstar_st = fstar[:,iOff-gw:iOff+gw+1]
            bstar_st = bstar[iOff2-gw:iOff2+gw+1]
            x_st = xGhost[  iOff-gw:iOff+gw+1] # x at the stencil for ui
            (Gl, Gr) = self.flux(u_st, x_st, funH.H(x_st, tloc), fstar_st, bstar_st, eqn)
            #fails += fail
            tend[:,i] = -(Gr - Gl)/dx
           #print ('fails at ', tend[:,i])
            
            if fail==1:
                print ('fails at ', x[i])
                tend[:,i] += eqn.S(u[:,i])*funH.Hx(x[i], tloc)
        if fails>0:
            print ("{}/{} stencils failed to find a steady state solution this timestep".format(fails, N))
        return tend
    
    def gf(self, u, x, Hx, H, eqn, initCond, funH, gw, dx, tloc):
        nvars = eqn.dim()
        N = len(x)-2*gw

        #fstar = np.zeros((nvars,max(N+2*gw,N+nsteps)))
        fstar = np.zeros((nvars, N+2*gw)) 

        if nsteps > gw :
            uloc = np.zeros((nvars,nsteps+N+gw))
            xloc = np.zeros((nsteps+N+gw))
            bstar = np.zeros((nsteps+N+gw))
            uloc[:,nsteps-gw:nsteps] =  u[:,0:gw] ### local extended u for the ode
            xloc[nsteps-gw:nsteps] =  x[0:gw] ### local extended u for the ode
            k=1
            for i in reversed(range(nsteps-gw)):
                #uloc[:,i] = np.exp(x[0]-k*dx)**2  #Ugly hack for convergence in steady case
                #uloc[:,i] = np.exp(x[0]-k*dx +0.1*np.sin(100*(x[0]-k*dx)))  #Ugly hack for convergence in stationary solution with oscillatory smooth H
                #uloc[:,i] = u[:,0]
                xloc[i] = x[0]-k*dx
                k +=1

            for j in range(nsteps-gw):
                uloc[:,j] = uloc[:,nsteps-gw] - (nsteps-gw-j)*(uloc[:,nsteps-gw+1]-uloc[:,nsteps-gw]) #extrapolation
                
#            uloc[:,:] = initCond.u0(xloc, funH.H(xloc, tloc))
            #return    
            uloc[:,nsteps:]=u[:,gw:]    
            xloc[nsteps:]=x[gw:]    
            #print (uloc)
        elif gw == nsteps:
            uloc = np.zeros((nvars,N+2*gw))
            xloc = np.zeros((N+2*gw))
            bstar = np.zeros((N+2*gw))
            uloc[:,:] =  u[:,:] ### initatilization of the multistep method
            xloc[:] =  x[:] ### initatilization of the multistep method
        else:
            uloc = np.zeros((nvars,N+nsteps+gw))
            xloc = np.zeros(N+nsteps+gw)
            bstar = np.zeros(N+nsteps+gw)
            iOff=gw-nsteps
            uloc[:,0:nsteps] =  u[:,gw-nsteps:nsteps+iOff] ### initatilization of the multistep method
            xloc[0:nsteps] =  x[gw-nsteps:nsteps+iOff] ### initatilization of the multistep method
            uloc[:,nsteps:] = u[:,gw:]
            xloc[nsteps:] = x[gw:]

#--------------------------ugly---only for discrete_ab 
#        if(nsteps >gw):
#            uloc = initCond.u0(xloc, funH.H(xloc, tloc))
#            uloc[:,nsteps:]=u[:,gw:]
#---------------------------------------------------------
        #fstar[:,0:nsteps] =  eqn.F(u[:,0]) ### initatilization of the multistep method
        fstar[:,0:gw] =  0#eqn.F(u[:,0:gw]) ### initatilization of the multistep method

        if nvars == 2:

            bstar[0:nsteps] =  funH.H(xloc[0:nsteps],tloc) ### initatilization of the multistep method


            for i in range(N+gw):   #-------------------this part is just for SW
                iOff = nsteps + i #+max(gw,nsteps) # i with offset for {fstar}Ghost
                sumSBx=odi.B_odeint(eqn, Hx, H, xloc, iOff, tloc)

                bstar[i+nsteps] = bstar[i+nsteps-1] + dx*sumSBx

 
        for i in range(N+gw):
            iOff = nsteps + i #+max(gw,nsteps) # i with offset for {fstar}Ghost
            sumSHx=odi.odeint(eqn, bstar, Hx, H, uloc, xloc, iOff, tloc)

            fstar[:,i+gw] = fstar[:,i+gw-1] + dx*sumSHx

        #if nsteps< 2*gw :
        #    fstar[:,N+nsteps:N+nsteps+(2*gw-nsteps)] = fstar[:,N+nsteps-1]



        return fstar, bstar

    def flux(self, u, x, H, fstar, bstar, eqn):
        nvars = eqn.dim()
        Grm = np.zeros(nvars)
        Grp = np.zeros(nvars)
        Glm = np.zeros(nvars)
        Glp = np.zeros(nvars)
        
        Urm = np.zeros(nvars)
        Urp = np.zeros(nvars)
        Ulm = np.zeros(nvars)
        Ulp = np.zeros(nvars)
        
        i = (u.shape[1]-1)/2
        i = int(i)
        
        lambda_min = u[1,i]/u[0,i] - np.sqrt(9.81*u[0,i])
        lambda_max = u[1,i]/u[0,i] - np.sqrt(9.81*u[0,i])
        sonicl = -1
        for k in range(1,2):
            if u[1,i-k]/u[0,i-k] - np.sqrt(9.81*u[0,i-k]) < lambda_min:
                lambda_min = u[1,i-k]/u[0,i-k] - np.sqrt(9.81*u[0,i-k])
            if u[1,i-k]/u[0,i-k] - np.sqrt(9.81*u[0,i-k]) > lambda_max:
                lambda_max = u[1,i-k]/u[0,i-k] - np.sqrt(9.81*u[0,i-k])
            
        if lambda_min*lambda_max < 0:
            sonicl = 1
            
        lambda_min = u[1,i]/u[0,i] - np.sqrt(9.81*u[0,i])
        lambda_max = u[1,i]/u[0,i] - np.sqrt(9.81*u[0,i])
        sonicr = -1
        for k in range(1,2):
            if u[1,i+k]/u[0,i+k] - np.sqrt(9.81*u[0,i+k]) < lambda_min:
                lambda_min = u[1,i+k]/u[0,i+k] - np.sqrt(9.81*u[0,i+k])
            if u[1,i+k]/u[0,i+k] - np.sqrt(9.81*u[0,i+k]) > lambda_max:
                lambda_max = u[1,i+k]/u[0,i+k] - np.sqrt(9.81*u[0,i+k])
            
        if lambda_min*lambda_max < 0:
            sonicr = 1

        if nvars == 2 and compute_source =='hydrostatic_reconstruction':
            phi = eqn.F_hr(u, bstar, H) - fstar
        else:
            phi = eqn.F(u) - fstar
  
        for var in range(nvars):
            Grm[var] = wr.wenorec(self.order, phi[var,1:-1]) # at i+1/2^-
            Grp[var] = wr.wenorec(self.order, phi[var,-1:1:-1]) # at i+1/2^+
            Glm[var] = wr.wenorec(self.order, phi[var,0:-2]) # at i-1/2^-
            Glp[var] = wr.wenorec(self.order, phi[var,-2:0:-1]) # at i-1/2^+
        #    if sonicr > 0 :
        #        Urm[var] = u[var,i] #wr.wenorec(self.order, u[var,1:-1]) # at i+1/2^-
        #        Urp[var] = u[var,i+1] #wr.wenorec(self.order, u[var,-1:1:-1]) # at i+1/2^+
        #    if sonicl > 0 :
        #        Ulm[var] = u[var,i-1]#wr.wenorec(self.order, u[var,0:-2]) # at i-1/2^-
        #        Ulp[var] = u[var,i] #wr.wenorec(self.order, u[var,-2:0:-1]) # at i-1/2^+
                
                
        if sonicr > 0 :
            Urm[0] = wr.wenorec(self.order, u[0,1:-1]) # at i+1/2^-
            Urp[0] = wr.wenorec(self.order, u[0,-1:1:-1]) # at i+1/2^+
        if sonicl > 0 :
            Ulm[0] = wr.wenorec(self.order, u[0,0:-2]) # at i-1/2^-
            Ulp[0] = wr.wenorec(self.order, u[0,-2:0:-1]) # at i-1/2^+
     
        if sonicr > 0 :
            Gr = Grm
        else :
            Gr = np.dot(eqn.Piplus(u[:,i], u[:,i+1]),Grm) + np.dot(eqn.Piminus(u[:,i], u[:,i+1]),Grp)
        if sonicl > 0 :
            Gl = Glm
        else :
            Gl = np.dot(eqn.Piplus(u[:,i-1], u[:,i]),Glm) + np.dot(eqn.Piminus(u[:,i-1], u[:,i]),Glp)
            
        
       # if sonicr > 0 :
       #     ar = np.sqrt( 0.5*9.81*( u[0,i] + u[0,i+1] ) ) + ( abs( u[1,i] ) + abs( u[0,i+1]  ) )/( u[0,i] + u[0,i+1] )
       #     Gr = Gr - ar*( Urp - Urm )
       # if sonicl > 0 :
       #     al = np.sqrt( 0.5*9.81*( u[0,i] + u[0,i-1] ) ) + ( abs( u[1,i] ) + abs( u[1,i-1]  ) )/( u[0,i] + u[0,i-1] )
       #     Gl = Gl - al*( Ulp - Ulm )
    
        return (Gl, Gr)
    
