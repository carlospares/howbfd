# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import numpy as np

class BoundaryCond:
    # Identifiers for BCs
    PERIODIC = 400
    IN_OUT = 401 # inflow (as in IC) on left, homogeneous Neumann on right
    LIN_EXTRAP = 402 # linear extrapolation (for spatial domain)
    FORCE_STEADY = 403 # fill left (r. right) ghost cells with steady state consistent with value at first (r. last) cell
    FORCE_STEADY_ARBITRARY = 404 # fill ghost cells with whatever initcond.STEADY says (used to be the default, but bad idea)
    FORCE_STEADY_INIT = 405 # force steady state which agrees with HConstr = H[0], uConstr = u0[:,0]
    WALL = 406
    INIT = 407
    HDOWNQUP = 408
    SUBCR = 409
    SUBCR_RE = 410

    def __init__(self, cf):
        self.bc = cf.boundary

    def expand_with_bcs(self, uNew, uOld, gw, eqn, initCond, funH, xGhost, tloc):
        """ Take the array of values and make a copy, augmented with BCs
        Input:
            uNew: array (assumed initialized) which will be written (size N+2*gw)
            uOld: array with values for inner cells (size N)
            gw: number of ghost cells: should be (q+1)/2 for WENO-q
            eqn: object of class Equation
            initCond: object of class InitCond
            funH: object of class FunH
            xGhost: array of values for x (of length N+2*gw, ie including ghost cells)
        Output:
            None (uNew updated in place)
        """
        (nvars,N) = uOld.shape
        uNew[:,gw:N+gw] = uOld[:,:]
        if self.bc==BoundaryCond.PERIODIC:
            uNew[:,:gw] = uOld[:,-gw:]
            uNew[:,-gw:] = uOld[:,:gw]
        elif self.bc==BoundaryCond.IN_OUT:
            #-----IN_OUT
##            uNew[:,:gw] = initCond.u0(xGhost[:gw], funH.H(xGhost[:gw]))
            uNew[:,:gw] =  eqn.steady(funH.H(xGhost[:gw], tloc), xGhost[:gw])
            uNew[:,-gw:] = uOld[:,-1:-1-gw:-1] # naively try to make derivative zero

        elif self.bc==BoundaryCond.SUBCR:
            #---transcritical 
            uNew[:,:gw] =  eqn.steady(funH.H(xGhost[:gw], tloc), xGhost[:gw])
            uNew[0,:gw] =  uOld[0,gw:0:-1]
##
            uNew[:,-gw:] =  eqn.steady(funH.H(xGhost[:gw], tloc), xGhost[-gw:])
            uNew[1,-gw:] = uOld[1,-1:-1-gw:-1]
            uNew[:,-gw:] = uOld[:,-1:-1-gw:-1]
            # ---subcritical
#            uNew[:,:gw] =  eqn.steady(funH.H(xGhost[:gw], tloc), xGhost[:gw])
#            uNew[0,:gw] =  uOld[0,gw:0:-1]

#            uNew[:,-gw:] =  eqn.steady(funH.H(xGhost[:gw], tloc), xGhost[-gw:])
#            uNew[1,-gw:] = uOld[1,-1:-1-gw:-1] # naively try to make derivative zero
            #-----subcritical reversed
        elif self.bc==BoundaryCond.SUBCR_RE:
            uNew[:,-gw:] =  eqn.steady(funH.H(xGhost[:gw], tloc), xGhost[:gw])
            uNew[0,-gw:] =  uOld[0,gw:0:-1]

            uNew[:,:gw] =  eqn.steady(funH.H(xGhost[:gw], tloc), xGhost[-gw:])
            uNew[1,:gw] =  uOld[1,-1:-1-gw:-1] # naively try to make derivative zero


        elif self.bc==BoundaryCond.LIN_EXTRAP:
            for j in range(gw):
                uNew[:,j] = uOld[:,0] - (gw-j)*(uOld[:,1]-uOld[:,0])
                uNew[:,N+gw+j] = uOld[:,-1] + (j+1)*(uOld[:,-1]-uOld[:,-2])
        elif self.bc==BoundaryCond.FORCE_STEADY:
            uNew[:,:gw] = eqn.steady_constraint(funH.H(xGhost[gw],tloc), uOld[:,0], funH.H(xGhost[:gw,tloc]), xGhost[:gw],np.zeros((nvars,gw)))
            uNew[:,-gw:] = eqn.steady_constraint(funH.H(xGhost[-gw-1,tloc]), uOld[:,-1], funH.H(xGhost[-gw:],tloc), xGhost[-gw:], np.zeros((nvars,gw)))
        elif self.bc==BoundaryCond.FORCE_STEADY_ARBITRARY:
            uNew[:,:gw] = eqn.steady(funH.H(xGhost[:gw],tloc), xGhost[:gw])
            uNew[:,-gw:] = eqn.steady(funH.H(xGhost[-gw:],tloc), xGhost[-gw:])
        elif self.bc==BoundaryCond.FORCE_STEADY_INIT:
            HConstr = funH.H(xGhost[gw],tloc)
            uConstr = initCond.u0(np.array([xGhost[gw]]), np.array([HConstr]))
            #import matplotlib.pyplot as plt
            #print eqn.steady_constraint(HConstr, uConstr, funH.H(xGhost))
            #plt.plot(eqn.steady_constraint(HConstr, uConstr, funH.H(xGhost))[0], '-*')
            #plt.plot(eqn.steady_constraint(HConstr, uConstr, funH.H(xGhost))[1], '-*')
            #plt.title("steady constraint")
            #plt.show()
            uNew[:,:gw] = eqn.steady_constraint(HConstr, uConstr, funH.H(xGhost[:gw],tloc),xGhost[:gw], np.zeros((nvars,gw)))
            uNew[:,-gw:] = eqn.steady_constraint(HConstr, uConstr, funH.H(xGhost[-gw:],tloc), xGhost[-gw:], np.zeros((nvars,gw)))
        elif self.bc==BoundaryCond.WALL:
#            print eqn
            uNew[0,-gw:] = uOld[0,-1:-1-gw:-1]
            uNew[0,:gw] = uOld[0, gw:0:-1]
            uNew[1,-gw:] = -uOld[1,-1:-1-gw:-1]
            uNew[1,:gw] = -uOld[1, gw:0:-1]
        elif self.bc==BoundaryCond.INIT:
            uNew[:,:gw] = initCond.u0(xGhost[:gw], funH.H(xGhost[:gw],tloc))
            uNew[:,-gw:] = initCond.u0(xGhost[-gw:], funH.H(xGhost[-gw:],tloc))
        elif self.bc == BoundaryCond.HDOWNQUP:
            uNew[0,-gw:] = uOld[0,-1:-1-gw:-1]
            uNew[0,:gw] = 2.
            uNew[1,-gw:] = 2.5
            uNew[1,:gw] = uOld[1, gw:0:-1]            

    def x_expand_with_bcs(self, xNew, xOld, gw):
        """ Take x and make a copy, augmented with BCs
        Input:
            xNew: array (assumed initialized) which will be written (size N+2*gw)
            xOld: array with values for inner cells (size N)
            gw: number of ghost cells: should be (q+1)/2 for WENO-q
        Output:
            None (xNew updated in place)
        """
        N = len(xOld)
        xNew[gw:N+gw] = xOld
        if self.bc==BoundaryCond.PERIODIC:
            xNew[:gw] = xOld[-gw:]
            xNew[-gw:] = xOld[:gw]
        else: # LIN_EXTRAP
            for j in range(gw):
                xNew[j] = xOld[0] - (gw-j)*(xOld[1]-xOld[0])
                xNew[N+gw+j] = xOld[-1] + (j+1)*(xOld[-1]-xOld[-2])
