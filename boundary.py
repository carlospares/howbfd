# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import numpy as np

class BoundaryCond:
    # Identifiers for BCs
    PERIODIC = 400
    IN_OUT = 401 # inflow (as in IC) on left, homogeneous Neumann on right
    LIN_EXTRAP = 402 # linear extrapolation (for spatial domain)
    FORCE_STEADY = 403 # make BCs be values for u0(x)

    def __init__(self, bc):
        self.bc = bc

    def expand_with_bcs(self, uNew, uOld, gw, eqn, initCond, funH, xGhost=0):
        """ Take the array of values and make a copy, augmented with BCs
        Input:
            uNew: array (assumed initialized) which will be written (size N+2*gw)
            uOld: array with values for inner cells (size N)
            gw: number of ghost cells: should be (q+1)/2 for WENO-q
            inflow: if appropriate, value to be written at inflow boundary
            xGhost: if appropriate, x at ghost cells (for FORCE_STEADY)
        Output:
            None (uNew updated in place)
        """
        (nvars,N) = uOld.shape
        uNew[:,gw:N+gw] = uOld[:,:]
        if self.bc==BoundaryCond.PERIODIC:
            uNew[:,:gw] = uOld[:,-gw:]
            uNew[:,-gw:] = uOld[:,:gw]
        elif self.bc==BoundaryCond.IN_OUT:
            uNew[:,:gw] = initCond.u0(xGhost[:gw], funH.H(xGhost[:gw]))
            uNew[:,-gw:] = uOld[:,-1:-1-gw:-1] # naively try to make derivative zero
        elif self.bc==BoundaryCond.LIN_EXTRAP:
            for j in range(gw):
                uNew[:,j] = uOld[:,0] - (gw-j)*(uOld[:,1]-uOld[:,0])
                uNew[:,N+gw+j] = uOld[:,-1] + (j+1)*(uOld[:,-1]-uOld[:,-2])
        elif self.bc==BoundaryCond.FORCE_STEADY:
            uNew[:,:gw] = eqn.steady(funH.H(xGhost[:gw]))
            uNew[:,-gw:] = eqn.steady(funH.H(xGhost[-gw:]))
#            print uNew[:,:gw], funH.H(xGhost[:gw]), uNew[:,-gw:], funH.H(xGhost[-gw:])

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