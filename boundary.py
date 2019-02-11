import numpy as np

class BoundaryCond:
    # Identifiers for BCs
    PERIODIC = 0
    IN_OUT = 1 # inflow on left, homogeneous Neumann on right
    LIN_EXTRAP = 2 # linear extrapolation (for spatial domain)
    FORCE_STEADY = 3 # make BCs be values for exp(x)

    def __init__(self, bc):
        self.bc = bc

    def expand_with_bcs(self, uNew, uOld, gw, inflow=0, xGhost=0):
        """ Take an array and make a copy, augmented with BCs
        Input:
            uNew: array (assumed initialized) which will be written (size N+2*gw)
            uOld: array with values for inner cells (size N)
            gw: number of ghost cells: should be (q-1)/2 for WENO-q
            bdry: type of boundary condition (see below)
            inflow: if appropriate, value to be written at inflow boundary
        Output:
            None (uNew updated in place)
        """
        N = len(uOld)
        uNew[gw:N+gw] = uOld
        if self.bc==BoundaryCond.PERIODIC:
            uNew[:gw] = uOld[-gw:]
            uNew[-gw:] = uOld[:gw]
        elif self.bc==BoundaryCond.IN_OUT:
            uNew[:gw] = inflow
            uNew[-gw:] = uOld[-1:-1-gw:-1]
        elif self.bc==BoundaryCond.LIN_EXTRAP:
            for j in range(gw):
                uNew[j] = uOld[0] - (j+1)*(uOld[1]-uOld[0])
                uNew[N+gw+j] = uOld[-1] + (j+1)*(uOld[-1]-uOld[-2])
        elif self.bc==BoundaryCond.FORCE_STEADY:
            uNew[:gw] = np.exp(xGhost[:gw])
            uNew[-gw:] = np.exp(xGhost[:gw])