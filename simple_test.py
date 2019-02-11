from __future__ import division
import numpy as np
import wenorec as wr
from initcond import InitCond
from equation import Equation
from boundary import BoundaryCond
import numflux as flux
import howbdf_io as wb_io

####################################
# Options                          #
####################################
init = InitCond.STEADY             # see initcond.py
perturb = InitCond.PERT_PATCH      # see initcond.py
equation = Equation.LINEAR         # see equation.py
bdry = BoundaryCond.FORCE_STEADY   # see boundary.py
order = 3                          # only 3 and 5 for now
N = 100                            # number of spatial points
M = 10000                          # number of timesteps
a = 0                              # left interval limit
b = 1                              # right interval limit
T = 10                             # end time
plot_every = 0.5                   # generate a plot at every (this many) seconds
show_plots = False                 # show plots or save only?
well_balanced = True               # is it well balanced? or basic WENO?
####################################

initCond = InitCond(init, perturb)
eqn = Equation(equation)
bdry = BoundaryCond(bdry)
xBdry = BoundaryCond(BoundaryCond.PERIODIC if bdry.bc==BoundaryCond.PERIODIC else BoundaryCond.LIN_EXTRAP)
dt = T/M  # always float division (imported from __future__)
interfaces = np.linspace(a,b,N+1) # we won't really use them
x = 0.5*(interfaces[1:] + interfaces[:-1]) # midpoints (for periodic BCs)
u = initCond.u0(x)
# print u
inflow_left = initCond.u0(a)
gw = int((order-1)/2)+1 # number of ghost cells
dx = x[1]-x[0]
# print x

xGhost = np.zeros(N+2*gw) # storage for x with ghost cells
uGhost = np.zeros(N+2*gw) # same for u
xBdry.expand_with_bcs(xGhost, x, gw) # add BCs to x
tend = np.zeros(N) # tend[i] = (d/dt)u_i

for n in range(M):
    # create expanded array for u with appropriate BCs:
    bdry.expand_with_bcs(uGhost, u, gw, inflow=inflow_left, xGhost=xGhost)  # apply BC to u
    for i in range(N):
        iOff = i+gw # i with offset for {u,x}Ghost
        u_st = uGhost[iOff-gw:iOff+gw+1] # u at the stencil for ui
        x_st = xGhost[iOff-gw:iOff+gw+1] # x at the stencil for ui
        if well_balanced:
            g_st = eqn.g(u[i], u_st, x[i], x_st) # g evaluated at the stencil
        else: # basic WENO
            g_st = eqn.F(u_st)
        Grm = wr.wenorec(order, g_st[1:-1]) # at i+1/2^-
        Grp = wr.wenorec(order, g_st[-1:1:-1]) # at i+1/2^+
        Glm = wr.wenorec(order, g_st[0:-2]) # at i-1/2^-
        Glp = wr.wenorec(order, g_st[-2:0:-1]) # at i-1/2^+
        (critL, critR) = eqn.upw_criterion(u_st)
        Gr = flux.upwind(Grm, Grp, critR)
        Gl = flux.upwind(Glm, Glp, critL)
        tend[i] = -(Gr - Gl)/dx;

        if not well_balanced: # WENO needs source terms
            tend[i] += eqn.SHx(u[i])

    u = u + dt*tend
    wb_io.plot_if_appropriate(x, u, n*dt, dt, plot_every=plot_every, show_plots=show_plots)