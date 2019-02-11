import numpy as np
from initcond import InitCond
from equation import Equation
from boundary import BoundaryCond
from numflux import Flux
from howbfd_io import IoManager

########################################
# Options                              #
########################################
init = InitCond.STEADY                 # see initcond.py
perturb = InitCond.PERT_GAUSS          # see initcond.py, PERT_NONE to omit
equation = Equation.BURGERS            # see equation.py
boundary = BoundaryCond.FORCE_STEADY   # see boundary.py
numflux = Flux.UPWIND                  # see numflux.py
order = 3                              # 3, 5, 7, 9, 11
well_balanced = True                   # is it well balanced? or basic WENO?
N = 200                                # number of spatial points
cfl = 0.9                              # cfl number to use for dt
a = -1                                 # left interval limit
b = 1                                  # right interval limit
T = 8                                  # end time
plot_every = 0.5                       # call io every (this many) seconds
show_plots = False                     # show plots?
save_plots = False                     # save plot images?
save_npys = True                       # save npy with solution snapshot?
########################################

initCond = InitCond(init, perturb)
eqn = Equation(equation)
bdry = BoundaryCond(boundary)
flux = Flux(numflux, order, well_balanced)

interfaces = np.linspace(a,b,N+1) # we won't really use them
x = 0.5*(interfaces[1:] + interfaces[:-1]) # midpoints (for periodic BCs)
u = initCond.u0(x)
inflow_left = initCond.u0(a)
gw = int((order-1)/2)+1 # number of ghost cells
dx = x[1]-x[0]


if bdry.bc==BoundaryCond.PERIODIC:
    xBdry = BoundaryCond(BoundaryCond.PERIODIC)
else:
    xBdry = BoundaryCond(BoundaryCond.LIN_EXTRAP)
xGhost = np.zeros(N+2*gw) # storage for x with ghost cells
uGhost = np.zeros(N+2*gw) # same for u
xBdry.expand_with_bcs(xGhost, x, gw) # add BCs to x
tend = np.zeros(N) # tend[i] = (d/dt)u_i
io_manager = IoManager(plot_every, T)

# identifies options used to run
tag = "i{}-{}e{}f{}b{}w{}n{}o{}_".format(init, perturb, equation, numflux,
                                   boundary, int(well_balanced), N, order)

t = 0
while t < T:
    dt = min(cfl*dx/eqn.max_vel(u), io_manager.get_next_plot_time() - t)
    # create expanded array for u with appropriate BCs:
    bdry.expand_with_bcs(uGhost, u, gw, inflow=inflow_left, xGhost=xGhost)  # apply BC to u
    for i in range(N):
        iOff = i+gw # i with offset for {u,x}Ghost
        u_st = uGhost[iOff-gw:iOff+gw+1] # u at the stencil for ui, size 2gw+1
        x_st = xGhost[iOff-gw:iOff+gw+1] # x at the stencil for ui
        (Gl, Gr) = flux.flux(u_st, x_st, eqn)

        tend[i] = -(Gr - Gl)/dx + (1-well_balanced)*eqn.SHx(u[i]);
    t += dt
    u = u + dt*tend
    io_manager.io_if_appropriate(xGhost, uGhost, t, show_plot=show_plots,
                                 tag=tag, save_plot=False, save_npy=True)