import numpy as np
from initcond import InitCond
from equation import Equation
from boundary import BoundaryCond
from numflux import Flux
from howbfd_io import IoManager

########################################
# Options                              #
########################################
equation = Equation.BURGERS            # see equation.py
init = InitCond.PWPOLY                 # see initcond.py
perturb = InitCond.PERT_NONE           # see initcond.py, PERT_NONE to omit
sw_h = Equation.SWE_H_PWPOLY           # see equation.py
boundary = BoundaryCond.FORCE_STEADY   # see boundary.py
numflux = Flux.RUSANOV                 # see numflux.py
order = 3                              # 3, 5, 7, 9, 11
well_balanced = True                   # is it well balanced? or basic WENO?
N = 100                                # number of spatial points
cfl = 0.5                              # cfl number to use for dt
a = 0                                  # left interval limit
b = 2                                  # right interval limit
T = 1                                  # end time
plot_every = 0.001                     # call io every (this many) seconds
show_plots = True                      # show plots?
save_plots = False                     # save plot images?
save_npys = False                      # save npy with solution snapshot?
########################################

bdry = BoundaryCond(boundary)
gw = int((order-1)/2)+1 # number of ghost cells
interfaces = np.linspace(a,b,N+1) # we won't really use them
x = 0.5*(interfaces[1:] + interfaces[:-1]) # midpoints (periodic BCs are well-def)
xGhost = np.zeros(N+2*gw) # storage for x with ghost cells
bdry.x_expand_with_bcs(xGhost, x, gw) # add BCs to x

eqn = Equation(equation, xGhost, sw_h)
initCond = InitCond(init, eqn, perturb)
flux = Flux(numflux, order, well_balanced)

u = initCond.u0(x) # value of u0 at midpoint of cells
nvars = eqn.dim()
dx = x[1]-x[0]

uGhost = np.zeros((nvars, N+2*gw)) # same for u
tend = np.zeros((nvars,N)) # tend[i] = (d/dt)u_i
io_manager = IoManager(plot_every, T, eqn)

# identifies options used to run
tag = io_manager.get_tag(init, perturb, equation, numflux,
                         boundary, int(well_balanced), N, order)

t = 0
while t < T:
    dt = min(cfl*dx/eqn.max_vel(u), io_manager.get_next_plot_time() - t)
    # create expanded array for u with appropriate BCs:
    bdry.expand_with_bcs(uGhost, u, gw, initCond, xGhost=xGhost)  # apply BC to u
    for i in range(N):
        iOff = i+gw # i with offset for {u,x}Ghost
        u_st = uGhost[:,iOff-gw:iOff+gw+1] # u at the stencil for ui, size 2gw+1
        x_st = xGhost[iOff-gw:iOff+gw+1] # x at the stencil for ui
        (Gl, Gr) = flux.flux(u_st, x_st, eqn, dt)
        tend[:,i] = -(Gr - Gl)/dx
        if not well_balanced:
            tend[:,i] += eqn.SHx(x[i], u[:,i])
    t += dt
    u = u + dt*tend
    io_manager.io_if_appropriate(x, u, t, show_plot=show_plots, tag=tag,
                                 save_plot=save_plots, save_npy=save_npys)

io_manager.statistics(x, u, initCond)