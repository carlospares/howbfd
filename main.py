import numpy as np
from initcond import InitCond
from equation import Equation
from boundary import BoundaryCond
from numflux import Flux
from howbfd_io import IoManager, parse_command_line, safe_name
import importlib
import os

### Get config file from command line, or load default:
quickConf = parse_command_line() # from howbfd_io
if quickConf == "":
    # import howbfd_config as cf
    cf = importlib.import_module("howbfd_config")
else:
    configfile = safe_name(quickConf[0]) # from howbfd_io
    cf = importlib.import_module(configfile)

### Set up problem:
bdry = BoundaryCond(cf.boundary)
gw = int((cf.order-1)/2)+1 # number of ghost cells
interfaces = np.linspace(cf.a,cf.b,cf.N+1) # we won't really use them
x = 0.5*(interfaces[1:] + interfaces[:-1]) # midpoints (periodic BCs are well-def)
xGhost = np.zeros(cf.N+2*gw) # storage for x with ghost cells
bdry.x_expand_with_bcs(xGhost, x, gw) # add BCs to x

eqn = Equation(cf.equation, xGhost, cf.sw_H, cf.sw_H_noise_factor)
initCond = InitCond(cf.init, eqn, cf.perturb_init)
flux = Flux(cf.numflux, cf.order, cf.well_balanced)

u = initCond.u0(x) # value of u0 at midpoint of cells
nvars = eqn.dim()
dx = x[1]-x[0]

uGhost = np.zeros((nvars, cf.N+2*gw)) # same for u
tend = np.zeros((nvars,cf.N)) # tend[i] = (d/dt)u_i
io_manager = IoManager(cf.plot_every, cf.T, eqn)

# identifies options used to run
tag = io_manager.get_tag(cf.init, cf.perturb_init, cf.equation, cf.numflux,
                         cf.boundary, int(cf.well_balanced), cf.N, cf.order)

### Main loop
t = 0
while t < cf.T:
    dt = min(cf.cfl*dx/eqn.max_vel(u), io_manager.get_next_plot_time() - t)
    # create expanded array for u with appropriate BCs:
    bdry.expand_with_bcs(uGhost, u, gw, initCond, xGhost=xGhost)  # apply BC to u
    for i in range(cf.N):
        iOff = i+gw # i with offset for {u,x}Ghost
        u_st = uGhost[:,iOff-gw:iOff+gw+1] # u at the stencil for ui, size 2gw+1
        x_st = xGhost[iOff-gw:iOff+gw+1] # x at the stencil for ui
        (Gl, Gr) = flux.flux(u_st, x_st, eqn, dt)
        tend[:,i] = -(Gr - Gl)/dx
        if not cf.well_balanced:
            tend[:,i] += eqn.SHx(x[i], u[:,i])
    t += dt
    u = u + dt*tend
    io_manager.io_if_appropriate(x, u, t, show_plot=cf.show_plots, tag=tag,
                                 save_plot=cf.save_plots, save_npy=cf.save_npys)

### Write some final statistics
io_manager.statistics(x, u, initCond)