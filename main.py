#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import numpy as np
from pylab import *
from initcond import InitCond
from eq_factory import equation_factory
from functionH import FunH
from boundary import BoundaryCond
from numflux import Flux
from timest import TimeStepping
from howbfd_io import IoManager, parse_command_line, safe_name
import importlib
import os

### Get config file from command line, or load default:
quickConf = parse_command_line() # from howbfd_io
configfile = safe_name(quickConf[0])
cf = importlib.import_module(configfile)


### Set up problem:
bdry = BoundaryCond(cf.boundary)
gw = int((cf.order-1)/2)+1 # number of ghost cells
interfaces = np.linspace(cf.a,cf.b,cf.N+1) # we won't really use them
x = 0.5*(interfaces[1:] + interfaces[:-1]) # midpoints (so periodic BCs are OK)
xGhost = np.zeros(cf.N+2*gw) # storage for x with ghost cells
bdry.x_expand_with_bcs(xGhost, x, gw) # add BCs to x

# equation_factory returns an object of the appropriate subclass of Equation
eqn = equation_factory(cf.equation)
funH = FunH(cf.funh, x, cf.H_noise_factor)
H = funH.H(x)
initCond = InitCond(cf.init, eqn, cf.perturb_init)
flux = Flux(cf.numflux, cf.order, cf.well_balanced, cf.is_conservative)
timest = TimeStepping(cf.timest)


u = initCond.u0(x, funH.H(x)) # value of u0 at midpoint of cells
nvars = eqn.dim()
dx = x[1]-x[0]

#uGhost = np.zeros((nvars, cf.N+2*gw)) # same for u
io_manager = IoManager(cf.plot_every, cf.T, eqn)

# identifies options used to run
tag = io_manager.get_tag(cf.init, cf.perturb_init, cf.equation, cf.numflux,
                         cf.boundary, int(cf.well_balanced), cf.N, cf.order)

t = 0
# Deal with possible plot at t=0
io_manager.io_if_appropriate(x, u, H, t, show_plot=cf.show_plots, tag=tag,
                             save_plot=cf.save_plots, save_npy=cf.save_npys)

### Main loop
while t < cf.T:
    dt = min(cf.cfl*dx/eqn.max_vel(u), io_manager.get_next_plot_time() - t)
    # create expanded array for u with appropriate BCs:
    u = timest.update(x,u, flux, bdry,funH,eqn,cf.well_balanced,cf.N, gw,nvars, dx,dt)
    t += dt
    io_manager.io_if_appropriate(x, u, H, t, show_plot=cf.show_plots, tag=tag,
                                 save_plot=cf.save_plots, save_npy=cf.save_npys)

### Write some final statistics
#plot(x, initCond.u0(x, funH.H(x))[0],'k',x, eqn.steady(funH.H(x))[0],'r')
#show()
io_manager.statistics(x, u, funH.H(x), eqn)