#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import numpy as np
import matplotlib.pyplot as plt
from initcond import InitCond
from eq_factory import equation_factory
from functionH import FunH
from boundary import BoundaryCond
from numflux import Flux
from timest import TimeStepping
from howbfd_io import IoManager, parse_command_line#, safe_name
import os

### Get config file from command line, or load default:
config = parse_command_line() # from howbdf_io, defaults to howbdf_config

### Set up problem:
bdry = BoundaryCond(config)
gw = int((config.order-1)/2)+1 # number of ghost cells

# equation_factory returns an object of the appropriate subclass of Equation
eqn = equation_factory(config)

initCond = InitCond(eqn, config)
flux = Flux(config)
timest = TimeStepping(config)
io_manager = IoManager(eqn, config)

# in case it's needed, prepare data structures for order of convergence
dxs = np.zeros(config.refinements+1)
errors = np.zeros(config.refinements+1)

for level in range(0, config.refinements+1):
    N = config.N * (2**level)
    print "Starting simulation with N={}...".format(N)
    interfaces = np.linspace(config.a,config.b,N+1) # we won't really use them
    x = 0.5*(interfaces[1:] + interfaces[:-1]) # midpoints (so periodic BCs are OK)
    xGhost = np.zeros(N+2*gw) # storage for x with ghost cells
    bdry.x_expand_with_bcs(xGhost, x, gw) # add BCs to x
    funH = FunH(x, config)
    H = funH.H(x)

    u = initCond.u0(x, H) # value of u0 at midpoint of cells
    dx = x[1]-x[0]

    t = 0
    # Deal with possible plot at t=0
    io_manager.io_if_appropriate(x, u, H, t, config)

    ### Main loop
    while t < config.T:
        dt = min(config.cfl*dx/eqn.max_vel(u), io_manager.get_next_plot_time() - t)
        u = timest.update(x, u, flux, bdry, funH, initCond, eqn, gw, dx, dt, config)
        t += dt
        io_manager.io_if_appropriate(x, u, H, t, config)
    
    # io_manager.statistics(x, u, funH.H(x), eqn)
    exact = eqn.exact(x, t, config)
    errors[level] = np.sum(np.abs( (exact[:,N/2:] - u[:,N/2:]) ))*dx
    # ^ ugly hack! Compute error only in right half of domain to avoid inflow from left
    print errors[level]
    dxs[level] = dx
    
    io_manager.reset_timer() # otherwise only level=0 will plot

if config.refinements > 0:
    io_manager.plot_eoc(dxs, errors, timest.order())