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
#import importlib
import os

### Get config file from command line, or load default:
cf = parse_command_line() # from howbdf_io

### Set up problem:
bdry = BoundaryCond(cf.boundary)
gw = int((cf.order-1)/2)+1 # number of ghost cells

# equation_factory returns an object of the appropriate subclass of Equation
eqn = equation_factory(cf.equation)

initCond = InitCond(cf.init, eqn, cf.perturb_init)
flux = Flux(cf.numflux, cf.order, cf.well_balanced, cf.is_conservative)
timest = TimeStepping(cf.timest)
nvars = eqn.dim()
io_manager = IoManager(cf.plot_every, cf.T, eqn)

dxs = np.zeros(cf.refinements+1)
errors = np.zeros(cf.refinements+1)

for level in range(0, cf.refinements+1):
    N = cf.N * (2**level)
    print "Starting simulation with N={}...".format(N)
    interfaces = np.linspace(cf.a,cf.b,N+1) # we won't really use them
    x = 0.5*(interfaces[1:] + interfaces[:-1]) # midpoints (so periodic BCs are OK)
    xGhost = np.zeros(N+2*gw) # storage for x with ghost cells
    bdry.x_expand_with_bcs(xGhost, x, gw) # add BCs to x

    funH = FunH(cf.funh, x, cf.H_noise_factor)
    H = funH.H(x)

    u = initCond.u0(x, funH.H(x)) # value of u0 at midpoint of cells
    dx = x[1]-x[0]
    
    # print x.shape
    # print eqn.exact(x,0,cf).shape
    # plt.plot(x, eqn.exact(x,0,cf).T)
    # plt.show()

    # identifies options used to run
    tag = io_manager.get_tag(cf.init, cf.perturb_init, cf.equation, cf.numflux,
                             cf.boundary, int(cf.well_balanced), N, cf.order)

    t = 0
    # Deal with possible plot at t=0
    io_manager.io_if_appropriate(x, u, H, t, cf, show_plot=cf.show_plots, tag=tag,
                                 save_plot=cf.save_plots, save_npy=cf.save_npys)

    ### Main loop
    while t < cf.T:
        dt = min(cf.cfl*dx/eqn.max_vel(u), io_manager.get_next_plot_time() - t)
        u = timest.update(x,u, flux, bdry,funH,initCond, eqn,cf.well_balanced, N, gw,nvars, dx,dt)
        t += dt
        io_manager.io_if_appropriate(x, u, H, t, cf, show_plot=cf.show_plots, tag=tag,
                                     save_plot=cf.save_plots, save_npy=cf.save_npys)

    ### Write some final statistics
    # plt.plot(x, initCond.u0(x, funH.H(x))[0],'k',x, eqn.steady(funH.H(x))[0],'r')
    # plt.show()
    
    io_manager.statistics(x, u, funH.H(x), eqn)
    io_manager.reset_timer()

    exact = eqn.exact(x, t, cf)
    errors[level] = np.sum(np.abs( (exact[:,N/2:] - u[:,N/2:]) ))*dx
    print errors[level]
    dxs[level] = dx

if cf.refinements > 0:
    plt.close()
    z = np.polyfit(np.log(dxs), np.log(errors), 1)
    z = z[0].round(2)
    plt.loglog(dxs, errors, 'x-', label="Best fit {}".format(z))
    plt.loglog(dxs, [ (h**cf.timest) * errors[0] / (dxs[0]**cf.timest) for h in dxs], '--', label='h^{}'.format(cf.timest))
    plt.grid()
    plt.legend()

    plt.show()