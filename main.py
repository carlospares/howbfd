#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import numpy as np
import matplotlib.pyplot as plt
from initcond import InitCond
from eq_factory import equation_factory
from nm_factory import nummeth_factory
from functionH import FunH
from boundary import BoundaryCond
from timest import TimeStepping
from howbfd_io import IoManager, parse_command_line#, safe_name
from time import clock
import os


### Get config file from command line, or load default:
config = parse_command_line() # from howbdf_io, defaults to howbdf_config

#print config.a, config.b, config.T,  config.plot_exact

### Set up problem:
bdry = BoundaryCond(config)
gw = int((config.order-1)/2)+1 # number of ghost cells


# equation_factory returns an object of the appropriate subclass of Equation
eqn = equation_factory(config)
nm = nummeth_factory(config)

initCond = InitCond(eqn, config)
timest = TimeStepping(config)
io_manager = IoManager(eqn, config)

# in case it's needed, prepare data structures for order of convergence
dxs = np.zeros(config.refinements+1)
errors = np.zeros(config.refinements+1)

tini = clock()
for level in range(0, config.refinements+1):
    N = config.N * (2**level)
#    N = config.N
    print "Starting simulation with N={}...".format(N)
    interfaces = np.linspace(config.a,config.b,N+1) # we won't really use them
    x = 0.5*(interfaces[1:] + interfaces[:-1]) # midpoints (so periodic BCs are OK)
    xGhost = np.zeros(N+2*gw) # storage for x with ghost cells
    bdry.x_expand_with_bcs(xGhost, x, gw) # add BCs to x
    funH = FunH(xGhost, config)
    t = 0.0
    H = funH.H(x, t)
    u = initCond.u0(x, H) # value of u0 at midpoint of cells
    dx = x[1]-x[0]
#    print '[', 0,',', np.sum(u[0])*dx, '],'
    # Deal with possible plot at t=0
    uin = u
    io_manager.io_if_appropriate(x, u, H, t, config)

    ### Main loop
    while t < config.T:
        dt = min(config.cfl*dx/eqn.max_vel(u), io_manager.get_next_plot_time() - t)
#        dt = min(dx**(5/3.),io_manager.get_next_plot_time() - t)
        u = timest.update(x, u, nm, bdry, funH, initCond, eqn, gw, dx, dt, config, t)
        t += dt
        io_manager.io_if_appropriate(x, u, H, t, config)
        #io_manager.io_if_appropriate(x, uin-u, H, t, config)
#        print '[', t,',', np.sum(u[0])*dx, '],' 
    
    # io_manager.statistics(x, u, funH.H(x), eqn)
    exact = eqn.exact(x, t, H, config)
#    errors[level] = np.sum(np.abs( (exact[:,N/4:3*N/4] - u[:,N/4:3*N/4]) ))*dx

    errors[level] = np.sum(np.abs(exact - u))*dx
    
    
    # ^ ugly hack! Compute error only in center of domain to avoid BCs
    print "Error at N={} is {}".format(N, errors[level])
    if level > 0:
        order = (np.log(errors[level-1]) -np.log(errors[level]))/np.log(2.)
        print "Order: "+str(order)
    dxs[level] = dx
    
    io_manager.reset_timer() # otherwise only level=0 will plot
tfin = clock()
print 'CPU Time: ' + str(tfin-tini)


for i in range(N):
    print x[i],uin[0,i],u[0,i]


if config.refinements > 0:
    io_manager.plot_eoc(dxs, errors, timest.order())
    


