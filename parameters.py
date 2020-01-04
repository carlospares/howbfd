# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from numflux import Flux
from timest import TimeStepping
import importlib


class Parameters:
    # For a detailed explanation, see howbfd_config
    
    def __init__(self, filename):
        cf = importlib.import_module(filename)
        
        # try:
            # self.equation = cf.equation
        # except AttributeError:
            # print "No equation defined in {}, assuming linear transport".format(filename)
            # self.equation = Equation.LINEAR
        
        # try:
            # self.init = cf.init
        # except AttributeError:
            # print "No initial data defined in {}, assuming steady state".format(filename)
            # self.equation = InitCond.STEADY
            
        # try:
            # self.funh = cf.funh
        # except AttributeError:
            # print "No H defined in {}, assuming identity".format(filename)
            # self.funh = FunH.IDENT
        
        # try:
            # self.boundary = cf.boundary
        # except AttributeError:
            # print "No BC defined in {}, assuming periodic".format(filename)
            # self.boundary = BoundaryCond.PERIODIC
            
        # TO DO:
        # Give default values to parameters as above (instead of letting program crash)?
        # Not sure if good idea, either makes "silent" choices for parameters or spams output
        self.equation = cf.equation
        self.init = cf.init
        self.funh = cf.funh
        self.boundary = cf.boundary
        self.H_noise_factor = cf.H_noise_factor
        self.perturb_init = cf.perturb_init
        self.nummeth = cf.nummeth
        self.timest = cf.timest
        self.order = cf.order
        self.N = cf.N
        self.cfl = cf.cfl
        self.a = cf.a
        self.b = cf.b 
        self.T = cf.T
        self.plot_every = cf.plot_every
        self.show_plots = cf.show_plots
        self.save_plots = cf.save_plots
        self.save_npys = cf.save_npys
        
        self.plot_exact = cf.plot_exact
        self.refinements = 0
        