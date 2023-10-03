# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

#import numpy as np
#import sys
#import matplotlib.pyplot as plt
#from scipy.interpolate import interp1d # UnivariateSpline
#from scipy.optimize import newton
#from Funciones_salto_estacionario import phi, phiu

class NumericalMethod:


    # Identifiers for equation
    UPWIND = 100 # upwind with projection
    UPWINDWB= 101   # well-balanced upwind
    UPWINDWBCONS = 110 # well-balanced upwind conservative (only conservative for supercritical problems)
    UPWINDWB1 = 112 # upwind well-balanced for only one solution
    UPWINDGF = 113 # upwind well-balanced global flux
    RUSANOVG = 102 # standard Rusanov with alpha computed globally
    RUSANOVGWB = 103 # Rusanov WB for everty stationary solution with global alpha (no conservative)# 
    RUSANOV = 104  # standard Rusanov with alpha computed locally
    RUSANOVWB = 105 # Rusanov WB for every stationary solution with local alpha (no conservative)
    RUSANOVGWB1 = 111 #Rusanov with  global alpha that  preserves ONE stationary solution (conservative)
    RUSANOVGWBCONS = 106 # A try of Rusanov with global alpha WB for every stationary solution and conservative for the mass (doesn't work very well)
    FLUXSPLIT = 107 # Flux-splitting for the shallow water  equations
    FLUXSPLITWB = 108 # Flux-splitting well-balanced for every stationary solution (no conservative)
    FLUXSPLITWBCONS = 109 # A try of Flux-splitting WB for every stationary solution and conservative (doesn't work very well)

#    SW = 702       # eq_sw.py
#    SW_WAR = 703



    ###########################################################################
    # Following functions MUST be overridden by subclasses for anything to work:
    ###########################################################################
    
    def tend(self, x, u, nm, bdry, funH, initCond, eqn, gw, dx, dt, cf):
        """ Tendences"""
        raise NotImplementedError

    def flux(self, u, x, H, eqn, fstar, dt=0.1):
        """Numerical flux"""
        raise NotImplementedError
   
    def gf(self, u, x, Hx, eqn, gw, dx):
        """Numerical flux"""
        raise NotImplementedError
