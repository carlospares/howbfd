# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import numpy as np
# from equation import Equation

class InitCond:
    # Identifiers for initial condition
    SIN = 500
    SHOCK = 501
    STEADY = 502 # use eqn-dependent steady state
    PWPOLY = 503
    FLAT = 504
    RAREFACTION = 505
    WATER_AT_REST = 506
    ORDER_TEST = 507

    # Identifiers for perturbation (if relevant)
    PERT_NONE = 600
    PERT_POLY = 601
    PERT_PATCH = 602
    PERT_SIN = 603
    PERT_GAUSS = 604
    PERT_MGAUSS = 605

    C = 0 # average of sine perturbation


    def __init__(self, eqn, cf):
        self.initCond = cf.init
        self.pert = cf.perturb_init
        self.eqn = eqn

    def u0(self, x, H):
        """ Initial condition """
        U0 = np.zeros((self.eqn.dim(), len(x)))
        if self.initCond==InitCond.SHOCK:
            U0[0] = 2.0*(x <= 0.5) + 1.0 * (x>0.5)
        elif self.initCond==InitCond.SIN:
            U0[0] = 1 + np.sin(2*np.pi*x)
        elif self.initCond==InitCond.STEADY:
            U0 = self.eqn.steady(H)
        elif self.initCond==InitCond.PWPOLY:
            U0[0] = 1. + (0.13+0.05*(x-10)*(x-10))*(x>=8)*(x<=12)+0.33*((x<8)+(x>12))
        elif self.initCond==InitCond.FLAT:
            U0[0] = np.ones_like(x)
        elif self.initCond==InitCond.RAREFACTION:
            U0[0] = 1.0*(x <= 0) + 2.0*(x>0)
        elif self.initCond==InitCond.WATER_AT_REST:
            U0[0] = 1 + H
        elif self.initCond==InitCond.ORDER_TEST:
            U0[0] = np.exp(-1./(1 -x**2))*(x < 1)*(x > -1)
#            U0[0] = (1- x**2)*(x < 1)*(x > -1)
            
        return U0 + self.perturbation(x)

    def perturbation(self, x):
        pert = np.zeros((self.eqn.dim(), len(x)))
        if self.pert == InitCond.PERT_SIN:
            pert[0] = 0.01*(self.C+np.sin(2*np.pi*x))
        elif self.pert == InitCond.PERT_POLY:
            pert[0] = 0.01*x*(1-x)
        elif self.pert == InitCond.PERT_PATCH:
            pert[0] = 1*(x>=-0.5)*(x<=-0.3)
        elif self.pert == InitCond.PERT_GAUSS:
            pert[0] = 0.3*np.exp(-200*(x+0.5)*(x+0.5))
        elif self.pert == InitCond.PERT_MGAUSS:
            pert[0] = -0.3*np.exp(-200*x*x)
        return pert