# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import numpy as np
# from equation import Equation

class InitCond:
    # Identifiers for initial condition
    SIN = 0
    SHOCK = 1
    STEADY = 2 # use eqn-dependent steady state
    PWPOLY = 3
    FLAT = 4

    # Identifiers for perturbation (if relevant)
    PERT_NONE = 0
    PERT_POLY = 1
    PERT_PATCH = 2
    PERT_SIN = 3
    PERT_GAUSS = 4
    PERT_MGAUSS = 5

    C = 0 # average of sine perturbation


    def __init__(self, ic, eqn, pert=PERT_NONE):
        self.initCond = ic
        self.pert = pert
        self.eqn = eqn

    def u0(self, x):
        """ Initial condition """
        U0 = np.zeros((self.eqn.dim(), len(x)))
        if self.initCond==InitCond.SHOCK:
            U0[0] = 1.0*(x <= 0.5) + 2.0 * (x>0.5)
        elif self.initCond==InitCond.SIN:
            U0[0] = 1 + np.sin(2*np.pi*x)
        elif self.initCond==InitCond.STEADY:
            U0 = self.eqn.steady(x)
        elif self.initCond==InitCond.PWPOLY:
            U0[0] = 1. + (0.13+0.05*(x-10)*(x-10))*(x>=8)*(x<=12)+0.33*((x<8)+(x>12))
        elif self.initCond==InitCond.FLAT:
            U0[0] = np.ones_like(x)
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