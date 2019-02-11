import numpy as np

class InitCond:
    # Identifiers for initial condition
    SIN = 0
    SHOCK = 1
    STEADY = 2

    # Identifiers for perturbation (if relevant)
    PERT_NONE = 0
    PERT_POLY = 1
    PERT_PATCH = 2
    PERT_SIN = 3
    PERT_GAUSS = 4

    C = 0 # average of sine perturbation


    def __init__(self, ic, pert=PERT_NONE):
        self.initCond = ic
        self.pert = pert

    def u0(self, x):
        """ Initial condition """
        if self.initCond==InitCond.SHOCK:
            U0 = 1.0*(x <= 0.5) + 2.0 * (x>0.5)
        elif self.initCond==InitCond.SIN:
            U0 = 1 + np.sin(2*np.pi*x)
        elif self.initCond==InitCond.STEADY:
            U0 = np.exp(x)
        return U0 + self.perturbation(x)

    def perturbation(self, x):
        if self.pert == InitCond.PERT_SIN:
            return 0.01*(self.C+np.sin(2*np.pi*x))
        elif self.pert == InitCond.PERT_POLY:
            return 0.01*x*(1-x)
        elif self.pert == InitCond.PERT_PATCH:
            return 1*(x>=0.6)*(x<=0.7)
        elif self.pert == InitCond.PERT_GAUSS:
            return 0.3*np.exp(-200*(x+0.5)*(x+0.5))
        else:
            return np.zeros(np.size(x))