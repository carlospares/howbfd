# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import numpy as np
from scipy.optimize import *
# from equation import Equation

class InitCond:
    # Identifiers for initial condition
    READ_FROM_FILE = 0
    SIN = 500
    SHOCK = 501
    STEADY = 502 # use eqn-dependent steady state
    PWPOLY = 503
    FLAT = 504
    RAREFACTION = 505
    WATER_AT_REST = 506
    ORDER_TEST = 507
    TWO_ST=508
    WATER_MASS = 509
    TWO_ST_SW = 510
    MMSburg  = 511
#    SW_TRANS = 510

    # Identifiers for perturbation (if relevant)
    PERT_NONE = 600
    PERT_POLY = 601
    PERT_PATCH = 602
    PERT_SIN = 603
    PERT_GAUSS = 604
    PERT_MGAUSS = 605
    PERT_WB = 607
    PERT_WM = 608

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
            U0 = self.eqn.steady(H,x)
        elif self.initCond==InitCond.PWPOLY:
            U0[0] = 1. + (0.13+0.05*(x-10)*(x-10))*(x>=8)*(x<=12)+0.33*((x<8)+(x>12))
        elif self.initCond==InitCond.FLAT:
            U0[0] = np.ones_like(x)
        elif self.initCond==InitCond.RAREFACTION:
            U0[0] = 1.0*(x <= 0) + 2.0*(x>0)
        elif self.initCond==InitCond.READ_FROM_FILE:
            xx = []
            N=len(x)
            if N == 25:
                file_in = open('initial_data/analytical_sw/initial_sub_25.dat', 'r')
            elif N== 50:
                file_in = open('initial_data/analytical_sw/initial_sub_50.dat', 'r')
            elif N== 100:
                file_in = open('initial_data/analytical_sw/initial_sub_100.dat', 'r')
            elif N== 200:
                file_in = open('initial_data/analytical_sw/initial_sub_200.dat', 'r')
            elif N== 400:
                file_in = open('initial_data/analytical_sw/initial_sub_400.dat', 'r')
            elif N== 800:
                file_in = open('initial_data/analytical_sw/initial_sub_800.dat', 'r')
            for y in file_in.read().split('\n'):
                #if y.isdigit():
                xx.append(float(y))
            U0[0] = xx#2 + H
            U0[1] = U0[1]+24.0
        elif self.initCond==InitCond.WATER_AT_REST:
            xx = []
            N=len(x)
            U0[0] = 2.0 + H
            U0[1] = 0.0
        elif self.initCond==InitCond.WATER_MASS:
            U0[0]= 1 + H + 1.*(x>9)*(x < 11)
        elif self.initCond==InitCond.ORDER_TEST:
            U0[0] = (x**6*(1 - 6*(x-1) + 21*(x-1)**2 - 56*(x-1)**3 + 126*(x-1)**4-252*(x-1)**5))*( x >0)*(x < 1) + 1.*(x >= 1)
            #            U0[0] =  (x**5*(1 - 5*(x-1) + 15*(x-1)**2 - 35*(x-1)**3 + 70*(x-1)**4))*( x >0)*(x < 1) + 1.*(x >= 1)
#            U0[0] =  (x**3 -3*x**3*(x-1) + 6*x**3*(x-1)**2)*(x > 0)*(x<1)+ 1.*(x>=1.)
#            U0[0] = (-2*x**3 + 3* x**2)*(x > 0)*(x<1)+ 1.*(x>=1.)
#            U0[0] = np.exp(-1./(1 -x**2))*(x < 1)*(x > -1)
#            U0[0] = (1- x**2)*(x < 1)*(x > -1)
        elif self.initCond==InitCond.TWO_ST:
            U0[0] = 4*np.exp(x)*(x<0)+ np.exp(x)*(x>=0)
        elif self.initCond==InitCond.TWO_ST_SW:
            U0 = self.eqn.twosteady(H,x)
#            U0 = self.trans(x,H)
        elif self.initCond==InitCond.MMSburg:
            U0 = np.exp(-(x-5.0)*(x-5.0))
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
            pert[0] = 0.1*np.exp(-200*(x+0.5)*(x+0.5))
        elif self.pert == InitCond.PERT_MGAUSS:
            pert[0] = -0.3*np.exp(-200*x*x)
        elif self.pert ==InitCond.PERT_WB:
            pert[0] = .02*(x<=-.3)*(x>=-.4)
        elif self.pert == InitCond.PERT_WM:
            pert[0] = .5*(x<7.)*(x > 5.)
        return pert
    
#    def trans(self,x, H):
#        L = len(x)
#        g = 9.812
#        q = 2.5
#        C = 17.56957396120237
#        U0 = np.zeros((2,L))
#        U0[1,:] = q
#        def f(h, H): #Esta es la función que igualamos a 0 y a la que le aplicamos el método de Newton o el método de la biyección.
#            y = h + q**2/(2*g*h**2) - H  -  C/g
#            return y
#        hmin = (q**2/g)**(1/3.) #Es donde se encuentra el mínimo de la función.
#
#        tol1 = 1e-4
#        tol2 = 1e-10
#        tol3 = 1e-8
#    #Cuando Hr es próximo al Hmin, la solución es próxima a hmin.
#        for i in range(L):
#            fhmin= f(hmin, H[i])
#            if (hmin< tol2): 
#                z = (hmin, hmin)
#            elif ((abs(fhmin) > tol1)):
#                Hi = H[i]
#                try:
#                    z = (bisect(f, 1e-20, hmin,Hi), newton(f, hmin*2, args = (Hi,))) #El primer estado es el supercrítico y el segundo es el subcrítico.
#                except:
#                    print 'fallo'
#                    z = (hmin,hmin)
#            else:
#                z = (hmin, hmin)
#    
#            if x[i] <= 0:
#                U0[0,i] = z[1]
#            else:
#                U0[0,i] = z[0]
#        return U0
