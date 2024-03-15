# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 20:23:32 2019

@author: Carlos Parés Madroñal
"""
from scipy.interpolate import InterpolatedUnivariateSpline
import numpy as np
from boundary import BoundaryCond

class FunH:
    FLAT = 300
    IDENT = 301
    DISC = 302
    BUMP = 303
    DISC_BOT= 304
    IDpSIN= 305
    BUMP2 = 306
    BUMPD = 307
    BUMPS = 310
    BUMPT = 311
    SLOPE = 308
    MMSburg   = 309

    SEED = 11235813 # seed for reproducibility
    def __init__(self, xGhost, cf):
        self.funH = cf.funh
        self.noise_amplit = cf.H_noise_factor
        if self.noise_amplit != 0:
            xNeeded = xGhost[cf.gw:-cf.gw] if cf.boundary == BoundaryCond.PERIODIC else xGhost
            # unless BCs are periodic, we need to generate H also for ghost points
            np.random.seed(self.SEED) # so we get consistent results
            Hnoise = self.noise_amplit*np.random.rand(len(xNeeded))
            self.Hnoiseinterp = InterpolatedUnivariateSpline(xNeeded, Hnoise, k=1) # faster than interp1d!


    def H(self, x, t):
        """ Initial condition """
        if self.funH==self.FLAT:
            H = 0.1*np.ones_like(x)
        elif self.funH==self.IDENT:
            H = np.copy(x)
        elif self.funH==self.DISC:
            H = .1*x*(x <= 0) + (.9 +x)*(x > 0)
        elif self.funH == self.BUMP:
            H = (0.13+0.05*(x-10)*(x-10))*(x>=8)*(x<=12)+0.33*((x<8)+(x>12))
        elif self.funH==self.DISC_BOT:
            H = (1 + x)*(x <= 1) + (2+x)*(x > 1)
        elif self.funH==self.IDpSIN:
            A = 0.1
            alfa = 100#25
            H = x+ A*np.sin(alfa*x)
        elif self.funH == self.BUMP2:
            H = -.25*(1 + np.cos(5*np.pi*x))*(x<.2)*(x>-.2)
        elif self.funH == self.BUMPD:
            H = -.25*(1 + np.cos(5*np.pi*(x+1.2)))*(x<-1)*(x>-1.4) + 1.*(x > 0)
        elif self.funH == self.BUMPS:
            H = -.05*np.sin(x-12.5)*np.exp(1-(x-12.5)*(x-12.5))
        elif self.funH == self.BUMPT:
            H = -.2*np.exp( 1-1./(1.-pow(abs(x-10)/5,2)) )*(abs(x-10)<5)
            #H = -(0.2-0.05*pow((x-10),2))*(x<12)*(x>8)
        elif self.funH == self.SLOPE:
            H = x + 11
        if self.noise_amplit != 0:
            H += self.Hnoiseinterp(x)
        elif self.funH==self.MMSburg:
            H = np.exp(-(x-5.0-t)*(x-5.0-t))
        return H
    
    def Hx(self, x, t):
        """ Initial condition """
        #if self.noise_amplit != 0:
        #    print "[ERROR] Tried to compute H_x but H has noise on, not smooth"
        #    raise NotImplementedError
            
        if self.funH==FunH.FLAT:
            Hx = np.zeros_like(x)
        elif self.funH==FunH.IDENT:
            Hx = np.ones_like(x)
        elif self.funH==FunH.DISC:
            Hx = .1*(x <= 0) + 1.*(x > 0)
        elif self.funH == self.BUMP:
            Hx= 0.1*(x-10)*(x>=8)*(x<=12)
        elif self.funH==self.DISC_BOT:
            Hx = np.ones_like(x)
        elif self.funH == self.IDpSIN:
            A = 0.1
            alfa = 100#25
            Hx = np.ones_like(x) + A*alfa*np.cos(alfa*x)
        elif self.funH == self.BUMP2:
            Hx =  1.25*np.pi*np.sin(5*np.pi*x)*(x<.2)*(x>-.2)
        elif self.funH == self.SLOPE:
            Hx = np.ones_like(x)
        elif self.funH == self.BUMPD:
            Hx =  1.25*np.pi*np.sin(5*np.pi*(x+1.2))*(x<-1)*(x>-1.4)
        elif self.funH == self.BUMPS:
            Hx =-0.05*np.cos(x-12.5)*np.exp(1-(x-12.5)*(x-12.5))+0.1*np.sin(x-12.5)*(x-12.5)*np.exp(1-(x-12.5)*(x-12.5))
        elif self.funH == self.MMSburg:
            Hx =  -2.0*(x-5-t)*np.exp(-(x-5.0-t)*(x-5.0-t))        
        return Hx

