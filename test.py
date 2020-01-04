# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 23:19:13 2019

@author: Usuario
"""

from pylab import *
g = 9.8
HC = 0.3607567072803111
hc = 1.40360489
qc = 4.67599032
#H = 0.27194503
H = HC + 3/2.*hstar - .5*qc**2/(g*hc**2) - hc 
hstar = (qc**2/g)**(1/3.)
h = linspace(1,2, 200)
phi = .5*qc**2/h**2 + g*h  - .5*qc**2/hc**2 - g*hc + g*HC
plot(h,phi, h, g*H*ones(size(h)), hstar, g*H, '*')