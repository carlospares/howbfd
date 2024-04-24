#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 06:37:19 2024

@author: carlos
"""

from  pylab import *
import numpy as np
import matplotlib.pyplot as plt

###################  AB3 continuous H

def AB3Dirac(a,b,N,f0, S, Hx, jumpH):
# resuelbe en [-L,L]; N tiene que ser par
    h = (b-a)/N
    x = linspace(a,b,N+1)
    f =zeros(N+1)
    beta = [23/12, -16/12, 5/12]
    f[0:3] = exp(x[0:3]**2)
    for k in range(3,N+1):
        f[k] = f[k-1] + h*(beta[0]*S(f[k-1])*Hx(x[k-1])+ beta[1]*S(f[k-2])*Hx(x[k-2])+ beta[2]*S(f[k-3])*Hx(x[k-3]))
        if abs(x[k])< h/4:
            f[k] = 0.5*(sqrt(2*f[k])*exp(jumpH))**2
    return (x,f)
    
def AM4smooth(a,b,N,f0, S, Hx, jumpH):
# resuelbe en [-L,L]; N tiene que ser par
    h = 1.0*(b-a)/N
    x = linspace(a,b,N+1)
    f =zeros(N+1)
    beta = [9./24., 19./24., -5./24., 1./24]
    f[0:3] = exp(x[0:3])
    for k in range(3,N+1):
        psi1 = 0.5*f[k-1]*f[k-1] + h*(beta[1]*S(f[k-1])*Hx(x[k-1])+ beta[2]*S(f[k-2])*Hx(x[k-2])+ beta[3]*S(f[k-3])*Hx(x[k-3]))
        den  = 0.5 - h*beta[0]*Hx(x[k])
        psi  = psi1/den
        f[k] = sqrt(psi)
    
    return (x,f)

def AB4Dirac(a,b,N,f0, S, Hx,jumpH):
    h = (b-a)/N
    x = linspace(a,b,N+1)
    f =zeros(N+1)
    beta = [55/24, -59/24, 37/24, -3/8]
    f[0:4] = exp(x[0:4]**2)
    for k in range(4,N+1):
        f[k] = f[k-1] + h*(beta[0]*S(f[k-1])*Hx(x[k-1])+ beta[1]*S(f[k-2])*Hx(x[k-2])+ beta[2]*S(f[k-3])*Hx(x[k-3]) + beta[3]*S(f[k-4])*Hx(x[k-4]))
        if abs(x[k])< h/4:
            f[k] = 0.5*(sqrt(2*f[k])*exp(jumpH))**2

    return (x,f)





def S(f):
    result = f*f
    return result

def Hx(x):
    return 1.0

def exacta(x, jumpH):
    return  exp(x)*(x<=0) + 0.5*(sqrt(2)*exp(jumpH))**2*exp(x)*(x>0)


a = -1
b = 1



mesh = [20 ,40,80,160,320, 640]

leg = ['N = '+str(N) for N in mesh]
leg.append('exacta')

figure('AM4 continuous')

print('AM4. H continuous')
print('#################')

jumpH = 0

for N in mesh:
    (x,f) = AM4smooth(a,b,N, 1, S,Hx,jumpH)
    plot(x,f)
    error = sum(abs(exacta(x,jumpH)-f))*1.0*(b-a)/N
    print('N= ',N,' error= ',error)
    if N!= mesh[0]:
        order = (log(errorold)- log(error))/log(2)
        print('order= ', order)
    errorold = error
    plot(x,f)
    legend(leg)

###################  AB3 continuous Discontinuous H




