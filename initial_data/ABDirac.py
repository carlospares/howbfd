#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 06:37:19 2024

@author: carlos
"""

from  pylab import *

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
    return 2*f

def Hx(x):
    return x

def exacta(x, jumpH):
    return  exp(x**2)*(x<=0) + 0.5*(sqrt(2)*exp(jumpH))**2*exp(x**2)*(x>0)


a = -1
b = 1



mesh = [20,40,80,160,320, 640]

leg = ['N = '+str(N) for N in mesh]
leg.append('exacta')

figure('AB3 continuous')

print('AB3. H continuous')
print('#################')

jumpH = 0

for N in mesh:
    (x,f) = AB3Dirac(a,b,N, 1, S,Hx,jumpH)
    plot(x,f)
    error = sum(abs(exacta(x,jumpH)-f))*(b-a)/N
    print('N= ',N,' error= ',error)
    if N!= mesh[0]:
        order = (log(errorold)- log(error))/log(2)
        print('order= ', order)
    errorold = error
plot(x,exacta(x,jumpH))
legend(leg)

###################  AB3 continuous Discontinuous H








figure('AB3 discontinuous')
print('AB3. H discontinuous')
print('#################')

jumpH = 1
for N in mesh:
    (x,f) = AB3Dirac(a,b,N, 1, S,Hx,jumpH)
    plot(x,f)
    error = sum(abs(exacta(x, jumpH)-f))*(b-a)/N
    print('N= ',N,' error= ',error)
    if N!= mesh[0]:
        order = (log(errorold)- log(error))/log(2)
        print('order= ', order)
    errorold = error
plot(x,exacta(x,jumpH))
legend(leg)

###################  AB4 continuous H





figure('AB4 continuous')
print('AB4. H continuous')
print('#################')

jumpH = 0

for N in mesh:
    (x,f) = AB4Dirac(-1,1,N, 1, S,Hx,jumpH)
    plot(x,f)
    error = sum(abs(exacta(x,jumpH)-f))*2/N
    print('N= ',N,' error= ',error)
    if N!= mesh[0]:
        order = (log(errorold)- log(error))/log(2)
        print('order= ', order)
    errorold = error
plot(x,exacta(x,jumpH))
legend(leg)

###################  AB4 Discontinuous H


figure('AB4 discontinuous')
print('AB4. H discontinuous')
print('#################')

jumpH = 1

for N in mesh:
    (x,f) = AB4Dirac(a,b,N, 1, S,Hx,jumpH)
    plot(x,f)
    error = sum(abs(exacta(x,jumpH)-f))*(b-a)/N
    print('N= ',N,' error= ',error)
    if N!= mesh[0]:
        order = (log(errorold)- log(error))/log(2)
        print('order= ', order)
    errorold = error
plot(x,exacta(x,jumpH))
legend(leg)

