#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import numpy as np
from equation import Equation
from functionH import FunH



def odeint(nsteps, multmeth, eqn, arg0, arg1, arg2, arg3, arg4, arg5):
    nvars = eqn.dim()
    if nsteps == 4 and multmeth == 'AM':
        if nvars == 1:
            return adamsmoulton4(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
        elif nvars ==2:
            return adamsmoulton4SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
    elif nsteps == 4 and multmeth == 'AB':
        if nvars == 1:
            return adamsbashforth4(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
        elif nvars ==2:
            return adamsbashforth4SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
    elif nsteps == 6 and multmeth == 'AM':
        if nvars == 1:
            return adamsmoulton6(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
        elif nvars ==2:
            return adamsmoulton6SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
    elif nsteps == 6 and multmeth == 'AB':
        if nvars == 1:
            return adamsbashforth6(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
        elif nvars ==2:
            return adamsbashforth6SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
    elif nsteps == 8 and multmeth == 'AM' :
        if nvars == 1:
            return adamsmoulton8(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
        elif nvars ==2:
            return adamsmoulton8SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
    elif nsteps == 8 and multmeth == 'AB' :
        if nvars == 1:
            return adamsbashforth8(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
        elif nvars ==2:
            return adamsbashforth8SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
    elif nsteps == 2 and multmeth == 'AM':
        if nvars == 1:
            return adamsmoulton2(eqn, arg0, arg1, arg2, arg3, arg4, arg5)# MARIO!!!!
        elif nvars ==2:
            return adamsmoulton2SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5)# MARIO!!!!
    elif nsteps == 3 and multmeth == 'AM':
        if nvars == 1:
            return adamsmoulton3(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
        elif nvars ==2:
            return adamsmoulton3SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
    elif nsteps == 2 and multmeth == 'AB':
        if nvars == 1:
            return adamsbashforth2(eqn, arg0, arg1, arg2, arg3, arg4, arg5)# MARIO!!!!
        elif nvars ==2:
            return adamsbashforth2SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5)# MARIO!!!!
    elif nsteps == 3 and multmeth == 'AB':
        if nvars == 1:
            return adamsbashforth3(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
        elif nvars ==2:
            return adamsbashforth3SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
        
def adamsbashforth2(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 2
    ab_coeff=[-1./2., 3./2]
    
    sumSHx = np.zeros(nvars)
    for j in [-2, -1]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps]*eqn.S(u[:,i+j])*Hx(x[i+j],t)
    return sumSHx
    
def adamsbashforth2SW(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 2
    ab_coeff=[-1./2., 3./2]
    
    ddx = x[i] - x[i-1]
    g = 9.812
    
    # Collect:
    # - stencil nodes
    # - eta values (REMARK: really specific to gravity source !!!!)
    # - bathymetry and bathymetry derivatives values
    xx  = np.zeros(nsteps)
    eta = np.zeros(nsteps)
    bb  = np.zeros(nsteps)
    for l in range(0,nsteps):
        xx[l] = x[i-nsteps+l+1]
        eta[l] = -H(x[i-nsteps+l],t)+u[0,i-nsteps+l]
        bb[l] = H(x[i-nsteps+l+1],t)
    
    
    Bx = np.zeros(nsteps)
    for q in range(0,nsteps):
        Bx[q] = 0.0
        LL = Lprime( nsteps, xx, x[i-nsteps+q] )
        for p in range(0,nsteps):
            Bx[q] = Bx[q] + LL[p]*bb[p]
        
    sumSHx = np.zeros(nvars)
    sumSHx[1] = 0.5*g*( bb[nsteps-1]*bb[nsteps-1] - bb[nsteps-2]*bb[nsteps-2]  )/ddx
    for j in [-2, -1]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps]*g*eta[j+nsteps]*Bx[j+nsteps]
        
    return sumSHx

def adamsbashforth3(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 3
    ab_coeff=[5./12., -16./12., 23./12.]

    sumSHx = np.zeros(nvars)
    for j in [-3, -2, -1]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps]*eqn.S(u[:,i+j])*Hx(x[i+j],t)

    return sumSHx
    
    
def adamsbashforth3SW(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 3
    ab_coeff=[5./12., -16./12., 23./12.]

    ddx = x[i] - x[i-1]
    g = 9.812
    
    # Collect:
    # - stencil nodes
    # - eta values (REMARK: really specific to gravity source !!!!)
    # - bathymetry and bathymetry derivatives values
    xx  = np.zeros(nsteps)
    eta = np.zeros(nsteps)
    bb  = np.zeros(nsteps)
    for l in range(0,nsteps):
        xx[l] = x[i-nsteps+l+1]
        eta[l] = -H(x[i-nsteps+l],t)+u[0,i-nsteps+l]
        bb[l] = H(x[i-nsteps+l+1],t)
    
    Bx = np.zeros(nsteps)
    for q in range(0,nsteps):
        Bx[q] = 0.0
        LL = Lprime( nsteps, xx, x[i-nsteps+q] )
        for p in range(0,nsteps):
            Bx[q] = Bx[q] + LL[p]*bb[p]
      
    # Compute integrated source
    sumSHx = np.zeros(nvars)
    sumSHx[1] = 0.5*g*( bb[nsteps-1]*bb[nsteps-1] - bb[nsteps-2]*bb[nsteps-2]  )/ddx
    for j in [-3, -2, -1]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps]*g*eta[j+nsteps]*Bx[j+nsteps]

    return sumSHx

def adamsbashforth4(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 4
    ab_coeff=[-9./24., 37./24., -59./24., 55./24.]

    sumSHx = np.zeros(nvars)
    for j in [-4, -3, -2, -1]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps]*eqn.S(u[:,i+j])*Hx(x[i+j],t)

    return sumSHx

def adamsbashforth4SW(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 4
    ab_coeff=[-9./24., 37./24., -59./24., 55./24.]

    ddx = x[i] - x[i-1]
    g = 9.812

    # Collect:
    # - stencil nodes
    # - eta values (REMARK: really specific to gravity source !!!!)
    # - bathymetry and bathymetry derivatives values
    xx  = np.zeros(nsteps)
    eta = np.zeros(nsteps)
    bb  = np.zeros(nsteps)
    for l in range(0,nsteps):
        xx[l] = x[i-nsteps+l+1]
        eta[l] = -H(x[i-nsteps+l],t)+u[0,i-nsteps+l]
        bb[l] = H(x[i-nsteps+l+1],t)


    Bx = np.zeros(nsteps)
    for q in range(0,nsteps):
        Bx[q] = 0.0
        LL = Lprime( nsteps, xx, x[i-nsteps+q] )
        for p in range(0,nsteps):
            Bx[q] = Bx[q] + LL[p]*bb[p]

    sumSHx = np.zeros(nvars)
    sumSHx[1] = 0.5*g*( bb[nsteps-1]*bb[nsteps-1] - bb[nsteps-2]*bb[nsteps-2]  )/ddx
    for j in [-4, -3, -2, -1]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps]*g*eta[j+nsteps]*Bx[j+nsteps]

    return sumSHx

def adamsbashforth6(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 6
    ab_coeff=[-475./1440., 2877./1440., -7298./1440., 9982./1440,  -7923./1440., 4277./1440.]

    sumSHx = np.zeros(nvars)
    for j in [-6, -5, -4, -3, -2, -1]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps]*eqn.S(u[:,i+j])*Hx(x[i+j], t)

    return sumSHx

def adamsbashforth6SW(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 6
    ab_coeff=[-475./1440., 2877./1440., -7298./1440., 9982./1440,  -7923./1440., 4277./1440.]

    ddx = x[i] - x[i-1]
    g = 9.812

    # Collect:
    # - stencil nodes
    # - eta values (REMARK: really specific to gravity source !!!!)
    # - bathymetry and bathymetry derivatives values
    xx  = np.zeros(nsteps)
    eta = np.zeros(nsteps)
    bb  = np.zeros(nsteps)
    for l in range(0,nsteps):
        xx[l] = x[i-nsteps+l+1]
        eta[l] = -H(x[i-nsteps+l],t)+u[0,i-nsteps+l]
        bb[l] = H(x[i-nsteps+l+1],t)


    Bx = np.zeros(nsteps)
    for q in range(0,nsteps):
        Bx[q] = 0.0
        LL = Lprime( nsteps, xx, x[i-nsteps+q] )
        for p in range(0,nsteps):
            Bx[q] = Bx[q] + LL[p]*bb[p]

    sumSHx = np.zeros(nvars)
    sumSHx[1] = 0.5*g*( bb[nsteps-1]*bb[nsteps-1] - bb[nsteps-2]*bb[nsteps-2]  )/ddx
    for j in [-6, -5, -4, -3, -2, -1]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps]*g*eta[j+nsteps]*Bx[j+nsteps]

    return sumSHx

def adamsbashforth8(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 8
    ab_coeff=[-36799./120960., 295767./120960., -1041723./120960., 2102243./120960.,  -2664477./120960., 2183877./120960., -1152169./120960., 434241./120960.]

    sumSHx = np.zeros(nvars)
    for j in [-8, -7, -6, -5, -4, -3, -2, -1]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps]*eqn.S(u[:,i+j])*Hx(x[i+j], t)

    return sumSHx

def adamsbashforth8SW(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 8
    ab_coeff=[-36799./120960., 295767./120960., -1041723./120960., 2102243./120960.,  -2664477./120960., 2183877./120960., -1152169./120960., 434241./120960.]

    ddx = x[i] - x[i-1]
    g = 9.812

    # Collect:
    # - stencil nodes
    # - eta values (REMARK: really specific to gravity source !!!!)
    # - bathymetry and bathymetry derivatives values
    xx  = np.zeros(nsteps)
    eta = np.zeros(nsteps)
    bb  = np.zeros(nsteps)
    for l in range(0,nsteps):
        xx[l] = x[i-nsteps+l+1]
        eta[l] = -H(x[i-nsteps+l],t)+u[0,i-nsteps+l]
        bb[l] = H(x[i-nsteps+l+1],t)


    Bx = np.zeros(nsteps)
    for q in range(0,nsteps):
        Bx[q] = 0.0
        LL = Lprime( nsteps, xx, x[i-nsteps+q] )
        for p in range(0,nsteps):
            Bx[q] = Bx[q] + LL[p]*bb[p]

    sumSHx = np.zeros(nvars)
    sumSHx[1] = 0.5*g*( bb[nsteps-1]*bb[nsteps-1] - bb[nsteps-2]*bb[nsteps-2]  )/ddx
    for j in [-8, -7, -6, -5, -4, -3, -2, -1]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps]*g*eta[j+nsteps]*Bx[j+nsteps]

    return sumSHx

def adamsmoulton2(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 2
    ab_coeff=[1./2., 1./2.]

    sumSHx = np.zeros(nvars)
    for j in [-1, 0]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps-1]*eqn.S(u[:,i+j])*Hx(x[i+j], t)

    return sumSHx

def adamsmoulton2SW(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 2
    ab_coeff=[1./2., 1./2.]

    # Collect stencil nodes
    xx = [ x[i-1], x[i] ]
    ddx = x[i] - x[i-1]
    g = 9.812
    # Collect eta values (REMARK: really specific to gravity source !!!!)
    eta = [ -H(x[i-1],t)+u[0,i-1],  -H(x[i],t)+u[0,i] ]
    
    # Collect bathymetry and bathymetry derivatives values
    bb = [ H(x[i-1],t), H(x[i],t)]
    
    Bx = np.zeros(nsteps)
    for q in range(0,nsteps):
        Bx[q] = 0.0
        LL = Lprime( nsteps, xx, x[i-nsteps+q+1] )
        for p in range(0,nsteps):
            Bx[q] = Bx[q] + LL[p]*bb[p]
        
    sumSHx = np.zeros(nvars)
    sumSHx[nvars-1] = 0.5*g*( bb[nsteps-1]*bb[nsteps-1] - bb[nsteps-2]*bb[nsteps-2]  )/ddx
    for j in [-2, -1]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps]*g*eta[j+nsteps]*Bx[j+nsteps]

    return sumSHx

    

def adamsmoulton3(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 3
    ab_coeff=[-1./12., 8./12., 5./12.]

    sumSHx = np.zeros(nvars)
    for j in [-2, -1, 0]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps-1]*eqn.S(u[:,i+j])*Hx(x[i+j], t)

    return sumSHx


def adamsmoulton3SW(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 3
    ab_coeff=[-1./12., 8./12., 5./12.]
    
    ddx = x[i] - x[i-1]
    g = 9.812

    # Collect:
    # - stencil nodes
    # - eta values (REMARK: really specific to gravity source !!!!)
    # - bathymetry and bathymetry derivatives values
    xx  = np.zeros(nsteps)
    eta = np.zeros(nsteps)
    bb  = np.zeros(nsteps)
    for l in range(0,nsteps):
        xx[l] = x[i-nsteps+l+1]
        eta[l] = -H(x[i-nsteps+l+1],t)+u[0,i-nsteps+l+1]
        bb[l] = H(x[i-nsteps+l+1],t)
    
    Bx = np.zeros(nsteps)
    for q in range(0,nsteps):
        Bx[q] = 0.0
        LL = Lprime( nsteps, xx, x[i-nsteps+q+1] )
        for p in range(0,nsteps):
            Bx[q] = Bx[q] + LL[p]*bb[p]
      
    # Compute integrated source
    
    sumSHx = np.zeros(nvars)
    sumSHx[nvars-1] = 0.5*g*( bb[nsteps-1]*bb[nsteps-1] - bb[nsteps-2]*bb[nsteps-2]  )/ddx
    for j in [-2, -1, 0]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps-1]*g*eta[j+nsteps-1]*Bx[j+nsteps-1]
        

    return sumSHx

def adamsmoulton4(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 4
    ab_coeff=[1./24., -5./24., 19./24., 9./24]

    sumSHx = np.zeros(nvars)
    for j in [-3, -2, -1, 0]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps-1]*eqn.S(u[:,i+j])*Hx(x[i+j], t)

    return sumSHx

def adamsmoulton4SW(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 4
    ab_coeff=[1./24., -5./24., 19./24., 9./24]

    ddx = x[i] - x[i-1]
    g = 9.812

    # Collect:
    # - stencil nodes
    # - eta values (REMARK: really specific to gravity source !!!!)
    # - bathymetry and bathymetry derivatives values
    xx  = np.zeros(nsteps)
    eta = np.zeros(nsteps)
    bb  = np.zeros(nsteps)
    for l in range(0,nsteps):
        xx[l] = x[i-nsteps+l+1]
        eta[l] = -H(x[i-nsteps+l+1],t)+u[0,i-nsteps+l+1]
        bb[l] = H(x[i-nsteps+l+1],t)

    Bx = np.zeros(nsteps)
    for q in range(0,nsteps):
        Bx[q] = 0.0
        LL = Lprime( nsteps, xx, x[i-nsteps+q+1] )
        for p in range(0,nsteps):
            Bx[q] = Bx[q] + LL[p]*bb[p]

    # Compute integrated source

    sumSHx = np.zeros(nvars)
    sumSHx[nvars-1] = 0.5*g*( bb[nsteps-1]*bb[nsteps-1] - bb[nsteps-2]*bb[nsteps-2]  )/ddx
    for j in [-3, -2, -1, 0]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps-1]*g*eta[j+nsteps-1]*Bx[j+nsteps-1]

    return sumSHx

def adamsmoulton6(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 6
    ab_coeff=[27./1440., -173./1440., 482./1440., -798./1440,  1427./1440., 475./1440.]

    sumSHx = np.zeros(nvars)
    for j in [-5, -4, -3, -2, -1, 0]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps-1]*eqn.S(u[:,i+j])*Hx(x[i+j], t)

    return sumSHx

def adamsmoulton6SW(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 6
    ab_coeff=[27./1440., -173./1440., 482./1440., -798./1440,  1427./1440., 475./1440.]


    ddx = x[i] - x[i-1]
    g = 9.812

    # Collect:
    # - stencil nodes
    # - eta values (REMARK: really specific to gravity source !!!!)
    # - bathymetry and bathymetry derivatives values
    xx  = np.zeros(nsteps)
    eta = np.zeros(nsteps)
    bb  = np.zeros(nsteps)
    for l in range(0,nsteps):
        xx[l] = x[i-nsteps+l+1]
        eta[l] = -H(x[i-nsteps+l+1],t)+u[0,i-nsteps+l+1]
        bb[l] = H(x[i-nsteps+l+1],t)

    Bx = np.zeros(nsteps)
    for q in range(0,nsteps):
        Bx[q] = 0.0
        LL = Lprime( nsteps, xx, x[i-nsteps+q+1] )
        for p in range(0,nsteps):
            Bx[q] = Bx[q] + LL[p]*bb[p]

    # Compute integrated source

    sumSHx = np.zeros(nvars)
    sumSHx[nvars-1] = 0.5*g*( bb[nsteps-1]*bb[nsteps-1] - bb[nsteps-2]*bb[nsteps-2]  )/ddx
    for j in [-5, -4, -3, -2, -1, 0]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps-1]*g*eta[j+nsteps-1]*Bx[j+nsteps-1]

    return sumSHx

def adamsmoulton8(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 8
    ab_coeff=[1375./120960., -11351./120960., 41499./120960.,  -88547./120960., 123133./120960., -121797./120960, 139849./120960., 36799/120960.]

    sumSHx = np.zeros(nvars)
    for j in [-7, -6, -5, -4, -3, -2, -1, 0]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps-1]*eqn.S(u[:,i+j])*Hx(x[i+j], t)

    return sumSHx

def adamsmoulton8SW(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 8
    ab_coeff=[1375./120960., -11351./120960., 41499./120960.,  -88547./120960., 123133./120960., -121797./120960, 139849./120960., 36799/120960.]
 
    ddx = x[i] - x[i-1]
    g = 9.812

    # Collect:
    # - stencil nodes
    # - eta values (REMARK: really specific to gravity source !!!!)
    # - bathymetry and bathymetry derivatives values
    xx  = np.zeros(nsteps)
    eta = np.zeros(nsteps)
    bb  = np.zeros(nsteps)
    for l in range(0,nsteps):
        xx[l] = x[i-nsteps+l+1]
        eta[l] = -H(x[i-nsteps+l+1],t)+u[0,i-nsteps+l+1]
        bb[l] = H(x[i-nsteps+l+1],t)

    Bx = np.zeros(nsteps)
    for q in range(0,nsteps):
        Bx[q] = 0.0
        LL = Lprime( nsteps, xx, x[i-nsteps+q+1] )
        for p in range(0,nsteps):
            Bx[q] = Bx[q] + LL[p]*bb[p]

    # Compute integrated source

    sumSHx = np.zeros(nvars)
    sumSHx[nvars-1] = 0.5*g*( bb[nsteps-1]*bb[nsteps-1] - bb[nsteps-2]*bb[nsteps-2]  )/ddx
    for j in [-7, -6, -5, -4, -3, -2, -1, 0]:
        sumSHx[nvars-1] += ab_coeff[j+nsteps-1]*g*eta[j+nsteps-1]*Bx[j+nsteps-1]

    return sumSHx

def Lprime(m, xx, y):

    LL = np.zeros(m)
    
    for l in range(0, m):
        den = 1.0
        for j in range(0, m):
            if j!=l:
                den = den*(xx[l]-xx[j])
                
        num = 0.0
        for i in range(0, m):
            if i!=l:
                prod = 1.0
                for j in range(0, m):
                    if (j!=l) and (j!=i):
                        prod = prod*(y-xx[j])
                num = num + prod
        LL[l] = num/den

    return LL
