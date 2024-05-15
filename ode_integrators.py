#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import numpy as np
from equation import Equation
from functionH import FunH
from howbfd_io import IoManager, parse_command_line

### Get config file from command line, or load default:
config = parse_command_line() # from howbdf_io, defaults to howbdf_config

def odeint(nsteps, multmeth, eqn, arg0, arg1, arg2, arg3, arg4, arg5):
    nvars = eqn.dim()
    if config.system == 'No':
        if multmeth == 'AM':
            if nsteps == 2 :
                return adamsmoulton2(eqn, arg0, arg1, arg2, arg3, arg4, arg5) # MARIO!!!!
            if nsteps == 3 :
                return adamsmoulton3(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
            if nsteps == 4 :
                return adamsmoulton4(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
            if nsteps == 6 :
                return adamsmoulton6(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
            if nsteps == 8 :
                return adamsmoulton8(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
        elif multmeth == 'AB':
            if nsteps == 2 :
                return adamsbashforth2(eqn, arg0, arg1, arg2, arg3, arg4, arg5) # MARIO!!!!
            if nsteps == 3 :
                return adamsbashforth3(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
            if nsteps == 4 :
                return adamsbashforth4(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
            if nsteps == 6 :
                return adamsbashforth6(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
            if nsteps == 8 :
                return adamsbashforth8(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
    elif config.system == 'SW':
        if multmeth == 'AM':
            if nsteps == 2 :
                return adamsmoulton2SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5) # MARIO!!!!
            if nsteps == 3 :
                return adamsmoulton3SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
            if nsteps == 4 :
                return adamsmoulton4SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
            if nsteps == 6 :
                return adamsmoulton6SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
            if nsteps == 8 :
                return adamsmoulton8SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
        elif multmeth == 'AB':
            if nsteps == 2 :
                return adamsbashforth2SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5) # MARIO!!!!
            if nsteps == 3 :
                return adamsbashforth3SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
            if nsteps == 4 :
                return adamsbashforth4SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
            if nsteps == 6 :
                return adamsbashforth6SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
            if nsteps == 8 :
                return adamsbashforth8SW(eqn, arg0, arg1, arg2, arg3, arg4, arg5)
        
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

def adamsmoulton4(eqn, Hx, H, u, x, i, t ):
    funH = FunH(x, config)
    d_index = None
    if config.funh == FunH.DISC:
        d_index=funH.find_disc(x,1.0) #check again for the threshold
        
    Y = funH.get_disc_points(x)
        
    dx = x[2] - x[1]
    nvars = eqn.dim()
    sumSHx = np.zeros(nvars)
    nsteps= 4
    ab_coeff=[1./24., -5./24., 19./24., 9./24]
    
    indicator = 'normal'
    if d_index != None :
        for num in d_index:
            if num + 1 == i:
                indicator = 'jump'
            elif i >= num+1+1 and i<=num + 1 + nsteps:
                indicator='AM2'

    if ( indicator == 'AM2'):
            sumSHx[nvars-1] += adamsmoulton2(eqn, Hx, H, u, x, i, t)
    elif (indicator == 'jump'):
        #print i, 'hello'
        #dH = H(Y[j]+ 0.0000000001, t ) - H(Y[j] - 0.0000000001, t )
        dH = H(x[i-1]+ 0.0000000001, t ) - H(x[i-1] - 0.0000000001, t ) #if the dicontinuity is on a mesh point
        if(abs(dH) <= 0.000001):
            dH = H(x[i-1]+  dx , t ) - H(x[i-1] , t ) #if the disc is on the face
            #dH = H(x[i-1]+ 0.5*dx + 0.0000000001, t ) - H(x[i-1]+ 0.5*dx - 0.0000000001, t ) #if the disc is on the face 
        
        delta = eqn.discH_jumpF( u[:,i-1], u[:,i], i, dH, x, t)
        sumSHx[nvars-1] += delta/dx

        #Left integration
        sumSHx[nvars-1] += eqn.S(u[:,i-1])*Hx(x[i-1],t)*0.5
    
        #Right integration
        sumSHx[nvars-1] += eqn.S(u[:,i])*Hx(x[i],t)*0.5
    else :
        for k in [-3, -2, -1, 0]:
            sumSHx[nvars-1] += ab_coeff[k+nsteps-1]*eqn.S(u[:,i+k])*Hx(x[i+k], t)

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
    funH = FunH(x, config)
    d_index = None
    if config.funh == FunH.DISC:
        d_index=funH.find_disc(x,1.0) #check again for the threshold

    dx = x[2] - x[1]
    nvars = eqn.dim()
    sumSHx = np.zeros(nvars)
    nsteps= 6
    ab_coeff=[27./1440., -173./1440., 482./1440., -798./1440,  1427./1440., 475./1440.]

    indicator = 'normal'
    if d_index != None :
        for num in d_index:
            if num + 1 == i:
                indicator = 'jump'
            elif i >= num+1+1 and i<=num + 1 + nsteps:
                indicator='AM2'

    if ( indicator == 'AM2'):
            sumSHx[nvars-1] += adamsmoulton2(eqn, Hx, H, u, x, i, t)
    elif (indicator == 'jump'):
        #print i, 'hello'
        dH = H(x[i-1]+ 0.0000000001, t ) - H(x[i-1] - 0.0000000001, t ) #if the dicontinuity is on a mesh point
        if(abs(dH) <= 0.000001):
            dH = H(x[i-1]+  dx , t ) - H(x[i-1] , t ) #if the disc is on the face
            #dH = H(x[i-1]+ 0.5*dx + 0.0000000001, t ) - H(x[i-1]+ 0.5*dx - 0.0000000001, t ) #if the disc is on the face

        delta = eqn.discH_jumpF( u[:,i-1], u[:,i], i, dH, x, t)
        sumSHx[nvars-1] += delta/dx

        #Left integration
        sumSHx[nvars-1] += eqn.S(u[:,i-1])*Hx(x[i-1],t)*0.5

        #Right integration
        sumSHx[nvars-1] += eqn.S(u[:,i])*Hx(x[i],t)*0.5
    else :
        for k in [-5, -4, -3, -2, -1, 0]:
            sumSHx[nvars-1] += ab_coeff[k+nsteps-1]*eqn.S(u[:,i+k])*Hx(x[i+k], t)

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
    funH = FunH(x, config)
    d_index = None
    if config.funh == FunH.DISC:
        d_index=funH.find_disc(x,1.0) #check again for the threshold

    dx = x[2] - x[1]
    nvars = eqn.dim()
    sumSHx = np.zeros(nvars)
    nsteps= 8
    ab_coeff=[1375./120960., -11351./120960., 41499./120960.,  -88547./120960., 123133./120960., -121797./120960, 139849./120960., 36799/120960.]


    indicator = 'normal'
    if d_index != None :
        for num in d_index:
            if num + 1 == i:
                indicator = 'jump'
            elif i >= num+1+1 and i<=num + 1 + nsteps:
                indicator='AM2'

    if ( indicator == 'AM2'):
            sumSHx[nvars-1] += adamsmoulton2(eqn, Hx, H, u, x, i, t)
    elif (indicator == 'jump'):
        #print i, 'hello'
        dH = H(x[i-1]+ 0.0000000001, t ) - H(x[i-1] - 0.0000000001, t ) #if the dicontinuity is on a mesh point
        if(abs(dH) <= 0.000001):
            dH = H(x[i-1]+  dx , t ) - H(x[i-1] , t ) #if the disc is on the face
            #dH = H(x[i-1]+ 0.5*dx + 0.0000000001, t ) - H(x[i-1]+ 0.5*dx - 0.0000000001, t ) #if the disc is on the face

        delta = eqn.discH_jumpF( u[:,i-1], u[:,i], i, dH, x, t)
        sumSHx[nvars-1] += delta/dx

        #Left integration
        sumSHx[nvars-1] += eqn.S(u[:,i-1])*Hx(x[i-1],t)*0.5

        #Right integration
        sumSHx[nvars-1] += eqn.S(u[:,i])*Hx(x[i],t)*0.5
    else :
        for k in [-7, -6, -5, -4, -3, -2, -1, 0]:
            sumSHx[nvars-1] += ab_coeff[k+nsteps-1]*eqn.S(u[:,i+k])*Hx(x[i+k], t)

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
    
def Lbasis(m, xx, y):

    LL = np.zeros(m)
    
    for l in range(0, m):
        den = 1.0
        for j in range(0, m):
            if j!=l:
                den = den*(xx[l]-xx[j])
                
        num = 1.0
        for j in range(0, m):
            if j!=l:
                num = num*(y-xx[j])
        
        LL[l] = num/den

    return LL
    
def disc_int(eqn, xi, xip1, m, xx, uu, Hx, t):
 
    Ix = 0
    dx = xip1-xi
    
    # 6 points Gauss-Legendre formula
     # Point 1
    s = 0.238619186083197 ;
    y = ( 1.0 - s )*0.5*xi + ( 1.0 + s )*0.5*xip1
    w = 0.5*0.467913934572691
    
    LL = Lbasis(m,xx,y)
    for q in range(0,m):
        Ix += w*dx*LL[q]*Hx(xx[q],t)*eqn.S(uu[:,q])

   
    # Point 2
    s = -0.238619186083197
    y = ( 1.0 - s )*0.5*xi + ( 1.0 + s )*0.5*xip1
    w = 0.5*0.467913934572691

    LL = Lbasis(m,xx,y)
    for q in range(0,m):
        Ix += w*dx*LL[q]*Hx(xx[q],t)*eqn.S(uu[:,q])
     
     # Point 3
    s = 0.661209386466265
    y = ( 1.0 - s )*0.5*xi + ( 1.0 + s )*0.5*xip1
    w = 0.5*0.360761573048139

    LL = Lbasis(m,xx,y)
    for q in range(0,m):
        Ix += w*dx*LL[q]*Hx(xx[q],t)*eqn.S(uu[:,q])
        
     # Point 4
    s = -0.661209386466265
    y = ( 1.0 - s )*0.5*xi + ( 1.0 + s )*0.5*xip1
    w = 0.5*0.360761573048139
    
    LL = Lbasis(m,xx,y)
    for q in range(0,m):
        Ix += w*dx*LL[q]*Hx(xx[q],t)*eqn.S(uu[:,q])

     # Point 5
    s = 0.932469514203152
    y = ( 1.0 - s )*0.5*xi + ( 1.0 + s )*0.5*xip1
    w =  0.5*0.171324492379170

    LL = Lbasis(m,xx,y)
    for q in range(0,m):
        Ix += w*dx*LL[q]*Hx(xx[q],t)*eqn.S(uu[:,q])
        
     # Point 6
    s = -0.932469514203152
    y = ( 1.0 - s )*0.5*xi + ( 1.0 + s )*0.5*xip1
    w =  0.5*0.171324492379170
    
    LL = Lbasis(m,xx,y)
    for q in range(0,m):
        Ix += w*dx*LL[q]*Hx(xx[q],t)*eqn.S(uu[:,q])
    
    return Ix
