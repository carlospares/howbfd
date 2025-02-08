#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import numpy as np
from functionH import FunH
from howbfd_io import IoManager, parse_command_line
from multistep_coefficients import AB_COEFFS
from multistep_coefficients import AM_COEFFS

### Get config file from command line, or load default:
config = parse_command_line() # from howbdf_io, defaults to howbdf_config
if config.ode =='AB':
    ab_coeff = AB_COEFFS[config.steps]

if config.ode =='AM':
    ab_coeff = AM_COEFFS[config.steps]

def B_odeint(eqn, Hx, H, x, i, t):
    nvars = eqn.dim()
    nsteps = config.steps
    if config.ode =='AB':
        sumSHx = 0.
        for j in range(-nsteps,0):
            sumSHx  += ab_coeff[j+nsteps]*( Hx(x[i+j],t) )

    if config.ode =='AM':
        sumSHx = 0.
        for j in range(-nsteps+1,1):
            sumSHx  += ab_coeff[j+nsteps-1]*( Hx(x[i+j],t) )

    return sumSHx

#------------------------------------------------------------------------------------------------------------------------------------

def odeint(eqn, arg0, arg1, arg2, arg3, arg4, arg5, arg6):
    nvars = eqn.dim()

    if config.system == 'No':
        if config.ode == 'AM':
                return adamsmoulton(eqn, arg1, arg2, arg3, arg4, arg5, arg6) # MARIO!!!!
        elif config.ode == 'AB':
                return adamsbashforth(eqn, arg1, arg2, arg3, arg4, arg5, arg6) 
    elif config.system == 'SW':
        if config.ode == 'AM':
            return adamsmoultonSW(eqn, arg0, arg1, arg2, arg3, arg4, arg5, arg6)
        elif config.ode == 'AB':
            return adamsbashforthSW(eqn, arg0, arg1, arg2, arg3, arg4, arg5, arg6) 

#------------------------------------------------------------------------------------------------------------------------------------

def adamsbashforth(eqn, Hx, H, u, x, i, t):
    funH = FunH(x, config)
    d_index = None
    if config.funh == FunH.DISC:
        d_index=funH.find_disc(x,1.0) #check again for the threshold
        
    Y = funH.get_disc_points(x)
        
    dx = x[2] - x[1]
    nvars = eqn.dim()
    sumSHx = np.zeros(nvars)
    nsteps= config.steps
    
    indicator = 'normal'
    if d_index != None :
        for num in d_index:
            if num + 1 == i:
                indicator = 'jump'
            elif i >= num+1+1 and i<=num + nsteps:
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
        sumSHx[nvars-1] += eqn.S(u[:,i-1])*Hx(x[i-1],t)*0.5 - eqn.sigma(u[:,i-1])*0.5
    
        #Right integration
        sumSHx[nvars-1] += eqn.S(u[:,i])*Hx(x[i],t)*0.5 - eqn.sigma(u[:,i])*0.5
    else :
        for k in range(-nsteps,0):
            sumSHx[nvars-1] += ab_coeff[k+nsteps]*( eqn.S(u[:,i+k])*Hx(x[i+k],t) - eqn.sigma(u[:,i+k]) )

#------------------------------------------------------------------------------------------------------------------------------------

def adamsbashforthSW(eqn, B, Hx, H, u, x, i, t):
    funH = FunH(x, config)
    d_index = None
    if config.funh == FunH.DISC or config.funh == FunH.STEP:
        d_index=funH.find_disc(x,1.0) #check again for the threshold
        
    Y = funH.get_disc_points(x)
        
    ddx = x[i] - x[i-1]
    dx=ddx

    nvars = eqn.dim()
    nsteps= config.steps
    
    indicator = 'normal'
    if d_index != None :
        for num in d_index:
            if num +1 == i:
                indicator = 'jump'
            elif i >= num+1 and i<=num + 1+ nsteps:
                indicator='AM2'

    compute_source = config.compute_source
    #print(indicator,d_index,i)

    g = 9.812

    # Collect:
    # - stencil nodes
    # - eta values (REMARK: really specific to gravity source !!!!)
    # - bathymetry and bathymetry derivatives values
    xx  = np.zeros(nsteps)
    eta = np.zeros(nsteps)
    bb  = np.zeros(nsteps)
    sig = np.zeros(nsteps)
    for l in range(0,nsteps):
        xx[l] = x[i-nsteps+l+1]
        FF = eqn.sigma(u[:,i-nsteps+l])
        sig[l]=FF[1]

    if ( indicator == 'AM2'):
        sumSHx=0.
        sumSHx  = adamsmoulton2SW(eqn, B, Hx, H, u, x, i, t)
        #print(indicator,x[i],i)
    elif (indicator == 'jump'):
        #print(indicator,x[i])
        sumSHx = 0.
        sumSHx = np.zeros(nvars)
        dH = H(x[i-1]+ 0.0000000001, t ) - H(x[i-1] - 0.0000000001, t ) #if the dicontinuity is on a mesh point
        if(abs(dH) <= 0.000001):
            dH = H(x[i-1]+  dx , t ) - H(x[i-1] , t ) #if the disc is on the face
            #dH = H(x[i-1]+ 0.5*dx + 0.0000000001, t ) - H(x[i-1]+ 0.5*dx - 0.0000000001, t ) #if the disc is on the face 
        
        delta = eqn.discH_jumpF( u[:,i-1], u[:,i], i, dH, x, t)
        sumSHx[nvars-1]  += delta/dx

        #Left integration
        tmp  = eqn.S(u[:,i-1])*Hx(x[i-1],t)*0.5 - eqn.sigma(u[:,i-1])*0.5
        sumSHx[nvars-1] += tmp[1]
    
        #Right integration
        tmp  = eqn.S(u[:,i])*Hx(x[i],t)*0.5 - eqn.sigma(u[:,i])*0.5
        sumSHx[nvars-1] += tmp[1]
    else:
        #print(indicator,x[i],i)

        if compute_source == 'analytic_source':
    
            sumSHx=0
            for j in range(-nsteps,0):
                sumSHx += ab_coeff[j+nsteps]*( eqn.S(u[:,i+j])*Hx(x[i+j],t)) #- sig[j+nsteps])

        elif compute_source == 'source_reconstruction':

            for l in range(0,nsteps):
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
            for j in range(-nsteps,0):
                sumSHx[nvars-1] += ab_coeff[j+nsteps]*( g*eta[j+nsteps]*Bx[j+nsteps] - sig[j+nsteps] )

        elif compute_source == 'hydrostatic_reconstruction':

            for l in range(0,nsteps):
                bb[l] = B[i-nsteps+l+1] #reconstructed topography- be carefull has the value in i and the i-(nsteps-1) nodes
                eta[l] = u[0,i-nsteps+l] # new version of keeping the lake at rest

            sumSHx = np.zeros(nvars)
            sumSHx[1] = 0.5*g*( bb[nsteps-1]*bb[nsteps-1] - bb[nsteps-2]*bb[nsteps-2]  )/ddx
            for j in range(-nsteps,0):
                FF = eqn.sigma(u[:,j+nsteps])
                sumSHx[nvars-1] += ab_coeff[j+nsteps]*( g*eta[j+nsteps]*Hx(x[j+i],t)- sig[j+nsteps]) #new version for keeping lake at rest solving for eta
        else :
            print(compute_source, 'This type of reconstruction does not exist')

    return sumSHx

#------------------------------------------------------------------------------------------------------------------------------------

def adamsmoulton(eqn, Hx, H, u, x, i, t ):
    funH = FunH(x, config)
    d_index = None
    if config.funh == FunH.DISC:
        d_index=funH.find_disc(x,1.0) #check again for the threshold
        
    Y = funH.get_disc_points(x)
        
    dx = x[2] - x[1]
    nvars = eqn.dim()
    sumSHx = np.zeros(nvars)
    nsteps= config.steps
    
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
            #dH = H(x[i-1]+  dx , t ) - H(x[i-1] , t ) #if the disc is on the face
            dH = H(x[i-1]+ 0.5*dx + 0.0000000001, t ) - H(x[i-1]+ 0.5*dx - 0.0000000001, t ) #if the disc is on the face 
        
        delta = eqn.discH_jumpF( u[:,i-1], u[:,i], i, dH, x, t)
        sumSHx[nvars-1] += delta/dx

        #Left integration
        sumSHx[nvars-1] += eqn.S(u[:,i-1])*Hx(x[i-1],t)*0.5 - eqn.sigma(u[:,i-1])*0.5
    
        #Right integration
        sumSHx[nvars-1] += eqn.S(u[:,i])*Hx(x[i],t)*0.5 - eqn.sigma(u[:,i])*0.5
    else :
        for k in range(-nsteps+1,1):
            sumSHx[nvars-1] += ab_coeff[k+nsteps-1]*( eqn.S(u[:,i+k])*Hx(x[i+k], t) - eqn.sigma(u[:,i+k]) )

    return sumSHx



#------------------------------------------------------------------------------------------------------------------------------------
def adamsmoultonSW(eqn, B, Hx, H, u, x, i, t):
    funH = FunH(x, config)
    d_index = None
    if config.funh == FunH.STEP:
        d_index=funH.find_disc(x,1.0) #check again for the threshold
        
    Y = funH.get_disc_points(x)
        
    dx = x[2] - x[1]
    ddx = x[i] - x[i-1]
    nvars = eqn.dim()
    sumSHx = np.zeros(nvars)
    nsteps= config.steps

    compute_source = config.compute_source
    
    indicator = 'normal'
    if d_index != None :
        for num in d_index:
            if num + 1 == i:
                indicator = 'jump'
            elif i >= num+1+1 and i<=num + 1 + nsteps:
                indicator='AM2'

    g = 9.812

    # Collect:
    # - stencil nodes
    # - eta values (REMARK: really specific to gravity source !!!!)
    # - bathymetry and bathymetry derivatives values
    xx  = np.zeros(nsteps)
    eta = np.zeros(nsteps)
    sig = np.zeros(nsteps)
    bb  = np.zeros(nsteps)
    for l in range(0,nsteps):
        xx[l] = x[i-nsteps+l+1]
        FF = eqn.sigma(u[:,i-nsteps+l+1])
        sig[l]=FF[1]


    if ( indicator == 'AM2'):
        #print(indicator,i)
        sumSHx=0.
        sumSHx = adamsmoulton2SW(eqn, B, Hx, H, u, x, i, t)
    elif (indicator == 'jump'):
        #print(indicator,i)
        sumSHx = np.zeros(nvars)
        dH = H(x[i-1]+ 0.0000000001, t ) - H(x[i-1] - 0.0000000001, t ) #if the dicontinuity is on a mesh point
        if(abs(dH) <= 0.000001):
            dH = H(x[i-1]+  dx , t ) - H(x[i-1] , t ) #if the disc is on the face
            #dH = H(x[i-1]+ 0.5*dx + 0.0000000001, t ) - H(x[i-1]+ 0.5*dx - 0.0000000001, t ) #if the disc is on the face 
        
        #print(indicator,x[i],x[i-1],dH,u[:,i-1],u[:,i])
        delta = eqn.discH_jumpF( u[:,i-1], u[:,i], i, dH, x, t)
        sumSHx[nvars-1]  += delta/dx

        #Left integration
        tmp  = eqn.S(u[:,i-1])*Hx(x[i-1],t)*0.5 - eqn.sigma(u[:,i-1])*0.5
        sumSHx[nvars-1] +=tmp[1]
    
        #Right integration
        tmp  = eqn.S(u[:,i])*Hx(x[i],t)*0.5 - eqn.sigma(u[:,i])*0.5
        sumSHx[nvars-1] +=tmp[1]
    else:    
    #------------------------------analytic source  
        if compute_source == 'analytic_source':

            sumSHx=0.
            for j in range(-nsteps+1,1):
                sumSHx += ab_coeff[j+nsteps-1]*( eqn.S(u[:,i+j])*Hx(x[i+j], t) - sig[j+nsteps-1])

        elif compute_source == 'source_reconstruction':

    #------------------------------source reconstruction 
            for l in range(0,nsteps):
                eta[l] = -H(x[i-nsteps+l+1],t)+u[0,i-nsteps+l+1]
                bb[l] = H(x[i-nsteps+l+1],t)

            Bx = np.zeros(nsteps)
            for q in range(0,nsteps):
                Bx[q] = 0.0
                LL = Lprime( nsteps, xx, x[i-nsteps+q+1] )
                for p in range(0,nsteps):
                    Bx[q] = Bx[q] + LL[p]*bb[p]

            sumSHx = np.zeros(nvars)
            sumSHx[nvars-1] = 0.5*g*( bb[nsteps-1]*bb[nsteps-1] - bb[nsteps-2]*bb[nsteps-2]  )/ddx
            for j in range(-nsteps+1,1):
                sumSHx[nvars-1] += ab_coeff[j+nsteps-1]*( g*eta[j+nsteps-1]*Bx[j+nsteps-1] - sig[j+nsteps-1] )

        elif compute_source == 'hydrostatic_reconstruction':

    #----------------------hydrostatic reconstruction

            for l in range(0,nsteps):
                bb[l] = B[i-nsteps+l+1] #reconstructed topography
                eta[l] = u[0,i-nsteps+l+1]

            # Compute integrated source
            sumSHx = np.zeros(nvars)
            sumSHx[nvars-1] = 0.5*g*( bb[nsteps-1]*bb[nsteps-1] - bb[nsteps-2]*bb[nsteps-2]  )/ddx
            for j in range(-nsteps+1,1):
                sumSHx[nvars-1] += ab_coeff[j+nsteps-1]*( g*eta[j+nsteps-1]*Hx(x[j+i],t))
        else:
            print(compute_source, 'This type of reconstruction does not exist')
            print()

    return sumSHx

#------------------------------------------------------------------------------------------------------------------------------------
#def adamsbashforth2(eqn, Hx, H, u, x, i, t):
#    nvars = eqn.dim()
#    nsteps= 2
#    ab_coeff=[-1./2., 3./2]
#    
#    sumSHx = np.zeros(nvars)
#    for j in [-2, -1]:
#        sumSHx[nvars-1] += ab_coeff[j+nsteps]*( eqn.S(u[:,i+j])*Hx(x[i+j],t) - eqn.sigma(u[:,i+j]) )
#    return sumSHx
    

#------------------------------------------------------------------------------------------------------------------------------------
def adamsmoulton2(eqn, Hx, H, u, x, i, t):
    nvars = eqn.dim()
    nsteps= 2
    ab_coeff2=[1./2., 1./2.]

    sumSHx = np.zeros(nvars)
    for j in [-1, 0]:
        sumSHx[nvars-1] += ab_coeff2[j+nsteps-1]*( eqn.S(u[:,i+j])*Hx(x[i+j], t) - eqn.sigma(u[:,i+j]) )

    return sumSHx
#------------------------------------------------------------------------------------------------------------------------------------
def adamsmoulton2SW(eqn, B, Hx, H, u, x, i, t):

    compute_source = config.compute_source
 
    nvars = eqn.dim()
    nsteps = 2
    ab_coeff2=[1./2., 1./2.]

    ddx = x[i] - x[i-1]
    g = 9.812

    # Collect:
    # - stencil nodes
    # - eta values (REMARK: really specific to gravity source !!!!)
    # - bathymetry and bathymetry derivatives values
    xx  = np.zeros(nsteps)
    eta = np.zeros(nsteps)
    sig = np.zeros(nsteps)
    bb  = np.zeros(nsteps)
    for l in range(0,nsteps):
        xx[l] = x[i-nsteps+l+1]
        FF = eqn.sigma(u[:,i-nsteps+l+1])
        sig[l]=FF[1]

#------------------------------analytic source  
    if compute_source == 'analytic_source':

        sumSHx=0.
        for j in range(-nsteps+1,1):
            sumSHx += ab_coeff2[j+nsteps-1]*( eqn.S(u[:,i+j])*Hx(x[i+j], t) - sig[j+nsteps-1])

    elif compute_source == 'source_reconstruction':

#------------------------------source reconstruction 
        for l in range(0,nsteps):
            eta[l] = -H(x[i-nsteps+l+1],t)+u[0,i-nsteps+l+1]
            bb[l] = H(x[i-nsteps+l+1],t)

        Bx = np.zeros(nsteps)
        for q in range(0,nsteps):
            Bx[q] = 0.0
            LL = Lprime( nsteps, xx, x[i-nsteps+q+1] )
            for p in range(0,nsteps):
                Bx[q] = Bx[q] + LL[p]*bb[p]

        sumSHx = np.zeros(nvars)
        sumSHx[nvars-1] = 0.5*g*( bb[nsteps-1]*bb[nsteps-1] - bb[nsteps-2]*bb[nsteps-2]  )/ddx
        for j in range(-nsteps+1,1):
            sumSHx[nvars-1] += ab_coeff2[j+nsteps-1]*( g*eta[j+nsteps-1]*Bx[j+nsteps-1] - sig[j+nsteps-1] )

    elif compute_source == 'hydrostatic_reconstruction':

#----------------------hydrostatic reconstruction

        for l in range(0,nsteps):
            bb[l] = B[i-nsteps+l+1] #reconstructed topography
            eta[l] = u[0,i-nsteps+l+1]

        # Compute integrated source
        sumSHx = np.zeros(nvars)
        sumSHx[nvars-1] = 0.5*g*( bb[nsteps-1]*bb[nsteps-1] - bb[nsteps-2]*bb[nsteps-2]  )/ddx
        for j in range(-nsteps+1,1):
            sumSHx[nvars-1] += ab_coeff2[j+nsteps-1]*( g*eta[j+nsteps-1]*Hx(x[j+i],t))
    else:
        print(compute_source, 'This type of reconstruction does not exist')
        print()

    return sumSHx
#-------------------------------------------------------------------------------------------------------------------

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
#-------------------------------------------------------------------------------------------------------------------
    
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
#-------------------------------------------------------------------------------------------------------------------
    
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
