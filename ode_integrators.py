#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import numpy as np
from equation import Equation
from functionH import FunH



def odeint(nsteps, multmeth, eqn, arg0, arg1, arg2, arg3, arg4):
    if nsteps == 4 and multmeth == 'AM':
        return adamsmoulton4(eqn, arg0, arg1, arg2, arg3, arg4)
    elif nsteps == 4 and multmeth == 'AB':
        return adamsbashforth4(eqn, arg0, arg1, arg2, arg3, arg4)
    elif nsteps == 6 and multmeth == 'AM':
        return adamsmoulton6(eqn, arg0, arg1, arg2, arg3, arg4)
    elif nsteps == 6 and multmeth == 'AB':
        return adamsbashforth6(eqn, arg0, arg1, arg2, arg3, arg4)
    elif nsteps == 8 and multmeth == 'AM' :
        return adamsmoulton8(eqn, arg0, arg1, arg2, arg3, arg4)
    elif nsteps == 8 and multmeth == 'AB' :
        return adamsbashforth8(eqn, arg0, arg1, arg2, arg3, arg4)

def adamsbashforth4(eqn, Hx, u, x, i, t):
    nsteps= 4
    ab_coeff=[-9./24., 37./24., -59./24., 55./24]

    sumSHx = 0.
    for j in [-4, -3, -2, -1]:
        sumSHx += ab_coeff[j+nsteps]*eqn.S(u[:,i+j])*Hx(x[i+j],t)

    return sumSHx

def adamsbashforth6(eqn, Hx, u, x, i, t):
    nsteps= 6
    ab_coeff=[-475./1440., 2877./1440., -7298./1440., 9982./1440,  -7923./1440., 4277./1440.]

    sumSHx = 0.
    for j in [-6, -5, -4, -3, -2, -1]:
        sumSHx += ab_coeff[j+nsteps]*eqn.S(u[:,i+j])*Hx(x[i+j], t)

    return sumSHx

def adamsbashforth8(eqn, Hx, u, x, i, t):
    nsteps= 8
    ab_coeff=[-36799/120960., 295767./120960., -1041723./120960., 2102243./120960.,  -2664477./120960., 2183877./120960., -1152169./120960., 434241./120960.]

    sumSHx = 0.
    for j in [-8, -7, -6, -5, -4, -3, -2, -1]:
        sumSHx += ab_coeff[j+nsteps]*eqn.S(u[:,i+j])*Hx(x[i+j], t)

    return sumSHx

def adamsmoulton4(eqn, Hx, u, x, i, t):
    nsteps= 4
    ab_coeff=[1./24., -5./24., 19./24., 9./24]

    sumSHx = 0.
    for j in [-3, -2, -1, 0]:
        sumSHx += ab_coeff[j+nsteps-1]*eqn.S(u[:,i+j])*Hx(x[i+j], t)

    return sumSHx

def adamsmoulton6(eqn, Hx, u, x, i, t):
    nsteps= 6
    ab_coeff=[27./1440., -173./1440., 482./1440., -798./1440,  1427./1440., 475./1440.]

    sumSHx = 0.
    for j in [-5, -4, -3, -2, -1, 0]:
        sumSHx += ab_coeff[j+nsteps-1]*eqn.S(u[:,i+j])*Hx(x[i+j], t)

    return sumSHx

def adamsmoulton8(eqn, Hx, u, x, i, t):
    nsteps= 8
    ab_coeff=[1375./120960., -11351./120960., 41499./120960.,  -88547./120960., 123133./120960., -12179./12090, 139849./120960., 36799/120960.]

    sumSHx = 0.
    for j in [-7, -6, -5, -4, -3, -2, -1, 0]:
        sumSHx += ab_coeff[j+nsteps-1]*eqn.S(u[:,i+j])*Hx(x[i+j], t)

    return sumSHx
