# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 19:23:29 2019

@author: Usuario
"""
a = -.5
b = 2
N = 200
x = linspace(a,b, 200)

tag3wb = 'i508-600e700H301f100b401w1n200o3_1.0.npy'
tag3nwb = 'i508-600e700H301f100b401w0n200o3_1.0.npy'
tag5wb = 'i508-600e700H301f100b401w1n200o5_1.0.npy'
tag5nwb = 'i508-600e700H301f100b401w0n200o5_1.0.npy'
u3wb = load(tag3wb)
u3nwb = load(tag3nwb)
u5wb = load(tag5wb)
u5nwb = load(tag5nwb)

plot(x,u3wb, x, u3nwb, x, u5wb, x, u5nwb)
show()
