# -*- coding: utf-8 -*-
#Ernesto Pimentel García

# -*- coding: utf-8 -*-
import time #No importamos todo con * porque puede coincidir el nombre de alguna función con otra función de otra importación.
#También podemos poner import time as tt, para llamar más rápido a dicho archivo.
from numpy import *
#numpy exporta muchos términos matemáticos
from matplotlib.pyplot import *
#Esta es para pintar.

from scipy.optimize import *

def feval(mf0,*args):
    """ Funcion auxiliar que permite la evaluacion de funciones
	definidas por el usuario"""
    return eval(mf0)(*args)

g = 9.81
heps = 1e-12

def phi(hl, ul, Hl, Hr): #Parte del estado WL= (hL, uL,HL) y calcula el h del estado al que se llega en el nivel Hr.
    def f(h): #Esta es la función que igualamos a 0 y a la que le aplicamos el método de Newton o el método de la biyección.
        y = Hl - Hr + ((ul**2)/(2*g))*((hl**2)/(h**2) - 1.) + h - hl
        return y
    hmin = ((ul**2)*(hl**2)/g)**(1/3.) #Es donde se encuentra el mínimo de la función.
    fhmin= f(hmin)
    tol1 = 1e-4
    tol2 = 1e-10
    tol3 = 1e-8
    #Cuando Hr es próximo al Hmin, la solución es próxima a hmin.
    if (abs(ul) < tol3): #Cuando la velocidad es próxima a 0:
        z = ((hl -Hl + Hr), (hl -Hl + Hr))
    elif (hmin< tol2): 
        z = (hmin, hmin)
    elif ((abs(fhmin) > tol1)):
        z = (bisect(f, 1e-20, hmin), newton(f, hmin*2)) #El primer estado es el supercrítico y el segundo es el subcrítico.
    else:
        z = (hmin, hmin)
    return z

def phiu(hl, ul, hphi): #para calcular la u a la que llega.
    if (hphi ==0):
        y=0
    else:
        y = (hl*ul)/hphi
    return y