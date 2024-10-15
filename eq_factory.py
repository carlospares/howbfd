    # -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

from equation import Equation
from eq_linear import LinearEquation
from eq_burgers import BurgersEquation
from eq_sw import SWEquation
from eq_euler_grav import EulerEquationGRAV
from eq_euler_sec import EulerEquationSEC
from eq_sw_war import SWEquationWAR

def equation_factory(cf):
    eqn = cf.equation
    if eqn == Equation.LINEAR:
        return LinearEquation()
    elif eqn == Equation.BURGERS:
        return BurgersEquation()
    elif eqn == Equation.SW:
        return SWEquation()
    elif eqn == Equation.SW_WAR:
        return SWEquationWAR()
    elif eqn == Equation.EulerGRAV:
        return EulerEquationGRAV()
    elif eqn == Equation.EulerSEC:
        return EulerEquationSEC()
    else:
        print ("[ERROR] Equation not recognized!!")
        print ("	...Returning generic Equation object, but nothing will work.")
        return Equation()
