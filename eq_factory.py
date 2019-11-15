# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

from equation import Equation
from eq_linear import LinearEquation
from eq_burgers import BurgersEquation
from eq_sw import SWEquation
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
    else:
        print "[ERROR] Equation not recognized!!"
        print "	...Returning generic Equation object, but nothing will work."
        return Equation()