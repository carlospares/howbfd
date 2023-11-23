# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

from nummeth import NumericalMethod
from nm_upwind import Upwind
from nm_upwind_wb import UpwindWB
from nm_upwind_wb_cons import UpwindWBCons
from nm_upwind_wb1 import UpwindWB1
from nm_upwind_gf import UpwindGF #added by Maria
from nm_rusanovg import RusanovG
from nm_rusanovg_wb import RusanovGWB
from nm_rusanov import Rusanov
from nm_rusanov_wb import RusanovWB
from nm_rusanovg_wb_cons import RusanovGWBCons
from nm_flsplit import FluxSplit
from nm_flsplit_wb import FluxSplitWB
from nm_flsplit_wb_cons import FluxSplitWBCons
#from nm_swupwind import SWUpwind
#from nm_swupwind_wb import SWUpwindWB
from nm_rusanovg_wb1 import RusanovGWB1
#from eq_burgers import BurgersEquation
#from eq_sw import SWEquation
#from eq_sw_war import SWEquationWAR

def nummeth_factory(cf):
    nm = cf.nummeth
    if nm == NumericalMethod.UPWIND:
        return Upwind(cf)
    elif nm == NumericalMethod.UPWINDWB:
        return UpwindWB(cf)
    elif nm == NumericalMethod.UPWINDWBCONS:
        return UpwindWBCons(cf)
    elif nm == NumericalMethod.UPWINDWB1:
        return UpwindWB1(cf)
    elif nm == NumericalMethod.UPWINDGF: #added by Maria
        return UpwindGF(cf)
    elif nm == NumericalMethod.RUSANOVG:
        return RusanovG(cf)
    elif nm == NumericalMethod.RUSANOVGWB:
        return RusanovGWB(cf)
    elif nm == NumericalMethod.RUSANOV:
        return Rusanov(cf)
    elif nm == NumericalMethod.RUSANOVWB:
        return  RusanovWB(cf)
    elif nm == NumericalMethod.RUSANOVGWBCONS:
        return RusanovGWBCons(cf)
    elif nm == NumericalMethod.FLUXSPLIT:
        return FluxSplit(cf)
    elif nm == NumericalMethod.FLUXSPLITWB:
        return FluxSplitWB(cf)
    elif nm == NumericalMethod.FLUXSPLITWBCONS:
        return FluxSplitWBCons(cf)
#    elif nm == NumericalMethod.SWUPWIND:
#        return  SWUpwind(cf)
#    elif nm == NumericalMethod.SWUPWINDWB:
#        return SWUpwindWB(cf)
    elif nm == NumericalMethod.RUSANOVGWB1:
        return RusanovGWB1(cf)
    else:
        print "[ERROR] Equation not recognized!!"
        print "	...Returning generic Numerical Method object, but nothing will work."
        return NumericalMethod()
