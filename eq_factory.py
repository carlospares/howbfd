from equation import Equation
from linear_equation import LinearEquation
from burgers_equation import BurgersEquation
from sw_equation import SWEquation

def equation_factory(eqn, x=[0], H=SWEquation.H_FLAT, noise_amplit=0):
	if eqn == Equation.LINEAR:
		return LinearEquation()
	elif eqn == Equation.BURGERS:
		return BurgersEquation()
	elif eqn == Equation.SW:
		return SWEquation(x, H, noise_amplit)
	else:
		print "[ERROR] Equation not recognized!!"
		print "	...Returning generic Equation object, but nothing will work."
		return Equation()