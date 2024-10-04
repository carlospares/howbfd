from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from nummeth import NumericalMethod
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.LINEAR
init = InitCond.TWO_ST
funh = FunH.IDENT
H_noise_factor = 0.0
perturb_init = InitCond.PERT_NONE
boundary = BoundaryCond.IN_OUT
nummeth = NumericalMethod.UPWINDGF
timest = TimeStepping.TVDRK3
order = 5
N = 100
cfl = 0.5
a = -.5
b = 2.
T = 1.
plot_every = .25
show_plots = True
save_plots = True
save_npys = True
plot_exact = False

print ("Loaded config/linear_upwind.py!")
