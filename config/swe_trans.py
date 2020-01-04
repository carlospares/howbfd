from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from nummeth import NumericalMethod
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.SW
init = InitCond.STEADY
funh = FunH.BUMPD
H_noise_factor = 0
perturb_init = InitCond.PERT_NONE
boundary = BoundaryCond.IN_OUT
nummeth = NumericalMethod.UPWINDWB
timest = TimeStepping.TVDRK3
order = 3
N = 200
a = -3
b = 3
cfl = 0.5
T = 4.
plot_exact = False
plot_every = 1.
show_plots = True
save_plots = False
save_npys = True



print "Loaded config/swe_rest.py!"