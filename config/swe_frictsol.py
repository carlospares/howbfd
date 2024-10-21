from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from nummeth import NumericalMethod
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.SW
init = InitCond.Frictsol
funh = FunH.FRICTSOL
H_noise_factor = 0
perturb_init = InitCond.PERT_NONE#PERT_WB#PERT_NONE
boundary = BoundaryCond.IN_OUT
nummeth = NumericalMethod.UPWINDGF
timest = TimeStepping.TVDRK3
order = 3
N = 160
a = 0
b = 1
cfl = 0.6
T = 50.
steps=6
ode='AM'
system='SW'
plot_exact = False
plot_every = 1
show_plots = True
save_plots = False
save_npys = True



print("Loaded config/swe_rest.py!")
