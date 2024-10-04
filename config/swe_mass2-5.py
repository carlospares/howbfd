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
funh = FunH.BUMP
H_noise_factor = 0
perturb_init = InitCond.PERT_WM
boundary = BoundaryCond.FORCE_STEADY
nummeth = NumericalMethod.RUSANOVG
timest = TimeStepping.TVDRK3
order = 5
N = 200
a = -10
b = 30
cfl = 0.3
T =.5
plot_every = 0.5
show_plots = True
save_plots = False
save_npys = True
plot_exact = False

print )"Loaded config/swe_rest.py!")
