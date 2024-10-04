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
perturb_init = InitCond.PERT_WB
boundary = BoundaryCond.FORCE_STEADY_INIT
nummeth = NumericalMethod.RUSANOVGWB1
timest = TimeStepping.TVDRK3
order = 3
N = 100
a = 0
b = 20
cfl = 0.5
T = 8
plot_every = 0.5
plot_exact= False
show_plots = True
save_plots = False
save_npys = False

print ("Loaded config/swe_pwpoly.py!")
