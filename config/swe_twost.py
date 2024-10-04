from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from nummeth import NumericalMethod
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.SW
init = InitCond.TWO_ST_SW
funh = FunH.SLOPE
H_noise_factor = 0
perturb_init = InitCond.PERT_NONE
boundary = BoundaryCond.INIT
nummeth = NumericalMethod.UPWINDWB
timest = TimeStepping.TVDRK3
order = 3
N = 200
a = -10
b = 10
cfl = 0.5
T =.25
plot_every = 0.25
show_plots = True
save_plots = False
save_npys = False
plot_exact = False

print ("Loaded config/swe_rest.py!")
