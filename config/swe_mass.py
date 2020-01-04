from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from nummeth import NumericalMethod
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.SW
init = InitCond.WATER_MASS
funh = FunH.SLOPE
H_noise_factor = 0
perturb_init = InitCond.PERT_NONE
boundary = BoundaryCond.WALL
nummeth = NumericalMethod.RUSANOVGWB1
timest = TimeStepping.EULER
order = 3
N = 200
a = -5
b = 25
cfl = 0.5
T =.5
plot_every = 0.5
show_plots = True
save_plots = False
save_npys = True
plot_exact = False

print "Loaded config/swe_rest.py!"