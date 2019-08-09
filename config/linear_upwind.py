from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from numflux import Flux
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.LINEAR
init = InitCond.SIN
funh = FunH.FLAT
H_noise_factor = 0.0
perturb_init = InitCond.PERT_NONE
boundary = BoundaryCond.IN_OUT
numflux = Flux.UPWIND
timest = TimeStepping.TVDRK2
order = 3
well_balanced = True
is_conservative = True
N = 100
cfl = 0.5
a = -1
b = 1
T = 15
plot_every = 2
show_plots = True
save_plots = False
save_npys = False

print "Loaded config/linear_upwind.py!"