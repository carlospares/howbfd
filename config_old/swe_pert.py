from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from numflux import Flux
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.SW
init = InitCond.STEADY
funh = FunH.BUMP2
H_noise_factor = 0
perturb_init = InitCond.PERT_NONE
boundary = BoundaryCond.INIT
numflux = Flux.UPWIND_SW
timest = TimeStepping.TVDRK3
order = 3
well_balanced = True
is_conservative = True
N = 150
a = -3
b = 3
cfl = 0.01
T = .05
plot_exact = False
plot_every = .01
show_plots = True
save_plots = False
save_npys = True



print "Loaded config/swe_rest.py!"