from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from numflux import Flux
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.BURGERS
init = InitCond.STEADY
funh = FunH.IDENT
H_noise_factor = 0.0
perturb_init = InitCond.PERT_NONE
boundary = BoundaryCond.FORCE_STEADY
numflux = Flux.RUSANOV
timest = TimeStepping.TVDRK2
order = 3
well_balanced = True
is_conservative = True
N = 100
cfl = 0.5
a = -1
b = 1
T = 2
plot_every = 0.5
show_plots = True
save_plots = False
save_npys = False

print "Loaded config/burgers_rusanov.py!"