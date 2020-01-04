from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from numflux import Flux
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.SW
init = InitCond.WATER_MASS
funh = FunH.BUMP2
H_noise_factor = 0
perturb_init = InitCond.PERT_NONE
boundary = BoundaryCond.WALL
numflux = Flux.RUSANOV
timest = TimeStepping.EULER
order = 3
well_balanced = False
is_conservative = True
N = 200
a = -5
b = 25
cfl = 0.5
T =2.
plot_every = 0.5
show_plots = True
save_plots = False
save_npys = True

print "Loaded config/swe_rest.py!"