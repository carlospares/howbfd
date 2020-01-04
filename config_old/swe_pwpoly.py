from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from numflux import Flux
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.SW
init = InitCond.PWPOLY
funh = FunH.BUMP
H_noise_factor = 0.1
perturb_init = InitCond.PERT_NONE
boundary = BoundaryCond.FORCE_STEADY_INIT
numflux = Flux.RUSANOV
timest = TimeStepping.EULER
order = 3
well_balanced = True
is_conservative = True
N = 100
a = 0
b = 20
cfl = 0.5
T = 8
plot_every = 0.5
show_plots = True
save_plots = False
save_npys = False

print "Loaded config/swe_pwpoly.py!"