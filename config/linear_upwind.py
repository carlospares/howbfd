from initcond import InitCond
from equation import Equation
from sw_equation import SWEquation
from boundary import BoundaryCond
from numflux import Flux

# For a detailed explanation, see howbfd_config

equation = Equation.LINEAR
init = InitCond.STEADY
perturb_init = InitCond.PERT_PATCH
boundary = BoundaryCond.FORCE_STEADY
numflux = Flux.UPWIND
order = 3
well_balanced = True
N = 100
cfl = 0.5
a = -1
b = 1
T = 40
plot_every = 2
show_plots = True
save_plots = False
save_npys = False
sw_H = SWEquation.H_PWPOLY # not relevant
sw_H_noise_factor = 0.1 # not relevant