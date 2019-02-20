from initcond import InitCond
from equation import Equation
from boundary import BoundaryCond
from numflux import Flux

# For a detailed explanation, see howbfd_config

equation = Equation.SWE_REST
init = InitCond.PWPOLY
sw_H = Equation.SWE_H_PWPOLY
sw_H_noise_factor = 0.1
perturb_init = InitCond.PERT_NONE
boundary = BoundaryCond.FORCE_STEADY
numflux = Flux.RUSANOV
order = 3
well_balanced = True
N = 100
cfl = 0.5
a = -1
b = 1
T = 8
plot_every = 0.5
show_plots = True
save_plots = False
save_npys = False