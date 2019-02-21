from initcond import InitCond
from equation import Equation
from eq_sw import SWEquation
from boundary import BoundaryCond
from numflux import Flux

# For a detailed explanation, see howbfd_config
equation = Equation.SW
init = InitCond.PWPOLY
perturb_init = InitCond.PERT_NONE
sw_H = SWEquation.H_PWPOLY
sw_H_noise_factor = 0.1
boundary = BoundaryCond.FORCE_STEADY
numflux = Flux.RUSANOV
order = 3
well_balanced = True
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