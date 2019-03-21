from initcond import InitCond
from equation import Equation
from eq_sw import SWEquation
from boundary import BoundaryCond
from numflux import Flux

# For a detailed explanation, see howbfd_config

equation = Equation.SW
init = InitCond.FLAT
sw_H = SWEquation.H_FLAT
sw_H_noise_factor = 0.0
perturb_init = InitCond.PERT_MGAUSS
boundary = BoundaryCond.FORCE_STEADY
numflux = Flux.RUSANOV
order = 3
well_balanced = True
is_conservative = True
N = 100
cfl = 0.5
a = -1
b = 1
T = 0.5
plot_every = 0.02
show_plots = True
save_plots = True
save_npys = False
print "Loaded config/swe_pond.py!"