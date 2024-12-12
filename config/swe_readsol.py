from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from nummeth import NumericalMethod
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.SW
init = InitCond.WATER_AT_REST#DISCRETE_AB#STEADY#READ_FROM_FILE#STEADY #WATER_AT_REST
funh = FunH.STEP#BUMPS#BUMP2#D# BUMP
H_noise_factor = 0
perturb_init = InitCond.PERT_NONE#PERT_RIEMANN#DISC#PERT_WB #PERT_NONE
boundary = BoundaryCond.SUBCR#FORCE_DISCRETE_STEADY_INIT#SUPER#_RE#SUBCR#SUBCR_RE#IN_OUT #FORCE_STEADY_INIT
nummeth = NumericalMethod.UPWINDGF#RUSANOVGWB #UPWINDGF
timest = TimeStepping.TVDRK3#EULER#TVDRK3
order = 3
N = 25
a = 0#-3 #0
b = 25 #3 #25
cfl = 0.6
T = 300.0
steps=4
ode='AB'
system='SW'
compute_source='analytic_source' #analytic_source, source_reconstruction,hydraustatic_reconstruction
plot_exact = False
plot_every = 300.0
show_plots = True
save_plots = False
save_npys = False



print ("Loaded config/swe_readsol.py!")
