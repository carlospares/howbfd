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
funh = FunH.STEP#BUMPS#STEP#BUMP2#D# BUMPS
H_noise_factor = 0
perturb_init = InitCond.PERT_NONE#DISC#GAUSS#DISC#PERT_RIEMANN#DISC#PERT_WB #PERT_NONE
boundary = BoundaryCond.WALL#FORCE_STEADY_INIT#FORCE_DISCRETE_STEADY_INIT#SUPER_RE#FORCE_STEADY_INIT#FORCE_DISCRETE_STEADY_INIT#SUBCR#SUPER#FORCE_DISCRETE_STEADY_INIT#SUPER#_RE#SUBCR#SUBCR_RE#IN_OUT #FORCE_STEADY_INIT
nummeth = NumericalMethod.UPWINDNCGF#RUSANOVGWB #UPWINDGF
timest = TimeStepping.TVDRK3
order = 7
N = 25
a = 0#-3 #0
b = 25#3 #25
cfl = 0.4
T = 1.0
steps=4
ode='AB'
system='SW'
compute_source='analytic_source'#'source_reconstruction'#'analytic_source'#'hydrostatic_reconstruction' 
plot_exact = False
plot_every = 1.5
show_plots = False
save_plots = False
save_npys = False



print ("Loaded config/swe_readsol.py!")
