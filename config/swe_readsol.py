from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from nummeth import NumericalMethod
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.SW
init = InitCond.STEADY#READ_FROM_FILE#STEADY #WATER_AT_REST
funh = FunH.BUMPS#STEP#BUMP2#D# BUMP
H_noise_factor = 0
perturb_init = InitCond.PERT_NONE#PERT_RIEMANN#DISC#PERT_WB #PERT_NONE
boundary = BoundaryCond.SUBCR#_RE#SUBCR#SUBCR_RE#IN_OUT #FORCE_STEADY_INIT
nummeth = NumericalMethod.UPWINDGF#RUSANOVGWB #UPWINDGF
timest = TimeStepping.TVDRK3#EULER#TVDRK3
order = 3
N = 25
a = 0#-3 #0
b = 25 #3 #25
cfl = 0.5
T = 30.
steps=8
ode='AM'
system='SW'
plot_exact = False
plot_every = 100.0
show_plots = False
save_plots = False
save_npys = False



print ("Loaded config/swe_readsol.py!")
