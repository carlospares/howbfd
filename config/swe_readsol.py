from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from nummeth import NumericalMethod
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.SW
init = InitCond.STEADY #READ_FROM_FILE#STEADY #WATER_AT_REST
funh = FunH.BUMP2#BUMPS#D# BUMP
H_noise_factor = 0
perturb_init = InitCond.PERT_NONE#DISC#NONE#PERT_WB #PERT_NONE
boundary = BoundaryCond.SUBCR#SUBCR_RE#IN_OUT #FORCE_STEADY_INIT
nummeth = NumericalMethod.UPWINDGF#RUSANOVGWB #UPWINDGF
timest = TimeStepping.TVDRK2#EULER#TVDRK3
order = 5
N = 401
a = 0#-3 #0
b = 25 #3 #25
cfl = 0.6
T = 150.0
steps=4
ode='AM'
system='SW'
plot_exact = False
plot_every = 25.0
show_plots = True
save_plots = False
save_npys = False



print ("Loaded config/swe_rest.py!")
