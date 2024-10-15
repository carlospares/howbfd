from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from nummeth import NumericalMethod
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.EulerSEC
init = InitCond.Eulersec#READ_FROM_FILE#STEADY #WATER_AT_REST
funh = FunH.EUL_SEC#BUMP2#BUMPS#D# BUMP
H_noise_factor = 0
perturb_init = InitCond.PERT_NONE#DISC#NONE#PERT_WB #PERT_NONE
boundary = BoundaryCond.SUBCR#SUBCR_RE#IN_OUT #FORCE_STEADY_INIT
nummeth = NumericalMethod.UPWINDGF#RUSANOVGWB #UPWINDGF
timest = TimeStepping.TVDRK3
order = 3
N = 50
a = 0#-3 #0
b = 25 #3 #25
cfl = 0.5
T = 250.0
steps=4
ode='AB'
system='SW'
plot_exact = False
plot_every = 1.0
show_plots = True
save_plots = False
save_npys = False



print ("Loaded config/swe_rest.py!")