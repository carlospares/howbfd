from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from nummeth import NumericalMethod
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.SW
init = InitCond.WATER_AT_REST#STEADY #WATER_AT_REST
funh = FunH.STEP#BUMPS#D# BUMP
H_noise_factor = 0
perturb_init = InitCond.PERT_DISC#GAUSS #PERT_NONE
boundary = BoundaryCond.WALL#IN_OUT #FORCE_STEADY_INIT
nummeth = NumericalMethod.UPWINDGF#RUSANOVGWB #UPWINDGF
timest = TimeStepping.TVDRK3
order = 3
N = 100
a = 0#-3 #0
b = 25 #3 #25
cfl = 0.5
T = 10.5#0.1
plot_exact = False
plot_every = 0.5
show_plots = True
save_plots = False
save_npys = True



print "Loaded config/swe_rest.py!"
