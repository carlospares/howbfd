from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from nummeth import NumericalMethod
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.SW
init = InitCond.WATER_AT_REST#STEADY
funh = FunH.BUMP#S#STEP#BUMPS#D# BUMP
H_noise_factor = 0
perturb_init = InitCond.PERT_NONE#PERT_DISC#GAUSS #PERT_NONE
boundary = BoundaryCond.WALL#IN_OUT #FORCE_STEADY_INIT
nummeth = NumericalMethod.UPWINDGF#RUSANOVGWB #UPWINDGF
timest = TimeStepping.TVDRK3
order = 3
N = 100
a = 0#-3 #0
b = 25 #3 #25
cfl = 0.9
T = 10.#0.1
steps = 4
ode = 'AM'
system = 'SW'
plot_exact = False
plot_every = 300.0
show_plots = True
save_plots = False
save_npys = False



print "Loaded config/swe_rest.py!"
