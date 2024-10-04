from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from nummeth import NumericalMethod
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.SW
init = InitCond.READ_FROM_FILE
funh = FunH.BUMPS
H_noise_factor = 0
perturb_init = InitCond.PERT_NONE#WB
boundary = BoundaryCond.IN_OUT
nummeth = NumericalMethod.UPWIND
timest = TimeStepping.TVDRK3
order = 3
N = 100
a = 0
b = 25
cfl = 0.5
T = 20.0
steps = 4
ode = 'AM'
system ='SW'
plot_exact = False
plot_every = 200.5
show_plots = True
save_plots = False
save_npys = False



print "Loaded config/swe_rest.py!"
