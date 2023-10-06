from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from nummeth import NumericalMethod
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.BURGERS
init = InitCond.STEADY
funh = FunH.IDENT
H_noise_factor = 0.0
perturb_init = InitCond.PERT_GAUSS
boundary = BoundaryCond.FORCE_STEADY_INIT
nummeth = NumericalMethod.UPWINDGF
timest = TimeStepping.TVDRK3
order = 3
N = 100
cfl = 0.5
a = -1
b = 1
T = 1.0#0.0037157669102204603
plot_every = 0.25
show_plots = True
save_plots = True
save_npys = False
plot_exact = True

print "Loaded config/burgers_upwind_gf.py!"
