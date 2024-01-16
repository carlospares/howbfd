from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from nummeth import NumericalMethod
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.BURGERS
init = InitCond.STEADY#MMSburg
funh = FunH.IDpSIN#MMSburg#IDENT
H_noise_factor = 0.0
perturb_init = InitCond.PERT_NONE#GAUSS
boundary = BoundaryCond.FORCE_STEADY_INIT
nummeth = NumericalMethod.UPWIND#RUSANOV #UPWIND
timest = TimeStepping.TVDRK3
order = 3
N = 30
cfl = 0.5
a = -1#0
b = 1#15
T = 1.0#2.0
plot_every = 0.1
show_plots = True
save_plots = False
save_npys = False
plot_exact = True

print "Loaded config/burgers_rusanov.py!"
