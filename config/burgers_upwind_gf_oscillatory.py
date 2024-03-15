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
funh = FunH.IDpSIN#IDENT#IDpSIN
H_noise_factor = 0.0
perturb_init = InitCond.PERT_NONE#DISC#NONE#PERT_GAUSS
boundary = BoundaryCond.FORCE_STEADY_INIT
nummeth = NumericalMethod.UPWIND#GF
timest = TimeStepping.TVDRK3
order = 5
N = 150
cfl = 0.5
a = -1
b = 1
T = 0.7#0.0037157669102204603
plot_every = 2.1
show_plots = True
save_plots = False
save_npys = False
plot_exact = True

print "Loaded config/burgers_upwind_gf.py!"
