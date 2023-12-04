from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from nummeth import NumericalMethod
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.BURGERS
init = InitCond.MMSburg#STEADY
funh = FunH.MMSburg#IDENT
H_noise_factor = 0.0
perturb_init = InitCond.PERT_NONE #PERT_GAUSS
boundary = BoundaryCond.FORCE_STEADY_INIT
nummeth = NumericalMethod.UPWIND #RUSANOV #UPWIND
timest = TimeStepping.TVDRK3
order = 9
N = 30
cfl = 0.01
a = 0#-1
b = 15#1
T = 2.0#8
plot_every = 2.5
show_plots = True
save_plots = False
save_npys = False
plot_exact = True

print "Loaded config/burgers_rusanov.py!"
