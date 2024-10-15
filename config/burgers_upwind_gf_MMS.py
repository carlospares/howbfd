from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from nummeth import NumericalMethod
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.BURGERS
init = InitCond.MMSburg
funh = FunH.MMSburg
H_noise_factor = 0.0
perturb_init = InitCond.PERT_NONE
boundary = BoundaryCond.FORCE_STEADY_INIT
nummeth = NumericalMethod.UPWINDGF
timest = TimeStepping.TVDRK3
order = 3
N = 30
cfl = 0.5
a = 0
b = 15
T = 2.0#0.0037157669102204603
steps=4
ode='AB'
system='No'
plot_every = 2.5
show_plots = True
save_plots = False
save_npys = False
plot_exact = True

print ("Loaded config/burgers_upwind_gf.py!")
