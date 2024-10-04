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
funh = FunH.DISC#PAR#DISC#IDENT
H_noise_factor = 0.0
perturb_init = InitCond.PERT_GAUSS#NONE
boundary = BoundaryCond.FORCE_STEADY_INIT
nummeth = NumericalMethod.UPWINDGF#WB#GF
timest = TimeStepping.TVDRK3
order = 3 
N = 20
cfl = 0.3
a = -1
b = 1
T = 0.8#0.0037157669102204603
steps=4
ode='AM'
system='No'
plot_every = 0.1
show_plots = True 
save_plots = False
save_npys = True
plot_exact = True

print "Loaded config/burgers_upwind_gf.py!"
