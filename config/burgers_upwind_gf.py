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
funh = FunH.IDENT#DISC#PAR#DISC#IDENT
H_noise_factor = 0.0
perturb_init = InitCond.PERT_NONE#GAUSS#NONE
boundary = BoundaryCond.FORCE_STEADY_INIT
nummeth = NumericalMethod.UPWIND#GF#WB#GF
timest = TimeStepping.TVDRK3
order = 7 
N = 20
a = -1
b = 1
cfl = 0.75
T = 5.0#0.2#0.0037157669102204603
steps=8
ode='AM'
system='No'
compute_source='analytic_source'#'hydrostatic_reconstruction'#'source_reconstruction'#'analytic_source'#'hydrostatic_reconstruction' 
plot_every = 10.05
show_plots = False 
save_plots = False
save_npys = False
plot_exact = False

print ("Loaded config/burgers_upwind_gf.py!")
