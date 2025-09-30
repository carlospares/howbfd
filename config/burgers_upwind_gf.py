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
nummeth = NumericalMethod.UPWINDNCGF#WB#GF
timest = TimeStepping.TVDRK3
order = 3 
N = 20
a = -1
b = 1
cfl = 0.5
T = 5.0 #0.2#0.0037157669102204603
steps = 4
ode = 'AM'
system = 'No'
compute_source= 'analytic_source'#'hydrostatic_reconstruction'#'source_reconstruction'#'analytic_source'#'hydrostatic_reconstruction' 
plot_every = 1.0
show_plots = True 
save_plots = False
save_npys = False
plot_exact = False

print ("Loaded config/burgers_upwind_gf.py!")
