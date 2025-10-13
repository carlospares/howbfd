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
nummeth = NumericalMethod.UPWINDNCGF
timest = TimeStepping.TVDRK3
order = 5
N = 30
cfl = 0.05
a = 0
b = 15
T = 1.0#0.0037157669102204603
steps=4
ode='AB'
system='No'
compute_source= 'analytic_source'#'hydrostatic_reconstruction'#'source_reconstruction'#'analytic_source'#'hydrostatic_reconstruction' 
plot_every = 2.5
show_plots = False
save_plots = False
save_npys = False
plot_exact = True

print ("Loaded config/burgers_upwind_gf.py!")
