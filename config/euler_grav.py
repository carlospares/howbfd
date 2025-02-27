from initcond import InitCond
from equation import Equation
from functionH import FunH
from eq_sw import SWEquation
from boundary import BoundaryCond
from nummeth import NumericalMethod
from timest import TimeStepping

# For a detailed explanation, see howbfd_config

equation = Equation.EulerGRAV
init = InitCond.Eulerisothermal#DISCRETE_AB#STEADY#READ_FROM_FILE#STEADY #WATER_AT_REST
funh = FunH.EUL_ISOTHERMAL
H_noise_factor = 0
perturb_init = InitCond.PERT_NONE#DISC#PERT_RIEMANN#DISC#PERT_WB #PERT_NONE
boundary = BoundaryCond.LIN_EXTRAP#IN_OUT #FORCE_STEADY_INIT
nummeth = NumericalMethod.RUSANOV#WB #UPWINDGF
timest = TimeStepping.TVDRK3#TVDRK3#EULER#TVDRK3
order = 3
N = 400
a = 0
b = 1
cfl = 0.4
T = 2.0 # 0.1644
steps=4
ode='AM'
system='SW'
compute_source='analytic_source'#'hydrostatic_reconstruction'#'source_reconstruction'#'analytic_source'#'hydrostatic_reconstruction' 
plot_exact = False
plot_every = 1000.0
show_plots = True
save_plots = False
save_npys = False



print ("Loaded config/euler_grav.py!")
