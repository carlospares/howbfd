from initcond import InitCond
from equation import Equation
from boundary import BoundaryCond
from numflux import Flux

########################################
# DEFAULT OPTIONS!
# To use one of the premade files, call
# python main.py -c <any file in /config/ >
########################################
# See detailed explanation below
equation = Equation.SWE_REST
init = InitCond.PWPOLY
perturb_init = InitCond.PERT_NONE
sw_H = Equation.SWE_H_PWPOLY
sw_H_noise_factor = 0.1
boundary = BoundaryCond.FORCE_STEADY
numflux = Flux.RUSANOV
order = 3
well_balanced = True
N = 100
a = 0
b = 20
cfl = 0.5
T = 8
plot_every = 0.5
show_plots = True
save_plots = False
save_npys = False


# equation: decides the equation to be solved. It is one of:
# Equation.
#	LINEAR: 1D linear transport, u_t + alpha u_x = u (alpha constant)
#	BURGERS: 1D Burgers equation, u_t + u u_x = u^2
#	SWE_REST: 1D shallow water equation, variables [h, q=hu]
# See equation.py

# init: decides the initial condition of the problem. It is one of:
# InitCond.
#	SIN: 1 + sin(2*pi*x) [if more than one variable, others will be zero]
#	SHOCK: 1.0*(x <= 0.5) + 2.0 * (x>0.5) [same]
#	STEADY: arbitrary steady solution for equation defined in equation.steady
#	PWPOLY: 1+(0.13+0.05*(x-10)**2) in [8, 12]; 0.33 otherwise

# perturb_init: chooses what perturbation to add to the initial condition.
# Meant for perturbations to InitCond.STEADY, but works for all others.
# If there is more than one variable, will only apply it to the first.
# It is one of:
# InitCond.
#	PERT_NONE: no perturbation
#	PERT_POLY: 0.01*x*(1-x)
#	PERT_PATCH: 1 in [-0.5, -0.3], 0 otherwise
#	PERT_SIN: 0.01*(C+ sin(2*pi*x)), C a constant config. in initcond.py
#	PERT_GAUSS: 0.3*exp(-200*(x+0.5)**2)

# sw_H and sw_H_noise_factor decide the bathymetry profile of SWE
# unused for other equations at the moment
# sw_H is one of:
# Equation.SWE_H_
#	FLAT: 0.1
#	PWPOLY: (0.13+0.05*(x-10)**2) in [8,12]; 0.33 otherwise
# sw_H_noise_factor adds a U[0,1] perturbation to H; so the end result is
# H = sw_H + noise_factor*U[0,1] (0 for no noise)

# boundary: sets the boundary conditions to use. It must be one of 
# BoundaryCond.
#	PERIODIC: ghost cells will repeat values at opposite end
#   IN_OUT: inflow (as in initial condition) on left, homogeneous Neumann on right
#   FORCE_STEADY: make BCs be values for u0(x); might break if u0 not steady state.
# x (the spatial coordinates) will be expanded with periodic BCs if PERIODIC,
# linearly extrapolated otherwise (to x[N]+dx, x[N]+2dx... at right end, " at left)

# numflux sets the numerical flux. It must be one of:
# Flux.
# 	UPWIND: upwind flux, SCALAR ONLY
#	RUSANOV: Rusanov for scalars or systems

# order: order of WENO reconstruction (one of 3, 5, 7, 9, 11)

# well_balanced:
#	True: use high order well-balanced version
#	False: use plain WENO reconstruction + source terms for SH_x (requires H_x implemented)

# Spatial parameters:
#	N: number of spatial points
#	[a,b]: interval to solve
# Points are cell centers; ie x[0] = a+0.5dx

# Temporal parameters:
#	cfl: cfl number
#	T: time horizon of the simulation

# IO options:
#	plot_every: Every (this many) seconds, do...
#	show_plots: ...display the plot
#	save_plots: ...save the plot in /figs
#	save_npys:  ...save numpy array in /npys