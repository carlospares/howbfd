import numpy as np
import matplotlib.pyplot as plt

# Loading data
w3am4 = np.loadtxt('Riemann/Weno3_AM4_N100_T0p7.dat')
w3am6 = np.loadtxt('Riemann/Weno3_AM6_N100_T0p7.dat')
w3am8 = np.loadtxt('Riemann/Weno3_AM8_N100_T0p7.dat')
w3up = np.loadtxt('Riemann/Weno3_Upwind_N100_T0p7.dat')
#ref = np.loadtxt('Riemann/reference_upwind_N5000_T0p7.dat')
ref = np.loadtxt('Riemann/reference_w5_am6_N800_T0p7.dat')
analytic4 = np.loadtxt('../Weno3_AM4/out100')
analytic6 = np.loadtxt('../Weno3_AM6/out100')
analytic8 = np.loadtxt('../Weno3_AM8/out100.txt')
#analytic_ref = np.loadtxt('../../../initial_data/analytical_sw/supercritical/initial_sup_5000.dat')
analytic_ref = np.loadtxt('../Weno7_AM8/out800.txt')

#w3am4 = np.loadtxt('DISC_0p0001_from_analytical/Weno3_AM4_N100_T1p0.txt')
#w3am6 = np.loadtxt('DISC_0p0001_from_analytical/Weno3_AM6_N100_T1p0.txt')
#w3am8 = np.loadtxt('DISC_0p0001_from_analytical/Weno3_AM8_N100_T1p0.txt')
#w3up = np.loadtxt('DISC_0p0001_from_analytical/Weno3_Upwind_N100_T1p0.txt')
#ref = np.loadtxt('DISC_0p0001_from_analytical/reference_w7_am8_T1p0.txt')
#analytic = np.loadtxt('../../../initial_data/analytical_sw/subcritical/initial_sub_100.dat')
#analytic_ref = np.loadtxt('../../../initial_data/analytical_sw/subcritical/initial_sub_800.dat')

# First figure
plt.figure()
#pert = w3am4[:, 1] - analytic
pert4 = w3am4[:, 1] - analytic4[:,3]
pert6 = w3am6[:, 1] - analytic6[:,3]
pert8 = w3am8[:, 1] - analytic8[:,3]
pert_ref = ref[:, 1] - analytic_ref[:,3]

# Plot the perturbation (difference between the numerical solution and analytical)
plt.plot(w3am4[:, 0], w3am4[:, 1] - w3am4[:, 5], 'b', linewidth=2, label='eta')
plt.plot(w3am4[:, 0], -w3am4[:, 5], 'k', linewidth=2,label='-H(x)')

plt.xlabel(r'$x$', fontsize=18)
plt.ylabel(r'$\eta$', fontsize=18)
# Zoom-in region
#xlim_zoom = [1.5, 2.5]
#ylim_zoom = [np.sin(1.5), np.sin(2.5)]

# Create the zoomed-in axes
#zoom_axes = plt.axes([0.5, 0.3, 0.35, 0.35])  # Position in normalized units [left, bottom, width, height]
#plt.box(on=True)

## Plot the zoomed-in region in the new axes
#plt.plot(w3am4[:, 0], pert4, linewidth=1.5)
#plt.xlabel('Zoomed x')
#plt.ylabel(r'$\eta$')
#plt.title('Perturbation')
#
## Set limits for the zoomed-in plot
#plt.xlim([6, 11])
## You can set the y-limits manually if needed
## plt.ylim(ylim_zoom)

plt.tight_layout()
#plt.show()
#plt.savefig('initial_subcritical_pert.png')

# Second figure
plt.figure()

# Subplot 1: comparison of AM4, AM6, AM8
plt.subplot(2, 1, 1)
plt.plot(w3up[:, 0], w3up[:, 3] - w3up[:, 1] + pert4,'--', linewidth=2, label='WENO3',color='slategray')
plt.plot(ref[:, 0], ref[:, 3] - ref[:, 1] + pert_ref, 'r', linewidth=2.5, label='reference')
plt.plot(w3am4[:, 0], w3am4[:, 3] - w3am4[:, 1] + pert4, '-ob', linewidth=2, label='GF-AM4',markersize=3)
plt.plot(w3am6[:, 0], w3am6[:, 3] - w3am6[:, 1] + pert6, '--k', linewidth=2, label='GF-AM6')
plt.plot(w3am8[:, 0], w3am8[:, 3] - w3am8[:, 1] + pert8, '-.c', linewidth=2, label='GF-AM8')

plt.legend(loc='lower left',fontsize=11,framealpha=0.5)
plt.xlabel(r'$x$', fontsize=18)
plt.ylabel(r'$\eta$', fontsize=18)
plt.grid()
plt.axis([0, 25, -0.5, 1.25])

# Subplot 2: Comparison of AM4, AM6, AM8 for 'q' values
plt.subplot(2, 1, 2)
plt.plot(w3up[:, 0], w3up[:, 4], '--', linewidth=2, label='WENO3',color='slategray')
plt.plot(ref[:, 0], ref[:, 4], 'r', linewidth=2.5, label='reference')
plt.plot(w3am4[:, 0], w3am4[:, 4], '-ob', linewidth=2, label='GF-AM4',markersize=3)
plt.plot(w3am6[:, 0], w3am6[:, 4], '--k', linewidth=2, label='GF-AM6')
plt.plot(w3am8[:, 0], w3am8[:, 4], '-.c', linewidth=2, label='GF-AM8')
plt.axis([0, 25, 16, 24.5])

plt.legend(loc='lower left',fontsize=11,framealpha=0.5)
plt.xlabel(r'$x$', fontsize=18)
plt.ylabel(r'$q$', fontsize=18)
plt.grid()

plt.tight_layout()
#plt.show()
#plt.savefig('subcritical_pert_AM.png')

## Third figure: Upwind scheme
#plt.figure()
#
## Subplot 1: eta comparison for Upwind scheme
#plt.subplot(2, 1, 1)
#plt.plot(w3up[:, 0], w3up[:, 3] - w3up[:, 1] + pert4, 'b', linewidth=2, label='Upwind')
#plt.legend(loc='upper right')
#plt.xlabel(r'$x$', fontsize=18)
#plt.ylabel(r'$\eta$', fontsize=18)
#plt.grid()
#
## Subplot 2: q values for Upwind scheme
#plt.subplot(2, 1, 2)
#plt.plot(w3up[:, 0], w3up[:, 4], 'b', linewidth=2, label='Upwind')
#plt.legend(loc='upper right')
#plt.xlabel(r'$x$', fontsize=18)
#plt.ylabel(r'$q$', fontsize=18)
#plt.grid()
#
#plt.tight_layout()
plt.show()
##plt.savefig('subcritical_pert_Upwind.png',dpi=100)

