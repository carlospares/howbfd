import numpy as np
import matplotlib.pyplot as plt

# Loading data
w3am4 = np.loadtxt('weno5_am6_N100_riemann_nm0p05_T1.txt')
w3up = np.loadtxt('weno5_upwind_N100_riemann_nm0p05_T1.txt')
#w3up400 = np.loadtxt('weno5_upwind_N400_pert_T1.txt')
ref = np.loadtxt('reference_w7_am8_N800_riemann_nm0p05_T1.txt')
analytic4 = np.loadtxt('weno5_am6_N100_steady.txt')
analyticup = np.loadtxt('weno5_upwind_N100_steady.txt')
#analyticup400 = np.loadtxt('weno5_upwind_N400_steady.txt')
analytic_ref = np.loadtxt('reference_w7_am8_N800_steady.txt')


# First figure
plt.figure()
pert4 = w3am4[:, 1] - analytic4[:,3]
#pert6 = w3am6[:, 1] - analytic6[:,3]
#pert8 = w3am8[:, 1] - analytic8[:,3]
pertup = w3up[:, 1] - analyticup[:,3]
#pertup400 = w3up400[:, 1] - analyticup400[:,1]#3
pert_ref = ref[:, 1] - analytic_ref[:,3]

# Plot the perturbation (difference between the numerical solution and analytical)
plt.plot(w3am4[:, 0], w3am4[:, 1] - w3am4[:, 5], 'b', linewidth=2, label='$\eta$')
plt.plot(w3am4[:, 0], -w3am4[:, 5], 'k', linewidth=2,label='$-H(x)$')
plt.legend(loc='upper right',fontsize=11,framealpha=0.5,ncol=2,bbox_to_anchor=(1,0.95))
plt.xlabel(r'$x$', fontsize=18)
plt.axis([0, 25, -0.1, 3.2])

# Zoom-in region
#xlim_zoom = [1.5, 2.5]
#ylim_zoom = [np.sin(1.5), np.sin(2.5)]

# Create the zoomed-in axes
#zoom_axes = plt.axes([0.5, 0.3, 0.35, 0.35])  # Position in normalized units [left, bottom, width, height]
#plt.box(on=True)

# Plot the zoomed-in region in the new axes
#plt.plot(w3am4[:, 0], pert4, linewidth=1.5)
##plt.xlabel('Zoomed x')
#plt.ylabel(r'$\eta$')
#plt.title('Perturbation')

## Set limits for the zoomed-in plot
#plt.xlim([6, 11])
# You can set the y-limits manually if needed
# plt.ylim(ylim_zoom)

plt.tight_layout()
#plt.show()
#plt.savefig('initial_subcritical_pert.png')

# Second figure
plt.figure()

# Subplot 1: comparison of AM4, AM6, AM8
plt.subplot(2, 1, 1)
plt.plot(w3up[:, 0], w3up[:, 3] - w3up[:, 1] + pertup,'--', linewidth=2, label='WENO5',color='slategray')
#plt.plot(w3up400[:, 0], w3up400[:, 3] - w3up400[:, 1] + pertup400,'-b', linewidth=2, label='WENO5 (N=400)')
plt.plot(ref[:, 0], ref[:, 3] - ref[:, 1] + pert_ref, 'r', linewidth=2.5, label='reference')
#plt.plot(w3am4[:, 0], w3am4[:, 3] - w3am4[:, 1] + pert4, '-ob', linewidth=2, label='GF-AM4',markersize=3)
plt.plot(w3am4[:, 0], w3am4[:, 3] - w3am4[:, 1] + pert4, 'k', linewidth=2, label='GF-AM6')
#plt.plot(w3am6[:, 0], w3am6[:, 3] - w3am6[:, 1] + pert6, '--k', linewidth=2, label='GF-AM6')
#plt.plot(w3am8[:, 0], w3am8[:, 3] - w3am8[:, 1] + pert8, '-.c', linewidth=2, label='GF-AM8')

plt.legend(loc='upper right',fontsize=11,framealpha=1.0,ncol=1)
plt.xlabel(r'$x$', fontsize=18)
plt.ylabel(r'$\eta-\eta_s$', fontsize=18)
plt.grid()
plt.axis([0, 25, -0.1, 1.25])

# Subplot 2: Comparison of AM4, AM6, AM8 for 'q' values
plt.subplot(2, 1, 2)
plt.plot(w3up[:, 0], w3up[:, 4]-w3up[:,2], '--', linewidth=2, label='WENO5',color='slategray')
#plt.plot(w3up400[:, 0], w3up400[:, 4], '-b', linewidth=2, label='WENO5 (N=400)')
plt.plot(ref[:, 0], ref[:, 4]-ref[:,2], 'r', linewidth=2.5, label='reference')
plt.plot(w3am4[:, 0], w3am4[:, 4]-w3am4[:,2], '-k', linewidth=2, label='GF-AM6')
#plt.plot(w3am6[:, 0], w3am6[:, 4], '--k', linewidth=2, label='GF-AM6')
#plt.plot(w3am8[:, 0], w3am8[:, 4], '-.c', linewidth=2, label='GF-AM8')
plt.axis([0, 25, -0.1, 3.])

#plt.legend(loc='upper right',fontsize=11,framealpha=0.5,ncol=2)
plt.xlabel(r'$x$', fontsize=18)
plt.ylabel(r'$q-q_s$', fontsize=18)
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

