import numpy as np
import matplotlib.pyplot as plt

# Loading data
#-------pertrubation disc-----
w3 =    np.loadtxt('Burgers/pert_disc_steady/pert_W3_T0p8.txt')
w3nwb = np.loadtxt('Burgers/pert_disc_steady/pert_W3NWB_T0p8.txt')
ref =   np.loadtxt('Burgers/pert_disc_steady/reference_W9_T0p8.txt')



# First figure

plt.figure()

plt.plot(w3[:, 0], w3[:, 1], 'b', linewidth=2, label='initial state')
plt.plot(w3[:, 0], -w3[:, 3], 'k', linewidth=2,label='H(x)')
plt.legend(loc='upper right',fontsize=11,framealpha=0.5)
plt.xlabel(r'$x$', fontsize=18)
plt.ylabel(r'$U$', fontsize=18)
plt.grid()

# Create the zoomed-in axes
#zoom_axes = plt.axes([0.5, 0.3, 0.35, 0.35])  # Position in normalized units [left, bottom, width, height]
#plt.box(on=True)

#pert = w3[:, 1] - st[:,3]
#pert_ref = ref[:, 1] - st_ref[:,3]
#pert_nwb = w3nwb[:, 1] - st_nwb[:,3]
## Plot the zoomed-in region in the new axes
##plt.plot(w3[:, 0], w3[:,1]-w3[:,5], linewidth=1.5)
#plt.plot(w3[:, 0], w3[:,3]-w3[:,1]+pert, linewidth=1.5)
#plt.xlabel('Zoomed x')
#plt.ylabel(r'$\eta$')
#plt.title('Perturbation')

# Set limits for the zoomed-in plot
#plt.xlim([6, 11])
## You can set the y-limits manually if needed
# plt.ylim(ylim_zoom)

#Second figure
plt.figure()
#pert4 = w3am4[:, 1] - analytic4[:,3]
#pert6 = w3am6[:, 1] - analytic6[:,3]
#pert8 = w3am8[:, 1] - analytic8[:,3]
#pertup = w3up[:, 1] - analyticup[:,3]
#pert_ref = ref[:, 1] - analytic_ref[:,3]

# Plot the perturbation (difference between the numerical solution and analytical)
plt.plot(ref[:, 0], ref[:, 2], '--r', linewidth=3, label='reference')
plt.plot(w3nwb[:, 0], w3nwb[:, 2] ,'--', linewidth=2, label='WENO3',color='slategray')
plt.plot(w3[:, 0], w3[:, 2], 'b', linewidth=2, label='WBWENO3')
#plt.plot(w3[:, 0], w3[:, 3], 'k', linewidth=2,label='H(x)')
plt.legend(loc='upper right',fontsize=11,framealpha=0.5)
plt.xlabel(r'$x$', fontsize=18)
plt.ylabel(r'$U$', fontsize=18)
plt.grid()
#plt.axis([0, 25, 2-0.001, 2.001])
#plt.axis([0, 25, -0.00004, 0.0001])


plt.tight_layout()
#plt.show()
#plt.savefig('initial_subcritical_pert.png')

# Third figure
#
## Subplot 1:
#plt.subplot(2, 1, 1)
#plt.plot(w3[:, 0], w3[:, 3]-w3[:,5], 'b', linewidth=2, label='surface')
#plt.plot(w3[:, 0], -w3[:, 5], 'k', linewidth=2,label='bottom')
##plt.plot(w3[:, 0], w3[:, 3] - w3[:, 1] + pertup,'--', linewidth=2, label='WENO3',color='slategray')
###plt.plot(ref[:, 0], ref[:, 3] - ref[:, 1] + pert_ref, 'r', linewidth=2.5, label='reference')
##plt.plot(w3am4[:, 0], w3am4[:, 3] - w3am4[:, 1] + pert4, '-ob', linewidth=2, label='GF-AM4',markersize=3)
##plt.plot(w3am6[:, 0], w3am6[:, 3] - w3am6[:, 1] + pert6, '--k', linewidth=2, label='GF-AM6')
##plt.plot(w3am8[:, 0], w3am8[:, 3] - w3am8[:, 1] + pert8, '-.c', linewidth=2, label='GF-AM8')
##
#plt.axis([0, 25, -0.1, 2.2])
##plt.legend(loc='upper right',fontsize=11,framealpha=0.5,bbox_to_anchor=(1, 0.9))
##bbox_to_anchor=(1, 0.9),   # move it slightly downward
#plt.grid()
#
## Subplot 2: 
#plt.subplot(2, 1, 2)
#plt.plot(w3[:, 0], w3[:, 4], 'r', linewidth=2,label='q')
##plt.plot(w3up[:, 0], w3up[:, 4], '--', linewidth=2, label='WENO3',color='slategray')
###plt.plot(ref[:, 0], ref[:, 4], 'r', linewidth=2.5, label='reference')
###plt.plot(w3am4[:, 0], w3am4[:, 4], '-ob', linewidth=2, label='GF-AM4',markersize=3)
##plt.plot(w3am6[:, 0], w3am6[:, 4], '--k', linewidth=2, label='GF-AM6')
###plt.plot(w3am8[:, 0], w3am8[:, 4], '-.c', linewidth=2, label='GF-AM8')
##plt.axis([0, 25, -23.995, -24.005])
#plt.axis([0, 25, 4.415, 4.425])
#
#plt.legend(loc='upper right',fontsize=11,framealpha=0.5)
#plt.xlabel(r'$x$', fontsize=18)
##plt.ylabel(r'$q$', fontsize=18)
#plt.grid()

##plt.tight_layout()
##plt.show()
##plt.savefig('subcritical_pert_AM.png')

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

