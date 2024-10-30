import numpy as np
import matplotlib.pyplot as plt

# Loading data
w3am4 = np.loadtxt('weno3_am4_T0p2_N100.txt')
w3am6 = np.loadtxt('weno3_am6_T0p2_N100.txt')
w3am8 = np.loadtxt('weno3_am8_T0p2_N100.txt')
ref = np.loadtxt('reference_weno3_am4_T0p2_N1000.txt')
w3up = np.loadtxt('weno3_upwind_T0p2_N100.txt')


# First figure
plt.figure()


plt.plot(ref[:, 0], ref[:, 2]   ,'-', linewidth=2, label='reference ',color='r')
plt.plot(w3up[:, 0], w3up[:, 2]   ,'--', linewidth=2, label='WENO3-nWB',color='slategray')
plt.plot(w3am4[:, 0], w3am4[:, 2] , '-ob', linewidth=2, label='WENO3GF-AM4 ',markersize=3)
plt.plot(w3am6[:, 0], w3am6[:, 2] , '--k', linewidth=2, label='WENO3GF-AM6 ')
plt.plot(w3am8[:, 0], w3am8[:, 2] , '-.c', linewidth=2, label='WENO3GF-AM8 ')

plt.legend(loc='upper left',fontsize=11,framealpha=0.5)
plt.xlabel(r'$x$', fontsize=18)
plt.ylabel(r'$u$', fontsize=18)
plt.grid()
#plt.axis([0, 30, -0.0001, 0.0001])

plt.tight_layout()
plt.show()

