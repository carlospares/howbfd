import numpy as np
import matplotlib.pyplot as plt

# Loading data
am4 = np.loadtxt('AM4.txt')
ab4 = np.loadtxt('AB4.txt')
am6 = np.loadtxt('AM6.txt')
ab6 = np.loadtxt('AB6.txt')
am8 = np.loadtxt('AM8.txt')
ab8 = np.loadtxt('AB8.txt')
w3nwb = np.loadtxt('W3NB.txt')

plt.figure()

# Subplot 1: comparison of AM4, AM6, AM8
plt.plot((w3nwb[:, 6]), (w3nwb[:, 4]) ,'-r', linewidth=2, label='W5-nWB')
plt.plot((am4[:, 6]), (am4[:, 4]) ,'-^', linewidth=2, label='W5GF-AM4',color='slategray')
plt.plot((ab4[:, 6]), (ab4[:, 4]) ,'-vk', linewidth=2, label='W5GF-AB4')
plt.plot((am6[:, 6]), (am6[:, 4]) ,'-*', linewidth=2, label='W5GF-AM6',color='cyan')
plt.plot((ab6[:, 6]), (ab6[:, 4]) ,'--', linewidth=2, label='W5GF-AB6',color='#FF5733')
plt.plot((am8[:, 6]), (am8[:, 4]) ,'-v', linewidth=2, label='W5GF-AM8',color='brown')
plt.plot((ab8[0:4, 6]), (ab8[0:4, 4]) ,'-o', linewidth=2, label='W5GF-AB8',color='magenta')
#plt.plot(ref[:, 0], ref[:, 3] - ref[:, 1] + pert_ref, 'r', linewidth=2.5, label='reference')
#plt.plot(w3am4[:, 0], w3am4[:, 3] - w3am4[:, 1] + pert4, '-ob', linewidth=2, label='GF-AB4-A',markersize=3)
#plt.plot(w3am4_an[:, 0], w3am4_an[:, 3] - w3am4_an[:, 1] + pert4b, '--k', linewidth=2, label='GF-AB4-B')
#plt.plot(w3_up_nwb[:, 0], w3_up_nwb[:, 3] - w3_up_nwb[:, 1] + pert4b, '--c', linewidth=2, label='WENO3')
#plt.plot(w3am8_an[:, 0], w3am8_an[:, 3] - w3am8_an[:, 1] + pert4b, '--g', linewidth=2, label='GF-AB8-B')
#plt.plot(w3am6[:, 0], w3am6[:, 3] - w3am6[:, 1] + pert6, '--k', linewidth=2, label='GF-AM6')
#plt.plot(w3am8[:, 0], w3am8[:, 3] - w3am8[:, 1] + pert8, '-.c', linewidth=2, label='GF-AM8')

plt.xscale('log')   # Logarithmic X-axis
plt.yscale('log')   # Logarithmic Y-axis

plt.legend(loc='lower left',fontsize=11,framealpha=1.0,ncol=1)
plt.ylabel(r'$ Error $', fontsize=18)
plt.xlabel(r'$ CPU\ time $', fontsize=18)
plt.grid()
#plt.axis([0, 25, -0.0015, 0.0015])

plt.tight_layout()
plt.show()

