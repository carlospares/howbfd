import numpy as np
import matplotlib.pyplot as plt

# Load the data from text files
w5_am3_n250 = np.loadtxt('w5_am3_n250.txt')
w5_am4_n250 = np.loadtxt('w5_am4_n250.txt')
w5_am6_n250 = np.loadtxt('w5_am6_n250.txt')
w5_upwind_n250 = np.loadtxt('w5_upwind_n250.txt')

w5_am3_n250_stationary = np.loadtxt('w5_am3_n250_stationary.txt')
w5_am4_n250_stationary = np.loadtxt('w5_am4_n250_stationary.txt')
w5_am6_n250_stationary = np.loadtxt('w5_am6_n250_stationary.txt')
w5_upwind_n250_stationary = np.loadtxt('w5_upwind_n250_stationary.txt')

# Plot the data for WENO5
plt.plot(w5_am3_n250[:, 0], w5_am3_n250[:, 2] - w5_am3_n250_stationary[:, 2], linewidth=1.5, label='WENO5 GF-AM3')
plt.plot(w5_am4_n250[:, 0], w5_am4_n250[:, 2] - w5_am4_n250_stationary[:, 2], linewidth=1.5, label='WENO5 GF-AM4')
plt.plot(w5_am6_n250[:, 0], w5_am6_n250[:, 2] - w5_am6_n250_stationary[:, 2], linewidth=1.5, label='WENO5 GF-AM6')
plt.plot(w5_upwind_n250[:, 0], w5_upwind_n250[:, 2] - w5_upwind_n250_stationary[:, 2], linewidth=1.5, label='WENO5')

# Set grid, labels, and axis limits
plt.grid(True)
plt.xlabel(r'$x $', fontsize=18)
plt.ylabel(r'$u-u_{in} $', fontsize=18)

# Set legend
plt.legend(loc='northwest', fontsize=10)

# Set axis limits
plt.axis([-1, 1, 0, 20 * 10**-3])

# Show the plot
plt.show()

# (Optional) Save the figure (uncomment if needed)
# str_filename = 'weno5_AM_error_pert0p005'
# plt.savefig(str_filename + '.pdf', format='pdf', dpi=300)

