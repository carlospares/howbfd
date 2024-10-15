import numpy as np
import matplotlib.pyplot as plt

# Load the data from text files
w3_am3_n150 = np.loadtxt('w3_am3_n150.txt')
w3_am4_n150 = np.loadtxt('w3_am4_n150.txt')
w3_am6_n150 = np.loadtxt('w3_am6_n150.txt')
w3_upwind_n150 = np.loadtxt('w3_upwind_n150.txt')

w3_am3_n150_stationary = np.loadtxt('w3_am3_n150_stationary.txt')
w3_am4_n150_stationary = np.loadtxt('w3_am4_n150_stationary.txt')
w3_am6_n150_stationary = np.loadtxt('w3_am6_n150_stationary.txt')
w3_upwind_n150_stationary = np.loadtxt('w3_upwind_n150_stationary.txt')

# Plot the WENO3 data
plt.plot(w3_am3_n150[:, 0], w3_am3_n150[:, 2] - w3_am3_n150_stationary[:, 2], linewidth=1.5, label='AM3')
plt.plot(w3_am4_n150[:, 0], w3_am4_n150[:, 2] - w3_am4_n150_stationary[:, 2], linewidth=1.5, label='AM4')
plt.plot(w3_am6_n150[:, 0], w3_am6_n150[:, 2] - w3_am6_n150_stationary[:, 2], linewidth=1.5, label='AM6')
plt.plot(w3_upwind_n150[:, 0], w3_upwind_n150[:, 2] - w3_upwind_n150_stationary[:, 2], linewidth=1.5, label='upwind no wb')

# Add grid, labels, and legend for WENO3
plt.grid(True)
plt.xlabel(r'$x $', fontsize=18)
plt.ylabel(r'$u-u_{in} $', fontsize=18)
#plt.legend(loc='northwest', fontsize=10)
plt.legend()
plt.show()

# (Optional) Save the figure (uncomment if needed)
# plt.savefig('weno3_AM_error_pert0p005.pdf', format='pdf', dpi=300)

# Plot the WENO5 data
w5_am3_n150 = np.loadtxt('w5_am3_n150.txt')
w5_am4_n150 = np.loadtxt('w5_am4_n150.txt')
w5_am6_n150 = np.loadtxt('w5_am6_n150.txt')
w5_upwind_n150 = np.loadtxt('w5_upwind_n150.txt')
reference_p = np.loadtxt('reference_p_w3.txt')  # Assuming you have this data as well

w5_am3_n150_stationary = np.loadtxt('w5_am3_n150_stationary.txt')
w5_am4_n150_stationary = np.loadtxt('w5_am4_n150_stationary.txt')
w5_am6_n150_stationary = np.loadtxt('w5_am6_n150_stationary.txt')
w5_upwind_n150_stationary = np.loadtxt('w5_upwind_n150_stationary.txt')

# Plot WENO5 comparison
plt.plot(w5_am3_n150[:, 0], w5_am3_n150[:, 2] - w5_am3_n150_stationary[:, 2], linewidth=1.5, label='WENO5 AM3')
plt.plot(w5_am4_n150[:, 0], w5_am4_n150[:, 2] - w5_am4_n150_stationary[:, 2], linewidth=1.5, label='WENO5 AM4')
plt.plot(w5_am6_n150[:, 0], w5_am6_n150[:, 2] - w5_am6_n150_stationary[:, 2], linewidth=1.5, label='WENO5 AM6')
plt.plot(w5_upwind_n150[:, 0], w5_upwind_n150[:, 2] - w5_upwind_n150_stationary[:, 2], linewidth=1.5, label='WENO5 Upwind')

# Plot reference
out = reference_p[:, 2] - reference_p[:, 1]
out[500:1500] = 0  # Assuming this applies only to indices 500:1500
plt.plot(reference_p[:, 0], out, 'k', linewidth=1.5, label='Reference')

# Add grid, labels, and legend for WENO5
plt.grid(True)
plt.xlabel(r'$x $', fontsize=18)
plt.ylabel(r'$u-u_{in} $', fontsize=18)
plt.legend(loc='northwest', fontsize=10)
plt.show()

# (Optional) Save the figure (uncomment if needed)
# plt.savefig('weno5_AM_error_pert0p005.pdf', format='pdf', dpi=300)

