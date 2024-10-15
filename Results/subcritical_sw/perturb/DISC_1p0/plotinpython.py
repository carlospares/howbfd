import numpy as np
import matplotlib.pyplot as plt

# Load the data from text files
w3_am4_t0p7 = np.loadtxt('w3_am4_t0p7.txt')
w3_am6_t0p7 = np.loadtxt('w3_am6_t0p7.txt')
w3_upwind_t0p7 = np.loadtxt('w3_upwind_t0p7.txt')

# Plot the data
plt.plot(w3_upwind_t0p7[:, 0], w3_upwind_t0p7[:, 3], '--', linewidth=1.5, label='Upwind no wb')
plt.plot(w3_am6_t0p7[:, 0], w3_am6_t0p7[:, 3], '--', linewidth=1.5, label='AM6')
plt.plot(w3_am4_t0p7[:, 0], w3_am4_t0p7[:, 3], '--', linewidth=1.5, label='AM4')

# Set labels and grid
plt.xlabel(r'$x $', fontsize=18)
plt.ylabel(r'$\eta $', fontsize=18)
plt.grid(True)

# Legend settings
plt.legend(fontsize=10)

# Save the figure as a PDF (optional, commented out as in the original code)
# str_filename = 'weno3_sup_DISCbig_t0p7_eta'
# plt.savefig(str_filename + '.pdf', format='pdf', dpi=300)

# Display the plot
plt.show()

