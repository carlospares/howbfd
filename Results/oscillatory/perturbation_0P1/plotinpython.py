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

# Plot the data
plt.plot(w3_am3_n150[:, 0], w3_am3_n150[:, 2] - w3_am3_n150_stationary[:, 2], linewidth=1.5, label='AM3')
plt.plot(w3_am4_n150[:, 0], w3_am4_n150[:, 2] - w3_am4_n150_stationary[:, 2], linewidth=1.5, label='AM4')
plt.plot(w3_am6_n150[:, 0], w3_am6_n150[:, 2] - w3_am6_n150_stationary[:, 2], linewidth=1.5, label='AM6')
plt.plot(w3_upwind_n150[:, 0], w3_upwind_n150[:, 2] - w3_upwind_n150_stationary[:, 2], linewidth=1.5, label='upwind no wb')

# Add grid, labels, and legend
plt.grid(True)
plt.xlabel(r'$x $', fontsize=18)
plt.ylabel(r'$u-u_{in} $', fontsize=18)

# Set legend and its font size
plt.legend(loc='northwest', fontsize=10)

# Optionally, set axis limits (uncomment if needed)
# plt.axis([-1, 1, -1, 3])

# Show the plot
plt.show()

# (Optional) Save the figure (uncomment if needed)
# str_filename = 'plot_output'
# plt.savefig(str_filename + '.pdf', format='pdf', dpi=300)

