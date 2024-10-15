import numpy as np
import matplotlib.pyplot as plt

# Load the data from text files
w3_am4_n100 = np.loadtxt('w3_am4_n100.txt')
w3_am8_n100 = np.loadtxt('w3_am8_n100.txt')
w3_upwind_n100 = np.loadtxt('w3_upwind_n100.txt')
w3_am4_n100_st = np.loadtxt('w3_am4_n100_st.txt')
w3_am8_n100_st = np.loadtxt('w3_am8_n100_st.txt')
w3_upwind_n100_st = np.loadtxt('w3_upwind_n100_st.txt')

# Create the first subplot
plt.subplot(2, 1, 1)
plt.plot(w3_upwind_n100[:, 0], w3_upwind_n100[:, 3] - w3_upwind_n100_st[:, 3], linewidth=1.5, label='WENO3')
plt.plot(w3_am8_n100[:, 0], w3_am8_n100[:, 3] - w3_am8_n100_st[:, 3], linewidth=1.5, label='GF-AM8')
plt.plot(w3_am4_n100[:, 0], w3_am4_n100[:, 3] - w3_am4_n100_st[:, 3], linewidth=1.5, label='GF-AM6')

# Set labels and grid for the first subplot
plt.xlabel(r'$x $', fontsize=18)
plt.ylabel(r'$h-h_{in} $', fontsize=18)
plt.grid(True)
plt.legend(fontsize=10)

# Create the second subplot
plt.subplot(2, 1, 2)
plt.plot(w3_am8_n100[:, 0], w3_am8_n100[:, 3] - w3_am8_n100_st[:, 3], linewidth=1.5, label='GF-AM8')
plt.plot(w3_am4_n100[:, 0], w3_am4_n100[:, 3] - w3_am4_n100_st[:, 3], linewidth=1.5, label='GF-AM6')

# Set labels and grid for the second subplot
plt.xlabel(r'$x $', fontsize=18)
plt.ylabel(r'$h-h_{in} $', fontsize=18)
plt.grid(True)
plt.legend(fontsize=10)

# Save the figure (this line is equivalent to the commented MATLAB savefig code)
str_filename = 'weno3_sub_DISCsmall_n100_eta'
plt.savefig(str_filename + '.pdf', format='pdf', dpi=300)

# Show the plot
plt.show()

