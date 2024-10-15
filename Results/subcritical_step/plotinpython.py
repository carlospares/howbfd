import numpy as np
import matplotlib.pyplot as plt

# Load the data from the text file
weno3_am4_per_200 = np.loadtxt('weno3_am4_per_200.txt')

# Plot the data
plt.plot(weno3_am4_per_200[:, 0], weno3_am4_per_200[:, 3] - weno3_am4_per_200[:, 1], linewidth=2)

# Set labels
plt.xlabel('x')
plt.ylabel('Difference (h-h_in)')

# Display the plot
plt.grid(True)
plt.show()
