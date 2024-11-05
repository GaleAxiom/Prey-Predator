import numpy as np
import matplotlib.pyplot as plt

# Load the data from solution.csv
data = np.loadtxt('output/frame_u_phi0.75_gamma0.4_022000.csv', delimiter=',')

# Plot the data using imshow
plt.imshow(data)

# Add a colorbar for reference
plt.colorbar()

# Show the plot
plt.show()