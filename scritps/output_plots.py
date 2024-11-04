import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os

# Load data from CSV
data = np.genfromtxt('output/message.txt', delimiter=',')
data = data[1:]  # Skip header if present

charaterizatiion = ['Homo = 1', 'Homo > 1', 'Cold Spots', 'Stripes', 'Hot Spots', 'Plates', 'Chaos']


int_to_char = {i: char for i, char in enumerate(charaterizatiion)}

# Sort data according to the rule data[:,0] * 1000000 + data[:,1]
data = data[np.argsort(data[:, 0] * 1000000 + data[:, 1])]
data[:,2] = data[:,2] - 1

# Extract coordinates and values
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

# Determine the grid shape
num_x = len(np.unique(x))
num_y = len(np.unique(y))

# Reshape z to be a 2D array with the correct dimensions
Z = z.reshape((num_x, num_y))

# Create a meshgrid for x and y
X, Y = np.meshgrid(np.unique(x), np.unique(y))

boundaries = np.arange(len(charaterizatiion) + 1) - 0.5
ticks = np.arange(len(charaterizatiion))

#coolwarm, nipy_spectral
cmap = plt.get_cmap('nipy_spectral', len(charaterizatiion))
sc = plt.pcolormesh(X, Y, Z.T, cmap=cmap, shading='auto')
norm = mcolors.BoundaryNorm(boundaries, cmap.N)
cbar = plt.colorbar(sc, boundaries=boundaries, ticks=ticks, norm=norm)
cbar.ax.set_yticklabels((charaterizatiion))
plt.xlabel(r'$\phi$')
plt.ylabel(r'$\gamma$')
plt.xticks(np.arange(0.5, 1.01, 0.1))
plt.yticks(np.arange(0.2, 0.75, 0.05))
plt.savefig('output/classification.svg')
plt.clf()  # Clear the figure for the next plot

plt.scatter(x, y, c=z, cmap=cmap, norm=norm)
plt.colorbar()
plt.savefig('output/classification_scatter.svg')


# Load another data file
data = np.genfromtxt('/Users/florentdistree/Documents/Uni/Irreversibility/output/frame_u_phi0.9_gamma0.54_600000.csv', delimiter=',')
plt.imshow(data, cmap='seismic', origin='lower')
plt.colorbar()
plt.savefig('output/plates.svg')