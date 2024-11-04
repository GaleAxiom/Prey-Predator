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

phi = np.linspace(np.min(data[:, 0]), np.max(data[:, 0]), num_x)
gamma_transcritical = phi/2
plt.plot(phi, gamma_transcritical, label=r'Transcritical', color='red', linewidth=1.4)

rho_r = 3.333
rho_d = 0.3
gamma_hopf = phi - 1/rho_r
plt.plot(phi[gamma_hopf > gamma_transcritical], gamma_hopf[gamma_hopf > gamma_transcritical], label=r'Hopf', color='cyan', linewidth=1.4)

gamma_turing = rho_d/rho_r * (1 - np.sqrt(2 + rho_r/rho_d * phi))**2
print(phi[(gamma_turing > gamma_hopf) & (gamma_turing < phi - rho_d/rho_r)])
plt.plot(phi[(gamma_turing > gamma_hopf) & (gamma_turing < phi - rho_d/rho_r)], \
        gamma_turing[(gamma_turing > gamma_hopf) & (gamma_turing < phi - rho_d/rho_r)], label=r'Turing', color='yellow', linewidth=1.4)

plt.legend()
plt.xlabel(r'$\phi$')
plt.ylabel(r'$\gamma$')
plt.xticks(np.arange(0.5, 1.01, 0.1))
plt.yticks(np.arange(0.2, 0.75, 0.05))
plt.tight_layout()
plt.savefig('output/classification.svg')
plt.clf()  # Clear the figure for the next plot
