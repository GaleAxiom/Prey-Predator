import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np
import glob
import os

# Create a list to store frames
frames = []

# Find all frame_u_*.csv files in the output directory
csv_files_u = sorted(glob.glob('output/frame_u_*.csv'))
csv_files_v = sorted(glob.glob('output/frame_v_*.csv'))

for csv_file in csv_files_u:
    # Load data from CSV
    u = np.loadtxt(csv_file, delimiter=',')
    
    # Plot u
    plt.imshow(u, cmap='hot', origin='lower', extent=[0, 1, 0, 1])
    plt.colorbar(label='Concentration')
    # Save the plot as an image
    image_path = csv_file.replace('.csv', '.png')
    plt.savefig(image_path)
    plt.clf()
    frames.append(imageio.imread(image_path))

# Save frames as a GIF
imageio.mimsave('GIFs/output_animation_u.gif', frames, duration=0.5)

frames = []

for csv_file in csv_files_v:
    # Load data from CSV
    u = np.loadtxt(csv_file, delimiter=',')
    
    # Plot u
    plt.imshow(u, cmap='hot', origin='lower', extent=[0, 1, 0, 1])
    plt.colorbar(label='Concentration')
    # Save the plot as an image
    image_path = csv_file.replace('.csv', '.png')
    plt.savefig(image_path)
    plt.clf()
    frames.append(imageio.imread(image_path))

# Save frames as a GIF
imageio.mimsave('GIFs/output_animation_v.gif', frames, duration=0.5)

# Remove the CSV and PNG files
for csv_file in csv_files_u:
    image_path = csv_file.replace('.csv', '.png')
    os.remove(csv_file)
    os.remove(image_path)

for csv_file in csv_files_v:
    image_path = csv_file.replace('.csv', '.png')
    os.remove(csv_file)
    os.remove(image_path)