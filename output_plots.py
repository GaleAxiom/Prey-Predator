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

# for csv_file in csv_files_u:
#     # Load data from CSV
#     u = np.loadtxt(csv_file, delimiter=',')
    
#     # Plot u
#     plt.imshow(u, cmap='seismic', origin='lower', extent=[0, 1, 0, 1])
#     plt.colorbar(label='Concentration')
#     # Save the plot as an image
#     image_path = csv_file.replace('.csv', '.png')
#     plt.savefig(image_path)
#     plt.clf()
#     frames.append(imageio.imread(image_path))

# # Save frames as a GIF
# imageio.mimsave('GIFs/output_animation_u.gif', frames, duration=1)

# frames = []

# for csv_file in csv_files_v:
#     # Load data from CSV
#     u = np.loadtxt(csv_file, delimiter=',')
    
#     # Plot u
#     plt.imshow(u, cmap='seismic', origin='lower', extent=[0, 1, 0, 1])
#     plt.colorbar(label='Concentration')
#     # Save the plot as an image
#     image_path = csv_file.replace('.csv', '.png')
#     plt.savefig(image_path)
#     plt.clf()
#     frames.append(imageio.imread(image_path))

# # Save frames as a GIF
# imageio.mimsave('GIFs/output_animation_v.gif', frames, duration=1)

for i in range(0, len(csv_files_u)):
    u = np.loadtxt(csv_files_u[i], delimiter=',')
    plt.imshow(u, cmap='seismic', origin='lower', extent=[0, 1, 0, 1])
    plt.colorbar(label='Concentration')
    image_path = csv_files_u[i].rstrip('.csv')
    plt.savefig(image_path + '.svg')
    plt.clf()

    # v = np.loadtxt(csv_files_v[i], delimiter=',')
    # plt.imshow(v, cmap='seismic', origin='lower', extent=[0, 1, 0, 1])
    # plt.colorbar(label='Concentration')
    # image_path = csv_files_v[i].rstrip('.csv')
    # plt.savefig(image_path + '.svg')
    # plt.clf()


# Remove the CSV files in the output directory
for csv_file in csv_files_u:
    os.remove(csv_file)