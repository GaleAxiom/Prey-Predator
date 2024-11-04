import os
import numpy as np
import imageio.v2 as imageio
import matplotlib.pyplot as plt


# load data from CSV files frames_u_*.csv
frames = []
for file in os.listdir('output/'):
    if file.startswith('frame_u_'):
        frames.append(np.genfromtxt(f'output/{file}', delimiter=','))
frames = np.array(frames)

# Convert CSV files to SVG images
for i, frame in enumerate(frames):
        plt.imshow(frame.T, cmap='seismic', origin='lower')
        plt.colorbar()
        plt.savefig(f'output/frame_{i}.png')
        plt.clf()  # Clear the figure for the next plot

# Convert SVG images to GIF
# Create a list to store the file paths of the SVG images
svg_files = [f'output/frame_{i}.png' for i in range(len(frames))]

# Create a list to store the file paths of the GIF frames
gif_frames = []

# Convert SVG images to GIF frames
for svg_file in svg_files:
    gif_frame = imageio.imread(svg_file)
    gif_frames.append(gif_frame)

# Save the GIF file
output_file = 'output/animation.gif'
imageio.mimsave(output_file, gif_frames, duration=0.5)

# # Remove the SVG images
# for svg_file in svg_files:
#     os.remove(svg_file)

# # Remove the CSV files
# for file in os.listdir('output/'):
#     if file.startswith('frame_u_'):
#         os.remove(f'output/{file}')