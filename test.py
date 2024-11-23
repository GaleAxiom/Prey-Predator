import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# # Load the data from solution.csv
# data_u = np.loadtxt('output/spatial_average_u_phi0.840000_gamma0.480000.csv', delimiter=',')
# data_v = np.loadtxt('output/spatial_average_v_phi0.840000_gamma0.480000.csv', delimiter=',')

# t, v = data_v[:,0], data_v[:,1]
# plt.plot(t, v)
# plt.plot(t, data_u[:,1], label='u')
# plt.xlabel('Time')
# plt.ylabel('Concentration')
# plt.legend()
# plt.show()


plt.rcParams.update({
    "axes.grid": True,
    "grid.color" : (0.5, 0.5, 0.5, 0.4),
    "figure.figsize": (10, 5),
    "figure.constrained_layout.use": True,
    "font.size": 28,
    "axes.labelsize": 28,
    "font.family": "Times New Roman",
    "mathtext.fontset": "cm",
    "axes.prop_cycle": plt.cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#bcbd22', '#17becf', '#7f7f7f', 'darkblue', 'gold', 'salmon', 'teal', '#e377c2', 'tan', 'darkred']),
    "lines.linewidth": 2
})

def plot_spatial_data(file_u, file_v, title):
    # Load the data from the specified CSV files
    data_u = np.loadtxt(file_u, delimiter=',')
    data_v = np.loadtxt(file_v, delimiter=',')
    t, u = data_u[:,0], data_u[:,1]
    t, v = data_v[:,0], data_v[:,1]

    # Normalize the time data for colormap
    norm = plt.Normalize(t.min(), t.max())
    cmap = cm.viridis

    # Plot the data using line segments with color mapping
    fig, ax = plt.subplots()
    for i in range(len(t) - 1):
        ax.plot(u[i:i+2], v[i:i+2], color=cmap(norm(t[i])))

    ax.set_xlabel(r'Prey concentration: $u$[-]')
    ax.set_ylabel(r'Predator concentration: $v$[-]')

    # Add a colorbar for reference
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label='time')
    # Save the plot
    plt.savefig('figures/' + str(title) + '.svg')
    # plt.show()
    

plot_spatial_data('output/spatial_average_u_chaos.csv', 'output/spatial_average_v_chaos.csv', 'chaos')
plot_spatial_data('output/spatial_average_u_cold_spots.csv', 'output/spatial_average_v_cold_spots.csv', 'cold_spots')
plot_spatial_data('output/spatial_average_u_hot_spots.csv', 'output/spatial_average_v_hot_spots.csv', 'hot spots')
plot_spatial_data('output/spatial_average_u_stripes.csv', 'output/spatial_average_v_stripes.csv', 'stripes')


#! Temporal plot
def plot_temporal_data(file_u, file_v, title):
    # Load the data from the specified CSV files
    data_u = np.loadtxt(file_u, delimiter=',')
    data_v = np.loadtxt(file_v, delimiter=',')
    t, u = data_u[:,0], data_u[:,1]
    t, v = data_v[:,0], data_v[:,1]

    # Plot the data
    fig, ax = plt.subplots()
    ax.plot(t, u, label='u')
    ax.plot(t, v, label='v')
    ax.set_xlabel(r'Time')
    ax.set_ylabel(r'Concentration')
    ax.legend()

    # Show the plot
    plt.savefig('figures/' + str(title) + '_temporal.svg')

#! use of average data
plot_temporal_data('output/spatial_average_u_chaos.csv', 'output/spatial_average_v_chaos.csv', 'chaos_avg')
plot_temporal_data('output/spatial_average_u_cold_spots.csv', 'output/spatial_average_v_cold_spots.csv', 'cold_spots_avg')
plot_temporal_data('output/spatial_average_u_hot_spots.csv', 'output/spatial_average_v_hot_spots.csv', 'hot spots_avg')
plot_temporal_data('output/spatial_average_u_stripes.csv', 'output/spatial_average_v_stripes.csv', 'stripes_avg')

plot_temporal_data('output/spatial_average_u_chaos.csv', 'output/spatial_average_v_chaos.csv', 'chaos_point')
plot_temporal_data('output/spatial_average_u_cold_spots.csv', 'output/spatial_average_v_cold_spots.csv', 'cold_spots_point')
plot_temporal_data('output/spatial_average_u_hot_spots.csv', 'output/spatial_average_v_hot_spots.csv', 'hot spots_point')
plot_temporal_data('output/spatial_average_u_stripes.csv', 'output/spatial_average_v_stripes.csv', 'stripes_point')



def plot_matrix_to_svg(matrix, title):
    plt.imshow(matrix)
    plt.colorbar()
    plt.set_cmap('seismic')

    plt.xlabel(r'X-axis')
    plt.ylabel(r'Y-axis')

    plt.savefig('figures/' + str(title) + '.svg')
    plt.clf()

# Example usage
u = np.loadtxt('output/frame_u_stripes.csv', delimiter=',')
v = np.loadtxt('output/frame_v_stripes.csv', delimiter=',')

plot_matrix_to_svg(u, 'u_stripes_frame')
plot_matrix_to_svg(v, 'v_stripes_frame')

u = np.loadtxt('output/frame_u_chaos.csv', delimiter=',')
v = np.loadtxt('output/frame_v_chaos.csv', delimiter=',')

plot_matrix_to_svg(u, 'u_chaos')
plot_matrix_to_svg(v, 'v_chaos')

u = np.loadtxt('output/frame_u_cold_spots.csv', delimiter=',')
v = np.loadtxt('output/frame_v_cold_spots.csv', delimiter=',')

plot_matrix_to_svg(u, 'u_cold_spots_frame')
plot_matrix_to_svg(v, 'v_cold_spots_frame')

u = np.loadtxt('output/frame_u_hot_spots.csv', delimiter=',')
v = np.loadtxt('output/frame_v_hot_spots.csv', delimiter=',')

plot_matrix_to_svg(u, 'u_hot_spots_frame')
plot_matrix_to_svg(v, 'v_hot_spots_frame')