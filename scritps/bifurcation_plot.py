import matplotlib.pyplot as plt

# Function to read data from a file and extract phi and gamma values
def read_data(file_path):
    data = open(file_path).read()
    lines = data.strip().split('\n')
    phi_values = []
    gamma_values = []
    for i in range(0, len(lines), 2):
        phi = float(lines[i].split(':')[1].strip())
        gamma = float(lines[i + 1].split(':')[1].strip())
        phi_values.append(phi)
        gamma_values.append(gamma)
    return phi_values, gamma_values

# Read data for each category
phi_values_stripes, gamma_values_stripes = read_data('output_stripes/value')
phi_values_hot, gamma_values_hot = read_data('output_hot_spots/value')
phi_values_cold, gamma_values_cold = read_data('output_cold_spots/value')
phi_values_homo_large, gamma_values_homo_large = read_data('output_homo_>1/value')
phi_values_chaos, gamma_values_chaos = read_data('output_chaos/value')

# Create scatter plot
plt.scatter(phi_values_stripes, gamma_values_stripes, color='blue', label='Stripes')
plt.scatter(phi_values_hot, gamma_values_hot, color='red', label='Hot')
plt.scatter(phi_values_cold, gamma_values_cold, color='green', label='Cold')
plt.scatter(phi_values_homo_large, gamma_values_homo_large, color='yellow', label='Homo > 1')
plt.scatter(phi_values_chaos, gamma_values_chaos, color='purple', label='Chaos')

# Add labels and title
plt.xlabel('Phi')
plt.ylabel('Gamma')
plt.title('Scatter Plot of Phi vs Gamma')
plt.legend()

# Show plot
plt.show()