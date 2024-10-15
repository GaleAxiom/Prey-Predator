#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <omp.h>

// Define parameters
const double L = 1.0;  // Length of the domain
const int N = 200;  // Number of grid points
const double T = 0.1;  // Total time
const double dt = 0.000001;  // Time step size
const double dx = L / (N - 1);

// Constants
const double rho_r = 0.3;
const double rho_d = 0.3;
const double U_s = 0.4 / 5;

std::vector<double> k_arr;
std::vector<double> W_arr;

// Function to compute derivatives with periodic boundary conditions
void compute_derivatives(const std::vector<std::vector<double> >& u, const std::vector<std::vector<double> >& v, std::vector<std::vector<double> >& du, std::vector<std::vector<double> >& dv, double k, double W) {
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            int ip = (i + 1) % N;
            int im = (i - 1 + N) % N;
            int jp = (j + 1) % N;
            int jm = (j - 1 + N) % N;

            double u_xx = (u[ip][j] + u[im][j] + u[i][jp] + u[i][jm] - 4 * u[i][j]) / (dx * dx);
            double u_ij = u[i][j];

            double v_xx = (v[ip][j] + v[im][j] + v[i][jp] + v[i][jm] - 4 * v[i][j]) / (dx * dx);
            double v_ij = v[i][j];

            double k_minus_U_s_u_ij = k - U_s * u_ij;
            double W_plus_U_s_u_ij = W + U_s * u_ij;

            du[i][j] = rho_r * ((k_minus_U_s_u_ij / (k - U_s)) * u_ij - (W_plus_U_s_u_ij / W_plus_U_s_u_ij) * u_ij * v_ij) + rho_d * u_xx;
            dv[i][j] = (W_plus_U_s_u_ij / W_plus_U_s_u_ij) * u_ij * v_ij - v_ij * v_ij + v_xx;
        }
    }
}

// Function to display the loading bar
void display_loading_bar(int step, int total_steps) {
    int bar_width = 70;
    float progress = (float)step / total_steps;
    int pos = bar_width * progress;

    std::cout << "[";
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

void save_to_csv(const std::vector<std::vector<double> >& data, const std::string& filename) {
    std::ofstream file(filename);
    file << std::fixed << std::setprecision(12);  // Set precision to 6 decimal places
    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i < row.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }
    file.close();
}

int main() {

    printf("Number of threads: %d\n", omp_get_max_threads());
    omp_set_num_threads(omp_get_max_threads());

    k_arr.push_back(0.2);
    W_arr.push_back(0.05);
    
    // Initialize concentration of preys and predators
    std::vector<std::vector<double> > u(N, std::vector<double>(N, 1.0));
    std::vector<std::vector<double> > v(N, std::vector<double>(N, 1.0));

    // Add small random perturbations
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            u[i][j] += ((double) rand() / RAND_MAX) * 0.1;
            v[i][j] += ((double) rand() / RAND_MAX) * 0.1;
        }
    }

    int num_steps = static_cast<int>(T / dt);

    const auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<std::vector<double> > du(N, std::vector<double>(N));
    std::vector<std::vector<double> > dv(N, std::vector<double>(N));

    std::vector<std::vector<double> > u_temp(N, std::vector<double>(N));
    std::vector<std::vector<double> > v_temp(N, std::vector<double>(N));

    std::vector<std::vector<double> > k1_u(N, std::vector<double>(N));
    std::vector<std::vector<double> > k1_v(N, std::vector<double>(N));
    std::vector<std::vector<double> > k2_u(N, std::vector<double>(N));
    std::vector<std::vector<double> > k2_v(N, std::vector<double>(N));
    std::vector<std::vector<double> > k3_u(N, std::vector<double>(N));
    std::vector<std::vector<double> > k3_v(N, std::vector<double>(N));
    std::vector<std::vector<double> > k4_u(N, std::vector<double>(N));
    std::vector<std::vector<double> > k4_v(N, std::vector<double>(N));
    std::vector<std::vector<double> > k5_u(N, std::vector<double>(N));
    std::vector<std::vector<double> > k5_v(N, std::vector<double>(N));
    std::vector<std::vector<double> > k6_u(N, std::vector<double>(N));
    std::vector<std::vector<double> > k6_v(N, std::vector<double>(N));

    for (double k : k_arr) {
        for (double W : W_arr) {
            for (int step = 0; step < num_steps; ++step) {
                // Update concentration fields using RK45

                compute_derivatives(u, v, k1_u, k1_v, k, W);

                #pragma omp parallel for
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        u_temp[i][j] = u[i][j] + dt / 4 * k1_u[i][j];
                        v_temp[i][j] = v[i][j] + dt / 4 * k1_v[i][j];
                    }
                }

                compute_derivatives(u_temp, v_temp, k2_u, k2_v, k, W);

                #pragma omp parallel for
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        u_temp[i][j] = u[i][j] + dt / 8 * k1_u[i][j] + dt / 8 * k2_u[i][j];
                        v_temp[i][j] = v[i][j] + dt / 8 * k1_v[i][j] + dt / 8 * k2_v[i][j];
                    }
                }

                compute_derivatives(u_temp, v_temp, k3_u, k3_v, k, W);

                #pragma omp parallel for
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        u_temp[i][j] = u[i][j] + dt / 2 * k3_u[i][j];
                        v_temp[i][j] = v[i][j] + dt / 2 * k3_v[i][j];
                    }
                }

                compute_derivatives(u_temp, v_temp, k4_u, k4_v, k, W);

                #pragma omp parallel for
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        u_temp[i][j] = u[i][j] + dt * (3 * k1_u[i][j] + 9 * k4_u[i][j]) / 16;
                        v_temp[i][j] = v[i][j] + dt * (3 * k1_v[i][j] + 9 * k4_v[i][j]) / 16;
                    }
                }

                compute_derivatives(u_temp, v_temp, k5_u, k5_v, k, W);

                #pragma omp parallel for
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        u_temp[i][j] = u[i][j] + dt * (-3 * k1_u[i][j] + 2 * k2_u[i][j] + 12 * k3_u[i][j] - 12 * k4_u[i][j] + 8 * k5_u[i][j]) / 7;
                        v_temp[i][j] = v[i][j] + dt * (-3 * k1_v[i][j] + 2 * k2_v[i][j] + 12 * k3_v[i][j] - 12 * k4_v[i][j] + 8 * k5_v[i][j]) / 7;
                    }
                }

                compute_derivatives(u_temp, v_temp, k6_u, k6_v, k, W);

                // Update u and v
                #pragma omp parallel for
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        u[i][j] += dt * (7 * k1_u[i][j] + 32 * k3_u[i][j] + 12 * k4_u[i][j] + 32 * k5_u[i][j] + 7 * k6_u[i][j]) / 90;
                        v[i][j] += dt * (7 * k1_v[i][j] + 32 * k3_v[i][j] + 12 * k4_v[i][j] + 32 * k5_v[i][j] + 7 * k6_v[i][j]) / 90;
                    }
                }
                
                if (step % 1000 == 0) {
                    //std::cout << "Step: " << step << " out of " << num_steps << std::endl;
                    // Display loading bar
                    display_loading_bar(step, num_steps);
                    std::ostringstream filename_u;
                    filename_u << "output/frame_u_" << std::setw(6) << std::setfill('0') << step << ".csv";
                    save_to_csv(u, filename_u.str());

                    std::ostringstream filename_v;
                    filename_v << "output/frame_v_" << std::setw(6) << std::setfill('0') << step << ".csv";
                    save_to_csv(v, filename_v.str());
                }
            }
        }
    }
    
    const auto end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s" << std::endl;
    return 0;
}