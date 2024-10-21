#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <omp.h>

// Define parameters
const double L = 100.0;  // Length of the domain
const int N = 200;  // Number of grid points
const double T = 300;  // Total time
const double dt = 0.0005;  // Time step size
const double dx = L / (N - 1);

// Constants
const double rho_r = 3.333;
const double rho_d = 0.3;
const double U_s = 0.4;

std::vector<double> k_arr;
std::vector<double> W_arr;

// Function to compute derivatives with periodic boundary conditions
void compute_derivatives(const std::vector<std::vector<double> >& u, const std::vector<std::vector<double> >& v, std::vector<std::vector<double> >& du, std::vector<std::vector<double> >& dv, double k, double W) {
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            int ip1 = (i + 1) % N;
            int im1 = (i - 1 + N) % N;
            int ip2 = (i + 2) % N;
            int im2 = (i - 2 + N) % N;
            int jp1 = (j + 1) % N;
            int jm1 = (j - 1 + N) % N;
            int jp2 = (j + 2) % N;
            int jm2 = (j - 2 + N) % N;

            double u_xx = (-u[ip2][j] + 16 * u[ip1][j] - 30 * u[i][j] + 16 * u[im1][j] - u[im2][j]) / (12 * dx * dx);
            double u_yy = (-u[i][jp2] + 16 * u[i][jp1] - 30 * u[i][j] + 16 * u[i][jm1] - u[i][jm2]) / (12 * dx * dx);
            double u_ij = u[i][j];

            double v_xx = (-v[ip2][j] + 16 * v[ip1][j] - 30 * v[i][j] + 16 * v[im1][j] - v[im2][j]) / (12 * dx * dx);
            double v_yy = (-v[i][jp2] + 16 * v[i][jp1] - 30 * v[i][j] + 16 * v[i][jm1] - v[i][jm2]) / (12 * dx * dx);
            double v_ij = v[i][j];

            double k_minus_U_s_u_ij = k - U_s * u_ij;
            double W_plus_U_s_u_ij = W + U_s * u_ij;

            du[i][j] = rho_r * ((k_minus_U_s_u_ij / (k - U_s)) * u_ij - ((W + U_s)/ W_plus_U_s_u_ij) * u_ij * v_ij) + rho_d * (u_xx + u_yy);
            dv[i][j] = ((W + U_s) / W_plus_U_s_u_ij) * u_ij * v_ij - v_ij * v_ij + v_xx + v_yy;
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

    k_arr.push_back(1.2);
    W_arr.push_back(0.15);
    
    // std::vector<double> W_list = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
    // W_arr.insert(W_arr.end(), W_list.begin(), W_list.end());

    // std::vector<double> k_list = {0.85, 1, 2, 2.5, 3, 3.5, 4};
    // k_arr.insert(k_arr.end(), k_list.begin(), k_list.end());

    for (double k : k_arr){
        std::cout << "phi :" << (k - 2* U_s)/ (k - U_s) << " : ";
        for (double W : W_arr){
            std::cout << "gamma : " << W/(W + U_s)  << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    
    // Initialize concentration of preys and predators
    std::vector<std::vector<double> > u(N, std::vector<double>(N, 1));
    std::vector<std::vector<double> > v(N, std::vector<double>(N, 1));

    // Add small random perturbations
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            u[i][j] += ((double) rand() / RAND_MAX) * 1e-5;
            v[i][j] += ((double) rand() / RAND_MAX) * 1e-5;
        }
    }

    int num_steps = static_cast<int>(T / dt);

    const auto start = std::chrono::high_resolution_clock::now();
    
    // Dormand–Prince coefficients
    const double a2 = 1.0 / 5.0;
    const double a3 = 3.0 / 10.0;
    const double a4 = 4.0 / 5.0;
    const double a5 = 8.0 / 9.0;
    const double b1 = 35.0 / 384.0;
    const double b2 = 0.0;
    const double b3 = 500.0 / 1113.0;
    const double b4 = 125.0 / 192.0;
    const double b5 = -2187.0 / 6784.0;
    const double b6 = 11.0 / 84.0;

    std::vector<std::vector<double>> k1_u(N, std::vector<double>(N));
    std::vector<std::vector<double>> k1_v(N, std::vector<double>(N));
    std::vector<std::vector<double>> k2_u(N, std::vector<double>(N));
    std::vector<std::vector<double>> k2_v(N, std::vector<double>(N));
    std::vector<std::vector<double>> k3_u(N, std::vector<double>(N));
    std::vector<std::vector<double>> k3_v(N, std::vector<double>(N));
    std::vector<std::vector<double>> k4_u(N, std::vector<double>(N));
    std::vector<std::vector<double>> k4_v(N, std::vector<double>(N));
    std::vector<std::vector<double>> k5_u(N, std::vector<double>(N));
    std::vector<std::vector<double>> k5_v(N, std::vector<double>(N));
    std::vector<std::vector<double>> k6_u(N, std::vector<double>(N));
    std::vector<std::vector<double>> k6_v(N, std::vector<double>(N));
    std::vector<std::vector<double>> k7_u(N, std::vector<double>(N));
    std::vector<std::vector<double>> k7_v(N, std::vector<double>(N));
    std::vector<std::vector<double>> u_temp(N, std::vector<double>(N));
    std::vector<std::vector<double>> v_temp(N, std::vector<double>(N));

    for (double k : k_arr) {
        for (double W : W_arr) {
            std::cout << "phi :" << (k - 2* U_s)/ (k - U_s) << std::endl;
            std::cout << "gamma : " << W/(W + U_s)  << std::endl;
            for (int step = 0; step < num_steps; ++step) {
                // Compute k1
                compute_derivatives(u, v, k1_u, k1_v, k, W);

                // Compute k2
                #pragma omp parallel for
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        u_temp[i][j] = u[i][j] + dt * a2 * k1_u[i][j];
                        v_temp[i][j] = v[i][j] + dt * a2 * k1_v[i][j];
                    }
                }
                compute_derivatives(u_temp, v_temp, k2_u, k2_v, k, W);

                // Compute k3
                #pragma omp parallel for
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        u_temp[i][j] = u[i][j] + dt * (3.0 / 40.0 * k1_u[i][j] + 9.0 / 40.0 * k2_u[i][j]);
                        v_temp[i][j] = v[i][j] + dt * (3.0 / 40.0 * k1_v[i][j] + 9.0 / 40.0 * k2_v[i][j]);
                    }
                }
                compute_derivatives(u_temp, v_temp, k3_u, k3_v, k, W);

                // Compute k4
                #pragma omp parallel for
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        u_temp[i][j] = u[i][j] + dt * (44.0 / 45.0 * k1_u[i][j] - 56.0 / 15.0 * k2_u[i][j] + 32.0 / 9.0 * k3_u[i][j]);
                        v_temp[i][j] = v[i][j] + dt * (44.0 / 45.0 * k1_v[i][j] - 56.0 / 15.0 * k2_v[i][j] + 32.0 / 9.0 * k3_v[i][j]);
                    }
                }
                compute_derivatives(u_temp, v_temp, k4_u, k4_v, k, W);

                // Compute k5
                #pragma omp parallel for
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        u_temp[i][j] = u[i][j] + dt * (19372.0 / 6561.0 * k1_u[i][j] - 25360.0 / 2187.0 * k2_u[i][j] + 64448.0 / 6561.0 * k3_u[i][j] - 212.0 / 729.0 * k4_u[i][j]);
                        v_temp[i][j] = v[i][j] + dt * (19372.0 / 6561.0 * k1_v[i][j] - 25360.0 / 2187.0 * k2_v[i][j] + 64448.0 / 6561.0 * k3_v[i][j] - 212.0 / 729.0 * k4_v[i][j]);
                    }
                }
                compute_derivatives(u_temp, v_temp, k5_u, k5_v, k, W);

                // Compute k6
                #pragma omp parallel for
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        u_temp[i][j] = u[i][j] + dt * (9017.0 / 3168.0 * k1_u[i][j] - 355.0 / 33.0 * k2_u[i][j] + 46732.0 / 5247.0 * k3_u[i][j] + 49.0 / 176.0 * k4_u[i][j] - 5103.0 / 18656.0 * k5_u[i][j]);
                        v_temp[i][j] = v[i][j] + dt * (9017.0 / 3168.0 * k1_v[i][j] - 355.0 / 33.0 * k2_v[i][j] + 46732.0 / 5247.0 * k3_v[i][j] + 49.0 / 176.0 * k4_v[i][j] - 5103.0 / 18656.0 * k5_v[i][j]);
                    }
                }
                compute_derivatives(u_temp, v_temp, k6_u, k6_v, k, W);

                // Compute k7
                #pragma omp parallel for
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        u_temp[i][j] = u[i][j] + dt * (35.0 / 384.0 * k1_u[i][j] + 500.0 / 1113.0 * k3_u[i][j] + 125.0 / 192.0 * k4_u[i][j] - 2187.0 / 6784.0 * k5_u[i][j] + 11.0 / 84.0 * k6_u[i][j]);
                        v_temp[i][j] = v[i][j] + dt * (35.0 / 384.0 * k1_v[i][j] + 500.0 / 1113.0 * k3_v[i][j] + 125.0 / 192.0 * k4_v[i][j] - 2187.0 / 6784.0 * k5_v[i][j] + 11.0 / 84.0 * k6_v[i][j]);
                    }
                }
                compute_derivatives(u_temp, v_temp, k7_u, k7_v, k, W);

                // Update u and v
                #pragma omp parallel for
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        u[i][j] += dt * (b1 * k1_u[i][j] + b2 * k2_u[i][j] + b3 * k3_u[i][j] + b4 * k4_u[i][j] + b5 * k5_u[i][j] + b6 * k6_u[i][j]);
                        v[i][j] += dt * (b1 * k1_v[i][j] + b2 * k2_v[i][j] + b3 * k3_v[i][j] + b4 * k4_v[i][j] + b5 * k5_v[i][j] + b6 * k6_v[i][j]);
                    }
                }

                display_loading_bar(step, num_steps);
                if (step == num_steps - 1) 
                {
                    std::ostringstream filename_u;
                    filename_u << "output/frame_u_k" << k << "_W" << W << "_" << std::setw(6) << std::setfill('0') << step << ".csv";
                    save_to_csv(u, filename_u.str());

                    std::ostringstream filename_v;
                    filename_v << "output/frame_v_k" << k << "_W" << W << "_" << std::setw(6) << std::setfill('0') << step << ".csv";
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