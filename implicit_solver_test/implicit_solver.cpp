#include <iostream>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <cmath>
#include <omp.h>
#include "Eigen/Dense" // Include Eigen library

// Function to display the loading bar
void display_loading_bar(int step, int total_steps, double variance) {
    int bar_width = 70;
    float progress = (float)step / total_steps;
    int pos = bar_width * progress;

    std::cout << "[";
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << "% Variance: " << variance << "\r";
    std::cout.flush();
}

int main(int argc, char* argv[]) {
    const double U_s = atof(argv[1]);
    const double phi = atof(argv[2]);
    const double gamma = atof(argv[3]);

    std::cout << "phi :" << phi << std::endl;
    std::cout << "gamma : " << gamma << std::endl;

    const double L = 100.0;  // Length of the domain
    const int N = 200; // Example size, replace with actual value
    const double T = 300.0; // Total time
    const double dt = 0.0005; // Time step

    const double rho_r = 3.333;
    const double rho_d = 0.3;
    const double dx = L / (N - 1);

    const auto k = (phi - 2.) / (phi - 1.) * U_s;
    const auto W = gamma * U_s / (1 - gamma);
    
    int num_steps = static_cast<int>(T / dt);

    const auto start = std::chrono::high_resolution_clock::now();

    std::cout << "Starting memory allocation..." << std::endl;
    std::vector<std::vector<double>> f_u(N, std::vector<double>(N));
    std::vector<std::vector<double>> f_v(N, std::vector<double>(N));
    std::vector<std::vector<double>> u(N, std::vector<double>(N, 1));
    std::vector<std::vector<double>> v(N, std::vector<double>(N, 1));
    std::vector<std::vector<double>> du(N, std::vector<double>(N));
    std::vector<std::vector<double>> dv(N, std::vector<double>(N));

    // Add small random perturbations
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            u[i][j] += ((double) rand() / RAND_MAX) * 1e-5;
            v[i][j] += ((double) rand() / RAND_MAX) * 1e-5;
        }
    }

    std::cout << "Starting simulation..." << std::endl;
    for (int step = 0; step < num_steps; ++step) {
        std::cout << "Step " << step << " of " << num_steps << std::endl;

        // Compute the right-hand side
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

        // Set up the linear system for u and v
        Eigen::MatrixXd A_u = Eigen::MatrixXd::Zero(N * N, N * N);
        Eigen::VectorXd b_u = Eigen::VectorXd::Zero(N * N);
        Eigen::MatrixXd A_v = Eigen::MatrixXd::Zero(N * N, N * N);
        Eigen::VectorXd b_v = Eigen::VectorXd::Zero(N * N);

        // Fill A_u and b_u for the current grid point
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                int idx = i * N + j;
                
                // Diagonal elements
                A_u(idx, idx) = 1.0 + dt * rho_d * (2.0 / (dx * dx) + 2.0 / (dx * dx));
                b_u(idx) = u[i][j] + dt * du[i][j];

                A_v(idx, idx) = 1.0 + dt * rho_d * (2.0 / (dx * dx) + 2.0 / (dx * dx));
                b_v(idx) = v[i][j] + dt * dv[i][j];

                // Off-diagonal elements for diffusion terms with periodic boundary conditions
                int ip1 = (i + 1) % N;
                int im1 = (i - 1 + N) % N;
                int jp1 = (j + 1) % N;
                int jm1 = (j - 1 + N) % N;

                A_u(idx, ip1 * N + j) = -dt * rho_d / (dx * dx);
                A_u(idx, im1 * N + j) = -dt * rho_d / (dx * dx);
                A_u(idx, i * N + jp1) = -dt * rho_d / (dx * dx);
                A_u(idx, i * N + jm1) = -dt * rho_d / (dx * dx);

                A_v(idx, ip1 * N + j) = -dt * rho_d / (dx * dx);
                A_v(idx, im1 * N + j) = -dt * rho_d / (dx * dx);
                A_v(idx, i * N + jp1) = -dt * rho_d / (dx * dx);
                A_v(idx, i * N + jm1) = -dt * rho_d / (dx * dx);
            }
        }

        // Solve the linear system using Eigen
        Eigen::VectorXd u_flat = A_u.colPivHouseholderQr().solve(b_u);
        Eigen::VectorXd v_flat = A_v.colPivHouseholderQr().solve(b_v);

        // Convert the flat solution back to 2D
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                int idx = i * N + j;
                u[i][j] = u_flat(idx);
                v[i][j] = v_flat(idx);
            }
        }

        // Compute variance (example, replace with actual computation)
        double variance = 0.0;
        #pragma omp parallel for reduction(+:variance)
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                variance += (u[i][j] - u[i][j]) * (u[i][j] - u[i][j]);
            }
        }
        variance /= (N * N);

        display_loading_bar(step, num_steps, variance);
    }

    const auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;

    return 0;
}