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
const double T = 500;  // Total time
const double dt = 0.0005;  // Time step size
const double dx = L / (N - 1);

// Constants
const double rho_r = 3.333;
const double rho_d = 0.3;

// Function to compute derivatives with periodic boundary conditions
void compute_derivatives(const std::vector<std::vector<double> >& u, const std::vector<std::vector<double> >& v, std::vector<std::vector<double> >& du, std::vector<std::vector<double> >& dv, double k, double W, double U_s) {
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

int main(int argc, char* argv[]) {

    printf("Number of threads: %d\n", omp_get_max_threads());
    omp_set_num_threads(omp_get_max_threads());

    const double U_s    = atof(argv[1]);
    const double phi    = atof(argv[2]);
    const double gamma  = atof(argv[3]);
    
    std::vector<std::vector<double> > u_init(N, std::vector<double>(N, 1));
    std::vector<std::vector<double> > v_init(N, std::vector<double>(N, 1));

    // Add small random perturbations
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            u_init[i][j] += ((double) rand() / RAND_MAX) * 1e-5;
            v_init[i][j] += ((double) rand() / RAND_MAX) * 1e-5;
        }
    }

    int num_steps = static_cast<int>(T / dt);

    const auto start = std::chrono::high_resolution_clock::now();

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

    std::vector<std::vector<double>> u(N, std::vector<double>(N));
    std::vector<std::vector<double>> v(N, std::vector<double>(N));


    const auto k = (phi - 2.) / (phi - 1.) * U_s;
    const auto W =  gamma * U_s / (1 - gamma);
    
    std::cout << "phi :" << (k - 2* U_s)/ (k - U_s) << std::endl;
    std::cout << "gamma : " << W/(W + U_s)  << std::endl;

    std::fill(u.begin(), u.end(), std::vector<double>(N, 0));
    std::fill(v.begin(), v.end(), std::vector<double>(N, 0));
    std::fill(k1_u.begin(), k1_u.end(), std::vector<double>(N, 0));
    std::fill(k1_v.begin(), k1_v.end(), std::vector<double>(N, 0));
    std::fill(k2_u.begin(), k2_u.end(), std::vector<double>(N, 0));
    std::fill(k2_v.begin(), k2_v.end(), std::vector<double>(N, 0));
    std::fill(k3_u.begin(), k3_u.end(), std::vector<double>(N, 0));
    std::fill(k3_v.begin(), k3_v.end(), std::vector<double>(N, 0));
    std::fill(k4_u.begin(), k4_u.end(), std::vector<double>(N, 0));
    std::fill(k4_v.begin(), k4_v.end(), std::vector<double>(N, 0));
    std::fill(k5_u.begin(), k5_u.end(), std::vector<double>(N, 0));
    std::fill(k5_v.begin(), k5_v.end(), std::vector<double>(N, 0));
    std::fill(k6_u.begin(), k6_u.end(), std::vector<double>(N, 0));
    std::fill(k6_v.begin(), k6_v.end(), std::vector<double>(N, 0));
    std::fill(k7_u.begin(), k7_u.end(), std::vector<double>(N, 0));
    std::fill(k7_v.begin(), k7_v.end(), std::vector<double>(N, 0));
    std::fill(u_temp.begin(), u_temp.end(), std::vector<double>(N, 0));
    std::fill(v_temp.begin(), v_temp.end(), std::vector<double>(N, 0));


    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            u[i][j] = u_init[i][j];
            v[i][j] = v_init[i][j];
        }
    }

    for (int step = 0; step < num_steps; ++step) {
        // Compute k1
        compute_derivatives(u, v, k1_u, k1_v, k, W, U_s);

        // Compute k2
        #pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                u_temp[i][j] = u[i][j] + dt * 1.0 / 5.0 * k1_u[i][j];
                v_temp[i][j] = v[i][j] + dt * 1.0 / 5.0 * k1_v[i][j];
            }
        }
        compute_derivatives(u_temp, v_temp, k2_u, k2_v, k, W, U_s);

        // Compute k3
        #pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                u_temp[i][j] = u[i][j] + dt * (3.0 / 40.0 * k1_u[i][j] + 9.0 / 40.0 * k2_u[i][j]);
                v_temp[i][j] = v[i][j] + dt * (3.0 / 40.0 * k1_v[i][j] + 9.0 / 40.0 * k2_v[i][j]);
            }
        }
        compute_derivatives(u_temp, v_temp, k3_u, k3_v, k, W, U_s);

        // Compute k4
        #pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                u_temp[i][j] = u[i][j] + dt * (44.0 / 45.0 * k1_u[i][j] - 56.0 / 15.0 * k2_u[i][j] + 32.0 / 9.0 * k3_u[i][j]);
                v_temp[i][j] = v[i][j] + dt * (44.0 / 45.0 * k1_v[i][j] - 56.0 / 15.0 * k2_v[i][j] + 32.0 / 9.0 * k3_v[i][j]);
            }
        }
        compute_derivatives(u_temp, v_temp, k4_u, k4_v, k, W, U_s);

        // Compute k5
        #pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                u_temp[i][j] = u[i][j] + dt * (19372.0 / 6561.0 * k1_u[i][j] - 25360.0 / 2187.0 * k2_u[i][j] + 64448.0 / 6561.0 * k3_u[i][j] - 212.0 / 729.0 * k4_u[i][j]);
                v_temp[i][j] = v[i][j] + dt * (19372.0 / 6561.0 * k1_v[i][j] - 25360.0 / 2187.0 * k2_v[i][j] + 64448.0 / 6561.0 * k3_v[i][j] - 212.0 / 729.0 * k4_v[i][j]);
            }
        }
        compute_derivatives(u_temp, v_temp, k5_u, k5_v, k, W, U_s);

        // Compute k6
        #pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                u_temp[i][j] = u[i][j] + dt * (9017.0 / 3168.0 * k1_u[i][j] - 355.0 / 33.0 * k2_u[i][j] + 46732.0 / 5247.0 * k3_u[i][j] + 49.0 / 176.0 * k4_u[i][j] - 5103.0 / 18656.0 * k5_u[i][j]);
                v_temp[i][j] = v[i][j] + dt * (9017.0 / 3168.0 * k1_v[i][j] - 355.0 / 33.0 * k2_v[i][j] + 46732.0 / 5247.0 * k3_v[i][j] + 49.0 / 176.0 * k4_v[i][j] - 5103.0 / 18656.0 * k5_v[i][j]);
            }
        }
        compute_derivatives(u_temp, v_temp, k6_u, k6_v, k, W, U_s);

        // Compute k7
        #pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                u_temp[i][j] = u[i][j] + dt * (35.0 / 384.0 * k1_u[i][j] + 500.0 / 1113.0 * k3_u[i][j] + 125.0 / 192.0 * k4_u[i][j] - 2187.0 / 6784.0 * k5_u[i][j] + 11.0 / 84.0 * k6_u[i][j]);
                v_temp[i][j] = v[i][j] + dt * (35.0 / 384.0 * k1_v[i][j] + 500.0 / 1113.0 * k3_v[i][j] + 125.0 / 192.0 * k4_v[i][j] - 2187.0 / 6784.0 * k5_v[i][j] + 11.0 / 84.0 * k6_v[i][j]);
            }
        }
        compute_derivatives(u_temp, v_temp, k7_u, k7_v, k, W, U_s);

        // Update u and v
        #pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                u[i][j] += dt * (35.0 / 384.0 * k1_u[i][j] + 0.0 * k2_u[i][j] + 500.0 / 1113.0 * k3_u[i][j] + 125.0 / 192.0 * k4_u[i][j] -2187.0 / 6784.0 * k5_u[i][j] + 11.0 / 84.0 * k6_u[i][j]);
                v[i][j] += dt * (35.0 / 384.0 * k1_v[i][j] + 0.0 * k2_v[i][j] + 500.0 / 1113.0 * k3_v[i][j] + 125.0 / 192.0 * k4_v[i][j] -2187.0 / 6784.0 * k5_v[i][j] + 11.0 / 84.0 * k6_v[i][j]);
            }
        }

        // Compute variance with respect to the center point
        double variance = 0.0;
        #pragma omp parallel for reduction(+:variance)
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                double diff = u[i][j] - u[N/2][N/2];
                variance += diff * diff;
            }
        }
        variance /= (N * N);

        if(variance < 1e-20) {
            break;
        }

        display_loading_bar(step, num_steps, variance);
    }

    std::ostringstream filename_u;
    filename_u << "output/frame_u_phi" << phi << "_gamma" << gamma << "_" << std::setw(6) << std::setfill('0') << num_steps << ".csv";
    save_to_csv(u, filename_u.str());

    std::ostringstream filename_v;
    filename_v << "output/frame_v_phi" << phi << "_gamma" << gamma << "_" << std::setw(6) << std::setfill('0') << num_steps << ".csv";
    save_to_csv(v, filename_v.str());

    const auto end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s" << std::endl;
    return 0;
}