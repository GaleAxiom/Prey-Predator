#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <omp.h>

// Function to solve the tridiagonal system using the Thomas algorithm
void thomas_algorithm(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, std::vector<double>& d, std::vector<double>& x) {
    int n = b.size();
    std::vector<double> c_prime(n, 0.0);
    std::vector<double> d_prime(n, 0.0);

    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        double m = 1.0 / (b[i] - a[i] * c_prime[i - 1]);
        c_prime[i] = c[i] * m;
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) * m;
    }

    x[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }
}

void save_to_csv(const std::vector<std::vector<double> >& data, const std::string& filename) {
    std::ofstream file(filename);
    file << std::fixed << std::setprecision(12);  // Set precision to 12 decimal places
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
    // Parameters
    const int nx = 100; // Number of spatial points in x direction
    const int ny = 100; // Number of spatial points in y direction
    const int nt = 1000000; // Number of time steps
    const double Lx = 1.0; // Length of the rod in x direction
    const double Ly = 1.0; // Length of the rod in y direction
    const double T = 100.0; // Total time

    // Parameters for the equations
    const double rho_d = 0.3; // Thermal diffusivity
    const double rho_r = 3.33; // Reaction rate
    const double U_s = 0.4; // Saturation concentration
    const double phi = 0.75; // Ratio of the reaction rates
    const double gamma = 0.4; // Ratio of the reaction rates

    const auto k = (phi - 2.) / (phi - 1.) * U_s;
    const auto W = gamma * U_s / (1 - gamma);

    const double dx = Lx / (nx - 1);
    const double dy = Ly / (ny - 1);
    const double dt = T / nt;
    const double rx_u = rho_d * dt / (12 * dx * dx);
    const double ry_u = rho_d * dt / (12 * dy * dy);
    const double rx_v = dt / (12 * dx * dx);
    const double ry_v = dt / (12 * dy * dy);

    // Initial condition
    std::vector<std::vector<double>> u(nx, std::vector<double>(ny, 1.0));
    std::vector<std::vector<double>> v(nx, std::vector<double>(ny, 1.0));
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            u[i][j] += ((double) rand() / RAND_MAX) * 1e-5;
            v[i][j] += ((double) rand() / RAND_MAX) * 1e-5;
        }
    }

    // Time-stepping loop
    for (int t = 0; t < nt; ++t) {
        // Iterative solver for non-linear equations
        const int max_iter = 10000;
        const double tol = 1e-4;
        bool converged = false;
        double max_error = 0.0;

        if (t % 1000 == 0) {
            // Save the results to a CSV file
            std::cout << "Saving frame " << t << std::endl;
            std::ostringstream filename_u;
            filename_u << "output/frame_u_phi" << phi << "_gamma" << gamma << "_" << std::setw(6) << std::setfill('0') << t << ".csv";
            save_to_csv(u, filename_u.str());
        }

        for (int iter = 0; iter < max_iter; ++iter) {
            max_error = 0.0;

            // Solve temperature in x-direction
            for (int j = 0; j < ny; ++j) {
                std::vector<double> a(nx - 1, -rx_u);
                std::vector<double> b(nx, 1.0 + 30 * rx_u);
                std::vector<double> c(nx - 1, -rx_u);
                std::vector<double> d(nx, 0.0);
                for (int i = 0; i < nx; ++i) {
                    int ip1 = (i + 1) % nx;
                    int ip2 = (i + 2) % nx;
                    int im1 = (i - 1 + nx) % nx;
                    int im2 = (i - 2 + nx) % nx;
                    d[i] = (1.0 - 30 * rx_u) * u[i][j] + rx_u * (16 * (u[ip1][j] + u[im1][j]) - (u[ip2][j] + u[im2][j]))
                            + dt * rho_r * ((k - U_s * u[i][j]) / (k - U_s) * u[i][j] - (W + U_s) / (W + U_s * u[i][j]) * u[i][j] * v[i][j]);
                }

                std::vector<double> temp_new(nx, 0.0);
                thomas_algorithm(a, b, c, d, temp_new);

                for (int i = 0; i < nx; ++i) {
                    double error = std::abs(temp_new[i] - u[i][j]);
                    max_error = std::max(max_error, error);
                    u[i][j] = temp_new[i];
                }
            }

            // Solve temperature in y-direction
            for (int i = 0; i < nx; ++i) {
                std::vector<double> a(ny - 1, -ry_u);
                std::vector<double> b(ny, 1.0 + 30 * ry_u);
                std::vector<double> c(ny - 1, -ry_u);
                std::vector<double> d(ny, 0.0);
                for (int j = 0; j < ny; ++j) {
                    int jp1 = (j + 1) % ny;
                    int jp2 = (j + 2) % ny;
                    int jm1 = (j - 1 + ny) % ny;
                    int jm2 = (j - 2 + ny) % ny;
                    d[j] = (1.0 - 30 * ry_u) * u[i][j] + ry_u * (16 * (u[i][jp1] + u[i][jm1]) - (u[i][jp2] + u[i][jm2]))
                            + dt * rho_r * ((k - U_s * u[i][j]) / (k - U_s) * u[i][j] - (W + U_s) / (W + U_s * u[i][j]) * u[i][j] * v[i][j]);
                }

                std::vector<double> temp_new(ny, 0.0);
                thomas_algorithm(a, b, c, d, temp_new);

                for (int j = 0; j < ny; ++j) {
                    double error = std::abs(temp_new[j] - u[i][j]);
                    max_error = std::max(max_error, error);
                    u[i][j] = temp_new[j];
                }
            }

            // Solve concentration in x-direction
            for (int j = 0; j < ny; ++j) {
                std::vector<double> a(nx - 1, -rx_v);
                std::vector<double> b(nx, 1.0 + 30 * rx_v);
                std::vector<double> c(nx - 1, -rx_v);
                std::vector<double> d(nx, 0.0);
                for (int i = 0; i < nx; ++i) {
                    int ip1 = (i + 1) % nx;
                    int ip2 = (i + 2) % nx;
                    int im1 = (i - 1 + nx) % nx;
                    int im2 = (i - 2 + nx) % nx;
                    d[i] = (1.0 - 30 * rx_v) * v[i][j] + rx_v * (16 * (v[ip1][j] + v[im1][j]) - (v[ip2][j] + v[im2][j]))
                            + dt * (W + U_s) / (W + U_s * u[i][j]) * u[i][j] * v[i][j] - v[i][j] * v[i][j];
                }

                std::vector<double> conc_new(nx, 0.0);
                thomas_algorithm(a, b, c, d, conc_new);

                for (int i = 0; i < nx; ++i) {
                    double error = std::abs(conc_new[i] - v[i][j]);
                    max_error = std::max(max_error, error);
                    v[i][j] = conc_new[i];
                }
            }

            // Solve concentration in y-direction
            for (int i = 0; i < nx; ++i) {
                std::vector<double> a(ny - 1, -ry_v);
                std::vector<double> b(ny, 1.0 + 30 * ry_v);
                std::vector<double> c(ny - 1, -ry_v);
                std::vector<double> d(ny, 0.0);
                for (int j = 0; j < ny; ++j) {
                    int jp1 = (j + 1) % ny;
                    int jp2 = (j + 2) % ny;
                    int jm1 = (j - 1 + ny) % ny;
                    int jm2 = (j - 2 + ny) % ny;
                    d[j] = (1.0 - 30 * ry_v) * v[i][j] + ry_v * (16 * (v[i][jp1] + v[i][jm1]) - (v[i][jp2] + v[i][jm2]))
                            + dt * (W + U_s) / (W + U_s * u[i][j]) * u[i][j] * v[i][j] - v[i][j] * v[i][j];
                }

                std::vector<double> conc_new(ny, 0.0);
                thomas_algorithm(a, b, c, d, conc_new);

                for (int j = 0; j < ny; ++j) {
                    double error = std::abs(conc_new[j] - v[i][j]);
                    max_error = std::max(max_error, error);
                    v[i][j] = conc_new[j];
                }
            }

            std::cout << "Time :" << t << ", Iteration: " << iter << ", Error: " << max_error << std::endl;
            if (max_error < tol) {
                converged = true;
                break;
            }
        }

        if (!converged) {
            std::cerr << "Warning: Solution did not converge within the maximum number of iterations." << std::endl;
        }
    }

    return 0;
}